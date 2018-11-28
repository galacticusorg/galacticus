!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module which implements calculations of stellar population luminosities in the AB magnitude system.

module Stellar_Population_Luminosities
  !% Implements calculations of stellar population luminosities in the AB magnitude system.
  use, intrinsic :: ISO_C_Binding
  use            :: FGSL                                  , only : fgsl_interp_accel, fgsl_function, fgsl_integration_workspace, FGSL_Integ_Gauss15
  use            :: Abundances_Structure
  use            :: ISO_Varying_String
  use            :: Locks
  use            :: Stellar_Population_Spectra_Postprocess
  implicit none
  private
  public :: Stellar_Population_Luminosity, Stellar_Population_Luminosity_Track

  type luminosityTable
     !% Structure for holding tables of simple stellar population luminosities.
     integer                                                            :: agesCount                          , metallicitiesCount
     integer         (c_size_t)                                         :: isTabulatedMaximum         =0
     logical                            , allocatable, dimension(:)     :: isTabulated
     double precision                   , allocatable, dimension(:)     :: age                                , metallicity
     double precision                   , allocatable, dimension(:,:,:) :: luminosity
     ! Interpolation structures.
     logical                                                            :: resetAge                   =.true. , resetMetallicity                   =.true.
     type            (fgsl_interp_accel)                                :: interpolationAcceleratorAge        , interpolationAcceleratorMetallicity
  end type luminosityTable

  ! Array of simple stellar population luminosity tables.
  type            (luminosityTable), allocatable, dimension(:) :: luminosityTables

  ! Module global variables used in integrations.
  double precision                                                          :: ageTabulate                                  , redshiftTabulate
  integer                                                                   :: filterIndexTabulate
  integer         (c_size_t                                  )              :: populationIDTabulate
  type            (abundances                                )              :: abundancesTabulate
  class           (stellarPopulationSpectraPostprocessorClass), allocatable :: stellarPopulationSpectraPostprocessorTabulate
  !$omp threadprivate(ageTabulate,redshiftTabulate,abundancesTabulate,filterIndexTabulate,populationIDTabulate,stellarPopulationSpectraPostprocessorTabulate)

  ! Flag indicating if this module has been initialized yet.
  logical                                                      :: moduleInitialized                                      =.false.

  ! Tolerance used in integrations.
  double precision                                             :: stellarPopulationLuminosityIntegrationToleranceRelative
  logical                                                      :: stellarPopulationLuminosityIntegrationToleranceDegrade

  ! Option controlling writing of luminosities to file.
  logical                                                      :: stellarPopulationLuminosityStoreToFile
  type            (varying_string )                            :: stellarPopulationLuminosityStoreDirectory
  
  ! Option controlling behavior when maximum age of stellar populations is exceeded.
  logical                                                      :: stellarPopulationLuminosityMaximumAgeExceededIsFatal

  ! Read/write lock used to control access to luminosity tables.
  type            (ompReadWriteLock)                           :: luminosityTableLock
 
contains

  subroutine Stellar_Population_Luminosity_Tabulate(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    !% Tabulate stellar population luminosity in the given filters.
    use, intrinsic :: ISO_C_Binding
    use            :: Numerical_Constants_Astronomical
    use            :: File_Utilities
    use            :: IO_HDF5
    use            :: String_Handling
    use            :: Input_Parameters
    use            :: Galacticus_Paths
    use            :: Galacticus_Display
    use            :: Galacticus_Error
    use            :: Instruments_Filters
    use            :: Numerical_Integration
    use            :: Memory_Management
    use            :: Stellar_Populations
    use            :: Stellar_Population_Spectra
    use            :: Input_Parameters
    use            :: String_Handling
    implicit none
    integer                                                    , intent(in   ), dimension(:    ) :: filterIndex                           , luminosityIndex
    double precision                                           , intent(in   ), dimension(:    ) :: redshift
    type            (stellarPopulationSpectraPostprocessorList), intent(in   ), dimension(:    ) :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                   ), intent(inout)                   :: stellarPopulation_
    type            (luminosityTable                          ), allocatable  , dimension(:    ) :: luminosityTablesTemporary
    double precision                                           , allocatable  , dimension(:,:,:) :: luminosityTemporary
    logical                                                    , allocatable  , dimension(:    ) :: isTabulatedTemporary
    double precision                                                          , dimension(2    ) :: wavelengthRange
    type            (lockDescriptor                           )                                  :: lockFileDescriptor
    class           (stellarPopulationSpectraClass            ), pointer                         :: stellarPopulationSpectra_
    integer         (c_size_t                                 )                                  :: iAge                                  , iLuminosity           , &
         &                                                                                          iMetallicity                          , jLuminosity           , &
         &                                                                                          populationID
    integer                                                                                      :: loopCountMaximum                      , loopCount             , &
         &                                                                                          errorStatus                           , luminosityIndexMaximum
    logical                                                                                      :: computeTable                          , calculateLuminosity
    double precision                                                                             :: toleranceRelative                     , normalization
    type            (fgsl_function                            )                                  :: integrandFunction
    type            (fgsl_integration_workspace               )                                  :: integrationWorkspace
    type            (varying_string                           )                                  :: message                               , luminositiesFileName  , &
         &                                                                                          descriptorString
    character       (len=16                                   )                                  :: datasetName                           , redshiftLabel         , &
         &                                                                                          label
    type            (hdf5Object                               )                                  :: luminositiesFile
    type            (inputParameters                          )                                  :: descriptor

    ! Determine if we have created space for this population yet.
    if (.not.moduleInitialized) then
       !$omp critical (Luminosity_Tables_Initialize)
       if (.not.moduleInitialized) then
          ! Read the parameter controlling integration tolerance.
          !# <inputParameter>
          !#   <name>stellarPopulationLuminosityIntegrationToleranceRelative</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>4.0d-3</defaultValue>
          !#   <description>The relative tolerance used when integrating the flux of stellar populations through filters.</description>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>stellarPopulationLuminosityIntegrationToleranceDegrade</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>If {\normalfont \ttfamily true}, automatically degrade the relative tolerance used when integrating the flux of stellar populations through filters to ensure convergence.</description>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          ! Read parameters controlling storing luminosities to file.
          !# <inputParameter>
          !#   <name>stellarPopulationLuminosityStoreToFile</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.true.</defaultValue>
          !#   <description>Specifies whether or not stellar populations luminosities (integrated under a filter) should be stored to file for rapid reuse.</description>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>stellarPopulationLuminosityStoreDirectory</name>
          !#   <defaultValue>galacticusPath(pathTypeDataDynamic)//'stellarPopulations'</defaultValue>
          !#   <attachedTo>module</attachedTo>
          !#   <description>
          !#    Specifies the directory to which stellar populations luminosities (integrated under a filter) should be stored to file for rapid reuse.
          !#   </description>
          !#   <type>string</type>
          !#   <cardinality>1</cardinality>
          !# </inputParameter>
          ! Read the parameter controlling behavior if maximum age of stellar populations are exceeded.
          !# <inputParameter>
          !#   <name>stellarPopulationLuminosityMaximumAgeExceededIsFatal</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.true.</defaultValue>
          !#   <description>Specifies whether or not exceeding the maximum available age of the stellar population is fatal.</description>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Initialize the luminosity table lock.
          luminosityTableLock=ompReadWriteLock()
          ! Flag that this module is now initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical (Luminosity_Tables_Initialize)
    end if

    ! Obtain a read lock on the luminosity tables.
    call luminosityTableLock%setRead()
    
    ! Allocate table storage. First test if the tables must be resized (or allocated).
    populationID=stellarPopulation_%uniqueID()
    if (.not.allocated(luminosityTables).or.size(luminosityTables) < populationID) then
       ! A resize or allocation is required. Obtain a write lock on the tables and then re-test if resizing is required (as
       ! another thread could potentially have alreaedy done this while we waited for the write lock).
       call luminosityTableLock%setWrite(haveReadLock=.true.)
       if (allocated(luminosityTables)) then
          if (size(luminosityTables) < populationID) then
             call Move_Alloc(luminosityTables,luminosityTablesTemporary)
             allocate(luminosityTables(populationID))
             luminosityTables(1:size(luminosityTablesTemporary))=luminosityTablesTemporary
             luminosityTables(size(luminosityTablesTemporary)+1:populationID)%isTabulatedMaximum=0
             deallocate(luminosityTablesTemporary)
             call Memory_Usage_Record(sizeof(luminosityTables(1)),blockCount=0)
          end if
       else
          allocate(luminosityTables(populationID))
          luminosityTables%isTabulatedMaximum=0
          call Memory_Usage_Record(sizeof(luminosityTables))
       end if
       call luminosityTableLock%unsetWrite(haveReadLock=.true.)
    end if
    ! Determine if we have tabulated luminosities for this luminosityIndex in this population yet.
    luminosityIndexMaximum=maxval(luminosityIndex)    
    if (.not.allocated(luminosityTables(populationID)%isTabulated) .or. luminosityTables(populationID)%isTabulatedMaximum < luminosityIndexMaximum) then
       ! Tabulation is required. Obtain a write lock on the luminosity tables and then retest if tabulation is required (as
       ! another thread may have tabulated while we waited for the write lock).
       call luminosityTableLock%setWrite(haveReadLock=.true.)
       luminosityIndexMaximum=maxval(luminosityIndex)    
       if (.not.allocated(luminosityTables(populationID)%isTabulated) .or. luminosityTables(populationID)%isTabulatedMaximum < luminosityIndexMaximum) then
          call File_Lock_Initialize(lockFileDescriptor)
          stellarPopulationSpectra_ => stellarPopulation_%spectra()
          luminosityIndexMaximum    =  maxval(luminosityIndex)
          do iLuminosity=1,size(luminosityIndex)
             if (allocated(luminosityTables(populationID)%isTabulated)) then
                if (size(luminosityTables(populationID)%isTabulated) >= luminosityIndex(iLuminosity)) then
                   computeTable=.not.luminosityTables(populationID)%isTabulated(luminosityIndex(iLuminosity))
                else
                   call Move_Alloc (luminosityTables(populationID)%isTabulated,isTabulatedTemporary)
                   call Move_Alloc (luminosityTables(populationID)%luminosity ,luminosityTemporary )
                   call allocateArray(luminosityTables(populationID)%isTabulated,[luminosityIndexMaximum])
                   call allocateArray(luminosityTables(populationID)%luminosity ,[luminosityIndexMaximum&
                        &,luminosityTables(populationID)%agesCount,luminosityTables(populationID)%metallicitiesCount])
                   luminosityTables(populationID)%isTabulated(1:size(isTabulatedTemporary)    )=isTabulatedTemporary
                   luminosityTables(populationID)%isTabulated(  size(isTabulatedTemporary)+1:luminosityIndexMaximum)=.false.
                   luminosityTables(populationID)%luminosity (1:size(isTabulatedTemporary),:,:)=luminosityTemporary
                   call deallocateArray(isTabulatedTemporary)
                   call deallocateArray(luminosityTemporary)
                   computeTable=.true.
                end if
             else
                call allocateArray(luminosityTables(populationID)%isTabulated,[luminosityIndexMaximum])
                luminosityTables(populationID)%isTabulated=.false.
                ! Since we have not yet tabulated any luminosities yet for this population, we need to get a list of suitable
                ! metallicities and ages at which to tabulate.
                call stellarPopulationSpectra_%tabulation(luminosityTables(populationID)%agesCount &
                     &,luminosityTables(populationID)%metallicitiesCount,luminosityTables(populationID)%age&
                     &,luminosityTables(populationID)%metallicity)
                where (luminosityTables(populationID)%metallicity > 0.0d0)
                   luminosityTables(populationID)%metallicity=log10(luminosityTables(populationID)%metallicity/metallicitySolar)
                elsewhere
                   luminosityTables(populationID)%metallicity=logMetallicityZero
                end where
                call allocateArray(luminosityTables(populationID)%luminosity,[luminosityIndexMaximum&
                     &,luminosityTables(populationID)%agesCount ,luminosityTables(populationID)%metallicitiesCount])
                computeTable=.true.
             end if

             ! If we haven't, do so now.
             if (computeTable) then

                ! Determine if we can read the required luminosity from file.
                calculateLuminosity=.true.
                if (stellarPopulationLuminosityStoreToFile) then
                   ! Construct name of the file to which this would be stored.
                   luminositiesFileName=stellarPopulationLuminosityStoreDirectory                                                                                                                // &
                        &               "/stellarLuminosities::filter:"                                                                                                                          // &
                        &               Filter_Name                                                                                                (                    filterIndex(iLuminosity))// &
                        &               "::postprocessor:"                                                                                                                                       // &
                        &               stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_%hashedDescriptor(includeSourceDigest=.true.                  )// &
                        &               "::population:"                                                                                                                                          // &
                        &               stellarPopulation_                                                                        %hashedDescriptor(includeSourceDigest=.true.                  )// &
                        &               ".hdf5"
                   if (File_Exists(luminositiesFileName)) then
                      ! Construct the dataset name.
                      write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                      datasetName="redshift"//adjustl(trim(redshiftLabel))
                      ! Open the file and check for the required dataset.
                      call hdf5Access%set()
                      call File_Lock(char(luminositiesFileName),lockFileDescriptor,lockIsShared=.true.)
                      call luminositiesFile%openFile(char(luminositiesFileName),readOnly=.true.)
                      if (luminositiesFile%hasDataset(trim(datasetName))) then
                         ! Read the dataset.
                         call luminositiesFile%readDatasetStatic(trim(datasetName),luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:))
                         ! We do not need to calculate this luminosity.
                         calculateLuminosity=.false.
                      end if
                      call luminositiesFile%close()
                      call File_Unlock(lockFileDescriptor)
                      call hdf5Access%unset()
                   end if
                end if

                ! Compute the luminosity if necessary.
                if (calculateLuminosity) then
                   ! Display a message and counter.
                   message=var_str('Tabulating stellar luminosities for stellar population #')//populationID//', luminosity '
                   write (redshiftLabel,'(f6.3)') redshift(iLuminosity)
                   message=message                                                                                     // &
                        &  Filter_Name                                          (filterIndex             (iLuminosity))// &
                        &  ":z"                                                                                        // &
                        &  trim(adjustl(redshiftLabel))                                                                // &
                        &  " "                                                                                         // &
                        &                                                                                 iLuminosity  // &
                        &  " of "                                                                                      // &
                        &  size(luminosityIndex)
                   call Galacticus_Display_Indent (message,verbosityWorking)
                   call Galacticus_Display_Counter(0,.true.,verbosityWorking)             
                   ! Get wavelength extent of the filter.
                   wavelengthRange=Filter_Extent(filterIndex(iLuminosity))
                   ! Integrate over the wavelength range.
                   filterIndexTabulate = filterIndex     (iLuminosity)
                   redshiftTabulate    = redshift        (iLuminosity)
                   populationIDTabulate= populationID
                   loopCountMaximum    =+luminosityTables(populationID)%metallicitiesCount &
                        &               *luminosityTables(populationID)%agesCount
                   loopCount           = 0
                   !$omp parallel private(iAge,iMetallicity,integrandFunction,integrationWorkspace,toleranceRelative,errorStatus) copyin(filterIndexTabulate,redshiftTabulate,populationIDTabulate)
                   allocate(stellarPopulationSpectraPostprocessorTabulate,mold=stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_)
                   call stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_%deepCopy(stellarPopulationSpectraPostprocessorTabulate)
                   !$omp do
                   do iAge=1,luminosityTables(populationID)%agesCount
                      ageTabulate=luminosityTables(populationID)%age(iAge)
                      do iMetallicity=1,luminosityTables(populationID)%metallicitiesCount
                         ! Update the counter.
                         !$omp atomic
                         loopCount=loopCount+1
                         call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountMaximum)),.false.,verbosityWorking)
                         call abundancesTabulate%metallicitySet(luminosityTables(populationID)%metallicity(iMetallicity) &
                              &,metallicityType=metallicityTypeLogarithmicByMassSolar)
                         toleranceRelative=stellarPopulationLuminosityIntegrationToleranceRelative
                         errorStatus      =errorStatusFail
                         do while (errorStatus /= errorStatusSuccess)
                            luminosityTables(populationID)%luminosity(                              &
                                 &                                luminosityIndex(iLuminosity), &
                                 &                                iAge                        , &
                                 &                                iMetallicity                  &
                                 &                               )                              &
                                 & =Integrate(                                                  &
                                 &            wavelengthRange(1)                              , &
                                 &            wavelengthRange(2)                              , &
                                 &            Filter_Luminosity_Integrand                     , &
                                 &            integrandFunction                               , &
                                 &            integrationWorkspace                            , &
                                 &            toleranceAbsolute          =0.0d0               , &
                                 &            toleranceRelative          =toleranceRelative   , &
                                 &            integrationRule            =FGSL_Integ_Gauss15  , &
                                 &            maxIntervals               =10000               , &
                                 &            errorStatus                =errorStatus           &
                                 &           )
                            call Integrate_Done(integrandFunction,integrationWorkspace)
                            if (errorStatus /= errorStatusSuccess) then
                               if (stellarPopulationLuminosityIntegrationToleranceDegrade.and.toleranceRelative < 1.0d0) then
                                  toleranceRelative=2.0d0*toleranceRelative
                                  write (label,'(e9.3)') 2.0d0*stellarPopulationLuminosityIntegrationToleranceRelative
                                  message=         "WARNING: increasing relative tolerance for stellar population luminosities to"          //char(10)
                                  message=message//trim(adjustl(label))//" and retrying integral"//{introspection:location}
                                  call Galacticus_Warn(message)
                               else if (stellarPopulationLuminosityIntegrationToleranceDegrade) then
                                  message="integration of stellar populations failed"
                                  call Galacticus_Error_Report(message//{introspection:location})
                               else
                                  write (label,'(e9.3)') 2.0d0*stellarPopulationLuminosityIntegrationToleranceRelative
                                  message=         "integration of stellar populations failed"                                              //char(10)
                                  message=message//"HELP: consider increasing the [stellarPopulationLuminosityIntegrationToleranceRelative]"//char(10)
                                  message=message//"      parameter to "//trim(adjustl(label))//" to reduce the integration tolerance"      //char(10)
                                  message=message//"      required if you can accept this lower accuracy."
                                  call Galacticus_Error_Report(message//{introspection:location})
                               end if
                            end if
                         end do
                      end do
                   end do
                   !$omp end do
                   deallocate(stellarPopulationSpectraPostprocessorTabulate)
                   !$omp end parallel
                   ! Clear the counter and write a completion message.
                   call Galacticus_Display_Counter_Clear(           verbosityWorking)
                   call Galacticus_Display_Unindent     ('finished',verbosityWorking)
                   ! Get the normalization by integrating a zeroth magnitude (AB) source through the filter.
                   normalization=Integrate(wavelengthRange(1),wavelengthRange(2),Filter_Luminosity_Integrand_AB &
                        &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative&
                        &=stellarPopulationLuminosityIntegrationToleranceRelative)
                   call Integrate_Done(integrandFunction,integrationWorkspace)
                   ! Normalize the luminosity.
                   luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:) &
                        &=luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:)/normalization
                   ! Store the luminosities to file.
                   if (stellarPopulationLuminosityStoreToFile) then
                      ! Construct the dataset name.
                      write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                      datasetName="redshift"//adjustl(trim(redshiftLabel))
                      ! Open the file.
                      descriptor=inputParameters()
                      call stellarPopulation_%descriptor(descriptor,includeMethod=.true.)
                      descriptorString=descriptor%serializeToString()
                      call hdf5Access%set()
                      call File_Lock(char(luminositiesFileName),lockFileDescriptor,lockIsShared=.false.)
                      call luminositiesFile%openFile      (char(luminositiesFileName)             )
                      if (.not.luminositiesFile%hasAttribute('parameters')) call luminositiesFile%writeAttribute(char(descriptorString    ),'parameters')
                      ! Write the dataset.
                      if (.not.luminositiesFile%hasDataset(trim(datasetName))) &
                           & call luminositiesFile%writeDataset(luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:),datasetName=trim(datasetName),commentText="Tabulated luminosities at redshift z="//adjustl(trim(redshiftLabel)))
                      ! Close the file.
                      call luminositiesFile%close()
                      call File_Unlock(lockFileDescriptor)
                      call hdf5Access%unset()
                   end if
                end if
                ! Flag that calculations have been performed for this filter.
                luminosityTables(populationID)%isTabulated(luminosityIndex(iLuminosity))=.true.
                if (luminosityIndex(iLuminosity) > luminosityTables(populationID)%isTabulatedMaximum) then
                   jLuminosity=luminosityTables(populationID)%isTabulatedMaximum
                   do while (jLuminosity < size(luminosityTables(populationID)%isTabulated) .and. luminosityTables(populationID)%isTabulated(min(jLuminosity+1,size(luminosityTables(populationID)%isTabulated))))
                      jLuminosity=jLuminosity+1
                   end do
                   luminosityTables(populationID)%isTabulatedMaximum=jLuminosity
                end if
             end if
          end do
       end if
       call luminosityTableLock%unsetWrite(haveReadLock=.true.)
    end if
    ! Release our read lock on the luminosity tables.
    call luminosityTableLock%unsetRead()
    return

  contains

    double precision function Filter_Luminosity_Integrand(wavelength)
      !% Integrand for the luminosity through a given filter.
      use Instruments_Filters
      implicit none
      double precision, intent(in   ) :: wavelength
      double precision                :: wavelengthRedshifted
      
      ! If this luminosity is for a redshifted spectrum, then we shift wavelength at which we sample the stellar population spectrum
      ! to be a factor of (1+z) smaller. We therefore integrate over the stellar SED at shorter wavelengths, since these will be
      ! shifted into the filter by z=0. Factor of 1/wavelength appears since we want to integrate F_nu (dnu / nu) and dnu =
      ! -c/lambda^2 dlambda. Note that we follow the convention of Hogg et al. (2002) and assume that the filter response gives the
      ! fraction of incident photons received by the detector at a given wavelength, multiplied by the relative photon response
      ! (which will be 1 for a photon-counting detector such as a CCD, or proportional to the photon energy for a
      ! bolometer/calorimeter type detector).
      wavelengthRedshifted=wavelength/(1.0d0+redshiftTabulate)
      Filter_Luminosity_Integrand=+Filter_Response                                         (filterIndexTabulate             ,wavelength          ) &
           &                      *stellarPopulationSpectra_                    %luminosity(abundancesTabulate  ,ageTabulate,wavelengthRedshifted) &
           &                      *stellarPopulationSpectraPostprocessorTabulate%multiplier(wavelengthRedshifted,ageTabulate,redshiftTabulate    ) &
           &                      /                                                                                          wavelength
      return
    end function Filter_Luminosity_Integrand
    
    double precision function Filter_Luminosity_Integrand_AB(wavelength)
      !% Integrand for the luminosity of a zeroth magnitude (AB) source through a given filter.
      use Instruments_Filters
      use Numerical_Constants_Astronomical
      implicit none
      double precision, intent(in   ) :: wavelength
      ! Luminosity of a zeroth magintude (AB) source in Solar luminosities per Hz.
      double precision, parameter     :: luminosityZeroPointABSolar=luminosityZeroPointAB/luminositySolar
      
      Filter_Luminosity_Integrand_AB=Filter_Response(filterIndexTabulate,wavelength)*luminosityZeroPointABSolar/wavelength
      return
    end function Filter_Luminosity_Integrand_AB
    
  end subroutine Stellar_Population_Luminosity_Tabulate
  
  function Stellar_Population_Luminosity(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,abundancesStellar,age,redshift)
    !% Returns the luminosity for a $1 M_\odot$ simple {\normalfont \ttfamily stellarPopulation\_} of given {\normalfont \ttfamily
    !% abundances} and {\normalfont \ttfamily age} and observed through the filter specified by {\normalfont \ttfamily
    !% filterIndex}.
    use, intrinsic :: ISO_C_Binding
    use            :: Galacticus_Error
    use            :: Numerical_Interpolation
    use            :: Stellar_Populations
    implicit none
    integer                                                    , intent(in   ), dimension(:)                     :: filterIndex                           , luminosityIndex
    double precision                                           , intent(in   ), dimension(:)                     :: age                                   , redshift
    type            (stellarPopulationSpectraPostprocessorList), intent(in   ), dimension(:)                     :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                   ), intent(inout)                                   :: stellarPopulation_
    type            (abundances                               ), intent(in   )                                   :: abundancesStellar
    double precision                                                          , dimension(size(luminosityIndex)) :: Stellar_Population_Luminosity
    double precision                                                          , dimension(0:1                  ) :: hAge                                  , hMetallicity
    integer         (c_size_t                                 )                                                  :: iAge                                  , iLuminosity             , &
         &                                                                                                          iMetallicity                          , jAge                    , &
         &                                                                                                          jMetallicity                          , populationID
    double precision                                                                                             :: ageLast                               , metallicity

    ! Tabulate the luminosities.
    call Stellar_Population_Luminosity_Tabulate(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    ! Obtain a read lock on the luminosity tables.
    call luminosityTableLock%setRead()
    ! Get interpolation in metallicity.
    populationID=stellarPopulation_%uniqueID()
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < luminosityTables(populationID)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > luminosityTables(populationID)%metallicity(luminosityTables(populationID)%metallicitiesCount)) then
       iMetallicity=luminosityTables(populationID)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       iMetallicity=Interpolate_Locate(luminosityTables(populationID)%metallicity &
            &,luminosityTables(populationID)%interpolationAcceleratorMetallicity,metallicity &
            &,luminosityTables(populationID)%resetMetallicity)
       hMetallicity=Interpolate_Linear_Generate_Factors(luminosityTables(populationID)%metallicity ,iMetallicity,metallicity)
    end if
    ! Do the interpolation.
    Stellar_Population_Luminosity(:)= 0.0d0
    ageLast                         =-1.0d0
    do iLuminosity=1,size(luminosityIndex)
       ! Only compute luminosities for entries with positive age (negative age implies that the luminosity required is for a
       ! population observed prior to the formation of this population).
       if (age(iLuminosity) >= 0.0d0) then
          ! Get interpolation in age if the age for this luminosity differs from the previous one.
          if (iLuminosity == 1 .or. age(iLuminosity) /= ageLast) then
             ! Check for out of range age.
             if (age(iLuminosity) > luminosityTables(populationID)%age(luminosityTables(populationID)%agesCount)) then
                if (stellarPopulationLuminosityMaximumAgeExceededIsFatal) then
                   call Galacticus_Error_Report('age exceeds the maximum tabulated'//{introspection:location})
                else
                   iAge=luminosityTables(populationID)%agesCount-1
                   hAge=[0.0d0,1.0d0]
                end if
             else
                iAge=Interpolate_Locate(luminosityTables(populationID)%age &
                     &,luminosityTables(populationID)%interpolationAcceleratorAge,age(iLuminosity),luminosityTables(populationID)%resetAge)
                hAge=Interpolate_Linear_Generate_Factors(luminosityTables(populationID)%age,iAge&
                     &,age(iLuminosity))
             end if
             ageLast=age(iLuminosity)
          end if
          do jAge=0,1
             do jMetallicity=0,1
                Stellar_Population_Luminosity(iLuminosity)=Stellar_Population_Luminosity(iLuminosity)&
                     &+luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),iAge +jAge,iMetallicity+jMetallicity)&
                     &*hAge(jAge)*hMetallicity(jMetallicity)
             end do
          end do
       end if
    end do
    call luminosityTableLock%unsetRead()
    ! Prevent interpolation from returning negative fluxes.
    Stellar_Population_Luminosity=max(Stellar_Population_Luminosity,0.0d0)
    return
  end function Stellar_Population_Luminosity
  
  subroutine Stellar_Population_Luminosity_Track(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,abundancesStellar,redshift,ages,luminosities)
    !% Returns the luminosity for a $1 M_\odot$ simple stellar population of given {\normalfont \ttfamily abundances} drawn from
    !% the given {\normalfont \ttfamily stellarPopulation} and observed through the filter specified by {\normalfont \ttfamily
    !% filterIndex}, for all available ages.
    use, intrinsic :: ISO_C_Binding
    use            :: Galacticus_Error
    use            :: Memory_Management
    use            :: Numerical_Interpolation
    use            :: Stellar_Populations
    implicit none
    integer                                                    , intent(in   ), dimension(:    )              :: filterIndex                           , luminosityIndex
    double precision                                           , intent(in   ), dimension(:    )              :: redshift
    type            (stellarPopulationSpectraPostprocessorList), intent(in   ), dimension(:)                  :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                   ), intent(inout)                                :: stellarPopulation_
    type            (abundances                               ), intent(in   )                                :: abundancesStellar
    double precision                                           , intent(  out), dimension(:    ), allocatable :: ages
    double precision                                           , intent(  out), dimension(:  ,:), allocatable :: luminosities
    double precision                                                          , dimension(0:1  )              :: hMetallicity
    integer         (c_size_t                                 )                                               :: iLuminosity                           , iMetallicity   , &
         &                                                                                                       jMetallicity                          , populationID
    double precision                                                                                          :: metallicity

    ! Obtain a read lock on the luminosity tables.
    call luminosityTableLock%setRead()
    ! Tabulate the luminosities.
    call Stellar_Population_Luminosity_Tabulate(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    ! Get interpolation in metallicity.
    populationID=stellarPopulation_%uniqueID()
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < luminosityTables(populationID)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > luminosityTables(populationID)%metallicity(luminosityTables(populationID)%metallicitiesCount)) then
       iMetallicity=luminosityTables(populationID)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       iMetallicity=Interpolate_Locate(luminosityTables(populationID)%metallicity &
            &,luminosityTables(populationID)%interpolationAcceleratorMetallicity,metallicity &
            &,luminosityTables(populationID)%resetMetallicity)
       hMetallicity=Interpolate_Linear_Generate_Factors(luminosityTables(populationID)%metallicity ,iMetallicity,metallicity)
    end if
    ! Allocate arrays for ages and luminosities.
    call allocateArray(ages        ,[luminosityTables(populationID)%agesCount                      ])
    call allocateArray(luminosities,[luminosityTables(populationID)%agesCount,size(luminosityIndex)])
    ! Assign ages.
    ages=luminosityTables(populationID)%age
    ! Do the interpolation.
    luminosities(:,:)=0.0d0
    do iLuminosity=1,size(luminosityIndex)
       do jMetallicity=0,1
          luminosities                                     (:,                iLuminosity                             )= &
               & +luminosities                             (:,                iLuminosity                             )  &
               & +luminosityTables(populationID)%luminosity(  luminosityIndex(iLuminosity),:,iMetallicity+jMetallicity)  &
               & *hMetallicity                             (                                              jMetallicity)
       end do
    end do
    ! Release the read lock on the luminosity tables.
    call luminosityTableLock%setRead()
    ! Prevent interpolation from returning negative fluxes.
    luminosities=max(luminosities,0.0d0)
    return
  end subroutine Stellar_Population_Luminosity_Track
    
end module Stellar_Population_Luminosities
