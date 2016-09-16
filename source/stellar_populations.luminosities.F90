!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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
  use FGSL
  use, intrinsic :: ISO_C_Binding
  use Abundances_Structure
  use ISO_Varying_String
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
  double precision                                             :: ageTabulate                                                    , redshiftTabulate
  integer                                                      :: filterIndexTabulate                                            , imfIndexTabulate, &
       &                                                          postprocessingChainIndexTabulate
  type            (abundances     )                            :: abundancesTabulate
  !$omp threadprivate(ageTabulate,redshiftTabulate,abundancesTabulate,filterIndexTabulate,imfIndexTabulate,postprocessingChainIndexTabulate)

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

contains

  subroutine Stellar_Population_Luminosity_Tabulate(luminosityIndex,filterIndex,postprocessingChainIndex,imfIndex,redshift)
    !% Tabulate stellar population luminosity in the given filters.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Astronomical
    use File_Utilities
    use IO_HDF5
    use String_Handling
    use Star_Formation_IMF
    use Input_Parameters
    use Galacticus_Input_Paths
    use Galacticus_Display
    use Galacticus_Error
    use Instruments_Filters
    use Numerical_Integration
    use Memory_Management
    use Stellar_Population_Spectra
    use Stellar_Population_Spectra_Postprocess
    use Input_Parameters2
    implicit none
    integer                                                                                       , intent(in   ) :: filterIndex                              (:), imfIndex                   , &
         &                                                                                                           luminosityIndex                          (:), postprocessingChainIndex(:)
    double precision                                                                              , intent(in   ) :: redshift                                 (:)
    type            (luminosityTable              ), allocatable, dimension(:)                                    :: luminosityTablesTemporary
    double precision                               , allocatable, dimension(:,:,:)                                :: luminosityTemporary
    logical                                        , allocatable, dimension(:)                                    :: isTabulatedTemporary
    double precision                                            , dimension(2)                                    :: wavelengthRange
    type            (lockDescriptor               )                                                               :: lockFileDescriptor
    class           (stellarPopulationSpectraClass), pointer                                                      :: stellarPopulationSpectra_
    integer         (c_size_t                     )                                                               :: iAge                                        , iLuminosity                , &
         &                                                                                                           iMetallicity                                , jLuminosity
    integer                                                                                                       :: loopCountMaximum                            , loopCount                  , &
         &                                                                                                           errorStatus                                 , luminosityIndexMaximum
    logical                                                                                                       :: computeTable                                , calculateLuminosity        , &
         &                                                                                                           stellarLuminositiesUniqueLabelConstructed
    double precision                                                                                              :: toleranceRelative                           , normalization
    type            (fgsl_function                )                                                               :: integrandFunction
    type            (fgsl_integration_workspace   )                                                               :: integrationWorkspace
    type            (varying_string               )                                                               :: message                                     , luminositiesFileName       , &
         &                                                                                                           stellarLuminositiesUniqueLabel
    character       (len=16                       )                                                               :: datasetName                                 , redshiftLabel              , &
         &                                                                                                           label
    type            (hdf5Object                   )                                                               :: luminositiesFile
    
    ! Determine if we have created space for this IMF yet.
    !$omp critical (Luminosity_Tables_Initialize)
    if (.not.moduleInitialized) then
       ! Read the parameter controlling integration tolerance.
       !@ <inputParameter>
       !@   <name>stellarPopulationLuminosityIntegrationToleranceRelative</name>
       !@   <defaultValue>$2 \times 10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The relative tolerance used when integrating the flux of stellar populations through filters.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationLuminosityIntegrationToleranceRelative',stellarPopulationLuminosityIntegrationToleranceRelative,defaultValue=2.0d-3)
       !@ <inputParameter>
       !@   <name>stellarPopulationLuminosityIntegrationToleranceDegrade</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    If {\normalfont \ttfamily true}, automatically degrade the relative tolerance used when integrating the flux of stellar populations through filters to ensure convergence.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationLuminosityIntegrationToleranceDegrade',stellarPopulationLuminosityIntegrationToleranceDegrade,defaultValue=.false.)
       ! Read parameters controlling storing luminosities to file.
       !@ <inputParameter>
       !@   <name>stellarPopulationLuminosityStoreToFile</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not stellar populations luminosities (integrated under a filter) should be stored to file for rapid reuse.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationLuminosityStoreToFile',stellarPopulationLuminosityStoreToFile,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>stellarPopulationLuminosityStoreDirectory</name>
       !@   <defaultValue>{\normalfont \ttfamily \$GALACTICUS\_ROOT\_V094/data/stellarPopulations}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies the directory to which stellar populations luminosities (integrated under a filter) should be stored to file for rapid reuse.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationLuminosityStoreDirectory',stellarPopulationLuminosityStoreDirectory,defaultValue=char(Galacticus_Input_Path()//"data/stellarPopulations"))
       ! Read the parameter controlling behavior if maximum age of stellar populations are exceeded.
       !@ <inputParameter>
       !@   <name>stellarPopulationLuminosityMaximumAgeExceededIsFatal</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not exceeding the maximum available age of the stellar population is fatal.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationLuminosityMaximumAgeExceededIsFatal',stellarPopulationLuminosityMaximumAgeExceededIsFatal,defaultValue=.true.)
       ! Flag that this module is now initialized.
       moduleInitialized=.true.
    end if

    if (allocated(luminosityTables)) then
       if (size(luminosityTables) < imfIndex) then
          call Move_Alloc(luminosityTables,luminosityTablesTemporary)
          allocate(luminosityTables(imfIndex))
          luminosityTables(1:size(luminosityTablesTemporary))=luminosityTablesTemporary
          luminosityTables(size(luminosityTablesTemporary)+1:imfIndex)%isTabulatedMaximum=0
          deallocate(luminosityTablesTemporary)
          call Memory_Usage_Record(sizeof(luminosityTables(1)),blockCount=0)
       end if
    else
       allocate(luminosityTables(imfIndex))
       luminosityTables%isTabulatedMaximum=0
       call Memory_Usage_Record(sizeof(luminosityTables))
    end if

    ! Determine if we have tabulated luminosities for this luminosityIndex in this IMF yet.
    luminosityIndexMaximum=maxval(luminosityIndex)    
    if (.not.allocated(luminosityTables(imfIndex)%isTabulated) .or. luminosityTables(imfIndex)%isTabulatedMaximum < luminosityIndexMaximum) then
       call File_Lock_Initialize(lockFileDescriptor)
       stellarLuminositiesUniqueLabelConstructed = .false.
       stellarPopulationSpectra_                 => null()
       luminosityIndexMaximum                    =  maxval(luminosityIndex)
       do iLuminosity=1,size(luminosityIndex)
          if (allocated(luminosityTables(imfIndex)%isTabulated)) then
             if (size(luminosityTables(imfIndex)%isTabulated) >= luminosityIndex(iLuminosity)) then
                computeTable=.not.luminosityTables(imfIndex)%isTabulated(luminosityIndex(iLuminosity))
             else
                call Move_Alloc (luminosityTables(imfIndex)%isTabulated,isTabulatedTemporary)
                call Move_Alloc (luminosityTables(imfIndex)%luminosity ,luminosityTemporary )
                call allocateArray(luminosityTables(imfIndex)%isTabulated,[luminosityIndexMaximum])
                call allocateArray(luminosityTables(imfIndex)%luminosity ,[luminosityIndexMaximum&
                     &,luminosityTables(imfIndex)%agesCount,luminosityTables(imfIndex)%metallicitiesCount])
                luminosityTables(imfIndex)%isTabulated(1:size(isTabulatedTemporary)    )=isTabulatedTemporary
                luminosityTables(imfIndex)%isTabulated(  size(isTabulatedTemporary)+1:luminosityIndexMaximum)=.false.
                luminosityTables(imfIndex)%luminosity (1:size(isTabulatedTemporary),:,:)=luminosityTemporary
                call deallocateArray(isTabulatedTemporary)
                call deallocateArray(luminosityTemporary)
                computeTable=.true.
             end if
          else
             call allocateArray(luminosityTables(imfIndex)%isTabulated,[luminosityIndexMaximum])
             luminosityTables(imfIndex)%isTabulated=.false.
             ! Since we have not yet tabulated any luminosities yet for this IMF, we need to get a list of suitable metallicities and
             ! ages at which to tabulate.
             if (.not.associated(stellarPopulationSpectra_)) stellarPopulationSpectra_ => stellarPopulationSpectra()
             call stellarPopulationSpectra_%tabulation(imfIndex,luminosityTables(imfIndex)%agesCount &
                  &,luminosityTables(imfIndex)%metallicitiesCount,luminosityTables(imfIndex)%age&
                  &,luminosityTables(imfIndex)%metallicity)
             where (luminosityTables(imfIndex)%metallicity > 0.0d0)
                luminosityTables(imfIndex)%metallicity=log10(luminosityTables(imfIndex)%metallicity/metallicitySolar)
             elsewhere
                luminosityTables(imfIndex)%metallicity=logMetallicityZero
             end where
             call allocateArray(luminosityTables(imfIndex)%luminosity,[luminosityIndexMaximum&
                  &,luminosityTables(imfIndex)%agesCount ,luminosityTables(imfIndex)%metallicitiesCount])
             computeTable=.true.
          end if

          ! If we haven't, do so now.
          if (computeTable) then

             ! Determine if we can read the required luminosity from file.
             calculateLuminosity=.true.
             if (stellarPopulationLuminosityStoreToFile) then
                ! Construct name of the file to which this would be stored.
                !# <uniqueLabel>
                !#  <function>Stellar_Population_Luminosities_Label</function>
                !#  <ignoreRegex>^stellarPopulationSpectraPostprocess.*</ignoreRegex>
                !#  <ignoreRegex>^starFormationImf.*</ignoreRegex>
                !#  <ignoreRegex>^imf.*</ignoreRegex>
                !#  <ignoreRegex>^stellarAstrophysics.*</ignoreRegex>
                !#  <ignoreRegex>^stellarFeedback.*</ignoreRegex>
                !#  <ignoreRegex>^stellarProperties.*</ignoreRegex>
                !#  <ignoreRegex>^stellarWinds.*</ignoreRegex>
                !#  <ignoreRegex>^stellarTracks.*</ignoreRegex>
                !#  <ignoreRegex>^supernovae.*Method$</ignoreRegex>
                !#  <ignore>supernovaEnergy</ignore>
                !#  <ignore>initialMassForSupernovaeTypeII</ignore>
                !#  <ignore>elementsToTrack</ignore>
                !#  <ignore>stellarPopulationLuminosityStoreDirectory</ignore>
                !# </uniqueLabel>
                if (.not.stellarLuminositiesUniqueLabelConstructed) then
                   stellarLuminositiesUniqueLabel=Stellar_Population_Luminosities_Label(includeSourceDigest=.true.,asHash=.true.)
                   stellarLuminositiesUniqueLabelConstructed=.true.
                end if
                luminositiesFileName=stellarPopulationLuminosityStoreDirectory                                                   // &
                     &               "/stellarLuminosities::IMF:"                                                                // &
                     &               IMF_Name                                             (imfIndex                             )// &
                     &               "::filter:"                                                                                 // &
                     &               Filter_Name                                          (filterIndex             (iLuminosity))// &
                     &               "::postprocessing:"                                                                         // &
                     &               Stellar_Population_Spectrum_Postprocess_Chain_Methods(postprocessingChainIndex(iLuminosity))// &
                     &               "::dependencies:"                                                                           // &
                     &               stellarLuminositiesUniqueLabel                                                              // &
                     &               ".hdf5"
                if (File_Exists(luminositiesFileName)) then
                   ! Construct the dataset name.
                   write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                   datasetName="redshift"//adjustl(trim(redshiftLabel))
                   ! Open the file and check for the required dataset.
                   !$omp critical (HDF5_Access)
                   call File_Lock(char(luminositiesFileName),lockFileDescriptor,lockIsShared=.true.)
                   call luminositiesFile%openFile(char(luminositiesFileName),readOnly=.true.)
                   if (luminositiesFile%hasDataset(trim(datasetName))) then
                      ! Read the dataset.
                      call luminositiesFile%readDatasetStatic(trim(datasetName),luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:))
                      ! We do not need to calculate this luminosity.
                      calculateLuminosity=.false.
                   end if
                   call luminositiesFile%close()
                   call File_Unlock(lockFileDescriptor)
                   !$omp end critical (HDF5_Access)
                end if
             end if

             ! Compute the luminosity if necessary.
             if (calculateLuminosity) then
                ! Display a message and counter.
                message='Tabulating stellar luminosities for '//char(IMF_Name(imfIndex))//' IMF, luminosity '
                write (redshiftLabel,'(f6.3)') redshift(iLuminosity)
                message=message                                                                                     // &
                     &  Filter_Name                                          (filterIndex             (iLuminosity))// &
                     &  ":"                                                                                         // &
                     &  Stellar_Population_Spectrum_Postprocess_Chain_Methods(postprocessingChainIndex(iLuminosity))// &
                     &  ":z"                                                                                        // &
                     &  trim(adjustl(redshiftLabel))                                                                // &
                     &  " "                                                                                         // &
                     &                                                                                 iLuminosity  // &
                     &  " of "                                                                                      // &
                     &  size(luminosityIndex)
                call Galacticus_Display_Indent (message,verbosityWorking)
                call Galacticus_Display_Counter(0,.true.,verbosityWorking)             
                ! Get stellar population spectra object if necessary.
                if (.not.associated(stellarPopulationSpectra_)) stellarPopulationSpectra_ => stellarPopulationSpectra()
                ! Get wavelength extent of the filter.
                wavelengthRange=Filter_Extent(filterIndex(iLuminosity))
                ! Integrate over the wavelength range.
                filterIndexTabulate             =filterIndex             (iLuminosity)
                postprocessingChainIndexTabulate=postprocessingChainIndex(iLuminosity)
                redshiftTabulate                =redshift                (iLuminosity)
                imfIndexTabulate                =imfIndex
                loopCountMaximum                =luminosityTables(imfIndex)%metallicitiesCount*luminosityTables(imfIndex)%agesCount
                loopCount                       =0
                !$omp parallel do private(iAge,iMetallicity,integrandFunction,integrationWorkspace,toleranceRelative,errorStatus) copyin(filterIndexTabulate,postprocessingChainIndexTabulate,redshiftTabulate,imfIndexTabulate)
                do iAge=1,luminosityTables(imfIndex)%agesCount
                   ageTabulate=luminosityTables(imfIndex)%age(iAge)
                   do iMetallicity=1,luminosityTables(imfIndex)%metallicitiesCount
                      ! Update the counter.
                      !$omp atomic
                      loopCount=loopCount+1
                      call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountMaximum)),.false.,verbosityWorking)
                      call abundancesTabulate%metallicitySet(luminosityTables(imfIndex)%metallicity(iMetallicity) &
                           &,metallicityType=metallicityTypeLogarithmicByMassSolar)
                      toleranceRelative=stellarPopulationLuminosityIntegrationToleranceRelative
                      errorStatus      =errorStatusFail
                      do while (errorStatus /= errorStatusSuccess)
                         luminosityTables(imfIndex)%luminosity(                              &
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
                               call Galacticus_Error_Report('Stellar_Population_Luminosity',message)
                            else
                               write (label,'(e9.3)') 2.0d0*stellarPopulationLuminosityIntegrationToleranceRelative
                               message=         "integration of stellar populations failed"                                              //char(10)
                               message=message//"HELP: consider increasing the [stellarPopulationLuminosityIntegrationToleranceRelative]"//char(10)
                               message=message//"      parameter to "//trim(adjustl(label))//" to reduce the integration tolerance"      //char(10)
                               message=message//"      required if your can accept this lower accuracy."
                               call Galacticus_Error_Report('Stellar_Population_Luminosity',message)
                            end if
                         end if
                      end do
                   end do
                end do
                !$omp end parallel do
                ! Clear the counter and write a completion message.
                call Galacticus_Display_Counter_Clear(           verbosityWorking)
                call Galacticus_Display_Unindent     ('finished',verbosityWorking)
                ! Get the normalization by integrating a zeroth magnitude (AB) source through the filter.
                normalization=Integrate(wavelengthRange(1),wavelengthRange(2),Filter_Luminosity_Integrand_AB &
                     &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative&
                     &=stellarPopulationLuminosityIntegrationToleranceRelative)
                call Integrate_Done(integrandFunction,integrationWorkspace)
                ! Normalize the luminosity.
                luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:) &
                     &=luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:)/normalization
                ! Store the luminosities to file.
                if (stellarPopulationLuminosityStoreToFile) then
                   ! Construct the dataset name.
                   write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                   datasetName="redshift"//adjustl(trim(redshiftLabel))
                   ! Open the file.
                   !$omp critical (HDF5_Access)
                   call File_Lock(char(luminositiesFileName),lockFileDescriptor,lockIsShared=.false.)
                   call luminositiesFile%openFile(char(luminositiesFileName))
                   ! Write the dataset.
                   if (.not.luminositiesFile%hasDataset(trim(datasetName))) &
                        & call luminositiesFile%writeDataset(luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:),datasetName=trim(datasetName),commentText="Tabulated luminosities at redshift z="//adjustl(trim(redshiftLabel)))
                   ! Close the file.
                   call luminositiesFile%close()
                   call File_Unlock(lockFileDescriptor)
                   !$omp end critical (HDF5_Access)
                end if
             end if
             ! Flag that calculations have been performed for this filter.
             luminosityTables(imfIndex)%isTabulated(luminosityIndex(iLuminosity))=.true.
             if (luminosityIndex(iLuminosity) > luminosityTables(imfIndex)%isTabulatedMaximum) then
                jLuminosity=luminosityTables(imfIndex)%isTabulatedMaximum
                do while (jLuminosity < size(luminosityTables(imfIndex)%isTabulated) .and. luminosityTables(imfIndex)%isTabulated(min(jLuminosity+1,size(luminosityTables(imfIndex)%isTabulated))))
                   jLuminosity=jLuminosity+1
                end do
                luminosityTables(imfIndex)%isTabulatedMaximum=jLuminosity
             end if
          end if
       end do
    end if
    !$omp end critical (Luminosity_Tables_Initialize)
    return

  contains

    double precision function Filter_Luminosity_Integrand(wavelength)
      !% Integrand for the luminosity through a given filter.
      use Stellar_Population_Spectra_Postprocess
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
      Filter_Luminosity_Integrand=Filter_Response(filterIndexTabulate,wavelength)*stellarPopulationSpectra_%luminosity(abundancesTabulate &
           &,ageTabulate,wavelengthRedshifted,imfIndexTabulate)*Stellar_Population_Spectrum_PostProcess(postprocessingChainIndexTabulate,wavelengthRedshifted,ageTabulate,redshiftTabulate)/wavelength
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
  
  function Stellar_Population_Luminosity(luminosityIndex,filterIndex,postprocessingChainIndex,imfIndex,abundancesStellar,age,redshift)
    !% Returns the luminosity for a $1 M_\odot$ simple stellar population of given {\normalfont \ttfamily abundances} and {\normalfont \ttfamily age} drawn from IMF
    !% specified by {\normalfont \ttfamily imfIndex} and observed through the filter specified by {\normalfont \ttfamily filterIndex}.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    use Numerical_Interpolation
    implicit none
    integer                                                                                    , intent(in   ) :: filterIndex                  (:), imfIndex                   , &
         &                                                                                                        luminosityIndex              (:), postprocessingChainIndex(:)
    double precision                                                                           , intent(in   ) :: age                          (:), redshift                (:)
    type            (abundances                )                                               , intent(in   ) :: abundancesStellar
    double precision                                         , dimension(size(luminosityIndex))                :: Stellar_Population_Luminosity
    double precision                                         , dimension(0:1)                                  :: hAge                            , hMetallicity
    integer         (c_size_t                  )                                                               :: iAge                            , iLuminosity                , &
         &                                                                                                        iMetallicity                    , jAge                       , &
         &                                                                                                        jMetallicity
    double precision                                                                                           :: ageLast                         , metallicity

    ! Tabulate the luminosities.
    call Stellar_Population_Luminosity_Tabulate(luminosityIndex,filterIndex,postprocessingChainIndex,imfIndex,redshift)
    ! Get interpolation in metallicity.
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < luminosityTables(imfIndex)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > luminosityTables(imfIndex)%metallicity(luminosityTables(imfIndex)%metallicitiesCount)) then
       iMetallicity=luminosityTables(imfIndex)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       iMetallicity=Interpolate_Locate(luminosityTables(imfIndex)%metallicity &
            &,luminosityTables(imfIndex)%interpolationAcceleratorMetallicity,metallicity &
            &,luminosityTables(imfIndex)%resetMetallicity)
       hMetallicity=Interpolate_Linear_Generate_Factors(luminosityTables(imfIndex)%metallicity ,iMetallicity,metallicity)
    end if

    ! Do the interpolation.
    Stellar_Population_Luminosity(:)= 0.0d0
    ageLast                         =-1.0d0
    !$omp critical (Luminosity_Tables_Initialize)
    do iLuminosity=1,size(luminosityIndex)
       ! Only compute luminosities for entries with positive age (negative age implies that the luminosity required is for a
       ! population observed prior to the formation of this population).
       if (age(iLuminosity) > 0.0d0) then
          ! Get interpolation in age if the age for this luminosity differs from the previous one.
          if (iLuminosity == 1 .or. age(iLuminosity) /= ageLast) then
             ! Check for out of range age.
             if (age(iLuminosity) > luminosityTables(imfIndex)%age(luminosityTables(imfIndex)%agesCount)) then
                if (stellarPopulationLuminosityMaximumAgeExceededIsFatal) then
                   call Galacticus_Error_Report('Stellar_Population_Luminosity','age exceeds the maximum tabulated')
                else
                   iAge=luminosityTables(imfIndex)%agesCount-1
                   hAge=[0.0d0,1.0d0]
                end if
             else
                iAge=Interpolate_Locate(luminosityTables(imfIndex)%age &
                     &,luminosityTables(imfIndex)%interpolationAcceleratorAge,age(iLuminosity),luminosityTables(imfIndex)%resetAge)
                hAge=Interpolate_Linear_Generate_Factors(luminosityTables(imfIndex)%age,iAge&
                     &,age(iLuminosity))
             end if
             ageLast=age(iLuminosity)
          end if
          do jAge=0,1
             do jMetallicity=0,1
                Stellar_Population_Luminosity(iLuminosity)=Stellar_Population_Luminosity(iLuminosity)&
                     &+luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),iAge +jAge,iMetallicity+jMetallicity)&
                     &*hAge(jAge)*hMetallicity(jMetallicity)
             end do
          end do
       end if
    end do
    !$omp end critical (Luminosity_Tables_Initialize)
    ! Prevent interpolation from returning negative fluxes.
    Stellar_Population_Luminosity=max(Stellar_Population_Luminosity,0.0d0)
    return
  end function Stellar_Population_Luminosity
  
  subroutine Stellar_Population_Luminosity_Track(luminosityIndex,filterIndex,postprocessingChainIndex,imfIndex,abundancesStellar,redshift,ages,luminosities)
    !% Returns the luminosity for a $1 M_\odot$ simple stellar population of given {\normalfont \ttfamily abundances} drawn from IMF
    !% specified by {\normalfont \ttfamily imfIndex} and observed through the filter specified by {\normalfont \ttfamily filterIndex}, for all available ages.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    use Memory_Management
    use Numerical_Interpolation
    implicit none
    integer                                                                    , intent(in   ) :: filterIndex      (:), imfIndex                   , &
         &                                                                                        luminosityIndex  (:), postprocessingChainIndex(:)
    double precision                                                           , intent(in   ) :: redshift         (:)
    type            (abundances                )                               , intent(in   ) :: abundancesStellar
    double precision                            , allocatable, dimension(:    ), intent(  out) :: ages
    double precision                            , allocatable, dimension(:  ,:), intent(  out) :: luminosities
    double precision                                         , dimension(0:1  )                :: hMetallicity
    integer         (c_size_t                  )                                               :: iLuminosity         , iMetallicity               , &
         &                                                                                        jMetallicity
    double precision                                                                           :: metallicity

    ! Tabulate the luminosities.
    call Stellar_Population_Luminosity_Tabulate(luminosityIndex,filterIndex,postprocessingChainIndex,imfIndex,redshift)
    ! Get interpolation in metallicity.
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < luminosityTables(imfIndex)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > luminosityTables(imfIndex)%metallicity(luminosityTables(imfIndex)%metallicitiesCount)) then
       iMetallicity=luminosityTables(imfIndex)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       iMetallicity=Interpolate_Locate(luminosityTables(imfIndex)%metallicity &
            &,luminosityTables(imfIndex)%interpolationAcceleratorMetallicity,metallicity &
            &,luminosityTables(imfIndex)%resetMetallicity)
       hMetallicity=Interpolate_Linear_Generate_Factors(luminosityTables(imfIndex)%metallicity ,iMetallicity,metallicity)
    end if
    ! Allocate arrays for ages and luminosities.
    call allocateArray(ages        ,[luminosityTables(imfIndex)%agesCount                      ])
    call allocateArray(luminosities,[luminosityTables(imfIndex)%agesCount,size(luminosityIndex)])
    ! Assign ages.
    ages=luminosityTables(imfIndex)%age
    ! Do the interpolation.
    luminosities(:,:)=0.0d0
    !$omp critical (Luminosity_Tables_Initialize)
    do iLuminosity=1,size(luminosityIndex)
       do jMetallicity=0,1
          luminosities                                 (:,                iLuminosity                             )= &
               & +luminosities                         (:,                iLuminosity                             )  &
               & +luminosityTables(imfIndex)%luminosity(  luminosityIndex(iLuminosity),:,iMetallicity+jMetallicity)  &
               & *hMetallicity                         (                                              jMetallicity)
       end do
    end do
    !$omp end critical (Luminosity_Tables_Initialize)
    ! Prevent interpolation from returning negative fluxes.
    luminosities=max(luminosities,0.0d0)
    return
  end subroutine Stellar_Population_Luminosity_Track
    
end module Stellar_Population_Luminosities
