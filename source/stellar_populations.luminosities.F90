!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  implicit none
  private
  public :: Stellar_Population_Luminosity

  type luminosityTable
     !% Structure for holding tables of simple stellar population luminosities.
     integer                                                            :: agesCount                         , metallicitiesCount
     logical                            , allocatable, dimension(:)     :: isTabulated
     double precision                   , allocatable, dimension(:)     :: age                               , metallicity
     double precision                   , allocatable, dimension(:,:,:) :: luminosity
     ! Interpolation structures.
     logical                                                            :: resetAge                   =.true., resetMetallicity                   =.true.
     type            (fgsl_interp_accel)                                :: interpolationAcceleratorAge       , interpolationAcceleratorMetallicity
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

  ! Optoin controlling writing of luminosities to file.
  logical                                                      :: stellarPopulationLuminosityStoreToFile

contains

  function Stellar_Population_Luminosity(luminosityIndex,filterIndex,postprocessingChainIndex,imfIndex,abundancesStellar,age,redshift)
    !% Returns the luminosity for a $1 M_\odot$ simple stellar population of given {\tt abundances} and {\tt age} drawn from IMF
    !% specified by {\tt imfIndex} and observed through the filter specified by {\tt filterIndex}.
    use Memory_Management
    use Stellar_Population_Spectra
    use Stellar_Population_Spectra_Postprocess
    use Instruments_Filters
    use Numerical_Integration
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    use Galacticus_Display
    use Galacticus_Input_Paths
    use Input_Parameters
    use Star_Formation_IMF
    use ISO_Varying_String
    use String_Handling
    use IO_HDF5
    use File_Utilities
    use MPI_Utilities
    implicit none
    integer                                                                                    , intent(in   ) :: filterIndex                  (:), imfIndex                   , &
         &                                                                                                        luminosityIndex              (:), postprocessingChainIndex(:)
    double precision                                                                           , intent(in   ) :: age                          (:), redshift                (:)
    type            (abundances                )                                               , intent(in   ) :: abundancesStellar
    double precision                                         , dimension(size(luminosityIndex))                :: Stellar_Population_Luminosity
    type            (luminosityTable           ), allocatable, dimension(:)                                    :: luminosityTablesTemporary
    double precision                            , allocatable, dimension(:,:,:)                                :: luminosityTemporary
    logical                                     , allocatable, dimension(:)                                    :: isTabulatedTemporary
    double precision                                         , dimension(2)                                    :: wavelengthRange
    double precision                                         , dimension(0:1)                                  :: hAge                            , hMetallicity
    integer                                                                                                    :: iAge                            , iLuminosity                , &
         &                                                                                                        iMetallicity                    , jAge                       , &
         &                                                                                                        jMetallicity                    , loopCount                  , &
         &                                                                                                        loopCountMaximum
    logical                                                                                                    :: computeTable                    , calculateLuminosity
    double precision                                                                                           :: ageLast                         , metallicity                , &
         &                                                                                                        normalization
    type            (c_ptr                     )                                                               :: parameterPointer
    type            (fgsl_function             )                                                               :: integrandFunction
    type            (fgsl_integration_workspace)                                                               :: integrationWorkspace
    type            (varying_string            )                                                               :: message                         , luminositiesFileName
    character       (len=16                    )                                                               :: datasetName                     , redshiftLabel
    type            (hdf5Object                )                                                               :: luminositiesFile

    ! Determine if we have created space for this IMF yet.
    !$omp critical (Luminosity_Tables_Initialize)
    if (.not.moduleInitialized) then
       ! Read the parameter controlling integration tolerance.
       !@ <inputParameter>
       !@   <name>stellarPopulationLuminosityIntegrationToleranceRelative</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The relative tolerance used when integrating the flux of stellar populations through filters.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationLuminosityIntegrationToleranceRelative',stellarPopulationLuminosityIntegrationToleranceRelative,defaultValue=1.0d-3)
       ! Read the parameter controlling storing luminosities to file.
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
       ! Flag that this module is now initialized.
       moduleInitialized=.true.
    end if

    if (allocated(luminosityTables)) then
       if (size(luminosityTables) < imfIndex) then
          call Move_Alloc(luminosityTables,luminosityTablesTemporary)
          allocate(luminosityTables(imfIndex))
          luminosityTables(1:size(luminosityTablesTemporary))=luminosityTablesTemporary
          deallocate(luminosityTablesTemporary)
          call Memory_Usage_Record(sizeof(luminosityTables(1)),blockCount=0)
       end if
    else
       allocate(luminosityTables(imfIndex))
       call Memory_Usage_Record(sizeof(luminosityTables))
    end if

    ! Determine if we have tabulated luminosities for this luminosityIndex in this IMF yet.
    do iLuminosity=1,size(luminosityIndex)
       if (allocated(luminosityTables(imfIndex)%isTabulated)) then
          if (size(luminosityTables(imfIndex)%isTabulated) >= luminosityIndex(iLuminosity)) then
             computeTable=.not.luminosityTables(imfIndex)%isTabulated(luminosityIndex(iLuminosity))
          else
             call Move_Alloc (luminosityTables(imfIndex)%isTabulated,isTabulatedTemporary)
             call Move_Alloc (luminosityTables(imfIndex)%luminosity ,luminosityTemporary )
             call Alloc_Array(luminosityTables(imfIndex)%isTabulated,[luminosityIndex(iLuminosity)])
             call Alloc_Array(luminosityTables(imfIndex)%luminosity ,[luminosityIndex(iLuminosity)&
                  &,luminosityTables(imfIndex)%agesCount,luminosityTables(imfIndex)%metallicitiesCount])
             luminosityTables(imfIndex)%isTabulated(1:size(isTabulatedTemporary)    )=isTabulatedTemporary
             luminosityTables(imfIndex)%isTabulated(  size(isTabulatedTemporary)+1:luminosityIndex(iLuminosity))=.false.
             luminosityTables(imfIndex)%luminosity (1:size(isTabulatedTemporary),:,:)=luminosityTemporary
             call Dealloc_Array(isTabulatedTemporary)
             call Dealloc_Array(luminosityTemporary)
             computeTable=.true.
          end if
       else
          call Alloc_Array(luminosityTables(imfIndex)%isTabulated,[luminosityIndex(iLuminosity)])
          luminosityTables(imfIndex)%isTabulated=.false.
          ! Since we have not yet tabulated any luminosities yet for this IMF, we need to get a list of suitable metallicities and
          ! ages at which to tabulate.
          call Stellar_Population_Spectrum_Tabulation(imfIndex,luminosityTables(imfIndex)%agesCount &
               &,luminosityTables(imfIndex)%metallicitiesCount,luminosityTables(imfIndex)%age&
               &,luminosityTables(imfIndex)%metallicity)
          where (luminosityTables(imfIndex)%metallicity > 0.0d0)
             luminosityTables(imfIndex)%metallicity=log10(luminosityTables(imfIndex)%metallicity/metallicitySolar)
          elsewhere
             luminosityTables(imfIndex)%metallicity=logMetallicityZero
          end where
          call Alloc_Array(luminosityTables(imfIndex)%luminosity,[luminosityIndex(iLuminosity)&
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
             !# </uniqueLabel>
             luminositiesFileName=Galacticus_Input_Path()                                                                                   // &
                  &               "data/stellarPopulations/stellarLuminosities::IMF:"                                                       // &
                  &               IMF_Name                                             (imfIndex                                           )// &
                  &               "::filter:"                                                                                               // &
                  &               Filter_Name                                          (filterIndex             (iLuminosity)              )// &
                  &               "::postprocessing:"                                                                                       // &
                  &               Stellar_Population_Spectrum_Postprocess_Chain_Methods(postprocessingChainIndex(iLuminosity)              )// &
                  &               "::dependencies:"                                                                                         // &
                  &               Stellar_Population_Luminosities_Label                (includeSourceDigest=.true.           ,asHash=.true.)// &
                  &               ".hdf5"
             if (File_Exists(luminositiesFileName)) then
                ! Construct the dataset name.
                write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                datasetName="redshift"//adjustl(trim(redshiftLabel))
                ! Open the file and check for the required dataset.
                !$omp critical (HDF5_Access)
                call luminositiesFile%openFile(char(luminositiesFileName),readOnly=.true.)
                if (luminositiesFile%hasDataset(trim(datasetName))) then
                   ! Read the dataset.
                   call luminositiesFile%readDatasetStatic(trim(datasetName),luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:))
                   ! We do not need to calculate this luminosity.
                   calculateLuminosity=.false.
                end if
                call luminositiesFile%close()
                !$omp end critical (HDF5_Access)
             end if
          end if
  
          ! Compute the luminosity if necessary.
          if (calculateLuminosity) then
             ! Display a message and counter.
             message='Tabulating stellar luminosities for '//char(IMF_Name(imfIndex))//' IMF, luminosity '
             message=message//iLuminosity//' of '//size(luminosityIndex)
             call Galacticus_Display_Indent (message,verbosityWorking)
             call Galacticus_Display_Counter(0,.true.,verbosityWorking)             
             ! Get wavelength extent of the filter.
             wavelengthRange=Filter_Extent(filterIndex(iLuminosity))
             ! Integrate over the wavelength range.
             filterIndexTabulate             =filterIndex             (iLuminosity)
             postprocessingChainIndexTabulate=postprocessingChainIndex(iLuminosity)
             redshiftTabulate                =redshift                (iLuminosity)
             imfIndexTabulate                =imfIndex
             loopCountMaximum                =luminosityTables(imfIndex)%metallicitiesCount*luminosityTables(imfIndex)%agesCount
             loopCount                       =0
             !$omp parallel do private(iAge,iMetallicity,integrandFunction,integrationWorkspace) copyin(filterIndexTabulate,postprocessingChainIndexTabulate,redshiftTabulate,imfIndexTabulate)
             do iAge=1,luminosityTables(imfIndex)%agesCount
                ageTabulate=luminosityTables(imfIndex)%age(iAge)
                do iMetallicity=1,luminosityTables(imfIndex)%metallicitiesCount
                   ! Update the counter.
                   !$omp atomic
                   loopCount=loopCount+1
                   call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountMaximum)),.false.,verbosityWorking)
                   call abundancesTabulate%metallicitySet(luminosityTables(imfIndex)%metallicity(iMetallicity) &
                        &,metallicityType=logarithmicByMassSolar)
                   luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),iAge,iMetallicity) &
                        &=Integrate(wavelengthRange(1),wavelengthRange(2),Filter_Luminosity_Integrand,parameterPointer &
                        &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative&
                        &=stellarPopulationLuminosityIntegrationToleranceRelative,integrationRule =FGSL_Integ_Gauss15,maxIntervals&
                        &=10000)
                   call Integrate_Done(integrandFunction,integrationWorkspace)
                end do
             end do
             !$omp end parallel do
             ! Clear the counter and write a completion message.
             call Galacticus_Display_Counter_Clear(           verbosityWorking)
             call Galacticus_Display_Unindent     ('finished',verbosityWorking)
             ! Get the normalization by integrating a zeroth magnitude (AB) source through the filter.
             normalization=Integrate(wavelengthRange(1),wavelengthRange(2),Filter_Luminosity_Integrand_AB,parameterPointer &
                  &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative&
                  &=stellarPopulationLuminosityIntegrationToleranceRelative)
             call Integrate_Done(integrandFunction,integrationWorkspace)
             ! Normalize the luminosity.
             luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:) &
                  &=luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:)/normalization
             ! Store the luminosities to file.
             if (stellarPopulationLuminosityStoreToFile.and.mpiSelf%isMaster()) then
                ! Construct the dataset name.
                write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                datasetName="redshift"//adjustl(trim(redshiftLabel))
                ! Open the file.
                !$omp critical (HDF5_Access)
                call luminositiesFile%openFile(char(luminositiesFileName))
                ! Write the dataset.
                call luminositiesFile%writeDataset(luminosityTables(imfIndex)%luminosity(luminosityIndex(iLuminosity),:,:),datasetName=trim(datasetName),commentText="Tabulated luminosities at redshift z="//adjustl(trim(redshiftLabel)))
                ! Close the file.
                call luminositiesFile%close()
                !$omp end critical (HDF5_Access)
             end if
          end if
          ! Flag that calculations have been performed for this filter.
          luminosityTables(imfIndex)%isTabulated(luminosityIndex(iLuminosity))=.true.
       end if
    end do

    ! Get interpolation in metallicity.
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=logarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < luminosityTables(imfIndex)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > luminosityTables(imfIndex)%metallicity(luminosityTables(imfIndex)%metallicitiesCount)) then
       iMetallicity=luminosityTables(imfIndex)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       iMetallicity=Interpolate_Locate(luminosityTables(imfIndex)%metallicitiesCount,luminosityTables(imfIndex)%metallicity &
            &,luminosityTables(imfIndex)%interpolationAcceleratorMetallicity,metallicity &
            &,luminosityTables(imfIndex)%resetMetallicity)
       hMetallicity=Interpolate_Linear_Generate_Factors(luminosityTables(imfIndex)%metallicitiesCount &
            &,luminosityTables(imfIndex)%metallicity ,iMetallicity,metallicity)
    end if

    ! Do the interpolation.
    Stellar_Population_Luminosity(:)= 0.0d0
    ageLast                         =-1.0d0
    do iLuminosity=1,size(luminosityIndex)
       ! Only compute luminosities for entries with positive age (negative age implies that the luminosity required is for a
       ! population observed prior to the formation of this population).
       if (age(iLuminosity) > 0.0d0) then
          ! Get interpolation in age if the age for this luminosity differs from the previous one.
          if (iLuminosity == 1 .or. age(iLuminosity) /= ageLast) then
             ! Check for out of range age.
             if (age(iLuminosity) > luminosityTables(imfIndex)%age(luminosityTables(imfIndex)%agesCount)) call&
                  & Galacticus_Error_Report('Stellar_Population_Luminosity','age exceeds the maximum tabulated')
             iAge=Interpolate_Locate(luminosityTables(imfIndex)%agesCount,luminosityTables(imfIndex)%age &
                  &,luminosityTables(imfIndex)%interpolationAcceleratorAge,age(iLuminosity),luminosityTables(imfIndex)%resetAge)
             hAge=Interpolate_Linear_Generate_Factors(luminosityTables(imfIndex)%agesCount,luminosityTables(imfIndex)%age,iAge&
                  &,age(iLuminosity))
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

  function Filter_Luminosity_Integrand(wavelength,parameterPointer) bind(c)
    !% Integrand for the luminosity through a given filter.
    use Stellar_Population_Spectra
    use Stellar_Population_Spectra_Postprocess
    use Instruments_Filters
    implicit none
    real            (kind=c_double)        :: Filter_Luminosity_Integrand
    real            (kind=c_double), value :: wavelength
    type            (c_ptr        ), value :: parameterPointer
    double precision                       :: wavelengthRedshifted

    ! If this luminosity is for a redshifted spectrum, then we shift wavelength at which we sample the stellar population spectrum
    ! to be a factor of (1+z) smaller. We therefore integrate over the stellar SED at shorter wavelengths, since these will be
    ! shifted into the filter by z=0. Factor of 1/wavelength appears since we want to integrate F_nu (dnu / nu) and dnu =
    ! -c/lambda^2 dlambda. Note that we follow the convention of Hogg et al. (2002) and assume that the filter response gives the
    ! fraction of incident photons received by the detector at a given wavelength, multiplied by the relative photon response
    ! (which will be 1 for a photon-counting detector such as a CCD, or proportional to the photon energy for a
    ! bolometer/calorimeter type detector).
    wavelengthRedshifted=wavelength/(1.0d0+redshiftTabulate)
    Filter_Luminosity_Integrand=Filter_Response(filterIndexTabulate,wavelength)*Stellar_Population_Spectrum(abundancesTabulate &
         &,ageTabulate,wavelengthRedshifted,imfIndexTabulate)*Stellar_Population_Spectrum_PostProcess(postprocessingChainIndexTabulate,wavelengthRedshifted,ageTabulate,redshiftTabulate)/wavelength
    return
  end function Filter_Luminosity_Integrand

  function Filter_Luminosity_Integrand_AB(wavelength,parameterPointer) bind(c)
    !% Integrand for the luminosity of a zeroth magnitude (AB) source through a given filter.
    use Instruments_Filters
    use Numerical_Constants_Astronomical
    implicit none
    real            (kind=c_double)            :: Filter_Luminosity_Integrand_AB
    real            (kind=c_double), value     :: wavelength
    type            (c_ptr        ), value     :: parameterPointer
    ! Luminosity of a zeroth magintude (AB) source in Solar luminosities per Hz.
    double precision               , parameter :: luminosityZeroPointABSolar    =luminosityZeroPointAB/luminositySolar

    Filter_Luminosity_Integrand_AB=Filter_Response(filterIndexTabulate,wavelength)*luminosityZeroPointABSolar/wavelength
    return
  end function Filter_Luminosity_Integrand_AB

end module Stellar_Population_Luminosities
