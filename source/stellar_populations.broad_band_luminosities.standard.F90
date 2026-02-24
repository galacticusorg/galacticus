!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  Implements the standard stellar populations broad band luminosities class.
  !!}
  
  use, intrinsic :: ISO_C_Binding                         , only : c_size_t
  use            :: Locks                                 , only : ompReadWriteLock
  use            :: Numerical_Interpolation               , only : interpolator
  use            :: Stellar_Population_Spectra            , only : stellarPopulationSpectraClass
  use            :: Stellar_Population_Spectra_Postprocess, only : stellarPopulationSpectraPostprocessorClass

  type luminosityTable
     !!{
     Structure for holding tables of simple stellar population luminosities.
     !!}
     integer                                                       :: agesCount           , metallicitiesCount
     integer         (c_size_t    )                                :: isTabulatedMaximum=0
     logical                       , allocatable, dimension(:    ) :: isTabulated
     double precision              , allocatable, dimension(:    ) :: age                 , metallicity
     double precision              , allocatable, dimension(:,:,:) :: luminosity
     type            (interpolator)                                :: interpolatorAge     , interpolatorMetallicity
  end type luminosityTable

  !![
  <stellarPopulationBroadBandLuminosities name="stellarPopulationBroadBandLuminositiesStandard">
   <description>The standard stellar populations broad band luminosities class.</description>
  </stellarPopulationBroadBandLuminosities>
  !!]
  type, extends(stellarPopulationBroadBandLuminositiesClass) :: stellarPopulationBroadBandLuminositiesStandard
     !!{
     The standard stellar population broad band luminosities class.
     !!}
     private
     type            (luminosityTable ), allocatable, dimension(:) :: luminosityTables
     double precision                                              :: integrationToleranceRelative
     logical                                                       :: integrationToleranceDegrade
     logical                                                       :: storeToFile
     type            (varying_string  )                            :: storeDirectory
     logical                                                       :: maximumAgeExceededIsFatal
     type            (ompReadWriteLock)                            :: luminosityTableLock
   contains
     !![
     <methods>
       <method description="Tabulate broad band stellar luminosities." method="tabulate" />
     </methods>
     !!]
     procedure :: luminosities     => standardLuminosities
     procedure :: luminosityTracks => standardLuminosityTracks
     procedure :: tabulate         => standardTabulate
  end type stellarPopulationBroadBandLuminositiesStandard

  interface stellarPopulationBroadBandLuminositiesStandard
     !!{
     Constructors for the \refClass{stellarPopulationBroadBandLuminositiesStandard} stellar population broad band luminosities class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface stellarPopulationBroadBandLuminositiesStandard

  ! Module scope variables used in integrations.
  double precision                                                      :: age_                                   , redshift_
  integer                                                               :: indexFilter_
  integer         (c_size_t                                  )          :: populationID_
  type            (abundances                                )          :: abundances_
  class           (stellarPopulationSpectraPostprocessorClass), pointer :: stellarPopulationSpectraPostprocessor__
  class           (stellarPopulationSpectraClass             ), pointer :: stellarPopulationSpectra__
  type            (interpolator                              ), pointer :: filterResponse_
  !$omp threadprivate(age_,redshift_,abundances_,indexFilter_,populationID_,stellarPopulationSpectraPostprocessor__,stellarPopulationSpectra__,filterResponse_)

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationBroadBandLuminositiesStandard} stellar population broad band luminosities class which takes a
    parameter set as input.
    !!}
    use :: Input_Paths     , only : inputPath      , pathTypeDataDynamic
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (stellarPopulationBroadBandLuminositiesStandard)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    double precision                                                                :: integrationToleranceRelative
    logical                                                                         :: storeToFile                 , integrationToleranceDegrade, &
         &                                                                             maximumAgeExceededIsFatal
    type            (varying_string                                )                :: storeDirectory

    !![
    <inputParameter>
      <name>integrationToleranceRelative</name>
      <defaultValue>4.0d-3</defaultValue>
      <description>The relative tolerance used when integrating the flux of stellar populations through filters.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>integrationToleranceDegrade</name>
      <defaultValue>.false.</defaultValue>
      <description>If {\normalfont \ttfamily true}, automatically degrade the relative tolerance used when integrating the flux of stellar populations through filters to ensure convergence.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>storeToFile</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not stellar populations luminosities (integrated under a filter) should be stored to file for rapid reuse.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>storeDirectory</name>
      <defaultValue>inputPath(pathTypeDataDynamic)//'stellarPopulations'</defaultValue>
      <description>Specifies the directory to which stellar populations luminosities (integrated under a filter) should be stored to file for rapid reuse.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>maximumAgeExceededIsFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not exceeding the maximum available age of the stellar population is fatal.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarPopulationBroadBandLuminositiesStandard(integrationToleranceRelative,integrationToleranceDegrade,maximumAgeExceededIsFatal,storeToFile,storeDirectory)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(integrationToleranceRelative,integrationToleranceDegrade,maximumAgeExceededIsFatal,storeToFile,storeDirectory) result(self)
    !!{
    Internal constructor for the \refClass{stellarPopulationBroadBandLuminositiesStandard} stellar population broad band luminosities class.
    !!}
    implicit none
    type            (stellarPopulationBroadBandLuminositiesStandard)                :: self
    double precision                                                , intent(in   ) :: integrationToleranceRelative
    logical                                                         , intent(in   ) :: storeToFile                 , integrationToleranceDegrade, &
         &                                                                             maximumAgeExceededIsFatal
    type            (varying_string                                ), intent(in   ) :: storeDirectory
    !![
    <constructorAssign variables="integrationToleranceRelative, integrationToleranceDegrade, maximumAgeExceededIsFatal,storeToFile,storeDirectory"/>
    !!]
    
    self%luminosityTableLock=ompReadWriteLock()
    return
  end function standardConstructorInternal

  function standardLuminosities(self,luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,abundancesStellar,age,redshift)
    !!{
    Returns the luminosity for a $1 M_\odot$ simple {\normalfont \ttfamily stellarPopulation\_} of given {\normalfont \ttfamily
    abundances} and {\normalfont \ttfamily age} and observed through the filter specified by {\normalfont \ttfamily
    filterIndex}.
    !!}
    use            :: Abundances_Structure, only : Abundances_Get_Metallicity, logMetallicityZero, metallicityTypeLogarithmicByMassSolar
    use            :: Error               , only : Error_Report
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: Stellar_Populations , only : stellarPopulationClass
    implicit none
    class           (stellarPopulationBroadBandLuminositiesStandard)                                  , intent(inout) :: self
    integer                                                         , dimension( :                   ), intent(in   ) :: luminosityIndex                       , filterIndex
    double precision                                                , dimension( :                   ), intent(in   ) :: age                                   , redshift
    type            (abundances                                    )                                  , intent(in   ) :: abundancesStellar
    type            (stellarPopulationSpectraPostprocessorList     ), dimension( :                   ), intent(inout) :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                        )                                  , intent(inout) :: stellarPopulation_
    double precision                                                , dimension(size(luminosityIndex))                :: standardLuminosities
    double precision                                                , dimension(0:1                  )                :: hAge                                  , hMetallicity
    integer         (c_size_t                                      )                                                  :: iAge                                  , iMetallicity   , &
         &                                                                                                               iLuminosityStart                      , iLuminosityEnd , &
         &                                                                                                               populationID
    double precision                                                                                                  :: metallicity

    ! Tabulate the luminosities.
    call self%tabulate(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    ! Obtain a read lock on the luminosity tables.
    call self%luminosityTableLock%setRead()
    ! Get interpolation in metallicity.
    populationID=stellarPopulation_%uniqueID()
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < self%luminosityTables(populationID)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > self%luminosityTables(populationID)%metallicity(self%luminosityTables(populationID)%metallicitiesCount)) then
       iMetallicity=self%luminosityTables(populationID)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       call self%luminosityTables(populationID)%interpolatorMetallicity%linearFactors(metallicity,iMetallicity,hMetallicity)
    end if
    ! Do the interpolations.
    iLuminosityStart=0
    do while (iLuminosityStart < size(luminosityIndex))
       iLuminosityStart=iLuminosityStart+1
       ! Find all contiguous luminosities at this age.
       iLuminosityEnd=iLuminosityStart
       do while (age(iLuminosityEnd) == age(iLuminosityStart) .or. (age(iLuminosityStart) < 0.0d0 .and. age(iLuminosityEnd) < 0.0d0))
          iLuminosityEnd=iLuminosityEnd+1
          if (iLuminosityEnd > size(luminosityIndex)) exit
       end do
       iLuminosityEnd=iLuminosityEnd-1
       ! Only compute luminosities for entries with positive age (negative age implies that the luminosity required is for a
       ! population observed prior to the formation of this population).
       if (age(iLuminosityStart) < 0.0d0) then
          standardLuminosities(iLuminosityStart:iLuminosityEnd)=0.0d0
       else
          ! Check for out of range age.
          if (age(iLuminosityStart) > self%luminosityTables(populationID)%age(self%luminosityTables(populationID)%agesCount)) then
             if (self%maximumAgeExceededIsFatal) then
                call Error_Report('age exceeds the maximum tabulated'//{introspection:location})
             else
                iAge=self%luminosityTables(populationID)%agesCount-1
                hAge=[0.0d0,1.0d0]
             end if
          else
             call self%luminosityTables(populationID)%interpolatorAge%linearFactors(age(iLuminosityStart),iAge,hAge)
          end if
          standardLuminosities(iLuminosityStart:iLuminosityEnd)=                                                                                                 &
               & +self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosityStart:iLuminosityEnd),iAge+0,iMetallicity+0)*hAge(0)*hMetallicity(0) &
               & +self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosityStart:iLuminosityEnd),iAge+0,iMetallicity+1)*hAge(0)*hMetallicity(1) &
               & +self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosityStart:iLuminosityEnd),iAge+1,iMetallicity+0)*hAge(1)*hMetallicity(0) &
               & +self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosityStart:iLuminosityEnd),iAge+1,iMetallicity+1)*hAge(1)*hMetallicity(1)
       end if
       iLuminosityStart=iLuminosityEnd
    end do
    call self%luminosityTableLock%unsetRead()
    ! Prevent interpolation from returning negative fluxes.
    standardLuminosities=max(standardLuminosities,0.0d0)
    return
  end function standardLuminosities

  subroutine standardLuminosityTracks(self,luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,abundancesStellar,redshift,ages,luminosities)
    !!{
    Returns the luminosity for a $1 M_\odot$ simple stellar population of given {\normalfont \ttfamily abundances} drawn from
    the given {\normalfont \ttfamily stellarPopulation} and observed through the filter specified by {\normalfont \ttfamily
    filterIndex}, for all available ages.
    !!}
    use            :: Abundances_Structure, only : Abundances_Get_Metallicity, logMetallicityZero, metallicityTypeLogarithmicByMassSolar
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    implicit none
    class           (stellarPopulationBroadBandLuminositiesStandard), intent(inout) :: self
    integer                                                         , intent(in   ), dimension( :   )              :: luminosityIndex                       , filterIndex
    double precision                                                , intent(in   ), dimension( :   )              :: redshift
    double precision                                                , intent(  out), dimension( :   ), allocatable :: ages
    double precision                                                , intent(  out), dimension( : ,:), allocatable :: luminosities
    type            (stellarPopulationSpectraPostprocessorList     ), intent(inout), dimension( :   )              :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                        ), intent(inout)                                :: stellarPopulation_
    type            (abundances                                    ), intent(in   )                                :: abundancesStellar
    double precision                                                               , dimension(0:1  )              :: hMetallicity
    integer         (c_size_t                                      )                                               :: iLuminosity                           , iMetallicity, &
         &                                                                                                            jMetallicity                          , populationID
    double precision                                                                                               :: metallicity

    ! Obtain a read lock on the luminosity tables.
    call self%luminosityTableLock%setRead()
    ! Tabulate the luminosities.
    call self%tabulate(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    ! Get interpolation in metallicity.
    populationID=stellarPopulation_%uniqueID()
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < self%luminosityTables(populationID)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > self%luminosityTables(populationID)%metallicity(self%luminosityTables(populationID)%metallicitiesCount)) then
       iMetallicity=self%luminosityTables(populationID)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       call self%luminosityTables(populationID)%interpolatorMetallicity%linearFactors(metallicity,iMetallicity,hMetallicity)
    end if
    ! Allocate arrays for ages and luminosities.
    allocate(ages        (self%luminosityTables(populationID)%agesCount                      ))
    allocate(luminosities(self%luminosityTables(populationID)%agesCount,size(luminosityIndex)))
    ! Assign ages.
    ages=self%luminosityTables(populationID)%age
    ! Do the interpolation.
    luminosities(:,:)=0.0d0
    do iLuminosity=1,size(luminosityIndex)
       do jMetallicity=0,1
          luminosities                                          (:,                iLuminosity                             )= &
               & +luminosities                                  (:,                iLuminosity                             )  &
               & +self%luminosityTables(populationID)%luminosity(  luminosityIndex(iLuminosity),:,iMetallicity+jMetallicity)  &
               & *hMetallicity                                  (                                              jMetallicity)
       end do
    end do
    ! Release the read lock on the luminosity tables.
    call self%luminosityTableLock%setRead()
    ! Prevent interpolation from returning negative fluxes.
    luminosities=max(luminosities,0.0d0)
    return
  end subroutine standardLuminosityTracks

  subroutine standardTabulate(self,luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    !!{
    Tabulate stellar population luminosity in the given filters.
    !!}
    use            :: Abundances_Structure            , only : logMetallicityZero          , metallicityTypeLogarithmicByMassSolar
    use            :: Display                         , only : displayCounter              , displayCounterClear                  , displayGreen            , displayIndent        , &
          &                                                    displayMagenta              , displayReset                         , displayUnindent         , verbosityLevelWorking
    use            :: File_Utilities                  , only : File_Exists                 , File_Lock                            , File_Unlock             , lockDescriptor       , &
         &                                                     Directory_Make                  , File_Path
    use            :: Error                           , only : Error_Report                , Warn                                 , errorStatusFail         , errorStatusSuccess
    use            :: HDF5_Access                     , only : hdf5Access
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: ISO_Varying_String              , only : assignment(=)               , char                                 , operator(//)            , var_str
    use            :: Input_Parameters                , only : inputParameters
    use            :: Instruments_Filters             , only : Filter_Extent               , Filter_Name                          , Filter_Response_Function
    use            :: Numerical_Constants_Astronomical, only : metallicitySolar
    use            :: Numerical_Integration           , only : GSL_Integ_Gauss15           , integrator
    use            :: String_Handling                 , only : operator(//)                , stringXMLFormat
    use            :: Table_Labels                    , only : extrapolationTypeExtrapolate
    implicit none
    class           (stellarPopulationBroadBandLuminositiesStandard), intent(inout)                   :: self
    integer                                                         , intent(in   ), dimension(:    ) :: filterIndex                                   , luminosityIndex
    double precision                                                , intent(in   ), dimension(:    ) :: redshift
    type            (stellarPopulationSpectraPostprocessorList     ), intent(inout), dimension(:    ) :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                        ), intent(inout)                   :: stellarPopulation_
    type            (luminosityTable                               ), allocatable  , dimension(:    ) :: luminosityTablesTemporary
    double precision                                                , allocatable  , dimension(:,:,:) :: luminosityTemporary
    logical                                                         , allocatable  , dimension(:    ) :: isTabulatedTemporary
    double precision                                                               , dimension(2    ) :: wavelengthRange
    type            (lockDescriptor                                )                                  :: lockFileDescriptor
    class           (stellarPopulationSpectraClass                 ), pointer                         :: stellarPopulationSpectra_
    class           (stellarPopulationSpectraPostprocessorClass    ), pointer                         :: stellarPopulationSpectraPostprocessorPrevious_
    type            (integrator                                    ), allocatable                     :: integratorAB_
    type            (integrator                                    ), allocatable  , save             :: integrator_
    !$omp threadprivate(integrator_)
    type            (inputParameters                               )               , save             :: descriptor
    !$omp threadprivate(descriptor)
    integer         (c_size_t                                      )                                  :: iAge                                          , iLuminosity                          , &
         &                                                                                               iMetallicity                                  , jLuminosity                          , &
         &                                                                                               populationID
    integer                                                                                           :: loopCountMaximum                              , loopCount                            , &
         &                                                                                               errorStatus                                   , luminosityIndexMaximum
    logical                                                                                           :: computeTable                                  , calculateLuminosity                  , &
         &                                                                                               stellarPopulationHashedDescriptorComputed     , copyDone
    double precision                                                                                  :: toleranceRelative                             , normalization
    type            (varying_string                                )               , save             :: message                                       , luminositiesFileName                 , &
         &                                                                                               descriptorString                              , stellarPopulationHashedDescriptor    , &
         &                                                                                               postprocessorHashedDescriptor
    !$omp threadprivate(message,luminositiesFileName,descriptorString,stellarPopulationHashedDescriptor,postprocessorHashedDescriptor)
    character       (len=16                                        )                                  :: datasetName                                   , redshiftLabel                        , &
         &                                                                                               label
    
    ! Obtain a read lock on the luminosity tables.
    call self%luminosityTableLock%setRead()
    ! Allocate table storage. First test if the tables must be resized (or allocated).
    populationID=stellarPopulation_%uniqueID()
    if (.not.allocated(self%luminosityTables).or.size(self%luminosityTables) < populationID) then
       ! A resize or allocation is required. Obtain a write lock on the tables and then re-test if resizing is required (as
       ! another thread could potentially have already done this while we waited for the write lock).
       call self%luminosityTableLock%setWrite(haveReadLock=.true.)
       if (allocated(self%luminosityTables)) then
          if (size(self%luminosityTables) < populationID) then
             call Move_Alloc(self%luminosityTables,luminosityTablesTemporary)
             allocate(self%luminosityTables(populationID))
             !![
	     <workaround type="gfortran" PR="46897" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=46897">
	       <seeAlso type="gfortran" PR="57696" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=57696"/>
	       <description>
		 Type-bound defined assignment not done because multiple part array references would occur in intermediate expressions.
	       </description>
	     !!]
             do iLuminosity=1,size(luminosityTablesTemporary)
                self%luminosityTables(iLuminosity)=luminosityTablesTemporary(iLuminosity)
             end do
             !![
	     </workaround>
             !!]
             self%luminosityTables(size(luminosityTablesTemporary)+1:populationID)%isTabulatedMaximum=0
             deallocate(luminosityTablesTemporary)
          end if
       else
          allocate(self%luminosityTables(populationID))
          self%luminosityTables%isTabulatedMaximum=0
       end if
       call self%luminosityTableLock%unsetWrite(haveReadLock=.true.)
    end if
    ! Determine if we have tabulated luminosities for this luminosityIndex in this population yet.
    luminosityIndexMaximum=maxval(luminosityIndex)
    if (.not.allocated(self%luminosityTables(populationID)%isTabulated) .or. self%luminosityTables(populationID)%isTabulatedMaximum < luminosityIndexMaximum) then
       ! Tabulation is required. Obtain a write lock on the luminosity tables and then retest if tabulation is required (as
       ! another thread may have tabulated while we waited for the write lock).
       call self%luminosityTableLock%setWrite(haveReadLock=.true.)
       luminosityIndexMaximum=maxval(luminosityIndex)
       if (.not.allocated(self%luminosityTables(populationID)%isTabulated) .or. self%luminosityTables(populationID)%isTabulatedMaximum < luminosityIndexMaximum) then
          !$omp critical(broadBandLuminositiesStandardComputeTable)
          stellarPopulationSpectra_                      => stellarPopulation_%spectra()
          luminosityIndexMaximum                         =  maxval(luminosityIndex)
          stellarPopulationHashedDescriptorComputed      =  .false.
          stellarPopulationSpectraPostprocessorPrevious_ => null()
          !$omp parallel private(iAge,iMetallicity,toleranceRelative,errorStatus,copyDone)
          copyDone=.false.
          do iLuminosity=1,size(luminosityIndex)
             !$omp single
             if (allocated(self%luminosityTables(populationID)%isTabulated)) then
                if (size(self%luminosityTables(populationID)%isTabulated) >= luminosityIndex(iLuminosity)) then
                   computeTable=.not.self%luminosityTables(populationID)%isTabulated(luminosityIndex(iLuminosity))
                else
                   call move_alloc(self%luminosityTables(populationID)%isTabulated,isTabulatedTemporary                                                                                                        )
                   call move_alloc(self%luminosityTables(populationID)%luminosity ,luminosityTemporary                                                                                                         )
                   allocate       (self%luminosityTables(populationID)%isTabulated(luminosityIndexMaximum                                                                                                     ))
                   allocate       (self%luminosityTables(populationID)%luminosity (luminosityIndexMaximum,self%luminosityTables(populationID)%agesCount,self%luminosityTables(populationID)%metallicitiesCount))
                   self%luminosityTables(populationID)%isTabulated(1:size(isTabulatedTemporary)    )=isTabulatedTemporary
                   self%luminosityTables(populationID)%isTabulated(  size(isTabulatedTemporary)+1:luminosityIndexMaximum)=.false.
                   self%luminosityTables(populationID)%luminosity (1:size(isTabulatedTemporary),:,:)=luminosityTemporary
                   deallocate(isTabulatedTemporary)
                   deallocate(luminosityTemporary)
                   computeTable=.true.
                end if
             else
                allocate(self%luminosityTables(populationID)%isTabulated(luminosityIndexMaximum))
                self%luminosityTables(populationID)%isTabulated=.false.
                ! Since we have not yet tabulated any luminosities yet for this population, we need to get a list of suitable
                ! metallicities and ages at which to tabulate.
                call stellarPopulationSpectra_%tabulation(self%luminosityTables(populationID)%agesCount,self%luminosityTables(populationID)%metallicitiesCount,self%luminosityTables(populationID)%age,self%luminosityTables(populationID)%metallicity)
                where (self%luminosityTables(populationID)%metallicity > 0.0d0)
                   self%luminosityTables(populationID)%metallicity=log10(self%luminosityTables(populationID)%metallicity/metallicitySolar)
                elsewhere
                   self%luminosityTables(populationID)%metallicity=logMetallicityZero
                end where
                allocate(self%luminosityTables(populationID)%luminosity(luminosityIndexMaximum,self%luminosityTables(populationID)%agesCount,self%luminosityTables(populationID)%metallicitiesCount))
                self%luminosityTables(populationID)%interpolatorAge        =interpolator(self%luminosityTables(populationID)%age        ,extrapolationType=extrapolationTypeExtrapolate)
                self%luminosityTables(populationID)%interpolatorMetallicity=interpolator(self%luminosityTables(populationID)%metallicity,extrapolationType=extrapolationTypeExtrapolate)
                computeTable=.true.
             end if
             !$omp end single
             ! If we haven't, do so now.
             if (computeTable) then

                !$omp single
                ! Determine if we can read the required luminosity from file.
                calculateLuminosity=.true.
                if (self%storeToFile) then
                   ! Construct name of the file to which this would be stored.
                   if (.not.stellarPopulationHashedDescriptorComputed) then
                      stellarPopulationHashedDescriptor        =stellarPopulation_%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)
                      stellarPopulationHashedDescriptorComputed=.true.
                   end if
                   if (.not.associated(stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_,stellarPopulationSpectraPostprocessorPrevious_)) then
                      postprocessorHashedDescriptor                  =  stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)
                      stellarPopulationSpectraPostprocessorPrevious_ => stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_
                   end if
                   luminositiesFileName=self%storeDirectory                  // &
                        &               "/stellarLuminosities::filter:"      // &
                        &               Filter_Name(filterIndex(iLuminosity))// &
                        &               "::postprocessor:"                   // &
                        &               postprocessorHashedDescriptor        // &
                        &               "::population:"                      // &
                        &               stellarPopulationHashedDescriptor    // &
                        &               ".hdf5"
                   if (File_Exists(luminositiesFileName)) then
                      ! Construct the dataset name.
                      write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                      datasetName="redshift"//adjustl(trim(redshiftLabel))
                      ! Open the file and check for the required dataset.
                      ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
                      call File_Lock(char(luminositiesFileName),lockFileDescriptor,lockIsShared=.true.)
                      block
                        type(hdf5Object) :: luminositiesFile
                        !$ call hdf5Access%set()
                        luminositiesFile=hdf5Object(char(luminositiesFileName),readOnly=.true.)
                        if (luminositiesFile%hasDataset(trim(datasetName))) then
                           ! Read the dataset.
                           call luminositiesFile%readDatasetStatic(trim(datasetName),self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:))
                           ! We do not need to calculate this luminosity.
                           calculateLuminosity=.false.
                        end if
                        !$ call hdf5Access%unset()
                      end block
                      call File_Unlock(lockFileDescriptor)
                   end if
                end if
                !$omp end single copyprivate(luminositiesFileName,stellarPopulationHashedDescriptor,postprocessorHashedDescriptor)
                ! Compute the luminosity if necessary.
                if (calculateLuminosity) then
                   !$omp single
                   ! Display a message and counter.
                   message=var_str('Tabulating stellar luminosities for stellar population #')//populationID//', luminosity '
                   write (redshiftLabel,'(f6.3)') redshift(iLuminosity)
                   message=message                                               // &
                        &  Filter_Name                 (filterIndex(iLuminosity))// &
                        &  ":z"                                                  // &
                        &  trim(adjustl(redshiftLabel))                          // &
                        &  " "                                                   // &
                        &                                           iLuminosity  // &
                        &  " of "                                                // &
                        &  size(luminosityIndex)
                   call displayIndent (message,verbosityLevelWorking)
                   call displayCounter(0,.true.,verbosityLevelWorking)
                   ! Get wavelength extent of the filter.
                   wavelengthRange=Filter_Extent(filterIndex(iLuminosity))
                   ! Integrate over the wavelength range.
                   populationID_   =populationID
                   indexFilter_    =filterIndex (iLuminosity)
                   redshift_       =redshift    (iLuminosity)
                   loopCountMaximum=+self%luminosityTables(populationID)%metallicitiesCount &
                        &           *self%luminosityTables(populationID)%agesCount
                   loopCount       = 0
                   !$omp end single copyprivate(indexFilter_,redshift_,populationID_)
                   if (.not.copyDone) then
                      copyDone=.true.
                      allocate(integrator_)
                      integrator_=integrator(integrandFilteredLuminosity,toleranceRelative=self%integrationToleranceRelative,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=10000_c_size_t)
                      allocate(stellarPopulationSpectra__             ,mold=stellarPopulationSpectra_                                                                 )
                      !$omp critical(broadBandLuminositiesDeepCopy)
                      !![
                      <deepCopyReset variables="stellarPopulationSpectra_"/>
                      <deepCopy source="stellarPopulationSpectra_" destination="stellarPopulationSpectra__"/>
                      <deepCopyFinalize variables="stellarPopulationSpectra__"/>
                      !!]
                      !$omp end critical(broadBandLuminositiesDeepCopy)
                   end if
                   ! The postprocessor can differ for each luminosity, so deep copy it always.
                   allocate(stellarPopulationSpectraPostprocessor__,mold=stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_)
                   !$omp critical(broadBandLuminositiesDeepCopy)
                   !![
                   <deepCopyReset variables="stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_"/>
                   <deepCopy source="stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_" destination="stellarPopulationSpectraPostprocessor__"/>
                   <deepCopyFinalize variables="stellarPopulationSpectraPostprocessor__"/>
                   !!]
                   !$omp end critical(broadBandLuminositiesDeepCopy)
                   ! Get the filter response function as an interpolator.
                   filterResponse_ => Filter_Response_Function(filterIndex(iLuminosity))
                   !$omp do schedule(dynamic)
                   do iAge=1,self%luminosityTables(populationID)%agesCount
                      age_=self%luminosityTables(populationID)%age(iAge)
                      do iMetallicity=1,self%luminosityTables(populationID)%metallicitiesCount
                         ! Update the counter.
                         !$omp atomic
                         loopCount=loopCount+1
                         call displayCounter(int(100.0d0*dble(loopCount)/dble(loopCountMaximum)),.false.,verbosityLevelWorking)
                         call abundances_%metallicitySet(self%luminosityTables(populationID)%metallicity(iMetallicity) &
                              &,metallicityType=metallicityTypeLogarithmicByMassSolar)
                         toleranceRelative=self%integrationToleranceRelative
                         errorStatus      =errorStatusFail
                         do while (errorStatus /= errorStatusSuccess)
                            call integrator_%toleranceSet(toleranceRelative=toleranceRelative)
                            self%luminosityTables(populationID)%luminosity(                              &
                                 &                                         luminosityIndex(iLuminosity), &
                                 &                                         iAge                        , &
                                 &                                         iMetallicity                  &
                                 &                                        )                              &
                                 & =integrator_%integrate(                           &
                                 &                               wavelengthRange(1), &
                                 &                               wavelengthRange(2), &
                                 &                        status=errorStatus         &
                                 &                       )
                            if (errorStatus /= errorStatusSuccess) then
                               if (self%integrationToleranceDegrade.and.toleranceRelative < 1.0d0) then
                                  toleranceRelative=2.0d0*toleranceRelative
                                  write (label,'(e9.3)') 2.0d0*self%integrationToleranceRelative
                                  message=         displayMagenta()//"WARNING:"//displayReset()//" increasing relative tolerance for stellar population luminosities to"    //char(10)
                                  message=message//trim(adjustl(label))//" and retrying integral"//{introspection:location}
                                  call Warn(message)
                               else if (self%integrationToleranceDegrade) then
                                  message="integration of stellar populations failed"
                                  call Error_Report(message//{introspection:location})
                               else
                                  write (label,'(e9.3)')       self%integrationToleranceRelative
                                  message=         "integration of stellar populations failed"                                                  //char(10)
                                  message=message//displayGreen()
                                  message=message//"HELP: "
                                  message=message//displayReset()
                                  message=message//      "consider increasing the integrationtolerance parameter from the current value of "
                                  message=message//trim(adjustl(label))
                                  write (label,'(e9.3)') 2.0d0*self%integrationToleranceRelative
                                  message=message//"      to "
                                  message=message//trim(adjustl(label))
                                  message=message//      " if you can accept this lower accuracy."                                               //char(10)//char(10)
                                  message=message//"      To do this, set the highlited option in your parameter file:"                          //char(10)//char(10)
                                  message=message//stringXMLFormat('<stellarPopulationBroadBandLuminosities value="standard">**B<integrationToleranceRelative value="'//trim(adjustl(label))//'"/>**C</stellarPopulationBroadBandLuminosities>',indentInitial=6)//char(10)//char(10)
                                  message=message//"      Alternatively you can allow tolerances to be automatically degraded where"             //char(10)
                                  message=message//"      needed to ensure convergence by setting the highlighted option in your parameter file:"//char(10)//char(10)
                                  message=message//stringXMLFormat('<stellarPopulationBroadBandLuminosities value="standard">**B<integrationToleranceDegrade value="true"/>**C</stellarPopulationBroadBandLuminosities>',indentInitial=6)//char(10)
                                 call Error_Report(message//{introspection:location})
                               end if
                            end if
                         end do
                      end do
                   end do
                   !$omp end do
                   !![
                   <objectDestructor name="stellarPopulationSpectraPostprocessor__"/>
                   !!]
                   !$omp single
                   ! Clear the counter and write a completion message.
                   call displayCounterClear(           verbosityLevelWorking)
                   call displayUnindent    ('finished',verbosityLevelWorking)
                   ! Get the normalization by integrating a zeroth magnitude (AB) source through the filter.
                   if (.not.allocated(integratorAB_)) then
                      allocate(integratorAB_)
                      integratorAB_=integrator(integrandFilteredLuminosityAB,toleranceRelative=self%integrationToleranceRelative)
                   end if
                   normalization=integratorAB_%integrate(wavelengthRange(1),wavelengthRange(2))
                   ! Normalize the luminosity.
                   self         %luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:) &
                        & =+self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:) &
                        &  /normalization
                   ! Store the luminosities to file.
                   if (self%storeToFile) then
                      ! Construct the dataset name.
                      write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                      datasetName="redshift"//adjustl(trim(redshiftLabel))
                      ! Open the file.
                      descriptor=inputParameters()
                      call stellarPopulation_%descriptor(descriptor,includeClass=.true.)
                      descriptorString=descriptor%serializeToString()
                      call descriptor%destroy()
                      ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
                      call Directory_Make(File_Path(luminositiesFileName)                                        )
                      call File_Lock     (          luminositiesFileName ,lockFileDescriptor,lockIsShared=.false.)
                      block
                        type(hdf5Object) :: luminositiesFile
                        !$ call hdf5Access%set()
                        luminositiesFile=hdf5Object(luminositiesFileName)
                        if (.not.luminositiesFile%hasAttribute('parameters')) call luminositiesFile%writeAttribute(char(descriptorString),'parameters')
                        ! Write the dataset.
                        if (.not.luminositiesFile%hasDataset(trim(datasetName))) &
                             & call luminositiesFile%writeDataset(self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:),datasetName=trim(datasetName),comment="Tabulated luminosities at redshift z="//adjustl(trim(redshiftLabel)))
                        !$ call hdf5Access%unset()
                      end block
                      call File_Unlock(lockFileDescriptor)
                   end if
                   !$omp end single
                   ! Destroy our filter response function.
                   deallocate(filterResponse_)
                end if
                !$omp single
                ! Flag that calculations have been performed for this filter.
                self%luminosityTables(populationID)%isTabulated(luminosityIndex(iLuminosity))=.true.
                if (luminosityIndex(iLuminosity) > self%luminosityTables(populationID)%isTabulatedMaximum) then
                   jLuminosity=self%luminosityTables(populationID)%isTabulatedMaximum
                   do while (jLuminosity < size(self%luminosityTables(populationID)%isTabulated) .and. self%luminosityTables(populationID)%isTabulated(min(jLuminosity+1,size(self%luminosityTables(populationID)%isTabulated))))
                      jLuminosity=jLuminosity+1
                   end do
                   self%luminosityTables(populationID)%isTabulatedMaximum=jLuminosity
                end if
                !$omp end single
             end if
             !$omp barrier
          end do
          !$omp single
          if (allocated(integratorAB_)) deallocate(integratorAB_)
          !$omp end single
          if (copyDone) then
             deallocate(integrator_)
             !![
             <objectDestructor name="stellarPopulationSpectra__"/>
             !!]
          end if
          !$omp end parallel
          !$omp end critical(broadBandLuminositiesStandardComputeTable)
       end if
       call self%luminosityTableLock%unsetWrite(haveReadLock=.true.)
    end if
    ! Release our read lock on the luminosity tables.
    call self%luminosityTableLock%unsetRead()
    return

  contains

    double precision function integrandFilteredLuminosity(wavelength)
      !!{
      Integrand for the luminosity through a given filter.
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength
      double precision                :: wavelengthRedshifted

      ! If this luminosity is for a redshifted spectrum, then we shift wavelength at which we sample the stellar population
      ! spectrum to be a factor of (1+z) smaller. We therefore integrate over the stellar SED at shorter wavelengths, since these
      ! will be shifted into the filter by z=0. The factor of 1/λ appears since we want to integrate F_ν (dν / ν) and dν =
      ! -c/λ² dλ. Note that we follow the convention of Hogg et al. (2002) and assume that the filter response gives the fraction
      ! of incident photons received by the detector at a given wavelength, multiplied by the relative photon response (which will
      ! be 1 for a photon-counting detector such as a CCD, or proportional to the photon energy for a bolometer/calorimeter type
      ! detector).
      wavelengthRedshifted=wavelength/(1.0d0+redshift_)
      integrandFilteredLuminosity=+filterResponse_                        %interpolate(wavelength                             ) &
           &                      *stellarPopulationSpectra__             %luminosity (abundances_  ,age_,wavelengthRedshifted) &
           &                      *stellarPopulationSpectraPostprocessor__%multiplier (wavelengthRedshifted,age_,redshift_    ) &
           &                      /                                                    wavelength
      return
    end function integrandFilteredLuminosity

    double precision function integrandFilteredLuminosityAB(wavelength)
      !!{
      Integrand for the luminosity of a zeroth magnitude (AB) source through a given filter.
      !!}
      use :: Numerical_Constants_Astronomical, only : luminositySolar, luminosityZeroPointAB
      implicit none
      double precision, intent(in   ) :: wavelength
      ! Luminosity of a zeroth magnitude (AB) source in Solar luminosities per Hz.
      double precision, parameter     :: luminosityZeroPointABSolar=luminosityZeroPointAB/luminositySolar

      integrandFilteredLuminosityAB=+filterResponse_           %interpolate(wavelength) &
           &                        *luminosityZeroPointABSolar                         &
           &                        /wavelength
      return
    end function integrandFilteredLuminosityAB

  end subroutine standardTabulate

