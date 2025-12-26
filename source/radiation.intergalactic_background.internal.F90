!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a class for intergalactic background light which computes the background internally.
  !!}

  use :: Atomic_Cross_Sections_Ionization_Photo, only : atomicCrossSectionIonizationPhotoClass
  use :: Cosmology_Functions                   , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                  , only : cosmologyParametersClass
  use :: Intergalactic_Medium_State            , only : intergalacticMediumStateClass
  use :: Kind_Numbers                          , only : kind_int8
  use :: Numerical_Interpolation               , only : interpolator
  use :: Output_Times                          , only : outputTimesClass
  use :: Star_Formation_Rates_Disks            , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids        , only : starFormationRateSpheroidsClass
  use :: Stellar_Population_Selectors          , only : stellarPopulationSelectorClass

  !![
  <radiationField name="radiationFieldIntergalacticBackgroundInternal">
   <description>A radiation field class for intergalactic background light with properties computed internally.</description>
   <stateStore>
     <stateStore variables="accretionDiskSpectra_" store="accretionDiskSpectraStateStore_" restore="accretionDiskSpectraStateRestore_" module="Functions_Global"/>
   </stateStore>
   <deepCopy>
     <ignore variables="accretionDiskSpectra_"/>
   </deepCopy>
  </radiationField>
  !!]
  type, extends(radiationFieldIntergalacticBackground) :: radiationFieldIntergalacticBackgroundInternal
     !!{
     A radiation field class for intergalactic background light with properties computed internally
     !!}
     private
     class           (cosmologyParametersClass              ), pointer                     :: cosmologyParameters_               => null()
     class           (cosmologyFunctionsClass               ), pointer                     :: cosmologyFunctions_                => null()
     class           (intergalacticMediumStateClass         ), pointer                     :: intergalacticMediumState_          => null()
     class           (atomicCrossSectionIonizationPhotoClass), pointer                     :: atomicCrossSectionIonizationPhoto_ => null()
     class           (*                                     ), pointer                     :: accretionDiskSpectra_              => null()
     class           (starFormationRateDisksClass           ), pointer                     :: starFormationRateDisks_            => null()
     class           (starFormationRateSpheroidsClass       ), pointer                     :: starFormationRateSpheroids_        => null()
     class           (stellarPopulationSelectorClass        ), pointer                     :: stellarPopulationSelector_         => null()
     class           (outputTimesClass                      ), pointer                     :: outputTimes_                       => null()
     integer                                                                               :: wavelengthCountPerDecade
     integer         (c_size_t                              )                              :: wavelengthCount
     double precision                                                                      :: wavelengthMinimum                           , wavelengthMaximum
     integer                                                                               :: timeCountPerDecade
     integer         (c_size_t                              )                              :: timeCount
     double precision                                                                      :: redshiftMinimum                             , redshiftMaximum
     double precision                                                                      :: timeMinimum                                 , timeMaximum
     double precision                                        , allocatable, dimension(:  ) :: wavelength                                  , redshift                       , &
          &                                                                                   time_                                       , crossSectionNeutralHydrogen    , &
          &                                                                                   crossSectionNeutralHelium                   , crossSectionSinglyIonizedHelium, &
          &                                                                                   spectrum
     double precision                                        , allocatable, dimension(:,:) :: emissivityODE                               , emissivity
     double precision                                                     , dimension(0:1) :: timeODE
     double precision                                                                      :: timeCurrent
     type            (interpolator                          )                              :: interpolatorWavelength                      , interpolatorTime
     class           (*                                     ), pointer                     :: statePrevious                      => null()
     double precision                                                                      :: timePrevious
     integer         (kind_int8                             )                              :: universeUniqueIDPrevious
     integer         (c_size_t                              )                              :: iTime
     double precision                                                     , dimension(0:1) :: hTime
   contains
     final     ::                      intergalacticBackgroundInternalDestructor
     procedure :: flux              => intergalacticBackgroundInternalFlux
     procedure :: time              => intergalacticBackgroundInternalTime
     procedure :: timeSet           => intergalacticBackgroundInternalTimeSet
     procedure :: timeDependentOnly => intergalacticBackgroundInternalTimeDependentOnly
     procedure :: autoHook          => intergalacticBackgroundInternalAutoHook
     procedure :: deepCopy          => intergalacticBackgroundInternalDeepCopy
     procedure :: deepCopyReset     => intergalacticBackgroundInternalDeepCopyReset
     procedure :: deepCopyFinalize  => intergalacticBackgroundInternalDeepCopyFinalize
  end type radiationFieldIntergalacticBackgroundInternal

  interface radiationFieldIntergalacticBackgroundInternal
     !!{
     Constructors for the \refClass{radiationFieldIntergalacticBackgroundInternal} radiation field class.
     !!}
     module procedure intergalacticBackgroundInternalConstructorParameters
     module procedure intergalacticBackgroundInternalConstructorInternal
  end interface radiationFieldIntergalacticBackgroundInternal

  type :: intergalacticBackgroundInternalState
     !!{
     Class used to store the state of the intergalactic background radiation field for the internal solver. This will be stored
     as an attribute of the universe object.
     !!}
     double precision, allocatable, dimension(:,:) :: flux
     double precision                              :: timePrevious, timeNext
  end type intergalacticBackgroundInternalState

  ! Module-scope pointer to self for ODE solving.
  class(radiationFieldIntergalacticBackgroundInternal), pointer :: self_
  !$omp threadprivate(self_)

contains

  function intergalacticBackgroundInternalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiationFieldIntergalacticBackgroundInternal} radiation field class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter                , inputParameters
    use :: Functions_Global, only : accretionDiskSpectraConstruct_, accretionDiskSpectraDestruct_
    implicit none
    type            (radiationFieldIntergalacticBackgroundInternal)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (cosmologyParametersClass                     ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                      ), pointer       :: cosmologyFunctions_
    class           (intergalacticMediumStateClass                ), pointer       :: intergalacticMediumState_
    class           (atomicCrossSectionIonizationPhotoClass       ), pointer       :: atomicCrossSectionIonizationPhoto_
    class           (*                                            ), pointer       :: accretionDiskSpectra_
    class           (starFormationRateDisksClass                  ), pointer       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass              ), pointer       :: starFormationRateSpheroids_
    class           (stellarPopulationSelectorClass               ), pointer       :: stellarPopulationSelector_
    class           (outputTimesClass                             ), pointer       :: outputTimes_
    integer                                                                        :: wavelengthCountPerDecade          , timeCountPerDecade
    double precision                                                               :: wavelengthMinimum                 , wavelengthMaximum , &
         &                                                                            redshiftMinimum                   , redshiftMaximum

    !![
    <inputParameter>
      <name>wavelengthCountPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of wavelength to use for calculations of the cosmic background radiation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMinimum</name>
      <defaultValue>100.0d0</defaultValue>
      <description>The minimum wavelength (in units of \AA) to use in calculations of the cosmic background radiation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMaximum</name>
      <defaultValue>100000.0d0</defaultValue>
      <description>The maximum wavelength (in units of \AA) to use in calculations of the cosmic background radiation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeCountPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of time to use for calculations of the cosmic background radiation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum redshift to use in calculations of the cosmic background radiation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <defaultValue>30.0d0</defaultValue>
      <description>The maximum redshift to use in calculations of the cosmic background radiation.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"               name="cosmologyParameters_"               source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                name="cosmologyFunctions_"                source="parameters"/>
    <objectBuilder class="intergalacticMediumState"          name="intergalacticMediumState_"          source="parameters"/>
    <objectBuilder class="atomicCrossSectionIonizationPhoto" name="atomicCrossSectionIonizationPhoto_" source="parameters"/>
    <objectBuilder class="starFormationRateDisks"            name="starFormationRateDisks_"            source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"        name="starFormationRateSpheroids_"        source="parameters"/>
    <objectBuilder class="stellarPopulationSelector"         name="stellarPopulationSelector_"         source="parameters"/>
    <objectBuilder class="outputTimes"                       name="outputTimes_"                       source="parameters"/>
    !!]
    call accretionDiskSpectraConstruct_(parameters,accretionDiskSpectra_)
    self=radiationFieldIntergalacticBackgroundInternal(wavelengthMinimum,wavelengthMaximum,wavelengthCountPerDecade,redshiftMinimum,redshiftMaximum,timeCountPerDecade,cosmologyParameters_,cosmologyFunctions_,intergalacticMediumState_,atomicCrossSectionIonizationPhoto_,accretionDiskSpectra_,starFormationRateDisks_,starFormationRateSpheroids_,stellarPopulationSelector_,outputTimes_)
    !![
    <inputParametersValidate source="parameters" extraAllowedNames="accretionDiskSpectra"/>
    <objectDestructor name="cosmologyParameters_"              />
    <objectDestructor name="cosmologyFunctions_"               />
    <objectDestructor name="intergalacticMediumState_"         />
    <objectDestructor name="atomicCrossSectionIonizationPhoto_"/>
    <objectDestructor name="starFormationRateDisks_"           />
    <objectDestructor name="starFormationRateSpheroids_"       />
    <objectDestructor name="stellarPopulationSelector_"        />
    <objectDestructor name="outputTimes_"                      />
    !!]
    if (associated(accretionDiskSpectra_)) call accretionDiskSpectraDestruct_(accretionDiskSpectra_)
    return
  end function intergalacticBackgroundInternalConstructorParameters

  function intergalacticBackgroundInternalConstructorInternal(wavelengthMinimum,wavelengthMaximum,wavelengthCountPerDecade,redshiftMinimum,redshiftMaximum,timeCountPerDecade,cosmologyParameters_,cosmologyFunctions_,intergalacticMediumState_,atomicCrossSectionIonizationPhoto_,accretionDiskSpectra_,starFormationRateDisks_,starFormationRateSpheroids_,stellarPopulationSelector_,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{radiationFieldIntergalacticBackgroundInternal} radiation field class.
    !!}
    use :: Numerical_Ranges, only : Make_Range          , rangeTypeLogarithmic
    use :: Table_Labels    , only : extrapolationTypeFix, extrapolationTypeZero
    implicit none
    type            (radiationFieldIntergalacticBackgroundInternal)                        :: self
    integer                                                        , intent(in   )         :: wavelengthCountPerDecade          , timeCountPerDecade
    double precision                                               , intent(in   )         :: wavelengthMinimum                 , wavelengthMaximum , &
         &                                                                                    redshiftMinimum                   , redshiftMaximum
    class           (cosmologyParametersClass                     ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass                      ), intent(in   ), target :: cosmologyFunctions_
    class           (intergalacticMediumStateClass                ), intent(in   ), target :: intergalacticMediumState_
    class           (atomicCrossSectionIonizationPhotoClass       ), intent(in   ), target :: atomicCrossSectionIonizationPhoto_
    class           (*                                            ), intent(in   ), target :: accretionDiskSpectra_
    class           (starFormationRateDisksClass                  ), intent(in   ), target :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass              ), intent(in   ), target :: starFormationRateSpheroids_
    class           (stellarPopulationSelectorClass               ), intent(in   ), target :: stellarPopulationSelector_
    class           (outputTimesClass                             ), intent(in   ), target :: outputTimes_
    integer         (c_size_t                                     )                        :: iTime                             , iWavelength
    !![
    <constructorAssign variables="wavelengthMinimum, wavelengthMaximum, wavelengthCountPerDecade, redshiftMinimum, redshiftMaximum, timeCountPerDecade, *cosmologyParameters_, *cosmologyFunctions_, *intergalacticMediumState_, *atomicCrossSectionIonizationPhoto_, *accretionDiskSpectra_, *starFormationRateDisks_, *starFormationRateSpheroids_, *stellarPopulationSelector_, *outputTimes_"/>
    !!]

    ! Build tables of wavelength and time for cosmic background radiation.
    self%timeMaximum=self%cosmologyFunctions_%cosmicTime                 (                      &
         &           self%cosmologyFunctions_%expansionFactorFromRedshift (                     &
         &                                                                 self%redshiftMinimum &
         &                                                                )                     &
         &                                                               )
    self%timeMinimum=self%cosmologyFunctions_%cosmicTime                 (                      &
         &           self%cosmologyFunctions_%expansionFactorFromRedshift (                     &
         &                                                                 self%redshiftMaximum &
         &                                                                )                     &
         &                                                               )
    self%wavelengthCount=+int(                                     &
         &                    +dble(self%wavelengthCountPerDecade) &
         &                    *log10(                              &
         &                           +self%wavelengthMaximum       &
         &                           /self%wavelengthMinimum       &
         &                          )                              &
         &                   )                                     &
         &               +1
    self%timeCount      =+int(                                     &
         &                    +dble(self%timeCountPerDecade      ) &
         &                    *log10(                              &
         &                           +self%timeMaximum             &
         &                           /self%timeMinimum             &
         &                          )                              &
         &                   )                                     &
         &               +1
    allocate(self%wavelength   (self%wavelengthCount                 ))
    allocate(self%spectrum     (self%wavelengthCount                 ))
    allocate(self%time_        (                       self%timeCount))
    allocate(self%redshift     (                       self%timeCount))
    allocate(self%emissivity   (self%wavelengthCount,  self%timeCount))
    allocate(self%emissivityODE(self%wavelengthCount,0:2             ))
    self%wavelength=Make_Range(                             &
         &                         self%wavelengthMinimum , &
         &                         self%wavelengthMaximum , &
         &                     int(self%wavelengthCount  ), &
         &                     rangeTypeLogarithmic         &
         &                    )
    self%time_     =Make_Range(                             &
         &                         self%timeMinimum       , &
         &                         self%timeMaximum       , &
         &                     int(self%timeCount        ), &
         &                     rangeTypeLogarithmic         &
         &                    )
    ! Convert times to redshifts.
    do iTime=1,self%timeCount
       self%redshift(iTime)                                                            &
            & =self%cosmologyFunctions_%redshiftFromExpansionFactor(                   &
            &  self%cosmologyFunctions_%expansionFactor             (                  &
            &                                                        self%time_(iTime) &
            &                                                       )                  &
            &                                                      )
    end do
    ! Initialize the background radiation to zero.
    self%spectrum  =0.0d0
    ! Initialize the emissivity to zero.
    self%emissivity=0.0d0
    ! Construct tables of photoionization cross-sections.
    allocate(self%crossSectionNeutralHydrogen    (self%wavelengthCount))
    allocate(self%crossSectionNeutralHelium      (self%wavelengthCount))
    allocate(self%crossSectionSinglyIonizedHelium(self%wavelengthCount))
    do iWavelength=1,self%wavelengthCount
       self%crossSectionNeutralHydrogen    (iWavelength)=self%atomicCrossSectionIonizationPhoto_%crossSection(1,1,1,self%wavelength(iWavelength))
       self%crossSectionNeutralHelium      (iWavelength)=self%atomicCrossSectionIonizationPhoto_%crossSection(2,1,1,self%wavelength(iWavelength))
       self%crossSectionSinglyIonizedHelium(iWavelength)=self%atomicCrossSectionIonizationPhoto_%crossSection(2,2,1,self%wavelength(iWavelength))
    end do
    ! Build interpolators.
    self%interpolatorWavelength=interpolator(self%wavelength,extrapolationType=extrapolationTypeZero)
    self%interpolatorTime      =interpolator(self%time_     ,extrapolationType=extrapolationTypeFix )
    ! Initialize state.
    self%statePrevious            => null()
    self%timePrevious             =  -1.0d0
    self%universeUniqueIDPrevious =  -1_kind_int8
    return
  end function intergalacticBackgroundInternalConstructorInternal

  subroutine intergalacticBackgroundInternalAutoHook(self)
    use :: Events_Hooks, only : universePreEvolveEventGlobal
    implicit none
    class(radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self

    ! Hook to universe pre-evolve events.
    !$omp masked
    call universePreEvolveEventGlobal%attach(self,intergalacticBackgroundInternalUniversePreEvolve,label='radiationFieldIntergalacticBackgroundInternal')
    !$omp end masked
    return
  end subroutine intergalacticBackgroundInternalAutoHook

  subroutine intergalacticBackgroundInternalDestructor(self)
    !!{
    Destructor for the \refClass{radiationFieldIntergalacticBackgroundInternal} radiation field class.
    !!}
    use :: Events_Hooks    , only : universePreEvolveEventGlobal
    use :: Functions_Global, only : accretionDiskSpectraDestruct_
    implicit none
    type(radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"              />
    <objectDestructor name="self%cosmologyFunctions_"               />
    <objectDestructor name="self%intergalacticMediumState_"         />
    <objectDestructor name="self%atomicCrossSectionIonizationPhoto_"/>
    <objectDestructor name="self%starFormationRateDisks_"           />
    <objectDestructor name="self%starFormationRateSpheroids_"       />
    <objectDestructor name="self%stellarPopulationSelector_"        />
    <objectDestructor name="self%outputTimes_"                      />
    !!]
    if (associated(self%accretionDiskSpectra_)) call accretionDiskSpectraDestruct_(self%accretionDiskSpectra_)
    !$omp masked
    if (universePreEvolveEventGlobal%isAttached(self,intergalacticBackgroundInternalUniversePreEvolve)) call universePreEvolveEventGlobal%detach(self,intergalacticBackgroundInternalUniversePreEvolve)
    !$omp end masked
    return
  end subroutine intergalacticBackgroundInternalDestructor

  double precision function intergalacticBackgroundInternalTime(self)
    !!{
    Return the epoch.
    !!}
    implicit none
    class(radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self

    intergalacticBackgroundInternalTime=self%timeCurrent
    return
  end function intergalacticBackgroundInternalTime

  subroutine intergalacticBackgroundInternalTimeSet(self,time)
    !!{
    Set the epoch.
    !!}
    implicit none
    class           (radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self
    double precision                                               , intent(in   ) :: time

    self%timeCurrent=time
    return
  end subroutine intergalacticBackgroundInternalTimeSet

 double precision function intergalacticBackgroundInternalFlux(self,wavelength,node)
    !!{
    Return the flux in the internally-computed intergalactic background.
    !!}
    use            :: Error        , only : Error_Report
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (radiationFieldIntergalacticBackgroundInternal), intent(inout)  :: self
    double precision                                               , intent(in   )  :: wavelength
    type            (treeNode                                     ), intent(inout)  :: node
    double precision                                               , dimension(0:1) :: hWavelength
    double precision                                               , parameter      :: timeTolerance=1.0d-3
    logical                                                                         :: lockHeld
    character       (len=16                                       )                 :: timeCurrent         , timeNext
    integer         (c_size_t                                     )                 :: iWavelength         , jWavelength, &
         &                                                                             jTime

    if (node%hostTree%hostUniverse%uniqueID /= self%universeUniqueIDPrevious .or. self%timeCurrent /= self%timePrevious) then
       lockHeld=node%hostTree%hostUniverse%lock%setNonBlocking()
       ! Get the state of the radiation field.
       self%universeUniqueIDPrevious =  node%hostTree%hostUniverse%uniqueID
       self%timePrevious             =  self                      %timeCurrent
       self%statePrevious            => node%hostTree%hostUniverse%attributes %value('radiationFieldIntergalacticBackgroundInternal')
       select type (state => self%statePrevious)
       type is (intergalacticBackgroundInternalState)
          ! Check that the time is within the applicable range.
          if (self%timeCurrent > state%timeNext*(1.0d0+timeTolerance)) then
             write (timeCurrent,'(e16.8)') self %timeCurrent
             write (timeNext   ,'(e16.8)') state%timeNext
             call Error_Report(                                                                       &
                  &            'time is out of range for intergalactic radiation field: '//char(10)// &
                  &            '   timeCurrent = Gyr'//adjustl(trim(timeCurrent))        //char(10)// &
                  &            '  >'                                                     //char(10)// &
                  &            '   timeNext    = Gyr'//adjustl(trim(timeNext   ))        //char(10)// &
                  &            {introspection:location}                                               &
                  &           )
          end if
          ! Find interpolation in array of times.
          call self%interpolatorTime%linearFactors(self%timeCurrent,self%iTime,self%hTime)
          if (self%timeCurrent > state%timePrevious) self%hTime=[1.0d0,0.0d0]
       class default
          call Error_Report('state has unknown type'//{introspection:location})
       end select
       if (lockHeld) call node%hostTree%hostUniverse%lock%unset()
    end if
    ! Find interpolation in the array of wavelengths.
    call self%interpolatorWavelength%linearFactors(     wavelength ,iWavelength,hWavelength)
    ! Interpolate in wavelength and time.
    intergalacticBackgroundInternalFlux=0.0d0
    select type (state => self%statePrevious)
    type is (intergalacticBackgroundInternalState)
       do jTime=0,1
          do jWavelength=0,1
             intergalacticBackgroundInternalFlux=+intergalacticBackgroundInternalFlux                                           &
                  &                              +     hWavelength                   (jWavelength                             ) &
                  &                              *self%hTime                         (                        jTime           ) &
                  &                              *state%flux                         (jWavelength+iWavelength,jTime+self%iTime)
          end do
       end do
       intergalacticBackgroundInternalFlux=max(intergalacticBackgroundInternalFlux,0.0d0)
    class default
       call Error_Report('state has unknown type'//{introspection:location})
    end select
    return
  end function intergalacticBackgroundInternalFlux

  subroutine intergalacticBackgroundInternalUniversePreEvolve(self,universe_)
    !!{
    Attach an initial event to the universe to cause the background radiation update function to be called.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : universe    , universeEvent
    implicit none
    class  (*                                   ), intent(inout), target :: self
    type   (universe                            ), intent(inout)         :: universe_
    type   (universeEvent                       ), pointer               :: event
    type   (intergalacticBackgroundInternalState), pointer               :: state
    logical                                                              :: lockHeld
    
    select type (self)
    class is (radiationFieldIntergalacticBackgroundInternal)
       ! If the universe object already has an "radiationFieldIntergalacticBackgroundInternal" attribute, then do not add a new
       ! event here - we want only one event per universe.
       if (.not.universe_%attributes%exists('radiationFieldIntergalacticBackgroundInternal')) then
          ! Create the first interrupt event in the universe object.
          event                       => universe_%createEvent( )
          event%time                  =  self     %time_      (1)
          event%creator               => self
          event%task                  => intergalacticBackgroundInternalUpdate
          lockHeld=universe_%lock%setNonBlocking()
          allocate(state                                          )
          allocate(state%flux(self%wavelengthCount,self%timeCount))
          state%timeNext    =self%time_(1)
          state%timePrevious=0.0d0
          state%flux        =0.0d0
          call universe_%attributes%set('radiationFieldIntergalacticBackgroundInternal',state)
          if (lockHeld) call universe_%lock%unset()
          self%universeUniqueIDPrevious =  -1_kind_int8
          self%timePrevious             =  -1.0d0
          self%statePrevious            => null()
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine intergalacticBackgroundInternalUniversePreEvolve

  logical function intergalacticBackgroundInternalUpdate(event,universe_) result (success)
    !!{
    Update the radiation background for a given universe.
    !!}
    use            :: Abundances_Structure        , only : abundances                       , max
    use            :: Arrays_Search               , only : searchArrayClosest
    use            :: Display                     , only : displayIndent                    , displayMessage          , displayUnindent
    use            :: Error                       , only : Error_Report
    use            :: Functions_Global            , only : accretionDiskSpectraSpectrumNode_
    use            :: Output_HDF5                 , only : outputFile
    use            :: Galacticus_Nodes            , only : defaultDiskComponent             , defaultSpheroidComponent, mergerTreeList , nodeComponentBasic, &
          &                                                nodeComponentDisk                , nodeComponentSpheroid   , treeNode       , universe          , &
          &                                                universeEvent
    use            :: HDF5_Access                 , only : hdf5Access
    use            :: IO_HDF5                     , only : hdf5Object
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: ISO_Varying_String          , only : varying_string
    use            :: Merger_Tree_Walkers         , only : mergerTreeWalkerAllNodes
    use            :: Numerical_Constants_Math    , only : Pi
    use            :: Numerical_Constants_Physical, only : plancksConstant                  , speedLight
    use            :: Numerical_Constants_Prefixes, only : centi
    use            :: Numerical_Constants_Units   , only : metersToAngstroms                , ergs
    use            :: Numerical_Integration       , only : integrator
    use            :: Numerical_ODE_Solvers       , only : odeSolver
    use            :: Stellar_Population_Spectra  , only : stellarPopulationSpectraClass
    use            :: Stellar_Populations         , only : stellarPopulationClass
    implicit none
    class           (universeEvent                       ), intent(inout) :: event
    type            (universe                            ), intent(inout) :: universe_
    type            (mergerTreeList                      ), pointer       :: forest
    type            (treeNode                            ), pointer       :: node
    class           (nodeComponentBasic                  ), pointer       :: basic
    class           (nodeComponentDisk                   ), pointer       :: disk
    class           (nodeComponentSpheroid               ), pointer       :: spheroid
    type            (universeEvent                       ), pointer       :: eventNew
    class           (stellarPopulationClass              ), pointer       :: stellarPopulationDisk_               , stellarPopulationSpheroid_
    class           (stellarPopulationSpectraClass       ), pointer       :: stellarPopulationSpectraDisk_        , stellarPopulationSpectraSpheroid_       , &
         &                                                                   stellarPopulationSpectra_
    double precision                                      , parameter     :: odeToleranceAbsolute         =1.0d-30, odeToleranceRelative             =1.0d-3
    double precision                                      , parameter     :: integrationToleranceAbsolute =1.0d-30, integrationToleranceRelative     =1.0d-3
    double precision                                      , parameter     :: timeTolerance                =1.0d-03
    type            (abundances                          ), target        :: gasAbundancesDisk                    , gasAbundancesSpheroid
    type            (abundances                          ), pointer       :: gasAbundances
    class           (*                                   ), pointer       :: state
    type            (odeSolver                           )                :: solver
    type            (integrator                          )                :: integrator_
    type            (mergerTreeWalkerAllNodes            )                :: treeWalker
    double precision                                                      :: starFormationRateDisk                , starFormationRateSpheroid               , &
         &                                                                   gasMassDisk                          , gasMassSpheroid                         , &
         &                                                                   ageEnd                               , ageStart                                , &
         &                                                                   stellarSpectrumDisk                  , stellarSpectrumSpheroid                 , &
         &                                                                   timeStart                            , timeEnd                                 , &
         &                                                                   treeTimeLatest                       , wavelength
    type            (varying_string                      )                :: message
    character       (len=6                               )                :: label
    type            (hdf5Object                          )                :: outputGroup                          , outputDataset
    integer         (c_size_t                            )                :: iTime                                , iWavelength                             , &
         &                                                                   iNow
    logical                                                               :: firstTime                            , lockHeld

    ! Guard on event creator class.
    select type (self => event%creator)
    class is (radiationFieldIntergalacticBackgroundInternal)
       ! Display message.
       write (label,'(f6.3)') event%time
       message="Evolving cosmic background radiation to time "//trim(label)//" Gyr"
       call displayIndent(message)
       ! Find the current timestep.
       iNow=searchArrayClosest(self%time_,event%time)
       ! Construct an integrator.
       integrator_=integrator(stellarSpectraConvolution,toleranceAbsolute=integrationToleranceAbsolute,toleranceRelative=integrationToleranceRelative)
       ! Iterate over all nodes.
       call displayMessage('Accumulating emissivity')
       treeTimeLatest=0.0d0
       forest => universe_%trees
       do while (associated(forest))
          treeWalker=mergerTreeWalkerAllNodes(forest%tree,spanForest=.true.)
          do while (treeWalker%next(node))
             basic          =>                    node %basic()
             treeTimeLatest =  max(treeTimeLatest,basic%time ())
             if (basic%time() == event%time) then
                ! Get the star formation rates and metallicities for this node.
                disk                      => node    %disk                            (    )
                spheroid                  => node    %spheroid                        (    )
                starFormationRateDisk     =  self    %starFormationRateDisks_    %rate(node)
                starFormationRateSpheroid =  self    %starFormationRateSpheroids_%rate(node)
                gasMassDisk               =  disk    %massGas                         (    )
                gasMassSpheroid           =  spheroid%massGas                         (    )
                gasAbundancesDisk         =  disk    %abundancesGas                   (    )
                gasAbundancesSpheroid     =  spheroid%abundancesGas                   (    )
                if (starFormationRateDisk     > 0.0d0) gasAbundancesDisk    =gasAbundancesDisk    /gasMassDisk
                if (starFormationRateSpheroid > 0.0d0) gasAbundancesSpheroid=gasAbundancesSpheroid/gasMassSpheroid
                if (starFormationRateDisk > 0.0d0 .or. starFormationRateSpheroid > 0.0d0) then
                   ! Find stellar spectra for disk and spheroid.
                   stellarPopulationDisk_            => self%stellarPopulationSelector_%select (starFormationRateDisk    ,gasAbundancesDisk    ,defaultDiskComponent    )
                   stellarPopulationSpheroid_        => self%stellarPopulationSelector_%select (starFormationRateSpheroid,gasAbundancesSpheroid,defaultSpheroidComponent)
                   stellarPopulationSpectraDisk_     =>      stellarPopulationDisk_    %spectra(                                                                        )
                   stellarPopulationSpectraSpheroid_ =>      stellarPopulationSpheroid_%spectra(                                                                        )
                   ! Find the duration of the current timestep.
                   ! Accumulate emissivity to each timestep.
                   firstTime=.true.
                   do iTime=1,self%timeCount
                      ! Skip times in the past.
                      if (self%time_(iTime) < event%time) cycle
                      ! Compute age of the currently forming population at this time.
                      ageEnd=self%time_(iTime)-event%time
                      if (iTime == 1) then
                         ageStart=                                  0.0d0
                      else
                         ageStart=max(self%time_(iTime-1)-event%time,0.0d0)
                      end if
                      ! Iterate over wavelength
                      do iWavelength=1,self%wavelengthCount
                         wavelength                =  self%wavelength(iWavelength)
                         stellarPopulationSpectra_ => stellarPopulationSpectraDisk_
                         gasAbundances             => gasAbundancesDisk
                         stellarSpectrumDisk       =  integrator_%integrate(ageStart,ageEnd)
                         stellarPopulationSpectra_ => stellarPopulationSpectraSpheroid_
                         gasAbundances             => gasAbundancesSpheroid
                         stellarSpectrumSpheroid   =  integrator_%integrate(ageStart,ageEnd)
                         self%emissivity         (iWavelength,iTime) &
                              & =+self%emissivity(iWavelength,iTime) &
                              &  +(                                  &
                              &    +stellarSpectrumDisk              &
                              &    *starFormationRateDisk            &
                              &    +stellarSpectrumSpheroid          &
                              &    *starFormationRateSpheroid        &
                              &   )                                  &
                              &  *node%hostTree%volumeWeight
                         ! Add AGN emission. This accumulates only to the the current time.
                         if (firstTime)                                                                                                                                  &
                              & self%emissivity  (iWavelength,iTime)=+self         %emissivity                       (                                iWavelength,iTime) &
                              &                                      +              accretionDiskSpectraSpectrumNode_(self%accretionDiskSpectra_,node, wavelength      ) &
                              &                                      *node%hostTree%volumeWeight
                      end do
                      firstTime=.false.
                   end do
                end if
             end if
          end do
          forest => forest%next
       end do
       ! Evolve the cosmic background radiation up to this timestep.
       state => universe_%attributes%value('radiationFieldIntergalacticBackgroundInternal')
       select type (state)
       type is (intergalacticBackgroundInternalState)
          if (iNow > 1) then
             call displayMessage('Solving cosmic background radiation evolution')
             self%timeODE      (  0:1)=self%time_     (  iNow-1:iNow)
             self%emissivityODE(:,0  )=self%emissivity(:,iNow-1     )
             self%emissivityODE(:,  1)=self%emissivity(:,       iNow)
             self_     => self
             timeStart =  self%time_(iNow-1)
             timeEnd   =  self%time_(iNow  )
             solver    =  odeSolver(self%wavelengthCount,intergalacticBackgroundInternalODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
             call solver%solve(timeStart,timeEnd,self%spectrum)
             ! Convert.
             lockHeld=universe_%lock%setNonBlocking()
             state%flux(:,iNow)=max(                                   &
                  &                 +plancksConstant                   &
                  &                 *speedLight                   **2  &
                  &                 *metersToAngstroms                 &
                  &                 /4.0d0                             &
                  &                 /Pi                                &
                  &                 /self%wavelength                   &
                  &                 *self%spectrum                     &
                  &                 *centi                        **2  &
                  &                 /ergs                            , &
                  &                 +0.0d0                             &
                  &                )
             if (lockHeld) call universe_%lock%unset()
          end if
          ! Add the next event to the universe.
          lockHeld=universe_%lock%setNonBlocking()
          state    %timePrevious=self%time_(iNow  )
          if (iNow < self%timeCount) then
             state%timeNext     =self%time_(iNow+1)
          else
             ! At the final step, set the next time to the current time plus a small tolerance to allow interpolation of the
             ! radiation field to slightly later times (which may occur due to numerical inaccuracies).
             state%timeNext     =self%time_(iNow  )*(1.0d0+timeTolerance)
          end if
          if (lockHeld) call universe_%lock%unset()
          self%universeUniqueIDPrevious =  -1_kind_int8
          self%timePrevious             =  -1.0d0
          self%statePrevious            => null()
          if     (                                                                            &
               &              iNow    <                        self             %timeCount    &
               &  .and.                                                                       &
               &   self%time_(iNow+1) < treeTimeLatest                                        &
               &  .and.                                                                       &
               &   self%time_(iNow+1) < self%outputTimes_%time(self%outputTimes_%count    ()) &
               & ) then
             eventNew         => universe_%createEvent(      )
             eventNew%time    =  self     %time_      (iNow+1)
             eventNew%creator => event    %creator
             eventNew%task    => intergalacticBackgroundInternalUpdate
          else
             ! Output the results to file.
             !$ call hdf5Access%set()
             outputGroup=outputFile%openGroup('backgroundRadiation','Cosmic background radiation data.')
             call    outputGroup  %writeDataset  (self%wavelength        ,'wavelength','Wavelength at which the background radiation is tabulated [Å].'    ,datasetReturned=outputDataset)
             call    outputDataset%writeAttribute(1.0d0/metersToAngstroms,'unitsInSI'                                                                                                    )
             call    outputGroup  %writeDataset  (self%redshift          ,'redshift'  ,'Redshift at which the background radiation is tabulated [].'       ,datasetReturned=outputDataset)
             call    outputDataset%writeAttribute(0.0d0                  ,'unitsInSI'                                                                                                    )
             call    outputGroup  %writeDataset  (state%flux             ,'flux'      ,'Flux is the cosmic background radiation [erg cm⁻² s⁻¹ Hz⁻¹ sr⁻¹].' ,datasetReturned=outputDataset)
             call    outputDataset%writeAttribute(ergs/centi**2          ,'unitsInSI'                                                                                                    )
             !$ call hdf5Access   %unset         (                                                                                                                                       )
          end if
       end select
       ! Display message.
       call displayUnindent('done')
       ! Return true since we've performed our task.
       success=.true.
    class default
       ! Incorrect event creator type.
       success=.false.
       call Error_Report('incorrect event creator class'//{introspection:location})
    end select
    return

  contains

    double precision function stellarSpectraConvolution(age)
      !!{
      Integrand for convolution of stellar spectra.
      !!}
      use :: Error, only : errorStatusInputDomain, errorStatusSuccess
      implicit none
      double precision, intent(in   ) :: age
      integer                         :: status

      stellarSpectraConvolution=stellarPopulationSpectra_%luminosity(               &
           &                                                         gasAbundances, &
           &                                                         age          , &
           &                                                         wavelength   , &
           &                                                         status         &
           &                                                        )
      if     (                                                                   &
           &   status /= errorStatusSuccess                                      &
           &  .and.                                                              &
           &   status /= errorStatusInputDomain                                  &
           & ) call Error_Report(                                                &
           &                     'stellar population spectrum function failed'// &
           &                     {introspection:location}                        &
           &                    )
      return
    end function stellarSpectraConvolution

  end function intergalacticBackgroundInternalUpdate

  integer function intergalacticBackgroundInternalODEs(time,spectrum,spectrumRateOfChange)
    !!{
    Evaluates the ODEs controlling the evolution of cosmic background radiation.
    !!}
    use :: Interface_GSL                   , only : GSL_Success
    use :: Numerical_Constants_Astronomical, only : gigaYear         , heliumByMassPrimordial, hydrogenByMassPrimordial, luminositySolar, &
          &                                         massSolar        , megaParsec
    use :: Numerical_Constants_Atomic      , only : atomicMassHelium , atomicMassHydrogen    , atomicMassUnit
    use :: Numerical_Constants_Physical    , only : plancksConstant  , speedLight
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Numerical_Constants_Units       , only : metersToAngstroms
    implicit none
    double precision, intent(in   )               :: time
    double precision, intent(in   ), dimension(:) :: spectrum
    double precision, intent(  out), dimension(:) :: spectrumRateOfChange
    double precision                              :: spectralGradient    (self_%wavelengthCount)
    double precision                              :: expansionFactor

    ! Get the expansion factor.
    expansionFactor=self_%cosmologyFunctions_%expansionFactor(time)
    ! Add source terms, linearly interpolating between timesteps. Convert from emissivity units [L☉ Hz⁻¹ Mpc⁻³] to
    ! background units [photons m⁻³ Hz⁻¹].
    spectrumRateOfChange(1:self_%wavelengthCount)                  &
         & =+(                                                     &
         &    +                          self_%emissivityODE(:,0)  &
         &    +(self_%emissivityODE(:,1)-self_%emissivityODE(:,0)) &
         &    *(            time        -self_%      timeODE(  0)) &
         &    /(self_%      timeODE(  1)-self_%      timeODE(  0)) &
         &   )                                                     &
         &  *self_%wavelength                                      &
         &  *luminositySolar                                       &
         &  *gigaYear                                              &
         &  /metersToAngstroms                                     &
         &  /plancksConstant                                       &
         &  /speedLight                                            &
         &  /megaParsec**3
    ! Add expansion dilution: -3 H(t)      n_ν
    spectrumRateOfChange(1:self_%wavelengthCount)=+spectrumRateOfChange(1:self_%wavelengthCount)            &
         &                                        -spectrum            (1:self_%wavelengthCount)            &
         &                                        *3.0d0                                                    &
         &                                        *self_%cosmologyFunctions_%expansionRate(expansionFactor)
    ! Add redshifting       : +  H(t) d(ν n_ν)/dν
    spectralGradient    (1                      )=-  spectrum          (1                      )
    spectralGradient    (2:self_%wavelengthCount)=-                                                self_%wavelength(2:self_%wavelengthCount  )**2 &
         &                                        *(                                                                                              &
         &                                          +spectrum          (2:self_%wavelengthCount  )/self_%wavelength(2:self_%wavelengthCount  )    &
         &                                          -spectrum          (1:self_%wavelengthCount-1)/self_%wavelength(1:self_%wavelengthCount-1)    &
         &                                         )                                                                                              &
         &                                        /(                                                                                              &
         &                                          +                                              self_%wavelength(2:self_%wavelengthCount  )    &
         &                                          -                                              self_%wavelength(1:self_%wavelengthCount-1)    &
         &                                         )
    spectrumRateOfChange(1:self_%wavelengthCount)=+spectrumRateOfChange(1:self_%wavelengthCount)            &
         &                                        +spectralGradient                                         &
         &                                        *self_%cosmologyFunctions_%expansionRate(expansionFactor)
    ! Absorption.
    where (spectrum(1:self_%wavelengthCount) > 0.0d0)
       spectrumRateOfChange         (1:self_%wavelengthCount)                                &
            & =+spectrumRateOfChange(1:self_%wavelengthCount)                                &
            &  -spectrum            (1:self_%wavelengthCount)                                &
            &  *gigaYear                                                                     &
            &  *massSolar                                                                    &
            &  /megaParsec                                                               **3 &
            &  *centi                                                                    **2 &
            &  *speedLight                                                                   &
            &  /atomicMassUnit                                                               &
            &  *(                                                                            &
            &    +self_%crossSectionNeutralHydrogen                                          &
            &    *hydrogenByMassPrimordial                                                   &
            &    *self_%intergalacticMediumState_      %neutralHydrogenFraction    (time)    &
            &    /atomicMassHydrogen                                                         &
            &    +self_%crossSectionNeutralHelium                                            &
            &    *heliumByMassPrimordial                                                     &
            &    *self_%intergalacticMediumState_      %neutralHeliumFraction      (time)    &
            &    /atomicMassHelium                                                           &
            &    +self_%crossSectionSinglyIonizedHelium                                      &
            &    *heliumByMassPrimordial                                                     &
            &    *self_%intergalacticMediumState_      %singlyIonizedHeliumFraction(time)    &
            &    /atomicMassHelium                                                           &
            &   )                                                                            &
            &  *  self_%cosmologyParameters_           %OmegaBaryon                (    )    &
            &  *  self_%cosmologyParameters_           %densityCritical            (    )    &
            &  /  self_%cosmologyFunctions_            %expansionFactor            (time)**3
    end where
    ! Return success.
    intergalacticBackgroundInternalODEs=GSL_Success
    return
  end function intergalacticBackgroundInternalODEs

  logical function intergalacticBackgroundInternalTimeDependentOnly(self)
    !!{
    Return false as, while this radiation field depends only on time, it can not be evaluated for arbitrary times.
    !!}
    implicit none
    class(radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self
    !$GLC attributes unused :: self

    intergalacticBackgroundInternalTimeDependentOnly=.false.
    return
  end function intergalacticBackgroundInternalTimeDependentOnly

  subroutine intergalacticBackgroundInternalDeepCopyReset(self)
    !!{
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : accretionDiskSpectraDeepCopyReset_
    implicit none
    class(radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self
    
    self%copiedSelf => null()
    if (associated(self%cosmologyParameters_              )) call self%cosmologyParameters_              %deepCopyReset()
    if (associated(self%cosmologyFunctions_               )) call self%cosmologyFunctions_               %deepCopyReset()
    if (associated(self%intergalacticMediumState_         )) call self%intergalacticMediumState_         %deepCopyReset()
    if (associated(self%atomicCrossSectionIonizationPhoto_)) call self%atomicCrossSectionIonizationPhoto_%deepCopyReset()
    if (associated(self%starFormationRateDisks_           )) call self%starFormationRateDisks_           %deepCopyReset()
    if (associated(self%starFormationRateSpheroids_       )) call self%starFormationRateSpheroids_       %deepCopyReset()
    if (associated(self%stellarPopulationSelector_        )) call self%stellarPopulationSelector_        %deepCopyReset()
    if (associated(self%outputTimes_                      )) call self%outputTimes_                      %deepCopyReset()
    if (associated(self%accretionDiskSpectra_             )) call accretionDiskSpectraDeepCopyReset_(self%accretionDiskSpectra_)
    return
  end subroutine intergalacticBackgroundInternalDeepCopyReset
  
  subroutine intergalacticBackgroundInternalDeepCopyFinalize(self)
    !!{
    Finalize a deep reset of the object.
    !!}
    use :: Functions_Global, only : accretionDiskSpectraDeepCopyFinalize_
    implicit none
    class(radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self
    
    if (associated(self%cosmologyParameters_              )) call self%cosmologyParameters_              %deepCopyFinalize()
    if (associated(self%cosmologyFunctions_               )) call self%cosmologyFunctions_               %deepCopyFinalize()
    if (associated(self%intergalacticMediumState_         )) call self%intergalacticMediumState_         %deepCopyFinalize()
    if (associated(self%atomicCrossSectionIonizationPhoto_)) call self%atomicCrossSectionIonizationPhoto_%deepCopyFinalize()
    if (associated(self%starFormationRateDisks_           )) call self%starFormationRateDisks_           %deepCopyFinalize()
    if (associated(self%starFormationRateSpheroids_       )) call self%starFormationRateSpheroids_       %deepCopyFinalize()
    if (associated(self%stellarPopulationSelector_        )) call self%stellarPopulationSelector_        %deepCopyFinalize()
    if (associated(self%outputTimes_                      )) call self%outputTimes_                      %deepCopyFinalize()
    if (associated(self%accretionDiskSpectra_             )) call accretionDiskSpectraDeepCopyFinalize_(self%accretionDiskSpectra_)
    return
  end subroutine intergalacticBackgroundInternalDeepCopyFinalize
  
  subroutine intergalacticBackgroundInternalDeepCopy(self,destination)
    !!{
    Perform a deep copy of the object.
    !!}
    use :: Functions_Global, only : accretionDiskSpectraDeepCopy_
    implicit none
    class(radiationFieldIntergalacticBackgroundInternal), intent(inout), target :: self
    class(radiationFieldClass                          ), intent(inout)         :: destination

    call self%deepCopy_(destination)
    select type (destination)
    type is (radiationFieldIntergalacticBackgroundInternal)
       nullify(destination%accretionDiskSpectra_)
       if (associated(self%accretionDiskSpectra_)) then
          allocate(destination%accretionDiskSpectra_,mold=self%accretionDiskSpectra_)
          call accretionDiskSpectraDeepCopy_(self%accretionDiskSpectra_,destination%accretionDiskSpectra_)
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): galacticstructure : [destination] : ')//loc(destination)//' : '//loc(destination%accretionDiskSpectra_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine intergalacticBackgroundInternalDeepCopy

