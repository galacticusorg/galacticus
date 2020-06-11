!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a class for intergalactic background light which computes the background internally.

  use :: Accretion_Disk_Spectra                , only : accretionDiskSpectraClass
  use :: Atomic_Cross_Sections_Ionization_Photo, only : atomicCrossSectionIonizationPhotoClass
  use :: Cosmology_Functions                   , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                  , only : cosmologyParametersClass
  use :: Intergalactic_Medium_State            , only : intergalacticMediumStateClass
  use :: Numerical_Interpolation               , only : interpolator
  use :: Output_Times                          , only : outputTimesClass
  use :: Star_Formation_Rates_Disks            , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids        , only : starFormationRateSpheroidsClass
  use :: Stellar_Population_Selectors          , only : stellarPopulationSelectorClass

  !# <radiationField name="radiationFieldIntergalacticBackgroundInternal">
  !#  <description>A radiation field class for intergalactic background light with properties computed internally.</description>
  !# </radiationField>
  type, extends(radiationFieldIntergalacticBackground) :: radiationFieldIntergalacticBackgroundInternal
     !% A radiation field class for intergalactic background light with properties computed internally
     private
     class           (cosmologyParametersClass              ), pointer                     :: cosmologyParameters_               => null()
     class           (cosmologyFunctionsClass               ), pointer                     :: cosmologyFunctions_                => null()
     class           (intergalacticMediumStateClass         ), pointer                     :: intergalacticMediumState_          => null()
     class           (atomicCrossSectionIonizationPhotoClass), pointer                     :: atomicCrossSectionIonizationPhoto_ => null()
     class           (accretionDiskSpectraClass             ), pointer                     :: accretionDiskSpectra_              => null()
     class           (starFormationRateDisksClass           ), pointer                     :: starFormationRateDisks_            => null()
     class           (starFormationRateSpheroidsClass       ), pointer                     :: starFormationRateSpheroids_        => null()
     class           (stellarPopulationSelectorClass        ), pointer                     :: stellarPopulationSelector_         => null()
     class           (outputTimesClass                      ), pointer                     :: outputTimes_                       => null()
     integer                                                                               :: wavelengthCountPerDecade                    , wavelengthCount
     double precision                                                                      :: wavelengthMinimum                           , wavelengthMaximum
     integer                                                                               :: timeCountPerDecade                          , timeCount
     double precision                                                                      :: redshiftMinimum                             , redshiftMaximum
     double precision                                                                      :: timeMinimum                                 , timeMaximum
     double precision                                        , allocatable, dimension(:  ) :: wavelength                                  , redshift                       , &
          &                                                                                   time                                        , crossSectionNeutralHydrogen    , &
          &                                                                                   crossSectionNeutralHelium                   , crossSectionSinglyIonizedHelium, &
          &                                                                                   spectrum
     double precision                                        , allocatable, dimension(:,:) :: emissivityODE                               , emissivity
     double precision                                                     , dimension(0:1) :: timeODE
     double precision                                                                      :: timeCurrent
     type            (interpolator                          )                              :: interpolatorWavelength                      , interpolatorTime
   contains
     final     ::             intergalacticBackgroundInternalDestructor
     procedure :: flux     => intergalacticBackgroundInternalFlux
     procedure :: timeSet  => intergalacticBackgroundInternalTimeSet
     procedure :: autoHook => intergalacticBackgroundInternalAutoHook
  end type radiationFieldIntergalacticBackgroundInternal

  interface radiationFieldIntergalacticBackgroundInternal
     !% Constructors for the {\normalfont \ttfamily intergalacticBackgroundInternal} radiation field class.
     module procedure intergalacticBackgroundInternalConstructorParameters
     module procedure intergalacticBackgroundInternalConstructorInternal
  end interface radiationFieldIntergalacticBackgroundInternal

  type :: intergalacticBackgroundInternalState
     !% Class used to store the state of the intergalactic background radiation field for the internal solver. This will be stored
     !% as an attribute of the universe object.
     double precision, allocatable, dimension(:,:) :: flux
     double precision                              :: timePrevious, timeNext
  end type intergalacticBackgroundInternalState

  ! Module-scope pointer to self for ODE solving.
  class(radiationFieldIntergalacticBackgroundInternal), pointer :: intergalacticBackgroundInternalSelf
  !$omp threadprivate(intergalacticBackgroundInternalSelf)

contains

  function intergalacticBackgroundInternalConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily intergalacticBackgroundInternal} radiation field class which takes a parameter list as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (radiationFieldIntergalacticBackgroundInternal)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (cosmologyParametersClass                     ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                      ), pointer       :: cosmologyFunctions_
    class           (intergalacticMediumStateClass                ), pointer       :: intergalacticMediumState_
    class           (atomicCrossSectionIonizationPhotoClass       ), pointer       :: atomicCrossSectionIonizationPhoto_
    class           (accretionDiskSpectraClass                    ), pointer       :: accretionDiskSpectra_
    class           (starFormationRateDisksClass                  ), pointer       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass              ), pointer       :: starFormationRateSpheroids_
    class           (stellarPopulationSelectorClass               ), pointer       :: stellarPopulationSelector_
    class           (outputTimesClass                             ), pointer       :: outputTimes_
    integer                                                                        :: wavelengthCountPerDecade          , timeCountPerDecade
    double precision                                                               :: wavelengthMinimum                 , wavelengthMaximum , &
         &                                                                            redshiftMinimum                   , redshiftMaximum

    !# <inputParameter>
    !#   <name>wavelengthCountPerDecade</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of wavelength to use for calculations of the cosmic background radiation.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>wavelengthMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>100.0d0</defaultValue>
    !#   <description>The minimum wavelength (in units of \AA) to use in calculations of the cosmic background radiation.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>wavelengthMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>100000.0d0</defaultValue>
    !#   <description>The maximum wavelength (in units of \AA) to use in calculations of the cosmic background radiation.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeCountPerDecade</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of bins per decade of time to use for calculations of tge cosmic background radiation.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshiftMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The minimum redshift to use in calculations of the cosmic background radiation.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshiftMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>30.0d0</defaultValue>
    !#   <description>The maximum redshift to use in calculations of the cosmic background radiation.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters"               name="cosmologyParameters_"               source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"                name="cosmologyFunctions_"                source="parameters"/>
    !# <objectBuilder class="intergalacticMediumState"          name="intergalacticMediumState_"          source="parameters"/>
    !# <objectBuilder class="atomicCrossSectionIonizationPhoto" name="atomicCrossSectionIonizationPhoto_" source="parameters"/>
    !# <objectBuilder class="accretionDiskSpectra"              name="accretionDiskSpectra_"              source="parameters"/>
    !# <objectBuilder class="starFormationRateDisks"            name="starFormationRateDisks_"            source="parameters"/>
    !# <objectBuilder class="starFormationRateSpheroids"        name="starFormationRateSpheroids_"        source="parameters"/>
    !# <objectBuilder class="stellarPopulationSelector"         name="stellarPopulationSelector_"         source="parameters"/>
    !# <objectBuilder class="outputTimes"                       name="outputTimes_"                       source="parameters"/>
    self=radiationFieldIntergalacticBackgroundInternal(wavelengthMinimum,wavelengthMaximum,wavelengthCountPerDecade,redshiftMinimum,redshiftMaximum,timeCountPerDecade,cosmologyParameters_,cosmologyFunctions_,intergalacticMediumState_,atomicCrossSectionIonizationPhoto_,accretionDiskSpectra_,starFormationRateDisks_,starFormationRateSpheroids_,stellarPopulationSelector_,outputTimes_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"              />
    !# <objectDestructor name="cosmologyFunctions_"               />
    !# <objectDestructor name="intergalacticMediumState_"         />
    !# <objectDestructor name="atomicCrossSectionIonizationPhoto_"/>
    !# <objectDestructor name="accretionDiskSpectra_"             />
    !# <objectDestructor name="starFormationRateDisks_"           />
    !# <objectDestructor name="starFormationRateSpheroids_"       />
    !# <objectDestructor name="stellarPopulationSelector_"        />
    !# <objectDestructor name="outputTimes_"                      />
    return
  end function intergalacticBackgroundInternalConstructorParameters

  function intergalacticBackgroundInternalConstructorInternal(wavelengthMinimum,wavelengthMaximum,wavelengthCountPerDecade,redshiftMinimum,redshiftMaximum,timeCountPerDecade,cosmologyParameters_,cosmologyFunctions_,intergalacticMediumState_,atomicCrossSectionIonizationPhoto_,accretionDiskSpectra_,starFormationRateDisks_,starFormationRateSpheroids_,stellarPopulationSelector_,outputTimes_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily intergalacticBackgroundInternal} radiation field class.
    use :: Memory_Management, only : allocateArray
    use :: Numerical_Ranges , only : Make_Range   , rangeTypeLogarithmic
    implicit none
    type            (radiationFieldIntergalacticBackgroundInternal)                        :: self
    integer                                                        , intent(in   )         :: wavelengthCountPerDecade          , timeCountPerDecade
    double precision                                               , intent(in   )         :: wavelengthMinimum                 , wavelengthMaximum , &
         &                                                                                    redshiftMinimum                   , redshiftMaximum
    class           (cosmologyParametersClass                     ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass                      ), intent(in   ), target :: cosmologyFunctions_
    class           (intergalacticMediumStateClass                ), intent(in   ), target :: intergalacticMediumState_
    class           (atomicCrossSectionIonizationPhotoClass       ), intent(in   ), target :: atomicCrossSectionIonizationPhoto_
    class           (accretionDiskSpectraClass                    ), intent(in   ), target :: accretionDiskSpectra_
    class           (starFormationRateDisksClass                  ), intent(in   ), target :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass              ), intent(in   ), target :: starFormationRateSpheroids_
    class           (stellarPopulationSelectorClass               ), intent(in   ), target :: stellarPopulationSelector_
    class           (outputTimesClass                             ), intent(in   ), target :: outputTimes_
    integer                                                                                :: iTime                             , iWavelength
    !# <constructorAssign variables="wavelengthMinimum, wavelengthMaximum, wavelengthCountPerDecade, redshiftMinimum, redshiftMaximum, timeCountPerDecade, *cosmologyParameters_, *cosmologyFunctions_, *intergalacticMediumState_, *atomicCrossSectionIonizationPhoto_, *accretionDiskSpectra_, *starFormationRateDisks_, *starFormationRateSpheroids_, *stellarPopulationSelector_, *outputTimes_"/>

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
    call allocateArray(self%wavelength   ,[self%wavelengthCount               ]                  )
    call allocateArray(self%spectrum     ,[self%wavelengthCount               ]                  )
    call allocateArray(self%time         ,[                     self%timeCount]                  )
    call allocateArray(self%redshift     ,[                     self%timeCount]                  )
    call allocateArray(self%emissivity   ,[self%wavelengthCount,self%timeCount]                  )
    call allocateArray(self%emissivityODE,[self%wavelengthCount,2             ],lowerBounds=[1,0])
    self%wavelength=Make_Range(                        &
         &                     self%wavelengthMinimum, &
         &                     self%wavelengthMaximum, &
         &                     self%wavelengthCount  , &
         &                     rangeTypeLogarithmic    &
         &                    )
    self%time      =Make_Range(                        &
         &                     self%timeMinimum      , &
         &                     self%timeMaximum      , &
         &                     self%timeCount        , &
         &                     rangeTypeLogarithmic    &
         &                    )
    ! Convert times to redshifts.
    do iTime=1,self%timeCount
       self%redshift(iTime)                                                           &
            & =self%cosmologyFunctions_%redshiftFromExpansionFactor(                  &
            &  self%cosmologyFunctions_%expansionFactor             (                 &
            &                                                        self%time(iTime) &
            &                                                       )                 &
            &                                                      )
    end do
    ! Initialize the background radiation to zero.
    self%spectrum  =0.0d0
    ! Initialize the emissivity to zero.
    self%emissivity=0.0d0
    ! Construct tables of photoionization cross-sections.
    call allocateArray(self%crossSectionNeutralHydrogen    ,[self%wavelengthCount])
    call allocateArray(self%crossSectionNeutralHelium      ,[self%wavelengthCount])
    call allocateArray(self%crossSectionSinglyIonizedHelium,[self%wavelengthCount])
    do iWavelength=1,self%wavelengthCount
       self%crossSectionNeutralHydrogen    (iWavelength)=self%atomicCrossSectionIonizationPhoto_%crossSection(1,1,1,self%wavelength(iWavelength))
       self%crossSectionNeutralHelium      (iWavelength)=self%atomicCrossSectionIonizationPhoto_%crossSection(2,1,1,self%wavelength(iWavelength))
       self%crossSectionSinglyIonizedHelium(iWavelength)=self%atomicCrossSectionIonizationPhoto_%crossSection(2,2,1,self%wavelength(iWavelength))
    end do
    ! Build interpolators.
    self%interpolatorWavelength=interpolator(self%wavelength)
    self%interpolatorTime      =interpolator(self%time      )
    return
  end function intergalacticBackgroundInternalConstructorInternal

  subroutine intergalacticBackgroundInternalAutoHook(self)
    use :: Events_Hooks, only : universePreEvolveEvent
    implicit none
    class(radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self

    ! Hook to universe pre-evolve events.
    !$omp master
    call universePreEvolveEvent%attach(self,intergalacticBackgroundInternalUniversePreEvolve)
    !$omp end master
    return
  end subroutine intergalacticBackgroundInternalAutoHook

  subroutine intergalacticBackgroundInternalDestructor(self)
    !% Destructor for the {\normalfont \ttfamily intergalacticBackgroundInternal} radiation field class.
    implicit none
    type(radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"              />
    !# <objectDestructor name="self%cosmologyFunctions_"               />
    !# <objectDestructor name="self%intergalacticMediumState_"         />
    !# <objectDestructor name="self%atomicCrossSectionIonizationPhoto_"/>
    !# <objectDestructor name="self%accretionDiskSpectra_"             />
    !# <objectDestructor name="self%starFormationRateDisks_"           />
    !# <objectDestructor name="self%starFormationRateSpheroids_"       />
    !# <objectDestructor name="self%stellarPopulationSelector_"        />
    !# <objectDestructor name="self%outputTimes_"                      />
    return
  end subroutine intergalacticBackgroundInternalDestructor

  subroutine intergalacticBackgroundInternalTimeSet(self,time)
    !% Set the epoch.
    implicit none
    class           (radiationFieldIntergalacticBackgroundInternal), intent(inout) :: self
    double precision                                               , intent(in   ) :: time

    self%timeCurrent=time
    return
  end subroutine intergalacticBackgroundInternalTimeSet

 double precision function intergalacticBackgroundInternalFlux(self,wavelength,node)
    !% Return the flux in the internally-computed intergalatic background.
    use            :: Galacticus_Error, only : Galacticus_Error_Report
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    class           (radiationFieldIntergalacticBackgroundInternal), intent(inout)  :: self
    double precision                                               , intent(in   )  :: wavelength
    type            (treeNode                                     ), intent(inout)  :: node
    double precision                                               , dimension(0:1) :: hWavelength         , hTime
    double precision                                               , parameter      :: timeTolerance=1.0d-3
    class           (*                                            ), pointer        :: state
    integer         (c_size_t                                     )                 :: iWavelength         , jWavelength, &
         &                                                                             iTime               , jTime

    !$omp critical (radiationFieldIntergalacticBackgroundInternalCritical)
    ! Get the state of the radiation field.
    state => node%hostTree%hostUniverse%attributes%value('radiationFieldIntergalacticBackgroundInternal')
    select type (state)
    type is (intergalacticBackgroundInternalState)
       ! Check that the time is within the applicable range.
       if (self%timeCurrent > state%timeNext*(1.0d0+timeTolerance)) call Galacticus_Error_Report('time is out of range'//{introspection:location})
       ! Find interpolation in the array of wavelengths.
       call self%interpolatorWavelength%linearFactors(     wavelength ,iWavelength,hWavelength)
       ! Find interpolation in array of times.
       call self%interpolatorTime      %linearFactors(self%timeCurrent,iTime      ,hTime      )
       if (self%timeCurrent > state%timePrevious) hTime=[1.0d0,0.0d0]
       ! Interpolate in wavelength and time.
       intergalacticBackgroundInternalFlux=0.0d0
       do jTime=0,1
          do jWavelength=0,1
             intergalacticBackgroundInternalFlux=+intergalacticBackgroundInternalFlux                                      &
                  &                              +hTime                              (                        jTime      ) &
                  &                              *hWavelength                        (jWavelength                        ) &
                  &                              *state%flux                         (jWavelength+iWavelength,jTime+iTime)
          end do
       end do
    class default
       intergalacticBackgroundInternalFlux=0.0d0
    end select
    !$omp end critical (radiationFieldIntergalacticBackgroundInternalCritical)
    intergalacticBackgroundInternalFlux=max(intergalacticBackgroundInternalFlux,0.0d0)
    return
  end function intergalacticBackgroundInternalFlux

  subroutine intergalacticBackgroundInternalUniversePreEvolve(self,universe_)
    !% Attach an initial event to the universe to cause the background radiation update function to be called.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : universe               , universeEvent
    implicit none
    class(*                                   ), intent(inout), target :: self
    type (universe                            ), intent(inout)         :: universe_
    type (universeEvent                       ), pointer               :: event
    type (intergalacticBackgroundInternalState), pointer               :: state

    select type (self)
    class is (radiationFieldIntergalacticBackgroundInternal)
       ! If the universe object already has an "radiationFieldIntergalacticBackgroundInternal" attribute, then do not add a new event here - we want only one event per universe.
       if (.not.universe_%attributes%exists('radiationFieldIntergalacticBackgroundInternal')) then
          ! Create the first interrupt event in the universe object.
          event                       => universe_%createEvent( )
          event%time                  =  self     %time       (1)
          event%creator               => self
          event%task                  => intergalacticBackgroundInternalUpdate
          !$omp critical (radiationFieldIntergalacticBackgroundInternalCritical)
          allocate(state                                          )
          allocate(state%flux(self%wavelengthCount,self%timeCount))
          state%timeNext    =self%time(1)
          state%timePrevious=0.0d0
          state%flux        =0.0d0
          call universe_%attributes%set('radiationFieldIntergalacticBackgroundInternal',state)
          !$omp end critical (radiationFieldIntergalacticBackgroundInternalCritical)
       end if
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine intergalacticBackgroundInternalUniversePreEvolve

  logical function intergalacticBackgroundInternalUpdate(event,universe_) result (success)
    !% Update the radiation background for a given universe.
    use            :: Abundances_Structure        , only : abundances                   , max
    use            :: Arrays_Search               , only : searchArrayClosest
    use            :: FODEIV2                     , only : fodeiv2_driver               , fodeiv2_system
    use            :: Galacticus_Display          , only : Galacticus_Display_Indent    , Galacticus_Display_Message, Galacticus_Display_Unindent
    use            :: Galacticus_Error            , only : Galacticus_Error_Report
    use            :: Galacticus_HDF5             , only : galacticusOutputFile
    use            :: Galacticus_Nodes            , only : defaultDiskComponent         , defaultSpheroidComponent  , mergerTreeList             , nodeComponentBasic, &
          &                                                nodeComponentDisk            , nodeComponentSpheroid     , treeNode                   , universe          , &
          &                                                universeEvent
    use            :: IO_HDF5                     , only : hdf5Access                   , hdf5Object
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: ISO_Varying_String          , only : varying_string
    use            :: Merger_Tree_Walkers         , only : mergerTreeWalkerAllNodes
    use            :: Numerical_Constants_Math    , only : Pi
    use            :: Numerical_Constants_Physical, only : plancksConstant              , speedLight
    use            :: Numerical_Constants_Prefixes, only : centi
    use            :: Numerical_Constants_Units   , only : angstromsPerMeter            , ergs
    use            :: Numerical_Integration       , only : integrator
    use            :: ODEIV2_Solver               , only : ODEIV2_Solve                 , ODEIV2_Solver_Free
    use            :: Stellar_Population_Spectra  , only : stellarPopulationSpectraClass
    use            :: Stellar_Populations         , only : stellarPopulationClass
    implicit none
    class           (universeEvent                       ), intent(in   ) :: event
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
    type            (abundances                          ), target        :: gasAbundancesDisk                    , gasAbundancesSpheroid
    type            (abundances                          ), pointer       :: gasAbundances
    class           (*                                   ), pointer       :: state
    type            (fodeiv2_system                      ), save          :: ode2System
    type            (fodeiv2_driver                      ), save          :: ode2Driver
    type            (integrator                          )                :: integrator_
    logical                                               , save          :: odeReset
    type            (mergerTreeWalkerAllNodes            )                :: treeWalker
    double precision                                                      :: starFormationRateDisk                , starFormationRateSpheroid               , &
         &                                                                   gasMassDisk                          , gasMassSpheroid                         , &
         &                                                                   ageEnd                               , ageStart                                , &
         &                                                                   stellarSpectrumDisk                  , stellarSpectrumSpheroid                 , &
         &                                                                   timeStart                            , timeEnd                                 , &
         &                                                                   treeTimeLatest                       , wavelength
    integer                                                               :: iTime                                , iWavelength
    type            (varying_string                      )                :: message
    character       (len=6                               )                :: label
    type            (hdf5Object                          )                :: outputGroup                          , outputDataset
    integer         (c_size_t                            )                :: iNow
    logical                                                               :: firstTime

    ! Guard on event creator class.
    select type (self => event%creator)
    class is (radiationFieldIntergalacticBackgroundInternal)
       ! Display message.
       write (label,'(f6.3)') event%time
       message="Evolving cosmic background radiation to time "//trim(label)//" Gyr"
       call Galacticus_Display_Indent(message)
       ! Find the current timestep.
       iNow=searchArrayClosest(self%time,event%time)
       ! Iterate over all nodes.
       call Galacticus_Display_Message('Accumulating emissivity')
       treeTimeLatest=0.0d0
       forest => universe_%trees
       do while (associated(forest))
          treeWalker=mergerTreeWalkerAllNodes(forest%tree,spanForest=.true.)
          do while (treeWalker%next(node))
             basic          =>                    node %basic()
             treeTimeLatest =  max(treeTimeLatest,basic%time ())
             if (basic%time() == event%time) then
                ! Get the star formation rates and metallicites for this node.
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
                      if (self%time(iTime) < event%time) cycle
                      ! Compute age of the currently forming population at this time.
                      ageEnd=self%time(iTime)-event%time
                      if (iTime == 1) then
                         ageStart=                                  0.0d0
                      else
                         ageStart=max(self%time(iTime-1)-event%time,0.0d0)
                      end if
                      ! Iterate over wavelength
                      integrator_=integrator(stellarSpectraConvolution,toleranceAbsolute=integrationToleranceAbsolute,toleranceRelative=integrationToleranceRelative)
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
                         if (firstTime)                                                                                               &
                              & self%emissivity  (iWavelength,iTime)=+self%emissivity                        (iWavelength,iTime     ) &
                              &                                      +self%accretionDiskSpectra_%spectrum    (node       ,wavelength) &
                              &                                      *node%hostTree             %volumeWeight
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
             call Galacticus_Display_Message('Solving cosmic background radiation evolution')
             self%timeODE      (  0:1)=self%time      (  iNow-1:iNow)
             self%emissivityODE(:,0  )=self%emissivity(:,iNow-1     )
             self%emissivityODE(:,  1)=self%emissivity(:,       iNow)
             intergalacticBackgroundInternalSelf => self
             odeReset=.true.
             timeStart=self%time(iNow-1)
             timeEnd  =self%time(iNow  )
             call ODEIV2_Solve(                                            &
                  &            ode2Driver                                , &
                  &            ode2System                                , &
                  &            timeStart                                 , &
                  &            timeEnd                                   , &
                  &            self%wavelengthCount                      , &
                  &            self%spectrum                             , &
                  &            intergalacticBackgroundInternalODEs       , &
                  &            odeToleranceAbsolute                      , &
                  &            odeToleranceRelative                      , &
                  &            reset=odeReset                              &
                  &           )
             call ODEIV2_Solver_Free(ode2Driver,ode2System)
             ! Convert
             !$omp critical (radiationFieldIntergalacticBackgroundInternalCritical)
             state%flux(:,iNow)=max(                                   &
                  &                 +plancksConstant                   &
                  &                 *speedLight                   **2  &
                  &                 *angstromsPerMeter                 &
                  &                 /4.0d0                             &
                  &                 /Pi                                &
                  &                 /self%wavelength                   &
                  &                 *self%spectrum                     &
                  &                 *centi                        **2  &
                  &                 /ergs                            , &
                  &                 +0.0d0                             &
                  &                )
             !$omp end critical (radiationFieldIntergalacticBackgroundInternalCritical)
          end if
          ! Add the next event to the universe.
          !$omp critical (radiationFieldIntergalacticBackgroundInternalCritical)
          state%timeNext=self%time(iNow+1)
          !$omp end critical (radiationFieldIntergalacticBackgroundInternalCritical)
          if     (                                                                           &
               &             iNow    <                        self             %timeCount    &
               &  .and.                                                                      &
               &   self%time(iNow+1) < treeTimeLatest                                        &
               &  .and.                                                                      &
               &   self%time(iNow+1) < self%outputTimes_%time(self%outputTimes_%count    ()) &
               & ) then
             !$omp critical (radiationFieldIntergalacticBackgroundInternalCritical)
             state%timePrevious=self%time(iNow)
             !$omp end critical (radiationFieldIntergalacticBackgroundInternalCritical)
             eventNew         => universe_%createEvent(      )
             eventNew%time    =  self     %time       (iNow+1)
             eventNew%creator => event    %creator
             eventNew%task    => intergalacticBackgroundInternalUpdate
          else
             ! Output the results to file.
             call hdf5Access%set()
             outputGroup=galacticusOutputFile%openGroup('backgroundRadiation','Cosmic background radiation data.')
             call outputGroup  %writeDataset  (self%wavelength        ,'wavelength','Wavelength at which the background radiation is tabulated [Å].'    ,datasetReturned=outputDataset)
             call outputDataset%writeAttribute(1.0d0/angstromsPerMeter,'unitsInSI'                                                                                                    )
             call outputDataset%close         (                                                                                                                                       )
             call outputGroup  %writeDataset  (self%redshift          ,'redshift'  ,'Redshift at which the background radiation is tabulated [].'       ,datasetReturned=outputDataset)
             call outputDataset%writeAttribute(0.0d0                  ,'unitsInSI'                                                                                                    )
             call outputDataset%close         (                                                                                                                                       )
             call outputGroup  %writeDataset  (state%flux             ,'flux'      ,'Flux is the cosmic background radiation [erg cm⁻² s⁻¹ Hz⁻¹ sr⁻¹].' ,datasetReturned=outputDataset)
             call outputDataset%writeAttribute(ergs/centi**2          ,'unitsInSI'                                                                                                    )
             call outputDataset%close         (                                                                                                                                       )
             call outputGroup  %close         (                                                                                                                                       )
             call hdf5Access   %unset         (                                                                                                                                       )
          end if
       end select
       ! Display message.
       call Galacticus_Display_Unindent('done')
       ! Return true since we've performed our task.
       success=.true.
    class default
       ! Incorrect event creator type.
       success=.false.
       call Galacticus_Error_Report('incorrect event creator class'//{introspection:location})
    end select
    return

  contains

    double precision function stellarSpectraConvolution(age)
      !% Integrand for convolution of stellar spectra.
      use :: Galacticus_Error, only : errorStatusInputDomain, errorStatusSuccess
      implicit none
      double precision, intent(in   ) :: age
      integer                         :: status

      stellarSpectraConvolution=stellarPopulationSpectra_%luminosity(               &
           &                                                         gasAbundances, &
           &                                                         age          , &
           &                                                         wavelength   , &
           &                                                         status         &
           &                                                        )
      if     (                                                                              &
           &   status /= errorStatusSuccess                                                 &
           &  .and.                                                                         &
           &   status /= errorStatusInputDomain                                             &
           & ) call Galacticus_Error_Report(                                                &
           &                                'stellar population spectrum function failed'// &
           &                                {introspection:location}                        &
           &                               )
      return
    end function stellarSpectraConvolution

  end function intergalacticBackgroundInternalUpdate

  integer function intergalacticBackgroundInternalODEs(time,spectrum,spectrumRateOfChange)
    !% Evaluates the ODEs controlling the evolution of cosmic background radiation.
    use :: Interface_GSL                   , only : GSL_Success
    use :: Numerical_Constants_Astronomical, only : gigaYear         , heliumByMassPrimordial, hydrogenByMassPrimordial, luminositySolar, &
          &                                         massSolar        , megaParsec
    use :: Numerical_Constants_Atomic      , only : atomicMassHelium , atomicMassHydrogen    , atomicMassUnit
    use :: Numerical_Constants_Physical    , only : plancksConstant  , speedLight
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Numerical_Constants_Units       , only : angstromsPerMeter
    implicit none
    double precision, intent(in   )               :: time
    double precision, intent(in   ), dimension(:) :: spectrum
    double precision, intent(  out), dimension(:) :: spectrumRateOfChange
    double precision                              :: spectralGradient    (intergalacticBackgroundInternalSelf%wavelengthCount)
    double precision                              :: expansionFactor

    ! Get the expansion factor.
    expansionFactor=intergalacticBackgroundInternalSelf%cosmologyFunctions_%expansionFactor(time)
    ! Add source terms, linearly interpolating between timesteps. Convert from emissivity units [L☉ Hz⁻¹ Mpc⁻³] to
    ! background units [photons m⁻³ Hz⁻¹].
    spectrumRateOfChange(1:intergalacticBackgroundInternalSelf%wavelengthCount)                                                &
         & =+(                                                                                                                 &
         &    +                                                        intergalacticBackgroundInternalSelf%emissivityODE(:,0)  &
         &    +(intergalacticBackgroundInternalSelf%emissivityODE(:,1)-intergalacticBackgroundInternalSelf%emissivityODE(:,0)) &
         &    *(                                          time        -intergalacticBackgroundInternalSelf%      timeODE(  0)) &
         &    /(intergalacticBackgroundInternalSelf%      timeODE(  1)-intergalacticBackgroundInternalSelf%      timeODE(  0)) &
         &   )                                                                                                                 &
         &  *intergalacticBackgroundInternalSelf%wavelength                                                                    &
         &  *luminositySolar                                                                                                   &
         &  *gigaYear                                                                                                          &
         &  /angstromsPerMeter                                                                                                 &
         &  /plancksConstant                                                                                                   &
         &  /speedLight                                                                                                        &
         &  /megaParsec**3
    ! Add expansion dilution: -3 H(t)      n_ν
    spectrumRateOfChange(1:intergalacticBackgroundInternalSelf%wavelengthCount)=+spectrumRateOfChange(1:intergalacticBackgroundInternalSelf%wavelengthCount)            &
         &                                                                      -spectrum            (1:intergalacticBackgroundInternalSelf%wavelengthCount)            &
         &                                                                      *3.0d0                                                                                  &
         &                                                                      *intergalacticBackgroundInternalSelf%cosmologyFunctions_%expansionRate(expansionFactor)
    ! Add redshifting       : +  H(t) d(ν n_ν)/dν
    spectralGradient    (1                                                    )=-  spectrum          (1                                                      )
    spectralGradient    (2:intergalacticBackgroundInternalSelf%wavelengthCount)=-                                                                              intergalacticBackgroundInternalSelf%wavelength(2:intergalacticBackgroundInternalSelf%wavelengthCount  )**2 &
         &                                                                      *(                                                                                                                                                                                        &
         &                                                                        +spectrum          (2:intergalacticBackgroundInternalSelf%wavelengthCount  )/intergalacticBackgroundInternalSelf%wavelength(2:intergalacticBackgroundInternalSelf%wavelengthCount  )    &
         &                                                                        -spectrum          (1:intergalacticBackgroundInternalSelf%wavelengthCount-1)/intergalacticBackgroundInternalSelf%wavelength(1:intergalacticBackgroundInternalSelf%wavelengthCount-1)    &
         &                                                                       )                                                                                                                                                                                        &
         &                                                                      /(                                                                                                                                                                                        &
         &                                                                        +                                                                            intergalacticBackgroundInternalSelf%wavelength(2:intergalacticBackgroundInternalSelf%wavelengthCount  )    &
         &                                                                        -                                                                            intergalacticBackgroundInternalSelf%wavelength(1:intergalacticBackgroundInternalSelf%wavelengthCount-1)    &
         &                                                                      )
    spectrumRateOfChange(1:intergalacticBackgroundInternalSelf%wavelengthCount)=+spectrumRateOfChange(1:intergalacticBackgroundInternalSelf%wavelengthCount  )                                                                                                            &
         &                                                                      +spectralGradient                                                                                                                                                                         &
         &                                                                      *intergalacticBackgroundInternalSelf%cosmologyFunctions_%expansionRate(expansionFactor)
    ! Absorption.
    where (spectrum(1:intergalacticBackgroundInternalSelf%wavelengthCount) > 0.0d0)
       spectrumRateOfChange         (1:intergalacticBackgroundInternalSelf%wavelengthCount)                                &
            & =+spectrumRateOfChange(1:intergalacticBackgroundInternalSelf%wavelengthCount)                                &
            &  -spectrum            (1:intergalacticBackgroundInternalSelf%wavelengthCount)                                &
            &  *gigaYear                                                                                                   &
            &  *massSolar                                                                                                  &
            &  /megaParsec                                                                                             **3 &
            &  *centi                                                                                                  **2 &
            &  *speedLight                                                                                                 &
            &  /atomicMassUnit                                                                                             &
            &  *(                                                                                                          &
            &    +intergalacticBackgroundInternalSelf%crossSectionNeutralHydrogen                                          &
            &    *hydrogenByMassPrimordial                                                                                 &
            &    *intergalacticBackgroundInternalSelf%intergalacticMediumState_      %neutralHydrogenFraction    (time)    &
            &    /atomicMassHydrogen                                                                                       &
            &    +intergalacticBackgroundInternalSelf%crossSectionNeutralHelium                                            &
            &    *heliumByMassPrimordial                                                                                   &
            &    *intergalacticBackgroundInternalSelf%intergalacticMediumState_      %neutralHeliumFraction      (time)    &
            &    /atomicMassHelium                                                                                         &
            &    +intergalacticBackgroundInternalSelf%crossSectionSinglyIonizedHelium                                      &
            &    *heliumByMassPrimordial                                                                                   &
            &    *intergalacticBackgroundInternalSelf%intergalacticMediumState_      %singlyIonizedHeliumFraction(time)    &
            &    /atomicMassHelium                                                                                         &
            &   )                                                                                                          &
            &  *  intergalacticBackgroundInternalSelf%cosmologyParameters_           %OmegaBaryon                (    )    &
            &  *  intergalacticBackgroundInternalSelf%cosmologyParameters_           %densityCritical            (    )    &
            &  /  intergalacticBackgroundInternalSelf%cosmologyFunctions_            %expansionFactor            (time)**3
    end where
    ! Return success.
    intergalacticBackgroundInternalODEs=GSL_Success
    return
  end function intergalacticBackgroundInternalODEs
