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

  !+ Contributions to this file made by: Daniel McAndrew.

  use :: Atomic_Cross_Sections_Ionization_Photo      , only : atomicCrossSectionIonizationPhoto      , atomicCrossSectionIonizationPhotoClass
  use :: Atomic_Ionization_Potentials                , only : atomicIonizationPotential              , atomicIonizationPotentialClass
  use :: Atomic_Radiation_Gaunt_Factors              , only : gauntFactor                            , gauntFactorClass
  use :: Atomic_Rates_Excitation_Collisional         , only : atomicExcitationRateCollisional        , atomicExcitationRateCollisionalClass
  use :: Atomic_Rates_Ionization_Collisional         , only : atomicIonizationRateCollisional        , atomicIonizationRateCollisionalClass
  use :: Atomic_Rates_Recombination_Dielectronic     , only : atomicRecombinationRateDielectronic    , atomicRecombinationRateDielectronicClass
  use :: Atomic_Rates_Recombination_Radiative        , only : atomicRecombinationRateRadiative       , atomicRecombinationRateRadiativeClass       , recombinationCaseB
  use :: Atomic_Rates_Recombination_Radiative_Cooling, only : atomicRecombinationRateRadiativeCooling, atomicRecombinationRateRadiativeCoolingClass
  use :: Cosmological_Density_Field                  , only : cosmologicalMassVariance               , cosmologicalMassVarianceClass
  use :: Cosmology_Functions                         , only : cosmologyFunctions                     , cosmologyFunctionsClass
  use :: Cosmology_Parameters                        , only : cosmologyParameters                    , cosmologyParametersClass                    , hubbleUnitsTime
  use :: Galacticus_Nodes                            , only : treeNode
  use :: Intergalactic_Medium_State                  , only : intergalacticMediumState               , intergalacticMediumStateClass               , intergalacticMediumStateInternal     , intergalacticMediumStateRecFast
  use :: Linear_Growth                               , only : linearGrowth                           , linearGrowthClass
  use :: Output_Times                                , only : outputTimes                            , outputTimesClass
  use :: Radiation_Fields                            , only : radiationField                         , radiationFieldClass                         , radiationFieldIntergalacticBackground

  !![
  <universeOperator name="universeOperatorIntergalacticMediumStateEvolve">
   <description>An operator on universes which attaches hooks to compute evolution of the intergalactic medium.</description>
  </universeOperator>
  !!]
  type, extends(universeOperatorClass) :: universeOperatorIntergalacticMediumStateEvolve
     !!{
     Implementation of a universeOperator which computes and outputs the power spectrum and related quantities.
     !!}
     private
     class           (outputTimesClass                            ), pointer                     :: outputTimes_                             => null()
     class           (cosmologyParametersClass                    ), pointer                     :: cosmologyParameters_                     => null()
     class           (cosmologyFunctionsClass                     ), pointer                     :: cosmologyFunctions_                      => null()
     class           (linearGrowthClass                           ), pointer                     :: linearGrowth_                            => null()
     class           (cosmologicalMassVarianceClass               ), pointer                     :: cosmologicalMassVariance_                => null()
     class           (radiationFieldClass                         ), pointer                     :: radiationField_                          => null()
     class           (gauntFactorClass                            ), pointer                     :: gauntFactor_                             => null()
     class           (atomicCrossSectionIonizationPhotoClass      ), pointer                     :: atomicCrossSectionIonizationPhoto_       => null()
     class           (atomicIonizationPotentialClass              ), pointer                     :: atomicIonizationPotential_               => null()
     class           (atomicRecombinationRateDielectronicClass    ), pointer                     :: atomicRecombinationRateDielectronic_     => null()
     class           (atomicRecombinationRateRadiativeClass       ), pointer                     :: atomicRecombinationRateRadiative_        => null()
     class           (atomicRecombinationRateRadiativeCoolingClass), pointer                     :: atomicRecombinationRateRadiativeCooling_ => null()
     class           (atomicIonizationRateCollisionalClass        ), pointer                     :: atomicIonizationRateCollisional_         => null()
     class           (atomicExcitationRateCollisionalClass        ), pointer                     :: atomicExcitationRateCollisional_         => null()
     class           (intergalacticMediumStateClass               ), pointer                     :: intergalacticMediumState_                => null()
     type            (intergalacticMediumStateRecFast             )                              :: intergalacticMediumStateInitial
     double precision                                              , allocatable, dimension(:  ) :: temperature                                       , massFiltering  , &
          &                                                                                         clumpingFactor                                    , opticalDepth
     integer                                                                                     :: timeCountPerDecade                                , timeCount
     double precision                                                                            :: redshiftMinimum                                   , redshiftMaximum
     double precision                                                                            :: timeMinimum                                       , timeMaximum
     double precision                                              , allocatable, dimension(:  ) :: time                                              , redshift
     double precision                                              , allocatable, dimension(:,:) :: densityHydrogen                                   , densityHelium  , &
          &                                                                                         massFilteringComposite
   contains
     !![
     <methods>
       <method description="Set the state of the IGM state class up to the given time index." method="stateSet" />
     </methods>
     !!]
     final     ::             intergalacticMediumStateEvolveDestructor
     procedure :: operate  => intergalacticMediumStateEvolveOperate
     procedure :: stateSet => intergalacticMediumStateEvolveStateSet
  end type universeOperatorIntergalacticMediumStateEvolve

  interface universeOperatorIntergalacticMediumStateEvolve
     !!{
     Constructors for the \refClass{universeOperatorIntergalacticMediumStateEvolve} universeOperator.
     !!}
     module procedure intergalacticMediumStateEvolveConstructorParameters
     module procedure intergalacticMediumStateEvolveConstructorInternal
  end interface universeOperatorIntergalacticMediumStateEvolve

  class(universeOperatorIntergalacticMediumStateEvolve), pointer :: self_
  type (treeNode                                      ), pointer :: node_
  !$omp threadprivate(self_,node_)

contains

  function intergalacticMediumStateEvolveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{universeOperatorIntergalacticMediumStateEvolve} universeOperator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (universeOperatorIntergalacticMediumStateEvolve)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (outputTimesClass                              ), pointer       :: outputTimes_
    class           (cosmologyParametersClass                      ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass                             ), pointer       :: linearGrowth_
    class           (cosmologicalMassVarianceClass                 ), pointer       :: cosmologicalMassVariance_
    class           (radiationFieldClass                           ), pointer       :: radiationField_
    class           (gauntFactorClass                              ), pointer       :: gauntFactor_
    class           (atomicCrossSectionIonizationPhotoClass        ), pointer       :: atomicCrossSectionIonizationPhoto_
    class           (atomicIonizationPotentialClass                ), pointer       :: atomicIonizationPotential_
    class           (atomicRecombinationRateDielectronicClass      ), pointer       :: atomicRecombinationRateDielectronic_
    class           (atomicRecombinationRateRadiativeClass         ), pointer       :: atomicRecombinationRateRadiative_
    class           (atomicRecombinationRateRadiativeCoolingClass  ), pointer       :: atomicRecombinationRateRadiativeCooling_
    class           (atomicIonizationRateCollisionalClass          ), pointer       :: atomicIonizationRateCollisional_
    class           (atomicExcitationRateCollisionalClass          ), pointer       :: atomicExcitationRateCollisional_
    class           (intergalacticMediumStateClass                 ), pointer       :: intergalacticMediumState_
    integer                                                                         :: timeCountPerDecade
    double precision                                                                :: redshiftMinimum                         , redshiftMaximum, &
         &                                                                             timeMinimum                             , timeMaximum

    !![
    <inputParameter>
      <name>timeCountPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of time to use for calculations of the properties of the universe.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum redshift to use in calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <defaultValue>400.0d0</defaultValue>
      <description>The maximum redshift to use in calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"                     name="cosmologyParameters_"                     source="parameters"                                                      />
    <objectBuilder class="cosmologyFunctions"                      name="cosmologyFunctions_"                      source="parameters"                                                      />
    <objectBuilder class="linearGrowth"                            name="linearGrowth_"                            source="parameters"                                                      />
    <objectBuilder class="cosmologicalMassVariance"                name="cosmologicalMassVariance_"                source="parameters"                                                      />
    <objectBuilder class="outputTimes"                             name="outputTimes_"                             source="parameters"                                                      />
    <objectBuilder class="gauntFactor"                             name="gauntFactor_"                             source="parameters"                                                      />
    <objectBuilder class="atomicCrossSectionIonizationPhoto"       name="atomicCrossSectionIonizationPhoto_"       source="parameters"                                                      />
    <objectBuilder class="atomicIonizationPotential"               name="atomicIonizationPotential_"               source="parameters"                                                      />
    <objectBuilder class="atomicRecombinationRateDielectronic"     name="atomicRecombinationRateDielectronic_"     source="parameters"                                                      />
    <objectBuilder class="atomicRecombinationRateRadiative"        name="atomicRecombinationRateRadiative_"        source="parameters"                                                      />
    <objectBuilder class="atomicRecombinationRateRadiativeCooling" name="atomicRecombinationRateRadiativeCooling_" source="parameters"                                                      />
    <objectBuilder class="atomicIonizationRateCollisional"         name="atomicIonizationRateCollisional_"         source="parameters"                                                      />
    <objectBuilder class="atomicExcitationRateCollisional"         name="atomicExcitationRateCollisional_"         source="parameters"                                                      />
    <objectBuilder class="intergalacticMediumState"                name="intergalacticMediumState_"                source="parameters"                                                      />
    <objectBuilder class="radiationField"                          name="radiationField_"                          source="parameters" parameterName="radiationFieldIntergalacticBackground"/>
    !!]
    timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    timeMaximum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
    self=universeOperatorIntergalacticMediumStateEvolve(timeMinimum,timeMaximum,timeCountPerDecade,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_,outputTimes_,gauntFactor_,atomicCrossSectionIonizationPhoto_,atomicIonizationPotential_,atomicRecombinationRateDielectronic_,atomicRecombinationRateRadiative_,atomicRecombinationRateRadiativeCooling_,atomicIonizationRateCollisional_,atomicExcitationRateCollisional_,intergalacticMediumState_,radiationField_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"                    />
    <objectDestructor name="cosmologyFunctions_"                     />
    <objectDestructor name="linearGrowth_"                           />
    <objectDestructor name="cosmologicalMassVariance_"               />
    <objectDestructor name="outputTimes_"                            />
    <objectDestructor name="gauntFactor_"                            />
    <objectDestructor name="atomicCrossSectionIonizationPhoto_"      />
    <objectDestructor name="atomicIonizationPotential_"              />
    <objectDestructor name="atomicRecombinationRateDielectronic_"    />
    <objectDestructor name="atomicRecombinationRateRadiativeCooling_"/>
    <objectDestructor name="atomicRecombinationRateRadiative_"       />
    <objectDestructor name="atomicIonizationRateCollisional_"        />
    <objectDestructor name="atomicExcitationRateCollisional_"        />
    <objectDestructor name="intergalacticMediumState_"               />
    <objectDestructor name="radiationField_"                         />
    !!]
    return
  end function intergalacticMediumStateEvolveConstructorParameters

  function intergalacticMediumStateEvolveConstructorInternal(timeMinimum,timeMaximum,timeCountPerDecade,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_,outputTimes_,gauntFactor_,atomicCrossSectionIonizationPhoto_,atomicIonizationPotential_,atomicRecombinationRateDielectronic_,atomicRecombinationRateRadiative_,atomicRecombinationRateRadiativeCooling_,atomicIonizationRateCollisional_,atomicExcitationRateCollisional_,intergalacticMediumState_,radiationField_) result(self)
    !!{
    Internal constructor for the \refClass{universeOperatorIntergalacticMediumStateEvolve} universeOperator class.
    !!}
    use            :: Error                                , only : Error_Report
    use, intrinsic :: ISO_C_Binding                        , only : c_size_t
    use            :: Intergalactic_Medium_Filtering_Masses, only : intergalacticMediumFilteringMassGnedin2000
    use            :: Numerical_Comparison                 , only : Values_Agree
    use            :: Numerical_Constants_Astronomical     , only : heliumByMassPrimordial                    , hydrogenByMassPrimordial, massSolar     , megaParsec
    use            :: Numerical_Constants_Atomic           , only : atomicMassHelium                          , atomicMassHydrogen      , atomicMassUnit
    use            :: Numerical_Ranges                     , only : Make_Range                                , rangeTypeLogarithmic
    implicit none
    type            (universeOperatorIntergalacticMediumStateEvolve)                        :: self
    class           (outputTimesClass                              ), intent(in   ), target :: outputTimes_
    class           (cosmologyParametersClass                      ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass                       ), intent(in   ), target :: cosmologyFunctions_
    class           (linearGrowthClass                             ), intent(in   ), target :: linearGrowth_
    class           (cosmologicalMassVarianceClass                 ), intent(in   ), target :: cosmologicalMassVariance_
    class           (radiationFieldClass                           ), intent(in   ), target :: radiationField_
    class           (gauntFactorClass                              ), intent(in   ), target :: gauntFactor_
    class           (atomicCrossSectionIonizationPhotoClass        ), intent(in   ), target :: atomicCrossSectionIonizationPhoto_
    class           (atomicIonizationPotentialClass                ), intent(in   ), target :: atomicIonizationPotential_
    class           (atomicRecombinationRateDielectronicClass      ), intent(in   ), target :: atomicRecombinationRateDielectronic_
    class           (atomicRecombinationRateRadiativeClass         ), intent(in   ), target :: atomicRecombinationRateRadiative_
    class           (atomicRecombinationRateRadiativeCoolingClass  ), intent(in   ), target :: atomicRecombinationRateRadiativeCooling_
    class           (atomicIonizationRateCollisionalClass          ), intent(in   ), target :: atomicIonizationRateCollisional_
    class           (atomicExcitationRateCollisionalClass          ), intent(in   ), target :: atomicExcitationRateCollisional_
    class           (intergalacticMediumStateClass                 ), intent(in   ), target :: intergalacticMediumState_
    type            (intergalacticMediumFilteringMassGnedin2000    )                        :: intergalacticMediumFilteringMass_
    integer                                                         , intent(in   )         :: timeCountPerDecade
    double precision                                                , intent(in   )         :: timeMinimum, timeMaximum
    double precision                                                , dimension(3)          :: massFilteringODEs
    integer                                                                                 :: iTime                                   , atomicNumber          , &
         &                                                                                     ionizationState
    double precision                                                                        :: atomicMass                              , density               , &
          &                                                                                    ionicFraction                           , massFractionPrimordial, &
          &                                                                                    massFilteringVarianceInitial
    !![
    <constructorAssign variables="timeMinimum, timeMaximum ,timeCountPerDecade, *cosmologyParameters_, *cosmologyFunctions_, *linearGrowth_, *cosmologicalMassVariance_, *outputTimes_, *gauntFactor_, *atomicCrossSectionIonizationPhoto_, *atomicIonizationPotential_, *atomicRecombinationRateDielectronic_, *atomicRecombinationRateRadiative_, *atomicRecombinationRateRadiativeCooling_, *atomicIonizationRateCollisional_, *atomicExcitationRateCollisional_, *intergalacticMediumState_, *radiationField_"/>
    !!]

    ! Get a RecFast intergalactic medium object for setting initial conditions.
    self%intergalacticMediumStateInitial=intergalacticMediumStateRecFast(cosmologyFunctions_,cosmologyParameters_)
    ! Build tables of properties and time for the temperature and ionization state densities in the IGM.
    self%timeCount=+int(                               &
         &              +dble(self%timeCountPerDecade) &
         &              *log10(                        &
         &                     +self%timeMaximum       &
         &                     /self%timeMinimum       &
         &                    )                        &
         &             )                               &
         &         +1
     ! Allocate arrays for all required IGM properties.
     allocate(self%temperature           (self%timeCount  ))
     allocate(self%densityHydrogen       (self%timeCount,2))
     allocate(self%densityHelium         (self%timeCount,3))
     allocate(self%time                  (self%timeCount  ))
     allocate(self%redshift              (self%timeCount  ))
     allocate(self%massFilteringComposite(self%timeCount,2))
     allocate(self%clumpingFactor        (self%timeCount  ))
     allocate(self%opticalDepth          (self%timeCount  ))
     allocate(self%massFiltering         (self%timeCount  ))
     ! Build grid of times.
     self%time=Make_Range(                      &
          &               self%timeMinimum    , &
          &               self%timeMaximum    , &
          &               self%timeCount      , &
          &               rangeTypeLogarithmic  &
          &              )
     ! Convert times to redshifts.
     do iTime=1,self%timeCount
        if (Values_Agree(self%time(iTime),self%outputTimes_%time(self%outputTimes_%count()),relTol=1.0d-6)) &
             & self%time(iTime)=self%outputTimes_%time(self%outputTimes_%count())
        self%redshift(iTime)                                                            &
             & =self%cosmologyFunctions_ %redshiftFromExpansionFactor(                  &
             &   self%cosmologyFunctions_%expansionFactor             (                 &
             &                                                         self%time(iTime) &
             &                                                        )                 &
             &                                                       )
     end do
     ! Initialize arrays to unphysical values.
     self%temperature           =-1.0d0
     self%densityHydrogen       =-1.0d0
     self%densityHelium         =-1.0d0
     self%massFilteringComposite=-1.0d0
     self%clumpingFactor        =-1.0d0
     self%opticalDepth          =-1.0d0
     self%massFiltering         =-1.0d0
     ! Initialize the temperature to that computed by RecFast at the initial time.
     self%temperature(1)=self%intergalacticMediumStateInitial%temperature(self%timeMinimum)
     ! Initialize the densities of ionic species using fractions from RecFast at the initial time.
     do atomicNumber=1,2
        select case (atomicNumber)
        case (1)
           massFractionPrimordial=hydrogenByMassPrimordial
           atomicMass            =atomicMassHydrogen
        case (2)
           massFractionPrimordial=heliumByMassPrimordial
           atomicMass            =atomicMassHelium
        end select
        do ionizationState=1,atomicNumber+1
           ionicFraction=0.0d0
           select case (atomicNumber)
           case (1)
              select case (ionizationState)
              case (1)
                 ionicFraction=self%intergalacticMediumStateInitial%neutralHydrogenFraction      (self%timeMinimum)
              case (2)
                 ionicFraction=self%intergalacticMediumStateInitial%singlyIonizedHydrogenFraction(self%timeMinimum)
              end select
           case (2)
              select case (ionizationState)
              case (1)
                 ionicFraction=self%intergalacticMediumStateInitial%neutralHeliumFraction        (self%timeMinimum)
              case (2)
                 ionicFraction=self%intergalacticMediumStateInitial%singlyIonizedHeliumFraction  (self%timeMinimum)
              case (3)
                 ionicFraction=self%intergalacticMediumStateInitial%doublyIonizedHeliumFraction  (self%timeMinimum)
              end select
           case default
              call Error_Report('unknown atomic number'//{introspection:location})
           end select
           density=+massFractionPrimordial                                    &
                &  *self%cosmologyParameters_%OmegaBaryon    (           )    &
                &  *self%cosmologyParameters_%densityCritical(           )    &
                &  /self%cosmologyFunctions_ %expansionFactor(timeMinimum)**3 &
                &  /atomicMassUnit                                            &
                &  /atomicMass                                                &
                &  *massSolar                                                 &
                &  /megaParsec                                            **3 &
                &  *ionicFraction
           select case (atomicNumber)
           case (1)
              self%densityHydrogen(1,ionizationState)=density
           case (2)
              self%densityHelium  (1,ionizationState)=density
           end select
        end do
     end do
     ! Set the composite variables used to solve for filtering mass.
     intergalacticMediumFilteringMass_=intergalacticMediumFilteringMassGnedin2000(.true.,self%cosmologyParameters_,self%cosmologyFunctions_,self%linearGrowth_,self%intergalacticMediumState_)
     call intergalacticMediumFilteringMass_%conditionsInitialODEs(self%timeMinimum,massFilteringODEs)
     self%massFilteringComposite(1,:)=massFilteringODEs(1:2)
     self%massFiltering         (1  )=massFilteringODEs(3  )
     ! Find the initial mass variance on the filtering mass scale.
     massFilteringVarianceInitial=self%cosmologicalMassVariance_%rootVariance(self%massFiltering(1),timeMinimum)
     ! Set the initial clumping factor.
     self%clumpingFactor(1)=1.0d0+massFilteringVarianceInitial**2
     ! Initialize optical depth.
     self%opticalDepth  (1)=0.0d0
     ! Initialize the IGM state object.
     call self%stateSet(1_c_size_t)
     return
   end function intergalacticMediumStateEvolveConstructorInternal

   subroutine intergalacticMediumStateEvolveDestructor(self)
     !!{
     Destructor for the \refClass{universeOperatorIntergalacticMediumStateEvolve} universeOperator class.
     !!}
     implicit none
     type(universeOperatorIntergalacticMediumStateEvolve), intent(inout) :: self

     !![
     <objectDestructor name="self%cosmologyParameters_"                    />
     <objectDestructor name="self%cosmologyFunctions_"                     />
     <objectDestructor name="self%linearGrowth_"                           />
     <objectDestructor name="self%cosmologicalMassVariance_"               />
     <objectDestructor name="self%outputTimes_"                            />
     <objectDestructor name="self%radiationField_"                         />
     <objectDestructor name="self%gauntFactor_"                            />
     <objectDestructor name="self%atomicCrossSectionIonizationPhoto_ "     />
     <objectDestructor name="self%atomicIonizationPotential_"              />
     <objectDestructor name="self%atomicRecombinationRateDielectronic_"    />
     <objectDestructor name="self%atomicRecombinationRateRadiative_"       />
     <objectDestructor name="self%atomicRecombinationRateRadiativeCooling_"/>
     <objectDestructor name="self%atomicIonizationRateCollisional_"        />
     <objectDestructor name="self%atomicExcitationRateCollisional_"        />
     <objectDestructor name="self%intergalacticMediumState_"               />
     !!]
     return
   end subroutine intergalacticMediumStateEvolveDestructor

   subroutine intergalacticMediumStateEvolveOperate(self,universe_)
     !!{
     Attach an initial event to a merger tree to cause the properties update function to be called.
     !!}
     use :: Galacticus_Nodes, only : universe, universeEvent
     implicit none
     class(universeOperatorIntergalacticMediumStateEvolve), intent(inout), target :: self
     type (universe                                      ), intent(inout)         :: universe_
     type (universeEvent                                 ), pointer               :: event

     ! Create the first interrupt event in the universe object.
     event         => universe_%createEvent( )
     event%time    =  self     %time       (1)
     event%creator => self
     event%task    => intergalacticMediumStateEvolveUpdate
     return
   end subroutine intergalacticMediumStateEvolveOperate

   logical function intergalacticMediumStateEvolveUpdate(event,universe_) result (success)
     !!{
     Update the properties for a given universe.
     !!}
     use            :: Arrays_Search           , only : searchArrayClosest
     use            :: Display                 , only : displayIndent     , displayMessage, displayUnindent
     use            :: Error                   , only : Error_Report
     use            :: Output_HDF5             , only : outputFile
     use            :: Galacticus_Nodes        , only : mergerTree        , mergerTreeList, nodeComponentBasic, treeNode, &
          &                                             universe          , universeEvent
     use            :: HDF5_Access             , only : hdf5Access
     use            :: IO_HDF5                 , only : hdf5Object
     use, intrinsic :: ISO_C_Binding           , only : c_size_t
     use            :: ISO_Varying_String      , only : varying_string
     use            :: Numerical_Constants_Math, only : Pi
     use            :: Numerical_ODE_Solvers   , only : odeSolver
     implicit none
     class           (universeEvent     ), intent(in   )            :: event
     type            (universe          ), intent(inout)            :: universe_
     type            (universeEvent     ), pointer                  :: newEvent
     double precision                    , parameter                :: odeToleranceAbsolute =1.0d-3     ,                 &
          &                                                            odeToleranceRelative =1.0d-3     ,                 &
          &                                                            timeToleranceRelative=1.0d-6
     type            (mergerTree        ), pointer                  :: tree
     type            (mergerTreeList    ), pointer                  :: forest
     type            (treeNode          ), pointer                  :: node
     class           (nodeComponentBasic), pointer                  :: basic
     integer         (c_size_t          ), parameter                :: propertyCount        =10_c_size_t
     double precision                    , parameter                :: massFilteringScale   =1.0d+4
     double precision                    , dimension(propertyCount) :: properties                       , propertyScales
     type            (odeSolver         )                           :: solver
     type            (varying_string    )                           :: message
     character       (len=6             )                           :: label
     type            (hdf5Object        )                           :: igmGroup                         , igmDataset
     integer         (c_size_t          )                           :: iNow
     double precision                                               :: treetimeLatest                   , timeCurrent    , &
          &                                                            timeMaximum

#ifdef USEMPI
     call Error_Report('intergalactic medium state evolver not implemented under MPI'//{introspection:location})
#endif
     select type (self => event%creator)
     class is (universeOperatorIntergalacticMediumStateEvolve)
        ! Display message.
        write (label,'(f6.3)') event%time
        message = "Evolving IGM properties to time "//trim(label)//" Gyr"
        call displayIndent(message)
        ! Find the current timestep.
        iNow = searchArrayClosest(self%time,event%time)
        ! Evolve the properties up to this timestep.
        if (iNow > 1) then
           ! Get required objects.
           self_ => self
           node_ => universe_%trees%tree%nodeBase
           ! Map properties to a contiguous array.
           properties( 1   )=self%temperature           (iNow-1    )
           properties( 2: 3)=self%densityHydrogen       (iNow-1,1:2)
           properties( 4: 6)=self%densityHelium         (iNow-1,1:3)
           properties( 7: 8)=self%massFilteringComposite(iNow-1,1:2)
           properties( 9   )=self%opticalDepth          (iNow-1    )
           properties(10   )=self%massFiltering         (iNow-1    )
           ! Set property scales.
           propertyScales( 1  )=        self%temperature           (iNow-1  )
           propertyScales( 2:3)=abs(sum(self%densityHydrogen       (iNow-1,:)))
           propertyScales( 4:6)=abs(sum(self%densityHelium         (iNow-1,:)))
           propertyScales( 7  )=+self%linearGrowth_%value(self%time(iNow-1  ))   &
                &               *(                                               &
                &                 +self%massFiltering(iNow-1)                    &
                &                 *3.0d0                                         &
                &                 /(                                             &
                &                   +4.0d0                                       &
                &                   *Pi                                          &
                &                   *self%cosmologyParameters_%OmegaMatter    () &
                &                   *self%cosmologyParameters_%densityCritical() &
                &                  )                                             &
                &                )**(2.0d0/3.0d0)                                &
                &               *Pi**2
           propertyScales( 8  )=+propertyScales(7)                                         &
                &               *self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime)
           propertyScales( 9  )=self%opticalDepth (iNow-1)
           propertyScales(10  )=self%massFiltering(iNow-1)
           ! Display message
           call displayMessage('Solving properties evolution')
           timeCurrent=self%time(iNow-1)
           solver     =odeSolver(propertyCount,intergalacticMediumStateEvolveODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative,scale=propertyScales)    
           call solver%solve(timeCurrent,self%time(iNow),properties)
           self%temperature           (iNow    )=max(properties( 1   ),0.0d0)
           self%densityHydrogen       (iNow,1:2)=max(properties( 2: 3),0.0d0)
           self%densityHelium         (iNow,1:3)=max(properties( 4: 6),0.0d0)
           self%massFilteringComposite(iNow,1:2)=    properties( 7: 8)
           self%opticalDepth          (iNow    )=max(properties( 9   ),0.0d0)
           self%massFiltering         (iNow    )=    properties(10   )
           ! Compute the filtering mass at this time.
           self%clumpingFactor        (iNow    )=+1.0d0                                                                                    &
                &                                +self%cosmologicalMassVariance_%rootVariance(self%massFiltering(iNow),self%time(iNow))**2
        end if
        ! Find the latest time across all trees in the universe.
        treeTimeLatest =  0.0d0
        forest         => universe_%trees
        do while (associated(forest))
           tree => forest%tree
           do while (associated(tree))
              node           => tree%nodeBase
              basic          => node%basic()
              treetimeLatest =  max(treetimeLatest,basic%time())
              tree           => tree%nextTree
           end do
           forest => forest%next
        end do
        ! Add the next event to the universe.
        timeMaximum=min(treetimeLatest,self%outputTimes_%time(self%outputTimes_%count()))
        if (iNow < self%timeCount .and. self%time(iNow+1) <= timeMaximum) then
           newEvent         => universe_%createEvent()
           newEvent%time    =  min(self%time(iNow+1),timeMaximum*(1.0d0-timeToleranceRelative))
           newEvent%creator => self
           newEvent%task    => intergalacticMediumStateEvolveUpdate
        else
           ! Output the results to file.
           !$ call hdf5Access%set()
           igmGroup=outputFile%openGroup('igmProperties', 'Properties of the intergalactic medium.')
           call igmGroup  %writeDataset  (self%redshift            ,'redshift'        ,'Redshift [].'                                  ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(0.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%temperature         ,'temperature'     ,'Temperature of the IGM [K].'                   ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%densityHydrogen(:,1),'densityHydrogen1','Density of H1 in the IGM [m⁻³].'               ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%densityHydrogen(:,2),'densityHydrogen2','Density of H2 in the IGM [m⁻³].'               ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%densityHelium  (:,1),'densityHelium1'  ,'Density of He1 in the IGM [m⁻³].'              ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%densityHelium  (:,2),'densityHelium2'  ,'Density of He2 in the IGM [m⁻³].'              ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%densityHelium  (:,3),'densityHelium3'  ,'Density of He3 in the IGM [m⁻³].'              ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%clumpingFactor      ,'clumpingFactor'  ,'Clumping factor in the IGM [].'                ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%opticalDepth        ,'opticalDepth'    ,'Electron scattering optical depth from z=0 [].',datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           call igmGroup  %writeDataset  (self%massFiltering       ,'filteringMass'   ,'Filtering mass in the IGM [M☉].'               ,datasetReturned=igmDataset)
           call igmDataset%writeAttribute(1.0d0                    ,'unitsInSI'                                                                                   )
           !$ call hdf5Access%unset()
        end if
        ! Store the past history to the default IGM state class.
        call self%stateSet(iNow)
        ! Display message.
        call displayUnindent('done')
        class default
        call Error_Report('incorrect class'//{introspection:location})
     end select
     ! Return true since we've performed our task.
     success=.true.
     return
   end function intergalacticMediumStateEvolveUpdate

   integer function intergalacticMediumStateEvolveODEs(time,properties,propertiesRateOfChange)
     !!{
     Evaluates the ODEs controlling the evolution temperature.
     !!}
     use :: Interface_GSL                        , only : GSL_Success
     use :: Intergalactic_Medium_Filtering_Masses, only : gnedin2000ODEs
     use :: Numerical_Constants_Astronomical     , only : gigaYear
     use :: Numerical_Constants_Atomic           , only : massHeliumAtom    , massHydrogenAtom
     use :: Numerical_Constants_Math             , only : Pi
     use :: Numerical_Constants_Physical         , only : boltzmannsConstant, electronMass     , electronRadius, fineStructure      , &
          &                                               plancksConstant   , radiationConstant, speedLight    , thomsonCrossSection
     use :: Numerical_Constants_Prefixes         , only : centi
     use :: Numerical_Constants_Units            , only : metersToAngstroms , electronVolt
     use :: Numerical_Integration                , only : integrator
     implicit none
     double precision                                          , intent(in  )                :: time
     double precision                                          , intent(in   ), dimension(:) :: properties
     double precision                                          , intent(  out), dimension(:) :: propertiesRateOfChange
     double precision                                          , parameter                   :: dielectronicRecombinationRateHeIEnergyLoss=40.74d0 ! electron volts.
     double precision                                          , parameter                   :: massFilteringMinimum                      =1.0d2
     double precision                                                         , dimension(2) :: densityHydrogen_                                  , massFilteringComposite_            , &
          &                                                                                     massFilteringCompositeRateOfChange
     double precision                                                         , dimension(3) :: densityHelium_                                    , massFilteringODEsRateOfChange      , &
          &                                                                                     massFilteringODEsProperties
     type            (integrator                              )                              :: integratorPhotoionization                         , integratorPhotoheating
     integer                                                                                 :: electronNumber                                    , atomicNumber                       , &
          &                                                                                     ionizationState                                   , shellNumber                        , &
          &                                                                                     photoionizationGroundIonizationState              , photoionizationGroundElectronNumber, &
          &                                                                                     iProperty
     double precision                                                                        :: temperature                                       , clumpingFactor                     , &
          &                                                                                     electronDensityRateOfChange                       , densityElectron                    , &
          &                                                                                     densityTotal                                      , ionizationPhotoRateFrom            , &
          &                                                                                     ionizationPhotoRateTo                             , opticalDepthRateOfChange           , &
          &                                                                                     massFilteringRateOfChange                         , opticalDepth                       , &
          &                                                                                     collisionIonizationRateFrom                       , collisionIonizationRateTo          , &
          &                                                                                     densityLowerIon                                   , densityUpperIon                    , &
          &                                                                                     densityThisIon                                    , recombinationDielectronicRateTo    , &
          &                                                                                     recombinationDielectronicRateFrom                 , recombinationRateTo                , &
          &                                                                                     recombinationRateFrom                             , wavelengthMinimum                  , &
          &                                                                                     wavelengthMaximum                                 , heatingRate                        , &
          &                                                                                     massParticleMean                                  , massFiltering_

     ! Extract properties from the contiguous array.
     temperature            =max(properties( 1   ),0.0d0               )
     densityHydrogen_       =    properties( 2: 3)
     densityHelium_         =    properties( 4: 6)
     massFilteringComposite_=    properties( 7: 8)
     opticalDepth           =    properties( 9   )
     massFiltering_         =max(properties(10   ),massFilteringMinimum)
     ! Compute total and electron densities.
     densityElectron        =+    max(      densityHydrogen_(2),0.0d0)  &
          &                  +    max(      densityHelium_  (2),0.0d0)  &
          &                  +    max(2.0d0*densityHelium_  (3),0.0d0)
     densityTotal           =+sum(max(      densityHydrogen_   ,0.0d0)) &
          &                  +sum(max(      densityHelium_     ,0.0d0)) &
          &                  +densityElectron
     ! Initialize heating rate to zero.
     heatingRate          =  0.0d0
     ! Compute rates of change of filtering mass composite parameters and optical depth.
     if (densityTotal > 0.0d0) then
        ! Find mean particle mass.
        massParticleMean=+(                                        &
             &             +massHydrogenAtom*sum(densityHydrogen_) &
             &             +massHeliumAtom  *sum(densityHelium_  ) &
             &             +electronMass    *    densityElectron   &
             &            )                                        &
             &           /densityTotal
        ! Evaluate optical depth term.
        opticalDepthRateOfChange=+speedLight                                                                               &
             &                   *thomsonCrossSection                                                                      &
             &                   *densityElectron                                                                          &
             &                   /self_%cosmologyFunctions_%expansionRate(self_%cosmologyFunctions_%expansionFactor(time)) &
             &                   *                                        self_%cosmologyFunctions_%expansionFactor(time)  &
             &                   *gigayear
     else
        massParticleMean        =0.0d0
        opticalDepthRateOfChange=0.0d0
     end if
     ! Evaluate the rates of change for the filtering mass variables.
     massFilteringODEsProperties       (1:2)=massFilteringComposite_
     massFilteringODEsProperties       (3  )=massFiltering_
     massFilteringODEsRateOfChange          =gnedin2000ODEs(self_%cosmologyParameters_,self_%cosmologyFunctions_,self_%linearGrowth_,time,massParticleMean,temperature,massFilteringODEsProperties)
     massFilteringCompositeRateOfChange     =massFilteringODEsRateOfChange(1:2)
     massFilteringRateOfChange              =massFilteringODEsRateOfChange(3  )
     ! Compute the clumping factor.
     clumpingFactor=+1.0d0                                                                &
          &         +self_%cosmologicalMassVariance_%rootVariance(massFiltering_,time)**2
     ! Build integrators.
     integratorPhotoionization=integrator(integrandPhotoionizationRate       ,toleranceRelative=1.0d-2)
     integratorPhotoheating   =integrator(integrandPhotoionizationHeatingRate,toleranceRelative=1.0d-3)
     ! Iterate over ionic species.
     iProperty=1 ! Counter for ionic species in properties array.
     do atomicNumber=1,2
        do ionizationState=1,atomicNumber+1
           ! Increment property array counter.
           iProperty=iProperty+1
           ! Determine electron number.
           electronNumber=atomicNumber+1-ionizationState
           ! Get density of this ionic state.
           densityThisIon    =max(properties(iProperty  ),0.0d0)
           ! Get density of upper ionic state (i.e. current ion minus one electron).
           if (ionizationState < atomicNumber+1) then
              densityUpperIon=max(properties(iProperty+1),0.0d0)
           else
              densityUpperIon=                            0.0d0
           end if
           ! Get density of lower ionic state (i.e. current ion plus one electron).
           if (ionizationState > 1) then
              densityLowerIon=max(properties(iProperty-1),0.0d0)
           else
              densityLowerIon=                            0.0d0
           end if
           ! Specify electron shell number. For H and He, this is always the 1s shell.
           shellNumber=1
           ! Compute collisional ionization rates from this ion.
           if (electronNumber  > 0) then
              collisionIonizationRateFrom=-self_%atomicIonizationRateCollisional_%rate     (atomicNumber,ionizationState  ,temperature) &
                   &                      *densityThisIon                                                                               &
                   &                      *densityElectron
              heatingRate                =+heatingRate                                                                                  &
                   &                      -self_%atomicIonizationPotential_      %potential(atomicNumber,electronNumber +1            ) &
                   &                      *electronVolt                                                                                 &
                   &                      *collisionIonizationRateFrom                                                                  &
                   &                      *gigaYear                                                                                     &
                   &                      *centi**3                                                                                     &
                   &                      *clumpingFactor
            else
              collisionIonizationRateFrom=+0.0d0
           end if
           ! Compute collisional ionization rates to this ion.
           if (ionizationState > 1) then
              collisionIonizationRateTo  =+self_%atomicIonizationRateCollisional_%rate     (atomicNumber,ionizationState-1,temperature) &
                   &                      *densityLowerIon                                                                              &
                   &                      *densityElectron
           else
              collisionIonizationRateTo  =+0.0d0
           end if
           ! Compute recombination rates from this ion.
           if (ionizationState > 1) then
              recombinationRateFrom      =-self_%atomicRecombinationRateRadiative_%rate(atomicNumber,ionizationState-1,temperature,recombinationCaseB) &
                   &                      *densityThisIon                                                                                              &
                   &                      *densityElectron
           else
              recombinationRateFrom      =+0.0d0
           end if
           ! Compute recombination rates to this ion.
           if (electronNumber  > 0) then
              recombinationRateTo        =+self_%atomicRecombinationRateRadiative_       %rate(atomicNumber,ionizationState,temperature,recombinationCaseB) &
                   &                      *densityUpperIon                                                                                                  &
                   &                      *densityElectron
              heatingRate                =+heatingRate                                                                                                      &
                   &                      -self_%atomicRecombinationRateRadiativeCooling_%rate(atomicNumber,ionizationState,temperature,recombinationCaseB) &
                   &                      *densityThisIon                                                                                                   &
                   &                      *densityElectron                                                                                                  &
                   &                      *gigaYear                                                                                                         &
                   &                      *centi**3                                                                                                         &
                   &                      *clumpingFactor                                                                                                   &
                   &                      *0.75d0                                                                                                           &
                   &                      *boltzmannsConstant                                                                                               &
                   &                      *temperature
           else
              recombinationRateTo        =+0.0d0
           end if
           ! Compute dielectronic recombination rates from this ion.
           if (ionizationState > 1) then
              recombinationDielectronicRateFrom=-self_%atomicRecombinationRateDielectronic_%rate(atomicNumber,electronNumber+1,temperature) &
                   &                            *densityThisIon                                                                             &
                   &                            *densityElectron
           else
              recombinationDielectronicRateFrom=+0.0d0
           end if
           ! Compute dielectronic recombination rates to this ion.
           if (electronNumber  > 0) then
              recombinationDielectronicRateTo  =+self_%atomicRecombinationRateDielectronic_%rate(atomicNumber,electronNumber  ,temperature) &
                   &                            *densityUpperIon                                                                            &
                   &                            *densityElectron
              heatingRate                      =+heatingRate                                                                                &
                   &                            -dielectronicRecombinationRateHeIEnergyLoss                                                 &
                   &                            *electronVolt                                                                               &
                   &                            *recombinationDielectronicRateFrom                                                          &
                   &                            *gigaYear                                                                                   &
                   &                            *centi**3                                                                                   &
                   &                            *clumpingFactor
           else
              recombinationDielectronicRateTo  =+0.0d0
           end if
           ! Set the epoch for the intergalactic background radiation field.
           call self_%radiationField_%timeSet(time)
           ! Compute rate of photoionizations from this ion.
           if (electronNumber  > 0) then
              ! Set the ground state for this photoionization calculation.
              photoionizationGroundIonizationState=ionizationState
              ! Set the minimum and maximum wavelengths for photoionization.
              wavelengthMinimum=+0.0d0
              wavelengthMaximum=+plancksConstant                                                         &
                   &            *speedLight                                                              &
                   &            /self_%atomicIonizationPotential_%potential(atomicNumber,electronNumber) &
                   &            /electronVolt                                                            &
                   &            *metersToAngstroms
              ! Integrate photoionizations over wavelength.
              ionizationPhotoRateFrom=-integratorPhotoionization%integrate(wavelengthMinimum,wavelengthMaximum) &
                   &                  *densityThisIon
           else
              ionizationPhotoRateFrom =+0.0d0
           end if
           if (ionizationState > 1) then
              ! Set the ground state for this photoionization calculation.
              photoionizationGroundIonizationState=ionizationState-1
              photoionizationGroundElectronNumber =electronNumber +1
              ! Set the minimum and maximum wavelengths for photoionization.
              wavelengthMinimum=0.0d0
              wavelengthMaximum=+plancksConstant                                                           &
                   &            *speedLight                                                                &
                   &            /self_%atomicIonizationPotential_%potential(atomicNumber,electronNumber+1) &
                   &            /electronVolt                                                              &
                   &            *metersToAngstroms
              ! Integrate photoionizations over wavelength.
              ionizationPhotoRateTo            =+integratorPhotoionization%integrate(wavelengthMinimum,wavelengthMaximum) &
                   &                            *densityLowerIon
              heatingRate                      =+heatingRate                                                              &
                   &                            +integratorPhotoheating   %integrate(wavelengthMinimum,wavelengthMaximum) &
                   &                            *densityLowerIon                                                          &
                   &                            *gigaYear
           else
              ionizationPhotoRateTo =+0.0d0
           end if
           ! Compute heating rate due to Bremsstrahlung.
           if (ionizationState > 1) then
              heatingRate=+heatingRate                                &
                   &      -16.0d0                                     &
                   &      / 3.0d0                                     &
                   &      *sqrt(                                      &
                   &            +2.0d0                                &
                   &            *Pi                                   &
                   &            /3.0d0                                &
                   &           )                                      &
                   &      *dble(ionizationState-1) **2                &
                   &      *densityThisIon                             &
                   &      *densityElectron                            &
                   &      *electronRadius          **3                &
                   &      *speedLight                                 &
                   &      /electronRadius                             &
                   &      *sqrt(                                      &
                   &            +electronMass                         &
                   &            *speedLight        **2                &
                   &            *boltzmannsConstant                   &
                   &            *temperature                          &
                   &           )                                      &
                   &      *fineStructure                              &
                   &      *self_%gauntFactor_%total(                  &
                   &                                atomicNumber    , &
                   &                                electronNumber+1, &
                   &                                temperature       &
                   &                               )                  &
                   &      *gigaYear                                   &
                   &      *clumpingFactor
           end if
           ! Add collisional excitation cooling rate.
           heatingRate=+heatingRate                                                        &
                &      -self_%atomicExcitationRateCollisional_%coolingRate(                &
                &                                                          atomicNumber  , &
                &                                                          electronNumber, &
                &                                                          temperature     &
                &                                                         )                &
                &      *clumpingFactor
           ! Compute net rate of change of density.
           propertiesRateOfChange(iProperty)=                                                   &
                ! Collisional ionization.
                & +(+collisionIonizationRateFrom      +collisionIonizationRateTo      )         &
                & *gigaYear                                                                     &
                & *centi**3                                                                     &
                & *clumpingFactor                                                               &
                ! Recombination.
                & +(+recombinationRateFrom            +recombinationRateTo            )         &
                & *gigaYear                                                                     &
                & *centi**3                                                                     &
                & *clumpingFactor                                                               &
                ! Dielectronic recombination.
                & +(+recombinationDielectronicRateFrom+recombinationDielectronicRateTo)         &
                & *gigaYear                                                                     &
                & *centi**3                                                                     &
                & *clumpingFactor                                                               &
                ! Photoionization.
                & +(+ionizationPhotoRateFrom          +ionizationPhotoRateTo          )         &
                & *gigaYear                                                                     &
                ! Cosmological expansion.
                & -3.0d0                                                                        &
                & *self_%cosmologyFunctions_%expansionRate(self_%cosmologyFunctions_%expansionFactor(time)) &
                & *densityThisIon
        end do
     end do
     ! Compute the rate of change of electron density.
     electronDensityRateOfChange=+propertiesRateOfChange(3) &
          &                      +propertiesRateOfChange(5) &
          &                      +propertiesRateOfChange(6)
     ! Compute rate of change of temperature due to cosmological expansion.
     propertiesRateOfChange(1)=-2.0d0                                           &
          &                    *self_%cosmologyFunctions_%expansionRate  (      &
          &                     self_%cosmologyFunctions_%expansionFactor (     &
          &                                                                time &
          &                                                               )     &
          &                                                              )      &
          &                    * temperature
     ! Compute rate of change of temperature due to atomic processes.
     if (densityTotal > 0.0d0)                                                            &
          & propertiesRateOfChange(1)=                                                    &
          &                    +propertiesRateOfChange(1)                                 &
          ! Accumulated atomic process heating rate.
          &                    +heatingRate                                               &
          &                    /1.5d0                                                     &
          &                    /boltzmannsConstant                                        &
          &                    /densityTotal                                              &
          ! CMB Compton scattering heating/cooling rate.
          &                   +speedLight                                                 &
          &                   *densityElectron                                            &
          &                   *4.0d0                                                      &
          &                   *thomsonCrossSection                                        &
          &                   *radiationConstant                                          &
          &                   *  self_%cosmologyFunctions_%temperatureCMBEpochal(time)**4 &
          &                   *(                                                          &
          &                     +self_%cosmologyFunctions_%temperatureCMBEpochal(time)    &
          &                     -temperature                                              &
          &                    )                                                          &
          &                   /electronMass                                               &
          &                   /speedLight                                             **2 &
          &                   /1.5d0                                                      &
          &                   /densityTotal                                               &
          &                   *gigaYear                                                   &
          ! Particle number rate of change.
          &    +electronDensityRateOfChange                                               &
          &    /densityTotal                                                              &
          &    +3.0d0                                                                     &
          &    *self_%cosmologyFunctions_%expansionRate  (                                &
          &     self_%cosmologyFunctions_%expansionFactor (                               &
          &                                                time                           &
          &                                               )                               &
          &                                              )
     ! Transfer rates of change to contiguous array.
     propertiesRateOfChange( 7: 8)=massFilteringCompositeRateOfChange
     propertiesRateOfChange( 9   )=opticalDepthRateOfChange
     propertiesRateOfChange(10   )=massFilteringRateOfChange
     ! Return success.
     intergalacticMediumStateEvolveODEs=GSL_Success

   contains

     double precision function integrandPhotoionizationRate(wavelength)
       !!{
       Integrand function used to compute the rate of photoionizations of an ionic species.
       !!}
       use :: Numerical_Constants_Units, only : ergs
       implicit none
       double precision, intent(in   ) :: wavelength
       double precision                :: photonDensity, photonFlux

       if (wavelength <= 0.0d0) then
          integrandPhotoionizationRate=0.0d0
       else
          photonFlux   =+self_%radiationField_%flux(wavelength,node_) &
               &        /centi**2                                     &
               &        *ergs
          photonDensity=+4.0d0                                        &
               &        *Pi                                           &
               &        *photonFlux                                   &
               &        /plancksConstant                              &
               &        /speedLight                                   &
               &        /wavelength
          integrandPhotoionizationRate=+speedLight                                                                                  &
               &                       *self_%atomicCrossSectionIonizationPhoto_%crossSection(                                      &
               &                                                                              atomicNumber                        , &
               &                                                                              photoionizationGroundIonizationState, &
               &                                                                              shellNumber                         , &
               &                                                                              wavelength                            &
               &                                                                             )                                      &
               &                       *centi**2                                                                                    &
               &                       *photonDensity
       end if
       return
     end function integrandPhotoionizationRate

     double precision function integrandPhotoionizationHeatingRate(wavelength)
       !!{
       Integrand function used to compute the rate of photoionization heating of an ionic species.
       !!}
       use :: Numerical_Constants_Units, only : ergs
       implicit none
       double precision, intent(in   ) :: wavelength
       double precision                :: photonDensity, photonFlux

       if (wavelength < 0.0d0) then
          integrandPhotoionizationHeatingRate=0.0d0
       else
          photonFlux   =+self_%radiationField_%flux(wavelength,node_) &
               &        /centi**2                                     &
               &        *ergs
          photonDensity=+4.0d0                                        &
               &        *Pi                                           &
               &        *photonFlux                                   &
               &        /plancksConstant                              &
               &        /speedLight                                   &
               &        /wavelength
          integrandPhotoionizationHeatingRate=+speedLight                                                                                  &
               &                              *self_%atomicCrossSectionIonizationPhoto_%crossSection(                                      &
               &                                                                                     atomicNumber                        , &
               &                                                                                     photoionizationGroundIonizationState, &
               &                                                                                     shellNumber                         , &
               &                                                                                     wavelength                            &
               &                                                                                    )                                      &
               &                              *centi**2                                                                                    &
               &                              *photonDensity                                                                               &
               &                              *(                                                                                           &
               &                                +plancksConstant                                                                           &
               &                                *speedLight                                                                                &
               &                                *metersToAngstroms                                                                         &
               &                                /wavelength                                                                                &
               &                                -self_%atomicIonizationPotential_      %potential   (                                      &
               &                                                                                     atomicNumber                        , &
               &                                                                                     photoionizationGroundElectronNumber   &
               &                                                                                    )                                      &
               &                                *electronVolt                                                                              &
               &                               )
       end if
       return
     end function integrandPhotoionizationHeatingRate

   end function intergalacticMediumStateEvolveODEs

   subroutine intergalacticMediumStateEvolveStateSet(self,iNow)
     use, intrinsic :: ISO_C_Binding, only : c_size_t
     implicit none
     class  (universeOperatorIntergalacticMediumStateEvolve), intent(inout) :: self
     integer(c_size_t                                      ), intent(in   ) :: iNow

     !![
     <eventHook name="intergalacticMediumStateEvolveUpdate">
      <interface>
        double precision, intent(in   ), dimension(:) :: time            , densityHydrogen1
        double precision, intent(in   ), dimension(:) :: densityHydrogen2, densityHelium1
        double precision, intent(in   ), dimension(:) :: densityHelium2  , densityHelium3
        double precision, intent(in   ), dimension(:) :: temperature     , massFiltering
      </interface>
      <callWith>self%time           (1:iNow  ), &amp;
        &amp;   self%densityHydrogen(1:iNow,1), &amp;
        &amp;   self%densityHydrogen(1:iNow,2), &amp;
        &amp;   self%densityHelium  (1:iNow,1), &amp;
        &amp;   self%densityHelium  (1:iNow,2), &amp;
        &amp;   self%densityHelium  (1:iNow,3), &amp;
        &amp;   self%temperature    (1:iNow  ), &amp;
        &amp;   self%massFiltering  (1:iNow  )</callWith>
     </eventHook>
     !!]
     return
   end subroutine intergalacticMediumStateEvolveStateSet
