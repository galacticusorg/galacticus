!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements the standard black hole node component.
!!}

module Node_Component_Black_Hole_Standard
  !!{
  Implement black hole tree node methods.
  !!}
  use :: Accretion_Disks                     , only : accretionDisksClass
  use :: Black_Hole_Binary_Initial_Separation, only : blackHoleBinaryInitialSeparationClass
  use :: Black_Hole_Binary_Mergers           , only : blackHoleBinaryMergerClass
  use :: Black_Hole_Binary_Recoil_Velocities , only : blackHoleBinaryRecoilClass
  use :: Black_Hole_Binary_Separations       , only : blackHoleBinarySeparationGrowthRateClass
  use :: Cooling_Radii                       , only : coolingRadiusClass
  use :: Cosmology_Parameters                , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleClass
  use :: Hot_Halo_Temperature_Profiles       , only : hotHaloTemperatureProfileClass
  implicit none
  private
  public :: Node_Component_Black_Hole_Standard_Rate_Compute       , Node_Component_Black_Hole_Standard_Scale_Set        , &
       &    Node_Component_Black_Hole_Standard_Thread_Uninitialize, Node_Component_Black_Hole_Standard_Output_Properties, &
       &    Node_Component_Black_Hole_Standard_Output_Names       , Node_Component_Black_Hole_Standard_Output_Count     , &
       &    Node_Component_Black_Hole_Standard_Output             , Node_Component_Black_Hole_Standard_Initialize       , &
       &    Node_Component_Black_Hole_Standard_Post_Evolve        , Node_Component_Black_Hole_Standard_Thread_Initialize

  !![
  <component>
   <class>blackHole</class>
   <name>standard</name>
   <isDefault>true</isDefault>
   <output instances="first"/>
   <properties>
    <property>
      <name>mass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <classDefault>defaultBlackHoleComponent%massSeed()</classDefault>
      <output unitsInSI="massSolar" comment="Mass of the black hole."/>
    </property>
    <property>
      <name>spin</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <getFunction>Node_Component_Black_Hole_Standard_Spin</getFunction>
      <classDefault>self%spinSeed()</classDefault>
      <output unitsInSI="0.0d0" comment="Spin of the black hole."/>
    </property>
    <property>
      <name>radialPosition</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>tripleInteractionTime</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>massSeed</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" isDeferred="get"  />
    </property>
    <property>
      <name>spinSeed</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <getFunction>Node_Component_Black_Hole_Standard_Seed_Spin</getFunction>
    </property>
    <property>
      <name>accretionRate</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>radiativeEfficiency</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
   </properties>
   <bindings>
    <binding method="enclosedMass" function="Node_Component_Black_Hole_Standard_Enclosed_Mass" bindsTo="component"/>
    <binding method="acceleration" function="Node_Component_Black_Hole_Standard_Acceleration"  bindsTo="component"/>
    <binding method="tidalTensor"  function="Node_Component_Black_Hole_Standard_Tidal_Tensor"  bindsTo="component"/>
   </bindings>
   <functions>objects.nodes.components.black_hole.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(cosmologyParametersClass                ), pointer :: cosmologyParameters_
  class(accretionDisksClass                     ), pointer :: accretionDisks_
  class(blackHoleBinaryRecoilClass              ), pointer :: blackHoleBinaryRecoil_
  class(blackHoleBinaryInitialSeparationClass   ), pointer :: blackHoleBinaryInitialSeparation_
  class(blackHoleBinaryMergerClass              ), pointer :: blackHoleBinaryMerger_
  class(blackHoleBinarySeparationGrowthRateClass), pointer :: blackHoleBinarySeparationGrowthRate_
  class(coolingRadiusClass                      ), pointer :: coolingRadius_
  class(hotHaloTemperatureProfileClass          ), pointer :: hotHaloTemperatureProfile_
  class(darkMatterHaloScaleClass                ), pointer :: darkMatterHaloScale_
  !$omp threadprivate(accretionDisks_,cosmologyParameters_,blackHoleBinaryRecoil_,blackHoleBinaryInitialSeparation_,blackHoleBinaryMerger_,blackHoleBinarySeparationGrowthRate_,coolingRadius_,hotHaloTemperatureProfile_,darkMatterHaloScale_)

  ! Accretion model parameters.
  ! Enhancement factors for the accretion rate.
  double precision :: bondiHoyleAccretionEnhancementHotHalo         , bondiHoyleAccretionEnhancementSpheroid
  ! Temperature of accreting gas.
  double precision :: bondiHoyleAccretionTemperatureSpheroid
  ! Control for hot mode only accretion.
  logical          :: bondiHoyleAccretionHotModeOnly

  ! Seed mass for black holes.
  double precision :: blackHoleSeedMass

  ! Feedback parameters.
  double precision :: blackHoleWindEfficiency                       , blackHoleRadioModeFeedbackEfficiency
  logical          :: blackHoleHeatsHotHalo                         , blackHoleWindEfficiencyScalesWithRadiativeEfficiency

  ! Output options.
  logical          :: blackHoleOutputAccretion
  logical          :: blackHoleOutputData
  logical          :: blackHoleOutputMergers

  ! Record of whether cold mode is explicitly tracked.
  logical          :: coldModeTracked

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Black_Hole_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Initialize(parameters_)
    !!{
    Initializes the standard black hole component module.
    !!}
    use :: Galacticus_Nodes, only : defaultHotHaloComponent, nodeComponentBlackHoleStandard
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type(inputParameters               ), intent(inout) :: parameters_
    type(nodeComponentBlackHoleStandard)                :: blackHoleStandardComponent

    ! Bind deferred functions.
    call blackHoleStandardComponent%massSeedFunction(Node_Component_Black_Hole_Standard_Seed_Mass)
    ! Get the seed mass
    !![
    <inputParameter>
      <name>blackHoleSeedMass</name>
      <source>parameters_</source>
      <defaultValue>100.0d0</defaultValue>
      <description>The mass of the seed black hole placed at the center of each newly formed galaxy.</description>
    </inputParameter>
    !!]
    ! Get accretion rate enhancement factors.
    !![
    <inputParameter>
      <name>bondiHoyleAccretionEnhancementSpheroid</name>
      <defaultValue>5.0d0</defaultValue>
      <description>The factor by which the Bondi-Hoyle accretion rate of spheroid gas onto black holes in enhanced.</description>
      <source>parameters_</source>
    </inputParameter>
    <inputParameter>
      <name>bondiHoyleAccretionEnhancementHotHalo</name>
      <defaultValue>6.0d0</defaultValue>
      <description>The factor by which the Bondi-Hoyle accretion rate of hot halo gas onto black holes in enhanced.</description>
      <source>parameters_</source>
    </inputParameter>
    <inputParameter>
      <name>bondiHoyleAccretionHotModeOnly</name>
      <defaultValue>.true.</defaultValue>
      <description>Determines whether accretion from the hot halo should only occur if the halo is in the hot accretion mode.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]

    ! Get temperature of accreting gas.
    !![
    <inputParameter>
      <name>bondiHoyleAccretionTemperatureSpheroid</name>
      <defaultValue>1.0d2</defaultValue>
      <description>The assumed temperature (in Kelvin) of gas in the spheroid when computing Bondi-Hoyle accretion rates onto black holes.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]

    ! Get wind efficiency and scaling.
    !![
    <inputParameter>
      <name>blackHoleWindEfficiency</name>
      <defaultValue>2.4d-3</defaultValue>
      <description>The efficiency of the black hole-driven wind: $L_\mathrm{wind} = \epsilon_\mathrm{wind} \dot{M}_\bullet \clight^2$.</description>
      <source>parameters_</source>
    </inputParameter>
    <inputParameter>
      <name>blackHoleWindEfficiencyScalesWithRadiativeEfficiency</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether the black hole wind efficiency should scale with the radiative efficiency of the accretion disk.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]

    ! Options controlling AGN feedback.
    !![
    <inputParameter>
      <name>blackHoleHeatsHotHalo</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not the black hole launched jets should heat the hot halo.</description>
      <source>parameters_</source>
    </inputParameter>
    <inputParameter>
      <name>blackHoleRadioModeFeedbackEfficiency</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Efficiency with which radio-mode feedback is coupled to the hot halo.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]

    ! Get options controlling output.
    !![
    <inputParameter>
      <name>blackHoleOutputAccretion</name>
      <defaultValue>.false.</defaultValue>
      <description>Determines whether or not accretion rates and jet powers will be output.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]

    ! Get options controlling output.
    !![
    <inputParameter>
      <name>blackHoleOutputData</name>
      <defaultValue>.false.</defaultValue>
      <description>Determines whether or not properties for all black holes (rather than just the central black hole) will be output.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]

    !![
    <inputParameter>
      <name>blackHoleOutputMergers</name>
      <defaultValue>.false.</defaultValue>
      <description>Determines whether or not properties of black hole mergers will be output.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]

    ! Check if cold mode is explicitly tracked.
    coldModeTracked=defaultHotHaloComponent%massColdIsGettable()
    ! Bind deferred functions.
    call blackHoleStandardComponent%      accretionRateFunction(Node_Component_Black_Hole_Standard_Accretion_Rate      )
    call blackHoleStandardComponent%radiativeEfficiencyFunction(Node_Component_Black_Hole_Standard_Radiative_Efficiency)
    return
  end subroutine Node_Component_Black_Hole_Standard_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Black_Hole_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node standard black hole module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent     , openMPThreadBindingAtLevel, dependencyRegEx, dependencyDirectionBefore
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    use :: Input_Parameters, only : inputParameter           , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies

    if (defaultBlackHoleComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionBefore,'^remnantStructure:')
       call satelliteMergerEvent%attach(defaultBlackHoleComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentBlackHoleStandard',dependencies=dependencies)
       !![
       <objectBuilder class="cosmologyParameters"                 name="cosmologyParameters_"                 source="parameters_"/>
       <objectBuilder class="accretionDisks"                      name="accretionDisks_"                      source="parameters_"/>
       <objectBuilder class="blackHoleBinaryRecoil"               name="blackHoleBinaryRecoil_"               source="parameters_"/>
       <objectBuilder class="blackHoleBinaryInitialSeparation"    name="blackHoleBinaryInitialSeparation_"    source="parameters_"/>
       <objectBuilder class="blackHoleBinaryMerger"               name="blackHoleBinaryMerger_"               source="parameters_"/>
       <objectBuilder class="blackHoleBinarySeparationGrowthRate" name="blackHoleBinarySeparationGrowthRate_" source="parameters_"/>
       <objectBuilder class="coolingRadius"                       name="coolingRadius_"                       source="parameters_"/>
       <objectBuilder class="hotHaloTemperatureProfile"           name="hotHaloTemperatureProfile_"           source="parameters_"/>
       <objectBuilder class="darkMatterHaloScale"                 name="darkMatterHaloScale_"                 source="parameters_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Black_Hole_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the tree node standard black hole module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    implicit none

    if (defaultBlackHoleComponent%standardIsActive()) then
       call satelliteMergerEvent%detach(defaultBlackHoleComponent,satelliteMerger)
       !![
       <objectDestructor name="cosmologyParameters_"                />
       <objectDestructor name="accretionDisks_"                     />
       <objectDestructor name="blackHoleBinaryRecoil_"              />
       <objectDestructor name="blackHoleBinaryInitialSeparation_"   />
       <objectDestructor name="blackHoleBinaryMerger_"              />
       <objectDestructor name="blackHoleBinarySeparationGrowthRate_"/>
       <objectDestructor name="coolingRadius_"                      />
       <objectDestructor name="hotHaloTemperatureProfile_"          />
       <objectDestructor name="darkMatterHaloScale_"                />
       !!]
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Thread_Uninitialize

  !![
  <rateComputeTask>
   <unitName>Node_Component_Black_Hole_Standard_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the black hole node mass rate of change.
    !!}
    use :: Galacticus_Nodes                , only : defaultBlackHoleComponent, interruptTask        , nodeComponentBasic  , nodeComponentBlackHole, &
          &                                         nodeComponentHotHalo     , nodeComponentSpheroid, propertyTypeInactive, treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear                 , megaParsec
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant       , speedLight
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    type            (treeNode              ), intent(inout)            :: node
    logical                                 , intent(inout)            :: interrupt
    procedure       (interruptTask         ), intent(inout), pointer   :: interruptProcedure
    integer                                 , intent(in   )            :: propertyType
    class           (nodeComponentBlackHole)               , pointer   :: blackHoleCentral                                                                                                                                   , blackHole
    class           (nodeComponentSpheroid )               , pointer   :: spheroid
    class           (nodeComponentHotHalo  )               , pointer   :: hotHalo
    class           (nodeComponentBasic    )               , pointer   :: basic
    double precision                                       , parameter :: windVelocity                =1.0d4                                                                                                                                                     !    Velocity of disk wind.
    double precision                                       , parameter :: ismTemperature              =1.0d4                                                                                                                                                     !    Temperature of the ISM.
    double precision                                       , parameter :: criticalDensityNormalization=2.0d0*massHydrogenAtom*speedLight**2*megaParsec/3.0d0/Pi/boltzmannsConstant/gigaYear/ismTemperature/kilo/windVelocity
    integer                                                            :: iInstance                                                                                                                                          , instanceCount
    double precision                                                   :: accretionRateHotHalo                                                                                                                               , accretionRateSpheroid                                          , &
         &                                                                criticalDensityRadius2                                                                                                                            , energyInputRate                                                , &
         &                                                                heatingRate                                                                                                                                       , jetEfficiency                                                  , &
         &                                                                massAccretionRate                                                                                                                                 , radiativeEfficiency                                            , &
         &                                                                restMassAccretionRate                                                                                                                             , spheroidDensityOverCriticalDensity                             , &
         &                                                                spheroidDensityRadius2                                                                                                                            , spheroidGasMass                                                , &
         &                                                                spheroidRadius                                                                                                                                    , windEfficiencyNet                                              , &
         &                                                                windFraction                                                                                                                                      , hotModeFraction

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    if (defaultBlackHoleComponent%standardIsActive()) then
       ! Get a count of the number of black holes associated with this node.
       instanceCount=node%blackHoleCount()
       ! Get the central black hole.
       blackHoleCentral => node%blackHole(instance=1)
       ! Get the basic, spheroid, and hot halo components.
       basic    => node%basic   ()
       spheroid => node%spheroid()
       hotHalo  => node%hotHalo ()
       ! Iterate over instances.
       do iInstance=1,max(instanceCount,1)
          ! Get the black hole.
          blackHole => node%blackHole(instance=iInstance)
          ! Find the rate of rest mass accretion onto the black hole.
          call Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(blackHole,accretionRateSpheroid,accretionRateHotHalo)
          restMassAccretionRate=accretionRateSpheroid+accretionRateHotHalo
          ! Finish if there is no accretion.
          if (restMassAccretionRate <= 0.0d0) cycle
          ! Find the radiative efficiency of the accretion.
          radiativeEfficiency=accretionDisks_%efficiencyRadiative(blackHole,restMassAccretionRate)
          ! Find the jet efficiency.
          if (restMassAccretionRate > 0.0d0) then
             jetEfficiency=accretionDisks_%powerJet(blackHole,restMassAccretionRate)/restMassAccretionRate/(speedLight&
                  &/kilo)**2
          else
             jetEfficiency=0.0d0
          end if
          ! Find the rate of increase in mass of the black hole.
          massAccretionRate=restMassAccretionRate*(1.0d0-radiativeEfficiency-jetEfficiency)
          ! If no black hole component currently exists and we have some accretion then interrupt and create a black hole.
          if (instanceCount == 0 .and. massAccretionRate /= 0.0d0) then
             interrupt=.true.
             interruptProcedure => Node_Component_Black_Hole_Standard_Create
             return
          end if
          ! Skip to the next black hole if this one has non-positive mass and a negative accretion rate.
          if (blackHole%mass() <= 0.0d0 .and. massAccretionRate < 0.0d0) cycle
          ! Add the accretion to the black hole.
          call blackHole%massRate       ( massAccretionRate                                 )
          ! Remove the accreted mass from the spheroid component.
          call spheroid %massGasSinkRate(-accretionRateSpheroid                             )
          ! Remove the accreted mass from the hot halo component.
          call hotHalo  %   massSinkRate(-accretionRateHotHalo ,interrupt,interruptProcedure)
          ! Set spin-up rate due to accretion.
          if (restMassAccretionRate > 0.0d0) call blackHole%spinRate(accretionDisks_%rateSpinUp(blackHole,restMassAccretionRate))
          ! Add heating to the hot halo component.
          if (blackHoleHeatsHotHalo) then
             ! Compute jet coupling efficiency based on whether halo is cooling quasistatically. Reduce this efficiency as the gas
             ! content in the halo drops below the cosmological mean.
             if (coldModeTracked) then
                hotModeFraction=1.0d0
             else
                hotModeFraction=Hot_Mode_Fraction(node)
             end if
             ! Get jet power.
             heatingRate=blackHoleRadioModeFeedbackEfficiency*jetEfficiency*restMassAccretionRate*(speedLight/kilo)**2
             ! Pipe this power to the hot halo.
             call hotHalo%heatSourceRate(heatingRate,interrupt,interruptProcedure)
          end if
          ! Add energy to the spheroid component.
          if (blackHoleWindEfficiency > 0.0d0) then
             spheroidGasMass=spheroid%massGas()
             if (spheroidGasMass > 0.0d0) then
                spheroidRadius=spheroid%radius()
                if (spheroidRadius > 0.0d0) then
                   spheroidDensityRadius2=3.0d0*spheroidGasMass/4.0d0/Pi/spheroidRadius
                   criticalDensityRadius2=criticalDensityNormalization*blackHoleWindEfficiency*restMassAccretionRate
                   ! Construct an interpolating factor such that the energy input from the wind drops to zero below half of the
                   ! critical density.
                   spheroidDensityOverCriticalDensity=spheroidDensityRadius2/criticalDensityRadius2-0.5d0
                   if (spheroidDensityOverCriticalDensity <= 0.0d0) then
                      ! No energy input below half of critical density.
                      windFraction=0.0d0
                   else if (spheroidDensityOverCriticalDensity >= 1.0d0) then
                      ! Full energy input above 1.5 times critical density.
                      windFraction=1.0d0
                   else
                      ! Smooth polynomial interpolating function between these limits.
                      windFraction=3.0d0*spheroidDensityOverCriticalDensity**2-2.0d0*spheroidDensityOverCriticalDensity**3
                   end if
                   ! Include scaling with radiative efficiency if requested,
                   windEfficiencyNet=windFraction*blackHoleWindEfficiency
                   if (blackHoleWindEfficiencyScalesWithRadiativeEfficiency) windEfficiencyNet=windEfficiencyNet*radiativeEfficiency
                   ! Compute the energy input and send it down the spheroid gas energy input pipe.
                   energyInputRate=windEfficiencyNet*restMassAccretionRate*(speedLight/kilo)**2
                   call spheroid%energyGasInputRate(energyInputRate)
                end if
             end if
          end if
       end do
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Rate_Compute

  !![
  <scaleSetTask>
   <unitName>Node_Component_Black_Hole_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent, nodeComponentBlackHole, nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    double precision                        , parameter              :: scaleMassRelative=1.0d-4
    double precision                        , parameter              :: scaleSizeRelative=1.0d-4
    double precision                        , parameter              :: scaleSizeAbsolute=1.0d-6
    double precision                        , parameter              :: scaleMassAbsolute=1.0d+0
    class           (nodeComponentSpheroid )               , pointer :: spheroid
    class           (nodeComponentBlackHole)               , pointer :: blackHole
    integer                                                          :: instance

    ! Determine if the standard implementation is active and at least one black hole exists.
    if (defaultBlackHoleComponent%standardIsActive().and.node%blackHoleCount() > 0) then
       ! Get the spheroid component.
       spheroid => node%spheroid()
       ! Loop over instances.
       do instance=1,node%blackHoleCount()
          ! Get the black hole.
          blackHole => node%blackHole(instance=instance)
          ! Set scale for mass.
          call blackHole%massScale(                                                       &
               &                   max(                                                   &
               &                                                 blackHole%massSeed   (), &
               &                       max(                                               &
               &                               scaleMassRelative*spheroid %massStellar(), &
               &                           max(                                           &
               &                               scaleMassAbsolute                        , &
               &                                                 blackHole%mass       ()  &
               &                              )                                           &
               &                          )                                               &
               &                      )                                                   &
               &                  )

          ! Set scale for spin.
          call blackHole%spinScale(1.0d0)

          ! Set scale for radius.
          call blackHole%radialPositionScale(                                                    &
               &                             maxval(                                             &
               &                                  [                                              &
               &                                   scaleSizeAbsolute,                            &
               &                                   scaleSizeRelative*spheroid %halfMassRadius(), &
               &                                                     blackHole%radialPosition()  &
               &                                  ]                                              &
               &                                   )                                             &
               &                            )
       end do
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Scale_Set

  subroutine satelliteMerger(self,node)
    !!{
    Merge any black hole associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode
    implicit none
    class           (*                     ), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    type            (treeNode              ), pointer       :: hostNode
    class           (nodeComponentBlackHole), pointer       :: blackHoleHostCentral, blackHole         , &
         &                                                     blackHolePrimary    , blackHoleSecondary
    integer                                                 :: instance
    double precision                                        :: blackHoleMassNew    , blackHoleSpinNew  , &
         &                                                     massBlackHole1      , massBlackHole2    , &
         &                                                     radiusInitial       , recoilVelocity    , &
         &                                                     spinBlackHole1      , spinBlackHole2
    !$GLC attributes unused :: self
    
    ! Find the node to merge with.
    hostNode => node%mergesWith()
    ! Find the initial radius of the satellite black hole in the remnant.
    radiusInitial=blackHoleBinaryInitialSeparation_%separationInitial(node,hostNode)
    ! If the separation is non-positive, assume that the black holes merge instantaneously.
    if (radiusInitial <= 0.0d0) then
       ! Get the central black hole of the host galaxy.
       blackHoleHostCentral => hostNode%blackHole(instance=1,autoCreate=.true.)
       ! Loop over all black holes in the satellite galaxy.
       do instance=1,node%blackHoleCount()
          ! Get the black hole.
          blackHole => node%blackHole(instance=instance)
          ! Compute the outcome of the merger.
          call blackHoleBinaryMerger_%merge(                             &
               &                            blackHole           %mass(), &
               &                            blackHoleHostCentral%mass(), &
               &                            blackHole           %spin(), &
               &                            blackHoleHostCentral%spin(), &
               &                            blackHoleMassNew           , &
               &                            blackHoleSpinNew             &
               &                           )
          ! Merge the black holes instantaneously.
          ! Check which black hole is more massive in order to compute an appropriate recoil velocity
          if (blackHoleHostCentral%mass() >= blackHole%mass()) then
             blackHolePrimary   => blackHoleHostCentral
             blackHoleSecondary => blackHole
          else
             blackHolePrimary   => blackHole
             blackHoleSecondary => blackHoleHostCentral
          end if
          massBlackHole1=blackHolePrimary  %mass()
          massBlackHole2=blackHoleSecondary%mass()
          spinBlackHole1=blackHolePrimary  %spin()
          spinBlackHole2=blackHoleSecondary%spin()
          ! Now calculate the recoil velocity of the binary black hole and check wether it escapes the galaxy.
          recoilVelocity=blackHoleBinaryRecoil_%velocity(blackHolePrimary,blackHoleSecondary)
          if (Node_Component_Black_Hole_Standard_Recoil_Escapes(node,recoilVelocity,radius=0.0d0,ignoreCentralBlackHole=.true.)) then
             blackHoleMassNew=0.0d0
             blackHoleSpinNew=0.0d0
          end if
          ! Move the black hole to the host.
          call Node_Component_Black_Hole_Standard_Output_Merger(node,massBlackHole1,massBlackHole2)
          call blackHoleHostCentral%massSet(blackHoleMassNew           )
          call blackHoleHostCentral%spinSet(blackHoleSpinNew           )
          ! Reset the satellite black hole to zero mass.
          call blackHole           %massSet(blackHole       %massSeed())
          call blackHole           %spinSet(blackHole       %spinSeed())
       end do
    else
       ! Adjust the radii of the black holes in the satellite galaxy.
       do instance=node%blackHoleCount(),1,-1
          blackHole => node%blackHole(instance=instance)
          call blackHole%radialPositionSet(radiusInitial)
          ! Declares them as not having interacted in a triple black hole interaction.
          call blackHole%tripleInteractionTimeSet(0.0d0)
          ! Remove this black hole if it has no mass.
          if (blackHole%mass() <= 0.0d0) call node%blackHoleRemove(instance)
       end do
       ! Move black holes from the satellite to the host.
       call node%blackHoleMove(hostNode)
    end if
    return
  end subroutine satelliteMerger

  logical function Node_Component_Black_Hole_Standard_Recoil_Escapes(node,recoilVelocity,radius,ignoreCentralBlackHole)
    !!{
    Return true if the given recoil velocity is sufficient to eject a black hole from the halo.
    !!}
    use :: Galactic_Structure_Options   , only : componentTypeBlackHole
    use :: Galactic_Structure_Potentials, only : Galactic_Structure_Potential
    use :: Galacticus_Nodes             , only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: recoilVelocity        , radius
    logical                   , intent(in   ) :: ignoreCentralBlackHole
    double precision                          :: potentialCentral      , potentialCentralSelf, &
         &                                       potentialHalo         , potentialHaloSelf

    ! Compute relevant potentials.
    potentialCentral=Galactic_Structure_Potential(node,radius                                                                      )
    potentialHalo   =Galactic_Structure_Potential(node,darkMatterHaloScale_%virialRadius(node)                                     )
    if (ignoreCentralBlackHole) then
       ! Compute potential of central black hole to be subtracted off of total value.
       potentialCentralSelf=Galactic_Structure_Potential(node,radius                                 ,componentType=componentTypeBlackHole)
       potentialHaloSelf   =Galactic_Structure_Potential(node,darkMatterHaloScale_%virialRadius(node),componentType=componentTypeBlackHole)
    else
       ! No correction for central black hole as it is to be included.
       potentialCentralSelf=0.0d0
       potentialHaloSelf   =0.0d0
    end if
    ! Evaluate the escape condition.
    Node_Component_Black_Hole_Standard_Recoil_Escapes= &
         &  +0.5d0*recoilVelocity      **2             &
         &  +      potentialCentral                    &
         &  -      potentialCentralSelf                &
         & >                                           &
         &  +      potentialHalo                       &
         &  -      potentialHaloSelf
    return
  end function Node_Component_Black_Hole_Standard_Recoil_Escapes

  subroutine Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(blackHole,accretionRateSpheroid,accretionRateHotHalo)
    !!{
    Returns the rate of mass accretion onto the black hole in {\normalfont \ttfamily node}.
    !!}
    use :: Black_Hole_Fundamentals         , only : Black_Hole_Eddington_Accretion_Rate
    use :: Bondi_Hoyle_Lyttleton_Accretion , only : Bondi_Hoyle_Lyttleton_Accretion_Radius, Bondi_Hoyle_Lyttleton_Accretion_Rate
    use :: Galactic_Structure_Densities    , only : Galactic_Structure_Density
    use :: Galactic_Structure_Options      , only : componentTypeColdHalo                 , componentTypeHotHalo                , componentTypeSpheroid, coordinateSystemCylindrical, &
          &                                         massTypeGaseous
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole                , nodeComponentHotHalo                , nodeComponentSpheroid, treeNode
    use :: Ideal_Gases_Thermodynamics      , only : Ideal_Gas_Jeans_Length                , Ideal_Gas_Sound_Speed
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr               , gigaYear                            , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (nodeComponentBlackHole), intent(inout)          :: blackHole
    double precision                        , intent(  out)          :: accretionRateHotHalo      , accretionRateSpheroid
    type            (treeNode              )               , pointer :: node
    class           (nodeComponentSpheroid )               , pointer :: spheroid
    class           (nodeComponentHotHalo  )               , pointer :: hotHalo
    double precision                        , parameter              :: gasDensityMinimum   =1.0d0                              ! Lowest gas density to consider when computing accretion rates onto black hole (in units of M_Solar/Mpc^3).
    double precision                                                 :: accretionRadius           , accretionRateMaximum    , &
         &                                                              blackHoleMass             , gasDensity              , &
         &                                                              hotHaloTemperature        , hotModeFraction         , &
         &                                                              jeansLength               , position             (3), &
         &                                                              radiativeEfficiency       , relativeVelocity        , &
         &                                                              coldModeFraction

    ! Get the host node.
    node => blackHole%host()
    ! Get black hole mass.
    blackHoleMass=blackHole%mass()
    ! Check black hole mass is positive.
    if (blackHoleMass > 0.0d0) then
       ! Compute the relative velocity of black hole and gas. We assume that relative motion arises only from the radial
       ! migration of the black hole.
       relativeVelocity=blackHoleBinarySeparationGrowthRate_%growthRate(blackHole)*Mpc_per_km_per_s_To_Gyr
       ! Contribution from spheroid:
       ! Get the accretion radius. We take this to be the larger of the Bondi-Hoyle radius and the current radius position of
       ! the black hole.
       accretionRadius=max(                                                                                              &
            &               Bondi_Hoyle_Lyttleton_Accretion_Radius(blackHoleMass,bondiHoyleAccretionTemperatureSpheroid) &
            &              ,blackHole%radialPosition()                                                      &
            &             )
       ! Set the position.
       position=[accretionRadius,0.0d0,0.0d0]
       ! Get density of gas at the galactic center.
       gasDensity=Galactic_Structure_Density(node,position,coordinateSystem=coordinateSystemCylindrical,componentType&
            &=componentTypeSpheroid,massType =massTypeGaseous)
       ! Check if we have a non-negligible gas density.
       if (gasDensity > gasDensityMinimum) then
          ! Get the spheroid component.
          spheroid => node%spheroid()
          ! Get the Jeans length scale.
          jeansLength=Ideal_Gas_Jeans_Length(bondiHoyleAccretionTemperatureSpheroid,gasDensity)
          ! Limit the smoothing scale to the scale of the spheroid.
          jeansLength=min(jeansLength,spheroid%radius())
          ! If the Jeans length exceeds the Bondi-Hoyle-Lyttleton accretion radius, then recompute gas density for a larger
          ! radius, as the gas should be smoothly distributed on scales below the Jeans length.
          if (jeansLength > accretionRadius) then
             ! Set the position.
             position=[jeansLength,0.0d0,0.0d0]
             ! Get density of gas at the galactic center.
             gasDensity=Galactic_Structure_Density(node,position,coordinateSystem=coordinateSystemCylindrical,componentType&
                  &=componentTypeSpheroid,massType =massTypeGaseous)
          end if
          ! Compute the accretion rate.
          accretionRateSpheroid=max(bondiHoyleAccretionEnhancementSpheroid*Bondi_Hoyle_Lyttleton_Accretion_Rate(blackHoleMass&
               &,gasDensity ,relativeVelocity,bondiHoyleAccretionTemperatureSpheroid),0.0d0)
          ! Get the radiative efficiency of the accretion.
          radiativeEfficiency=accretionDisks_%efficiencyRadiative(blackHole,accretionRateSpheroid)
          ! Limit the accretion rate to the Eddington limit.
          if (radiativeEfficiency > 0.0d0) accretionRateSpheroid=min(accretionRateSpheroid&
               &,Black_Hole_Eddington_Accretion_Rate(blackHole) /radiativeEfficiency)
       else
          ! Gas density is negative - set zero accretion rate.
          accretionRateSpheroid=0.0d0
       end if
       ! Contribution from hot halo:
       ! Get the hot halo component.
       hotHalo => node%hotHalo()
       ! Get halo gas temperature.
       hotHaloTemperature=hotHaloTemperatureProfile_%temperature(node,radius=0.0d0)
       ! Get the accretion radius.
       accretionRadius=Bondi_Hoyle_Lyttleton_Accretion_Radius(blackHoleMass,hotHaloTemperature)
       accretionRadius=min(accretionRadius,hotHalo%outerRadius())
       ! Set the position.
       position=[accretionRadius,0.0d0,0.0d0]
       ! Find the fraction of gas in the halo which is in the hot mode. Set this to unity if hot/cold mode is not to be
       ! considered.
       select case (bondiHoyleAccretionHotModeOnly)
       case (.true.)
          if (coldModeTracked) then
             hotModeFraction=1.0d0
          else
             hotModeFraction=Hot_Mode_Fraction(node)
          end if
          coldModeFraction=0.0d0
       case (.false.)
          hotModeFraction=1.0d0
          if (coldModeTracked) then
             coldModeFraction=1.0d0
          else
             coldModeFraction=0.0d0
          end if
       end select
       ! Get density of gas at the galactic center - scaled by the fraction in the hot accretion mode.
       gasDensity=                                                                                 &
            &      hotModeFraction                                                                 &
            &     *Galactic_Structure_Density(                                                     &
            &                                 node                                           , &
            &                                 position                                           , &
            &                                 coordinateSystem=coordinateSystemCylindrical       , &
            &                                 componentType   =componentTypeHotHalo              , &
            &                                 massType        =massTypeGaseous                     &
            &                                )
       if (coldModeTracked.and.coldModeFraction > 0.0d0)                                           &
            & gasDensity=                                                                          &
            &             gasDensity                                                               &
            &            +coldModeFraction                                                         &
            &            *Galactic_Structure_Density(                                              &
            &                                        node                                    , &
            &                                        position                                    , &
            &                                        coordinateSystem=coordinateSystemCylindrical, &
            &                                        componentType   =componentTypeColdHalo      , &
            &                                        massType        =massTypeGaseous              &
            &                                       )
       ! Check if we have a non-zero gas density.
       if (gasDensity > gasDensityMinimum) then
          ! Compute the accretion rate.
          accretionRateHotHalo=max(bondiHoyleAccretionEnhancementHotHalo*Bondi_Hoyle_Lyttleton_Accretion_Rate(blackHoleMass&
               &,gasDensity,relativeVelocity,hotHaloTemperature,accretionRadius),0.0d0)
          ! Limit the accretion rate to the total mass of the hot halo, divided by the sound crossing time.
          accretionRateMaximum=max(hotHalo%mass()/(hotHalo%outerRadius()/(kilo*gigaYear/megaParsec)&
               &/Ideal_Gas_Sound_Speed(hotHaloTemperature)),0.0d0)
          accretionRateHotHalo=min(accretionRateHotHalo,accretionRateMaximum)
          ! Get the radiative efficiency of the accretion.
          radiativeEfficiency=accretionDisks_%efficiencyRadiative(blackHole,accretionRateHotHalo)
          ! Limit the accretion rate to the Eddington limit.
          if (radiativeEfficiency > 0.0d0) accretionRateHotHalo=min(accretionRateHotHalo&
               &,Black_Hole_Eddington_Accretion_Rate(blackHole)/radiativeEfficiency)
       else
          ! No gas density, so zero accretion rate.
          accretionRateHotHalo=0.0d0
       end if
    else
       accretionRateSpheroid=0.0d0
       accretionRateHotHalo =0.0d0
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Mass_Accretion_Rate

  subroutine Node_Component_Black_Hole_Standard_Create(node)
    !!{
    Creates a black hole component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode
    implicit none
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentBlackHole)               , pointer :: blackHole

    ! Create the black hole.
    blackHole => node%blackHole(autoCreate=.true.)
    ! Set to the seed mass.
    call blackHole%          massSet(blackHole%massSeed())
    call blackHole%          spinSet(blackHole%spinSeed())
    call blackHole%radialPositionSet(               0.0d0)
    return
  end subroutine Node_Component_Black_Hole_Standard_Create

  !![
  <mergerTreeOutputNames>
   <unitName>Node_Component_Black_Hole_Standard_Output_Names</unitName>
   <sortName>Node_Component_Black_Hole_Standard_Output</sortName>
  </mergerTreeOutputNames>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Output_Names(node,integerProperty,integerProperties,doubleProperty,doubleProperties,time)
    !!{
    Set names of black hole properties to be written to the \glc\ output file.
    !!}
    use :: Galacticus_Nodes                  , only : treeNode
    use :: Numerical_Constants_Astronomical  , only : gigaYear             , massSolar
    use :: Numerical_Constants_Prefixes      , only : kilo
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
    implicit none
    type            (treeNode             )              , intent(inout) :: node
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleProperty   , integerProperty
    type            (outputPropertyInteger), dimension(:), intent(inout) :: integerProperties
    type            (outputPropertyDouble ), dimension(:), intent(inout) :: doubleProperties
    !$GLC attributes unused :: time

    if (Node_Component_Black_Hole_Standard_Matches(node)) then
       integerProperty=integerProperty+1
       integerProperties(integerProperty)%name     ='blackHoleCount'
       integerProperties(integerProperty)%comment  ='Number of super-massive black holes in the galaxy.'
       integerProperties(integerProperty)%unitsInSI=0.0d0
       if (blackHoleOutputAccretion) then
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%name     ='blackHoleAccretionRate'
          doubleProperties(doubleProperty)%comment  ='Rest-mass accretion rate onto the black hole.'
          doubleProperties(doubleProperty)%unitsInSI=massSolar/gigaYear
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%name     ='blackHoleJetPower'
          doubleProperties(doubleProperty)%comment  ='Power of the black hole-driven jet.'
          doubleProperties(doubleProperty)%unitsInSI=massSolar*kilo**2/gigaYear
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%name     ='blackHoleRadiativeEfficiency'
          doubleProperties(doubleProperty)%comment  ='The radiative efficiency of the black hole accretion system.'
          doubleProperties(doubleProperty)%unitsInSI=0.0d0
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Names

  !![
  <mergerTreeOutputPropertyCount>
   <unitName>Node_Component_Black_Hole_Standard_Output_Count</unitName>
   <sortName>Node_Component_Black_Hole_Standard_Output</sortName>
  </mergerTreeOutputPropertyCount>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Output_Count(node,integerPropertyCount,doublePropertyCount,time)
    !!{
    Account for the number of black hole properties to be written to the the \glc\ output file.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount  , integerPropertyCount
    integer                   , parameter     :: extraPropertyCount =3
    !$GLC attributes unused :: time

    if (Node_Component_Black_Hole_Standard_Matches(node)) then
       integerPropertyCount=integerPropertyCount+1
       if (blackHoleOutputAccretion) doublePropertyCount=doublePropertyCount+extraPropertyCount
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Count

  !![
  <mergerTreeOutputTask>
   <unitName>Node_Component_Black_Hole_Standard_Output</unitName>
   <sortName>Node_Component_Black_Hole_Standard_Output</sortName>
  </mergerTreeOutputTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Output(node,integerProperty,integerBufferCount,integerProperties,doubleProperty,doubleBufferCount,doubleProperties,time,instance)
    !!{
    Store black hole properties in the \glc\ output file buffers.
    !!}
    use :: Galacticus_Nodes                  , only : nodeComponentBlackHole, treeNode
    use :: Kind_Numbers                      , only : kind_int8
    use :: Multi_Counters                    , only : multiCounter
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger , outputPropertyDouble
    implicit none
    double precision                        , intent(in   )               :: time
    type            (treeNode              ), intent(inout)               :: node
    integer                                 , intent(inout)               :: doubleBufferCount    , doubleProperty       , &
         &                                                                   integerBufferCount   , integerProperty
    type            (outputPropertyInteger ), intent(inout), dimension(:) :: integerProperties
    type            (outputPropertyDouble  ), intent(inout), dimension(:) :: doubleProperties
    type            (multiCounter          ), intent(inout)               :: instance
    class           (nodeComponentBlackHole)               , pointer      :: blackHole
    double precision                                                      :: accretionRateHotHalo , accretionRateSpheroid, &
         &                                                                   restMassAccretionRate
    !$GLC attributes unused :: time, instance

    if (Node_Component_Black_Hole_Standard_Matches(node)) then
       ! Store the properties.
       if (blackHoleOutputAccretion) then
          ! Get the black hole component.
          blackHole => node%blackHole(instance=1)
          ! Get the rest mass accretion rate.
          call Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(blackHole,accretionRateSpheroid,accretionRateHotHalo)
          restMassAccretionRate=accretionRateSpheroid+accretionRateHotHalo
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%scalar(doubleBufferCount)=restMassAccretionRate
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%scalar(doubleBufferCount)=accretionDisks_%powerJet           (blackHole,restMassAccretionRate)
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%scalar(doubleBufferCount)=accretionDisks_%efficiencyRadiative(blackHole,restMassAccretionRate)
       end if
       ! Count number of black holes associated with this galaxy.
       integerProperty=integerProperty+1
       integerProperties(integerProperty)%scalar(integerBufferCount)=node%blackHoleCount()
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Output

  logical function Node_Component_Black_Hole_Standard_Matches(node)
    !!{
    Return true if the black hole component of {\normalfont \ttfamily node} is a match to the standard implementation.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent, nodeComponentBlackHole, nodeComponentBlackHoleStandard, treeNode
    implicit none
    type (treeNode              ), intent(inout) :: node
    class(nodeComponentBlackHole), pointer       :: blackHole

    ! Get the black hole component.
    blackHole => node%blackHole()
    ! Ensure that it is of the standard class.
    Node_Component_Black_Hole_Standard_Matches=.false.
    select type (blackHole)
    class is (nodeComponentBlackHoleStandard)
       Node_Component_Black_Hole_Standard_Matches=.true.
    type  is (nodeComponentBlackHole        )
       Node_Component_Black_Hole_Standard_Matches=defaultBlackHoleComponent%standardIsActive()
    end select
    return
  end function Node_Component_Black_Hole_Standard_Matches

  subroutine Node_Component_Black_Hole_Standard_Output_Merger(node,massBlackHole1,massBlackHole2)
    !!{
    Outputs properties of merging black holes.
    !!}
    use :: Galacticus_HDF5 , only : galacticusOutputFile
    use :: Galacticus_Nodes, only : nodeComponentBasic  , treeNode
    use :: IO_HDF5         , only : hdf5Access          , hdf5Object
    implicit none
    type            (treeNode          ), intent(inout) :: node
    double precision                    , intent(in   ) :: massBlackHole1    , massBlackHole2
    class           (nodeComponentBasic), pointer       :: basic
    type            (hdf5Object        )                :: mergersGroup

    ! Exit if merger data is not to be output.
    if (.not.blackHoleOutputMergers) return

    ! Ignore mergers with zero mass black holes.
    if (massBlackHole2 <= 0.0d0    ) return

    ! Get the basic component.
    basic => node%basic()

    ! Open the group to which black hole mergers should be written.
    !$ call hdf5Access%set()
    mergersGroup=galacticusOutputFile%openGroup("blackHoleMergers","Black hole mergers data.")
    ! Append to the datasets.
    call mergersGroup%writeDataset([massBlackHole1            ],"massBlackHole1","Mass of the first merging black hole." ,appendTo=.true.)
    call mergersGroup%writeDataset([massBlackHole2            ],"massBlackHole2","Mass of the second merging black hole.",appendTo=.true.)
    call mergersGroup%writeDataset([basic%time()              ],"timeOfMerger"  ,"The time of the black hole merger."    ,appendTo=.true.)
    call mergersGroup%writeDataset([node%hostTree%volumeWeight],"volumeWeight"  ,"The weight for the black hole merger." ,appendTo=.true.)
    ! Close the group.
    call mergersGroup%close()
    !$ call hdf5Access%unset()
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Merger

  !![
  <mergerTreeExtraOutputTask>
   <unitName>Node_Component_Black_Hole_Standard_Output_Properties</unitName>
  </mergerTreeExtraOutputTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Output_Properties(node,iOutput,treeIndex,nodePassesFilter)
    !!{
    Output properties for all black holes in {\normalfont \ttfamily node}.
    !!}
    use            :: Galacticus_HDF5   , only : galacticusOutputFile
    use            :: Galacticus_Nodes  , only : nodeComponentBlackHole, treeNode
    use            :: IO_HDF5           , only : hdf5Access            , hdf5Object
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : assignment(=)         , char           , varying_string
    use            :: Kind_Numbers      , only : kind_int8
    use            :: Memory_Management , only : allocateArray         , deallocateArray
    use            :: String_Handling   , only : operator(//)
    implicit none
    type            (treeNode              ), intent(inout), pointer      :: node
    integer         (kind=kind_int8        ), intent(in   )               :: treeIndex
    integer         (c_size_t              ), intent(in   )               :: iOutput
    logical                                 , intent(in   )               :: nodePassesFilter
    class           (nodeComponentBlackHole)               , pointer      :: blackHole
    integer         (kind=kind_int8        ), allocatable  , dimension(:) :: mergerTreeIndex     , nodeIndex
    double precision                        , allocatable  , dimension(:) :: mass                , massAccretionRate    , radiativeEfficiency, &
         &                                                                   radius              , spin                 , timescale
    double precision                                                      :: accretionRateHotHalo, accretionRateSpheroid
    integer                                                               :: blackHoleCount      , instance
    type            (hdf5Object            )                              :: blackHolesGroup     , outputGroup
    type            (varying_string        )                              :: groupName

    ! If black hole output was requested , output their properties.
    if (nodePassesFilter .and. blackHoleOutputData) then
       ! Get a count of the number of black holes present.
       blackHoleCount=node%blackHoleCount()
       ! Open the output group.
       call hdf5Access%set()
       blackHolesGroup=galacticusOutputFile%openGroup("blackHole","Black hole data.")
       groupName="Output"
       groupName=groupName//iOutput
       outputGroup=blackHolesGroup%openGroup(char(groupName),"Properties of black holes for all trees at each output.")
       call hdf5Access%unset()
       ! Allocate array to store profile.
       call allocateArray(radius             ,[blackHoleCount])
       call allocateArray(spin               ,[blackHoleCount])
       call allocateArray(mass               ,[blackHoleCount])
       call allocateArray(timescale          ,[blackHoleCount])
       call allocateArray(massAccretionRate  ,[blackHoleCount])
       call allocateArray(radiativeEfficiency,[blackHoleCount])
       call allocateArray(nodeIndex          ,[blackHoleCount])
       call allocateArray(mergerTreeIndex    ,[blackHoleCount])
       ! Construct arrays of black hole properties.
       do instance=1,blackHoleCount
          blackHole => node%blackHole(instance=instance)
          call  Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(blackHole,accretionRateSpheroid&
               &,accretionRateHotHalo)
          mass               (instance)=blackHole%mass()
          spin               (instance)=blackHole%spin()
          radius             (instance)=blackHole%radialPosition()
          massAccretionRate  (instance)=accretionRateSpheroid+accretionRateHotHalo
          radiativeEfficiency(instance)=accretionDisks_%efficiencyRadiative(blackHole,massAccretionRate(instance))
          nodeIndex          (instance)=node%index()
          mergerTreeIndex    (instance)=treeIndex
          if (instance > 1) then
             if (blackHoleBinarySeparationGrowthRate_%growthRate(blackHole) /= 0.0d0 )then
                timescale(instance)=-blackHole                           %radialPosition(         ) &
                     &              /blackHoleBinarySeparationGrowthRate_%growthRate    (blackHole)
             else
                timescale(instance)=0.0d0
             end if
          else
             timescale   (instance)=0.0d0
          end if
       end do
       ! Write dataset to the group, first the arrays containing all data.
       call hdf5Access%set()
       call outputGroup%writeDataset(mass,               "mass"               ,"The black hole masses.",                appendTo=.true.)
       call outputGroup%writeDataset(spin,               "spin"               ,"The black hole spins.",                 appendTo=.true.)
       call outputGroup%writeDataset(radius,             "radius"             ,"The black hole radial positions.",      appendTo=.true.)
       call outputGroup%writeDataset(timescale,          "timescale"          ,"The black hole timescales for merger.", appendTo=.true.)
       call outputGroup%writeDataset(radiativeEfficiency,"radiativeEfficiency","The black hole radiative efficiencies.",appendTo=.true.)
       call outputGroup%writeDataset(massAccretionRate,  "accretionRate"      ,"The black hole accretion rates.",       appendTo=.true.)
       call outputGroup%writeDataset(nodeIndex,          "nodeIndex"          ,"The black hole host galaxy inices.",    appendTo=.true.)
       call outputGroup%writeDataset(mergerTreeIndex,    "mergerTreeIndex"    ,"The black hole merger tree indices.",   appendTo=.true.)

       ! Deallocatate profile arrays.
       call deallocateArray(mass               )
       call deallocateArray(spin               )
       call deallocateArray(radius             )
       call deallocateArray(timescale          )
       call deallocateArray(radiativeEfficiency)
       call deallocateArray(massAccretionRate  )
       call deallocateArray(nodeIndex          )
       call deallocateArray(mergerTreeIndex    )
       call outputGroup    %close()
       call blackHolesGroup%close()
       call hdf5Access%unset()
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Properties

  double precision function Hot_Mode_Fraction(node)
    !!{
    A simple interpolating function which is used as a measure of the fraction of a halo which is in the hot accretion mode.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , parameter     :: coolingRadiusFractionalTransitionMinimum=0.9d0
    double precision          , parameter     :: coolingRadiusFractionalTransitionMaximum=1.0d0
    double precision                          :: coolingRadiusFractional                       , x

    coolingRadiusFractional=+coolingRadius_      %      radius(node) &
         &                  /darkMatterHaloScale_%virialRadius(node)
    if      (coolingRadiusFractional < coolingRadiusFractionalTransitionMinimum) then
       Hot_Mode_Fraction=1.0d0
    else if (coolingRadiusFractional > coolingRadiusFractionalTransitionMaximum) then
       Hot_Mode_Fraction=0.0d0
    else
       x=      (coolingRadiusFractional                 -coolingRadiusFractionalTransitionMinimum) &
            & /(coolingRadiusFractionalTransitionMaximum-coolingRadiusFractionalTransitionMinimum)
       Hot_Mode_Fraction=x**2*(2.0d0*x-3.0d0)+1.0d0
    end if
    return
  end function Hot_Mode_Fraction

  double precision function Node_Component_Black_Hole_Standard_Accretion_Rate(self)
    !!{
    Return the rest mass accretion rate onto a standard black hole.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHoleStandard
    implicit none
    class           (nodeComponentBlackHoleStandard), intent(inout) :: self
    double precision                                                :: accretionRateSpheroid, accretionRateHotHalo

    call Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(self,accretionRateSpheroid,accretionRateHotHalo)
    Node_Component_Black_Hole_Standard_Accretion_Rate=accretionRateSpheroid+accretionRateHotHalo
    return
  end function Node_Component_Black_Hole_Standard_Accretion_Rate

  double precision function Node_Component_Black_Hole_Standard_Radiative_Efficiency(self)
    !!{
    Return the radiative efficiency of a standard black hole.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHoleStandard
    implicit none
    class(nodeComponentBlackHoleStandard), intent(inout) :: self

    Node_Component_Black_Hole_Standard_Radiative_Efficiency=accretionDisks_%efficiencyRadiative(self,self%accretionRate())
    return
  end function Node_Component_Black_Hole_Standard_Radiative_Efficiency

  !![
  <postStepTask>
  <unitName>Node_Component_Black_Hole_Standard_Post_Evolve</unitName>
  </postStepTask>
  !!]
  subroutine Node_Component_Black_Hole_Standard_Post_Evolve(node,status)
    !!{
    Keep black hole spin in physical range.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent, nodeComponentBlackHole, treeNode
    use :: Interface_GSL   , only : GSL_Failure
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    integer                                 , intent(inout)          :: status
    class           (nodeComponentBlackHole)               , pointer :: blackHole
    double precision                        , parameter              :: spinMaximum=0.9999d0
    integer                                                          :: i                   , instanceCount
    double precision                                                 :: spin

    ! Check if the standard component is active.
    if (defaultBlackHoleComponent%standardIsActive()) then
       ! Get a count of the number of black holes associated with this node.
       instanceCount=node%blackHoleCount()
       if (instanceCount > 0) then
          ! Iterate over instances.
          do i=1,instanceCount
             ! Get the black hole component.
             blackHole => node%blackHole(instance=i)
             if (blackHole%spin() > spinMaximum .or. blackHole%spin() < 0.0d0) then
                status    =GSL_Failure
                spin      =max(min(blackHole%spin(),spinMaximum),0.0d0)
                call blackHole%spinSet(spin)
             end if
          end do
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Post_Evolve

  double precision function Node_Component_Black_Hole_Standard_Seed_Mass(self)
    !!{
    Return the seed mass for standard black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHoleStandard
    implicit none
    class(nodeComponentBlackHoleStandard), intent(inout) :: self
    !$GLC attributes unused :: self
    
    Node_Component_Black_Hole_Standard_Seed_Mass=blackHoleSeedMass
    return
  end function Node_Component_Black_Hole_Standard_Seed_Mass

end module Node_Component_Black_Hole_Standard
