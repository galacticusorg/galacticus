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
Contains a module which implements the simple black hole node component.
!!}

module Node_Component_Black_Hole_Simple
  !!{
  Implements the simple black hole node component.
  !!}
  use :: Black_Hole_Binary_Mergers     , only : blackHoleBinaryMergerClass
  use :: Cooling_Radii                 , only : coolingRadiusClass
  use :: Dark_Matter_Halo_Scales       , only : darkMatterHaloScaleClass
  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass
  implicit none
  private
  public :: Node_Component_Black_Hole_Simple_Initialize         , Node_Component_Black_Hole_Simple_Scale_Set        , &
       &    Node_Component_Black_Hole_Simple_Thread_Uninitialize, Node_Component_Black_Hole_Simple_Output_Names     , &
       &    Node_Component_Black_Hole_Simple_Output_Count       , Node_Component_Black_Hole_Simple_Output           , &
       &    Node_Component_Black_Hole_Simple_Rate_Compute       , Node_Component_Black_Hole_Simple_Thread_Initialize

  !![
  <component>
   <class>blackHole</class>
   <name>simple</name>
   <isDefault>false</isDefault>
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
      <name>massSeed</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" isDeferred="get" />
    </property>
   </properties>
   <bindings>
      <binding method="enclosedMass" function="Node_Component_Black_Hole_Simple_Enclosed_Mass" bindsTo="component" />
      <binding method="acceleration" function="Node_Component_Black_Hole_Simple_Acceleration"  bindsTo="component" />
      <binding method="tidalTensor"  function="Node_Component_Black_Hole_Simple_Tidal_Tensor"  bindsTo="component" />
   </bindings>
   <functions>objects.nodes.components.black_hole.simple.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass       ), pointer :: darkMatterHaloScale_
  class(coolingRadiusClass             ), pointer :: coolingRadius_
  class(blackHoleBinaryMergerClass     ), pointer :: blackHoleBinaryMerger_
  class(starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_
  !$omp threadprivate(darkMatterHaloScale_,coolingRadius_,blackHoleBinaryMerger_,starFormationRateSpheroids_)

  ! Seed mass for black holes.
  double precision :: blackHoleSeedMass

  ! Feedback parameters.
  double precision :: blackHoleHeatingEfficiency                   , blackHoleWindEfficiency , &
       &              blackHoleToSpheroidStellarGrowthRatio
  logical          :: blackHoleHeatsHotHalo

  ! Output options.
  logical          :: blackHoleOutputAccretion

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Black_Hole_Simple_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Initialize(parameters_)
    !!{
    Initializes the simple black hole node component module.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHoleSimple
    use :: Input_Parameters, only : inputParameter              , inputParameters
    implicit none
    type(inputParameters             ), intent(inout) :: parameters_
    type(nodeComponentBlackHoleSimple)                :: blackHoleSimple

    ! Bind deferred functions.
    call blackHoleSimple%massSeedFunction(Node_Component_Black_Hole_Simple_Seed_Mass)
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
      <name>blackHoleToSpheroidStellarGrowthRatio</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The ratio of the rates of black hole growth and spheroid stellar mass growth.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]
    ! Options controlling AGN feedback.
    !![
    <inputParameter>
      <name>blackHoleHeatsHotHalo</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not the black hole should heat the hot halo.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]
    if (blackHoleHeatsHotHalo) then
       !![
       <inputParameter>
         <name>blackHoleHeatingEfficiency</name>
         <defaultValue>1.0d-3</defaultValue>
         <description>The efficiency with which accretion onto a black hole heats the hot halo.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]
    else
       blackHoleHeatingEfficiency=0.0d0
    end if
    ! Get options controlling winds.
    !![
    <inputParameter>
      <name>blackHoleWindEfficiency</name>
      <defaultValue>2.2157d-3</defaultValue>
      <description>The efficiency of the black hole accretion-driven wind.</description>
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
    return
  end subroutine Node_Component_Black_Hole_Simple_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Black_Hole_Simple_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node random spin module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent     , openMPThreadBindingAtLevel, dependencyRegEx, dependencyDirectionAfter
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    use :: Input_Parameters, only : inputParameter           , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies

    if (defaultBlackHoleComponent%simpleIsActive()) then
       !![
       <objectBuilder class="darkMatterHaloScale"        name="darkMatterHaloScale_"        source="parameters_"/>
       <objectBuilder class="coolingRadius"              name="coolingRadius_"              source="parameters_"/>
       <objectBuilder class="blackHoleBinaryMerger"      name="blackHoleBinaryMerger_"      source="parameters_"/>
       <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters_"/>
       !!]
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent%attach(defaultBlackHoleComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentBlackHoleSimple',dependencies=dependencies)
   end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Black_Hole_Simple_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Thread_Uninitialize()
    !!{
    Uninitializes the tree node random spin module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent
    implicit none

    if (defaultBlackHoleComponent%simpleIsActive()) then
       !![
       <objectDestructor name="darkMatterHaloScale_"       />
       <objectDestructor name="coolingRadius_"             />
       <objectDestructor name="blackHoleBinaryMerger_"     />
       <objectDestructor name="starFormationRateSpheroids_"/>
       !!]
       call satelliteMergerEvent%detach(defaultBlackHoleComponent,satelliteMerger)
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Thread_Uninitialize

  !![
  <scaleSetTask>
   <unitName>Node_Component_Black_Hole_Simple_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole   , nodeComponentBlackHoleSimple, nodeComponentSpheroid, treeNode, &
         &                          defaultBlackHoleComponent
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    class           (nodeComponentBlackHole)               , pointer :: blackHole
    class           (nodeComponentSpheroid )               , pointer :: spheroid
    double precision                        , parameter              :: massScaleAbsolute=1.0d+0

    ! Check if we are the default method.
    if (.not.defaultBlackHoleComponent%simpleIsActive()) return
    ! Get the black hole component.
    blackHole => node%blackHole()
    ! Ensure that it is of the standard class.
    select type (blackHole)
    class is (nodeComponentBlackHoleSimple)
       ! Get the spheroid component.
       spheroid => node%spheroid()
       ! Set scale for mass.
       call blackHole%massScale(                                                                                  &
            &                   max(                                                                              &
            &                                                                     blackHole%massSeed          (), &
            &                       max(                                                                          &
            &                               blackHoleToSpheroidStellarGrowthRatio*spheroid %massStellar       (), &
            &                           max(                                                                      &
            &                                                                                massScaleAbsolute  , &
            &                                                                     blackHole%mass              ()  &
            &                              )                                                                      &
            &                          )                                                                          &
            &                      )                                                                              &
            &                  )
    end select
    return
  end subroutine Node_Component_Black_Hole_Simple_Scale_Set

  !![
  <rateComputeTask>
   <unitName>Node_Component_Black_Hole_Simple_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the black hole mass rate of change.
    !!}
    use :: Galacticus_Nodes            , only : defaultBlackHoleComponent, interruptTask        , nodeComponentBlackHole, nodeComponentBlackHoleSimple, &
          &                                     nodeComponentHotHalo     , nodeComponentSpheroid, propertyTypeInactive  , treeNode
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    type            (treeNode                ), intent(inout)          :: node
    logical                                   , intent(inout)          :: interrupt
    procedure       (interruptTask           ), intent(inout), pointer :: interruptProcedure
    integer                                   , intent(in   )          :: propertyType
    class           (nodeComponentBlackHole  )               , pointer :: blackHole
    class           (nodeComponentSpheroid   )               , pointer :: spheroid
    class           (nodeComponentHotHalo    )               , pointer :: hotHalo
    double precision                          , parameter              :: coolingRadiusFractionalTransitionMinimum=0.9d0
    double precision                          , parameter              :: coolingRadiusFractionalTransitionMaximum=1.0d0
    double precision                                                   :: coolingRadiusFractional                       , couplingEfficiency   , &
         &                                                                energyInputRate                               , heatingRate          , &
         &                                                                massAccretionRate                             , restMassAccretionRate, &
         &                                                                x

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    if (defaultBlackHoleComponent%simpleIsActive()) then

       ! Get the spheroid component.
       spheroid => node%spheroid()

       ! Find the rate of rest mass accretion onto the black hole.
       restMassAccretionRate=blackHoleToSpheroidStellarGrowthRatio*starFormationRateSpheroids_%rate(node)

       ! Finish if there is no accretion.
       if (restMassAccretionRate <= 0.0d0) return

       ! Find the rate of increase in mass of the black hole.
       massAccretionRate=restMassAccretionRate*max((1.0d0-blackHoleHeatingEfficiency-blackHoleWindEfficiency),0.0d0)

       ! Get the black hole component.
       blackHole => node%blackHole()

       ! Detect black hole component type.
       select type (blackHole)
       type is (nodeComponentBlackHole)
          ! Generic type - interrupt and create a simple black hole if accretion rate is non-zero.
          if (massAccretionRate /= 0.0d0) then
             interrupt=.true.
             interruptProcedure => Node_Component_Black_Hole_Simple_Create
          end if
          return
       class is (nodeComponentBlackHoleSimple)
          ! Simple type - continue processing.
          call blackHole%massRate       (     massAccretionRate)
          ! Remove the accreted mass from the spheroid component.
          call spheroid %massGasSinkRate(-restMassAccretionRate)
          ! Add heating to the hot halo component.
          if (blackHoleHeatsHotHalo) then
             ! Compute jet coupling efficiency based on whether halo is cooling quasistatically.
             coolingRadiusFractional=+coolingRadius_      %      radius(node) &
                  &                  /darkMatterHaloScale_%virialRadius(node)
             if      (coolingRadiusFractional < coolingRadiusFractionalTransitionMinimum) then
                couplingEfficiency=1.0d0
             else if (coolingRadiusFractional > coolingRadiusFractionalTransitionMaximum) then
                couplingEfficiency=0.0d0
             else
                x=      (coolingRadiusFractional                 -coolingRadiusFractionalTransitionMinimum) &
                     & /(coolingRadiusFractionalTransitionMaximum-coolingRadiusFractionalTransitionMinimum)
                couplingEfficiency=x**2*(2.0d0*x-3.0d0)+1.0d0
             end if
             ! Compute the heating rate.
             heatingRate=couplingEfficiency*blackHoleHeatingEfficiency*restMassAccretionRate*(speedLight/kilo)**2
             ! Pipe this power to the hot halo.
             hotHalo => node%hotHalo()
             call hotHalo%heatSourceRate(heatingRate,interrupt,interruptProcedure)
          end if
          ! Add energy to the spheroid component.
          if (blackHoleWindEfficiency > 0.0d0) then
             ! Compute the energy input and send it down the spheroid gas energy input pipe.
             energyInputRate=blackHoleWindEfficiency*restMassAccretionRate*(speedLight/kilo)**2
             call spheroid%energyGasInputRate(energyInputRate)
          end if
       end select
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Rate_Compute

  subroutine satelliteMerger(self,node)
    !!{
    Merge (instantaneously) any simple black hole associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode
    implicit none
    class           (*                     ), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    type            (treeNode              ), pointer       :: nodeHost
    class           (nodeComponentBlackHole), pointer       :: blackHoleHost   , blackHole
    double precision                                        :: massBlackHoleNew, spinBlackHoleNew
    !$GLC attributes unused :: self
    
    ! Find the node to merge with.
    nodeHost      => node    %mergesWith(                 )
    ! Get the black holes.
    blackHole     => node    %blackHole (autoCreate=.true.)
    blackHoleHost => nodeHost%blackHole (autoCreate=.true.)
    ! Compute the effects of the merger.
    call blackHoleBinaryMerger_%merge(blackHole    %mass(), &
         &                            blackHoleHost%mass(), &
         &                            0.0d0               , &
         &                            0.0d0               , &
         &                            massBlackHoleNew    , &
         &                            spinBlackHoleNew      &
         &                           )
    ! Move the black hole to the host.
    call blackHoleHost%massSet(massBlackHoleNew)
    call blackHole    %massSet(           0.0d0)
    return
  end subroutine satelliteMerger

  subroutine Node_Component_Black_Hole_Simple_Create(node)
    !!{
    Creates a simple black hole component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode
    implicit none
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentBlackHole)               , pointer :: blackHole

    ! Create the component.
    blackHole => node%blackHole(autoCreate=.true.)
    ! Set the seed mass.
    call blackHole%massSet(blackHoleSeedMass)
    return
  end subroutine Node_Component_Black_Hole_Simple_Create

  !![
  <mergerTreeOutputNames>
   <unitName>Node_Component_Black_Hole_Simple_Output_Names</unitName>
   <sortName>Node_Component_Black_Hole_Simple_Output</sortName>
  </mergerTreeOutputNames>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Output_Names(node,integerProperty,integerProperties,doubleProperty,doubleProperties,time)
    !!{
    Set names of black hole properties to be written to the \glc\ output file.
    !!}
    use :: Galacticus_Nodes                  , only : treeNode
    use :: Numerical_Constants_Astronomical  , only : gigaYear             , massSolar
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
    implicit none
    type            (treeNode)                           , intent(inout) :: node
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleProperty   , integerProperty
    type            (outputPropertyInteger), dimension(:), intent(inout) :: integerProperties
    type            (outputPropertyDouble ), dimension(:), intent(inout) :: doubleProperties
    !$GLC attributes unused :: time, integerProperty, integerProperties

    ! Ensure that the black hole component is of the simple class.
    if (Node_Component_Black_Hole_Simple_Matches(node)) then
       if (blackHoleOutputAccretion) then
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%name     ='blackHoleAccretionRate'
          doubleProperties(doubleProperty)%comment  ='Rest-mass accretion rate onto the black hole.'
          doubleProperties(doubleProperty)%unitsInSI=massSolar/gigaYear
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Output_Names

  !![
  <mergerTreeOutputPropertyCount>
   <unitName>Node_Component_Black_Hole_Simple_Output_Count</unitName>
   <sortName>Node_Component_Black_Hole_Simple_Output</sortName>
  </mergerTreeOutputPropertyCount>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Output_Count(node,integerPropertyCount,doublePropertyCount,time)
    !!{
    Account for the number of black hole properties to be written to the the \glc\ output file.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount  , integerPropertyCount
    integer                   , parameter     :: extraPropertyCount =1
    !$GLC attributes unused :: time, integerPropertyCount

    ! Ensure that the black hole component is of the simple class.
    if (Node_Component_Black_Hole_Simple_Matches(node)) then
       if (blackHoleOutputAccretion) doublePropertyCount=doublePropertyCount+extraPropertyCount
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Output_Count

  !![
  <mergerTreeOutputTask>
   <unitName>Node_Component_Black_Hole_Simple_Output</unitName>
   <sortName>Node_Component_Black_Hole_Simple_Output</sortName>
  </mergerTreeOutputTask>
  !!]
  subroutine Node_Component_Black_Hole_Simple_Output(node,integerProperty,integerBufferCount,integerProperties,doubleProperty,doubleBufferCount,doubleProperties,time,instance)
    !!{
    Store black hole properties in the \glc\ output file buffers.
    !!}
    use :: Galacticus_Nodes                  , only : nodeComponentBlackHole, nodeComponentSpheroid, treeNode
    use :: Kind_Numbers                      , only : kind_int8
    use :: Multi_Counters                    , only : multiCounter
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger , outputPropertyDouble
    implicit none
    double precision                        , intent(in   )               :: time
    type            (treeNode              ), intent(inout)               :: node
    integer                                 , intent(inout)               :: doubleBufferCount    , doubleProperty , &
         &                                                                   integerBufferCount   , integerProperty
    type            (outputPropertyInteger ), intent(inout), dimension(:) :: integerProperties
    type            (outputPropertyDouble  ), intent(inout), dimension(:) :: doubleProperties
    type            (multiCounter          ), intent(inout)               :: instance
    class           (nodeComponentBlackHole)               , pointer      :: blackHole
    double precision                                                      :: restMassAccretionRate
    !$GLC attributes unused :: time, integerProperty, integerBufferCount, integerProperties, instance

    ! Ensure that the black hole component is of the simple class.
    if (Node_Component_Black_Hole_Simple_Matches(node)) then
       ! Get the black hole component.
       blackHole => node%blackHole()
       ! Store the properties.
       if (blackHoleOutputAccretion) then
          ! Get the rest mass accretion rate.
          restMassAccretionRate=blackHoleToSpheroidStellarGrowthRatio*starFormationRateSpheroids_%rate(node)
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%scalar(doubleBufferCount)=restMassAccretionRate
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Output

  logical function Node_Component_Black_Hole_Simple_Matches(node)
    !!{
    Return true if the black hole component of {\normalfont \ttfamily node} is a match to the simple implementation.
    !!}
    use :: Galacticus_Nodes, only : defaultBlackHoleComponent, nodeComponentBlackHole, nodeComponentBlackHoleSimple, treeNode
    implicit none
    type (treeNode              ), intent(inout) :: node
    class(nodeComponentBlackHole), pointer       :: blackHole

    ! Get the black hole component.
    blackHole => node%blackHole()
    ! Ensure that it is of the simple class.
    Node_Component_Black_Hole_Simple_Matches=.false.
    select type (blackHole)
    class is (nodeComponentBlackHoleSimple)
       Node_Component_Black_Hole_Simple_Matches=.true.
    type  is (nodeComponentBlackHole       )
       Node_Component_Black_Hole_Simple_Matches=defaultBlackHoleComponent%simpleIsActive()
    end select
    return
  end function Node_Component_Black_Hole_Simple_Matches

  double precision function Node_Component_Black_Hole_Simple_Seed_Mass(self)
    !!{
    Return the seed mass for simple black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHoleSimple
    implicit none
    class(nodeComponentBlackHoleSimple), intent(inout) :: self
    !$GLC attributes unused :: self
    
    Node_Component_Black_Hole_Simple_Seed_Mass=blackHoleSeedMass
    return
  end function Node_Component_Black_Hole_Simple_Seed_Mass

end module Node_Component_Black_Hole_Simple
