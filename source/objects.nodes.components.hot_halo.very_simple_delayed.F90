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
Contains a module which implements an extension to the very simple hot halo node component by including an outflowed reservoir
with delayed reincorporation.
!!}

module Node_Component_Hot_Halo_VS_Delayed
  !!{
  Implements an extension to the very simple hot halo node component by including an outflowed reservoir
  with delayed reincorporation.
  !!}
  implicit none
  private
  public :: Node_Component_Hot_Halo_VS_Delayed_Initialize, Node_Component_Hot_Halo_VS_Delayed_Scale_Set

  !![
  <component>
   <class>hotHalo</class>
   <name>verySimpleDelayed</name>
   <extends>
    <class>hotHalo</class>
    <name>verySimple</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>outflowedMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
    </property>
    <property>
      <name>outflowedAbundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the outflowed phase of the hot halo."/>
    </property>
   </properties>
   <bindings>
    <binding method="massBaryonic" function="Node_Component_Hot_Halo_Very_Simple_Delayed_Mass_Baryonic"/>
   </bindings>
   <functions>objects.nodes.components.hot_halo.very_simple_delayed.bound_functions.inc</functions>
  </component>
  !!]

  ! Options controlling the numerical implementation.
  double precision :: scaleRelativeMass

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Hot_Halo_VS_Delayed_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_VS_Delayed_Initialize(parameters)
    !!{
    Initializes the very simple hot halo component module.
    !!}
    use :: Galacticus_Nodes, only : defaultHotHaloComponent, nodeComponentHotHaloVerySimpleDelayed
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type(inputParameters                      ), intent(inout) :: parameters
    type(nodeComponentHotHaloVerySimpleDelayed)                :: hotHalo
    type(inputParameters                      )                :: subParameters

    !$omp critical (Node_Component_Hot_Halo_Very_Simple_Delayed_Initialize)
    if (defaultHotHaloComponent%verySimpleDelayedIsActive()) then
       ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
       call hotHalo%      outflowingMassRateFunction(Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate      )
       call hotHalo%outflowingAbundancesRateFunction(Node_Component_Hot_Halo_VS_Delayed_Outflowing_Abundances_Rate)
       ! Find our parameters.
       subParameters=parameters%subParameters('componentHotHalo')
       ! Read parameters controlling the physical implementation.
       !![
       <inputParameter>
         <name>scaleRelativeMass</name>
         <defaultValue>1.0d-2</defaultValue>
         <description>The mass scale, relative to the total mass of the node, below which calculations in the delayed very simple hot halo component are allowed to become inaccurate.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
    end if
    !$omp end critical (Node_Component_Hot_Halo_Very_Simple_Delayed_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Initialize

  subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas from a galaxy and deposit it into very simple hot halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (                    ), intent(inout), optional, pointer :: interruptProcedure
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Funnel the outflowing gas into the outflowed reservoir.
    call self%outflowedMassRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate

  subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Abundances_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas abundances from a galaxy and deposit them into very simple hot halo.
    !!}
    use :: Abundances_Structure, only : abundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo
    implicit none
    class    (nodeComponentHotHalo), intent(inout)                    :: self
    type     (abundances          ), intent(in   )                    :: rate
    logical                        , intent(inout), optional          :: interrupt
    procedure(                    ), intent(inout), optional, pointer :: interruptProcedure
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Funnel the outflowing gas abundances into the outflowed reservoir.
    call self%outflowedAbundancesRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Abundances_Rate

  !![
  <scaleSetTask>
   <unitName>Node_Component_Hot_Halo_VS_Delayed_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Hot_Halo_VS_Delayed_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure, only : unitAbundances
    use :: Galacticus_Nodes    , only : nodeComponentBasic     , nodeComponentHotHalo, nodeComponentHotHaloVerySimpleDelayed, treeNode, &
         &                              defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                                               :: massVirial

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%verySimpleDelayedIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of our class.
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! The the basic component.
       basic      => node %basic()
       ! Get virial properties.
       massVirial =  basic%mass ()
       ! Set the scale.
       call hotHalo%outflowedMassScale      (               massVirial*scaleRelativeMass)
       call hotHalo%outflowedAbundancesScale(unitAbundances*massVirial*scaleRelativeMass)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Scale_Set

end module Node_Component_Hot_Halo_VS_Delayed
