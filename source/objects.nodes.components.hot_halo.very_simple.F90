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
Contains a module which implements a very simple hot halo node component.
!!}

module Node_Component_Hot_Halo_Very_Simple
  !!{
  Implements a very simple hot halo node component.
  !!}
  implicit none
  private
  public :: Node_Component_Hot_Halo_Very_Simple_Initialize, Node_Component_Hot_Halo_Very_Simple_Scale_Set

  !![
  <component>
   <class>hotHalo</class>
   <name>verySimple</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>mass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
    </property>
    <property>
      <name>abundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the hot halo."/>
    </property>
    <property>
      <name>unaccretedMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of gas that failed to accrete into the hot halo."/>
    </property>
    <property>
      <name>outflowingMass</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>outflowingAbundances</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>abundances</type>
      <rank>0</rank>
    </property>
    <property>
      <name>outerRadius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <output unitsInSI="megaParsec" comment="Outer radius of the hot halo."/>
    </property>
   </properties>
   <bindings>
    <binding method="massBaryonic" function="Node_Component_Hot_Halo_Very_Simple_Mass_Baryonic"/>
   </bindings>
   <functions>objects.nodes.components.hot_halo.very_simple.bound_functions.inc</functions>
  </component>
  !!]

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Initialize(parameters)
    !!{
    Initializes the very simple hot halo component module.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Galacticus_Nodes, only : nodeComponentHotHaloVerySimple
    implicit none
    type(inputParameters               ), intent(inout) :: parameters
    type(nodeComponentHotHaloVerySimple)                :: hotHalo
    !$GLC attributes unused :: parameters
    
    ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
    call hotHalo%      outflowingMassRateFunction(Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate      )
    call hotHalo%outflowingAbundancesRateFunction(Node_Component_Hot_Halo_Very_Simple_Outflowing_Abundances_Rate)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Initialize

  subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
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

    ! Funnel the outflow gas into the hot halo.
    call self%massRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate

  subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Abundances_Rate(self,rate,interrupt,interruptProcedure)
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

    ! Funnel the outflow gas abundances into the hot halo.
    call self%abundancesRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Abundances_Rate

  !![
  <scaleSetTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure, only : unitAbundances
    use :: Galacticus_Nodes    , only : nodeComponentBasic     , nodeComponentHotHalo, nodeComponentHotHaloVerySimple, treeNode, &
         &                              defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    double precision                      , parameter              :: scaleMassRelative=1.0d-2
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                                               :: massVirial

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%verySimpleIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of the very simple class.
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimple)
       ! The the basic component.
       basic => node%basic()
       ! Get virial properties.
       massVirial=basic%mass()
       ! Set the scale.
       call hotHalo%          massScale(               massVirial*scaleMassRelative)
       call hotHalo%unaccretedMassScale(               massVirial*scaleMassRelative)
       call hotHalo%    abundancesScale(unitAbundances*massVirial*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Scale_Set

end module Node_Component_Hot_Halo_Very_Simple
