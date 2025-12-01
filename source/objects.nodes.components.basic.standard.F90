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
Contains a module with the standard implementation of basic tree node methods.
!!}

module Node_Component_Basic_Standard
  !!{
  The standard implementation of basic tree node methods.
  !!}
  implicit none
  private
  public :: Node_Component_Basic_Standard_Scale_Set, Node_Component_Basic_Standard_Plausibility

  !![
  <component>
   <class>basic</class>
   <name>standard</name>
   <isDefault>true</isDefault>
   <properties>
    <property>
      <name>mass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Total mass of the node, assuming universal baryon fraction."/>
    </property>
    <property>
      <name>time</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>timeLastIsolated</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <getFunction>BasicStandardTimeLastIsolated</getFunction>
      <output unitsInSI="gigaYear" comment="Time at which node was last an isolated halo."/>
    </property>
    <property>
      <name>accretionRate</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
   <functions>objects.nodes.components.basic.standard.bound_functions.inc</functions>
  </component>
  !!]

contains

  !![
  <scaleSetTask>
   <unitName>Node_Component_Basic_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Basic_Standard_Scale_Set(node)
    !!{
    Set scales for properties in the standard implementation of the basic component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentBasicStandard, treeNode
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    double precision                    , parameter              :: timeScale        =1.0d-3
    double precision                    , parameter              :: scaleMassRelative=1.0d-6
    class           (nodeComponentBasic)               , pointer :: basic

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard class.
    select type (basic)
    class is (nodeComponentBasicStandard)
       ! Set scale for time.
       call basic%timeScale(timeScale                     )
       ! Set scale for mass.
       call basic%massScale(basic%mass()*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Basic_Standard_Scale_Set

  !![
  <radiusSolverPlausibility>
   <unitName>Node_Component_Basic_Standard_Plausibility</unitName>
  </radiusSolverPlausibility>
  !!]
  subroutine Node_Component_Basic_Standard_Plausibility(node)
    !!{
    Determines whether the basic is physically plausible. Require the mass and time to be positive.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentBasicStandard, treeNode
    implicit none
    type   (treeNode          ), intent(inout) :: node
    class  (nodeComponentBasic), pointer       :: basic

    basic => node%basic()
    select type (basic)
    class is (nodeComponentBasicStandard)
       if (basic%mass() <= 0.0d0 .or. basic%time() <= 0.0d0) then
          node%isPhysicallyPlausible=.false.
          node%isSolvable           =.false.
       end if
    end select
    return
  end subroutine Node_Component_Basic_Standard_Plausibility

end module Node_Component_Basic_Standard
