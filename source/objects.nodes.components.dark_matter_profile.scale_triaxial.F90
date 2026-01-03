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
Contains a module which implements a dark matter profile method that provides a scale radius and a shape parameter.
!!}

module Node_Component_Dark_Matter_Profile_Scale_Triaxial
  !!{
  Implements a dark matter profile method that provides a scale radius and a shape parameter.
  !!}
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Triaxial_Scale_Set

  !![
  <component>
   <class>darkMatterProfile</class>
   <name>scaleTriaxial</name>
   <isDefault>false</isDefault>
   <extends>
    <class>darkMatterProfile</class>
    <name>scale</name>
   </extends>
   <properties>
    <property>
      <name>axisRatios</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true"/>
      <output labels="[X,Y,Z]" unitsInSI="0.0d0" comment="Triaxial axis ratios the dark matter profile."/>
      <classDefault>[-1.0d0,-1.0d0,-1.0d0]</classDefault>
    </property>
   </properties>
  </component>
  !!]

contains

  !![
  <scaleSetTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Triaxial_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Triaxial_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScaleTriaxial, treeNode
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile

    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure it is of the scale+triaxial class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileScaleTriaxial)
       ! Set scale for the scale radius.
       call darkMatterProfile%axisRatiosScale([1.0d0,1.0d0,1.0d0])
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Triaxial_Scale_Set

end module Node_Component_Dark_Matter_Profile_Scale_Triaxial
