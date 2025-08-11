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
Contains a module which implements the non-central black hole node component.
!!}

module Node_Component_Black_Hole_Noncentral
  !!{
  Implement non-central black hole tree node methods.
  !!}
  implicit none
  private

  !![
  <component>
   <class>blackHole</class>
   <name>nonCentral</name>
   <extends>
    <class>blackHole</class>
    <name>standard</name>
   </extends>
   <isDefault>false</isDefault>
   <output instances="first"/>
   <properties>
    <property>
      <name>radialPosition</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
   </properties>
  </component>
  !!]

end module Node_Component_Black_Hole_Noncentral
