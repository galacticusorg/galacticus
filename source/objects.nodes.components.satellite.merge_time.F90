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
Contains a module which implements a satellite orbit component that simply tracks the time of merging.
!!}

module Node_Component_Satellite_Merge_Time
  !!{
  Implements a satellite orbit component that simply tracks the time of merging.
  !!}
  implicit none
  private

  !![
  <component>
   <class>satellite</class>
   <name>mergeTime</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>timeOfMerging</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <classDefault>huge(0.0d0)</classDefault>
    </property>
   </properties>
  </component>
  !!]

end module Node_Component_Satellite_Merge_Time
