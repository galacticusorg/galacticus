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
Contains a module which provides a class that implements active masses for star formation.
!!}

module Star_Formation_Active_Masses
  !!{
  Provides a class that implements calculations of active masses for star formation.
  !!}
  use :: Galacticus_Nodes, only : nodeComponent
  private

  !![
  <functionClass>
   <name>starFormationActiveMass</name>
   <descriptiveName>Active Masses for Star Formation</descriptiveName>
   <description>Class providing models of the actively star-forming gas mass in a galactic component —
    the mass (in $\mathrm{M}_\odot$) of gas that is eligible to form stars, which may differ from the total
    ISM mass depending on the star formation model. For example, only molecular gas or gas above a
    threshold density may be considered active. This active mass is passed to the star formation
    rate law to determine the overall rate of star formation in the disk or spheroid, and is the
    key coupling between the gas reservoir and the stellar mass growth rate.</description>
   <default>totalISM</default>
   <method name="massActive" >
    <description>Returns the mass (in $\mathrm{M}_\odot$) of gas which is actively star forming in the provided \mono{component}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(nodeComponent), intent(inout) :: component</argument>
   </method>
  </functionClass>
  !!]

end module Star_Formation_Active_Masses
