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
Contains a module which implements a class for calculations of stellar winds.
!!}

module Stellar_Astrophysics_Winds
  !!{
  Implements a class for calculations of stellar winds.
  !!}
  implicit none
  private

  !![
  <functionClass>
   <name>stellarWinds</name>
   <descriptiveName>Stellar Winds</descriptiveName>
   <description>
    Class providing models of stellar winds.
   </description>
   <default>leitherer1992</default>
   <method name="rateMassLoss" >
    <description>Return the mass loss rate (in $M_\odot$/Gyr) from stars of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily age} and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: initialMass, age, metallicity</argument>
   </method>
   <method name="velocityTerminal" >
    <description>Return the terminal velocity (in km/s) of winds from stars of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily age} and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: initialMass, age, metallicity</argument>
   </method>
  </functionClass>
  !!]

end module Stellar_Astrophysics_Winds
