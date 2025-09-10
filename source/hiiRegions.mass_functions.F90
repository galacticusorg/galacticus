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

!+    Contributions to this file made by: Sachi Weerasooriya

!!{
Contains a module which provides a class that implements HII region mass functions.
!!}

module HII_Region_Mass_Functions
  !!{
  Provides a class that implements HII region mass functions.
  !!}
  private

  !![
  <functionClass>
   <name>hiiRegionMassFunction</name>
   <descriptiveName>HII region mass functions</descriptiveName>
   <description>Class providing models of mass functions for HII regions.</description>
   <default>rosolowsky2021</default>
   <method name="cumulativeDistributionFunction">
    <description>Returns the cumulative distribution of the HII region mass function between a minimum and maximum mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: massMinimum, massMaximum</argument>
   </method>
   <method name="cumulativeMass">
    <description>Returns the cumulative mass from the HII region mass function between a minimum and maximum mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: massMinimum, massMaximum</argument>
   </method>
  </functionClass>
  !!]

end module HII_Region_Mass_Functions
