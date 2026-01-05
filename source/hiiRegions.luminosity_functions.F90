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

!+    Contributions to this file made by: Sachi Weerasooriya

!!{
Contains a module which provides a class that implements HII region luminosity functions.
!!}

module HII_Region_Luminosity_Functions
  !!{
  Provides a class that implements HII region luminosity functions.
  !!}
  private

  !![
  <functionClass>
   <name>hiiRegionLuminosityFunction</name>
   <descriptiveName>HII region luminosity function.</descriptiveName>
   <description>Class providing models of luminosity function for emission line.</description>
   <default>powerLaw</default>
   <method name="cumulativeDistributionFunction">
    <description>Returns the cumulative distribution of the HII region luminosity function between a minimum and maximum $Q_\mathrm{H}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum</argument>
   </method>
   <method name="cumulativeLuminosity">
    <description>Returns the cumulative luminosity from the HII region luminosity function between a minimum and maximum $Q_\mathrm{H}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum</argument>
   </method>
  </functionClass>
  !!]

end module HII_Region_Luminosity_Functions
