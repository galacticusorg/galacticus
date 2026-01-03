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
Provides a class that implements the distribution of HII region densities.
!!}

module HII_Region_Density_Distributions
  !!{
  Provides a class that implements calculations for hydrogen density distribution in a HII region.
  !!}

  private

  !![
  <functionClass>
   <name>hiiRegionDensityDistribution</name>
   <descriptiveName>HII region density distribution</descriptiveName>
   <description>
    Class providing models for the distribution of HII region hydrogen density.
   </description>
   <default>logNormal</default>
   <method name="cumulativeDensityDistribution">
    <description>Return the cumulative distribution of HII region hydrogen density between a minimum and maximum $n_\mathrm{H}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: densityHydrogenMinimum, densityHydrogenMaximum</argument>
   </method>
  </functionClass>
  !!]

end module HII_Region_Density_Distributions

