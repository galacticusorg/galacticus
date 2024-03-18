!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which provides a class that implements the distribution of luminosity function.
!!}

module hii_Region_Density_Distributions
  !!{
  Provides a class that implements calculations for hydrogen density distribution in a HII region.
  !!}

  private

  !![
  <functionClass>
   <name>hiiRegionDensityDistribution</name>
   <descriptiveName>HII region density distribution.</descriptiveName>
   <description>
    Class providing models for density distribution of hydrogen.
   </description>
   <default>densityDistribution</default>
   <method name="densityFunction" >
    <description>Return the density distribution of HII regions in galaxies. </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: densityHydrogen </argument>
   </method>
  </functionClass>
  !!]

end module hii_Region_Density_Distributions

