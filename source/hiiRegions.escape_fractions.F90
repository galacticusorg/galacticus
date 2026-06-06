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
Contains a module which provides a class that implements calculations of escape fractions from HII regions.
!!}
  
module HII_Region_Escape_Fraction
  !!{
  Provides a class that implements calculations of escape fractions from HII regions.
  !!}
  private
  
  !![
  <functionClass>
   <name>hiiRegionEscapeFraction</name>
   <descriptiveName>HII Regions Escape Fractions</descriptiveName>
   <description>
    Class providing models of the escape fraction $f_\mathrm{esc}$ of hydrogen-ionizing (Lyman-continuum)
    photons from individual H\textsc{ii} regions---the fraction of photons produced by massive stars that
    escape the H\textsc{ii} region into the surrounding \gls{ism} or \gls{igm} rather than being absorbed
    by gas within the region. The escape fraction is returned as a function of H\textsc{ii} region age,
    reflecting the time evolution as stellar winds and supernovae clear channels through the surrounding
    gas.
   </description>
   <default>fixed</default>
   <method name="escapeFraction" >
    <description>Returns the escape fraction for hydrogen-ionizing photons.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: ageHIIRegion </argument>
   </method>
  </functionClass>
  !!]

end module HII_Region_Escape_Fraction
