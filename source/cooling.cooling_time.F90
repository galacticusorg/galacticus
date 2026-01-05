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
Contains a module that implements calculations of the cooling time.
!!}

module Cooling_Times
  !!{
  Implements calculations of the cooling time.
  !!}
  use :: Abundances_Structure         , only : abundances
  use :: Chemical_Abundances_Structure, only : chemicalAbundances
  use :: Galacticus_Nodes             , only : treeNode
  use :: Radiation_Fields             , only : radiationFieldClass
  implicit none
  private

  !![
  <functionClass>
   <name>coolingTime</name>
   <descriptiveName>Cooling times.</descriptiveName>
   <description>
    Class providing models of the cooling time for gas in the hot atmosphere surrounding a galaxy.
   </description>
   <default>simple</default>
   <method name="time" >
    <description>Returns the cooling time for gas in the hot atmosphere surrounding the galaxy in units of Gyr.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode           ), intent(inout) :: node                      </argument>
    <argument>double precision                     , intent(in   ) :: temperature      , density</argument>
    <argument>type            (abundances         ), intent(in   ) :: gasAbundances             </argument>
    <argument>type            (chemicalAbundances ), intent(in   ) :: chemicalDensities         </argument>
    <argument>class           (radiationFieldClass), intent(inout) :: radiation                 </argument>
   </method>
   <method name="gradientDensityLogarithmic" >
    <description>Returns the logarithmic derivative of cooling time with respect to density for gas in the hot atmosphere surrounding the galaxy.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode           ), intent(inout) :: node                      </argument>
    <argument>double precision                     , intent(in   ) :: temperature      , density</argument>
    <argument>type            (abundances         ), intent(in   ) :: gasAbundances             </argument>
    <argument>type            (chemicalAbundances ), intent(in   ) :: chemicalDensities         </argument>
    <argument>class           (radiationFieldClass), intent(inout) :: radiation                 </argument>
   </method>
   <method name="gradientTemperatureLogarithmic" >
    <description>Returns the logarithmic derivative of cooling time with respect to temperature for gas in the hot atmosphere surrounding the galaxy.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode           ), intent(inout) :: node                      </argument>
    <argument>double precision                     , intent(in   ) :: temperature      , density</argument>
    <argument>type            (abundances         ), intent(in   ) :: gasAbundances             </argument>
    <argument>type            (chemicalAbundances ), intent(in   ) :: chemicalDensities         </argument>
    <argument>class           (radiationFieldClass), intent(inout) :: radiation                 </argument>
   </method>
  </functionClass>
  !!]

end module Cooling_Times
