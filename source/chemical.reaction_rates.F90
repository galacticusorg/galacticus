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
Contains a module that implements a class providing calculations of chemical reaction rates.
!!}

module Chemical_Reaction_Rates
  !!{
  Provides a class implementing chemical reaction rates.
  !!}
  use :: Chemical_Abundances_Structure, only : chemicalAbundances
  use :: Galacticus_Nodes             , only : treeNode
  use :: Radiation_Fields             , only : radiationFieldClass
  private

  !![
  <functionClass>
   <name>chemicalReactionRate</name>
   <descriptiveName>Chemical Reaction Rates</descriptiveName>
   <description>
    Class providing chemical reaction rates.
   </description>
   <default>zero</default>
   <method name="rates" >
    <description>Return the collisional excitation cooling rate , in units of J/m$^3$/s, for ion of given {\normalfont \ttfamily atomicNumber} and {\normalfont \ttfamily electronNumber} at temperature {\normalfont \ttfamily T} (in Kelvin).</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision                     , intent(in   ) :: lengthColumn   , temperature</argument>
    <argument>type            (chemicalAbundances ), intent(in   ) :: chemicalDensity             </argument>
    <argument>double precision                     , intent(in   ) :: factorClumping              </argument>
    <argument>class           (radiationFieldClass), intent(inout) :: radiation                   </argument>
    <argument>type            (chemicalAbundances ), intent(inout) :: chemicalRates               </argument>
    <argument>type            (treeNode           ), intent(inout) :: node                        </argument>
   </method>
  </functionClass>
  !!]

end module Chemical_Reaction_Rates
