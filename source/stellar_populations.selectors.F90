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
Contains a module which implements a class for selecting stellar populations.
!!}

module Stellar_Population_Selectors
  !!{
  Implements a class for selecting stellar populations.
  !!}
  use :: Abundances_Structure, only : abundances
  use :: Galacticus_Nodes    , only : nodeComponent
  use :: Stellar_Populations , only : stellarPopulationClass
  implicit none
  private

  !![
  <functionClass>
   <name>stellarPopulationSelector</name>
   <descriptiveName>Stellar Population Selectors</descriptiveName>
   <description>Class providing selectors for stellar populations---algorithms that choose the appropriate
    \refClass{stellarPopulationClass} object to associate with a star formation event given the current
    star formation rate, the component metallicity, and the galaxy component. Different implementations
    may always return the same (fixed) stellar population, or may select a metal-poor population at low
    metallicities to represent Population~III stars or a starburst population at high star formation
    rates, allowing heterogeneous stellar populations within a single galaxy.</description>
   <default>fixed</default>
   <method name="select" >
    <description>Return the appropriate \refClass{stellarPopulationClass} object for a star formation event given the current star formation rate, the elemental abundances of the fuel, and the galaxy component type (disk, spheroid, etc.).</description>
    <type>class(stellarPopulationClass)</type>
    <pass>yes</pass>
    <argument>double precision               , intent(in   ) :: rateStarFormation</argument>
    <argument>type            (abundances   ), intent(in   ) :: abundances_      </argument>
    <argument>class           (nodeComponent), intent(in   ) :: component        </argument>
   </method>
   <method name="isStarFormationRateDependent" >
    <description>Return true if the selection of stellar population is dependent on star formation rate.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Stellar_Population_Selectors
