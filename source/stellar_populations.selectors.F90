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
   <description>
    Class providing selectors for stellar populations.
   </description>
   <default>fixed</default>
   <method name="select" >
    <description>Return a stellar population.</description>
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
