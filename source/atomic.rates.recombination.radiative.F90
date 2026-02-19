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

!!{
Contains a module which provides a class implementing atomic radiative recombination rates.
!!}

module Atomic_Rates_Recombination_Radiative
  !!{
  Provides a class implementing radiative recombination rates.
  !!}
  private

  !![
  <enumeration>
   <name>recombinationCase</name>
   <description>Enumeration of radiative recombination cases (A or B).</description>
   <visibility>public</visibility>
   <indexing>-2</indexing>
   <entry label="a"/>
   <entry label="b"/>
   <entry label="0"/>
   <entry label="1"/>
  </enumeration>
  !!]

  !![
  <functionClass>
   <name>atomicRecombinationRateRadiative</name>
   <descriptiveName>Atomic Radiative Recombination</descriptiveName>
   <description>Class providing atomic radiative recombination rates.</description>
   <default>verner1996</default>
   <method name="rate" >
    <description>Returns the radiative recombination rate in units of cm$^3$ s$^{-1}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer                                           , intent(in   )           :: atomicNumber, ionizationState</argument>
    <argument>double precision                                  , intent(in   )           :: temperature                  </argument>
    <argument>type            (enumerationRecombinationCaseType), intent(in   ), optional :: level                        </argument>
   </method>
  </functionClass>
  !!]

end module Atomic_Rates_Recombination_Radiative
