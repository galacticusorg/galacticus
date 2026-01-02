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
Contains a module which implements a class for calculations of Type Ia supernovae.
!!}

module Supernovae_Type_Ia
  !!{
  Implements a class for calculations of Type Ia supernovae.
  !!}
  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionClass
  implicit none
  private

  !![
  <functionClass>
   <name>supernovaeTypeIa</name>
   <descriptiveName>Supernovae Type Ia</descriptiveName>
   <description>
    Class providing models of supernovae type Ia, including the cumulative number occurring and metal yield.
   </description>
   <default>nagashima2005</default>
   <method name="massInitialRange" >
    <description>Return the range of initial stellar masses that contribute to the Type Ia population.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class           (initialMassFunctionClass), intent(inout) :: initialMassFunction_                    </argument>
    <argument>double precision                          , intent(in   ) :: age                 , metallicity       </argument>
    <argument>double precision                          , intent(  out) :: massInitialMinimum  , massInitialMaximum</argument>
   </method>
   <method name="number" >
    <description>Return the cumulative number of Type Ia supernovae from a stellar population of the given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily age}, and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class           (initialMassFunctionClass), intent(inout), target :: initialMassFunction_                  </argument>
    <argument>double precision                          , intent(in   )         :: initialMass         , age, metallicity</argument>
   </method>
   <method name="yield" >
    <description>Return the cumulative yield from Type Ia supernoave from a stellar population of the given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily age}, and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class           (initialMassFunctionClass), intent(inout)           :: initialMassFunction_                  </argument>
    <argument>double precision                          , intent(in   )           :: initialMass         , age, metallicity</argument>
    <argument>integer                                   , intent(in   ), optional :: atomIndex                             </argument>
   </method>
  </functionClass>
  !!]

end module Supernovae_Type_Ia
