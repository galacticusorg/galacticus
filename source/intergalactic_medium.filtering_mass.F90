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
Contains a module which provides a class for filtering masses.
!!}

module Intergalactic_Medium_Filtering_Masses
  !!{
  Provides a class for filtering masses.
  !!}
  private

  !![
  <functionClass>
   <name>intergalacticMediumFilteringMass</name>
   <descriptiveName>Intergalactic Medium Filtering Mass</descriptiveName>
   <description>Class providing intergalactic medium filtering mass.</description>
   <default>gnedin2000</default>
   <method name="massFiltering" >
    <description>Return the filtering mass at the given {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="massFilteringRateOfChange" >
    <description>Return the rate of change of the filtering mass at the given {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time</argument>
   </method>
   <method name="fractionBaryons" >
    <description>Return the fraction of baryons accreted into a halo of the given {\normalfont \ttfamily mass} at the {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass, time</argument>
   </method>
   <method name="fractionBaryonsGradientMass" >
    <description>Return the gradient with respect to mass of the fraction of baryons accreted into a halo of the given {\normalfont \ttfamily mass} at the {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass, time</argument>
   </method>
   <method name="fractionBaryonsRateOfChange" >
    <description>Return the rate of change of the fraction of baryons accreted into a halo of the given {\normalfont \ttfamily mass} at the {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass, time</argument>
   </method>
  </functionClass>
  !!]

end module Intergalactic_Medium_Filtering_Masses
