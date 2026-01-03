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
Contains a module which implements a class for calculations of stellar astrophysics.
!!}

module Stellar_Astrophysics
  !!{
  Implements a class for calculations of stellar astrophysics.
  !!}
  private

  !![
  <functionClass>
   <name>stellarAstrophysics</name>
   <descriptiveName>Stellar Astrophysics</descriptiveName>
   <description>
    Class providing models of stellar astrophysics including recycled mass, metal yield, and lifetime as a function of initial
    properties.
   </description>
   <default>file</default>
   <method name="massInitial" >
    <description>Returns the initial mass of a star of given {\normalfont \ttfamily lifetime} and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: lifetime, metallicity</argument>
   </method>
   <method name="massEjected" >
    <description>Returns the mass ejected by a star of given {\normalfont \ttfamily massInitial} and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: massInitial,metallicity</argument>
   </method>
   <method name="massYield" >
    <description>Returns the metal mass yielded by a star of given {\normalfont \ttfamily massInitial} and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: massInitial, metallicity</argument>
    <argument>integer         , intent(in   ), optional :: atomIndex</argument>
   </method>
   <method name="lifetime" >
    <description>Returns the lifetime of a star of given {\normalfont \ttfamily massInitial} and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: massInitial, metallicity</argument>
   </method>
  </functionClass>
  !!]

end module Stellar_Astrophysics
