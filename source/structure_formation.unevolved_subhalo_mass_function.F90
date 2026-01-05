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
Contains a module which provides a class that implements subhalo mass functions.
!!}

module Unevolved_Subhalo_Mass_Functions
  implicit none
  private

  !![
  <functionClass>
   <name>unevolvedSubhaloMassFunction</name>
   <descriptiveName>Unevolved Subhalo Mass Function</descriptiveName>
   <description>Class providing unevolved subhalo mass functions.</description>
   <default>giocoli2008</default>
   <method name="differential" >
    <description>Return the differential unevolved subhalo mass function per halo for {\normalfont \ttfamily mass} [$M_\odot$] subhalos in {\normalfont \ttfamily massHost} [$M_\odot$] hosts at {\normalfont \ttfamily time} [Gyr].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time, mass, massHost</argument>
   </method>
   <method name="integrated" >
    <description>Return the unevolved subhalo mass function per host at {\normalfont \ttfamily time} [Gyr] in hosts of mass {\normalfont \ttfamily massHost} [$M_\odot$] integrated between {\normalfont \ttfamily massLow} and {\normalfont \ttfamily massHigh} [$M_\odot$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time, massLow, massHigh, massHost</argument>
   </method>
  </functionClass>
  !!]

end module Unevolved_Subhalo_Mass_Functions
