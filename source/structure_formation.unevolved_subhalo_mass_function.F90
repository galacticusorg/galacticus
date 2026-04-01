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
   <description>Class providing unevolved subhalo mass functions $\mathrm{d}N/\mathrm{d}m$---the mass spectrum of
    subhalos at the time they first fell into their host, before tidal stripping has reduced their mass.
    The unevolved subhalo mass function represents the initial conditions for subhalo evolution and is predicted from
    the extended Press-Schechter formalism or calibrated to N-body simulations. Both differential
    (per unit mass per host halo) and integrated (between two mass limits) forms are provided, as
    functions of the subhalo mass, host mass, and cosmic time.</description>
   <default>giocoli2008</default>
   <method name="differential" >
    <description>Return the differential unevolved subhalo mass function per halo for \mono{mass} [$\mathrm{M}_\odot$] subhalos in \mono{massHost} [$\mathrm{M}_\odot$] hosts at \mono{time} [Gyr].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time, mass, massHost</argument>
   </method>
   <method name="integrated" >
    <description>Return the unevolved subhalo mass function per host at \mono{time} [Gyr] in hosts of mass \mono{massHost} [$\mathrm{M}_\odot$] integrated between \mono{massLow} and \mono{massHigh} [$\mathrm{M}_\odot$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: time, massLow, massHigh, massHost</argument>
   </method>
  </functionClass>
  !!]

end module Unevolved_Subhalo_Mass_Functions
