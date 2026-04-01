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

!+ Contributions to this file made by: Andrew Benson, Daniel McAndrew.

!!{
Contains a module that implements a dielectronic recombination rate class.
!!}

module Atomic_Rates_Recombination_Dielectronic
  !!{
  Implements a dielectronic recombination rates class.
  !!}
  private

  !![
  <functionClass>
   <name>atomicRecombinationRateDielectronic</name>
   <descriptiveName>Atomic Dielectronic Recombination Rates</descriptiveName>
   <description>Class providing dielectronic recombination rate coefficients (in cm$^3$ s$^{-1}$) for ions of given
    atomic number and electron number at a specified electron temperature. Dielectronic recombination proceeds via
    simultaneous capture of a free electron and excitation of a bound electron, followed by radiative stabilization;
    it dominates over radiative recombination at temperatures above $\sim 10^5$ K for many ions and is therefore a
    critical process in the ionization balance and cooling of hot collisionally ionized gas in the \gls{igm} and
    \gls{cgm}.</description>
   <default>arnaud1985</default>
   <method name="rate" >
    <description>Return the dielectroninc recombination rate (in units of cm$^3$ s$^{-1}$) for the ion of given \mono{atomicNumber} and \mono{electronNumber} at the given \mono{temperature} (in Kelvin).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer         , intent(in   ) :: atomicNumber, electronNumber</argument>
    <argument>double precision, intent(in   ) :: temperature</argument>
   </method>
  </functionClass>
  !!]

end module Atomic_Rates_Recombination_Dielectronic
