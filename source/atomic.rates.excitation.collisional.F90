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

!+ Contributions to this file made by: Daniel McAndrew.

!!{
Contains a module which provides a class implementing atomic collisional excitation rates.
!!}

module Atomic_Rates_Excitation_Collisional
  !!{
  Provides a class implementing atomic collisional excitation rates.
  !!}
  private

  !![
  <functionClass>
   <name>atomicExcitationRateCollisional</name>
   <descriptiveName>Atomic Collisional Excitation</descriptiveName>
   <description>Class providing atomic collisional excitation rates.</description>
   <default>scholzWalters1991</default>
   <method name="coolingRate" >
    <description>Return the collisional excitation cooling rate , in units of J m$^3$ s$^{-1}$, for ion of given {\normalfont \ttfamily atomicNumber} and {\normalfont \ttfamily electronNumber} at temperature {\normalfont \ttfamily T} (in Kelvin).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer         , intent(in   ) :: atomicNumber, electronNumber</argument>
    <argument>double precision, intent(in   ) :: temperature                 </argument>
   </method>
  </functionClass>
  !!]

end module Atomic_Rates_Excitation_Collisional
