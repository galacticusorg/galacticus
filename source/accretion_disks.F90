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
Contains a module which provides a class that implements accretion disks.
!!}

module Accretion_Disks
  !!{
  Provides a class that implements accretion disks.
  !!}
  use :: Galacticus_Nodes, only : nodeComponentBlackHole
  private

  !![
  <functionClass>
   <name>accretionDisks</name>
   <descriptiveName>Accretion disks</descriptiveName>
   <description>
    A class implementing circumnuclear accretion disks. Circumnuclear accretion disks surrounding supermassive black holes at
    the centers of galaxies influence the evolution of both the black hole (via accretion rates of mass and angular momentum
    and possibly by extracting rotational energy from the black hole) and the surrounding galaxy if they lead to energetic
    outflows (e.g. jets) from the nuclear region.
   </description>
   <default>shakuraSunyaev</default>
   <method name="efficiencyRadiative" >
    <description>Returns the radiative efficiency of the accretion disk.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class           (nodeComponentBlackHole), intent(inout) :: blackHole</argument>
    <argument>double precision                        , intent(in   ) :: accretionRateMass</argument>
   </method>
   <method name="powerJet" >
    <description>Returns the power of the jet launched by the accretion disk in units of $M_\odot$ (km/s)$^2$ Gyr$^{-1}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class           (nodeComponentBlackHole), intent(inout) :: blackHole</argument>
    <argument>double precision                        , intent(in   ) :: accretionRateMass</argument>
   </method>
   <method name="rateSpinUp" >
    <description>Returns the spin-up rate of the black hole due to accretion from the accretion disk.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class           (nodeComponentBlackHole), intent(inout) :: blackHole</argument>
    <argument>double precision                        , intent(in   ) :: accretionRateMass</argument>
   </method>
  </functionClass>
  !!]

end module Accretion_Disks
