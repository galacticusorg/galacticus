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
Contains a module which implements a class for calculation of stellar tracks.
!!}

module Stellar_Astrophysics_Tracks
  !!{
  Implements a class for stellar tracks.
  !!}
  implicit none
  private

  !![
  <functionClass>
   <name>stellarTracks</name>
   <descriptiveName>Stellar Tracks</descriptiveName>
   <description>
    Class providing models of stellar tracks in luminosity and effective temperature as a function of initial mass,
    metallicity, and age.
   </description>
   <default>file</default>
   <method name="luminosity" >
    <description>Returns the bolometric luminosity of a star of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily metallicity} and {\normalfont \ttfamily age}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: initialMass, metallicity, age</argument>
   </method>
   <method name="temperatureEffective" >
    <description>Returns the effective temperature of a star of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily metallicity} and {\normalfont \ttfamily age}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: initialMass, metallicity, age</argument>
   </method>
  </functionClass>
  !!]

end module Stellar_Astrophysics_Tracks
