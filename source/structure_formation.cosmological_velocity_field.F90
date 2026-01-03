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
Contains a module which provides a class that implements models of the cosmological velocity field.
!!}

module Cosmological_Velocity_Field
  !!{
  Provides a class that implements models of the cosmological velocity field.
  !!}
  private

  !![
  <functionClass>
   <name>cosmologicalVelocityField</name>
   <descriptiveName>Cosmological Velocity Field</descriptiveName>
   <description>
    A class providing models of the cosmological velocity field.
   </description>
   <default>filteredPower</default>
   <method name="velocityRadialMeanPairwise">
    <description>Return the mean radial velocity (averaged over all positions; in km/s) at a given {\normalfont \ttfamily separation} (in units of Mpc) and {\normalfont \ttfamily time} (in units of Gyr). If {\normalfont \ttfamily includeHubbleFlow} is {\normalfont \ttfamily true} then the Hubble flow is included, otherwise only the peculiar component of the mean radial velocity is computed.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: separation, time</argument>
    <argument>logical         , intent(in   ) :: includeHubbleFlow</argument>
   </method>
   <method name="velocityDispersion1D">
    <description>Return the 1-D dispersion of the velocity field (in units of km/s) on a scale corresponding to the given {\normalfont \ttfamily mass} (in units of $\mathrm{M}_\odot$) and the given {\normalfont \ttfamily time} (in units of Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass, time</argument>
   </method>
  <method name="velocityDispersion1DHaloPairwise">
    <description>Return the 1-D dispersion of the velocity field (in units of km/s) for pairs of halos of the given {\normalfont \ttfamily mass1} and {\normalfont \ttfamily mass2} (in units of $\mathrm{M}_\odot$) at the given {\normalfont \ttfamily separation} (in units of Mpc) and the given {\normalfont \ttfamily time} (in units of Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass1, mass2, separation, time</argument>
   </method>
  </functionClass>
  !!]

end module Cosmological_Velocity_Field
