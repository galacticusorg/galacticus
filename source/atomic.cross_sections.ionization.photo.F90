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
Contains a module that implements the atomic photo-ionization cross-section class.
!!}

module Atomic_Cross_Sections_Ionization_Photo
  !!{
  Implements the atomic photo-ionization cross-section class.
  !!}
  private

  !![
  <functionClass>
   <name>atomicCrossSectionIonizationPhoto</name>
   <descriptiveName>Atomic cross sections for photo-ionization.</descriptiveName>
   <description>Class providing atomic cross sections for photo-ionization.</description>
   <default>verner</default>
   <method name="crossSection" >
    <description>Returns the cross-section for photoionization (in units of cm$^2$) for a given atom in a given ionization state at the specified {\normalfont \ttfamily wavelength} (given in units of \AA).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer         , intent(in   ) :: atomicNumber, ionizationState, shellNumber</argument>
    <argument>double precision, intent(in   ) :: wavelength</argument>
   </method>
  </functionClass>
  !!]

end module Atomic_Cross_Sections_Ionization_Photo
