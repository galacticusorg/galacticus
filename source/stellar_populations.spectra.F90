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

!+    Contributions to this file made by:  Alex Merson.

!!{
Contains a module that provides a class implementing stellar population spectra.
!!}

module Stellar_Population_Spectra
  !!{
  Provides a class implementing stellar population spectra.
  !!}
  use :: Abundances_Structure, only : abundances
  private

  !![
  <functionClass>
   <name>stellarPopulationSpectra</name>
   <descriptiveName>Stellar Population Spectra</descriptiveName>
   <description>
    Class providing stellar population spectra, which are used to construct integrated spectra of galaxies.
   </description>
   <default>FSPS</default>
   <method name="luminosity" >
    <description>Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population, composition {\normalfont \ttfamily abundances}, of the given {\normalfont \ttfamily age} (in Gyr), at the specified {\normalfont \ttfamily wavelength} (in Angstroms).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (abundances), intent(in   )           :: abundancesStellar</argument>
    <argument>double precision            , intent(in   )           :: age              , wavelength</argument>
    <argument>integer                     , intent(  out), optional :: status</argument>
   </method>
   <method name="tabulation" >
    <description>Return a tabulation of ages and metallicities at which stellar spectra should be tabulated.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer                                    , intent(  out) :: agesCount, metallicitiesCount</argument>
    <argument>double precision, allocatable, dimension(:), intent(  out) :: ages     , metallicity</argument>
   </method>
   <method name="wavelengths" >
    <description>Return a tabulation of wavelengths at which stellar spectra are defined.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer                                    , intent(  out) :: wavelengthsCount</argument>
    <argument>double precision, allocatable, dimension(:), intent(  out) :: wavelengths</argument>
   </method>
   <method name="wavelengthInterval" >
    <description>At a given wavelength, return the wavelength interval in the tabulation of wavelength for which stellar spectra are defined.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavelength</argument>
   </method>
  </functionClass>
  !!]

end module Stellar_Population_Spectra
