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
Contains a module which implements a class of modifiers of the power spectrum for the halo model.
!!}

module Halo_Model_Power_Spectrum_Modifiers
  !!{
  Implements a class of modifiers of the power spectrum for the halo model.
  !!}
  private

  !![
  <enumeration>
   <name>haloModelTerm</name>
   <description>Enumeration of terms in the halo model.</description>
   <visibility>public</visibility>
   <entry label="oneHalo"/>
   <entry label="twoHalo"/>
  </enumeration>
  !!]

  !![
  <functionClass>
   <name>haloModelPowerSpectrumModifier</name>
   <descriptiveName>Halo Model Power Spectrum Modifier</descriptiveName>
   <description>
    Class providing modifiers to the one-halo and two-halo power spectrum terms in halo model galaxy
    clustering calculations. The \mono{modify} method accepts the raw halo-model power spectrum at
    each wavenumber and modifies it in place for a specified term (one-halo or two-halo), allowing
    corrections such as baryonic feedback suppression, assembly bias, or scale-dependent effects to
    be applied without altering the underlying halo occupation distribution. Implementations may be
    mass-dependent (for per-halo corrections) or act globally on the full power spectrum array.
   </description>
   <default>identity</default>
   <method name="modify" >
    <description>Modify the power spectra in the halo model of clustering.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision                              , intent(in   ), dimension(:  )           :: wavenumber</argument>
    <argument>type            (enumerationHaloModelTermType), intent(in   )                           :: term</argument>
    <argument>double precision                              , intent(inout), dimension(:  )           :: powerSpectrum</argument>
    <argument>double precision                              , intent(inout), dimension(:,:), optional :: powerSpectrumCovariance</argument>
    <argument>double precision                              , intent(in   )                , optional :: mass</argument>
   </method>
  </functionClass>
  !!]

end module Halo_Model_Power_Spectrum_Modifiers
