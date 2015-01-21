!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements empirical models of conditional mass functions.

module Halo_Model_Power_Spectrum_Modifiers
  !% Implements empirical models of conditional mass functions.
  use ISO_Varying_String
  !# <include directive="haloModelPowerSpectrumModifier" type="functionModules" >
  include 'haloModelPowerSpectrumModifier.functionModules.inc'
  !# </include>
  implicit none
  private

  ! Labels for terms in the halo model.
  integer, public, parameter :: termOneHalo=1
  integer, public, parameter :: termTwoHalo=2

  !# <include directive="haloModelPowerSpectrumModifier" type="function" >
  !#  <descriptiveName>Halo Model Power Spectrum Modifier</descriptiveName>
  !#  <description>Object providing modifiers to the power spectra in the halo model.</description>
  !#  <default>null</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="modify" >
  !#   <description>Modify the power spectra in the halo model of clustering.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), dimension(:  )           :: wavenumber</argument>
  !#   <argument>integer         , intent(in   )                           :: term</argument>
  !#   <argument>double precision, intent(inout), dimension(:  )           :: powerSpectrum</argument>
  !#   <argument>double precision, intent(inout), dimension(:,:), optional :: powerSpectrumCovariance</argument>
  !#   <argument>double precision, intent(in   )                , optional :: mass</argument>
  !#  </method>
  include 'haloModelPowerSpectrumModifier.type.inc'
  !# </include>

end module Halo_Model_Power_Spectrum_Modifiers
