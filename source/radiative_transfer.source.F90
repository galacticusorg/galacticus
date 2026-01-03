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
Contains a module which provides photon sources for radiative transfer calculations.
!!}

module Radiative_Transfer_Sources
  !!{
  Provides a class that implements photon sources for radiative transfer calculations.
  !!}
  use :: Radiative_Transfer_Photon_Packet, only : radiativeTransferPhotonPacketClass
  private

  !![
  <functionClass>
   <name>radiativeTransferSource</name>
   <descriptiveName>Radiative Transfer Photon Sources</descriptiveName>
   <description>Class providing photon sources for radiative transfer calculations.</description>
   <default>point</default>
   <method name="initializePhotonPacket" >
    <description>Initialize a photon packet from the source.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
   </method>
   <method name="spectrum" >
    <description>Return the spectrum (in units of $L_\odot$ \AA$^{-1}$) of the source at the given {\normalfont \ttfamily wavelength}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: wavelength</argument>
    <argument>integer         , intent(in   ), optional :: sourceType</argument>
   </method>
   <method name="luminosity" >
    <description>Return the luminosity (in units of $L_\odot$) of the source at the given wavelength range.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: wavelengthMinimum, wavelengthMaximum</argument>
    <argument>integer         , intent(in   ), optional :: sourceType</argument>
   </method>
   <method name="sourceTypeCount" >
    <description>Return the number of source types.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="sourceTypeName" >
    <description>Return the name of the requested of source type.</description>
    <type>type(varying_string)</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: sourceType</argument>
   </method>
  </functionClass>
  !!]
  
end module Radiative_Transfer_Sources
