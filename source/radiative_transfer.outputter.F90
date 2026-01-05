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
Contains a module which provides photon outputters for radiative transfer calculations.
!!}

module Radiative_Transfer_Outputters
  !!{
  Provides a class that implements outputters for radiative transfer calculations.
  !!}
  use :: IO_HDF5                         , only : hdf5Object
  use :: Radiative_Transfer_Photon_Packet, only : radiativeTransferPhotonPacketClass
  use :: Radiative_Transfer_Sources      , only : radiativeTransferSourceClass
  private

  !![
  <functionClass>
   <name>radiativeTransferOutputter</name>
   <descriptiveName>Radiative Transfer Outputters</descriptiveName>
   <description>Class providing outputters for radiative transfer calculations.</description>
   <default>null</default>
   <method name="reset" >
    <description>Reset the outputter.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
    </code>
   </method>
   <method name="sourceProperties" >
    <description>Output properties of the source distribution.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferSourceClass), intent(inout) :: radiativeTransferSource_</argument>
    <argument>type (hdf5Object                  ), intent(inout) :: outputGroup</argument>
    <code>
     !$GLC attributes unused :: self, radiativeTransferSource_, outputGroup
    </code>
   </method>
   <method name="photonPacketEscapes" >
    <description>Handle cases where a photon packet escapes the computational domain.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
    <code>
     !$GLC attributes unused :: self, photonPacket
    </code>
   </method>
   <method name="finalize" >
    <description>Finalize any calculations output.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
    </code>
   </method>
   <method name="output" >
    <description>Perform final output.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type (hdf5Object), intent(inout) :: outputGroup</argument>
    <code>
     !$GLC attributes unused :: self, outputGroup
    </code>
   </method>
  </functionClass>
  !!]

end module Radiative_Transfer_Outputters
