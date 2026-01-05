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
Contains a module which provides matter types for radiative transfer calculations.
!!}

module Radiative_Transfer_Matters
  !!{
  Provides a class that implements matter types for radiative transfer calculations.
  !!}
  use, intrinsic :: ISO_C_Binding                          , only : c_size_t
  use            :: ISO_Varying_String                     , only : varying_string
  use            :: Computational_Domain_Volume_Integrators, only : computationalDomainVolumeIntegratorClass
  use            :: Radiative_Transfer_Photon_Packet       , only : radiativeTransferPhotonPacketClass
  private
  public :: radiativeTransferPropertiesMatter

  type, abstract :: radiativeTransferPropertiesMatter
     !!{
     Parent class for properties of matter used in radiative transfer calculations.
     !!}
   contains
  end type radiativeTransferPropertiesMatter
  
  !![
  <functionClass>
   <name>radiativeTransferMatter</name>
   <descriptiveName>Radiative Transfer Matter</descriptiveName>
   <description>Class providing matter types for radiative transfer calculations.</description>
   <default>atomic</default>
   <method name="propertyClass" >
    <description>Return the property class to be used.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter), intent(inout), allocatable :: properties</argument>
   </method>
   <method name="populateDomain" >
    <description>Populate a domain cell with matter.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter       ), intent(inout) :: properties</argument>
    <argument>class(computationalDomainVolumeIntegratorClass), intent(inout) :: integrator</argument>
    <argument>logical                                        , intent(in   ) :: onProcess</argument>
   </method>
   <method name="broadcastDomain" >
    <description>Broadcast populated domain to all MPI processes.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer                                   , intent(in   ) :: sendFromProcess</argument>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
   </method>
   <method name="reset" >
    <description>Reset the matter prior to a new iteration.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter       ), intent(inout) :: properties</argument>
   </method>
   <method name="absorptionCoefficient" >
    <description>Return the absorption coefficient of the matter.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter ), intent(inout) :: properties</argument>
    <argument>class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
   </method>
   <method name="accumulatePhotonPacket" >
    <description>Return the absorption coefficient of the matter.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class           (radiativeTransferPropertiesMatter ), intent(inout) :: properties</argument>
    <argument>class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
    <argument>double precision                                    , intent(in   ) :: absorptionCoefficient, lengthTraversed</argument>
   </method>
   <method name="interactWithPhotonPacket" >
    <description>Interact with a photon packet.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter ), intent(inout) :: properties</argument>
    <argument>class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
   </method>
   <method name="stateSolve" >
    <description>Solve for the state of the matter.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout), target   :: properties</argument>
    <argument>integer                                   , intent(  out), optional :: status</argument>
   </method>
   <method name="convergenceMeasure" >
    <description>Return a convergence measure for the matter.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
   </method>
   <method name="outputProperty" >
    <description>Return the value of the output property.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
    <argument>integer(c_size_t                         ), intent(in   ) :: output</argument>
   </method>
   <method name="countOutputs" >
    <description>Return the count of the number of properties to output.</description>
    <type>integer(c_size_t)</type>
    <pass>yes</pass>
   </method>
   <method name="outputName" >
    <description>Return the name of the indexed output property.</description>
    <type>type(varying_string)</type>
    <pass>yes</pass>
    <argument>integer(c_size_t), intent(in   ) :: output</argument>
   </method>
   <method name="accumulationReduction" >
    <description>Perform reduction across MPI processes of accumulated properties.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
   </method>
   <method name="broadcastState" >
    <description>Broadcast domain state to all MPI processes.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer                                   , intent(in   ) :: sendFromProcess</argument>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
   </method>
  </functionClass>
  !!]
  
end module Radiative_Transfer_Matters
