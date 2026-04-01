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
   <description>Class providing matter types for Monte Carlo radiative transfer calculations---the physical
    description of the gas and dust that photon packets interact with as they propagate through the
    computational domain. Methods populate domain cells with matter properties (density, temperature,
    ionization fractions), compute the absorption coefficient at each packet wavelength, accumulate
    absorbed energy from traversing packets, and update the matter state (e.g.\ photoionization
    equilibrium) after each iteration. The default implementation models atomic hydrogen gas.</description>
   <default>atomic</default>
   <method name="propertyClass" >
    <description>Allocate and return the concrete \refClass{radiativeTransferPropertiesMatter} object used to store per-cell matter properties (density, temperature, ionization fractions) in the computational domain.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter), intent(inout), allocatable :: properties</argument>
   </method>
   <method name="populateDomain" >
    <description>Populate the given domain cell with matter properties (density, composition, initial temperature) by querying the volume integrator, setting up the local state for subsequent photon-packet propagation.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter       ), intent(inout) :: properties</argument>
    <argument>class(computationalDomainVolumeIntegratorClass), intent(inout) :: integrator</argument>
    <argument>logical                                        , intent(in   ) :: onProcess</argument>
   </method>
   <method name="broadcastDomain" >
    <description>Broadcast the populated domain cell properties from the specified MPI process to all other processes, ensuring all ranks share a consistent initial matter state before photon propagation begins.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer                                   , intent(in   ) :: sendFromProcess</argument>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
   </method>
   <method name="reset" >
    <description>Reset the accumulated radiation-field tallies (absorbed energy, photon counts) stored in the cell properties to zero in preparation for a new Monte Carlo iteration of photon propagation.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter       ), intent(inout) :: properties</argument>
   </method>
   <method name="absorptionCoefficient" >
    <description>Return the absorption coefficient (in units of inverse length) of the matter in the given cell at the wavelength of the photon packet, used to compute the mean free path during packet propagation.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter ), intent(inout) :: properties</argument>
    <argument>class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
   </method>
   <method name="accumulatePhotonPacket" >
    <description>Accumulate the energy deposited by the photon packet into the cell properties, recording the absorbed luminosity and packet-crossing statistics for the current Monte Carlo iteration.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class           (radiativeTransferPropertiesMatter ), intent(inout) :: properties</argument>
    <argument>class           (radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
    <argument>double precision                                    , intent(in   ) :: absorptionCoefficient, lengthTraversed</argument>
   </method>
   <method name="interactWithPhotonPacket" >
    <description>Attempt a physical interaction (e.g.\ absorption, scattering, or re-emission) between the matter and the photon packet, updating the packet properties and returning true if an interaction occurred.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter ), intent(inout) :: properties</argument>
    <argument>class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
   </method>
   <method name="stateSolve" >
    <description>Solve for the equilibrium matter state (e.g.\ photoionization balance, dust temperature) in the given cell based on the radiation field accumulated during the current iteration, updating the cell properties in place.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout), target   :: properties</argument>
    <argument>integer                                   , intent(  out), optional :: status</argument>
   </method>
   <method name="convergenceMeasure" >
    <description>Return a scalar convergence measure (e.g.\ the fractional change in ionization fraction or temperature) for the given cell between the current and previous iterations, used to assess global radiative transfer convergence.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
   </method>
   <method name="outputProperty" >
    <description>Return the value of the indexed output property (e.g.\ ionization fraction, temperature, or density) for the given cell, as written to the output HDF5 file after convergence.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
    <argument>integer(c_size_t                         ), intent(in   ) :: output</argument>
   </method>
   <method name="countOutputs" >
    <description>Return the total number of scalar output properties provided by this matter implementation, used to allocate output arrays before iterating over individual property names and values.</description>
    <type>integer(c_size_t)</type>
    <pass>yes</pass>
   </method>
   <method name="outputName" >
    <description>Return the name of the output property at the given index, used as the dataset name when writing per-cell matter properties to the HDF5 output file.</description>
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
    <description>Broadcast the solved matter state (ionization fractions, temperatures, and other cell properties) from the specified MPI process to all other processes after the state-solve step.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer                                   , intent(in   ) :: sendFromProcess</argument>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
   </method>
  </functionClass>
  !!]
  
end module Radiative_Transfer_Matters
