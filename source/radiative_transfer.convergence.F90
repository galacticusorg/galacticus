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
Contains a module which provides convergence criteria for radiative transfer calculations.
!!}

module Radiative_Transfer_Convergences
  !!{
  Provides a class that implements convergence criteria for radiative transfer calculations.
  !!}
  use :: Radiative_Transfer_Matters      , only : radiativeTransferMatterClass      , radiativeTransferPropertiesMatter
  use :: Radiative_Transfer_Photon_Packet, only : radiativeTransferPhotonPacketClass
  private
  
  !![
  <enumeration>
   <name>statusCell</name>
   <description>Specifies whether a domain cell is the first cell, last cell, or an intermediate cell in a convergence test sweep, allowing implementations to perform special initialization or finalization actions at the domain boundaries.</description>
   <visibility>public</visibility>
   <entry label="first"/>
   <entry label="last" />
   <entry label="other"/>
  </enumeration>
  !!]
  
  !![
  <functionClass>
   <name>radiativeTransferConvergence</name>
   <descriptiveName>Radiative Transfer Convergence Criteria</descriptiveName>
   <description>Class providing convergence criteria for Monte Carlo radiative transfer calculations---diagnostics
    that assess whether the iterative solution (radiation field, ionization state, or dust temperature)
    has converged to a stable solution after sufficient photon packets have been propagated. Methods
    process escaping photon packets, test individual domain cells for convergence, and determine whether
    the entire grid has reached the required accuracy. Multiple iterations of photon propagation may
    be needed before convergence is declared and the final radiation field is recorded.</description>
   <default>always</default>
   <method name="photonPacketEscapes" >
    <description>Process a photon packet that has escaped the computational domain, accumulating any statistics needed by the convergence criterion before the packet is discarded.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(radiativeTransferPhotonPacketClass), intent(inout) :: photonPacket</argument>
   </method>
   <method name="testConvergence" >
    <description>Test whether a single domain cell has converged by comparing its current radiation field or matter state against the previous iteration, setting the \mono{converged} flag accordingly.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class  (radiativeTransferMatterClass     ), intent(inout) :: radiativeTransferMatter_</argument>
    <argument>class  (radiativeTransferPropertiesMatter), intent(inout) :: properties</argument>
    <argument>type   (enumerationStatusCellType        ), intent(in   ) :: statusCell</argument>
    <argument>logical                                   , intent(  out) :: converged</argument>
   </method>
  </functionClass>
  !!]

end module Radiative_Transfer_Convergences
