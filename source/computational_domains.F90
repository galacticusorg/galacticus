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
Contains a module which provides computational domains.
!!}

module Computational_Domains
  !!{
  Provides a class that implements computational domains.
  !!}
  use            :: IO_HDF5                         , only : hdf5Object
  use, intrinsic :: ISO_C_Binding                   , only : c_size_t
  use            :: Radiative_Transfer_Photon_Packet, only : radiativeTransferPhotonPacketClass
  private
  public :: domainIterator

  type, abstract :: domainIterator
     !!{
     Parent class for iterators over computational domains.
     !!}
   contains
     procedure(domainIteratorNext), deferred :: next
  end type domainIterator
  
  !![
  <functionClass>
   <name>computationalDomain</name>
   <descriptiveName>Computational Domains</descriptiveName>
   <description>Class providing computational domains.</description>
   <default>cartesian3D</default>
   <method name="iterator" >
    <description>Return an iterator which can be used to iterate over the domain.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(domainIterator), intent(inout), allocatable :: iterator</argument>
   </method>
   <method name="initialize" >
    <description>Initialize the computational domain.</description>
    <type>void</type>
    <pass>yes</pass>
   </method>
   <method name="reset" >
    <description>Reset computational domain prior to a new iteration.</description>
    <type>void</type>
    <pass>yes</pass>
   </method>
   <method name="indicesFromPosition" >
    <description>Return a set of indices identifying the domain cell containing a position. For points outside the domain values of {\normalfont \ttfamily -huge(0)} will be returned.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision                       , dimension(3), intent(in   ) :: position</argument>
    <argument>integer         (c_size_t), allocatable, dimension(:), intent(inout) :: indices</argument>
   </method>
   <method name="absorptionCoefficient" >
    <description>Return the absorption coefficient in the domain cell.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class  (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket</argument>
    <argument>integer(c_size_t                          ), dimension(:), intent(in   ) :: indices</argument>
   </method>
   <method name="lengthToCellBoundary" >
    <description>Return the length until a photon packet hits a cell boundary. Also return the indices of that neighboring cell.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class           (radiativeTransferPhotonPacketClass)                           , intent(inout) :: photonPacket</argument>
    <argument>integer         (c_size_t                          )             , dimension(:), intent(in   ) :: indices</argument>
    <argument>integer         (c_size_t                          ), allocatable, dimension(:), intent(inout) :: indicesNeighbor</argument>
    <argument>double precision                                                 , dimension(3), intent(  out) :: positionBoundary</argument>
   </method>
   <method name="accumulatePhotonPacket" >
    <description>Accumulate ``absorptions'' from the photon packet as it traverses a cell of the computational domain.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class           (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket</argument>
    <argument>integer         (c_size_t                          ), dimension(:), intent(in   ) :: indices</argument>
    <argument>double precision                                                  , intent(in   ) :: absorptionCoefficient, lengthTraversed</argument>
   </method>
   <method name="interactWithPhotonPacket" >
    <description>Interact the photon packet with matter in the computational domain, possibly modifying its properties. Return true if the photon packet is still alive, false if it is absorbed.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>class           (radiativeTransferPhotonPacketClass)              , intent(inout) :: photonPacket</argument>
    <argument>integer         (c_size_t                          ), dimension(:), intent(inout) :: indices</argument>
   </method>
   <method name="stateSolve" >
    <description>Solve for the state of matter in the computational domain.</description>
    <type>void</type>
    <pass>yes</pass>
   </method>
   <method name="converged" >
    <description>Return true if the computational domain is converged.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="output" >
    <description>Output the computational domain.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(hdf5Object), intent(inout) :: outputGroup</argument>
   </method>
  </functionClass>
  !!]
  
  abstract interface
     logical function domainIteratorNext(self)
       !!{
       Function which steps a computational domain iterator to the next iteration, returning false if no more iterations
       remain.
       !!}
       import domainIterator
       class(domainIterator), intent(inout) :: self
     end function domainIteratorNext
  end interface

end module Computational_Domains
