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

module Conditional_Mass_Functions
  !% Implements empirical models of conditional mass functions.
  use ISO_Varying_String
  !# <include directive="conditionalMassFunction" type="functionModules" >
  include 'conditionalMassFunction.functionModules.inc'
  !# </include>
  implicit none
  private

  integer, parameter, public :: haloModelGalaxyTypeAll      =0
  integer, parameter, public :: haloModelGalaxyTypeCentral  =1
  integer, parameter, public :: haloModelGalaxyTypeSatellite=2

  !# <include directive="conditionalMassFunction" type="function" >
  !#  <descriptiveName>Conditional Mass Function</descriptiveName>
  !#  <description>Object providing empirical models of conditional mass functions.</description>
  !#  <default>behroozi2010</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="massFunction" >
  !#   <description>Return the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm halo})$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: massHalo  , mass</argument>
  !#   <argument>integer         , intent(in   ), optional :: galaxyType</argument>
  !#  </method>
  !#  <method name="massFunctionVariance" >
  !#   <description>Return the variance in the cumulative conditional mass function, $\langle N(M_\star|M_{\rm halo}) \rangle \equiv \phi(M_\star|M_{\rm halo})$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: massHalo,massLow,massHigh</argument>
  !#  </method>
  include 'conditionalMassFunction.type.inc'
  !# </include>

end module Conditional_Mass_Functions
