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

!% Contains a module which implements gravitational lensing from large scale structure.

module Gravitational_Lensing
  !% Implements gravitational lensing from large scale structure.
  use Tables
  use ISO_Varying_String
  use, intrinsic :: ISO_C_Binding
  !# <include directive="gravitationalLensing" type="functionModules" >
  include 'gravitationalLensing.functionModules.inc'
  !# </include>
  implicit none
  private

  !# <include directive="gravitationalLensing" type="function" >
  !#  <descriptiveName>Gravitational Lensing</descriptiveName>
  !#  <description>Object providing gravitational lensing probabilities due to large scale structure.</description>
  !#  <default>takahashi2011</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="magnificationPDF" >
  !#   <description>Returns the differential probability function for magnification.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: magnification, redshift, scaleSource</argument>
  !#  </method>
  !#  <method name="magnificationCDF" >
  !#   <description>Returns the cumulative probability function for magnification.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: magnification, redshift, scaleSource</argument>
  !#  </method>
  include 'gravitationalLensing.type.inc'
  !# </include>

end module Gravitational_Lensing
