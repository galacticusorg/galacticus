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

!% Contains a module which implements a class implementing satellite orbital parameters at virial radius crossing.

module Virial_Orbits
  !% Implements a class implementing satellite orbital parameters at virial radius crossing.
  use ISO_Varying_String
  use Galacticus_Nodes
  use Virial_Density_Contrast
  use Kepler_Orbits
  !# <include directive="virialOrbit" type="functionModules" >
  include 'virialOrbit.functionModules.inc'
  !# </include>
  private

  !# <include directive="virialOrbit" type="function" >
  !#  <descriptiveName>Virial Orbits</descriptiveName>
  !#  <description>Class providing orbital parameters of satellite halos at the time of merging.</description>
  !#  <default>benson2005</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>yes</stateful>
  !#  <method name="orbit" >
  !#   <description>Returns an orbit object.</description>
  !#   <type>type(keplerOrbit)</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node, host</argument>
  !#   <argument>logical       , intent(in   )          :: acceptUnboundOrbits</argument>
  !#  </method>
  !#  <method name="densityContrastDefinition" >
  !#   <description>Returns a {\tt virialDensityContrast} object describing the virial density contrast used to define this orbit class.</description>
  !#   <type>class(virialDensityContrastClass)</type>
  !#   <pass>yes</pass>
  !#  </method>
   include 'virialOrbit.type.inc'
  !# </include>

end module Virial_Orbits
