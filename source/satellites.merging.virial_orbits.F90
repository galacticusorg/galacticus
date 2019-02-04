!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which provides a class implementing satellite orbital parameters at virial radius crossing.

module Virial_Orbits
  !% Provides a class implementing satellite orbital parameters at virial radius crossing.
  use Galacticus_Nodes       , only : treeNode
  use Virial_Density_Contrast
  use Kepler_Orbits
  private

  !# <functionClass>
  !#  <name>virialOrbit</name>
  !#  <descriptiveName>Virial Orbits</descriptiveName>
  !#  <description>Class providing orbital parameters of satellite halos at the time of merging.</description>
  !#  <default>benson2005</default>
  !#  <method name="orbit" >
  !#   <description>Returns an orbit object.</description>
  !#   <type>type(keplerOrbit)</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>type(treeNode), intent(inout) :: node               , host</argument>
  !#   <argument>logical       , intent(in   ) :: acceptUnboundOrbits</argument>
  !#  </method>
  !#  <method name="densityContrastDefinition" >
  !#   <description>Returns a {\normalfont \ttfamily virialDensityContrast} object describing the virial density contrast used to define this orbit class.</description>
  !#   <type>class(virialDensityContrastClass)</type>
  !#   <pass>yes</pass>
  !#  </method>
  !# </functionClass>

end module Virial_Orbits
