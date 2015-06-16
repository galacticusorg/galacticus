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

!% Contains a module that provides a class that implements modification of halo mass samples used in merger tree building.

module Merger_Trees_Mass_Function_Sampling_Modifiers
  !% Provides a class that implements modification of halo mass samples used in merger tree building.
  use ISO_Varying_String
  !# <include directive="massFunctionSamplingModifier" type="functionModules" >
  include 'massFunctionSamplingModifier.functionModules.inc'
  !# </include>
  private

  !# <include directive="massFunctionSamplingModifier" type="function" >
  !#  <descriptiveName>Mass Function Sampling Modifier</descriptiveName>
  !#  <description>Class providing modification of halo mass samples as used in merger tree building.</description>
  !#  <default>null</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="modify" >
  !#   <description>Modify a list of halo masses to be used in merger tree building.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, allocatable, dimension(:), intent(inout) :: treeHaloMass</argument>
  !#   <argument>double precision                           , intent(in   ) :: treeBaseTime</argument>
  !#  </method>
  include 'massFunctionSamplingModifier.type.inc'
  !# </include>

end module Merger_Trees_Mass_Function_Sampling_Modifiers
