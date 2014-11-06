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

!% Contains a module which provides an object that implements concentrations of dark matter halo profiles.

module Dark_Matter_Profiles_Concentration
  !% Provides an object that implements concentrations of dark matter halo profiles.
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  use               Galacticus_Nodes
  use               Virial_Density_Contrast
  use               Dark_Matter_Profiles
  use               FGSL
  !# <include directive="darkMatterProfileConcentration" type="functionModules" >
  include 'darkMatterProfileConcentration.functionModules.inc'
  !# </include>
  private

  !# <include directive="darkMatterProfileConcentration" type="function" >
  !#  <descriptiveName>Dark Matter Profile Concentrations</descriptiveName>
  !#  <description>Object providing dark matter profile concentrations.</description>
  !#  <default>gao2008</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>yes</stateful>
  !#  <method name="concentration" >
  !#   <description>Returns the concentration parameter for the given {\tt node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#  </method>
  !#  <method name="densityContrastDefinition" >
  !#   <description>Returns a {\tt virialDensityContrast} object describing the virial density contrast used to define this concentration.</description>
  !#   <type>class(virialDensityContrastClass)</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="darkMatterProfileDefinition" >
  !#   <description>Returns a {\tt darkMatterProfile} object describing the dark matter density profile used to define this concentration.</description>
  !#   <type>class(darkMatterProfileClass)</type>
  !#   <pass>yes</pass>
  !#  </method>
  include 'darkMatterProfileConcentration.type.inc'
  !# </include>
  
end module Dark_Matter_Profiles_Concentration
