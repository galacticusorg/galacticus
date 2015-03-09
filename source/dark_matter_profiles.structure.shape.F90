!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

module Dark_Matter_Profiles_Shape
  !% Provides an object that implements shape parameters of dark matter halo profiles.
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  use               Galacticus_Nodes
  !# <include directive="darkMatterProfileShape" type="functionModules" >
  include 'darkMatterProfileShape.functionModules.inc'
  !# </include>
  private

  !# <include directive="darkMatterProfileShape" type="function" >
  !#  <descriptiveName>Dark Matter Profile Shapes</descriptiveName>
  !#  <description>Object providing dark matter profile shape parameters.</description>
  !#  <default>gao2008</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>no</stateful>
  !#  <method name="shape" >
  !#   <description>Returns the shape parameter for the given {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#  </method>
  include 'darkMatterProfileShape.type.inc'
  !# </include>
  
end module Dark_Matter_Profiles_Shape
