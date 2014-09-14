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

!% Contains a module which provides a class that implements reincorportation of outflowed mass into the hot halo.

module Hot_Halo_Outflows_Reincorporations
  !% Provides a class that implements reincorportation of outflowed mass into the hot halo.
  use ISO_Varying_String
  use Galacticus_Nodes
  !# <include directive="hotHaloOutflowReincorporation" type="functionModules" >
  include 'hotHaloOutflowReincorporation.functionModules.inc'
  !# </include>
  private

  !# <include directive="hotHaloOutflowReincorporation" type="function" >
  !#  <descriptiveName>Hot Halo Outflow Reincorporation</descriptiveName>
  !#  <description>Class implementing reincorportation of outflowed mass into the hot halo.</description>
  !#  <default>velocityMaximumScaling</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="rate" >
  !#   <description>Return the rate at which outflowed mass is being reincorporated into the hot halo.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#  </method>
  include 'hotHaloOutflowReincorporation.type.inc'
  !# </include>

end module Hot_Halo_Outflows_Reincorporations
