!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which provides a class for calculations of dark matter halo mass accretion histories.

module Dark_Matter_Halo_Mass_Accretion_Histories
  !% Provides a class for calculations of dark matter halo mass accretion histories.
  use :: Galacticus_Nodes, only : treeNode
  private

  !# <functionClass>
  !#  <name>darkMatterHaloMassAccretionHistory</name>
  !#  <descriptiveName>Dark Matter Halo Mass Accretion Histories</descriptiveName>
  !#  <description>
  !#   Class providing dark matter halo mass accretion histories.
  !#  </description>
  !#  <default>wechsler2002</default>
  !#  <method name="time">
  !#   <description>Returns the time at which the given halo mass was reached.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: mass</argument>
  !#  </method>
  !#  <method name="massAccretionRate">
  !#   <description>Returns the mass accretion rate at the specified time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: time</argument>
  !#  </method>
  !# </functionClass>

end module Dark_Matter_Halo_Mass_Accretion_Histories
