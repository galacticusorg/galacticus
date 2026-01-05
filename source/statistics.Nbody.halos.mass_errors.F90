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
Contains a module which provides a class that implements errors on dark matter halo masses in
N-body simulations.
!!}

module Statistics_NBody_Halo_Mass_Errors
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>nbodyHaloMassError</name>
   <descriptiveName>N-body Halo Mass Errors</descriptiveName>
   <description>Class providing models of errors on N-body halo masses.</description>
   <default>null</default>
   <method name="errorFractional" >
    <description>Return the fractional error on the mass of an N-body halo corresponding to the given {\normalfont \ttfamily \gls{node}}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="correlation" >
    <description>Return the correlation in the  error on the mass of a pair of N-body halos corresponding to the given {\normalfont \ttfamily node1} and {\normalfont \ttfamily node2}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node1, node2</argument>
   </method>
   <method name="errorZeroAlways" >
    <description>Return {\normalfont \ttfamily true} if the mass error is always zero for any halo.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     nbodyHaloMassErrorErrorZeroAlways=.false.
    </code>
   </method>
  </functionClass>
  !!]

end module Statistics_NBody_Halo_Mass_Errors
