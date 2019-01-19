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

!% Contains a module which implements a class for calculations of the initial radius in the dark matter halo for use when solving
!% for galactic structure.

module Galactic_Structure_Initial_Radii
  !% Implements a class for calculations of the initial radius in the dark matter halo for use when solving for galactic
  !% structure.
  use Galacticus_Nodes
  implicit none
  private

  !# <functionClass>
  !#  <name>galacticStructureRadiiInitial</name>
  !#  <descriptiveName>Initial Mass Functions</descriptiveName>
  !#  <description>Class providing stellar initial mass functions.</description>
  !#  <default>gnedin2004</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <calculationReset>yes</calculationReset>
  !#  <method name="radius" >
  !#   <description>Find the initial radius in the dark matter halo of {\normalfont \ttfamily node} corresponding to the given final {\normalfont \ttfamily radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !#  <method name="radiusDerivative" >
  !#   <description>Find the derivative of the initial radius in the dark matter halo of {\normalfont \ttfamily node} with respect to the final radius corresponding to the given final {\normalfont \ttfamily radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !# </functionClass>

end module Galactic_Structure_Initial_Radii
