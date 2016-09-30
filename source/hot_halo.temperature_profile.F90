!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which provides an object that implements hot halo temperature profiles.

module Hot_Halo_Temperature_Profiles
  !% Provides an object that implements hot halo temperature profiles.
  use ISO_Varying_String
  use Mass_Distributions
  use Galacticus_Nodes
  private

  !# <include directive="hotHaloTemperatureProfile" type="function" >
  !#  <descriptiveName>Hot Halo Temperature profiles</descriptiveName>
  !#  <description>Class implementing hot halo temperarture profiles.</description>
  !#  <default>virial</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="temperature" >
  !#   <description>Return the temperature of the hot halo at the given {\normalfont \ttfamily radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !#  <method name="temperatureLogSlope" >
  !#   <description>Return the logarithmic slope of the temperature of the hot halo at the given {\normalfont \ttfamily radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: radius</argument>
  !#  </method>
  include 'hotHaloTemperatureProfile.type.inc'
  !# </include>

end module Hot_Halo_Temperature_Profiles
