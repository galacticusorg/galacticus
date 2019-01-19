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

!% Contains a module which provides a class that implements ram pressure stripping in spheroids.

module Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids
  !% Provides a class that implements calculations of ram pressure stripping in spheroids.
  use Galacticus_Nodes, only : treeNode
  
  !# <functionClass>
  !#  <name>ramPressureStrippingSpheroids</name>
  !#  <descriptiveName>Ram pressure stripping in spheroids</descriptiveName>
  !#  <description>Class providing models of ram pressure stripping in spheroids.</description>
  !#  <default>null</default>
  !#  <method name="rateMassLoss" >
  !#   <description>Returns the rate of mass loss (in $M_\odot$~Gyr$^{-1}$) due to ram pressure stripping of the spheroid component of {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#  </method>
  !# </functionClass>
  
end module Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids
