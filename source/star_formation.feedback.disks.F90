!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which provides an object that implements cosmological parameters.

module Star_Formation_Feedback_Disks
  !% Provides an object that implements calculations of feedback from star formation in disks.
  use Galacticus_Nodes
  
  !# <functionClass>
  !#  <name>starFormationFeedbackDisks</name>
  !#  <descriptiveName>Feedback from star formation in disks</descriptiveName>
  !#  <description>Class providing models of feedback from star formation in disks.</description>
  !#  <default>powerLaw</default>
  !#  <method name="outflowRate" >
  !#   <description>Returns the outflow rate due to star formation in the disk component of {\normalfont \ttfamily node} in units of $M_\odot/$Gyr.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: rateEnergyInput, rateStarFormation</argument>
  !#  </method>
  !# </functionClass>

end module Star_Formation_Feedback_Disks
