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

!% Contains a module with the standard implementation of N-body component method.

module Node_Component_NBody_Standard
  !% The standard implementation of N-body component method.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_NBody_Standard_Promote

  !# <component>
  !#  <class>nBody</class>
  !#  <name>standard</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>velocityMaximum</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="kilo" comment="Maximum velocity of the halo rotation curve."/>
  !#   </property>
  !#   <property>
  !#     <name>velocityDispersion</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="kilo" comment="Velocity dispersion of the dark matter halo."/>
  !#   </property>
  !#   <property>
  !#     <name>particleCount</name>
  !#     <type>longInteger</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="0.0d0" comment="Number of N-body particles in this halo."/>
  !#   </property>
  !#  </properties>
  !# </component>

contains

   !# <nodePromotionTask>
   !#  <unitName>Node_Component_NBody_Standard_Promote</unitName>
   !# </nodePromotionTask>
   subroutine Node_Component_NBody_Standard_Promote(thisNode)
     !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the properties of {\tt thisNode}
     !% to be those of its parent.
     implicit none
     type (treeNode          ), intent(inout), pointer :: thisNode
     type (treeNode          )               , pointer :: parentNode
     class(nodeComponentNBody)               , pointer :: parentNBody, thisNBody

     ! Get the nBody component.
     thisNBody => thisNode%nBody()
     ! Ensure that it is of the standard class.
     select type (thisNBody)
     class is (nodeComponentNBodyStandard)
        ! Get the parent node and its nBody component.
        parentNode  => thisNode  %parent
        parentNBody => parentNode%nBody()
        ! Copy properties from the parent.
        call thisNBody%   velocityMaximumSet(parentNBody%velocityMaximum   ())
        call thisNBody%velocityDispersionSet(parentNBody%velocityDispersion())
        call thisNBody%     particleCountSet(parentNBody%particleCount     ())
     end select
     return
   end subroutine Node_Component_NBody_Standard_Promote

end module Node_Component_NBody_Standard
