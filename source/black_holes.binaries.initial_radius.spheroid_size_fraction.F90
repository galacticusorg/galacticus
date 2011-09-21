!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements a black hole binary initial separation which is a fixed fraction of the scale radius of the
!% larger of the host and satellite spheroids.

module Black_Hole_Binary_Initial_Radii_Spheroid_Size
  !% Implements a black hole binary initial separation which is a fixed fraction of the scale radius of the larger of the host and
  !% satellite spheroids.
  implicit none
  private
  public :: Black_Hole_Binary_Initial_Radii_Spheroid_Size_Initialize

  ! The fraction of the spheroid radius at which merging black holes should be placed.
  double precision :: blackHoleInitialRadiusSpheroidRadiusRatio

contains

  !# <blackHoleBinaryInitialRadiiMethod>
  !#  <unitName>Black_Hole_Binary_Initial_Radii_Spheroid_Size_Initialize</unitName>
  !# </blackHoleBinaryInitialRadiiMethod>
  subroutine Black_Hole_Binary_Initial_Radii_Spheroid_Size_Initialize(blackHoleBinaryInitialRadiiMethod,Black_Hole_Binary_Initial_Radius_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: blackHoleBinaryInitialRadiiMethod
    procedure(double precision), pointer, intent(inout) :: Black_Hole_Binary_Initial_Radius_Get
    
    if (blackHoleBinaryInitialRadiiMethod == 'spheroidRadiusFraction') then
       Black_Hole_Binary_Initial_Radius_Get => Black_Hole_Binary_Initial_Radius_Spheroid_Size
       !@ <inputParameter>
       !@   <name>blackHoleInitialRadiusSpheroidRadiusRatio</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The fraction of the spheroid radius at which merging black holes will be initially placed.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('blackHoleInitialRadiusSpheroidRadiusRatio',blackHoleInitialRadiusSpheroidRadiusRatio,defaultValue=0.0d0)
    end if
    return
  end subroutine Black_Hole_Binary_Initial_Radii_Spheroid_Size_Initialize

  double precision function Black_Hole_Binary_Initial_Radius_Spheroid_Size(thisNode,hostNode)
    !% Returns an initial separation for a binary black holes that is a fixed fraction of the scale radius of the larger of the
    !% host and satellite spheroids.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode,hostNode

    Black_Hole_Binary_Initial_Radius_Spheroid_Size=blackHoleInitialRadiusSpheroidRadiusRatio*max(                                     &
         &                                                                                        Tree_Node_Spheroid_Radius(thisNode) &
         &                                                                                       ,Tree_Node_Spheroid_Radius(hostNode) &
         &                                                                                      )
    return
  end function Black_Hole_Binary_Initial_Radius_Spheroid_Size

end module Black_Hole_Binary_Initial_Radii_Spheroid_Size
