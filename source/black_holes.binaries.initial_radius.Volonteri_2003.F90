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


!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements a black hole binary initial separation which follows that of \cite{volonteri_assembly_2003}.

module Black_Hole_Binary_Initial_Radii_Volonteri_2003
  !% Implements a black hole binary initial separation which follows that of \cite{volonteri_assembly_2003}.
  implicit none
  private
  public :: Black_Hole_Binary_Initial_Radii_Volonteri_2003_Initialize

contains

  !# <blackHoleBinaryInitialRadiiMethod>
  !#  <unitName>Black_Hole_Binary_Initial_Radii_Volonteri_2003_Initialize</unitName>
  !# </blackHoleBinaryInitialRadiiMethod>
  subroutine Black_Hole_Binary_Initial_Radii_Volonteri_2003_Initialize(blackHoleBinaryInitialRadiiMethod,Black_Hole_Binary_Initial_Radius_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: blackHoleBinaryInitialRadiiMethod
    procedure(double precision), pointer, intent(inout) :: Black_Hole_Binary_Initial_Radius_Get
    
    if (blackHoleBinaryInitialRadiiMethod == 'Volonteri 2003') Black_Hole_Binary_Initial_Radius_Get => Black_Hole_Binary_Initial_Radius_Volonteri_2003
    return
  end subroutine Black_Hole_Binary_Initial_Radii_Volonteri_2003_Initialize

  double precision function Black_Hole_Binary_Initial_Radius_Volonteri_2003(thisNode,hostNode)
    !% Returns an initial separation for binary black holes using the method of \cite{volonteri_assembly_2003}, with the assumption that
    !% the local velocity dispersion is approximately the dark matter halo virial velocity.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode,hostNode

    Black_Hole_Binary_Initial_Radius_Volonteri_2003=gravitationalConstantGalacticus*(       Tree_Node_Black_Hole_Mass(thisNode)    &
       &                                                                           +        Tree_Node_Black_Hole_Mass(hostNode)    &
       &                                                                            )                                              &
       &                                                             /(2.0d0       * Dark_Matter_Halo_Virial_Velocity(hostNode)**2 &
       &                                                              )

    return
  end function Black_Hole_Binary_Initial_Radius_Volonteri_2003

end module Black_Hole_Binary_Initial_Radii_Volonteri_2003
