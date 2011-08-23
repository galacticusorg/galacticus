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

!% Contains a module which implements a black hole binary initial separation based on tidal disruption of the satellite galaxy.

module Black_Hole_Binary_Initial_Radii_Tidal_Radius
  !% Implements a black hole binary initial separation based on tidal disruption of the satellite galaxy.
  use Tree_Nodes
  implicit none
  private
  public :: Black_Hole_Binary_Initial_Radii_Tidal_Radius_Initialize

  ! Variables used in root finding.
  double precision         :: radiusHalfMass,massHalf
  type(treeNode),  pointer :: activeNode
  !$omp threadprivate(radiusHalfMass,massHalf,activeNode)

contains

  !# <blackHoleBinaryInitialRadiiMethod>
  !#  <unitName>Black_Hole_Binary_Initial_Radii_Tidal_Radius_Initialize</unitName>
  !# </blackHoleBinaryInitialRadiiMethod>
  subroutine Black_Hole_Binary_Initial_Radii_Tidal_Radius_Initialize(blackHoleBinaryInitialRadiiMethod,Black_Hole_Binary_Initial_Radius_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: blackHoleBinaryInitialRadiiMethod
    procedure(double precision), pointer, intent(inout) :: Black_Hole_Binary_Initial_Radius_Get
    
    if (blackHoleBinaryInitialRadiiMethod == 'tidal radius') Black_Hole_Binary_Initial_Radius_Get => Black_Hole_Binary_Initial_Radius_Tidal_Radius
    return
  end subroutine Black_Hole_Binary_Initial_Radii_Tidal_Radius_Initialize

  double precision function Black_Hole_Binary_Initial_Radius_Tidal_Radius(thisNode,hostNode)
    !% Returns an initial separation for a binary black holes through tidal disruption.
    use, intrinsic :: ISO_C_Binding
    use ISO_Varying_String
    use Galactic_Structure_Options
    use Galacticus_Error
    use Root_Finder
    use FGSL
    use Galactic_Structure_Enclosed_Masses
    use Numerical_Constants_Math
    use Dark_Matter_Halo_Scales
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),          intent(inout), pointer  :: thisNode, hostNode
    type(fgsl_function),     save                    :: rootFunction
    type(fgsl_root_fsolver), save                    :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision                                 :: radiusMinimum,radiusMaximum
    type(c_ptr)                                      :: parameterPointer

    ! Assume zero separation by default.
    Black_Hole_Binary_Initial_Radius_Tidal_Radius=0.0d0

    ! If the primary black hole has zero mass (i.e. has been ejected), then return immediately.
    if (Tree_Node_Black_Hole_Mass(thisNode) <= 0.0d0) return

    ! Get the half-mass radius of the satellite galaxy.
    radiusHalfMass=Galactic_Structure_Radius_Enclosing_Mass(thisNode,fractionalMass=0.5d0,massType=massTypeGalactic)

    ! Get the mass within the half-mass radius.
    massHalf=Galactic_Structure_Enclosed_Mass(thisNode,radiusHalfMass,massType=massTypeGalactic)

    ! Solve for the radius around the host at which the satellite gets disrupted.
    activeNode => hostNode
    radiusMinimum=Galactic_Structure_Radius_Enclosing_Mass(hostNode,fractionalMass=0.5d0,massType=massTypeGalactic)
    do while (Tidal_Radius_Root(radiusMinimum,parameterPointer) <= 0.0d0)
       radiusMinimum=0.5d0*radiusMinimum
    end do
    radiusMaximum=Galactic_Structure_Radius_Enclosing_Mass(hostNode,fractionalMass=0.5d0,massType=massTypeGalactic)
    do while (Tidal_Radius_Root(radiusMaximum,parameterPointer) >= 0.0d0)
       radiusMaximum=2.0d0*radiusMaximum
    end do
    Black_Hole_Binary_Initial_Radius_Tidal_Radius=Root_Find(radiusMinimum,radiusMaximum,Tidal_Radius_Root,parameterPointer&
         &,rootFunction,rootFunctionSolver,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    return
  end function Black_Hole_Binary_Initial_Radius_Tidal_Radius

  function Tidal_Radius_Root(radius,parameterPointer) bind(c)
    !% Root function used in solving for the radius of tidal disruption of a satellite galaxy.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: Tidal_Radius_Root

    ! Evaluate the root function.
    Tidal_Radius_Root= Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeGalactic)/radius        **3 &
       &              -massHalf                                                                     /radiusHalfMass**3
    return
  end function Tidal_Radius_Root

end module Black_Hole_Binary_Initial_Radii_Tidal_Radius
