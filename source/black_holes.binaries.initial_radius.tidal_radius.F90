!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements a black hole binary initial separation based on tidal disruption of the satellite galaxy.

module Black_Hole_Binary_Initial_Radii_Tidal_Radius
  !% Implements a black hole binary initial separation based on tidal disruption of the satellite galaxy.
  use Galacticus_Nodes
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
    type     (varying_string  ),          intent(in   ) :: blackHoleBinaryInitialRadiiMethod
    procedure(double precision), pointer, intent(inout) :: Black_Hole_Binary_Initial_Radius_Get
    
    if (blackHoleBinaryInitialRadiiMethod == 'tidalRadius') Black_Hole_Binary_Initial_Radius_Get => Black_Hole_Binary_Initial_Radius_Tidal_Radius
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
    type            (treeNode              ), intent(inout), pointer :: thisNode, hostNode
    class           (nodeComponentBlackHole),                pointer :: thisBlackHoleComponent
    type            (fgsl_function         ), save                   :: rootFunction
    type            (fgsl_root_fsolver     ), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision                                                 :: radiusMinimum,radiusMaximum
    type            (c_ptr                 )                         :: parameterPointer

    ! Assume zero separation by default.
    Black_Hole_Binary_Initial_Radius_Tidal_Radius=0.0d0
    ! Get the black hole component.
    thisBlackHoleComponent => thisNode%blackHole(instance=1)
    ! If the primary black hole has zero mass (i.e. has been ejected), then return immediately.
    if (thisBlackHoleComponent%mass() <= 0.0d0) return
    ! Get the half-mass radius of the satellite galaxy.
    radiusHalfMass=Galactic_Structure_Radius_Enclosing_Mass(thisNode,fractionalMass=0.5d0,massType=massTypeGalactic)
    ! Get the mass within the half-mass radius.
    massHalf=Galactic_Structure_Enclosed_Mass(thisNode,radiusHalfMass,massType=massTypeGalactic)
    ! Return zero radius for massless galaxy.
    if (radiusHalfMass <= 0.0d0 .or. massHalf <= 0.0d0) return
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
    !$omp critical (Black_Hole_Binary_Initial_Radius_Tidal_Radius_Root)
    Black_Hole_Binary_Initial_Radius_Tidal_Radius=Root_Find(radiusMinimum,radiusMaximum,Tidal_Radius_Root,parameterPointer&
         &,rootFunction,rootFunctionSolver,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    !$omp end critical (Black_Hole_Binary_Initial_Radius_Tidal_Radius_Root)
    return
  end function Black_Hole_Binary_Initial_Radius_Tidal_Radius

  function Tidal_Radius_Root(radius,parameterPointer) bind(c)
    !% Root function used in solving for the radius of tidal disruption of a satellite galaxy.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double), value :: radius
    type(c_ptr   ), value :: parameterPointer
    real(c_double)        :: Tidal_Radius_Root

    ! Evaluate the root function.
    Tidal_Radius_Root= Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeGalactic)/massHalf &
       &              -(radius/radiusHalfMass)**3
    return
  end function Tidal_Radius_Root

end module Black_Hole_Binary_Initial_Radii_Tidal_Radius
