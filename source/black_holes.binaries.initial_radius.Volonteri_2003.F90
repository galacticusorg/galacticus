!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
    
    if (blackHoleBinaryInitialRadiiMethod == 'Volonteri2003') Black_Hole_Binary_Initial_Radius_Get => Black_Hole_Binary_Initial_Radius_Volonteri_2003
    return
  end subroutine Black_Hole_Binary_Initial_Radii_Volonteri_2003_Initialize

  double precision function Black_Hole_Binary_Initial_Radius_Volonteri_2003(thisNode,hostNode)
    !% Returns an initial separation for binary black holes using the method of \cite{volonteri_assembly_2003}, with the assumption that
    !% the local velocity dispersion is approximately the dark matter halo virial velocity.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    implicit none
    type (treeNode              ), intent(inout), pointer :: thisNode,hostNode
    class(nodeComponentBlackHole),                pointer :: thisBlackHoleComponent,hostBlackHoleComponent

    ! Get the black hole components.
    thisBlackHoleComponent => thisNode%blackHole()
    hostBlackHoleComponent => hostNode%blackHole()
    ! Compute the initial separation.
    Black_Hole_Binary_Initial_Radius_Volonteri_2003=gravitationalConstantGalacticus                &
         &                                          *(                                             &
         &                                             thisBlackHoleComponent%mass()               &
         &                                            +hostBlackHoleComponent%mass()               &
         &                                           )                                             &
         &                                          /2.0d0                                         &
         &                                          /Dark_Matter_Halo_Virial_Velocity(hostNode)**2

    return
  end function Black_Hole_Binary_Initial_Radius_Volonteri_2003

end module Black_Hole_Binary_Initial_Radii_Volonteri_2003
