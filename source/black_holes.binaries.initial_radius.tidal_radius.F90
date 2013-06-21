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

!% Contains a module which implements a black hole binary initial separation based on tidal disruption of the satellite galaxy.

module Black_Hole_Binary_Initial_Radii_Tidal_Radius
  !% Implements a black hole binary initial separation based on tidal disruption of the satellite galaxy.
  use Galacticus_Nodes
  implicit none
  private
  public :: Black_Hole_Binary_Initial_Radii_Tidal_Radius_Initialize

  ! Variables used in root finding.
  double precision                    :: massHalf  , radiusHalfMass
  type            (treeNode), pointer :: activeNode
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
    type     (varying_string                               ), intent(in   )          :: blackHoleBinaryInitialRadiiMethod
    procedure(Black_Hole_Binary_Initial_Radius_Tidal_Radius), intent(inout), pointer :: Black_Hole_Binary_Initial_Radius_Get

    if (blackHoleBinaryInitialRadiiMethod == 'tidalRadius') Black_Hole_Binary_Initial_Radius_Get => Black_Hole_Binary_Initial_Radius_Tidal_Radius
    return
  end subroutine Black_Hole_Binary_Initial_Radii_Tidal_Radius_Initialize

  double precision function Black_Hole_Binary_Initial_Radius_Tidal_Radius(thisNode,hostNode)
    !% Returns an initial separation for a binary black holes through tidal disruption.
    use ISO_Varying_String
    use Galactic_Structure_Options
    use Galacticus_Error
    use Root_Finder
    use Galactic_Structure_Enclosed_Masses
    use Numerical_Constants_Math
    use Dark_Matter_Halo_Scales
    use Galacticus_Display
    use String_Handling
    implicit none
    type (treeNode              ), intent(inout), pointer :: hostNode              , thisNode
    class(nodeComponentBlackHole)               , pointer :: thisBlackHoleComponent
    type (rootFinder            ), save                   :: finder
    !$omp threadprivate(finder)
    ! Assume zero separation by default.
    Black_Hole_Binary_Initial_Radius_Tidal_Radius=0.0d0
    ! Get the black hole component.
    thisBlackHoleComponent => thisNode%blackHole(instance=1)
    ! If the primary black hole has zero mass (i.e. has been ejected), then return immediately.
    if (thisBlackHoleComponent%mass() <= 0.0d0) return
    ! Get the half-mass radius of the satellite galaxy.
    radiusHalfMass=Galactic_Structure_Radius_Enclosing_Mass(thisNode,fractionalMass=0.5d0,massType=massTypeGalactic)
    ! Get the mass within the half-mass radius.
    massHalf      =Galactic_Structure_Enclosed_Mass        (thisNode,radiusHalfMass      ,massType=massTypeGalactic)
    ! Return zero radius for massless galaxy.
    if (radiusHalfMass <= 0.0d0 .or. massHalf <= 0.0d0) return
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(Tidal_Radius_Root                               )
       call finder%rangeExpand (                                                             &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
       call finder%tolerance   (toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    end if
    ! Solve for the radius around the host at which the satellite gets disrupted.
    activeNode => hostNode
    Black_Hole_Binary_Initial_Radius_Tidal_Radius=finder%find                              &
         &  (                                                                              &
         &   rootGuess=Galactic_Structure_Radius_Enclosing_Mass(                           &
         &                                                      hostNode                 , &
         &                                                      fractionalMass=0.5d0     , &
         &                                                      massType=massTypeGalactic  &
         &                                                     )                           &
         &  )
    return
  end function Black_Hole_Binary_Initial_Radius_Tidal_Radius

  double precision function Tidal_Radius_Root(radius)
    !% Root function used in solving for the radius of tidal disruption of a satellite galaxy.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    double precision, intent(in   ) :: radius

    ! Evaluate the root function.
    Tidal_Radius_Root= Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeGalactic)/massHalf &
       &              -(radius/radiusHalfMass)**3
    return
  end function Tidal_Radius_Root

end module Black_Hole_Binary_Initial_Radii_Tidal_Radius
