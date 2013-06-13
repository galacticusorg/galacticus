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

!% Contains a module which implements galactic structure tasks for the standard black hole node component.

module Node_Component_Black_Hole_Standard_Structure_Tasks
  !% Implements galactic structure tasks for the standard black hole tree node component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Black_Hole_Standard_Rotation_Curve         , &
       &    Node_Component_Black_Hole_Standard_Potential    , Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient
  
contains

  !# <rotationCurveTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Rotation_Curve</unitName>
  !# </rotationCurveTask>
  double precision function Node_Component_Black_Hole_Standard_Rotation_Curve(thisNode,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve for the central black hole. Assumes a point mass black hole with a Keplerian rotation curve,
    !% \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer  :: thisNode
    integer                                 , intent(in   )           :: componentType,massType
    double precision                        , intent(in   )           :: radius
    logical                                 , intent(in   ), optional :: haloLoaded
    class           (nodeComponentBlackHole),                pointer  :: thisBlackHoleComponent
    double precision                                                  :: componentMass

    ! Set to zero by default.
    Node_Component_Black_Hole_Standard_Rotation_Curve=0.0d0
    ! Get the black hole component and check that it is of the standard class.
    thisBlackHoleComponent => thisNode%blackHole(instance=1)
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleStandard)
       ! Check if the radius exceeds the gravitational radius.
       if (radius > Black_Hole_Gravitational_Radius(thisBlackHoleComponent)) then
          ! Radius is larger than the gravitational radius - compute the rotation speed.
          componentMass=thisBlackHoleComponent%enclosedMass(radius,componentType,massType,weightByMass &
               &,weightIndexNull,haloLoaded)
          if (componentMass > 0.0d0) Node_Component_Black_Hole_Standard_Rotation_Curve=sqrt(gravitationalConstantGalacticus&
               &*componentMass/radius)
       else 
          ! Radius is less than the gravitational radius - return the speed of light.
          Node_Component_Black_Hole_Standard_Rotation_Curve=speedLight*milli
       end if
    end select
    return
  end function Node_Component_Black_Hole_Standard_Rotation_Curve

  !# <potentialTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Potential</unitName>
  !# </potentialTask>
  double precision function Node_Component_Black_Hole_Standard_Potential(thisNode,radius,componentType,massType,haloLoaded)
    !% Compute the gravitational potential due to a black hole.
    use Galacticus_Nodes
    use Numerical_Constants_Physical 
    use Galactic_Structure_Options
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer  :: thisNode
    integer                                 , intent(in   )           :: componentType,massType
    double precision                        , intent(in   )           :: radius
    logical                                 , intent(in   ), optional :: haloLoaded
    class           (nodeComponentBlackHole),                pointer  :: thisBlackHoleComponent
    double precision                                                  :: componentMass

    Node_Component_Black_Hole_Standard_Potential=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    ! Get the black hole component and check that it is of the standard class.
    thisBlackHoleComponent => thisNode%blackHole(instance=1)
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleStandard)
       if (Black_Hole_Gravitational_Radius(thisBlackHoleComponent) <=0.0d0) return
       ! Computes the potential - limit the radius to the gravitational radius to avoid divergent potentials.
       componentMass=thisBlackHoleComponent%enclosedMass(radius,componentType,massType,weightByMass&
            &,weightIndexNull,haloLoaded)
       Node_Component_Black_Hole_Standard_Potential=-gravitationalConstantGalacticus*componentMass/max(radius &
            &,Black_Hole_Gravitational_Radius(thisBlackHoleComponent))
    end select
    return
  end function Node_Component_Black_Hole_Standard_Potential

  !# <rotationCurveGradientTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient</unitName>
  !# </rotationCurveGradientTask>
  double precision function Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient(thisNode,radius,componentType &
       &,massType,haloLoaded)
    !% Computes the rotation curve gradient for the central black hole. Assumes a point mass black hole with a Keplerian 
    !% rotation curve, \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer  :: thisNode
    integer                                 , intent(in   )           :: componentType,massType
    double precision                        , intent(in   )           :: radius
    logical                                 , intent(in   ), optional :: haloLoaded
    class           (nodeComponentBlackHole),                pointer  :: thisBlackHoleComponent
    double precision                                                  :: componentMass

    ! Set to zero by default.
    Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (      radius        <=            0.0d0                                              ) return
    ! Get the black hole component and check that it is of the standard class.
    thisBlackHoleComponent => thisNode%blackHole(instance=1)
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleStandard)
       componentMass=thisBlackHoleComponent%enclosedMass(radius,componentType,massType,weightByMass,weightIndexNull&
            &,haloLoaded)
       if (componentMass == 0.0d0) return       
       if (radius > Black_Hole_Gravitational_Radius(thisBlackHoleComponent)) then
          Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient=     &
               &                         -gravitationalConstantGalacticus &
               &                         *componentMass                   &
               &                         /radius**2
       else
          Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient=0.0d0
       end if
    end select
    return
  end function Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient
  
end module Node_Component_Black_Hole_Standard_Structure_Tasks
