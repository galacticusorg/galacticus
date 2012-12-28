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
  public :: Node_Component_Black_Hole_Standard_Enclosed_Mass, Node_Component_Black_Hole_Standard_Rotation_Curve         , &
       &    Node_Component_Black_Hole_Standard_Potential    , Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient
  
contains

  !# <enclosedMassTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Enclosed_Mass</unitName>
  !# </enclosedMassTask>
  subroutine Node_Component_Black_Hole_Standard_Enclosed_Mass(thisNode,radius,massType,componentType,weightBy,weightIndex&
       &,componentMass)
    !% Computes the mass within a given radius for a central black hole. Black hole is treated as a point mass.
    use Galactic_Structure_Options
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(in   )          :: massType,componentType,weightBy,weightIndex
    double precision                        , intent(in   )          :: radius
    class           (nodeComponentBlackHole),                pointer :: thisBlackHoleComponent
    double precision                        , intent(  out)          :: componentMass

    ! Set zero enclosed mass by default.
    componentMass=0.0d0
    ! Return the black hole mass only if massType and componentType are of black hole type.
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    ! Get the black hole component and check that it is of the standard class.
    thisBlackHoleComponent => thisNode%blackHole(instance=1)
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleStandard)
       if (radius >= 0.0d0) then
          select case (weightBy)
          case (weightByMass)
             componentMass=thisBlackHoleComponent%mass()
          end select
       end if
    end select
    return
  end subroutine Node_Component_Black_Hole_Standard_Enclosed_Mass

  !# <rotationCurveTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Rotation_Curve</unitName>
  !# </rotationCurveTask>
  subroutine Node_Component_Black_Hole_Standard_Rotation_Curve(thisNode,radius,massType,componentType,componentVelocity)
    !% Computes the rotation curve for the central black hole. Assumes a point mass black hole with a Keplerian rotation curve,
    !% \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(in   )          :: massType,componentType
    double precision                        , intent(in   )          :: radius
    double precision                        , intent(  out)          :: componentVelocity
    class           (nodeComponentBlackHole),                pointer :: thisBlackHoleComponent
    double precision                                                 :: componentMass

    ! Set to zero by default.
    componentVelocity=0.0d0
    ! Get the black hole component and check that it is of the standard class.
    thisBlackHoleComponent => thisNode%blackHole(instance=1)
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleStandard)
       ! Check if the radius exceeds the gravitational radius.
       if (radius > Black_Hole_Gravitational_Radius(thisBlackHoleComponent)) then
          ! Radius is larger than the gravitational radius - compute the rotation speed.
          call Node_Component_Black_Hole_Standard_Enclosed_Mass(thisNode,radius,massType,componentType,weightByMass&
               &,weightIndexNull,componentMass)
          if (componentMass > 0.0d0) componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass/radius)
       else 
          ! Radius is less than the gravitational radius - return the speed of light.
          componentVelocity=speedLight*milli
       end if
    end select
    return
  end subroutine Node_Component_Black_Hole_Standard_Rotation_Curve

  !# <potentialTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Potential</unitName>
  !# </potentialTask>
  subroutine Node_Component_Black_Hole_Standard_Potential(thisNode,radius,componentType,massType,componentPotential)
    !% Compute the gravitational potential due to a black hole.
    use Galacticus_Nodes
    use Numerical_Constants_Physical 
    use Galactic_Structure_Options
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(in   )          :: componentType,massType
    double precision                        , intent(in   )          :: radius
    double precision                        , intent(  out)          :: componentPotential
    class           (nodeComponentBlackHole),                pointer :: thisBlackHoleComponent
    double precision                                                 :: componentMass

    componentPotential=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    ! Get the black hole component and check that it is of the standard class.
    thisBlackHoleComponent => thisNode%blackHole(instance=1)
    if (Black_Hole_Gravitational_Radius(thisBlackHoleComponent) <=0.0d0) return
    ! Computes the potential - limit the radius to the gravitational radius to avoid divergent potentials.
    call Node_Component_Black_Hole_Standard_Enclosed_Mass(thisNode,radius,massType,componentType,weightByMass,weightIndexNull &
         &,componentMass)
    componentPotential=-gravitationalConstantGalacticus*componentMass/max(radius&
         &,Black_Hole_Gravitational_Radius(thisBlackHoleComponent))
    return
  end subroutine Node_Component_Black_Hole_Standard_Potential

  !# <rotationCurveGradientTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient</unitName>
  !# </rotationCurveGradientTask>
  subroutine Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient(thisNode,radius,massType,componentType&
       &,componentRotationCurveGradient)
    !% Computes the rotation curve gradient for the central black hole. Assumes a point mass black hole with a Keplerian 
    !% rotation curve, \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(in   )          :: massType,componentType
    double precision                        , intent(in   )          :: radius
    double precision                        , intent(  out)          :: componentRotationCurveGradient
    class           (nodeComponentBlackHole),                pointer :: thisBlackHoleComponent
    double precision                                                 :: componentMass

    ! Set to zero by default.
    componentRotationCurveGradient=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (      radius        <=            0.0d0                                              ) return
    call Node_Component_Black_Hole_Standard_Enclosed_Mass(thisNode,radius,massType,componentType,weightByMass,weightIndexNull&
         &,componentMass)
    if (componentMass == 0.0d0) return
    if (radius > Black_Hole_Gravitational_Radius(thisBlackHoleComponent)) then
       componentRotationCurveGradient=-gravitationalConstantGalacticus &
            &                         *componentMass                   &
            &                         /radius**2
    else
       componentRotationCurveGradient=0.0d0
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient
  
end module Node_Component_Black_Hole_Standard_Structure_Tasks
