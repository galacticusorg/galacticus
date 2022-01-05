!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!!{
Contains a module which implements galactic structure tasks for the standard black hole node component.
!!}

module Node_Component_Black_Hole_Standard_Structure_Tasks
  !!{
  Implements galactic structure tasks for the standard black hole tree node component.
  !!}
  implicit none
  private
  public :: Node_Component_Black_Hole_Standard_Rotation_Curve         , Node_Component_Black_Hole_Standard_Potential, &
       &    Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient

contains

  !![
  <rotationCurveTask>
   <unitName>Node_Component_Black_Hole_Standard_Rotation_Curve</unitName>
  </rotationCurveTask>
  !!]
  double precision function Node_Component_Black_Hole_Standard_Rotation_Curve(node,radius,componentType,massType)
    !!{
    Computes the rotation curve for the central black hole. Assumes a point mass black hole with a Keplerian rotation curve,
    \emph{except} that the rotation speed is limited to never exceed the speed of light.
    !!}
    use :: Black_Hole_Fundamentals         , only : Black_Hole_Gravitational_Radius
    use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole         , nodeComponentBlackHoleStandard, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Prefixes    , only : milli
    implicit none
    type            (treeNode              ), intent(inout)           :: node
    integer                                 , intent(in   )           :: componentType, massType
    double precision                        , intent(in   )           :: radius
    class           (nodeComponentBlackHole)               , pointer  :: blackHole
    double precision                                                  :: componentMass

    ! Set to zero by default.
    Node_Component_Black_Hole_Standard_Rotation_Curve=0.0d0
    ! Get the black hole component and check that it is of the standard class.
    blackHole => node%blackHole(instance=1)
    select type (blackHole)
    class is (nodeComponentBlackHoleStandard)
       ! Check if the radius exceeds the gravitational radius.
       if (radius > Black_Hole_Gravitational_Radius(blackHole)) then
          ! Radius is larger than the gravitational radius - compute the rotation speed.
          componentMass=blackHole%enclosedMass(radius,componentType,massType,weightByMass &
               &,weightIndexNull)
          if (componentMass > 0.0d0) Node_Component_Black_Hole_Standard_Rotation_Curve=sqrt(gravitationalConstantGalacticus&
               &*componentMass/radius)
       else
          ! Radius is less than the gravitational radius - return the speed of light.
          Node_Component_Black_Hole_Standard_Rotation_Curve=speedLight*milli
       end if
    end select
    return
  end function Node_Component_Black_Hole_Standard_Rotation_Curve

  !![
  <potentialTask>
   <unitName>Node_Component_Black_Hole_Standard_Potential</unitName>
  </potentialTask>
  !!]
  double precision function Node_Component_Black_Hole_Standard_Potential(node,radius,componentType,massType,status)
    !!{
    Compute the gravitational potential due to a black hole.
    !!}
    use :: Black_Hole_Fundamentals         , only : Black_Hole_Gravitational_Radius
    use :: Galactic_Structure_Options      , only : componentTypeAll               , componentTypeBlackHole        , massTypeAll, massTypeBlackHole, &
          &                                         weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole         , nodeComponentBlackHoleStandard, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode              ), intent(inout)           :: node
    integer                                 , intent(in   )           :: componentType, massType
    double precision                        , intent(in   )           :: radius
    integer                                 , intent(inout), optional :: status
    class           (nodeComponentBlackHole)               , pointer  :: blackHole
    double precision                                                  :: componentMass
    !$GLC attributes unused :: status

    Node_Component_Black_Hole_Standard_Potential=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    ! Get the black hole component and check that it is of the standard class.
    blackHole => node%blackHole(instance=1)
    select type (blackHole)
    class is (nodeComponentBlackHoleStandard)
       if (Black_Hole_Gravitational_Radius(blackHole) <=0.0d0) return
       ! Computes the potential - limit the radius to the gravitational radius to avoid divergent potentials.
       componentMass=blackHole%enclosedMass(radius,componentType,massType,weightByMass&
            &,weightIndexNull)
       Node_Component_Black_Hole_Standard_Potential=-gravitationalConstantGalacticus*componentMass/max(radius &
            &,Black_Hole_Gravitational_Radius(blackHole))
    end select
    return
  end function Node_Component_Black_Hole_Standard_Potential

  !![
  <rotationCurveGradientTask>
   <unitName>Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient</unitName>
  </rotationCurveGradientTask>
  !!]
  double precision function Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient(node,radius,componentType &
       &,massType)
    !!{
    Computes the rotation curve gradient for the central black hole. Assumes a point mass black hole with a Keplerian
    rotation curve, \emph{except} that the rotation speed is limited to never exceed the speed of light.
    !!}
    use :: Black_Hole_Fundamentals         , only : Black_Hole_Gravitational_Radius
    use :: Galactic_Structure_Options      , only : componentTypeAll               , componentTypeBlackHole        , massTypeAll, massTypeBlackHole, &
          &                                         weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole         , nodeComponentBlackHoleStandard, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode              ), intent(inout)           :: node
    integer                                 , intent(in   )           :: componentType, massType
    double precision                        , intent(in   )           :: radius
    class           (nodeComponentBlackHole)               , pointer  :: blackHole
    double precision                                                  :: componentMass

    ! Set to zero by default.
    Node_Component_Black_Hole_Standard_Rotation_Curve_Gradient=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (      radius        <=            0.0d0                                              ) return
    ! Get the black hole component and check that it is of the standard class.
    blackHole => node%blackHole(instance=1)
    select type (blackHole)
    class is (nodeComponentBlackHoleStandard)
       componentMass=blackHole%enclosedMass(radius,componentType,massType,weightByMass,weightIndexNull)
       if (componentMass == 0.0d0) return
       if (radius > Black_Hole_Gravitational_Radius(blackHole)) then
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
