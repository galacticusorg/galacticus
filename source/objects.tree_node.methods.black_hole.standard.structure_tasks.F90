!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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

!% Contains a module which implements galactic structure tasks for the standard black hole tree node component.

module Tree_Node_Methods_Black_Hole_Structure_Tasks
  !% Implements galactic structure tasks for the standard black hole tree node component.
  use Tree_Nodes
  use Tree_Node_Methods_Black_Hole_Data
  use Components
  implicit none
  private
  public :: Black_Hole_Rotation_Curve_Standard,Black_Hole_Enclosed_Mass_Standard,Black_Hole_Potential_Standard,Black_Hole_Rotation_Curve_Gradient_Standard
  
contains

  !# <enclosedMassTask>
  !#  <unitName>Black_Hole_Enclosed_Mass_Standard</unitName>
  !# </enclosedMassTask>
  subroutine Black_Hole_Enclosed_Mass_Standard(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for a central black hole. Black hole is treated as a point mass.
    use Galactic_Structure_Options
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType,weightBy,weightIndex
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass

    ! Set zero enclosed mass by default.
    componentMass=0.0d0

    ! Return the black hole mass only if massType and componentType are of black hole type.
    if (.not.methodSelected                                                                  ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (.not.thisNode%componentExists(componentIndex)                                        ) return

    ! Compute the enclosed mass.
    if (radius >= 0.0d0) then
       select case (weightBy)
       case (weightByMass)
          componentMass=Tree_Node_Black_Hole_Mass(thisNode,instance=1)
       end select
    end if
    return
  end subroutine Black_Hole_Enclosed_Mass_Standard

  !# <rotationCurveTask>
  !#  <unitName>Black_Hole_Rotation_Curve_Standard</unitName>
  !# </rotationCurveTask>
  subroutine Black_Hole_Rotation_Curve_Standard(thisNode,radius,massType,componentType,componentVelocity)
    !% Computes the rotation curve for the central black hole. Assumes a point mass black hole with a Keplerian rotation curve,
    !% \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Tree_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes

    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentVelocity
    double precision                         :: componentMass

    ! Set to zero by default.
    componentVelocity=0.0d0

    ! Compute if a black hole is present.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then       
       ! Check if the radius exceeds the gravitational radius. Do the calculation here rather than calling the gravitational
       ! radius function since we want to enforce the calculation to always use the first black hole instance. <v0.9.1> This ugly
       ! solution should be solved by passing a black hole object directly to the gravitational radius function.
       if (radius > gravitationalConstantGalacticus*Tree_Node_Black_Hole_Mass(thisNode,instance=1)/(milli*speedLight)**2) then
          ! Radius is larger than the gravitational radius - compute the rotation speed.
          call Black_Hole_Enclosed_Mass_Standard(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
          if (componentMass > 0.0d0) componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass/radius)
       else 
          ! Radius is less than the gravitational radius - return the speed of light.
          componentVelocity=speedLight*milli
       end if
    end if
    return
  end subroutine Black_Hole_Rotation_Curve_Standard

  !# <potentialTask>
  !#  <unitName>Black_Hole_Potential_Standard</unitName>
  !# </potentialTask>
  subroutine Black_Hole_Potential_Standard(thisNode,radius,componentType,massType,componentPotential)
    !% Compute the gravitational potential due to a black hole.
    use Tree_Nodes
    use Numerical_Constants_Physical 
    use Galactic_Structure_Options
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: componentType,massType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentPotential
    double precision                         :: componentMass

    componentPotential=0.0d0
    if (.not.methodSelected                                                                  ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (.not.thisNode%componentExists(componentIndex)                                        ) return
    if (Black_Hole_Gravitational_Radius(thisNode) <=0.0d0) return

    ! Computes the potential - limit the radius to the gravitational radius to avoid divergent potentials.
    call Black_Hole_Enclosed_Mass_Standard(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
    componentPotential=-gravitationalConstantGalacticus                       &
          &            *componentMass                                         &
          &            /max(radius,Black_Hole_Gravitational_Radius(thisNode))
    return
  end subroutine Black_Hole_Potential_Standard

  !# <rotationCurveGradientTask>
  !#  <unitName>Black_Hole_Rotation_Curve_Gradient_Standard</unitName>
  !# </rotationCurveGradientTask>
  subroutine Black_Hole_Rotation_Curve_Gradient_Standard(thisNode,radius,massType,componentType,componentRotationCurveGradient)
    !% Computes the rotation curve gradient for the central black hole. Assumes a point mass black hole with a Keplerian 
    !% rotation curve, \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Tree_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentRotationCurveGradient
    double precision                         :: componentMass

    ! Set to zero by default.
    componentRotationCurveGradient=0.0d0
    if (.not.methodSelected                                                                  ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (.not.thisNode%componentExists(componentIndex)                                        ) return
    if (radius <= 0.0d0) return
    call Black_Hole_Enclosed_Mass_Standard(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
    if (componentMass ==0.0d0 ) return
    if (radius > gravitationalConstantGalacticus*componentMass/(milli*speedLight)**2) then
       componentRotationCurveGradient=-gravitationalConstantGalacticus &
            &                         *componentMass                   &
            &                         /radius**2
    else
       componentRotationCurveGradient=0.0d0
    end if
    return
  end subroutine Black_Hole_Rotation_Curve_Gradient_Standard
  
end module Tree_Node_Methods_Black_Hole_Structure_Tasks
