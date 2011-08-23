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

!% Contains a module which implements a black hole binary recoil velocity which follows the formulae in
!% \cite{campanelli_large_2007}, derived from post-Newtonian evaluations.

module Black_Hole_Binary_Recoil_Velocities_Standard
  !% Implements a black hole binary recoil velocity which follows the formulae in \cite{campanelli_large_2007}, derived from
  !% post-Newtonian evaluations.
  implicit none
  private
  public :: Black_Hole_Binary_Recoil_Velocity_Standard_Initialize

contains

  !# <blackHoleBinaryRecoilVelocityMethod>
  !#  <unitName>Black_Hole_Binary_Recoil_Velocity_Standard_Initialize</unitName>
  !# </blackHoleBinaryRecoilVelocityMethod>
  subroutine Black_Hole_Binary_Recoil_Velocity_Standard_Initialize(blackHoleBinaryRecoilVelocityMethod,Black_Hole_Binary_Recoil_Velocity_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: blackHoleBinaryRecoilVelocityMethod
    procedure(double precision), pointer, intent(inout) :: Black_Hole_Binary_Recoil_Velocity_Get

    if (blackHoleBinaryRecoilVelocityMethod == 'Campanelli 2007') Black_Hole_Binary_Recoil_Velocity_Get => Black_Hole_Binary_Recoil_Velocity_Standard
    return
  end subroutine Black_Hole_Binary_Recoil_Velocity_Standard_Initialize

  double precision function Black_Hole_Binary_Recoil_Velocity_Standard(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)
    !% Returns the recoil velocity of a black hole binary, accounting for the parallel and perpendicular velocities, plus that of
    !% the binary's center of mass. Constants used are retrieved from the articles by: \cite{koppitz_recoil_2007} for $H=(7.3\pm
    !% 0.3)10^3$~km/s, \cite{gonzalez_maximum_2007} for $A=1.2 \times 10^4$~km/s $B=-0.93$, \cite{gonzalez_supermassive_2007} for $K
    !% \cos(\delta\theta)=(6,-5.3)10^4$~km/s and $K=(6.0\pm 0.1)10^4$~km/s.
    use FGSL
    use Tree_Nodes
    use Pseudo_Random
    use Numerical_Constants_Math
    implicit none
    type(fgsl_rng),   save      :: pseudoSequenceObject
    logical,          save      :: pseudoSequenceReset=.true.
    double precision, intent(in):: massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2
    double precision, parameter :: H=7.3d3, K=6.0d4, A=1.2d4, B=-0.93
    double precision            :: velocityRecoil,velocityOrthogonal,velocityParallel,velocityMass,q    
    double precision            :: alpha1Orthogonal,alpha2Orthogonal,alpha1Parallel,alpha2Parallel
    double precision            :: theta, phi

    ! If either black hole has non-positive mass, recoil velocity must be zero.
    if (massBlackHole1 <= 0.0d0 .or. massBlackHole2 <= 0.0d0) then
       Black_Hole_Binary_Recoil_Velocity_Standard=0.0d0
       return
    end if

    ! Select angular coordinates at random for the spin of the black hole.
    phi  =     2.0d0*Pi*Pseudo_Random_Get(pseudoSequenceObject,reset=pseudoSequenceReset)
    theta=acos(2.0d0*   Pseudo_Random_Get(pseudoSequenceObject,reset=pseudoSequenceReset)-1.0d0)

    ! Compute the mass ratio of the two black holes.
    q=massBlackHole1/massBlackHole2

    ! Compute alpha (and components), the angular momentum of the black hole per unit mass. This is equal to the spin scalar
    ! (since we don't track the direction of spin).
    alpha1Orthogonal=0.0d0
    alpha1Parallel  =spinBlackHole1
    alpha2Orthogonal=spinBlackHole2*sin(theta)*cos(phi)
    alpha2Parallel  =spinBlackHole2*cos(theta)

    ! Compute velocity kicks in each direction.
    velocityOrthogonal=H*(q**2      /(1.0d0+q)**5)*(alpha2Parallel  -q*alpha1Parallel  )
    velocityParallel  =K*(q**2      /(1.0d0+q)**5)*(alpha2Orthogonal-q*alpha1Orthogonal)
    velocityMass      =A*(q**2*(1-q)/(1.0d0+q)**5)*(1.0d0+B*(q/(1.0d0+q)**2))

    ! Compute the net recoil velocity.
    Black_Hole_Binary_Recoil_Velocity_Standard=dsqrt(velocityParallel**2+velocityMass**2+velocityOrthogonal**2)     
    return
  end function Black_Hole_Binary_Recoil_Velocity_Standard

end module Black_Hole_Binary_Recoil_Velocities_Standard
