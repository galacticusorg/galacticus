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

!% Contains a module which implements a black hole binary recoil velocity which follows the formulae in
!% \cite{campanelli_large_2007}, derived from post-Newtonian evaluations.

module Black_Hole_Binary_Recoil_Velocities_Standard
  !% Implements a black hole binary recoil velocity which follows the formulae in \cite{campanelli_large_2007}, derived from
  !% post-Newtonian evaluations.
  use FGSL
  implicit none
  private
  public :: Black_Hole_Binary_Recoil_Velocity_Standard_Initialize, Black_Hole_Binary_Recoil_Velocity_Standard_Snapshot,&
       & Black_Hole_Binary_Recoil_Velocity_Standard_State_Store, Black_Hole_Binary_Recoil_Velocity_Standard_State_Retrieve

  ! Random sequence objects.
  type(fgsl_rng) :: pseudoSequenceObject,clonedPseudoSequenceObject
  logical        :: pseudoSequenceReset=.true.,pseudoSequenceResetSnapshot
  !$omp threadprivate(pseudoSequenceReset,pseudoSequenceObject,pseudoSequenceResetSnapshot,clonedPseudoSequenceObject)

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

    if (blackHoleBinaryRecoilVelocityMethod == 'Campanelli2007') Black_Hole_Binary_Recoil_Velocity_Get => Black_Hole_Binary_Recoil_Velocity_Standard
    return
  end subroutine Black_Hole_Binary_Recoil_Velocity_Standard_Initialize

  double precision function Black_Hole_Binary_Recoil_Velocity_Standard(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)
    !% Returns the recoil velocity of a black hole binary, accounting for the parallel and perpendicular velocities, plus that of
    !% the binary's center of mass. Constants used are retrieved from the articles by: \cite{koppitz_recoil_2007} for $H=(7.3\pm
    !% 0.3)10^3$~km/s, \cite{gonzalez_maximum_2007} for $A=1.2 \times 10^4$~km/s $B=-0.93$, \cite{gonzalez_supermassive_2007} for $K
    !% \cos(\delta\theta)=(6,-5.3)10^4$~km/s and $K=(6.0\pm 0.1)10^4$~km/s.
    use Tree_Nodes
    use Pseudo_Random
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in):: massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2
    double precision, parameter :: H=7.3d3, K=6.0d4, A=1.2d4, B=-0.93
    double precision            :: velocityOrthogonal,velocityParallel,velocityMass,q    
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

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Black_Hole_Binary_Recoil_Velocity_Standard_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Black_Hole_Binary_Recoil_Velocity_Standard_Snapshot
    !% Store a snapshot of the random number generator internal state.
    implicit none

    if (.not.pseudoSequenceReset) clonedPseudoSequenceObject=FGSL_Rng_Clone(pseudoSequenceObject)
    pseudoSequenceResetSnapshot=pseudoSequenceReset
    return
  end subroutine Black_Hole_Binary_Recoil_Velocity_Standard_Snapshot
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Black_Hole_Binary_Recoil_Velocity_Standard_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Black_Hole_Binary_Recoil_Velocity_Standard_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) pseudoSequenceResetSnapshot
    if (.not.pseudoSequenceResetSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Black_Hole_Binary_Recoil_Velocity_Standard_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Black_Hole_Binary_Recoil_Velocity_Standard_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Black_Hole_Binary_Recoil_Velocity_Standard_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) pseudoSequenceReset
    if (.not.pseudoSequenceReset) call Pseudo_Random_Retrieve(pseudoSequenceObject,fgslStateFile)
    return
  end subroutine Black_Hole_Binary_Recoil_Velocity_Standard_State_Retrieve
    
end module Black_Hole_Binary_Recoil_Velocities_Standard
