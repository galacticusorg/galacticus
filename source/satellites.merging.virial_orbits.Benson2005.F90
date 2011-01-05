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


!% Contains a module which implements the \cite{benson_orbital_2005} orbital parameter distribution for merging subhalos.

module Virial_Orbits_Benson2005
  !% Implements the \cite{benson_orbital_2005} orbital parameter distribution for merging subhalos.
  use FGSL
  private
  public :: Virial_Orbital_Parameters_Benson2005_Initialize, Virial_Orbital_Parameters_Benson2005_Snapshot,&
       & Virial_Orbital_Parameters_Benson2005_State_Store, Virial_Orbital_Parameters_Benson2005_State_Retrieve

  type(fgsl_rng) :: pseudoSequenceObject,clonedPseudoSequenceObject
  logical        :: resetSequence=.true.,resetSequenceSnapshot
  !$omp threadprivate(pseudoSequenceObject,resetSequence)
  
contains

  !# <virialOrbitsMethod>
  !#  <unitName>Virial_Orbital_Parameters_Benson2005_Initialize</unitName>
  !# </virialOrbitsMethod>
  subroutine Virial_Orbital_Parameters_Benson2005_Initialize(virialOrbitsMethod,Virial_Orbital_Parameters_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: virialOrbitsMethod
    procedure(),          pointer, intent(inout) :: Virial_Orbital_Parameters_Get
    
    if (virialOrbitsMethod.eq.'Benson2005') Virial_Orbital_Parameters_Get => Virial_Orbital_Parameters_Benson2005
    return
  end subroutine Virial_Orbital_Parameters_Benson2005_Initialize

  subroutine Virial_Orbital_Parameters_Benson2005(thisNode,acceptUnboundOrbits,velocityRadial,velocityTangential,angularMomentum&
       &,orbitalEnergy,eccentricity,semimajorAxis)
    !% Return orbital velocities of a satellite selected at random from the fitting function found by \cite{benson_orbital_2005}.
    use Pseudo_Random
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    logical,          intent(in)              :: acceptUnboundOrbits
    double precision, intent(out),   optional :: velocityRadial,velocityTangential,angularMomentum,orbitalEnergy,eccentricity&
         &,semimajorAxis
    type(treeNode),                  pointer  :: hostNode
    double precision, parameter               :: pMax=1.96797d0, velocityMax=3.0d0
    double precision, parameter               :: a(9)=(/0.390052d+01,0.247973d+01,0.102373d+02,0.683922d+00,0.353953d+00&
         &,0.107716d+01 ,0.509837d+00,0.206204d+00,0.314641d+00/)
    double precision, parameter               :: EPS_BOUND=1.0d-4 ! Tolerence to ensure that orbits are sufficiently bound.
    double precision                          :: orbitalA,orbitalB,b1,b2,distributionFunction,uniformRandom&
         &,angularMomentumInternal ,energyInternal,velocityRadialInternal,velocityTangentialInternal,velocityScale,radialScale&
         &,massFactor
    logical                                   :: foundOrbit

    ! Determine the mass parameter.
    hostNode => thisNode%parentNode
    massFactor=1.0d0+Tree_Node_Mass(thisNode)/Tree_Node_Mass(hostNode)
    ! Select an orbit.
    foundOrbit=.false.
    do while(.not.foundOrbit)
       ! Select potential radial and tangential velocities.
       velocityRadialInternal    =Pseudo_Random_Get(pseudoSequenceObject,resetSequence)*velocityMax
       velocityTangentialInternal=Pseudo_Random_Get(pseudoSequenceObject,resetSequence)*velocityMax
       ! Evaluate distribution function for these parameters.
       b1=a(3)*dexp(-a(4)*(velocityTangentialInternal-a(5))**2)
       b2=a(6)*dexp(-a(7)*(velocityTangentialInternal-a(8))**2)
       distributionFunction=a(1)*velocityTangentialInternal*dexp(-a(2)*((velocityTangentialInternal-a(9))**2.0))*dexp(-b1&
            &*(velocityRadialInternal-b2)**2)
       if (distributionFunction > pMax) call Galacticus_Error_Report('Virial_Orbital_Parameters_Benson2005','distribution&
            & function exceeds expected peak value')
       uniformRandom=pMax*Pseudo_Random_Get(pseudoSequenceObject,resetSequence)
       if (uniformRandom <= distributionFunction) then
          foundOrbit=.true.
          ! If requested, check that the orbit is bound. We require it to have E<-EPS_BOUND to ensure that it is sufficiently
          ! bound that later rounding errors will not make it appear unbound.
          if (.not.acceptUnboundOrbits) then
             angularMomentumInternal=velocityTangentialInternal/massFactor
             energyInternal=-1.0d0+0.5d0*(velocityRadialInternal**2+velocityTangentialInternal**2)/massFactor
             foundOrbit=(energyInternal < -EPS_BOUND)
          end if
       end if
    end do
    angularMomentumInternal=velocityTangentialInternal/massFactor
    orbitalA=1.0d0-massFactor/(velocityTangentialInternal**2)
    orbitalB=-velocityRadialInternal/velocityTangentialInternal
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale=Dark_Matter_Halo_Virial_Radius(hostNode)
    if (present(velocityRadial))     velocityRadial    =velocityRadialInternal    *velocityScale
    if (present(velocityTangential)) velocityTangential=velocityTangentialInternal*velocityScale
    if (present(angularMomentum))    angularMomentum   =angularMomentumInternal*velocityScale*radialScale
    if (present(orbitalEnergy  ))    orbitalEnergy     =(-1.0d0+0.5d0*(velocityRadialInternal**2+velocityTangentialInternal**2) &
         &/massFactor)*velocityScale**2
    if (present(eccentricity))       eccentricity      =dsqrt(orbitalA**2+orbitalB**2)*(velocityTangentialInternal**2)/massFactor
    if (present(semimajorAxis))      semimajorAxis     =radialScale*massFactor/(2.0*massFactor-velocityRadialInternal**2 &
         &-velocityTangentialInternal**2)
    return
  end subroutine Virial_Orbital_Parameters_Benson2005

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Virial_Orbital_Parameters_Benson2005_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Virial_Orbital_Parameters_Benson2005_Snapshot
    !% Store a snapshot of the random number generator internal state.
    implicit none

    if (.not.resetSequence) clonedPseudoSequenceObject=FGSL_Rng_Clone(pseudoSequenceObject)
    resetSequenceSnapshot=resetSequence
    return
  end subroutine Virial_Orbital_Parameters_Benson2005_Snapshot
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Virial_Orbital_Parameters_Benson2005_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Virial_Orbital_Parameters_Benson2005_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) resetSequenceSnapshot
    if (.not.resetSequenceSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Virial_Orbital_Parameters_Benson2005_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Virial_Orbital_Parameters_Benson2005_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Virial_Orbital_Parameters_Benson2005_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) resetSequence
    if (.not.resetSequence) call Pseudo_Random_Retrieve(pseudoSequenceObject,fgslStateFile)
    return
  end subroutine Virial_Orbital_Parameters_Benson2005_State_Retrieve
  
end module Virial_Orbits_Benson2005
