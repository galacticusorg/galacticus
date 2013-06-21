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

!% Contains a module which implements the \cite{benson_orbital_2005} orbital parameter distribution for merging subhalos.

module Virial_Orbits_Benson2005
  !% Implements the \cite{benson_orbital_2005} orbital parameter distribution for merging subhalos.
  use FGSL
  implicit none
  private
  public :: Virial_Orbital_Parameters_Benson2005_Initialize, Virial_Orbital_Parameters_Benson2005_Snapshot,&
       & Virial_Orbital_Parameters_Benson2005_State_Store, Virial_Orbital_Parameters_Benson2005_State_Retrieve

  type   (fgsl_rng) :: clonedPseudoSequenceObject       , pseudoSequenceObject  
  logical           :: resetSequence             =.true., resetSequenceSnapshot 
  !$omp threadprivate(pseudoSequenceObject,resetSequence,clonedPseudoSequenceObject,resetSequenceSnapshot)
contains

  !# <virialOrbitsMethod>
  !#  <unitName>Virial_Orbital_Parameters_Benson2005_Initialize</unitName>
  !# </virialOrbitsMethod>
  subroutine Virial_Orbital_Parameters_Benson2005_Initialize(virialOrbitsMethod,Virial_Orbital_Parameters_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Kepler_Orbits
    implicit none
    type     (varying_string                      ), intent(in   )          :: virialOrbitsMethod            
    procedure(Virial_Orbital_Parameters_Benson2005), intent(inout), pointer :: Virial_Orbital_Parameters_Get 
    
    if (virialOrbitsMethod == 'Benson2005') Virial_Orbital_Parameters_Get => Virial_Orbital_Parameters_Benson2005
    return
  end subroutine Virial_Orbital_Parameters_Benson2005_Initialize

  function Virial_Orbital_Parameters_Benson2005(thisNode,hostNode,acceptUnboundOrbits) result (thisOrbit)
    !% Return orbital parameters of a satellite selected at random from the fitting function found by \cite{benson_orbital_2005}.
    use Pseudo_Random
    use Kepler_Orbits
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    implicit none
    type            (keplerOrbit       )                                    :: thisOrbit                                                                                                                                                                                                                                         
    type            (treeNode          )           , intent(inout), pointer :: hostNode                                                                                                                                        , thisNode                                                                                        
    logical                                        , intent(in   )          :: acceptUnboundOrbits                                                                                                                                                                                                                               
    class           (nodeComponentBasic)                          , pointer :: hostBasicComponent                                                                                                                              , thisBasicComponent                                                                              
    double precision                    , parameter                         :: pMax                   =1.96797d0                                                                                                               , velocityMax                                                                           =3.0d0    
    double precision                    , parameter                         :: a                   (9)=(/0.390052d+01,0.247973d+01,0.102373d+02,0.683922d+00,0.353953d+00,0.107716d+01,0.509837d+00,0.206204d+00,0.314641d+00/)                                                                                                  
    double precision                    , parameter                         :: EPS_BOUND              =1.0d-4                                                                                                                                               !   Tolerence to ensure that orbits are sufficiently bound.          
    double precision                                                        :: b1                                                                                                                                              , b2                                                                                          , & 
         &                                                                     distributionFunction                                                                                                                            , energyInternal                                                                              , & 
         &                                                                     uniformRandom                                                                                                                                   , velocityRadialInternal                                                                      , & 
         &                                                                     velocityScale                                                                                                                                   , velocityTangentialInternal                                                                      
    logical                                                                 :: foundOrbit                                                                                                                                                                                                                                        
    
    ! Reset the orbit.
    call thisOrbit%reset()
    ! Set masses and radius of the orbit.
    thisBasicComponent => thisNode%basic()
    hostBasicComponent => hostNode%basic()
    call thisOrbit%massesSet(thisBasicComponent%mass(),hostBasicComponent%mass())
    call thisOrbit%radiusSet(Dark_Matter_Halo_Virial_Radius(hostNode))

    ! Select an orbit.
    foundOrbit=.false.
    do while(.not.foundOrbit)
       ! Select potential radial and tangential velocities.
       velocityRadialInternal    =Pseudo_Random_Get(pseudoSequenceObject,resetSequence)*velocityMax
       velocityTangentialInternal=Pseudo_Random_Get(pseudoSequenceObject,resetSequence)*velocityMax
       ! Evaluate distribution function for these parameters.
       b1=a(3)*exp(-a(4)*(velocityTangentialInternal-a(5))**2)
       b2=a(6)*exp(-a(7)*(velocityTangentialInternal-a(8))**2)
       distributionFunction=a(1)*velocityTangentialInternal*exp(-a(2)*((velocityTangentialInternal-a(9))**2.0))*exp(-b1&
            &*(velocityRadialInternal-b2)**2)
       if (distributionFunction > pMax) call Galacticus_Error_Report('Virial_Orbital_Parameters_Benson2005','distribution&
            & function exceeds expected peak value')
       uniformRandom=pMax*Pseudo_Random_Get(pseudoSequenceObject,resetSequence)
       if (uniformRandom <= distributionFunction) then
          foundOrbit=.true.
          ! If requested, check that the orbit is bound. We require it to have E<-EPS_BOUND to ensure that it is sufficiently
          ! bound that later rounding errors will not make it appear unbound.
          if (.not.acceptUnboundOrbits) then
             energyInternal=-1.0d0+0.5d0*(velocityRadialInternal**2+velocityTangentialInternal**2)*thisOrbit%specificReducedMass()
             foundOrbit=(energyInternal < -EPS_BOUND)
          end if
       end if
    end do
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    call thisOrbit%velocityRadialSet    (velocityRadialInternal    *velocityScale)
    call thisOrbit%velocityTangentialSet(velocityTangentialInternal*velocityScale)
    return
  end function Virial_Orbital_Parameters_Benson2005

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
    integer           , intent(in   ) :: stateFile     
    type   (fgsl_file), intent(in   ) :: fgslStateFile 
    
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
    integer           , intent(in   ) :: stateFile     
    type   (fgsl_file), intent(in   ) :: fgslStateFile 
    
    read (stateFile) resetSequence
    if (.not.resetSequence) call Pseudo_Random_Retrieve(pseudoSequenceObject,fgslStateFile)
    return
  end subroutine Virial_Orbital_Parameters_Benson2005_State_Retrieve
  
end module Virial_Orbits_Benson2005
