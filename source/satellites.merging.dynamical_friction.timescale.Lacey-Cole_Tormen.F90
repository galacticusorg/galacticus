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

!% Contains a module which implements calculations of satellite merging times using the \cite{lacey_merger_1993} method with a
!% parameterization of orbital parameters designed to fit the results of \cite{tormen_rise_1997} as described by 
!% \cite{cole_hierarchical_2000}.

module Dynamical_Friction_Lacey_Cole_Tormen
  !% Implements calculations of satellite merging times using the \cite{lacey_merger_1993} method with a parameterization of
  !% orbital parameters designed to fit the results of \cite{tormen_rise_1997} as described by \cite{cole_hierarchical_2000}.
  use FGSL
  implicit none
  private
  public :: Satellite_Time_Until_Merging_Lacey_Cole_Tormen_Initialize,&
       & Satellite_Time_Until_Merging_Lacey_Cole_Tormen_State_Retrieve,&
       & Satellite_Time_Until_Merging_Lacey_Cole_Tormen_State_Store, Satellite_Time_Until_Merging_Lacey_Cole_Tormen_Snapshot

  ! Random number objects
  type(fgsl_rng) :: randomSequenceObject,clonedPseudoSequenceObject
  logical        :: resetRandomSequence=.true.,resetRandomSequenceSnapshot
  !$omp threadprivate(resetRandomSequence,randomSequenceObject,clonedPseudoSequenceObject,resetRandomSequenceSnapshot)

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_Lacey_Cole_Tormen_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_Lacey_Cole_Tormen_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: satelliteMergingMethod
    procedure(double precision), pointer, intent(inout) :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'Lacey-Cole+Tormen') Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_Lacey_Cole_Tormen
    return
  end subroutine Satellite_Time_Until_Merging_Lacey_Cole_Tormen_Initialize

  double precision function Satellite_Time_Until_Merging_Lacey_Cole_Tormen(thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{lacey_merger_1993} method.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Dynamical_Friction_Timescale_Utilities
    use Kepler_Orbits
    use Gaussian_Random
    implicit none
    type (treeNode          ), pointer, intent(inout) :: thisNode
    type (keplerOrbit       ),          intent(inout) :: thisOrbit
    type (treeNode          ), pointer                :: hostNode
    class(nodeComponentBasic), pointer                :: thisBasicComponent,hostBasicComponent
    double precision,          parameter              :: inverseTwoB1=1.169335453d0 ! 1/2/B(1).
    double precision,          parameter              :: orbitalFactorDistributionSigma= 0.26d0 ! Cole et al. (2000).
    double precision,          parameter              :: orbitalFactorDistributionMean =-0.14d0 ! Cole et al. (2000).
    double precision                                  :: massRatio,randomDeviate,log10OrbitalFactor,orbitalFactor

    ! Find the host node.
    hostNode => thisNode%parent
    ! Compute the orbital factor - selected at random from a lognormal distribution.
    randomDeviate=Gaussian_Random_Get(randomSequenceObject,orbitalFactorDistributionSigma,resetRandomSequence)
    log10OrbitalFactor=orbitalFactorDistributionMean+randomDeviate
    orbitalFactor=10.0d0**log10OrbitalFactor
    ! Compute mass ratio.
    thisBasicComponent => thisNode%basic()
    hostBasicComponent => hostNode%basic()
    massRatio=hostBasicComponent%mass()/thisBasicComponent%mass()
    ! Check for a greater than unity mass ratio.
    if (massRatio <= 1.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       Satellite_Time_Until_Merging_Lacey_Cole_Tormen=0.0d0
    else
       ! Compute dynamical friction timescale.
       Satellite_Time_Until_Merging_Lacey_Cole_Tormen=Dynamical_Friction_Timescale_Multiplier()*orbitalFactor&
            &*Dark_Matter_Halo_Dynamical_Timescale(hostNode)*inverseTwoB1*massRatio /dlog(massRatio)
    end if
  return
  end function Satellite_Time_Until_Merging_Lacey_Cole_Tormen

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Satellite_Time_Until_Merging_Lacey_Cole_Tormen_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Satellite_Time_Until_Merging_Lacey_Cole_Tormen_Snapshot
    !% Store a snapshot of the random number generator internal state.
    implicit none

    if (.not.resetRandomSequence) clonedPseudoSequenceObject=FGSL_Rng_Clone(randomSequenceObject)
    resetRandomSequenceSnapshot=resetRandomSequence
    return
  end subroutine Satellite_Time_Until_Merging_Lacey_Cole_Tormen_Snapshot
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Satellite_Time_Until_Merging_Lacey_Cole_Tormen_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Satellite_Time_Until_Merging_Lacey_Cole_Tormen_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) resetRandomSequenceSnapshot
    if (.not.resetRandomSequenceSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Satellite_Time_Until_Merging_Lacey_Cole_Tormen_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Satellite_Time_Until_Merging_Lacey_Cole_Tormen_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Satellite_Time_Until_Merging_Lacey_Cole_Tormen_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) resetRandomSequence
    if (.not.resetRandomSequence) call Pseudo_Random_Retrieve(randomSequenceObject,fgslStateFile)
    return
  end subroutine Satellite_Time_Until_Merging_Lacey_Cole_Tormen_State_Retrieve
    
end module Dynamical_Friction_Lacey_Cole_Tormen
