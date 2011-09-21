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
  !$omp threadprivate(resetRandomSequence,randomSequenceObject)

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
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Dynamical_Friction_Timescale_Utilities
    use Kepler_Orbits_Structure
    use Gaussian_Random
    implicit none
    type(treeNode),    pointer, intent(inout) :: thisNode
    type(keplerOrbit),          intent(inout) :: thisOrbit
    type(treeNode),    pointer                :: hostNode
    double precision,  parameter              :: inverseTwoB1=1.169335453d0 ! 1/2/B(1).
    double precision,  parameter              :: orbitalFactorDistributionSigma= 0.26d0 ! Cole et al. (2000).
    double precision,  parameter              :: orbitalFactorDistributionMean =-0.14d0 ! Cole et al. (2000).
    double precision                          :: massRatio,randomDeviate,log10OrbitalFactor,orbitalFactor

    ! Find the host node.
    hostNode => thisNode%parentNode
    ! Compute the orbital factor - selected at random from a lognormal distribution.
    randomDeviate=Gaussian_Random_Get(randomSequenceObject,orbitalFactorDistributionSigma,resetRandomSequence)
    log10OrbitalFactor=orbitalFactorDistributionMean+randomDeviate
    orbitalFactor=10.0d0**log10OrbitalFactor
    ! Compute mass ratio.
    massRatio=Tree_Node_Mass(hostNode)/Tree_Node_Mass(thisNode)
    ! Compute dynamical friction timescale.
    Satellite_Time_Until_Merging_Lacey_Cole_Tormen=Dynamical_Friction_Timescale_Multiplier()*orbitalFactor&
         &*Dark_Matter_Halo_Dynamical_Timescale(hostNode)*inverseTwoB1*massRatio /dlog(massRatio)
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
