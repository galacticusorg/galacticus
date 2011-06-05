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


!% Contains a module which implements the \cite{wetzel_orbits_2010} orbital parameter distribution for merging subhalos.

module Virial_Orbits_Wetzel2010
  !% Implements the \cite{wetzel_orbits_2010} orbital parameter distribution for merging subhalos.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  private
  public :: Virial_Orbital_Parameters_Wetzel2010_Initialize, Virial_Orbital_Parameters_Wetzel2010_Snapshot,&
       & Virial_Orbital_Parameters_Wetzel2010_State_Store, Virial_Orbital_Parameters_Wetzel2010_State_Retrieve

  type(fgsl_rng) :: pseudoSequenceObject,clonedPseudoSequenceObject
  logical        :: resetSequence=.true.,resetSequenceSnapshot
  !$omp threadprivate(pseudoSequenceObject,resetSequence)

  ! Table of the cumulative distribution for the pericentric radius.
  integer,          parameter                 :: pericentricRadiusPointsPerDecade=10
  integer                                     :: pericentricRadiusCount
  double precision, parameter                 :: pericentricRadiusMinimum=1.0d-6, pericentricRadiusMaximum=1.0d2
  double precision, allocatable, dimension(:) :: pericentricRadiusTableRadius,pericentricRadiusTableCumulativeProbability

  ! Parameters of the fitting functions.
  double precision, parameter                 :: circularityP1= 0.0d0,circularityAlpha1=0.242d0,circularityBeta1= 2.360d0,circularityGamma1=0.108d0, circularityGamma2=1.05d0
  double precision, parameter                 :: pericenterP1 =-4.0d0,pericenterAlpha1 =0.450d0,pericenterBeta1 =-0.395d0,pericenterGamma1 =0.109d0, pericenterGamma2 =0.85d0
  double precision, parameter                 :: r1Minimum    =0.05d0, c1Maximum=9.999999d0

  ! Global variables used in root finding.
  double precision                            :: uniformDeviate,C0,C1
  !$omp threadprivate(uniformDeviate,C0,C1)

contains

  !# <virialOrbitsMethod>
  !#  <unitName>Virial_Orbital_Parameters_Wetzel2010_Initialize</unitName>
  !# </virialOrbitsMethod>
  subroutine Virial_Orbital_Parameters_Wetzel2010_Initialize(virialOrbitsMethod,Virial_Orbital_Parameters_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Memory_Management
    use Numerical_Ranges
    use Hypergeometric_Functions
    use Kepler_Orbits_Structure
    implicit none
    type(varying_string),                  intent(in)    :: virialOrbitsMethod
    procedure(type(keplerOrbit)), pointer, intent(inout) :: Virial_Orbital_Parameters_Get
    integer                                              :: iRadius
    double precision                                     :: x,xGamma2
    
    if (virialOrbitsMethod == 'Wetzel2010') then
       ! Set procedure pointer to our orbital parameter function.
       Virial_Orbital_Parameters_Get => Virial_Orbital_Parameters_Wetzel2010

       ! Construct a look-up table for the pericentric radius distribution.
       ! Determine number of points to use in the tabulation.
       pericentricRadiusCount=int(dlog10(pericentricRadiusMaximum/pericentricRadiusMinimum)*dble(pericentricRadiusPointsPerDecade))+1
       ! Allocate space for the table.
       call Alloc_Array(pericentricRadiusTableRadius               ,[pericentricRadiusCount])
       call Alloc_Array(pericentricRadiusTableCumulativeProbability,[pericentricRadiusCount])
       ! Construct a range of radii.
       pericentricRadiusTableRadius=Make_Range(pericentricRadiusMinimum,pericentricRadiusMaximum,pericentricRadiusCount,rangeType=rangeTypeLogarithmic)
       ! For each radius, compute the cumulative probability.
       do iRadius=1,pericentricRadiusCount
          x      =pericentricRadiusTableRadius(iRadius)
          xGamma2=x**pericenterGamma2
          pericentricRadiusTableCumulativeProbability(iRadius)=dexp(-xGamma2)*x*(                                                                                                          &
               &  pericenterGamma2*(1.0d0+pericenterGamma2*(1.0d0+xGamma2))*Hypergeometric_1F1([2.0d0],[(1.0d0+3.0d0*pericenterGamma2)/pericenterGamma2],xGamma2)/(1.0d0+pericenterGamma2) &
               & +                                                          Hypergeometric_1F1([1.0d0],[(1.0d0+3.0d0*pericenterGamma2)/pericenterGamma2],xGamma2)*(1.0d0+pericenterGamma2) &
               &                                                                )
       end do
       ! Normalize to unit probability.
       pericentricRadiusTableCumulativeProbability=pericentricRadiusTableCumulativeProbability/pericentricRadiusTableCumulativeProbability(pericentricRadiusCount)
    end if
    return
  end subroutine Virial_Orbital_Parameters_Wetzel2010_Initialize

 function Virial_Orbital_Parameters_Wetzel2010(thisNode,hostNode,acceptUnboundOrbits) result (thisOrbit)
    !% Return orbital velocities of a satellite selected at random from the fitting function found by \cite{wetzel_orbits_2010}.
    use Pseudo_Random
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    use Critical_Overdensity
    use Numerical_Interpolation
    use Root_Finder
    use Cosmology_Functions
    use FGSL
    use Kepler_Orbits_Structure
    implicit none
    type(keplerOrbit)                                :: thisOrbit
    type(treeNode),          intent(inout), pointer  :: thisNode,hostNode
    logical,                 intent(in)              :: acceptUnboundOrbits
    double precision,        parameter               :: toleranceAbsolute =0.0d0, toleranceRelative =1.0d-2
    double precision,        parameter               :: circularityMinimum=0.0d0, circularityMaximum=1.0d0
    double precision,        parameter               :: redshiftMaximum   =5.0d0, expansionFactorMinimum=1.0d0/(1.0d0+redshiftMaximum)
    type(fgsl_interp),       save                    :: interpolationObject
    type(fgsl_interp_accel), save                    :: interpolationAccelerator
    logical,                 save                    :: interpolationReset=.true.
    !$omp threadprivate(interpolationObject,interpolationAccelerator,interpolationReset)
    type(fgsl_function),     save                    :: rootFunction
    type(fgsl_root_fsolver), save                    :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    type(c_ptr)                                      :: parameterPointer
    double precision                                 :: g1,R1,timeNode,massCharacteristic,expansionFactor&
         &,pericentricRadius,apocentricRadius,probabilityTotal,circularity,eccentricityInternal,radialScale
    logical                                          :: foundOrbit

    ! Reset the orbit.
    call thisOrbit%reset()
    ! Set masses and radius of the orbit.
    call thisOrbit%massesSet(Tree_Node_Mass(thisNode),Tree_Node_Mass(hostNode))
    call thisOrbit%radiusSet(Dark_Matter_Halo_Virial_Radius(hostNode))

    ! Get the time at which this node exists.
    timeNode=Tree_Node_Time(thisNode)

    ! Get the expansion factor.
    expansionFactor=Expansion_Factor(timeNode)

    ! Limit the expansion factor to the smallest value considered by Wetzel.
    if (expansionFactor < expansionFactorMinimum) then
       expansionFactor=              expansionFactorMinimum
       timeNode       =Cosmology_Age(expansionFactorMinimum)
    end if

    ! Get the characteristic mass, M*.
    massCharacteristic=Critical_Overdensity_Collapsing_Mass(timeNode)

    ! Compute parameter of the circularity fitting function. We limit C1 to a given maximum - the fit is not explored in this
    ! regime and without the truncation we get problems evaluating hypergeometric functions.
    g1=(1.0d0/expansionFactor)**circularityP1
    C1=min(circularityAlpha1*(1.0d0+circularityBeta1*(g1*thisOrbit%hostMass()/massCharacteristic)**circularityGamma1),c1Maximum)
    C0=1.0d0
    probabilityTotal=Circularity_Cumulative_Probability(circularityMaximum)
    C0=1.0d0/probabilityTotal

    ! Compute parameter of the pericentric distance fitting function. Since the fit for R1 can lead to negative pericentric
    ! distances in some cases we force R1 to always be above a specified minimum.
    g1=(1.0d0/expansionFactor)**pericenterP1
    R1=max(pericenterAlpha1*(1.0d0+pericenterBeta1*(g1*thisOrbit%hostMass()/massCharacteristic)**pericenterGamma1),r1Minimum)

    ! Search for an orbit.
    foundOrbit=.false.
    do while (.not.foundOrbit)
       
       ! Compute pericentric radius by inversion in table.
       uniformDeviate=Pseudo_Random_Get(pseudoSequenceObject,resetSequence)
       pericentricRadius=R1*Interpolate(pericentricRadiusCount,pericentricRadiusTableCumulativeProbability&
            &,pericentricRadiusTableRadius ,interpolationObject,interpolationAccelerator,uniformDeviate ,reset=interpolationReset)
       
       ! Compute circularity by root finding in the cumulative probability distribution.
       uniformDeviate=Pseudo_Random_Get(pseudoSequenceObject,resetSequence)
       circularity=Root_Find(circularityMinimum,circularityMaximum,Circularity_Root,parameterPointer,rootFunction &
            &,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
       
       ! Check that this is an orbit which actually reaches the virial radius.
       eccentricityInternal=dsqrt(1.0-circularity**2)
       apocentricRadius    =pericentricRadius*(1.0d0+eccentricityInternal)/(1.0d0-eccentricityInternal)
       foundOrbit=apocentricRadius >= 1.0d0 .and. pericentricRadius <= 1.0d0
    end do
    
    ! Get length scale for this orbit.
    radialScale  =Dark_Matter_Halo_Virial_Radius(hostNode)
    
    ! Set eccentricity and periapsis.
    call thisOrbit%eccentricitySet    (dsqrt(1.0-circularity**2)    )
    call thisOrbit%radiusPericenterSet(pericentricRadius*radialScale)
    
    return
  end function Virial_Orbital_Parameters_Wetzel2010

  function Circularity_Root(circularity,parameterPointer) bind(c)
    !% Function used in finding the circularity corresponding to a given cumulative probability.
    real(c_double)          :: Circularity_Root
    real(c_double), value   :: circularity
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: cumulativeProbability

    cumulativeProbability=Circularity_Cumulative_Probability(circularity)
    Circularity_Root=cumulativeProbability-uniformDeviate
    return
  end function Circularity_Root

  double precision function Circularity_Cumulative_Probability(circularity)
    !% The cumulative probability distribution for orbital circularity.
    use Hypergeometric_Functions
    implicit none
    double precision, intent(in) :: circularity

    Circularity_Cumulative_Probability=C0*(circularity**(circularityGamma2+1.0d0))*Hypergeometric_2F1([-C1,1.0d0&
         &+circularityGamma2],[2.0d0+circularityGamma2],circularity)/(circularityGamma2+1.0d0)
    return
  end function Circularity_Cumulative_Probability

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Virial_Orbital_Parameters_Wetzel2010_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Virial_Orbital_Parameters_Wetzel2010_Snapshot
    !% Store a snapshot of the random number generator internal state.
    implicit none

    if (.not.resetSequence) clonedPseudoSequenceObject=FGSL_Rng_Clone(pseudoSequenceObject)
    resetSequenceSnapshot=resetSequence
    return
  end subroutine Virial_Orbital_Parameters_Wetzel2010_Snapshot
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Virial_Orbital_Parameters_Wetzel2010_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Virial_Orbital_Parameters_Wetzel2010_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) resetSequenceSnapshot
    if (.not.resetSequenceSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Virial_Orbital_Parameters_Wetzel2010_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Virial_Orbital_Parameters_Wetzel2010_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Virial_Orbital_Parameters_Wetzel2010_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) resetSequence
    if (.not.resetSequence) call Pseudo_Random_Retrieve(pseudoSequenceObject,fgslStateFile)
    return
  end subroutine Virial_Orbital_Parameters_Wetzel2010_State_Retrieve
  
end module Virial_Orbits_Wetzel2010
