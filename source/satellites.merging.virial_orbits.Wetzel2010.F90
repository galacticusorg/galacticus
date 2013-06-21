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

!% Contains a module which implements the \cite{wetzel_orbits_2010} orbital parameter distribution for merging subhalos.

module Virial_Orbits_Wetzel2010
  !% Implements the \cite{wetzel_orbits_2010} orbital parameter distribution for merging subhalos.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Virial_Orbital_Parameters_Wetzel2010_Initialize, Virial_Orbital_Parameters_Wetzel2010_Snapshot,&
       & Virial_Orbital_Parameters_Wetzel2010_State_Store, Virial_Orbital_Parameters_Wetzel2010_State_Retrieve

  type            (fgsl_rng)                            :: clonedPseudoSequenceObject                            , pseudoSequenceObject                     
  logical                                               :: resetSequence                              =.true.    , resetSequenceSnapshot                    
  !$omp threadprivate(pseudoSequenceObject,resetSequence,clonedPseudoSequenceObject,resetSequenceSnapshot)
  ! Table of the cumulative distribution for the pericentric radius.
  integer                   , parameter                 :: pericentricRadiusPointsPerDecade           =10                                                   
  integer                                               :: pericentricRadiusCount                                                                           
  double precision          , parameter                 :: pericentricRadiusMaximum                   =1.0d2     , pericentricRadiusMinimum    =1.0d-6      
  double precision          , allocatable, dimension(:) :: pericentricRadiusTableCumulativeProbability           , pericentricRadiusTableRadius             
  
  ! Parameters of the fitting functions.
  double precision          , parameter                 :: circularityAlpha1                          =0.242d0   , circularityBeta1            =2.360d0 , & 
       &                                                   circularityGamma1                          =0.108d0   , circularityGamma2           =1.05d0  , & 
       &                                                   circularityP1                              =0.0d0                                                
  double precision          , parameter                 :: pericenterAlpha1                           =0.450d0   , pericenterBeta1             =-0.395d0, & 
       &                                                   pericenterGamma1                           =0.109d0   , pericenterGamma2            =0.85d0  , & 
       &                                                   pericenterP1                               =-4.0d0                                               
  double precision          , parameter                 :: c1Maximum                                  =9.999999d0, r1Minimum                   =0.05d0      
  
  ! Global variables used in root finding.
  double precision                                      :: C0                                                    , C1                                   , & 
       &                                                   uniformDeviate                                                                                   
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
    use Kepler_Orbits
    implicit none
    type            (varying_string                      ), intent(in   )          :: virialOrbitsMethod                     
    procedure       (Virial_Orbital_Parameters_Wetzel2010), intent(inout), pointer :: Virial_Orbital_Parameters_Get          
    integer                                                                        :: iRadius                                
    double precision                                                               :: x                            , xGamma2 
    
    if (virialOrbitsMethod == 'Wetzel2010') then
       ! Set procedure pointer to our orbital parameter function.
       Virial_Orbital_Parameters_Get => Virial_Orbital_Parameters_Wetzel2010

       ! Construct a look-up table for the pericentric radius distribution.
       ! Determine number of points to use in the tabulation.
       pericentricRadiusCount=int(log10(pericentricRadiusMaximum/pericentricRadiusMinimum)*dble(pericentricRadiusPointsPerDecade))+1
       ! Allocate space for the table.
       call Alloc_Array(pericentricRadiusTableRadius               ,[pericentricRadiusCount])
       call Alloc_Array(pericentricRadiusTableCumulativeProbability,[pericentricRadiusCount])
       ! Construct a range of radii.
       pericentricRadiusTableRadius=Make_Range(pericentricRadiusMinimum,pericentricRadiusMaximum,pericentricRadiusCount,rangeType=rangeTypeLogarithmic)
       ! For each radius, compute the cumulative probability.
       do iRadius=1,pericentricRadiusCount
          x      =pericentricRadiusTableRadius(iRadius)
          xGamma2=x**pericenterGamma2
          pericentricRadiusTableCumulativeProbability(iRadius)=exp(-xGamma2)*x*(                                                                                                          &
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
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Critical_Overdensity
    use Numerical_Interpolation
    use Root_Finder
    use Cosmology_Functions
    use Kepler_Orbits
    implicit none
    type            (keplerOrbit       )                                          :: thisOrbit                                                                              
    type            (treeNode          )                 , intent(inout), pointer :: hostNode                                              , thisNode                       
    logical                                              , intent(in   )          :: acceptUnboundOrbits                                                                    
    class           (nodeComponentBasic)                                , pointer :: hostBasicComponent                                    , thisBasicComponent             
    double precision                          , parameter                         :: toleranceAbsolute       =0.0d0                        , toleranceRelative   =1.0d-2    
    double precision                          , parameter                         :: circularityMaximum      =1.0d0                        , circularityMinimum  =0.0d0     
    double precision                          , parameter                         :: expansionFactorMinimum  =1.0d0/(1.0d0+redshiftMaximum), redshiftMaximum     =5.0d0     
    type            (fgsl_interp       ), save                                    :: interpolationObject                                                                    
    type            (fgsl_interp_accel ), save                                    :: interpolationAccelerator                                                               
    logical                             , save                                    :: interpolationReset      =.true.                                                        
    !$omp threadprivate(interpolationObject,interpolationAccelerator,interpolationReset)
    type            (rootFinder        ), save                                    :: finder                                                                                 
    !$omp threadprivate(finder)
    double precision                                                              :: R1                                                    , apocentricRadius           , & 
         &                                                                           circularity                                           , eccentricityInternal       , & 
         &                                                                           expansionFactor                                       , g1                         , & 
         &                                                                           massCharacteristic                                    , pericentricRadius          , & 
         &                                                                           probabilityTotal                                      , radialScale                , & 
         &                                                                           timeNode                                                                               
    logical                                                                       :: foundOrbit                                                                             
    
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(Circularity_Root                   )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
    end if

    ! Reset the orbit.
    call thisOrbit%reset()
    ! Set masses and radius of the orbit.
    thisBasicComponent => thisNode%basic()
    hostBasicComponent => hostNode%basic()
    call thisOrbit%massesSet(thisBasicComponent%mass(),hostBasicComponent%mass())
    call thisOrbit%radiusSet(Dark_Matter_Halo_Virial_Radius(hostNode))

    ! Get the time at which this node exists.
    timeNode=thisBasicComponent%time()

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
       circularity=finder%find(rootRange=[circularityMinimum,circularityMaximum])
       ! Check that this is an orbit which actually reaches the virial radius.
       eccentricityInternal=sqrt(1.0-circularity**2)
       apocentricRadius    =pericentricRadius*(1.0d0+eccentricityInternal)/(1.0d0-eccentricityInternal)
       foundOrbit=apocentricRadius >= 1.0d0 .and. pericentricRadius <= 1.0d0
    end do
    
    ! Get length scale for this orbit.
    radialScale  =Dark_Matter_Halo_Virial_Radius(hostNode)
    
    ! Set eccentricity and periapsis.
    call thisOrbit%eccentricitySet    (sqrt(1.0-circularity**2)    )
    call thisOrbit%radiusPericenterSet(pericentricRadius*radialScale)
    
    return
  end function Virial_Orbital_Parameters_Wetzel2010

  double precision function Circularity_Root(circularity)
    !% Function used in finding the circularity corresponding to a given cumulative probability.
    double precision, intent(in   ) :: circularity           
    double precision                :: cumulativeProbability 
    
    cumulativeProbability=Circularity_Cumulative_Probability(circularity)
    Circularity_Root=cumulativeProbability-uniformDeviate
    return
  end function Circularity_Root

  double precision function Circularity_Cumulative_Probability(circularity)
    !% The cumulative probability distribution for orbital circularity.
    use Hypergeometric_Functions
    implicit none
    double precision, intent(in   ) :: circularity 
    
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
    use Pseudo_Random
    implicit none
    integer           , intent(in   ) :: stateFile     
    type   (fgsl_file), intent(in   ) :: fgslStateFile 
    
    write (stateFile) resetSequenceSnapshot
    if (.not.resetSequenceSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Virial_Orbital_Parameters_Wetzel2010_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Virial_Orbital_Parameters_Wetzel2010_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Virial_Orbital_Parameters_Wetzel2010_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    integer           , intent(in   ) :: stateFile     
    type   (fgsl_file), intent(in   ) :: fgslStateFile 
    
    read (stateFile) resetSequence
    if (.not.resetSequence) call Pseudo_Random_Retrieve(pseudoSequenceObject,fgslStateFile)
    return
  end subroutine Virial_Orbital_Parameters_Wetzel2010_State_Retrieve
  
end module Virial_Orbits_Wetzel2010
