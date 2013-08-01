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

!% Contains a module which implements calculations of satellite merging times using the \cite{jiang_fitting_2008} method.

module Dynamical_Friction_Jiang2008
  !% Implements calculations of satellite merging times using the \cite{jiang_fitting_2008} method.
  use FGSL
  implicit none
  private
  public :: Satellite_Time_Until_Merging_Jiang2008_Initialize,Satellite_Time_Until_Merging_Jiang2008_Snapshot&
       &,Satellite_Time_Until_Merging_Jiang2008_State_Store,Satellite_Time_Until_Merging_Jiang2008_State_Retrieve

  ! Scatter (in log(T_merge)) to add to the merger times.
  double precision :: satelliteMergingJiang2008Scatter
  ! Random number objects
  type(fgsl_rng) :: randomSequenceObject,clonedPseudoSequenceObject
  logical        :: resetRandomSequence=.true.,resetRandomSequenceSnapshot
  !$omp threadprivate(resetRandomSequence,randomSequenceObject,clonedPseudoSequenceObject,resetRandomSequenceSnapshot)

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_Jiang2008_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_Jiang2008_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use Galacticus_Nodes
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type     (varying_string                        ), intent(in   )          :: satelliteMergingMethod
    procedure(Satellite_Time_Until_Merging_Jiang2008), intent(inout), pointer :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'Jiang2008') then
       ! Set the function pointer to our implementation.
       Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_Jiang2008
       ! Get input parameters.
          !@ <inputParameter>
          !@   <name>satelliteMergingJiang2008Scatter</name>
          !@   <defaultValue>$0$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to add random scatter to the dynamical friction timescales in the {\tt Jiang2008} satellite merging time implementation.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteMergingJiang2008Scatter',satelliteMergingJiang2008Scatter,defaultValue=0.0d0)
          ! Check that required properties are gettable.
          if (.not.defaultBasicComponent%massIsGettable()) call Galacticus_Error_Report('Satellite_Time_Until_Merging_Jiang2008_Initialize','this method requires that the "mass" property of the basic component be gettable')
    end if
    return
  end subroutine Satellite_Time_Until_Merging_Jiang2008_Initialize

  double precision function Satellite_Time_Until_Merging_Jiang2008(thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{jiang_fitting_2008} method.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles
    use Dynamical_Friction_Timescale_Utilities
    use Kepler_Orbits
    use Satellite_Orbits
    use Gaussian_Random
    implicit none
    type            (treeNode          )           , intent(inout), pointer :: thisNode
    type            (keplerOrbit       )           , intent(inout)          :: thisOrbit
    type            (treeNode          )                          , pointer :: hostNode
    class           (nodeComponentBasic)                          , pointer :: hostBasicComponent                   , thisBasicComponent
    logical                             , parameter                         :: acceptUnboundOrbits          =.false.
    double precision                    , parameter                         :: C                            =0.43d0 , a                 =0.94d0, &  !   Fitting parameters from Jiang's paper.
         &                                                                     b                            =0.60d0 , d                 =0.60d0
    double precision                                                        :: equivalentCircularOrbitRadius        , massRatio                , &
         &                                                                     orbitalCircularity                   , radialScale              , &
         &                                                                     velocityScale                        , randomDeviate                                                                                        

    ! Find the host node.
    hostNode => thisNode%parent
    ! Get the equivalent circular orbit.
    equivalentCircularOrbitRadius=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(hostNode,thisOrbit)
    ! Get velocity scale.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale  =Dark_Matter_Halo_Virial_Radius  (hostNode)
    ! Compute orbital circularity.
    orbitalCircularity=thisOrbit%angularMomentum()/equivalentCircularOrbitRadius/Dark_Matter_Profile_Circular_Velocity(hostNode&
         &,equivalentCircularOrbitRadius)
    ! Compute mass ratio (mass in host [not including satellite] divided by mass in satellite).
    thisBasicComponent => thisNode%basic()
    hostBasicComponent => hostNode%basic()
    massRatio=hostBasicComponent%mass()/thisBasicComponent%mass()-1.0d0
    ! Check for a non-zero mass ratio.
    if (massRatio <= 0.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       Satellite_Time_Until_Merging_Jiang2008=0.0d0
    else
       ! Compute dynamical friction timescale.
       Satellite_Time_Until_Merging_Jiang2008=Dynamical_Friction_Timescale_Multiplier()&
            &*Dark_Matter_Halo_Dynamical_Timescale(hostNode)*sqrt(equivalentCircularOrbitRadius/radialScale)*((a&
            &*(orbitalCircularity**b)+d)/2.0/C)*massRatio/log(1.0d0+massRatio)
       ! Add scatter if necessary.
       if (satelliteMergingJiang2008Scatter > 0.0d0) then
          randomDeviate=Gaussian_Random_Get(randomSequenceObject,satelliteMergingJiang2008Scatter,resetRandomSequence)
          Satellite_Time_Until_Merging_Jiang2008=Satellite_Time_Until_Merging_Jiang2008*exp(randomDeviate)
       end if
    end if
    return
  end function Satellite_Time_Until_Merging_Jiang2008

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Satellite_Time_Until_Merging_Jiang2008_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Satellite_Time_Until_Merging_Jiang2008_Snapshot
    !% Store a snapshot of the random number generator internal state.
    implicit none

    if (.not.resetRandomSequence) clonedPseudoSequenceObject=FGSL_Rng_Clone(randomSequenceObject)
    resetRandomSequenceSnapshot=resetRandomSequence
    return
  end subroutine Satellite_Time_Until_Merging_Jiang2008_Snapshot
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Satellite_Time_Until_Merging_Jiang2008_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Satellite_Time_Until_Merging_Jiang2008_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) resetRandomSequenceSnapshot
    if (.not.resetRandomSequenceSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Satellite_Time_Until_Merging_Jiang2008_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Satellite_Time_Until_Merging_Jiang2008_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Satellite_Time_Until_Merging_Jiang2008_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) resetRandomSequence
    if (.not.resetRandomSequence) call Pseudo_Random_Retrieve(randomSequenceObject,fgslStateFile)
    return
  end subroutine Satellite_Time_Until_Merging_Jiang2008_State_Retrieve
    
end module Dynamical_Friction_Jiang2008
