!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Implements calculations of satellite merging times using the \cite{jiang_fitting_2008} method.
  use FGSL
  
  !# <satelliteMergingTimescales name="satelliteMergingTimescalesJiang2008">
  !#  <description>Computes the merging timescale using the method of \cite{jiang_fitting_2008}.</description>
  !# </satelliteMergingTimescales>

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesJiang2008
     !% A class implementing the \cite{jiang_fitting_2008} method for satellite merging timescales.
     private
     ! Scatter (in log(T_merge)) to add to the merger times.
     double precision           :: scatter
     ! Random number objects
     type            (fgsl_rng) :: randomSequenceObject, clonedPseudoSequenceObject
     logical                    :: resetRandomSequence , resetRandomSequenceSnapshot
   contains
     final     ::                     jiang2008Destructor
     procedure :: stateStore       => jiang2008StateStore
     procedure :: stateRestore     => jiang2008StateRestore
     procedure :: stateSnapshot    => jiang2008StateSnapshot
     procedure :: timeUntilMerging => jiang2008TimeUntilMerging
  end type satelliteMergingTimescalesJiang2008

  interface satelliteMergingTimescalesJiang2008
     !% Constructors for the \cite{jiang_fitting_2008} merging timescale class.
     module procedure jiang2008DefaultConstructor
     module procedure jiang2008Constructor
  end interface satelliteMergingTimescalesJiang2008

  ! Initialization state.
  logical          :: jiang2008Initialized=.false.

  ! Default scatter.
  double precision :: satelliteMergingJiang2008Scatter

contains

  function jiang2008DefaultConstructor()
    !% Default constructor for the \cite{jiang_fitting_2008} merging timescale class.
    use Galacticus_Display
    use Input_Parameters
    use Galacticus_Nodes
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(satelliteMergingTimescalesJiang2008) :: jiang2008DefaultConstructor

    if (.not.jiang2008Initialized) then
       !$omp critical (satelliteMergingTimescalesJiang2008Initialize)
       if (.not.jiang2008Initialized) then
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
          if (.not.defaultBasicComponent%massIsGettable()) call Galacticus_Error_Report('jiang2008DefaultConstructor','this method requires that the "mass" property of the basic component be gettable')
          ! Record that we are now initialized.
          jiang2008Initialized=.true.
       end if
       !$omp end critical (satelliteMergingTimescalesJiang2008Initialize)
    end if
    jiang2008DefaultConstructor=jiang2008Constructor(satelliteMergingJiang2008Scatter)
   return
  end function jiang2008DefaultConstructor

  function jiang2008Constructor(scatter)
    !% Constructor for the \cite{jiang_fitting_2008} merging timescale class.
    implicit none
    type            (satelliteMergingTimescalesJiang2008)                :: jiang2008Constructor
    double precision                                     , intent(in   ) :: scatter

    jiang2008Constructor%resetRandomSequence=.true.
    jiang2008Constructor%scatter            =scatter
    return
  end function jiang2008Constructor

  subroutine jiang2008Destructor(self)
    !% Destructor for the \cite{jiang_fitting_2008} merging timescale class.
     use Gaussian_Random
     implicit none
     type(satelliteMergingTimescalesJiang2008), intent(inout) :: self

     ! Destroy the random number object.
     if (self%resetRandomSequence        ) call Gaussian_Random_Free(self%randomSequenceObject      )
     if (self%resetRandomSequenceSnapshot) call Gaussian_Random_Free(self%clonedPseudoSequenceObject)
    return
  end subroutine jiang2008Destructor

  double precision function jiang2008TimeUntilMerging(self,thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{jiang_fitting_2008} method.
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles
    use Dynamical_Friction_Timescale_Utilities
    use Satellite_Orbits
    use Gaussian_Random
    use Galacticus_Error
    implicit none
    class           (satelliteMergingTimescalesJiang2008)           , intent(inout)          :: self
    type            (treeNode                           )           , intent(inout), pointer :: thisNode
    type            (keplerOrbit                        )           , intent(inout)          :: thisOrbit
    type            (treeNode                           )                          , pointer :: hostNode
    class           (nodeComponentBasic                 )                          , pointer :: hostBasic                            , thisBasic
    logical                                              , parameter                         :: acceptUnboundOrbits          =.false.
    double precision                                     , parameter                         :: timeInfinite                 =1.0d30

    double precision                                     , parameter                         :: C                            =0.43d0 , a                 =0.94d0, &  !   Fitting parameters from Jiang's paper.
         &                                                                                      b                            =0.60d0 , d                 =0.60d0
    integer                                                                                  :: errorCode
    double precision                                                                         :: equivalentCircularOrbitRadius        , massRatio                , &
         &                                                                                      orbitalCircularity                   , radialScale              , &
         &                                                                                      velocityScale                        , randomDeviate

    ! Find the host node.
    hostNode => thisNode%parent
    ! Get the equivalent circular orbit.
    equivalentCircularOrbitRadius=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(hostNode,thisOrbit,errorCode)
    ! Check error codes.
    select case (errorCode)
    case (errorCodeOrbitUnbound     )
       jiang2008TimeUntilMerging=timeInfinite
       return
    case (errorCodeNoEquivalentOrbit)
       ! Circularity is not defined. Assume instantaneous merging.
       jiang2008TimeUntilMerging=0.0d0
       return
    case (errorCodeSuccess          )
    case default
       call Galacticus_Error_Report('jiang2008TimeUntilMerging','unrecognized error code')
    end select
    ! Get velocity scale.
    velocityScale=Dark_Matter_Halo_Virial_Velocity(hostNode)
    radialScale  =Dark_Matter_Halo_Virial_Radius  (hostNode)
    ! Compute orbital circularity.
    orbitalCircularity= thisOrbit%angularMomentum()                                                   &
         &             /equivalentCircularOrbitRadius                                                 &
         &             /Dark_Matter_Profile_Circular_Velocity(hostNode,equivalentCircularOrbitRadius)
    ! Compute mass ratio (mass in host [not including satellite] divided by mass in satellite).
    thisBasic => thisNode%basic()
    hostBasic => hostNode%basic()
    massRatio=hostBasic%mass()/thisBasic%mass()-1.0d0
    ! Check for a non-zero mass ratio.
    if (massRatio <= 0.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       jiang2008TimeUntilMerging=0.0d0
    else
       ! Compute dynamical friction timescale.
       jiang2008TimeUntilMerging       &
            & =Dynamical_Friction_Timescale_Multiplier(        ) &
            & *Dark_Matter_Halo_Dynamical_Timescale   (hostNode) &
            & *sqrt(equivalentCircularOrbitRadius/radialScale)   &
            & *((a*(orbitalCircularity**b)+d)/2.0d0/C)           &
            & *          massRatio                               &
            & /log(1.0d0+massRatio)
       ! Add scatter if necessary.
       if (self%scatter > 0.0d0) then
          randomDeviate=Gaussian_Random_Get(self%randomSequenceObject,self%scatter,self%resetRandomSequence)
          jiang2008TimeUntilMerging=jiang2008TimeUntilMerging*exp(randomDeviate)
       end if
    end if
    return
  end function jiang2008TimeUntilMerging

  subroutine jiang2008StateSnapshot(self)
    !% Store a snapshot of the random number generator internal state.
    implicit none
    class(satelliteMergingTimescalesJiang2008), intent(inout) :: self

    if (.not.self%resetRandomSequence) self%clonedPseudoSequenceObject=FGSL_Rng_Clone(self%randomSequenceObject)
    self%resetRandomSequenceSnapshot=self%resetRandomSequence
    return
  end subroutine jiang2008StateSnapshot

  subroutine jiang2008StateStore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    class  (satelliteMergingTimescalesJiang2008), intent(inout) :: self
    integer                                     , intent(in   ) :: stateFile
    type   (fgsl_file                          ), intent(in   ) :: fgslStateFile
    
    write (stateFile) self%resetRandomSequenceSnapshot
    if (.not.self%resetRandomSequenceSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine jiang2008StateStore
  
  subroutine jiang2008StateRestore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    class  (satelliteMergingTimescalesJiang2008), intent(inout) :: self
    integer                                     , intent(in   ) :: stateFile
    type   (fgsl_file                          ), intent(in   ) :: fgslStateFile
    
    read (stateFile) self%resetRandomSequence
    if (.not.self%resetRandomSequence) call Pseudo_Random_Retrieve(self%randomSequenceObject,fgslStateFile)
    return
  end subroutine jiang2008StateRestore
