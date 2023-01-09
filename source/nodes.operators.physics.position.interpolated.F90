!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !!{
  Implements a node operator class that interpolates positions of nodes using the approach of \cite{merson_lightcone_2013}.
  !!}
  
  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <nodeOperator name="nodeOperatorPositionInterpolated">
   <description>
    A node operator class that interpolates positions of nodes using the approach of \cite{merson_lightcone_2013}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorPositionInterpolated
     !!{
     A node operator class that interpolates positions of nodes using the approach of \cite{merson_lightcone_2013}.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: lengthBox
     logical                                            :: isPeriodic
     integer                                            :: timeMaximumID                , coefficientsID
   contains
     !![
     <methods>
       <method method="computeInterpolation" description="Compute the interpolation for the given node."/>
     </methods>
     !!]
     final     ::                                        positionInterpolatedDestructor
     procedure :: differentialEvolutionAnalytics      => positionInterpolatedDifferentialEvolutionAnalytics
     procedure :: differentialEvolutionSolveAnalytics => positionInterpolatedDifferentialEvolutionSolveAnalytics
     procedure :: differentialEvolution               => positionInterpolatedDifferentialEvolution
     procedure :: differentialEvolutionPost           => positionInterpolatedDifferentialEvolutionPost
     procedure :: nodePromote                         => positionInterpolatedNodePromote
     procedure :: nodesMerge                          => positionInterpolatedNodesMerge
     procedure :: nodeInitialize                      => positionInterpolatedNodeInitialize
     procedure :: autoHook                            => positionInterpolatedAutoHook
     procedure :: computeInterpolation                => positionInterpolatedComputeInterpolation
  end type nodeOperatorPositionInterpolated
  
  interface nodeOperatorPositionInterpolated
     !!{
     Constructors for the {\normalfont \ttfamily positionInterpolated} node operator class.
     !!}
     module procedure positionInterpolatedConstructorParameters
     module procedure positionInterpolatedConstructorInternal
  end interface nodeOperatorPositionInterpolated

  !! Submodule-scope pointer to self.
  class(nodeOperatorPositionInterpolated), pointer :: self_
  !$omp threadprivate(self_)
  
contains
  
  function positionInterpolatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily positionInterpolated} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorPositionInterpolated)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    double precision                                                  :: lengthBox

    !![
    <inputParameter>
      <name>lengthBox</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The periodic length of the positions. For non-periodic positions, a value of zero should be given.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nodeOperatorPositionInterpolated(lengthBox,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function positionInterpolatedConstructorParameters

  function positionInterpolatedConstructorInternal(lengthBox,cosmologyFunctions_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily positionInterpolated} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorPositionInterpolated)                        :: self
    double precision                                  , intent(in   )         :: lengthBox
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="lengthBox, *cosmologyFunctions_"/>
    !!]

    self%isPeriodic=lengthBox > 0.0d0
    !![
    <addMetaProperty component="position" name="positionInterpolatedTimeMaximum"  id="self%timeMaximumID"  rank="0" isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="position" name="positionInterpolatedCoefficients" id="self%coefficientsID" rank="1" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function positionInterpolatedConstructorInternal

  subroutine positionInterpolatedAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent, openMPThreadBindingAtLevel
    implicit none
    class(nodeOperatorPositionInterpolated), intent(inout) :: self
    
    call satelliteHostChangeEvent%attach(self,satelliteHostChange,openMPThreadBindingAtLevel,label='positionInterpolated')
    return
  end subroutine positionInterpolatedAutoHook

  subroutine positionInterpolatedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily positionInterpolated} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent
    implicit none
    type(nodeOperatorPositionInterpolated), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    if (satelliteHostChangeEvent%isAttached(self,satelliteHostChange )) call satelliteHostChangeEvent%detach(self,satelliteHostChange)
    return
  end subroutine positionInterpolatedDestructor

  subroutine positionInterpolatedDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition
    implicit none
    class(nodeOperatorPositionInterpolated), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node
    class(nodeComponentPosition           ), pointer       :: position

    position => node%position()
    call position%positionAnalytic()
    call position%velocityAnalytic()
    return
  end subroutine positionInterpolatedDifferentialEvolutionAnalytics

 subroutine positionInterpolatedDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Interrupt evolution to update position interpolation.
    !!}
    use :: Galacticus_Nodes, only : interruptTask, nodeComponentBasic, nodeComponentPosition
    implicit none
    class    (nodeOperatorPositionInterpolated), intent(inout), target  :: self
    type     (treeNode                        ), intent(inout)          :: node
    logical                                    , intent(inout)          :: interrupt
    procedure(interruptTask                   ), intent(inout), pointer :: functionInterrupt
    integer                                    , intent(in   )          :: propertyType
    class    (nodeComponentBasic              )               , pointer :: basic
    class    (nodeComponentPosition           )               , pointer :: position
    !$GLC attributes unused :: propertyType

    basic    => node%basic   ()      
    position => node%position()      
    if (basic%time() >= position%floatRank0MetaPropertyGet(self%timeMaximumID)) then
       interrupt         =  .true.
       self_             => self
       functionInterrupt => positionInterpolatedComputeInterpolation_
    end if
    return
  end subroutine positionInterpolatedDifferentialEvolution

  subroutine positionInterpolatedDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Compute the interpolated position and velocity of the node.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentPosition  , nodeComponentBasic
    use :: Numerical_Interpolation         , only : interpolator
    use :: Histories                       , only : history
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    implicit none
    class           (nodeOperatorPositionInterpolated), intent(inout)                  :: self
    type            (treeNode                        ), intent(inout)                  :: node
    double precision                                  , intent(in   )                  :: time
    class           (nodeComponentBasic              ), pointer                        :: basic                       , basicParent
    class           (nodeComponentPosition           ), pointer                        :: position                    , positionParent
    double precision                                  , dimension(3     )              :: position_                   , velocity_
    double precision                                  , dimension(4,3   )              :: coefficientsCubic
    double precision                                  , dimension(    20)              :: coefficientsSpiral
    double precision                                  , dimension(  2, 2)              :: coefficientsAngle           , coefficientsLogRadius
    double precision                                  , dimension(2,2, 3)              :: vectorInPlaneNormal
    double precision                                  , dimension(:     ), allocatable :: coefficients
    integer                                           , parameter                      :: historyBeginPosition=1      , historyEndPosition   =3, &
         &                                                                                historyBeginVelocity=4      , historyEndVelocity   =6
    double precision                                  , parameter                      :: toleranceTime       =1.0d-2
    double precision                                                                   :: lengthBox                   , angle                  , &
         &                                                                                logRadius                   , timeParent
    integer         (c_size_t                        )                                 :: i
    logical                                                                            :: usingHistory
    type            (interpolator                    )                                 :: interpolator_
    type            (history                         )                                 :: history_

    ! Check that the time is within the allowable range.
    position => node%position()
    basic    => node %basic  ()
    coefficients=position%floatRank1MetaPropertyGet(self%coefficientsID)
    ! If the maximum interpolation time is exceeded, do not attempt to set the position. In this case we are waiting for an
    ! interrupt to differential evolution to trigger a recompute of the interpolation.
    if (size(coefficients) > 0 .and. time > position%floatRank0MetaPropertyGet(self%timeMaximumID)) return
    ! Use the size of the interpolation coefficient array to determine the type of interpolation. A more elegant/robust solution
    ! would be preferable.
    if (size(coefficients) == 0) then
       ! Interpolation coefficients have not yet been computed or are not defined. Use the non-interpolated position.
       usingHistory=.false.
       history_    =position%positionHistory()
       if (history_%exists()) then
          if (history_%time(1) <= basic%time()) then
             interpolator_=interpolator        (history_%time                                                  )
             i            =interpolator_%locate(basic   %time(),closest=.true.                                 )
             position_    =history_     %data  (i              ,        historyBeginPosition:historyEndPosition)
             velocity_    =history_     %data  (i              ,        historyBeginVelocity:historyEndVelocity)
             usingHistory =.true.
          end if
       end if
       if (.not.usingHistory) then
          position_=position%position()
          velocity_=position%velocity()
       end if
    else if (size(coefficients) == 12) then
       ! Using cubic polynomial interpolation.
       coefficientsCubic=reshape(coefficients,[4,3])
       do i=1,3
          position_(i)=+coefficientsCubic(1,i)*time**3 &
               &       +coefficientsCubic(2,i)*time**2 &
               &       +coefficientsCubic(3,i)*time    &
               &       +coefficientsCubic(4,i)
          velocity_(i)=+3.0d0*coefficientsCubic(1,i)*time**2 &
               &      +2.0d0*coefficientsCubic(2,i)*time    &
               &      +      coefficientsCubic(3,i)
       end do
       ! Convert from comoving back to physical position/velocity, and to km/s.
       position_=+position_                                      &
            &    *self%cosmologyFunctions_%expansionFactor(time)
       velocity_=+velocity_                                      &
            &    *self%cosmologyFunctions_%expansionFactor(time) &
            &    *Mpc_per_km_per_s_To_Gyr 
    else if (size(coefficients) == 20) then
       ! Use logarithmic spiral interoplation in physical position.
       coefficientsSpiral   = coefficients
       vectorInPlaneNormal  = reshape(coefficientsSpiral( 1:12),[2,2,3])
       coefficientsAngle    = reshape(coefficientsSpiral(13:16),[2,2  ])
       coefficientsLogRadius= reshape(coefficientsSpiral(17:20),[2,2  ])
       angle                = coefficientsAngle    (1,1)+coefficientsAngle    (1,2)*time
       logRadius            = coefficientsLogRadius(1,1)+coefficientsLogRadius(1,2)*time
       position_            =+(                                       &
            &                  +vectorInPlaneNormal(1,1,:)*cos(angle) &
            &                  +vectorInPlaneNormal(1,2,:)*sin(angle) &
            &                 )                                       &
            &                *exp(logRadius)
       velocity_            =+(                                       &
            &                  +vectorInPlaneNormal(2,1,:)*cos(angle) &
            &                  +vectorInPlaneNormal(2,2,:)*sin(angle) &
            &                 )                                       &
            &                *exp(logRadius)
       ! Our interpolation is relative to the host halo center - add that on now.
       basicParent    => node       %parent%basic   ()
       positionParent => node       %parent%position()
       timeParent     =  basicParent%time           ()
       call basicParent%timeSet(time      )
       position_=+               position_   &
            &    +positionParent%position ()
       velocity_=+               velocity_   &
            &    +positionParent%velocity ()
       call basicParent%timeSet(timeParent)
    end if
    ! Handle periodic boundaries.
    if (self%isPeriodic) then
       lengthBox=self%lengthBox*self%cosmologyFunctions_%expansionFactor(time)
       do i=1,3
          do while (position_(i) <  0.0d0    )
             position_(i)=position_(i)+lengthBox
          end do
          do while (position_(i) >= lengthBox)
             position_(i)=position_(i)-lengthBox
          end do
       end do
    end if
    return
  end subroutine positionInterpolatedDifferentialEvolutionSolveAnalytics

  subroutine positionInterpolatedDifferentialEvolutionPost(self,node)
    !!{
    Trigger interpolation recalculation after an evolution step.
    !!}
    implicit none
    class(nodeOperatorPositionInterpolated), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node

    call self%computeInterpolation(node)
    return
  end subroutine positionInterpolatedDifferentialEvolutionPost

  subroutine positionInterpolatedNodesMerge(self,node)
    !!{
    Trigger interpolation recalculation after a node merger.
    !!}
    implicit none
    class(nodeOperatorPositionInterpolated), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node

    call self%computeInterpolation(node)
    return
  end subroutine positionInterpolatedNodesMerge

  subroutine positionInterpolatedNodeInitialize(self,node)
    !!{
    Trigger interpolation calculation at node initialization
    !!}
    implicit none
    class(nodeOperatorPositionInterpolated), intent(inout)         :: self
    type (treeNode                        ), intent(inout), target :: node

    call self%computeInterpolation(node)
    return
  end subroutine positionInterpolatedNodeInitialize

  subroutine satelliteHostChange(self,node)
    !!{
    Handle cases where a satellite switches host node.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node

    select type (self)
    class is (nodeOperatorPositionInterpolated)
       call self%computeInterpolation(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteHostChange

  subroutine positionInterpolatedNodePromote(self,node)
    !!{
    Promote the node to its parent by copying position interpolation information from the parent node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition
    implicit none
    class(nodeOperatorPositionInterpolated), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node
    class(nodeComponentPosition           ), pointer       :: positionParent, position
    
    position       => node       %position()
    positionParent => node%parent%position()
    call position%floatRank0MetaPropertySet(self% timeMaximumID,positionParent%floatRank0MetaPropertyGet(self% timeMaximumID))
    call position%floatRank1MetaPropertySet(self%coefficientsID,positionParent%floatRank1MetaPropertyGet(self%coefficientsID))
    return
  end subroutine positionInterpolatedNodePromote

  subroutine positionInterpolatedComputeInterpolation_(node)
    !!{
    Interrupt function to recompute interpolation.
    !!}
    type(treeNode), intent(inout), target :: node
    
    call self_%computeInterpolation(node)
    return
  end subroutine positionInterpolatedComputeInterpolation_
    
  subroutine positionInterpolatedComputeInterpolation(self,node)
    !!{
    Compute interpolation coefficients for positions. The approach here follows that of \cite{merson_lightcone_2013}. For halos
    which are not orbiting within a host halo during the current interval, interpolation is via a cubic polynomial in each
    Cartesian coordinate of \emph{comoving} position matched to the position and velocity at the initial and final times. For
    halos which are orbiting within a halo halo during the current interval, interpolation is via a logarithmic spiral in the
    relative physical position of the halo and the center of the host halo.
    !!}
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Error                           , only : Error_Report    
    use            :: Galacticus_Nodes                , only : nodeComponentBasic       , nodeComponentPosition  , treeNode                          , nodeEvent                   , &
         &                                                     nodeEventSubhaloPromotion, nodeEventBranchJump    , nodeEventSubhaloPromotionIntertree, nodeEventBranchJumpIntertree
    use            :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    use            :: Linear_Algebra                  , only : vector                   , matrix                 , assignment(=)
    use            :: Histories                       , only : history
    use            :: Numerical_Comparison            , only : Values_Differ            , Values_Agree           , Values_Less_Than
    use            :: Numerical_Interpolation         , only : interpolator
    use            :: Vectors                         , only : Vector_Product           , Vector_Magnitude
    implicit none
    class(nodeOperatorPositionInterpolated), intent(inout)     :: self
    type (treeNode                        ), intent(inout)     :: node
    class           (nodeComponentPosition), pointer           :: position                     , positionParent       , &
         &                                                        positionGrandParent
    class           (nodeComponentBasic   ), pointer           :: basic                        , basicParent          , &
         &                                                        basicGrandParent
    type            (treeNode             ), pointer           :: nodeParent                   , nodeGrandparent      , &
         &                                                        nodeJump
    class           (nodeEvent            ), pointer           :: event
    double precision                       , dimension(  2, 3) :: positionComoving             , velocityComoving     , &
         &                                                        positionRelative
    double precision                       , dimension(2,2, 3) :: vectorInPlaneNormal
    double precision                       , dimension(  2   ) :: time                         , expansionFactor
    double precision                       , dimension(  2, 2) :: coefficientsAngle            , coefficientsLogRadius
    double precision                       , dimension(     3) :: vectorNormal
    double precision                       , dimension(  4, 3) :: coefficientsCubic
    double precision                       , dimension(    20) :: coefficientsSpiral
    double precision                       , dimension(     0) :: coefficientsNull
    type            (vector               ), allocatable       :: coordinates                  , coefficients
    type            (matrix               ), allocatable       :: terms
    double precision                       , parameter         :: separationTiny        =1.0d-6
    type            (history              )                    :: positionHistory
    integer         (c_size_t             )                    :: i                            , j                    , &
         &                                                        k
    type            (interpolator         )                    :: interpolator_
    double precision                                           :: timeParent                   , timeGrandparent
    logical                                                    :: useSpiralInterpolation       , searchEvent

    ! Get the position component of this node in which we will store the interpolation coefficients.
    position => node%position()
    ! If the node has no parent, there is no interpolation to do, so return.
    if (.not.associated(node%parent)) then
       ! Set a null interpolation - this will caused the fixed-at-snapshot position to be used instead.
       call position%floatRank1MetaPropertySet(self%coefficientsID,coefficientsNull)
       return
    end if
    ! If the node is a satellite, exists at, or after its parent, and the parent has no parent, there is no interpolation to do, so return.
    if (node%isSatellite()) then
       basic       => node       %basic()
       basicParent => node%parent%basic()
       if (Values_Less_Than(basicParent%time(),basic%time(),relTol=1.0d-2)) then
          nodeParent => node%parent
          do while (nodeParent%isSatellite())
             nodeParent => nodeParent%parent
          end do
          if (.not.associated(nodeParent%parent)) then
             ! Set a null interpolation - this will caused the fixed-at-snapshot position to be used instead.
             call position%floatRank1MetaPropertySet(self%coefficientsID,coefficientsNull)
             return
          end if
       end if
    end if
    ! Determine whether to use cubic polynomial or logarithmic spiral interpolation.
    positionHistory =  position%positionHistory()
    basic           => node    %basic          ()
    nodeJump        => null                    ()
    ! First check if we have a position history.
    if (positionHistory%exists()) then
       ! A position history exists. Determine if this node is already a satellite in a host halo - if it is then (by default) we
       ! will use logarithmic spiral interpolation. If the node is not yet a subhalo in a host, then we use cubic polynomial
       ! interpolation instead.
       useSpiralInterpolation=node%isSatellite()
       ! This is to decide if we should look for an event - in cases where we're still utilizing position history, but there's no
       ! grandparent that we can trace to the time where we need it.
       searchEvent  =.false.
       interpolator_=interpolator        (positionHistory%time  )
       i            =interpolator_%locate(basic          %time())
       if(basic%time() >= positionHistory%time(i) .and. size(positionHistory%time) > i) then
          nodeParent  => node      %parent
          basicParent => nodeParent%basic ()
          nodeGrandparent => nodeParent%parent
          if (associated(nodeGrandparent)) then
             basicGrandparent => nodeGrandparent%basic()
             do while (Values_Less_Than(basicGrandparent%time(),positionHistory%time(2),relTol=1.0d-2))
                if (.not.associated(nodeGrandparent%parent)) then
                   searchEvent=.true.
                   exit
                end if
                nodeGrandparent  => nodeGrandparent%parent
                basicGrandparent => nodeGrandparent%basic ()
             end do
          end if
       end if
       ! Check if the current time of the node matches or exceeds the final time in its position history.
       if   (                                                                                 &
          &                  basic%time() >= positionHistory%time(size(positionHistory%time)) &
          &  .or.                                                                             &
          &   (                                                                               &
          &     Values_Agree(basic%time(),positionHistory%time(1),relTol=1.0d-2)              &
          &    .and.                                                                          &
          &     size        (             positionHistory%time                 ) == 1         &
          &   )                                                                               &
          &  .or.                                                                             &
          &   searchEvent                                                                     &
          & ) then
          ! Our node exists at a time beyond the end of its position history. We must check if there is some event attached to
          ! this node which connects it to another node (which we can then use for the future position).
          event => node%event
          do while (associated(event))
             ! Look for a handled event type. Subhalo promotions and branch jumps are handled.
             select type (event)
             type is (nodeEventSubhaloPromotion)
                nodeJump => event%node
                exit
             type is (nodeEventBranchJump      )
                nodeJump => event%node
                exit
             type is (nodeEventSubhaloPromotionInterTree)
                call Error_Report('inter-tree subhalo promotions are not supported'//{introspection:location})
             type is (nodeEventBranchJumpInterTree)
                call Error_Report('inter-tree branch jumps are not supported'      //{introspection:location})
             end select
             event => event%next
          end do
          if (associated(nodeJump)) then
             ! An event which causes our node to jump to some other non-hosted node or another branch was found. We assume cubic
             ! polynomial interpolation between the initial and final locations.
             useSpiralInterpolation=.false.
          else
             ! No history remains, and no event exists. This is the end of the life of this node, so we do not need to compute any
             ! interpolation.
             call position%floatRank1MetaPropertySet(self%coefficientsID,coefficientsNull)
             call position%floatRank0MetaPropertySet(self%timeMaximumID,huge(0.0d0)      )
             return
          end if
       else
          ! The current time of the node does not exceed the time in its position history. We can use logarithmic spiral
          ! interpolation if: a) this is a satellite node, and b) the current time matches or exceeds the first time in the node's
          ! position history.
          useSpiralInterpolation=useSpiralInterpolation .and. basic%time() >= positionHistory%time(1)
          ! Check if satellite has been orphanized into a host which does not exist at the time of the position history. If it
          ! has, we can not use logarithmic spiral interpolation as we do not have a position for the host halo.
          nodeParent  => node      %parent       
          basicParent => nodeParent%basic ()
          if (node%isSatellite() .and. Values_Less_Than(basic%time(),basicParent%time(),relTol=1.0d-2)) useSpiralInterpolation=.false.
       end if
    else
       ! No position history is available - this must be an isolated node and so we must use cubic polynomial interpolation.
       useSpiralInterpolation=.false.
    end if
    ! Nullify the position interpolation for this node before attempting to compute the new interpolation. This ensures that we do
    ! not attempt to interpolate the position of this node outside the allowed range of the previous interpolation.
    call position%floatRank1MetaPropertySet(self%coefficientsID,coefficientsNull)
    ! Branch on the type of interpolation we are to compute.
    if (.not.useSpiralInterpolation) then
       ! Cubic polynomial interpolation.
       !! Determine the intial and final position to use for this time interval.
       if (associated(nodeJump)) then
          ! There is a node event which will jump this node to another node. Use the position of that node as the final position.
          basicParent              => nodeJump             %basic   (                     )
          positionParent           => nodeJump             %position(                     )
          positionComoving   (1,:) =  position             %position(                     )
          velocityComoving   (1,:) =  position             %velocity(                     )
          time               (1  ) =  basic                %time    (                     )
          time               (2  ) =  basicParent          %time    (                     )
          positionComoving   (2,:) =  positionParent       %position(                     )
          velocityComoving   (2,:) =  positionParent       %velocity(                     )
       else if (node%isPrimaryProgenitor()) then
          ! The node is the primary progenitor, so use the position of its parent as the final position.
          basicParent              => node          %parent%basic   (                     )
          positionParent           => node          %parent%position(                     )
          time               (1  ) =  basic                %time    (                     )
          time               (2  ) =  basicParent          %time    (                     )
          positionComoving   (1,:) =  position             %position(                     )
          positionComoving   (2,:) =  positionParent       %position(                     )
          velocityComoving   (1,:) =  position             %velocity(                     )
          velocityComoving   (2,:) =  positionParent       %velocity(                     )
       else if (positionHistory%exists()) then
          ! This node will become a subhalo at the end of this time interval - use its position history to determine the final
          ! location.
          if (Values_Less_Than(basic%time(),positionHistory%time(1),relTol=1.0d-2)) then
             ! The current time is significantly less than the first time in the position history. Therefore, use the current
             ! time/position as our starting point, and the first time/position in the history and our ending point.
             time            (1  ) =  basic                 %time    (                    )
             time            (2  ) =  positionHistory       %time    (1                   )
             positionComoving(1,:) =  position              %position(                    )
             positionComoving(2,:) =  positionHistory       %data    (1  ,1:3             )
             velocityComoving(1,:) =  position              %velocity(                    )
             velocityComoving(2,:) =  positionHistory       %data    (1  ,4:6             )
          else
             ! The current time is at or after the first time in the position history. Therefore, locate the point in the position
             ! history that spans the current time.
             interpolator_         =  interpolator                 (positionHistory%time  )
             i                     =  interpolator_         %locate(basic          %time())
             time            (1  ) =  positionHistory       %time  (i                     )
             time            (2  ) =  positionHistory       %time  (i+1                   )
             positionComoving(1,:) =  positionHistory       %data  (i  ,1:3               )
             positionComoving(2,:) =  positionHistory       %data  (i+1,1:3               )
             velocityComoving(1,:) =  positionHistory       %data  (i  ,4:6               )
             velocityComoving(2,:) =  positionHistory       %data  (i+1,4:6               )

          end if
       else
          ! We have no means to determine the start/end points for this interval. This node must have no future, and so we can
          ! safely set a null interpolation infinitely far into the future.
          call position%floatRank1MetaPropertySet(self%coefficientsID,coefficientsNull)
          call position%floatRank0MetaPropertySet(self%timeMaximumID ,huge(0.0d0)     )
          return
       end if
       ! Compute the interpolation only if there is a non-zero time interval over which to interpolate. Halos about to experience
       ! a node promotion or branch jump event can exist at precisely the time of the halo they are about to become.
       if (time(1) < time(2)) then
          ! Get comoving position/velocity at start/end times.
          expansionFactor (1  )=                      self%cosmologyFunctions_%expansionFactor(time(1))
          expansionFactor (2  )=                      self%cosmologyFunctions_%expansionFactor(time(2))
          positionComoving(1,:)=positionComoving(1,:)/                         expansionFactor(     1 )
          positionComoving(2,:)=positionComoving(2,:)/                         expansionFactor(     2 )
          velocityComoving(1,:)=velocityComoving(1,:)/                         expansionFactor(     1 )/Mpc_per_km_per_s_To_Gyr
          velocityComoving(2,:)=velocityComoving(2,:)/                         expansionFactor(     2 )/Mpc_per_km_per_s_To_Gyr
          ! Handle periodic positions.
          if (self%isPeriodic) then
             do j=1,3
                if (positionComoving(2,j) > positionComoving(1,j)+0.5d0*self%lengthBox) positionComoving(2,j)=positionComoving(2,j)-self%lengthBox
                if (positionComoving(2,j) < positionComoving(1,j)-0.5d0*self%lengthBox) positionComoving(2,j)=positionComoving(2,j)+self%lengthBox
             end do
          end if
          ! Solve for the interpolation coefficients in each Cartesian axis.
          do i=1,3
             allocate(terms       )
             allocate(coordinates )
             allocate(coefficients)
             coordinates=vector(                                                                                                             &
                  &                               [                                                                                          &
                  &                                 positionComoving(1,i),positionComoving(2,i),velocityComoving(1,i),velocityComoving(2,i)  &
                  &                               ]                                                                                          &
                  &             )
             terms      =matrix(                                                                                                             &
                  &             transpose(                                                                                                   &
                  &                       reshape(                                                                                           &
                  &                               [                                                                                          &
                  &                                      time(1)**3      ,      time(1)**2     ,time(1)              ,1.0d0                , &
                  &                                      time(2)**3      ,      time(2)**2     ,time(2)              ,1.0d0                , &
                  &                                3.0d0*time(1)**2      ,2.0d0*time(1)        ,1.0d0                ,0.0d0                , &
                  &                                3.0d0*time(2)**2      ,2.0d0*time(2)        ,1.0d0                ,0.0d0                  &
                  &                               ]                                                                                        , &
                  &                               [4,4]                                                                                      &
                  &                              )                                                                                           &
                  &                      )                                                                                                   &
                  &            )
             coefficients          =terms%linearSystemSolve(coordinates )
             coefficientsCubic(:,i)=                        coefficients
             deallocate(terms       )
             deallocate(coordinates )
             deallocate(coefficients)
          end do
          ! Store the computed interpolation coefficients and the time to which they are valid.
          call position%floatRank1MetaPropertySet(self%coefficientsID,reshape(coefficientsCubic,[12]))
          call position%floatRank0MetaPropertySet(self%timeMaximumID ,   time(                    2 ))
       end if
    else
       ! Compute logarithmic spiral interpolation.
       basic           => node    %basic          ()
       position        => node    %position       ()
       positionHistory =  position%positionHistory()
       if (.not.positionHistory%exists()) call Error_Report('position history does not exist'//{introspection:location})
       ! Find the interval in the position history which spans the current epoch.
       interpolator_   =interpolator          (positionHistory%time  )
       i               =interpolator_  %locate(basic          %time())
       time         (1)=positionHistory%time  (i                     )
       time         (2)=positionHistory%time  (i+1                   )
       ! We want to find a parent (and grandparent) whose existance spans these times.
       nodeParent  => node      %parent       
       basicParent => nodeParent%basic ()
       do while (Values_Less_Than(basicParent%time(),time(1),relTol=1.0d-2))
          nodeParent  => nodeParent%parent
          basicParent => nodeParent%basic ()
       end do
       ! If the parent is not (approximately) at our initial time, back up to the prior parent (if one exists).
       if (Values_Differ(basicParent%time(),time(1),relTol=1.0d-2) .and. associated(nodeParent%firstChild)) then
          nodeParent  => nodeParent%firstChild
          basicParent => nodeParent%basic     ()
       end if
       ! Find a grand-parent which exists at (or after) the final time.
       nodeGrandparent => nodeParent%parent
       if (associated(nodeGrandparent)) then
          basicGrandparent => nodeGrandparent%basic()
          do while (Values_Less_Than(basicGrandparent%time(),time(2),relTol=1.0d-2))
             nodeGrandparent => nodeGrandparent%parent
             if (associated(nodeGrandparent)) then
                basicGrandparent => nodeGrandparent%basic()
             else
                ! No grandparent halo exists at the final time. Use a null interpolation to ensure that
                ! the fixed-at-snapshot position will be used for the node. Set the maximum time for the interpolation to infinity to
                ! avoid an infinite loop of attempting to compute this interpolation.
                call position%floatRank1MetaPropertySet(self%coefficientsID,coefficientsNull)
                call position%floatRank0MetaPropertySet(self%timeMaximumID ,huge(0.0d0)     )
                return
             end if
          end do
          ! Back up to the previous grandparent. This grandparent therefore exists up until our final time.
          nodeGrandparent  => nodeGrandparent%firstChild
          basicGrandparent => nodeGrandparent%basic     ()
          ! Calculate interpolating factors for positions and velocities.
          positionParent      => nodeParent      %position()
          positionGrandparent => nodeGrandparent %position()
          timeParent          =  basicParent     %time    ()
          timeGrandparent     =  basicGrandparent%time    ()
          ! Iterate over position and velocity.
          do k=1,2
             ! Find displacement vectors (in physical coordinates) from the host center. Adjust the parent/grandparent nodes to
             ! the initial/final times of our interval such that their positions will be interpolated to the correct time.
             select case (k)
             case (1)
                ! Position.
                call basicParent     %timeSet(time(1))
                positionRelative(1,:)=positionHistory%data(i  ,1:3)-positionParent     %position()
                call basicGrandparent%timeSet(time(2))
                positionRelative(2,:)=positionHistory%data(i+1,1:3)-positionGrandParent%position()
             case (2)
                ! Velocity
                call basicParent     %timeSet(time(1))
                positionRelative(1,:)=positionHistory%data(i  ,4:6)-positionParent     %velocity()
                call basicGrandparent%timeSet(time(2))
                positionRelative(2,:)=positionHistory%data(i+1,4:6)-positionGrandParent%velocity()
             end select
             ! Handle periodic positions.
             if (k == 1 .and. self%isPeriodic) then
                expansionFactor(1)=self%cosmologyFunctions_%expansionFactor(time(1))
                expansionFactor(2)=self%cosmologyFunctions_%expansionFactor(time(2))
                do j=1,3
                   if (positionRelative(2,j)/expansionFactor(2) > positionRelative(1,j)/expansionFactor(1)+0.5d0*self%lengthBox) positionRelative(2,j)=positionRelative(2,j)-self%lengthBox*expansionFactor(2)
                   if (positionRelative(2,j)/expansionFactor(2) < positionRelative(1,j)/expansionFactor(1)-0.5d0*self%lengthBox) positionRelative(2,j)=positionRelative(2,j)+self%lengthBox*expansionFactor(2)
                end do
             end if
             ! Handle zero relative positions.
             if (all(positionRelative(2,:) == 0.0d0)) then
                ! Final relative position is zero.
                if (all(positionRelative(1,:) == 0.0d0)) then
                   ! Initial position is zero also, use an arbitrary direction.
                   positionRelative(2,:)=separationTiny*[1.0d0,0.0d0,0.0d0]
                else
                   positionRelative(2,:)=separationTiny*positionRelative(1,:)/Vector_Magnitude(positionRelative(1,:))
                end if
             end if
             if (all(positionRelative(1,:) == 0.0d0)) then
                ! Initial relative position is zero.
                if (all(positionRelative(2,:) == 0.0d0)) then
                   ! Final position is zero also, use an arbitrary direction.
                   positionRelative(1,:)=separationTiny*[1.0d0,0.0d0,0.0d0]
                else
                   positionRelative(1,:)=separationTiny*positionRelative(2,:)/Vector_Magnitude(positionRelative(2,:))
                end if
             end if
             ! Construct the normal vector to the orbital plane, plus a vector in the plane which is normal to our initial time
             ! relative position.
             vectorNormal              =+Vector_Product  (positionRelative   (  1,:),positionRelative(2,:))
             vectorNormal              =+                 vectorNormal                                      &
                  &                     /Vector_Magnitude(vectorNormal                                    )
             vectorInPlaneNormal(k,1,:)=+                 positionRelative   (  1,:)                        &
                  &                     /Vector_Magnitude(positionRelative   (  1,:)                      )
             vectorInPlaneNormal(k,2,:)=+Vector_Product  (vectorInPlaneNormal(k,1,:),vectorNormal         )
             vectorInPlaneNormal(k,2,:)=+                 vectorInPlaneNormal(k,2,:)                        &
                  &                     /Vector_Magnitude(vectorInPlaneNormal(k,2,:)                      )
             !! The in plane normal vector is defined only up to +/- inversion. Choose the option which minimizes the angle between
             !! it and the relative position vector at the final time. This is how we enforce the choice of Merson et al. (2013) to
             !! always assumes that the motion of the particle is over the minimal span of angle.
             if     (                                                                &
                  &   Dot_Product(positionRelative(2,:),-vectorInPlaneNormal(k,2,:)) &
                  &  >                                                               &
                  &   Dot_Product(positionRelative(2,:),+vectorInPlaneNormal(k,2,:)) &
                  & ) vectorInPlaneNormal(k,2,:)=-vectorInPlaneNormal(k,2,:)
             !! Find the linear fit coefficients to angle.
             coefficientsAngle    (k,2)=+acos(                                                                              &
                  &                           min(                                                                          &
                  &                               +1.0d0                                                                  , &
                  &                               max(                                                                      &
                  &                                   -1.0d0                                                              , &
                  &                                   +Dot_Product     (positionRelative(2,:),vectorInPlaneNormal(k,1,:))   &
                  &                                   /Vector_Magnitude(positionRelative(2,:)                           )   &
                  &                                  )                                                                      &
                  &                               )                                                                         &
                  &                          )                                                                              &
                  &                     /    (                                                                              &
                  &                           +time(2)                                                                      &
                  &                           -time(1)                                                                      &
                  &                          )
             coefficientsAngle    (k,1)=-time(1) &
                  &                     *coefficientsAngle(k,2)
             !! Find the linear fit coefficients to log-radius.
             coefficientsLogRadius(k,2)=+(                                              &
                  &                       +log(Vector_Magnitude(positionRelative(2,:))) &
                  &                       -log(Vector_Magnitude(positionRelative(1,:))) &
                  &                      )                                              &
                  &                     /(                                              &
                  &                       +time(2) &
                  &                       -time(1) &
                  &                      )
             coefficientsLogRadius(k,1)=-time(1)                                        &
                  &                     *coefficientsLogRadius(k,2)                     &
                  &                     +log(Vector_Magnitude(positionRelative(1,:)))             
          end do
          ! Reset parent/grandparent node times.
          call basicParent     %timeSet(timeParent     )
          call basicGrandparent%timeSet(timeGrandparent)
          ! Store the computed interpolation coefficients.
          coefficientsSpiral( 1:12)=reshape(vectorInPlaneNormal  ,[12])
          coefficientsSpiral(13:16)=reshape(coefficientsAngle    ,[ 4])
          coefficientsSpiral(17:20)=reshape(coefficientsLogRadius,[ 4])
          call position%floatRank1MetaPropertySet(self%coefficientsID,coefficientsSpiral   )
          call position%floatRank0MetaPropertySet(self%timeMaximumID ,time              (2))
       else
          ! No grandparent halo exists. This can occur in a tree which ceases to exist. Use a null interpolation to ensure that
          ! the fixed-at-snapshot position will be used for the node. Set the maximum time for the interpolation to infinity to
          ! avoid an infinite loop of attempting to compute this interpolation.
          call position%floatRank1MetaPropertySet(self%coefficientsID,coefficientsNull)
          call position%floatRank0MetaPropertySet(self%timeMaximumID ,huge(0.0d0)     )
       end if
    end if
    return
  end subroutine positionInterpolatedComputeInterpolation
