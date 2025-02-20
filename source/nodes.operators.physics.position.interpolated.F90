!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  use :: Galacticus_Nodes   , only : treeNode
  
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
     logical                                            :: isPeriodic                   , wrapPeriodic
     integer                                            :: coefficientsID
   contains
     !![
     <methods>
       <method method="interpolate" description="Perform the interpolation."/>
     </methods>
     !!]
     final     ::                                        positionInterpolatedDestructor
     procedure :: differentialEvolutionAnalytics      => positionInterpolatedDifferentialEvolutionAnalytics
     procedure :: predeterminedSolveAnalytics         => positionInterpolatedPredeterminedSolveAnalytics
     procedure :: differentialEvolutionSolveAnalytics => positionInterpolatedDifferentialEvolutionSolveAnalytics
     procedure :: nodeTreeInitialize                  => positionInterpolatedNodeTreeInitialize
     procedure :: interpolate                         => positionInterpolatedInterpolate
  end type nodeOperatorPositionInterpolated
  
  interface nodeOperatorPositionInterpolated
     !!{
     Constructors for the {\normalfont \ttfamily positionInterpolated} node operator class.
     !!}
     module procedure positionInterpolatedConstructorParameters
     module procedure positionInterpolatedConstructorInternal
  end interface nodeOperatorPositionInterpolated

  type :: nodeTrace
     !!{
     Type used for tracing nodes through their positional history. 
     !!}
     double precision                     :: time
     type            (treeNode ), pointer :: node        => null(), nodeHost => null()
     type            (nodeTrace), pointer :: next        => null()
     logical                              :: isSatellite
     integer         (c_size_t )          :: iHistory
  end type nodeTrace

  ! Count of coefficients for different interpolation types.
  integer(c_size_t), parameter :: countCoefficientsCubicPolynomial  =12_c_size_t
  integer(c_size_t), parameter :: countCoefficientsLogarithmicSpiral=20_c_size_t
  
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
    logical                                                           :: wrapPeriodic

    !![
    <inputParameter>
      <name>lengthBox</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The periodic length of the positions. For non-periodic positions, a value of zero should be given.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wrapPeriodic</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, interpolated positions that lie outside of the periodic box will be wrapped back into the box.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nodeOperatorPositionInterpolated(lengthBox,wrapPeriodic,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function positionInterpolatedConstructorParameters

  function positionInterpolatedConstructorInternal(lengthBox,wrapPeriodic,cosmologyFunctions_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily positionInterpolated} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorPositionInterpolated)                        :: self
    double precision                                  , intent(in   )         :: lengthBox
    logical                                           , intent(in   )         :: wrapPeriodic
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="lengthBox, wrapPeriodic, *cosmologyFunctions_"/>
    !!]

    self%isPeriodic=lengthBox > 0.0d0
    !![
    <addMetaProperty component="position" name="positionInterpolatedCoefficients" id="self%coefficientsID" rank="1" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function positionInterpolatedConstructorInternal

  subroutine positionInterpolatedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily positionInterpolated} node operator class.
    !!}
    implicit none
    type(nodeOperatorPositionInterpolated), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
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

  subroutine positionInterpolatedDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Compute the interpolated position and velocity of the node.
    !!}
    implicit none
    class           (nodeOperatorPositionInterpolated), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    double precision                                  , intent(in   ) :: time
   
    ! Call the function to perform the interpolation.
    call self%interpolate(node,time)
    return
  end subroutine positionInterpolatedDifferentialEvolutionSolveAnalytics
  
  subroutine positionInterpolatedPredeterminedSolveAnalytics(self,node,time)
    !!{
    Compute the interpolated position and velocity of the node.
    !!}
    implicit none
    class           (nodeOperatorPositionInterpolated), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    double precision                                  , intent(in   ) :: time
    
    ! Call the function to perform the interpolation.
    call self%interpolate(node,time)
    return
  end subroutine positionInterpolatedPredeterminedSolveAnalytics
  
  subroutine positionInterpolatedInterpolate(self,node,time)
    !!{
    Compute the interpolated position and velocity of the node.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentPosition, nodeComponentBasic
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    class           (nodeOperatorPositionInterpolated), intent(inout)                  :: self
    type            (treeNode                        ), intent(inout)                  :: node
    double precision                                  , intent(in   )                  :: time
    class           (nodeComponentPosition           ), pointer                        :: position
    double precision                                  , dimension(3     )              :: position_          , velocity_
    double precision                                  , dimension(4,3   )              :: coefficientsCubic
    double precision                                  , dimension(    20)              :: coefficientsSpiral
    double precision                                  , dimension(  2, 2)              :: coefficientsAngle  , coefficientsLogRadius
    double precision                                  , dimension(2,2, 3)              :: vectorInPlaneNormal
    double precision                                  , dimension(2     )              :: angle              , logRadius
    double precision                                  , dimension(:     ), allocatable :: coefficients
    integer         (c_size_t                        )                                 :: countTrace         , iTrace               , &
         &                                                                                offset
    integer                                                                            :: i
    logical                                                                            :: isSpiral
    double precision                                                                   :: lengthBox
    !$GLC attributes initialized :: coefficients
    
    ! Extract all interpolation coefficients.
    position     => node    %position                 (                   )
    coefficients =  position%floatRank1MetaPropertyGet(self%coefficientsID)
    ! Determine the number of steps in the interpolation.
    countTrace  =size(coefficients)/(2_c_size_t+countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)
    ! Find the appropriate time.
    iTrace=1_c_size_t
    do while (iTrace < countTrace .and. time >= coefficients(iTrace+1_c_size_t))
       iTrace=iTrace+1_c_size_t
    end do
    ! Extract cubic polynomial coefficients for this time.
    isSpiral         =coefficients(countTrace+iTrace) > 0.0d0
    offset           =2_c_size_t*countTrace+(countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)*(iTrace-1_c_size_t)
    coefficientsCubic=reshape(coefficients(offset+1_c_size_t:offset+countCoefficientsCubicPolynomial),[4,3])
    ! Interpolate the cubic polynomial.
    do i=1,3
       position_(i)=+      coefficientsCubic(1,i)*time**3 &
            &       +      coefficientsCubic(2,i)*time**2 &
            &       +      coefficientsCubic(3,i)*time    &
            &       +      coefficientsCubic(4,i)
       velocity_(i)=+3.0d0*coefficientsCubic(1,i)*time**2 &
            &       +2.0d0*coefficientsCubic(2,i)*time    &
            &       +      coefficientsCubic(3,i)
    end do
    ! Convert from comoving back to physical position/velocity, and to km/s.
    position_=+position_                                      &
         &    *self%cosmologyFunctions_%expansionFactor(time)
    velocity_=+velocity_                                      &
         &    *self%cosmologyFunctions_%expansionFactor(time) &
         &    *MpcPerKmPerSToGyr 
    if (isSpiral) then
      ! Add on logarithmic spiral interpolation in physical position for satellite nodes.
      offset               =offset+countCoefficientsCubicPolynomial
      coefficientsSpiral   =coefficients(offset+1_c_size_t:offset+countCoefficientsLogarithmicSpiral)
      vectorInPlaneNormal  =reshape(coefficientsSpiral( 1:12),[2,2,3])
      coefficientsAngle    =reshape(coefficientsSpiral(13:16),[2,2  ])
      coefficientsLogRadius=reshape(coefficientsSpiral(17:20),[2,2  ])
      angle                =coefficientsAngle    (:,1)+coefficientsAngle    (:,2)*time
      logRadius            =coefficientsLogRadius(:,1)+coefficientsLogRadius(:,2)*time
      position_            =+position_                                  &
           &                +(                                          &
           &                  +vectorInPlaneNormal(1,1,:)*cos(angle(1)) &
           &                  +vectorInPlaneNormal(1,2,:)*sin(angle(1)) &
           &                 )                                          &
           &                *exp(logRadius(1))
      velocity_            =+velocity_                                  &
           &                +(                                          &
           &                  +vectorInPlaneNormal(2,1,:)*cos(angle(2)) &
           &                  +vectorInPlaneNormal(2,2,:)*sin(angle(2)) &
           &                 )                                          &
           &                *exp(logRadius(2))
   end if
    ! Handle periodic boundaries.
    if (self%isPeriodic .and. self%wrapPeriodic) then
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
    ! Set the position and velocity.
    call position%positionSet(position_)
    call position%velocitySet(velocity_)
    return
  end subroutine positionInterpolatedInterpolate

  subroutine positionInterpolatedNodeTreeInitialize(self,node)
    !!{
    Trigger interpolation calculation at node initialization
    !!}
    use :: Error                       , only : Error_Report
    use :: Galacticus_Nodes            , only : nodeComponentBasic       , nodeComponentPosition, nodeComponentSatellite            , nodeEvent                   , &
         &                                      nodeEventSubhaloPromotion, nodeEventBranchJump  , nodeEventSubhaloPromotionIntertree, nodeEventBranchJumpIntertree
    use :: Histories                   , only : history
    use :: Satellite_Merging_Timescales, only : satelliteMergeTimeInfinite
    use :: String_Handling             , only : operator(//)
    use :: ISO_Varying_String          , only : var_str
    implicit none
    class           (nodeOperatorPositionInterpolated), intent(inout)  , target      :: self
    type            (treeNode                        ), intent(inout)  , target      :: node
    type            (treeNode                        )                 , pointer     :: nodeDescendant      , nodeHost                  , &
         &                                                                              nodeMergeTarget     , nodeMergeTargetPrevious
    class           (nodeComponentBasic              )                 , pointer     :: basic               , basicHost                 , &
         &                                                                              basicMergeTarget
    class           (nodeComponentPosition           )                 , pointer     :: position            , positionMergeTarget
    class           (nodeComponentSatellite          )                 , pointer     :: satellite
    class           (nodeEvent                       )                 , pointer     :: event, eventPrior
    type            (nodeTrace                       )                 , pointer     :: traceTail           , traceHead 
    double precision                                  , dimension(   3)              :: positionReference
    double precision                                  , dimension(4 ,3)              :: coefficientsCubic
    double precision                                  , dimension(20  )              :: coefficientsSpiral
    double precision                                  , dimension( :  ), allocatable :: coefficients
    character       (len=12                          )                               :: label
    type            (history                         )                               :: positionHistory     , positionHistoryMergeTarget
    double precision                                                                 :: time                , timeMerger                , &
         &                                                                              timeBranchJump      , timeSubhaloPromotion      , &
         &                                                                              timeJumpLatest
    integer         (c_size_t                        )                               :: iHistory            , countHistory              , &
         &                                                                              iTrace              , countTrace                , &
         &                                                                              offset              , i
    logical                                                                          :: isSatellite         , haveHistory               , &
         &                                                                              haveMerger          , haveBranchJump            , &
         &                                                                              haveSubhaloPromotion, isInitialSatellite        , &
         &                                                                              isHost

    ! Consider only branch tips.
    if (associated(node%firstChild)) return
    ! Look for subhalo promotion events. Halos that are apparently at the tip of a branch, but have a subhalo promotion associated
    ! with their initial time are actually the target of the promotion. We do not want to process them here. Instead, the subhalo
    ! that promotes to them will be processed, and will trace through this subhalo promotion.
    haveSubhaloPromotion =  .false.
    event                => node%event
    do while (associated(event))
       ! Look for a handled event type. Subhalo promotions and branch jumps are handled.
       select type (event)
       type is (nodeEventSubhaloPromotion)
          haveSubhaloPromotion=.true.
          timeSubhaloPromotion=event%time
          exit
       end select
       event => event%next
    end do
    if (haveSubhaloPromotion) then
       basic => node%basic()
       ! Check if the time of this promotion event coincides with the time of this node. If it does, we do not process this node.
       if (basic%time() == event%time) return
    end if
    ! Initialize the trace data structure.
    traceHead  => null()
    traceTail  => null()
    countTrace =  0_c_size_t
    ! Find the initial node and host.
    nodeDescendant => node
    nodeHost       => node
    ! If the initial node is a satellite, find its isolated host.
    do while (nodeHost%isSatellite())
       nodeHost => nodeHost%parent
    end do
    isSatellite       =.not.associated(nodeDescendant,nodeHost)
    isInitialSatellite=isSatellite
    ! Set the initial time.
    basic => nodeDescendant%basic()
    time  =  basic         %time ()
    ! Initialize history status.
    iHistory    =-1_c_size_t
    countHistory=+0_c_size_t
    ! Initialize merger status.
    haveMerger      =  .false.
    timeMerger      =  huge(0.0d0)
    nodeMergeTarget => null()    
    ! Follow the node through the tree.
    do while (associated(nodeDescendant))
       ! Determine if any position history is available for the current descendant.
       basic           => nodeDescendant %basic          ()
       basicHost       => nodeHost       %basic          ()
       position        => nodeDescendant %position       ()
       positionHistory =  position       %positionHistory()
       haveHistory     =  positionHistory%exists         ()
       if (haveHistory) then
          countHistory=size(positionHistory%time)
          if (isInitialSatellite) then
             ! For nodes that are initially satellites and have a position history, start them on the first step of that history.
             iHistory          =+1_c_size_t
          end if
       else
          iHistory    =-1_c_size_t
          countHistory=+0_c_size_t
       end if
       ! Find any future events associated with this descendant.
       call findEvents(presentOnly=.false.)
       ! Determine any merger time for this node.
       satellite => nodeDescendant%satellite()
       if (.not.haveMerger .and. satellite%timeOfMerging() < satelliteMergeTimeInfinite) then
          haveMerger=.true.
          timeMerger=satellite%timeOfMerging()
          ! Use any cloned copy of this node (if available) to determine the merge target.
          if (associated(nodeDescendant%parent) .and. nodeDescendant%parent%index() == nodeDescendant%index()) then
             nodeMergeTarget => nodeDescendant%parent%mergeTarget
          else
             nodeMergeTarget => nodeDescendant%mergeTarget
          end if
          ! If we have a merge target available, follow it. 
          if (associated(nodeMergeTarget)) then
             ! Look for a clone - if found, use the clone as our target as it may have a position history.
             if     (                                                           &
                  &   associated(nodeMergeTarget%parent)                        &
                  &  .and.                                                      &
                  &   nodeMergeTarget%parent%index() == nodeMergeTarget%index() &
                  & )                                                           &
                  & nodeMergeTarget => nodeMergeTarget%parent
             basicMergeTarget           => nodeMergeTarget    %basic          ()
             positionMergeTarget        => nodeMergeTarget    %position       ()
             positionHistoryMergeTarget =  positionMergeTarget%positionHistory()
             if (positionHistoryMergeTarget%exists()) then
                if      (timeMerger > positionHistoryMergeTarget%time(size(positionHistoryMergeTarget%time))) then
                   timeMerger=positionHistoryMergeTarget%time(size(positionHistoryMergeTarget%time))
                else if (timeMerger < positionHistoryMergeTarget%time(1                                    )) then
                   timeMerger=basicMergeTarget          %time(                                     )
                else
                   i=size(positionHistoryMergeTarget%time)
                   do while (positionHistoryMergeTarget%time(i) > timeMerger)
                      i=i-1
                   end do
                   timeMerger=positionHistoryMergeTarget%time(i)
                end if
             else
                if (timeMerger > basicMergeTarget%time()) timeMerger=basicMergeTarget%time()
             end if
          else
             ! Trace the host here to find the descendant just prior to the merge time. We only look at hosts, so no need to
             ! consider branch jumps or subhalo promotions.
             nodeMergeTarget         => nodeHost
             basicMergeTarget        => nodeMergeTarget%basic()
             nodeMergeTargetPrevious => null()
             do while (basicMergeTarget%time() <= timeMerger .and. associated(nodeMergeTarget%parent))
                nodeMergeTargetPrevious => nodeMergeTarget
                nodeMergeTarget         => nodeMergeTarget%parent
                basicMergeTarget        => nodeMergeTarget%basic ()
             end do
             ! Back up to the prior node if the current one exceeds the merger time - we need node positions to coincide from the
             ! time at which they merge (or earlier).
             if (basicMergeTarget%time() > timeMerger) then
                nodeMergeTarget  => nodeMergeTargetPrevious
                basicMergeTarget => nodeMergeTarget        %basic()
             end if
             timeMerger=basicMergeTarget%time()
          end if
       end if
       ! Orphanize an initial satellite that has no history.
       if (isInitialSatellite.and..not.haveHistory) &
            & nodeDescendant => nodeHost
       isInitialSatellite=.false.
       ! Determine satellite status for the current descendant.
       isSatellite=.not.associated(nodeDescendant,nodeHost)
       if (.not.isSatellite) iHistory=-1_c_size_t
       ! Append this step to the trace.
       countTrace=countTrace+1_c_size_t
       if (associated(traceTail)) then
          allocate(traceTail%next)
          traceTail => traceTail%next
       else
          allocate(traceHead     )
          traceTail => traceHead
       end if
       traceTail%time        =  time
       traceTail%isSatellite =  isSatellite
       traceTail%iHistory    =  iHistory
       traceTail%node        => nodeDescendant
       traceTail%nodeHost    => nodeHost
       traceTail%next        => null()
       ! Move to the next step.
       !! Case where we want to move to the parent node.
       if     (                                                        &
            &                    nodeDescendant%isPrimaryProgenitor()  & ! 1. Node is primary progenitor - so we can move directly to its descendant.
            &  .or.                                                    &
            &    .not.associated(nodeDescendant%parent               ) & ! 2. Node has no descendant - we have reached the end of the tree - move to that null
            &  .or.                                                    & !    descendant to terminate the tree walk.
            &    .not.haveHistory                                      & ! 3. Node is not the primary progenitor, but has no position history. It becomes an orphan
            &  .or.                                                    & !    in its parent, so move our pointer to that parent such that we will follow its position.
            &   (     haveHistory .and. iHistory == countHistory)      & ! 4. We have reached the end of the node's position history.  It becomes an orphan
            & ) then                                                     !    in its parent, so move our pointer to that parent such that we will follow its position.
          !! Check for a subhalo promotion.
          if (haveSubhaloPromotion) then
             ! We have a subhalo promotion - move our pointer to that node.
             nodeDescendant => event         %node
             nodeHost       => nodeDescendant        ! This is always an isolated halo by construction, so must be self-hosting.
          else
             nodeDescendant => nodeHost      %parent ! Always move to the host's parent - this allows us to correctly handle satellites that have been orphanized.
             nodeHost       => nodeDescendant        ! This is always an isolated halo by construction, so must be self-hosting.
          end if
          ! Update the time if the descendant exists.
          if (associated(nodeDescendant)) then
             basic => nodeDescendant%basic()
             time  =  basic         %time ()
          end if
       !! Case where we want to follow a position history.
       else if (haveHistory) then
          ! We have a position history to use.
          if (iHistory < 0_c_size_t) then
             ! This is our first use of the history - find the initial step in history that we want to use.
             iHistory=1_c_size_t
             do while (iHistory < countHistory .and. time >= positionHistory%time(iHistory))
                iHistory=iHistory+1_c_size_t
             end do
          else
             ! This is a subsequent step in the history, simply move to the next step.
             iHistory=iHistory+1_c_size_t
          end if
          ! Check if we have not yet reached the end of history.
          if (iHistory <= countHistory) then
             ! Update the time to this new history step.
             time=positionHistory%time(iHistory)
             ! Move the host pointer.
             basicHost => nodeHost%basic()
             do while (associated(nodeHost) .and. basicHost%time() < time)
                nodeHost => nodeHost%parent
                if (associated(nodeHost)) basicHost => nodeHost%basic()
             end do
             ! If the host exists after our current descendant, back up to the previous host.
             if (associated(nodeHost) .and. basicHost%time() > time) nodeHost => nodeHost%firstChild
          else
             ! The end of the history was reached.
             !! Check for a subhalo promotion.
             if (haveSubhaloPromotion) then
                ! We have a subhalo promotion - move our pointer to that node.
                nodeDescendant => event         %node
                nodeHost       => nodeDescendant      ! This is always an isolated halo by construction, so must be self-hosting.
             else
                ! This descendant becomes orphanized in its parent, so move our pointer to that parent such that we will follow
                ! its position.
                nodeDescendant => nodeHost      %parent
                nodeHost       => nodeDescendant        ! This is always an isolated halo by construction, so must be self-hosting.
             end if
             ! Update the time if the descendant exists.
             if (associated(nodeDescendant)) then
                basic => nodeDescendant%basic()
                time  =  basic         %time ()
             end if
          end if

       end if
       ! Check for a merger.
       if (haveMerger .and. time >= timeMerger) then
          ! Move to the host, unless the merger time is before the current time, in which case we need to back up to the host's
          ! progenitor.
          if (associated(nodeMergeTarget)) then
             ! A merge target is available - move to that node.
             nodeDescendant   => nodeMergeTarget
             ! Determine if this descendant is its own host.
             !! If it is the primary progenitor, then it must be its own host.
             isHost=nodeDescendant%isPrimaryProgenitor().or..not.associated(nodeDescendant%parent)
             !! Otherwise, it is its own host if the merger happens prior to it becoming a satellite in its parent.
             if (.not.isHost) then
                basicHost => nodeDescendant%parent%basic()
                isHost    =  timeMerger < basicHost%time()
             end if
             ! If we are merging onto some subhalo, but we merge at a time after it is orphanized - then we need to move to its host instead
             if (isHost) then
                nodeHost => nodeDescendant
             else
                ! The descendant is not its own host - trace down until we find the host at the time of merging.
                nodeHost => nodeDescendant%parent
                ! Find any branch jump and follow the latest of these before the merging time,
                eventPrior     => nodeDescendant%event
                timeJumpLatest =  -huge(0.0d0)
                do while (associated(eventPrior))
                   select type (eventPrior)
                   type is (nodeEventBranchJump)
                      if (eventPrior%time < timeMerger .and. eventPrior%time > timeJumpLatest .and. associated(eventPrior%task)) then
                         nodeHost       => eventPrior%node
                         timeJumpLatest =  eventPrior%time                         
                      end if
                   end select
                   eventPrior => eventPrior%next
                end do
                ! Find the isolated host.
                do while (nodeHost%isSatellite())
                   nodeHost => nodeHost%parent
                end do
                basicHost => nodeHost%basic()
                do while (basicHost%time() < timeMerger)
                   nodeHost  => nodeHost%parent
                   basicHost => nodeHost%basic ()
                end do
                ! We may be merging into this descendant at some mid-point in its position history - seek that position now.
                position        => nodeDescendant%position       ()
                positionHistory =  position      %positionHistory()
                if (positionHistory%exists()) then
                   iHistory    =1_c_size_t
                   countHistory=size(positionHistory%time)
                   do while (iHistory < countHistory .and. timeMerger > positionHistory%time(iHistory))
                      iHistory=iHistory+1_c_size_t
                   end do
                end if
             end if
          else
             ! Node merge target is not available - assume merging with the current host.
             if (time == timeMerger) then
                nodeDescendant => nodeHost
             else
                nodeDescendant => nodeHost%firstChild
             end if
             nodeHost   => nodeDescendant
             basic      => nodeDescendant%basic()
             timeMerger =  basic         %time ()
          end if
          time            =  timeMerger
          timeMerger      =  huge(0.0d0)
          haveMerger      =  .false.
          nodeMergeTarget => null()
          ! We have merged - forget any pre-existing branch jumps.
          haveBranchJump=.false.          
       end if
       ! Check for any event associated with this node if it is a satellite. In a case where node "A" merges with node "B" at
       ! some time, t, and node "B" has a branch jump event at that same time, t, we want to follow the host node through the
       ! branch jump right away.
       if (associated(nodeDescendant) .and. .not.associated(nodeDescendant,nodeHost)) call findEvents(presentOnly=.true.)
       ! Check for a branch jump.
       if (haveBranchJump .and. time == timeBranchJump) then
          if (associated(nodeDescendant,nodeHost)) nodeDescendant => event%node ! Orphanized, so descendant always follows the host.
          nodeHost => event%node
       end if
    end do
    ! Allocate storage for all coefficients.    
    !! NOTE: This is currently inefficient - we allocate enough space for both logarithmic spiral and cubic polynomial
    !! coefficients for each timestep, but for non-satellites we only need the cubic polynomial coefficients.
    allocate(coefficients((countTrace-1_c_size_t)*(2_c_size_t+countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)))
    coefficients=-huge(0.0d0)
    ! Extract the final position to use as a reference point when handling periodic boundaries.
    position => traceTail%node%position()
    if (traceTail%isSatellite) then
       positionHistory  = position                           %positionHistory(                      )
       positionReference=+positionHistory                    %data           (traceTail%iHistory,1:3) &
            &            /self           %cosmologyFunctions_%expansionFactor(traceTail%time        )
    else
       positionReference=+position                           %position       (                      ) &
            &            /self           %cosmologyFunctions_%expansionFactor(traceTail%time        )
    end if
    ! Walk through the trace and compute interpolations.
    traceTail => traceHead
    iTrace    =  0_c_size_t
    do while (associated(traceTail))
       ! Determine the interpolation coefficients for this step.
       if (associated(traceTail%next)) then ! Skip the final step as we can't interpolate into the future.
          ! Catch any duplicated times.
          if (traceTail%time == traceTail%next%time) then
             write (label,'(e12.6)') traceTail%time
             call Error_Report(var_str('duplicated time (')//trim(adjustl(label))//' Gyr) in position trace for node '//node%index()//{introspection:location})
          end if
          ! Store the time.
          coefficients(1_c_size_t+iTrace)=traceTail%time
          ! Determine the type of interpolation to use.
          if (traceTail%isSatellite .and. traceTail%next%isSatellite) then
             ! Logarithmic spiral interpolation - the node is a satellite at this step, and at the next step.
             coefficients(1_c_size_t+(countTrace-1_c_size_t)+iTrace)=1.0d0
             !! Get the logarithmic spiral interpolation for the host halo.
             coefficientsSpiral=computeCoefficientsLogarithmicSpiral(traceTail                )
             !! Get the cubic polynomial interpolation for the host halo.
             coefficientsCubic =computeCoefficientsCubicPolynomial  (traceTail,useHost=.true. )
             !! Store the coefficients.
             offset                                                                   =2_c_size_t*(countTrace-1_c_size_t)+(countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)*iTrace
             coefficients(offset+1_c_size_t:offset+countCoefficientsCubicPolynomial  )=reshape(coefficientsCubic ,[countCoefficientsCubicPolynomial  ])
             offset                                                                   =offset+countCoefficientsCubicPolynomial
             coefficients(offset+1_c_size_t:offset+countCoefficientsLogarithmicSpiral)=reshape(coefficientsSpiral,[countCoefficientsLogarithmicSpiral])
          else
             ! Polynomial interpolation - whenever the node is not a satellite at this and the next step.
             coefficients(1_c_size_t+(countTrace-1_c_size_t)+iTrace)=0.0d0
             !! Get the cubic polynomial coefficients.
             coefficientsCubic=computeCoefficientsCubicPolynomial   (traceTail,useHost=.false.)
             !! Store the coefficients.
             offset                                                                 =2_c_size_t*(countTrace-1_c_size_t)+(countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)*iTrace
             coefficients(offset+1_c_size_t:offset+countCoefficientsCubicPolynomial)=reshape(coefficientsCubic,[countCoefficientsCubicPolynomial])
          end if
       end if
       ! Move to the next step, cleaning up our list as we go.
       traceHead => traceTail%next
       deallocate(traceTail)
       traceTail => traceHead
       iTrace    =  iTrace+1_c_size_t
    end do
    position => node%position()
    call position%floatRank1MetaPropertySet(self%coefficientsID,coefficients)
    return

  contains

    subroutine findEvents(presentOnly)
      !!{
      Locate events associated with this node.
      !!}
      implicit none
      logical                    , intent(in   ) :: presentOnly
      class           (nodeEvent), pointer       :: event_
      double precision                           :: timeEvent

      ! Find any future events associated with this descendant.
      haveBranchJump       =  .false.
      haveSubhaloPromotion =  .false.
      timeEvent            =  huge(0.0d0)
      timeBranchJump       =  huge(0.0d0)
      timeSubhaloPromotion =  huge(0.0d0)
      event                => null(     )
      event_               => nodeDescendant%event
      do while (associated(event_))
         ! Only consider events that occur earlier than the currently found event time. Also, ignore the paired event in the
         ! receiving halo.
         if (event_%time < timeEvent .and. associated(event_%task)) then
            ! Look for a handled event type in the future. Subhalo promotions and branch jumps are handled.
            select type (event_)
            type is (nodeEventSubhaloPromotion)
               if ((.not. presentOnly .and. event_%time > time) .or. (presentOnly .and. event_%time == time)) then
                  haveSubhaloPromotion =  .true.
                  timeSubhaloPromotion =  event_%time
                  event                => event_
                  timeEvent            =  event %time
               end if
            type is (nodeEventBranchJump      )
               if ((.not. presentOnly .and. event_%time > time) .or. (presentOnly .and. event_%time == time)) then
                  haveBranchJump =  .true.
                  timeBranchJump =  event_%time
                  event          => event_
                  timeEvent      =  event %time
               end if
            type is (nodeEventSubhaloPromotionInterTree)
               call Error_Report('inter-tree subhalo promotions are not supported'//{introspection:location})
            type is (nodeEventBranchJumpInterTree)
               call Error_Report('inter-tree branch jumps are not supported'      //{introspection:location})
            end select
         end if
         event_ => event_%next
      end do
      return
    end subroutine findEvents
    
    function computeCoefficientsLogarithmicSpiral(trace) result(coefficientsLogarithmicSpiral)
      !!{
      Compute coefficients of a logarithmic spiral interpolation for position and velocity.
      !!}
      use :: Galacticus_Nodes                , only : nodeComponentPosition
      use :: Histories                       , only : history
      use :: Vectors                         , only : Vector_Product       , Vector_Magnitude
      implicit none
      double precision                       , dimension(20)             :: coefficientsLogarithmicSpiral
      type            (nodeTrace            ), intent(in   )   , target  :: trace
      type            (nodeTrace            )                  , pointer :: trace_
      class           (nodeComponentPosition)                  , pointer :: position                            , positionHost      
      double precision                       , dimension(  2  )          :: time
      double precision                       , dimension(  2,3)          :: positionRelative
      double precision                       , dimension(  2,2)          :: coefficientsAngle                   , coefficientsLogRadius
      double precision                       , dimension(    3)          :: positionSatellite_                  , positionHost_        , &
           &                                                                vectorNormal
      double precision                       , dimension(2,2,3)          :: vectorInPlaneNormal
      double precision                       , parameter                 :: separationTiny               =1.0d-6
      type            (history              )                            :: positionHistory
      integer                                                            :: i                                   , j                    , &
           &                                                                k
      double precision                                                   :: expansionFactor

      ! Iterate over position and velocity.
      do k=1,2
         ! Extract data at start and end times.
         do i=1,2
            ! Select the point in the trace to use.
            select case (i)
            case (1)
               trace_ => trace
            case (2)
               trace_ => trace%next
            end select
            time(i)=trace_%time
            ! Validate that the node is a satellite.
            if (.not.trace_%isSatellite) call Error_Report('expected a satellite node'//{introspection:location})
            ! Extract positional components and history.
            position        => trace_  %node    %position       ()
            positionHost    => trace_  %nodeHost%position       ()
            positionHistory =  position         %positionHistory()
            ! Find displacement vectors (in physical coordinates) from the host center.
            select case (k)
            case (1)
               ! Position.
               positionSatellite_=positionHistory%data    (trace_%iHistory,1:3)
               positionHost_     =positionHost   %position(                   )
               ! Handle periodic positions.
               if (self%isPeriodic) then
                  expansionFactor=self%cosmologyFunctions_%expansionFactor(trace_%time)
                  do j=1,3
                     if (positionSatellite_(j)/expansionFactor > positionReference(j)+0.5d0*self%lengthBox) positionSatellite_(j)=positionSatellite_(j)-self%lengthBox*expansionFactor
                     if (positionSatellite_(j)/expansionFactor < positionReference(j)-0.5d0*self%lengthBox) positionSatellite_(j)=positionSatellite_(j)+self%lengthBox*expansionFactor
                     if (positionHost_     (j)/expansionFactor > positionReference(j)+0.5d0*self%lengthBox) positionHost_     (j)=positionHost_     (j)-self%lengthBox*expansionFactor
                     if (positionHost_     (j)/expansionFactor < positionReference(j)-0.5d0*self%lengthBox) positionHost_     (j)=positionHost_     (j)+self%lengthBox*expansionFactor
                  end do
               end if
               positionRelative(i,:)=positionSatellite_                          -positionHost_
            case (2)
               ! Velocity
               positionRelative(i,:)=positionHistory   %data(trace_%iHistory,4:6)-positionHost %velocity()
            end select
         end do
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
         coefficientsAngle    (k,1)=-time(1)                &
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
      ! Store the computed interpolation coefficients.
      coefficientsLogarithmicSpiral( 1:12)=reshape(vectorInPlaneNormal  ,[12])
      coefficientsLogarithmicSpiral(13:16)=reshape(coefficientsAngle    ,[ 4])
      coefficientsLogarithmicSpiral(17:20)=reshape(coefficientsLogRadius,[ 4])
      return
    end function computeCoefficientsLogarithmicSpiral
    
    function computeCoefficientsCubicPolynomial(trace,useHost) result(coefficientsCubicPolynomial)
      !!{
      Compute coefficients of a cubic polynomial interpolation for position and velocity.
      !!}
      use :: Galacticus_Nodes                , only : nodeComponentPosition
      use :: Histories                       , only : history
      use :: Linear_Algebra                  , only : vector               , matrix, assignment(=)
      use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
      use :: Numerical_Comparison            , only : Values_Agree
      implicit none
      double precision                                  , dimension(4,3)              :: coefficientsCubicPolynomial
      type            (nodeTrace                       ), intent(in   ) , target      :: trace
      logical                                           , intent(in   )               :: useHost
      class           (nodeComponentPosition           )                , pointer     :: position
      type            (nodeTrace                       )                , pointer     :: trace_
      double precision                                  , dimension(2  )              :: times
      double precision                                  , dimension(2,3)              :: positionComoving            , velocityComoving, &
           &                                                                             positionPhysical            , velocityPhysical
      type            (vector                          )                , allocatable :: coordinates                 , coefficients
      type            (matrix                          )                , allocatable :: terms
      integer                                                                         :: i                           , j
      type            (history                         )                              :: positionHistory
      logical                                                                         :: isClone
      
      do i=1,2
         select case (i)
         case (1)
            trace_ => trace
         case (2)
            trace_ => trace%next
         end select
         times(i)=trace_%time
         if (useHost) then
            position => trace_%nodeHost%position()
         else
            position => trace_%node    %position()
         end if
         if (.not.useHost .and. trace_%isSatellite) then
            positionHistory      =position       %positionHistory(                   )
            positionPhysical(i,:)=positionHistory%data           (trace_%iHistory,1:3)
            velocityPhysical(i,:)=positionHistory%data           (trace_%iHistory,4:6)
         else
            positionPhysical(i,:)=position       %position       (                   )
            velocityPhysical(i,:)=position       %velocity       (                   )
         end if
         ! Convert from physical to comoving coordinates, and, for velocities, from km/s to Mpc/Gyr.
         positionComoving(i,:)=positionPhysical(i,:)/self%cosmologyFunctions_%expansionFactor(times(i))
         velocityComoving(i,:)=velocityPhysical(i,:)/self%cosmologyFunctions_%expansionFactor(times(i))/MpcPerKmPerSToGyr
      end do
      ! Handle periodic positions.
      if (self%isPeriodic) then
         do i=1,2
            do j=1,3
               if (positionComoving(i,j) > positionReference(j)+0.5d0*self%lengthBox) positionComoving(i,j)=positionComoving(i,j)-self%lengthBox
               if (positionComoving(i,j) < positionReference(j)-0.5d0*self%lengthBox) positionComoving(i,j)=positionComoving(i,j)+self%lengthBox
            end do
         end do
      end if
      ! Solve for the interpolation coefficients in each Cartesian axis.
      isClone=Values_Agree(times(1),times(2),relTol=2.0d-9)
      do i=1,3
         if (isClone) then
            ! For clones use a simple linear interpolation which should be sufficiently accurate over the tiny offset in times between the clone and original node.
            coefficientsCubicPolynomial(:,i)=[                                                                                      &
                 &                            0.0d0                                                                               , &
                 &                            0.0d0                                                                               , &
                 &                            (-positionComoving(1,i)         +positionComoving(2,i)         )/(times(2)-times(1)), &
                 &                            (+positionComoving(1,i)*times(2)-positionComoving(2,i)*times(1))/(times(2)-times(1))  &
                 &                           ]
         else
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
                 &                                      times(1)**3     ,      times(1)**2    ,times(1)             ,1.0d0                , &
                 &                                      times(2)**3     ,      times(2)**2    ,times(2)             ,1.0d0                , &
                 &                                3.0d0*times(1)**2     ,2.0d0*times(1)       ,1.0d0                ,0.0d0                , &
                 &                                3.0d0*times(2)**2     ,2.0d0*times(2)       ,1.0d0                ,0.0d0                  &
                 &                               ]                                                                                        , &
                 &                               [4,4]                                                                                      &
                 &                              )                                                                                           &
                 &                      )                                                                                                   &
                 &            )
            coefficients                    =terms%linearSystemSolve(coordinates )
            coefficientsCubicPolynomial(:,i)=                        coefficients
            deallocate(terms       )
            deallocate(coordinates )
            deallocate(coefficients)
         end if
      end do
      return
    end function computeCoefficientsCubicPolynomial
      
  end subroutine positionInterpolatedNodeTreeInitialize
