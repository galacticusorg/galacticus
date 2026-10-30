!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a node operator class that interpolates positions of nodes using the approach of :cite:t:`merson_lightcone_2013`.
  !!}
  
  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Galacticus_Nodes   , only : treeNode
  use :: Kind_Numbers       , only : kind_int8
  
  !![
  <enumeration docformat="rst">
   <name>orphanTrajectory</name>
   <description>
   Options for the trajectory followed by an orphaned satellite once its position history is exhausted.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <entry label="host" />
   <entry label="decay"/>
  </enumeration>
  !!]

  !![
  <nodeOperator name="nodeOperatorPositionInterpolated" docformat="rst">
   <description>
   A node operator class that interpolates the 3D comoving positions of nodes between snapshots using the approach of :cite:t:`merson_lightcone_2013`, enabling construction of lightcones from discrete simulation outputs. ``lengthBox`` sets the comoving box side length; ``wrapPeriodic`` controls periodic boundary treatment; ``nodeIndicesReport`` lists node indices for which detailed interpolation information should be reported during processing.

   ``orphanTrajectory`` controls the trajectory of a satellite once its position history is exhausted (i.e. once it becomes an orphan). For ``host`` (the default) the orphan moves to the center of its host halo over the following snapshot interval, and then follows the host. For ``decay`` the orphan instead follows a decaying orbit about its host which reaches the host center at the time of merging. The orbit begins at the final position of the satellite relative to its host, in the plane defined by its final relative position and velocity, and is that of a satellite on a circular orbit in a singular isothermal sphere decaying under dynamical friction :cite:p:`binney_galactic_2008`: the square of the orbital radius decreases linearly in time, while the tangential speed is held fixed at its final value. This orbit is followed between successive steps of the host; over the final step before merging, position and velocity are instead interpolated (with a cubic polynomial) from the orbit to the node into which the orphan merges, which avoids the divergence of the radial speed of the orbit at the time of merging. Position is therefore continuous at the time at which the satellite becomes an orphan, and at the time of merging, but velocity in general is not continuous at the time at which the satellite becomes an orphan. The ``host`` trajectory is used for any orphan which has no finite time of merging, which has no position history at all, or whose host halo itself becomes a satellite before the time of merging (from that time until the orphan merges).
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorPositionInterpolated
     !!{RST
     A node operator class that interpolates positions of nodes using the approach of :cite:t:`merson_lightcone_2013`.
     !!}
     private
     class           (cosmologyFunctionsClass          ), pointer                   :: cosmologyFunctions_ => null()
     integer         (kind_int8                        ), allocatable, dimension(:) :: nodeIndicesReport
     double precision                                                               :: lengthBox
     logical                                                                        :: isPeriodic                   , wrapPeriodic
     integer                                                                        :: coefficientsID
     type            (enumerationOrphanTrajectoryType  )                            :: orphanTrajectory
   contains
     !![
     <methods docformat="rst">
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
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorPositionInterpolated` node operator class.
     !!}
     module procedure positionInterpolatedConstructorParameters
     module procedure positionInterpolatedConstructorInternal
  end interface nodeOperatorPositionInterpolated

  type :: nodeTrace
     !!{RST
     Type used for tracing nodes through their positional history.
     !!}
     double precision                                 :: time
     type            (treeNode ), pointer             :: node             => null() , nodeHost => null()
     type            (nodeTrace), pointer             :: next             => null()
     logical                                          :: isSatellite
     integer         (c_size_t )                      :: iHistory
     logical                                          :: isOrphanDecaying =  .false.
     double precision           , dimension(10)       :: orbitOrphan
  end type nodeTrace

  ! Count of coefficients for different interpolation types.
  integer         (c_size_t), parameter :: countCoefficientsCubicPolynomial  =12_c_size_t
  integer         (c_size_t), parameter :: countCoefficientsLogarithmicSpiral=20_c_size_t
  integer         (c_size_t), parameter :: countCoefficientsOrphanDecay      =10_c_size_t

  ! Labels for the type of interpolation used for each step of the trace.
  double precision          , parameter :: interpolationCubicPolynomial      =0.0d0
  double precision          , parameter :: interpolationLogarithmicSpiral    =1.0d0
  double precision          , parameter :: interpolationOrphanDecay          =2.0d0
  
contains
  
  function positionInterpolatedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorPositionInterpolated` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorPositionInterpolated)                            :: self
    type            (inputParameters                 ), intent(inout)             :: parameters
    class           (cosmologyFunctionsClass         ), pointer                   :: cosmologyFunctions_
    integer         (kind_int8                       ), allocatable, dimension(:) :: nodeIndicesReport
    double precision                                                              :: lengthBox
    logical                                                                       :: wrapPeriodic
    type            (varying_string                  )                            :: orphanTrajectory

    allocate(nodeIndicesReport(parameters%count('nodeIndicesReport',zeroIfNotPresent=.true.)))
    if (parameters%isPresent('nodeIndicesReport')) then
       !![
       <inputParameter docformat="rst">
	 <name>nodeIndicesReport</name>
	 <description>
	 A list of node indices for which reporting should be performed.
	 </description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter docformat="rst">
      <name>lengthBox</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The periodic length of the positions. For non-periodic positions, a value of zero should be given.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wrapPeriodic</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, interpolated positions that lie outside of the periodic box will be wrapped back into the box.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>orphanTrajectory</name>
      <defaultValue>var_str('host')</defaultValue>
      <description>
      The trajectory followed by a satellite once its position history is exhausted. "``host``" moves the orphan to the center of its host halo over the following snapshot interval. "``decay``" follows a decaying orbit about the host which reaches the host center at the time of merging.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nodeOperatorPositionInterpolated(lengthBox,wrapPeriodic,nodeIndicesReport,enumerationOrphanTrajectoryEncode(char(orphanTrajectory),includesPrefix=.false.),cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function positionInterpolatedConstructorParameters

  function positionInterpolatedConstructorInternal(lengthBox,wrapPeriodic,nodeIndicesReport,orphanTrajectory,cosmologyFunctions_) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorPositionInterpolated` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorPositionInterpolated)                              :: self
    double precision                                  , intent(in   )               :: lengthBox
    logical                                           , intent(in   )               :: wrapPeriodic
    class           (cosmologyFunctionsClass         ), intent(in   ), target       :: cosmologyFunctions_
    integer         (kind_int8                       ), intent(in   ), dimension(:) :: nodeIndicesReport
    type            (enumerationOrphanTrajectoryType ), intent(in   )               :: orphanTrajectory
    !![
    <constructorAssign variables="lengthBox, wrapPeriodic, nodeIndicesReport, orphanTrajectory, *cosmologyFunctions_"/>
    !!]

    self%isPeriodic=lengthBox > 0.0d0
    !![
    <addMetaProperty component="position" name="positionInterpolatedCoefficients" id="self%coefficientsID" rank="1" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function positionInterpolatedConstructorInternal

  subroutine positionInterpolatedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorPositionInterpolated` node operator class.
    !!}
    implicit none
    type(nodeOperatorPositionInterpolated), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine positionInterpolatedDestructor

  subroutine positionInterpolatedDifferentialEvolutionAnalytics(self,node)
    !!{RST
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
    !!{RST
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
    !!{RST
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
    !!{RST
    Compute the interpolated position and velocity of the node.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentPosition, nodeComponentBasic
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    use :: Display                         , only : displayMessage       , displayIndent     , displayUnindent
    implicit none
    class           (nodeOperatorPositionInterpolated), intent(inout)              :: self
    type            (treeNode                        ), intent(inout)              :: node
    double precision                                  , intent(in   )              :: time
    class           (nodeComponentPosition           ), pointer                    :: position
    double precision                                  , dimension(3     )          :: position_          , velocity_
    double precision                                  , dimension(4,3   )          :: coefficientsCubic
    double precision                                  , dimension(    20)          :: coefficientsSpiral
    double precision                                  , dimension(  2, 2)          :: coefficientsAngle  , coefficientsLogRadius
    double precision                                  , dimension(2,2, 3)          :: vectorInPlaneNormal
    double precision                                  , dimension(2     )          :: angle              , logRadius
    double precision                                  , dimension(:     ), pointer :: coefficients
    integer         (c_size_t                        )                             :: countTrace         , iTrace               , &
         &                                                                            offset
    integer                                                                        :: i
    logical                                                                        :: isSpiral           , isOrphanDecay        , &
         &                                                                            report
    double precision                                  , dimension(3     )          :: positionOrphan     , velocityOrphan
    double precision                                                               :: lengthBox
    character       (len=1024                        )                             :: label
    !$GLC attributes initialized :: coefficients

    report=any(node%index() == self%nodeIndicesReport)
    if (report) then
       write (label,'(i12," | ",e12.6)') node%index(),time
       call displayIndent("Position interpolation for node | time: "//trim(adjustl(label)))
    end if
    ! Extract all interpolation coefficients.
    position     => node    %position                          (                   )
    coefficients => position%floatRank1MetaPropertyGetReference(self%coefficientsID)
    ! Determine the number of steps in the interpolation.
    countTrace  =size(coefficients)/(2_c_size_t+countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)
    ! Find the appropriate time.
    iTrace=1_c_size_t
    do while (iTrace < countTrace .and. time >= coefficients(iTrace+1_c_size_t))
       iTrace=iTrace+1_c_size_t
    end do
    if (report) then
       write (label,'(i4,1x,"of",1x,i4," at time ",e12.6)') iTrace,countTrace,coefficients(iTrace)
       call displayMessage("Index: "//trim(adjustl(label)))
    end if
    ! Extract cubic polynomial coefficients for this time.
    isSpiral         =coefficients(countTrace+iTrace) == interpolationLogarithmicSpiral
    isOrphanDecay    =coefficients(countTrace+iTrace) == interpolationOrphanDecay
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
    if (report) then
       write (label,'(12(e12.6,1x))                ') coefficientsCubic
       call displayMessage("Cubic polynomial coefficients:        "//trim(adjustl(label)))
       write (label,'(3(e12.6,1x),"| ",3(e12.6,1x))') position_,velocity_
       call displayMessage("Cubic polynomial position|velocity:   "//trim(adjustl(label)))
    end if
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
       if (report) then
          write (label,'(20(e12.6,1x))'                ) coefficientsSpiral
          call displayMessage("Logarithmic spiral coefficients:      "//trim(adjustl(label)))
          write (label,'(3(e12.6,1x),"| ",3(e12.6,1x))') position_,velocity_
          call displayMessage("Logarithmic spiral position|velocity: "//trim(adjustl(label)))
       end if
    else if (isOrphanDecay) then
       ! Add on the offset of an orphan on its decaying orbit, in physical position.
       offset   =offset+countCoefficientsCubicPolynomial
       call positionInterpolatedOrphanOffset(coefficients(offset+1_c_size_t:offset+countCoefficientsOrphanDecay),time,positionOrphan,velocityOrphan)
       position_=+position_      &
            &    +positionOrphan
       velocity_=+velocity_      &
            &    +velocityOrphan
       if (report) then
          write (label,'(10(e12.6,1x))'                ) coefficients(offset+1_c_size_t:offset+countCoefficientsOrphanDecay)
          call displayMessage("Orphan decaying orbit coefficients:   "   //trim(adjustl(label)))
          write (label,'(3(e12.6,1x),"| ",3(e12.6,1x))') position_,velocity_
          call displayMessage("Orphan decaying orbit position|velocity: "//trim(adjustl(label)))
       end if
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
    if (report) call displayUnindent("")
    return
  end subroutine positionInterpolatedInterpolate

  subroutine positionInterpolatedNodeTreeInitialize(self,node)
    !!{RST
    Trigger interpolation calculation at node initialization
    !!}
    use :: Error                       , only : Error_Report
    use :: Display                     , only : displayIndent            , displayUnindent      , displayMessage
    use :: Galacticus_Nodes            , only : nodeComponentBasic       , nodeComponentPosition, nodeComponentSatellite, nodeEvent, &
         &                                      nodeEventSubhaloPromotion, nodeEventBranchJump
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
    class           (nodeEvent                       )                 , pointer     :: event               , eventPrior
    type            (nodeTrace                       )                 , pointer     :: traceTail           , traceHead 
    double precision                                  , dimension(   3)              :: positionReference
    double precision                                  , dimension(4 ,3)              :: coefficientsCubic
    double precision                                  , dimension(20  )              :: coefficientsSpiral
    double precision                                  , dimension(10  )              :: orbitOrphan
    double precision                                  , dimension( :  ), allocatable :: coefficients
    character       (len=128                         )                               :: label
    character       (len=512                         )                               :: message
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
         &                                                                              isHost              , reporting                 , &
         &                                                                              isOrphanDecaying    , isOrphanBegun             , &
         &                                                                              isOrphanReverted    , isOrphanDecayEnded

    ! Consider only branch tips.
    if (associated(node%firstChild)) return
    reporting=any(node%index() == self%nodeIndicesReport)
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
    ! Set the initial time.
    basic => nodeDescendant%basic()
    time  =  basic         %time ()
    ! Apply any branch jump that occurs at the initial time of this branch tip. A branch jump redirects the host (and, for a
    ! self-hosting node, the descendant) onto the branch into which the node jumps. The main tree walk does not catch a jump at
    ! the initial time: findEvents(presentOnly=.false.) at the top of the loop finds only jumps in the future, while the
    ! present-time check (findEvents(presentOnly=.true.)) is reached only after the time has been advanced to the next history
    ! step. We therefore apply it explicitly here, so that satellite positions are interpolated relative to the correct host
    ! from the very first step - consistent with how the merger logic follows such jumps when another branch merges onto this node.
    eventPrior => node%event
    do while (associated(eventPrior))
       select type (eventPrior)
       type is (nodeEventBranchJump)
          if (eventPrior%time == time .and. associated(eventPrior%task)) then
             if (associated(nodeDescendant,nodeHost)) nodeDescendant => eventPrior%node
             nodeHost => eventPrior%node
          end if
       end select
       eventPrior => eventPrior%next
    end do
    isSatellite       =.not.associated(nodeDescendant,nodeHost)
    isInitialSatellite=isSatellite
    ! Initialize history status.
    iHistory    =-1_c_size_t
    countHistory=+0_c_size_t
    ! Initialize merger status.
    haveMerger      =  .false.
    timeMerger      =  huge(0.0d0)
    nodeMergeTarget => null()    
    ! Initialize orphan status.
    isOrphanDecaying  =.false.
    isOrphanReverted  =.false.
    isOrphanDecayEnded=.false.
    orbitOrphan       =0.0d0
    ! Follow the node through the tree.
    if (reporting) then
       write (label,'(i12)') node%index()
       call displayIndent('Tracing node: '//trim(adjustl(label)))
       call displayMessage('nodeIndex    hostIndex    time         timeHost     isSatellite iHistory countHistory')
       call displayMessage('_____________________________________________________________________________________')
    end if
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
          if (reporting) then
             write (label,'(e12.6)') timeMerger
             call displayMessage("merger - time: "//trim(adjustl(label)))
          end if
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
       ! Append this step to the trace - unless an orphan has just reverted from a decaying orbit to following its host, in which
       ! case the final step in the trace has already been updated to describe this step.
       if (.not.isOrphanReverted) then
          countTrace=countTrace+1_c_size_t
          if (associated(traceTail)) then
             allocate(traceTail%next)
             traceTail => traceTail%next
          else
             allocate(traceHead     )
             traceTail => traceHead
          end if
          if (reporting) then
             write (label,'(i12,1x,i12,1x,e12.6,1x,e12.6,1x,l1,1x,i4,1x,i4)') nodeDescendant%index(),nodeHost%index(),time,basicHost%time(),isSatellite,iHistory,countHistory
             call displayMessage(trim(adjustl(label)))
          end if
          traceTail%time        =  time
          traceTail%isSatellite =  isSatellite
          traceTail%iHistory    =  iHistory
          traceTail%node        => nodeDescendant
          traceTail%nodeHost    => nodeHost
          traceTail%next        => null()
          if (isOrphanDecaying) then
             traceTail%isOrphanDecaying=.true.
             traceTail%orbitOrphan     =orbitOrphan
          end if
       end if
       isOrphanReverted=.false.
       ! Move to the next step.
       !! Case where we are following the decaying orbit of an orphan.
       if (isOrphanDecaying) then
          call orphanDecayStep()
          if (isOrphanReverted) cycle
          ! Update the time if the orphan has reverted to following its host.
          if (.not.isOrphanDecaying .and. associated(nodeDescendant)) then
             basic => nodeDescendant%basic()
             time  =  basic         %time ()
          end if
       !! Case where we want to move to the parent node.
       else if (                                                       &
            &                    nodeDescendant%isPrimaryProgenitor()  & ! 1. Node is primary progenitor - so we can move directly to its descendant.
            &  .or.                                                    &
            &    .not.associated(nodeDescendant%parent               ) & ! 2. Node has no descendant - we have reached the end of the tree - move to that null
            &  .or.                                                    & !    descendant to terminate the tree walk.
            &    .not.haveHistory                                      & ! 3. Node is not the primary progenitor, but has no position history. It becomes an orphan
            &  .or.                                                    & !    in its parent, so move our pointer to that parent such that we will follow its position.
            &   (     haveHistory .and. iHistory == countHistory)      & ! 4. We have reached the end of the node's position history.  It becomes an orphan
            & ) then                                                     !    in its parent, so move our pointer to that parent such that we will follow its position.
          !! Check for a subhalo promotion.
          if (reporting) then
             write (label,'(l1)') haveSubhaloPromotion
             call displayMessage("trace to parent - promotion?: "//trim(adjustl(label)))
          end if
          if (haveSubhaloPromotion) then
             ! We have a subhalo promotion - move our pointer to that node.
             nodeDescendant    => event         %node
             nodeHost          => nodeDescendant               ! This is always an isolated halo by construction, so must be self-hosting.
             ! The promoted halo survives as an independent central, so cancel any pending (stale) galaxy merger.
             haveMerger         =  .false.
             timeMerger         =  huge(0.0d0)
             nodeMergeTarget    => null()
             isOrphanDecayEnded =  .false.
          else
             ! A satellite which has reached the end of its position history may become an orphan following a decaying orbit.
             isOrphanBegun=.false.
             if (.not.associated(nodeDescendant,nodeHost) .and. haveHistory) then
                if (iHistory == countHistory) isOrphanBegun=orphanDecayBegin()
             end if
             if (isOrphanBegun) then
                call orphanDecayStep()
                if (isOrphanReverted) cycle
             else
                ! Always move to the host's parent - this allows us to correctly handle satellites that have been orphanized.
                if (.not.associated(nodeDescendant,nodeHost) .and. associated(nodeHost%parent) .and. nodeHost%parent%index() == nodeHost%index()) then
                   nodeDescendant => nodeHost      %parent%parent ! Parent is a clone (and descendant is a satellite), so skip over it.
                else
                   nodeDescendant => nodeHost      %parent        ! Parent is not a clone, use directly.
                end if
                nodeHost          => nodeDescendant               ! This is always an isolated halo by construction, so must be self-hosting.
             end if
          end if
          ! Update the time if the descendant exists.
          if (associated(nodeDescendant) .and. .not.isOrphanDecaying) then
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
          if (reporting) then
             write (label,'(l1)') iHistory == countHistory
             call displayMessage("trace history - last?: "//trim(adjustl(label)))
          end if
          if (iHistory <= countHistory) then
             ! Update the time to this new history step.
             time=positionHistory%time(iHistory)
             ! Move the host pointer.
             basicHost => nodeHost%basic()
             do while (associated(nodeHost) .and. basicHost%time() < time)
                nodeHost => nodeHost%parent
                if (associated(nodeHost)) then
                   basicHost => nodeHost%basic()
                   if (reporting) then
                      write (label,'(2(e12.6,1x))') time,basicHost%time()
                      call displayMessage("trace host - time, timeHost: "//trim(adjustl(label)))
                   end if
                else
                   if (reporting) then
                      write (label,'(1(e12.6,1x),"???")') time
                      call displayMessage("trace host - time, timeHost: "//trim(adjustl(label)))
                   end if
                end if
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
                ! The promoted halo survives as an independent central, so cancel any pending (stale) galaxy merger.
                haveMerger         =  .false.
                timeMerger         =  huge(0.0d0)
                nodeMergeTarget    => null()
                isOrphanDecayEnded =  .false.
             else
                if (orphanDecayBegin()) then
                   ! The satellite becomes an orphan following a decaying orbit. Step back to the final entry in the history, which
                   ! is the final position of the satellite.
                   iHistory=countHistory
                   call orphanDecayStep()
                   if (isOrphanReverted) cycle
                else
                   ! This descendant becomes orphanized in its parent, so move our pointer to that parent such that we will follow
                   ! its position.
                   nodeDescendant => nodeHost      %parent
                   nodeHost       => nodeDescendant        ! This is always an isolated halo by construction, so must be self-hosting.
                end if
             end if
             ! Update the time if the descendant exists.
             if (associated(nodeDescendant) .and. .not.isOrphanDecaying) then
                basic => nodeDescendant%basic()
                time  =  basic         %time ()
             end if
          end if
       end if
       ! Check for a merger.
       if (haveMerger .and. time >= timeMerger) then
          ! Any decaying orbit of an orphan ends at the merger.
          isOrphanDecaying=.false.
          ! Move to the host, unless the merger time is before the current time, in which case we need to back up to the host's
          ! progenitor.
          if (reporting) then
             if (associated(nodeMergeTarget)) then
                write (label,'(l1,a2,i12,a1)') associated(nodeMergeTarget),' [',nodeMergeTarget%index(),']'
             else
                write (label,'(l1)'          ) associated(nodeMergeTarget)
             end if
             call displayMessage("trace merger - target?: "//trim(adjustl(label)))
          end if
          if (associated(nodeMergeTarget)) then
             ! A merge target is available - move to that node.
             nodeDescendant => nodeMergeTarget
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
          time               =  timeMerger
          timeMerger         =  huge(0.0d0)
          haveMerger         =  .false.
          nodeMergeTarget    => null()
          isOrphanDecayEnded =  .false.
          ! We have merged - forget any pre-existing branch jumps.
          haveBranchJump=.false.          
       end if
       ! Check for any event associated with this node if it is a satellite. In a case where node "A" merges with node "B" at
       ! some time, t, and node "B" has a branch jump event at that same time, t, we want to follow the host node through the
       ! branch jump right away.
       if (associated(nodeDescendant) .and. .not.associated(nodeDescendant,nodeHost)) call findEvents(presentOnly=.true.)
       ! Check for a branch jump.
       if (reporting) then
          write (label,'(l1,1x,2(e12.6,1x))') haveBranchJump,time,timeBranchJump
          call displayMessage("trace events - branch jump?, time, timeJump: "//trim(adjustl(label)))
       end if
       if (haveBranchJump .and. time == timeBranchJump) then
          if (associated(nodeDescendant,nodeHost) .or. isOrphanDecaying) nodeDescendant => event%node ! Orphanized, so descendant always follows the host.
          nodeHost         => event%node
          isOrphanDecaying =  .false.
       end if
    end do
    if (reporting) call displayUnindent("")
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
    if (reporting) then
       write (label,'(i12)') node%index()
       call displayIndent('Trace report for node: '//trim(adjustl(label)))
       call displayMessage('nodeIndex    hostIndex    time     isSatellite indexHistory')
       call displayMessage('___________________________________________________________')
    end if
    traceTail => traceHead
    iTrace    =  0_c_size_t
    do while (associated(traceTail))
       if (reporting) then
          write (message,'(i12,1x,i12,1x,f8.5,1x,l1,10x,i4)') traceTail%node%index(),traceTail%nodeHost%index(),traceTail%time,traceTail%isSatellite,traceTail%iHistory
          call displayMessage(message)
       end if
       ! Determine the interpolation coefficients for this step.
       if (associated(traceTail%next)) then ! Skip the final step as we can't interpolate into the future.
          ! Catch any duplicated times.
          if (traceTail%time == traceTail%next%time) then
             if (reporting) then
                write (message,'(i12,1x,i12,1x,f8.5,1x,l1,10x,i4)') traceTail%next%node%index(),traceTail%next%nodeHost%index(),traceTail%next%time,traceTail%next%isSatellite,traceTail%next%iHistory
                call displayMessage(message)
             end if
             write (label,'(e12.6)') traceTail%time
             call Error_Report(var_str('duplicated time (')//trim(adjustl(label))//' Gyr) in position trace for node '//node%index()//{introspection:location})
          end if
          ! Store the time.
          coefficients(1_c_size_t+iTrace)=traceTail%time
          ! Determine the type of interpolation to use.
          if (traceTail%isOrphanDecaying .and. traceTail%next%isOrphanDecaying) then
             ! Decaying orbit of an orphan - the node is following its decaying orbit at this step, and at the next step.
             coefficients(1_c_size_t+(countTrace-1_c_size_t)+iTrace)=interpolationOrphanDecay
             !! Get the cubic polynomial interpolation for the host halo.
             coefficientsCubic =computeCoefficientsCubicPolynomial(traceTail,useHost=.true. )
             !! Store the coefficients.
             offset                                                                 =2_c_size_t*(countTrace-1_c_size_t)+(countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)*iTrace
             coefficients(offset+1_c_size_t:offset+countCoefficientsCubicPolynomial)=reshape(coefficientsCubic,[countCoefficientsCubicPolynomial])
             offset                                                                 =offset+countCoefficientsCubicPolynomial
             coefficients(offset+1_c_size_t:offset+countCoefficientsOrphanDecay    )=traceTail%orbitOrphan
             if (reporting) then
                call displayIndent("orphan")
                write (message,'(10(e12.6,1x))') traceTail%orbitOrphan
                call displayMessage("   orphan: "//trim(message))
                write (message,'(12(e12.6,1x))') reshape(coefficientsCubic ,[countCoefficientsCubicPolynomial  ])
                call displayMessage("     host: "//trim(message))
                call displayUnindent("")
             end if
          else if (traceTail%isSatellite .and. traceTail%next%isSatellite .and. .not.traceTail%isOrphanDecaying) then
             ! Logarithmic spiral interpolation - the node is a satellite at this step, and at the next step.
             coefficients(1_c_size_t+(countTrace-1_c_size_t)+iTrace)=interpolationLogarithmicSpiral
             !! Get the logarithmic spiral interpolation for the host halo.
             coefficientsSpiral=computeCoefficientsLogarithmicSpiral(traceTail                )
             !! Get the cubic polynomial interpolation for the host halo.
             coefficientsCubic =computeCoefficientsCubicPolynomial  (traceTail,useHost=.true. )
             !! Store the coefficients.
             offset                                                                   =2_c_size_t*(countTrace-1_c_size_t)+(countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)*iTrace
             coefficients(offset+1_c_size_t:offset+countCoefficientsCubicPolynomial  )=reshape(coefficientsCubic ,[countCoefficientsCubicPolynomial  ])
             offset                                                                   =offset+countCoefficientsCubicPolynomial
             coefficients(offset+1_c_size_t:offset+countCoefficientsLogarithmicSpiral)=reshape(coefficientsSpiral,[countCoefficientsLogarithmicSpiral])
             if (reporting) then
                call displayIndent("satellite")
                write (message,'(20(e12.6,1x))') reshape(coefficientsSpiral,[countCoefficientsLogarithmicSpiral])
                call displayMessage("satellite: "//trim(message))
                write (message,'(12(e12.6,1x))') reshape(coefficientsCubic ,[countCoefficientsCubicPolynomial  ])
                call displayMessage("     host: "//trim(message))
                call displayUnindent("")
             end if
          else
             ! Polynomial interpolation - whenever the node is not a satellite at this and the next step, or is ending a decaying
             ! orbit as an orphan.
             coefficients(1_c_size_t+(countTrace-1_c_size_t)+iTrace)=interpolationCubicPolynomial
             !! Get the cubic polynomial coefficients.
             coefficientsCubic=computeCoefficientsCubicPolynomial   (traceTail,useHost=.false.)
             !! Store the coefficients.
             offset                                                                 =2_c_size_t*(countTrace-1_c_size_t)+(countCoefficientsCubicPolynomial+countCoefficientsLogarithmicSpiral)*iTrace
             coefficients(offset+1_c_size_t:offset+countCoefficientsCubicPolynomial)=reshape(coefficientsCubic,[countCoefficientsCubicPolynomial])
             if (reporting) then
                call displayIndent("non-satellite")
                write (message,'(12(e12.6,1x))') reshape(coefficientsCubic ,[countCoefficientsCubicPolynomial  ])
                call displayMessage("     self: "//trim(message))
                call displayUnindent("")
             end if
          end if
       end if
       ! Move to the next step, cleaning up our list as we go.
       traceHead => traceTail%next
       deallocate(traceTail)
       traceTail => traceHead
       iTrace    =  iTrace+1_c_size_t
    end do
    if (reporting) call displayUnindent('done')
    position => node%position()
    call position%floatRank1MetaPropertySet(self%coefficientsID,coefficients)
    return

  contains

    subroutine findEvents(presentOnly)
      !!{RST
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
            end select
         end if
         event_ => event_%next
      end do
      return
    end subroutine findEvents

    logical function orphanDecayBegin()
      !!{RST
      Begin a decaying orbit for a satellite which has reached the final entry in its position history (which must be the final
      step in the trace), if the ``decay`` orphan trajectory is selected and the orbit can be constructed. Returns true if the
      decaying orbit was begun.
      !!}
      use :: Vectors, only : Vector_Product, Vector_Magnitude
      implicit none
      class           (nodeComponentPosition), pointer      :: positionSatellite, positionHost_
      class           (nodeComponentBasic   ), pointer      :: basicHost_
      type            (history              )               :: positionHistory_
      double precision                       , dimension(3) :: positionRelative , velocityRelative   , &
           &                                                   vectorRadial     , vectorTangential   , &
           &                                                   vectorNormal
      double precision                                      :: expansionFactor  , expansionFactorHost, &
           &                                                   lengthBox_       , momentumAngular
      integer                                               :: j

      orphanDecayBegin=.false.
      ! A decaying orbit requires that it be selected, that the satellite has a finite time of merging after the final entry in
      ! its position history, and that it is a satellite at that final entry. A decaying orbit which has ended because the host
      ! became a satellite is not resumed - the orphan then follows its host until it merges.
      if (self%orphanTrajectory /= orphanTrajectoryDecay) return
      if (     isOrphanDecayEnded                       ) return
      if (.not.haveMerger                               ) return
      if (.not.traceTail%isSatellite                    ) return
      if (     timeMerger       <= traceTail%time       ) return
      ! Find the final position and velocity of the satellite relative to its host. As for the logarithmic spiral interpolation,
      ! the physical position of the host is translated to that at the expansion factor of the satellite.
      positionSatellite   => traceTail%node    %position       ()
      positionHost_       => traceTail%nodeHost%position       ()
      basicHost_          => traceTail%nodeHost%basic          ()
      positionHistory_    =  positionSatellite %positionHistory()
      expansionFactor     =  self%cosmologyFunctions_%expansionFactor(traceTail %time  )
      expansionFactorHost =  self%cosmologyFunctions_%expansionFactor(basicHost_%time())
      positionRelative    = +positionHistory_%data    (traceTail%iHistory,1:3) &
           &                -positionHost_   %position(                      ) &
           &                *expansionFactor                                   &
           &                /expansionFactorHost
      velocityRelative    = +positionHistory_%data    (traceTail%iHistory,4:6) &
           &                -positionHost_   %velocity(                      )
      ! Handle periodic positions by taking the nearest periodic image of the satellite.
      if (self%isPeriodic) then
         lengthBox_=self%lengthBox*expansionFactor
         do j=1,3
            positionRelative(j)=positionRelative(j)-lengthBox_*dble(nint(positionRelative(j)/lengthBox_))
         end do
      end if
      ! A satellite located precisely at the center of its host has no defined orbit.
      if (Vector_Magnitude(positionRelative) <= 0.0d0) return
      ! Construct unit vectors in the radial direction, and in the tangential direction in the orbital plane.
      vectorRadial   =+positionRelative                                    &
           &          /Vector_Magnitude(positionRelative                 )
      vectorNormal   =+Vector_Product  (positionRelative,velocityRelative)
      momentumAngular=+Vector_Magnitude(vectorNormal                     )
      if (momentumAngular > 0.0d0) then
         vectorNormal=+vectorNormal    &
              &       /momentumAngular
      else
         ! The orbit is purely radial, so the orbital plane is undefined. Choose any plane containing the radial vector, using the
         ! Cartesian axis which is most nearly perpendicular to it.
         vectorNormal                                 =0.0d0
         vectorNormal(minloc(abs(vectorRadial),dim=1))=1.0d0
         vectorNormal                                 =+                 Vector_Product(vectorRadial,vectorNormal)  &
              &                                        /Vector_Magnitude(Vector_Product(vectorRadial,vectorNormal))
      end if
      vectorTangential=Vector_Product(vectorNormal,vectorRadial)
      ! Store the parameters of the decaying orbit.
      orbitOrphan( 1: 3)=vectorRadial
      orbitOrphan( 4: 6)=vectorTangential
      orbitOrphan( 7   )=Vector_Magnitude(positionRelative                 )
      orbitOrphan( 8   )=traceTail%time
      orbitOrphan( 9   )=timeMerger
      orbitOrphan(10   )=Dot_Product     (velocityRelative,vectorTangential)
      ! The final step in the trace is the start of the decaying orbit.
      isOrphanDecaying          =.true.
      traceTail%isOrphanDecaying=.true.
      traceTail%orbitOrphan     =orbitOrphan
      orphanDecayBegin          =.true.
      if (reporting) then
         write (message,'(10(e12.6,1x))') orbitOrphan
         call displayMessage("orphan - begin decaying orbit: "//trim(message))
      end if
      return
    end function orphanDecayBegin

    subroutine orphanDecayStep()
      !!{RST
      Advance an orphan on its decaying orbit to the next step of its host. If the end of the tree is reached the trace ends on the
      decaying orbit. If the host itself becomes a satellite the decaying orbit ends, and the orphan reverts to following its host
      (as for the ``host`` trajectory) until it merges: the final step in the trace is then updated to place the orphan at the
      center of its host, ``isOrphanReverted`` is set to indicate that the trace walk must continue from that step, and
      ``isOrphanDecayEnded`` is set to prevent a further decaying orbit from beginning when the position history of the host is
      exhausted.
      !!}
      implicit none
      class(nodeComponentBasic), pointer :: basicHost_

      do while (.true.)
         if (.not.associated(nodeHost%parent)) then
            ! The end of the tree has been reached - the trace ends on the decaying orbit.
            isOrphanDecaying =  .false.
            nodeDescendant   => null()
            if (reporting) call displayMessage("orphan - end decaying orbit at end of tree")
            exit
         else if (nodeHost%isPrimaryProgenitor()) then
            ! Move to the next step of the host, skipping any which do not lie after the final step in the trace.
            nodeHost   => nodeHost  %parent
            basicHost_ => nodeHost  %basic ()
            time       =  basicHost_%time  ()
            if (time > traceTail%time) exit
         else if (traceTail%time > traceTail%orbitOrphan(8)) then
            ! The host is itself becoming a satellite. The orphan reverts to following its host, as for the "host" trajectory. The
            ! final step in the trace is updated to place the orphan at the center of its host (the interpolation over the
            ! preceding step then carries the orphan smoothly from its decaying orbit to the host center), and the trace walk
            ! continues from that step, following the host.
            traceTail%node             => nodeHost
            traceTail%isSatellite      =  .false.
            traceTail%isOrphanDecaying =  .false.
            traceTail%iHistory         =  -1_c_size_t
            nodeDescendant             => nodeHost
            iHistory                   =  -1_c_size_t
            isOrphanDecaying           =  .false.
            isOrphanReverted           =  .true.
            isOrphanDecayEnded         =  .true.
            if (reporting) call displayMessage("orphan - end decaying orbit as host becomes a satellite")
            exit
         else
            ! The host is itself becoming a satellite at the first step of the decaying orbit. The orphan reverts to following its
            ! host, exactly as for the "host" trajectory.
            isOrphanDecaying   =  .false.
            isOrphanDecayEnded =  .true.
            nodeDescendant     => nodeHost%parent
            nodeHost           => nodeDescendant
            if (reporting) call displayMessage("orphan - end decaying orbit as host becomes a satellite at its first step")
            exit
         end if
      end do
      return
    end subroutine orphanDecayStep
    
    function computeCoefficientsLogarithmicSpiral(trace) result(coefficientsLogarithmicSpiral)
      !!{RST
      Compute coefficients of a logarithmic spiral interpolation for position and velocity.
      !!}
      use :: Galacticus_Nodes                , only : nodeComponentPosition, nodeComponentBasic
      use :: Histories                       , only : history
      use :: Vectors                         , only : Vector_Product       , Vector_Magnitude
      implicit none
      double precision                       , dimension(20)             :: coefficientsLogarithmicSpiral
      type            (nodeTrace            ), intent(in   )   , target  :: trace
      type            (nodeTrace            )                  , pointer :: trace_
      class           (nodeComponentPosition)                  , pointer :: position                            , positionHost      
      class           (nodeComponentBasic   )                  , pointer :: basicHost
      double precision                       , dimension(  2  )          :: time
      double precision                       , dimension(  2,3)          :: positionRelative                    , positionSatelliteComoving, &
           &                                                                positionHostComoving
      double precision                       , dimension(  2,2)          :: coefficientsAngle                   , coefficientsLogRadius
      double precision                       , dimension(    3)          :: positionSatellite_                  , positionHost_            , &
           &                                                                vectorNormal
      double precision                       , dimension(2,2,3)          :: vectorInPlaneNormal
      double precision                       , parameter                 :: separationTiny               =1.0d-6
      type            (history              )                            :: positionHistory
      integer                                                            :: i                                   , j                        , &
           &                                                                k
      double precision                                                   :: expansionFactor                     , expansionFactorHost

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
            basicHost       => trace_  %nodeHost%basic          ()
            positionHistory =  position         %positionHistory()
            ! Find displacement vectors (in physical coordinates) from the host center.
            select case (k)
            case (1)
               ! Position.
               positionSatellite_=positionHistory%data    (trace_%iHistory,1:3)
               positionHost_     =positionHost   %position(                   )
               ! Get expansion factors.
               expansionFactor    =self%cosmologyFunctions_%expansionFactor(trace_   %time  )
               expansionFactorHost=self%cosmologyFunctions_%expansionFactor(basicHost%time())
               ! Handle periodic positions.
               if (self%isPeriodic) then
                  do j=1,3
                     if (positionSatellite_(j)/expansionFactor     > positionReference(j)+0.5d0*self%lengthBox) positionSatellite_(j)=positionSatellite_(j)-self%lengthBox*expansionFactor
                     if (positionSatellite_(j)/expansionFactor     < positionReference(j)-0.5d0*self%lengthBox) positionSatellite_(j)=positionSatellite_(j)+self%lengthBox*expansionFactor
                     if (positionHost_     (j)/expansionFactorHost > positionReference(j)+0.5d0*self%lengthBox) positionHost_     (j)=positionHost_     (j)-self%lengthBox*expansionFactorHost
                     if (positionHost_     (j)/expansionFactorHost < positionReference(j)-0.5d0*self%lengthBox) positionHost_     (j)=positionHost_     (j)+self%lengthBox*expansionFactorHost
                  end do
               end if
               ! Record comoving positions for reporting.
               if (reporting) then
                  positionSatelliteComoving(i,:)=positionSatellite_/expansionFactor
                  positionHostComoving     (i,:)=positionHost_     /expansionFactorHost
               end if
               ! If the host is not at the same expansion factor as the satellite we must translate its physical position to what it would be at the expansion factor of the satellite.
               positionHost_=positionHost_*expansionFactor/expansionFactorHost
               ! Compute the relative physical position.
               positionRelative         (i,:)=positionSatellite_                          -positionHost_
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
      ! Report.
      if (reporting) then
         call displayIndent("logarithmic polynomial")
         call displayIndent("satellite")
         write (message,'(3(e12.6,1x))') positionSatelliteComoving(1,:)
         call displayMessage("start: "//trim(message))
         write (message,'(3(e12.6,1x))') positionSatelliteComoving(2,:)
         call displayMessage("  end: "//trim(message))
         call displayUnindent("")
         call displayIndent("     host")
         write (message,'(3(e12.6,1x))') positionHostComoving     (1,:)
         call displayMessage("start: "//trim(message))
         write (message,'(3(e12.6,1x))') positionHostComoving     (2,:)
         call displayMessage("  end: "//trim(message))
         call displayUnindent("")
         call displayUnindent("")
      end if
      ! Store the computed interpolation coefficients.
      coefficientsLogarithmicSpiral( 1:12)=reshape(vectorInPlaneNormal  ,[12])
      coefficientsLogarithmicSpiral(13:16)=reshape(coefficientsAngle    ,[ 4])
      coefficientsLogarithmicSpiral(17:20)=reshape(coefficientsLogRadius,[ 4])
      return
    end function computeCoefficientsLogarithmicSpiral
    
    function computeCoefficientsCubicPolynomial(trace,useHost) result(coefficientsCubicPolynomial)
      !!{RST
      Compute coefficients of a cubic polynomial interpolation for position and velocity.
      !!}
      use :: Galacticus_Nodes                , only : nodeComponentPosition, nodeComponentBasic
      use :: Histories                       , only : history
      use :: Linear_Algebra                  , only : vector               , matrix, assignment(=)
      use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
      use :: Numerical_Comparison            , only : Values_Agree
      implicit none
      double precision                                  , dimension(4,3)              :: coefficientsCubicPolynomial
      type            (nodeTrace                       ), intent(in   ) , target      :: trace
      logical                                           , intent(in   )               :: useHost
      class           (nodeComponentBasic              )                , pointer     :: basic
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
      double precision                                                                :: expansionFactor
      character       (len=512                         )                              :: message

      do i=1,2
         select case (i)
         case (1)
            trace_ => trace
         case (2)
            trace_ => trace%next
         end select
         times(i)=trace_%time
         if (useHost) then
            position        => trace_%nodeHost           %position       (              )
            basic           => trace_%nodeHost           %basic          (              )
            expansionFactor =  self  %cosmologyFunctions_%expansionFactor(basic%time ( ))
         else
            position        => trace_%node               %position       (              )
            expansionFactor =  self  %cosmologyFunctions_%expansionFactor(      times(i))
         end if
         if (.not.useHost .and. trace_%isOrphanDecaying .and. trace_%time > trace_%orbitOrphan(8)) then
            ! An orphan on its decaying orbit (after the final entry in its position history) - its position is that of the host
            ! plus its offset on the decaying orbit.
            position             =>  trace_%nodeHost%position()
            call positionInterpolatedOrphanOffset(trace_%orbitOrphan,times(i),positionPhysical(i,:),velocityPhysical(i,:))
            positionPhysical(i,:)=  +positionPhysical(i,:) &
                 &                  +position%position()
            velocityPhysical(i,:)=  +velocityPhysical(i,:) &
                 &                  +position%velocity()
         else if (.not.useHost .and. trace_%isSatellite) then
            positionHistory      =position       %positionHistory(                   )
            positionPhysical(i,:)=positionHistory%data           (trace_%iHistory,1:3)
            velocityPhysical(i,:)=positionHistory%data           (trace_%iHistory,4:6)
         else
            positionPhysical(i,:)=position       %position       (                   )
            velocityPhysical(i,:)=position       %velocity       (                   )
         end if
         ! Convert from physical to comoving coordinates, and, for velocities, from km/s to Mpc/Gyr.
         positionComoving(i,:)=positionPhysical(i,:)/expansionFactor
         velocityComoving(i,:)=velocityPhysical(i,:)/expansionFactor/MpcPerKmPerSToGyr
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
      ! Report.
      if (reporting) then
         call displayIndent("cubic polynomial")
         write (message,'(3(e12.6,1x))') positionComoving(1,:)
         call displayMessage("start: "//trim(message))
         write (message,'(3(e12.6,1x))') positionComoving(2,:)
         call displayMessage("  end: "//trim(message))
         call displayUnindent("")
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

  subroutine positionInterpolatedOrphanOffset(orbit,time,position,velocity)
    !!{RST
    Compute the offset from its host (in physical position, and in velocity) of an orphan on a decaying orbit. The orbit is that of
    a satellite on a circular orbit in a singular isothermal sphere, decaying under dynamical friction :cite:p:`binney_galactic_2008`,
    such that the square of the radius decreases linearly in time to zero at the time of merging, while the tangential speed is
    constant. The ``orbit`` array contains the initial radial and tangential unit vectors, the initial radius, the initial time, the
    time of merging, and the tangential speed.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    double precision, intent(in   ), dimension(10) :: orbit
    double precision, intent(in   )                :: time
    double precision, intent(  out), dimension( 3) :: position         , velocity
    double precision                               :: fractionRemaining, duration, &
         &                                            radius           , angle   , &
         &                                            speedRadial

    ! Find the fraction of the decay time remaining.
    duration         =+orbit(9)-orbit(8)
    fractionRemaining=+max(0.0d0,min(1.0d0,(orbit(9)-time)/duration))
    if (fractionRemaining <= 0.0d0) then
       ! The orbit has decayed completely.
       position=0.0d0
       velocity=0.0d0
       return
    end if
    ! Compute the radius, angle, and radial speed. The angle is the integral over time of the angular speed, (tangential
    ! speed)/(radius). Speeds are converted from km/s to Mpc/Gyr where necessary.
    radius     =+orbit(7)                  &
         &      *sqrt(fractionRemaining)
    angle      =+2.0d0                     &
         &      *orbit(10)                 &
         &      /MpcPerKmPerSToGyr         &
         &      *duration                  &
         &      /orbit(7)                  &
         &      *(                         &
         &        +1.0d0                   &
         &        -sqrt(fractionRemaining) &
         &       )
    speedRadial=-0.5d0                     &
         &      *orbit(7)                  &
         &      /duration                  &
         &      /sqrt(fractionRemaining)   &
         &      *MpcPerKmPerSToGyr
    ! Construct the position and velocity.
    position=+radius                       &
         &   *(                            &
         &     +orbit(1:3)*cos(angle)      &
         &     +orbit(4:6)*sin(angle)      &
         &    )
    velocity=+speedRadial                  &
         &   *(                            &
         &     +orbit(1:3)*cos(angle)      &
         &     +orbit(4:6)*sin(angle)      &
         &    )                            &
         &   +orbit(10)                    &
         &   *(                            &
         &     -orbit(1:3)*sin(angle)      &
         &     +orbit(4:6)*cos(angle)      &
         &    )
    return
  end subroutine positionInterpolatedOrphanOffset
