!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a preset position component in which positions are interpolated using the approach of
\cite{merson_lightcone_2013}.
!!}

module Node_Component_Position_Preset_Interpolated
  !!{
  Implements a preset position component in which positions are interpolated using the approach of
  \cite{merson_lightcone_2013}.
  !!}
  implicit none
  private
  public :: threadInitialize                                   , threadUninitialize                               , &
       &    computeInterpolation                               , rateCompute                                      , &
       &    initialize                                         , nodeComponentPositionPresetInterpolatedStateStore, &
       &    nodeComponentPositionPresetInterpolatedStateRestore
  
  !![
  <component>
   <class>position</class>
   <name>presetInterpolated</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>position</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <getFunction bindsTo="component">PositionPresetInterpolatedPosition</getFunction>
      <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Position of the node (in physical coordinates)."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <getFunction bindsTo="component">PositionPresetInterpolatedVelocity</getFunction>
      <output labels="[X,Y,Z]" unitsInSI="kilo" comment="Velocity of the node (in physical coordinates)."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>positionHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>interpolationCoefficients</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>interpolationTimeMaximum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
   <functions>objects.nodes.components.position.preset.interpolated.bound_functions.inc</functions>
  </component>
  !!]

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine initialize(parameters_)
    !!{
    Initializes the ``preset-interpolated'' position component module.
    !!}
    use :: Galacticus_Nodes                                , only : defaultPositionComponent
    use :: Input_Parameters                                , only : inputParameter          , inputParameters
    use :: Node_Component_Position_Preset_Interpolated_Data, only : positionPresetBoxLength , isPeriodic
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    
    !$omp critical (nodeComponentPositionPresetInterpolatedInitialize)
    if (defaultPositionComponent%presetInterpolatedIsActive()) then
       !![
       <inputParameter>
         <name>positionPresetBoxLength</name>
         <defaultValue>0.0d0</defaultValue>
         <description>The periodic length of the positions. For non-periodic positions, a value of zero should be given.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]
       isPeriodic=positionPresetBoxLength > 0.0d0
    end if
    !$omp end critical (nodeComponentPositionPresetInterpolatedInitialize)
    return
  end subroutine initialize
  
  !![
  <nodeComponentThreadInitializationTask>
   <unitName>threadInitialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine threadInitialize(parameters_)
    !!{
    Initializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks                                    , only : nodePromotionEvent      , postEvolveEvent, openMPThreadBindingAtLevel
    use :: Galacticus_Nodes                                , only : defaultPositionComponent
    use :: Input_Parameters                                , only : inputParameters
    use :: Node_Component_Position_Preset_Interpolated_Data, only : cosmologyFunctions_
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultPositionComponent%presetInterpolatedIsActive()) then
       !![
       <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters_"/>
       !!]
       call nodePromotionEvent%attach(defaultPositionComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentPositionPresetInterpolated')
       call    postEvolveEvent%attach(defaultPositionComponent,postEvolve   ,openMPThreadBindingAtLevel,label='nodeComponentPositionPresetInterpolated')
    end if
    return
  end subroutine threadInitialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>threadUninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine threadUninitialize()
    !!{
    Uninitializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks                                    , only : nodePromotionEvent      , postEvolveEvent
    use :: Galacticus_Nodes                                , only : defaultPositionComponent
    use :: Node_Component_Position_Preset_Interpolated_Data, only : cosmologyFunctions_
    implicit none

    if (defaultPositionComponent%presetInterpolatedIsActive()) then
       !![
       <objectDestructor name="cosmologyFunctions_"/>
       !!]
       if (nodePromotionEvent%isAttached(defaultPositionComponent,nodePromotion)) call nodePromotionEvent%detach(defaultPositionComponent,nodePromotion)
       if (   postEvolveEvent%isAttached(defaultPositionComponent,postEvolve   )) call    postEvolveEvent%detach(defaultPositionComponent,postEvolve   )
    end if
    return
  end subroutine threadUninitialize

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, update the position of {\normalfont \ttfamily
    node} to that of the parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition, nodeComponentPositionPresetInterpolated, treeNode
    implicit none
    class(*                    ), intent(inout)          :: self
    type (treeNode             ), intent(inout), target  :: node
    class(nodeComponentPosition)               , pointer :: positionParent, position
    !$GLC attributes unused :: self
    
    position       => node       %position()
    positionParent => node%parent%position()
    select type (positionParent)
    class is (nodeComponentPositionPresetInterpolated)
       call position%                 positionSet(positionParent%position                 ())
       call position%                 velocitySet(positionParent%velocity                 ())
       call position%          positionHistorySet(positionParent%positionHistory          ())
       call position%interpolationCoefficientsSet(positionParent%interpolationCoefficients())
       call position% interpolationTimeMaximumSet(positionParent%interpolationTimeMaximum ())
    end select
    return
  end subroutine nodePromotion

  subroutine postEvolve(self,node)
    !!{
    Trigger interpoltion recalculation after an evolution step.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentPositionPresetInterpolated, treeNode
    implicit none
    class(*       ), intent(inout)          :: self
    type (treeNode), intent(inout), target  :: node

    select type (self)
    class is (nodeComponentPositionPresetInterpolated)
       call computeInterpolation(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine postEvolve

  !![
  <nodeMergerTask>
   <unitName>computeInterpolation</unitName>
  </nodeMergerTask>
  <satelliteHostChangeTask>
   <unitName>computeInterpolation</unitName>
  </satelliteHostChangeTask>
  <mergerTreeInitializeTask>
   <unitName>computeInterpolation</unitName>
  </mergerTreeInitializeTask>
  !!]
  subroutine computeInterpolation(node)
    !!{
    Compute interpolation coefficients for positions. The approach here follows that of \cite{merson_lightcone_2013}. For halos
    which are not orbiting within a host halo during the current interval, interpolation is via a cubic polynomial in each
    Cartesian coordinate of \emph{comoving} position matched to the position and velocity at the initial and final times. For
    halos which are orbiting within a halo halo during the current interval, interpolation is via a logarithmic spiral in the
    relative physical position of the halo and the center of the host halo.
    !!}
    use, intrinsic :: ISO_C_Binding                                   , only : c_size_t
    use            :: Error                                           , only : Error_Report    
    use            :: Galacticus_Nodes                                , only : nodeComponentBasic       , nodeComponentPosition  , treeNode                          , nodeEvent                   , &
         &                                                                     nodeEventSubhaloPromotion, nodeEventBranchJump    , nodeEventSubhaloPromotionIntertree, nodeEventBranchJumpIntertree, &
         &                                                                     defaultPositionComponent
    use            :: Node_Component_Position_Preset_Interpolated_Data, only : cosmologyFunctions_      , positionPresetBoxLength, isPeriodic
    use            :: Numerical_Constants_Astronomical                , only : Mpc_per_km_per_s_To_Gyr
    use            :: Linear_Algebra                                  , only : vector                   , matrix                 , assignment(=)
    use            :: Histories                                       , only : history
    use            :: Numerical_Comparison                            , only : Values_Differ            , Values_Agree           , Values_Less_Than
    use            :: Numerical_Interpolation                         , only : interpolator
    use            :: Vectors                                         , only : Vector_Product           , Vector_Magnitude
    implicit none
    type            (treeNode             ), intent(inout)    , target  :: node
    class           (nodeComponentPosition)                   , pointer :: position                     , positionParent       , &
         &                                                                 positionGrandParent
    class           (nodeComponentBasic   )                   , pointer :: basic                        , basicParent          , &
         &                                                                 basicGrandParent
    type            (treeNode             )                   , pointer :: nodeParent                   , nodeGrandparent      , &
         &                                                                 nodeJump
    class           (nodeEvent            )                   , pointer :: event
    double precision                       , dimension(  2, 3)          :: positionComoving             , velocityComoving     , &
         &                                                                 positionRelative
    double precision                       , dimension(2,2, 3)          :: vectorInPlaneNormal
    double precision                       , dimension(  2   )          :: time                         , expansionFactor
    double precision                       , dimension(  2, 2)          :: coefficientsAngle            , coefficientsLogRadius
    double precision                       , dimension(     3)          :: vectorNormal
    double precision                       , dimension(  4, 3)          :: coefficientsCubic
    double precision                       , dimension(    20)          :: coefficientsSpiral
    double precision                       , dimension(     0)          :: coefficientsNull
    type            (vector               ), allocatable                :: coordinates                  , coefficients
    type            (matrix               ), allocatable                :: terms
    double precision                       , parameter                  :: separationTiny        =1.0d-6
    type            (history              )                             :: positionHistory
    integer         (c_size_t             )                             :: i                            , j                    , &
         &                                                                 k
    type            (interpolator         )                             :: interpolator_
    double precision                                                    :: timeParent                   , timeGrandparent
    logical                                                             :: useSpiralInterpolation       , searchEvent

    ! If this component is not active return immediately.
    if (.not.defaultPositionComponent%presetInterpolatedIsActive()) return
    ! Get the position component of this node in which we will store the interpolation coefficients.
    position => node%position()      
    ! If the node has no parent, there is no interpolation to do, so return.
    if (.not.associated(node%parent)) then
       ! Set a null interpolation - this will caused the fixed-at-snapshot position to be used instead.
       call position%interpolationCoefficientsSet(coefficientsNull)
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
             call position%interpolationCoefficientsSet(coefficientsNull)
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
             call position%interpolationCoefficientsSet(coefficientsNull)
             call position%interpolationTimeMaximumSet (huge(0.0d0))
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
    call position%interpolationCoefficientsSet(coefficientsNull)
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
          call position%interpolationCoefficientsSet(coefficientsNull)
          call position%interpolationTimeMaximumSet (huge(0.0d0))
          return
       end if
       ! Compute the interpolation only if there is a non-zero time interval over which to interpolate. Halos about to experience
       ! a node promotion or branch jump event can exist at precisely the time of the halo they are about to become.
       if (time(1) < time(2)) then
          ! Get comoving position/velocity at start/end times.
          expansionFactor (1  )=                      cosmologyFunctions_%expansionFactor(time(1))
          expansionFactor (2  )=                      cosmologyFunctions_%expansionFactor(time(2))
          positionComoving(1,:)=positionComoving(1,:)/                    expansionFactor(     1 )
          positionComoving(2,:)=positionComoving(2,:)/                    expansionFactor(     2 )
          velocityComoving(1,:)=velocityComoving(1,:)/                    expansionFactor(     1 )/Mpc_per_km_per_s_To_Gyr
          velocityComoving(2,:)=velocityComoving(2,:)/                    expansionFactor(     2 )/Mpc_per_km_per_s_To_Gyr
          ! Handle periodic positions.
          if (isPeriodic) then
             do j=1,3
                if (positionComoving(2,j) > positionComoving(1,j)+0.5d0*positionPresetBoxLength) positionComoving(2,j)=positionComoving(2,j)-positionPresetBoxLength
                if (positionComoving(2,j) < positionComoving(1,j)-0.5d0*positionPresetBoxLength) positionComoving(2,j)=positionComoving(2,j)+positionPresetBoxLength
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
          call position%interpolationCoefficientsSet(reshape(coefficientsCubic,[12]))
          call position%interpolationTimeMaximumSet (   time(                    2 ))
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
                call position%interpolationCoefficientsSet(coefficientsNull)
                call position%interpolationTimeMaximumSet (huge(0.0d0))
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
             if (k == 1 .and. isPeriodic) then
                expansionFactor(1)=cosmologyFunctions_%expansionFactor(time(1))
                expansionFactor(2)=cosmologyFunctions_%expansionFactor(time(2))
                do j=1,3
                   if (positionRelative(2,j)/expansionFactor(2) > positionRelative(1,j)/expansionFactor(1)+0.5d0*positionPresetBoxLength) positionRelative(2,j)=positionRelative(2,j)-positionPresetBoxLength*expansionFactor(2)
                   if (positionRelative(2,j)/expansionFactor(2) < positionRelative(1,j)/expansionFactor(1)-0.5d0*positionPresetBoxLength) positionRelative(2,j)=positionRelative(2,j)+positionPresetBoxLength*expansionFactor(2)
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
             coefficientsLogRadius(k,1)=-time(1) &
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
          call position%interpolationCoefficientsSet(coefficientsSpiral   )
          call position%interpolationTimeMaximumSet (time              (2))
       else
          ! No grandparent halo exists. This can occur in a tree which ceases to exist. Use a null interpolation to ensure that
          ! the fixed-at-snapshot position will be used for the node. Set the maximum time for the interpolation to infinity to
          ! avoid an infinite loop of attempting to compute this interpolation.
          call position%interpolationCoefficientsSet(coefficientsNull)
          call position%interpolationTimeMaximumSet (huge(0.0d0))
       end if
    end if
    return
  end subroutine computeInterpolation

  !![
  <rateComputeTask>
   <unitName>rateCompute</unitName>
  </rateComputeTask>
  !!]
  subroutine rateCompute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Interrupt evolution to update position interpolation.
    !!}
    use :: Galacticus_Nodes, only : defaultPositionComponent, interruptTask, nodeComponentBasic, nodeComponentPosition, &
         &                          treeNode
    implicit none
    type     (treeNode             ), intent(inout)          :: node
    logical                         , intent(inout)          :: interrupt
    procedure(interruptTask        ), intent(inout), pointer :: interruptProcedure
    integer                         , intent(in   )          :: propertyType
    class    (nodeComponentBasic   )               , pointer :: basic
    class    (nodeComponentPosition)               , pointer :: position
    !$GLC attributes unused :: propertyType

    ! If this component is not active return immediately.
    if (.not.defaultPositionComponent%presetInterpolatedIsActive()) return
    ! If the node has no parent, there is no interpolation to do, so return.
    basic    => node%basic   ()      
    position => node%position()      
    if (basic%time() >= position%interpolationTimeMaximum()) then
       interrupt          =  .true.
       interruptProcedure => computeInterpolation
    end if
    return
  end subroutine rateCompute

  !![
  <stateStoreTask>
   <unitName>nodeComponentPositionPresetInterpolatedStateStore</unitName>
  </stateStoreTask>
  !!]
  subroutine nodeComponentPositionPresetInterpolatedStateStore(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display                                         , only : displayMessage     , verbosityLevelInfo
    use            :: Node_Component_Position_Preset_Interpolated_Data, only : cosmologyFunctions_
    use, intrinsic :: ISO_C_Binding                                   , only : c_ptr              , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentPosition -> presetInterpolated',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="cosmologyFunctions_"/>
    !!]
    return
  end subroutine nodeComponentPositionPresetInterpolatedStateStore

  !![
  <stateRetrieveTask>
   <unitName>nodeComponentPositionPresetInterpolatedStateRestore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine nodeComponentPositionPresetInterpolatedStateRestore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display                                         , only : displayMessage     , verbosityLevelInfo
    use            :: Node_Component_Position_Preset_Interpolated_Data, only : cosmologyFunctions_
    use, intrinsic :: ISO_C_Binding                                   , only : c_ptr              , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentPosition -> presetInterpolated',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="cosmologyFunctions_"/>
    !!]
    return
  end subroutine nodeComponentPositionPresetInterpolatedStateRestore

end module Node_Component_Position_Preset_Interpolated
