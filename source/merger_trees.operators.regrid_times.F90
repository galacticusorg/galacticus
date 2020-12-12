!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Contains a module which implements a merger tree operator which restructures the tree onto a fixed grid of timesteps.

  use :: Cosmological_Density_Field, only : criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Linear_Growth             , only : linearGrowthClass

  !# <mergerTreeOperator name="mergerTreeOperatorRegridTimes">
  !#  <description>
  !#   A merger tree operator class which will interpolate the merger tree structure onto a new array of timesteps. The timestep
  !#   array is specified via the parameters:
  !#   \begin{description}
  !#   \item[{\normalfont \ttfamily [snapshotSpacing]}] The spacing of the timesteps. Five options are available: {\normalfont
  !#   \ttfamily linear} will space timesteps uniformly in expansion factor, {\normalfont \ttfamily logarithmic} will space
  !#   timesteps uniformly in the logarithm of expansion factor, {\normalfont \ttfamily logCriticalDensity} will space timesteps
  !#   uniformly in the logarithm of critical density, $\delta_\mathrm{c}$, {\normalfont \ttfamily millennium} will use times
  !#   corresponding to the redshifts of snapshots in the Millennium Simulation database, while {\normalfont \ttfamily read} will
  !#   use times corresponding to the redshifts specified in the {\normalfont \ttfamily mergerTreeRegridRedshifts} parameter;
  !#   \item[{\normalfont \ttfamily [expansionFactorStart]}] The smallest expansion factor in the array (ignored for {\normalfont
  !#   \ttfamily millennium} and {\normalfont \ttfamily read} spacings);
  !#   \item[{\normalfont \ttfamily [expansionFactorEnd]}] The largest expansion factor in the array (ignored for {\normalfont
  !#   \ttfamily millennium} and {\normalfont \ttfamily read} spacings);
  !#   \item[{\normalfont \ttfamily [regridCount]}] The number of timesteps in the array;
  !#   \end{description}
  !#   Along each branch of the tree, new halos are inserted at times corresponding to the times in the resulting array. The
  !#   masses of these nodes are linearly interpolated between the existing nodes on the branch. Once these new nodes have been
  !#   added, all other nodes are removed from the tree\footnote{The base node of the tree is never removed, even if it does not
  !#   lie on one of the times in the constructed array.} The processing is useful to construct representations of trees as they
  !#   would be if only sparse time sampling were available. As such, it is useful for exploring how the number of snapshots in
  !#   merger trees extracted from N-body simulations\index{merger tree!N-body} affects the properties of galaxies that form in
  !#   them.
  !#  </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorRegridTimes
     !% A merger tree operator class which restructures the tree onto a fixed grid of timesteps.
     private
     class           (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_  => null()
     class           (criticalOverdensityClass), pointer                   :: criticalOverdensity_ => null()
     class           (linearGrowthClass       ), pointer                   :: linearGrowth_        => null()
     logical                                                               :: dumpTrees
     double precision                                                      :: snapTolerance
     double precision                          , allocatable, dimension(:) :: timeGrid
   contains
     final     ::                        regridTimesDestructor
     procedure :: operatePreEvolution => regridTimesOperatePreEvolution
  end type mergerTreeOperatorRegridTimes

  interface mergerTreeOperatorRegridTimes
     !% Constructors for the regrid times merger tree operator class.
     module procedure regridTimesConstructorParameters
     module procedure regridTimesConstructorInternal
  end interface mergerTreeOperatorRegridTimes

  ! Enumeration for snapshot spacing.
  !# <enumeration>
  !#  <name>snapshotSpacing</name>
  !#  <description>Specifies the spacing of snapshots when regridding merger trees.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <entry label="linear"                 />
  !#  <entry label="logarithmic"            />
  !#  <entry label="logCriticalOverdensity" />
  !#  <entry label="millennium"             />
  !#  <entry label="list"                   />
  !# </enumeration>

contains

  function regridTimesConstructorParameters(parameters) result(self)
    !% Constructor for the regrid times merger tree operator class which takes a parameter set as input.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (mergerTreeOperatorRegridTimes)                              :: self
    type            (inputParameters              ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass      ), pointer                     :: cosmologyFunctions_
    class           (criticalOverdensityClass     ), pointer                     :: criticalOverdensity_
    class           (linearGrowthClass            ), pointer                     :: linearGrowth_
    double precision                               , allocatable  , dimension(:) :: snapshotTimes
    logical                                                                      :: dumpTrees
    integer                                                                      :: regridCount          , snapshotSpacing   , &
         &                                                                          iTime
    double precision                                                             :: expansionFactorStart , expansionFactorEnd, &
         &                                                                          snapTolerance
    type            (varying_string               )                              :: snapshotSpacingText

    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="criticalOverdensity" name="criticalOverdensity_" source="parameters"/>
    !# <objectBuilder class="linearGrowth"        name="linearGrowth_"        source="parameters"/>
    !# <inputParameter>
    !#   <name>dumpTrees</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not to dump merger trees as they are regridded.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>snapTolerance</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The fractional tolerance used in deciding if a node should be snapped to a time on the grid.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>regridCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>100</defaultValue>
    !#   <description>Number of points in time to use when regridding merger trees.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>expansionFactorStart</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>Starting expansion factor to use when regridding merger trees.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>expansionFactorEnd</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>Ending expansion factor to use when regridding merger trees.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>snapshotSpacing</name>
    !#   <variable>snapshotSpacingText</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('logarithmic')</defaultValue>
    !#   <description>Type of spacing to use in merger tree regridding (linear or logarithmic).</description>
    !# </inputParameter>
    ! Find the spacing type to be used.
    snapshotSpacing=enumerationSnapshotSpacingEncode(char(snapshotSpacingText),includesPrefix=.false.)
    ! Read redshifts if necessary.
    if (snapshotSpacing == snapshotSpacingList) then
       allocate(snapshotTimes(parameters%count('snapshotRedshifts')))
       if (size(snapshotTimes) /= regridCount) call Galacticus_Error_Report('mismatch between [regridCount] and size of [snapshotRedshifts]'//{introspection:location})
       !# <inputParameter>
       !#   <name>snapshotRedshifts</name>
       !#   <variable>snapshotTimes</variable>
       !#   <source>parameters</source>
       !#   <description>The redshifts at which merger trees are to regridded when the {\normalfont \ttfamily [snapshotSpacing]}$=${\normalfont \ttfamily list} option is selected.</description>
       !# </inputParameter>
       do iTime=1,size(snapshotTimes)
          snapshotTimes(iTime)=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(snapshotTimes(iTime)))
       end do
    end if
    ! Build the instance.
    self=mergerTreeOperatorRegridTimes(snapTolerance,regridCount,expansionFactorStart,expansionFactorEnd,snapshotSpacing,dumpTrees,snapshotTimes,cosmologyFunctions_,criticalOverdensity_,linearGrowth_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    !# <objectDestructor name="criticalOverdensity_"/>
    !# <objectDestructor name="linearGrowth_"       />
    return
  end function regridTimesConstructorParameters

  function regridTimesConstructorInternal(snapTolerance,regridCount,expansionFactorStart,expansionFactorEnd,snapshotSpacing,dumpTrees,snapshotTimes,cosmologyFunctions_,criticalOverdensity_,linearGrowth_) result(self)
    !% Internal constructor for the regrid times merger tree operator class.
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Memory_Management, only : allocateArray
    use :: Numerical_Ranges , only : Make_Range             , rangeTypeLinear, rangeTypeLogarithmic
    use :: Sorting          , only : sort
    implicit none
    type            (mergerTreeOperatorRegridTimes)                                        :: self
    integer                                        , intent(in   )                         :: regridCount
    double precision                               , intent(in   )                         :: expansionFactorStart, expansionFactorEnd, &
         &                                                                                    snapTolerance
    integer                                        , intent(in   )                         :: snapshotSpacing
    logical                                        , intent(in   )                         :: dumpTrees
    double precision                               , intent(in   ), optional, dimension(:) :: snapshotTimes
    class           (cosmologyFunctionsClass      ), intent(in   ), target                 :: cosmologyFunctions_
    class           (criticalOverdensityClass     ), intent(in   ), target                 :: criticalOverdensity_
    class           (linearGrowthClass            ), intent(in   ), target                 :: linearGrowth_
    integer                                                                                :: iTime
    !# <constructorAssign variables="dumpTrees, snapTolerance, *cosmologyFunctions_, *criticalOverdensity_, *linearGrowth_"/>

    ! Validate arguments.
    if (regridCount < 2) call Galacticus_Error_Report('regridCount > 2 is required'//{introspection:location})
    ! Construct array of grid expansion factors.
    call allocateArray(self%timeGrid,[regridCount])
    select case (snapshotSpacing)
    case (snapshotSpacingLinear                )
       self%timeGrid=Make_Range(expansionFactorStart,expansionFactorEnd,regridCount,rangeTypeLinear     )
       ! Convert expansion factors to time.
       do iTime=1,regridCount
          self%timeGrid(iTime)=self%cosmologyFunctions_%cosmicTime(self%timeGrid(iTime))
       end do
    case (snapshotSpacingLogarithmic           )
       self%timeGrid=Make_Range(expansionFactorStart,expansionFactorEnd,regridCount,rangeTypeLogarithmic)
       ! Convert expansion factors to time.
       do iTime=1,regridCount
          self%timeGrid(iTime)=self%cosmologyFunctions_%cosmicTime(self%timeGrid(iTime))
       end do
    case (snapshotSpacingLogCriticalOverdensity)
       ! Build a logarithmic grid in critical overdensity.
       self%timeGrid=Make_Range(                                                                                             &
            &                   +self%criticalOverdensity_%value(self%cosmologyFunctions_%cosmicTime(expansionFactorStart))  &
            &                   /self%linearGrowth_       %value(self%cosmologyFunctions_%cosmicTime(expansionFactorStart)), &
            &                   +self%criticalOverdensity_%value(self%cosmologyFunctions_%cosmicTime(expansionFactorEnd  ))  &
            &                   /self%linearGrowth_       %value(self%cosmologyFunctions_%cosmicTime(expansionFactorEnd  )), &
            &                    regridCount                                                                               , &
            &                    rangeTypeLogarithmic                                                                        &
            &                  )
       ! Convert critical overdensity to time.
       do iTime=1,regridCount
          self%timeGrid(iTime)=self%criticalOverdensity_%timeOfCollapse(self%timeGrid(iTime))
       end do
    case (snapshotSpacingMillennium            )
       ! Use the timesteps used in the original Millennium Simulation as reported by Croton et al.  (2006; MNRAS; 365;
       ! 11; http://adsabs.harvard.edu/abs/2006MNRAS.365...11C). Note that we specifically use the redshifts for these
       ! snapshots as reported by the Millennium Database using query: "select z from Snapshots..MR".
       ! Check for consistent number of timesteps.
       if (regridCount /= 60) call Galacticus_Error_Report('"millennium" grid spacing requires exactly 60 timesteps'//{introspection:location})
       ! Convert redshifts to time.
       self%timeGrid=[                                                                         &
            &         19.915688d0,18.243723d0,16.724525d0,15.343073d0,14.085914d0,12.940780d0, &
            &         11.896569d0,10.943864d0,10.073462d0, 9.277915d0, 8.549912d0, 7.883204d0, &
            &          7.272188d0, 6.711586d0, 6.196833d0, 5.723864d0, 5.288833d0, 4.888449d0, &
            &          4.519556d0, 4.179469d0, 3.865683d0, 3.575905d0, 3.308098d0, 3.060419d0, &
            &          2.831182d0, 2.618862d0, 2.422044d0, 2.239486d0, 2.070027d0, 1.912633d0, &
            &          1.766336d0, 1.630271d0, 1.503636d0, 1.385718d0, 1.275846d0, 1.173417d0, &
            &          1.077875d0, 0.988708d0, 0.905463d0, 0.827699d0, 0.755036d0, 0.687109d0, &
            &          0.623590d0, 0.564177d0, 0.508591d0, 0.456577d0, 0.407899d0, 0.362340d0, &
            &          0.319703d0, 0.279802d0, 0.242469d0, 0.207549d0, 0.174898d0, 0.144383d0, &
            &          0.115883d0, 0.089288d0, 0.064493d0, 0.041403d0, 0.019933d0, 0.000000d0  &
            &        ]
       do iTime=1,regridCount
          self%timeGrid(iTime)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%timeGrid(iTime)))
       end do
    case (snapshotSpacingList                  )
       if (.not.present(snapshotTimes)) call Galacticus_Error_Report('"list" grid spacing requires a list of snapshot times be supplied'//{introspection:location})
       self%timeGrid=snapshotTimes
       call sort(self%timeGrid)
    end select
    return
  end function regridTimesConstructorInternal

  subroutine regridTimesDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorRegridTimes), intent(inout) :: self

    !# <objectDestructor name="self%linearGrowth_"       />
    !# <objectDestructor name="self%criticalOverdensity_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    return
  end subroutine regridTimesDestructor

  subroutine regridTimesOperatePreEvolution(self,tree)
    !% Perform a regrid times operation on a merger tree.
    use            :: Galacticus_Display     , only : Galacticus_Display_Indent, Galacticus_Display_Unindent  , Galacticus_Display_Message, verbosityWorking
    use            :: Galacticus_Error       , only : Galacticus_Error_Report  , Galacticus_Warn
    use            :: Galacticus_Nodes       , only : mergerTree               , nodeComponentBasic           , nodeComponentSatellite    , nodeEvent       , &
          &                                           treeNode                 , treeNodeList
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: ISO_Varying_String     , only : var_str
    use            :: Kind_Numbers           , only : kind_int8
    use            :: Merger_Tree_Walkers    , only : mergerTreeWalkerAllNodes , mergerTreeWalkerIsolatedNodes
    use            :: Merger_Trees_Dump      , only : Merger_Tree_Dump
    use            :: Numerical_Comparison   , only : Values_Agree
    use            :: Numerical_Interpolation, only : interpolator
    use            :: String_Handling        , only : operator(//)
    implicit none
    class           (mergerTreeOperatorRegridTimes), intent(inout), target                :: self
    type            (mergerTree                   ), intent(inout), target                :: tree
    type            (treeNode                     )                             , pointer :: nodeChild                       , mergee     , &
         &                                                                                   nodeSibling                     , node
    type            (treeNodeList                 ), allocatable  , dimension(:)          :: newNodes
    integer         (kind=kind_int8               ), allocatable  , dimension(:)          :: highlightNodes
    class           (nodeComponentBasic           )                             , pointer :: basicChild                      , basicParent, &
         &                                                                                   basic
    class           (nodeComponentSatellite       )                             , pointer :: mergeeSatellite
    type            (mergerTree                   )                             , pointer :: currentTree
    class           (nodeEvent                    )                             , pointer :: event                           , pairedEvent
    type            (mergerTreeWalkerAllNodes     )                                       :: treeWalkerAllNodes
    type            (mergerTreeWalkerIsolatedNodes)                                       :: treeWalkerIsolatedNodes
    logical                                                                               :: mergeTargetWarningIssued=.false.
    type            (interpolator                 )                                       :: interpolator_
    integer         (c_size_t                     )                                       :: iNow                            , iParent    , &
         &                                                                                   iTime                           , countNodes
    integer                                                                               :: allocErr
    double precision                                                                      :: massNow                         , massParent , &
         &                                                                                   timeNow                         , timeParent
    integer         (kind=kind_int8               )                                       :: firstNewNode                    , nodeIndex

    ! Build an interpolator.
    interpolator_=interpolator(self%timeGrid)
    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       call Galacticus_Display_Indent(var_str('Regridding tree ')//currentTree%index,verbosityWorking)
       ! Dump the unprocessed tree if required.
       if (self%dumpTrees) call Merger_Tree_Dump(                              &
            &                                    currentTree                 , &
            &                                    backgroundColor    ='white' , &
            &                                    nodeColor          ='black' , &
            &                                    highlightColor     ='black' , &
            &                                    edgeColor          ='black' , &
            &                                    nodeStyle          ='solid' , &
            &                                    highlightStyle     ='filled', &
            &                                    edgeStyle          ='solid' , &
            &                                    labelNodes         =.false. , &
            &                                    scaleNodesByLogMass=.true.  , &
            &                                    edgeLengthsToTimes =.true.    &
            &                                   )
       ! Iterate through to tree to:
       !  a) Find the current maximum node index in the tree, and;
       !  b) Snap halos to snapshot times if requested, and;
       !  c) Count the number of nodes in the tree.
       nodeIndex =0_kind_int8
       countNodes=0_c_size_t
       treeWalkerAllNodes=mergerTreeWalkerAllNodes(currentTree)
       do while (treeWalkerAllNodes%next(node))
          nodeIndex=max(nodeIndex,node%index())
          ! Count node.
          countNodes=countNodes+1_c_size_t
          ! Check for merge targets being set - these are not supported under the regridding transformation, so issue a warning.
          if (associated(node%mergeTarget).and..not.mergeTargetWarningIssued) then
             !$omp critical (mergeTargetWarning)
             if (.not.mergeTargetWarningIssued) then
                call Galacticus_Warn(                                                                                &
                     &               'WARNING: nodes in this tree have merge targets set'               //char(10)// &
                     &               '         this is not supported by the regridding operator'        //char(10)// &
                     &               '         your tree may crash or deadlock'                         //char(10)// &
                     &               '         to avoid this problem do not preset merge targets, e.g. '//char(10)// &
                     &               '           <mergerTreeReadPresetMergerNodes value="false"/>'                   &
                     &              )
                mergeTargetWarningIssued=.true.
             end if
             !$omp end critical (mergeTargetWarning)
          end if
          ! Check if this node can be snapped to a grid time.
          if (self%snapTolerance > 0.0d0) then
             ! Get the basic component.
             basic   => node %basic()
             ! Get the time for this node.
             timeNow =  basic%time ()
             ! Find the closest time in the new time grid.
             iNow    =  interpolator_%locate(timeNow,closest=.true.)
             ! Test how close the node is to this time.
             if (Values_Agree(timeNow,self%timeGrid(iNow),relTol=self%snapTolerance)) then
                ! Adjust the time of the node.
                call basic%timeSet(self%timeGrid(iNow))
                ! Check for mergees with their merge times set to match the time of this node.
                mergee => node%firstMergee
                do while (associated(mergee))
                   mergeeSatellite => mergee         %satellite    ()
                   ! Get the merge time for this mergee.
                   timeNow         =  mergeeSatellite%timeOfMerging()
                   ! Find the closest time in the new time grid.
                   iNow    =  interpolator_%locate(timeNow,closest=.true.)
                   if (Values_Agree(timeNow,self%timeGrid(iNow),relTol=self%snapTolerance)) &
                        & call mergeeSatellite%timeOfMergingSet(self%timeGrid(iNow))
                   mergee => mergee%siblingMergee
                end do
                ! Check for events with their event times set to match the time of this node.
                event => node%event
                do while (associated(event))
                   ! Get the merge time for this event.
                   timeNow=event%time
                   ! Find the closest time in the new time grid.
                   iNow   =interpolator_%locate(timeNow,closest=.true.)
                   if (Values_Agree(timeNow,self%timeGrid(iNow),relTol=self%snapTolerance)) then
                      event%time=self%timeGrid(iNow)
                      if (associated(event%node)) then
                         pairedEvent => event%node%event
                         do while (associated(pairedEvent))
                            if (pairedEvent%ID == event%ID) then
                               pairedEvent%time=self%timeGrid(iNow)
                               exit
                            end if
                            pairedEvent => pairedEvent%next
                         end do
                      end if
                   end if
                   event => event%next
                end do
             end if
          end if
       end do
       firstNewNode=nodeIndex+1
       call Galacticus_Display_Message(var_str('Tree contains ')//countNodes//' nodes prior to regridding',verbosityWorking)
       ! Walk the tree, locating branches which intersect grid times.
       treeWalkerIsolatedNodes=mergerTreeWalkerIsolatedNodes(currentTree)
       do while (treeWalkerIsolatedNodes%next(node))
          ! Skip this node if it is the root node.
          if (associated(node%parent)) then
             basic       => node       %basic()
             basicParent => node%parent%basic()
             ! Get the time of this node and its parent.
             timeNow   =basic  %time()
             timeParent=basicParent%time()
             ! Get masses of these halos.
             massNow   =basic  %mass()
             massParent=basicParent%mass()
             if (node%isPrimaryProgenitor()) then
                ! Remove the mass in any non-primary progenitors - we don't want to include
                ! their mass in the estimated mass growth rate of this node.
                nodeChild => node%parent%firstChild%sibling
                do while (associated(nodeChild))
                   basicChild => nodeChild%basic()
                   massParent =  massParent-basicChild%mass()
                   nodeChild  => nodeChild%sibling
                end do
                ! Do not let the parent mass decrease along the branch.
                massParent=max(massParent,massNow)
             else
                ! Halo is not the primary progenitor of its parent. Assume that its mass does
                ! not grow further.
                massParent=massNow
             end if
             ! Locate these times in the list of grid times.
             iNow   =interpolator_%locate(timeNow   )
             iParent=interpolator_%locate(timeParent)
             ! For nodes existing precisely at a grid time, ignore this grid point. (These are,
             ! typically, nodes which have been created at these points.)
             if (timeParent == self%timeGrid(iParent)) iParent=iParent-1
             ! If the branch from node to parent spans one or more grid times, insert new nodes
             ! at those points.
             if (iParent > iNow) then
                ! Create new nodes.
                allocate(newNodes(iParent-iNow),stat=allocErr)
                if (allocErr/=0) call Galacticus_Error_Report('unable to allocate new nodes'//{introspection:location})
                do iTime=iNow+1,iParent
                   nodeIndex=nodeIndex+1_kind_int8
                   newNodes(iTime-iNow)%node => treeNode(hostTree=currentTree)
                   call newNodes(iTime-iNow)%node%indexSet(nodeIndex)
                end do
                ! Assign node properties and build links.
                do iTime=iNow+1,iParent
                   ! Assign a time and a mass
                   basic => newNodes(iTime-iNow)%node%basic(autoCreate=.true.)
                   call basic%timeSet(                              self%timeGrid(iTime)                              )
                   call basic%massSet(massNow+(massParent-massNow)*(self%timeGrid(iTime)-timeNow)/(timeParent-timeNow))
                   ! Link to child node.
                   if (iTime > iNow+1 ) newNodes(iTime-iNow)%node%firstChild => newNodes(iTime-iNow-1)%node
                   ! Link to parent node.
                   if (iTime < iParent) newNodes(iTime-iNow)%node%parent     => newNodes(iTime-iNow+1)%node
                end do
                ! Link final node to the parent.
                newNodes(iParent-iNow)%node%parent  => node%parent
                ! Link final node sibling to current node sibling.
                newNodes(iParent-iNow)%node%sibling => node%sibling
                ! Link the parent to the final node.
                if (node%isPrimaryProgenitor()) then
                   ! Node is the main progenitor of its parent, so simply replace it with the
                   ! final node in our list.
                   node%parent%firstChild  => newNodes(iParent-iNow)%node
                else
                   ! Node is not the main progenitor of its parent, so find the child node that
                   ! has it as a sibling.
                   nodeChild => node%parent%firstChild
                   do while (.not.associated(nodeChild%sibling,node))
                      nodeChild => nodeChild%sibling
                   end do
                   nodeChild%sibling => newNodes(iParent-iNow)%node
                end if
                ! Link the child of the first node to the node being processed.
                newNodes(1)%node%firstChild  => node
                ! Nullify any sibling of the node being processed.
                node%sibling => null()
                ! Link the parent of the node being processed to the first node of the list.
                node%parent  => newNodes(1)%node
                ! Erase the node list.
                deallocate(newNodes)
             end if
          end if
       end do
       ! Dump the intermediate tree if required.
       if (self%dumpTrees) then
          allocate(highlightNodes(nodeIndex-firstNewNode+2))
          highlightNodes(1)=currentTree%baseNode%index()
          do nodeIndex=1,nodeIndex-firstNewNode+1
             highlightNodes(nodeIndex+1)=firstNewNode+nodeIndex-1
          end do
          call Merger_Tree_Dump(                                    &
               &                currentTree                       , &
               &                highlightNodes     =highlightNodes, &
               &                backgroundColor    ='white'       , &
               &                nodeColor          ='black'       , &
               &                highlightColor     ='black'       , &
               &                edgeColor          ='black'       , &
               &                nodeStyle          ='solid'       , &
               &                highlightStyle     ='filled'      , &
               &                edgeStyle          ='dotted'      , &
               &                labelNodes         =.false.       , &
               &                scaleNodesByLogMass=.true.        , &
               &                edgeLengthsToTimes =.true.          &
               &               )
          deallocate(highlightNodes)
       end if
       ! Walk the tree removing nodes not at grid times.
       countNodes              =0_c_size_t
       treeWalkerIsolatedNodes=mergerTreeWalkerIsolatedNodes(currentTree)
       do while (treeWalkerIsolatedNodes%next(node))
          basic => node%basic()
          ! Get the time for this node.
          timeNow=basic%time()
          ! Find the closest time in the new time grid.
          iNow   =interpolator_%locate(timeNow,closest=.true.)
          ! If this node does not lie precisely on the grid then remove it.
          if (associated(node%parent) .and. timeNow /= self%timeGrid(iNow)) then
             if (node%isPrimaryProgenitor()) then
                ! Handle primary progenitor nodes.
                if (associated(node%firstChild)) then
                   ! Handle primary progenitors with children
                   nodeChild => node%firstChild
                   ! Assign all children a parent that is the parent of the current node.
                   do while (associated(nodeChild))
                      nodeChild%parent => node %parent
                      if (.not.associated(nodeChild%sibling)) then
                         nodeChild%sibling => node%sibling
                         nodeChild             => null()
                      else
                         nodeChild             => nodeChild%sibling
                      end if
                   end do
                   ! Assign the current node's parent a child that is the child of the current node.
                   node%parent%firstChild => node%firstChild
                else
                   ! Handle primary nodes with no children - simply make the parent's main
                   ! progenitor the sibling of the current node.
                   node%parent%firstChild => node%sibling
                end if
             else
                ! Handle non-primary nodes.
                if (associated(node%firstChild)) then
                   ! Handle non-primary nodes with children.
                   ! Assign all children a parent that is the parent of the current node.
                   nodeChild => node%firstChild
                   do while (associated(nodeChild))
                      nodeChild%parent => node %parent
                      if (.not.associated(nodeChild%sibling)) then
                         nodeChild%sibling => node%sibling
                         nodeChild => null()
                      else
                         nodeChild            => nodeChild%sibling
                      end if
                   end do
                   ! Find which sibling points the current node and link in the children of the
                   ! current node.
                   nodeSibling => node%parent%firstChild
                   do while (.not.associated(nodeSibling%sibling,node))
                      nodeSibling => nodeSibling%sibling
                   end do
                   nodeSibling%sibling => node%firstChild
                else
                   ! Handle non-primary nodes with no children - just snip it out of the sibling list.
                   nodeSibling => node%parent%firstChild
                   do while (.not.associated(nodeSibling%sibling,node))
                      nodeSibling => nodeSibling%sibling
                   end do
                   nodeSibling%sibling => node%sibling
                end if
             end if
             ! Destroy the node.
             call node%destroy()
             deallocate(node)
             call treeWalkerIsolatedNodes%previous(node)
          else
             countNodes=countNodes+1_c_size_t
          end if
       end do
       call Galacticus_Display_Message(var_str('Tree contains ')//countNodes//' nodes after regridding',verbosityWorking)
       ! Dump the processed tree if required.
       if (self%dumpTrees) call Merger_Tree_Dump(                              &
            &                                    currentTree                 , &
            &                                    backgroundColor    ='white' , &
            &                                    nodeColor          ='black' , &
            &                                    highlightColor     ='black' , &
            &                                    edgeColor          ='black' , &
            &                                    nodeStyle          ='solid' , &
            &                                    highlightStyle     ='filled', &
            &                                    edgeStyle          ='solid' , &
            &                                    labelNodes         =.false. , &
            &                                    scaleNodesByLogMass=.true.  , &
            &                                    edgeLengthsToTimes =.true.    &
            &                                   )
       call Galacticus_Display_Unindent('Done',verbosityWorking)
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine regridTimesOperatePreEvolution
