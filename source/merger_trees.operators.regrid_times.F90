!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !# <mergerTreeOperator name="mergerTreeOperatorRegridTimes">
  !#  <description>Provides a merger tree operator which restructures the tree onto a fixed grid of timesteps.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorRegridTimes
     !% A merger tree operator class which restructures the tree onto a fixed grid of timesteps.
     private
     logical                                     :: dumpTrees
     double precision, allocatable, dimension(:) :: timeGrid
   contains
     final     ::            regridTimesDestructor
     procedure :: operate => regridTimesOperate
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

  function regridTimesConstructorParameters(parameters)
    !% Constructor for the regrid times merger tree operator class which takes a parameter set as input.
    use Cosmology_Functions
    implicit none
    type            (mergerTreeOperatorRegridTimes)                              :: regridTimesConstructorParameters
    type            (inputParameters              ), intent(in   )               :: parameters
    class           (cosmologyFunctionsClass      ), pointer                     :: cosmologyFunctions_
    double precision                               , allocatable  , dimension(:) :: snapshotTimes
    logical                                                                      :: dumpTrees
    integer                                                                      :: regridCount                     , snapshotSpacing   , &
         &                                                                          iTime
    double precision                                                             :: expansionFactorStart            , expansionFactorEnd
    type            (varying_string               )                              :: snapshotSpacingText
    !# <inputParameterList label="allowedParameterNames" />
        
    !# <inputParameter>
    !#   <name>dumpTrees</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not to dump merger trees as they are regridded.</description>
    !#   <type>boolean</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>regridCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>100</defaultValue>
    !#   <description>Number of points in time to use when regridding merger trees.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>expansionFactorStart</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <description>Starting expansion factor to use when regridding merger trees.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>expansionFactorEnd</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>Ending expansion factor to use when regridding merger trees.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>snapshotSpacing</name>
    !#   <variable>snapshotSpacingText</variable>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('logarithmic')</defaultValue>
    !#   <description>Type of spacing to use in merger tree regridding (linear or logarithmic).</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    ! Find the spacing type to be used.
    snapshotSpacing=enumerationSnapshotSpacingEncode(char(snapshotSpacingText),includesPrefix=.false.)
    ! Read redshifts if necessary.
    if (snapshotSpacing == snapshotSpacingList) then
       !# <inputParameter>
       !#   <name>snapshotRedshifts</name>
       !#   <variable>snapshotTimes</variable>
       !#   <source>parameters</source>
       !#   <description>The redshifts at which merger trees are to regridded when the {\normalfont \ttfamily [snapshotSpacing]}$=${\normalfont \ttfamily list} option is selected.</description>
       !#   <type>real</type>
       !#   <cardinality>0..</cardinality>
       !# </inputParameter>
       cosmologyFunctions_ => cosmologyFunctions()
       do iTime=1,size(snapshotTimes)
          snapshotTimes(iTime)=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(snapshotTimes(iTime)))
       end do
    end if
    ! Build the instance.
    regridTimesConstructorParameters=regridTimesConstructorInternal(regridCount,expansionFactorStart,expansionFactorEnd,snapshotSpacing,dumpTrees)
    return
  end function regridTimesConstructorParameters

  function regridTimesConstructorInternal(regridCount,expansionFactorStart,expansionFactorEnd,snapshotSpacing,dumpTrees,snapshotTimes)
    !% Internal constructor for the regrid times merger tree operator class.
    use Cosmology_Functions
    use Critical_Overdensities
    use Numerical_Ranges
    use Galacticus_Error
    use Memory_Management
    use Sort
    implicit none
    type            (mergerTreeOperatorRegridTimes)                                        :: regridTimesConstructorInternal
    integer                                        , intent(in   )                         :: regridCount
    double precision                               , intent(in   )                         :: expansionFactorStart          , expansionFactorEnd
    integer                                        , intent(in   )                         :: snapshotSpacing
    logical                                        , intent(in   )                         :: dumpTrees
    double precision                               , intent(in   ), optional, dimension(:) :: snapshotTimes
    class           (cosmologyFunctionsClass      ), pointer                               :: cosmologyFunctions_
    class           (criticalOverdensityClass     ), pointer                               :: criticalOverdensity_
    integer                                                                                :: iTime
    
    ! Validate arguments.
    if (regridCount < 2) call Galacticus_Error_Report('regridTimesConstructorInternal','regridCount > 2 is required')
    ! Store options.
    regridTimesConstructorInternal%dumpTrees=dumpTrees
    ! Construct array of grid expansion factors.
    call Alloc_Array(regridTimesConstructorInternal%timeGrid,[regridCount])
    cosmologyFunctions_  => cosmologyFunctions ()
    criticalOverdensity_ => criticalOverdensity()
    select case (snapshotSpacing)
    case (snapshotSpacingLinear                )
       regridTimesConstructorInternal%timeGrid=Make_Range(expansionFactorStart,expansionFactorEnd,regridCount,rangeTypeLinear     )
       ! Convert expansion factors to time.
       do iTime=1,regridCount
          regridTimesConstructorInternal%timeGrid(iTime)=cosmologyFunctions_%cosmicTime(regridTimesConstructorInternal%timeGrid(iTime))
       end do
    case (snapshotSpacingLogarithmic           )
       regridTimesConstructorInternal%timeGrid=Make_Range(expansionFactorStart,expansionFactorEnd,regridCount,rangeTypeLogarithmic)
       ! Convert expansion factors to time.
       do iTime=1,regridCount
          regridTimesConstructorInternal%timeGrid(iTime)=cosmologyFunctions_%cosmicTime(regridTimesConstructorInternal%timeGrid(iTime))
       end do
    case (snapshotSpacingLogCriticalOverdensity)
       ! Build a logarithmic grid in critical overdensity.
       regridTimesConstructorInternal%timeGrid                                                              &
            & =Make_Range(                                                                                  &
            &              criticalOverdensity_%value(cosmologyFunctions_%cosmicTime(expansionFactorStart)) &
            &             ,criticalOverdensity_%value(cosmologyFunctions_%cosmicTime(expansionFactorEnd  )) &
            &             ,regridCount                                                                      &
            &             ,rangeTypeLogarithmic                                                             &
            &            )
       ! Convert critical overdensity to time.
       do iTime=1,regridCount
          regridTimesConstructorInternal%timeGrid(iTime)=criticalOverdensity_%timeOfCollapse(regridTimesConstructorInternal%timeGrid(iTime))
       end do
    case (snapshotSpacingMillennium            )
       ! Use the timesteps used in the original Millennium Simulation as reported by Croton et al.  (2006; MNRAS; 365;
       ! 11; http://adsabs.harvard.edu/abs/2006MNRAS.365...11C). Note that we specifically use the redshifts for these
       ! snapshots as reported by the Millennium Database using query: "select z from Snapshots..MR".
       ! Check for consistent number of timesteps.
       if (regridCount /= 60) call Galacticus_Error_Report('regridTimesConstructorInternal','"millennium" grid spacing requires exactly 60 timesteps')
       ! Convert expansion factors to time.
       regridTimesConstructorInternal%timeGrid=[                                                                         &
            &                                   19.915688d0,18.243723d0,16.724525d0,15.343073d0,14.085914d0,12.940780d0, &
            &                                   11.896569d0,10.943864d0,10.073462d0, 9.277915d0, 8.549912d0, 7.883204d0, &
            &                                    7.272188d0, 6.711586d0, 6.196833d0, 5.723864d0, 5.288833d0, 4.888449d0, &
            &                                    4.519556d0, 4.179469d0, 3.865683d0, 3.575905d0, 3.308098d0, 3.060419d0, &
            &                                    2.831182d0, 2.618862d0, 2.422044d0, 2.239486d0, 2.070027d0, 1.912633d0, &
            &                                    1.766336d0, 1.630271d0, 1.503636d0, 1.385718d0, 1.275846d0, 1.173417d0, &
            &                                    1.077875d0, 0.988708d0, 0.905463d0, 0.827699d0, 0.755036d0, 0.687109d0, &
            &                                    0.623590d0, 0.564177d0, 0.508591d0, 0.456577d0, 0.407899d0, 0.362340d0, &
            &                                    0.319703d0, 0.279802d0, 0.242469d0, 0.207549d0, 0.174898d0, 0.144383d0, &
            &                                    0.115883d0, 0.089288d0, 0.064493d0, 0.041403d0, 0.019933d0, 0.000000d0  &
            &                                  ]
       do iTime=1,regridCount
          regridTimesConstructorInternal%timeGrid(iTime)=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(regridTimesConstructorInternal%timeGrid(iTime)))
       end do
    case (snapshotSpacingList                  )
       if (.not.present(snapshotTimes)) call Galacticus_Error_Report('regridTimesConstructorInternal','"list" grid spacing requires a list of snapshot times be supplied')
       regridTimesConstructorInternal%timeGrid=snapshotTimes
       call Sort_Do(regridTimesConstructorInternal%timeGrid)
    end select
    return
  end function regridTimesConstructorInternal

  elemental subroutine regridTimesDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorRegridTimes), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine regridTimesDestructor

  subroutine regridTimesOperate(self,tree)
    !% Perform a regrid times operation on a merger tree.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Nodes
    use               Galacticus_Error
    use               FGSL
    use               Numerical_Interpolation
    use               Kind_Numbers
    use               Merger_Trees_Dump
    implicit none
    class           (mergerTreeOperatorRegridTimes), intent(inout)                        :: self
    type            (mergerTree                   ), intent(inout), target                :: tree
    type            (treeNode                     )                             , pointer :: nodeChild               , nodeNext   , &
         &                                                                                   nodeSibling             , node
    type            (treeNodeList                 ), allocatable  , dimension(:)          :: newNodes
    integer         (kind=kind_int8               ), allocatable  , dimension(:)          :: highlightNodes
    class           (nodeComponentBasic           )                             , pointer :: basicChild              , basicParent, &
         &                                                                                   basic
    type            (mergerTree                   )                             , pointer :: currentTree
    type            (fgsl_interp_accel            )                                       :: interpolationAccelerator
    logical                                                                               :: interpolationReset
    integer         (c_size_t                     )                                       :: iNow                    , iParent    , &
         &                                                                                   iTime
    integer                                                                               :: allocErr
    double precision                                                                      :: massNow                 , massParent , &
         &                                                                                   timeNow                 , timeParent
    integer         (kind=kind_int8               )                                       :: firstNewNode            , nodeIndex

    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       ! Dump the unprocessed tree if required.
       if (self%dumpTrees) call Merger_Tree_Dump(                              &
            &                                    currentTree%index,            &
            &                                    currentTree%baseNode        , &
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
       ! Ensure interpolation accelerator gets reset.
       interpolationReset=.true.
       ! Find the current maximum node index in the tree.
       nodeIndex=0_kind_int8
       node => currentTree%baseNode
       do while (associated(node))
          nodeIndex=max(nodeIndex,node%index())
          call node%walkTree(node)
       end do
       firstNewNode=nodeIndex+1
       ! Walk the tree, locating branches which intersect grid times.
       node => currentTree%baseNode
       do while (associated(node))
          basic => node%basic()
          ! Skip this node if it is the root node.
          if (associated(node%parent)) then
             basicParent => node%parent%basic()             
             ! Get the time of this node and its parent.
             timeNow   =basic  %time()
             timeParent=basicParent%time()             
             ! Get masses of these halos.
             massNow   =basic  %mass()
             massParent=basicParent%mass()
             if (node%isPrimaryProgenitor()) then
                ! Remove the mass in any non-primary progenitors - we don't want to include their mass in the estimated mass
                ! growth rate of this node.
                nodeChild => node%parent%firstChild%sibling
                do while (associated(nodeChild))
                   basicChild => nodeChild%basic()
                   massParent          =  massParent-basicChild%mass()
                   nodeChild           => nodeChild%sibling
                end do
             else
                ! Halo is not the primary progenitor of its parent. Assume that its mass does not grow further.
                massParent=massNow
             end if
             ! Locate these times in the list of grid times.
             iNow   =Interpolate_Locate(self%timeGrid,interpolationAccelerator,timeNow   ,reset=interpolationReset)
             iParent=Interpolate_Locate(self%timeGrid,interpolationAccelerator,timeParent,reset=interpolationReset)
             ! For nodes existing precisely at a grid time, ignore this grid point. (These are, typically, nodes which have been created at these points.)
             if (timeParent == self%timeGrid(iParent)) iParent=iParent-1
             ! If the branch from node to parent spans one or more grid times, insert new nodes at those points.
             if (iParent > iNow) then
                ! Create new nodes.
                allocate(newNodes(iParent-iNow),stat=allocErr)
                if (allocErr/=0) call Galacticus_Error_Report('regridTimesOperate','unable to allocate new nodes')
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
                   ! Node is the main progenitor of its parent, so simply replace it with the final node in our list.
                   node%parent%firstChild  => newNodes(iParent-iNow)%node
                else
                   ! Node is not the main progenitor of its parent, so find the child node that has it as a sibling.
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
          ! Step to the next node.
          call node%walkTree(node)
       end do       
       ! Dump the intermediate tree if required.
       if (self%dumpTrees) then
          allocate(highlightNodes(nodeIndex-firstNewNode+2))
          highlightNodes(1)=currentTree%baseNode%index()
          do nodeIndex=1,nodeIndex-firstNewNode+1
             highlightNodes(nodeIndex+1)=firstNewNode+nodeIndex-1
          end do
          call Merger_Tree_Dump(                                    &
               &                currentTree%index,                  &
               &                currentTree%baseNode              , &
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
       node => currentTree%baseNode
       do while (associated(node))
          basic => node%basic()
          ! Record the next node to walk to.
          call node%walkTree(nodeNext)
          ! Get the time for this node.
          timeNow=basic%time()
          ! Find the closest time in the new time grid.
          iNow   =Interpolate_Locate(self%timeGrid,interpolationAccelerator,timeNow,reset=interpolationReset,closest=.true.)
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
                   ! Handle primary nodes with no children - simply make the parents main progenitor the sibling of the current node.
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
                   ! Find which sibling points the current node and link in the children of the current node.
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
          end if
          ! Step to the next node.
          node => nodeNext
       end do       
       ! Clean up interpolation objects.
       call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)
       ! Dump the processed tree if required.
       if (self%dumpTrees) call Merger_Tree_Dump(                               &
            &                                    currentTree%index,             &
            &                                    currentTree%baseNode         , &
            &                                    backgroundColor     ='white' , &
            &                                    nodeColor           ='black' , &
            &                                    highlightColor      ='black' , &
            &                                    edgeColor           ='black' , &
            &                                    nodeStyle           ='solid' , &
            &                                    highlightStyle      ='filled', &
            &                                    edgeStyle           ='solid' , &
            &                                    labelNodes          =.false. , &
            &                                    scaleNodesByLogMass =.true.  , &
            &                                    edgeLengthsToTimes  =.true.    &
            &                                   )
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine regridTimesOperate
