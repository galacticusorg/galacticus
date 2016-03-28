!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% Contains a module which implements a merger tree operator which prunes branches to end at a fixed time.

  !# <mergerTreeOperator name="mergerTreeOperatorPruneByTime">
  !#  <description>Provides a merger tree operator which prunes branches to end at a fixed time.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneByTime
     !% A merger tree operator class which prunes branches to end at a fixed time.
     private
     double precision :: massMinimum , massMaximum, &
          &              timeEarliest
   contains
     final     ::            pruneByTimeDestructor
     procedure :: operate => pruneByTimeOperate
  end type mergerTreeOperatorPruneByTime
  
  interface mergerTreeOperatorPruneByTime
     !% Constructors for the prune-by-time merger tree operator class.
     module procedure pruneByTimeConstructorParameters
     module procedure pruneByTimeConstructorInternal
  end interface mergerTreeOperatorPruneByTime

contains

  function pruneByTimeConstructorParameters(parameters)
    !% Constructor for the prune-by-time merger tree operator class which takes a parameter set as input.
    use Cosmology_Functions
    implicit none
    type (mergerTreeOperatorPruneByTime)                :: pruneByTimeConstructorParameters
    type (inputParameters              ), intent(in   ) :: parameters
    class(cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    !# <inputParameterList label="allowedParameterNames" />
        
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>redshiftEarliest</name>
    !#   <source>parameters</source>
    !#   <variable>pruneByTimeConstructorParameters%timeEarliest</variable>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>Redshift at which to truncate merger tree branches.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMinimum</name>
    !#   <source>parameters</source>
    !#   <variable>pruneByTimeConstructorParameters%massMinimum</variable>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>Minimum mass for which to consider merger tree branches for truncation.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMaximum</name>
    !#   <source>parameters</source>
    !#   <variable>pruneByTimeConstructorParameters%massMaximum</variable>
    !#   <defaultValue>huge(0.0d0)</defaultValue>
    !#   <description>Maximum mass for which to consider merger tree branches for truncation.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    cosmologyFunctions_ => cosmologyFunctions()
    pruneByTimeConstructorParameters%timeEarliest             &
         & =cosmologyFunctions_ %cosmicTime(                  &
         &   cosmologyFunctions_%expansionFactorFromRedshift( &
         &    pruneByTimeConstructorParameters%timeEarliest   &
         &                                                  ) &
         &                                 )
    return
  end function pruneByTimeConstructorParameters

  function pruneByTimeConstructorInternal(timeEarliest,massMinimum,massMaximum)
    !% Internal constructor for the prune-by-time merger tree operator class.
    implicit none
    type            (mergerTreeOperatorPruneByTime)                :: pruneByTimeConstructorInternal
    double precision                               , intent(in   ) :: massMinimum                   , massMaximum, &
         &                                                            timeEarliest

    pruneByTimeConstructorInternal%timeEarliest=timeEarliest
    pruneByTimeConstructorInternal%massMinimum =massMinimum
    pruneByTimeConstructorInternal%massMaximum =massMaximum
    return
  end function pruneByTimeConstructorInternal

  elemental subroutine pruneByTimeDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorPruneByTime), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine pruneByTimeDestructor

  subroutine pruneByTimeOperate(self,tree)
    !% Perform a prune-by-time operation on a merger tree.
    use Merger_Trees_Pruning_Utilities
    implicit none
    class           (mergerTreeOperatorPruneByTime), intent(inout)          :: self
    type            (mergerTree                   ), intent(inout), target  :: tree
    type            (treeNode                     )               , pointer :: node              , nodeNext    , &
         &                                                                     nodeNew           , nodeChild
    class           (nodeComponentBasic           )               , pointer :: basic             , basicParent , &
         &                                                                     basicChild
    type            (mergerTree                   )               , pointer :: currentTree
    double precision                                                        :: massNow           , massParent  , &
         &                                                                     timeNow           , timeParent  , &
         &                                                                     massAtTimeEarliest
    
    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       ! Walk the tree, locating branches which cross the earliest allowed time.
       node => currentTree%baseNode
       do while (associated(node))
          ! Find the node to walk to next.
          nodeNext => node%walkTree()          
          ! Skip this node if it is the root node.
          if (associated(node%parent)) then
             ! Get basic components.
             basic       => node       %basic()
             basicParent => node%parent%basic()             
             ! Get the time of this node and its parent.
             timeNow   =basic  %time()
             timeParent=basicParent%time()             
             ! If the branch from node to parent spans the earliest time, insert a new node at that time.
             if (timeParent > self%timeEarliest .and. timeNow < self%timeEarliest) then
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
                ! Determine mass at truncation time.
                massAtTimeEarliest=massNow+(massParent-massNow)*(self%timeEarliest-timeNow)/(timeParent-timeNow)
                ! Check if mass is within range.
                if     (                                        &
                     &   massAtTimeEarliest >= self%massMinimum &
                     &  .and.                                   &
                     &   massAtTimeEarliest <= self%massMaximum &
                     & ) then
                   ! Update the node to walk to next
                   if (associated(node%sibling)) then
                      nodeNext => node%sibling
                   else
                      nodeNext => node%parent
                   end if
                   ! Create new node.
                   nodeNew => treeNode(hostTree=currentTree)
                   call nodeNew%indexSet(node%index())
                   ! Assign a time and a mass
                   basic => nodeNew%basic(autoCreate=.true.)
                   call basic%timeSet(self%timeEarliest )
                   call basic%massSet(massAtTimeEarliest)
                   ! No child node.
                   nodeNew%firstChild => null()
                   ! Link to parent node.
                   nodeNew%parent     => node%parent
                   ! Link  sibling to current node sibling.
                   nodeNew%sibling    => node%sibling
                   ! Link the parent if necessary.
                   if (node%isPrimaryProgenitor()) then
                      ! Node is the main progenitor of its parent, so simply replace it with the final node in our list.
                      node%parent%firstChild  => nodeNew
                   else
                      ! Node is not the main progenitor of its parent, so find the child node that has it as a sibling.
                      nodeChild => node%parent%firstChild
                      do while (.not.associated(nodeChild%sibling,node))
                         nodeChild => nodeChild%sibling
                      end do
                      nodeChild%sibling => nodeNew
                   end if
                   ! Clean the branch.
                   call Merger_Tree_Prune_Clean_Branch(node)
                   ! Destroy the branch.
                   call currentTree%destroyBranch(node)
                end if
             end if
          end if
          ! Step to the next node.
          node => nodeNext
       end do
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine pruneByTimeOperate
