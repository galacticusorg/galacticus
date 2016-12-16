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

!% Contains a module which implements a prune-by-lightcone operator on merger trees.

  use Geometry_Lightcones
  use Satellite_Oprhan_Distributions
  
  !# <mergerTreeOperator name="mergerTreeOperatorPruneLightcone" defaultThreadPrivate="yes">
  !#  <description>Provides a pruning-by-lightcone operator on merger trees. Trees which have no nodes which lie within the lightcone are completely pruned away.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneLightcone
     !% A pruning-by-mass merger tree operator class.
     private
     class  (geometryLightconeClass          ), pointer :: geometryLightcone_
     class  (satelliteOrphanDistributionClass), pointer :: satelliteOrphanDistribution_
     logical                                            :: bufferIsolatedHalos         , positionHistoryAvailable
   contains     
     final     ::            pruneLightconeDestructor
     procedure :: operate => pruneLightconeOperate
  end type mergerTreeOperatorPruneLightcone

  interface mergerTreeOperatorPruneLightcone
     !% Constructors for the pruning-by-lightcone merger tree operator class.
     module procedure pruneLightconeConstructorParameters
     module procedure pruneLightconeConstructorInternal
  end interface mergerTreeOperatorPruneLightcone

contains

  function pruneLightconeConstructorParameters(parameters) result(self)
    !% Constructor for the prune-by-lightcone merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorPruneLightcone)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>bufferIsolatedHalos</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>self%bufferIsolatedHalos</variable>
    !#   <description>If true, intersection of a tree with the lightcone will be determined using the positions non-isolated (a.k.a. ``satellite'') halos, and of isolated halos (a.k.a ``centrals'') with a buffer region (with radius equal to the extent of the orphan satellite distribution---see \S\ref{sec:methodsSatelliteOrphanDistribution}) placed around each such halo, and any intersection of that region with the lightcone is sufficient to prevent pruning of the tree. If this parameter is {\normalfont \ttfamily false} then (unbuffered) positions of all halos are used for determining intersection with the lightcone---this requires complete (i.e. throughout the extent of their existance) knowledge of non-isolated halos prior to application of this operator.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="geometryLightcone"           name="self%geometryLightcone_"           source="parameters"/>
    !# <objectBuilder class="satelliteOrphanDistribution" name="self%satelliteOrphanDistribution_" source="parameters"/>
    call pruneLightconeValidate(self)
    return
  end function pruneLightconeConstructorParameters

  function pruneLightconeConstructorInternal(geometryLightcone_,satelliteOrphanDistribution_,bufferIsolatedHalos) result(self)
    !% Internal constructor for the prune-by-lightcone merger tree operator class.
    implicit none
    type   (mergerTreeOperatorPruneLightcone)                        :: self
    class  (geometryLightconeClass          ), intent(in   ), target :: geometryLightcone_
    class  (satelliteOrphanDistributionClass), intent(in   ), target :: satelliteOrphanDistribution_
    logical                                  , intent(in   )         :: bufferIsolatedHalos
    !# <constructorAssign variables="*geometryLightcone_, *satelliteOrphanDistribution_, bufferIsolatedHalos"/>

    call pruneLightconeValidate(self)
    return
  end function pruneLightconeConstructorInternal

  subroutine pruneLightconeDestructor(self)
    !% Destructor for the lightcone merger tree operator function class.
    implicit none
    type(mergerTreeOperatorPruneLightcone), intent(inout) :: self
    
    !# <objectDestructor name="self%geometryLightcone_"          />
    !# <objectDestructor name="self%satelliteOrphanDistribution_"/>
    return
  end subroutine pruneLightconeDestructor

  subroutine pruneLightconeValidate(self)
    !% Validate the lightcone pruning operator.
    use Galacticus_Error
    use Array_Utilities
    implicit none
    class(mergerTreeOperatorPruneLightcone), intent(inout) :: self

    ! If buffering is applied to isolated halos, then satellite halos are to be checked for intersection only for as long as their
    ! position is known. In this case, check if satellite position history can be obtained, and also require that satellite time
    ! of merging can be both read and written.
    if (self%bufferIsolatedHalos) then
       self%positionHistoryAvailable=defaultPositionComponent%positionHistoryIsGettable()
       if     (                                                                                                                    &
            &  .not.(                                                                                                              &
            &         defaultSatelliteComponent%timeOfMergingIsGettable()                                                          &
            &        .and.                                                                                                         &
            &         defaultSatelliteComponent%timeOfMergingIsSettable()                                                          &
            &       )                                                                                                              &
            & )                                                                                                                    &
            & call Galacticus_Error_Report                                                                                         &
            &      (                                                                                                               &
            &       'pruneLightconeValidate'                                                                                     , &
            &       'buffering isolated halos requires that the position history property of the position component be gettable'// &
            &       Galacticus_Component_List(                                                                                     &
            &                                 'satellite'                                                                       ,  &
            &                                   defaultSatelliteComponent%timeOfMergingAttributeMatch(requireGettable=.true.)      &
            &                                  .intersection.                                                                      &
            &                                   defaultSatelliteComponent%timeOfMergingAttributeMatch(requireSettable=.true.)      &
            &                                 )                                                                                    &
            &      )
    end if
    return
  end subroutine pruneLightconeValidate
  
  subroutine pruneLightconeOperate(self,tree)
    !% Perform a prune-by-lightcone operation on a merger tree.
    use Merger_Trees_Pruning_Utilities
    use Histories
    implicit none
    class           (mergerTreeOperatorPruneLightcone), intent(inout)         :: self
    type            (mergerTree                      ), intent(inout), target :: tree
    type            (treeNode                        ), pointer               :: node               , nodeNext            , &
         &                                                                       nodeParent
    type            (mergerTree                      ), pointer               :: treeCurrent        , treeNext
    class           (nodeComponentBasic              ), pointer               :: basic
    class           (nodeComponentPosition           ), pointer               :: position
    class           (nodeComponentSatellite          ), pointer               :: satellite
    type            (history                         )                        :: positionHistory
    double precision                                                          :: radiusBuffer       , timeOfMergingCurrent
    logical                                                                   :: intersectsLightcone

    ! Set buffer size to zero by default.
    radiusBuffer=0.0d0
    ! Iterate over trees.
    treeCurrent => tree
    nodeParent  => null()
    do while (associated(treeCurrent))
       ! Get root node of the tree.
       node  => treeCurrent%baseNode
       do while (associated(node))
          ! If buffering isolated halos, set a suitable buffer radius. For isolated halos this is the maximum extent of their
          ! orphan satellite population. For non-isolated halos, no buffer is required. In addition, for non-isolated halos, limit
          ! the time over which they will be checked for intersection with the lightcone to the maximum time for which they have a
          ! defined position.
          if (self%bufferIsolatedHalos) then
             ! Store current parent of the node.
             nodeParent => node%parent
             ! Handle satellites and centrals.
             if (node%isSatellite()) then
                ! No buffer is required for a satellite node.
                radiusBuffer=0.0d0
                ! Find the extent of the position history known for this node. Limit the merging time to the final time for which
                ! position is known (or the current time if no position history is available). Record the current time of merging
                ! so that it can be reset after testing for intersection.
                position             => node    %position       ()
                basic                => node    %basic          ()
                satellite            => node    %satellite      ()
                positionHistory      =  position%positionHistory()
                timeOfMergingCurrent =  satellite%timeOfMerging ()
                if (positionHistory%exists()) then
                   call satellite%timeOfMergingSet(positionHistory%time(size(positionHistory%time)))
                else
                   call satellite%timeOfMergingSet(basic          %time(                          ))
                end if
             else
                ! For a central - set a suitable buffer.
                radiusBuffer=self%satelliteOrphanDistribution_%extent(node)
                ! Also temporarily decouple from the tree.
                node%parent => null()
             end if
          end if
          ! Test for intersection with the lightcone.
          intersectsLightcone=self%geometryLightcone_%isInLightcone(node,atPresentEpoch=.false.,radiusBuffer=radiusBuffer)
          ! Reset the time of merging if it was adjusted above.
          if (self%bufferIsolatedHalos) then
             ! Recouple to the tree.
             node%parent => nodeParent
             ! Reset time of merging.
             if (node%isSatellite()) call satellite%timeOfMergingSet(timeOfMergingCurrent)
          end if
          ! If intersection with lightcone was detected then this tree can not be pruned - return immediately.
          if (intersectsLightcone) return
          ! Walk to the next node in the tree.
          node => node%walkTree()
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Entire forest is outside lightcone. Destroy all but the base node in the first tree. (Leaving just the base node makes the
    ! tree inert - i.e. it can not do anything.) Destroy any additional trees in the forest.
    node => tree%baseNode%firstChild
    do while (associated(node))
       nodeNext => node%sibling
       call Merger_Tree_Prune_Clean_Branch(node)
       call node%destroyBranch()
       deallocate(node)
       node => nodeNext
    end do
    treeCurrent          => tree%nextTree
    tree       %nextTree => null()
    do while (associated(treeCurrent))
       treeNext => treeCurrent%nextTree
       ! Destroy the tree.
       call treeCurrent%destroy()
       deallocate(treeCurrent)
       ! Move to the next tree.
       treeCurrent => treeNext
    end do
    return
  end subroutine pruneLightconeOperate
