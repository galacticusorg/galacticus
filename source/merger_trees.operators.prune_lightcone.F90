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

  !!{
  Implements a prune-by-lightcone operator on merger trees.
  !!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneLightcone">
   <description>
     Provides a pruning-by-lightcone operator on merger trees, \emph{intended for use with the newer approach, in which galaxies
     are evolved to precisely the time of lightcone crossing}. For the older approach (in which galaxies were evolved to one of a
     fixed set of snapshots, and then output into the section of the lightcone corresponding to that snapshot) see the
     \refClass{mergerTreeOperatorPruneLightconeSnapshots} class.

     Trees which have no nodes which lie within the lightcone are completely pruned away. If the parameter \mono{[splitTrees]} is
     set to \mono{true} then any parts of a merger tree which does intersect the lightcone that exist after the latest time at which a
     constituent node of the tree intersects the lightcone will be pruned away also (possibly causing the tree to be split into
     multiple trees in a forest).
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorPruneLightconeSnapshots) :: mergerTreeOperatorPruneLightcone
     !!{
     A pruning-by-lightcone merger tree operator class.
     !!}
     private
   contains
     procedure :: processNode   => pruneLightconeProcessNode
     procedure :: timeIntersect => pruneLightconeTimeIntersect
  end type mergerTreeOperatorPruneLightcone

  interface mergerTreeOperatorPruneLightcone
     !!{
     Constructors for the pruning-by-lightcone merger tree operator class.
     !!}
     module procedure pruneLightconeConstructorParameters
     module procedure pruneLightconeConstructorInternal
  end interface mergerTreeOperatorPruneLightcone

contains

  function pruneLightconeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-by-lightcone merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeOperatorPruneLightcone)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (geometryLightconeClass          ), pointer       :: geometryLightcone_
    class (outputTimesClass                 ), pointer       :: outputTimes_
    logical                                                  :: splitTrees

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>splitTrees</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, prune away any nodes of the tree that are not needed to determine evolution up to the latest time at which a node is present inside the lightcone. This typically leads to a tree splitting into a forest of trees.</description>
    </inputParameter>
    <objectBuilder class="geometryLightcone" name="geometryLightcone_" source="parameters"/>
    <objectBuilder class="outputTimes"       name="outputTimes_"       source="parameters"/>
    !!]
    self=mergerTreeOperatorPruneLightcone(geometryLightcone_,outputTimes_,splitTrees)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="geometryLightcone_"/>
    <objectDestructor name="outputTimes_"      />
    !!]
    return
  end function pruneLightconeConstructorParameters

  function pruneLightconeConstructorInternal(geometryLightcone_,outputTimes_,splitTrees) result(self)
    !!{
    Internal constructor for the prune-by-lightcone merger tree operator class.
    !!}
    implicit none
    type   (mergerTreeOperatorPruneLightcone)                        :: self
    class  (geometryLightconeClass          ), intent(in   ), target :: geometryLightcone_
    class  (outputTimesClass                ), intent(in   ), target :: outputTimes_
    logical                                  , intent(in   )         :: splitTrees
    !![
    <constructorAssign variables="*geometryLightcone_, *outputTimes_, splitTrees"/>
    !!]

    self%bufferIsolatedHalos=.false.
    return
  end function pruneLightconeConstructorInternal

  logical function pruneLightconeProcessNode(self,node) result(processNode)
    !!{
    Return true if the given node should be processed for intersection with the lightcone.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeEvent, nodeEventSubhaloPromotion
    implicit none
    class           (mergerTreeOperatorPruneLightcone), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    class           (nodeComponentBasic              ), pointer       :: basic
    class           (nodeEvent                       ), pointer       :: event
    logical                                                           :: haveSubhaloPromotion
    double precision                                                  :: timeSubhaloPromotion

    processNode=.true.
    ! Do not process non-leaf nodes of the tree.
    if (associated(node%firstChild)) processNode=.false.
    ! Look for subhalo promotion events. Halos that are apparently at the tip of a branch, but have a subhalo promotion associated
    ! with their initial time are actually the target of the promotion. We do not want to process them here. Instead, the subhalo
    ! that promotes to them will be processed, and will trace through this subhalo promotion.
    haveSubhaloPromotion =  .false.
    event                => node%event
    do while (associated(event))
       ! Look for a handled event type. Subhalo promotions are handled.
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
       if (basic%time() == timeSubhaloPromotion) processNode=.false.
    end if
    return
  end function pruneLightconeProcessNode
  
  double precision function pruneLightconeTimeIntersect(self,node) result(timeIntersect)
    !!{
    Return the latest time at which the given node intersects the lightcone.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, mergerTree
    implicit none
    class           (mergerTreeOperatorPruneLightcone), intent(inout)               :: self
    type            (treeNode                        ), intent(inout)               :: node
    type            (mergerTree                      ), pointer                     :: tree
    class           (nodeComponentBasic              ), pointer                     :: basic
    double precision                                  , allocatable  , dimension(:) :: timesCrossing
    double precision                                                                :: timeEnd                , timeLightconeCrossing
    logical                                                                         :: nodeIntersectsLightcone

    ! Find the latest time in this forest.
    timeEnd =  0.0d0
    tree    => node%hostTree%firstTree
    do while (associated(tree))
       basic   => tree%nodeBase%basic()
       timeEnd =  max(timeEnd,basic%time())
       tree    => tree%nextTree
    end do
    ! Test for intersection with the lightcone.
    basic                   => node                   %basic                (                                       )
    timeLightconeCrossing   =  self%geometryLightcone_%timeLightconeCrossing(node,basic%time(),timeEnd,timesCrossing)
    nodeIntersectsLightcone =   allocated(timesCrossing)     &
         &                     .and.                         &
         &                      size     (timesCrossing) > 0
    if (nodeIntersectsLightcone) then
       timeIntersect=maxval(timesCrossing)
    else
       timeIntersect=0.0d0
    end if
    return
  end function pruneLightconeTimeIntersect
  
