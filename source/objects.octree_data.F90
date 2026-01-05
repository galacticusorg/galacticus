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
Contains a module which implements an object to store octree data structure.
!!}

module Octree_Data_Structure
  !!{
  Implements an object to store octree data structure.
  !!}
  use, intrinsic :: ISO_C_Binding , only : c_size_t
  implicit none
  private
  public :: octreeData

  !![
  <generic identifier="Type">
   <instance label="Double"      intrinsic="double precision              "/>
   <instance label="DoubleRank1" intrinsic="double precision, dimension(:)"/>
  </generic>
  !!]

  type :: octreeData
     !!{
     Type to give an octree data structure.
     !!}
     private
     type(octreeNode), pointer :: rootNode => null()
   contains
     !![
     <methods>
      <method description="Build an octree given the particle coordinates and weights." method="build"/>
      <method description="Destroy the octree object." method="destroy"/>
      <method description="Make a copy of the octree data structure." method="copy"/>
      <method description="Add a particle to the octree." method="addParticle"/>
      <method description="Remove a particle from the octree." method="removeParticle"/>
      <method description="Traverse the octree and compute the quantities needed." method="traverseCompute"/>
     </methods>
     !!]
     procedure :: build             => Build_Octree
     procedure :: destroy           => Destroy_Octree
     procedure :: copy              => Copy_Octree
     procedure :: addParticle       => Add_Particle_To_Octree
     procedure :: removeParticle    => Remove_Particle_From_Octree
     procedure ::                      Traverse_Octree_Compute{Type¦label}
     generic   :: traverseCompute   => Traverse_Octree_Compute{Type¦label}
  end type octreeData

  type :: octreeNodePointer
     !!{
     Pointer to an octreeNode object.
     !!}
     private
     type(octreeNode), pointer :: p => null()
  end type octreeNodePointer

  type :: octreeNode
     !!{
     Type to give a node in the octree.
     !!}
     private
     type            (octreeNodePointer), dimension(8) :: children
     type            (octreeNode       ), pointer      :: parent       => null()
     double precision                   , dimension(3) :: centerOfMass =  0.0d0  , boxBottomLeft
     double precision                                  :: nodeWeight   =  0.0d0  , boxWidth
     integer         (c_size_t         )               :: particleCount=  0      , depth        , &
          &                                               childIndex   =  0
  end type octreeNode

  abstract interface
     subroutine evaluate{Type¦label}(value,centerOfMass,nodeWeight,relativePosition,separation,separationSquare)
       {Type¦intrinsic}              , intent(inout) :: value
       double precision, dimension(3), intent(in   ) :: centerOfMass    , relativePosition
       double precision              , intent(in   ) :: nodeWeight      , separation      , &
            &                                           separationSquare
     end subroutine evaluate{Type¦label}
  end interface

contains

  subroutine Build_Octree(self,coordinates,weights)
    !!{
    Build an octree given the coordinates and weights of particles using the \cite{barnes_hierarchical_1986} algorithm.
    !!}
    implicit none
    class           (octreeData),                 intent(inout) :: self
    double precision            , dimension(:,:), intent(in   ) :: coordinates
    double precision            , dimension(:  ), intent(in   ) :: weights
    double precision            , dimension(3  )                :: boxBottomLeft           , coordinateMin, &
         &                                                         coordinateMax
    double precision                                            :: boxWidth                , extension
    double precision            , parameter                     :: extensionFraction=1.0d-1
    integer         (c_size_t  )                                :: particleCount           , i

    particleCount=size(weights)
    ! Determine the boundaries of the box that encloses all the particles.
    coordinateMin=minval(coordinates,dim=2)
    coordinateMax=maxval(coordinates,dim=2)
    boxWidth     =maxval(coordinateMax-coordinateMin)
    ! Extend the box boundaries.
    extension    =extensionFraction*boxWidth
    boxWidth     =boxWidth     +      extension
    boxBottomLeft=coordinateMin-0.5d0*extension
    ! Initialize the root node.
    call Create_Node(self%rootNode,boxBottomLeft,boxWidth,0_c_size_t)
    ! Loop over particles.
    do i=1, particleCount
       call self%addParticle(coordinates(:,i),weights(i))
    end do
    return
  end subroutine Build_Octree

  subroutine Destroy_Octree(self)
    !!{
    Destroy the octree.
    !!}
    implicit none
    class(octreeData), intent(inout) :: self

    ! Destroy the root node and all its child nodes.
    call Destroy_Node(self%rootNode)
    return
  end subroutine Destroy_Octree

  subroutine Copy_Octree(self,selfCopy)
    !!{
    Make a copy of the octree.
    !!}
    implicit none
    class(octreeData), intent(in   ) :: self
    class(octreeData), intent(inout) :: selfCopy

    ! Copy the whole octree structure.
    call Copy_Node(self%rootNode,selfCopy%rootNode)
    return
  end subroutine Copy_Octree

  subroutine Add_Particle_To_Octree(self,coordinate,weight)
    !!{
    Add a particle to the octree.
    !!}
    use :: Display            , only : displayMessage, verbosityLevelWarn
    use :: ISO_Varying_String , only : varying_string, assignment(=)     , operator(//)
    implicit none
    class           (octreeData    ),               intent(inout) :: self
    double precision                , dimension(3), intent(in   ) :: coordinate
    double precision                              , intent(in   ) :: weight
    type            (octreeNode    ),               pointer       :: workNode
    double precision                , dimension(3)                :: centerOfMass
    double precision                                              :: nodeWeight
    integer         (c_size_t      )                              :: childIndex
    type            (varying_string)                              :: message

    workNode => self%rootNode
    if     (                                                             &
         &   weight > 0.0d0                                              &
         &  .and.                                                        &
         &   all(coordinate >  workNode%boxBottomLeft                  ) &
         &  .and.                                                        &
         &   all(coordinate <= workNode%boxBottomLeft+workNode%boxWidth) &
         & ) then
       do while(associated(workNode))
          if (workNode%particleCount == 0) then
             ! If it reaches the leaf node, set the properties of the leaf node and walk to the next particle.
             workNode%particleCount =  1
             workNode%nodeWeight    =  weight
             workNode%centerOfMass  =  coordinate
             workNode               => null()
          else
             nodeWeight  =workNode%nodeWeight
             centerOfMass=workNode%centerOfMass
             ! If the current node has no child nodes, create its first child node based on the properties of the current node.
             if (workNode%particleCount == 1) then
                ! Find the index of the first child node.
                childIndex=getChildIndex(workNode,centerOfMass)
                call Create_Child_Node(workNode,childIndex)
             end if
             ! Add the new particle.
             workNode%particleCount =  workNode%particleCount+1
             workNode%nodeWeight    =  nodeWeight+weight
             workNode%centerOfMass  =  (nodeWeight*centerOfMass+weight*coordinate)/workNode%nodeWeight
             ! Find the index of child node that the new particle belongs to.
             childIndex=getChildIndex(workNode,coordinate)
             ! Create the child node if it does not exist.
             if (.not.associated(workNode%children(childIndex)%p)) then
                call Create_Child_Node(workNode,childIndex)
             end if
             ! Move to the child node.
             workNode               => workNode%children(childIndex)%p
          end if
       end do
    else
       if (weight <= 0.0d0) then
          message='Warning: Particle has nonpositive weight. It will not be added to the octree.'//char(10)
       else
          message='Warning: Particle is outside the root box of the current octree. It will not be added to the octree. To add the particle, the octree needs to be rebuilt with a larger root box.'//char(10)
       end if
       call displayMessage(message,verbosityLevelWarn)
    end if
    return
  end subroutine Add_Particle_To_Octree

  subroutine Remove_Particle_From_Octree(self,coordinate,weight)
    !!{
    Remove a particle from the octree.
    !!}
    use :: Display            , only : displayMessage, verbosityLevelWarn
    use :: ISO_Varying_String , only : varying_string, assignment(=)     , operator(//)
    implicit none
    class           (octreeData    ),               intent(inout) :: self
    double precision                , dimension(3), intent(in   ) :: coordinate
    double precision                              , intent(in   ) :: weight
    type            (octreeNode    ), pointer                     :: workNode    , leafNode, &
         &                                                           parentNode
    double precision                , dimension(3)                :: centerOfMass
    double precision                                              :: nodeWeight
    integer         (c_size_t      )                              :: childIndex  , i
    type            (varying_string)                              :: message

    ! Find the leaf node that the particle belongs to.
    leafNode => getLeafNode(self%rootNode,coordinate,weight)
    if (.not. associated(leafNode)) then
       message='Warning: Particle not found in the octree. No changes will be made.'//char(10)
       call displayMessage(message,verbosityLevelWarn)
    else
       workNode => leafNode
       do while(associated(workNode))
          if (workNode%particleCount == 1) then
             workNode%particleCount=0
             workNode%nodeWeight   =0.0d0
             workNode%centerOfMass =0.0d0
             parentNode => workNode%parent
             if (associated(parentNode)) then
                childIndex=workNode%childIndex
                call Destroy_Node(parentNode%children(childIndex)%p)
             end if
             ! Move to the parent node.
             workNode   => parentNode
          else
             nodeWeight  =workNode%nodeWeight
             centerOfMass=workNode%centerOfMass
             workNode%particleCount=workNode%particleCount-1
             workNode%nodeWeight   =nodeWeight-weight
             workNode%centerOfMass =(nodeWeight*centerOfMass-weight*coordinate)/workNode%nodeWeight
             ! If the current node has only one particle left, destroy all its existed child nodes.
             if (workNode%particleCount == 1) then
                do i=1,8
                   if (associated(workNode%children(i)%p)) call Destroy_Node(workNode%children(i)%p)
                end do
             end if
             ! Move to the parent node.
             workNode   => workNode%parent
          end if
       end do
    end if
    return
  end subroutine Remove_Particle_From_Octree

  subroutine Traverse_Octree_Compute{Type¦label}(self,coordinate,weight,thetaTolerance,value,evaluateFunction,isExternalParticle)
    !!{
    Octree traversal using the \cite{barnes_hierarchical_1986} algorithm.
    !!}
    implicit none
    class           (octreeData          ),               intent(in   ) :: self
    double precision                      , dimension(3), intent(in   ) :: coordinate
    double precision                                    , intent(in   ) :: weight                  , thetaTolerance
    {Type¦intrinsic}                                    , intent(inout) :: value
    procedure       (evaluate{Type¦label})                              :: evaluateFunction
    logical                               , optional    , intent(in   ) :: isExternalParticle
    type            (octreeNode          ), pointer                     :: workNode                , siblingNode
    double precision                      , dimension(3)                :: centerOfMass            , relativePosition
    double precision                                                    :: nodeWeight              , separation                        , &
         &                                                                 separationSquare        , theta
    double precision                      , parameter                   :: thetaLarge=10.0d0
    integer         (c_size_t            )                              :: childIndex
    logical                                                             :: isExternalParticleActual

    ! Check whether the particle whose property (acceleration, potential,...) is to be computed is an external particle,
    !  i.e. not contained in the current octree.
    isExternalParticleActual=.false.
    if (present(isExternalParticle)) isExternalParticleActual=isExternalParticle
    ! Start the octree traversal.
    workNode => self%rootNode
    value = 0.0d0
    do while(associated(workNode))
       ! Skip empty node.
       if (workNode%nodeWeight <= 0.0d0) then
          siblingNode => getNonEmptySibling(workNode)
          do while(.not.associated(siblingNode) .and. associated(workNode%parent))
             workNode    => workNode%parent
             siblingNode => getNonEmptySibling(workNode)
          end do
          workNode    => siblingNode
          cycle
       end if
       nodeWeight  =workNode%nodeWeight
       centerOfMass=workNode%centerOfMass
       ! Separation between the particle and the current node.
       relativePosition=coordinate-workNode%centerOfMass
       separationSquare=sum (relativePosition**2)
       separation      =sqrt(separationSquare   )
       ! Compute the opening angle.
       if (separation > 0.0d0) then
          theta=workNode%boxWidth/separation
       else
          theta=thetaLarge
       end if
       if (theta < thetaTolerance .or. workNode%particleCount <= 1) then
          if (.not.isExternalParticleActual) then
             ! Check whether the current particle belongs to the node.
             if (all(coordinate > workNode%boxBottomLeft) .and. all(coordinate <= workNode%boxBottomLeft+workNode%boxWidth)) then
                nodeWeight=workNode%nodeWeight-weight
                if (nodeWeight > 0.0d0) then
                   centerOfMass    =(workNode%nodeWeight*workNode%centerOfMass-weight*coordinate)/nodeWeight
                   relativePosition=coordinate-centerOfMass
                   separationSquare=sum (relativePosition**2)
                   separation      =sqrt(separationSquare)
                else
                   centerOfMass    =0.0d0
                   relativePosition=0.0d0
                   separationSquare=0.0d0
                   separation      =0.0d0
                end if
             end if
          end if
          ! Do necessary computations.
          call evaluateFunction(value,centerOfMass,nodeWeight,relativePosition,separation,separationSquare)
          siblingNode => getNonEmptySibling(workNode)
          do while(.not.associated(siblingNode) .and. associated(workNode%parent))
             workNode    => workNode%parent
             siblingNode => getNonEmptySibling(workNode)
          end do
          workNode    => siblingNode
       else
          ! Find the first non-empty child node.
          childIndex =1
          do while((.not.associated(workNode%children(childIndex)%p)) .and. childIndex < 8)
             childIndex =  childIndex+1
          end do
          workNode   => workNode%children(childIndex)%p
       end if
    end do
    return
  end subroutine Traverse_Octree_Compute{Type¦label}

  subroutine Create_Node(node,boxBottomLeft,boxWidth,depth)
    !!{
    Create a node in the octree.
    !!}
    implicit none
    type            (octreeNode), pointer                       :: node
    double precision            , dimension(3), intent(in   )   :: boxBottomLeft
    double precision                          , intent(in   )   :: boxWidth
    integer         (c_size_t  )              , intent(in   )   :: depth
    integer         (c_size_t  )                                :: i

    if (associated(node)) call Destroy_Node(node)
    allocate(node)
    node%boxBottomLeft=boxBottomLeft
    node%boxWidth     =boxWidth
    node%depth        =depth
    do i=1,8
       node%children(i)%p => null()
    end do
    return
  end subroutine Create_Node

  recursive subroutine Copy_Node(node,nodeCopy)
    !!{
    Make a copy of a node and all its child nodes.
    !!}
    implicit none
    type   (octreeNode), pointer, intent(in   ) :: node
    type   (octreeNode), pointer                :: nodeCopy
    integer(c_size_t  )                         :: i

    nodeCopy => null()
    if (associated(node)) then
       allocate(nodeCopy)
       ! Copy the node properties.
       nodeCopy%boxBottomLeft=  node%boxBottomLeft
       nodeCopy%boxWidth     =  node%boxWidth
       nodeCopy%depth        =  node%depth
       nodeCopy%centerOfMass =  node%centerOfMass
       nodeCopy%nodeWeight   =  node%nodeWeight
       nodeCopy%particleCount=  node%particleCount
       nodeCopy%childIndex   =  node%childIndex
       nodeCopy%parent       => null()
       ! Copy the child nodes.
       do i=1,8
          call Copy_Node(node%children(i)%p,nodeCopy%children(i)%p)
          if (associated(nodeCopy%children(i)%p)) then
             nodeCopy%children(i)%p%parent => nodeCopy
          end if
       end do
    end if
    return
  end subroutine Copy_Node

  recursive subroutine Destroy_Node(node)
    !!{
    Destroy a node and all its child nodes.
    !!}
    implicit none
    type   (octreeNode), pointer :: node
    integer(c_size_t  )          :: i

    do i=1,8
       if (associated(node%children(i)%p)) then
          call Destroy_Node(node%children(i)%p)
          node%children(i)%p => null()
       end if
    end do
    node%parent => null()
    deallocate(node)
    node        => null()
    return
  end subroutine Destroy_Node

  subroutine Create_Child_Node(node,childIndex)
    !!{
    Create a child node with the given child index.
    !!}
    implicit none
    type            (octreeNode), pointer                       :: node
    integer         (c_size_t  )                , intent(in   ) :: childIndex
    double precision            , dimension(3  )                :: childBoxBottomLeft
    double precision                                            :: childBoxWidth
    double precision            , dimension(3,8), parameter     :: convertTable=reshape([0.0d0,0.0d0,0.0d0,  &
         &                                                                               0.5d0,0.0d0,0.0d0,  &
         &                                                                               0.0d0,0.5d0,0.0d0,  &
         &                                                                               0.5d0,0.5d0,0.0d0,  &
         &                                                                               0.0d0,0.0d0,0.5d0,  &
         &                                                                               0.5d0,0.0d0,0.5d0,  &
         &                                                                               0.0d0,0.5d0,0.5d0,  &
         &                                                                               0.5d0,0.5d0,0.5d0], &
         &                                                                               shape=[3,8])

    ! Find the position and width of the child box given the child index. The indices of the 8 child nodes
    ! are determined following the Z-order definition.
    childBoxBottomLeft=node%boxBottomLeft+convertTable(:,childIndex)*node%boxWidth
    childBoxWidth     =0.5d0*node%boxWidth
    ! Create the child node.
    call Create_Node(node%children(childIndex)%p,childBoxBottomLeft,childBoxWidth,node%depth+1)
    node%children(childIndex)%p%parent     => node
    node%children(childIndex)%p%childIndex =  childIndex
    ! If it is the first child node, copy the property of the parent node to the new child node.
    if (node%particleCount == 1) then
       node%children(childIndex)%p%particleCount=1
       node%children(childIndex)%p%nodeWeight   =node%nodeWeight
       node%children(childIndex)%p%centerOfMass =node%centerOfMass
    end if
    return
  end subroutine Create_Child_Node

  function getChildIndex(node,coordinate)
    !!{
    Return the index of the child node that a particle belongs to.
    !!}
    use :: Display            , only : displayMessage, verbosityLevelWarn
    use :: ISO_Varying_String , only : varying_string, assignment(=)     , operator(//)
    implicit none
    integer         (c_size_t      )                                  :: getChildIndex
    type            (octreeNode    ), pointer         , intent(in   ) :: node
    double precision                , dimension(3    ), intent(in   ) :: coordinate
    integer         (c_size_t      ), dimension(3    )                :: coordinateIndex
    integer         (c_size_t      ), dimension(2,2,2), parameter     :: convertTable=reshape([1,2,3,4,5,6,7,8], &
         &                                                                                     shape=[2,2,2])
    character       (len=80        )                                  :: valueString
    type            (varying_string)                                  :: message

    ! Find the child index given the coordinate. The indices of the 8 child nodes
    ! are determined following the Z-order definition.
    coordinateIndex=ceiling(2.0d0*(coordinate-node%boxBottomLeft)/node%boxWidth)
    if (minval(coordinateIndex) < 1 .or. maxval(coordinateIndex) > 2) then
       ! Report a warning if the particle is considered to be outside the current node box due to roundoff errors.
       message='Warning: Particel is considered to be outside the current node box. This may be caused by runoff errors.'//char(10)
       write (valueString,'(e26.17e3,e26.17e3,e26.17e3)') coordinate
       message=message//'  Particle coordinate: '//trim(valueString)//char(10)
       write (valueString,'(e26.17e3,e26.17e3,e26.17e3)') node%boxBottomLeft
       message=message//'  Box bottom left    : '//trim(valueString)//char(10)
       write (valueString,'(e26.17e3,e26.17e3,e26.17e3)') node%boxBottomLeft+node%boxWidth
       message=message//'  Box top right      : '//trim(valueString)//char(10)
       call displayMessage(message,verbosityLevelWarn)
       ! Assign the particle to the nearest node.
       where(coordinateIndex < 1)
          coordinateIndex = 1
       else where (coordinateIndex > 2)
          coordinateIndex = 2
       end where
    end if
    getChildIndex=convertTable(coordinateIndex(1),coordinateIndex(2),coordinateIndex(3))
    return
  end function getChildIndex

  function getNonEmptySibling(node) result(sibling)
    !!{
    Return the next non-empty sibling node.
    !!}
    implicit none
    type   (octreeNode), pointer                :: sibling
    type   (octreeNode), pointer, intent(in   ) :: node
    type   (octreeNode), pointer                :: workNode
    integer(c_size_t  )                         :: childIndex, siblingIndex

    childIndex =  node%childIndex
    workNode   => node%parent
    if (childIndex > 0 .and. childIndex < 8) then
       siblingIndex=childIndex+1
       do while((.not. associated(workNode%children(siblingIndex)%p)) .and. siblingIndex < 8)
          siblingIndex=siblingIndex+1
       end do
       sibling => workNode%children(siblingIndex)%p
    else
       ! If no non-empty sibling node is found, return a null pointer.
       sibling => null()
    end if
    return
  end function getNonEmptySibling

  function getLeafNode(rootNode,coordinate,weight) result(leafNode)
    !!{
    Return the leaf node that a particle belongs to.
    !!}
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    type            (octreeNode    ), pointer                     :: leafNode
    type            (octreeNode    ), pointer     , intent(in   ) :: rootNode
    double precision                , dimension(3), intent(in   ) :: coordinate
    double precision                              , intent(in   ) :: weight
    type            (octreeNode    ), pointer                     :: workNode
    integer         (c_size_t      )                              :: childIndex

    workNode => rootNode
    leafNode => null()
    if     (                                                             &
         &   weight > 0.0d0                                              &
         &  .and.                                                        &
         &   all(coordinate >  workNode%boxBottomLeft                  ) &
         &  .and.                                                        &
         &   all(coordinate <= workNode%boxBottomLeft+workNode%boxWidth) &
         & ) then
       ! If the particle has nonzero weight and is inside the root box, search for
       ! the leaf node it belongs to.
       do while(associated(workNode))
          if (workNode%particleCount == 1) then
             leafNode => workNode
             workNode => null()
          else
             childIndex =  getChildIndex(workNode,coordinate)
             ! Move to the child node.
             workNode   => workNode%children(childIndex)%p
          end if
       end do
       ! Check whether the particle coordinate and weight matches with the leaf node properties.
       if (associated(leafNode)) then
          if (                                                                    &
               any(Values_Differ(leafNode%centerOfMass,coordinate,relTol=1.0d-6)) &
              .or.                                                                &
                   Values_Differ(leafNode%nodeWeight  ,weight    ,relTol=1.0d-6)  &
             ) then
             leafNode => null()
          end if
       end if
    end if
    return
  end function getLeafNode

end module Octree_Data_Structure
