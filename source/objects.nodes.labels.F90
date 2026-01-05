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
Contains a module which provides functionality for labeling nodes.
!!}

module Nodes_Labels
  !!{
  Provides functionality for labeling nodes.
  !!}
  use :: ISO_Varying_String, only : varying_string
  use :: Kind_Numbers      , only : kind_int8
  private
  public :: nodeLabelRegister    , nodeLabelSet  , nodeLabelNames, nodeLabelList     , &
       &    nodeLabelDescriptions, nodeLabelUnset, nodeLabelCount, nodeLabelIsPresent

  integer                                            :: nodeLabelsID=-1
  type   (varying_string), allocatable, dimension(:) :: labels_        , descriptions_

contains

  integer function nodeLabelRegister(label,description) result(labelID)
    !!{
    Register an allowed label for nodes, and return an ID for that label.
    !!}
    use :: Error             , only : Error_Report
    use :: Events_Hooks      , only : nodePromotionEventGlobal, openMPThreadBindingNone
    use :: ISO_Varying_String, only : operator(==)            , assignment(=)
    implicit none
    character(len=*         ), intent(in   )               :: label
    character(len=*         ), intent(in   ), optional     :: description
    type     (varying_string), allocatable  , dimension(:) :: labelsTmp  , descriptionsTmp
    integer                                                :: i

    labelID=-1
    !$omp critical(lockNodeLabels)
    if (allocated(labels_)) then
       do i=1,size(labels_)
          if (labels_(i) == label) then
             labelID=i
             if (descriptions_(i) == '' .and. present(description)) &
                  & descriptions_(i)=description
             exit
          end if
       end do
    end if
    if (labelID == -1) then
       if (allocated(labels_)) then
          call move_alloc(labels_      ,labelsTmp      )
          call move_alloc(descriptions_,descriptionsTmp)
          allocate(labels_      (size(labelsTmp      )+1))
          allocate(descriptions_(size(descriptionsTmp)+1))
          labels_      (1:size(labelsTmp      ))=labelsTmp
          descriptions_(1:size(descriptionsTmp))=descriptionsTmp
       else
          allocate(labels_      (                      1))
          allocate(descriptions_(                      1))
       end if
       labelID                  =size(labels_)
       labels_         (labelID)=label
       if (present(description)) then
          descriptions_(labelID)=description
       else
          descriptions_(labelID)=''
       end if
    end if
    if (size(labels_) > bit_size(0_kind_int8)) call Error_Report('insufficient storage space for labels'//{introspection:location})
    if (nodeLabelsID < 0) then
       !![
       <addMetaProperty component="basic" name="nodeLabels" type="longInteger" isCreator="yes" id="nodeLabelsID"/>
       !!]
       call nodePromotionEventGlobal%attach(nodeLabelsID,mergeLabels,openMPThreadBindingNone,label='mergeLabels')
    end if
    !$omp end critical(lockNodeLabels)
    return
  end function nodeLabelRegister

  subroutine mergeLabels(ID,node)
    !!{
    Merge labels of nodes on promotion.
    !!}
    use :: Galacticus_Nodes, only : treeNode    , nodeComponentBasic
    implicit none
    class  (*                 ), intent(inout)          :: ID
    type   (treeNode          ), intent(inout), target  :: node
    integer(kind_int8         )                         :: label        , labelParent, &
         &                                                 labelCombined
    class  (nodeComponentBasic)               , pointer :: basic        , basicParent
    !$GLC attributes unused :: ID

    basic         => node              %basic                          (            )
    basicParent   => node       %parent%basic                          (            )
    label         =  basic             %longIntegerRank0MetaPropertyGet(nodeLabelsID)
    labelParent   =  basicParent       %longIntegerRank0MetaPropertyGet(nodeLabelsID)
    labelCombined =  ior(label,labelParent)
    call basic%longIntegerRank0MetaPropertySet(nodeLabelsID,labelCombined)
    return
  end subroutine mergeLabels

  integer function nodeLabelCount()
    !!{
    Return a count of the number of node labels.
    !!}
    implicit none
    
    !$omp critical(lockNodeLabels)
    if (allocated(labels_)) then
       nodeLabelCount=size(labels_)
    else
       nodeLabelCount=0
    end if
    !$omp end critical(lockNodeLabels)
    return
  end function nodeLabelCount
  
  subroutine nodeLabelUnset(labelID,node)
    !!{
    Unset a label on a node.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : treeNode    , nodeComponentBasic
    implicit none
    integer                    , intent(in   ) :: labelID
    type   (treeNode          ), intent(inout) :: node
    class  (nodeComponentBasic), pointer       :: basic

    !$omp critical(lockNodeLabels)
    if (.not.allocated(labels_) .or. labelID > size(labels_)) call Error_Report('label ID is out of range'//{introspection:location})
    !$omp end critical(lockNodeLabels)
    basic => node%basic(autoCreate=.true.)
    call basic%longIntegerRank0MetaPropertySet(                                                                     &
         &                                                                                 nodeLabelsID           , &
         &                                     ibclr(                                                               &
         &                                           basic%longIntegerRank0MetaPropertyGet(nodeLabelsID),labelID-1  &
         &                                          )                                                               &
         &                                    )
    return
  end subroutine nodeLabelUnset
  
  subroutine nodeLabelSet(labelID,node)
    !!{
    Set a label on a node.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : treeNode    , nodeComponentBasic
    implicit none
    integer                    , intent(in   ) :: labelID
    type   (treeNode          ), intent(inout) :: node
    class  (nodeComponentBasic), pointer       :: basic

    !$omp critical(lockNodeLabels)
    if (.not.allocated(labels_) .or. labelID > size(labels_)) call Error_Report('label ID is out of range'//{introspection:location})
    !$omp end critical(lockNodeLabels)
    basic => node%basic(autoCreate=.true.)
    call basic%longIntegerRank0MetaPropertySet(                                                                     &
         &                                                                                 nodeLabelsID           , &
         &                                     ibset(                                                               &
         &                                           basic%longIntegerRank0MetaPropertyGet(nodeLabelsID),labelID-1  &
         &                                          )                                                               &
         &                                    )
    return
  end subroutine nodeLabelSet
  
  logical function nodeLabelIsPresent(labelID,node)
    !!{
    Return true if the specified label is present in the node.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : treeNode    , nodeComponentBasic
    implicit none
    integer                    , intent(in   ) :: labelID
    type   (treeNode          ), intent(inout) :: node
    class  (nodeComponentBasic), pointer       :: basic
    integer                    , parameter     :: bitTrue=ibset(0,0)
    integer(kind_int8         )                :: label_ 

    !$omp critical(lockNodeLabels)
    if (.not.allocated(labels_) .or. labelID > size(labels_)) call Error_Report('label ID is out of range'//{introspection:location})
    !$omp end critical(lockNodeLabels)
    basic              => node %basic                          (autoCreate=.true.      )
    label_             =  basic%longIntegerRank0MetaPropertyGet(           nodeLabelsID)
    nodeLabelIsPresent =  ibits(label_,labelID-1,1) == bitTrue
    return
  end function nodeLabelIsPresent

  subroutine nodeLabelNames(labels)
    !!{
    Return a list of label names.
    !!}
    implicit none
    type(varying_string), intent(inout), allocatable, dimension(:) :: labels

    !$omp critical(lockNodeLabels)
    if (allocated(labels_)) then
       if (allocated(labels)) then
          if (size(labels) /= size(labels_)) then
             deallocate(labels               )
             allocate  (labels(size(labels_)))
          end if
       else
          allocate     (labels(size(labels_)))
       end if
       labels=labels_
    else
       if (allocated(labels)) deallocate(labels)
    end if
    !$omp end critical(lockNodeLabels)
    return
  end subroutine nodeLabelNames

  subroutine nodeLabelDescriptions(descriptions)
    !!{
    Return a list of label descriptions.
    !!}
    implicit none
    type(varying_string), intent(inout), allocatable, dimension(:) :: descriptions

    !$omp critical(lockNodeLabels)
    if (allocated(descriptions_)) then
       if (allocated(descriptions)) then
          if (size(descriptions) /= size(descriptions_)) then
             deallocate(descriptions                    )
             allocate  (descriptions(size(descriptions)))
          end if
       else
          allocate     (descriptions(size(descriptions)))
       end if
       descriptions=descriptions_
    else
       if (allocated(descriptions)) deallocate(descriptions)
    end if
    !$omp end critical(lockNodeLabels)
    return
  end subroutine nodeLabelDescriptions

  subroutine nodeLabelList(node,labels)
    !!{
    Return a list of labels for the provided node. Labels are returned as an array, with value 0 or 1 to indicate that the label
    is not or is set.
    !!}
    use :: Galacticus_Nodes, only : treeNode, nodeComponentBasic
    implicit none
    type   (treeNode          ), intent(inout)                            :: node
    integer(kind_int8         ), intent(inout), allocatable, dimension(:) :: labels
    class  (nodeComponentBasic), pointer                                  :: basic
    integer(kind_int8         )                                           :: label_            , i
    integer                    , parameter                                :: bitTrue=ibset(0,0)

    !$omp critical(lockNodeLabels)
    if (allocated(labels_)) then
       if (allocated(labels)) then
          if (size(labels) /= size(labels_)) then
             deallocate(labels               )
             allocate  (labels(size(labels_)))
          end if
       else
          allocate     (labels(size(labels_)))
       end if
       labels  =  0
       basic  => node %basic                          (autoCreate=.true.      )
       label_ =  basic%longIntegerRank0MetaPropertyGet(           nodeLabelsID)
       do i=1,size(labels_)
          if (ibits(label_,i-1,1) == bitTrue) labels(i)=1
       end do
    else
       if (allocated(labels)) deallocate(labels)
    end if
    !$omp end critical(lockNodeLabels)
    return
  end subroutine nodeLabelList

end module Nodes_Labels
