!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which handles outputting of tree index data to the \glc\ output file.

module Galacticus_Output_Trees_Tree_Weights
  !% Handles outputting of tree index data to the \glc\ output file.
  use Galacticus_Nodes, only : treeNode
  implicit none
  private
  public :: Galacticus_Output_Tree_Tree_Weights, Galacticus_Output_Tree_Tree_Weights_Property_Count, Galacticus_Output_Tree_Tree_Weights_Names

  ! Number of properties.
  integer, parameter :: weightPropertyCount   =1

  ! Flag indicating whether or not tree weight is to be output.
  logical            :: outputTreeWeights

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputTreeInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Weights_Initialize()
    !% Initializes the module by determining whether or not tree weight data should be output.
    use Input_Parameters
    use Galacticus_Error
    implicit none

    if (.not.outputTreeInitialized) then
       !$omp critical(Galacticus_Output_Tree_Weights_Initialize)
       if (.not.outputTreeInitialized) then
          !# <inputParameter>
          !#   <name>outputTreeWeights</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not tree weights should be included in the output.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Flag that module is now initialized.
          outputTreeInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Weights_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Weights_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Tree_Weights_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Tree_Weights</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Tree_Weights_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of tree weight properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode)              , intent(inout) :: thisNode
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    !GCC$ attributes unused :: thisNode, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI, time
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Weights_Initialize()

    if (outputTreeWeights) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='mergerTreeWeight'
       doublePropertyComments(doubleProperty)='Tree weight for this node.'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Tree_Weights_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Tree_Weights_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Tree_Weights</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Tree_Weights_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout) :: thisNode
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: thisNode, time, integerPropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Weights_Initialize()

    if (outputTreeWeights) doublePropertyCount=doublePropertyCount+weightPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Tree_Weights_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Tree_Weights</unitName>
  !#  <sortName>Galacticus_Output_Tree_Tree_Weights</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Tree_Weights(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store link properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Multi_Counters
    implicit none
    double precision                , intent(in   ) :: time
    type            (treeNode      ), intent(inout) :: thisNode
    integer                         , intent(inout) :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                             integerProperty
    integer         (kind=kind_int8), intent(inout) :: integerBuffer    (:,:)
    double precision                , intent(inout) :: doubleBuffer     (:,:)
    type            (multiCounter  ), intent(inout) :: instance
    !GCC$ attributes unused :: integerBufferCount, integerProperty, integerBuffer, time, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Weights_Initialize()

    if (outputTreeWeights) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=thisNode%hostTree%volumeWeight
    end if
    return
  end subroutine Galacticus_Output_Tree_Tree_Weights

end module Galacticus_Output_Trees_Tree_Weights
