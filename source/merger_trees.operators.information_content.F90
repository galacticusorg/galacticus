!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a merger tree operator which computes the cladistic
  information content \cite{thorley_information_1998} of merger trees.
  !!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t
  use            :: Kind_Numbers , only : kind_int8

  !![
  <mergerTreeOperator name="mergerTreeOperatorInformationContent">
   <description>
    A merger tree operator which computes the cladistic information content
    \cite{thorley_information_1998} of merger trees. This is output to a group in the output
    file with name specified by the {\normalfont \ttfamily [outputGroupName]} parameter. Two
    datasets are written to this group: {\normalfont \ttfamily treeIndex} which gives the
    index of each tree, and {\normalfont \ttfamily informationContent} which gives the
    cladistic information content in units of bits.
  </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorInformationContent
     !!{
     A merger tree operator class which computes the cladistic information content
     \cite{thorley_information_1998} of merger trees.
     !!}
     private
     type            (varying_string)                            :: outputGroupName
     integer         (c_size_t      )                            :: treeCount
     integer         (kind=kind_int8), dimension(:), allocatable :: treeIndex
     double precision                , dimension(:), allocatable :: informationContent
   contains
     procedure :: operatePreEvolution => informationContentOperatePreEvolution
     procedure :: finalize            => informationContentFinalize
  end type mergerTreeOperatorInformationContent

  interface mergerTreeOperatorInformationContent
     !!{
     Constructors for the prune-hierarchy merger tree operator class.
     !!}
     module procedure informationContentConstructorParameters
     module procedure informationContentConstructorInternal
  end interface mergerTreeOperatorInformationContent

contains

  function informationContentConstructorParameters(parameters) result(self)
    !!{
    Constructor for the information content merger tree operator class which takes a parameter set as input.
    !!}
    implicit none
    type(mergerTreeOperatorInformationContent)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    type(varying_string                      )                :: outputGroupName

    !![
    <inputParameter>
      <name>outputGroupName</name>
      <source>parameters</source>
      <defaultValue>var_str('treeInformationContent')</defaultValue>
      <description>The name of an \gls{hdf5} group to which tree information content should be written.</description>
    </inputParameter>
    !!]
    self=mergerTreeOperatorInformationContent(outputGroupName)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function informationContentConstructorParameters

  function informationContentConstructorInternal(outputGroupName) result(self)
    !!{
    Internal constructor for the information content merger tree operator class.
    !!}
    implicit none
    type(mergerTreeOperatorInformationContent)                :: self
    type(varying_string                      ), intent(in   ) :: outputGroupName
    !![
    <constructorAssign variables="outputGroupName"/>
    !!]
    
    self%treeCount=0_c_size_t
   return
  end function informationContentConstructorInternal

  subroutine informationContentOperatePreEvolution(self,tree)
    !!{
    Perform a information content operation on a merger tree.
    !!}
    use :: Factorials         , only : Logarithmic_Double_Factorial
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    class           (mergerTreeOperatorInformationContent), intent(inout), target       :: self
    type            (mergerTree                          ), intent(inout), target       :: tree
    type            (mergerTree                          )               , pointer      :: treeCurrent
    type            (treeNode                            )               , pointer      :: nodeChild                    , node
    integer         (c_size_t                            )               , parameter    :: treeCountIncrement      =1000
    integer         (kind=kind_int8                      ), allocatable  , dimension(:) :: treeIndexTmp
    type            (mergerTreeWalkerIsolatedNodes       )                              :: treeWalker
    double precision                                      , allocatable  , dimension(:) :: informationContentTmp
    integer                                                                             :: childCount                   , leafCount
    double precision                                                                    :: logPermittedBifurcations     , logPossibleBifurcations, &
         &                                                                                 informationContent

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Walk the tree, counting the number of leaves and accumulated the log of the number of permitted bifurcations.
       leafCount                =  0
       logPermittedBifurcations =  0.0d0
       treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
       do while (treeWalker%next(node))
          ! Check for leaf nodes.
          if (.not.associated(node%firstChild)) then
             ! Increment leaf node counter.
             leafCount=leafCount+1
          else
             ! Count child nodes.
             childCount =  0
             nodeChild  => node%firstChild
             do while (associated(nodeChild))
                childCount =  childCount       +1
                nodeChild  => nodeChild%sibling
             end do
             ! Increment the number of permitted bifurcations based on this number of children.
             logPermittedBifurcations=logPermittedBifurcations+Logarithmic_Double_Factorial(2*childCount-3)
          end if
       end do
       ! Compute logarithm of the possible bifurcations.
       logPossibleBifurcations=Logarithmic_Double_Factorial(2*leafCount-3)
       ! Compute the CIC in bits.
       informationContent=(logPossibleBifurcations-logPermittedBifurcations)/log(2.0d0)
       ! Ensure arrays are large enough to store this tree.
       if (.not.allocated(self%treeIndex)) then
          allocate(self%treeIndex         (treeCountIncrement))
          allocate(self%informationContent(treeCountIncrement))
       else if (self%treeCount >= size(self%treeIndex)) then
          call move_alloc(self%treeIndex         ,     treeIndexTmp                             )
          call move_alloc(self%informationContent,     informationContentTmp                    )
          allocate       (self%treeIndex         (size(treeIndexTmp         )+treeCountIncrement))
          allocate       (self%informationContent(size(informationContentTmp)+treeCountIncrement))
          self%treeIndex         (1:size(treeIndexTmp         ))=treeIndexTmp
          self%informationContent(1:size(informationContentTmp))=informationContentTmp
          deallocate(treeIndexTmp         )
          deallocate(informationContentTmp)
       end if
       ! Store the information content.
       self%treeCount                         =self              %treeCount+1
       self%treeIndex         (self%treeCount)=treeCurrent       %index
       self%informationContent(self%treeCount)=informationContent
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine informationContentOperatePreEvolution

  subroutine informationContentFinalize(self)
    !!{
    Outputs tree information content function.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5       , only : hsize_t
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class  (mergerTreeOperatorInformationContent), intent(inout) :: self
    integer(hsize_t                             ), parameter     :: chunkSize              =1024_hsize_t
    type   (hdf5Object                          )                :: informationContentGroup             , dataset
    logical                                                      :: preexisting

    ! Check if we have data to output.
    if (allocated(self%treeIndex)) then
       !$ call hdf5Access%set()
       ! Output information content information.
       informationContentGroup   =outputFile             %openGroup     (char(self%outputGroupName),'Cladistic information content of trees.'                                       )
       preexisting               =informationContentGroup%hasDataset    (                           'treeIndex'                                                                     )
       call                       informationContentGroup%writeDataset  (self%treeIndex            ,'treeIndex'                                 ,appendTo=.true.,chunkSize=chunkSize)
       call                       informationContentGroup%writeDataset  (self%informationContent   ,'informationContent',datasetReturned=dataset,appendTo=.true.,chunkSize=chunkSize)
       if (.not.preexisting) call dataset                %writeAttribute('bits'                    ,"units"                                                                         )
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine informationContentFinalize
