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

  !% Contains a module which implements a merger tree operator which computes the cladistic
  !% information content \cite{thorley_information_1998} of merger trees.

  use, intrinsic :: ISO_C_Binding
  
  !# <mergerTreeOperator name="mergerTreeOperatorInformationContent">
  !#  <description>
  !#   A merger tree operator which computes the cladistic information content
  !#   \cite{thorley_information_1998} of merger trees. This is output to a group in the output
  !#   file with name specified by the {\normalfont \ttfamily [outputGroupName]} parameter. Two
  !#   datasets are written to this group: {\normalfont \ttfamily treeIndex} which gives the
  !#   index of each tree, and {\normalfont \ttfamily informationContent} which gives the
  !#   cladistic information content in units of bits.
  !# </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorInformationContent
     !% A merger tree operator class which computes the cladistic information content
     !% \cite{thorley_information_1998} of merger trees.
     private
     type            (varying_string)                            :: outputGroupName
     integer         (c_size_t      )                            :: treeCount
     integer         (kind=kind_int8), dimension(:), allocatable :: treeIndex
     double precision                , dimension(:), allocatable :: informationContent
   contains
     final     ::             informtionContentDestructor
     procedure :: operate  => informtionContentOperate
     procedure :: finalize => informationContentFinalize
  end type mergerTreeOperatorInformationContent
  
  interface mergerTreeOperatorInformationContent
     !% Constructors for the prune-hierarchy merger tree operator class.
     module procedure informtionContentConstructorParameters
     module procedure informtionContentConstructorInternal
  end interface mergerTreeOperatorInformationContent

contains

  function informtionContentConstructorParameters(parameters)
    !% Constructor for the information content merger tree operator class which takes a parameter set as input.
    implicit none
    type   (mergerTreeOperatorInformationContent)                :: informtionContentConstructorParameters
    type   (inputParameters                     ), intent(in   ) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
        
    !# <inputParameter>
    !#   <name>outputGroupName</name>
    !#   <source>parameters</source>
    !#   <variable>informtionContentConstructorParameters%outputGroupName</variable>
    !#   <defaultValue>var_str('treeInformationContent')</defaultValue>
    !#   <description>The name of an \gls{hdf5} group to which tree information content should be written.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    informtionContentConstructorParameters%treeCount=0_c_size_t
    return
  end function informtionContentConstructorParameters

  function informtionContentConstructorInternal(outputGroupName)
    !% Internal constructor for the information content merger tree operator class.
    use Galacticus_Error
    implicit none
    type(mergerTreeOperatorInformationContent)                :: informtionContentConstructorInternal
    type(varying_string                      ), intent(in   ) :: outputGroupName

    informtionContentConstructorInternal%outputGroupName=outputGroupName
    informtionContentConstructorInternal%treeCount      =0_c_size_t
   return
  end function informtionContentConstructorInternal

  elemental subroutine informtionContentDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorInformationContent), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine informtionContentDestructor

  subroutine informtionContentOperate(self,tree)
    !% Perform a information content operation on a merger tree.
    use Factorials
    use Merger_Trees_Pruning_Utilities
    use Memory_Management
    implicit none
    class           (mergerTreeOperatorInformationContent), intent(inout)               :: self
    type            (mergerTree                          ), intent(inout), target       :: tree
    type            (mergerTree                          )               , pointer      :: treeCurrent
    type            (treeNode                            )               , pointer      :: nodeChild                    , node
    integer         (c_size_t                            )               , parameter    :: treeCountIncrement      =1000
    integer         (kind=kind_int8                      ), allocatable  , dimension(:) :: treeIndexTmp
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
       node                     => treeCurrent%baseNode
       do while (associated(node))
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
          node => node%walkTree()
       end do
       ! Compute logarithm of the possible bifurcations.
       logPossibleBifurcations=Logarithmic_Double_Factorial(2*leafCount-3)
       ! Compute the CIC in bits.
       informationContent=(logPossibleBifurcations-logPermittedBifurcations)/log(2.0d0)
       ! Ensure arrays are large enough to store this tree.
       if (.not.allocated(self%treeIndex)) then
          call Alloc_Array(self%treeIndex         ,[treeCountIncrement])
          call Alloc_Array(self%informationContent,[treeCountIncrement])
       else if (self%treeCount >= size(self%treeIndex)) then
          call Move_Alloc (self%treeIndex         ,      treeIndexTmp                             )
          call Move_Alloc (self%informationContent,      informationContentTmp                    )
          call Alloc_Array(self%treeIndex         ,shape(treeIndexTmp         )+treeCountIncrement)
          call Alloc_Array(self%informationContent,shape(informationContentTmp)+treeCountIncrement)
          self%treeIndex         (1:size(treeIndexTmp         ))=treeIndexTmp
          self%informationContent(1:size(informationContentTmp))=informationContentTmp
          call Dealloc_Array(treeIndexTmp         )
          call Dealloc_Array(informationContentTmp)
       end if
       ! Store the information content.
       self%treeCount                         =self              %treeCount+1
       self%treeIndex         (self%treeCount)=treeCurrent       %index
       self%informationContent(self%treeCount)=informationContent
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine informtionContentOperate

  subroutine informationContentFinalize(self)
    !% Outputs tree information content function.
    use IO_HDF5
    use Galacticus_HDF5
    implicit none
    class(mergerTreeOperatorInformationContent), intent(inout) :: self
    type (hdf5Object                          )                :: informationContentGroup, dataset

    ! Check if we have data to output.
    if (allocated(self%treeIndex)) then
       !$omp critical(HDF5_Access)
       ! Output information content information.
       informationContentGroup=galacticusOutputFile%openGroup(char(self%outputGroupName),'Cladistic information content of trees.')
       call informationContentGroup%writeDataset  (self%treeIndex         ,'treeIndex'                                 )
       call informationContentGroup%writeDataset  (self%informationContent,'informationContent',datasetReturned=dataset)
       call dataset                %writeAttribute('bits'                 ,"units"                                     )
       call dataset                %close         (                                                                    )
       call informationContentGroup%close         (                                                                    )
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine informationContentFinalize
