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

  !% Implements a merger tree constructor class which constructs a merger tree given a full specification in XML.

  use FoX_DOM

  !# <mergerTreeConstructor name="mergerTreeConstructorFullySpecified">
  !#  <description>Merger tree constructor class which constructs a merger tree given a full specification in XML.</description>
  !#  <deepCopy>
  !#    <increment variables="document%copyCount" atomic="yes"/>
  !#  </deepCopy>
  !# </mergerTreeConstructor>
  type, extends(mergerTreeConstructorClass) :: mergerTreeConstructorFullySpecified
     !% A class implementing merger tree construction from a full specification of the tree in XML.
     private
     type   (varying_string   )          :: fileName
     type   (documentContainer), pointer :: document  => null()
     type   (nodeList         ), pointer :: trees     => null()
     integer(c_size_t         )          :: treeCount
   contains
     final     ::              fullySpecifiedDestructor
     procedure :: construct => fullySpecifiedConstruct
  end type mergerTreeConstructorFullySpecified

  type :: documentContainer
     !% A container for XML docoment.
     private
     type   (node             ), pointer :: doc       => null()
     integer                             :: copyCount
  end type documentContainer

  interface mergerTreeConstructorFullySpecified
     !% Constructors for the {\normalfont \ttfamily fullySpecified} merger tree constructor class.
     module procedure fullySpecifiedConstructorParameters
     module procedure fullySpecifiedConstructorInternal
  end interface mergerTreeConstructorFullySpecified

contains
  
  function fullySpecifiedConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily fullySpecified} merger tree operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(mergerTreeConstructorFullySpecified)                :: self
    type(inputParameters                    ), intent(inout) :: parameters
    type(varying_string                     )                :: fileName

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the file containing the merger tree specification.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    self=mergerTreeConstructorFullySpecified(fileName)
    !# <inputParametersValidate source="parameters"/>
    return
  end function fullySpecifiedConstructorParameters

  function fullySpecifiedConstructorInternal(fileName) result(self)
    !% Internal constructor for the {\normalfont \ttfamily fullySpecified} merger tree operator class.
    use Galacticus_Error
    implicit none
    type(mergerTreeConstructorFullySpecified)                :: self
    type(varying_string                     ), intent(in   ) :: fileName
    integer                                                  :: ioErr
    !# <constructorAssign variables="fileName"/>

    !$omp critical (FoX_DOM_Access)
    if (.not.associated(self%document)) allocate(self%document)
    ! Parse the merger tree file.
    self%document%doc => parseFile(char(self%fileName),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('unable to read or parse fully-specified merger tree file'//{introspection:location})
    self%document%copyCount = 1
    ! Get the list of trees.
    self%trees => getElementsByTagname(self%document%doc,"tree")
    ! Count the number of trees.
    self%treeCount=getLength(self%trees)
    if (self%treeCount <= 0) call Galacticus_Error_Report('no trees were specified'//{introspection:location})
    !$omp end critical (FoX_DOM_Access)
    return    
  end function fullySpecifiedConstructorInternal

  subroutine fullySpecifiedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily fullySpecified} merger tree constructor class.
    implicit none
    type(mergerTreeConstructorFullySpecified), intent(inout) :: self

    ! Reduce the count of document copies.
    !$omp atomic
    self%document%copyCount=self%document%copyCount-1
    ! Destroy the XML document only if the count of document copies decreases to 0.
    if (self%document%copyCount == 0) then
       !$omp critical (FoX_DOM_Access)
       call destroy(self%document%doc)
       if (associated(self%document)) deallocate(self%document)
       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine fullySpecifiedDestructor
  
  function fullySpecifiedConstruct(self,treeNumber) result(tree)
    !% Construct a fully-specified merger tree.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Nodes   , only : treeNodeList
    use               FoX_DOM
    use               Kind_Numbers
    use               Galacticus_Error
    use               Memory_Management
    use               Galacticus_Display
    use               Pseudo_Random
    implicit none
    type            (mergerTree                         ), pointer                     :: tree
    class           (mergerTreeConstructorFullySpecified), intent(inout)               :: self
    integer         (c_size_t                           ), intent(in   )               :: treeNumber
    type            (treeNodeList                       ), allocatable  , dimension(:) :: nodeArray
    type            (node                               ), pointer                     :: treeDefinition, nodeDefinition
    type            (nodeList                           ), pointer                     :: nodes
    integer                                                                            :: i             , nodeCount
    integer         (kind_int8                          )                              :: indexValue
    double precision                                                                   :: uniformRandom

    ! Read one tree.
    if (treeNumber > 0_c_size_t .and. treeNumber <= self%treeCount) then
       !$omp critical (FoX_DOM_Access)
       ! Select one tree.
       treeDefinition => item(self%trees,int(treeNumber-1))
       ! Get the list of nodes in this tree.
       nodes => getElementsByTagname(treeDefinition,"node")
       nodeCount=getLength(nodes)
       if (nodeCount <= 0) call Galacticus_Error_Report('no nodes were specified'//{introspection:location})
       ! Create the tree.
       allocate(tree)
       ! Create an array of nodes.
       allocate(nodeArray(nodeCount))
       call Memory_Usage_Record(sizeof(nodeArray))
       ! Iterate over nodes.
       do i=1,nodeCount
          ! Create the node.
          nodeArray(i)%node => treeNode(hostTree=tree)
          ! Get the node definition.
          nodeDefinition => item(nodes,i-1)
          ! Assign an index to the node.
          call nodeArray(i)%node%indexSet(indexNode(nodeDefinition,'index'))
       end do
       !$omp end critical (FoX_DOM_Access)
       ! Initialize the tree root to null.
       tree%index            =  int(treeNumber)
       tree%initializedUntil =  0.0d0
       tree%firstTree        => tree
       tree%baseNode         => null()
       tree%event            => null()
       call tree%properties%initialize()
       ! Restart the random number sequence.
       tree%randomNumberGenerator=pseudoRandom()
       uniformRandom=tree%randomNumberGenerator%uniformSample(ompThreadOffset=.false.,mpiRankOffset=.false.,incrementSeed=int(tree%index))
       ! Begin writing report.
       call Galacticus_Display_Indent('Initial conditions of fully-specified tree',verbosityInfo)
       ! Iterate over nodes.
       do i=1,nodeCount
          ! Get the node definition.
          !$omp critical (FoX_DOM_Access)
          nodeDefinition => item(nodes,i-1)
          ! Build parent pointers.
          indexValue=indexNode(nodeDefinition,'parent'        )
          nodeArray(i)%node%parent         => nodeLookup(nodeArray,indexValue)
          ! Build child pointers.
          indexValue=indexNode(nodeDefinition,'firstChild'    )
          nodeArray(i)%node%firstChild     => nodeLookup(nodeArray,indexValue)
          ! Build sibling pointers.
          indexValue=indexNode(nodeDefinition,'sibling'       )
          nodeArray(i)%node%sibling        => nodeLookup(nodeArray,indexValue)
          ! Build satellite pointers.
          indexValue=indexNode(nodeDefinition,'firstSatellite',required=.false.)
          nodeArray(i)%node%firstSatellite => nodeLookup(nodeArray,indexValue)
          !$omp end critical (FoX_DOM_Access)
          ! Assign the tree root node if this node has no parent.
          if (.not.associated(nodeArray(i)%node%parent)) then
             if (associated(tree%baseNode)) call Galacticus_Error_Report('multiple root nodes found in the tree'//{introspection:location})
             tree%baseNode => nodeArray(i)%node
          end if
          ! Build components.
          call nodeArray(i)%node%componentBuilder(nodeDefinition)
          ! Dump the node.
          if (Galacticus_Verbosity_Level() > verbosityInfo) call nodeArray(i)%node%serializeASCII()
       end do
       ! Finish writing report.
       call Galacticus_Display_Unindent('done',verbosityInfo)
       ! Destroy the node array.
       deallocate(nodeArray)
       ! Check that we found a root node.
       if (.not.associated(tree%baseNode)) call Galacticus_Error_Report('no root node was found'//{introspection:location})
    else
       nullify(tree)
    end if
    return

  contains

    function indexNode(nodeDefinition,indexType,required)
      !% Extract and return an index from a node definition as used when constructing fully-specified merger trees.
      use FoX_Dom
      use Galacticus_Error
      use Kind_Numbers
      implicit none
      integer  (kind=kind_int8)                          :: indexNode
      type     (node          ), intent(in   ), pointer  :: nodeDefinition
      character(len=*         ), intent(in   )           :: indexType
      logical                  , intent(in   ), optional :: required
      type     (nodeList      )               , pointer  :: indexElements
      type     (node          )               , pointer  :: indexElement
      integer                                            :: indexValue
      !# <optionalArgument name="required" defaultsTo=".true." />

      ! Find all matching tags.
      indexElements => getElementsByTagname(nodeDefinition,indexType)
      if (getLength(indexElements) > 1) call Galacticus_Error_Report('multiple indices specified'//{introspection:location})
      if (getLength(indexElements) < 1) then
         if (required_) then
            call Galacticus_Error_Report('required index not specified'//{introspection:location})
         else
            indexNode=-1
            return
         end if
      end if
      ! Get the index element.
      indexElement => item(indexElements,0)
      ! Extract the value.
      call extractDataContent(indexElement,indexValue)
      ! Transfer to function result.
      indexNode=indexValue
      return
    end function indexNode

    function nodeLookup(nodeArray,indexValue) result (node)
      !% Find the position of a node in the {\normalfont \ttfamily nodeArray} array given its {\normalfont \ttfamily indexValue}.
      use Galacticus_Nodes, only : treeNode, treeNodeList
      use Kind_Numbers
      use Galacticus_Error
      implicit none
      type   (treeNode      ), pointer                     :: node
      type   (treeNodeList  ), dimension(:), intent(in   ) :: nodeArray
      integer(kind=kind_int8)              , intent(in   ) :: indexValue
      integer                                              :: i

      node => null()    
      if (indexValue < 0_kind_int8) return
      do i=1,size(nodeArray)
         if (nodeArray(i)%node%index() == indexValue) then
            node => nodeArray(i)%node
            exit
         end if
      end do
      if (.not.associated(node)) call Galacticus_Error_Report('unable to find requested node'//{introspection:location})
      return
    end function nodeLookup

  end function fullySpecifiedConstruct
