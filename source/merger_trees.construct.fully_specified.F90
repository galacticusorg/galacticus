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

! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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
  Implements a merger tree constructor class which constructs a merger tree given a full specification in XML.
  !!}

  use :: FoX_DOM                 , only : node
  use :: IO_XML                  , only : xmlNodeList
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  use :: Merger_Tree_Seeds       , only : mergerTreeSeedsClass

  !![
  <mergerTreeConstructor name="mergerTreeConstructorFullySpecified">
   <description>
    A merger tree constructor class which constructs a merger tree given a full specification in XML. This class will construct
    a merger tree, and set properties of components in each node, using a description read from an XML document. The document
    is specified via the {\normalfont \ttfamily [fileName]} input parameter.
    
    The tree specification document looks as follows:
    \begin{verbatim}
    &lt;!-- Simple initial conditions test case --&gt;
    &lt;initialConditions&gt;
    
      &lt;node&gt;
        &lt;index&gt;2&lt;/index&gt;
        &lt;parent&gt;1&lt;/parent&gt;
        &lt;firstChild&gt;-1&lt;/firstChild&gt;
        &lt;sibling&gt;-1&lt;/sibling&gt;
        &lt;basic&gt;
          &lt;time&gt;1.0&lt;/time&gt;
          &lt;timeLastIsolated&gt;1.0&lt;/timeLastIsolated&gt;
          &lt;mass&gt;1.0e12&lt;/mass&gt;
          &lt;accretionRate&gt;7.9365079e9&lt;/accretionRate&gt;
        &lt;/basic&gt;
        &lt;spin&gt;
          &lt;spin&gt;0.1&lt;/spin&gt;
        &lt;/spin&gt;
        &lt;disk&gt;
          &lt;massGas&gt;1.0e10&lt;/massGas&gt;
          &lt;angularMomentum&gt;1.0e10&lt;/angularMomentum&gt;
          &lt;abundancesGas&gt;
    	&lt;metals&gt;1.0e9&lt;/metals&gt;
    	&lt;Fe&gt;1.0e9&lt;/Fe&gt;
          &lt;/abundancesGas&gt;
        &lt;/disk&gt;
      &lt;/node&gt;
    
      &lt;node&gt;
        &lt;index&gt;1&lt;/index&gt;
        &lt;parent&gt;-1&lt;/parent&gt;
        &lt;firstChild&gt;2&lt;/firstChild&gt;
        &lt;sibling&gt;-1&lt;/sibling&gt;
        &lt;basic&gt;
          &lt;time&gt;13.8&lt;/time&gt;
          &lt;timeLastIsolated&gt;13.8&lt;/timeLastIsolated&gt;
          &lt;mass&gt;1.1e12&lt;/mass&gt;
          &lt;accretionRate&gt;7.8125e9&lt;/accretionRate&gt;
        &lt;/basic&gt;
        &lt;position&gt;
          &lt;position&gt;1.23&lt;/position&gt;
          &lt;position&gt;6.31&lt;/position&gt;
          &lt;position&gt;3.59&lt;/position&gt;
        &lt;/position&gt;
      &lt;/node&gt;
    
    &lt;/initialConditions&gt;
    \end{verbatim}
    The document consists of a set of {\normalfont \ttfamily node} elements, each of which defines a single node in the merger
    tree. Each {\normalfont \ttfamily node} element must specify the {\normalfont \ttfamily index} of the node, along with the
    index of the node's {\normalfont \ttfamily parent}, {\normalfont \ttfamily firstChild}, and {\normalfont \ttfamily
    sibling}.
    
    Each {\normalfont \ttfamily node} element may contain elements which specify the properties of a component in the node. For
    example, a {\normalfont \ttfamily basic} element will specify properties of the ``basic'' component. If multiple elements
    for a given component type are present, then multiple instances of that component will be created in the node.
    
    Within a component definition element scalar properties are set using an element with the same name as that property
    (e.g. {\normalfont \ttfamily mass} in the {\normalfont \ttfamily basic} components in the above example). Rank-1 properties
    are set using a list of elements with the same name as the property (e.g. {\normalfont \ttfamily position} in the
    {\normalfont \ttfamily position} component in the above example).
    
    For composite properties (e.g. abundances), the specification element should contain sub-elements that specify each
    property of the composite. Currently only the {\normalfont \ttfamily abundances} object supports specification in this way,
    as detailed below:
    \begin{description}
     \item [{\normalfont \ttfamily abundances}] (See {\normalfont \ttfamily abundancesGas} in the above example.) The total
     metal content is specified via a {\normalfont \ttfamily metals} element. If other elements are being tracked, their
     content is specified via an element with the short-name of the element (e.g. {\normalfont \ttfamily Fe} for iron).
    \end{description}
   </description>
   <deepCopy>
     <increment variables="document%copyCount" atomic="yes"/>
   </deepCopy>
   <runTimeFileDependencies paths="fileName"/>
  </mergerTreeConstructor>
  !!]
  type, extends(mergerTreeConstructorClass) :: mergerTreeConstructorFullySpecified
     !!{
     A class implementing merger tree construction from a full specification of the tree in XML.
     !!}
     private
     class  (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
     class  (mergerTreeSeedsClass      ), pointer                   :: mergerTreeSeeds_       => null()
     type   (varying_string            )                            :: fileName
     type   (documentContainer         ), pointer                   :: document               => null()
     type   (xmlNodeList               ), allocatable, dimension(:) :: trees
     integer(c_size_t                  )                            :: treeCount
   contains
     final     ::              fullySpecifiedDestructor
     procedure :: construct => fullySpecifiedConstruct
  end type mergerTreeConstructorFullySpecified

  type :: documentContainer
     !!{
     A container for XML document.
     !!}
     private
     type   (node), pointer :: doc       => null()
     integer                :: copyCount =  0
  end type documentContainer

  interface mergerTreeConstructorFullySpecified
     !!{
     Constructors for the \refClass{mergerTreeConstructorFullySpecified} merger tree constructor class.
     !!}
     module procedure fullySpecifiedConstructorParameters
     module procedure fullySpecifiedConstructorInternal
  end interface mergerTreeConstructorFullySpecified

contains

  function fullySpecifiedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeConstructorFullySpecified} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeConstructorFullySpecified)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(randomNumberGeneratorClass         ), pointer       :: randomNumberGenerator_
    class(mergerTreeSeedsClass               ), pointer       :: mergerTreeSeeds_
    type (varying_string                     )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file containing the merger tree specification.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    <objectBuilder class="mergerTreeSeeds"       name="mergerTreeSeeds_"       source="parameters"/>
    !!]
    self=mergerTreeConstructorFullySpecified(fileName,randomNumberGenerator_,mergerTreeSeeds_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    <objectDestructor name="mergerTreeSeeds_"      />
    !!]
    return
  end function fullySpecifiedConstructorParameters

  function fullySpecifiedConstructorInternal(fileName,randomNumberGenerator_,mergerTreeSeeds_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeConstructorFullySpecified} merger tree operator class.
    !!}
    use :: FoX_DOM           , only : parseFile
    use :: Error             , only : Error_Report
    use :: IO_XML            , only : XML_Get_Elements_By_Tag_Name
    use :: File_Utilities    , only : File_Exists
    use :: Display           , only : displayGreen                , displayReset
    use :: ISO_Varying_String, only : varying_string              , var_str     , operator(//)
    implicit none
    type   (mergerTreeConstructorFullySpecified)                        :: self
    type   (varying_string                     ), intent(in   )         :: fileName
    class  (randomNumberGeneratorClass         ), intent(in   ), target :: randomNumberGenerator_
    class  (mergerTreeSeedsClass               ), intent(in   ), target :: mergerTreeSeeds_
    integer                                                             :: ioErr
    type   (varying_string                     )                        :: message
    !![
    <constructorAssign variables="fileName, *randomNumberGenerator_, *mergerTreeSeeds_"/>
    !!]

    !$omp critical (FoX_DOM_Access)
    if (.not.associated(self%document)) allocate(self%document)
    ! Parse the merger tree file.
    self%document%doc => parseFile(char(self%fileName),iostat=ioErr)
    if (ioErr /= 0) then
       message=var_str("unable to read or parse fully-specified merger tree file ")//"'"//self%fileName//"'"
       if (File_Exists(self%fileName)) then
          message=message//char(10)//displayGreen()//"HELP:"//displayReset()//" check that the XML in this file is valid (e.g. `xmllint --noout "//self%fileName//"` will display any XML errors"
       else
          message=message//" - file does not exist"
       end if
       call Error_Report(message//{introspection:location})
    end if
    self%document%copyCount = 1
    ! Get the list of trees.
    call XML_Get_Elements_By_Tag_Name(self%document%doc,"tree",self%trees)
    ! Count the number of trees.
    self%treeCount=size(self%trees)
    if (self%treeCount <= 0) call Error_Report('no trees were specified'//{introspection:location})
    !$omp end critical (FoX_DOM_Access)
    return
  end function fullySpecifiedConstructorInternal

  subroutine fullySpecifiedDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeConstructorFullySpecified} merger tree constructor class.
    !!}
    use :: FoX_DOM, only : destroy
    implicit none
    type(mergerTreeConstructorFullySpecified), intent(inout) :: self

    ! Reduce the count of document copies.
    if (associated(self%document)) then
       !$omp atomic
       self%document%copyCount=self%document%copyCount-1
       ! Destroy the XML document only if the count of document copies decreases to 0.
       if (self%document%copyCount == 0) then
          !$omp critical (FoX_DOM_Access)
          call destroy(self%document%doc)
          if (associated(self%document)) deallocate(self%document)
          !$omp end critical (FoX_DOM_Access)
       end if
    end if
    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    <objectDestructor name="self%mergerTreeSeeds_"      />
    !!]
    return
  end subroutine fullySpecifiedDestructor

  function fullySpecifiedConstruct(self,treeNumber,finished) result(tree)
    !!{
    Construct a fully-specified merger tree.
    !!}
    use            :: Display         , only : displayIndent               , displayUnindent, displayVerbosity, verbosityLevelInfo
    use            :: FoX_DOM         , only : node
    use            :: Error           , only : Error_Report
    use            :: Galacticus_Nodes, only : mergerTree                  , treeNode       , treeNodeList
    use            :: IO_XML          , only : XML_Get_Elements_By_Tag_Name
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Kind_Numbers    , only : kind_int8
    implicit none
    type   (mergerTree                         ), pointer                     :: tree
    class  (mergerTreeConstructorFullySpecified), intent(inout)               :: self
    integer(c_size_t                           ), intent(in   )               :: treeNumber
    logical                                     , intent(  out)               :: finished
    type   (treeNodeList                       ), allocatable  , dimension(:) :: nodeArray
    type   (node                               ), pointer                     :: treeDefinition, nodeDefinition
    type   (xmlNodeList                        ), allocatable  , dimension(:) :: nodes
    integer                                                                   :: i             , nodeCount
    integer(kind_int8                          )                              :: indexValue

    ! Read one tree.
    if (treeNumber > 0_c_size_t .and. treeNumber <= self%treeCount) then
       !$omp critical (FoX_DOM_Access)
       ! Select one tree.
       treeDefinition => self%trees(int(treeNumber-1))%element
       ! Get the list of nodes in this tree.
       call XML_Get_Elements_By_Tag_Name(treeDefinition,"node",nodes)
       nodeCount=size(nodes)
       if (nodeCount <= 0) call Error_Report('no nodes were specified'//{introspection:location})
       ! Create the tree.
       allocate(tree)
       ! Create an array of nodes.
       allocate(nodeArray(nodeCount))
       ! Iterate over nodes.
       do i=1,nodeCount
          ! Create the node.
          nodeArray(i)%node => treeNode(hostTree=tree)
          ! Get the node definition.
          nodeDefinition => nodes(i-1)%element
          ! Assign an index to the node.
          call nodeArray(i)%node%indexSet(indexNode(nodeDefinition,'index'))
       end do
       !$omp end critical (FoX_DOM_Access)
       ! Initialize the tree root to null.
       tree%index             =  int(treeNumber)
       tree%volumeWeight      =  1.0d0
       tree%initializedUntil  =  0.0d0
       tree%isTreeInitialized =  .false.
       tree%firstTree         => tree
       tree%nodeBase          => null()
       tree%event             => null()
       call tree%properties%initialize()
       ! Restart the random number sequence.
       allocate(tree%randomNumberGenerator_,mold=self%randomNumberGenerator_)
       !$omp critical(mergerTreeConstructFullySpecifiedDeepCopyReset)
       !![
       <deepCopyReset variables="self%randomNumberGenerator_"/>
       <deepCopy source="self%randomNumberGenerator_" destination="tree%randomNumberGenerator_"/>
       <deepCopyFinalize variables="tree%randomNumberGenerator_"/>
       !!]
       !$omp end critical(mergerTreeConstructFullySpecifiedDeepCopyReset)
       call self                 %randomSequenceNonDeterministicWarn(tree)
       call self%mergerTreeSeeds_%set                               (tree)
       ! Begin writing report.
       call displayIndent('Initial conditions of fully-specified tree',verbosityLevelInfo)
       ! Iterate over nodes.
       do i=1,nodeCount
          ! Get the node definition.
          !$omp critical (FoX_DOM_Access)
          nodeDefinition => nodes(i-1)%element
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
             if (associated(tree%nodeBase)) call Error_Report('multiple root nodes found in the tree'//{introspection:location})
             tree%nodeBase => nodeArray(i)%node
          end if
          ! Build components.
          call nodeArray(i)%node%componentBuilder(nodeDefinition)
          ! Dump the node.
          if (displayVerbosity() > verbosityLevelInfo) call nodeArray(i)%node%serializeASCII()
       end do
       ! Finish writing report.
       call displayUnindent('done',verbosityLevelInfo)
       ! Destroy the node array.
       deallocate(nodeArray)
       ! Check that we found a root node.
       if (.not.associated(tree%nodeBase)) call Error_Report('no root node was found'//{introspection:location})
    else
       nullify(tree)
    end if
    finished=.not.associated(tree)
    return

  contains

    function indexNode(nodeDefinition,indexType,required)
      !!{
      Extract and return an index from a node definition as used when constructing fully-specified merger trees.
      !!}
      use :: FoX_Dom     , only : node                        , extractDataContent
      use :: Kind_Numbers, only : kind_int8
      use :: IO_XML      , only : XML_Get_Elements_By_Tag_Name
      use :: Error       , only : Error_Report
      implicit none
      integer  (kind=kind_int8)                             :: indexNode
      type     (node          ), intent(in   ), pointer     :: nodeDefinition
      character(len=*         ), intent(in   )              :: indexType
      logical                  , intent(in   ), optional    :: required
      type     (xmlNodeList   ), dimension(:) , allocatable :: indexElements
      type     (node          )               , pointer     :: indexElement
      integer                                               :: indexValue
      !![
      <optionalArgument name="required" defaultsTo=".true." />
      !!]
      
      ! Find all matching tags.
      call XML_Get_Elements_By_Tag_Name(nodeDefinition,indexType,indexElements)
      if (size(indexElements) > 1) call Error_Report('multiple indices specified'//{introspection:location})
      if (size(indexElements) < 1) then
         if (required_) then
            call Error_Report('required index not specified'//{introspection:location})
         else
            indexNode=-1
            return
         end if
      end if
      ! Get the index element.
      indexElement => indexElements(0)%element
      ! Extract the value.
      call extractDataContent(indexElement,indexValue)
      ! Transfer to function result.
      indexNode=indexValue
      return
    end function indexNode

    function nodeLookup(nodeArray,indexValue) result (node)
      !!{
      Find the position of a node in the {\normalfont \ttfamily nodeArray} array given its {\normalfont \ttfamily indexValue}.
      !!}
      use :: Error           , only : Error_Report
      use :: Galacticus_Nodes, only : treeNode    , treeNodeList
      use :: Kind_Numbers    , only : kind_int8
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
      if (.not.associated(node)) call Error_Report('unable to find requested node'//{introspection:location})
      return
    end function nodeLookup

  end function fullySpecifiedConstruct
