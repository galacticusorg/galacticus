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
Contains a module which implements reading of parameters from an XML data file.
!!}

! Specify an explicit dependence on the git2.o object file.
!: $(BUILDPATH)/git2.o

module Input_Parameters
  !!{
  Implements reading of parameters from an XML file.
  !!}
  use, intrinsic :: ISO_C_Binding     , only : c_char        , c_int
  use            :: FoX_dom           , only : node
  use            :: Function_Classes  , only : functionClass
  use            :: IO_HDF5           , only : hdf5Object
  use            :: ISO_Varying_String, only : varying_string
  use            :: Kind_Numbers      , only : kind_int8
  use            :: String_Handling   , only : char
  use            :: Hashes            , only : integerHash
  use            :: Locks             , only : ompLock
  private
  public :: inputParameters, inputParameter, inputParameterList
  
  !![
  <generic identifier="Type">
   <instance label="Logical"        intrinsic="logical"                                          outputConverter="regEx¦(.*)¦char($1)¦"/>
   <instance label="Integer"        intrinsic="integer"                                          outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="Long"           intrinsic="integer         (kind=kind_int8)"                 outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="Double"         intrinsic="double precision"                                 outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="Character"      intrinsic="character       (len=*         )"                 outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="VarStr"         intrinsic="type            (varying_string)"                 outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="LogicalRank1"   intrinsic="logical                         , dimension(:)"   outputConverter="regEx¦(.*)¦char($1)¦"/>
   <instance label="IntegerRank1"   intrinsic="integer                         , dimension(:)"   outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="DoubleRank1"    intrinsic="double precision                , dimension(:)"   outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="LongRank1"      intrinsic="integer         (kind=kind_int8), dimension(:)"   outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="CharacterRank1" intrinsic="character       (len=*         ), dimension(:)"   outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="VarStrRank1"    intrinsic="type            (varying_string), dimension(:)"   outputConverter="regEx¦(.*)¦$1¦"      />
   <instance label="DoubleRank2"    intrinsic="double precision                , dimension(:,:)" outputConverter="regEx¦(.*)¦$1¦"      />
  </generic>
  !!]

  !![
  <generic identifier="cType">
   <instance label="Double" intrinsic="real   (kind=c_double)" />
   <instance label="Long"   intrinsic="integer(kind=c_long  )" />
  </generic>
  !!]

  type :: genericObjectList
     !!{
     A list-type for unlimited polymorphic pointers.
     !!}
     private
     class(functionClass), pointer :: object => null()
  end type genericObjectList

  type :: inputParameter
     !!{
     A class to handle input parameters for \glc.
     !!}
     private
     type   (node             ), pointer                     :: content         => null()
     type   (inputParameter   ), pointer    , public         :: parent          => null() , firstChild => null() , &
          &                                                     sibling         => null() , referenced => null()
     type   (genericObjectList), allocatable, dimension(:,:) :: objects
     type   (varying_string   )                              :: contentOriginal
     logical                                                 :: created         =  .false., removed    =  .false., &
          &                                                     evaluated       =  .false., active     =  .true. , &
          &                                                     activeEvaluated =  .false.
   contains
     !![
     <methods>
       <method description="Return true if this is a valid parameter node." method="isParameter" />
       <method description="Return true the object corresponding to this parameter has been created." method="objectCreated" />
       <method description="Return a pointer to the object corresponding to this parameter." method="objectGet" />
       <method description="Set a pointer to the object corresponding to this parameter." method="objectSet" />
       <method description="Reset all objects in this parameter and any sub-parameters." method="reset" />
       <method description="Set the value of this parameter." method="set" />
       <method description="Return the value of this parameter in a simple textual context." method="get" />
       <method description="Destroy this parameter and all subparameters." method="destroy" />
     </methods>
     !!]
     procedure :: isParameter   => inputParameterIsParameter
     procedure :: objectCreated => inputParameterObjectCreated
     procedure :: objectGet     => inputParameterObjectGet
     procedure :: objectSet     => inputParameterObjectSet
     procedure :: destroy       => inputParameterDestroy
     procedure :: reset         => inputParameterReset
     procedure ::                  inputParameterSetDouble
     procedure ::                  inputParameterSetVarStr
     generic   :: set           => inputParameterSetDouble
     generic   :: set           => inputParameterSetVarStr
     procedure :: get           => inputParameterGet
  end type inputParameter

  !![
  <deepCopyActions class="inputParameters">
   <inputParameters>
    <methodCall method="lockReinitialize"/>
   </inputParameters>
  </deepCopyActions>
  !!]
    
  type :: inputParameters
     private
     type   (node           ), pointer, public :: document               => null()
     type   (node           ), pointer         :: rootNode               => null()
     type   (hdf5Object     )                  :: outputParameters                 , outputParametersContainer
     type   (inputParameter ), pointer, public :: parameters             => null()
     type   (inputParameters), pointer, public :: parent                 => null() , original                  => null()
     logical                                   :: outputParametersCopied =  .false., outputParametersTemporary = .false., &
          &                                       isNull                 =  .false., strict                    = .false.
     type   (integerHash    ), allocatable     :: warnedDefaults
     type   (ompLock        ), pointer         :: lock                   => null()
   contains
     !![
     <methods>
       <method description="Build a tree of {\normalfont \ttfamily inputParameter} objects from the structure of an XML parameter file." method="buildTree" />
       <method description="Resolve references in the tree of {\normalfont \ttfamily inputParameter} objects." method="resolveReferences" />
       <method description="Evaluate conditionals in the tree of {\normalfont \ttfamily inputParameter} objects." method="evaluateConditionals" />
       <method description="Return the HDF5 group to which this parameters content will be written." method="parametersGroup" />
       <method description="Open an output group for parameters in the given HDF5 object." method="parametersGroupOpen" />
       <method description="Copy the HDF5 output group for parameters from another parameters object." method="parametersGroupCopy" />
       <method description="Check that a given parameter name is a valid name, aborting if not." method="validateName" />
       <method description="Check that parameters are valid and, optionally, check if they match expected names in the provided list." method="checkParameters" />
       <method description="Return the XML node containing the named parameter." method="node" />
       <method description="Return true if the named parameter is present in the set." method="isPresent" />
       <method description="Return a count of the number copies of the named parameter. If the parameter is not present, this function aborts, unless {\normalfont \ttfamily zeroIfNotPresent} is set to {\normalfont \ttfamily true}, in which case a result of 0 is returned." method="copiesCount" />
       <method description="Return a count of the number of values in the named parameter. If the parameter is not present, this function aborts, unless {\normalfont \ttfamily zeroIfNotPresent} is set to {\normalfont \ttfamily true}, in which case a result of 0 is returned." method="count" />
       <method description="Return the set of subparameters of the named parameter." method="subParameters" />
       <method description="Return the parent parameters given the path to a parameter" method="findParent"/>
       <method description="Return the value of a parameter specified by name or XML node. A default value can be specified only if the parameter is specified by name. Supported types include rank-0 and rank-1 logicals, integers, long integers, doubles, characters, and varying strings." method="value" />
       <method description="Serialize input parameters to a string." method="serializeToString" />
       <method description="Serialize input parameters to an XML file." method="serializeToXML" />
       <method description="Add a parameter." method="addParameter" />
       <method description="Reset all objects in this parameter set." method="reset" />
       <method description="Destroy the parameters document." method="destroy" />
       <method description="Return the path to this parameters object." method="path" />
       <method description="Reinitialize lock." method="lockReinitialize" />
     </methods>
     !!]
     final     ::                         inputParametersFinalize
     procedure :: buildTree            => inputParametersBuildTree
     procedure :: resolveReferences    => inputParametersResolveReferences
     procedure :: evaluateConditionals => inputParametersEvaluateConditionals
     procedure :: destroy              => inputParametersDestroy
     procedure :: parametersGroup      => inputParametersParametersGroup
     procedure :: parametersGroupOpen  => inputParametersParametersGroupOpen
     procedure :: parametersGroupCopy  => inputParametersParametersGroupCopy
     procedure :: validateName         => inputParametersValidateName
     procedure :: checkParameters      => inputParametersCheckParameters
     procedure :: node                 => inputParametersNode
     procedure :: isPresent            => inputParametersIsPresent
     procedure :: copiesCount          => inputParametersCopiesCount
     procedure :: count                => inputParametersCount
     procedure :: subParameters        => inputParametersSubParameters
     procedure :: findParent           => inputParametersFindParent
     procedure ::                         inputParametersValueName{Type¦label}
     procedure ::                         inputParametersValueNode{Type¦label}
     generic   :: value                => inputParametersValueName{Type¦label}
     generic   :: value                => inputParametersValueNode{Type¦label}
     procedure :: serializeToString    => inputParametersSerializeToString
     procedure :: serializeToXML       => inputParametersSerializeToXML
     procedure :: addParameter         => inputParametersAddParameter
     procedure :: reset                => inputParametersReset
     procedure :: path                 => inputParametersPath
     procedure :: lockReinitialize     => inputParametersLockReinitialize
  end type inputParameters

  interface inputParameters
     !!{
     Constructors for the {\normalfont \ttfamily inputParameters} class.
     !!}
     module procedure inputParametersConstructorVarStr
     module procedure inputParametersConstructorFileChar
     module procedure inputParametersConstructorNode
     module procedure inputParametersConstructorCopy
     module procedure inputParametersConstructorNull
  end interface inputParameters

  ! Define a type to hold lists of parameters (and values) prior to output.
  type :: inputParameterList
     !!{
     A class to hold lists of parameters (and values) prior to output.
     !!}
     integer                                            :: count
     type   (varying_string), allocatable, dimension(:) :: name , value
   contains
     !![
     <methods>
       <method description="Serialize a list of input parameters to an XML document." method="serializeToXML"/>
       <method description="Add a parameter and value to the list."                   method="add"           />
     </methods>
     !!]
     final     ::                      inputParameterListDestructor
     procedure :: add               => inputParameterListAdd
     procedure :: serializeToXML    => inputParameterListSerializeToXML
  end type inputParameterList

  interface inputParameterList
     !!{
     Constructors for {\normalfont \ttfamily inputParameterList} objects.
     !!}
     module procedure inputParameterListConstructor
  end interface inputParameterList

  !![
  <enumeration>
   <name>inputParameterErrorStatus</name>
   <description>Error status codes used by the input parameters module.</description>
   <entry label="success"        />
   <entry label="notPresent"     />
   <entry label="parse"          />
   <entry label="emptyValue"     />
   <entry label="ambiguousValue" />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>inputParameterType</name>
   <description>Types for input parameters.</description>
   <entry label="double" />
   <entry label="integer"/>
   <entry label="text"   />
  </enumeration>
  !!]

  ! Maximum length allowed for parameter entries.
  integer, parameter :: parameterLengthMaximum=1024

  ! Interface to the (auto-generated) knownParameterNames() function.
  interface
     subroutine knownParameterNames(names)
       import varying_string
       type(varying_string), dimension(:), allocatable, intent(inout) :: names
     end subroutine knownParameterNames
  end interface

#ifdef GIT2AVAIL
  interface
     function gitDescendantOf(repoPath,commitHash,ancestorHash) bind(c,name='gitDescendantOf')
       !!{
       Template for a C function that returns whether a commit is an ancestor of another commit.
       !!}
       import c_char, c_int
       integer  (c_int )                :: gitDescendantOf
       character(c_char)                :: repoPath
       character(c_char), dimension(41) :: commitHash     , ancestorHash
     end function gitDescendantOf
  end interface
#endif

contains

  function inputParametersConstructorNull() result(self)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class creating a null instance.
    !!}
    use :: FoX_dom, only : createDocument, getDocumentElement, getImplementation, setLiveNodeLists
    implicit none
    type(inputParameters) :: self

    allocate(self%warnedDefaults)
    allocate(self%lock          )
    !$omp critical (FoX_DOM_Access)
    self%document       => createDocument    (                                  &
         &                                    getImplementation()             , &
         &                                    qualifiedName      ="parameters", &
         &                                    docType            =null()        &
         &                                 )
    self%rootNode       => getDocumentElement(self%document)
    call setLiveNodeLists(self%document,.false.)
    !$omp end critical (FoX_DOM_Access)
    self%parameters     => null              (             )
    self%warnedDefaults =  integerHash       (             )
    self%lock           =  ompLock           (             )
    self%isNull         = .true.
   return
  end function inputParametersConstructorNull

  function inputParametersConstructorVarStr(xmlString,allowedParameterNames,outputParametersGroup,noOutput,changeFiles) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from an XML file
    specified as a variable length string.
    !!}
    use :: FoX_dom           , only : node
    use :: ISO_Varying_String, only : char                             , extract    , operator(==)
    use :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name
    use :: FoX_dom           , only : node                             , parseString
    use :: ISO_Varying_String, only : char                             , extract    , operator(==)
    implicit none
    type     (inputParameters)                                           :: self
    type     (varying_string    )              , intent(in   )           :: xmlString
    type     (varying_string    ), dimension(:), intent(in   ), optional :: allowedParameterNames, changeFiles
    type     (hdf5Object        ), target      , intent(in   ), optional :: outputParametersGroup
    logical                                    , intent(in   ), optional :: noOutput
    type     (node              ), pointer                               :: doc                  , parameterNode

    ! Check if we have been passed XML or a file name.
    if (extract(xmlString,1,1) == "<") then
       ! Parse the string.
       !$omp critical (FoX_DOM_Access)
       doc           => parseString(char(xmlString))
       parameterNode => XML_Get_First_Element_By_Tag_Name(              &
            &                                             doc         , &
            &                                             'parameters'  &
            &                                            ) 
       !$omp end critical (FoX_DOM_Access)
       self=inputParameters(                       &
            &               parameterNode        , &
            &               allowedParameterNames, &
            &               outputParametersGroup, &
            &               noOutput               &
            &              )
    else
       self=inputParameters(                       &
            &               char(xmlString)      , &
            &               allowedParameterNames, &
            &               outputParametersGroup, &
            &               noOutput             , &
            &               changeFiles            &
            &              )
    end if
    return
  end function inputParametersConstructorVarStr

  function inputParametersConstructorFileChar(fileName,allowedParameterNames,outputParametersGroup,noOutput,changeFiles) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from an XML file
    specified as a character variable.
    !!}
    use :: Display           , only : displayGreen                     , displayReset
    use :: File_Utilities    , only : File_Exists
    use :: FoX_dom           , only : node                             , getAttribute, setAttribute          , getParentNode, &
         &                            removeChild                      , getNodeName , hasAttribute          , appendChild  , &
         &                            importNode                       , insertBefore, getNextSibling        , destroy      , &
         &                            getExceptionCode                 , inException , DOMException          , cloneNode    , &
         &                            getNodeType                      , ELEMENT_NODE
    use :: Error             , only : Error_Report
    use :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Parse   , XML_Get_Child_Elements, xmlNodeList  , &
         &                            XML_Path_Exists
    use :: ISO_Varying_String, only : trim                             , char        , assignment(=)         , operator(//) , &
         &                            operator(==)                     , index       , extract
    implicit none
    type     (inputParameters)                                        :: self
    character(len=*          )              , intent(in   )           :: fileName
    type     (varying_string ), dimension(:), intent(in   ), optional :: allowedParameterNames, changeFiles
    type     (hdf5Object     ), target      , intent(in   ), optional :: outputParametersGroup
    logical                                 , intent(in   ), optional :: noOutput
    type     (xmlNodeList    ), dimension(:), allocatable             :: childNodes           , newNodes
    type     (node           ), pointer                               :: doc                  , parameterNode    , &
         &                                                               changesDoc           , childNode        , &
         &                                                               changeNode           , changeNodeParent , &
         &                                                               changesNode          , importedNode     , &
         &                                                               newNode              , changeSiblingNode, &
         &                                                               targetNode           , clonedNode
    integer                                                           :: errorStatus          , i                , &
         &                                                               j                    , k
    logical                                                           :: append
    type     (varying_string )                                        :: changePath           , targetPath       , &
         &                                                               valueUpdated
    character(len=32         )                                        :: changeType

    ! Check that the file exists.
    if (.not.File_Exists(fileName)) call Error_Report("parameter file '"//trim(fileName)//"' does not exist"//{introspection:location})
    ! Open and parse the data file.
    !$omp critical (FoX_DOM_Access)
    doc => XML_Parse(fileName,iostat=errorStatus)
    if (errorStatus /= 0) then
       if (File_Exists(fileName)) then
          call Error_Report(                                                                                                                                                                &
               &            'Unable to parse parameter file: "'//trim(fileName)//'"'//char(10)//                                                                                            &
               &            displayGreen()//"HELP:"//displayReset()//" check that the XML in this file is valid (e.g. `xmllint --noout "//trim(fileName)//"` will display any XML errors"// &
               &            {introspection:location}                                                                                                                                        &
               &           )
       else
          call Error_Report(                                                                                                                                                                &
               &            'Unable to find parameter file: "' //trim(fileName)//'"'//                                                                                                      &
               &            {introspection:location}                                                                                                                                        &
               &           )
       end if
    end if
    parameterNode => XML_Get_First_Element_By_Tag_Name(              &
         &                                             doc         , &
         &                                             'parameters'  &
         &                                            )
    !$omp end critical (FoX_DOM_Access)
    ! Apply changes from any changes files.
    if (present(changeFiles).and.size(changeFiles) > 0) then
       do i=1,size(changeFiles)
          !$omp critical (FoX_DOM_Access)
          changesDoc =>  XML_Parse(char(changeFiles(i)),iostat=errorStatus)
          !$omp end critical (FoX_DOM_Access)
          if (errorStatus /= 0) then
             if (File_Exists(changeFiles(i))) then
                call Error_Report(                                                                                                                                                                      &
                     &            'Unable to parse parameter changes file: "'//trim(changeFiles(i))//'"'//char(10)//                                                                                    &
                     &            displayGreen()//"HELP:"//displayReset()//" check that the XML in this file is valid (e.g. `xmllint --noout "//trim(changeFiles(i))//"` will display any XML errors"// &
                     &            {introspection:location}                                                                                                                                              &
                     &           )
             else
                call Error_Report(                                                                                                                                                                      &
                     &            'Unable to find parameter file: "' //trim(changeFiles(i))//'"'//                                                                                                      &
                     &            {introspection:location}                                                                                                                                              &
                     &           )
             end if
          end if
          ! Get the root node.
          !$omp critical (FoX_DOM_Access)
          changesNode => XML_Get_First_Element_By_Tag_Name(            &
               &                                           changesDoc, &
               &                                           'changes'   &
               &                                          )
          ! Process all `change` nodes.
          call XML_Get_Child_Elements(changesNode,childNodes)
          do j=0,size(childNodes)-1
             childNode => childNodes(j)%element
             ! Skip non-`change` nodes.
             if (.not.getNodeName(childNode) == "change") cycle
             ! Validate that the node has the required attributes.
             if (.not.hasAttribute(childNode,"type")) call Error_Report('`change` element must have the `type` attribute'//{introspection:location})
             if (.not.hasAttribute(childNode,"path")) call Error_Report('`change` element must have the `path` attribute'//{introspection:location})
             ! Extract the change type.
             changeType=getAttribute(childNode,"type")
             ! Find the node in the parameters document to be changed.
             changePath=getAttribute(childNode,"path")
             if (XML_Path_Exists(parameterNode,char(changePath))) then
                if (changePath == "") then
                   changeNode =>                                   parameterNode
                else
                   changeNode => XML_Get_First_Element_By_Tag_Name(parameterNode,char(changePath),directChildrenOnly=.true.)
                end if
                if (trim(changeType) == "replaceOrAppend") changeType="replace"
             else if (trim(changeType) == "replaceOrAppend") then
                ! replaceOrAppend change, but path does not exist - therefore we switch to appending to the parent.
                changePath =  extract(changePath,1,index(changePath,"/",back=.true.)-1)
                changeNode => XML_Get_First_Element_By_Tag_Name(parameterNode,char(changePath),directChildrenOnly=.true.)
                changeType =  "append"
             else
                call Error_Report("path '"//trim(changePath)//"' does not exist"//{introspection:location})
             end if
             ! Process each type of change.
             select case (trim(changeType))
             case ("remove")
                ! Remove the identified node.
                changeNodeParent => getParentNode(                 changeNode)
                changeNode       => removeChild  (changeNodeParent,changeNode)
             case ("update")
                ! Update the value of the identified node.
                if (.not.hasAttribute(changeNode,"value" )) call Error_Report('can not update the `value` in a parameter that has no `value`'//{introspection:location})
                if (.not.hasAttribute(childNode ,"value" )) call Error_Report('`change` element must have the `value` attribute'             //{introspection:location})
                if (     hasAttribute(childNode ,"append")) then
                   append=getAttribute(childNode,"append") == "true"
                else
                   append=.false.
                end if
                if (append) then
                   valueUpdated=getAttribute(changeNode,"value")
                else
                   valueUpdated=""
                end if
                valueUpdated=valueUpdated//getAttribute(childNode,"value")
                call setAttribute(changeNode,"value",char(valueUpdated))
             case ("append")
                ! Append new parameters.
                call XML_Get_Child_Elements(childNode,newNodes)
                do k=0,size(newNodes)-1
                   newNode      => newNodes(k)%element
                   importedNode => importNode (doc       ,newNode     ,deep=.true.)
                   importedNode => appendChild(changeNode,importedNode            )
                end do
             case ("insertBefore","insertAfter","replace")
                ! Insert new parameters before/after the identified node, or replace that node.
                changeNodeParent => getParentNode(changeNode)
                call XML_Get_Child_Elements(childNode,newNodes)
                do k=0,size(newNodes)-1
                   newNode      => newNodes(k)%element
                   importedNode => importNode(doc,newNode,deep=.true.)
                   if (trim(changeType) == "insertAfter") then
                      changeSiblingNode => getNextSibling(changeNode)
                      if (associated(changeSiblingNode)) then
                         importedNode => insertBefore(changeNodeParent,importedNode,changeSiblingNode)
                      else
                         importedNode => appendChild (changeNodeParent,importedNode                  )
                      end if
                   else
                      importedNode    => insertBefore(changeNodeParent,importedNode,changeNode       )
                   end if
                end do
                if (trim(changeType) == "replace") changeNode => removeChild(changeNodeParent,changeNode)
             case ("replaceWith")
                ! Replace a parameter with another parameter.
                if (.not.hasAttribute(childNode,"target")) call Error_Report('`change` element with `type="replaceWith"` must have the `target` attribute'//{introspection:location})
                targetPath=getAttribute(childNode,"target")
                if (.not.XML_Path_Exists(parameterNode,char(targetPath))) call Error_Report("target path '"//trim(targetPath)//"' does not exist"//{introspection:location})
                targetNode       => XML_Get_First_Element_By_Tag_Name(parameterNode,char(targetPath),directChildrenOnly=.true.)
                clonedNode       => cloneNode    (targetNode                            ,deep=.true.)
                changeNodeParent => getParentNode(                            changeNode            )
                clonedNode       => insertBefore (changeNodeParent,clonedNode,changeNode            )
                changeNode       => removeChild  (changeNodeParent           ,changeNode            )
             case ("encapsulate")
                ! Encapsulate the identified node within the provided content.
                !! First insert all new content before the change node.
                changeNodeParent => getParentNode(changeNode)
                call XML_Get_Child_Elements(childNode,newNodes)
                clonedNode => null()
                do k=0,size(newNodes)-1
                   newNode      => newNodes(k)%element
                   importedNode => importNode(doc,newNode,deep=.true.)
                   importedNode => insertBefore(changeNodeParent,importedNode,changeNode)
                   ! Keep a pointer to the first node of the new content - this is where the change node will be encapsulated.
                   if (.not.associated(clonedNode) .and. getNodeType(importedNode) == ELEMENT_NODE) clonedNode => importedNode
                end do
                ! Remove the change node from the document.
                changeNode => removeChild(changeNodeParent,changeNode)
                ! Reinsert the change node into the encapsulating node.
                changeNode => appendChild(clonedNode      ,changeNode)
             case default
                call Error_Report("unknown change type `"//trim(changeType)//"`"//{introspection:location})
             end select
          end do
          call destroy(changesDoc)
          !$omp end critical (FoX_DOM_Access)
       end do
    end if
    ! Construct the parameter tree from the XML document.
    self=inputParameters(                                &
         &                        parameterNode        , &
         &                        allowedParameterNames, &
         &                        outputParametersGroup, &
         &                        noOutput             , &
         &               fileName=fileName               &
         &              )
    return
  end function inputParametersConstructorFileChar

  function inputParametersConstructorCopy(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from an existing parameters object.
    !!}
    implicit none
    type(inputParameters)                :: self
    type(inputParameters), intent(in   ) :: parameters

    self            =  inputParameters(parameters%rootNode  ,noOutput=.true.,noBuild=.true.)
    self%parameters =>                 parameters%parameters
    self%parent     =>                 parameters%parent
    self%original   =>                 parameters%original       
    if (allocated(parameters%warnedDefaults)) then
       if (allocated(self%warnedDefaults)) deallocate(self%warnedDefaults)
       allocate(self%warnedDefaults)
       self%warnedDefaults=parameters%warnedDefaults
    end if
    if (associated(parameters%lock)) then
       deallocate(self%lock)
       allocate  (self%lock)
       self%lock=ompLock()
    end if
    return
  end function inputParametersConstructorCopy

  function inputParametersConstructorNode(parametersNode,allowedParameterNames,outputParametersGroup,noOutput,noBuild,fileName) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from a FoX node.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_null_char
    use            :: Display           , only : displayGreen                     , displayMessage  , displayMagenta  , displayReset  , &
         &                                       verbosityLevelSilent
    use            :: File_Utilities    , only : File_Name_Temporary              , File_Remove
    use            :: FoX_dom           , only : getOwnerDocument                 , node            , setLiveNodeLists  , getTextContent, &
         &                                       hasAttribute                     , getAttributeNode, extractDataContent
    use            :: Error             , only : Error_Report
#ifdef GIT2AVAIL
    use            :: Input_Paths       , only : pathTypeExec                     , inputPath
    use            :: Output_Versioning , only : Version
#else
    use            :: Error             , only : Warn
#endif
    use            :: ISO_Varying_String, only : assignment(=)                    , char           , operator(//)      , operator(/=)   , &
         &                                       var_str
    use            :: String_Handling   , only : String_Strip,String_C_To_Fortran , operator(//)
    use            :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Path_Exists
    use            :: Display           , only : displayMessage
    use            :: IO_HDF5           , only : ioHDF5AccessInitialize
    use            :: HDF5_Access       , only : hdf5Access
    implicit none
    type     (inputParameters)                                        :: self
    type     (node           ), pointer     , intent(in   )           :: parametersNode
    type     (varying_string ), dimension(:), intent(in   ), optional :: allowedParameterNames
    character(len=*          )              , intent(in   ), optional :: fileName
    type     (hdf5Object     ), target      , intent(in   ), optional :: outputParametersGroup
    logical                                 , intent(in   ), optional :: noOutput                   , noBuild
    type     (varying_string )                                        :: message
    type     (node           ), pointer                               :: lastModifiedNode           , revisionNode       , &
         &                                                               strictNode
    logical                                                           :: hasRevision                , hasStrict
    character(len=41         )                                        :: commitHashParameters
#ifdef GIT2AVAIL
    integer  (c_int          ), dimension(:), allocatable             :: isAncestorOfParameters
    integer  (c_int          )                                        :: isAncestorOfSelf
    character(len=41         )                                        :: commitHashSelf
    character(len=42         )                                        :: commitHashSelf_
    integer                                                           :: i
#endif
    type     (varying_string ), dimension(:), allocatable  , save     :: allowedParameterNamesGlobal
    !$omp threadprivate(allowedParameterNamesGlobal)
    !![
    <optionalArgument name="noOutput" defaultsTo=".false." />
    <optionalArgument name="noBuild"  defaultsTo=".false." />
    !!]
#include "os.inc"

    allocate(self%warnedDefaults)
    allocate(self%lock          )
    self%isNull         =  .false.
    self%rootNode       =>             parametersNode
    self%parent         => null       (              )
    self%warnedDefaults =  integerHash(              )
    self%lock           =  ompLock    (              )
    !$omp critical (FoX_DOM_Access)
    self%document         => getOwnerDocument(parametersNode)
    call setLiveNodeLists(self%document,.false.)
    !$omp end critical (FoX_DOM_Access)
    if (.not.noBuild_) then
       allocate(self%parameters)
       self%parameters%content    => null()
       self%parameters%parent     => null()
       self%parameters%firstChild => null()
       self%parameters%sibling    => null()
       self%parameters%referenced => null()
       !$omp critical (FoX_DOM_Access)
       call self%buildTree           (self%parameters,parametersNode)
       call self%resolveReferences   (                              )
       !$omp end critical (FoX_DOM_Access)
       call self%evaluateConditionals(                              )
    end if
    ! Set a pointer to HDF5 object to which to write parameters.
    if (present(outputParametersGroup)) then
       !$ call hdf5Access%  set()
       self%outputParameters         =outputParametersGroup%openGroup('Parameters',attributesCompactMaxiumum=0)
       self%outputParametersCopied   =.false.
       self%outputParametersTemporary=.false.
       !$ call hdf5Access%unset()
    else if (.not.noOutput_) then
       ! The HDF5 access lock may not yet have been initialized. Ensure it is before using it.
       call ioHDF5AccessInitialize()
       !$ call hdf5Access%  set()
       call self%outputParametersContainer%openFile(                                      &
            &                                       char(                                 &
            &                                            File_Name_Temporary(             &
            &                                                                'glcTmpPar', &
#ifdef __APPLE__
            &                                                                '/tmp'       &
#else
            &                                                                '/dev/shm'   &
#endif
            &                                                               )             &
            &                                            )                                &
            &                                      )
       self%outputParameters         =self%outputParametersContainer%openGroup('Parameters',attributesCompactMaxiumum=0)
       self%outputParametersCopied   =.false.
       self%outputParametersTemporary=.true.
       !$ call hdf5Access%unset()
       call File_Remove(self%outputParametersContainer%name())
    end if
    ! Get allowed parameter names.
    if (.not.allocated(allowedParameterNamesGlobal)) &
         & call knownParameterNames(allowedParameterNamesGlobal)
    ! Check for migration information.
    if (present(fileName)) then
       if (XML_Path_Exists(self%rootNode,"lastModified")) then
          ! Look for a "lastModified" element in the parameter file.
          !$omp critical (FoX_DOM_Access)
          lastModifiedNode => XML_Get_First_Element_By_Tag_Name(self%rootNode        ,'lastModified')
          hasRevision      =  hasAttribute                     (     lastModifiedNode,'revision'    )
          hasStrict        =  hasAttribute                     (     lastModifiedNode,'strict'      )
          if (hasRevision) then
             revisionNode         => getAttributeNode(lastModifiedNode,'revision')
             commitHashParameters =  getTextContent  (revisionNode               )//c_null_char
          end if
          if (hasStrict  ) then
             strictNode           => getAttributeNode(lastModifiedNode,'strict'  )
             call extractDataContent(strictNode,self%strict)
          end if
          !$omp end critical (FoX_DOM_Access)
#ifdef GIT2AVAIL
          if (hasRevision) then
             ! A revision was available in the parameter file.
             !! Build an array of known migration commit hashes.
             !![
	     <parameterMigration/>
             !!]
             !! Extract the commit hash at which Galacticus was built.
             call Version(commitHashSelf_)
             commitHashSelf=trim(commitHashSelf_)//c_null_char
             !! Iterate over known migration commit hashes and check if they are ancestors.
             allocate(isAncestorOfParameters(size(commitHash)))
             do i=1,size(commitHash)
                isAncestorOfParameters(i)=gitDescendantOf(char(inputPath(pathTypeExec))//c_null_char,commitHashParameters,commitHash(i))
             end do
             if (any(isAncestorOfParameters /= 0_c_int .and. isAncestorOfParameters /= 1_c_int)) then
                message=var_str("parameter file revision check failed (#1; error code; ")//maxval(isAncestorOfParameters)//")"
                if (self%strict) then
                   call Error_Report(message//{introspection:location})
                else
                   call displayMessage(displayMagenta()//"WARNING: "//displayReset()//message)
                end if
             else if (any(isAncestorOfParameters == 0)) then
                ! Parameter file is missing migrations - issue a warning.
                message="parameter file may be missing important parameter updates - consider updating by running:"//char(10)//char(10)//"              ./scripts/aux/parametersMigrate.pl "//trim(fileName)//" newParameterFile.xml"
                if (self%strict) then
                   call Error_Report(message//{introspection:location})
                else
                   call displayMessage(displayMagenta()//"WARNING: "//displayReset()//message//char(10),verbosityLevelSilent)
                end if
             end if
             isAncestorOfSelf=gitDescendantOf(char(inputPath(pathTypeExec))//c_null_char,commitHashSelf,commitHashParameters)
             if (isAncestorOfSelf /= 0_c_int .and. isAncestorOfSelf /= 1_c_int) then
                message=var_str("parameter file revision check failed (#2; error code: ")//isAncestorOfSelf//")"
                if (self%strict) then
                   call Error_Report(message//{introspection:location})
                else
                   call displayMessage(displayMagenta()//"WARNING: "//displayReset()//message)
                end if
             else if (isAncestorOfSelf == 0_c_int) then
                ! Parameters are more recent than the executable - issue a warning.
                message="parameter file revision is newer than this executable - consider updating your copy of Galacticus"
                if (self%strict) then
                   call Error_Report(message//{introspection:location})
                else
                   call displayMessage(displayMagenta()//"WARNING: "//displayReset()//message,verbosityLevelSilent)
                end if
             end if
          end if
#else
          message="can not check if parameter file is up to date (`libgit` is not available)"
          if (self%strict) then
             call Error_Report(message//{introspection:location})
          else
             call Warn(displayMagenta()//"WARNING: "//displayReset()//message)
          end if
#endif
       end if
    end if
    ! Check parameters.
    call self%checkParameters(allowedParameterNamesGlobal=allowedParameterNamesGlobal,allowedParameterNames=allowedParameterNames)
    return
  end function inputParametersConstructorNode

  recursive subroutine inputParametersBuildTree(self,parentParameter,parametersNode)
    !!{
    Build a tree representation of the input parameter file.
    !!}
    use :: FoX_DOM, only : ELEMENT_NODE          , getNodeType, node
    use :: IO_XML , only : XML_Get_Child_Elements, xmlNodeList
    implicit none
    class  (inputParameters), intent(inout)              :: self
    type   (inputParameter ), intent(inout), pointer     :: parentParameter
    type   (node           ), intent(in   ), pointer     :: parametersNode
    type   (xmlNodeList    ), dimension(:) , allocatable :: childNodes
    type   (node           )               , pointer     :: childNode
    type   (inputParameter )               , pointer     :: currentParameter
    integer                                              :: i

    call XML_Get_Child_Elements(parametersNode,childNodes)
    currentParameter => null()
    do i=0,size(childNodes)-1
       childNode => childNodes(i)%element
       if (getNodeType(childNode) == ELEMENT_NODE) then
          if (associated(currentParameter)) then
             allocate(currentParameter%sibling)
             currentParameter => currentParameter%sibling
          else
             allocate(parentParameter%firstChild)
             currentParameter => parentParameter%firstChild
          end if
          currentParameter%content    => childNode
          currentParameter%parent     => parentParameter
          currentParameter%firstChild => null()
          currentParameter%sibling    => null()
          currentParameter%referenced => null()
          call self%buildTree(currentParameter,childNode)
       end if
    end do
    return
  end subroutine inputParametersBuildTree

  subroutine inputParametersResolveReferences(self)
    !!{
    Resolve references in a parameter tree.
    !!}
    use :: FoX_dom           , only : DOMException , ELEMENT_NODE  , getAttributeNode, getNodeName, &
          &                           getNodeType  , getTextContent, hasAttribute    , inException
    use :: ISO_Varying_String, only : assignment(=), operator(==)  , char
    use :: Error             , only : Error_Report
    implicit none
    class  (inputParameters), intent(inout) :: self
    type   (inputParameter ), pointer       :: currentParameter, referencedParameter
    type   (node           ), pointer       :: identifierNode  , identifierReferenceNode
    type   (varying_string )                :: identifier      , identifierReference
    type   (DOMException   )                :: exception
    logical                                 :: found

    ! Begin walking the parameter tree.
    currentParameter => inputParametersWalkTree(self%parameters)
    do while (associated(currentParameter))
       ! Find parameters which reference another parameter.
       if     (                                                                &
            &   getNodeType (currentParameter%content        ) == ELEMENT_NODE &
            &  .and.                                                           &
            &   hasAttribute(currentParameter%content,'idRef')                 &
            & ) then
          ! Search for a parameter with the referenced ID and the same name.
          found               =  .false.
          referencedParameter => inputParametersWalkTree(self%parameters)
          do while (associated(referencedParameter))
             ! If found, set a pointer to this other parameter which will be later dereferenced for parameter extraction.
             if     (                                                                                         &
                  &   getNodeType (referencedParameter%content     ) == ELEMENT_NODE                          &
                  &  .and.                                                                                    &
                  &   hasAttribute(referencedParameter%content,'id')                                          &
                  &  .and.                                                                                    &
                  &   getNodeName (referencedParameter%content     ) == getNodeName(currentParameter%content) &
                  & ) then
                identifierNode          => getAttributeNode(referencedParameter    %content,   'id'     )
                identifierReferenceNode => getAttributeNode(currentParameter       %content,   'idRef'  )
                identifier              =  getTextContent  (identifierNode                 ,ex=exception)
                if (inException(exception)) call Error_Report('unable to parse identifier'//{introspection:location})
                identifierReference     =  getTextContent  (identifierReferenceNode        ,ex=exception)
                if (inException(exception)) call Error_Report('unable to parse identifier'//{introspection:location})
                if (identifier == identifierReference) then
                   currentParameter%referenced => referencedParameter
                   found=.true.
                   exit
                end if
             end if
             ! Walk to next node.
             referencedParameter => inputParametersWalkTree(referencedParameter)
          end do
          if (.not.found) then
             identifierReferenceNode => getAttributeNode(currentParameter       %content,   'idRef'  )
             identifierReference     =  getTextContent  (identifierReferenceNode        ,ex=exception)
             if (inException(exception)) call Error_Report('unable to parse identifier'//{introspection:location})
             call Error_Report('unable to find referenced parameter "'//char(identifierReference)//'"'//{introspection:location})
          end if
       end if
       ! Walk to next node.
       currentParameter => inputParametersWalkTree(currentParameter)
    end do
    return
  end subroutine inputParametersResolveReferences

  subroutine inputParametersEvaluateConditionals(self)
    !!{
    Evaluate conditionals in a parameter tree.
    !!}
    use :: FoX_dom           , only : DOMException      , ELEMENT_NODE      , getAttributeNode, inException, &
         &                            getNodeType       , getTextContent    , hasAttribute    , node       , &
         &                            getNodeName
    use :: ISO_Varying_String, only : assignment(=)     , operator(==)      , char            , extract    , &
         &                            index             , adjustl           , trim
    use :: Error             , only : Error_Report
    use :: String_Handling   , only : String_Count_Words, String_Split_Words
    use :: IO_XML            , only : XML_Path_Exists
    implicit none
    class    (inputParameters           ), intent(inout)              :: self
    type     (inputParameter            ), pointer                    :: currentParameter
    type     (node                      ), pointer                    :: conditionNode 
    type     (inputParameter            ), pointer                    :: parameterTest
    type     (node                      ), pointer                    :: node_
    character(len=parameterLengthMaximum), dimension(:) , allocatable :: parameterNames
    type     (varying_string            )                             :: condition
    type     (DOMException              )                             :: exception
    character(len=parameterLengthMaximum)                             :: parameterName   , parameterLeafName, &
         &                                                               valueTest       , valueParameter
    logical                                                           :: operatorEquals  , matches          , &
         &                                                               allEvaluated    , didEvaluate
    integer                                                           :: countNames      , i

    ! Iterate until all parameters are evaluated.
    didEvaluate=.true.
    do while (didEvaluate)
       allEvaluated=.true.
       didEvaluate =.false.
       ! Begin walking the parameter tree.
       currentParameter => inputParametersWalkTree(self%parameters)
       do while (associated(currentParameter))
          ! Find parameters with conditionals.
          if     (                                                                 &
               &   getNodeType (currentParameter%content         ) == ELEMENT_NODE &
               &  .and.                                                            &
               &   hasAttribute(currentParameter%content,'active')                 &
               & ) then
             ! Extract the condition.
             conditionNode => getAttributeNode(currentParameter%content,   'active' )
             condition     =  getTextContent  (conditionNode           ,ex=exception)
             if (inException(exception)) call Error_Report('unable to parse conditional'//{introspection:location})
             ! Parse the condition.
             !! Extract the name of the parameter being conditioned upon.
             if (extract(condition,1,1) == "[" .and. index(condition,"]") > 2) then
                parameterName=extract(condition,2,index(condition,"]")-1)
                condition    =adjustl(extract(condition,index(condition,"]")+1))
             else
                call Error_Report("unable to parse parameter name in conditional"//{introspection:location})
             end if
             !! Extract the operator (currently only `==` and `!=` are supported).
             if      (extract(condition,1,2) == "==") then
                operatorEquals=.true.
             else if (extract(condition,1,2) == "!=") then
                operatorEquals=.false.
             else
                operatorEquals=.true.
                call Error_Report("unable to parse operator in conditional"//{introspection:location})
             end if
             condition=adjustl(extract(condition,3))
             !! Extract the string being compared to.
             valueTest=trim(condition)
             ! Get the value of the parameter.
             countNames=String_Count_Words(parameterName,":")
             allocate(parameterNames(countNames))
             call String_Split_Words(parameterNames,parameterName,":")
             parameterLeafName=parameterNames(countNames)
             parameterTest => currentParameter
             if (trim(parameterNames(1)) /= "." .and. trim(parameterNames(1)) /= "..") then
                ! Path is absolute - move to the root parameter.
                do while (associated(parameterTest%parent))
                   parameterTest => parameterTest%parent
                end do
             end if
             do i=1,countNames
                if (trim(parameterNames(i)) == ".") then
                   ! Self - no need to move.
                else if (trim(parameterNames(i)) == "..") then
                   ! Move to the parent parameter.
                   if (.not.associated(parameterTest%parent)) call Error_Report('no parent parameter exists'//{introspection:location})
                   parameterTest => parameterTest%parent
                else
                   ! Move to the named parameter.
                   parameterTest => parameterTest%firstChild
                   do while (associated(parameterTest))
                      node_ => parameterTest%content
                      !$omp critical (FoX_DOM_Access)
                      matches= parameterTest%active                          &
                           &  .and.                                          &
                           &   getNodeType(node_) == ELEMENT_NODE            &
                           &  .and.                                          &
                           &   (                                             &
                           &        hasAttribute(node_,'value')              &
                           &    .or.                                         &
                           &     XML_Path_Exists(node_,"value")              &
                           &    .or.                                         &
                           &        hasAttribute(node_,"idRef")              &
                           &   )                                             &
                           &  .and.                                          &
                           &   trim(parameterNames(i)) == getNodeName(node_)
                      !$omp end critical (FoX_DOM_Access)
                      if (matches) exit
                      parameterTest => parameterTest%sibling
                   end do
                   if (.not.associated(parameterTest)) call Error_Report('no child parameter exists'//{introspection:location})
                   if (associated(parameterTest%referenced)) parameterTest => parameterTest%referenced
                end if
             end do
             if (parameterTest%activeEvaluated) then
                valueParameter=parameterTest%get()
                ! Perform the test and set parameter active state.
                currentParameter%active         =(trim(valueParameter) == trim(valueTest)) .eqv. operatorEquals
                currentParameter%activeEvaluated=.true.
                didEvaluate                     =.true.
             else
                allEvaluated                    =.false.
             end if
             deallocate(parameterNames)
          else
             currentParameter%activeEvaluated=.true.
             didEvaluate                     =.true.
          end if
          ! Walk to next node.
          currentParameter => inputParametersWalkTree(currentParameter)
       end do
       if (.not.didEvaluate) call Error_Report('failed to evaluate parameter active statuses'//{introspection:location})
       if (allEvaluated) exit
    end do
    return
  end subroutine inputParametersEvaluateConditionals

  function inputParametersWalkTree(currentNode) result(nextNode)
    !!{
    Perform a depth-first walk of a parameter tree.
    !!}
    implicit none
    type(inputParameter), pointer                :: nextNode
    type(inputParameter), pointer, intent(in   ) :: currentNode

    nextNode => currentNode
    if (.not.associated(nextNode%parent)) then
       do while (associated(nextNode%firstChild))
          nextNode => nextNode%firstChild
       end do
       if (associated(nextNode,currentNode)) nullify(nextNode)
    else
       if (associated(nextNode%sibling)) then
          nextNode => nextNode%sibling
          do while (associated(nextNode%firstChild))
             nextNode => nextNode%firstChild
          end do
       else
          nextNode => nextNode%parent
          if (.not.associated(nextNode%parent)) nextNode => null()
       end if
    end if
    return
  end function inputParametersWalkTree

  subroutine inputParametersDestroy(self)
    !!{
    Destructor for the {\normalfont \ttfamily inputParameters} class.
    !!}
    use :: FoX_DOM, only : destroy
    implicit none    
    class(inputParameters), intent(inout) :: self

    ! Destroy the parameters document. Note that we do not use a finalizer for input parameters. This could destroy part of a
    ! document which was still pointed to from elsewhere, leaving a dangling pointer. Instead, destruction only occurs when
    ! explicitly requested.
    !$omp critical (FoX_DOM_Access)
    call destroy(self%document)
    !$omp end critical (FoX_DOM_Access)
    nullify(self%document)
    if (associated(self%parameters)) then
       call self%parameters%destroy()
       deallocate(self%parameters)
    end if
    call inputParametersFinalize(self)
    return
  end subroutine inputParametersDestroy

  subroutine inputParametersFinalize(self)
    !!{
    Finalizer for the {\normalfont \ttfamily inputParameters} class.
    !!}
    use :: FoX_dom           , only : destroy
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : char
    implicit none
    type(inputParameters), intent(inout) :: self
    type(varying_string )                :: fileNameTemporary

    if (self%isNull) then
       !$omp critical (FoX_DOM_Access)
       call destroy(self%document)
       !$omp end critical (FoX_DOM_Access)
    end if
    nullify(self%document  )
    nullify(self%rootNode  )
    nullify(self%parameters)
    nullify(self%parent    )
    if (allocated (self%warnedDefaults)) deallocate(self%warnedDefaults)
    if (associated(self%lock          )) deallocate(self%lock          )
    !$ call hdf5Access%set()
    if (self%outputParameters%isOpen().and..not.self%outputParametersCopied) then
       if (self%outputParametersTemporary) then
          ! Close and remove the temporary parameters file.
          fileNameTemporary=self%outputParametersContainer%name()
          call self%outputParameters         %close  ()
          call self%outputParametersContainer%close  ()
          call self%outputParameters         %destroy()
          call self%outputParametersContainer%destroy()
       else
          ! Simply close our parameters group.
          call self%outputParameters%close  ()
          call self%outputParameters%destroy()
       end if
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine inputParametersFinalize

  recursive subroutine inputParameterDestroy(self)
    !!{
    Destructor for the {\normalfont \ttfamily inputParameter} class.
    !!}
    class(inputParameter), intent(inout) :: self
    type (inputParameter), pointer       :: child, childNext

    ! We do not destroy the XML node content here as it may be used elsewhere.
    ! Destroy all children - this will trigger recursive destruction of all grandchildren, etc.
    child => self%firstChild
    do while (associated(child))
       childNext => child%sibling
       call child%destroy()
       deallocate(child)
       child => childNext
    end do
    nullify(self%firstChild)
    call self%reset(children=.false.)
    return
  end subroutine inputParameterDestroy

  logical function inputParameterIsParameter(self)
    !!{
    Return true if this is a valid parameter.
    !!}
    use :: FoX_dom, only : ELEMENT_NODE   , getNodeType, hasAttribute
    use :: IO_XML , only : XML_Path_Exists
    implicit none
    class(inputParameter), intent(in   ) :: self

    !$omp critical (FoX_DOM_Access)
    if (associated(self%content) .and. getNodeType(self%content) == ELEMENT_NODE) then
       inputParameterIsParameter=    hasAttribute(self%content,'value') &
            &                    .or.                                   &
            &                     XML_Path_Exists(self%content,'value') &
            &                    .or.                                   &
            &                        hasAttribute(self%content,"idRef")
    else
       inputParameterIsParameter=.false.
    end if
    !$omp end critical (FoX_DOM_Access)
    return
  end function inputParameterIsParameter

  logical function inputParameterObjectCreated(self)
    !!{
    Return true if the specified instance of the object associated with this parameter has been created.
    !!}
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_Get_Level, OMP_In_Parallel
    implicit none
    class  (inputParameter), intent(in   ) :: self
    integer                                :: instance, level

    inputParameterObjectCreated=.false.
    !$omp critical (inputParameterObjects)
    if (allocated(self%objects)) then
       !$ if (OMP_In_Parallel()) then
       !$    level   =OMP_Get_Level     ()
       !$    instance=OMP_Get_Thread_Num()+1
       !$ else
             level   =                     0
             instance=                     0
       !$ end if
       if (level <= ubound(self%objects,dim=1))                                           &
            & inputParameterObjectCreated=associated(self%objects(level,instance)%object)
    end if
    !$omp end critical (inputParameterObjects)
    return
  end function inputParameterObjectCreated

  function inputParameterObjectGet(self)
    !!{
    Return a pointer to the object associated with this parameter.
    !!}
    use    :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_Get_Level, OMP_In_Parallel
    implicit none
    class  (functionClass ), pointer       :: inputParameterObjectGet
    class  (inputParameter), intent(in   ) :: self
    integer                                :: instance               , level

    !$omp critical (inputParameterObjects)
    if (allocated(self%objects)) then
       !$ if (OMP_In_Parallel()) then
       !$    level   =OMP_Get_Level     ()
       !$    instance=OMP_Get_Thread_Num()+1
       !$ else
             level   =                     0
             instance=                     0
       !$ end if
       if (level <= ubound(self%objects,dim=1)) then
          inputParameterObjectGet => self%objects(level,instance)%object
       else
          call Error_Report('object not allocated for this parameter'//{introspection:location})
       end if
    else
       call Error_Report('object not allocated for this parameter'//{introspection:location})
    end if
    !$omp end critical (inputParameterObjects)
    return
  end function inputParameterObjectGet

  subroutine inputParameterObjectSet(self,object)
    !!{
    Set a pointer to the object associated with this parameter.
    !!}
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_Get_Level, OMP_Get_Max_Threads, OMP_In_Parallel
    implicit none
    class  (inputParameter   ), intent(inout)               :: self
    class  (functionClass    ), intent(in   ) , target      :: object
    type   (genericObjectList), dimension(:,:), allocatable :: objects
    integer                                                 :: instance, level

    !$omp critical (inputParameterObjects)
    !$ if (OMP_In_Parallel()) then
    !$    level   =OMP_Get_Level     ()
    !$    instance=OMP_Get_Thread_Num()+1
    !$ else
          level   =                     0
          instance=                     0
    !$ end if
    if (.not.allocated(self%objects)) then
       allocate(self%objects(0,0:1))
       !$ deallocate(self%objects)
       !$ allocate(self%objects(0:max(1,level),0:OMP_Get_Max_Threads()))
    else
       if (level > ubound(self%objects,dim=1)) then
          call move_alloc(self%objects,objects)
          allocate(self%objects(0:level,0:1))
          !$ deallocate(self%objects)
          !$ allocate(self%objects(0:level,0:OMP_Get_Max_Threads()))
          self%objects(0:ubound(objects,dim=1),:)=objects
          deallocate(objects)
       end if
    end if
    if (.not.associated(self%objects(level,instance)%object)) then      
       self%objects(level,instance)%object => object
       !![
       <referenceCountIncrement owner="self%objects(level,instance)" object="object"/>
       !!]
    end if
    !$omp end critical (inputParameterObjects)
    return
  end subroutine inputParameterObjectSet

  recursive subroutine inputParameterReset(self,children,evaluations)
    !!{
    Reset objects associated with this parameter and any sub-parameters.
    !!}
    use :: iso_varying_string, only : char
    use    :: FoX_DOM, only : destroy, getNodeName
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_Get_Level, OMP_In_Parallel
    implicit none
    class  (inputParameter), intent(inout), target   :: self
    logical                , intent(in   ), optional :: children, evaluations
    type   (inputParameter), pointer                 :: child   , childNext
    integer                                          :: instance, level
    !![
    <optionalArgument name="children"    defaultsTo=".true." />
    <optionalArgument name="evaluations" defaultsTo=".false."/>
    !!]

    ! Destroy any objects associated with this parameter node.
    if (allocated(self%objects)) then
       !$omp critical (inputParameterObjects)
       !$ if (OMP_In_Parallel()) then
       !$    level   =OMP_Get_Level     ()
       !$    instance=OMP_Get_Thread_Num()+1
       !$ else
             level   =                     0
             instance=                     0
       !$ end if
       if (level <= ubound(self%objects,dim=1)) then
          if (associated(self%objects(level,instance)%object)) then
             !![
             <objectDestructor name="self%objects(level,instance)%object"/>
             !!]
          end if
       end if
       !$omp end critical (inputParameterObjects)
    end if
    ! For evaluated parameters, restore them to their unevaluated state.
    if (evaluations_) then
       if (self%evaluated) then
          self%evaluated=.false.
          call self%set(self%contentOriginal)
       end if
    end if
    ! Reset parameter active state.
    self%activeEvaluated=.false.
    self%active         =.true.
    ! Clean up children if requested.
    if (children_) then
       ! Remove if this parameter was created.
       if (self%created) self%removed=.true.
       ! Reset children if requested.
       child => self%firstChild
       do while (associated(child))
          childNext => child%sibling
          call child%reset(children,evaluations)
          child => childNext
       end do
    end if
    return
  end subroutine inputParameterReset

  subroutine inputParameterSetDouble(self,value)
    !!{
    Set the value of a parameter.
    !!}
    use :: FoX_DOM, only : setAttribute
    implicit none
    class           (inputParameter), intent(inout) :: self
    double precision                , intent(in   ) :: value
    character       (len=24        )                :: valueText

    write (valueText,'(e24.16)') value
    !$omp critical (FoX_DOM_Access)
    call setAttribute(self%content,"value",trim(valueText))
    !$omp end critical (FoX_DOM_Access)
    return
  end subroutine inputParameterSetDouble

  subroutine inputParameterSetVarStr(self,value)
    !!{
    Set the value of a parameter.
    !!}
    use :: FoX_DOM           , only : setAttribute
    use :: ISO_Varying_String, only : char
    implicit none
    class    (inputParameter), intent(inout) :: self
    type     (varying_string), intent(in   ) :: value

    !$omp critical (FoX_DOM_Access)
    call setAttribute(self%content,"value",char(value))
    !$omp end critical (FoX_DOM_Access)
    return
  end subroutine inputParameterSetVarStr

  function inputParameterGet(self)
    !!{
    Get the value of a parameter.
    !!}
    use :: FoX_dom           , only : DOMException , getAttributeNode, getNodeName   , hasAttribute, &
          &                           inException  , node            , getTextContent
    use :: ISO_Varying_String, only : assignment(=)
    use :: Error             , only : Error_Report
    implicit none
    type   (varying_string)                :: inputParameterGet
    class  (inputParameter), intent(inout) :: self
    type   (node          ), pointer       :: valueElement
    type   (DOMException  )                :: exception
    logical                                :: hasValueAttribute

    !$omp critical (FoX_DOM_Access)
    hasValueAttribute=hasAttribute(self%content,'value')
    if (hasValueAttribute) then
       valueElement      => getAttributeNode(self%content,"value"     )
       inputParameterGet =  getTextContent  (valueElement,ex=exception)
       if (inException(exception)) call Error_Report(                                 &
            &                                        'unable to parse parameter [' // &
            &                                         getNodeName   (self%content) // &
            &                                        ']'                           // &
            &                                         {introspection:location}        &
            &                                       )
    else
       call Error_Report('no parameter value present'//{introspection:location})
    end if
    !$omp end critical (FoX_DOM_Access)
    return
  end function inputParameterGet

  subroutine inputParametersCheckParameters(self,allowedParameterNamesGlobal,allowedParameterNames,allowedMultiParameterNames)
    use :: Error              , only : Error_Report
    use :: Display            , only : displayIndent              , displayMagenta  , displayMessage                 , displayReset        , &
          &                            displayUnindent            , displayVerbosity, enumerationVerbosityLevelEncode, verbosityLevelSilent
    use :: FoX_dom            , only : DOMException               , destroy         , extractDataContent             , getAttributeNode    , &
          &                            getNodeName                , hasAttribute    , inException                    , node                , &
          &                            getParentNode
    use :: ISO_Varying_String , only : assignment(=)              , char            , operator(//)                   , operator(==)
    use :: Regular_Expressions, only : regEx
    use :: String_Handling    , only : String_Levenshtein_Distance
    implicit none
    class    (inputParameters                         )              , intent(inout)                    :: self
    type     (varying_string                          ), dimension(:), intent(in   ), optional, target  :: allowedParameterNamesGlobal, allowedParameterNames, &
         &                                                                                                 allowedMultiParameterNames
    type     (node                                    )                                       , pointer :: node_                      , ignoreWarningsNode   , &
         &                                                                                                 node__
    type     (inputParameter                          )                                       , pointer :: currentParameter
    type     (varying_string                          ), dimension(:)                         , pointer :: allowedParameterNames_
    type     (regEx                                   ), save                                           :: regEx_
    !$omp threadprivate(regEx_)
    logical                                                                                             :: warningsFound              , parameterMatched     , &
         &                                                                                                 verbose                    , ignoreWarnings       , &
         &                                                                                                 isException                , hasAttribute_        , &
         &                                                                                                 haveAllowedNames
    type     (enumerationInputParameterErrorStatusType)                                                 :: errorStatus
    integer                                                                                             :: allowedParametersCount     , status               , &
         &                                                                                                 distance                   , distanceMinimum      , &
         &                                                                                                 i                          , j
    character(len=1024                                )                                                 :: parameterValue
    character(len=1024                                )                                                 :: unknownName                , allowedParameterName , &
         &                                                                                                 parameterNameGuess         , unknownNamePath
    type     (varying_string                          )                                                 :: message                    , verbosityLevel
    type     (integerHash                             )                                                 :: parameterNamesSeen
    type     (DOMException                            )                                                 :: exception
    
    ! Determine whether we should be verbose.
    verbose=displayVerbosity() > verbosityLevelSilent
    if (self%isPresent('verbosityLevel')) then
       call self%value('verbosityLevel',verbosityLevel)
       verbose=enumerationVerbosityLevelEncode(char(verbosityLevel),includesPrefix=.false.) > verbosityLevelSilent
    end if
    ! Validate parameters.
    warningsFound     =.false.
    parameterNamesSeen=integerHash()
    if (associated(self%parameters)) then
       currentParameter => self%parameters%firstChild
       do while (associated(currentParameter))
          if (currentParameter%isParameter() .and. currentParameter%active) then
             node_ => currentParameter%content
             ! Attempt to read the parameter value.
             call self%value(currentParameter,parameterValue,errorStatus,writeOutput=.false.,evaluate=.false.)
             ! Determine if warnings should be ignored for this parameter.
             ignoreWarnings=.false.
             !$omp critical (FoX_DOM_Access)
             hasAttribute_=hasAttribute(node_,'ignoreWarnings')
             !$omp end critical (FoX_DOM_Access)
             if (hasAttribute_) then
                !$omp critical (FoX_DOM_Access)
                ignoreWarningsNode => getAttributeNode(node_,'ignoreWarnings')
                call extractDataContent(ignoreWarningsNode,ignoreWarnings,iostat=status,ex=exception)
                isException=inException(exception)
                unknownName=getNodeName(node_)
                !$omp end critical (FoX_DOM_Access)
                if (isException .or. status /= 0) &
                     & call Error_Report("unable to parse attribute 'ignoreWarnings' in parameter ["//trim(unknownName)//"]"//{introspection:location})
             end if
             ! Check for a match with allowed parameter names.
             haveAllowedNames=present(allowedParameterNamesGlobal).or.present(allowedParameterNames)
             parameterMatched=.not.haveAllowedNames
             allowedParametersCount=0
             do i=1,2
                select case (i)
                case (1)
                   if (.not.present(allowedParameterNamesGlobal)) cycle
                   allowedParameterNames_ => allowedParameterNamesGlobal
                case (2)
                   if (.not.present(allowedParameterNames)) cycle
                   allowedParameterNames_ => allowedParameterNames
                end select
                do j=1,size(allowedParameterNames_)
                   allowedParameterName=allowedParameterNames_(j)
                   if (allowedParameterName(1:6) == "regEx:") then
                      regEx_=regEx(allowedParameterName(7:len_trim(allowedParameterName)))
                      !$omp critical (FoX_DOM_Access)
                      parameterMatched=regEx_%matches(getNodeName(node_))
                      !$omp end critical (FoX_DOM_Access)
                      call regEx_%destroy()
                   else
                      !$omp critical (FoX_DOM_Access)
                      parameterMatched=(getNodeName(node_) == trim(allowedParameterName))
                      !$omp end critical (FoX_DOM_Access)
                   end if
                   if (parameterMatched) exit
                end do
             end do
             ! Report on warnings.
             if     (                                                   &
                  &   (                                                 &
                  &     errorStatus /= inputParameterErrorStatusSuccess &
                  &    .or.                                             &
                  &     .not.parameterMatched                           &
                  &   )                                                 &
                  &  .and.                                              &
                  &   .not.ignoreWarnings                               &
                  &  .and.                                              &
                  &   .not.warningsFound                                &
                  & ) then
                if (verbose) call displayIndent(displayMagenta()//'WARNING:'//displayReset()//' problems found with input parameters:')
                warningsFound=.true.
             end if
             if (errorStatus /= inputParameterErrorStatusSuccess .and. .not.ignoreWarnings .and. verbose) then
                !$omp critical (FoX_DOM_Access)
                select case (errorStatus%ID)
                case (inputParameterErrorStatusEmptyValue    %ID)
                   message='empty value for parameter ['    //getNodeName(node_)//']'
                case (inputParameterErrorStatusAmbiguousValue%ID)
                   message='ambiguous value for parameter ['//getNodeName(node_)//']'
                end select
                !$omp end critical (FoX_DOM_Access)
                call displayMessage(message)
             end if
             if (haveAllowedNames .and. .not.parameterMatched .and. .not.ignoreWarnings .and. verbose) then
                node__          => getParentNode(node_)
                unknownNamePath =  ""
                !$omp critical (FoX_DOM_Access)
                unknownName=getNodeName(node_)
                unknownNamePath=""
                do while (associated(node__))
                   if (getNodeName(node__) /= "#document")                                     &
                        & unknownNamePath =  getNodeName  (node__)//"/"//trim(unknownNamePath)
                   node__                 => getParentNode(node__)
                end do
                !$omp end critical (FoX_DOM_Access)
                distanceMinimum=-1
                do i=1,2
                   select case (i)
                   case (1)
                      if (.not.present(allowedParameterNamesGlobal)) cycle
                      allowedParameterNames_ => allowedParameterNamesGlobal
                   case (2)
                      if (.not.present(allowedParameterNames)) cycle
                      allowedParameterNames_ => allowedParameterNames
                   end select
                   do j=1,size(allowedParameterNames_)
                      allowedParameterName=allowedParameterNames_(j)
                      if (allowedParameterName(1:6) == "regEx:") cycle
                      distance=String_Levenshtein_Distance(trim(unknownName),trim(allowedParameterName))
                      if (distance < distanceMinimum .or. 0 > distanceMinimum) then
                         distanceMinimum   =distance
                         parameterNameGuess=allowedParameterName
                      end if
                   end do
                end do
                if (verbose) then
                   message='unrecognized parameter ['//trim(unknownName)//' in '//trim(unknownNamePath)//']'
                   if (distanceMinimum >= 0) message=message//' (did you mean ['//trim(parameterNameGuess)//']?)'
                   call displayMessage(message)
                end if
             end if
             ! Check for duplicated parameters.
             !$omp critical (FoX_DOM_Access)
             unknownName=getNodeName(node_)
             !$omp end critical (FoX_DOM_Access)
             if (parameterNamesSeen%exists(unknownName)) then
                parameterMatched=.false.
                if (present(allowedMultiParameterNames)) &
                     & parameterMatched=any(unknownName == allowedMultiParameterNames)
                if (.not.parameterMatched .and. .not.ignoreWarnings) then
                   if (.not.warningsFound.and.verbose) call displayIndent(displayMagenta()//'WARNING:'//displayReset()//' problems found with input parameters:')
                   warningsFound=.true.
                   if (verbose) then
                      message='multiple copies of parameter ['//getNodeName(node_)//'] present - only the first will be utilized'
                      call displayMessage(message)
                   end if
                end if
             else
                call parameterNamesSeen%set(unknownName,1)
             end if
          end if
          currentParameter => currentParameter%sibling
       end do
    end if
    if (warningsFound .and. verbose) call displayUnindent('')
    if (warningsFound .and. self%strict) call Error_Report('warnings found and strict compliance requested'//{introspection:location})
    return
  end subroutine inputParametersCheckParameters

  function inputParametersParametersGroup(self) result(parametersGroup)
    !!{
    Return the HDF5 group to which this parameters content will be written.
    !!}
    implicit none
    type(hdf5Object      )                :: parametersGroup
    class(inputParameters), intent(inout) :: self

    if (self%outputParameters%isOpen()) call self%outputParameters%deepCopy(parametersGroup)
    return
  end function inputParametersParametersGroup

  subroutine inputParametersParametersGroupOpen(self,outputGroup)
    !!{
    Open an output group for parameters in the given HDF5 object.
    !!}
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : char
    implicit none
    class(inputParameters), intent(inout) :: self
    type (hdf5Object     ), intent(inout) :: outputGroup
    type (varying_string )                :: fileNameTemporary

    !$ call hdf5Access%set()
    if (self%outputParameters%isOpen().and.self%outputParametersTemporary) then
       ! Parameters have been written to a temporary file. Copy them to our new group.
       call self%outputParametersContainer%copy('Parameters',outputGroup)
       fileNameTemporary=self%outputParametersContainer%name()
       call self%outputParameters         %close()
       call self%outputParametersContainer%close()
       self%outputParameters=outputGroup%openGroup('Parameters',attributesCompactMaxiumum=0)
       self%outputParametersTemporary=.false.
    else
       if (self%outputParameters%isOpen()) call self%outputParameters%close()
       self%outputParameters      =outputGroup%openGroup('Parameters',attributesCompactMaxiumum=0)
    end if
    !$ call hdf5Access%unset()
    self%outputParametersCopied=.false.
    return
  end subroutine inputParametersParametersGroupOpen

  subroutine inputParametersParametersGroupCopy(self,inputParameters_)
    !!{
    Copy an output group for parameters in the given HDF5 object.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class(inputParameters), intent(inout) :: self
    class(inputParameters), intent(in   ) :: inputParameters_

    !$ call hdf5Access%set()
    self%outputParameters      =inputParameters_%outputParameters
    self%outputParametersCopied=.true.
    !$ call hdf5Access%unset()
    return
  end subroutine inputParametersParametersGroupCopy

  subroutine inputParametersValidateName(self,parameterName)
    !!{
    Validate a parameter name.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class    (inputParameters), intent(in   ) :: self
    character(len=*          ), intent(in   ) :: parameterName
    !$GLC attributes unused :: self

    if (trim(parameterName) == "value") call Error_Report('"value" is not a valid parameter name'//{introspection:location})
    return
  end subroutine inputParametersValidateName

  function inputParametersNode(self,parameterName,requireValue,copyInstance,writeOutput)
    !!{
    Return the node containing the parameter.
    !!}
    use :: FoX_DOM           , only : ELEMENT_NODE   , getNodeName   , getNodeType     , hasAttribute, &
          &                           node           , getTextContent, getAttributeNode
    use :: Error             , only : Error_Report
    use :: IO_XML            , only : XML_Path_Exists
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : assignment(=)  , char
    implicit none
    type     (inputParameter ), pointer                 :: inputParametersNode
    class    (inputParameters), intent(inout)           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue       , writeOutput
    integer                   , intent(in   ), optional :: copyInstance
    type     (node           ), pointer                 :: node_              , identifierNode
    integer                                             :: skipInstances
    type     (varying_string )                          :: attributeName
    !![
    <optionalArgument name="requireValue" defaultsTo=".true."/>
    <optionalArgument name="writeOutput"  defaultsTo=".true."/>
    <optionalArgument name="copyInstance" defaultsTo="1"     />
    !!]

    call self%validateName(parameterName)
    !$omp critical (FoX_DOM_Access)
    inputParametersNode => self%parameters%firstChild
    skipInstances=copyInstance_-1
    do while (associated(inputParametersNode))
       if (.not.inputParametersNode%removed.and.inputParametersNode%active) then
          node_ => inputParametersNode%content
          if (getNodeType(node_) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(node_)) then
             if     (                                &
                  &  .not.requireValue_              &
                  &  .or.                            &
                  &      hasAttribute(node_,'value') &
                  &  .or.                            &
                  &   XML_Path_Exists(node_,"value") &
                  &  .or.                            &
                  &      hasAttribute(node_,"idRef") &
                  & ) then
                if (skipInstances > 0) then
                   skipInstances=skipInstances-1
                else
                   exit
                end if
             end if
          end if
       end if
       inputParametersNode => inputParametersNode%sibling
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.associated(inputParametersNode)) call Error_Report('parameter node ['//trim(parameterName)//'] not found'//{introspection:location})
    if (associated(inputParametersNode%referenced)) then
       ! We have a reference to another parameter.
       !$omp critical (FoX_DOM_Access)
       attributeName  =  getNodeName     (inputParametersNode%content        )
       identifierNode => getAttributeNode(inputParametersNode%content,'idRef')
       !$omp end critical (FoX_DOM_Access)
       !$ call hdf5Access%set()
       if (self%outputParameters%isOpen()) then
          if (.not.self%outputParameters%hasAttribute(char(attributeName))) then
             call self%outputParameters%writeAttribute("{idRef:"//getTextContent(identifierNode)//"}",char(attributeName))
          end if
       end if
       !$ call hdf5Access%unset()
       ! Return the referenced parameter.
       inputParametersNode => inputParametersNode%referenced
    end if
    return
  end function inputParametersNode

  logical function inputParametersIsPresent(self,parameterName,requireValue,searchInParents)
    !!{
    Return true if the specified parameter is present.
    !!}
    use :: FoX_dom, only : ELEMENT_NODE   , getNodeName, getNodeType, hasAttribute, &
          &                node
    use :: IO_XML , only : XML_Path_Exists
    implicit none
    class    (inputParameters), intent(in   )           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue    , searchInParents
    type     (node           ), pointer                 :: node_
    type     (inputParameter ), pointer                 :: currentParameter, currentParent
    !![
    <optionalArgument name="requireValue"    defaultsTo=".true." />
    <optionalArgument name="searchInParents" defaultsTo=".false."/>
    !!]

    call self%validateName(parameterName)
    inputParametersIsPresent=.false.
    if (.not.associated(self%parameters)) return
    !$omp critical (FoX_DOM_Access)
    currentParent => self%parameters
    do while (associated(currentParent))
       currentParameter => currentParent%firstChild
       do while (associated(currentParameter))
          if (.not.currentParameter%removed.and.currentParameter%active) then
             node_ => currentParameter%content
             if (getNodeType(node_) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(node_)) then
                if     (                                &
                     &  .not.requireValue_              &
                     &  .or.                            &
                     &      hasAttribute(node_,'value') &
                     &  .or.                            &
                     &   XML_Path_Exists(node_,"value") &
                     &  .or.                            &
                     &      hasAttribute(node_,"idRef") &
                     & ) then
                   inputParametersIsPresent=.true.
                   exit
                end if
             end if
          end if
          currentParameter => currentParameter%sibling
       end do
       if (searchInParents_ .and. .not. inputParametersIsPresent) then
          currentParent => currentParent%parent
       else
          currentParent => null()
       end if
    end do
    !$omp end critical (FoX_DOM_Access)
    return
  end function inputParametersIsPresent

  integer function inputParametersCopiesCount(self,parameterName,requireValue,zeroIfNotPresent)
    !!{
    Return true if the specified parameter is present.
    !!}
    use :: FoX_dom, only : ELEMENT_NODE   , getNodeName, getNodeType, hasAttribute, &
          &                node
    use :: Error  , only : Error_Report
    use :: IO_XML , only : XML_Path_Exists
    implicit none
    class    (inputParameters), intent(in   )           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue    , zeroIfNotPresent
    type     (node           ), pointer                 :: node_
    type     (inputParameter ), pointer                 :: currentParameter
    !![
    <optionalArgument name="zeroIfNotPresent" defaultsTo=".false."/>
    <optionalArgument name="requireValue"     defaultsTo=".true." />
    !!]

    call self%validateName(parameterName)
    if (self%isPresent(parameterName,requireValue_)) then
       inputParametersCopiesCount=0
       !$omp critical (FoX_DOM_Access)
       currentParameter => self%parameters%firstChild
       do while (associated(currentParameter))
          if (.not.currentParameter%removed.and.currentParameter%active) then
             node_ => currentParameter%content
             if (getNodeType(node_) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(node_)) then
                if     (                                  &
                     &   .not.requireValue_               &
                     &  .or.                              &
                     &   (                                &
                     &        hasAttribute(node_,'value') &
                     &    .or.                            &
                     &     XML_Path_Exists(node_,"value") &
                     &    .or.                            &
                     &        hasAttribute(node_,"idRef") &
                     &   )                                &
                     & )                                  &
                     & inputParametersCopiesCount=inputParametersCopiesCount+1
             end if
          end if
          currentParameter => currentParameter%sibling
      end do
       !$omp end critical (FoX_DOM_Access)
    else if (zeroIfNotPresent_) then
       inputParametersCopiesCount=0
    else
       inputParametersCopiesCount=0
       call Error_Report('parameter ['//parameterName//'] is not present'//{introspection:location})
    end if
    return
  end function inputParametersCopiesCount

  integer function inputParametersCount(self,parameterName,zeroIfNotPresent)
    !!{
    Return a count of the number of values in a parameter.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    use :: String_Handling   , only : String_Count_Words
    implicit none
    class    (inputParameters), intent(inout)           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: zeroIfNotPresent
    type     (varying_string )                          :: parameterText
    !![
    <optionalArgument name="zeroIfNotPresent" defaultsTo=".false." />
    !!]

    if (self%isPresent(parameterName)) then
       call self%value(parameterName,parameterText,writeOutput=.false.,evaluate=.false.)
       inputParametersCount=String_Count_Words(char(parameterText))
    else
       if (zeroIfNotPresent_) then
          inputParametersCount=0
       else
          inputParametersCount=0
          call Error_Report('parameter ['//parameterName//'] is not present'//{introspection:location})
       end if
    end if
    return
  end function inputParametersCount

  subroutine inputParametersFindParent(self,parameterPath,parent,parameterName)
    !!{
    Return the parent containing the specified parameter path.
    !!}
    use :: Error          , only : Error_Report
    use :: String_Handling, only : String_Count_Words, String_Split_Words
    implicit none
    class    (inputParameters           ), intent(in   ), target      :: self
    character(len=*                     ), intent(in   )              :: parameterPath
    type     (inputParameters           ), intent(  out), pointer     :: parent
    character(len=parameterLengthMaximum), intent(  out)              :: parameterName
    type     (inputParameters           )               , pointer     :: rootParameters   , subParameters, &
         &                                                               subParametersNext
    character(len=parameterLengthMaximum), dimension(:) , allocatable :: parameterNames
    integer                                                           :: countNames       , i
    
    countNames=String_Count_Words(parameterPath,":")
    allocate(parameterNames(countNames))
    call String_Split_Words(parameterNames,parameterPath,":")
    parameterName  =  parameterNames(countNames)
    subParameters  => null          (          )
    if (trim(parameterNames(1)) /= "." .and. trim(parameterNames(1)) /= "..") then
       ! Path is absolute - start from the original parameters.
       if (associated(self%original)) then
          rootParameters => self%original
       else
          rootParameters => self
       end if
    else
       ! Path is relative - simply start from the current parameter.
       rootParameters => self
    end if
    do i=1,countNames-1
       if      (trim(parameterNames(i)) == "." ) then
          ! Self - no need to move.
          if (i == 1) then
             allocate(subParameters)
             subParameters=inputParameters(rootParameters)
          end if
       else if (trim(parameterNames(i)) == "..") then
          ! Move to the parent parameter.
          if (i == 1) then
             if (.not.associated(rootParameters%parent)) call Error_Report('no parent parameter exists'//{introspection:location})
             allocate(subParameters)
             subParameters    =inputParameters(rootParameters%parent)
          else
             if (.not.associated( subParameters%parent)) call Error_Report('no parent parameter exists'//{introspection:location})
             allocate(subParametersNext)
             subParametersNext=inputParameters(subParameters %parent)
             deallocate(subParameters)
             subParameters => subParametersNext
          end if
       else
          ! Move to the named parameter.
          if (i == 1) then
             allocate(subParameters)
             subParameters    =rootParameters%subParameters(trim(parameterNames(i)),requireValue=.false.)
          else
             allocate(subParametersNext)
             subParametersNext=subParameters%subParameters(trim(parameterNames(i)),requireValue=.false.)
             deallocate(subParameters)
             subParameters => subParametersNext
          end if
       end if
    end do
    allocate(parent)
    if (countNames == 1) then
       parent=inputParameters(rootParameters)
    else
       parent=inputParameters( subParameters)
       if (associated(subParameters)) deallocate(subParameters)
    end if
    return
  end subroutine inputParametersFindParent

  function inputParametersSubParameters(self,parameterName,requireValue,requirePresent,copyInstance)
    !!{
    Return sub-parameters of the specified parameter.
    !!}
    use :: FoX_dom           , only : node
    use :: Error             , only : Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : ioHDF5AccessInitialize
    use :: ISO_Varying_String, only : assignment(=)         , char, operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    type     (inputParameters)                          :: inputParametersSubParameters
    class    (inputParameters), intent(inout), target   :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue                , requirePresent
    integer                   , intent(in   ), optional :: copyInstance
    type     (inputParameter ), pointer                 :: parameterNode
    integer                                             :: copyCount                   , copyInstance_
    type     (varying_string )                          :: groupName
    !![
    <optionalArgument name="requirePresent" defaultsTo=".true." />
    !!]

    ! The HDF5 access lock may not yet have been initialized. Ensure it is before using it.
    call ioHDF5AccessInitialize()
    if (.not.self%isPresent(parameterName,requireValue)) then
       if (requirePresent_) then
          call Error_Report('parameter ['//trim(parameterName)//'] not found'//{introspection:location})
       else
          inputParametersSubParameters=inputParameters()
       end if
       copyCount                                  =  1
    else
       copyCount                                  =  self%copiesCount(parameterName        ,requireValue=requireValue                          )
       parameterNode                              => self%node       (parameterName        ,requireValue=requireValue,copyInstance=copyInstance)
       if (associated(parameterNode%referenced)) parameterNode => parameterNode%referenced
       inputParametersSubParameters               =  inputParameters (parameterNode%content,noOutput    =.true.      ,noBuild     =.true.      )
       inputParametersSubParameters%parameters    => parameterNode
    end if
    inputParametersSubParameters%parent           => self
    if (associated(self%original)) then
       inputParametersSubParameters%original => self%original
    else
       inputParametersSubParameters%original => self
    end if
    !$ call hdf5Access%set()
    if (self%outputParameters%isOpen()) then
       groupName=parameterName
       if (copyCount > 1) then
          if (present(copyInstance)) then
             copyInstance_=copyInstance
          else
             copyInstance_=1
          end if
          groupName=groupName//"["//copyInstance//"]"
       end if
       inputParametersSubParameters%outputParameters=self%outputParameters%openGroup(char(groupName))
    end if
    !$ call hdf5Access%unset()
    return
  end function inputParametersSubParameters

  function inputParametersPath(self)
    !!{
    Return the path to the given parameters.
    !!}
    use :: FoX_dom           , only : getNodeName
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    implicit none
    type (varying_string )                        :: inputParametersPath
    class(inputParameters), intent(inout), target :: self
    type (inputParameter ), pointer               :: parameterNode
    
    parameterNode       => self%parameters
    inputParametersPath =  ""
    do while (associated(parameterNode))
       if (associated(parameterNode%content)) inputParametersPath=inputParametersPath//getNodeName(parameterNode%content)//"/"
       parameterNode => parameterNode%parent
    end do
    return
  end function inputParametersPath

  recursive subroutine inputParametersValueName{Type¦label}(self,parameterName,parameterValue,defaultValue,errorStatus,writeOutput,copyInstance,evaluate)
    !!{
    Return the value of the parameter specified by name.
    !!}
    use :: Error             , only : Error_Report, Warn
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : char
    use :: MPI_Utilities     , only : mpiSelf
    implicit none
    class           (inputParameters                         ), intent(inout), target   :: self
    character       (len=*                                   ), intent(in   )           :: parameterName
    {Type¦intrinsic}                                          , intent(  out)           :: parameterValue
    {Type¦intrinsic}                                          , intent(in   ), optional :: defaultValue
    type            (enumerationInputParameterErrorStatusType), intent(  out), optional :: errorStatus
    integer                                                   , intent(in   ), optional :: copyInstance
    logical                                                   , intent(in   ), optional :: writeOutput   , evaluate
    type            (inputParameters                         ), pointer                 :: parametersRoot
    type            (inputParameter                          ), pointer                 :: parameterNode
    type            (varying_string                          )                          :: parameterPath
    !![
    <optionalArgument name="writeOutput" defaultsTo=".true." />
    !!]

    if (self%isPresent(parameterName)) then
       parameterNode => self%node(parameterName,copyInstance=copyInstance,writeOutput=writeOutput)
       call self%value(parameterNode,parameterValue,errorStatus,writeOutput,evaluate)
    else if (present(defaultValue)) then
       parametersRoot => self
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       parameterPath=self%path()
       call parametersRoot%lock%set()
       if (.not.parametersRoot%warnedDefaults%exists(parameterPath)) then
          if (mpiSelf%isMaster()) call Warn("Using default value for parameter '["//char(parameterPath)//parameterName//"]'")
          call parametersRoot%warnedDefaults%set(parameterPath,1)
       end if
       call parametersRoot%lock%unset()
       parameterValue=defaultValue
       ! Write the parameter file to an HDF5 object.
       if (self%outputParameters%isOpen().and.writeOutput_) then
          !$ call hdf5Access%set()
          if (.not.self%outputParameters%hasAttribute(parameterName)) call self%outputParameters%writeAttribute({Type¦outputConverter¦parameterValue},parameterName)
          !$ call hdf5Access%unset()
       end if
    else if (present(errorStatus )) then
       errorStatus   =inputParameterErrorStatusNotPresent
    else
       call Error_Report('parameter ['//parameterName//'] not present and no default given'//{introspection:location})
    end if
    return
  end subroutine inputParametersValueName{Type¦label}

  recursive subroutine inputParametersValueNode{Type¦label}(self,parameterNode,parameterValue,errorStatus,writeOutput,evaluate)
    !!{
    Return the value of the specified parameter.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_int64_t                        , c_size_t
    use            :: FoX_dom           , only : DOMException                     , getAttributeNode  , getNodeName   , hasAttribute      , &
          &                                      inException                      , node              , getTextContent, extractDataContent
    use            :: Error             , only : Error_Report
    use            :: ISO_Varying_String, only : assignment(=)                    , char              , operator(//)  , operator(==)      , &
          &                                      trim                             , var_str
    use            :: String_Handling   , only : String_Count_Words               , String_Split_Words, operator(//)
    use            :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Path_Exists
    use            :: HDF5_Access       , only : hdf5Access
    implicit none
    class           (inputParameters                         ), intent(inout), target   :: self
    type            (inputParameter                          ), intent(inout), target   :: parameterNode
    {Type¦intrinsic}                                          , intent(  out)           :: parameterValue
    type            (enumerationInputParameterErrorStatusType), intent(  out), optional :: errorStatus
    logical                                                   , intent(in   ), optional :: writeOutput        , evaluate
#ifdef MATHEVALAVAIL
    integer         (c_int64_t                               )                          :: evaluator
    ! Declarations of GNU libmatheval procedures used.
    integer         (c_int64_t                               ), external                :: Evaluator_Create_
    double precision                                          , external                :: Evaluator_Evaluate_
    external                                                                            :: Evaluator_Destroy_
#endif
    type            (inputParameter                          )               , pointer  :: sibling
    type            (node                                    )               , pointer  :: valueElement       , identifierNode
    type            (inputParameters                         )               , pointer  :: parentParameters
    type            (DOMException                            )                          :: exception
    integer                                                                             :: copyInstance       , copyCount       , &
         &                                                                                 status
    logical                                                                             :: hasValueAttribute  , hasValueElement , &
         &                                                                                 isException        , isPresent       , &
         &                                                                                 isDouble           , isText          , &
         &                                                                                 haveDefault
    character       (len=parameterLengthMaximum              )                          :: expression         , parameterName   , &
         &                                                                                 workText           , content         , &
         &                                                                                 workValueText      , formatSpecifier , &
         &                                                                                 parameterLeafName  , defaultValue
    type            (varying_string                          )                          :: attributeName      , nodeName        , &
         &                                                                                 siblingName
    double precision                                                                    :: workValueDouble
    integer         (c_size_t                                )                          :: workValueInteger
    type            (enumerationInputParameterTypeType       )                          :: parameterType
    {Type¦match¦^Long.*¦character(len=parameterLengthMaximum) :: parameterText¦}
    {Type¦match¦^(Character|VarStr)Rank1$¦type(varying_string) :: parameterText¦}
    !![
    <optionalArgument name="writeOutput" defaultsTo=".true." />
    <optionalArgument name="evaluate"    defaultsTo=".true." />
    !!]

    if (present(errorStatus)) errorStatus=inputParameterErrorStatusSuccess
    !$omp critical (FoX_DOM_Access)
    hasValueAttribute =  hasAttribute   (parameterNode%content,'value')
    hasValueElement   =  XML_Path_Exists(parameterNode%content,'value')
    !$omp end critical (FoX_DOM_Access)
    if (hasValueAttribute .and. hasValueElement) then
       if (present(errorStatus)) then
          errorStatus=inputParameterErrorStatusAmbiguousValue
       else
          call Error_Report('ambiguous value attribute and element'//{introspection:location})
       end if
    else if (hasValueAttribute .or. hasValueElement) then
       !$omp critical (FoX_DOM_Access)
       if (hasValueAttribute) then
          valueElement => getAttributeNode                 (parameterNode%content,"value"                          )
       else
          valueElement => XML_Get_First_Element_By_Tag_Name(parameterNode%content,"value",directChildrenOnly=.true.)
       end if
       content=getTextContent(valueElement)
       !$omp end critical (FoX_DOM_Access)
       if (trim(content) == "") then
          if (present(errorStatus)) then
             errorStatus=inputParameterErrorStatusEmptyValue
          else
             !$omp critical (FoX_DOM_Access)
             nodeName=getNodeName(parameterNode%content)
             !$omp end critical (FoX_DOM_Access)
             call Error_Report(                               &
                  &            'empty value in parameter ['// &
                  &             nodeName                   // &
                  &            ']'                         // &
                  &            {introspection:location}       &
                  &           )
          end if
       else
          ! Evaluate if an expression.
          !$omp critical (FoX_DOM_Access)
          expression=getTextContent(valueElement)
          !$omp end critical (FoX_DOM_Access)
          if (expression(1:1) == "=" .and. evaluate_) then
             {Type¦match¦^(Double|Character|VarStr)$¦if (.true.) then¦if (.false.) then}
                {Type¦match¦^Double$¦isDouble=.true.¦isDouble=.false.}
                {Type¦match¦^(Character|VarStr)¦isText=.true.¦isText=.false.}
                ! This is an expression, and we have a scalar, floating point or text type - it can be evaluated.             
                !! Mark this parameter as being evaluated and store its original content. This allows the parameter to be reset to
                !! its original (unevaluated) state if necessary.
                parameterNode%evaluated      =.true.
                parameterNode%contentOriginal=expression
                !! Remove the initial "=".
                expression=expression(2:len_trim(expression))
                !! Replace other parameter values inside the parameter.
                do while (index(expression,"[") /= 0)
                   parameterName=expression(index(expression,"[")+1:index(expression,"]")-1)
                   parameterType=inputParameterTypeDouble
                   if (isText) then
                      if (index(parameterName,"|") > 0) then
                         formatSpecifier=parameterName(1:index(parameterName,"|")-1                        )
                         parameterName  =parameterName(  index(parameterName,"|")+1:len_trim(parameterName))
                         if (formatSpecifier(1:1) /= "%") call Error_Report('unrecognized format specifier'//{introspection:location})
                         select case (formatSpecifier(len_trim(formatSpecifier):len_trim(formatSpecifier)))
                         case ("s"    )
                            parameterType=inputParameterTypeText
                            formatSpecifier="(a)"
                         case ("d"    )
                            parameterType=inputParameterTypeInteger
                            formatSpecifier="(i"//formatSpecifier(2:len_trim(formatSpecifier)-1)//")"
                         case ("f","e")
                            parameterType=inputParameterTypeDouble
                            formatSpecifier="("//formatSpecifier(len_trim(formatSpecifier):len_trim(formatSpecifier))//formatSpecifier(2:len_trim(formatSpecifier)-1)//")"
                         case default
                            call Error_Report('unrecognized format specifier'//{introspection:location})
                         end select
                      else
                         call Error_Report('inserted parameters must have a format specifier'//{introspection:location})
                      end if
                   else if (isDouble) then
                      haveDefault=index(parameterName,"|") > 0
                      if (haveDefault) then
                         defaultValue =parameterName(  index(parameterName,"|")+1:len_trim(parameterName))
                         parameterName=parameterName(1:index(parameterName,"|")-1                        )
                      end if
                   end if
                   !! Find the named parameter's parent and extract the value from it.
                   call self%findParent(parameterName,parentParameters,parameterLeafName)
                   isPresent=parentParameters%isPresent(trim(parameterLeafName))
                   if (isPresent) then
                      if      (parameterType == inputParameterTypeDouble ) then
                         call parentParameters%value(trim(parameterLeafName),workValueDouble )
                      else if (parameterType == inputParameterTypeInteger) then
                         call parentParameters%value(trim(parameterLeafName),workValueInteger)
                      else if (parameterType == inputParameterTypeText   ) then
                         call parentParameters%value(trim(parameterLeafName),workValueText   )
                      end if
                      deallocate(parentParameters)
                   else if (isDouble .and. haveDefault) then
                      ! The parameter does not exist, but a default value is available - use that default.
                      read (defaultValue,*) workValueDouble
                      deallocate(parentParameters)
                   else
                      !$omp critical (FoX_DOM_Access)
                      expression=getTextContent(valueElement)
                      !$omp end critical (FoX_DOM_Access)
                      call Error_Report('parameter `'//trim(parameterName)//'` referenced in expression `'//trim(expression)//'` does not exist'//{introspection:location})
                   end if
                   if (isDouble) then
                      write (workText,'(e24.16)') workValueDouble
                   else if (isText) then
                      if      (parameterType == inputParameterTypeDouble ) then
                         write (worktext,formatSpecifier) workValueDouble
                      else if (parameterType == inputParameterTypeInteger) then
                         write (worktext,formatSpecifier) workValueInteger
                      else if (parameterType == inputParameterTypeText   ) then
                         write (worktext,formatSpecifier) workValueText
                      end if
                   end if
                   expression=expression(1:index(expression,"[")-1)//trim(adjustl(workText))//expression(index(expression,"]")+1:len_trim(expression))
                end do
                !! Evaluate the expression.
                if (isDouble) then
#ifdef MATHEVALAVAIL
                   evaluator=Evaluator_Create_(trim(expression))
                   if (evaluator == 0) call Error_Report("failed to parse expression '"//trim(expression)//"' - see https://galacticusorg.github.io/libmatheval/doc/evaluator_005fcreate.html for supported operators and functions"//{introspection:location})
                   workValueDouble=Evaluator_Evaluate_(evaluator,0,"",0.0d0)
                   call Evaluator_Destroy_(evaluator)
                   call parameterNode%set(workValueDouble)
                   !$omp critical (FoX_DOM_Access)
                   valueElement => getAttributeNode(parameterNode%content,"value")
                   !$omp end critical (FoX_DOM_Access)
#else
                   call Error_Report('derived parameters require libmatheval, but it is not installed'//{introspection:location})
#endif
                else if (isText) then
                   call parameterNode%set(var_str(trim(expression)))
                   !$omp critical (FoX_DOM_Access)
                   valueElement => getAttributeNode(parameterNode%content,"value")
                   !$omp end critical (FoX_DOM_Access)
                end if
             end if
          end if
          ! Extract the value.          
          status=0
          !$omp critical (FoX_DOM_Access)
          {Type¦match¦^(?!(Long|VarStr|Character))¦call extractDataContent(valueElement,parameterValue,iostat=status,ex=exception)¦}
          {Type¦match¦^(Long.*|(Character|VarStr)Rank1)$¦parameterText=getTextContent(valueElement,ex=exception)¦}
          {Type¦match¦^(VarStr|Character)$¦parameterValue=getTextContent(valueElement,ex=exception)¦}
          isException=inException(exception)
          !$omp end critical (FoX_DOM_Access)
          if (isException .or. status /= 0) then
             if (present(errorStatus)) then
                errorStatus=inputParameterErrorStatusParse
             else
                !$omp critical (FoX_DOM_Access) 
                attributeName=getNodeName   (parameterNode%content)
                content      =getTextContent(valueElement         )
                !$omp end critical (FoX_DOM_Access) 
                call Error_Report(                                &
                     &            'unable to parse parameter ['// &
                     &             attributeName               // &
                     &            ']='                         // &
                     &             content                     // &
                     &             {introspection:location}       &
                     &           )
             end if
          else
             ! Convert type as necessary.
             {Type¦match¦^Long.*¦read (parameterText,*) parameterValue¦}
             {Type¦match¦^(Character|VarStr)Rank1$¦call String_Split_Words(parameterValue,char(parameterText))¦}
             ! Write the parameter file to an HDF5 object.
             if (self%outputParameters%isOpen().and.writeOutput_) then
                !$omp critical (FoX_DOM_Access)
                attributeName=getNodeName(parameterNode%content)
                !$omp end critical (FoX_DOM_Access)
                copyCount    =self%copiesCount(char(attributeName))
                if (copyCount > 1) then
                   copyInstance =  copyCount
                   sibling      => parameterNode%sibling
                   do while (associated(sibling))
                      siblingName=getNodeName(sibling%content)
                      if (siblingName == attributeName) &
                           & copyInstance=copyInstance-1
                      sibling      => sibling%sibling
                   end do
                   attributeName=attributeName//"["//copyInstance//"]"
                end if
                if (hasAttribute(parameterNode%content,'id')) then
                   identifierNode => getAttributeNode(parameterNode%content,'id')
                   attributeName  =  attributeName//"{id:"//getTextContent(identifierNode)//"}"
                end if
                !$ call hdf5Access%set()
                if (.not.self%outputParameters%hasAttribute(char(attributeName))) call self%outputParameters%writeAttribute({Type¦outputConverter¦parameterValue},char(attributeName))
                !$ call hdf5Access%unset()
             end if
          end if
       end if
    end if
    return
  end subroutine inputParametersValueNode{Type¦label}

  function inputParameterListConstructor() result(self)
    !!{
    Construct an {\normalfont \ttfamily inputParameterList} object.
    !!}
    implicit none
    type(inputParameterList) :: self

    self%count=0
    return
  end function inputParameterListConstructor

  subroutine inputParameterListDestructor(self)
    !!{
    Destroy an {\normalfont \ttfamily inputParameterList} object.
    !!}
    implicit none
    type(inputParameterList), intent(inout) :: self

    if (allocated(self%name )) deallocate(self%name )
    if (allocated(self%value)) deallocate(self%value)
    return
  end subroutine inputParameterListDestructor

  subroutine inputParameterListAdd(self,name,value)
    !!{
    Add a parameter to a list of input parameters to an XML document.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class    (inputParameterList), intent(inout)               :: self
    character(len=*             ), intent(in   )               :: name             , value
    type     (varying_string    ), allocatable  , dimension(:) :: nameTemporary    , valueTemporary
    integer                      , parameter                   :: sizeIncrement =10

    if (.not.allocated(self%name)) then
       allocate(self%name (sizeIncrement))
       allocate(self%value(sizeIncrement))
       self%count=0
    else if (self%count == size(self%name)) then
       call Move_Alloc(self%name , nameTemporary)
       call Move_Alloc(self%value,valueTemporary)
       allocate(self%name (size( nameTemporary)+sizeIncrement))
       allocate(self%value(size(valueTemporary)+sizeIncrement))
       self%name (1:size( nameTemporary))= nameTemporary
       self%value(1:size(valueTemporary))=valueTemporary
       deallocate( nameTemporary)
       deallocate(valueTemporary)
    end if
    self%count=self%count+1
    self%name (self%count)=name
    self%value(self%count)=value
    return
  end subroutine inputParameterListAdd

  subroutine inputParameterListSerializeToXML(self,parameterDoc)
    !!{
    Serialize a list of input parameters to an XML document.
    !!}
    use :: FoX_wXML          , only : xml_AddAttribute, xml_EndElement, xml_NewElement, xmlf_t
    use :: ISO_Varying_String, only : char
    class (inputParameterList), intent(in   ) :: self
    type  (xmlf_t            ), intent(inout) :: parameterDoc
    integer                                   :: i

    if (allocated(self%name)) then
       do i=1,self%count
          call xml_NewElement  (parameterDoc,        char(self%name (i)))
          call xml_AddAttribute(parameterDoc,"value",char(self%value(i)))
          call xml_EndElement  (parameterDoc,        char(self%name (i)))
       end do
    end if
    return
  end subroutine inputParameterListSerializeToXML

  recursive function inputParametersSerializeToString(self,hashed)
    !!{
    Serialize input parameters to a string.
    !!}
    use :: FoX_dom             , only : getNodeName  , node
    use :: Hashes_Cryptographic, only : Hash_MD5
    use :: ISO_Varying_String  , only : assignment(=), char, operator(//)
    implicit none
    type   (varying_string                          )                          :: inputParametersSerializeToString
    class  (inputParameters                         ), intent(inout)           :: self
    logical                                          , intent(in   ), optional :: hashed
    type   (node                                    ), pointer                 :: node_
    type   (inputParameter                          ), pointer                 :: currentParameter
    type   (inputParameters                         ), allocatable             :: subParameters
    type   (enumerationInputParameterErrorStatusType)                          :: errorStatus
    type   (varying_string                          )                          :: parameterValue                  , nodeName
    logical                                                                    :: firstParameter
    !![
    <optionalArgument name="hashed" defaultsTo=".false." />
    !!]

    inputParametersSerializeToString=""
    firstParameter=.true.
    currentParameter => self%parameters%firstChild
    do while (associated(currentParameter))
       if (currentParameter%isParameter()) then
          if (.not.firstParameter) inputParametersSerializeToString=inputParametersSerializeToString//"_"
          call self%value(currentParameter,parameterValue,errorStatus,writeOutput=.false.)
          node_ => currentParameter%content
          !$omp critical (FoX_DOM_Access)
          nodeName=getNodeName(node_)
          !$omp end critical (FoX_DOM_Access)
          inputParametersSerializeToString=inputParametersSerializeToString// &
               &                           nodeName                        // &
               &                           ":"                             // &
               &                           adjustl(char(parameterValue))
          !$omp critical (FoX_DOM_Access)
          nodeName=getNodeName(node_)
          !$omp end critical (FoX_DOM_Access)
          allocate(subParameters)
          subParameters=self%subParameters(char(nodeName))
          if (associated(subParameters%parameters%firstChild)) inputParametersSerializeToString=inputParametersSerializeToString // &
               &                                                                                "{"                              // &
               &                                                                                subParameters%serializeToString()// &
               &                                                                                "}"
          deallocate(subParameters)
          firstParameter=.false.
       end if
       currentParameter => currentParameter%sibling
    end do
    if (hashed_) inputParametersSerializeToString=Hash_MD5(inputParametersSerializeToString)
    return
  end function inputParametersSerializeToString

  subroutine inputParametersSerializeToXML(self,parameterFile)
    !!{
    Serialize input parameters to an XML file.
    !!}
    use :: FoX_DOM           , only : serialize
    use :: ISO_Varying_String, only : char
    implicit none
    class(inputParameters), intent(in   ) :: self
    type (varying_string ), intent(in   ) :: parameterFile

    !$omp critical (FoX_DOM_Access)
    call serialize(self%document,char(parameterFile))
    !$omp end critical (FoX_DOM_Access)
    return
  end subroutine inputParametersSerializeToXML

  subroutine inputParametersAddParameter(self,parameterName,parameterValue,writeOutput)
    !!{
    Add a parameter to the set.
    !!}
    use :: FoX_DOM    , only : ELEMENT_NODE, appendChild, createElementNS, getNamespaceURI, &
          &                    getNodeName , getNodeType, hasAttribute   , node           , &
          &                    setAttribute
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class    (inputParameters), intent(inout)           :: self
    character(len=*          ), intent(in   )           :: parameterName   , parameterValue
    logical                   , intent(in   ), optional :: writeOutput
    type     (node           ), pointer                 :: parameterNode   , dummy
    type     (inputParameter ), pointer                 :: currentParameter
    logical                                             :: previouslyAdded
    !![
    <optionalArgument name="writeOutput" defaultsTo=".true."/>
    !!]
    
    ! Check if the parameter was previously added but then removed, in which case it still exists but is just flagged as removed.
    previouslyAdded=.false.
    if (associated(self%parameters)) then
       !$omp critical(FoX_DOM_Access)
       currentParameter => self%parameters%firstChild
       do while (associated(currentParameter))
          parameterNode => currentParameter%content
          if (getNodeType(parameterNode) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(parameterNode)) then
             if (.not.hasAttribute(parameterNode,'id')) then
                previouslyAdded=.true.
                exit
             end if
          end if
          currentParameter => currentParameter%sibling
       end do
       !$omp end critical(FoX_DOM_Access)
    end if
    if (previouslyAdded) then
       !$omp critical(FoX_DOM_Access)
       call setAttribute(parameterNode,"value",trim(parameterValue))
       !$omp end critical(FoX_DOM_Access)
    else
       !$omp critical(FoX_DOM_Access)
       parameterNode   => createElementNS(self%document,getNamespaceURI(self%document),parameterName)
       call setAttribute(parameterNode,"value",trim(parameterValue))
       dummy           => appendChild  (self%rootNode,parameterNode)
       !$omp end critical(FoX_DOM_Access)
       ! Write the parameter file to an HDF5 object.
       if (self%outputParameters%isOpen().and.writeOutput_) then
          !$ call hdf5Access%set()
          if (.not.self%outputParameters%hasAttribute(parameterName)) call self%outputParameters%writeAttribute(parameterValue,parameterName)
          !$ call hdf5Access%unset()
       end if
       if (associated(self%parameters)) then
          if (associated(self%parameters%firstChild)) then
             currentParameter => self%parameters%firstChild
             do while (associated(currentParameter%sibling))
                currentParameter => currentParameter%sibling
             end do
             allocate(currentParameter%sibling)
             currentParameter => currentParameter%sibling
          else
             allocate(self%parameters%firstChild)
             currentParameter => self%parameters%firstChild
          end if
       else
          allocate(self%parameters           )
          allocate(self%parameters%firstChild)
          currentParameter => self%parameters%firstChild
       end if
       currentParameter%content    => parameterNode
       currentParameter%parent     => self%parameters
       currentParameter%firstChild => null()
       currentParameter%sibling    => null()
       currentParameter%referenced => null()
       currentParameter%created    =  .true.
       currentParameter%evaluated  =  .false.
    end if
    currentParameter%removed=.false.
    currentParameter%active =.true.
    return
  end subroutine inputParametersAddParameter

  subroutine inputParametersReset(self)
    !!{
    Reset all objects in a parameter set.
    !!}
    implicit none
    class(inputParameters), intent(inout) :: self

    call self%parameters%reset(evaluations=.true.)
    return
  end subroutine inputParametersReset

  subroutine inputParametersLockReinitialize(self)
    !!{
    Reinitialize the OpenMP lock.
    !!}
    implicit none
    class(inputParameters), intent(inout) :: self

   if (associated(self%lock)) then
       nullify(self%lock)
       allocate(self%lock)
       self%lock=ompLock()
    end if
    return
  end subroutine inputParametersLockReinitialize

end module Input_Parameters
