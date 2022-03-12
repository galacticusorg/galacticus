!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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

module Input_Parameters
  !!{
  Implements reading of parameters from an XML file.
  !!}
  use :: FoX_dom           , only : node
  use :: Function_Classes  , only : functionClass
  use :: IO_HDF5           , only : hdf5Object
  use :: ISO_Varying_String, only : varying_string
  use :: Kind_Numbers      , only : kind_int8
  use :: String_Handling   , only : char
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
     type   (node             ), pointer                   :: content
     type   (inputParameter   ), pointer    , public       :: parent                 , firstChild        , &
          &                                                   sibling                , referenced
     type   (genericObjectList), allocatable, dimension(:) :: objects
     type   (varying_string   )                            :: contentOriginal
     logical                                               :: created        =.false., removed   =.false.,&
          &                                                   evaluated      =.false.
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

  type :: inputParameters
     private
     type   (node           ), pointer, public :: document               => null()
     type   (node           ), pointer         :: rootNode               => null()
     type   (hdf5Object     )                  :: outputParameters                 , outputParametersContainer
     type   (inputParameter ), pointer         :: parameters             => null()
     type   (inputParameters), pointer, public :: parent                 => null()
     logical                                   :: outputParametersCopied =  .false., outputParametersTemporary=.false., &
          &                                       isNull                 =  .false.
        contains
     !![
     <methods>
       <method description="Build a tree of {\normalfont \ttfamily inputParameter} objects from the structure of an XML parameter file." method="buildTree" />
       <method description="Resolve references in the tree of {\normalfont \ttfamily inputParameter} objects." method="resolveReferences" />
       <method description="Open an output group for parameters in the given HDF5 object." method="parametersGroupOpen" />
       <method description="Copy the HDF5 output group for parameters from another parameters object." method="parametersGroupCopy" />
       <method description="Check that a given parameter name is a valid name, aborting if not." method="validateName" />
       <method description="Check that parameters are valid and, optionally, check if they match expected names in the provided list." method="checkParameters" />
       <method description="Return the XML node containing the named parameter." method="node" />
       <method description="Return true if the named parameter is present in the set." method="isPresent" />
       <method description="Return a count of the number copies of the named parameter. If the parameter is not present, this function aborts, unless {\normalfont \ttfamily zeroIfNotPresent} is set to {\normalfont \ttfamily true}, in which case a result of 0 is returned." method="copiesCount" />
       <method description="Return a count of the number of values in the named parameter. If the parameter is not present, this function aborts, unless {\normalfont \ttfamily zeroIfNotPresent} is set to {\normalfont \ttfamily true}, in which case a result of 0 is returned." method="count" />
       <method description="Return the set of subparameters of the named parameter." method="subParameters" />
       <method description="Return the value of a parameter specified by name or XML node. A default value can be specified only if the parameter is specified by name. Supported types include rank-0 and rank-1 logicals, integers, long integers, doubles, characters, and varying strings." method="value" />
       <method description="Serialize input parameters to a string." method="serializeToString" />
       <method description="Serialize input parameters to an XML file." method="serializeToXML" />
       <method description="Add a parameter." method="addParameter" />
       <method description="Reset all objects in this parameter set." method="reset" />
       <method description="Destroy the parameters document." method="destroy" />
     </methods>
     !!]
     final     ::                        inputParametersFinalize
     procedure :: buildTree           => inputParametersBuildTree
     procedure :: resolveReferences   => inputParametersResolveReferences
     procedure :: destroy             => inputParametersDestroy
     procedure :: parametersGroupOpen => inputParametersParametersGroupOpen
     procedure :: parametersGroupCopy => inputParametersParametersGroupCopy
     procedure :: validateName        => inputParametersValidateName
     procedure :: checkParameters     => inputParametersCheckParameters
     procedure :: node                => inputParametersNode
     procedure :: isPresent           => inputParametersIsPresent
     procedure :: copiesCount         => inputParametersCopiesCount
     procedure :: count               => inputParametersCount
     procedure :: subParameters       => inputParametersSubParameters
     procedure ::                        inputParametersValueName{Type¦label}
     procedure ::                        inputParametersValueNode{Type¦label}
     generic   :: value               => inputParametersValueName{Type¦label}
     generic   :: value               => inputParametersValueNode{Type¦label}
     procedure :: serializeToString   => inputParametersSerializeToString
     procedure :: serializeToXML      => inputParametersSerializeToXML
     procedure :: addParameter        => inputParametersAddParameter
     procedure :: reset               => inputParametersReset
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
       <method description="Serialize a list of input parameters to an XML document." method="serializeToXML" />
       <method description="Add a parameter and value to the list." method="add" />
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

  ! Maximum length allowed for parameter entries.
  integer, parameter :: parameterLengthMaximum=1024

  ! Interface to the (auto-generated) knownParameterNames() function.
  interface
     subroutine knownParameterNames(names)
       import varying_string
       type(varying_string), dimension(:), allocatable, intent(inout) :: names
     end subroutine knownParameterNames
  end interface

contains

  function inputParametersConstructorNull()
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class creating a null instance.
    !!}
    use :: FoX_dom, only : createDocument, getDocumentElement, getImplementation, setLiveNodeLists
    implicit none
    type(inputParameters) :: inputParametersConstructorNull

    inputParametersConstructorNull%document   => createDocument    (                                  &
         &                                                          getImplementation()             , &
         &                                                          qualifiedName      ="parameters", &
         &                                                          docType            =null()        &
         &                                                         )
    inputParametersConstructorNull%rootNode   => getDocumentElement(inputParametersConstructorNull%document)
    inputParametersConstructorNull%parameters => null()
    inputParametersConstructorNull%isNull     = .true.
    !$omp critical (FoX_DOM_Access)
    call setLiveNodeLists(inputParametersConstructorNull%document,.false.)
    !$omp end critical (FoX_DOM_Access)
   return
  end function inputParametersConstructorNull

  function inputParametersConstructorVarStr(xmlString,allowedParameterNames,outputParametersGroup,noOutput)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from an XML file
    specified as a variable length string.
    !!}
    use :: FoX_dom           , only : node
    use :: ISO_Varying_String, only : char, extract, operator(==)
    use :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, parseString => parseStringTS
    use :: FoX_dom           , only : node
    use :: ISO_Varying_String, only : char, extract, operator(==)
    implicit none
    type     (inputParameters)                                           :: inputParametersConstructorVarStr
    type     (varying_string    )              , intent(in   )           :: xmlString
    character(len=*             ), dimension(:), intent(in   ), optional :: allowedParameterNames
    type     (hdf5Object        ), target      , intent(in   ), optional :: outputParametersGroup
    logical                                    , intent(in   ), optional :: noOutput
    type     (node              ), pointer                               :: parameterNode

    ! Check if we have been passed XML or a file name.
    if (extract(xmlString,1,1) == "<") then
       ! Parse the string.
       !$omp critical (FoX_DOM_Access)
       parameterNode => parseString(char(xmlString))
       !$omp end critical (FoX_DOM_Access)
       inputParametersConstructorVarStr=inputParametersConstructorNode      (                                                 &
            &                                                                XML_Get_First_Element_By_Tag_Name(               &
            &                                                                                                  parameterNode, &
            &                                                                                                  'parameters'   &
            &                                                                                                 )             , &
            &                                                                allowedParameterNames                          , &
            &                                                                outputParametersGroup                          , &
            &                                                                noOutput                                         &
            &                                                               )
    else
       inputParametersConstructorVarStr=inputParametersConstructorFileVarStr(                                                 &
            &                                                                xmlString                                      , &
            &                                                                allowedParameterNames                          , &
            &                                                                outputParametersGroup                          , &
            &                                                                noOutput                                         &
            &                                                               )
    end if
    return
  end function inputParametersConstructorVarStr

  function inputParametersConstructorFileVarStr(fileName,allowedParameterNames,outputParametersGroup,noOutput)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from an XML file
    specified as a variable length string.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    type     (inputParameters)                                           :: inputParametersConstructorFileVarStr
    type     (varying_string    )              , intent(in   )           :: fileName
    character(len=*             ), dimension(:), intent(in   ), optional :: allowedParameterNames
    type     (hdf5Object        ), target      , intent(in   ), optional :: outputParametersGroup
    logical                                    , intent(in   ), optional :: noOutput

    inputParametersConstructorFileVarStr=inputParametersConstructorFileChar(                       &
         &                                                                  char(fileName)       , &
         &                                                                  allowedParameterNames, &
         &                                                                  outputParametersGroup, &
         &                                                                  noOutput               &
         &                                                                 )
    return
  end function inputParametersConstructorFileVarStr

  function inputParametersConstructorFileChar(fileName,allowedParameterNames,outputParametersGroup,noOutput)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from an XML file
    specified as a character variable.
    !!}
    use :: File_Utilities  , only : File_Exists
    use :: FoX_dom         , only : node
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: IO_XML          , only : XML_Get_First_Element_By_Tag_Name, XML_Parse
    implicit none
    type     (inputParameters)                                        :: inputParametersConstructorFileChar
    character(len=*          )              , intent(in   )           :: fileName
    character(len=*          ), dimension(:), intent(in   ), optional :: allowedParameterNames
    type     (hdf5Object     ), target      , intent(in   ), optional :: outputParametersGroup
    logical                                 , intent(in   ), optional :: noOutput
    type     (node           ), pointer                               :: parameterNode
    integer                                                           :: errorStatus

    ! Check that the file exists.
    if (.not.File_Exists(fileName)) call Galacticus_Error_Report("parameter file '"//trim(fileName)//"' does not exist"//{introspection:location})
    ! Open and parse the data file.
    !$omp critical (FoX_DOM_Access)
    parameterNode => XML_Parse(fileName,iostat=errorStatus)
    if (errorStatus /= 0) then
       if (File_Exists(fileName)) then
          call Galacticus_Error_Report('Unable to parse parameter file: "'//trim(fileName)//'"'//{introspection:location})
       else
          call Galacticus_Error_Report('Unable to find parameter file: "' //trim(fileName)//'"'//{introspection:location})
       end if
    end if
    !$omp end critical (FoX_DOM_Access)
    inputParametersConstructorFileChar=inputParametersConstructorNode(                                                 &
         &                                                            XML_Get_First_Element_By_Tag_Name(               &
         &                                                                                              parameterNode, &
         &                                                                                              'parameters'   &
         &                                                                                             )             , &
         &                                                            allowedParameterNames                          , &
         &                                                            outputParametersGroup                          , &
         &                                                            noOutput                                         &
         &                                                           )
    return
  end function inputParametersConstructorFileChar

  function inputParametersConstructorCopy(parameters)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from an existing parameters object.
    !!}
    implicit none
    type(inputParameters)                :: inputParametersConstructorCopy
    type(inputParameters), intent(in   ) :: parameters

    inputParametersConstructorCopy            =  inputParameters(parameters%rootNode  ,noOutput=.true.,noBuild=.true.)
    inputParametersConstructorCopy%parameters =>                 parameters%parameters
    inputParametersConstructorCopy%parent     =>                 parameters%parent
    return
  end function inputParametersConstructorCopy

  function inputParametersConstructorNode(parametersNode,allowedParameterNames,outputParametersGroup,noOutput,noBuild)
    !!{
    Constructor for the {\normalfont \ttfamily inputParameters} class from an FoX node.
    !!}
    use :: Display           , only : displayMessage                   , displayGreen   , displayReset
    use :: File_Utilities    , only : File_Name_Temporary
    use :: FoX_dom           , only : getOwnerDocument                 , node           , setLiveNodeLists
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)                    , char           , operator(//)    , operator(/=)
    use :: String_Handling   , only : String_Strip
    use :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Path_Exists, getTextContent => getTextContentTS
    use :: Display           , only : displayMessage
    use :: File_Utilities    , only : File_Name_Temporary
    use :: FoX_dom           , only : getOwnerDocument                 , node           , setLiveNodeLists
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)                    , char           , operator(//)    , operator(/=)
    use :: String_Handling   , only : String_Strip
    use :: IO_HDF5           , only : ioHDF5AccessInitialize
    use :: HDF5_Access       , only : hdf5Access
    implicit none
    type     (inputParameters)                                        :: inputParametersConstructorNode
    type     (node           ), pointer     , intent(in   )           :: parametersNode
    character(len=*          ), dimension(:), intent(in   ), optional :: allowedParameterNames
    type     (hdf5Object     ), target      , intent(in   ), optional :: outputParametersGroup
    logical                                 , intent(in   ), optional :: noOutput                      , noBuild
    type     (node           ), pointer                               :: versionElement
    type     (varying_string ), dimension(:), allocatable             :: allowedParameterNamesCombined , allowedParameterNamesTmp
    integer                                                           :: allowedParameterFromFileCount , allowedParameterCount
    character(len=  10       )                                        :: versionLabel
    type     (varying_string )                                        :: message
    !![
    <optionalArgument name="noOutput" defaultsTo=".false." />
    <optionalArgument name="noBuild"  defaultsTo=".false." />
    !!]
#include "os.inc"
    
    inputParametersConstructorNode%isNull   =  .false.
    inputParametersConstructorNode%document => getOwnerDocument(parametersNode)
    inputParametersConstructorNode%rootNode =>                  parametersNode
    inputParametersConstructorNode%parent   => null            (              )
    !$omp critical (FoX_DOM_Access)
    call setLiveNodeLists(inputParametersConstructorNode%document,.false.)
    if (.not.noBuild_) then
       allocate(inputParametersConstructorNode%parameters)
       inputParametersConstructorNode%parameters%content    => null()
       inputParametersConstructorNode%parameters%parent     => null()
       inputParametersConstructorNode%parameters%firstChild => null()
       inputParametersConstructorNode%parameters%sibling    => null()
       inputParametersConstructorNode%parameters%referenced => null()
       call inputParametersConstructorNode%buildTree        (inputParametersConstructorNode%parameters,parametersNode)
       call inputParametersConstructorNode%resolveReferences(                                                        )
    end if
    !$omp end critical (FoX_DOM_Access)
    ! Set a pointer to HDF5 object to which to write parameters.
    if (present(outputParametersGroup)) then
       !$ call hdf5Access%  set()
       inputParametersConstructorNode%outputParameters         =outputParametersGroup%openGroup('Parameters')
       inputParametersConstructorNode%outputParametersCopied   =.false.
       inputParametersConstructorNode%outputParametersTemporary=.false.
       !$ call hdf5Access%unset()
    else if (.not.noOutput_) then
       ! The HDF5 access lock may not yet have been initialized. Ensure it is before using it.
       call ioHDF5AccessInitialize()
       !$ call hdf5Access%  set()
       call inputParametersConstructorNode%outputParametersContainer%openFile(                                      &
            &                                                                 char(                                 &
            &                                                                      File_Name_Temporary(             &
            &                                                                                          'glcTmpPar', &
#ifdef __APPLE__
            &                                                                                          '/tmp'       &
#else
            &                                                                                          '/dev/shm'   &
#endif
            &                                                                                         )             &
            &                                                                      )                                &
            &                                                                )
       inputParametersConstructorNode%outputParameters         =inputParametersConstructorNode%outputParametersContainer%openGroup('Parameters')
       inputParametersConstructorNode%outputParametersCopied   =.false.
       inputParametersConstructorNode%outputParametersTemporary=.true.
       !$ call hdf5Access%unset()
    end if
    ! Get allowed parameter names.
    call knownParameterNames(allowedParameterNamesCombined)
    allowedParameterFromFileCount=size(allowedParameterNamesCombined)
    ! Add in parameter names explicitly listed.
    if (present(allowedParameterNames)) then
       allowedParameterCount=size(allowedParameterNames)
       if (allocated(allowedParameterNamesCombined)) then
          call Move_Alloc(allowedParameterNamesCombined,allowedParameterNamesTmp)
          allocate(allowedParameterNamesCombined(size(allowedParameterNamesTmp)+size(allowedParameterNames)))
          allowedParameterNamesCombined(                                                     &
               &                                                      1                    : &
               &                        allowedParameterFromFileCount                        &
               &                       )=allowedParameterNamesTmp
          allowedParameterNamesCombined(                                                     &
               &                        allowedParameterFromFileCount+1                    : &
               &                        allowedParameterFromFileCount+allowedParameterCount  &
               &                       )=allowedParameterNames
          deallocate(allowedParameterNamesTmp)
       else
          allocate(allowedParameterNamesCombined(size(allowedParameterNames)))
          allowedParameterNamesCombined=allowedParameterNames
       end if
    end if
    if (.not.allocated(allowedParameterNamesCombined)) allocate(allowedParameterNamesCombined(0))
    ! Check for version information.
    !$omp critical (FoX_DOM_Access)
    if (XML_Path_Exists(inputParametersConstructorNode%rootNode,"version")) then
       versionElement => XML_Get_First_Element_By_Tag_Name(inputParametersConstructorNode%rootNode,"version")
       versionLabel=getTextContent(versionElement)
       if (String_Strip(versionLabel) /= "0.9.4") then
          message=displayGreen()//"HELP:"//displayReset()                           // &
               &  " Parameter file appears to be for version "                      // &
               &  String_Strip(versionLabel)                              //char(10)// &
               &  "      Consider using: scripts/aux/parametersMigrate.pl"          // &
               &  " oldParameters.xml"                                              // &
               &  " newParameters.xml"                                    //char(10)// &
               &  "      to migrate your parameter file."
          call displayMessage(message)
       end if
    end if
    !$omp end critical (FoX_DOM_Access)
    ! Check parameters.
    call inputParametersConstructorNode%checkParameters(allowedParameterNamesCombined)
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
    Build a tree representation of the input parameter file.
    !!}
    use :: FoX_dom           , only : ELEMENT_NODE           , getNodeName  , getNodeType     , hasAttribute  , &
         &                            DOMException           , inException  , getAttributeNode, getTextContent
    use :: ISO_Varying_String, only : operator(==)           , assignment(=)
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    implicit none
    class(inputParameters), intent(inout) :: self
    type (inputParameter ), pointer       :: currentParameter, referencedParameter
    type (node           ), pointer       :: identifierNode  , identifierReferenceNode
    type (varying_string )                :: identifier      , identifierReference
    type (DOMException   )                :: exception

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
                if (inException(exception)) call Galacticus_Error_Report('unable to parse identifier'//{introspection:location})
                identifierReference     =  getTextContent  (identifierReferenceNode        ,ex=exception)
                if (inException(exception)) call Galacticus_Error_Report('unable to parse identifier'//{introspection:location})
                if (identifier == identifierReference) currentParameter%referenced => referencedParameter
             end if
             ! Walk to next node.
             referencedParameter => inputParametersWalkTree(referencedParameter)
          end do
       end if
       ! Walk to next node.
       currentParameter => inputParametersWalkTree(currentParameter)
    end do
    return
  end subroutine inputParametersResolveReferences

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
    use :: File_Utilities    , only : File_Remove
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
    !$ call hdf5Access%set()
    if (self%outputParameters%isOpen().and..not.self%outputParametersCopied) then
       if (self%outputParametersTemporary) then
          ! Close and remove the temporary parameters file.
          fileNameTemporary=self%outputParametersContainer%name()
          call self%outputParameters         %close  ()
          call self%outputParametersContainer%close  ()
          call self%outputParameters         %destroy()
          call self%outputParametersContainer%destroy()
          call File_Remove(char(fileNameTemporary))
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
       inputParameterIsParameter=                        &
            &    .not.hasAttribute(self%content,'id'   ) &
            &   .and.                                    &
            &    (                                       &
            &         hasAttribute(self%content,'value') &
            &     .or.                                   &
            &      XML_Path_Exists(self%content,'value') &
            &     .or.                                   &
            &         hasAttribute(self%content,"idRef") &
            &    )
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
    !$ use :: OMP_Lib, only : OMP_Get_Ancestor_Thread_Num, OMP_In_Parallel
    implicit none
    class  (inputParameter), intent(in   ) :: self
    integer                                :: instance

    !$omp critical (inputParameterObjects)
    if (allocated(self%objects)) then
       !$ if (OMP_In_Parallel()) then
       !$    instance=OMP_Get_Ancestor_Thread_Num(1)+1
       !$ else
             instance=                               0
       !$ end if
       inputParameterObjectCreated=associated(self%objects(instance)%object)
    else
       inputParameterObjectCreated=.false.
    end if
    !$omp end critical (inputParameterObjects)
    return
  end function inputParameterObjectCreated

  function inputParameterObjectGet(self)
    !!{
    Return a pointer to the object associated with this parameter.
    !!}
    use    :: Galacticus_Error, only : Galacticus_Error_Report
    !$ use :: OMP_Lib         , only : OMP_Get_Ancestor_Thread_Num, OMP_In_Parallel
    implicit none
    class  (functionClass ), pointer       :: inputParameterObjectGet
    class  (inputParameter), intent(in   ) :: self
    integer                                :: instance

    !$omp critical (inputParameterObjects)
    if (allocated(self%objects)) then
       !$ if (OMP_In_Parallel()) then
       !$    instance=OMP_Get_Ancestor_Thread_Num(1)+1
       !$ else
             instance=                               0
       !$ end if
       inputParameterObjectGet => self%objects(instance)%object
    else
       call Galacticus_Error_Report('object not allocated for this parameter'//{introspection:location})
    end if
    !$omp end critical (inputParameterObjects)
    return
  end function inputParameterObjectGet

  subroutine inputParameterObjectSet(self,object)
    !!{
    Set a pointer to the object associated with this parameter.
    !!}
    !$ use :: OMP_Lib, only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Max_Threads, OMP_In_Parallel
    implicit none
    class  (inputParameter), intent(inout)         :: self
    class  (functionClass ), intent(in   ), target :: object
    integer                                        :: instance

    !$omp critical (inputParameterObjects)
    !$ if (OMP_In_Parallel()) then
    !$    instance=OMP_Get_Ancestor_Thread_Num(1)+1
    !$ else
          instance=                               0
    !$ end if
    if (.not.allocated(self%objects)) then
       allocate(self%objects(0:1))
       !$ deallocate(self%objects)
       !$ allocate(self%objects(0:OMP_Get_Max_Threads()))
    end if
    self%objects(instance)%object => object
    !![
    <referenceCountIncrement owner="self%objects(instance)" object="object"/>
    !!]
    !$omp end critical (inputParameterObjects)
    return
  end subroutine inputParameterObjectSet

  recursive subroutine inputParameterReset(self,children,evaluations)
    !!{
    Reset objects associated with this parameter and any sub-parameters.
    !!}
    use :: FoX_DOM, only : destroy
    implicit none
    class  (inputParameter), intent(inout), target   :: self
    logical                , intent(in   ), optional :: children, evaluations
    type   (inputParameter), pointer                 :: child   , childNext
    integer                                          :: i
    !![
    <optionalArgument name="children"    defaultsTo=".true." />
    <optionalArgument name="evaluations" defaultsTo=".false."/>
    !!]

    ! Destroy any objects associated with this parameter node.
    if (allocated(self%objects)) then
       do i=lbound(self%objects,dim=1),ubound(self%objects,dim=1)
          if (associated(self%objects(i)%object)) then
             !![
             <objectDestructor name="self%objects(i)%object"/>
             !!]
          end if
       end do
    end if
    ! For evaluated parameters, restore them to their unevaluated state.
    if (evaluations_) then
       if (self%evaluated) then
          self%evaluated=.false.
          call self%set(self%contentOriginal)
       end if
    end if
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
    call setAttribute(self%content,"value",trim(valueText))
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

    call setAttribute(self%content,"value",char(value))
    return
  end subroutine inputParameterSetVarStr

  function inputParameterGet(self)
    !!{
    Get the value of a parameter.
    !!}
    use :: FoX_dom           , only : DOMException           , getAttributeNode, getNodeName, hasAttribute, &
          &                           inException            , node
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    use :: IO_XML            , only : getTextContent          => getTextContentTS
    use :: FoX_dom           , only : DOMException           , getAttributeNode, getNodeName, hasAttribute, &
          &                           inException            , node
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
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
       if (inException(exception)) call Galacticus_Error_Report(                                 &
            &                                                   'unable to parse parameter [' // &
            &                                                    getNodeName   (self%content) // &
            &                                                   ']'                           // &
            &                                                    {introspection:location}        &
            &                                                  )
    else
       call Galacticus_Error_Report('no parameter value present'//{introspection:location})
    end if
    !$omp end critical (FoX_DOM_Access)
    return
  end function inputParameterGet

  subroutine inputParametersCheckParameters(self,allowedParameterNames,allowedMultiParameterNames)
    use :: Galacticus_Error   , only : Galacticus_Error_Report
    use :: Display            , only : displayIndent                  , displayMessage      , displayUnindent, displayVerbosity, &
          &                            enumerationVerbosityLevelEncode, verbosityLevelSilent, displayMagenta , displayReset
    use :: FoX_dom            , only : destroy                        , getNodeName         , node           , hasAttribute    , &
         &                             getAttributeNode               , extractDataContent  , inException    , DOMException
    use :: ISO_Varying_String , only : assignment(=)                  , operator(//)        , char           , operator(==)
    use :: Regular_Expressions, only : regEx
    use :: Hashes             , only : integerHash
    use :: String_Handling    , only : String_Levenshtein_Distance
    implicit none
    class    (inputParameters)              , intent(inout)           :: self
    type     (varying_string ), dimension(:), intent(in   ), optional :: allowedParameterNames
    type     (varying_string ), dimension(:), intent(in   ), optional :: allowedMultiParameterNames
    type     (node           ), pointer                               :: node_                     , ignoreWarningsNode
    type     (inputParameter ), pointer                               :: currentParameter
    type     (regEx          ), save                                  :: regEx_
    !$omp threadprivate(regEx_)
    logical                                                           :: warningsFound             , parameterMatched    , &
         &                                                               verbose                   , ignoreWarnings      , &
         &                                                               isException
    integer                                                           :: allowedParametersCount    , errorStatus         , &
         &                                                               distance                  , distanceMinimum     , &
         &                                                               j
    character(len=1024       )                                        :: parameterValue
    character(len=1024       )                                        :: unknownName               , allowedParameterName, &
         &                                                               parameterNameGuess
    type     (varying_string )                                        :: message                   , verbosityLevel
    type     (integerHash    )                                        :: parameterNamesSeen
    type     (DOMException   )                                        :: exception

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
          if (currentParameter%isParameter()) then
             node_ => currentParameter%content
             ! Attempt to read the parameter value.
             call self%value(currentParameter,parameterValue,errorStatus,writeOutput=.false.)
             ! Determine if warnings should be ignored for this parameter.
             ignoreWarnings=.false.
             if (hasAttribute(node_,'ignoreWarnings')) then
                ignoreWarningsNode => getAttributeNode(node_,'ignoreWarnings')
                call extractDataContent(ignoreWarningsNode,ignoreWarnings,iostat=errorStatus,ex=exception)
                isException=inException(exception)
                if (isException .or. errorStatus /= 0) &
                     & call Galacticus_Error_Report("unable to parse attribute 'ignoreWarnings' in parameter ["//getNodeName(node_)//"]"//{introspection:location})
             end if
             ! Check for a match with allowed parameter names.
             allowedParametersCount=0
             if (present(allowedParameterNames)) allowedParametersCount=size(allowedParameterNames)
             if (allowedParametersCount > 0) then
                parameterMatched=.false.
                j=1
                do while (.not.parameterMatched .and. j <= allowedParametersCount)
                   allowedParameterName=allowedParameterNames(j)
                   if (allowedParameterName(1:6) == "regEx:") then
                      regEx_=regEx(allowedParameterName(7:len_trim(allowedParameterName)))
                      parameterMatched=regEx_%matches(getNodeName(node_))
                      call regEx_%destroy()
                   else
                      parameterMatched=(getNodeName(node_) == trim(allowedParameterName))
                   end if
                   j=j+1
                end do
             else
                parameterMatched=.true.
             end if
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
                select case (errorStatus)
                case (inputParameterErrorStatusEmptyValue    )
                   message='empty value for parameter ['    //getNodeName(node_)//']'
                case (inputParameterErrorStatusAmbiguousValue)
                   message='ambiguous value for parameter ['//getNodeName(node_)//']'
                end select
                !$omp end critical (FoX_DOM_Access)
                call displayMessage(message)
             end if
             if (allowedParametersCount > 0 .and. .not.parameterMatched .and. .not.ignoreWarnings .and. verbose) then
                !$omp critical (FoX_DOM_Access)
                unknownName    =getNodeName(node_)
                !$omp end critical (FoX_DOM_Access)
                distanceMinimum=-1
                do j=1,allowedParametersCount
                   allowedParameterName=allowedParameterNames(j)
                   if (allowedParameterName(1:6) == "regEx:") cycle
                   distance=String_Levenshtein_Distance(trim(unknownName),trim(allowedParameterName))
                   if (distance < distanceMinimum .or. 0 > distanceMinimum) then
                      distanceMinimum   =distance
                      parameterNameGuess=allowedParameterName
                   end if
                end do
                if (verbose) then
                   message='unrecognized parameter ['//trim(unknownName)//']'
                   if (distanceMinimum >= 0) message=message//' (did you mean ['//trim(parameterNameGuess)//']?)'
                   call displayMessage(message)
                end if
             end if
             ! Check for duplicated parameters.
             !$omp critical (FoX_DOM_Access)
             if (parameterNamesSeen%exists(getNodeName(node_))) then
                parameterMatched=.false.
                if (present(allowedMultiParameterNames)) &
                     & parameterMatched=any(getNodeName(node_) == allowedMultiParameterNames)
                if (.not.parameterMatched .and. .not.ignoreWarnings) then
                   if (.not.warningsFound.and.verbose) call displayIndent(displayMagenta()//'WARNING:'//displayReset()//' problems found with input parameters:')
                   warningsFound=.true.
                   if (verbose) then
                      message='multiple copies of parameter ['//getNodeName(node_)//'] present - only the first will be utilized'
                      call displayMessage(message)
                   end if
                end if
             else
                call parameterNamesSeen%set(getNodeName(node_),1)
             end if
             !$omp end critical (FoX_DOM_Access)
          end if
          currentParameter => currentParameter%sibling
       end do
    end if
    if (warningsFound .and. verbose) call displayUnindent('')
    return
  end subroutine inputParametersCheckParameters

  subroutine inputParametersParametersGroupOpen(self,outputGroup)
    !!{
    Open an output group for parameters in the given HDF5 object.
    !!}
    use :: File_Utilities    , only : File_Remove
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
       call File_Remove(char(fileNameTemporary))
       self%outputParameters=outputGroup%openGroup('Parameters')
       self%outputParametersTemporary=.false.
    else
       if (self%outputParameters%isOpen()) call self%outputParameters%close()
       self%outputParameters      =outputGroup%openGroup('Parameters')
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
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class    (inputParameters), intent(in   ) :: self
    character(len=*          ), intent(in   ) :: parameterName
    !$GLC attributes unused :: self

    if (trim(parameterName) == "value") call Galacticus_Error_Report('"value" is not a valid parameter name'//{introspection:location})
    return
  end subroutine inputParametersValidateName

  function inputParametersNode(self,parameterName,requireValue,copyInstance)
    !!{
    Return the node containing the parameter.
    !!}
    use :: FoX_dom         , only : ELEMENT_NODE           , getNodeName, getNodeType, hasAttribute, &
          &                         node
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: IO_XML          , only : XML_Path_Exists
    implicit none
    type     (inputParameter ), pointer                 :: inputParametersNode
    class    (inputParameters), intent(in   )           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue
    integer                   , intent(in   ), optional :: copyInstance
    type     (node           ), pointer                 :: node_
    integer                                             :: skipInstances
    !![
    <optionalArgument name="requireValue" defaultsTo=".true." />
    <optionalArgument name="copyInstance" defaultsTo="1"      />
    !!]

    call self%validateName(parameterName)
    !$omp critical (FoX_DOM_Access)
    inputParametersNode => self%parameters%firstChild
    skipInstances=copyInstance_-1
    do while (associated(inputParametersNode))
       if (.not.inputParametersNode%removed) then
          node_ => inputParametersNode%content
          if (getNodeType(node_) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(node_)) then
             if     (                                  &
                  &   .not.hasAttribute(node_,'id'   ) &
                  &  .and.                             &
                  &   (                                &
                  &    .not.requireValue_              &
                  &    .or.                            &
                  &        hasAttribute(node_,'value') &
                  &    .or.                            &
                  &     XML_Path_Exists(node_,"value") &
                  &    .or.                            &
                  &        hasAttribute(node_,"idRef") &
                  &   )                                &
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
    if (.not.associated(inputParametersNode)) call Galacticus_Error_Report('parameter node ['//trim(parameterName)//'] not found'//{introspection:location})
    if (associated(inputParametersNode%referenced)) inputParametersNode => inputParametersNode%referenced
    return
  end function inputParametersNode

  logical function inputParametersIsPresent(self,parameterName,requireValue)
    !!{
    Return true if the specified parameter is present.
    !!}
    use :: FoX_dom, only : ELEMENT_NODE   , getNodeName, getNodeType, hasAttribute, &
          &                node
    use :: IO_XML , only : XML_Path_Exists
    implicit none
    class    (inputParameters), intent(in   )           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue
    type     (node           ), pointer                 :: node_
    type     (inputParameter ), pointer                 :: currentParameter
    !![
    <optionalArgument name="requireValue" defaultsTo=".true." />
    !!]

    call self%validateName(parameterName)
    inputParametersIsPresent=.false.
    if (.not.associated(self%parameters)) return
    !$omp critical (FoX_DOM_Access)
    currentParameter => self%parameters%firstChild
    do while (associated(currentParameter))
       if (.not.currentParameter%removed) then
          node_ => currentParameter%content
          if (getNodeType(node_) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(node_)) then
             if     (                                  &
                  &   .not.hasAttribute(node_,'id'   ) &
                  &  .and.                             &
                  &   (                                &
                  &    .not.requireValue_              &
                  &    .or.                            &
                  &        hasAttribute(node_,'value') &
                  &    .or.                            &
                  &     XML_Path_Exists(node_,"value") &
                  &    .or.                            &
                  &        hasAttribute(node_,"idRef") &
                  &   )                                &
                  & ) then
                inputParametersIsPresent=.true.
                exit
             end if
          end if
       end if
       currentParameter => currentParameter%sibling
    end do
    !$omp end critical (FoX_DOM_Access)
    return
  end function inputParametersIsPresent

  integer function inputParametersCopiesCount(self,parameterName,requireValue,zeroIfNotPresent)
    !!{
    Return true if the specified parameter is present.
    !!}
    use :: FoX_dom         , only : ELEMENT_NODE           , getNodeName, getNodeType, hasAttribute, &
          &                         node
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: IO_XML          , only : XML_Path_Exists
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
          if (.not.currentParameter%removed) then
             node_ => currentParameter%content
             if (getNodeType(node_) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(node_)) then
                if     (                                    &
                     &   .not.requireValue_                 &
                     &  .or.                                &
                     &   (                                  &
                     &     .not.hasAttribute(node_,'id'   ) &
                     &    .and.                             &
                     &     (                                &
                     &          hasAttribute(node_,'value') &
                     &      .or.                            &
                     &       XML_Path_Exists(node_,"value") &
                     &      .or.                            &
                     &          hasAttribute(node_,"idRef") &
                     &     )                                &
                     &   )                                  &
                     & )                                    &
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
       call Galacticus_Error_Report('parameter ['//parameterName//'] is not present'//{introspection:location})
    end if
    return
  end function inputParametersCopiesCount

  integer function inputParametersCount(self,parameterName,zeroIfNotPresent)
    !!{
    Return a count of the number of values in a parameter.
    !!}
    use :: Galacticus_Error  , only : Galacticus_Error_Report
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
       call self%value(parameterName,parameterText,writeOutput=.false.)
       inputParametersCount=String_Count_Words(char(parameterText))
    else
       if (zeroIfNotPresent_) then
          inputParametersCount=0
       else
          inputParametersCount=0
          call Galacticus_Error_Report('parameter ['//parameterName//'] is not present'//{introspection:location})
       end if
    end if
    return
  end function inputParametersCount

  function inputParametersSubParameters(self,parameterName,requireValue,requirePresent,copyInstance)
    !!{
    Return sub-parameters of the specified parameter.
    !!}
    use :: FoX_dom           , only : node
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: ISO_Varying_String, only : assignment(=)          , char, operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    type     (inputParameters)                          :: inputParametersSubParameters
    class    (inputParameters), intent(in   ), target   :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue                , requirePresent
    integer                   , intent(in   ), optional :: copyInstance
    type     (inputParameter ), pointer                 :: parameterNode
    integer                                             :: copyCount
    type     (varying_string )                          :: groupName

    !![
    <optionalArgument name="requirePresent" defaultsTo=".true." />
    !!]

    if (.not.self%isPresent(parameterName,requireValue)) then
       if (requirePresent_) then
          call Galacticus_Error_Report('parameter ['//trim(parameterName)//'] not found'//{introspection:location})
       else
          inputParametersSubParameters=inputParameters()
       end if
       copyCount                               =  1
    else
       copyCount                               =  self%copiesCount(parameterName        ,requireValue=requireValue                          )
       parameterNode                           => self%node       (parameterName        ,requireValue=requireValue,copyInstance=copyInstance)
       inputParametersSubParameters            =  inputParameters (parameterNode%content,noOutput    =.true.      ,noBuild     =.true.      )
       inputParametersSubParameters%parameters => parameterNode
    end if
    inputParametersSubParameters%parent => self
    !$ call hdf5Access%set()
    if (self%outputParameters%isOpen()) then
       groupName=parameterName
       if (copyCount > 1) groupName=groupName//"["//copyInstance//"]"
       inputParametersSubParameters%outputParameters=self%outputParameters%openGroup(char(groupName))
    end if
    !$ call hdf5Access%unset()
    return
  end function inputParametersSubParameters

  recursive subroutine inputParametersValueName{Type¦label}(self,parameterName,parameterValue,defaultValue,errorStatus,writeOutput,copyInstance)
    !!{
    Return the value of the parameter specified by name.
    !!}
    use :: FoX_dom         , only : hasAttribute           , node
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: HDF5_Access     , only : hdf5Access
    implicit none
    class           (inputParameters), intent(inout)           :: self
    character       (len=*          ), intent(in   )           :: parameterName
    {Type¦intrinsic}                 , intent(  out)           :: parameterValue
    {Type¦intrinsic}                 , intent(in   ), optional :: defaultValue
    integer                          , intent(  out), optional :: errorStatus
    integer                          , intent(in   ), optional :: copyInstance
    logical                          , intent(in   ), optional :: writeOutput
    type            (inputParameter ), pointer                 :: parameterNode
    !![
    <optionalArgument name="writeOutput" defaultsTo=".true." />
    !!]

    if (self%isPresent(parameterName)) then
       parameterNode => self%node(parameterName,copyInstance=copyInstance)
       call self%value(parameterNode,parameterValue,errorStatus,writeOutput)
    else if (present(defaultValue)) then
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
       call Galacticus_Error_Report('parameter ['//parameterName//'] not present and no default given'//{introspection:location})
    end if
    return
  end subroutine inputParametersValueName{Type¦label}

  recursive subroutine inputParametersValueNode{Type¦label}(self,parameterNode,parameterValue,errorStatus,writeOutput)
    !!{
    Return the value of the specified parameter.
    !!}
    use :: FoX_dom           , only : DOMException                     , getAttributeNode  , getNodeName , hasAttribute, &
          &                           inException                      , node
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Path_Exists
    use :: ISO_Varying_String, only : assignment(=)                    , char              , operator(//), operator(==), &
          &                           trim
    use :: String_Handling   , only : String_Count_Words               , String_Split_Words, operator(//)
    use :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Path_Exists   , getTextContent  => getTextContentTS, extractDataContent => extractDataContentTS
    use :: FoX_dom           , only : DOMException                     , getAttributeNode  , getNodeName , hasAttribute, &
          &                           inException                      , node
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Path_Exists
    use :: ISO_Varying_String, only : assignment(=)                    , char              , operator(//), operator(==), &
          &                           trim
    use :: String_Handling   , only : String_Count_Words               , String_Split_Words, operator(//)
    implicit none
    class           (inputParameters           ), intent(inout), target      :: self
    type            (inputParameter            ), intent(inout), target      :: parameterNode
    {Type¦intrinsic}                            , intent(  out)              :: parameterValue
    integer                                     , intent(  out), optional    :: errorStatus
    logical                                     , intent(in   ), optional    :: writeOutput
#ifdef MATHEVALAVAIL
    integer         (kind_int8                 )                             :: evaluator
    ! Declarations of GNU libmatheval procedures used.
    integer         (kind_int8                 ), external                   :: Evaluator_Create_
    double precision                            , external                   :: Evaluator_Evaluate_
    external                                                                 :: Evaluator_Destroy_
#endif
    type            (inputParameter            )               , pointer     :: sibling
    type            (node                      )               , pointer     :: valueElement
    type            (inputParameters           )               , pointer     :: rootParameters     , subParameters  , &
         &                                                                      subParametersNext
    character       (len=parameterLengthMaximum), dimension(:) , allocatable :: parameterNames
    type            (DOMException              )                             :: exception
    integer                                                                  :: status             , i              , &
         &                                                                      countNames         , copyCount      , &
         &                                                                      copyInstance
    logical                                                                  :: hasValueAttribute  , hasValueElement, &
         &                                                                      isException
    character       (len=parameterLengthMaximum)                             :: expression         , parameterName  , &
         &                                                                      workText           , content
    type            (varying_string            )                             :: attributeName
    double precision                                                         :: workValue
    {Type¦match¦^Long.*¦character(len=parameterLengthMaximum) :: parameterText¦}
    {Type¦match¦^(Character|VarStr)Rank1$¦type(varying_string) :: parameterText¦}
    !![
    <optionalArgument name="writeOutput" defaultsTo=".true." />
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
          call Galacticus_Error_Report('ambiguous value attribute and element'//{introspection:location})
       end if
    else if (hasValueAttribute .or. hasValueElement) then
       !$omp critical (FoX_DOM_Access)
       if (hasValueAttribute) then
          valueElement => getAttributeNode                 (parameterNode%content,"value"                          )
       else
          valueElement => XML_Get_First_Element_By_Tag_Name(parameterNode%content,"value",directChildrenOnly=.true.)
       end if
       !$omp end critical (FoX_DOM_Access)
       if (trim(getTextContent(valueElement)) == "") then
          if (present(errorStatus)) then
             errorStatus=inputParameterErrorStatusEmptyValue
          else
             call Galacticus_Error_Report(                                       &
                  &                       'empty value in parameter ['        // &
                  &                        getNodeName(parameterNode%content) // &
                  &                       ']'                                 // &
                  &                       {introspection:location}               &
                  &                      )
          end if
       else
          ! Evaluate if an expression.
          !$omp critical (FoX_DOM_Access)
          expression=getTextContent(valueElement)
          !$omp end critical (FoX_DOM_Access)
          if (expression(1:1) == "=") then
             {Type¦match¦^Double$¦if (.true.) then¦if (.false.) then}
                ! This is an expression, and we have a scalar, floating point type - it can be evaluated.             
                !! Mark this parameter as being evaluated and store its original content. This allows the parameter to be reset to
                !! its original (unevaluated) state if necessary.
                parameterNode%evaluated      =.true.
                parameterNode%contentOriginal=expression
                !! Remove the initial "=".
                expression=expression(2:len_trim(expression))
                !! Replace other parameter values inside the parameter.
                do while (index(expression,"[") /= 0)
                   parameterName=expression(index(expression,"[")+1:index(expression,"]")-1)
                   countNames=String_Count_Words(parameterName,":")
                   allocate(parameterNames(countNames))
                   call String_Split_Words(parameterNames,parameterName,":")
                   rootParameters => self
                   do while (associated(rootParameters%parent))
                      rootParameters => rootParameters%parent
                   end do
                   do i=1,countNames-1
                      if (i == 1) then
                         allocate(subParameters)
                         subParameters    =rootParameters%subParameters(trim(parameterNames(i)),requireValue=.false.)
                      else
                         allocate(subParametersNext)
                         subParametersNext=subParameters %subParameters(trim(parameterNames(i)),requireValue=.false.)
                         deallocate(subParameters)
                         subParameters => subParametersNext
                      end if
                   end do
                   if (countNames == 1) then
                      call rootParameters%value(trim(parameterNames(countNames)),workValue)
                   else
                      call subParameters %value(trim(parameterNames(countNames)),workValue)
                      deallocate(subParameters)
                   end if
                   write (workText,'(e24.16)') workValue
                   deallocate(parameterNames)
                   expression=expression(1:index(expression,"[")-1)//trim(adjustl(workText))//expression(index(expression,"]")+1:len_trim(expression))
                end do
                !! Evaluate the expression.
#ifdef MATHEVALAVAIL
                evaluator=Evaluator_Create_(trim(expression))
                workValue=Evaluator_Evaluate_(evaluator,0,"",0.0d0)
                call Evaluator_Destroy_(evaluator)
                !$omp critical (FoX_DOM_Access)
                call parameterNode%set(workValue)
                valueElement => getAttributeNode                 (parameterNode%content,"value"                          )
                !$omp end critical (FoX_DOM_Access)
#else
                call Galacticus_Error_Report('derived parameters require libmatheval, but it is not installed'//{introspection:location})
#endif
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
                call Galacticus_Error_Report(                                &
                     &                       'unable to parse parameter ['// &
                     &                        attributeName               // &
                     &                       ']='                         // &
                     &                        content                     // &
                     &                        {introspection:location}       &
                     &                      )
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
                   sibling      => parameterNode
                   do while (associated(sibling%sibling))
                      copyInstance =  copyInstance-1
                      sibling      => sibling%sibling
                   end do
                   attributeName=attributeName//"["//copyInstance//"]"                
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

  function inputParameterListConstructor()
    !!{
    Construct an {\normalfont \ttfamily inputParameterList} object.
    !!}
    implicit none
    type(inputParameterList) :: inputParameterListConstructor

    inputParameterListConstructor%count=0
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
    type   (varying_string )                          :: inputParametersSerializeToString
    class  (inputParameters), intent(inout)           :: self
    logical                 , intent(in   ), optional :: hashed
    type   (node           ), pointer                 :: node_
    type   (inputParameter ), pointer                 :: currentParameter
    integer                                           :: errorStatus
    type   (varying_string )                          :: parameterValue
    logical                                           :: firstParameter
    type   (inputParameters)                          :: subParameters
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
          inputParametersSerializeToString=inputParametersSerializeToString// &
               &                           getNodeName(node_)              // &
               &                           ":"                             // &
               &                           adjustl(char(parameterValue))
          subParameters=self%subParameters(getNodeName(node_))
          if (associated(subParameters%parameters%firstChild)) inputParametersSerializeToString=inputParametersSerializeToString // &
               &                                                                                "{"                              // &
               &                                                                                subParameters%serializeToString()// &
               &                                                                                "}"
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

    call serialize(self%document,char(parameterFile))
    return
  end subroutine inputParametersSerializeToXML

  subroutine inputParametersAddParameter(self,parameterName,parameterValue)
    !!{
    Add a parameter to the set.
    !!}
    use :: FoX_dom, only : ELEMENT_NODE, appendChild, createElementNS, getNamespaceURI, &
          &                getNodeName , getNodeType, hasAttribute   , node           , &
          &                setAttribute
         implicit none
    class    (inputParameters), intent(inout) :: self
    character(len=*          ), intent(in   ) :: parameterName   , parameterValue
    type     (node           ), pointer       :: parameterNode   , dummy
    type     (inputParameter ), pointer       :: currentParameter
    logical                                   :: previouslyAdded
    
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

end module Input_Parameters
