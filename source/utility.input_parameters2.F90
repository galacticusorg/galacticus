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

!% Contains a module which implements reading of parameters from an XML data file.

!: $(BUILDPATH)/utility.input_parameters2.C.o

module Input_Parameters2
  !% Implements reading of parameters from an XML file.
  use Hashes_Cryptographic
  use Galacticus_Versioning
  use Galacticus_Build
  use ISO_Varying_String
  use String_Handling
  use Kind_Numbers
  use FoX_dom
  use IO_XML
  use IO_HDF5
  private
  public :: inputParameters, inputParameter, inputParameterList, globalParameters  

  ! Include public specifiers for functions that will generate unique labels for modules.
  include 'utility.input_parameters.unique_labels.visibilities.inc'

  !# <generic identifier="Type">
  !#  <instance label="Logical"        intrinsic="logical"                                        outputConverter="regEx¦(.*)¦char($1)¦"/>
  !#  <instance label="Integer"        intrinsic="integer"                                        outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="Long"           intrinsic="integer         (kind=kind_int8)"               outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="Double"         intrinsic="double precision"                               outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="Character"      intrinsic="character       (len=*         )"               outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="VarStr"         intrinsic="type            (varying_string)"               outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="LogicalRank1"   intrinsic="logical                         , dimension(:)" outputConverter="regEx¦(.*)¦char($1)¦"/>
  !#  <instance label="IntegerRank1"   intrinsic="integer                         , dimension(:)" outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="DoubleRank1"    intrinsic="double precision                , dimension(:)" outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="LongRank1"      intrinsic="integer         (kind=kind_int8), dimension(:)" outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="CharacterRank1" intrinsic="character       (len=*         ), dimension(:)" outputConverter="regEx¦(.*)¦$1¦"      />
  !#  <instance label="VarStrRank1"    intrinsic="type            (varying_string), dimension(:)" outputConverter="regEx¦(.*)¦$1¦"      />
  !# </generic>

  !# <generic identifier="cType">
  !#  <instance label="Double" intrinsic="real   (kind=c_double)" />
  !#  <instance label="Long"   intrinsic="integer(kind=c_long  )" />
  !# </generic>

  type :: genericObjectList
     !% A list-type for unlimited polymorphic pointers.
     private
     class(*), pointer :: object
  end type genericObjectList

  type :: inputParameter
     !% A class to handle input parameters for \glc.
     private
     type   (node             ), pointer                   :: content
     type   (inputParameter   ), pointer    , public       :: parent       , firstChild, sibling
     type   (genericObjectList), allocatable, dimension(:) :: objects
   contains
     !@ <objectMethods>
     !@   <object>inputParameter</object>
     !@   <objectMethod>
     !@     <method>isParameter</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if this is a valid parameter node.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>objectCreated</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true the object corresponding to this parameter has been created.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>objectGet</method>
     !@     <type>\textcolor{red}{\textless *class(*)\textgreater}</type>
     !@     <arguments></arguments>
     !@     <description>Return a pointer to the object corresponding to this parameter.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>objectSet</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless *class(*)\textgreater} object\argin, \logicalzero\ isShared\argin</arguments>
     !@     <description>Set a pointer to the object corresponding to this parameter.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: isParameter   => inputParameterIsParameter
     procedure :: objectCreated => inputParameterObjectCreated
     procedure :: objectGet     => inputParameterObjectGet
     procedure :: objectSet     => inputParameterObjectSet
  end type inputParameter

  type :: inputParameters
     private
     type   (node           ), pointer, public :: document
     type   (node           ), pointer         :: rootNode
     type   (hdf5Object     )                  :: outputParameters
     type   (inputParameter ), pointer         :: parameters
     type   (inputParameters), pointer, public :: parent
     logical                                   :: global
   contains
     !@ <objectMethods>
     !@   <object>inputParameters</object>
     !@   <objectMethod>
     !@     <method>buildTree</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(inputParameter)\textgreater} *parentParameter\arginout, \textcolor{red}{\textless type(node)\textgreater} *parametersNode\argin</arguments>
     !@     <description>Build a tree of {\normalfont \ttfamily inputParameter} objects from the structure of an XML parameter file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>markGlobal</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Mark an input parameter set as the global set.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isGlobal</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if these parameters are the global set.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>parametersGroupOpen</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(hdf5Object)\textgreater} outputGroup\arginout</arguments>
     !@     <description>Open an output group for parameters in the given HDF5 object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>validateName</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} parameterName\argin</arguments>
     !@     <description>Check that a given parameter name is a valid name, aborting if not.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>checkParameters</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(varying\_string)[:]\textgreater} [allowedParameterNames]\argin</arguments>
     !@     <description>Check that parameters are valid and, optionally, check if they match expected names in the provided list.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>node</method>
     !@     <type>\textcolor{red}{\textless type(node)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} parameterName\argin</arguments>
     !@     <description>Return the XML node containing the named parameter.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isPresent</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} parameterName\argin</arguments>
     !@     <description>Return true if the named parameter is present in the set.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>copiesCount</method>
     !@     <type>\intzero</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} parameterName\argin, \logicalzero [zeroIfNotPresent]\argin</arguments>
     !@     <description>Return a count of the number copies of the named parameter. If the parameter is not present, this function aborts, unless {\normalfont \ttfamily zeroIfNotPresent} is set to {\normalfont \ttfamily true}, in which case a result of 0 is returned.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>count</method>
     !@     <type>\intzero</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} parameterName\argin, \logicalzero [zeroIfNotPresent]\argin</arguments>
     !@     <description>Return a count of the number of values in the named parameter. If the parameter is not present, this function aborts, unless {\normalfont \ttfamily zeroIfNotPresent} is set to {\normalfont \ttfamily true}, in which case a result of 0 is returned.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>subParameters</method>
     !@     <type>\textcolor{red}{\textless type(inputParameters)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} parameterName\argin</arguments>
     !@     <description>Return the set of subparameters of the named parameter.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>value</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}({\textless *\textgreater} parameterName|{\textless node\textgreater} parameterNode)\argin, parameterValue, \textcolor{red}({\textless *\textgreater} [defaultValue]|)\argin, \enumInputParameterErrorStatus [errorStatus]\argout, \logicalzero [writeOutput]\argin</arguments>
     !@     <description>Return the value of a parameter specified by name or XML node. A default value can be specified only if the parameter is specified by name. Supported types include rank-0 and rank-1 logicals, integers, long integers, doubles, characters, and varying strings.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serializeToString</method>
     !@     <type>\textcolor{red}{\textless type(varying\_string)\textgreater}</type>
     !@     <arguments></arguments>
     !@     <description>Serialize input parameters to a string.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serializeToXML</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(varying\_string)\textgreater} parameterFile\argin</arguments>
     !@     <description>Serialize input parameters to an XML file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>addParameter</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} parameterName\argin, \textcolor{red}{\textless character(len=*)\textgreater} parameterValue\argin</arguments>
     !@     <description>Add a parameter.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Destroy the parameters document.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                        inputParametersFinalize
     procedure :: buildTree           => inputParametersBuildTree
     procedure :: destroy             => inputParametersDestroy
     procedure :: markGlobal          => inputParametersMarkGlobal
     procedure :: isGlobal            => inputParametersIsGlobal
     procedure :: parametersGroupOpen => inputParametersParametersGroupOpen
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
  end type inputParameters

  interface inputParameters
     !% Constructors for the {\normalfont \ttfamily inputParameters} class.
     module procedure inputParametersConstructorFileChar
     module procedure inputParametersConstructorFileVarStr
     module procedure inputParametersConstructorNode
     module procedure inputParametersConstructorNull
  end interface inputParameters

  ! Define a type to hold lists of parameters (and values) prior to output.
  type :: inputParameterList
     !% A class to hold lists of parameters (and values) prior to output.
     integer                                            :: count
     type   (varying_string), allocatable, dimension(:) :: name , value
   contains
     !@ <objectMethods>
     !@   <object>inputParameterList</object>
     !@   <objectMethod>
     !@     <method>serializeToXML</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Serialize a list of input parameters to an XML document.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>add</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} name\argin, \textcolor{red}{\textless character(len=*)\textgreater} value\argin</arguments>
     !@     <description>Add a parameter and value to the list.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                      inputParameterListDestructor    
     procedure :: add               => inputParameterListAdd
     procedure :: serializeToXML    => inputParameterListSerializeToXML
  end type inputParameterList

  interface inputParameterList
     !% Constructors for {\normalfont \ttfamily inputParameterList} objects.
     module procedure inputParameterListConstructor
  end interface inputParameterList
  
  !# <enumeration>
  !#  <name>inputParameterErrorStatus</name>
  !#  <description>Error status codes used by the input parameters module.</description>
  !#  <entry label="success"        />
  !#  <entry label="notPresent"     />
  !#  <entry label="parse"          />
  !#  <entry label="emptyValue"     />
  !#  <entry label="ambiguousValue" />
  !# </enumeration>

  ! Pointer to the global input parameters.
  type   (inputParameters), pointer   :: globalParameters               => null()

  ! Thread-global state recording whether objects being built are threadprivate or public.
  logical                 , public    :: parametersObjectBuildIsPrivate  
  !$omp threadprivate(parametersObjectBuildIsPrivate)
  
  ! Maximum length allowed for parameter entries.
  integer                 , parameter :: parameterLengthMaximum         =  1024
  
contains

  function inputParametersConstructorNull()
    !% Constructor for the {\normalfont \ttfamily inputParameters} class creating a null instance.
    implicit none
    type(inputParameters) :: inputParametersConstructorNull

    inputParametersConstructorNull%document   => createDocument    (                                  &
         &                                                          getImplementation()             , &
         &                                                          qualifiedName      ="parameters", &
         &                                                          docType            =null()        &
         &                                                         )
    inputParametersConstructorNull%rootNode   => getDocumentElement(inputParametersConstructorNull%document)
    inputParametersConstructorNull%parameters => null()
    inputParametersConstructorNull%global     = .false.
    call setLiveNodeLists(inputParametersConstructorNull%document,.true.)
    return
  end function inputParametersConstructorNull
  
  
  function inputParametersConstructorFileVarStr(fileName,allowedParameterNames,allowedParametersFile,outputParametersGroup)
    !% Constructor for the {\normalfont \ttfamily inputParameters} class from an XML file
    !% specified as a variable length string.
    implicit none
    type     (inputParameters)                                           :: inputParametersConstructorFileVarStr
    type     (varying_string    )              , intent(in   )           :: fileName
    character(len=*             ), dimension(:), intent(in   ), optional :: allowedParameterNames
    character(len=*             )              , intent(in   ), optional :: allowedParametersFile
    type     (hdf5Object        ), target      , intent(in   ), optional :: outputParametersGroup
    
    inputParametersConstructorFileVarStr=inputParametersConstructorFileChar(                       &
         &                                                                  char(fileName)       , &
         &                                                                  allowedParameterNames, &
         &                                                                  allowedParametersFile, &
         &                                                                  outputParametersGroup  &
         &                                                                 )
    return
  end function inputParametersConstructorFileVarStr
  
  function inputParametersConstructorFileChar(fileName,allowedParameterNames,allowedParametersFile,outputParametersGroup)
    !% Constructor for the {\normalfont \ttfamily inputParameters} class from an XML file
    !% specified as a character variable.
    use Galacticus_Error
    implicit none
    type     (inputParameters)                                        :: inputParametersConstructorFileChar
    character(len=*          )              , intent(in   )           :: fileName
    character(len=*          ), dimension(:), intent(in   ), optional :: allowedParameterNames
    character(len=*          )              , intent(in   ), optional :: allowedParametersFile
    type     (hdf5Object     ), target      , intent(in   ), optional :: outputParametersGroup
    type     (node           ), pointer                               :: parameterNode
    integer                                                           :: errorStatus
    
    ! Open and parse the data file.
    !$omp critical (FoX_DOM_Access)
    parameterNode => parseFile(fileName,iostat=errorStatus)
    if (errorStatus /= 0) call Galacticus_Error_Report('inputParametersConstructorFileChar','Unable to find or parse parameter file')
    !$omp end critical (FoX_DOM_Access)
    inputParametersConstructorFileChar=inputParametersConstructorNode(                                                 &
         &                                                            XML_Get_First_Element_By_Tag_Name(               &
         &                                                                                              parameterNode, &
         &                                                                                              'parameters'   &
         &                                                                                             )             , &
         &                                                            allowedParameterNames                          , &
         &                                                            allowedParametersFile                          , &
         &                                                            outputParametersGroup                            &
         &                                                           )
    return
  end function inputParametersConstructorFileChar

  function inputParametersConstructorNode(parametersNode,allowedParameterNames,allowedParametersFile,outputParametersGroup)
    !% Constructor for the {\normalfont \ttfamily inputParameters} class from an FoX node.
    use Galacticus_Error
    use Galacticus_Display
    use Galacticus_Input_Paths
    use File_Utilities
    implicit none
    type     (inputParameters)                                        :: inputParametersConstructorNode
    type     (node           ), pointer     , intent(in   )           :: parametersNode
    character(len=*          ), dimension(:), intent(in   ), optional :: allowedParameterNames
    character(len=*          )              , intent(in   ), optional :: allowedParametersFile
    type     (hdf5Object     ), target      , intent(in   ), optional :: outputParametersGroup
    type     (node           ), pointer                               :: thisNode                      , allowedParameterDoc    , &
         &                                                               versionElement
    type     (nodeList       ), pointer                               :: allowedParameterList
    type     (varying_string ), dimension(:), allocatable             :: allowedParameterNamesCombined, allowedParameterNamesTmp
    integer                                                           :: i                            , errorStatus             , &
         &                                                               allowedParameterFromFileCount, allowedParameterCount
    character(len=  10       )                                        :: versionLabel
    type     (varying_string )                                        :: message

    inputParametersConstructorNode%global   =  .false.
    inputParametersConstructorNode%document => getOwnerDocument(parametersNode)    
    inputParametersConstructorNode%rootNode =>                  parametersNode
    inputParametersConstructorNode%parent   => null            (              )
    !$omp critical (FoX_DOM_Access)
    allocate(inputParametersConstructorNode%parameters)
    inputParametersConstructorNode%parameters%firstChild => null()
    call inputParametersConstructorNode%buildTree(inputParametersConstructorNode%parameters,parametersNode)    
    !$omp end critical (FoX_DOM_Access)
    ! Set a pointer to HDF5 object to which to write parameters.
    if (present(outputParametersGroup)) inputParametersConstructorNode%outputParameters=outputParametersGroup%openGroup('Parameters')
    ! Parse allowed parameters file if available.
    allowedParameterFromFileCount=0
    if (present(allowedParametersFile)) then
       ! Check if the file exists.
       if (File_Exists(char(Galacticus_Input_Path())//BUILDPATH//'/'//allowedParametersFile)) then
          !$omp critical (FoX_DOM_Access)
          ! Parse the file.
          allowedParameterDoc => parseFile(char(Galacticus_Input_Path())//BUILDPATH//'/'//allowedParametersFile,iostat=errorStatus)
          if (errorStatus /= 0) call Galacticus_Error_Report('inputParametersConstructorNode','Unable to parse allowed parameters file')
          ! Extract allowed parameter names to array.
          allowedParameterList => getElementsByTagname(allowedParameterDoc,"parameter")
          allowedParameterFromFileCount=getLength(allowedParameterList)
          allocate(allowedParameterNamesCombined(allowedParameterFromFileCount+allowedParameterCount))
          do i=0,allowedParameterFromFileCount-1
             thisNode => item(allowedParameterList,i)
             allowedParameterNamesCombined(i+1)=getTextContent(thisNode)
          end do
          ! Destroy the allowed parameter names document.
          call destroy(allowedParameterDoc)
          !$omp end critical (FoX_DOM_Access)
      else
          call Galacticus_Display_Message("Allowed parameter file '"//allowedParametersFile//"' is missing - incorrect parameters will not be detected",verbosityWarn)
       end if
    end if
    ! Add in parameter names explicitly listed.
    if (present(allowedParameterNames)) then
       if (allocated(allowedParameterNamesCombined)) then
          call Move_Alloc(allowedParameterNamesCombined,allowedParameterNamesTmp)
          allocate(allowedParameterNamesCombined(size(allowedParameterNamesTmp)+size(allowedParameterNames)))
          allowedParameterNamesCombined(                                                     &
               &                                                      1                     :&
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
          message="HELP: Parameter file appears to be for version "                 // &
               &  String_Strip(versionLabel)                              //char(10)// &
               &  "      Consider using: scripts/aux/parametersMigrate.pl"          // &
               &  " oldParameters.xml"                                              // &
               &  " newParameters.xml"                                    //char(10)// &
               &  "      to migrate your parameter file."
          call Galacticus_Display_Message(message)
       end if
    end if
    !$omp end critical (FoX_DOM_Access)
    ! Check parameters.
    call inputParametersConstructorNode%checkParameters(allowedParameterNamesCombined)
    return
  end function inputParametersConstructorNode
  
  recursive subroutine inputParametersBuildTree(self,parentParameter,parametersNode)
    !% Build a tree representation of the input parameter file.
    implicit none
    class  (inputParameters), intent(inout)          :: self
    type   (inputParameter ), intent(inout), pointer :: parentParameter
    type   (node           ), intent(in   ), pointer :: parametersNode
    type   (nodeList       )               , pointer :: childNodes
    type   (node           )               , pointer :: childNode
    type   (inputParameter )               , pointer :: currentParameter
    integer                                          :: i
    
    childNodes       => getChildNodes(parametersNode)
    currentParameter => null()
    do i=0,getLength(childNodes)-1
       childNode => item(childNodes,i)
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
          call self%buildTree(currentParameter,childNode)
       end if
    end do
    return
  end subroutine inputParametersBuildTree

  subroutine inputParametersDestroy(self)
    !% Destructor for the {\normalfont \ttfamily inputParameters} class.
    class(inputParameters), intent(inout) :: self

    ! Destroy the parameters document. Note that we do not use a finalizer for input parameters. This could destroy part of a
    ! document which was still pointed to from elsewhere, leaving a dangling pointer. Instead, destruction only occurs when
    ! explicitly requested.
    !$omp critical (FoX_DOM_Access)
    call destroy(self%document)
    !$omp end critical (FoX_DOM_Access)
    call inputParametersFinalize(self)
    return
  end subroutine inputParametersDestroy
 
  subroutine inputParametersFinalize(self)
    !% Finalizer for the {\normalfont \ttfamily inputParameters} class.
    type(inputParameters), intent(inout) :: self

    nullify(self%document  )
    nullify(self%rootNode  )
    nullify(self%parameters)
    nullify(self%parent    )
    !$omp critical(HDF5_Access)
    if (self%outputParameters%isOpen()) call self%outputParameters%close()
    !$omp end critical(HDF5_Access)
    return
  end subroutine inputParametersFinalize
  
  logical function inputParameterIsParameter(self)
    !% Return true if this is a valid parameter.
    implicit none
    class(inputParameter), intent(in   )         :: self
    
    !$omp critical (FoX_DOM_Access)
    if (associated(self%content) .and. getNodeType(self%content) == ELEMENT_NODE) then
       inputParameterIsParameter=                &
            &   hasAttribute   (self%content,'value') &
            &  .or.                               &
            &   XML_Path_Exists(self%content,'value') 
    else
       inputParameterIsParameter=.false.
    end if
    !$omp end critical (FoX_DOM_Access)
    return
  end function inputParameterIsParameter
  
  logical function inputParameterObjectCreated(self)
    !% Return true if the specified instance of the object associated with this parameter has been created.
    !$ use OMP_Lib
    implicit none
    class  (inputParameter), intent(in   ) :: self
    integer                                :: instance
    
    !$omp critical (inputParameterObjects)
    if (allocated(self%objects)) then
       !$ if (parametersObjectBuildIsPrivate) then
       !$    instance=OMP_Get_Ancestor_Thread_Num(0)+1
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
    !% Return a pointer to the object associated with this parameter.
    !$ use OMP_Lib
    use Galacticus_Error
    implicit none
    class  (*             ), pointer       :: inputParameterObjectGet
    class  (inputParameter), intent(in   ) :: self
    integer                                :: instance
    
    !$omp critical (inputParameterObjects)
    if (allocated(self%objects)) then
       !$ if (parametersObjectBuildIsPrivate) then
       !$    instance=OMP_Get_Ancestor_Thread_Num(0)+1
       !$ else
             instance=                               0
       !$ end if
       inputParameterObjectGet => self%objects(instance)%object
    else
       call Galacticus_Error_Report('inputParameterObjectGet','object not allocated for this parameter')
    end if
    !$omp end critical (inputParameterObjects)
    return
  end function inputParameterObjectGet
  
  subroutine inputParameterObjectSet(self,object)
    !% Set a pointer to the object associated with this parameter.
    !$ use OMP_Lib
    implicit none
    class  (inputParameter), intent(inout)         :: self
    class  (*             ), intent(in   ), target :: object
    integer                                        :: instance
    
    !$omp critical (inputParameterObjects)
    !$ if (parametersObjectBuildIsPrivate) then
    !$    instance=OMP_Get_Ancestor_Thread_Num(0)+1
    !$ else
          instance=                               0 
    !$ end if
    if (.not.allocated(self%objects)) allocate(self%objects(0:OMP_Get_Max_Threads()))
    self%objects(instance)%object => object
    !$omp end critical (inputParameterObjects)
    return
  end subroutine inputParameterObjectSet
  
  subroutine inputParametersCheckParameters(self,allowedParameterNames)
    !$ use OMP_Lib
    use Regular_Expressions
    use String_Handling
    implicit none
    class    (inputParameters)              , intent(inout)           :: self
    type     (varying_string ), dimension(:), intent(in   ), optional :: allowedParameterNames
    type     (node           ), pointer                               :: thisNode
    type     (inputParameter ), pointer                               :: currentParameter
    type     (regEx          ), save                                  :: thisRegEx
    !$omp threadprivate(thisRegEx)
    logical                                                           :: warningsFound                , parameterMatched
    integer                                                           :: i                            , j                   , &
         &                                                               allowedParametersCount       , errorStatus         , &
         &                                                               distance                     , distanceMinimum     , &
         &                                                               allowedParameterFromFileCount, nodeCount
    character(len=1024       )                                        :: parameterValue
    character(len=1024       )                                        :: unknownName                  , allowedParameterName, &
         &                                                               parameterNameGuess

    
    ! Return if there are no parameters to check.
    ! Validate parameters.
    warningsFound=.false.
    if (associated(self%parameters)) then
       currentParameter => self%parameters%firstChild
       do while (associated(currentParameter))
          if (currentParameter%isParameter()) then
             thisNode => currentParameter%content
             ! Attempt to read the parameter value.
             call self%value(currentParameter,parameterValue,errorStatus,writeOutput=.false.)
             ! Check for a match with allowed parameter names.
             allowedParametersCount=0
             if (present(allowedParameterNames)) allowedParametersCount=size(allowedParameterNames)
             if (allowedParametersCount > 0) then
                parameterMatched=.false.
                j=1
                do while (.not.parameterMatched .and. j <= allowedParametersCount)
                   allowedParameterName=allowedParameterNames(j)
                   if (allowedParameterName(1:6) == "regEx:") then
                      thisRegEx=regEx(allowedParameterName(7:len_trim(allowedParameterName)))
                      parameterMatched=thisRegEx%matches(getNodeName(thisNode))
                      call thisRegEx%destroy()
                   else
                      parameterMatched=(getNodeName(thisNode) == trim(allowedParameterName))
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
                  &   .not.warningsFound                                &
                  & ) then
                !$ if (omp_in_parallel()) then
                !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
                !$ else
                !$    write (0,'(a2,a2,$)') "MM",": "
                !$ end if
                write (0,'(a)') '-> WARNING: problems found with input parameters:'
                warningsFound=.true.
             end if
             if (errorStatus /= inputParameterErrorStatusSuccess) then
                !$ if (omp_in_parallel()) then
                !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
                !$ else
                !$    write (0,'(a2,a2,$)') "MM",": "
                !$ end if
                !$omp critical (FoX_DOM_Access)
                select case (errorStatus)
                case (inputParameterErrorStatusEmptyValue    )
                   write (0,'(3a)') '    empty value for parameter ['    ,getNodeName(thisNode),']'
                case (inputParameterErrorStatusAmbiguousValue)
                   write (0,'(3a)') '    ambiguous value for parameter [',getNodeName(thisNode),']'
                end select
                !$omp end critical (FoX_DOM_Access)
             end if
             if (allowedParametersCount+allowedParameterFromFileCount > 0 .and. .not.parameterMatched) then
                !$ if (omp_in_parallel()) then
                !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
                !$ else
                !$    write (0,'(a2,a2,$)') "MM",": "
                !$ end if
                !$omp critical (FoX_DOM_Access)
                unknownName    =getNodeName(thisNode)
                !$omp end critical (FoX_DOM_Access)
                distanceMinimum=-1
                do j=1,allowedParametersCount+allowedParameterFromFileCount
                   allowedParameterName=allowedParameterNames(j)
                   if (allowedParameterName(1:6) == "regEx:") cycle
                   distance=String_Levenshtein_Distance(trim(unknownName),trim(allowedParameterName))
                   if (distance < distanceMinimum .or. 0 > distanceMinimum) then
                      distanceMinimum   =distance
                      parameterNameGuess=allowedParameterName
                   end if
                end do
                if (distanceMinimum < 0) then
                   write (0,'(3a)') '    unrecognised parameter [',trim(unknownName),']'
                else
                   write (0,'(5a)') '    unrecognised parameter [',trim(unknownName),'] (did you mean [',trim(parameterNameGuess),']?)'
                end if
             end if
          end if
          currentParameter => currentParameter%sibling   
       end do
    end if
    if (warningsFound) then
       !$ if (omp_in_parallel()) then
       !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
       !$ else
       !$    write (0,'(a2,a2,$)') "MM",": "
       !$ end if
       write (0,'(a)') '<-'
    end if
    return
  end subroutine inputParametersCheckParameters  

  subroutine inputParametersMarkGlobal(self)
    !% Mark an {\normalfont \ttfamily inputParameters} object as the global input parameters.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    implicit none
    class(inputParameters), intent(inout), target :: self
    type (c_ptr          )                        :: globalParametersC
    interface
       subroutine inputParametersSetGlobalC(globalParametersC) bind(c,name="inputParametersSetGlobalC")
         import c_ptr
         type(c_ptr), value :: globalParametersC
       end subroutine inputParametersSetGlobalC
    end interface

    ! Check that global parameters have not yet been assigned.
    if (associated(globalParameters)) call Galacticus_Error_Report('inputParametersMarkGlobal','global parameters cannot be reassigned')
    ! Mark as global.
    self%global=.true.
    ! Set global parameters pointer.
    globalParameters => self
    ! Set global parameters pointer in the C interface also.
    globalParametersC=c_loc(globalParameters)
    call inputParametersSetGlobalC(globalParametersC)
    return
  end subroutine inputParametersMarkGlobal
  
  logical function inputParametersIsGlobal(self)
    !% Return true if an {\normalfont \ttfamily inputParameters} object is the global input parameter set.
    implicit none
    class(inputParameters), intent(in   ), target :: self

    inputParametersIsGlobal=self%global
    return
  end function inputParametersIsGlobal

  subroutine inputParametersParametersGroupOpen(self,outputGroup)
    !% Open an output group for parameters in the given HDF5 object.
    implicit none
    class(inputParameters), intent(inout) :: self
    type (hdf5Object     ), intent(inout) :: outputGroup

    !$omp critical(HDF5_Access)
    if (self%outputParameters%isOpen()) call self%outputParameters%close()
    self%outputParameters=outputGroup%openGroup('Parameters')
    !$omp end critical(HDF5_Access)
    return
  end subroutine inputParametersParametersGroupOpen

  subroutine inputParametersValidateName(self,parameterName)
    !% Validate a parameter name.
    use Galacticus_Error
    implicit none
    class    (inputParameters), intent(in   ) :: self
    character(len=*          ), intent(in   ) :: parameterName

    if (trim(parameterName) == "value") call Galacticus_Error_Report('inputParametersValidateName','"value" is not a valid parameter name')
    return
  end subroutine inputParametersValidateName

  function inputParametersNode(self,parameterName,requireValue,copyInstance)
    !% Return the node containing the parameter.
    use Galacticus_Error
    implicit none
    type     (inputParameter ), pointer                 :: inputParametersNode
    class    (inputParameters), intent(in   )           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue
    integer                   , intent(in   ), optional :: copyInstance
    type     (node           ), pointer                 :: thisNode
    integer                                             :: i                  , skipInstances
    !# <optionalArgument name="requireValue" defaultsTo=".true." />
    !# <optionalArgument name="copyInstance" defaultsTo="1"      />

    call self%validateName(parameterName)
    !$omp critical (FoX_DOM_Access)
    inputParametersNode => self%parameters%firstChild
    skipInstances=copyInstance_-1
    do while (associated(inputParametersNode))
       thisNode => inputParametersNode%content
       if (getNodeType(thisNode) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(thisNode)) then
          if (.not.requireValue_ .or. hasAttribute(thisNode,'value') .or. XML_Path_Exists(thisNode,"value")) then
             if (skipInstances > 0) then
                skipInstances=skipInstances-1
             else
                exit
             end if
          end if
       end if
       inputParametersNode => inputParametersNode%sibling
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.associated(inputParametersNode)) call Galacticus_Error_Report('inputParametersNode','parameter node ['//trim(parameterName)//'] not found')
    return
  end function inputParametersNode
  
  logical function inputParametersIsPresent(self,parameterName,requireValue)
    !% Return true if the specified parameter is present.
    implicit none
    class    (inputParameters), intent(in   )           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue
    type     (node           ), pointer                 :: thisNode
    type     (inputParameter ), pointer                 :: currentParameter
    integer                                             :: i
    !# <optionalArgument name="requireValue" defaultsTo=".true." />

    call self%validateName(parameterName)
    inputParametersIsPresent=.false.
    if (.not.associated(self%parameters)) return
    !$omp critical (FoX_DOM_Access)
    currentParameter => self%parameters%firstChild
    do while (associated(currentParameter))
       thisNode => currentParameter%content
       if (getNodeType(thisNode) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(thisNode)) then
          if (.not.requireValue_ .or. hasAttribute(thisNode,'value') .or. XML_Path_Exists(thisNode,"value")) then
             inputParametersIsPresent=.true.
             exit
          end if
       end if
       currentParameter => currentParameter%sibling
    end do
    !$omp end critical (FoX_DOM_Access)
    return
  end function inputParametersIsPresent

  integer function inputParametersCopiesCount(self,parameterName,zeroIfNotPresent)
    !% Return true if the specified parameter is present.
    use Galacticus_Error
    implicit none
    class    (inputParameters), intent(in   )           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: zeroIfNotPresent
    type     (node           ), pointer                 :: thisNode
    type     (inputParameter ), pointer                 :: currentParameter
    integer                                             :: i
    !# <optionalArgument name="zeroIfNotPresent" defaultsTo=".false." />
    
    call self%validateName(parameterName)
    if (self%isPresent(parameterName)) then
       inputParametersCopiesCount=0
       !$omp critical (FoX_DOM_Access)
       currentParameter => self%parameters%firstChild
       do while (associated(currentParameter))
          thisNode => currentParameter%content
          if (getNodeType(thisNode) == ELEMENT_NODE .and. trim(parameterName) == getNodeName(thisNode)) then
             if (hasAttribute(thisNode,'value') .or. XML_Path_Exists(thisNode,"value")) &
                  & inputParametersCopiesCount=inputParametersCopiesCount+1
          end if
          currentParameter => currentParameter%sibling
      end do
       !$omp end critical (FoX_DOM_Access)
    else if (zeroIfNotPresent_) then
       inputParametersCopiesCount=0
    else
       call Galacticus_Error_Report('inputParametersCopiesCount','parameter ['//parameterName//'] is not present')  
    end if
    return
  end function inputParametersCopiesCount

  integer function inputParametersCount(self,parameterName,zeroIfNotPresent)
    !% Return a count of the number of values in a parameter.
    use Galacticus_Error
    implicit none
    class    (inputParameters), intent(inout)           :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: zeroIfNotPresent
    type     (varying_string )                          :: parameterText
    !# <optionalArgument name="zeroIfNotPresent" defaultsTo=".false." />

    if (self%isPresent(parameterName)) then
       call self%value(parameterName,parameterText,writeOutput=.false.)    
       inputParametersCount=String_Count_Words(char(parameterText))
    else
       if (zeroIfNotPresent_) then
          inputParametersCount=0
       else
          call Galacticus_Error_Report('inputParametersCount','parameter ['//parameterName//'] is not present')
       end if
    end if
    return
  end function inputParametersCount
  
  function inputParametersSubParameters(self,parameterName,requireValue,requirePresent,copyInstance)
    !% Return sub-parameters of the specified parameter.
    use Galacticus_Error
    implicit none
    type     (inputParameters)                          :: inputParametersSubParameters
    class    (inputParameters), intent(in   ), target   :: self
    character(len=*          ), intent(in   )           :: parameterName
    logical                   , intent(in   ), optional :: requireValue                , requirePresent
    integer                   , intent(in   ), optional :: copyInstance
    type     (inputParameter ), pointer                 :: parameterNode
    !# <optionalArgument name="requirePresent" defaultsTo=".true." />

    if (.not.self%isPresent(parameterName,requireValue)) then
       if (requirePresent_) then
          call Galacticus_Error_Report('inputParametersSubParameters','parameter not found')
       else
          inputParametersSubParameters=inputParameters()
       end if
    else
       parameterNode                => self%node      (parameterName        ,requireValue=requireValue,copyInstance=copyInstance)
       inputParametersSubParameters =  inputParameters(parameterNode%content                                                    )
    end if
    inputParametersSubParameters%parent => self
    inputParametersSubParameters%global =  self%global
    !$omp critical(HDF5_Access)
    if (self%outputParameters%isOpen()) inputParametersSubParameters%outputParameters=self%outputParameters%openGroup(trim(parameterName))
    !$omp end critical(HDF5_Access)
    return
  end function inputParametersSubParameters

  subroutine inputParametersValueName{Type¦label}(self,parameterName,parameterValue,defaultValue,errorStatus,writeOutput,copyInstance)
    !% Return the value of the parameter specified by name.
    use Galacticus_Error
    implicit none
    class           (inputParameters), intent(inout)           :: self
    character       (len=*          ), intent(in   )           :: parameterName
    {Type¦intrinsic}                 , intent(  out)           :: parameterValue
    {Type¦intrinsic}                 , intent(in   ), optional :: defaultValue
    integer                          , intent(  out), optional :: errorStatus
    integer                          , intent(in   ), optional :: copyInstance
    logical                          , intent(in   ), optional :: writeOutput
    type            (inputParameter ), pointer                 :: parameterNode
    !# <optionalArgument name="writeOutput" defaultsTo=".true." />

    if (self%isPresent(parameterName)) then
       parameterNode => self%node(parameterName,copyInstance=copyInstance)
       call self%value(parameterNode,parameterValue,errorStatus,writeOutput)
    else if (present(defaultValue)) then
       parameterValue=defaultValue
       ! Write the parameter file to an HDF5 object.
       if (self%outputParameters%isOpen().and.writeOutput_) then
          !$omp critical(HDF5_Access)
         if (.not.self%outputParameters%hasAttribute(parameterName)) call self%outputParameters%writeAttribute({Type¦outputConverter¦parameterValue},parameterName)
          !$omp end critical(HDF5_Access)
       end if
    else if (present(errorStatus )) then
       errorStatus   =inputParameterErrorStatusNotPresent
    else
       call Galacticus_Error_Report('inputParametersValueName{Type¦label}','parameter ['//parameterName//'] not present and no default given')
    end if
    return
  end subroutine inputParametersValueName{Type¦label}

  subroutine inputParametersValueNode{Type¦label}(self,parameterNode,parameterValue,errorStatus,writeOutput)
    !% Return the value of the specified parameter.
    use Galacticus_Error
    implicit none
    class           (inputParameters), intent(inout)           :: self
    type            (inputParameter ), intent(in   )           :: parameterNode
    {Type¦intrinsic}                 , intent(  out)           :: parameterValue
    integer                          , intent(  out), optional :: errorStatus
    logical                          , intent(in   ), optional :: writeOutput
    type            (node           )               , pointer  :: valueElement
    type            (DOMException   )                          :: exception
    integer                                                    :: status
    logical                                                    :: hasValueAttribute, hasValueElement
    {Type¦match¦^Long.*¦character(len=parameterLengthMaximum) :: parameterText¦}    
    {Type¦match¦^(Character|VarStr)Rank1$¦type(varying_string) :: parameterText¦}    
    !# <optionalArgument name="writeOutput" defaultsTo=".true." />

    if (present(errorStatus)) errorStatus=inputParameterErrorStatusSuccess
    !$omp critical (FoX_DOM_Access)
    hasValueAttribute =  hasAttribute   (parameterNode%content,'value')
    hasValueElement   =  XML_Path_Exists(parameterNode%content,'value')
    if (hasValueAttribute .and. hasValueElement) then
       if (present(errorStatus)) then
          errorStatus=inputParameterErrorStatusAmbiguousValue
       else
          call Galacticus_Error_Report('inputParametersValueNode{Type¦label}','ambiguous value attribute and element')
       end if
    else if (hasValueAttribute .or. hasValueElement) then
       if (hasValueAttribute) then
          valueElement => getAttributeNode                 (parameterNode%content,"value"                          )
       else
          valueElement => XML_Get_First_Element_By_Tag_Name(parameterNode%content,"value",directChildrenOnly=.true.)
       end if
       if (trim(getTextContent(valueElement)) == "") then
          if (present(errorStatus)) then
             errorStatus=inputParameterErrorStatusEmptyValue
          else
             call Galacticus_Error_Report(                                          &
                  &                       'inputParametersValueNode{Type¦label}' ,  &
                  &                       'empty value in parameter ['           // &
                  &                        getNodeName(parameterNode%content)            // &
                  &                       ']'                                       &
                  &                      )
          end if
       else
          status=0
          {Type¦match¦^(?!(Long|VarStr|Character))¦call extractDataContent(valueElement,parameterValue,iostat=status,ex=exception)¦}
          {Type¦match¦^(Long.*|(Character|VarStr)Rank1)$¦parameterText=getTextContent(valueElement,ex=exception)¦}
          {Type¦match¦^(VarStr|Character)$¦parameterValue=getTextContent(valueElement,ex=exception)¦}
          if (inException(exception) .or. status /= 0) then
             if (present(errorStatus)) then
                errorStatus=inputParameterErrorStatusParse
             else
                call Galacticus_Error_Report(                                          &
                     &                       'inputParametersValueNode{Type¦label}' ,  &
                     &                       'unable to parse parameter ['          // &
                     &                        getNodeName   (parameterNode%content)         // &
                     &                       ']='                                   // &
                     &                        getTextContent(valueElement )            &
                     &                      )
             end if
          else
             ! Convert type as necessary.
             {Type¦match¦^Long.*¦read (parameterText,*) parameterValue¦}
             {Type¦match¦^(Character|VarStr)Rank1$¦call String_Split_Words(parameterValue,char(parameterText))¦}
             ! Write the parameter file to an HDF5 object.
             if (self%outputParameters%isOpen().and.writeOutput_) then
                !$omp critical(HDF5_Access)
                if (.not.self%outputParameters%hasAttribute(getNodeName(parameterNode%content))) call self%outputParameters%writeAttribute({Type¦outputConverter¦parameterValue},getNodeName(parameterNode%content))
                !$omp end critical(HDF5_Access)
             end if
          end if
       end if
    end if
    !$omp end critical (FoX_DOM_Access)
    return
  end subroutine inputParametersValueNode{Type¦label}

  subroutine cInputParametersFinalize(cSelf) bind(c,name="inputParametersFinalize")
    !% A C-callable wrapper function to finalize an {\normalfont \ttfamily inputParameters} object.
    use, intrinsic :: ISO_C_Binding
    implicit none
    type(c_ptr          ), value   :: cSelf
    type(inputParameters), pointer :: self
    
    call c_f_pointer(cSelf,self)
    call inputParametersFinalize(self)
    return
  end subroutine cInputParametersFinalize

  subroutine cInputParametersSubParameters(cSelf,parameterNameLength,parameterName,subParameters) bind(c,name="inputParametersSubParameters")
    !% A C-callable wrapper function to obtain subparameters from an {\normalfont \ttfamily inputParameters} object.
    use, intrinsic :: ISO_C_Binding
    use               ISO_Varying_String
    implicit none
    type     (c_ptr          ), value                          :: cSelf
    integer  (kind=c_int     ), value                          :: parameterNameLength
    type     (inputParameters), pointer                        :: self
    type     (inputParameters), pointer                        :: subParametersF
    character(kind=c_char    ), dimension(parameterNameLength) :: parameterName
    type     (varying_string )                                 :: parameterNameF
    type     (c_ptr          )                                 :: subParameters
    
    parameterNameF=String_C_to_Fortran(parameterName)
    call c_f_pointer(cSelf,self)
    allocate(subParametersF)
    subParametersF=self%subParameters(char(parameterNameF))
    subParameters=c_loc(subParametersF)
    return
  end subroutine cInputParametersSubParameters

  subroutine cInputParametersValueName{cType¦label}(cSelf,parameterNameLength,parameterName,parameterValue,defaultValue,errorStatus,writeOutput) bind(c,name="inputParametersValueName{cType¦label}")
    !% C-callable interface to input parameter get-by-name functions.
    use, intrinsic :: ISO_C_Binding
    use               ISO_Varying_String
    implicit none
    type             (c_ptr          ), value                                    :: cSelf
    integer          (kind=c_int     ), value                                    :: parameterNameLength
    character        (kind=c_char    ), dimension(parameterNameLength)           :: parameterName
    {cType¦intrinsic}                 , intent(  out)                            :: parameterValue
    {cType¦intrinsic}                 , intent(in   )                            :: defaultValue
    integer          (kind=c_int     ), intent(  out)                            :: errorStatus
    logical          (kind=c_bool    ), intent(in   )                 , optional :: writeOutput
    type             (inputParameters), pointer                                  :: self
    type             (varying_string )                                           :: parameterNameF
    logical                                                                      :: writeOutput_

    parameterNameF=String_C_to_Fortran(parameterName)
    call c_f_pointer(cSelf,self)
    writeOutput_=.true.
    if (present(writeOutput)) writeOutput_=logical(writeOutput)
    call self%value(char(parameterNameF),parameterValue,defaultValue,errorStatus,writeOutput_)
    return
  end subroutine cInputParametersValueName{cType¦label}

  function inputParameterListConstructor()
    !% Construct an {\normalfont \ttfamily inputParameterList} object.
    implicit none
    type(inputParameterList) :: inputParameterListConstructor

    return
  end function inputParameterListConstructor

  subroutine inputParameterListDestructor(self)
    !% Destroy an {\normalfont \ttfamily inputParameterList} object.
    implicit none
    type(inputParameterList), intent(inout) :: self

    if (allocated(self%name )) deallocate(self%name )
    if (allocated(self%value)) deallocate(self%value)
    return
  end subroutine inputParameterListDestructor

  subroutine inputParameterListAdd(self,name,value)
    !% Add a parameter to a list of input parameters to an XML document.
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
    !% Serialize a list of input parameters to an XML document.
    use FoX_wXML
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
    !% Serialize input parameters to a string.
    use Hashes_Cryptographic
    implicit none
    type   (varying_string )                          :: inputParametersSerializeToString
    class  (inputParameters), intent(inout)           :: self
    logical                 , intent(in   ), optional :: hashed 
    type   (node           ), pointer                 :: thisNode
    type   (inputParameter ), pointer                 :: currentParameter
    integer                                           :: i                               , errorStatus
    type   (varying_string )                          :: parameterValue
    logical                                           :: firstParameter
    type   (inputParameters)                          :: subParameters
    !# <optionalArgument name="hashed" defaultsTo=".false." />

    inputParametersSerializeToString=""
    firstParameter=.true.
    currentParameter => self%parameters%firstChild
    do while (associated(currentParameter))
       if (currentParameter%isParameter()) then
          if (.not.firstParameter) inputParametersSerializeToString=inputParametersSerializeToString//"_"
          call self%value(currentParameter,parameterValue,errorStatus,writeOutput=.false.)
          thisNode => currentParameter%content
          inputParametersSerializeToString=inputParametersSerializeToString// &
               &                           getNodeName(thisNode)           // &
               &                           ":"                             // &
               &                           adjustl(trim(parameterValue))
          subParameters=self%subParameters(getNodeName(thisNode))
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
    !% Serialize input parameters to an XML file.
    implicit none
    class(inputParameters), intent(in   ) :: self
    type (varying_string ), intent(in   ) :: parameterFile
    
    call serialize(self%document,char(parameterFile))
    return
  end subroutine inputParametersSerializeToXML
    
  subroutine inputParametersAddParameter(self,parameterName,parameterValue)
    !% Add a parameter to the set.
    implicit none
    class    (inputParameters), intent(inout) :: self
    character(len=*          ), intent(in   ) :: parameterName   , parameterValue
    type     (node           ), pointer       :: parameterNode   , dummy
    type     (inputParameter ), pointer       :: currentParameter

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
    return
  end subroutine inputParametersAddParameter
  
  ! Include functions that generate unique labels for modules.
  include 'utility.input_parameters.unique_labels.inc'

end module Input_Parameters2
