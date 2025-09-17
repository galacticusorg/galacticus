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
Contains a module which implements various utility functions for extracting data from XML files.
!!}

module IO_XML
  !!{
  Implements various utility functions for extracting data from XML files.
  !!}
  use :: FoX_dom           , only : node
  use :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: XML_Extrapolation_Element_Decode , XML_Array_Read                , XML_Array_Read_Static       , &
       &    XML_Get_First_Element_By_Tag_Name, XML_Count_Elements_By_Tag_Name, XML_Path_Exists             , &
       &    XML_Extract_Text                 , XML_Parse                     , XML_Get_Elements_By_Tag_Name, &
       &    xmlNodeList                      , XML_Get_Child_Elements

  ! Interface for array reading functions.
  interface XML_Array_Read
     module procedure XML_Array_Read_One_Column
     module procedure XML_Array_Read_Two_Column
     module procedure XML_List_Array_Read_One_Column
  end interface XML_Array_Read
  interface XML_Array_Read_Static
     module procedure XML_Array_Read_Static_One_Column
     module procedure XML_List_Double_Array_Read_Static_One_Column
     module procedure XML_List_Integer_Array_Read_Static_One_Column
     module procedure XML_List_Character_Array_Read_Static_One_Column
  end interface XML_Array_Read_Static

  ! Interface for parsing function.
  interface XML_Parse
     module procedure XML_Parse_VarStr
     module procedure XML_Parse_Char
  end interface XML_Parse
  
  type :: xincludeNode
     !!{
     Type used while resolving XInclude references during XML parsing.
     !!}
     type(node          ), pointer :: nodeParent => null(), nodeXInclude => null()
     type(varying_string)          :: fileName            , xPath
  end type xincludeNode

  type :: xincludeNodeList
     !!{
     Type used while resolving XInclude references during XML parsing.
     !!}
     type(xmlNodeList), allocatable, dimension(:) :: nodes
  end type xincludeNodeList

  type :: xmlNodeList
     !!{
     Type used to provide lists of XML nodes.
     !!}
     type(node), pointer :: element => null()
  end type xmlNodeList
  
contains

  function XML_Extract_Text(xmlElement)
    !!{
    Extract the text from an XML element and return as a variable length string.
    !!}
    use :: FoX_dom           , only : getTextContent, node
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type(varying_string)                         :: XML_Extract_Text
    type(node          ), intent(in   ), pointer :: xmlElement

    XML_Extract_Text=getTextContent(xmlElement)
    return
  end function XML_Extract_Text

  subroutine XML_Array_Read_Static_One_Column(xmlElement,arrayElementName,column1)
    !!{
    Read one column of data from an array of XML elements.
    !!}
    use :: FoX_dom, only : extractDataContent, getElementsByTagName, node
    implicit none
    type            (node       )              , intent(in   ), pointer :: xmlElement
    character       (len=*      )              , intent(in   )          :: arrayElementName
    double precision             , dimension(:), intent(inout)          :: column1
    type            (node       )                             , pointer :: arrayElement
    type            (xmlNodeList), dimension(:), allocatable            :: arrayElements
    double precision             , dimension(1)                         :: dataValues
    integer                                                             :: i

    call XML_Get_Elements_By_Tag_Name(xmlElement,arrayElementName,arrayElements)
    do i=1,size(arrayElements)
       arrayElement => arrayElements(i-1)%element
       call extractDataContent(arrayElement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_Array_Read_Static_One_Column

  subroutine XML_Array_Read_One_Column(xmlElement,arrayElementName,column1)
    !!{
    Read one column of data from an array of XML elements.
    !!}
    use :: FoX_dom          , only : extractDataContent, getElementsByTagName, node
    implicit none
    type            (node       )                           , intent(in   ), pointer :: xmlElement
    character       (len=*      )                           , intent(in   )          :: arrayElementName
    double precision             , allocatable, dimension(:), intent(inout)          :: column1
    type            (node       )                                          , pointer :: arrayElement
    type            (xmlNodeList), allocatable, dimension(:)                         :: arrayElements
    double precision                          , dimension(1)                         :: dataValues
    integer                                                                          :: i

    call XML_Get_Elements_By_Tag_Name(xmlElement,arrayElementName,arrayElements)
    if (allocated(column1)) deallocate(column1)
    allocate(column1(size(arrayElements)))
    do i=1,size(arrayElements)
       arrayElement => arrayElements(i-1)%element
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_Array_Read_One_Column

  subroutine XML_Array_Read_Two_Column(xmlElement,arrayElementName,column1,column2)
    !!{
    Read two columns of data from an array of XML elements.
    !!}
    use :: FoX_dom          , only : extractDataContent, getElementsByTagName, node
    implicit none
    type            (node       )                           , intent(in   ), pointer :: xmlElement
    character       (len=*      )                           , intent(in   )          :: arrayElementName
    double precision             , allocatable, dimension(:), intent(inout)          :: column1         , column2
    type            (node       )                                          , pointer :: arrayElement
    type            (xmlNodeList), allocatable, dimension(:)                         :: arrayElements
    double precision                          , dimension(2)                         :: dataValues
    integer                                                                          :: i

    call XML_Get_Elements_By_Tag_Name(xmlElement,arrayElementName,arrayElements)
    if (allocated(column1)) deallocate(column1)
    if (allocated(column2)) deallocate(column2)
    allocate(column1(size(arrayElements)))
    allocate(column2(size(arrayElements)))
    do i=1,size(arrayElements)
       arrayElement => arrayElements(i-1)%element
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
       column2(i)=dataValues(2)
    end do
    return
  end subroutine XML_Array_Read_Two_Column

  subroutine XML_List_Array_Read_One_Column(xmlElements,arrayElementName,column1)
    !!{
    Read one column of data from an array of XML elements.
    !!}
    use :: FoX_dom          , only : extractDataContent, node
    implicit none
    type            (xmlNodeList)             , dimension(0:), intent(in   ) :: xmlElements
    character       (len=*      )                            , intent(in   ) :: arrayElementName
    double precision             , allocatable, dimension(: ), intent(inout) :: column1
    type            (node       ), pointer                                   :: arrayElement
    double precision                          , dimension(1 )                :: dataValues
    integer                                                                  :: i

    if (allocated(column1)) deallocate(column1)
    allocate(column1(size(xmlElements)))
    do i=1,size(xmlElements)
       arrayElement => XML_Get_First_Element_By_Tag_Name(xmlElements(i-1)%element,arrayElementName)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_List_Array_Read_One_Column

  subroutine XML_List_Double_Array_Read_Static_One_Column(xmlElements,arrayElementName,column1)
    !!{
    Read one column of integer data from an array of XML elements.
    !!}
    use :: FoX_dom, only : extractDataContent, node
    implicit none
    type            (xmlNodeList), dimension(0:), intent(in   )          :: xmlElements
    character       (len=*      )               , intent(in   )          :: arrayElementName
    double precision             , dimension(: ), intent(inout)          :: column1
    type            (node       )                              , pointer :: arrayElement
    double precision             , dimension(1 )                         :: dataValues
    integer                                                              :: i

    do i=1,size(xmlElements)
       arrayElement => XML_Get_First_Element_By_Tag_Name(xmlElements(i-1)%element,arrayElementName)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_List_Double_Array_Read_Static_One_Column

  subroutine XML_List_Integer_Array_Read_Static_One_Column(xmlElements,arrayElementName,column1)
    !!{
    Read one column of integer data from an array of XML elements.
    !!}
    use :: FoX_dom, only : extractDataContent, node
    implicit none
    type     (xmlNodeList), dimension(0:), intent(in   )          :: xmlElements
    character(len=*      )               , intent(in   )          :: arrayElementName
    integer               , dimension(: ), intent(inout)          :: column1
    type     (node       )                              , pointer :: arrayElement
    integer               , dimension(1 )                         :: dataValues
    integer                                                       :: i

    do i=1,size(xmlElements)
       arrayElement => XML_Get_First_Element_By_Tag_Name(xmlElements(i-1)%element,arrayElementName)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_List_Integer_Array_Read_Static_One_Column

  subroutine XML_List_Character_Array_Read_Static_One_Column(xmlElements,arrayElementName,column1)
    !!{
    Read one column of character data from an array of XML elements.
    !!}
    use :: FoX_dom, only : extractDataContent, node
    implicit none
    type     (xmlNodeList     ), dimension(0:), intent(in   )          :: xmlElements
    character(len=*           )               , intent(in   )          :: arrayElementName
    character(len=*           ), dimension(: ), intent(inout)          :: column1
    type     (node            )                              , pointer :: arrayElement
    character(len=len(column1)), dimension(1 )                         :: dataValues
    integer                                                            :: i

    do i=1,size(xmlElements)
       arrayElement => XML_Get_First_Element_By_Tag_Name(xmlElements(i-1)%element,arrayElementName)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_List_Character_Array_Read_Static_One_Column
  
  subroutine XML_Get_Child_Elements(xmlElement,elements)
    !!{
    Return a list of pointers to all child nodes.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: FoX_DOM      , only : getFirstChild, getNextSibling, hasChildNodes, node
    implicit none
    type   (xmlNodeList), intent(inout), allocatable, dimension(:) :: elements
    type   (node       ), intent(in   ), pointer                   :: xmlElement
    type   (node       )               , pointer                   :: childNode
    integer(c_size_t   )                                           :: countElements
    
    countElements=0_c_size_t
    if (hasChildNodes(xmlElement)) then
       childNode => getFirstChild(xmlElement)
       do while (associated(childNode))
          countElements=countElements+1_c_size_t
          childNode => getNextSibling(childNode)
       end do
    end if
    if (allocated(elements)) deallocate(elements)
    allocate(elements(0:countElements-1))
    if (hasChildNodes(xmlElement)) then
       countElements =  0_c_size_t
       childNode     => getFirstChild(xmlElement)
       do while (associated(childNode))
          elements(countElements)%element => childNode
          countElements=countElements+1_c_size_t
          childNode => getNextSibling(childNode)
       end do
    end if
    return
  end subroutine XML_Get_Child_Elements
  
  recursive subroutine XML_Get_Elements_By_Tag_Name(xmlElement,tagName,elements)
    !!{
    Return a list of pointers to all nodes matching a given {\normalfont \ttfamily tagName}.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: FoX_DOM      , only : Element_Node, getFirstChild, getNextSibling, getNodeName , &
          &                                 getNodeType , hasChildNodes, node          , getAttribute
    implicit none
    type     (xmlNodeList     ), intent(inout), allocatable, dimension(:) :: elements
    integer  (c_size_t        )                                           :: countElements , offset
    type     (node            ), intent(in   ), pointer                   :: xmlElement
    character(len=*           ), intent(in   )                            :: tagName
    type     (node            )               , pointer                   :: childNode
    type     (xmlNodeList     )               , allocatable, dimension(:) :: childElements
    logical                                                               :: matchAll      , matches
    character(len=len(tagName))                                           :: tagName_      , attributeName, &
         &                                                                   attributeValue

    countElements=XML_Count_Elements_By_Tag_Name(xmlElement,tagName)    
    offset       =0_c_size_t
    if (allocated(elements)) deallocate(elements)
    allocate(elements(0:countElements-1))
    if (hasChildNodes(xmlElement)) then
      if (index(tagName,"[@") > 0) then
         tagName_      =tagName(                    1:index(tagName,"[@")-1)
         attributeName =tagName(index(tagName,"[@")+2:index(tagName,"=" )-1)
         attributeValue=tagName(index(tagName,"=" )+2:index(tagName,"]" )-2)
      else
          tagName_      =tagName
          attributeName =""
          attributeValue=""
       end if
       matchAll  =  trim(tagName_) == "*"
       childNode => getFirstChild(xmlElement)
       do while (associated(childNode))
          call XML_Get_Elements_By_Tag_Name(childNode,tagName,childElements)
          elements(offset:offset+size(childElements)-1)=childElements
          offset=offset+size(childElements)
          deallocate(childElements)
          if (getNodeType(childNode) == Element_Node .and. (matchAll .or. getNodeName(childNode) == trim(tagName_))) then
             matches=attributeName == "" .or. getAttribute(childNode,trim(attributeName)) == trim(attributeValue)
             if (matches) then
                elements(offset)%element => childNode
                offset=offset+1_c_size_t
             end if
          end if
          childNode => getNextSibling(childNode)
       end do
    end if
    return
  end subroutine XML_Get_Elements_By_Tag_Name
  
  recursive function XML_Count_Elements_By_Tag_Name(xmlElement,tagName) result(countElements)
    !!{
    Return a count of all nodes matching a given {\normalfont \ttfamily tagName}.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: FoX_DOM      , only : Element_Node, getFirstChild, getNextSibling, getNodeName, &
          &                                 getNodeType , hasChildNodes, node          , getAttribute
    implicit none
    integer  (c_size_t        )                         :: countElements
    type     (node            ), intent(in   ), pointer :: xmlElement
    character(len=*           ), intent(in   )          :: tagName
    type     (node            )               , pointer :: childNode
    logical                                             :: matchAll
    character(len=len(tagName))                         :: tagName_      , attributeName, &
         &                                                 attributeValue

    countElements=0_c_size_t
    if (hasChildNodes(xmlElement)) then
       if (index(tagName,"[@") > 0) then
          tagName_      =tagName(                    1:index(tagName,"[@")-1)
          attributeName =tagName(index(tagName,"[@")+2:index(tagName,"=" )-1)
          attributeValue=tagName(index(tagName,"=" )+2:index(tagName,"]" )-2)
       else
          tagName_      =tagName
          attributeName =""
          attributeValue=""
       end if
       matchAll  =  trim(tagName_) == "*"
       childNode => getFirstChild(xmlElement)
       do while (associated(childNode))
          countElements=countElements+XML_Count_Elements_By_Tag_Name(childNode,tagName)
          if (getNodeType(childNode) == Element_Node .and. (matchAll .or. getNodeName(childNode) == trim(tagName_))) then
             if (attributeName == "" .or. getAttribute(childNode,trim(attributeName)) == trim(attributeValue)) countElements=countElements+1_c_size_t
          end if
          childNode => getNextSibling(childNode)
       end do
    end if
    return
  end function XML_Count_Elements_By_Tag_Name

  function XML_Get_First_Element_By_Tag_Name(xmlElement,tagName,directChildrenOnly) result(element)
    !!{
    Return a pointer to the first node in an XML node that matches the given {\normalfont \ttfamily tagName}.
    !!}
    use :: FoX_dom, only : getParentNode, node
    use :: Error  , only : Error_Report
    implicit none
    type     (node            )               , pointer      :: element
    type     (node            ), intent(in   ), pointer      :: xmlElement
    character(len=*           ), intent(in   )               :: tagName
    logical                    , intent(in   ), optional     :: directChildrenOnly
    type     (xmlNodeList     ), allocatable  , dimension(:) :: elementList
    type     (node            )               , pointer      :: parent
    character(len=len(tagName))                              :: currentTagName                   , path
    integer                                                  :: pathPosition                     , i
    logical                                                  :: found
    !![
    <optionalArgument name="directChildrenOnly" defaultsTo=".false."/>
    !!]

    ! Validate.
    if (index(tagName,"//") > 0                               ) call Error_Report('XPath `//` operator is not supported'                       //{introspection:location})
    if (index(tagName,"/" ) > 0 .and. .not.directChildrenOnly_) call Error_Report('only direct children supported when using XPath expressions'//{introspection:location})
    ! Find element.
    element => xmlElement
    path=tagName
    do while (path /= "")
       pathPosition=index(path,"/")
       if (pathPosition == 0) then
          currentTagName=path
          path          =""
       else
          currentTagName=path(             1:pathPosition-1)
          path          =path(pathPosition+1:len(path)     )
       endif
       call XML_Get_Elements_By_Tag_Name(element,currentTagName,elementList)
       if (size(elementList) < 1) then
          call Error_Report('no elements match tag name "'//trim(currentTagName)//'"'//{introspection:location})
       else
          if (directChildrenOnly_) then
             found=.false.
             do i=0,size(elementList)-1
                parent => getParentNode(elementList(i)%element)
                if (associated(parent,element)) then
                   found=.true.
                   element => elementList(i)%element
                   exit
                end if
             end do
             if (.not.found) call Error_Report('no direct child elements match tag name "'//trim(currentTagName)//'"'//{introspection:location})
          else
             element => elementList(0)%element
          end if
       end if
    end do
    return
  end function XML_Get_First_Element_By_Tag_Name

  logical function XML_Path_Exists(xmlElement,path)
    !!{
    Return true if the supplied {\normalfont \ttfamily path} exists in the supplied {\normalfont \ttfamily xmlElement}.
    !!}
    use :: FoX_dom, only : ELEMENT_NODE , getElementsByTagName, getLength, getNodeType, &
          &                getParentNode, node
    implicit none
    type     (node         ), intent(in   ), pointer      :: xmlElement
    character(len=*        ), intent(in   )               :: path
    type     (node         )               , pointer      :: element     , child         , &
         &                                                   parent
    character(len=len(path))                              :: currentPath , currentTagName
    integer                                               :: pathPosition, i
    type     (xmlNodeList  ), allocatable  , dimension(:) :: elementList

    XML_Path_Exists =  .true.
    element         => xmlElement
    currentPath     =  path
    do while (currentPath /= "")
       pathPosition=index(currentPath,"/")
       if (pathPosition == 0) then
          currentTagName=currentPath
          currentPath   =""
       else
          currentTagName=currentPath(             1:pathPosition-1)
          currentPath   =currentPath(pathPosition+1:len(path)     )
       endif
       call XML_Get_Elements_By_Tag_Name(element,currentTagName,elementList)
       if (size(elementList) < 1) then
          XML_Path_Exists=.false.
          return
       else
          XML_Path_Exists=.false.
          do i=0,size(elementList)-1
             child  => elementList  (i    )%element
             parent => getParentNode(child)
             if (getNodeType(child) == ELEMENT_NODE .and. associated(parent,element)) then
                element => child
                XML_Path_Exists=.true.
                exit
             end if
          end do
          if (.not.XML_Path_Exists) return
       end if
    end do
    return
  end function XML_Path_Exists

  subroutine XML_Extrapolation_Element_Decode(extrapolationElement,limitType,extrapolationMethod,allowedMethods)
    !!{
    Extracts information from a standard XML {\normalfont \ttfamily extrapolationElement}. Optionally a set of {\normalfont \ttfamily allowedMethods} can be
    specified---if the extracted method does not match one of these an error is issued.
    !!}
    use :: FoX_dom     , only : extractDataContent                , node
    use :: Error       , only : Error_Report
    use :: Table_Labels, only : enumerationExtrapolationTypeEncode, enumerationExtrapolationTypeType
    implicit none
    type     (node                            )              , intent(in   ), pointer     :: extrapolationElement
    character(len=*                           )              , intent(  out)              :: limitType
    type     (enumerationExtrapolationTypeType)              , intent(  out)              :: extrapolationMethod
    type     (enumerationExtrapolationTypeType), dimension(:), intent(in   ), optional    :: allowedMethods
    type     (node                            )                             , pointer     :: limitElement        , methodElement
    type     (xmlNodeList                     ), dimension(:)               , allocatable :: elementList
    character(len=32                          )                                           :: methodType

    ! Extract the limit type.
    call XML_Get_Elements_By_Tag_Name(extrapolationElement,"limit",elementList)
    if (size(elementList) /= 1) call Error_Report('extrapolation element must contain exactly one limit element'//{introspection:location})
    limitElement => elementList(0)%element
    call extractDataContent(limitElement,limitType)
    ! Extract the method type.
    call XML_Get_Elements_By_Tag_Name(extrapolationElement,"method",elementList)
    if (size(elementList) /= 1) call Error_Report('extrapolation element must contain exactly one method element'//{introspection:location})
    methodElement => elementList(0)%element
    call extractDataContent(methodElement,methodType)
    extrapolationMethod=enumerationExtrapolationTypeEncode(trim(methodType),includesPrefix=.false.)
    ! Validate the method type.
    if (present(allowedMethods)) then
       if (all(allowedMethods /= extrapolationMethod)) call Error_Report('unallowed extrapolation method'//{introspection:location})
    end if
    return
  end subroutine XML_Extrapolation_Element_Decode
  
  function XML_Parse_VarStr(fileName,iostat,ex,fileNameCurrent) result(document)
    !!{
    Parse an XML document, automatically resolve XInclude references.
    !!}
    use :: FoX_dom           , only : DOMException, node
    use :: ISO_Varying_String, only : char
    implicit none
    type   (node          ), pointer                 :: document
    type   (varying_string), intent(in   )           :: fileName
    type   (varying_string), intent(  out), optional :: fileNameCurrent
    integer                , intent(inout), optional :: iostat
    type   (DOMException  ), intent(  out), optional :: ex

    document => XML_Parse(char(fileName),iostat,ex,fileNameCurrent)
    return
  end function XML_Parse_VarStr
  
  function XML_Parse_Char(fileName,iostat,ex,fileNameCurrent) result(document)
    !!{
    Parse an XML document, automatically resolve XInclude references.
    !!}
    use :: File_Utilities    , only : File_Exists  , File_Name         , File_Path    , File_Name_Expand
    use :: FoX_dom           , only : DOMException , ELEMENT_NODE      , destroy      , getAttribute    , &
          &                           getChildNodes, getDocumentElement, getFirstChild, getNextSibling  , &
          &                           getNodeName  , getNodeType       , getParentNode, hasAttribute    , &
          &                           hasChildNodes, importNode        , insertBefore , node            , &
          &                           parseFile    , removeChild       , replaceChild , setLiveNodeLists, &
          &                           setAttribute , inException
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=), char              , extract      , len             , &
          &                           operator(//) , operator(==)      , operator(/=)
    implicit none
    type     (node            ), pointer                     :: document           , nodeNew       , &
         &                                                      nodeCurrent        , nodeParent    , &
         &                                                      nodeXInclude       , nodeImported  , &
         &                                                      nodeInsert         , nodeNext
    character(len=*           ), intent(in   )               :: fileName
    type     (varying_string  ), intent(  out), optional     :: fileNameCurrent
    integer                    , intent(inout), optional     :: iostat
    type     (DOMException    ), intent(  out), optional     :: ex
    type     (xmlNodeList     ), allocatable  , dimension(:) :: nodesCurrent
    type     (xincludeNode    ), allocatable  , dimension(:) :: stack              , stackTmp
    type     (xincludeNodeList), allocatable  , dimension(:) :: stackList          , stackListTmp
    integer                    , parameter                   :: stackExpandCount=10
    integer                                                  :: stackCount         , stackListCount, &
         &                                                      i                  , countElements
    type     (varying_string  )                              :: filePath           , fileLeaf      , &
         &                                                      nameInsert         , fileNameFull  , &
         &                                                      fileName_
    logical                                                  :: allElements

    ! Expand the file name.
    fileName_=File_Name_Expand(fileName)
    ! Extract the path and leaf name to our document.
    filePath=File_Path(fileName_)
    fileLeaf=File_Name(fileName_)
    ! Initialize the XInclude reference stack.
    allocate(stack(stackExpandCount))
    stackCount                     =  1
    document                       => null()
    stack(stackCount)%nodeParent   => null()
    stack(stackCount)%nodeXInclude => null()
    stack(stackCount)%fileName     =  fileLeaf
    stack(stackCount)%xPath        =  ""
    ! Initialize the nodeList stack.
    allocate(stackList(stackExpandCount))
    ! Process the document.
    allElements=.false.
    do while (stackCount > 0)
       ! Construct the full filename.
       if (extract(stack(stackCount)%fileName,1,1) == "/") then
          fileNameFull=          stack(stackCount)%fileName
       else
          fileNameFull=filePath//stack(stackCount)%fileName
       end if
       if (present(fileNameCurrent)) fileNameCurrent=fileNameFull
       ! Check that file exists.
       if (.not.File_Exists(fileNameFull)) call Error_Report('file "'//char(fileNameFull)//'" does not exist'//{introspection:location})
       ! Parse the document.
       nodeNew => parseFile(char(fileNameFull),iostat=iostat,ex=ex)
       if (present(iostat).and.iostat /= 0.or.present(ex).and.inException(ex)) return
       ! Paths in any xi:include elements are relative to the file they are defined in. We must update this to be relative to our base parameter file.
       call XML_Get_Elements_By_Tag_Name(nodeNew,"xi:include",nodesCurrent)
       do i=0,size(nodesCurrent)-1
          nodeCurrent => nodesCurrent(i)%element
          if (getNodeName(nodeCurrent) == "xi:include") then
             if (.not.hasAttribute(nodeCurrent,"href")) call Error_Report("missing 'href' in XInclude"//{introspection:location})
             fileNameFull=getAttribute(nodeCurrent,"href")
             if (extract(fileNameFull,1,1) /= "/") then
                fileNameFull=File_Path(stack(stackCount)%fileName)//getAttribute(nodeCurrent,"href")
                call setAttribute(nodeCurrent,"href",char(fileNameFull))
             end if
          end if
       end do
       ! Process the file into our combined document.
       nodeParent   => stack(stackCount)%nodeParent
       nodeXInclude => stack(stackCount)%nodeXInclude
       if (stack(stackCount)%xPath == "") then
          nodeInsert => getDocumentElement(nodeNew)
          nameInsert =  ""
       else
          if (XML_Path_Exists(nodeNew,char(stack(stackCount)%xPath))) then
             nodeInsert  => XML_Get_First_Element_By_Tag_Name(nodeNew,char(stack(stackCount)%xPath),directChildrenOnly=.true.)
             nameInsert  =  getNodeName  (nodeInsert)
             nodeInsert  => getParentNode(nodeInsert)
             allElements =  extract(stack(stackCount)%xPath,len(stack(stackCount)%xPath),len(stack(stackCount)%xPath)) == "*"
          else
             call Error_Report("XPath '"//stack(stackCount)%xPath//"' not found"//{introspection:location})
          end if
       end if
       ! We drop the entire stack after popping off just one element. This is because when we insert a new node into the document
       ! it will change all of the pointers, so we must rescan it to handle any additional xi:include elements.
       stackCount=0
       ! Insert this document.
       if (associated(nodeParent)) then
          ! Insert the newly parsed document into the parent node.
          if (nameInsert == "") then
             nodeImported => importNode  (document  ,nodeInsert  ,.true.      )
             call destroy(nodeNew)
             nodeNew      => replaceChild(nodeParent,nodeImported,nodeXInclude)
          else
             countElements =  0
             nodeCurrent   => getFirstChild(nodeInsert)
             do while (associated(nodeCurrent))
                if ((allElements .and. getNodeType(nodeCurrent) == ELEMENT_NODE) .or. getNodeName(nodeCurrent) == nameInsert) countElements=countElements+1
                nodeCurrent => getNextSibling(nodeCurrent)
             end do
             i           =  0
             nodeCurrent => getFirstChild(nodeInsert)
             do while (associated(nodeCurrent))
                nodeNext => getNextSibling(nodeCurrent)
                if ((getNodeType(nodeCurrent) /= ELEMENT_NODE .and. i > 0 .and. i < countElements) .or. (allElements .and. getNodeType(nodeCurrent) == ELEMENT_NODE) .or. getNodeName(nodeCurrent) == nameInsert) then
                   if (getNodeName(nodeCurrent) == nameInsert) i=i+1
                   nodeImported => importNode  (document,nodeCurrent,.true.)
                   nodeCurrent  => insertBefore(getParentNode(nodeXInclude),nodeImported,nodeXInclude)
                end if
                nodeCurrent => nodeNext
             end do
             call destroy(nodeNew)
             nodeNew      => getParentNode(        nodeXInclude) ! Reprocess from the parent to ensure we capture all newly added nodes.
             nodeXInclude => removeChild  (nodeNew,nodeXInclude)
          end if
       else
          ! No node, this is therefore the base document.
          document => nodeNew
          call setLiveNodeLists(document,.false.)
       end if
       ! Search for any XIncludes - rescan the entire document as we can handle inserting only one xi:include element at a time.
       stackListCount=0
       if (hasChildNodes(document)) then
          stackListCount=1
          call XML_Get_Child_Elements(document,stackList(stackListCount)%nodes)
       end if
       do while (stackListCount > 0)
          nodesCurrent  =stackList(stackListCount)%nodes
          stackListCount=          stackListCount       -1
          do i=0,size(nodesCurrent)-1
             nodeCurrent => nodesCurrent(i)%element
             if (getNodeName(nodeCurrent) == "xi:include") then
                if (.not.hasAttribute(nodeCurrent,"href")) call Error_Report("missing 'href' in XInclude"//{introspection:location})
                if (stackCount == size(stack)) then
                   call Move_Alloc(stack,stackTmp)
                   allocate(stack(stackCount+stackExpandCount))
                   stack(1:stackCount)=stackTmp
                   deallocate(stackTmp)
                end if
                stackCount                     =  stackCount+1
                stack(stackCount)%fileName     =  getAttribute (nodeCurrent,"href")
                stack(stackCount)%nodeParent   => getParentNode(nodeCurrent       )
                stack(stackCount)%nodeXInclude =>               nodeCurrent
                if (hasAttribute(nodeCurrent,"xpointer")) then
                   stack(stackCount)%xPath=getAttribute(nodeCurrent,"xpointer")
                   if     (                                                                                                           &
                        &   extract(stack(stackCount)%xPath,                          1 ,                          9 ) == "xpointer(" &
                        &  .and.                                                                                                      &
                        &   extract(stack(stackCount)%xPath,len(stack(stackCount)%xPath),len(stack(stackCount)%xPath)) == ")"         &
                        & ) then
                      stack(stackCount)%xPath=extract(stack(stackCount)%xPath,10,len(stack(stackCount)%xPath)-1)
                   else
                      call Error_Report("malformed XPath in XPointer: '"//stack(stackCount)%xPath//"'"//{introspection:location})
                   end if
                else
                   stack(stackCount)%xPath=""
                end if
             end if
             ! Add any child nodes to the nodeList stack.
             if (hasChildNodes(nodeCurrent)) then
                if (stackListCount == size(stackList)) then
                   call Move_Alloc(stackList,stackListTmp)
                   allocate(stackList(stackListCount+stackExpandCount))
                   stackList(1:stackListCount)=stackListTmp
                   deallocate(stackListTmp)
                end if
                stackListCount=stackListCount+1
                call XML_Get_Child_Elements(nodeCurrent,stackList(stackListCount)%nodes)
             end if
          end do
       end do
    end do
    call setLiveNodeLists(document,.true.)  
    return
  end function XML_Parse_Char

end module IO_XML
