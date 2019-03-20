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

!% Contains a module which implements various utility functions for extracting data from XML files.

module IO_XML
  !% Implements various utility functions for extracting data from XML files.
  use FoX_dom
  use ISO_Varying_String
  implicit none
  private
  public :: XML_Extrapolation_Element_Decode , XML_Array_Read  , XML_Array_Read_Static, &
       &    XML_Get_First_Element_By_Tag_Name, XML_Array_Length, XML_Path_Exists      , &
       &    XML_Extract_Text                 , XML_Parse

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

  type :: xincludeNode
     !% Type used while resolving XInclude references during XML parsing.
     type(node          ), pointer :: nodeParent, nodeXInclude
     type(varying_string)          :: fileName  , xPath
  end type xincludeNode
  
  type :: xincludeNodeList
     !% Type used while resolving XInclude references during XML parsing.
     type(nodeList), pointer :: nodes
  end type xincludeNodeList
  
contains

  function XML_Extract_Text(xmlElement)
    !% Extract the text from an XML element and return as a variable length string.
    implicit none
    type(varying_string)                         :: XML_Extract_Text
    type(node          ), intent(in   ), pointer :: xmlElement

    XML_Extract_Text=getTextContent(xmlElement)
    return
  end function XML_Extract_Text

  integer function XML_Array_Length(xmlElement,arrayElementName)
    !% Return the length of an array of XML elements.
    implicit none
    type     (node    ), intent(in   ), pointer :: xmlElement
    character(len=*   ), intent(in   )          :: arrayElementName
    type     (nodeList)               , pointer :: arrayElements

    arrayElements => getElementsByTagName(xmlElement,arrayElementName)
    XML_Array_Length=getLength(arrayElements)
    return
  end function XML_Array_Length

  subroutine XML_Array_Read_Static_One_Column(xmlElement,arrayElementName,column1)
    !% Read one column of data from an array of XML elements.
    implicit none
    type            (node    )              , intent(in   ), pointer :: xmlElement
    character       (len=*   )              , intent(in   )          :: arrayElementName
    double precision          , dimension(:), intent(inout)          :: column1
    type            (node    )                             , pointer :: arrayElement
    type            (nodeList)                             , pointer :: arrayElements
    double precision          , dimension(1)                         :: dataValues
    integer                                                          :: i

    arrayElements => getElementsByTagName(xmlElement,arrayElementName)
    do i=1,getLength(arrayElements)
       arrayElement => item(arrayElements,i-1)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_Array_Read_Static_One_Column

  subroutine XML_Array_Read_One_Column(xmlElement,arrayElementName,column1)
    !% Read one column of data from an array of XML elements.
    use Memory_Management
    implicit none
    type            (node    )                           , intent(in   ), pointer :: xmlElement
    character       (len=*   )                           , intent(in   )          :: arrayElementName
    double precision          , allocatable, dimension(:), intent(inout)          :: column1
    type            (node    )                                          , pointer :: arrayElement
    type            (nodeList)                                          , pointer :: arrayElements
    double precision                       , dimension(1)                         :: dataValues
    integer                                                                       :: i

    arrayElements => getElementsByTagName(xmlElement,arrayElementName)
    call allocateArray(column1,[getLength(arrayElements)])
    do i=1,getLength(arrayElements)
       arrayElement => item(arrayElements,i-1)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_Array_Read_One_Column

  subroutine XML_Array_Read_Two_Column(xmlElement,arrayElementName,column1,column2)
    !% Read two columns of data from an array of XML elements.
    use Memory_Management
    implicit none
    type            (node    )                           , intent(in   ), pointer :: xmlElement
    character       (len=*   )                           , intent(in   )          :: arrayElementName
    double precision          , allocatable, dimension(:), intent(inout)          :: column1         , column2
    type            (node    )                                          , pointer :: arrayElement
    type            (nodeList)                                          , pointer :: arrayElements
    double precision                       , dimension(2)                         :: dataValues
    integer                                                                       :: i

    arrayElements => getElementsByTagName(xmlElement,arrayElementName)
    call allocateArray(column1,[getLength(arrayElements)])
    call allocateArray(column2,[getLength(arrayElements)])
    do i=1,getLength(arrayElements)
       arrayElement => item(arrayElements,i-1)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
       column2(i)=dataValues(2)
    end do
    return
  end subroutine XML_Array_Read_Two_Column

  subroutine XML_List_Array_Read_One_Column(xmlElements,arrayElementName,column1)
    !% Read one column of data from an array of XML elements.
    use Memory_Management
    implicit none
    type            (nodeList)                           , intent(in   ), pointer :: xmlElements
    character       (len=*   )                           , intent(in   )          :: arrayElementName
    double precision          , allocatable, dimension(:), intent(inout)          :: column1
    type            (node    )                                          , pointer :: arrayElement    , xmlElement
    double precision                       , dimension(1)                         :: dataValues
    integer                                                                       :: i

    call allocateArray(column1,[getLength(xmlElements)])
    do i=1,getLength(xmlElements)
       xmlElement   => item                             (xmlElements,i-1             )
       arrayElement => XML_Get_First_Element_By_Tag_Name(xmlElement ,arrayElementName)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_List_Array_Read_One_Column

  subroutine XML_List_Double_Array_Read_Static_One_Column(xmlElements,arrayElementName,column1)
    !% Read one column of integer data from an array of XML elements.
    implicit none
    type            (nodeList)              , intent(in   ), pointer :: xmlElements
    character       (len=*   )              , intent(in   )          :: arrayElementName
    double precision          , dimension(:), intent(inout)          :: column1
    type            (node    )                             , pointer :: arrayElement    , xmlElement
    double precision          , dimension(1)                         :: dataValues
    integer                                                          :: i

    do i=1,getLength(xmlElements)
       xmlElement   => item                             (xmlElements,i-1             )
       arrayElement => XML_Get_First_Element_By_Tag_Name(xmlElement ,arrayElementName)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_List_Double_Array_Read_Static_One_Column

  subroutine XML_List_Integer_Array_Read_Static_One_Column(xmlElements,arrayElementName,column1)
    !% Read one column of integer data from an array of XML elements.
    implicit none
    type     (nodeList)              , intent(in   ), pointer :: xmlElements
    character(len=*   )              , intent(in   )          :: arrayElementName
    integer            , dimension(:), intent(inout)          :: column1
    type     (node    )                             , pointer :: arrayElement    , xmlElement
    integer            , dimension(1)                         :: dataValues
    integer                                                   :: i

    do i=1,getLength(xmlElements)
       xmlElement   => item                             (xmlElements,i-1             )
       arrayElement => XML_Get_First_Element_By_Tag_Name(xmlElement ,arrayElementName)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_List_Integer_Array_Read_Static_One_Column

  subroutine XML_List_Character_Array_Read_Static_One_Column(xmlElements,arrayElementName,column1)
    !% Read one column of character data from an array of XML elements.
    implicit none
    type     (nodeList        )              , intent(in   ), pointer :: xmlElements
    character(len=*           )              , intent(in   )          :: arrayElementName
    character(len=*           ), dimension(:), intent(inout)          :: column1
    type     (node            )                             , pointer :: arrayElement    , xmlElement
    character(len=len(column1)), dimension(1)                         :: dataValues
    integer                                                           :: i

    do i=1,getLength(xmlElements)
       xmlElement   => item                             (xmlElements,i-1             )
       arrayElement => XML_Get_First_Element_By_Tag_Name(xmlElement ,arrayElementName)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_List_Character_Array_Read_Static_One_Column

  function XML_Get_First_Element_By_Tag_Name(xmlElement,tagName,directChildrenOnly)
    !% Return a pointer to the first node in an XML node that matches the given {\normalfont \ttfamily tagName}.
    use Galacticus_Error
    implicit none
    type     (node            )               , pointer  :: XML_Get_First_Element_By_Tag_Name
    type     (node            ), intent(in   ), pointer  :: xmlElement
    character(len=*           ), intent(in   )           :: tagName
    logical                    , intent(in   ), optional :: directChildrenOnly
    type     (nodeList        )               , pointer  :: elementList
    type     (node            )               , pointer  :: parent
    character(len=len(tagName))                          :: currentTagName                   , path
    integer                                              :: pathPosition                     , i
    logical                                              :: directChildrenOnlyActual
    
    ! Set default options.
    directChildrenOnlyActual=.false.
    if (present(directChildrenOnly)) directChildrenOnlyActual=directChildrenOnly
    ! Find element.
    XML_Get_First_Element_By_Tag_Name => xmlElement
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
       elementList => getElementsByTagName(XML_Get_First_Element_By_Tag_Name,currentTagName)
       if (getLength(elementList) < 1) then
          call Galacticus_Error_Report('no elements match tag name "'//trim(currentTagName)//'"'//{introspection:location})
       else
          if (directChildrenOnlyActual) then
             do i=0,getLength(elementList)-1
                parent => getParentNode(item(elementList,i))
                if (associated(parent,XML_Get_First_Element_By_Tag_Name)) then
                   XML_Get_First_Element_By_Tag_Name => item(elementList,i)
                   exit
                end if
             end do
          else
             XML_Get_First_Element_By_Tag_Name => item(elementList,0)
          end if
       end if
    end do
    return
  end function XML_Get_First_Element_By_Tag_Name

  logical function XML_Path_Exists(xmlElement,path)
    !% Return true if the supplied {\normalfont \ttfamily path} exists in the supplied {\normalfont \ttfamily xmlElement}.
    implicit none
    type     (node         ), intent(in   ), pointer :: xmlElement
    character(len=*        ), intent(in   )          :: path
    type     (nodeList     )               , pointer :: elementList
    type     (node         )               , pointer :: element       , child         , &
         &                                              parent
    character(len=len(path))                         :: currentPath   , currentTagName
    integer                                          :: pathPosition  , i

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
       elementList => getElementsByTagName(element,currentTagName)
       if (getLength(elementList) < 1) then
          XML_Path_Exists=.false.
          return
       else
          XML_Path_Exists=.false.
          do i=0,getLength(elementList)-1
             child  => item         (elementList,i)
             parent => getParentNode(child        )
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
    !% Extracts information from a standard XML {\normalfont \ttfamily extrapolationElement}. Optionally a set of {\normalfont \ttfamily allowedMethods} can be
    !% specified---if the extracted method does not match one of these an error is issued.
    use Galacticus_Error
    use Table_Labels
    implicit none
    type     (Node    )              , intent(in   ), pointer  :: extrapolationElement
    character(len=*   )              , intent(  out)           :: limitType
    integer                          , intent(  out)           :: extrapolationMethod
    integer            , dimension(:), intent(in   ), optional :: allowedMethods
    type     (Node    )                             , pointer  :: limitElement        , methodElement
    type     (NodeList)                             , pointer  :: elementList
    character(len=32  )                                        :: methodType

    ! Extract the limit type.
    elementList => getElementsByTagname(extrapolationElement,"limit")
    if (getLength(elementList) /= 1) call Galacticus_Error_Report('extrapolation element must contain exactly one limit element'//{introspection:location})
    limitElement => item(elementList,0)
    call extractDataContent(limitElement,limitType)

    ! Extract the method type.
    elementList => getElementsByTagname(extrapolationElement,"method")
    if (getLength(elementList) /= 1) call Galacticus_Error_Report('extrapolation element must contain exactly one method element'//{introspection:location})
    methodElement => item(elementList,0)
    call extractDataContent(methodElement,methodType)
    extrapolationMethod=enumerationExtrapolationTypeEncode(trim(methodType),includesPrefix=.false.)
    ! Validate the method type.
    if (present(allowedMethods)) then
       if (all(allowedMethods /= extrapolationMethod)) call Galacticus_Error_Report('unallowed extrapolation method'//{introspection:location})
    end if

    return
  end subroutine XML_Extrapolation_Element_Decode

  function XML_Parse(fileName,iostat) result(document)
    !% Parse an XML document, automatically resolve XInclude references.
    use Galacticus_Error
    use File_Utilities
    implicit none
    type     (node            ), pointer                     :: document           , nodeNew       , &
         &                                                      nodeCurrent        , nodeParent    , &
         &                                                      nodeXInclude       , nodeImported  , &
         &                                                      nodeInsert         , nodeNext
    character(len=*           ), intent(in   )               :: fileName
    integer                    , intent(inout), optional     :: iostat
    type     (nodeList        ), pointer                     :: nodesCurrent
    type     (xincludeNode    ), allocatable  , dimension(:) :: stack              , stackTmp
    type     (xincludeNodeList), allocatable  , dimension(:) :: stackList          , stackListTmp
    integer                    , parameter                   :: stackExpandCount=10
    integer                                                  :: stackCount         , stackListCount, &
         &                                                      i                  , countElements
    type     (varying_string  )                              :: filePath           , fileLeaf      , &
         &                                                      nameInsert

    ! Extract the path and leaf name to our document.
    filePath=File_Path(fileName)
    fileLeaf=File_Name(fileName)
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
    do while (stackCount > 0)
       ! Parse the document.
       nodeNew      => parseFile(char(filePath//stack(stackCount)%fileName    ),iostat=iostat)
       nodeParent   =>                          stack(stackCount)%nodeParent
       nodeXInclude =>                          stack(stackCount)%nodeXInclude
       if (stack(stackCount)%xPath == "") then
          nodeInsert => getDocumentElement(nodeNew)
          nameInsert =  ""
       else
          if (XML_Path_Exists(nodeNew,char(stack(stackCount)%xPath))) then
             nodeInsert => XML_Get_First_Element_By_Tag_Name(nodeNew,char(stack(stackCount)%xPath),directChildrenOnly=.true.)
             nameInsert =  getNodeName  (nodeInsert)
             nodeInsert => getParentNode(nodeInsert)
          else
             call Galacticus_Error_Report("XPath '"//stack(stackCount)%xPath//"' not found"//{introspection:location})
          end if
       end if
       if (present(iostat).and.iostat /= 0) return
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
                if (getNodeName(nodeCurrent) == nameInsert) countElements=countElements+1
                nodeCurrent => getNextSibling(nodeCurrent)
             end do
             i           =  0
             nodeCurrent => getFirstChild(nodeInsert)
             do while (associated(nodeCurrent))
                nodeNext => getNextSibling(nodeCurrent)
                if ((getNodeType(nodeCurrent) /= ELEMENT_NODE .and. i > 0 .and. i < countElements) .or. getNodeName(nodeCurrent) == nameInsert) then
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
       end if
       ! Search for any XIncludes - rescan the entire document as we can handle inserting only one xi:include element at a time.
       stackListCount=0
       if (hasChildNodes(document)) then
          stackListCount                       =  1
          stackList     (stackListCount)%nodes => getChildNodes(document)
       end if
       do while (stackListCount > 0)
          nodesCurrent   => stackList(stackListCount)%nodes
          stackListCount =  stackListCount-1
          do i=0,getLength(nodesCurrent)-1
             nodeCurrent => item(nodesCurrent,i)
             if (getNodeName(nodeCurrent) == "xi:include") then
                if (.not.hasAttribute(nodeCurrent,"href")) call Galacticus_Error_Report("missing 'href' in XInclude"//{introspection:location})
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
                      call Galacticus_Error_Report("malformed XPath in XPointer: '"//stack(stackCount)%xPath//"'"//{introspection:location})
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
                stackListCount                       =  stackListCount+1
                stackList     (stackListCount)%nodes => getChildNodes(nodeCurrent)
             end if
          end do
       end do
    end do
    return
  end function XML_Parse

end module IO_XML
