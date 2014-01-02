!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  implicit none
  private
  public :: XML_Extrapolation_Element_Decode , XML_Array_Read  , XML_Array_Read_Static, &
       &    XML_Get_First_Element_By_Tag_Name, XML_Array_Length, XML_Path_Exists

  ! Labels for extrapolation methods.
  integer, parameter, public :: extrapolateFixed=1, extrapolatePowerLaw=2, extrapolateZero=0

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

contains

  integer function XML_Array_Length(xmlElement,arrayElementName)
    !% Return the length of an array of XML elements.
    use FoX_dom
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
    use FoX_dom
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
    use FoX_dom
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
    call Alloc_Array(column1,[getLength(arrayElements)])
    do i=1,getLength(arrayElements)
       arrayElement => item(arrayElements,i-1)
       call extractDataContent(arrayELement,dataValues)
       column1(i)=dataValues(1)
    end do
    return
  end subroutine XML_Array_Read_One_Column

  subroutine XML_Array_Read_Two_Column(xmlElement,arrayElementName,column1,column2)
    !% Read two columns of data from an array of XML elements.
    use FoX_dom
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
    call Alloc_Array(column1,[getLength(arrayElements)])
    call Alloc_Array(column2,[getLength(arrayElements)])
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
    use FoX_dom
    use Memory_Management
    implicit none
    type            (nodeList)                           , intent(in   ), pointer :: xmlElements
    character       (len=*   )                           , intent(in   )          :: arrayElementName
    double precision          , allocatable, dimension(:), intent(inout)          :: column1
    type            (node    )                                          , pointer :: arrayElement    , xmlElement
    double precision                       , dimension(1)                         :: dataValues
    integer                                                                       :: i

    call Alloc_Array(column1,[getLength(xmlElements)])
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
    use FoX_dom
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
    use FoX_dom
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
    use FoX_dom
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

  function XML_Get_First_Element_By_Tag_Name(xmlElement,tagName)
    !% Return a pointer to the first node in an XML node that matches the given {\tt tagName}.
    use FoX_dom
    use Galacticus_Error
    implicit none
    type     (node            )               , pointer :: XML_Get_First_Element_By_Tag_Name
    type     (node            ), intent(in   ), pointer :: xmlElement
    character(len=*           ), intent(in   )          :: tagName
    type     (nodeList        )               , pointer :: elementList
    character(len=len(tagName))                         :: currentTagName                   , path
    integer                                             :: pathPosition

    XML_Get_First_Element_By_Tag_Name => xmlElement
    path=tagName
    do while (path /= "")
       pathPosition=index(path,"/")
       if (pathPosition == 0) then
          currentTagName=path
          path          =""
       else
          currentTagName=path(             1:          pathPosition-1)
          path          =path(pathPosition+1:len(path)-pathPosition  )
       endif
       elementList => getElementsByTagName(XML_Get_First_Element_By_Tag_Name,currentTagName)
       if (getLength(elementList) < 1) then
          call Galacticus_Error_Report('XML_Get_First_Element_By_Tag_Name','no elements match tag name')
       else
          XML_Get_First_Element_By_Tag_Name => item(elementList,0)
       end if
    end do
    return
  end function XML_Get_First_Element_By_Tag_Name

  logical function XML_Path_Exists(xmlElement,path)
    !% Return true if the supplied {\tt path} exists in the supplied {\tt xmlElement}.
    use FoX_dom
    implicit none
    type     (node         ), intent(in   ), pointer :: xmlElement
    character(len=*        ), intent(in   )          :: path
    type     (nodeList     )               , pointer :: elementList
    type     (node         )               , pointer :: element
    character(len=len(path))                         :: currentPath , currentTagName
    integer                                          :: pathPosition

    XML_Path_Exists =  .true.
    element         => xmlElement
    currentPath     =  path
    do while (currentPath /= "")
       pathPosition=index(currentPath,"/")
       if (pathPosition == 0) then
          currentTagName=currentPath
          currentPath   =""
       else
          currentTagName=currentPath(             1:          pathPosition-1)
          currentPath   =currentPath(pathPosition+1:len(path)-pathPosition  )
       endif
       elementList => getElementsByTagName(element,currentTagName)
       if (getLength(elementList) < 1) then
          XML_Path_Exists=.false.
          return
       else
          element => item(elementList,0)
       end if
    end do
    return
  end function XML_Path_Exists

  subroutine XML_Extrapolation_Element_Decode(extrapolationElement,limitType,extrapolationMethod,allowedMethods)
    !% Extracts information from a standard XML {\tt extrapolationElement}. Optionally a set of {\tt allowedMethods} can be
    !% specified---if the extracted method does not match one of these an error is issued.
    use Galacticus_Error
    use FoX_dom
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
    if (getLength(elementList) /= 1) call Galacticus_Error_Report('Extrapolation_Element_Decode','extrapolation element must contain exactly one limit element')
    limitElement => item(elementList,0)
    call extractDataContent(limitElement,limitType)

    ! Extract the method type.
    elementList => getElementsByTagname(extrapolationElement,"method")
    if (getLength(elementList) /= 1) call Galacticus_Error_Report('Extrapolation_Element_Decode','extrapolation element must contain exactly one method element')
    methodElement => item(elementList,0)
    call extractDataContent(methodElement,methodType)
    select case (trim(methodType))
    case ('zero')
       extrapolationMethod=extrapolateZero
    case ('fixed')
       extrapolationMethod=extrapolateFixed
    case ('power law')
       extrapolationMethod=extrapolatePowerLaw
    case default
       call Galacticus_Error_Report('Extrapolation_Element_Decode','unrecognized extrapolation method')
    end select

    ! Validate the method type.
    if (present(allowedMethods)) then
       if (all(allowedMethods /= extrapolationMethod)) call Galacticus_Error_Report('Extrapolation_Element_Decode','unallowed extrapolation method')
    end if

    return
  end subroutine XML_Extrapolation_Element_Decode

end module IO_XML
