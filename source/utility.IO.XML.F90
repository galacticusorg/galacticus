!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements various utility functions for extracting data from XML files.

module IO_XML
  !% Implements various utility functions for extracting data from XML files.
  implicit none
  private
  public :: XML_Extrapolation_Element_Decode

  ! Labels for extrapolation methods.
  integer, parameter, public :: extrapolateZero=0, extrapolateFixed=1, extrapolatePowerLaw=2

contains
  
  subroutine XML_Extrapolation_Element_Decode(extrapolationElement,limitType,extrapolationMethod,allowedMethods)
    !% Extracts information from a standard XML {\tt extrapolationElement}. Optionally a set of {\tt allowedMethods} can be
    !% specified---if the extracted method does not match one of these an error is issued.
    use Galacticus_Error
    use FoX_dom
    implicit none
    type(Node),       pointer, intent(in)               :: extrapolationElement
    character(len=*),          intent(out)              :: limitType
    integer,                   intent(out)              :: extrapolationMethod
    integer,          optional,intent(in), dimension(:) :: allowedMethods
    type(Node),       pointer                           :: limitElement,methodElement
    type(NodeList),   pointer                           :: elementList
    character(len=32)                                   :: methodType

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
