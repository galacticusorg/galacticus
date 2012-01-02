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


!% Contains a module which implements searching of ordered arrays.

module Arrays_Search
  !% Implements searching of ordered arrays.
  implicit none
  private
  public :: Search_Array, Search_Array_For_Closest

  interface Search_Array
     !% Generic interface for array searching routines.
     module procedure Search_Array_Double
     module procedure Search_Array_VarString
     module procedure Search_Array_Integer8
  end interface Search_Array

contains

  integer function Search_Array_Double(arrayToSearch,valueToFind)
    !% Searches an array, $x=(${\tt arrayToSearch}$)$, for value, $v(=${\tt valueToFind}$)$, to find the index $i$ such that $x(i) \le v < x(i+1)$.
    use FGSL
    implicit none
    double precision, intent(in), dimension(:) :: arrayToSearch
    double precision, intent(in)               :: valueToFind
    
    ! Call the FGSL routine to do the search.
    Search_Array_Double=FGSL_Interp_BSearch(arrayToSearch,valueToFind,int(lbound(arrayToSearch,dim=1),kind=fgsl_size_t)-1,int(ubound(arrayToSearch,dim=1),kind=fgsl_size_t)-1)
    return
  end function Search_Array_Double
  
  integer function Search_Array_Integer8(arrayToSearch,valueToFind)
    !% Searches a long integer array, $x=(${\tt arrayToSearch}$)$, for value, $v(=${\tt valueToFind}$)$, to find the index $i$ such that $x(i) \le v < x(i+1)$.
    use Kind_Numbers
    implicit none
    integer(kind=kind_int8), intent(in), dimension(:) :: arrayToSearch
    integer(kind=kind_int8), intent(in)               :: valueToFind
    integer                                           :: jLower,jMidpoint,jUpper
    logical                                           :: isInside

    isInside=.true.
    ! Check whether valueToFind is outside range of arrayToSearch().
    if (arrayToSearch(size(arrayToSearch)) >= arrayToSearch(1)) then ! arrayToSearch() is in ascending order.
       if      (valueToFind < arrayToSearch(1                  )) then
          isInside=.false.
          Search_Array_Integer8=0
       else if (valueToFind > arrayToSearch(size(arrayToSearch))) then
          isInside=.false.
          Search_Array_Integer8=size(arrayToSearch)+1
       end if
    else ! arrayToSearch() is in descending order.
       if      (valueToFind > arrayToSearch(1                  )) then
          isInside=.false.
          Search_Array_Integer8=0
       else if (valueToFind < arrayToSearch(size(arrayToSearch))) then
          isInside=.false.
          Search_Array_Integer8=size(arrayToSearch)+1
       end if
    end if
    ! Binary search in array if valueToFind lies within array range.
    if (isInside) then
       jLower=0
       jUpper=size(arrayToSearch)+1
       do while (jUpper-jLower > 1)
          jMidpoint=(jUpper+jLower)/2
          if ((arrayToSearch(size(arrayToSearch)) >= arrayToSearch(1)) .eqv. (valueToFind >= arrayToSearch(jMidpoint))) then
             jLower=jMidpoint
          else
             jUpper=jMidpoint
          endif
       end do
       Search_Array_Integer8=jLower
    end if
    return
  end function Search_Array_Integer8
  
  integer function Search_Array_VarString(arrayToSearch,valueToFind)
    !% Searches an array, $x=(${\tt arrayToSearch}$)$, for value, $v(=${\tt valueToFind}$)$, to find the index $i$ such that $x(i)
    !% = v$. With this algorithm, if multiple elements of $x()$ have the same value, then the largest value of $i$ for which
    !% $x(i)=v$ occurs will be returned.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(in), dimension(:) :: arrayToSearch
    type(varying_string), intent(in)               :: valueToFind   
    integer                                        :: jLower,jMidpoint,jUpper
    logical                                        :: isInside

    isInside=.true.
    ! Check whether valueToFind is outside range of arrayToSearch().
    if (arrayToSearch(size(arrayToSearch)) >= arrayToSearch(1)) then ! arrayToSearch() is in ascending order.
       if      (valueToFind < arrayToSearch(1                  )) then
          isInside=.false.
          Search_Array_VarString=0
       else if (valueToFind > arrayToSearch(size(arrayToSearch))) then
          isInside=.false.
          Search_Array_VarString=size(arrayToSearch)+1
       end if
    else ! arrayToSearch() is in descending order.
       if      (valueToFind > arrayToSearch(1                  )) then
          isInside=.false.
          Search_Array_VarString=0
       else if (valueToFind < arrayToSearch(size(arrayToSearch))) then
          isInside=.false.
          Search_Array_VarString=size(arrayToSearch)+1
       end if
    end if
    ! Binary search in array if valueToFind lies within array range.
    if (isInside) then
       jLower=0
       jUpper=size(arrayToSearch)+1
       do while (jUpper-jLower > 1)
          jMidpoint=(jUpper+jLower)/2
          if ((arrayToSearch(size(arrayToSearch)) >= arrayToSearch(1)) .eqv. (valueToFind >= arrayToSearch(jMidpoint))) then
             jLower=jMidpoint
          else
             jUpper=jMidpoint
          endif
       end do
       Search_Array_VarString=jLower
    end if
    return
  end function Search_Array_VarString
  
  integer function Search_Array_For_Closest(arrayToSearch,valueToFind,tolerance)
    !% Searches an array, $x=(${\tt arrayToSearch}$)$, for the entry closest to value, $v(=${\tt valueToFind}$)$ and returns the
    !% index of that element in the array. Optionally, a tolerance may be specified within which the two values must match.
    use FGSL
    use Galacticus_Error
    use Numerical_Comparison
    implicit none
    double precision, intent(in), dimension(:) :: arrayToSearch
    double precision, intent(in)               :: valueToFind
    double precision, intent(in), optional     :: tolerance
    
    ! For a single element array, just return the only option.
    if (size(arrayToSearch,dim=1) <= 1) then
       Search_Array_For_Closest=lbound(arrayToSearch,dim=1)
       return
    end if

    ! Call the FGSL routine to do the search.
    if      (valueToFind < arrayToSearch(lbound(arrayToSearch,dim=1))) then
       Search_Array_For_Closest=lbound(arrayToSearch,dim=1)
    else if (valueToFind > arrayToSearch(ubound(arrayToSearch,dim=1))) then
       Search_Array_For_Closest=ubound(arrayToSearch,dim=1)
    else
       Search_Array_For_Closest=FGSL_Interp_BSearch(arrayToSearch,valueToFind,int(lbound(arrayToSearch,dim=1),kind=fgsl_size_t)-1&
            &,int(ubound(arrayToSearch,dim=1),kind=fgsl_size_t)-1)
       if (dabs(valueToFind-arrayToSearch(Search_Array_For_Closest)) >= dabs(valueToFind-arrayToSearch(Search_Array_For_Closest&
            &+1))) Search_Array_For_Closest=Search_Array_For_Closest+1
    end if
    if (present(tolerance)) then
       if (.not.Values_Agree(valueToFind,arrayToSearch(Search_Array_For_Closest),absTol=tolerance)) &
            & call Galacticus_Error_Report('Search_Array_For_Closest','no match found within specified tolerance')
    end if
    return
  end function Search_Array_For_Closest
  
end module Arrays_Search

