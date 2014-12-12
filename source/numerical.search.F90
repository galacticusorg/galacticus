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

!% Contains a module which implements searching of ordered arrays.

module Arrays_Search
  !% Implements searching of ordered arrays.
  use, intrinsic :: ISO_C_Binding
  use Kind_Numbers
  implicit none
  private
  public :: Search_Array, Search_Array_For_Closest, Search_Indexed

  interface Search_Array
     !% Generic interface for array searching routines.
     module procedure Search_Array_Double
     module procedure Search_Array_VarString
     module procedure Search_Array_Integer8
  end interface Search_Array

  interface Search_Indexed
     !% Generic interface for array searching routines using indexing.
     module procedure Search_Indexed_Integer8
  end interface Search_Indexed

contains

  function Search_Array_Double(arrayToSearch,valueToFind)
    !% Searches an array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, for value, $v(=${\normalfont \ttfamily valueToFind}$)$, to find the index $i$ such that $x(i) \le v < x(i+1)$.
    use FGSL
    implicit none
    integer         (c_size_t)                              :: Search_Array_Double
    double precision          , dimension(:), intent(in   ) :: arrayToSearch
    double precision                        , intent(in   ) :: valueToFind

    ! Call the FGSL routine to do the search.
    Search_Array_Double=FGSL_Interp_BSearch(arrayToSearch,valueToFind,int(lbound(arrayToSearch,dim=1),kind=fgsl_size_t)-1,int(ubound(arrayToSearch,dim=1),kind=fgsl_size_t)-1)
    return
  end function Search_Array_Double

   function Search_Array_Integer8(arrayToSearch,valueToFind)
    !% Searches a long integer array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, for value, $v(=${\normalfont \ttfamily valueToFind}$)$, to find the index $i$ such that $x(i) \le v < x(i+1)$.
    use Kind_Numbers
    implicit none
    integer(c_size_t      )                              :: Search_Array_Integer8
    integer(kind=kind_int8), dimension(:), intent(in   ) :: arrayToSearch
    integer(kind=kind_int8)              , intent(in   ) :: valueToFind
    integer(c_size_t      )                              :: jLower       , jMidpoint, jUpper
    logical                                              :: isInside

    isInside=.true.
    ! Check whether valueToFind is outside range of arrayToSearch().
    if (arrayToSearch(size(arrayToSearch,kind=c_size_t)) >= arrayToSearch(1)) then ! arrayToSearch() is in ascending order.
       if      (valueToFind < arrayToSearch(1                  )) then
          isInside=.false.
          Search_Array_Integer8=0
       else if (valueToFind > arrayToSearch(size(arrayToSearch,kind=c_size_t))) then
          isInside=.false.
          Search_Array_Integer8=size(arrayToSearch,kind=c_size_t)+1
       end if
    else ! arrayToSearch() is in descending order.
       if      (valueToFind > arrayToSearch(1                  )) then
          isInside=.false.
          Search_Array_Integer8=0
       else if (valueToFind < arrayToSearch(size(arrayToSearch,kind=c_size_t))) then
          isInside=.false.
          Search_Array_Integer8=size(arrayToSearch,kind=c_size_t)+1
       end if
    end if
    ! Binary search in array if valueToFind lies within array range.
    if (isInside) then
       jLower=0
       jUpper=size(arrayToSearch,kind=c_size_t)+1
       do while (jUpper-jLower > 1)
          jMidpoint=(jUpper+jLower)/2
          if ((arrayToSearch(size(arrayToSearch,kind=c_size_t)) >= arrayToSearch(1)) .eqv. (valueToFind >= arrayToSearch(jMidpoint))) then
             jLower=jMidpoint
          else
             jUpper=jMidpoint
          endif
       end do
       Search_Array_Integer8=jLower
    end if
    return
  end function Search_Array_Integer8

  function Search_Array_VarString(arrayToSearch,valueToFind)
    !% Searches an array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, for value, $v(=${\normalfont \ttfamily valueToFind}$)$, to find the index $i$ such that $x(i)
    !% = v$. With this algorithm, if multiple elements of $x()$ have the same value, then the largest value of $i$ for which
    !% $x(i)=v$ occurs will be returned.
    use ISO_Varying_String
    implicit none
    integer(c_size_t      )                              :: Search_Array_VarString
    type   (varying_string), dimension(:), intent(in   ) :: arrayToSearch
    type   (varying_string)              , intent(in   ) :: valueToFind
    integer(c_size_t      )                              :: jLower       , jMidpoint, jUpper
    logical                                              :: isInside

    isInside=.true.
    ! Check whether valueToFind is outside range of arrayToSearch().
    if (arrayToSearch(size(arrayToSearch,kind=c_size_t)) >= arrayToSearch(1)) then ! arrayToSearch() is in ascending order.
       if      (valueToFind < arrayToSearch(1                  )) then
          isInside=.false.
          Search_Array_VarString=0
       else if (valueToFind > arrayToSearch(size(arrayToSearch,kind=c_size_t))) then
          isInside=.false.
          Search_Array_VarString=size(arrayToSearch,kind=c_size_t)+1
       end if
    else ! arrayToSearch() is in descending order.
       if      (valueToFind > arrayToSearch(1                  )) then
          isInside=.false.
          Search_Array_VarString=0
       else if (valueToFind < arrayToSearch(size(arrayToSearch,kind=c_size_t))) then
          isInside=.false.
          Search_Array_VarString=size(arrayToSearch,kind=c_size_t)+1
       end if
    end if
    ! Binary search in array if valueToFind lies within array range.
    if (isInside) then
       jLower=0
       jUpper=size(arrayToSearch,kind=c_size_t)+1
       do while (jUpper-jLower > 1)
          jMidpoint=(jUpper+jLower)/2
          if ((arrayToSearch(size(arrayToSearch,kind=c_size_t)) >= arrayToSearch(1)) .eqv. (valueToFind >= arrayToSearch(jMidpoint))) then
             jLower=jMidpoint
          else
             jUpper=jMidpoint
          endif
       end do
       Search_Array_VarString=jLower
    end if
    return
  end function Search_Array_VarString

   function Search_Array_For_Closest(arrayToSearch,valueToFind,tolerance)
    !% Searches an array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, for the entry closest to value, $v(=${\normalfont \ttfamily valueToFind}$)$ and returns the
    !% index of that element in the array. Optionally, a tolerance may be specified within which the two values must match.
    use FGSL
    use Galacticus_Error
    use Numerical_Comparison
    implicit none
    integer         (c_size_t)                                        :: Search_Array_For_Closest
    double precision          , dimension(:), intent(in   )           :: arrayToSearch
    double precision                        , intent(in   )           :: valueToFind
    double precision                        , intent(in   ), optional :: tolerance

    ! For a single element array, just return the only option.
    if (size(arrayToSearch,dim=1,kind=c_size_t) <= 1) then
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
       if (abs(valueToFind-arrayToSearch(Search_Array_For_Closest)) >= abs(valueToFind-arrayToSearch(Search_Array_For_Closest&
            &+1))) Search_Array_For_Closest=Search_Array_For_Closest+1
    end if
    if (present(tolerance)) then
       if (.not.Values_Agree(valueToFind,arrayToSearch(Search_Array_For_Closest),absTol=tolerance)) &
            & call Galacticus_Error_Report('Search_Array_For_Closest','no match found within specified tolerance')
    end if
    return
  end function Search_Array_For_Closest

  function Search_Indexed_Integer8(arrayToSearch,arrayIndex,valueToFind)
    !% Searches a long integer array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, which is rank ordered when indexed by {\normalfont \ttfamily arrayIndex}, for value, $v(=${\normalfont \ttfamily valueToFind}$)$, to find the index $i$ such that $x(i) \le v < x(i+1)$.
    use Kind_Numbers
    use, intrinsic :: ISO_C_Binding
    implicit none
    integer(kind=kind_int8) :: Search_Indexed_Integer8
    integer(kind=kind_int8), dimension(:), intent(in   ) :: arrayToSearch
    integer(kind=c_size_t ), dimension(:), intent(in   ) :: arrayIndex
    integer(kind=kind_int8)              , intent(in   ) :: valueToFind
    integer(kind=kind_int8)                              :: jLower       , jMidpoint, jUpper
    logical                                              :: isInside

    isInside=.true.
    ! Check whether valueToFind is outside range of arrayToSearch().


if (arrayIndex(1) <= 0 .or. arrayIndex(1) > size(arrayToSearch,kind=c_size_t)) then
write (0,*) "problem #1 ",arrayIndex(1),size(arrayToSearch,kind=c_size_t)
write (404,*) "problem #1 ",arrayIndex(1),size(arrayToSearch,kind=c_size_t)
call flush(0)
call flush(404)
end if 
if (size(arrayToSearch,kind=c_size_t) > size(arrayIndex,kind=c_size_t)) then
write (0,*) "problem #2 ",size(arrayToSearch,kind=c_size_t),size(arrayIndex,kind=c_size_t)
write (404,*) "problem #2 ",size(arrayToSearch,kind=c_size_t),size(arrayIndex,kind=c_size_t)
call flush(0)
call flush(404)
end if
if (arrayIndex(size(arrayToSearch,kind=c_size_t)) <= 0 .or. arrayIndex(size(arrayToSearch,kind=c_size_t)) > size(arrayToSearch,kind=c_size_t)) then
write (0,*) "problem #3 ",arrayIndex(size(arrayToSearch,kind=c_size_t)),size(arrayToSearch,kind=c_size_t)
write (404,*) "problem #3 ",arrayIndex(size(arrayToSearch,kind=c_size_t)),size(arrayToSearch,kind=c_size_t)
call flush(0)
call flush(404)
end if

    if (arrayToSearch(arrayIndex(size(arrayToSearch,kind=c_size_t))) >= arrayToSearch(arrayIndex(1))) then ! arrayToSearch() is in ascending order.
       if      (valueToFind < arrayToSearch(arrayIndex(1                  ))) then
          isInside=.false.
          Search_Indexed_Integer8=0
       else if (valueToFind > arrayToSearch(arrayIndex(size(arrayToSearch,kind=c_size_t)))) then
          isInside=.false.
          Search_Indexed_Integer8=size(arrayToSearch,kind=c_size_t)+1
       end if
    else ! arrayToSearch() is in descending order.
       if      (valueToFind > arrayToSearch(arrayIndex(1                  ))) then
          isInside=.false.
          Search_Indexed_Integer8=0
       else if (valueToFind < arrayToSearch(arrayIndex(size(arrayToSearch,kind=c_size_t)))) then
          isInside=.false.
          Search_Indexed_Integer8=size(arrayToSearch,kind=c_size_t)+1
       end if
    end if
    ! Binary search in array if valueToFind lies within array range.
    if (isInside) then
       jLower=0
       jUpper=size(arrayToSearch,kind=c_size_t)+1
       do while (jUpper-jLower > 1)
          jMidpoint=(jUpper+jLower)/2
          if ((arrayToSearch(arrayIndex(size(arrayToSearch,kind=c_size_t))) >= arrayToSearch(arrayIndex(1))) .eqv. (valueToFind >= arrayToSearch(arrayIndex(jMidpoint)))) then
             jLower=jMidpoint
          else
             jUpper=jMidpoint
          endif
       end do
       Search_Indexed_Integer8=jLower
    end if
    return
  end function Search_Indexed_Integer8

end module Arrays_Search

