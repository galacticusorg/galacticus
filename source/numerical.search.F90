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
Contains a module which implements searching of ordered arrays.
!!}

! Add dependency on GSL library.
!; gsl

module Arrays_Search
  !!{
  Implements searching of ordered arrays.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t, c_double
  implicit none
  private
  public :: searchArray, searchArrayClosest, searchIndexed

  !![
  <generic identifier="Type">
   <instance label="integer8" intrinsic="integer(kind_int8)"  />
   <instance label="varstr"   intrinsic="type(varying_string)"/>
  </generic>
  !!]

  interface searchArray
     !!{
     Generic interface for array searching routines.
     !!}
     module procedure searchArrayDouble
     module procedure searchArray{Type¦label}
  end interface searchArray

  interface searchIndexed
     !!{
     Generic interface for array searching routines using indexing.
     !!}
     module procedure searchIndexedInteger8
  end interface searchIndexed

  interface
     function gsl_interp_bsearch(x_array,x,index_lo,index_hi) bind(c,name='gsl_interp_bsearch')
       !!{
       Template for the GSL binary search function.
       !!}
       import
       integer(c_size_t)                              :: gsl_interp_bsearch
       real   (c_double), dimension(*), intent(in   ) :: x_array
       real   (c_double), value                       :: x
       integer(c_size_t), value                       :: index_lo          , index_hi
     end function gsl_interp_bsearch
  end interface
  
contains

  function searchArrayDouble(arrayToSearch,valueToFind)
    !!{
    Searches an array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, for value, $v(=${\normalfont \ttfamily valueToFind}$)$,
    to find the index $i$ such that $x(i) \le v < x(i+1)$.
    !!}
    implicit none
    integer         (c_size_t)                              :: searchArrayDouble
    double precision          , dimension(:), intent(in   ) :: arrayToSearch
    double precision                        , intent(in   ) :: valueToFind
    
    ! Account for the 0/1 array indexing difference between Fortran and C here.
    searchArrayDouble=+GSL_Interp_BSearch(                                                           &
         &                                arrayToSearch                                            , &
         &                                valueToFind                                              , &
         &                                int(lbound(arrayToSearch,dim=1),kind=c_size_t)-1_c_size_t, &
         &                                int(ubound(arrayToSearch,dim=1),kind=c_size_t)-1_c_size_t  &
         &                               )                                                           &
         &              +1_c_size_t
    return
  end function searchArrayDouble

  function searchArray{Type¦label}(arrayToSearch,valueToFind)
    !!{
    Searches an array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, for value, $v(=${\normalfont \ttfamily valueToFind}$)$,
    to find the index $i$ such that $x(i) \le v < x(i+1)$.
    !!}
    {Type¦match¦^(varstr)$¦  use :: ISO_Varying_String, only : varying_string, operator(<), operator(>), operator(<=), operator(>=)¦}
    {Type¦match¦^(integer8)$¦use :: Kind_Numbers      , only : kind_int8¦}
    implicit none
    integer         (c_size_t)                              :: searchArray{Type¦label}
    {Type¦intrinsic}          , dimension(:), intent(in   ) :: arrayToSearch
    {Type¦intrinsic}                        , intent(in   ) :: valueToFind
    integer         (c_size_t)                              :: jLower                 , jMidpoint, &
         &                                                     jUpper
    logical                                                 :: isInside

    ! Check whether valueToFind is outside range of arrayToSearch().
    isInside=.true.
    if (arrayToSearch(size(arrayToSearch,kind=c_size_t)) >= arrayToSearch(1)) then ! arrayToSearch() is in ascending order.
       if      (valueToFind < arrayToSearch(1                                )) then
          isInside=.false.
          searchArray{Type¦label}=                                 0_c_size_t
       else if (valueToFind > arrayToSearch(size(arrayToSearch,kind=c_size_t))) then
          isInside=.false.
          searchArray{Type¦label}=size(arrayToSearch,kind=c_size_t)+1_c_size_t
       end if
    else ! arrayToSearch() is in descending order.
       if      (valueToFind > arrayToSearch(1                                )) then
          isInside=.false.
          searchArray{Type¦label}=                                 +0_c_size_t
       else if (valueToFind < arrayToSearch(size(arrayToSearch,kind=c_size_t))) then
          isInside=.false.
          searchArray{Type¦label}=size(arrayToSearch,kind=c_size_t)+1_c_size_t
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
       searchArray{Type¦label}=jLower
    end if
    return
  end function searchArray{Type¦label}

  function searchArrayClosest(arrayToSearch,valueToFind,tolerance,status)
    !!{
    Searches an array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, for the entry closest to value, $v(=${\normalfont
    \ttfamily valueToFind}$)$ and returns the index of that element in the array. Optionally, a tolerance may be specified
    within which the two values must match.
    !!}
    use :: Error               , only : Error_Report, errorStatusFail, errorStatusSuccess
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    integer         (c_size_t)                                        :: searchArrayClosest
    double precision          , dimension(:), intent(in   )           :: arrayToSearch
    double precision                        , intent(in   )           :: valueToFind
    double precision                        , intent(in   ), optional :: tolerance
    integer                                 , intent(  out), optional :: status

    if (present(status)) status=errorStatusSuccess
    ! For a single element array, just return the only option.
    if (size(arrayToSearch,dim=1,kind=c_size_t) <= 1) then
       searchArrayClosest=lbound(arrayToSearch,dim=1)
       return
    end if
    ! Search for the closest entry.
    if      (valueToFind < arrayToSearch(lbound(arrayToSearch,dim=1))) then
       searchArrayClosest=lbound(arrayToSearch,dim=1)
    else if (valueToFind > arrayToSearch(ubound(arrayToSearch,dim=1))) then
       searchArrayClosest=ubound(arrayToSearch,dim=1)
    else
       searchArrayClosest=+GSL_Interp_BSearch(                                                           &
            &                                 arrayToSearch                                            , &
            &                                 valueToFind                                              , &
            &                                 int(lbound(arrayToSearch,dim=1),kind=c_size_t)-1_c_size_t, &
            &                                 int(ubound(arrayToSearch,dim=1),kind=c_size_t)-1_c_size_t  &
            &                                )                                                           &
            &              +1_c_size_t
       if     (                                                               &
            &   abs(valueToFind-arrayToSearch(searchArrayClosest           )) &
            &  >=                                                             &
            &   abs(valueToFind-arrayToSearch(searchArrayClosest+1_c_size_t)) &
            & )                                                               &
            & searchArrayClosest=+searchArrayClosest                          &
            &                    +1_c_size_t
    end if
    if (present(tolerance)) then
       if (.not.Values_Agree(valueToFind,arrayToSearch(searchArrayClosest),absTol=tolerance)) then
          if (present(status)) then
             status=errorStatusFail
          else
             call Error_Report('no match found within specified tolerance'//{introspection:location})
          end if
       end if
    end if
    return
  end function searchArrayClosest

  function searchIndexedInteger8(arrayToSearch,arrayIndex,valueToFind)
    !!{
    Searches a long integer array, $x=(${\normalfont \ttfamily arrayToSearch}$)$, which is rank ordered when indexed by
    {\normalfont \ttfamily arrayIndex}, for value, $v(=${\normalfont \ttfamily valueToFind}$)$, to find the index $i$ such that
    $x(i) \le v < x(i+1)$.
    !!}
    use :: Kind_Numbers , only : kind_int8
    implicit none
    integer(kind=kind_int8)                              :: searchIndexedInteger8
    integer(kind=kind_int8), dimension(:), intent(in   ) :: arrayToSearch
    integer(kind=c_size_t ), dimension(:), intent(in   ) :: arrayIndex
    integer(kind=kind_int8)              , intent(in   ) :: valueToFind
    integer(kind=kind_int8)                              :: jLower               , jMidpoint, &
         &                                                  jUpper
    logical                                              :: isInside

    ! Check whether valueToFind is outside range of arrayToSearch().
    isInside=.true.
    if (arrayToSearch(arrayIndex(size(arrayToSearch,kind=c_size_t))) >= arrayToSearch(arrayIndex(1))) then ! arrayToSearch() is in ascending order.
       if      (valueToFind < arrayToSearch(arrayIndex(1                                ))) then
          isInside             =.false.
          searchIndexedInteger8=                                 +0_c_size_t
       else if (valueToFind > arrayToSearch(arrayIndex(size(arrayToSearch,kind=c_size_t)))) then
          isInside             =.false.
          searchIndexedInteger8=size(arrayToSearch,kind=c_size_t)+1_c_size_t
       end if
    else ! arrayToSearch() is in descending order.
       if      (valueToFind > arrayToSearch(arrayIndex(1                                ))) then
          isInside             =.false.
          searchIndexedInteger8=                                 +0_c_size_t
       else if (valueToFind < arrayToSearch(arrayIndex(size(arrayToSearch,kind=c_size_t)))) then
          isInside             =.false.
          searchIndexedInteger8=size(arrayToSearch,kind=c_size_t)+1_c_size_t
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
       searchIndexedInteger8=jLower
    end if
    return
  end function searchIndexedInteger8

end module Arrays_Search

