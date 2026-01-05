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
Contains a module which implements topological sorting.
!!}

module Sorting_Topological
  !!{
  Implements topological sorting.
  !!}
  implicit none
  private
  public :: Sort_Topological

contains

  subroutine Sort_Topological(countObjects,countDependencies,dependencies,order,countOrdered,status)
    !!{
    Topological sorting function. Based on the example from \href{https://www.rosettacode.org/wiki/Topological_sort\#Modern_Fortran}{Rosetta Code}. Arguments are:
    \begin{description}
     \item[{\normalfont \ttfamily countObjects}\argin] the number of objects to be sorted;
     \item[{\normalfont \ttfamily countDependencies}\argin] the number of dependencies;
     \item[{\normalfont \ttfamily dependencies}\argin] an array of dependencies, such that {\normalfont \ttfamily dependencies(:,1)} depends on {\normalfont \ttfamily dependencies(:,2)};
     \item[{\normalfont \ttfamily order}\argout] an array giving the order of the objects after sorting;
     \item[{\normalfont \ttfamily countOrdered}\argout] a count of the objects which were ordered by the sort, such that {\normalfont \ttfamily order(1:countOrdered)} contains the ordered objects, while the remainder of {\normalfont \ttfamily order()} contains objects that were unordered (i.e. had no dependencies).
    \end{description}
    The unordered objects are those for which no solution is available---i.e. the graph is not acyclic. So, if {\normalfont \ttfamily order}$<${\normalfont \ttfamily countObjects} then one or more circular dependencies existed in the graph.
    !!}
    use :: Error, only : Error_Report, errorStatusFail, errorStatusSuccess
    implicit none
    integer                                , intent(in   ) :: countObjects
    integer                                , intent(in   ) :: countDependencies
    integer, dimension(countDependencies,2), intent(in   ) :: dependencies
    integer, dimension(countObjects       ), intent(  out) :: order
    integer                                , intent(  out) :: countOrdered
    integer, optional                      , intent(  out) :: status
    integer                                                :: dependencyLeft   , dependencyRight, &
         &                                                    positionLeft     , positionRight  , &
         &                                                    i                , j              , &
         &                                                    k
    integer, dimension(countObjects       )                :: position

    do i=1,countObjects
       order   (i)=i
       position(i)=i
    end do
    k=1
    do
       j=k
       k=countObjects+1
       do i=1,countDependencies
          dependencyLeft =dependencies(i,1)
          dependencyRight=dependencies(i,2)
          positionLeft   =position(dependencyLeft )
          positionRight  =position(dependencyRight)
          if (dependencyLeft == dependencyRight .or. positionLeft >= k .or. positionLeft < j .or. positionRight < j) cycle
          k=k-1
          position(order(k)      )=  positionLeft
          position(dependencyLeft)=      k
          order   (  positionLeft)=order(k)
          order   (      k       )=dependencyLeft
       end do
       if (k <= j) exit
    end do
    countOrdered=j-1
    if (present(status)) status=errorStatusSuccess
    if (countOrdered < countObjects) then
       if (present(status)) then
          status=errorStatusFail
       else
          call Error_Report('circular dependencies'//{introspection:location})
       end if
    end if
    return
  end subroutine Sort_Topological

end module Sorting_Topological
