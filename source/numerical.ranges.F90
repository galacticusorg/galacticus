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
Contains a module which implements construction of numerical ranges.
!!}

module Numerical_Ranges
  !!{
  Implements construction of numerical ranges.
  !!}
  implicit none
  private
  public :: Make_Range

  ! Parameters to specify type of range required.
  integer, parameter, public :: rangeTypeLinear=0, rangeTypeLogarithmic=1, rangeTypeUndefined=-1

contains

  recursive function Make_Range(rangeMinimum,rangeMaximum,rangeNumber,rangeType,rangeBinned) result (rangeValues)
    !!{
    Builds a numerical range between {\normalfont \ttfamily rangeMinimum} and {\normalfont
    \ttfamily rangeMaximum} using {\normalfont \ttfamily rangeNumber} points and spacing as
    specified by {\normalfont \ttfamily rangeType} (defaulting to linear spacing if no
    {\normalfont \ttfamily rangeType} is given).
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision, intent(in   )           :: rangeMaximum                  , rangeMinimum
    integer         , intent(in   )           :: rangeNumber
    integer         , intent(in   ), optional :: rangeType
    logical         , intent(in   ), optional :: rangeBinned
    double precision                          :: rangeValues      (rangeNumber)
    integer                                   :: iRange                        , rangeTypeActual
    logical                                   :: rangeBinnedActual

    ! Find what type of range is required.
    if (present(rangeType)) then
       rangeTypeActual=rangeType
    else
       rangeTypeActual=rangeTypeLinear
    end if
    ! Determine if a binned range is required.
    if (present(rangeBinned)) then
       rangeBinnedActual=rangeBinned
    else
       rangeBinnedActual=.false.
    end if
    ! Build the range.
    select case (rangeTypeActual)
    case (rangeTypeLinear)
       ! Determine if we are returning a regular range or bin centers.
       if (rangeBinnedActual) then
          ! Build a linear binned range returning bin centers.
          forall(iRange=1:rangeNumber)
             rangeValues(iRange)=rangeMinimum+(rangeMaximum-rangeMinimum)*(dble(iRange-1)+0.5d0)/dble(rangeNumber)
          end forall
       else
          ! Check that the rangeNumber is valid.
          if (rangeNumber <= 1) call Error_Report('number of points in range must exceed 1'//{introspection:location})
          ! Build a linear range.
          forall(iRange=1:rangeNumber)
             rangeValues(iRange)=rangeMinimum+(rangeMaximum-rangeMinimum)*dble(iRange-1)/dble(rangeNumber-1)
          end forall
          ! Ensure upper limit is precise.
          rangeValues(rangeNumber)=rangeMaximum
       end if
    case (rangeTypeLogarithmic)
       ! Call ourself with logged limits and then exponentiate the result.
       rangeValues=exp(Make_Range(log(rangeMinimum),log(rangeMaximum),rangeNumber,rangeTypeLinear,rangeBinnedActual))
    case default
       call Error_Report('range type is unrecognized'//{introspection:location})
    end select
    return
  end function Make_Range

end module Numerical_Ranges
