!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements construction of numerical ranges.

module Numerical_Ranges
  !% Implements construction of numerical ranges.
  implicit none
  private
  public :: Make_Range

  ! Parameters to specify type of range required.
  integer, public, parameter :: rangeTypeUndefined=-1, rangeTypeLinear=0, rangeTypeLogarithmic=1
  
contains

  recursive function Make_Range(rangeMinimum,rangeMaximum,rangeNumber,rangeType) result (rangeValues)
    !% Builds a numerical range between {\tt rangeMinimum} and {\tt rangeMaximum} using {\tt rangeNumber} points and spacing as
    !% specified by {\tt rangeType} (defaulting to linear spacing if no {\tt rangeType} is given).
    use Galacticus_Error
    implicit none
    double precision, intent(in)           :: rangeMinimum,rangeMaximum
    integer,          intent(in)           :: rangeNumber
    integer,          intent(in), optional :: rangeType
    double precision                       :: rangeValues(rangeNumber)
    integer                                :: rangeTypeActual,iRange

    ! Find what type of range is required.
    if (present(rangeType)) then
       rangeTypeActual=rangeType
    else
       rangeTypeActual=rangeTypeLinear
    end if
    ! Build the range.
    select case (rangeTypeActual)
    case (rangeTypeLinear)
       ! Check that the rangeNumber is valid.
       if (rangeNumber <= 1) call Galacticus_Error_Report('Make_Range','number of points in range must exceed 1')
       ! Build a linear range.
       forall(iRange=1:rangeNumber)
          rangeValues(iRange)=rangeMinimum+(rangeMaximum-rangeMinimum)*dble(iRange-1)/dble(rangeNumber-1)
       end forall
    case (rangeTypeLogarithmic)
       ! Call ourself with logged limits and then exponentiate the result.
       rangeValues=exp(Make_Range(log(rangeMinimum),log(rangeMaximum),rangeNumber,rangeTypeLinear))
    case default
       call Galacticus_Error_Report('Make_Range','range type is unrecognized')
    end select
    return
  end function Make_Range

end module Numerical_Ranges
