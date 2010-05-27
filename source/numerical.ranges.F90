!% Contains a module which implements construction of numerical ranges.

module Numerical_Ranges
  !% Implements construction of numerical ranges.
  private
  public :: Make_Range

  ! Parameters to specify type of range required.
  integer, public, parameter :: rangeTypeLinear=0, rangeTypeLogarithmic=1
  
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
