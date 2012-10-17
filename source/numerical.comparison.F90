!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements comparisons of values.

module Numerical_Comparison
  !% Implements comparisons of values.
  implicit none
  private
  public :: Values_Differ, Values_Agree

  interface Values_Differ
     module procedure Values_Differ_Real
     module procedure Values_Differ_Double
  end interface Values_Differ

  interface Values_Agree
     module procedure Values_Agree_Real
     module procedure Values_Agree_Double
  end interface Values_Agree

contains

  logical function Values_Differ_Real(value1,value2,absTol,relTol)
    !% Returns true if {\tt value1} and {\tt value2} differ by more than {\tt absTol} in absolute terms, or {\tt relTol} in
    !% relative terms.
    implicit none
    real, intent(in)           :: value1,value2
    real, intent(in), optional :: absTol,relTol
    
    Values_Differ_Real=.false.
    if (present(absTol)) Values_Differ_Real=(abs(value1-value2) > absTol)
    if (present(relTol)) Values_Differ_Real=Values_Differ_Real.or.(abs(value1-value2) > 0.5d0*abs(value1+value2)*relTol)
    if (.not.(present(absTol).or.present(relTol))) Values_Differ_Real=(value1 /= value2)
    return
  end function Values_Differ_Real

  logical function Values_Differ_Double(value1,value2,absTol,relTol)
    !% Returns true if {\tt value1} and {\tt value2} differ by more than {\tt absTol} in absolute terms, or {\tt relTol} in
    !% relative terms.
    implicit none
    double precision, intent(in)           :: value1,value2
    double precision, intent(in), optional :: absTol,relTol
    
    Values_Differ_Double=.false.
    if (present(absTol)) Values_Differ_Double=(abs(value1-value2) > absTol)
    if (present(relTol)) Values_Differ_Double=Values_Differ_Double.or.(abs(value1-value2) > 0.5d0*abs(value1+value2)*relTol)
    if (.not.(present(absTol).or.present(relTol))) Values_Differ_Double=(value1 /= value2)
    return
  end function Values_Differ_Double

  logical function Values_Agree_Real(value1,value2,absTol,relTol)
    !% Returns true if {\tt value1} and {\tt value2} agree to within {\tt absTol} in absolute terms, or {\tt relTol} in
    !% relative terms.
    implicit none
    real,   intent(in)           :: value1,value2
    real,   intent(in), optional :: absTol,relTol
    logical                      :: agreeAbsolutely,agreeRelatively

    if (.not.(present(absTol).or.present(relTol))) then
       Values_Agree_Real=(value1 == value2)
       return
    end if
    if (present(absTol)) then
       agreeAbsolutely=(abs(value1-value2) <= absTol)
    else
       agreeAbsolutely=.true.
    end if
    if (present(relTol)) then
       agreeRelatively=(abs(value1-value2) <= 0.5d0*abs(value1+value2)*relTol)
    else
       agreeRelatively=.true.
    end if
    Values_Agree_Real=    (present(absTol).and.agreeAbsolutely) &
         &              .or.(present(relTol).and.agreeRelatively)
    return
  end function Values_Agree_Real
  
  logical function Values_Agree_Double(value1,value2,absTol,relTol)
    !% Returns true if {\tt value1} and {\tt value2} agree to within {\tt absTol} in absolute terms, or {\tt relTol} in
    !% relative terms.
    implicit none
    double precision, intent(in)           :: value1,value2
    double precision, intent(in), optional :: absTol,relTol
    logical                                :: agreeAbsolutely,agreeRelatively

    if (.not.(present(absTol).or.present(relTol))) then
       Values_Agree_Double=(value1 == value2)
       return
    end if
    if (present(absTol)) then
       agreeAbsolutely=(abs(value1-value2) <= absTol)
    else
       agreeAbsolutely=.true.
    end if
    if (present(relTol)) then
       agreeRelatively=(abs(value1-value2) <= 0.5d0*abs(value1+value2)*relTol)
    else
       agreeRelatively=.true.
    end if
    Values_Agree_Double=    (present(absTol).and.agreeAbsolutely) &
         &              .or.(present(relTol).and.agreeRelatively)
    return
  end function Values_Agree_Double
  
end module Numerical_Comparison
