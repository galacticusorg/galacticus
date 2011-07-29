!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements comparisons of values.

module Numerical_Comparison
  !% Implements comparisons of values.
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
