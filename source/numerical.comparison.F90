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
Contains a module which implements comparisons of values.
!!}

module Numerical_Comparison
  !!{
  Implements comparisons of values.
  !!}
  implicit none
  private
  public :: Values_Differ, Values_Agree, Values_Less_Than

  interface Values_Differ
     module procedure Values_Differ_Real
     module procedure Values_Differ_Double
     module procedure Values_Differ_Double_Complex
  end interface Values_Differ

  interface Values_Agree
     module procedure Values_Agree_Real
     module procedure Values_Agree_Double
     module procedure Values_Agree_Double_Complex
  end interface Values_Agree

  interface Values_Less_Than
     module procedure Values_Less_Than_Double
  end interface Values_Less_Than
  
contains

  elemental logical function Values_Differ_Real(value1,value2,absTol,relTol)
    !!{
    Returns true if {\normalfont \ttfamily value1} and {\normalfont \ttfamily value2} differ by more than {\normalfont \ttfamily absTol} in absolute terms, or {\normalfont \ttfamily relTol} in
    relative terms.
    !!}
    implicit none
    real, intent(in   )           :: value1, value2
    real, intent(in   ), optional :: absTol, relTol

    Values_Differ_Real=.false.
    if (present(absTol)) Values_Differ_Real=(abs(value1-value2) > absTol)
    if (present(relTol)) Values_Differ_Real=Values_Differ_Real.or.(abs(value1-value2) > 0.5d0*abs(value1+value2)*relTol)
    if (.not.(present(absTol).or.present(relTol))) Values_Differ_Real=(value1 /= value2)
    return
  end function Values_Differ_Real

  elemental logical function Values_Differ_Double(value1,value2,absTol,relTol)
    !!{
    Returns true if {\normalfont \ttfamily value1} and {\normalfont \ttfamily value2} differ by more than {\normalfont \ttfamily absTol} in absolute terms, or {\normalfont \ttfamily relTol} in
    relative terms.
    !!}
    implicit none
    double precision, intent(in   )           :: value1, value2
    double precision, intent(in   ), optional :: absTol, relTol

    Values_Differ_Double=.false.
    if (present(absTol)) Values_Differ_Double=(abs(value1-value2) > absTol)
    if (present(relTol)) Values_Differ_Double=Values_Differ_Double.or.(abs(value1-value2) > 0.5d0*abs(value1+value2)*relTol)
    if (.not.(present(absTol).or.present(relTol))) Values_Differ_Double=(value1 /= value2)
    return
  end function Values_Differ_Double

  elemental logical function Values_Differ_Double_Complex(value1,value2,absTol,relTol)
    !!{
    Returns true if {\normalfont \ttfamily value1} and {\normalfont \ttfamily value2} differ by more than {\normalfont \ttfamily absTol} in absolute terms, or {\normalfont \ttfamily relTol} in
    relative terms.
    !!}
    implicit none
    double complex, intent(in   )           :: value1, value2
    double complex, intent(in   ), optional :: absTol, relTol

    Values_Differ_Double_Complex=.false.
    if (      present(absTol)                    ) Values_Differ_Double_Complex= abs(real(value1  - value2)) >                                real(absTol) &
         &                                                                      .or.                                                                       &
         &                                                                       abs(imag(value1  - value2)) >                                imag(absTol)
    if (                         present(relTol) ) Values_Differ_Double_Complex= Values_Differ_Double_Complex                                              &
         &                                                                      .or.                                                                       &
         &                                                                       abs(real(value1  - value2)) > 0.5d0*abs(real(value1+value2))*real(relTol) &
         &                                                                      .or.                                                                       &
         &                                                                       abs(imag(value1  - value2)) > 0.5d0*abs(imag(value1+value2))*imag(relTol)
    if (.not.(present(absTol).or.present(relTol))) Values_Differ_Double_Complex=          value1 /= value2
    return
  end function Values_Differ_Double_Complex

  elemental logical function Values_Agree_Real(value1,value2,absTol,relTol)
    !!{
    Returns true if {\normalfont \ttfamily value1} and {\normalfont \ttfamily value2} agree to within {\normalfont \ttfamily absTol} in absolute terms, or {\normalfont \ttfamily relTol} in
    relative terms.
    !!}
    implicit none
    real   , intent(in   )           :: value1         , value2
    real   , intent(in   ), optional :: absTol         , relTol
    logical                          :: agreeAbsolutely, agreeRelatively

    if (value1 == value2) then
       Values_Agree_Real=.true.
    else
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
            &            .or.(present(relTol).and.agreeRelatively)
    end if
    return
  end function Values_Agree_Real

  elemental logical function Values_Agree_Double(value1,value2,absTol,relTol)
    !!{
    Returns true if {\normalfont \ttfamily value1} and {\normalfont \ttfamily value2} agree to within {\normalfont \ttfamily absTol} in absolute terms, or {\normalfont \ttfamily relTol} in
    relative terms.
    !!}
    implicit none
    double precision, intent(in   )           :: value1         , value2
    double precision, intent(in   ), optional :: absTol         , relTol
    logical                                   :: agreeAbsolutely, agreeRelatively

    if (value1 == value2) then
       Values_Agree_Double=.true.
    else
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
    end if
    return
  end function Values_Agree_Double

  elemental logical function Values_Agree_Double_Complex(value1,value2,absTol,relTol)
    !!{
    Returns true if {\normalfont \ttfamily value1} and {\normalfont \ttfamily value2} agree to within {\normalfont \ttfamily absTol} in absolute terms, or {\normalfont \ttfamily relTol} in
    relative terms.
    !!}
    implicit none
    double complex, intent(in   )           :: value1         , value2
    double complex, intent(in   ), optional :: absTol         , relTol
    logical                                 :: agreeAbsolutely, agreeRelatively

    if (value1 == value2) then
       Values_Agree_Double_Complex=.true.
    else
       if (present(absTol)) then
          agreeAbsolutely= abs(real(value1-value2)) <=                                real(absTol) &
               &          .and.                                                                    &
               &           abs(imag(value1-value2)) <=                                imag(absTol)
       else
          agreeAbsolutely=.true.
       end if
       if (present(relTol)) then
          agreeRelatively= abs(real(value1-value2)) <= 0.5d0*abs(real(value1+value2))*real(relTol) &
               &          .and.                                                                    &
               &           abs(imag(value1-value2)) <= 0.5d0*abs(imag(value1+value2))*imag(relTol)
       else
          agreeRelatively=.true.
       end if
       Values_Agree_Double_Complex=    (present(absTol).and.agreeAbsolutely) &
            &                      .or.(present(relTol).and.agreeRelatively)
    end if
    return
  end function Values_Agree_Double_Complex

  logical function Values_Less_Than_Double(value1,value2,absTol,relTol)
    !!{
    Returns true if {\normalfont \ttfamily value1} is significantly less than {\normalfont \ttfamily value2}, with tolerance
    {\normalfont \ttfamily absTol} in absolute terms, or {\normalfont \ttfamily relTol} in relative terms.
    !!}
    implicit none
    double precision, intent(in   )           :: value1         , value2
    double precision, intent(in   ), optional :: absTol         , relTol
    logical                                   :: lessThanAbsolutely, lessThanRelatively

    if (value1 >= value2) then
       Values_Less_Than_Double=.false.
    else
       if (present(absTol)) then
          lessThanAbsolutely=(value2-value1 >= absTol)
       else
          lessThanAbsolutely=.true.
       end if
       if (present(relTol)) then
          lessThanRelatively=(value2-value1 >= 0.5d0*abs(value1+value2)*relTol)
       else
          lessThanRelatively=.true.
       end if
       Values_Less_Than_Double=    (present(absTol).and.lessThanAbsolutely) &
            &                  .or.(present(relTol).and.lessThanRelatively)
    end if
    return
  end function Values_Less_Than_Double

end module Numerical_Comparison
