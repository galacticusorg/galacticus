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

!% Contains a module which implements calculations of gamma functions.

module Gamma_Functions
  !% Implements calculations of gamma functions.
  implicit none
  private
  public :: Gamma_Function_Incomplete,Gamma_Function_Incomplete_Complementary,Gamma_Function_Logarithmic,Gamma_Function&
       &,Inverse_Gamma_Function_Incomplete,Inverse_Gamma_Function_Incomplete_Complementary

contains

  double precision function Gamma_Function_Incomplete(exponent,argument)
    !% Computes the incomplete gamma function.
    use FGSL
    implicit none
    double precision, intent(in   ) :: argument, exponent

    Gamma_Function_Incomplete=FGSL_SF_Gamma_Inc_Q(exponent,argument)
    return
  end function Gamma_Function_Incomplete

  double precision function Gamma_Function_Incomplete_Complementary(exponent,argument)
    !% Computes the complementary incomplete gamma function.
    use FGSL
    implicit none
    double precision, intent(in   ) :: argument, exponent

    Gamma_Function_Incomplete_Complementary=FGSL_SF_Gamma_Inc_P(exponent,argument)
    return
  end function Gamma_Function_Incomplete_Complementary

  double precision function Gamma_Function(argument)
    !% Computes the gamma function.
    implicit none
    double precision, intent(in   ) :: argument

    Gamma_Function=exp(Gamma_Function_Logarithmic(argument))
    return
  end function Gamma_Function

  double precision function Gamma_Function_Logarithmic(argument)
    !% Computes the logarithm of the gamma function.
    use FGSL
    implicit none
    double precision, intent(in   ) :: argument

    Gamma_Function_Logarithmic=FGSL_SF_lnGamma(argument)
    return
  end function Gamma_Function_Logarithmic

  double precision function Inverse_Gamma_Function_Incomplete_Complementary(a,P)
    !% Returns the inverse of the incomplete function. That is, it returns $x$ given $P(a,x)$.
    use ISO_Varying_String
    use Incomplete_Gamma
    use Galacticus_Error
    implicit none
    double precision                , intent(in   ) :: P         , a
    integer                                         :: errorState
    double precision                                :: Q
    type            (varying_string)                :: message

    Q=1.0d0-P
    call GamInv(a,Inverse_Gamma_Function_Incomplete_Complementary,-1.0d0,P,Q,errorState)
    if (errorState < 0) then
       select case (errorState)
       case (-2)
          message='input error: a <= 0'
       case (-3)
          message='no solution obtained: A/a is too large'
       case (-4)
          message='input error: P or Q < 0 or P+Q != 1'
       case (-6)
          message='20 iterations were performed - this should not happen'
       case (-7)
          message='iteration failed'
       case (-8)
          message='accuracy lost'
       end select
       call Galacticus_Error_Report('Inverse_Gamma_Function_Incomplete_Complementary',message)
    end if
    return
  end function Inverse_Gamma_Function_Incomplete_Complementary

  double precision function Inverse_Gamma_Function_Incomplete(a,Q)
    !% Returns the inverse of the incomplete function. That is, it returns $x$ given $Q(a,x)$.
    implicit none
    double precision, intent(in   ) :: Q, a
    double precision                :: P

    P=1.0d0-Q
    Inverse_Gamma_Function_Incomplete=Inverse_Gamma_Function_Incomplete_Complementary(a,P)
    return
  end function Inverse_Gamma_Function_Incomplete

end module Gamma_Functions
