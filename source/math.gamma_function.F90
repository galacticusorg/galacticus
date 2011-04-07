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


!% Contains a module which implements calculations of gamma functions.

module Gamma_Functions
  !% Implements calculations of gamma functions.
  private
  public :: Gamma_Function_Incomplete,Gamma_Function_Incomplete_Complementary,Gamma_Function_Logarithmic,Gamma_Function&
       &,Inverse_Gamma_Function_Incomplete,Inverse_Gamma_Function_Incomplete_Complementary

contains

  double precision function Gamma_Function_Incomplete(exponent,argument)
    !% Computes the incomplete gamma function.
    use FGSL
    implicit none
    double precision, intent(in) :: exponent,argument

    Gamma_Function_Incomplete=FGSL_SF_Gamma_Inc_Q(exponent,argument)
    return
  end function Gamma_Function_Incomplete
  
  double precision function Gamma_Function_Incomplete_Complementary(exponent,argument)
    !% Computes the complementary incomplete gamma function.
    use FGSL
    implicit none
    double precision, intent(in) :: exponent,argument

    Gamma_Function_Incomplete_Complementary=FGSL_SF_Gamma_Inc_P(exponent,argument)
    return
  end function Gamma_Function_Incomplete_Complementary
  
  double precision function Gamma_Function(argument)
    !% Computes the gamma function.
    use FGSL
    implicit none
    double precision, intent(in) :: argument

    Gamma_Function=dexp(Gamma_Function_Logarithmic(argument))
    return
  end function Gamma_Function
  
  double precision function Gamma_Function_Logarithmic(argument)
    !% Computes the logarithm of the gamma function.
    use FGSL
    implicit none
    double precision, intent(in) :: argument

    Gamma_Function_Logarithmic=FGSL_SF_lnGamma(argument)
    return
  end function Gamma_Function_Logarithmic

  double precision function Inverse_Gamma_Function_Incomplete_Complementary(a,P)
    !% Returns the inverse of the incomplete function. That is, it returns $x$ given $P(a,x)$.
    use ISO_Varying_String
    use Incomplete_Gamma
    use Galacticus_Error
    implicit none
    double precision,     intent(in) :: a,P
    integer                          :: errorState
    double precision                 :: Q
    type(varying_string)             :: message

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
    double precision,     intent(in) :: a,Q
    double precision                 :: P

    P=1.0d0-Q
    Inverse_Gamma_Function_Incomplete=Inverse_Gamma_Function_Incomplete_Complementary(a,P)
    return
  end function Inverse_Gamma_Function_Incomplete
  
end module Gamma_Functions
