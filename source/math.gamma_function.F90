!% Contains a module which implements calculations of gamma functions.

module Gamma_Functions
  !% Implements calculations of gamma functions.
  private
  public :: Gamma_Function_Incomplete_Complementary

contains

  double precision function Gamma_Function_Incomplete_Complementary(exponent,argument)
    !% Computes the complementary incomplete gamma function.
    use FGSL
    implicit none
    double precision, intent(in) :: exponent,argument

    Gamma_Function_Incomplete_Complementary=FGSL_SF_Gamma_Inc_P(exponent,argument)
    return
  end function Gamma_Function_Incomplete_Complementary
  
end module Gamma_Functions
