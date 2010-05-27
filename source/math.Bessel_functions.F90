!% Contains a module which implements calculations of Bessel functions.

module Bessel_Functions
  !% Implements calculations of Bessel functions.
  use FGSL
  private
  public :: Bessel_Function_K0, Bessel_Function_K1, Bessel_Function_I0, Bessel_Function_I1

contains

  double precision function Bessel_Function_K0(argument)
    !% Computes the $K_0$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    Bessel_Function_K0=FGSL_SF_Bessel_Kc0(argument)
    return
  end function Bessel_Function_K0
  
  double precision function Bessel_Function_K1(argument)
    !% Computes the $K_1$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    Bessel_Function_K1=FGSL_SF_Bessel_Kc1(argument)
    return
  end function Bessel_Function_K1
  
  double precision function Bessel_Function_I0(argument)
    !% Computes the $I_0$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    Bessel_Function_I0=FGSL_SF_Bessel_Ic0(argument)
    return
  end function Bessel_Function_I0
  
  double precision function Bessel_Function_I1(argument)
    !% Computes the $I_1$ Bessel function.
    implicit none
    double precision, intent(in) :: argument

    Bessel_Function_I1=FGSL_SF_Bessel_Ic1(argument)
    return
  end function Bessel_Function_I1
  
end module Bessel_Functions
