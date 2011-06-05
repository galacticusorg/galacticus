!% Contains a module which implements hypergeometric functions.

module Hypergeometric_Functions
  !% Implements hypergeometric functions.
  use FGSL
  private
  public :: Hypergeometric_2F1
  
contains
  
  double precision function Hypergeometric_2F1(a,b,x)
    !% Evaluate the $_2F_1(a_1,a_2;b_1;x)$ hypergeometric function.
    implicit none
    double precision, intent(in) :: a(2),b(1),x

    Hypergeometric_2F1=FGSL_SF_Hyperg_2F1(a(1),a(2),b(1),x)
    return
  end function Hypergeometric_2F1
  
end module Hypergeometric_Functions
