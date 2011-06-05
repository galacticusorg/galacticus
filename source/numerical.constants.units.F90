!% Contains a module of useful unit conversions.

module Numerical_Constants_Units
  !% Contains various useful unit conversions.
  use Numerical_Constants_Astronomical
  use Numerical_Constants_Prefixes
  public
  
  ! Ergs in Joules.
  double precision, parameter :: ergs=1.0d-7

  ! Conversion from Mpc/(km/s) to Gyr.
  double precision, parameter :: Mpc_per_km_per_s_To_Gyr=megaParsec/kilo/gigaYear

end module Numerical_Constants_Units
