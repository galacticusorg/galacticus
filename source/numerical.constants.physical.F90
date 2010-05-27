!% Contains a module of useful physical constants.

module Numerical_Constants_Physical
  !% Contains various useful physical constants.
  use FGSL
  use Numerical_Constants_Prefixes
  public

  ! Speed of light (m/s).
  double precision, parameter :: speedLight=FGSL_CONST_MKSA_SPEED_OF_LIGHT
  
  ! Newton's gravitational constant (in Galacticus' M_Solar, Mpc, km/s unit system).
  double precision, parameter :: gravitationalConstantGalacticus=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT&
       &*FGSL_CONST_MKSA_SOLAR_MASS/(kilo**2)/FGSL_CONST_MKSA_PARSEC/mega

  ! Newton's gravitational constant (in SI units).
  double precision, parameter :: gravitationalConstant=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT

  ! Radiation constant (in units of J/m^3/K^4).
  double precision, parameter :: radiationConstant=4.0d0*FGSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT/FGSL_CONST_MKSA_SPEED_OF_LIGHT

  ! Boltzmann's constant (in units of J/K).
  double precision, parameter :: boltzmannsConstant=FGSL_CONST_MKSA_BOLTZMANN

  ! Thomson cross section (in units of m^2).
  double precision, parameter :: thomsonCrossSection=FGSL_CONST_MKSA_THOMSON_CROSS_SECTION

end module Numerical_Constants_Physical
