!% Contains a module of useful atomic constants.

module Numerical_Constants_Atomic
  !% Contains various useful atomic constants.
  use FGSL
  public

  ! Atomic mass unit (in kg).
  double precision, parameter :: atomicMassUnit    =FGSL_CONST_MKSA_UNIFIED_ATOMIC_MASS

  ! Atomic masses.
  double precision, parameter :: atomicMassHydrogen=1.007825d0
  double precision, parameter :: atomicMassHelium  =4.002602d0
  
  ! Mass of hydrogen atom (in kg).
  double precision, parameter :: massHydrogenAtom=atomicMassHydrogen*atomicMassUnit

end module Numerical_Constants_Atomic
