!% Contains a module which implements thermodynamic properties of ideal gases.

module Ideal_Gases_Thermodynamics
  !% Implements thermodynamic properties of ideal gases.
  private
  public :: Ideal_Gas_Sound_Speed, Ideal_Gas_Jeans_Length

contains

  double precision function Ideal_Gas_Jeans_Length(temperature,density)
    !% Return the Jeans length (in Mpc) for gas of given temperature and density).
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in) :: temperature,density

    Ideal_gas_Jeans_Length=Ideal_Gas_Sound_Speed(temperature)/dsqrt(gravitationalConstantGalacticus*density)
    return
  end function Ideal_Gas_Jeans_Length

  double precision function Ideal_Gas_Sound_Speed(temperature,meanAtomicMass)
    !% Return the sound speed (in km/s) for an ideal gas of given {\tt temperature} and (optionally) {\tt meanAtomicMass}.
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Atomic
    implicit none
    double precision, intent(in)           :: temperature
    double precision, intent(in), optional :: meanAtomicMass
    double precision                       :: meanAtomicMassActual
    
    ! Determine what mean atomic mass to use.
    if (present(meanAtomicMass)) then
       meanAtomicMassActual=meanAtomicMass
    else
       meanAtomicMassActual=meanAtomicMassPrimordial
    end if
    
    ! Compute the sound speed.
    Ideal_Gas_Sound_Speed=dsqrt(5.0d0*boltzmannsConstant*temperature/3.0d0/meanAtomicMassActual/atomicMassUnit)/kilo
    return
  end function Ideal_Gas_Sound_Speed
  
end module Ideal_Gases_Thermodynamics
