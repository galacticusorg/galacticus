!% Contains a module which implements calculations of Bondi-Hoyle-Lyttleton accretion (see \citealt{edgar_review_2004}).

module Bondi_Hoyle_Lyttleton_Accretion
  !% Implements calculations of Bondi-Hoyle-Lyttleton accretion (see \citealt{edgar_review_2004}).
  private
  public :: Bondi_Hoyle_Lyttleton_Accretion_Rate, Bondi_Hoyle_Lyttleton_Accretion_Radius

  ! Flag indicating if module is initialized.
  logical :: bondiHoyleAccretionInitialized=.false.

contains

  double precision function Bondi_Hoyle_Lyttleton_Accretion_Rate(mass,density,velocity,temperature)
    !% Computes the Bondi-Hoyle-Lyttleton accretion rate (in $M_\odot$ Gyr$^{-1}$; \citealt{edgar_review_2004}).
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Ideal_Gases_Thermodynamics
    implicit none
    double precision, intent(in) :: mass,density,velocity,temperature
    double precision             :: soundSpeed

    ! Compute the sound speed.
    soundSpeed=Ideal_Gas_Sound_Speed(temperature)

    ! Compute the accretion rate.
    Bondi_Hoyle_Lyttleton_Accretion_Rate=(kilo*gigaYear/megaParsec)*4.0d0*Pi*(gravitationalConstantGalacticus*mass)**2*density &
         &/(soundSpeed**2+velocity**2)**1.5d0

    return
  end function Bondi_Hoyle_Lyttleton_Accretion_Rate

  double precision function Bondi_Hoyle_Lyttleton_Accretion_Radius(mass,temperature)
    !% Computes the Bondi-Hoyle-Lyttleton accretion radius (in Mpc; \citealt{edgar_review_2004}).
    use Numerical_Constants_Physical
    use Ideal_Gases_Thermodynamics
    implicit none
    double precision, intent(in) :: mass,temperature
    double precision             :: soundSpeed

    ! Compute the sound speed.
    soundSpeed=Ideal_Gas_Sound_Speed(temperature)

    ! Compute the accretion radius.
    Bondi_Hoyle_Lyttleton_Accretion_Radius=gravitationalConstantGalacticus*mass/soundSpeed**2

    return
  end function Bondi_Hoyle_Lyttleton_Accretion_Radius

end module Bondi_Hoyle_Lyttleton_Accretion
