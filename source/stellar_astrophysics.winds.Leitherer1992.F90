!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements a calculation of winds from stellar populations using the fitting formulae of \cite{leitherer_deposition_1992}.

module Stellar_Astrophysics_Winds_Leitherer1992
  !% Implements a calculation of winds from stellar populations using the fitting formulae of \cite{leitherer_deposition_1992}.
  use Numerical_Constants_Astronomical
  private
  public :: Stellar_Winds_Leitherer1992_Initialize

  ! Minimum metallicity to which we trust Leitherer et al.'s metallicity scaling.
  double precision, parameter :: metallicityMinimum=1.0d-4*metallicitySolar

contains

  !# <stellarWindsMethod>
  !#  <unitName>Stellar_Winds_Leitherer1992_Initialize</unitName>
  !# </stellarWindsMethod>
  subroutine Stellar_Winds_Leitherer1992_Initialize(stellarWindsMethod,Stellar_Winds_Mass_Loss_Rate_Get,Stellar_Winds_Terminal_Velocity_Get)
    !% Initialize the ``Leitherer1992'' stellar winds module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string), intent(in) :: stellarWindsMethod
    procedure(),          pointer    :: Stellar_Winds_Mass_Loss_Rate_Get,Stellar_Winds_Terminal_Velocity_Get

    if (stellarWindsMethod == 'Leitherer1992') then
       ! Set procedure pointers.
       Stellar_Winds_Mass_Loss_Rate_Get    => Stellar_Winds_Mass_Loss_Rate_Leitherer1992
       Stellar_Winds_Terminal_Velocity_Get => Stellar_Winds_Terminal_Velocity_Leitherer1992
    end if
    return
  end subroutine Stellar_Winds_Leitherer1992_Initialize

  double precision function Stellar_Winds_Mass_Loss_Rate_Leitherer1992(initialMass,age,metallicity)
    !% Compute the mass loss rate (in $M_\odot$/Gyr) from a star of given {\tt initialMass}, {\tt age} and {\tt metallicity} using
    !% the fitting formula of \cite{leitherer_deposition_1992}.
    use Stellar_Astrophysics_Tracks
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity
    double precision             :: stellarLuminosity,stellarEffectiveTemperature

    ! Get luminosity and effective temperature of the star.
    stellarLuminosity          =Stellar_Luminosity           (initialMass,metallicity,age)
    stellarEffectiveTemperature=Stellar_Effective_Temperature(initialMass,metallicity,age)

    ! Compute mass loss rate using fitting formula. (Initial constant 9 converts Leitherer's mass loss rate from per year to per Gyr.)
    if (stellarLuminosity > 0.0d0 .and. stellarEffectiveTemperature > 0.0d0) then
       Stellar_Winds_Mass_Loss_Rate_Leitherer1992=10.0d0**(9.0d0-24.06d0+2.45d0*dlog10(stellarLuminosity)-1.10d0&
            &*dlog10(initialMass) +1.31d0*dlog10(stellarEffectiveTemperature)+0.80d0*dlog10(max(metallicity,metallicityMinimum)&
            &/metallicitySolar))
    else
       Stellar_Winds_Mass_Loss_Rate_Leitherer1992=0.0d0
    end if
    return
  end function Stellar_Winds_Mass_Loss_Rate_Leitherer1992

  double precision function Stellar_Winds_Terminal_Velocity_Leitherer1992(initialMass,age,metallicity)
    !% Compute the terminal velocity (in km/s) from a star of given {\tt initialMass}, {\tt age} and {\tt metallicity} using
    !% the fitting formula of \cite{leitherer_deposition_1992}.
    use Stellar_Astrophysics_Tracks
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity
    double precision             :: stellarLuminosity,stellarEffectiveTemperature

    ! Get luminosity and effective temperature of the star.
    stellarLuminosity          =Stellar_Luminosity           (initialMass,metallicity,age)
    stellarEffectiveTemperature=Stellar_Effective_Temperature(initialMass,metallicity,age)

    ! Compute mass loss rate using fitting formula.
    if (stellarLuminosity > 0.0d0 .and. stellarEffectiveTemperature > 0.0d0) then
       Stellar_Winds_Terminal_Velocity_Leitherer1992=10.0d0**(1.23d0-0.30d0*dlog10(stellarLuminosity)+0.55d0 *dlog10(initialMass)&
            &+0.64d0*dlog10(stellarEffectiveTemperature)+0.12d0*dlog10(max(metallicity,metallicityMinimum)/metallicitySolar))
    else
       Stellar_Winds_Terminal_Velocity_Leitherer1992=0.0d0
    end if
    return
  end function Stellar_Winds_Terminal_Velocity_Leitherer1992

end module Stellar_Astrophysics_Winds_Leitherer1992
