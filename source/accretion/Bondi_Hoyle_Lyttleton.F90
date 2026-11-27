!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which implements calculations of Bondi-Hoyle-Lyttleton accretion (see :cite:author:`edgar_review_2004` :cite:year:`edgar_review_2004`).
!!}

module Bondi_Hoyle_Lyttleton_Accretion
  !!{RST
  Implements calculations of Bondi-Hoyle-Lyttleton accretion (see :cite:author:`edgar_review_2004` :cite:year:`edgar_review_2004`).
  !!}
  implicit none
  private
  public :: Bondi_Hoyle_Lyttleton_Accretion_Rate, Bondi_Hoyle_Lyttleton_Accretion_Radius

contains

  double precision function Bondi_Hoyle_Lyttleton_Accretion_Rate(mass,density,velocity,temperature,radius)
    !!{RST
    Computes the Bondi-Hoyle-Lyttleton accretion rate (in :math:`\mathrm{M}_\odot` Gyr\ :math:`^{-1}`; :cite:author:`edgar_review_2004` :cite:year:`edgar_review_2004`).
    !!}
    use :: Ideal_Gases_Thermodynamics      , only : Ideal_Gas_Sound_Speed
    use :: Numerical_Constants_Astronomical, only : gigaYear             , megaParsec, gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision, intent(in   )           :: density   , mass, temperature, velocity
    double precision, intent(in   ), optional :: radius
    double precision                          :: soundSpeed

    ! Compute the sound speed.
    soundSpeed=Ideal_Gas_Sound_Speed(temperature)

    ! Compute the accretion rate.
    if (present(radius)) then
       Bondi_Hoyle_Lyttleton_Accretion_Rate=(kilo*gigaYear/megaParsec)*4.0d0*Pi*radius**2*density*sqrt(soundSpeed**2+velocity**2)
    else
       ! Assume that the accretion radius is equal to the Bondi-Hoyle-Lyttleton radius.
       Bondi_Hoyle_Lyttleton_Accretion_Rate=(kilo*gigaYear/megaParsec)*4.0d0*Pi*(gravitationalConstant_internal*mass)**2*density &
            &/(soundSpeed**2+velocity**2)**1.5d0
    end if
    return
  end function Bondi_Hoyle_Lyttleton_Accretion_Rate

  double precision function Bondi_Hoyle_Lyttleton_Accretion_Radius(mass,temperature,velocity) result(radiusAccretion)
    !!{RST
    Computes the Bondi-Hoyle-Lyttleton accretion radius, :math:`\mathrm{G} M/(c_\mathrm{s}^2+v^2)` (in Mpc;
    :cite:author:`edgar_review_2004` :cite:year:`edgar_review_2004`), where :math:`v` is the (optional) velocity of the gas
    relative to the accretor.
    !!}
    use :: Ideal_Gases_Thermodynamics      , only : Ideal_Gas_Sound_Speed
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    double precision, intent(in   )           :: mass           , temperature
    double precision, intent(in   ), optional :: velocity
    double precision                          :: velocitySquared

    velocitySquared=0.0d0
    if (temperature > 0.0d0) velocitySquared=velocitySquared+Ideal_Gas_Sound_Speed(temperature)**2
    if (present(velocity)  ) velocitySquared=velocitySquared+                      velocity    **2
    if (velocitySquared > 0.0d0) then
       radiusAccretion=gravitationalConstant_internal*mass/velocitySquared
    else
       radiusAccretion=huge(0.0d0)
    end if
    return
  end function Bondi_Hoyle_Lyttleton_Accretion_Radius

end module Bondi_Hoyle_Lyttleton_Accretion
