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

  !!{RST
  Contains a module which provides utility variables for the :cite:t:`schive_understanding_2014` soliton density profile.
  !!}

  module Mass_Distribution_Soliton_Schive2014
    !!{RST
    Provides utility variables for the :cite:t:`schive_understanding_2014` soliton density profile.
    !!}
    implicit none
    private

    ! Coefficient of the dimensionless radius in the soliton profile.
    double precision, parameter, public  :: coefficientCore              =0.091d0 ! Schive et al. (2014; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S; equation 3).

    !+ The closed-form expression for coefficientMassCore below was derived with assistance from Claude (Anthropic).

    ! Coefficient relating the solitonic core mass to the central density and core radius, such that
    !
    !  M(<r_c) = coefficientMassCore * π * ρ_c * r_c³,
    !
    ! where M(<r_c) is the mass enclosed within the core radius of the Schive et al. (2014; equation 3;
    ! https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S) soliton density profile,
    !
    !  ρ(r) = ρ_c / [ 1 + coefficientCore (r/r_c)² ]⁸.
    !
    ! This is simply the r=r_c case of the analytic solution evaluated in solitonMassEnclosedBySphere(), i.e.
    !
    !  coefficientMassCore = [ 3465 tan-¹(√c) + √c P(c) / (1+c)⁷ ] / ( 53760 c^(3/2) ),
    !
    ! with c=coefficientCore and P(c) the polynomial below. For c=0.091 this evaluates to 0.888026837983102.
    double precision, parameter, private :: coefficientMassCorePolynomial=+  3465.0d0*coefficientCore**6 &
         &                                                                + 23100.0d0*coefficientCore**5 &
         &                                                                + 65373.0d0*coefficientCore**4 &
         &                                                                +101376.0d0*coefficientCore**3 &
         &                                                                + 92323.0d0*coefficientCore**2 &
         &                                                                + 48580.0d0*coefficientCore    &
         &                                                                -  3465.0d0
    double precision, parameter, public  :: coefficientMassCore          =+(                                      &
         &                                                                  +3465.0d0*atan(sqrt(coefficientCore)) &
         &                                                                  +sqrt(coefficientCore)                &
         &                                                                  *coefficientMassCorePolynomial        &
         &                                                                  /(1.0d0+coefficientCore)**7           &
         &                                                                 )                                      &
         &                                                                /(                                      &
         &                                                                  +53760.0d0                            &
         &                                                                  *coefficientCore**1.5d0               &
         &                                                                 )

  end module Mass_Distribution_Soliton_Schive2014
