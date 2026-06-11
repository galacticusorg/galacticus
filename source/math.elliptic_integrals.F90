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
Contains a module which implements exponential integrals.
!!}

! Add dependency on GSL library.
!; gsl

module Elliptic_Integrals
  !!{RST
  Implements exponential integrals.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double, c_int
  implicit none
  private
  public :: Elliptic_Integral_K           , Elliptic_Integral_E, Elliptic_Integral_Pi, Incomplete_Elliptic_Integral_E, &
       &    Incomplete_Elliptic_Integral_F

  interface
     function gsl_sf_ellint_Kcomp(k,mode) bind(c,name='gsl_sf_ellint_Kcomp')
       !!{RST
       Template for the GSL :math:`K(k)` elliptic integral function
       !!}
       import c_double, c_int
       real   (c_double)        :: gsl_sf_ellint_Kcomp
       real   (c_double), value :: k
       integer(c_int   ), value :: mode
     end function gsl_sf_ellint_Kcomp
     function gsl_sf_ellint_Ecomp(k,mode) bind(c,name='gsl_sf_ellint_Ecomp')
       !!{RST
       Template for the GSL :math:`E(k)` elliptic integral function
       !!}
       import c_double, c_int
       real   (c_double)        :: gsl_sf_ellint_Ecomp
       real   (c_double), value :: k
       integer(c_int   ), value :: mode
     end function gsl_sf_ellint_Ecomp
     function gsl_sf_ellint_Pcomp(k,n,mode) bind(c,name='gsl_sf_ellint_Pcomp')
       !!{RST
       Template for the GSL :math:`\Pi(k,n)` elliptic integral function
       !!}
       import c_double, c_int
       real   (c_double)        :: gsl_sf_ellint_Pcomp
       real   (c_double), value :: k                  , n
       integer(c_int   ), value :: mode
     end function gsl_sf_ellint_Pcomp
     function gsl_sf_ellint_E(phi,k,mode) bind(c,name='gsl_sf_ellint_E')
       !!{RST
       Template for the GSL :math:`E(\phi,k)` incomplete elliptic integral function
       !!}
       import c_double, c_int
       real   (c_double)        :: gsl_sf_ellint_E
       real   (c_double), value :: phi            , k
       integer(c_int   ), value :: mode
     end function gsl_sf_ellint_E
     function gsl_sf_ellint_F(phi,k,mode) bind(c,name='gsl_sf_ellint_F')
       !!{RST
       Template for the GSL :math:`F(\phi,k)` incomplete elliptic integral function
       !!}
       import c_double, c_int
       real   (c_double)        :: gsl_sf_ellint_F
       real   (c_double), value :: phi            , k
       integer(c_int   ), value :: mode
     end function gsl_sf_ellint_F
  end interface
  
contains

  double precision function Elliptic_Integral_K(m)
    !!{RST
    Evaluate the :math:`K(m)` elliptical integral. Note that we use the :cite:t:`abramowitz_handbook_1970` notation here, while GSL uses the :cite:t:`carlson_computing_1979` notation. The conversion from :math:`m` to :math:`k=\sqrt{m}` is performed in the call to the GSL function.
    !!}
    use :: Interface_GSL, only : GSL_Prec_Double
    implicit none
    double precision, intent(in   ) :: m

    Elliptic_Integral_K=gsl_sf_ellint_Kcomp(sqrt(m),GSL_Prec_Double)
    return
  end function Elliptic_Integral_K

  double precision function Elliptic_Integral_E(m)
    !!{RST
    Evaluate the :math:`E(m)` elliptical integral. Note that we use the :cite:t:`abramowitz_handbook_1970` notation here, while GSL uses the :cite:t:`carlson_computing_1979` notation. The conversion from :math:`m` to :math:`k=\sqrt{m}` is performed in the call to the GSL function.
    !!}
    use :: Interface_GSL, only : GSL_Prec_Double
    implicit none
    double precision, intent(in   ) :: m

    Elliptic_Integral_E=gsl_sf_ellint_Ecomp(sqrt(m),GSL_Prec_Double)
    return
  end function Elliptic_Integral_E

  double precision function Elliptic_Integral_Pi(m,n)
    !!{RST
    Evaluate the :math:`\Pi(m,n)` elliptical integral. Note that we use the :cite:t:`abramowitz_handbook_1970` notation here, while GSL uses the :cite:t:`carlson_computing_1979` notation. The conversion from :math:`m` to :math:`k=\sqrt{m}` and :math:`n \rightarrow -n` is performed in the call to the GSL function.
    !!}
    use :: Interface_GSL, only : GSL_Prec_Double
    implicit none
    double precision, intent(in   ) :: m, n

    Elliptic_Integral_Pi=gsl_sf_ellint_Pcomp(sqrt(m),-n,GSL_Prec_Double)
    return
  end function Elliptic_Integral_Pi

  double precision function Incomplete_Elliptic_Integral_E(phi,m)
    !!{RST
    Evaluate the :math:`E(\phi,m)` elliptical integral. Note that we use the :cite:t:`abramowitz_handbook_1970` notation here, while GSL uses the :cite:t:`carlson_computing_1979` notation. The conversion from :math:`m` to :math:`k=\sqrt{m}` is performed in the call to the GSL function.
    !!}
    use :: Interface_GSL, only : GSL_Prec_Double
    implicit none
    double precision, intent(in   ) :: phi, m

    Incomplete_Elliptic_Integral_E=gsl_sf_ellint_E(phi,sqrt(m),GSL_Prec_Double)
    return
  end function Incomplete_Elliptic_Integral_E

  double precision function Incomplete_Elliptic_Integral_F(phi,m)
    !!{RST
    Evaluate the :math:`F(\phi,m)` elliptical integral. Note that we use the :cite:t:`abramowitz_handbook_1970` notation here, while GSL uses the :cite:t:`carlson_computing_1979` notation. The conversion from :math:`m` to :math:`k=\sqrt{m}` is performed in the call to the GSL function.
    !!}
    use :: Interface_GSL, only : GSL_Prec_Double
    implicit none
    double precision, intent(in   ) :: phi, m

    Incomplete_Elliptic_Integral_F=gsl_sf_ellint_F(phi,sqrt(m),GSL_Prec_Double)
    return
  end function Incomplete_Elliptic_Integral_F

end module Elliptic_Integrals
