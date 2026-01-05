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

!!{
Contains a module which implements polylogarithm functions.
!!}

! Add dependency on GSL library.
!; gsl

module Polylogarithms
  !!{
  Implements polylogarithm functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double
  implicit none
  private
  public :: Polylogarithm_2, Polylogarithm_3

  interface
     function gsl_sf_fermi_dirac_1(x) bind(c,name='gsl_sf_fermi_dirac_1')
       !!{
       Template for the GSL Fermi-Dirac integral of index 1 C function.
       !!}
       import
       real(c_double)        :: gsl_sf_fermi_dirac_1
       real(c_double), value :: x
     end function gsl_sf_fermi_dirac_1
  end interface

  interface
     function gsl_sf_fermi_dirac_2(x) bind(c,name='gsl_sf_fermi_dirac_2')
       !!{
       Template for the GSL Fermi-Dirac integral of index 2 C function.
       !!}
       import
       real(c_double)        :: gsl_sf_fermi_dirac_2
       real(c_double), value :: x
     end function gsl_sf_fermi_dirac_2
  end interface

contains

  double precision function Polylogarithm_2(x)
    !!{
    Evaluate the polylogarithm function of order 2.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    ! Use the relation between polylogarithms and the Fermi-Dirac integral
    ! (e.g. https://en.wikipedia.org/wiki/Polylogarithm#Relationship_to_other_functions).
    Polylogarithm_2=-GSL_SF_Fermi_Dirac_1(log(-x))
    return
  end function Polylogarithm_2

  double precision function Polylogarithm_3(x)
    !!{
    Evaluate the polylogarithm function of order 3.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    ! Use the relation between polylogarithms and the Fermi-Dirac integral
    ! (e.g. https://en.wikipedia.org/wiki/Polylogarithm#Relationship_to_other_functions).
    Polylogarithm_3=-GSL_SF_Fermi_Dirac_2(log(-x))
    return
  end function Polylogarithm_3

end module Polylogarithms
