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
  Implements a dust extinction curve read from a tabulation.
  !!}

  use :: Tables, only : table1DGeneric

  !![
  <dustExtinctionCurve name="dustExtinctionCurveTabulated" abstract="yes" docformat="rst">
   <description>
   An abstract dust extinction curve given by a tabulation of :math:`k(\lambda)` against inverse wavelength,
   :math:`x=1/\lambda` in inverse microns---the variable in which such curves are conventionally published.

   Unlike the analytic curves, a tabulated curve is normalized by *its own* interpolated value at the :math:`V`-band
   effective wavelength, so it returns exactly unity there by construction.
   </description>
  </dustExtinctionCurve>
  !!]
  ! Note: declared as a concrete type despite `abstract="yes"` in the directive above. The directive keeps the class
  ! out of the set which can be built from a parameter file, while leaving the Fortran type concrete so that it may
  ! carry a finalizer -- a final subroutine's `self` cannot be declared `type(...)` for an abstract type. This mirrors
  ! `stellarSpectraDustAttenuationTabulated`, the class replaced here.
  type, extends(dustExtinctionCurveClass) :: dustExtinctionCurveTabulated
     !!{RST
     A dust extinction curve given by a tabulation.
     !!}
     private
     type(table1DGeneric) :: attenuationTable
   contains
     final     ::                        tabulatedDestructor
     procedure :: attenuationRelative => tabulatedAttenuationRelative
     procedure :: wavelengthRange     => tabulatedWavelengthRange
  end type dustExtinctionCurveTabulated

contains

  subroutine tabulatedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustExtinctionCurveTabulated` dust extinction curve class.
    !!}
    implicit none
    type(dustExtinctionCurveTabulated), intent(inout) :: self

    call self%attenuationTable%destroy()
    return
  end subroutine tabulatedDestructor

  double precision function tabulatedAttenuationRelative(self,wavelength) result(attenuationRelative)
    !!{RST
    Return the relative dust opacity for a tabulated extinction curve.
    !!}
    use :: Numerical_Constants_Units, only : micronsToAngstroms
    implicit none
    class           (dustExtinctionCurveTabulated), intent(inout) :: self
    double precision                              , intent(in   ) :: wavelength
    double precision                              , parameter     :: xV        =1.0d0/(wavelengthVBand/micronsToAngstroms)
    double precision                                              :: x

    x                  =+1.0d0                                 &
         &              /(                                     &
         &                +wavelength                          &
         &                /micronsToAngstroms                  &
         &               )
    attenuationRelative=+self%attenuationTable%interpolate(x ) &
         &              /self%attenuationTable%interpolate(xV)
    return
  end function tabulatedAttenuationRelative

  subroutine tabulatedWavelengthRange(self,wavelengthMinimum,wavelengthMaximum)
    !!{RST
    Return the range of wavelengths spanned by the tabulation. Note that these curves are constructed with
    extrapolating tables, so they return a value outside this range rather than zero---but that value is an
    extrapolation of the published data, not a fit to it.
    !!}
    use :: Numerical_Constants_Units, only : micronsToAngstroms
    implicit none
    class           (dustExtinctionCurveTabulated), intent(inout) :: self
    double precision                              , intent(  out) :: wavelengthMinimum, wavelengthMaximum

    wavelengthMinimum=micronsToAngstroms/self%attenuationTable%x(-1)
    wavelengthMaximum=micronsToAngstroms/self%attenuationTable%x(+1)
    return
  end subroutine tabulatedWavelengthRange
