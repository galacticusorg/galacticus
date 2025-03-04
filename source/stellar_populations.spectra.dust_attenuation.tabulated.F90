!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements calculations of attenuation of stellar spectra using a tabulation.
  !!}

  use :: Tables, only : table1DGeneric

  !![
  <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationTabulated" abstract="yes">
   <description>Returns the dust attenuation of stellar spectra from a tabulated relation.</description>
  </stellarSpectraDustAttenuation>
  !!]

  type, extends(stellarSpectraDustAttenuationClass) :: stellarSpectraDustAttenuationTabulated
     !!{
     A class implementing calculations of attenuation of stellar spectra using a tabulated relation.
     !!}
     private
     type(table1DGeneric) :: attenuationTable
   contains
     final     ::                tabulatedDestructor
     procedure :: attenuation => tabulatedAttenuation
  end type stellarSpectraDustAttenuationTabulated

contains

  subroutine tabulatedDestructor(self)
    !!{
    Destructor for the ``tabulated'' stellar spectra dust attenuation class.
    !!}
    implicit none
    type(stellarSpectraDustAttenuationTabulated), intent(inout) :: self

    call self%attenuationTable%destroy()
    return
  end subroutine tabulatedDestructor

  double precision function tabulatedAttenuation(self,wavelength,age,vBandAttenuation)
    !!{
    Return attenuation of stellar spectra according to the model of \cite{gordon_quantitative_2003}.
    !!}
    use :: Numerical_Constants_Units, only : micronsToAngstroms
    implicit none
    class           (stellarSpectraDustAttenuationTabulated), intent(inout) :: self
    double precision                                        , intent(in   ) :: wavelength                                                    , age, &
         &                                                                     vBandAttenuation
    double precision                                        , parameter     :: xV              =1.0d0/(5504.61227375652d0/micronsToAngstroms)
    double precision                                                        :: x
    !$GLC attributes unused :: age

    x=1.0d0/(wavelength/micronsToAngstroms)
    tabulatedAttenuation=vBandAttenuation*self%attenuationTable%interpolate(x)/self%attenuationTable%interpolate(xV)
    return
  end function tabulatedAttenuation
