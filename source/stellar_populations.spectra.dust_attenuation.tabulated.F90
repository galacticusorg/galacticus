!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Implements calculations of attenuation of stellar spectra using a tabulation.
  
  !# <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationTabulated">
  !#  <description>Returns the dust attenuation of stellar spectra from a tabulated relation.</description>
  !#  <abstract>yes</abstract>
  !# </stellarSpectraDustAttenuation>

  use Tables

  type, extends(stellarSpectraDustAttenuationClass) :: stellarSpectraDustAttenuationTabulated
     !% A class implementing calculations of attenuation of stellar spectra using a tabulated relation.
     private
     type(table1DGeneric) :: attenuationTable
   contains
     final     ::                tabulatedDestructor
     procedure :: attenuation => tabulatedAttenuation
  end type stellarSpectraDustAttenuationTabulated

contains

  subroutine tabulatedDestructor(self)
    !% Destructor for the ``tabulated'' stellar spectra dust attenuation class.
    use Input_Parameters
    implicit none
    type(stellarSpectraDustAttenuationTabulated), intent(inout) :: self

    call self%attenuationTable%destroy()
    return
  end subroutine tabulatedDestructor

  double precision function tabulatedAttenuation(self,wavelength,age,vBandAttenuation)
    !% Return attenuation of stellar spectra according to the model of \cite{gordon_quantitative_2003}.
    use Numerical_Constants_Units
    implicit none
    class           (stellarSpectraDustAttenuationTabulated), intent(inout) :: self
    double precision                                        , intent(in   ) :: wavelength      , age, &
         &                                                                     vBandAttenuation
    double precision                                                        :: x

    x=1.0d0/(wavelength/angstromsPerMicron)
    tabulatedAttenuation=vBandAttenuation*self%attenuationTable%interpolate(x)
    return
  end function tabulatedAttenuation
