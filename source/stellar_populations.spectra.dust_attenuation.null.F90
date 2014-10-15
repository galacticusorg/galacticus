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

  !% Implements calculations of null dust attenuation of stellar spectra.
  
  !# <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationNull">
  !#  <description>Returns a null dust attenuation.</description>
  !# </stellarSpectraDustAttenuation>

  type, extends(stellarSpectraDustAttenuationClass) :: stellarSpectraDustAttenuationNull
     !% A class implementing null dust attenuation of stellar spectra.
     private
   contains
     procedure :: attenuation => nullAttenuation
  end type stellarSpectraDustAttenuationNull

  interface stellarSpectraDustAttenuationNull
     !% Constructors for the ``null'' stellar spectra dust attenuation class.
     module procedure nullDefaultConstructor
  end interface stellarSpectraDustAttenuationNull

contains

  function nullDefaultConstructor()
    !% Default constructor for the ``null'' stellar spectra dust attenuatio class.
    implicit none
    type(stellarSpectraDustAttenuationNull) :: nullDefaultConstructor

    return
  end function nullDefaultConstructor

  double precision function nullAttenuation(self,wavelength,age,vBandAttenuation)
    !% Return a null attenuation.
    implicit none
    class           (stellarSpectraDustAttenuationNull), intent(inout) :: self
    double precision                                   , intent(in   ) :: wavelength      , age, &
         &                                                                vBandAttenuation

    nullAttenuation=0.0d0
    return
  end function nullAttenuation
