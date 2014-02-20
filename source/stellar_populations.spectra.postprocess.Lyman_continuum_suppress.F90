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

  !% An implementation of a spectrum postprocessor that suppresses the Lyman continuum.

  !# <spectraPostprocessor name="spectraPostprocessorLycSuppress">
  !#  <description>Suppress the Lyman continuum in stellar populations.</description>
  !# </spectraPostprocessor>

  type, extends(spectraPostprocessorClass) :: spectraPostprocessorLycSuppress
     !% An lycSuppress spectrum postprocessor.
     private
   contains
     procedure :: apply => lycSuppressApply
  end type spectraPostprocessorLycSuppress

  interface spectraPostprocessorLycSuppress
     !% Constructors for the {\tt lycSuppress} spectrum postprocessor class.
     module procedure lycSuppressDefaultConstructor
  end interface spectraPostprocessorLycSuppress

contains

  function lycSuppressDefaultConstructor()
    !% Default constructor for the {\tt lycSuppress} spectrum postprocessor class.
    implicit none
    type(spectraPostprocessorLycSuppress), target :: lycSuppressDefaultConstructor
    
    return
  end function lycSuppressDefaultConstructor

  subroutine lycSuppressApply(self,wavelength,age,redshift,modifier)
    !% Suppress the Lyman continuum in a spectrum.
    use Numerical_Constants_Atomic
    implicit none
    class           (spectraPostprocessorLycSuppress), intent(inout) :: self
    double precision                                 , intent(in   ) :: age     , redshift, wavelength
    double precision                                 , intent(inout) :: modifier

    if (wavelength < lymanSeriesLimitWavelengthHydrogen) modifier=0.0d0
    return
  end subroutine lycSuppressApply
