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

  !% An implementation of a spectrum postprocessor that does nothing.

  !# <spectraPostprocessor name="spectraPostprocessorIdentity">
  !#  <description>Performs an indentity postprocessing of spectra.</description>
  !# </spectraPostprocessor>

  type, extends(spectraPostprocessorClass) :: spectraPostprocessorIdentity
     !% An identity spectrum postprocessor.
     private
   contains
     procedure :: apply => identityApply
  end type spectraPostprocessorIdentity

  interface spectraPostprocessorIdentity
     !% Constructors for the identity spectrum postprocessor class.
     module procedure identityDefaultConstructor
  end interface spectraPostprocessorIdentity

contains

  function identityDefaultConstructor()
    !% Default constructor for the identity spectrum postprocessor class.
    implicit none
    type(spectraPostprocessorIdentity), target :: identityDefaultConstructor

    return
  end function identityDefaultConstructor

  subroutine identityApply(self,wavelength,age,redshift,modifier)
    !% Perform an identity postprocessing on a spectrum.
    implicit none
    class           (spectraPostprocessorIdentity), intent(inout) :: self
    double precision                              , intent(in   ) :: age     , redshift, wavelength
    double precision                              , intent(inout) :: modifier

    modifier=modifier
    return
  end subroutine identityApply
