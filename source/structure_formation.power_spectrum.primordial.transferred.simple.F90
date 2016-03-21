!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% A simple transferred primordial power spectrum class.

  use Power_Spectra_Primordial
  use Transfer_Functions, only : transferFunction, transferFunctionClass
  
  !# <powerSpectrumPrimordialTransferred name="powerSpectrumPrimordialTransferredSimple" defaultThreadPrivate="yes">
  !#  <description>Implements a simple transferred primordial power spectrum.</description>
  !# </powerSpectrumPrimordialTransferred>
  type, extends(powerSpectrumPrimordialTransferredClass) :: powerSpectrumPrimordialTransferredSimple
     !% A simple transferred primordial power spectrum class.
     private
     class(transferFunctionClass       ), pointer :: transferFunction_        => null()
     class(powerSpectrumPrimordialClass), pointer :: powerSpectrumPrimordial_ => null()
   contains
     final     ::                          simpleDestructor
     procedure :: power                 => simplePower
     procedure :: logarithmicDerivative => simpleLogarithmicDerivative
     procedure :: descriptor            => simpleDescriptor
  end type powerSpectrumPrimordialTransferredSimple

  interface powerSpectrumPrimordialTransferredSimple
     !% Constructors for the ``simple'' transferred primordial power spectrum class.
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface powerSpectrumPrimordialTransferredSimple

contains

  function simpleConstructorParameters(parameters)
    !% Constructor for the ``simple'' transferred primordial power spectrum class which takes a
    !% parameter set as input.
    use Input_Parameters2
    implicit none
    type(powerSpectrumPrimordialTransferredSimple)                :: simpleConstructorParameters
    type(inputParameters                         ), intent(inout) :: parameters

    !# <objectBuilder class="powerSpectrumPrimordial" name="simpleConstructorParameters%powerSpectrumPrimordial_" source="parameters"/>
    !# <objectBuilder class="transferFunction"        name="simpleConstructorParameters%transferFunction_"        source="parameters"/>
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(powerSpectrumPrimordial_,transferFunction_)
    !% Internal constructor for the ``simple'' transferred primordial power spectrum class.
    implicit none
    type (powerSpectrumPrimordialTransferredSimple)                        :: simpleConstructorInternal
    class(powerSpectrumPrimordialClass            ), intent(in   ), target :: powerSpectrumPrimordial_
    class(transferFunctionClass                   ), intent(in   ), target :: transferFunction_

    simpleConstructorInternal%powerSpectrumPrimordial_ => powerSpectrumPrimordial_
    simpleConstructorInternal%transferFunction_        => transferFunction_
    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !% Destructor for the ``simple'' transferred primordial power spectrum class.
    implicit none
    type(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self

    !# <objectDestructor name="self%transferFunction_"       />
    !# <objectDestructor name="self%powerSpectrumPrimordial_"/>
   return
  end subroutine simpleDestructor

  double precision function simplePower(self,wavenumber)
    !% Return the transferred primordial power spectrum at the given {\normalfont \ttfamily
    !% wavenumber}.
    implicit none
    class           (powerSpectrumPrimordialTransferredSimple), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber

    simplePower=+self%transferFunction_       %value(wavenumber)**2 &
         &      *self%powerSpectrumPrimordial_%power(wavenumber)
    return
  end function simplePower

  double precision function simpleLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transferred primordial power spectrum at the
    !% given {\normalfont \ttfamily wavenumber}.
    implicit none
    class           (powerSpectrumPrimordialTransferredSimple), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber

    simpleLogarithmicDerivative=+2.0d0                                                           &
         &                      *self%transferFunction_       %logarithmicDerivative(wavenumber) &
         &                      +self%powerSpectrumPrimordial_%logarithmicDerivative(wavenumber)
    return
  end function simpleLogarithmicDerivative

  subroutine simpleDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self
    type (inputParameters                         ), intent(inout) :: descriptor
    type (inputParameters                         )                :: subParameters

    call descriptor%addParameter("powerSpectrumPrimordialTransferredMethod","simple")
    subParameters=descriptor%subparameters("powerSpectrumPrimordialTransferredMethod")
    call self%transferFunction_       %descriptor(subParameters)
    call self%powerSpectrumPrimordial_%descriptor(subParameters)
    return
  end subroutine simpleDescriptor
