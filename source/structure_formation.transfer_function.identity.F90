!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements an identity transfer function class.

  !# <transferFunction name="transferFunctionIdentity">
  !#  <description>Provides an identity transfer function.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionIdentity
     !% A identity transfer function class.
     private
   contains
     final     ::                          identityDestructor
     procedure :: value                 => identityValue
     procedure :: logarithmicDerivative => identityLogarithmicDerivative
     procedure :: halfModeMass          => identityHalfModeMass
  end type transferFunctionIdentity

  interface transferFunctionIdentity
     !% Constructors for the identity transfer function class.
     module procedure identityConstructorParameters
     module procedure identityConstructorInternal
  end interface transferFunctionIdentity

contains

  function identityConstructorParameters(parameters)
    !% Constructor for the identity transfer function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(transferFunctionIdentity)                :: identityConstructorParameters
    type(inputParameters         ), intent(in   ) :: parameters
    
    return
  end function identityConstructorParameters

  function identityConstructorInternal()
    !% Internal constructor for the identity transfer function class.
    implicit none
    type(transferFunctionIdentity) :: identityConstructorInternal

    return
  end function identityConstructorInternal

  elemental subroutine identityDestructor(self)
    !% Destructor for the identity transfer function class.
    implicit none
    type(transferFunctionIdentity), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine identityDestructor

  double precision function identityValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionIdentity), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    identityValue=1.0d0
    return
  end function identityValue

  double precision function identityLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionIdentity), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    identityLogarithmicDerivative=0.0d0
    return
  end function identityLogarithmicDerivative

  double precision function identityHalfModeMass(self)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function. Not supported in this implementation.
    use Galacticus_Error
    implicit none
    class(transferFunctionIdentity), intent(inout) :: self

    call Galacticus_Error_Report('identityHalfModeMass','not supported by this implementation')
    return
  end function identityHalfModeMass
