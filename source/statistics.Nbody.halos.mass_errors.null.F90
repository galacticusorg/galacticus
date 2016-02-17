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

!% Contains a module which implements a null N-body dark matter halo mass error class.

  !# <nbodyHaloMassError name="nbodyHaloMassErrorNull">
  !#  <description>A null N-body dark matter halo mass error class. Errors are always zero.</description>
  !# </nbodyHaloMassError>
  type, extends(nbodyHaloMassErrorClass) :: nbodyHaloMassErrorNull
     !% A null N-body halo mass error class.
     private
    contains
     final     ::                    nullDestructor
     procedure :: errorFractional => nullErrorFractional
  end type nbodyHaloMassErrorNull

  interface nbodyHaloMassErrorNull
     !% Constructors for the {\normalfont \ttfamily null} N-body halo mass error class.
     module procedure nbodyHaloMassErrorNullParameters
     module procedure nbodyHaloMassErrorNullInternal
  end interface nbodyHaloMassErrorNull

contains

  function nbodyHaloMassErrorNullParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily null} N-body halo mass error class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(nbodyHaloMassErrorNull)                :: nbodyHaloMassErrorNullParameters
    type(inputParameters       ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)
    return
  end function nbodyHaloMassErrorNullParameters

  function nbodyHaloMassErrorNullInternal()
    !% Internal constructor for the {\normalfont \ttfamily null} N-body halo mass error class.
    implicit none
    type(nbodyHaloMassErrorNull) :: nbodyHaloMassErrorNullInternal

    ! Nothing to do.
    return
  end function nbodyHaloMassErrorNullInternal

  subroutine nullDestructor(self)
    !% Destructor for the {\normalfont \ttfamily null} N-body halo mass error class.
    implicit none
    type(nbodyHaloMassErrorNull), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine nullDestructor

  double precision function nullErrorFractional(self,node)
    !% Return the fractional error on the mass of an N-body halo.
    implicit none
    class(nbodyHaloMassErrorNull), intent(inout)          :: self
    type (treeNode              ), intent(inout), pointer :: node

    nullErrorFractional=0.0d0
    return
  end function nullErrorFractional
  
