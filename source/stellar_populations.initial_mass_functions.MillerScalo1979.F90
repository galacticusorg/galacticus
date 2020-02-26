!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a stellar initial mass function class for the \cite{miller_initial_1979} \gls{imf}.

  !# <initialMassFunction name="initialMassFunctionMillerScalo1979">
  !#  <description>A stellar initial mass function class for the \cite{miller_initial_1979} \gls{imf}.</description>
  !# </initialMassFunction>
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionMillerScalo1979
     !% A stellar initial mass function class for the \cite{miller_initial_1979} \gls{imf}.
     private
   contains
     procedure :: label => millerScalo1979Label
  end type initialMassFunctionMillerScalo1979

  interface initialMassFunctionMillerScalo1979
     !% Constructors for the {\normalfont \ttfamily millerScalo1979} initial mass function class.
     module procedure millerScalo1979ConstructorParameters
     module procedure millerScalo1979ConstructorInternal
  end interface initialMassFunctionMillerScalo1979

contains

  function millerScalo1979ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily millerScalo1979} initial mass function class which takes a parameter list as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionMillerScalo1979)                :: self
    type(inputParameters                   ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=initialMassFunctionMillerScalo1979()
    return
  end function millerScalo1979ConstructorParameters

  function millerScalo1979ConstructorInternal() result(self)
    !% Internal constructor for the {\normalfont \ttfamily millerScalo1979} initial mass function.
    implicit none
    type(initialMassFunctionMillerScalo1979):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                                    &
         &                                                                         mass    =[+0.10d0,+1.00d0,+2.00d0,+1.00d1,+1.25d2], &
         &                                                                         exponent=[-1.25d0,-2.00d0,-2.30d0,-3.30d0        ]  &
         &                                                                        )
    return
  end function millerScalo1979ConstructorInternal

  function millerScalo1979Label(self)
    !% Return a label for this \gls{imf}.
    implicit none
    class(initialMassFunctionMillerScalo1979), intent(inout) :: self
    type (varying_string                    )                :: millerScalo1979Label
    !GCC$ attributes unused :: self

    millerScalo1979Label="MillerScalo1979"
    return
  end function millerScalo1979Label
