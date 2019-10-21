!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a stellar initial mass function class for the \cite{kennicutt_rate_1983} \gls{imf}.

  !# <initialMassFunction name="initialMassFunctionKennicutt1983">
  !#  <description>A stellar initial mass function class for the \cite{kennicutt_rate_1983} \gls{imf}.</description>
  !# </initialMassFunction>
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionKennicutt1983
     !% A stellar initial mass function class for the \cite{kennicutt_rate_1983} \gls{imf}.
     private
   contains
     procedure :: label => kennicutt1983Label
  end type initialMassFunctionKennicutt1983

  interface initialMassFunctionKennicutt1983
     !% Constructors for the {\normalfont \ttfamily kennicutt1983} initial mass function class.
     module procedure kennicutt1983ConstructorParameters
     module procedure kennicutt1983ConstructorInternal
  end interface initialMassFunctionKennicutt1983

contains

  function kennicutt1983ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily kennicutt1983} initial mass function class which takes a parameter list as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionKennicutt1983)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=initialMassFunctionKennicutt1983()
    return
  end function kennicutt1983ConstructorParameters

  function kennicutt1983ConstructorInternal() result(self)
    !% Internal constructor for the {\normalfont \ttfamily kennicutt1983} initial mass function.
    implicit none
    type(initialMassFunctionKennicutt1983):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                            &
         &                                                                         mass    =[+0.10d0,+1.00d0,+2.00d0,+1.25d2], &
         &                                                                         exponent=[-1.25d0,-2.00d0,-2.30d0        ]  &
         &                                                                        )
    return
  end function kennicutt1983ConstructorInternal

  function kennicutt1983Label(self)
    !% Return a label for this \gls{imf}.
    implicit none
    class(initialMassFunctionKennicutt1983), intent(inout) :: self
    type (varying_string                  )                :: kennicutt1983Label
    !GCC$ attributes unused :: self

    kennicutt1983Label="Kennicutt1983"
    return
  end function kennicutt1983Label
