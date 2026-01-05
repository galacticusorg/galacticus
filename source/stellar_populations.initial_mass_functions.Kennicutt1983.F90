!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implements a stellar initial mass function class for the \cite{kennicutt_rate_1983} \gls{imf}.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionKennicutt1983">
   <description>
    A stellar initial mass function class for the \cite{kennicutt_rate_1983} \gls{imf}:
    \begin{equation}
     \phi(M) \propto \left\{ \begin{array}{ll}
     M^{-1.25} &amp; \hbox{ for } 0.10M_\odot &lt; M &lt; 1.00M_\odot \\
     M^{-2.00} &amp; \hbox{ for } 1.00M_\odot &lt; M &lt; 2.00M_\odot \\
     M^{-2.30} &amp; \hbox{ for } 2.00M_\odot &lt; M &lt; 125M_\odot \\
     0 &amp; \hbox {otherwise.} \end{array} \right.
    \end{equation}
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionKennicutt1983
     !!{
     A stellar initial mass function class for the \cite{kennicutt_rate_1983} \gls{imf}.
     !!}
     private
   contains
     procedure :: label => kennicutt1983Label
  end type initialMassFunctionKennicutt1983

  interface initialMassFunctionKennicutt1983
     !!{
     Constructors for the \refClass{initialMassFunctionKennicutt1983} initial mass function class.
     !!}
     module procedure kennicutt1983ConstructorParameters
     module procedure kennicutt1983ConstructorInternal
  end interface initialMassFunctionKennicutt1983

contains

  function kennicutt1983ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{initialMassFunctionKennicutt1983} initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionKennicutt1983)                :: self
    type(inputParameters                 ), intent(inout) :: parameters

    self=initialMassFunctionKennicutt1983()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function kennicutt1983ConstructorParameters

  function kennicutt1983ConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{initialMassFunctionKennicutt1983} initial mass function.
    !!}
    implicit none
    type(initialMassFunctionKennicutt1983):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                            &
         &                                                                         mass    =[+0.10d0,+1.00d0,+2.00d0,+1.25d2], &
         &                                                                         exponent=[-1.25d0,-2.00d0,-2.30d0        ]  &
         &                                                                        )
    return
  end function kennicutt1983ConstructorInternal

  function kennicutt1983Label(self)
    !!{
    Return a label for this \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionKennicutt1983), intent(inout) :: self
    type (varying_string                  )                :: kennicutt1983Label
    !$GLC attributes unused :: self

    kennicutt1983Label="Kennicutt1983"
    return
  end function kennicutt1983Label
