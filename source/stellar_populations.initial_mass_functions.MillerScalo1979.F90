!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a stellar initial mass function class for the \cite{miller_initial_1979} \gls{imf}.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionMillerScalo1979">
   <description>
    A stellar initial mass function class for the \cite{miller_initial_1979} \gls{imf}:
    \begin{equation}
     \phi(M) \propto \left\{ \begin{array}{ll}
     M^{-1.25} &amp; \hbox{ for } 0.10M_\odot &lt; M &lt; 1.00M_\odot \\
     M^{-2.00} &amp; \hbox{ for } 1.00M_\odot &lt; M &lt; 2.00M_\odot \\
     M^{-2.30} &amp; \hbox{ for } 2.00M_\odot &lt; M &lt; 10.0M_\odot \\
     M^{-3.30} &amp; \hbox{ for } 10.0M_\odot &lt; M &lt; 125M_\odot \\
     0 &amp; \hbox {otherwise.} \end{array} \right.
    \end{equation}
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionMillerScalo1979
     !!{
     A stellar initial mass function class for the \cite{miller_initial_1979} \gls{imf}.
     !!}
     private
   contains
     procedure :: label => millerScalo1979Label
  end type initialMassFunctionMillerScalo1979

  interface initialMassFunctionMillerScalo1979
     !!{
     Constructors for the \refClass{initialMassFunctionMillerScalo1979} initial mass function class.
     !!}
     module procedure millerScalo1979ConstructorParameters
     module procedure millerScalo1979ConstructorInternal
  end interface initialMassFunctionMillerScalo1979

contains

  function millerScalo1979ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{initialMassFunctionMillerScalo1979} initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionMillerScalo1979)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=initialMassFunctionMillerScalo1979()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function millerScalo1979ConstructorParameters

  function millerScalo1979ConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{initialMassFunctionMillerScalo1979} initial mass function.
    !!}
    implicit none
    type(initialMassFunctionMillerScalo1979):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                                    &
         &                                                                         mass    =[+0.10d0,+1.00d0,+2.00d0,+1.00d1,+1.25d2], &
         &                                                                         exponent=[-1.25d0,-2.00d0,-2.30d0,-3.30d0        ]  &
         &                                                                        )
    return
  end function millerScalo1979ConstructorInternal

  function millerScalo1979Label(self)
    !!{
    Return a label for this \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionMillerScalo1979), intent(inout) :: self
    type (varying_string                    )                :: millerScalo1979Label
    !$GLC attributes unused :: self

    millerScalo1979Label="MillerScalo1979"
    return
  end function millerScalo1979Label
