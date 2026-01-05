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
  Implements a stellar initial mass function class for the \cite{scalo_stellar_1986} \gls{imf}.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionScalo1986">
   <description>
    A stellar initial mass function class for the \cite{scalo_stellar_1986} \gls{imf}:
    \begin{equation}
     \phi(M) \propto \left\{ \begin{array}{ll}
     M^{+1.60} &amp; \hbox{ for } 0.10M_\odot &lt; M &lt; 0.18M_\odot \\
     M^{-1.01} &amp; \hbox{ for } 0.18M_\odot &lt; M &lt; 0.42M_\odot \\
     M^{-2.75} &amp; \hbox{ for } 0.42M_\odot &lt; M &lt; 0.62M_\odot \\
     M^{-2.08} &amp; \hbox{ for } 0.62M_\odot &lt; M &lt; 1.18M_\odot \\
     M^{-3.50} &amp; \hbox{ for } 1.18M_\odot &lt; M &lt; 3.50M_\odot \\
     M^{-2.63} &amp; \hbox{ for } 3.50M_\odot &lt; M &lt; 125M_\odot \\
     0 &amp; \hbox {otherwise.} \end{array} \right.
    \end{equation}
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionScalo1986
     !!{
     A stellar initial mass function class for the \cite{scalo_stellar_1986} \gls{imf}.
     !!}
     private
   contains
     procedure :: label => scalo1986Label
  end type initialMassFunctionScalo1986

  interface initialMassFunctionScalo1986
     !!{
     Constructors for the \refClass{initialMassFunctionScalo1986} initial mass function class.
     !!}
     module procedure scalo1986ConstructorParameters
     module procedure scalo1986ConstructorInternal
  end interface initialMassFunctionScalo1986

contains

  function scalo1986ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{initialMassFunctionScalo1986} initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionScalo1986)                :: self
    type(inputParameters             ), intent(inout) :: parameters

    self=initialMassFunctionScalo1986()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function scalo1986ConstructorParameters

  function scalo1986ConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{initialMassFunctionScalo1986} initial mass function.
    !!}
    implicit none
    type(initialMassFunctionScalo1986):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                                                    &
         &                                                                         mass    =[+0.10d0,+0.18d0,+0.42d0,+0.62d0,+1.18d0,+3.50d0,+1.25d2], &
         &                                                                         exponent=[+1.60d0,-1.01d0,-2.75d0,-2.08d0,-3.50d0,-2.63d0        ]  &
         &                                                                        )
    return
  end function scalo1986ConstructorInternal

  function scalo1986Label(self)
    !!{
    Return a label for this \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionScalo1986), intent(inout) :: self
    type (varying_string              )                :: scalo1986Label
    !$GLC attributes unused :: self

    scalo1986Label="Scalo1986"
    return
  end function scalo1986Label
