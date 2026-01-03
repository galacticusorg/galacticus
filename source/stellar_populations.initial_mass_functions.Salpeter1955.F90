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
  Implements a stellar initial mass function class for the \cite{salpeter_luminosity_1955} \gls{imf}.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionSalpeter1955">
   <description>
  !!]

  !![
    A stellar initial mass function class for the \cite{salpeter_luminosity_1955} \gls{imf} defined as:
    \begin{equation}
     \phi(M) \propto \left\{ \begin{array}{ll} M^{-2.35} &amp; \hbox{ for } 0.1M_\odot &lt; M &lt; 125M_\odot \\ 0 &amp; \hbox
     {otherwise.} \end{array} \right.
    \end{equation}
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionSalpeter1955
     !!{
     A stellar initial mass function class for the \cite{salpeter_luminosity_1955} \gls{imf}.
     !!}
     private
   contains
     procedure :: label => salpeter1955Label
  end type initialMassFunctionSalpeter1955

  interface initialMassFunctionSalpeter1955
     !!{
     Constructors for the \refClass{initialMassFunctionSalpeter1955} initial mass function class.
     !!}
     module procedure salpeter1955ConstructorParameters
     module procedure salpeter1955ConstructorInternal
  end interface initialMassFunctionSalpeter1955

contains

  function salpeter1955ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{initialMassFunctionSalpeter1955} initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionSalpeter1955)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=initialMassFunctionSalpeter1955()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function salpeter1955ConstructorParameters

  function salpeter1955ConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{initialMassFunctionSalpeter1955} initial mass function.
    !!}
    implicit none
    type(initialMassFunctionSalpeter1955):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                            &
         &                                                                         mass    =[+0.10d0,+1.25d2], &
         &                                                                         exponent=[-2.35d0        ]  &
         &                                                                        )
    return
  end function salpeter1955ConstructorInternal

  function salpeter1955Label(self)
    !!{
    Return a label for this \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionSalpeter1955), intent(inout) :: self
    type (varying_string                 )                :: salpeter1955Label
    !$GLC attributes unused :: self

    salpeter1955Label="Salpeter1955"
    return
  end function salpeter1955Label
