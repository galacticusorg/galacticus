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
  Implements a stellar initial mass function class for the top-heavy stellar initial mass function from \cite{baugh_can_2005}.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionBaugh2005TopHeavy">
   <description>
    A stellar initial mass function class for the top-heavy stellar initial mass function from \cite{baugh_can_2005}:
    The {\normalfont \ttfamily Baugh2005TopHeavy} \gls{imf} is defined by \citep{baugh_can_2005}:
    \begin{equation}
     \phi(M) \propto 
     M^{-1} \hbox{ for } 0.15M_\odot &lt; M &lt; 125M_\odot
    \end{equation}
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionBaugh2005TopHeavy
     !!{
     A stellar initial mass function class for the top-heavy stellar initial mass function from \cite{baugh_can_2005}.
     !!}
     private
   contains
     procedure :: label => baugh2005TopHeavyLabel
  end type initialMassFunctionBaugh2005TopHeavy

  interface initialMassFunctionBaugh2005TopHeavy
     !!{
     Constructors for the \refClass{initialMassFunctionBaugh2005TopHeavy} initial mass function class.
     !!}
     module procedure baugh2005TopHeavyConstructorParameters
     module procedure baugh2005TopHeavyConstructorInternal
  end interface initialMassFunctionBaugh2005TopHeavy

contains

  function baugh2005TopHeavyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{initialMassFunctionBaugh2005TopHeavy} initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionBaugh2005TopHeavy)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=initialMassFunctionBaugh2005TopHeavy()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function baugh2005TopHeavyConstructorParameters

  function baugh2005TopHeavyConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{initialMassFunctionBaugh2005TopHeavy} initial mass function.
    !!}
    implicit none
    type(initialMassFunctionBaugh2005TopHeavy):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                            &
         &                                                                         mass    =[+0.15d0,+1.25d2], &
         &                                                                         exponent=[-1.00d0        ]  &
         &                                                                        )
    return
  end function baugh2005TopHeavyConstructorInternal

  function baugh2005TopHeavyLabel(self)
    !!{
    Return a label for this \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionBaugh2005TopHeavy), intent(inout) :: self
    type (varying_string                      )                :: baugh2005TopHeavyLabel
    !$GLC attributes unused :: self

    baugh2005TopHeavyLabel="Baugh2005TopHeavy"
    return
  end function baugh2005TopHeavyLabel
