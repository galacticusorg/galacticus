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
  Implements a stellar initial mass function class for the \cite{kroupa_variation_2001} \gls{imf}.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionKroupa2001">
   <description>
    A stellar initial mass function class for the \cite{kroupa_variation_2001} \gls{imf}.
    \begin{equation}
     \phi(M) \propto \left\{ \begin{array}{ll}
     M^{-0.3} &amp; \hbox{ for } 0.01M_\odot &lt; M &lt; 0.08M_\odot \\
     M^{-1.8} &amp; \hbox{ for } 0.08M_\odot &lt; M &lt; 0.5M_\odot \\
     M^{-2.7} &amp; \hbox{ for } 0.5M_\odot &lt; M &lt; 1M_\odot \\
     M^{-2.3} &amp; \hbox{ for } 1M_\odot &lt; M &lt; 125M_\odot \\
    0 &amp; \hbox {otherwise.} \end{array} \right.
    \end{equation}
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionKroupa2001
     !!{
     A stellar initial mass function class for the \cite{kroupa_variation_2001} \gls{imf}.
     !!}
     private
   contains
     procedure :: label => kroupa2001Label
  end type initialMassFunctionKroupa2001

  interface initialMassFunctionKroupa2001
     !!{
     Constructors for the \refClass{initialMassFunctionKroupa2001} initial mass function class.
     !!}
     module procedure kroupa2001ConstructorParameters
     module procedure kroupa2001ConstructorInternal
  end interface initialMassFunctionKroupa2001

contains

  function kroupa2001ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{initialMassFunctionKroupa2001} initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionKroupa2001)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=initialMassFunctionKroupa2001()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function kroupa2001ConstructorParameters

  function kroupa2001ConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{initialMassFunctionKroupa2001} initial mass function.
    !!}
    implicit none
    type(initialMassFunctionKroupa2001):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                                    &
         &                                                                         mass    =[+0.01d0,+0.08d0,+0.50d0,+1.00d0,+1.25d2], &
         &                                                                         exponent=[-0.30d0,-1.80d0,-2.70d0,-2.30d0        ]  &
         &                                                                        )
    return
  end function kroupa2001ConstructorInternal

  function kroupa2001Label(self)
    !!{
    Return a label for this \gls{imf}.
    !!}
    implicit none
    class(initialMassFunctionKroupa2001), intent(inout) :: self
    type (varying_string               )                :: kroupa2001Label
    !$GLC attributes unused :: self

    kroupa2001Label="Kroupa2001"
    return
  end function kroupa2001Label
