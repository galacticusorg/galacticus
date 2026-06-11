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

  !!{RST
  Implements a stellar initial mass function class used by the `BPASS <https://bpass.auckland.ac.nz/>`_ library.
  !!}

  !![
  <initialMassFunction name="initialMassFunctionBPASS" docformat="rst">
   <description>
   A stellar initial mass function class used by the `BPASS &lt;https://bpass.auckland.ac.nz/&gt;`_ library:

   .. math::

       \phi(M) \propto \left\{ \begin{array}{ll}
       M^{-1.30} &amp; \hbox{ for } 0.1\mathrm{M}_\odot &lt; M &lt; 0.5\mathrm{M}_\odot \\
       M^{-2.35} &amp; \hbox{ for } 1\mathrm{M}_\odot &lt; M &lt; 120\mathrm{M}_\odot \\
      0 &amp; \hbox {otherwise.} \end{array} \right.
   </description>
  </initialMassFunction>
  !!]
  type, extends(initialMassFunctionPiecewisePowerLaw) :: initialMassFunctionBPASS
     !!{RST
     A stellar initial mass function class used by the `BPASS <https://bpass.auckland.ac.nz/>`_ library.
     !!}
     private
   contains
     procedure :: label => bpassLabel
  end type initialMassFunctionBPASS

  interface initialMassFunctionBPASS
     !!{RST
     Constructors for the ``initialMassFunctionBPASS`` initial mass function class.
     !!}
     module procedure bpassConstructorParameters
     module procedure bpassConstructorInternal
  end interface initialMassFunctionBPASS

contains

  function bpassConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``initialMassFunctionBPASS`` initial mass function class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(initialMassFunctionBPASS)                :: self
    type(inputParameters         ), intent(inout) :: parameters

    self=initialMassFunctionBPASS()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function bpassConstructorParameters

  function bpassConstructorInternal() result(self)
    !!{RST
    Internal constructor for the ``initialMassFunctionBPASS`` initial mass function.
    !!}
    implicit none
    type(initialMassFunctionBPASS):: self

    self%initialMassFunctionPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                    &
         &                                                                         mass    =[+0.10d0,+0.50d0,+1.20d2], &
         &                                                                         exponent=[-1.30d0,-2.35d0        ]  &
         &                                                                        )
    return
  end function bpassConstructorInternal

  function bpassLabel(self)
    !!{RST
    Return a label for this :term:`IMF`.
    !!}
    implicit none
    class(initialMassFunctionBPASS), intent(inout) :: self
    type (varying_string          )                :: bpassLabel
    !$GLC attributes unused :: self

    bpassLabel="BPASS"
    return
  end function bpassLabel
