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
  Implements a stellar feedback class which assumes instantaneous injection of energy at the canonical rate. Primarily intended
  for testing purposes. 
  !!}

  !![
  <stellarFeedback name="stellarFeedbackInstantCanonical">
   <description>
    A stellar feedback class which assumes that energy input from a stellar population occurs instantly, and at the canonical rate.
   </description>
  </stellarFeedback>
  !!]
  type, extends(stellarFeedbackClass) :: stellarFeedbackInstantCanonical
     !!{
     A stellar feedback class which assumes instantaneous injection of energy at the canonical rate. Primarily intended
     for testing purposes.
     !!}
     private
   contains
     procedure :: energyInputCumulative => instantCanonicalEnergyInputCumulative
  end type stellarFeedbackInstantCanonical

  interface stellarFeedbackInstantCanonical
     !!{
     Constructors for the \refClass{stellarFeedbackInstantCanonical} stellar feedback class.
     !!}
     module procedure instantCanonicalConstructorParameters
  end interface stellarFeedbackInstantCanonical

contains

  function instantCanonicalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarFeedbackInstantCanonical} stellar feedback class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(stellarFeedbackInstantCanonical)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=stellarFeedbackInstantCanonical()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function instantCanonicalConstructorParameters

  double precision function instantCanonicalEnergyInputCumulative(self,initialMassFunction_,initialMass,age,metallicity)
    !!{
    Compute the cumulative energy input from a star of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily age} and {\normalfont \ttfamily metallicity}.
    !!}
    implicit none
    class           (stellarFeedbackInstantCanonical), intent(inout), target :: self
    class           (initialMassFunctionClass       ), intent(inout)         :: initialMassFunction_
    double precision                                 , intent(in   )         :: age                 , initialMass, &
         &                                                                      metallicity
    !$GLC attributes unused :: self, initialMassFunction_, initialMass, age, metallicity
    
    instantCanonicalEnergyInputCumulative=feedbackEnergyInputAtInfinityCanonical
    return
  end function instantCanonicalEnergyInputCumulative
