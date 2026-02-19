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
  Implements a Population III supernovae class with no population III supernovae.
  !!}
  
  !![
  <supernovaePopulationIII name="supernovaePopulationIIIZero">
   <description>
    A Population III supernovae class that has zero population III supernovae.
   </description>
  </supernovaePopulationIII>
  !!]
  type, extends(supernovaePopulationIIIClass) :: supernovaePopulationIIIZero
     !!{
     A Population III supernovae class that has zero population III supernovae.
     !!}
     private
   contains
     procedure :: energyCumulative => zeroEnergyCumulative
  end type supernovaePopulationIIIZero

  interface supernovaePopulationIIIZero
     !!{
     Constructors for the \refClass{supernovaePopulationIIIZero} Population III supernovae class.
     !!}
     module procedure zeroConstructorParameters
  end interface supernovaePopulationIIIZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{supernovaePopulationIIIZero} Population III supernovae class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(supernovaePopulationIIIZero)                :: self
    type(inputParameters            ), intent(inout) :: parameters

    self=supernovaePopulationIIIZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroEnergyCumulative(self,initialMass,age,metallicity) result(energy)
    !!{
    Compute the cumulative energy input from Population III supernovae - there are zero population III supernovae in this model,
    so zero energy.    
    !!}
    implicit none
    class           (supernovaePopulationIIIZero), intent(inout) :: self
    double precision                             , intent(in   ) :: age        , initialMass, &
         &                                                          metallicity
    !$GLC attributes unused :: self, initialMass, age, metallicity
    
    energy=0.0d0
    return
  end function zeroEnergyCumulative
