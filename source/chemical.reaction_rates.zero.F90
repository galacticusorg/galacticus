!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% An implementation of calculations of chemical reaction rates which assumes zero rates.
  
  !# <chemicalReactionRate name="chemicalReactionRateZero">
  !#  <description>A chemical reaction rate class in which all rates are zero.</description>
  !# </chemicalReactionRate>
  type, extends(chemicalReactionRateClass) :: chemicalReactionRateZero
     !% A chemical reaction rate class in which all rates are zero.
     private
   contains
     procedure :: rates => zeroRates
  end type chemicalReactionRateZero

  interface chemicalReactionRateZero
     !% Constructors for the {\normalfont \ttfamily zero} chemical reaction rates class.
     module procedure zeroConstructorParameters
  end interface chemicalReactionRateZero
  
contains

  function zeroConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily zero} chemical reaction rates class which takes a parameter set as
    !% input.
    use Input_Parameters
    implicit none
    type(chemicalReactionRateZero)                :: self
    type(inputParameters         ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=chemicalReactionRateZero()
    return
  end function zeroConstructorParameters

  subroutine zeroRates(self,temperature,chemicalDensity,radiation,chemicalRates)
    !% Return zero rates of chemical reactions.
    implicit none
    class           (chemicalReactionRateZero), intent(inout) :: self
    type            (chemicalAbundances      ), intent(in   ) :: chemicalDensity
    double precision                          , intent(in   ) :: temperature
    type            (radiationStructure      ), intent(in   ) :: radiation
    type            (chemicalAbundances      ), intent(inout) :: chemicalRates
    !GCC$ attributes unused :: self, chemicalDensity, temperature, radiation
    
    call chemicalRates%reset()
    return
  end subroutine zeroRates
