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
  Implements an atomic dielectronic recombination class which assumes zero rate.
  !!}

  !![
  <atomicRecombinationRateDielectronic name="atomicRecombinationRateDielectronicZero" docformat="rst">
   <description>
   Implements an atomic dielectronic recombination class which assumes zero rate.
   </description>
  </atomicRecombinationRateDielectronic>
  !!]
  type, extends(atomicRecombinationRateDielectronicClass) :: atomicRecombinationRateDielectronicZero
     !!{RST
     Implements an atomic dielectronic recombination class which assumes zero rate.
     !!}
     private
   contains
     procedure :: rate => zeroRate
  end type atomicRecombinationRateDielectronicZero

  interface atomicRecombinationRateDielectronicZero
     !!{RST
     Constructors for the ``atomicRecombinationRateDielectronicZero`` atomic dielectronic recombination rate class.
     !!}
     module procedure zeroConstructorParameters
  end interface atomicRecombinationRateDielectronicZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``atomicRecombinationRateDielectronicZero`` atomic dielectronic recombination rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(atomicRecombinationRateDielectronicZero)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self=atomicRecombinationRateDielectronicZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroRate(self,atomicNumber,electronNumber,temperature)
    !!{RST
    Returns a zero rate of dielectronic recombination.
    !!}
    implicit none
    class           (atomicRecombinationRateDielectronicZero), intent(inout) :: self
    double precision                                         , intent(in   ) :: temperature
    integer                                                  , intent(in   ) :: atomicNumber, electronNumber
    !$GLC attributes unused :: self, atomicNumber, electronNumber, temperature

    zeroRate=0.0d0
    return
  end function zeroRate
