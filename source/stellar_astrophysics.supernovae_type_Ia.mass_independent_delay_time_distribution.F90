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
  Implements a supernovae type Ia class for delay time distributions that are independent of progenitor mass.
  !!}
  
  !![
  <supernovaeTypeIa name="supernovaeTypeIaMassIndependentDTD" abstract="yes">
   <description>
    A supernovae type Ia class for delay time distributions that are independent of progenitor mass.
   </description>
  </supernovaeTypeIa>
  !!]
  type, abstract, extends(supernovaeTypeIaFixedYield) :: supernovaeTypeIaMassIndependentDTD
     !!{
     A supernovae type Ia class for delay time distributions that are independent of progenitor mass.
     !!}
     private
   contains    
     !![
     <methods>
       <method method="numberCumulative" description="Return the cumulative number of type Ia supernovae per Solar mass of stars formed at a given population age and metallicity."/>
     </methods>
     !!]
     procedure                                     :: massInitialRange => massIndependentDTDMassInitialRange
     procedure                                     :: number           => massIndependentDTDNumber
     procedure(numberCumulativeTemplate), deferred :: numberCumulative
  end type supernovaeTypeIaMassIndependentDTD

  abstract interface
     double precision function numberCumulativeTemplate(self,age,metallicity)
       !!{
       Interface for cumulative number of Type Ia SNe.
       !!}
       import supernovaeTypeIaMassIndependentDTD
       class           (supernovaeTypeIaMassIndependentDTD), intent(inout), target :: self
       double precision                                    , intent(in   )         :: age , metallicity
     end function numberCumulativeTemplate
  end interface

contains

  subroutine massIndependentDTDMassInitialRange(self,initialMassFunction_,age,metallicity,massInitialMinimum,massInitialMaximum)
    !!{
    Return the range of initial stellar masses contributing to the Type Ia population.
    !!}
    implicit none
    class           (supernovaeTypeIaMassIndependentDTD), intent(inout) :: self
    class           (initialMassFunctionClass          ), intent(inout) :: initialMassFunction_
    double precision                                    , intent(in   ) :: age                 , metallicity
    double precision                                    , intent(  out) :: massInitialMinimum  , massInitialMaximum
    !$GLC attributes unused :: self, age, metallicity
    
    massInitialMinimum=initialMassFunction_%massMinimum()
    massInitialMaximum=initialMassFunction_%massMaximum()
    return
  end subroutine massIndependentDTDMassInitialRange
  
  double precision function massIndependentDTDNumber(self,initialMassFunction_,initialMass,age,metallicity) result(number)
    !!{
    Compute the cumulative number of Type Ia supernovae originating per unit interval of secondary star mass with given
    {\normalfont \ttfamily initialMass} and {\normalfont \ttfamily metallicity} after a time {\normalfont \ttfamily age}. Here we
    assume that the total number of Type Ias is specified independent of secondary star mass.
    !!}
    implicit none
    class           (supernovaeTypeIaMassIndependentDTD), intent(inout), target :: self
    class           (initialMassFunctionClass          ), intent(inout), target :: initialMassFunction_
    double precision                                    , intent(in   )         :: age                 , initialMass, &
         &                                                                         metallicity
    !$GLC attributes unused :: initialMass

    number=self%numberCumulative(age,metallicity)/initialMassFunction_%numberCumulative(initialMassFunction_%massMinimum(),initialMassFunction_%massMaximum())
    return
  end function massIndependentDTDNumber
  
