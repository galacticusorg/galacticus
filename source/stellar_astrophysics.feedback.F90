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
Contains a module which implements a class that performs calculations of stellar feedback.
!!}

module Stellar_Feedback
  !!{
  Implements a class that performs calculations of stellar feedback.
  !!}
  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionClass
  implicit none
  private

  !![
  <functionClass>
   <name>stellarFeedback</name>
   <descriptiveName>Stellar Feedback</descriptiveName>
   <description>
    Class providing models of stellar feedback.
   </description>
   <default>standard</default>
   <method name="energyInputCumulative" >
    <description>Return the cumulative energy input from a stellar population of the given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily age}, and {\normalfont \ttfamily metallicity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class           (initialMassFunctionClass), intent(inout) :: initialMassFunction_                  </argument>
    <argument>double precision                          , intent(in   ) :: initialMass         , age, metallicity</argument>
   </method>
  </functionClass>
  !!]

  ! Canonical value of the total energy input from a single stellar population of 1 M☉ after infinite time. All feedback
  ! calculations which don't specifically use the energy input should be scaled to this value if they want to have the correct
  ! time and IMF dependencies. Value was computed for a Salpeter IMF. Units are M☉ (km/s)².
  double precision, parameter, public :: feedbackEnergyInputAtInfinityCanonical=4.517d5

end module Stellar_Feedback
