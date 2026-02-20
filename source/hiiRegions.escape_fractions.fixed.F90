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

  !+ Contributions to this file made by Sachi Weerasooriya

  !!{
  Implementation of a fixed escape fraction from HII regions.
  !!}

  !![
  <hiiRegionEscapeFraction name="hiiRegionEscapeFractionFixed">
   <description>
    Computes the escape fraction of hydrogen ionizing photons from HII regions. A fixed escape fraction of
    $f_\mathrm{esc}${\normalfont \ttfamily [escapeFraction]} is assumed for HII regions with ages less than
    $\tau_\mathrm{limit}=${\normalfont \ttfamily {ageLimit}}.
   </description>
  </hiiRegionEscapeFraction>
  !!]
  type, extends(hiiRegionEscapeFractionClass) :: hiiRegionEscapeFractionFixed
     !!{
     Implementation of a fixed escape fraction from HII regions.
     !!}
     private
     double precision :: escapeFraction_, ageLimit
   contains
     procedure :: escapeFraction => escapeFractionFixed
  end type hiiRegionEscapeFractionFixed

  interface hiiRegionEscapeFractionFixed
     !!{
     Constructors for the hiiRegionEscapeFractionFixed for HII region luminosity function class.
     !!}
     module procedure escapeFractionFixedConstructorParameters
     module procedure escapeFractionFixedConstructorInternal
  end interface hiiRegionEscapeFractionFixed

contains

  function escapeFractionFixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hiiRegionEscapeFractionFixed} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionEscapeFractionFixed)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: escapeFraction, ageLimit
    !![
    <inputParameter>
      <name>escapeFraction</name>
      <defaultValue>0.006d0</defaultValue>
      <description> Escape fraction of ionizing photons from young HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ageLimit</name>
      <defaultValue>0.03d0</defaultValue>
      <description>The age beyond which all ionizing photons are assumed to escape from HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionEscapeFractionFixed(escapeFraction,ageLimit)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function escapeFractionFixedConstructorParameters

  function escapeFractionFixedConstructorInternal(escapeFraction_,ageLimit) result(self)
    !!{
    Internal constructor for the \refClass{hiiRegionEscapeFractionFixed} 
    !!}
    
    implicit none
    type            (hiiRegionEscapeFractionFixed)                :: self
    double precision                              , intent(in   ) :: escapeFraction_, ageLimit
    !![
    <constructorAssign variables="escapeFraction_, ageLimit"/>
    !!]

    return
  end function escapeFractionFixedConstructorInternal

  double precision function escapeFractionFixed(self,ageHIIRegion) result(escapeFraction)
    !!{
    Computes the escape fraction.
    !!}
    implicit none
    class           (hiiRegionEscapeFractionFixed), intent(inout) :: self
    double precision                              , intent(in   ) :: ageHIIRegion

    if (ageHIIRegion >= self%ageLimit) then
      escapeFraction=1.0d0
    else
      escapeFraction=self%escapeFraction_
    end if
    return
  end function escapeFractionFixed
