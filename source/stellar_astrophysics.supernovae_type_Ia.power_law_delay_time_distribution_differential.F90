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
  Implements a supernovae type Ia class with a power-law delay time distribution.
  !!}

  !![
  <supernovaeTypeIa name="supernovaeTypeIaPowerLawDTDDifferential">
   <description>
    A supernovae type Ia class with a power-law delay time distribution. This class implements the differential form of the delay time distribution and is intended primarily for testing.
   </description>
  </supernovaeTypeIa>
  !!]
  type, extends(supernovaeTypeIaDifferentialDTD) :: supernovaeTypeIaPowerLawDTDDifferential
     !!{
     A supernovae type Ia class with a power-law delay time distribution.
     !!}
     private
     double precision :: timeMinimum  , exponent, &
          &              normalization
   contains    
     procedure :: numberDifferential => powerLawDTDDifferentialNumberDifferential
  end type supernovaeTypeIaPowerLawDTDDifferential

  interface supernovaeTypeIaPowerLawDTDDifferential
     !!{
     Constructors for the \refClass{supernovaeTypeIaPowerLawDTDDifferential} supernovae type Ia class.
     !!}
     module procedure powerLawDTDDifferentialConstructorParameters
     module procedure powerLawDTDDifferentialConstructorInternal
  end interface supernovaeTypeIaPowerLawDTDDifferential

contains

  function powerLawDTDDifferentialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{supernovaeTypeIaPowerLawDTDDifferential} supernovae type Ia class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (supernovaeTypeIaPowerLawDTDDifferential)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    double precision                                                         :: timeMinimum  , exponent, &
          &                                                                     normalization

    !![
    <inputParameter>
      <name>timeMinimum</name>
      <source>parameters</source>
      <defaultValue>40.0d-3</defaultValue>
      <defaultSource>\citep[][field galaxy DTD]{freundlich_delay_2012}</defaultSource>
      <description>The minimum time before which the delay time distribution is zero.</description>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <source>parameters</source>
      <defaultValue>-1.07d0</defaultValue>
      <defaultSource>\citep[][field galaxy DTD]{freundlich_delay_2012}</defaultSource>
      <description>The exponent $\alpha$ in the delay time distribution.</description>
    </inputParameter>
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <defaultValue>0.21d-3</defaultValue>
      <defaultSource>\citep[][field galaxy DTD]{freundlich_delay_2012}</defaultSource>
      <description>The normalization $R_1$ of the delay time distribution at 1~Gyr in units of Gyr$^{-1}\,\mathrm{M}_\odot^{-1}$.</description>
    </inputParameter>
    !!]
    self=supernovaeTypeIaPowerLawDTDDifferential(timeMinimum,exponent,normalization)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawDTDDifferentialConstructorParameters

  function powerLawDTDDifferentialConstructorInternal(timeMinimum,exponent,normalization) result(self)
    !!{
    Internal constructor for the \refClass{supernovaeTypeIaPowerLawDTDDifferential} supernovae type Ia class.
    !!}
    implicit none
    type            (supernovaeTypeIaPowerLawDTDDifferential)                :: self
    double precision                                         , intent(in   ) :: timeMinimum  , exponent, &
          &                                                                     normalization
    !![
    <constructorAssign variables="timeMinimum, exponent, normalization"/>
    !!]
    
    self%initialized=.false.
    return
  end function powerLawDTDDifferentialConstructorInternal
  
  double precision function powerLawDTDDifferentialTimeDelayMinimum(self) result(time)
    !!{
    Compute the minimum time in the type Ia supernovae delay time distribution.
    !!}
    implicit none
    class(supernovaeTypeIaPowerLawDTDDifferential), intent(inout) :: self

    time=self%timeMinimum
    return
  end function powerLawDTDDifferentialTimeDelayMinimum
  
  double precision function powerLawDTDDifferentialNumberDifferential(self,age,metallicity) result(number)
    !!{
    Compute the differential number of Type Ia SNe assuming a power-law delay time distribution.
    !!}
    implicit none
    class           (supernovaeTypeIaPowerLawDTDDifferential), intent(inout) :: self
    double precision                                         , intent(in   ) :: age , metallicity
    !$GLC attributes unused :: metallicity
    
    if (age > self%timeMinimum) then
       number=+self%normalization                &
            & *     age          **self%exponent
    else
       number=+0.0d0
    end if
    return
  end function powerLawDTDDifferentialNumberDifferential
  
