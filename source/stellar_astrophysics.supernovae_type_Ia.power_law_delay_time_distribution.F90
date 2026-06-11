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
  Implements a supernovae type Ia class with a power-law delay time distribution.
  !!}

  !![
  <supernovaeTypeIa name="supernovaeTypeIaPowerLawDTD" docformat="rst">
   <description>
   A supernovae type Ia class with a power-law delay time distribution.
   </description>
  </supernovaeTypeIa>
  !!]
  type, extends(supernovaeTypeIaMassIndependentDTD) :: supernovaeTypeIaPowerLawDTD
     !!{RST
     A supernovae type Ia class with a power-law delay time distribution.
     !!}
     private
     double precision :: timeMinimum  , exponent, &
          &              normalization
   contains    
     procedure :: numberCumulative => powerLawDTDNumberCumulative
  end type supernovaeTypeIaPowerLawDTD

  interface supernovaeTypeIaPowerLawDTD
     !!{RST
     Constructors for the ``supernovaeTypeIaPowerLawDTD`` supernovae type Ia class.
     !!}
     module procedure powerLawDTDConstructorParameters
     module procedure powerLawDTDConstructorInternal
  end interface supernovaeTypeIaPowerLawDTD

contains

  function powerLawDTDConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``supernovaeTypeIaPowerLawDTD`` supernovae type Ia class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (supernovaeTypeIaPowerLawDTD)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: timeMinimum  , exponent, &
          &                                                         normalization

    !![
    <inputParameter docformat="rst">
      <name>timeMinimum</name>
      <source>parameters</source>
      <defaultValue>40.0d-3</defaultValue>
      <defaultSource>
      :cite:p:`freundlich_delay_2021`
      </defaultSource>
      <description>
      The minimum time before which the delay time distribution is zero.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponent</name>
      <source>parameters</source>
      <defaultValue>-1.07d0</defaultValue>
      <defaultSource>
      :cite:p:`freundlich_delay_2021`
      </defaultSource>
      <description>
      The power-law index :math:`\alpha` in the Type Ia supernova delay time distribution :math:`\mathrm{d}N/\mathrm{d}t \propto t^\alpha`, describing how the rate declines with time after star formation; the default value of :math:`-1.07` follows :cite:t:`freundlich_delay_2021`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>normalization</name>
      <source>parameters</source>
      <defaultValue>0.21d-3</defaultValue>
      <defaultSource>
      :cite:p:`freundlich_delay_2021`
      </defaultSource>
      <description>
      The normalization :math:`R_1` of the delay time distribution at 1 Gyr in units of Gyr\ :math:`^{-1}\,\mathrm{M}_\odot^{-1}`.
      </description>
    </inputParameter>
    !!]
    self=supernovaeTypeIaPowerLawDTD(timeMinimum,exponent,normalization)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawDTDConstructorParameters

  function powerLawDTDConstructorInternal(timeMinimum,exponent,normalization) result(self)
    !!{RST
    Internal constructor for the ``supernovaeTypeIaPowerLawDTD`` supernovae type Ia class.
    !!}
    implicit none
    type            (supernovaeTypeIaPowerLawDTD)                :: self
    double precision                             , intent(in   ) :: timeMinimum  , exponent, &
          &                                                         normalization
    !![
    <constructorAssign variables="timeMinimum, exponent, normalization"/>
    !!]
    
    self%initialized=.false.
    return
  end function powerLawDTDConstructorInternal
  
  double precision function powerLawDTDNumberCumulative(self,age,metallicity) result(number)
    !!{RST
    Compute the cumulative number of Type Ia SNe assuming a power-law delay time distribution.
    !!}
    implicit none
    class           (supernovaeTypeIaPowerLawDTD), intent(inout), target :: self
    double precision                             , intent(in   )         :: age , metallicity
    !$GLC attributes unused :: metallicity
    
    if (age > self%timeMinimum) then
       number =+  self%normalization                        &
            &  /                      (1.0d0+self%exponent) &
            &  *(                                           &
            &    +     age          **(1.0d0+self%exponent) &
            &    -self%timeMinimum  **(1.0d0+self%exponent) &
            &   )
    else
       number=+0.0d0
    end if
    return
  end function powerLawDTDNumberCumulative
  
