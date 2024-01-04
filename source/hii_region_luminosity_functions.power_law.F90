!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of a power law luminosity function which scales with a given exponent.
  !!}
  !![
  <hiiRegionLuminosityFunction name="hiiRegionLuminosityFunctionPowerLaw">
   <description>
    A luminosity function class in which the luminosity scales with given exponent.
    \begin{equation}
     \phi(Q_H) = Q_H^{-\alpha+1} / (1-\alpha),
    \end{equation}
    
   Where $\alpha$ is the exponent and $Q_H$ is the rate of photon production rate.
   </description>
  </hiiRegionLuminosityFunction>
  !!]
  type, extends(hiiRegionLuminosityFunctionClass) :: hiiRegionLuminosityFunctionPowerLaw
     !!{
     Implementation of a luminosity function that scales to a power law
     !!}
     private
     double precision :: rateHydrogenIonizingPhotonsMinimum , &
          &              rateHydrogenIonizingPhotonsMaximum, exponent
     
   contains
     procedure :: cumulativeLuminosity => powerLawCumulativeLuminosity
  end type hiiRegionLuminosityFunctionPowerLaw

  interface hiiRegionLuminosityFunctionPowerLaw
     !!{
     Constructors for the hiiRegionLuminosityFunctionPowerLaw for star HII region luminosity function class.
     !!}
     module procedure powerLawCumulativeLuminosityConstructorParameters
     module procedure powerLawCumulativeLuminosityConstructorInternal
  end interface hiiRegionLuminosityFunctionPowerLaw

contains

  function powerLawCumulativeLuminosityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily dynamicalTime} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionLuminosityFunctionPowerLaw)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum,exponent
    !![
    <inputParameter>
      <name>exponent</name>
      <defaultValue>1.80d0</defaultValue>
      <description> Exponent of the differential luminosity function. </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionLuminosityFunctionPowerLaw(rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum,exponent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawCumulativeLuminosityConstructorParameters

  function powerLawCumulativeLuminosityConstructorInternal(rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum,exponent) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily emissionLinePowerLaw} luminosity function class.
    !!}
    
    implicit none
    type            (hiiRegionLuminosityFunctionPowerLaw)               :: self
    double precision                                     , intent(in   ) :: rateHydrogenIonizingPhotonsMinimum, &
         &                                                                  rateHydrogenIonizingPhotonsMaximum,exponent
    !![
    <constructorAssign variables=" rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum, exponent"/>
    !!]
    !! Assign values to the object's components
    self%rateHydrogenIonizingPhotonsMinimum = rateHydrogenIonizingPhotonsMinimum
    self%rateHydrogenIonizingPhotonsMaximum = rateHydrogenIonizingPhotonsMaximum
    self%exponent = exponent
    
    
    return
  end function powerLawCumulativeLuminosityConstructorInternal

  double precision function powerLawCumulativeLuminosity(self,rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum) result(luminosityFunctionIntegrated)
    !!{
    Returns the number of galaxies per QH based on the power law 
    \begin{equation}
    \Phi(Q_{max})-\Phi(Q_{min}) = (Q_{max}**(1-\alpha) - Q_{min}**(1-\alpha)) / (1-\alpha)},
    \end{equation}
    !!}

    implicit none
    class  (hiiRegionLuminosityFunctionPowerLaw), intent(inout) :: self
    double precision, intent(in   ) :: rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum
    
    ! Compute the unnormalized cumulative number of HII regions
    luminosityFunctionIntegrated  = (rateHydrogenIonizingPhotonsMaximum**(1.0d0-self%exponent)        &
       &     - rateHydrogenIonizingPhotonsMinimum**(-self%exponent+1.0d0))                            &
       &     / (1.0d0-self%exponent)
    
    return
  end function powerLawCumulativeLuminosity
