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

!+    Contributions to this file made by: Sachi Weerasooriya

  !!{
  Implementation of a power law luminosity function for HII regions.
  !!}

  !![
  <hiiRegionLuminosityFunction name="hiiRegionLuminosityFunctionPowerLaw">
   <description>
    An HII region luminosity function class in which the luminosity function is given by:
    \begin{equation}
     \phi(Q_H) \propto \left\{ \begin{array}{ll}  Q_\mathrm{H}^{-\alpha} &amp; \hbox{ if } Q_\mathrm{H,min} &lt; Q_\mathrm{H} &lt; Q_\mathrm{H,max} \\ 0 &amp; \hbox{ otherwise} \end{array} \right. ,
    \end{equation}
    Where $Q_H$ is the rate of photon production rate, $Q_\mathrm{H,min}=${\normalfont \ttfamily
    [rateHydrogenIonizingPhotonsMinimum]} and $Q_\mathrm{H,max}=${\normalfont \ttfamily [rateHydrogenIonizingPhotonsMaximum]} and
    the minimum and maximum HII region luminosities respectively, and $\alpha=${\normalfont \ttfamily [exponent]}.
   </description>
  </hiiRegionLuminosityFunction>
  !!]
  type, extends(hiiRegionLuminosityFunctionClass) :: hiiRegionLuminosityFunctionPowerLaw
     !!{
     Implementation of a power law HII region luminosity function.
     !!}
     private
     double precision :: rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum, &
          &              exponent                          , normalization
     
   contains
     procedure :: cumulativeDistributionFunction => powerLawCumulativeDistributionFunction
     procedure :: cumulativeLuminosity           => powerLawCumulativeLuminosity
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
    Constructor for the \refClass{hiiRegionLuminosityFunctionPowerLaw} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionLuminosityFunctionPowerLaw)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum, &
         &                                                                  exponent
    !![
    <inputParameter>
      <name>exponent</name>
      <defaultValue>1.73d0</defaultValue>
      <defaultSource>\citep{santoro_phangs-muse_2022}</defaultSource>
      <description>Exponent of the differential luminosity function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>rateHydrogenIonizingPhotonsMinimum</name>
      <defaultValue>1.0d48</defaultValue>
      <defaultSource>(\citealt{santoro_phangs-muse_2022}; approximate)</defaultSource>
      <description>Minimum luminosity of HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>rateHydrogenIonizingPhotonsMaximum</name>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>Maximum luminosity of HII regions.</description>
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
    Internal constructor for the \refClass{hiiRegionLuminosityFunctionPowerLaw} luminosity function class.
    !!}
    
    implicit none
    type            (hiiRegionLuminosityFunctionPowerLaw)                :: self
    double precision                                     , intent(in   ) :: rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum, &
         &                                                                  exponent
    !![
    <constructorAssign variables=" rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum, exponent"/>
    !!]

    ! Normalize to unit probability.
    self%normalization=+(                                                                &
         &               +self%rateHydrogenIonizingPhotonsMaximum**(1.0d0-self%exponent) &
         &               -self%rateHydrogenIonizingPhotonsMinimum**(1.0d0-self%exponent) &
         &              )                                                                &
         &             /                                           (1.0d0-self%exponent)
    return
  end function powerLawCumulativeLuminosityConstructorInternal

  double precision function powerLawCumulativeDistributionFunction(self,rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum) result(distributionFunction)
    !!{
    Returns the fraction of HII regions in the given range of luminosity.
    !!}
    implicit none
    class           (hiiRegionLuminosityFunctionPowerLaw), intent(inout) :: self
    double precision                                     , intent(in   ) :: rateHydrogenIonizingPhotonsMinimum , rateHydrogenIonizingPhotonsMaximum
    double precision                                                     :: rateHydrogenIonizingPhotonsMinimum_, rateHydrogenIonizingPhotonsMaximum_
    
    rateHydrogenIonizingPhotonsMinimum_=max(rateHydrogenIonizingPhotonsMinimum,self%rateHydrogenIonizingPhotonsMinimum)
    rateHydrogenIonizingPhotonsMaximum_=min(rateHydrogenIonizingPhotonsMaximum,self%rateHydrogenIonizingPhotonsMaximum)
    if (rateHydrogenIonizingPhotonsMaximum_ > rateHydrogenIonizingPhotonsMinimum_) then
       distributionFunction=+(                                                            &
            &                 +rateHydrogenIonizingPhotonsMaximum_**(1.0d0-self%exponent) & 
            &                 -rateHydrogenIonizingPhotonsMinimum_**(1.0d0-self%exponent) &
            &                )                                                            &
            &               /                                       (1.0d0-self%exponent) &
            &               /self%normalization
    else
       distributionFunction=+0.0d0
    end if
    return
  end function powerLawCumulativeDistributionFunction

  double precision function powerLawCumulativeLuminosity(self,rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum) result(luminosity)
    !!{
    Returns the total luminosity, $Q_\mathrm{H}$ of HII regions in the given range of luminosity.
    !!}
    implicit none
    class           (hiiRegionLuminosityFunctionPowerLaw), intent(inout) :: self
    double precision                                     , intent(in   ) :: rateHydrogenIonizingPhotonsMinimum , rateHydrogenIonizingPhotonsMaximum
    double precision                                                     :: rateHydrogenIonizingPhotonsMinimum_, rateHydrogenIonizingPhotonsMaximum_
    
    rateHydrogenIonizingPhotonsMinimum_=max(rateHydrogenIonizingPhotonsMinimum,self%rateHydrogenIonizingPhotonsMinimum)
    rateHydrogenIonizingPhotonsMaximum_=min(rateHydrogenIonizingPhotonsMaximum,self%rateHydrogenIonizingPhotonsMaximum)
    if (rateHydrogenIonizingPhotonsMaximum_ > rateHydrogenIonizingPhotonsMinimum_) then
       luminosity=+(                                                            &
            &       +rateHydrogenIonizingPhotonsMaximum_**(2.0d0-self%exponent) & 
            &       -rateHydrogenIonizingPhotonsMinimum_**(2.0d0-self%exponent) &
            &      )                                                            &
            &     /                                       (2.0d0-self%exponent) &
            &     /self%normalization
    else
       luminosity=+0.0d0
    end if
    return
  end function powerLawCumulativeLuminosity
