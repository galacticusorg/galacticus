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

!+    Contributions to this file made by: Sachi Weerasooriya, Andrew Benson

  !!{
  Implementation of a power law mass function for HII regions.
  !!}

  !![
  <hiiRegionMassFunction name="hiiRegionMassFunctionRosolowsky2021">
   <description>
    An HII region luminosity function class in which the GMC mass function is given by Rosolowsky et al. 2021:
    \begin{equation}
      \left(\frac{m}{M}\right)^{\alpha} \exp\!\left(-\frac{m}{M}\right)
    \end{equation}
    Where $M$ is the mass of GMC in, $M_\mathrm{min}=${\normalfont \ttfamily
    [massMinimum]} and $M_\mathrm{max}=${\normalfont \ttfamily [massMaximum]} and
    the minimum and maximum GMC masses respectively, and $\alpha=${\normalfont \ttfamily [exponent]}.
   </description>
  </hiiRegionLuminosityFunction>
  !!]
  type, extends(hiiRegionMassFunctionClass) :: hiiRegionMassFunctionRosolowsky2021
     !!{
     Implementation of a power law HII region mass function.
     !!}
     private
     double precision :: massMinimum, massMaximum, &
          &              exponent                          , normalization
     
   contains
     procedure :: cumulativeMasssFunction => powerLawCumulativeMassFunction
     procedure :: cumulativeMass           => powerLawCumulativeMass
  end type hiiRegionMassFunctionRosolowsky2021

  interface hiiRegionMassFunctionRosolowsky2021
     !!{
     Constructors for the hiiRegionMassFunctionRosolowsky2021 for star HII region mass function class.
     !!}
     module procedure powerLawCumulativeMassConstructorParameters
     module procedure powerLawCumulativeMassConstructorInternal
  end interface hiiRegionMassFunctionRosolowsky2021

contains

  function powerLawCumulativeMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hiiRegionMassFunctionRosolowsky2021} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionMassFunctionRosolowsky2021)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: massMinimum, massMaximum, massCutoff, &
         &                                                                  exponent
    !![
    <inputParameter>
      <name>exponent</name>
      <defaultValue>-1.4d0</defaultValue>
      <defaultSource>\citep{rosolowsky_2021}</defaultSource>
      <description>Exponent of the differential luminosity function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>epsilon</name>
      <defaultValue>1.0d-2</defaultValue>
      <defaultSource>\citep{rosolowsky_2021}</defaultSource>
      <description>Exponent of the differential luminosity function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d5</defaultValue>
      <defaultSource>(\citealt{rosolowsky_2021_2022}; approximate)</defaultSource>
      <description>Minimum mass of GMC.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d8</defaultValue>
      <description>Maximum mass of GMC.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massCutOff</name>
      <defaultValue>2.6d6</defaultValue>
      <description>Cut off mass of HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionMassFunctionRosolowsky2021(massMinimum,massMaximum,massCutOff,exponent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawCumulativeMassConstructorParameters

  function powerLawCumulativeMassConstructorInternal(massMinimum,massMaximum,massCutOff,exponent) result(self)
    !!{
    Internal constructor for the \refClass{hiiRegionMassFunctionRosolowsky2021} mass function class.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    type            (hiiRegionMassFunctionRosolowsky2021)                :: self
    double precision                                     , intent(in   ) :: massMinimum, massMaximum, massCutOff, &
         &                                                                  exponent, massStellarMinimum,massStellarMaximum
    double precision                                                     :: xMinimum, xMaximum
    !![
    <constructorAssign variables=" massMinimum, massMaximum, massCutOff, exponent"/>
    !!]
    
    self%massStellarMinimum=massMinimum*epsilon
    self%massStellarMaximum=massMaximum*epsilon
    self%massStellarCutOff =massCutOff *epsilon

    xMinimum=self%massStellarMinimum/self%massStellarCutOff
    xMaximum=self%massStellarMaximum/self%massStellarCutOff

    ! Normalize to unit probability.
    self%normalization=+(                                                                       &
    &                     +Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),xMinimum) &
    &                     -Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),xMaximum) &
    &                    )
    return
  end function powerLawCumulativeMassConstructorInternal

  double precision function powerLawCumulativeMassFunction(self,massMinimum,massMaximum) result(massFunction)
    !!{
    Returns the fraction of HII regions in the given range of masses.
    !!}
    implicit none
    class           (hiiRegionMassFunctionRosolowsky2021), intent(inout) :: self
    double precision                                     , intent(in   ) :: massMinimum , massMaximum
    double precision                                                     :: xMinimum, xMaximum
    
    xMinimum=max(massStellarMinimum/self%massStellarCutOff,self%massStellarMinimum/self%massStellarCutOff)
    xMaximum=min(massStellarMaximum/self%massStellarCutOff,self%massStellarMaximum/self%massStellarCutOff)
    if (xMaximum > xMinimum_) then
       massFunction=+(                                                                                    &
            &               +Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),xMaximum)     &
            &               -Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),xMinimum)     &
            &                )                                                                  &
            &               /self%normalization
    else
       massFunction=+0.0d0
    end if
    return
  end function powerLawCumulativeMassFunction

  double precision function powerLawCumulativeMass(self,massMinimum,massMaximum) result(massHIIRegion)
    !!{
    Returns the total mass, $M_*$ of HII region in the given range of mass.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    class           (hiiRegionMassFunctionRosolowsky2021), intent(inout) :: self
    double precision                                     , intent(in   ) :: massMinimum , massMaximum
    double precision                                                     :: xMinimum, xMaximum
    
    xMinimum=max(massMinimum/self%massCutOff,self%massMinimum/self%massCutOff)
    xMaximum=min(massMaximum/self%massCutOff,self%massMaximum/self%massCutOff)
    if (xMaximum > xMinimum) then
       massHIIRegion =      (                                                                         &
            &               +Gamma_Function_Incomplete_Unnormalized((self%exponent+2.0d0),xMaximum)     &
            &               -Gamma_Function_Incomplete_Unnormalized((self%exponent+2.0d0),xMinimum)     &
            &               )                                                                  &
            &               /self%normalization
    else
       massHIIRegion=+0.0d0
    end if
    return
  end function powerLawCumulativeMass
