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
  Implementation of a power law mass function for HII regions.
  !!}

  !![
  <hiiRegionMassFunction name="hiiRegionMassFunctionPowerLaw">
   <description>
    An HII region luminosity function class in which the mass function is given by:
    \begin{equation}
     \phi(M_*) \propto \left\{ \begin{array}{ll}  Q_\mathrm{H}^{-\alpha} &amp; \hbox{ if } Q_\mathrm{H,min} &lt; Q_\mathrm{H} &lt; Q_\mathrm{H,max} \\ 0 &amp; \hbox{ otherwise} \end{array} \right. ,
    \end{equation}
    Where $M_*$ is the mass of gas in HII region, $M_\mathrm{min}=${\normalfont \ttfamily
    [massMinimum]} and $M_\mathrm{max}=${\normalfont \ttfamily [massMaximum]} and
    the minimum and maximum HII region masses respectively, and $\alpha=${\normalfont \ttfamily [exponent]}.
   </description>
  </hiiRegionLuminosityFunction>
  !!]
  type, extends(hiiRegionMassFunctionClass) :: hiiRegionMassFunctionPowerLaw
     !!{
     Implementation of a power law HII region mass function.
     !!}
     private
     double precision :: massMinimum, massMaximum, &
          &              exponent                          , normalization
     
   contains
     procedure :: cumulativeMasssFunction => powerLawCumulativeMassFunction
     procedure :: cumulativeMass           => powerLawCumulativeMass
  end type hiiRegionMassFunctionPowerLaw

  interface hiiRegionMassFunctionPowerLaw
     !!{
     Constructors for the hiiRegionMassFunctionPowerLaw for star HII region mass function class.
     !!}
     module procedure powerLawCumulativeMassConstructorParameters
     module procedure powerLawCumulativeMassConstructorInternal
  end interface hiiRegionMassFunctionPowerLaw

contains

  function powerLawCumulativeMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hiiRegionMassFunctionPowerLaw} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionMassFunctionPowerLaw)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: massMinimum, massMaximum, massCutoff, &
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
      <name>massMinimum</name>
      <defaultValue>1.0d48</defaultValue>
      <defaultSource>(\citealt{santoro_phangs-muse_2022}; approximate)</defaultSource>
      <description>Minimum mass of HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>Maximum mass of HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massCutOff</name>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>Cut off mass of HII regions.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionMassFunctionPowerLaw(massMinimum,massMaximum,massCutOff,exponent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawCumulativeMassConstructorParameters

  function powerLawCumulativeMassConstructorInternal(massMinimum,massMaximum,massCutOff,exponent) result(self)
    !!{
    Internal constructor for the \refClass{hiiRegionMassFunctionPowerLaw} luminosity function class.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    type            (hiiRegionMassFunctionPowerLaw)                :: self
    double precision                                     , intent(in   ) :: massMinimum, massMaximum, massCutOff, &
         &                                                                  exponent
    double precision                                                     :: massMinimumScaled, massMaximumScaled
    !![
    <constructorAssign variables=" massMinimum, massMaximum, massCutOff, exponent"/>
    !!]
    
    massMinimumScaled=self%massMinimum/self%massCutOff
    massMaximumScaled=self%massMaximum/self%massCutOff
    ! Normalize to unit probability.
    self%normalization=+(                                                                       &
    &                     +Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),massMinimumScaled) &
    &                     -Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),massMaximumScaled) &
    &                    )*self%massCutOff
    return
  end function powerLawCumulativeMassConstructorInternal

  double precision function powerLawCumulativeMassFunction(self,massMinimum,massMaximum) result(massFunction)
    !!{
    Returns the fraction of HII regions in the given range of masses.
    !!}
    implicit none
    class           (hiiRegionMassFunctionPowerLaw), intent(inout) :: self
    double precision                                     , intent(in   ) :: massMinimum , massMaximum
    double precision                                                     :: massMinimumScaled_, massMaximumScaled_
    
    massMinimumScaled_=max(massMinimum/self%massCutOff,self%massMinimum/self%massCutOff)
    massMaximumScaled_=min(massMaximum/self%massCutOff,self%massMaximum/self%massCutOff)
    if (massMaximumScaled_ > massMinimumScaled_) then
       massFunction=+(                                                                                    &
            &                +(massMaximumScaled_**self%exponent)*exp(-massMaximumScaled_)   &
            &               -(massMinimumScaled_**self%exponent)*exp(-massMinimumScaled_)     &
            &                )                                                                &
            &               /self%normalization
    else
       massFunction=+0.0d0
    end if
    return
  end function powerLawCumulativeMassFunction

  double precision function powerLawCumulativeMass(self,massMinimum,massMaximum) result(massHIIRegion)
    !!{
    Returns the total mass, $M_*$ of HII regions in the given range of mass.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    class           (hiiRegionMassFunctionPowerLaw), intent(inout) :: self
    double precision                                     , intent(in   ) :: massMinimum , massMaximum
    double precision                                                     :: massMinimumScaled_, massMaximumScaled_
    
    massMinimumScaled_=max(massMinimum/self%massCutOff,self%massMinimum/self%massCutOff)
    massMaximumScaled_=min(massMaximum/self%massCutOff,self%massMaximum/self%massCutOff)
    if (massMaximumScaled_ > massMinimumScaled_) then
       massHIIRegion=+(                                                                                    &
            &               +Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),massMinimumScaled_)     &
            &               -Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),massMaximumScaled_)     &
            &                )                                                                  &
            &               /self%normalization
    else
       massHIIRegion=+0.0d0
    end if
    return
  end function powerLawCumulativeMass
