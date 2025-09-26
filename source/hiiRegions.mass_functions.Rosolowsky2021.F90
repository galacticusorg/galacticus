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
  Implementation of a mass function for HII regions following the model of \cite{rosolowsky_giant_2021}.
  !!}

  !![
  <hiiRegionMassFunction name="hiiRegionMassFunctionRosolowsky2021">
   <description>
     An HII region stellar mass function class in which giant molecular clouds are assumed to turn a fixed fraction
     $\epsilon=${\normalfont \ttfamily [epsilon]} of their mass into stars, and the giant molecular cloud mass function is given
     by \cite{rosolowsky_giant_2021}:
     \begin{equation}
     \phi(M) \propto \left(\frac{m}{M}\right)^{\alpha} \exp\!\left(-\frac{m}{M}\right),
     \end{equation}
     where $M$ is the mass of giant molecular cloud, $M_\mathrm{min}=${\normalfont \ttfamily [massMinimum]} and
     $M_\mathrm{max}=${\normalfont \ttfamily [massMaximum]} and the minimum and maximum giant molecular cloud masses respectively,
     and $\alpha=${\normalfont \ttfamily [exponent]}.
   </description>
  </hiiRegionMassFunction>
  !!]
  type, extends(hiiRegionMassFunctionClass) :: hiiRegionMassFunctionRosolowsky2021
     !!{
     Implementation of a mass function for HII regions following the model of \cite{rosolowsky_giant_2021}.
     !!}
     private
     double precision :: massMinimum       , massMaximum       , &
          &              massStellarMinimum, massStellarMaximum, &
          &              massCutOff        , massStellarCutOff , &
          &              epsilon           , exponent          , &
          &              normalization
     
   contains
     procedure :: cumulativeDistributionFunction => rosolowsky2021CumulativeDistributionFunction
     procedure :: cumulativeMass                 => rosolowsky2021CumulativeMass
  end type hiiRegionMassFunctionRosolowsky2021

  interface hiiRegionMassFunctionRosolowsky2021
     !!{
     Constructors for the \refClass{hiiRegionMassFunctionRosolowsky2021} HII region stellar mass function class.
     !!}
     module procedure rosolowsky2021CumulativeMassConstructorParameters
     module procedure rosolowsky2021CumulativeMassConstructorInternal
  end interface hiiRegionMassFunctionRosolowsky2021

contains

  function rosolowsky2021CumulativeMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hiiRegionMassFunctionRosolowsky2021} timescale for star formation class which takes a parameter
    set as input.    
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (hiiRegionMassFunctionRosolowsky2021)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: massMinimum, massMaximum, &
         &                                                                  massCutoff , exponent   , &
         &                                                                  epsilon
    !![
    <inputParameter>
      <name>epsilon</name>
      <defaultValue>3.0d-2</defaultValue>
      <defaultSource>\citep[][figure 2]{grudic_nature_2019}</defaultSource>
      <description>Fraction of giant molecular cloud mass turned into stars.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <defaultValue>-1.2d0</defaultValue>
      <defaultSource>\citep[][table 5; ``disk'' sample]{rosolowsky_2021}</defaultSource>
      <description>Exponent of the differential mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d4</defaultValue>
      <defaultSource>\citep{fukui_molecular_2010}</defaultSource>
      <description>Minimum mass of giant molecular clouds.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d8</defaultValue>
      <description>Maximum mass ofgiant molecular clouds.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massCutOff</name>
      <defaultValue>4.7d6</defaultValue>
      <defaultSource>\citep[][table 5; ``disk'' sample]{rosolowsky_2021}</defaultSource>
      <description>Cut off mass of giant molecular clouds.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionMassFunctionRosolowsky2021(massMinimum,massMaximum,massCutOff,exponent,epsilon)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function rosolowsky2021CumulativeMassConstructorParameters

  function rosolowsky2021CumulativeMassConstructorInternal(massMinimum,massMaximum,massCutOff,exponent,epsilon) result(self)
    !!{
    Internal constructor for the \refClass{hiiRegionMassFunctionRosolowsky2021} HII region stellar mass function class.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    type            (hiiRegionMassFunctionRosolowsky2021)                :: self
    double precision                                     , intent(in   ) :: massMinimum, massMaximum, &
         &                                                                  massCutOff , exponent   , &
         &                                                                  epsilon
    double precision                                                     :: xMinimum   , xMaximum
    !![
    <constructorAssign variables="massMinimum, massMaximum, massCutOff, exponent, epsilon"/>
    !!]

    ! Convert from GMC masses to stellar masses.
    self%massStellarMinimum=massMinimum*epsilon
    self%massStellarMaximum=massMaximum*epsilon
    self%massStellarCutOff =massCutOff *epsilon
    ! Normalize to unit probability.
    xMinimum          =+self%massStellarMinimum/self%massStellarCutOff
    xMaximum          =+self%massStellarMaximum/self%massStellarCutOff
    self%normalization=+Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),xMinimum) &
         &             -Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),xMaximum)
    return
  end function rosolowsky2021CumulativeMassConstructorInternal

  double precision function rosolowsky2021CumulativeDistributionFunction(self,massMinimum,massMaximum) result(massFunction)
    !!{
    Returns the fraction of HII regions in the given range of stellar mass.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    class           (hiiRegionMassFunctionRosolowsky2021), intent(inout) :: self
    double precision                                     , intent(in   ) :: massMinimum, massMaximum
    double precision                                                     :: xMinimum   , xMaximum

    ! Compute limiting, scaled masses.
    xMinimum=max(massMinimum/self%massStellarCutOff,self%massStellarMinimum/self%massStellarCutOff)
    xMaximum=min(massMaximum/self%massStellarCutOff,self%massStellarMaximum/self%massStellarCutOff)
    ! Evaluate the mass function.
    if (xMaximum > xMinimum) then
       massFunction=+(                                                                        &
            &         +Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),xMinimum) &
            &         -Gamma_Function_Incomplete_Unnormalized((self%exponent+1.0d0),xMaximum) &
            &        )                                                                        &
            &       /self%normalization
    else
       massFunction=+0.0d0
    end if
    return
  end function rosolowsky2021CumulativeDistributionFunction
  
  double precision function rosolowsky2021CumulativeMass(self,massMinimum,massMaximum) result(massHIIRegion)
    !!{
    Returns the total mass, $M_\star$, of HII region in the given range of mass.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Unnormalized
    implicit none
    class           (hiiRegionMassFunctionRosolowsky2021), intent(inout) :: self
    double precision                                     , intent(in   ) :: massMinimum, massMaximum
    double precision                                                     :: xMinimum   , xMaximum
    
    ! Compute limiting, scaled masses.
    xMinimum=max(massMinimum/self%massStellarCutOff,self%massStellarMinimum/self%massStellarCutOff)
    xMaximum=min(massMaximum/self%massStellarCutOff,self%massStellarMaximum/self%massStellarCutOff)
    ! Evaluate the mass.
    if (xMaximum > xMinimum) then
       massHIIRegion=+self%massStellarCutOff                                                   &
            &        *(                                                                        &
            &          +Gamma_Function_Incomplete_Unnormalized((self%exponent+2.0d0),xMinimum) &
            &          -Gamma_Function_Incomplete_Unnormalized((self%exponent+2.0d0),xMaximum) &
            &         )                                                                        &
            &        /self%normalization
    else
       massHIIRegion=+0.0d0
    end if
    return
  end function rosolowsky2021CumulativeMass
