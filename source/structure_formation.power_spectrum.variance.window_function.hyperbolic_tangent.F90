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

!+    Contributions to this file made by: Xiaolong Du

  !!{
  Implements a hyperbolic tangent power spectrum window function class.
  !!}
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionHyperbolicTangent">
   <description>
    A hyperbolic tangent window function for filtering of power spectra. The window function is given by:
    \begin{equation}
     W(k) = \tanh \left[ \frac{1}{(k R/c)^{\beta}} \right],
    \end{equation}
    where $R$ is related to $M$ via the standard relation $M = \frac{4\pi}{3}\bar\rho_m R^3$.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionHyperbolicTangent
     !!{
     A hyperbolic tangent power spectrum window function class.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: c                             , beta
   contains
     final     ::                      hyperbolicTangentDestructor
     procedure :: value             => hyperbolicTangentValue
     procedure :: wavenumberMaximum => hyperbolicTangentWavenumberMaximum
  end type powerSpectrumWindowFunctionHyperbolicTangent

  interface powerSpectrumWindowFunctionHyperbolicTangent
     !!{
     Constructors for the \refClass{powerSpectrumWindowFunctionHyperbolicTangent} power spectrum window function class.
     !!}
     module procedure hyperbolicTangentConstructorParameters
     module procedure hyperbolicTangentConstructorInternal
  end interface powerSpectrumWindowFunctionHyperbolicTangent

contains

  function hyperbolicTangentConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{powerSpectrumWindowFunctionHyperbolicTangent} power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionHyperbolicTangent)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (cosmologyParametersClass                    ), pointer       :: cosmologyParameters_
    double precision                                                              :: c                   , beta

    !![
    <inputParameter>
      <name>c</name>
      <source>parameters</source>
      <defaultValue>1.85d0</defaultValue>
      <description>The parameter $c$ in the hyperbolic tangent power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>4.41d0</defaultValue>
      <description>The parameter $\beta$ in the tangent power power spectrum window function.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=hyperbolicTangentConstructorInternal(cosmologyParameters_,c,beta)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function hyperbolicTangentConstructorParameters

  function hyperbolicTangentConstructorInternal(cosmologyParameters_,c,beta) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumWindowFunctionHyperbolicTangent} power spectrum window function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (powerSpectrumWindowFunctionHyperbolicTangent)                        :: self
    class           (cosmologyParametersClass                    ), target, intent(in   ) :: cosmologyParameters_
    double precision                                                      , intent(in   ) :: c                   , beta
    !![
    <constructorAssign variables="*cosmologyParameters_, c, beta"/>
    !!]

    return
  end function hyperbolicTangentConstructorInternal

  subroutine hyperbolicTangentDestructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumWindowFunctionHyperbolicTangent} power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionHyperbolicTangent), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine hyperbolicTangentDestructor

  double precision function hyperbolicTangentValue(self,wavenumber,smoothingMass,time)
    !!{
    Hyperbolic tangent window function used in computing the variance of the power spectrum.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumWindowFunctionHyperbolicTangent), intent(inout) :: self
    double precision                                              , intent(in   ) :: smoothingMass, wavenumber, &
         &                                                                           time
    double precision                                                              :: radius
    !$GLC attributes unused :: time

    radius =+(                                             &
         &    +3.0d0                                       &
         &    /4.0d0                                       &
         &    /Pi                                          &
         &    *smoothingMass                               &
         &    /self%cosmologyParameters_%OmegaMatter    () &
         &    /self%cosmologyParameters_%densityCritical() &
         &   )**(1.0d0/3.0d0)
    if (wavenumber <= 0.0d0) then
       hyperbolicTangentValue=+0.0d0
    else
       hyperbolicTangentValue=+tanh(              &
            &                       +1.0d0        &
            &                       /(            &
            &                         +wavenumber &
            &                         *radius     &
            &                         /self %c    &
            &                        )**self%beta &
            &                      )
    end if
    return
  end function hyperbolicTangentValue

  double precision function hyperbolicTangentWavenumberMaximum(self,smoothingMass)
    !!{
    Sets maximum wavenumber to effectively infinity (really large number).
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionHyperbolicTangent), intent(inout) :: self
    double precision                                              , intent(in   ) :: smoothingMass
    double precision                                              , parameter     :: wavenumberLarge=huge(1.0d0) ! Effective infinity.
    !$GLC attributes unused :: self, smoothingMass

    hyperbolicTangentWavenumberMaximum=wavenumberLarge
    return
  end function hyperbolicTangentWavenumberMaximum
