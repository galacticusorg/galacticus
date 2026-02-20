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

!!k
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!!{
Implements a \cite{giocoli_population_2008} unevolved dark matter subhalo mass function class.
!!}

  !![
  <unevolvedSubhaloMassFunction name="unevolvedSubhaloMassFunctionGiocoli2008">
   <description>The halo mass function is computed from the function given by \cite{giocoli_population_2008}.</description>
  </unevolvedSubhaloMassFunction>
  !!]
  type, extends(unevolvedSubhaloMassFunctionClass) :: unevolvedSubhaloMassFunctionGiocoli2008
     !!{
     An unevolved subhalo mass function class using the model of \cite{giocoli_population_2008}.
     !!}
     private
     double precision :: normalization, exponent
    contains
     procedure :: differential => giocoli2008Differential
     procedure :: integrated   => giocoli2008Integrated
  end type unevolvedSubhaloMassFunctionGiocoli2008

  interface unevolvedSubhaloMassFunctionGiocoli2008
     !!{
     Constructors for the \refClass{unevolvedSubhaloMassFunctionGiocoli2008} halo mass function class.
     !!}
     module procedure giocoli2008ConstructorParameters
     module procedure giocoli2008ConstructorInternal
  end interface unevolvedSubhaloMassFunctionGiocoli2008

  ! Coefficient in the exponential term.
  double precision, parameter :: coefficientExponential=6.283d0

contains

  function giocoli2008ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{unevolvedSubhaloMassFunctionGiocoli2008} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (unevolvedSubhaloMassFunctionGiocoli2008)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    double precision                                                         :: normalization, exponent

    !![
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <defaultValue>0.21d0</defaultValue>
      <defaultSource>\cite{giocoli_population_2008}</defaultSource>
      <description>The parameter $N_0$ in the \cite{giocoli_population_2008} unevolved subhalo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <source>parameters</source>
      <defaultValue>0.8d0</defaultValue>
      <defaultSource>\cite{giocoli_population_2008}</defaultSource>
      <description>The parameter $\alpha$ in the \cite{giocoli_population_2008} unevolved subhalo mass function fit.</description>
    </inputParameter>
    !!]
    self=unevolvedSubhaloMassFunctionGiocoli2008(normalization,exponent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
   return
  end function giocoli2008ConstructorParameters

  function giocoli2008ConstructorInternal(normalization,exponent) result(self)
    !!{
    Internal constructor for the \refClass{unevolvedSubhaloMassFunctionGiocoli2008} halo mass function class.
    !!}
    implicit none
    type            (unevolvedSubhaloMassFunctionGiocoli2008)                :: self
    double precision                                         , intent(in   ) :: normalization, exponent
    !![
    <constructorAssign variables="normalization, exponent"/>
    !!]

    return
  end function giocoli2008ConstructorInternal

  double precision function giocoli2008Differential(self,time,mass,massHost)
    !!{
    Return the differential unevolved subhalo mass function at the given time and mass.
    !!}
    implicit none
    class           (unevolvedSubhaloMassFunctionGiocoli2008), intent(inout) :: self
    double precision                                         , intent(in   ) :: time    , mass, &
         &                                                                      massHost
    double precision                                                         :: x
    !$GLC attributes unused :: time

    x                      =+mass                        &
         &                  /massHost                    &
         &                  /self%exponent
    giocoli2008Differential=+self%normalization          &
         &                  /self%exponent               &
         &                  /massHost                    &
         &                  /x**(1.0d0+self%exponent)    &
         &                  /exp(                        &
         &                       -coefficientExponential &
         &                       *x**3                   &
         &                      )
    return
  end function giocoli2008Differential

  double precision function giocoli2008Integrated(self,time,massLow,massHigh,massHost)
    !!{
    Return the integrated unevolved subhalo mass function at the given time and mass.
    !!}
    implicit none
    class           (unevolvedSubhaloMassFunctionGiocoli2008), intent(inout) :: self
    double precision                                         , intent(in   ) :: time    , massLow , &
         &                                                                      massHigh, massHost
    double precision                                                         :: xLow    , xHigh
    !$GLC attributes unused :: time

    xLow                 =+massLow                                                       &
         &                /massHost                                                      &
         &                /self%exponent
    xHigh                =+massHigh                                                      &
         &                /massHost                                                      &
         &                /self%exponent
    giocoli2008Integrated=+self%normalization                                            &
         &                /3.0d0                                                         &
         &                *(                                                             &
         &                  +(coefficientExponential*xLow **3)**(self%exponent/3.0d0)    &
         &                  /                        xLow     ** self%exponent           &
         &                  *gammaIncomplete(        xLow )                              &
         &                  -(coefficientExponential*xHigh**3)**(self%exponent/3.0d0)    &
         &                  /                        xHigh    ** self%exponent           &
         &                  *gammaIncomplete(        xHigh)                              &
         &                 )
    return

  contains

    double precision function gammaIncomplete(x)
      !!{
      Evaluate the incomplete gamma function, possibly for a negative exponent.
      !!}
      use :: Error          , only : Error_Report
      use :: Gamma_Functions, only : Gamma_Function, Gamma_Function_Incomplete
      implicit none
      double precision, intent(in   ) :: x

      if      (self%exponent < 0.0d0) then
         gammaIncomplete=+Gamma_Function_Incomplete(                                               &
              &                                           -self%exponent/3.0d0                   , &
              &                                     +coefficientExponential*x**3                   &
              &                                    )                                               &
              &          *Gamma_Function           (                                               &
              &                                           -self%exponent/3.0d0                     &
              &                                    )
      else if (self%exponent < 3.0d0) then
         gammaIncomplete=+Gamma_Function_Incomplete(                                               &
              &                                     +1.0d0-self%exponent/3.0d0                   , &
              &                                     +coefficientExponential*x**3                   &
              &                                    )                                               &
              &          *Gamma_Function           (                                               &
              &                                     +1.0d0-self%exponent/3.0d0                     &
              &                                    )                                               &
              &          *exp                      (                                               &
              &                                     -coefficientExponential*x**3                   &
              &                                    )                                               &
              &          /                           coefficientExponential**(self%exponent/3.0d0) &
              &          /                           x                     ** self%exponent
      else
         gammaIncomplete=0.0d0
         call Error_Report('exponent out of range'//{introspection:location})
      end if
      return
    end function gammaIncomplete

  end function giocoli2008Integrated

