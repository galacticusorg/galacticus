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

  !!{
  Implementation of a 1D gamma distribution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DGamma">
   <description>A 1D gamma distribution function class.</description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DGamma
     !!{
     Implementation of a 1D gamma distribution function.
     !!}
     private
     logical          :: limitLowerExists, limitUpperExists
     double precision :: limitLower      , limitUpper      , &
          &              cdfAtLowerLimit , cdfAtUpperLimit , &
          &              shape           , rate
   contains
     procedure :: density    => gammaDensity
     procedure :: cumulative => gammaCumulative
     procedure :: inverse    => gammaInverse
  end type distributionFunction1DGamma

  interface distributionFunction1DGamma
     !!{
     Constructors for the {\normalfont \ttfamily gamma} 1D distribution function class.
     !!}
     module procedure gammaConstructorParameters
     module procedure gammaConstructorInternal
  end interface distributionFunction1DGamma

contains

  function gammaConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily gamma} 1D distribution function class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DGamma)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass ), pointer       :: randomNumberGenerator_
    double precision                                             :: shape                 , rate      , &
         &                                                          limitLower            , limitUpper

    !![
    <inputParameter>
      <name>shape</name>
      <description>The shape parameter of the gamma distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>rate</name>
      <description>The rate parameter of the gamma distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>limitLower</name>
      <description>The lower limit of the gamma distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>limitUpper</name>
      <description>The upper limit of the gamma distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DGamma(shape,rate,randomNumberGenerator_,limitLower,limitUpper)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gammaConstructorParameters

  function gammaConstructorInternal(shape,rate,randomNumberGenerator_,limitLower,limitUpper) result(self)
    !!{
    Constructor for {\normalfont \ttfamily gamma} 1D distribution function class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (distributionFunction1DGamma)                                  :: self
    double precision                             , intent(in   )                   :: shape                 , rate
    double precision                             , intent(in   ), optional         :: limitLower            , limitUpper
    class           (randomNumberGeneratorClass ), intent(in   ), optional, target :: randomNumberGenerator_
    !![
    <constructorAssign variables="shape, rate, *randomNumberGenerator_"/>
    !!]

    if (rate <= 0.0d0 .or. shape <= 0.0d0) call Error_Report('rate>0 and shape>0 are required'//{introspection:location})
    self%limitLowerExists=.false.
    self%limitUpperExists=.false.
    self%cdfAtLowerLimit =0.0d0
    self%cdfAtUpperLimit =1.0d0
    if (present(limitLower)) then
       if (limitLower < 0.0d0) call Error_Report('limitLower≥0 is required'//{introspection:location})
       self%limitLower     =limitLower
       self%cdfAtLowerLimit=self%cumulative(limitLower)
    else
       self%cdfAtLowerLimit=0.0d0
    end if
    if (present(limitUpper)) then
       if (limitUpper < 0.0d0) call Error_Report('limitUpper≥0 is required'//{introspection:location})
       self%limitUpper     =limitUpper
       self%cdfAtUpperLimit=self%cumulative(limitUpper)
    else
       self%cdfAtUpperLimit=1.0d0
    end if
    self%limitLowerExists=present(limitLower)
    self%limitUpperExists=present(limitUpper)
    return
  end function gammaConstructorInternal

  double precision function gammaDensity(self,x)
    !!{
    Return the density of a Gamma distribution.
    !!}
    use :: Gamma_Functions, only : Gamma_Function
    implicit none
    class           (distributionFunction1DGamma), intent(inout) :: self
    double precision                             , intent(in   ) :: x

    if     (                                                   &
         &                                x < 0.0d0            &
         &  .or.                                               &
         &   (self%limitLowerExists .and. x < self%limitLower) &
         &  .or.                                               &
         &   (self%limitUpperExists .and. x > self%limitUpper) &
         & ) then
       gammaDensity=0.0d0
    else
       gammaDensity=+x**(self%shape-1.0d0)      &
            &       *exp(-x*self%rate)          &
            &       *self%rate**self%shape      &
            &       /Gamma_Function(self%shape) &
            &       /(                          &
            &         +self%cdfAtUpperLimit     &
            &         -self%cdfAtLowerLimit     &
            &        )
    end if
    return
  end function gammaDensity

  double precision function gammaCumulative(self,x)
    !!{
    Return the cumulative probability of a Gamma distribution.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete_Complementary
    implicit none
    class           (distributionFunction1DGamma), intent(inout) :: self
    double precision                             , intent(in   ) :: x

    if      (                            x <= 0.0d0          ) then
       gammaCumulative=0.0d0
    else if (self%limitLowerExists .and. x <  self%limitLower) then
       gammaCumulative=0.0d0
    else if (self%limitUpperExists .and. x >  self%limitUpper) then
       gammaCumulative=1.0d0
    else
       gammaCumulative=min(                                                            &
            &              1.0d0                                                     , &
            &              max(                                                        &
            &                   0.0d0                                                , &
            &                   (                                                      &
            &                    +Gamma_Function_Incomplete_Complementary(             &
            &                                                              self%shape, &
            &                                                              self%rate   &
            &                                                             *x           &
            &                                                            )             &
            &                    -self%cdfAtLowerLimit                                 &
            &                   )                                                      &
            &                  /(                                                      &
            &                    +self%cdfAtUpperLimit                                 &
            &                    -self%cdfAtLowerLimit                                 &
            &                   )                                                      &
            &                 )                                                        &
            &             )
    end if
    return
  end function gammaCumulative

  double precision function gammaInverse(self,p)
    !!{
    Return the inverse of a Gamma distribution.
    !!}
    use :: Error          , only : Error_Report
    use :: Gamma_Functions, only : Inverse_Gamma_Function_Incomplete_Complementary
    implicit none
    class           (distributionFunction1DGamma), intent(inout), target :: self
    double precision                             , intent(in   )         :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    gammaInverse=+Inverse_Gamma_Function_Incomplete_Complementary(                         &
         &                                                           self%shape          , &
         &                                                        +p                       &
         &                                                        *(                       &
         &                                                          +self%cdfAtUpperLimit  &
         &                                                          -self%cdfAtLowerLimit  &
         &                                                         )                       &
         &                                                        +  self%cdfAtLowerLimit  &
         &                                                       )                         &
         &       /                                                   self%rate
    return
  end function gammaInverse
