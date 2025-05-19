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
  Implementation of a peak-background split density 1D distribution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DVoight">
   <description>A 1D distribution function class for Voight profiles.</description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DVoight
     !!{
     Implementation of a voight 1D distribution function.
     !!}
     private
     logical          :: limitLowerExists, limitUpperExists
     double precision :: limitLower      , limitUpper      , &
          &              cdfAtLowerLimit , cdfAtUpperLimit , &
          &              gamma           , mu              , &
          &              sigma
 contains
     procedure :: density    => voightDensity
     procedure :: cumulative => voightCumulative
     procedure :: minimum    => voightMinimum
     procedure :: maximum    => voightMaximum
  end type distributionFunction1DVoight

  interface distributionFunction1DVoight
     !!{
     Constructors for the {\normalfont \ttfamily voight} 1D distribution function class.
     !!}
     module procedure voightConstructorParameters
     module procedure voightConstructorInternal
  end interface distributionFunction1DVoight

contains

  function voightConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily voight} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DVoight)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass  ), pointer       :: randomNumberGenerator_
    double precision                                              :: gamma                 , mu       , &
         &                                                           sigma                 ,limitLower, &
         &                                                           limitUpper

    !![
    <inputParameter>
      <name>gamma</name>
      <description>The parameter $\gamma$ of the Voight function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mu</name>
      <description>The parameter $\gamma$ of the Voight function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sigma</name>
      <description>The parameter $\sigma$ of the Voight function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>limitLower</name>
      <description>The lower limit of the Voight function.</description>
      <defaultValue>-huge(0.0d0)</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>limitUpper</name>
      <description>The upper limit of the Voight function.</description>
      <defaultValue>+huge(0.0d0)</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DVoight(gamma,mu,sigma,limitLower,limitUpper,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function voightConstructorParameters

  function voightConstructorInternal(gamma,mu,sigma,limitLower,limitUpper,randomNumberGenerator_) result(self)
    !!{
    Constructor for {\normalfont \ttfamily voight} 1D distribution function class.
    !!}
    type            (distributionFunction1DVoight)                                  :: self
    double precision                              , intent(in   )                   :: gamma                 , mu          , &
         &                                                                             sigma
    class           (randomNumberGeneratorClass  ), intent(in   ), optional, target :: randomNumberGenerator_
    double precision                              , intent(in   ), optional         :: limitLower            , limitUpper
    double precision                                                                :: cdfLower              , cdfUpper
    !![
    <constructorAssign variables="gamma, mu, sigma, *randomNumberGenerator_"/>
    !!]

    self    %limitLowerExists=.false.
    self    %limitUpperExists=.false.
    self    %cdfAtLowerLimit =0.0d0
    self    %cdfAtUpperLimit =1.0d0
    cdfLower                 =0.0d0
    cdfUpper                 =0.0d0
    if (present(limitLower)) cdfLower=self%cumulative(limitLower)
    if (present(limitUpper)) cdfUpper=self%cumulative(limitUpper)
    self%limitLowerExists=present(limitLower)
    self%limitUpperExists=present(limitUpper)
    if (self%limitLowerExists) then
       self%limitLower     =limitLower
       self%cdfAtLowerLimit=cdfLower
    end if
    if (self%limitUpperExists) then
       self%limitUpper     =limitUpper
       self%cdfAtUpperLimit=cdfUpper
    end if
    return
  end function voightConstructorInternal

  double precision function voightMaximum(self)
    !!{
    Return the maximum extent a Voight distribution.
    !!}
    implicit none
    class(distributionFunction1DVoight), intent(inout) :: self

    if (self%limitUpperExists) then
       voightMaximum=self%limitUpper
    else
       voightMaximum=+huge(0.0d0)
    end if
    return
  end function voightMaximum
  
  double precision function voightMinimum(self)
    !!{
    Return the minimum extent a Voight distribution.
    !!}
    implicit none
    class(distributionFunction1DVoight), intent(inout) :: self

    if (self%limitLowerExists) then
       voightMinimum=self%limitLower
    else
       voightMinimum=-huge(0.0d0)
    end if
    return
  end function voightMinimum
  
  double precision function voightDensity(self,x)
    !!{
    Return the density of a Voight distribution.
    !!}
    use :: Error_Functions         , only : Faddeeva
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DVoight), intent(inout) :: self
    double precision                              , intent(in   ) :: x
    double precision                                              :: x0
    double complex                                                :: z   , w

    if     (                                                   &
         &   (self%limitLowerExists .and. x < self%limitLower) &
         &  .or.                                               &
         &   (self%limitUpperExists .and. x > self%limitUpper) &
         & ) then
       voightDensity=0.0d0
    else
       ! Compute the value of x relative to the mean of the Gaussian component.
       x0=x-self%mu
       ! Evaluate the Faddeeva function at w(z).
       z =dcmplx(x0,self%gamma)/sqrt(2.0d0)/self%sigma
       w =Faddeeva(z)
       ! Compute the density.
       voightDensity=real(w)/sqrt(2.0d0*Pi)/self%sigma/(self%cdfAtUpperLimit-self%cdfAtLowerLimit)
    end if
    return
  end function voightDensity

  double precision function voightCumulative(self,x)
    !!{
    Return the cumulative probability of a Voight distribution.
    !!}
    use :: Error_Functions         , only : Error_Function
    use :: Hypergeometric_Functions, only : Hypergeometric_pFq
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DVoight), intent(inout)               :: self
    double precision                              , intent(in   )               :: x
    double complex                                , parameter    , dimension(2) :: a           =[dcmplx(1.0d0,0.0d0),dcmplx(1.0d0,0.0d0)]
    double complex                                , parameter    , dimension(2) :: b           =[dcmplx(1.5d0,0.0d0),dcmplx(2.0d0,0.0d0)]
    double precision                              , parameter                   :: widthMaximum=100.0d0
    double precision                                                            :: x0
    double complex                                                              :: z

    if      (self%limitLowerExists .and. x < self%limitLower) then
       voightCumulative=0.0d0
    else if (self%limitUpperExists .and. x > self%limitUpper) then
       voightCumulative=1.0d0
    else
       ! Compute the value of x relative to the mean of the Gaussian component.
       x0=x-self%mu
       ! Evaluate z.
       z =dcmplx(x0,self%gamma)/sqrt(2.0d0)/self%sigma
       ! Evaluate the cumulative distribution.
       voightCumulative=                             &
            & +(                                     &
            &   +real(                               &
            &         +0.5d0                         &
            &         +0.5d0                         &
            &         *Error_Function(z)             &
            &         +dcmplx(0.0d0,1.0d0)           &
            &         *z**2                          &
            &         /Pi                            &
            &         *Hypergeometric_pFq(a,b,-z**2) &
            &        )                               &
            &   -self%cdfAtLowerLimit                &
            &  )                                     &
            & /(                                     &
            &   +self%cdfAtUpperLimit                &
            &   -self%cdfAtLowerLimit                &
            &  )
    end if
    return
  end function voightCumulative
