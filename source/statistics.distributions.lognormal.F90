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
  Implementation of a log-normal 1D distribution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DLogNormal">
   <description>A normal 1D distribution function class.</description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DNormal) :: distributionFunction1DLogNormal
     !!{
     Implementation of a log-normal 1D distribution function:
     \begin{equation}
       p(x) = \frac{1}{x \sigma \sqrt{2\pi}} \exp\left( - \frac{[\log x - \log x_0]^2}{2\sigma^2} \right).
     \end{equation}
     The parameters of this distribution can be specified either directly via $x_0$ and $\sigma$, or indirectly via the mean and
     variance $\exp(\log x_0 + \sigma^2/2)$ and $[\exp(\sigma^2)-1] \exp(2\log x_0 + \sigma^2)$.     
     !!}
     private
     double precision :: x0, sigma
   contains
     procedure :: density    => logNormalDensity
     procedure :: cumulative => logNormalCumulative
     procedure :: inverse    => logNormalInverse
     procedure :: minimum    => logNormalMinimum
     procedure :: maximum    => logNormalMaximum
  end type distributionFunction1DLogNormal

  interface distributionFunction1DLogNormal
     !!{
     Constructors for the \refClass{distributionFunction1DLogNormal} 1D distribution function class.
     !!}
     module procedure logNormalConstructorParameters
     module procedure logNormalConstructorInternal
  end interface distributionFunction1DLogNormal

contains

  function logNormalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DLogNormal} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DLogNormal)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass     ), pointer       :: randomNumberGenerator_
    double precision                                                 :: mean                  , variance  , &
         &                                                              x0                    , sigma     , &
         &                                                              limitLower            , limitUpper

    if (parameters%isPresent('mean')) then
       if (.not.parameters%isPresent('variance')) call Error_Report('must specify both [mean] and [variance], or [x0] and [sigma]'//{introspection:location})
       if (parameters%isPresent('x0').or.parameters%isPresent('sigma')) call Error_Report('must specify either [mean] and [variance], or [x0] and [sigma], not both'//{introspection:location})
       !![
       <inputParameter>
	 <name>mean</name>
	 <description>The mean of the log-normal distribution (\emph{ignoring} any imposed upper and lower limits).</description>
	 <source>parameters</source>
       </inputParameter>
       <inputParameter>
	 <name>variance</name>
	 <description>The variance of the log-normal distribution (\emph{ignoring} any imposed upper and lower limits).</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    else if (parameters%isPresent('x0')) then
       if (.not.parameters%isPresent('sigma')) call Error_Report('must specify both [mean] and [variance], or [x0] and [sigma]'//{introspection:location})
       if (parameters%isPresent('mean').or.parameters%isPresent('variance')) call Error_Report('must specify either [mean] and [variance], or [x0] and [sigma], not both'//{introspection:location})
       !![
       <inputParameter>
	 <name>x0</name>
	 <description>The parameter $x_0$ of the log-normal distribution.</description>
	 <source>parameters</source>
       </inputParameter>
       <inputParameter>
	 <name>sigma</name>
	 <description>The parameter $\sigma$ of the log-normal distribution.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    else
       call Error_Report('must specify both [mean] and [variance], or [x0] and [sigma]'//{introspection:location})
    end if
    !![
    <inputParameter>
      <name>limitLower</name>
      <description>The lower limit of the normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>limitUpper</name>
      <description>The upper limit of the normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    <conditionalCall>
      <call>
	self=distributionFunction1DLogNormal(limitLower=limitLower,limitUpper=limitUpper,randomNumberGenerator_=randomNumberGenerator_{conditions})
      </call>
      <argument name="mean"     value="mean"     parameterPresent="parameters"/>
      <argument name="variance" value="variance" parameterPresent="parameters"/>
      <argument name="x0"       value="x0"       parameterPresent="parameters"/>
      <argument name="sigma"    value="sigma"    parameterPresent="parameters"/>
    </conditionalCall>
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function logNormalConstructorParameters

  function logNormalConstructorInternal(mean,variance,x0,sigma,limitLower,limitUpper,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DLogNormal} 1D distribution function class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (distributionFunction1DLogNormal)                                  :: self
    double precision                                 , intent(in   ), optional         :: mean                  , variance      , &
         &                                                                                x0                    , sigma         , &
         &                                                                                limitLower            , limitUpper
    class           (randomNumberGeneratorClass     ), intent(in   ), optional, target :: randomNumberGenerator_
    double precision                                                                   :: meanNormal            , varianceNormal

    if (present(mean)) then
       if (.not.present(variance)) call Error_Report('must specify both `mean` and `variance`, or `x0` and `sigma`'//{introspection:location})
       if (present(x0).or.present(sigma)) call Error_Report('must specify either `mean` and `variance`, or `x0` and `sigma`, not both'//{introspection:location})
       varianceNormal=+log(               &
            &              +(             &
            &                +variance    &
            &                /mean    **2 &
            &               )             &
            &              +1.0d0         &
            &             )
       meanNormal    =+log(               &
            &              +mean          &
            &             )               &
            &         -    varianceNormal &
            &         /2.0d0
    else if (present(x0)) then
       if (.not.present(sigma)) call Error_Report('must specify both `mean` and `variance`, or `x0` and `sigma`'//{introspection:location})
       if (present(mean).or.present(variance)) call Error_Report('must specify either `mean` and `variance`, or `x0` and `sigma`, not both'//{introspection:location})
       varianceNormal=sigma**2
       meanNormal    =log(x0)
    else
       call Error_Report('must specify both `mean` and `variance`, or `x0` and `sigma]'//{introspection:location})
    end if
    if (present(limitLower)) then
       if (present(limitUpper)) then
          self%distributionFunction1DNormal=distributionFunction1DNormal(meanNormal,varianceNormal,limitLower=log(limitLower),limitUpper=log(limitUpper),randomNumberGenerator_=randomNumberGenerator_)
       else
          self%distributionFunction1DNormal=distributionFunction1DNormal(meanNormal,varianceNormal,limitLower=log(limitLower)                           ,randomNumberGenerator_=randomNumberGenerator_)
       end if
    else
       if (present(limitUpper)) then
          self%distributionFunction1DNormal=distributionFunction1DNormal(meanNormal,varianceNormal                           ,limitUpper=log(limitUpper),randomNumberGenerator_=randomNumberGenerator_)
       else
          self%distributionFunction1DNormal=distributionFunction1DNormal(meanNormal,varianceNormal                                                      ,randomNumberGenerator_=randomNumberGenerator_)
       end if
    end if
    return
  end function logNormalConstructorInternal

  double precision function logNormalMinimum(self)
    !!{
    Return the minimum possible value of a uniform distribution.
    !!}
    implicit none
    class(distributionFunction1DLogNormal), intent(inout) :: self

    logNormalMinimum=exp(self%distributionFunction1DNormal%minimum())
    return
  end function logNormalMinimum

  double precision function logNormalMaximum(self)
    !!{
    Return the maximum possible value of a uniform distribution.
    !!}
    implicit none
    class(distributionFunction1DLogNormal), intent(inout) :: self

    logNormalMaximum=exp(self%distributionFunction1DNormal%maximum())
    return
  end function logNormalMaximum

  double precision function logNormalDensity(self,x)
    !!{
    Return the density of a normal distribution.
    !!}
    implicit none
    class           (distributionFunction1DLogNormal), intent(inout) :: self
    double precision                                 , intent(in   ) :: x

    if (x > 0.0d0) then
       logNormalDensity=+self%distributionFunction1DNormal%density(log(x)) &
            &           /x
    else
       logNormalDensity=+0.0d0
    end if
    return
  end function logNormalDensity

  double precision function logNormalCumulative(self,x)
    !!{
    Return the cumulative probability of a normal distribution.
    !!}
    implicit none
    class           (distributionFunction1DLogNormal), intent(inout) :: self
    double precision                                 , intent(in   ) :: x

    if (x <= 0.0d0) then
       logNormalCumulative=0.0d0
    else
       logNormalCumulative=self%distributionFunction1DNormal%cumulative(log(x))
    end if
    return
  end function logNormalCumulative

  double precision function logNormalInverse(self,p)
    !!{
    Return the inverse of a normal distribution.
    !!}
    implicit none
    class           (distributionFunction1DLogNormal), intent(inout), target :: self
    double precision                                 , intent(in   )         :: p

    logNormalInverse=exp(self%distributionFunction1DNormal%inverse(p))
    return
  end function logNormalInverse
