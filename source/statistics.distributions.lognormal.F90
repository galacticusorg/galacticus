!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  Implementation of a normal 1D distibution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DLogNormal">
   <description>A normal 1D distribution function class.</description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DNormal) :: distributionFunction1DLogNormal
     !!{
     Implementation of a normal 1D distibution function.
     !!}
     private
   contains
     procedure :: density    => logNormalDensity
     procedure :: cumulative => logNormalCumulative
     procedure :: inverse    => logNormalInverse
     procedure :: minimum    => logNormalMinimum
     procedure :: maximum    => logNormalMaximum
  end type distributionFunction1DLogNormal

  interface distributionFunction1DLogNormal
     !!{
     Constructors for the {\normalfont \ttfamily normal} 1D distribution function class.
     !!}
     module procedure logNormalConstructorParameters
     module procedure logNormalConstructorInternal
  end interface distributionFunction1DLogNormal

contains

  function logNormalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily normal} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DLogNormal)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass     ), pointer       :: randomNumberGenerator_
    double precision                                                 :: mean                  , variance  , &
         &                                                              limitLower            , limitUpper

    !![
    <inputParameter>
      <name>mean</name>
      <description>The mean of the log-normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>variance</name>
      <description>The variance of the log-normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
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
    !!]
    self=distributionFunction1DLogNormal(mean,variance,limitLower,limitUpper,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function logNormalConstructorParameters

  function logNormalConstructorInternal(mean,variance,limitLower,limitUpper,randomNumberGenerator_) result(self)
    !!{
    Constructor for ``normal'' 1D distribution function class.
    !!}
    type            (distributionFunction1DLogNormal)                                  :: self
    double precision                                 , intent(in   )                   :: mean                  , variance
    class           (randomNumberGeneratorClass     ), intent(in   ), optional, target :: randomNumberGenerator_
    double precision                                 , intent(in   ), optional         :: limitLower            , limitUpper
    double precision                                                                   :: meanNormal            , varianceNormal

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

    logNormalDensity=+self%distributionFunction1DNormal%density(log(x)) &
         &           /x
    return
  end function logNormalDensity

  double precision function logNormalCumulative(self,x)
    !!{
    Return the cumulative probability of a normal distribution.
    !!}
    implicit none
    class           (distributionFunction1DLogNormal), intent(inout) :: self
    double precision                                 , intent(in   ) :: x

    logNormalCumulative=self%distributionFunction1DNormal%cumulative(log(x))
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
