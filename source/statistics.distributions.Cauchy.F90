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

  !!{
  Implementation of a 1D Cauchy distribution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DCauchy">
   <description>
    A Cauchy distribution:
    \begin{equation}
     P(x) \propto \left[1+{x-x_0\over\gamma}\right]^{-1}.
    \end{equation}
    Specified using:
    \begin{description}
    \item[{\normalfont \ttfamily [median]}] The median, $x_0$;
    \item[{\normalfont \ttfamily [scale]}] The scale, $\gamma$;
    \end{description}
   </description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DCauchy
     !!{
     Implementation of a 1D Cauchy distribution function.
     !!}
     private
     double precision :: median, scale
   contains
     procedure :: density    => cauchyDensity
     procedure :: cumulative => cauchyCumulative
     procedure :: inverse    => cauchyInverse
  end type distributionFunction1DCauchy

  interface distributionFunction1DCauchy
     !!{
     Constructors for the \refClass{distributionFunction1DCauchy} 1D distribution function class.
     !!}
     module procedure cauchyConstructorParameters
     module procedure cauchyConstructorInternal
     module procedure cauchyConstructorProbability
  end interface distributionFunction1DCauchy

contains

  function cauchyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DCauchy} 1D distribution function class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DCauchy)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass  ), pointer       :: randomNumberGenerator_
    double precision                                              :: median                , scale

    !![
    <inputParameter>
      <name>median</name>
      <description>The median of the Cauchy distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scale</name>
      <description>The scale parameter of the Cauchy distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DCauchy(median,scale,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cauchyConstructorParameters

  function cauchyConstructorInternal(median,scale,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DCauchy} 1D distribution function class.
    !!}
    type            (distributionFunction1DCauchy)                                  :: self
    double precision                              , intent(in   )                   :: median                , scale
    class           (randomNumberGeneratorClass  ), intent(in   ), target, optional :: randomNumberGenerator_
   !![
   <constructorAssign variables="median, scale, *randomNumberGenerator_"/>
   !!]

    return
  end function cauchyConstructorInternal

  function cauchyConstructorProbability(median,limit,probability,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DCauchy} 1D distribution function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    type            (distributionFunction1DCauchy)                :: self
    double precision                              , intent(in   ) :: median                , limit, &
         &                                                           probability
    class           (randomNumberGeneratorClass  ), intent(in   ) :: randomNumberGenerator_

    self=distributionFunction1DCauchy(median,limit/tan(0.5d0*Pi*(1.0d0-probability)),randomNumberGenerator_)
    return
  end function cauchyConstructorProbability

  double precision function cauchyDensity(self,x)
    !!{
    Return the density of a Cauchy distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DCauchy), intent(inout) :: self
    double precision                              , intent(in   ) :: x

    cauchyDensity=1.0d0/Pi/self%scale/(1.0d0+((x-self%median)/self%scale)**2)
    return
  end function cauchyDensity

  double precision function cauchyCumulative(self,x)
    !!{
    Return the cumulative probability of a Cauchy distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DCauchy), intent(inout) :: self
    double precision                              , intent(in   ) :: x

    cauchyCumulative=0.5d0+atan((x-self%median)/self%scale)/Pi
    return
  end function cauchyCumulative

  double precision function cauchyInverse(self,p)
    !!{
    Return the inverse of a Cauchy distribution.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DCauchy), intent(inout), target :: self
    double precision                              , intent(in   )         :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    cauchyInverse=self%median+self%scale*tan(Pi*(p-0.5d0))
    return
  end function cauchyInverse
