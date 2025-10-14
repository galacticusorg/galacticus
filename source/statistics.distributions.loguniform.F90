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
  Implementation of a 1D distribution function which is uniform in the logarithm of the variable.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DLogUniform">
   <description>
    A distribution uniform in the logarithm of $x$ over a finite range
    \begin{equation}
     P(x) \propto \left\{ \begin{array}{ll} x^{-1} &amp; \hbox{ if } x_\mathrm{l} \leq x \leq x_\mathrm{u} \\ 0 &amp; \hbox{ otherwise.}  \end{array} \right.
    \end{equation}
    Specified using:
    \begin{description}
    \item[{\normalfont \ttfamily [minimum]}] The lower limit of the range, $x_\mathrm{l}$;
    \item[{\normalfont \ttfamily [maximum]}] The upper limit of the range, $x_\mathrm{u}$.
    \end{description}
   </description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DLogUniform
     !!{
     Implementation of a 1D distribution function which is uniform in the logarithm of the variable.
     !!}
     private
     double precision :: limitLower, limitUpper
   contains
     procedure :: density    => logUniformDensity
     procedure :: cumulative => logUniformCumulative
     procedure :: inverse    => logUniformInverse
     procedure :: minimum    => logUniformMinimum
     procedure :: maximum    => logUniformMaximum
  end type distributionFunction1DLogUniform

  interface distributionFunction1DLogUniform
     !!{
     Constructors for the \refClass{distributionFunction1DLogUniform} 1D distribution function class.
     !!}
     module procedure logUniformConstructorParameters
     module procedure logUniformConstructorInternal
  end interface distributionFunction1DLogUniform

contains

  function logUniformConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DLogUniform} 1D distribution function class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DLogUniform)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass      ), pointer       :: randomNumberGenerator_
    double precision                                                  :: limitLower            , limitUpper

    !![
    <inputParameter>
      <name>limitLower</name>
      <description>The lower limit of the log-uniform distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>limitUpper</name>
      <description>The upper limit of the log-uniform distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DLogUniform(limitLower,limitUpper,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function logUniformConstructorParameters

  function logUniformConstructorInternal(limitLower,limitUpper,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DLogUniform} 1D distribution function class.
    !!}
    use :: Error, only : Error_Report
    type            (distributionFunction1DLogUniform)                                  :: self
    double precision                                  , intent(in   )                   :: limitLower            , limitUpper
    class           (randomNumberGeneratorClass      ), intent(in   ), target, optional :: randomNumberGenerator_
    !![
    <constructorAssign variables="limitLower, limitUpper, *randomNumberGenerator_"/>
    !!]

    if (limitLower <= 0.0d0     ) call Error_Report('`limitLower` > 0 is required'           //{introspection:location})
    if (limitUpper <= 0.0d0     ) call Error_Report('`limitUpper` > 0 is required'           //{introspection:location})
    if (limitLower >= limitUpper) call Error_Report('`limitLower` < `limitUpper` is required'//{introspection:location})
    return
  end function logUniformConstructorInternal

  double precision function logUniformDensity(self,x)
    !!{
    Return the density of a uniform distribution.
    !!}
    implicit none
    class           (distributionFunction1DLogUniform), intent(inout) :: self
    double precision                                  , intent(in   ) :: x

    if (x < self%limitLower .or. x > self%limitUpper) then
       logUniformDensity=0.0d0
    else
       logUniformDensity=1.0d0/x/(log(self%limitUpper)-log(self%limitLower))
    end if
    return
  end function logUniformDensity

  double precision function logUniformCumulative(self,x)
    !!{
    Return the cumulative probability of a uniform distribution.
    !!}
    implicit none
    class           (distributionFunction1DLogUniform), intent(inout) :: self
    double precision                                  , intent(in   ) :: x

    if      (x < self%limitLower) then
       logUniformCumulative=0.0d0
    else if (x > self%limitUpper) then
       logUniformCumulative=1.0d0
    else
       logUniformCumulative=(log(x)-log(self%limitLower))/(log(self%limitUpper)-log(self%limitLower))
    end if
    return
  end function logUniformCumulative

  double precision function logUniformInverse(self,p)
    !!{
    Return the inverse of a uniform distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunction1DLogUniform), intent(inout), target :: self
    double precision                                  , intent(in   )         :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    logUniformInverse=exp(log(self%limitLower)+p*(log(self%limitUpper)-log(self%limitLower)))
    return
  end function logUniformInverse

  double precision function logUniformMinimum(self)
    !!{
    Return the minimum possible value of a logUniform distribution.
    !!}
    implicit none
    class(distributionFunction1DLogUniform), intent(inout) :: self

    logUniformMinimum=self%limitLower
    return
  end function logUniformMinimum

  double precision function logUniformMaximum(self)
    !!{
    Return the maximum possible value of a logUniform distribution.
    !!}
    implicit none
    class(distributionFunction1DLogUniform), intent(inout) :: self

    logUniformMaximum=self%limitUpper
    return
  end function logUniformMaximum
