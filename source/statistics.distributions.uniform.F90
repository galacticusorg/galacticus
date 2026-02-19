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
  Implementation of a uniform 1D distribution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DUniform">
   <description>
    A uniform distribution over a finite range
    \begin{equation}
    P(x) \propto \left\{ \begin{array}{ll} 1 &amp; \hbox{ if } x_\mathrm{l} \leq x \leq x_\mathrm{u} \\ 0 &amp; \hbox{ otherwise.}  \end{array} \right.
    \end{equation}
    Specified using:
    \begin{description}
    \item[{\normalfont \ttfamily [minimum]}] The lower limit of the range, $x_\mathrm{l}$;
    \item[{\normalfont \ttfamily [maximum]}] The upper limit of the range, $x_\mathrm{u}$.
    \end{description}
   </description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DUniform
     !!{
     Implementation of a uniform 1D distribution function.
     !!}
     private
     double precision :: limitLower, limitUpper
   contains
     procedure :: density    => uniformDensity
     procedure :: cumulative => uniformCumulative
     procedure :: inverse    => uniformInverse
     procedure :: minimum    => uniformMinimum
     procedure :: maximum    => uniformMaximum
  end type distributionFunction1DUniform

  interface distributionFunction1DUniform
     !!{
     Constructors for the \refClass{distributionFunction1DUniform} 1D distribution function class.
     !!}
     module procedure uniformConstructorParameters
     module procedure uniformConstructorInternal
  end interface distributionFunction1DUniform

contains

  function uniformConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DUniform} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DUniform)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass   ), pointer       :: randomNumberGenerator_
    double precision                                               :: limitLower            , limitUpper

    !![
    <inputParameter>
      <name>limitLower</name>
      <description>The lower limit of the uniform distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>limitUpper</name>
      <description>The upper limit of the uniform distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DUniform(limitLower,limitUpper,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function uniformConstructorParameters

  function uniformConstructorInternal(limitLower,limitUpper,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DUniform} 1D distribution function class.
    !!}
    type            (distributionFunction1DUniform)                                  :: self
    double precision                               , intent(in   )                   :: limitLower            , limitUpper
    class           (randomNumberGeneratorClass   ), intent(in   ), target, optional :: randomNumberGenerator_
    !![
    <constructorAssign variables="limitLower, limitUpper, *randomNumberGenerator_"/>
    !!]

    return
  end function uniformConstructorInternal

  double precision function uniformMinimum(self)
    !!{
    Return the minimum possible value of a uniform distribution.
    !!}
    implicit none
    class(distributionFunction1DUniform), intent(inout) :: self

    uniformMinimum=self%limitLower
    return
  end function uniformMinimum

  double precision function uniformMaximum(self)
    !!{
    Return the maximum possible value of a uniform distribution.
    !!}
    implicit none
    class(distributionFunction1DUniform), intent(inout) :: self

    uniformMaximum=self%limitUpper
    return
  end function uniformMaximum

  double precision function uniformDensity(self,x)
    !!{
    Return the density of a uniform distribution.
    !!}
    implicit none
    class           (distributionFunction1DUniform), intent(inout) :: self
    double precision                               , intent(in   ) :: x

    if (x < self%limitLower .or. x > self%limitUpper) then
       uniformDensity=0.0d0
    else
       uniformDensity=1.0d0/(self%limitUpper-self%limitLower)
    end if
    return
  end function uniformDensity

  double precision function uniformCumulative(self,x)
    !!{
    Return the cumulative probability of a uniform distribution.
    !!}
    implicit none
    class           (distributionFunction1DUniform), intent(inout) :: self
    double precision                               , intent(in   ) :: x

    if      (x < self%limitLower) then
       uniformCumulative=0.0d0
    else if (x > self%limitUpper) then
       uniformCumulative=1.0d0
    else
       uniformCumulative=(x-self%limitLower)/(self%limitUpper-self%limitLower)
    end if
    return
  end function uniformCumulative

  double precision function uniformInverse(self,p)
    !!{
    Return the inverse of a uniform distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunction1DUniform), intent(inout), target :: self
    double precision                               , intent(in   )         :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    uniformInverse=self%limitLower+p*(self%limitUpper-self%limitLower)
    return
  end function uniformInverse
