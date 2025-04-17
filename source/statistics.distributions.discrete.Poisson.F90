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
  Implementation of a Poisson 1D discrete distribution function.
  !!}

  !![
  <distributionFunctionDiscrete1D name="distributionFunctionDiscrete1DPoisson">
   <description>A Poisson 1D discrete distribution function class.</description>
  </distributionFunctionDiscrete1D>
  !!]
  type, extends(distributionFunctionDiscrete1DClass) :: distributionFunctionDiscrete1DPoisson
     !!{
     Implementation of a Poisson 1D discrete distribution function.
     !!}
     private
     double precision :: mean
   contains
     procedure :: mass            => poissonMass
     procedure :: massLogarithmic => poissonMassLogarithmic
     procedure :: cumulative      => poissonCumulative
     procedure :: inverse         => poissonInverse
     procedure :: minimum         => poissonMinimum
     procedure :: maximum         => poissonMaximum
  end type distributionFunctionDiscrete1DPoisson

  interface distributionFunctionDiscrete1DPoisson
     !!{
     Constructors for the {\normalfont \ttfamily poisson} 1D discrete distribution function class.
     !!}
     module procedure poissonConstructorParameters
     module procedure poissonConstructorInternal
  end interface distributionFunctionDiscrete1DPoisson

contains

  function poissonConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily poisson} 1D discrete distribution function class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunctionDiscrete1DPoisson)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass           ), pointer       :: randomNumberGenerator_
    double precision                                                       :: mean

    !![
    <inputParameter>
      <name>mean</name>
      <description>The mean of the distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunctionDiscrete1DPoisson(mean,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function poissonConstructorParameters

  function poissonConstructorInternal(mean,randomNumberGenerator_) result(self)
    !!{
    Constructor for {\normalfont \ttfamily poisson} 1D distribution function class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (distributionFunctionDiscrete1DPoisson)                                  :: self
    double precision                                       , intent(in   )                   :: mean
    class           (randomNumberGeneratorClass           ), intent(in   ), target, optional :: randomNumberGenerator_
    !![
    <constructorAssign variables="mean, *randomNumberGenerator_"/>
    !!]

    if (mean <= 0.0d0) call Error_Report('λ ∈ (0,∞]'//{introspection:location})
    return
  end function poissonConstructorInternal

  double precision function poissonMass(self,x) result(mass)
    !!{
    Return the mass of a Poisson discrete distribution.
    !!}
    implicit none
    class  (distributionFunctionDiscrete1DPoisson), intent(inout) :: self
    integer                                       , intent(in   ) :: x

    mass=exp(self%massLogarithmic(x))
    return
  end function poissonMass

  double precision function poissonMassLogarithmic(self,x) result(massLogarithmic)
    !!{
    Return the logarithmic mass of a Poisson discrete distribution.
    !!}
    use :: Error     , only : Error_Report
    use :: Factorials, only : Logarithmic_Factorial
    implicit none
    class  (distributionFunctionDiscrete1DPoisson), intent(inout) :: self
    integer                                       , intent(in   ) :: x

    if (x < 0) call Error_Report('k ∈ [0,∞]'//{introspection:location})
    massLogarithmic=+dble                 (x)*log(self%mean) &
         &          -                             self%mean  &
         &          -Logarithmic_Factorial(x)
    return
  end function poissonMassLogarithmic

  double precision function poissonCumulative(self,x) result(cumulativeDistribution)
    !!{
    Return the cumulative probability of a Poisson discrete distribution.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete
    implicit none
    class  (distributionFunctionDiscrete1DPoisson), intent(inout) :: self
    integer                                       , intent(in   ) :: x

    cumulativeDistribution=Gamma_Function_Incomplete(dble(x+1),self%mean)
    return
  end function poissonCumulative

  integer function poissonInverse(self,p) result(count)
    !!{
    Return the inverse of a Poisson discrete distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunctionDiscrete1DPoisson), intent(inout) :: self
    double precision                                       , intent(in   ) :: p
    !$GLC attributes unused :: self, p

    count=0
    call Error_Report('inverse function is not implemented'//{introspection:location})
    return
  end function poissonInverse

  integer function poissonMinimum(self) result(minimum)
    !!{
    Return the minimum possible value in a Poisson discrete distribution.
    !!}
    implicit none
    class(distributionFunctionDiscrete1DPoisson), intent(inout) :: self
    !$GLC attributes unused :: self

    minimum=0
    return
  end function poissonMinimum

  integer function poissonMaximum(self) result(maximum)
    !!{
    Return the maximum possible value in a Poisson discrete distribution.
    !!}
    implicit none
    class(distributionFunctionDiscrete1DPoisson), intent(inout) :: self
    !$GLC attributes unused :: self

    maximum=huge(0)
    return
  end function poissonMaximum
