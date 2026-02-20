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
  Implementation of a degree-3 non-central $\chi^2$ 1D distribution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DNonCentralChiDegree3">
   <description>
    A non-central $\chi^2$ distribution with 3 degrees of freedom.
   </description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DNonCentralChiDegree3
     !!{
     Implementation of a non-central $\chi^2$ distribution function with 3 degrees of freedom.
     !!}
     private
     double precision :: lambda
     double precision :: xMinimum, xMaximum
   contains
     procedure :: density    => nonCentralChiSquaredDegree3Density
     procedure :: cumulative => nonCentralChiSquaredDegree3Cumulative
     procedure :: minimum    => nonCentralChiSquaredDegree3Minimum
     procedure :: maximum    => nonCentralChiSquaredDegree3Maximum
  end type distributionFunction1DNonCentralChiDegree3

  interface distributionFunction1DNonCentralChiDegree3
     !!{
     Constructors for the \refClass{distributionFunction1DNonCentralChiDegree3} 1D distribution function class.
     !!}
     module procedure nonCentralChiSquaredDegree3ConstructorParameters
     module procedure nonCentralChiSquaredDegree3ConstructorInternal
  end interface distributionFunction1DNonCentralChiDegree3

contains

  function nonCentralChiSquaredDegree3ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DNonCentralChiDegree3} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DNonCentralChiDegree3)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass                ), pointer       :: randomNumberGenerator_
    double precision                                                            :: lambda

    !![
    <inputParameter>
      <name>lambda</name>
      <description>Non centrality parameter</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DNonCentralChiDegree3(lambda,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nonCentralChiSquaredDegree3ConstructorParameters

  function nonCentralChiSquaredDegree3ConstructorInternal(lambda,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DNonCentralChiDegree3} 1D distribution function class.
    !!}
    implicit none
    type            (distributionFunction1DNonCentralChiDegree3)                                  :: self
    double precision                                            , intent(in   )                   :: lambda
    class           (randomNumberGeneratorClass                ), intent(in   ), optional, target :: randomNumberGenerator_
    double precision                                            , parameter                        :: offsetMaximum        =10.d0
    !![
    <constructorAssign variables="lambda, *randomNumberGenerator_"/>
    !!]

    ! Compute the x values beyond which we approximate the PDF as zero.
    self%xMinimum=(max(0.0d0,sqrt(self%lambda)-offsetMaximum))**2
    self%xMaximum=(          sqrt(self%lambda)+offsetMaximum )**2
    return
  end function nonCentralChiSquaredDegree3ConstructorInternal

  double precision function nonCentralChiSquaredDegree3Minimum(self)
    !!{
    Return the minimum possible value of a degree-3 non central $\chi^2$ distribution.
    !!}
    implicit none
    class(distributionFunction1DNonCentralChiDegree3), intent(inout) :: self

    nonCentralChiSquaredDegree3Minimum=0.0d0
    return
  end function nonCentralChiSquaredDegree3Minimum

  double precision function nonCentralChiSquaredDegree3Maximum(self)
    !!{
    Return the maximum possible value of a degree-3 non central $\chi^2$ distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(distributionFunction1DNonCentralChiDegree3), intent(inout) :: self

    nonCentralChiSquaredDegree3Maximum=0.0d0
    call Error_Report('no maximum exists'//{introspection:location})
    return
  end function nonCentralChiSquaredDegree3Maximum

  double precision function nonCentralChiSquaredDegree3Density(self,x)
    !!{
    Return the density of a degree-3 non central $\chi^2$ distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Bessel_Functions        , only : Bessel_Function_In
    implicit none
    class           (distributionFunction1DNonCentralChiDegree3), intent(inout) :: self
    double precision                                            , intent(in   ) :: x

    if (x < self%xMinimum .or. x > self%xMaximum) then
       nonCentralChiSquaredDegree3Density=+0.0d0
    else
       nonCentralChiSquaredDegree3Density=+0.5d0                                  &
            &                             *exp(                                   &
            &                                  -(                                 &
            &                                    +     x                          &
            &                                    +self%lambda                     &
            &                                   )                                 &
            &                                  /2.0d0                             &
            &                                 )                                   &
            &                             *(                                      &
            &                               +     x                               &
            &                               /self%lambda                          &
            &                              )**0.25d0                              &
            &                             *Bessel_Function_In(                    &
            &                                                 +0.5d0            , &
            &                                                 +sqrt(              &
            &                                                       +     x       &
            &                                                       *self%lambda  &
            &                                                      )              &
            &                                                )
    end if
    return
  end function nonCentralChiSquaredDegree3Density

  double precision function nonCentralChiSquaredDegree3Cumulative(self,x)
    !!{
    Return the cumulative probability of a degree-3 non central $\chi^2$ distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DNonCentralChiDegree3), intent(inout) :: self
    double precision                                            , intent(in   ) :: x

    if      (x < self%xMinimum) then
       nonCentralChiSquaredDegree3Cumulative=+0.0d0
    else if (x > self%xMaximum) then
       nonCentralChiSquaredDegree3Cumulative=+1.0d0
    else
       ! Evaluate the CDF. Note that the sinh() term (see,
       ! e.g. https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution#Cumulative_distribution_function) is expressed in
       ! terms of exponentials and combined with the other exponential term here. This helps avoid floating point errors when
       ! evaluating the CDF for very large arguments/non-centrality parameters.
       nonCentralChiSquaredDegree3Cumulative=+1.0d0                                                 &
            &                                -(                                                     &
            &                                  +0.5d0*erfc((sqrt(x)-sqrt(self%lambda))/sqrt(2.0d0)) &
            &                                  +0.5d0*erfc((sqrt(x)+sqrt(self%lambda))/sqrt(2.0d0)) &
            &                                  +sqrt(+2.0d0/Pi)                                     &
            &                                  /       sqrt(+  self%lambda)                         &
            &                                  *(                                                   &
            &                                    +exp(                                              &
            &                                         -    (+x+self%lambda)/2.0d0                   &
            &                                         +sqrt(+x*self%lambda)                         &
            &                                        )                                              &
            &                                    -exp(                                              &
            &                                         -    (+x+self%lambda)/2.0d0                   &
            &                                         -sqrt(+x*self%lambda)                         &
            &                                        )                                              &
            &                                   )                                                   &
            &                                  /2.0d0                                               &
            &                   )
    end if              
    return
  end function nonCentralChiSquaredDegree3Cumulative
