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
  Implementation of a negative exponential density 1D distribution function.
  !!}

  !![
  <distributionFunction1D name="distributionFunction1DNegativeExponential">
   <description>A negative exponential 1D distribution function class.</description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DNegativeExponential
     !!{
     Implementation of a negative exponential 1D distribution function.
     !!}
     private
     double precision :: rate
   contains
     procedure :: density    => negativeExponentialDensity
     procedure :: cumulative => negativeExponentialCumulative
     procedure :: inverse    => negativeExponentialInverse
  end type distributionFunction1DNegativeExponential

  interface distributionFunction1DNegativeExponential
     !!{
     Constructors for the \refClass{distributionFunction1DNegativeExponential} 1D distribution function class.
     !!}
     module procedure negativeExponentialConstructorParameters
     module procedure negativeExponentialConstructorInternal
  end interface distributionFunction1DNegativeExponential

contains

  function negativeExponentialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DNegativeExponential} 1D distribution function class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DNegativeExponential)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass               ), pointer       :: randomNumberGenerator_
    double precision                                                           :: rate

    !![
    <inputParameter>
      <name>rate</name>
      <description>The rate parameter of the negative exponential distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DNegativeExponential(rate,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function negativeExponentialConstructorParameters

  function negativeExponentialConstructorInternal(rate,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DNegativeExponential} 1D distribution function class.
    !!}
    type            (distributionFunction1DNegativeExponential)                                  :: self
    double precision                                           , intent(in   )                   :: rate
    class           (randomNumberGeneratorClass               ), intent(in   ), target, optional :: randomNumberGenerator_
    !![
    <constructorAssign variables="rate, *randomNumberGenerator_"/>
    !!]

    return
  end function negativeExponentialConstructorInternal

  double precision function negativeExponentialDensity(self,x)
    !!{
    Return the density of a negative exponential distribution.
    !!}
    implicit none
    class           (distributionFunction1DNegativeExponential), intent(inout) :: self
    double precision                                           , intent(in   ) :: x

    if (x < 0.0d0) then
       negativeExponentialDensity=0.0d0
    else
       negativeExponentialDensity=self%rate*exp(-self%rate*x)
    end if
    return
  end function negativeExponentialDensity

  double precision function negativeExponentialCumulative(self,x)
    !!{
    Return the cumulative probability of a negative exponential distribution.
    !!}
    implicit none
    class           (distributionFunction1DNegativeExponential), intent(inout) :: self
    double precision                                           , intent(in   ) :: x

    if (x < 0.0d0) then
       negativeExponentialCumulative=0.0d0
    else
       negativeExponentialCumulative=1.0d0-exp(-self%rate*x)
    end if
    return
  end function negativeExponentialCumulative

  double precision function negativeExponentialInverse(self,p)
    !!{
    Return the inverse of a negative exponential distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunction1DNegativeExponential), intent(inout), target :: self
    double precision                                           , intent(in   )         :: p
    double precision                                           , parameter             :: pTiny=1.0d-6

    if      (                                         &
         &    p < 0.0d0                               &
         &   .or.                                     &
         &    p > 1.0d0                               &
         &  ) then
       negativeExponentialInverse=+0.0d0
       call Error_Report(                             &
            &            'probability out of range'// &
            &            {introspection:location}     &
            &           )
    else if (                                         &
         &    p > pTiny                               &
         &  ) then
       negativeExponentialInverse=-log(               &
            &                          +1.0d0         &
            &                          -p             &
            &                         )               &
            &                     /self%rate
    else
       negativeExponentialInverse=+(                  &
            &                       +p**1/1.0d0       &
            &                       +p**2/2.0d0       &
            &                       +p**3/3.0d0       &
            &                       +p**4/4.0d0       &
            &                      )                  &
            &                     /self%rate
    end if
    return
  end function negativeExponentialInverse
