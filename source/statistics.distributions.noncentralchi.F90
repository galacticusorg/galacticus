!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  <distributionFunction1D name="distributionFunction1DNonCentralChi">
   <description>
    A non-central chi^2 distribution with 3 degrees of freedom.
   </description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DNonCentralChi
     !!{
     Implementation of a non-central chi^2 1D distibution function.
     !!}
     private
     double precision :: lambda
   contains
     procedure :: density    => nonCentralChiDensity
     procedure :: cumulative => nonCentralChiCumulative
     procedure :: minimum    => nonCentralChiMinimum
     procedure :: maximum    => nonCentralChiMaximum
  end type distributionFunction1DNonCentralChi

  interface distributionFunction1DNonCentralChi
     !!{
     Constructors for the {\normalfont \ttfamily non-central chi^2} 1D distribution function class.
     !!}
     module procedure nonCentralChiConstructorParameters
     module procedure nonCentralChiConstructorInternal
  end interface distributionFunction1DNonCentralChi

contains

  function nonCentralChiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily normal} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DNonCentralChi)         :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass  ), pointer       :: randomNumberGenerator_
    double precision                                              :: lambda
    !![
    <inputParameter>
      <name>lambda</name>
      <description>Non centrality parameter</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DNonCentralChi(lambda,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nonCentralChiConstructorParameters

  function nonCentralChiConstructorInternal(lambda,randomNumberGenerator_) result(self)
    !!{
    Constructor for non central chi^2 1D distribution function class.
    !!}
    use :: Error_Functions, only : Error_Function
    implicit none
    type            (distributionFunction1DNonCentralChi)                           :: self
    double precision                              , intent(in   )                   :: lambda
    class           (randomNumberGeneratorClass  ), intent(in   ), optional, target :: randomNumberGenerator_
    !![
    <constructorAssign variables="lambda, *randomNumberGenerator_"/>
    !!]
    return
  end function nonCentralChiConstructorInternal

  double precision function nonCentralChiMinimum(self)
    !!{
    Return the minimum possible value of a non central chi^2 distribution.
    !!}
    implicit none
    class(distributionFunction1DNonCentralChi), intent(inout) :: self
    nonCentralChiMinimum=0.0d0
    return
  end function nonCentralChiMinimum

  double precision function nonCentralChiMaximum(self)
    !!{
    Return the maximum possible value of a non central chi^2 distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(distributionFunction1DNonCentralChi), intent(inout) :: self
    nonCentralChiMaximum=0.0d0
    call Error_Report('no maximum exists'//{introspection:location})
    return
  end function nonCentralChiMaximum

  double precision function nonCentralChiDensity(self,x)
    !!{
    Return the density of a non central chi^2 distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Bessel_Functions, only : Bessel_Function_In
    implicit none
    class           (distributionFunction1DNonCentralChi), intent(inout) :: self
    double precision                                     , intent(in   ) :: x

    if (x < 0.0d0) then
       nonCentralChiDensity=0.0d0
    else
       nonCentralChiDensity=+0.5d0*exp(-(x + self%lambda)/2.0d0)&
                           &*(x/self%lambda)**0.25d0            &
                           &*Bessel_Function_In(sqrt(x*self%lambda), 0.5d0)
    end if
    return
  end function nonCentralChiDensity

  double precision function nonCentralChiCumulative(self,x)
    !!{
    Return the cumulative probability of a non central chi^2 distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DNonCentralChi), intent(inout) :: self
    double precision                                     , intent(in   ) :: x

    if (x < 0.0d0) then
       nonCentralChiCumulative=0.0d0
    else
       nonCentralChiCumulative=+1.0d0 -(&
                              &+0.5d0*erfc((sqrt(x) - sqrt(self%lambda))/sqrt(2.0d0))     &
                              &+0.5d0*erfc((sqrt(x) + sqrt(self%lambda))/sqrt(2.0d0))     &
                              &+sqrt(2.0d0/Pi)*sinh(sqrt(x*self%lambda))/sqrt(self%lambda)&
                              &*exp(-(x + self%lambda)/2.0d0)                             &
                              &)
    end if
    return
  end function nonCentralChiCumulative
