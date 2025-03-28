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
  Implementation of a beta density 1D distribution function.
  !!}

  use :: Root_Finder, only : rootFinder

  !![
  <distributionFunction1D name="distributionFunction1DBeta">
   <description>A beta 1D distribution function class.</description>
  </distributionFunction1D>
  !!]
    type, extends(distributionFunction1DClass) :: distributionFunction1DBeta
     !!{
     Implementation of a beta 1D distribution function.
     !!}
     private
     double precision             :: alpha , beta
     type            (rootFinder) :: finder
   contains
     procedure :: density    => betaDensity
     procedure :: cumulative => betaCumulative
     procedure :: inverse    => betaInverse
     procedure :: minimum    => betaMinimum
     procedure :: maximum    => betaMaximum
  end type distributionFunction1DBeta

  interface distributionFunction1DBeta
     !!{
     Constructors for the {\normalfont \ttfamily beta} 1D distribution function class.
     !!}
     module procedure betaConstructorParameters
     module procedure betaConstructorInternal
  end interface distributionFunction1DBeta

  ! Module-scope variables used in finding the inverse of the distribution.
  class           (distributionFunction1DBeta), pointer :: self_
  double precision                                      :: probabilityCumulative
  !$omp threadprivate(self_,probabilityCumulative)

contains

  function betaConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily beta} 1D distribution function class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DBeta)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    double precision                                            :: alpha                 , beta

    !![
    <inputParameter>
      <name>alpha</name>
      <description>The $\alpha$ parameter of the beta distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <description>The $\beta$ parameter of the beta distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DBeta(alpha,beta,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function betaConstructorParameters

  function betaConstructorInternal(alpha,beta,randomNumberGenerator_) result(self)
    !!{
    Constructor for {\normalfont \ttfamily beta} 1D distribution function class.
    !!}
    type            (distributionFunction1DBeta)                                  :: self
    double precision                            , intent(in   )                   :: alpha                        , beta
    class           (randomNumberGeneratorClass), intent(in   ), target, optional :: randomNumberGenerator_
    double precision                            , parameter                       :: toleranceAbsolute     =1.0d-6, toleranceRelative=0.0d+0
    !![
    <constructorAssign variables="alpha, beta, *randomNumberGenerator_"/>
    !!]
    
    self%finder=rootFinder(                                     &
         &                 rootFunction     =betaRoot         , &
         &                 toleranceAbsolute=toleranceAbsolute, &
         &                 toleranceRelative=toleranceRelative  &
         &                )
    return
  end function betaConstructorInternal

  double precision function betaDensity(self,x)
    !!{
    Return the density of a beta distribution.
    !!}
    use :: Beta_Functions, only : Beta_Function
    implicit none
    class           (distributionFunction1DBeta), intent(inout) :: self
    double precision                            , intent(in   ) :: x

    if (x < 0.0d0 .or. x > 1.0d0) then
       betaDensity=0.0d0
    else
       betaDensity=+       x **(self%alpha-1.0d0)       &
            &      *(1.0d0-x)**(self%beta -1.0d0)       &
            &      /Beta_Function(self%alpha,self%beta)
    end if
    return
  end function betaDensity

  double precision function betaCumulative(self,x)
    !!{
    Return the cumulative probability of a beta distribution.
    !!}
    use :: Beta_Functions, only : Beta_Function_Incomplete_Normalized
    implicit none
    class           (distributionFunction1DBeta), intent(inout) :: self
    double precision                            , intent(in   ) :: x

    if      (x < 0.0d0) then
       betaCumulative=0.0d0
    else if (x > 1.0d0) then
       betaCumulative=1.0d0
    else
       betaCumulative=Beta_Function_Incomplete_Normalized(self%alpha,self%beta,x)
    end if
    return
  end function betaCumulative

  double precision function betaInverse(self,p)
    !!{
    Return the inverse of a beta distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunction1DBeta), intent(inout), target :: self
    double precision                            , intent(in   )         :: p

    if      (                                         &
         &    p < 0.0d0                               &
         &   .or.                                     &
         &    p > 1.0d0                               &
         &  ) then
       betaInverse=+0.0d0
       call Error_Report(                             &
            &            'probability out of range'// &
            &            {introspection:location}     &
            &           )
    else
       self_                 => self
       probabilityCumulative =  p
       betaInverse           =  self%finder%find(rootRange=[0.0d0,1.0d0])
    end if
    return
  end function betaInverse

  double precision function betaRoot(x)
    !!{
    Root function used in finding the inverse of the beta distribution.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    betaRoot=self_%cumulative(x)-probabilityCumulative
    return
  end function betaRoot

  double precision function betaMinimum(self)
    !!{
    Return the minimum value of a beta distribution.
    !!}
    implicit none
    class(distributionFunction1DBeta), intent(inout) :: self
    !$GLC attributes unused :: self

    betaMinimum=0.0d0
    return
  end function betaMinimum

  double precision function betaMaximum(self)
    !!{
    Return the minimum value of a beta distribution.
    !!}
    implicit none
    class(distributionFunction1DBeta), intent(inout) :: self
    !$GLC attributes unused :: self

    betaMaximum=1.0d0
    return
  end function betaMaximum
