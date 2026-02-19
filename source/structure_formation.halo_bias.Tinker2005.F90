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
  Implementation of radial dependence of halo bias using the fitting function of \citep{tinker_mass--light_2005}.
  !!}

  use :: Correlation_Functions_Two_Point, only : correlationFunctionTwoPointClass

  !![
  <darkMatterHaloBias name="darkMatterHaloBiasTinker2005">
   <description>
   A dark matter halo bias class which applies the radial dependence fitting function of \citep{tinker_mass--light_2005} to
   another bias class.
   </description>
  </darkMatterHaloBias>
  !!]
  type, extends(darkMatterHaloBiasClass) :: darkMatterHaloBiasTinker2005
     !!{
     Implementation of radial dependence of halo bias using the fitting function of \citep{tinker_mass--light_2005}.
     !!}
     private
     class(darkMatterHaloBiasClass         ), pointer :: darkMatterHaloBias_                   => null()
     class(correlationFunctionTwoPointClass), pointer :: correlationFunctionTwoPointNonLinear_ => null()
   contains
     final     ::               tinker2005Destructor
     procedure :: biasByMass => tinker2005BiasByMass
  end type darkMatterHaloBiasTinker2005

  interface darkMatterHaloBiasTinker2005
     !!{
     Constructors for the \refClass{darkMatterHaloBiasTinker2005} dark matter halo bias class.
     !!}
     module procedure tinker2005ConstructorParameters
     module procedure tinker2005ConstructorInternal
  end interface darkMatterHaloBiasTinker2005

contains

  function tinker2005ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloBiasTinker2005} dark matter halo mass bias which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterHaloBiasTinker2005    )                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(darkMatterHaloBiasClass         ), pointer       :: darkMatterHaloBias_
    class(correlationFunctionTwoPointClass), pointer       :: correlationFunctionTwoPointNonLinear_

    !![
    <objectBuilder class="darkMatterHaloBias"          name="darkMatterHaloBias_"                                                                        source="parameters"/>
    <objectBuilder class="correlationFunctionTwoPoint" name="correlationFunctionTwoPointNonlinear_" parameterName="correlationFunctionTwoPointNonlinear" source="parameters"/>
    !!]
    self=darkMatterHaloBiasTinker2005(darkMatterHaloBias_,correlationFunctionTwoPointNonLinear_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloBias_"/>
    !!]
    return
  end function tinker2005ConstructorParameters

  function tinker2005ConstructorInternal(darkMatterHaloBias_,correlationFunctionTwoPointNonLinear_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloBiasTinker2005} dark matter halo bias class.
    !!}
    implicit none
    type (darkMatterHaloBiasTinker2005    )                        :: self
    class(darkMatterHaloBiasClass         ), intent(in   ), target :: darkMatterHaloBias_
    class(correlationFunctionTwoPointClass), intent(in   ), target :: correlationFunctionTwoPointNonLinear_
    !![
    <constructorAssign variables="*correlationFunctionTwoPointNonLinear_, *darkMatterHaloBias_"/>
    !!]

    return
  end function tinker2005ConstructorInternal

  subroutine tinker2005Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloBiasTinker2005} dark matter halo bias class.
    !!}
    implicit none
    type(darkMatterHaloBiasTinker2005), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloBias_"                  />
    <objectDestructor name="self%correlationFunctionTwoPointNonlinear_"/>
    !!]
    return
  end subroutine tinker2005Destructor

  double precision function tinker2005BiasByMass(self,mass,time,radius)
    !!{
    Returns the bias of a dark matter halo given the mass and time.
    !!}
    implicit none
    class           (darkMatterHaloBiasTinker2005), intent(inout)           :: self
    double precision                              , intent(in   )           :: mass               , time
    double precision                              , intent(in   ), optional :: radius
    double precision                                                        :: correlationFunction, dependenceRadial

    ! Get the large scale bias.
    tinker2005BiasByMass=self%darkMatterHaloBias_%bias(mass,time)
    ! Apply radial dependence.
    if (present(radius)) then
       correlationFunction=+self%correlationFunctionTwoPointNonlinear_%correlation(radius,time)
       ! Fitting function from eqn. (B7) of Tinker et al. (2005; ApJ; 631; 41). 
       dependenceRadial   =+(1.0d0+1.17d0*correlationFunction)**1.49d0 &
            &              /(1.0d0+0.69d0*correlationFunction)**2.09d0
       ! Multiply the large scale by by the square-root of the radial dependence. The fitting function for the radial dependence,
       ! ζ(r), is defined such that:
       !
       !  b²(M,r)=b²(M) ζ(r).
       !
       ! Here we are computing just one power of bias, so we apply the square root.
       tinker2005BiasByMass=tinker2005BiasByMass*sqrt(dependenceRadial)
    end if
    return
  end function tinker2005BiasByMass
