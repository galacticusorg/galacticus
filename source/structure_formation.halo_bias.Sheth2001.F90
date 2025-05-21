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
  Implementation of halo bias using the algorithm of \cite{sheth_ellipsoidal_2001}.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <darkMatterHaloBias name="darkMatterHaloBiasSheth2001">
   <description>
    A dark matter halo mass bias class utilizing the algorithm of \cite{sheth_ellipsoidal_2001}.
   </description>
  </darkMatterHaloBias>
  !!]
  type, extends(darkMatterHaloBiasClass) :: darkMatterHaloBiasSheth2001
     !!{
     Implementation of a dark matter halo mass utilizing the algorithm of \cite{sheth_ellipsoidal_2001}.
     !!}
     private
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
   contains
     final     ::               sheth2001Destructor
     procedure :: biasByMass => sheth2001BiasByMass
  end type darkMatterHaloBiasSheth2001

  interface darkMatterHaloBiasSheth2001
     !!{
     Constructors for the \refClass{darkMatterHaloBiasSheth2001} dark matter halo bias class.
     !!}
     module procedure sheth2001ConstructorParameters
     module procedure sheth2001ConstructorInternal
  end interface darkMatterHaloBiasSheth2001

contains

  function sheth2001ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloBiasSheth2001} dark matter halo mass bias which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(darkMatterHaloBiasSheth2001)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    class(criticalOverdensityClass       ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterHaloBiasSheth2001(criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function sheth2001ConstructorParameters

  function sheth2001ConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloBiasSheth2001} dark matter halo bias class.
    !!}
    implicit none
    type (darkMatterHaloBiasSheth2001  )                        :: self
    class(criticalOverdensityClass     ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="*criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function sheth2001ConstructorInternal

  subroutine sheth2001Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloBiasSheth2001} dark matter halo bias class.
    !!}
    implicit none
    type(darkMatterHaloBiasSheth2001), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine sheth2001Destructor

  double precision function sheth2001BiasByMass(self,mass,time,radius)
    !!{
    Returns the bias of a dark matter halo given the mass and time.
    !!}
    implicit none
    class           (darkMatterHaloBiasSheth2001), intent(inout)           :: self
    double precision                             , intent(in   )           :: mass                 , time
    double precision                             , intent(in   ), optional :: radius
    double precision                             , parameter               :: a            =0.707d0, b    =0.5d0, &
         &                                                                    c            =0.600d0
    double precision                                                       :: deltaCritical        , sigma      , &
         &                                                                    nu
    !$GLC attributes unused :: radius

    ! Get critical overdensity for collapse and root-variance, then compute peak height parameter, nu.
    deltaCritical=self%criticalOverdensity_     %value       (time=time,mass=mass)
    sigma        =self%cosmologicalMassVariance_%rootVariance(time=time,mass=mass)
    nu           =+deltaCritical                                                   &
         &        /sigma
    ! Compute halo bias.
    sheth2001BiasByMass=+1.0d0                    &
         &              +(                        &
         &                +sqrt(a)                &
         &                *     a                 &
         &                *   nu**2               &
         &                +sqrt(a)                &
         &                *  b                    &
         &                *  (a*nu**2)**(1.0d0-c) &
         &                -  (a*nu**2)**(     +c) &
         &                /(                      &
         &                  +(a*nu**2)**(     +c) &
         &                  +b                    &
         &                  *(1.0d0-c      )      &
         &                  *(1.0d0-c/2.0d0)      &
         &                 )                      &
         &               )                        &
         &              /sqrt(a)                  &
         &              /deltaCritical
    return
  end function sheth2001BiasByMass
