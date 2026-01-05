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
  Implementation of halo bias using the Press-Schechter algorithm \citep{cole_biased_1989}.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <darkMatterHaloBias name="darkMatterHaloBiasPressSchechter">
   <description>
    A dark matter halo mass bias class consistent with the halo mass function of \cite{press_formation_1974} (see
    \citep{mo_analytic_1996}).
   </description>
  </darkMatterHaloBias>
  !!]
  type, extends(darkMatterHaloBiasClass) :: darkMatterHaloBiasPressSchechter
     !!{
     Implementation of a dark matter halo mass utilizing the Press-Schechter algorithm \citep{cole_biased_1989}.
     !!}
     private
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
   contains
     final     ::               pressSchechterDestructor
     procedure :: biasByMass => pressSchechterBiasByMass
  end type darkMatterHaloBiasPressSchechter

  interface darkMatterHaloBiasPressSchechter
     !!{
     Constructors for the \refClass{darkMatterHaloBiasPressSchechter} dark matter halo bias class.
     !!}
     module procedure pressSchechterConstructorParameters
     module procedure pressSchechterConstructorInternal
  end interface darkMatterHaloBiasPressSchechter

contains

  function pressSchechterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloBiasPressSchechter} dark matter halo mass bias which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(darkMatterHaloBiasPressSchechter)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    class(criticalOverdensityClass       ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterHaloBiasPressSchechter(criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function pressSchechterConstructorParameters

  function pressSchechterConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloBiasPressSchechter} dark matter halo bias class.
    !!}
    implicit none
    type (darkMatterHaloBiasPressSchechter)                        :: self
    class(criticalOverdensityClass        ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass   ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="*criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function pressSchechterConstructorInternal

  subroutine pressSchechterDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloBiasPressSchechter} dark matter halo bias class.
    !!}
    implicit none
    type(darkMatterHaloBiasPressSchechter), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine pressSchechterDestructor

  double precision function pressSchechterBiasByMass(self,mass,time,radius)
    !!{
    Returns the bias of a dark matter halo given the mass and time.
    !!}
    implicit none
    class           (darkMatterHaloBiasPressSchechter), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass         , time
    double precision                                  , intent(in   ), optional :: radius
    double precision                                                            :: deltaCritical, sigma, &
         &                                                                         nu
    !$GLC attributes unused :: radius
    
    ! Get critical overdensity for collapse and root-variance, then compute peak height parameter, nu.
    deltaCritical=+self%criticalOverdensity_     %value       (time=time,mass=mass)
    sigma        =+self%cosmologicalMassVariance_%rootVariance(time=time,mass=mass)
    nu           =+deltaCritical                                                    &
         &        /sigma
    ! Compute halo bias.
    pressSchechterBiasByMass=1.0d0+(nu**2-1.0d0)/deltaCritical
    return
  end function pressSchechterBiasByMass
