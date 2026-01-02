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
  Implementation of halo bias using the algorithm of \cite{tinker_large_2010}.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Virial_Density_Contrast   , only : virialDensityContrastClass

  !![
  <darkMatterHaloBias name="darkMatterHaloBiasTinker2010">
   <description>
    A dark matter halo mass bias class utilizing the algorithm of \cite{tinker_large_2010}. The bias is computed at the
    appropriate virial overdensity (see \refPhysics{virialDensityContrast}).
   </description>
  </darkMatterHaloBias>
  !!]
  type, extends(darkMatterHaloBiasClass) :: darkMatterHaloBiasTinker2010
     !!{
     Implementation of a dark matter halo mass utilizing the algorithm of \cite{tinker_large_2010}.
     !!}
     private
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (virialDensityContrastClass   ), pointer :: virialDensityContrast_    => null()
     double precision                                         :: timePrevious                       , massPrevious, &
          &                                                      lowerA                             , upperA      , &
          &                                                      upperC
   contains
     final     ::               tinker2010Destructor
     procedure :: biasByMass => tinker2010BiasByMass
  end type darkMatterHaloBiasTinker2010

  interface darkMatterHaloBiasTinker2010
     !!{
     Constructors for the \refClass{darkMatterHaloBiasTinker2010} dark matter halo bias class.
     !!}
     module procedure tinker2010ConstructorParameters
     module procedure tinker2010ConstructorInternal
  end interface darkMatterHaloBiasTinker2010

contains

  function tinker2010ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloBiasTinker2010} dark matter halo mass bias which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(darkMatterHaloBiasTinker2010  )                :: self
    type(inputParameters               ), intent(inout) :: parameters
    class(criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class(virialDensityContrastClass   ), pointer       :: virialDensityContrast_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="virialDensityContrast"    name="virialDensityContrast_"   source="parameters"/>
    !!]
    self=darkMatterHaloBiasTinker2010(criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="virialDensityContrast_"   />
    !!]
    return
  end function tinker2010ConstructorParameters

  function tinker2010ConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloBiasTinker2010} dark matter halo bias class.
    !!}
    implicit none
    type (darkMatterHaloBiasTinker2010 )                        :: self
    class(criticalOverdensityClass     ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass), intent(in   ), target :: cosmologicalMassVariance_
    class(virialDensityContrastClass   ), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="*criticalOverdensity_, *cosmologicalMassVariance_, *virialDensityContrast_"/>
    !!]

    self%timePrevious=-1.0d0
    self%massPrevious=-1.0d0
    return
  end function tinker2010ConstructorInternal

  subroutine tinker2010Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloBiasTinker2010} dark matter halo bias class.
    !!}
    implicit none
    type(darkMatterHaloBiasTinker2010), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%virialDensityContrast_"   />
    !!]
    return
  end subroutine tinker2010Destructor

  double precision function tinker2010BiasByMass(self,mass,time,radius)
    !!{
    Returns the bias of a dark matter halo given the mass and time.
    !!}
    implicit none
    class           (darkMatterHaloBiasTinker2010), intent(inout)           :: self
    double precision                              , intent(in   )           :: mass                 , time
    double precision                              , intent(in   ), optional :: radius
    double precision                              , parameter               :: lowerB       =1.500d0, lowerC=2.400d0, &
         &                                                                     upperB       =0.183d0
    double precision                                                        :: deltaCritical        , nu            , &
         &                                                                     sigma                , y
    !$GLC attributes unused :: radius

    ! Get critical overdensity for collapse and root-variance, then compute peak height parameter, nu.
    deltaCritical=self%criticalOverdensity_     %value       (time=time,mass=mass)
    sigma        =self%cosmologicalMassVariance_%rootVariance(time=time,mass=mass)
    if (sigma <= 0.0d0) then
       tinker2010BiasByMass=+1.0d0
    else
       nu           =+deltaCritical                                                   &
            &        /sigma
       ! Update fitting parameters if the time has changed.
       if (time /= self%timePrevious .or. mass /= self%massPrevious) then
          ! Store the new time and mass.
          self%timePrevious=time
          self%massPrevious=mass
          ! Compute logarithm of halo density contrast
          y     =log10(self%virialDensityContrast_%densityContrast(mass,time))
          ! Compute parameters as a function of halo overdensity (from Table 2 of Tinker et al. 2010)
          self%upperA=+1.000d0          +0.24d0*y*exp(-(4.0d0/y)**4)
          self%lowerA=-0.880d0+0.440d0*y
          self%upperC=+0.019d0+0.107d0*y+0.19d0  *exp(-(4.0d0/y)**4)
       end if
       ! Compute halo bias (equation 6 of Tinker et al. 2010).
       tinker2010BiasByMass=+1.0d0                        &
            &               -                 self%upperA &
            &               *  nu           **self%lowerA &
            &               /(                            &
            &                 +nu           **self%lowerA &
            &                 +deltaCritical**self%lowerA &
            &                )                            &
            &               +                      upperB &
            &               *  nu           **     lowerB &
            &               +                 self%upperC &
            &               *  nu           **     lowerC
    end if
    return
  end function tinker2010BiasByMass
