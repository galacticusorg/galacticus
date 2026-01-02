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
  Implements a \cite{reed_halo_2007} dark matter halo mass function class.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <haloMassFunction name="haloMassFunctionReed2007">
   <description>The halo mass function is computed from the function given by \cite{reed_halo_2007}.</description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionReed2007
     !!{
     A halo mass function class using the fitting function of \cite{reed_halo_2007}.
     !!}
     private
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
   contains
     final     ::                  reed2007Destructor
     procedure :: differential  => reed2007Differential
  end type haloMassFunctionReed2007

  interface haloMassFunctionReed2007
     !!{
     Constructors for the \refClass{haloMassFunctionReed2007} halo mass function class.
     !!}
     module procedure reed2007ConstructorParameters
     module procedure reed2007ConstructorInternal
  end interface haloMassFunctionReed2007

contains

  function reed2007ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionReed2007} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (haloMassFunctionReed2007     )                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class(criticalOverdensityClass     ), pointer       :: criticalOverdensity_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    !!]
    self=haloMassFunctionReed2007(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    !!]
    return
  end function reed2007ConstructorParameters

  function reed2007ConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionReed2007} halo mass function class.
    !!}
    implicit none
    type (haloMassFunctionReed2007     )                        :: self
    class(cosmologyParametersClass     ), target, intent(in   ) :: cosmologyParameters_
    class(cosmologicalMassVarianceClass), target, intent(in   ) :: cosmologicalMassVariance_
    class(criticalOverdensityClass     ), target, intent(in   ) :: criticalOverdensity_
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologicalMassVariance_, *criticalOverdensity_"/>
    !!]

    return
  end function reed2007ConstructorInternal

  subroutine reed2007Destructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionReed2007} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionReed2007), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"      />
    <objectDestructor name="self%cosmologicalMassVariance_" />
    <objectDestructor name="self%criticalOverdensity_"      />
    !!]
    return
  end subroutine reed2007Destructor

  double precision function reed2007Differential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (haloMassFunctionReed2007), intent(inout), target   :: self
    double precision                          , intent(in   )           :: time                               , mass
    type            (treeNode                ), intent(inout), optional :: node
    ! Parameter values from Reed et al. (2007), text after equations (11) and (12).
    double precision                          , parameter               :: c                          =1.080d0, a         =0.764d0/c, &
         &                                                                 normalization              =0.310d0, p         =0.300d0
    double precision                                                    :: rootVariance                       , peakHeight          , &
         &                                                                 peakHeightScaled                   , alpha               , &
         &                                                                 G1                                 , G2                  , &
         &                                                                 powerSpectrumSlopeEffective

    ! Determine the mass variance. If zero, return zero mass function.
    rootVariance=self%cosmologicalMassVariance_%rootVariance(mass,time)
    if (rootVariance <= 0.0d0) then
       reed2007Differential=0.0d0
       return
    end if
    ! Compute the mass function.
    peakHeight                 =+self%criticalOverdensity_%value(time=time,mass=mass,node=node) &
         &                      /     rootVariance
    alpha                      =self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time)
    !! Scaled peak height defined by Reed et al. (2007; test after equation 9).
    peakHeightScaled           =+sqrt(                                                    &
         &                            +c                                                  &
         &                            *a                                                  &
         &                           )                                                    &
         &                      *peakHeight
    !! Reed et al. (2007; equation 12).
    G1                         =+exp(-0.50d0*(log(peakHeightScaled)-0.788d0)**2/0.6d0**2)
    G2                         =+exp(-0.50d0*(log(peakHeightScaled)-1.138d0)**2/0.2d0**2)
    !! Reed et al. (2007; equation 13).
    powerSpectrumSlopeEffective=-6.0d0                                                    &
         &                      *alpha                                                    &
         &                      -3.0d0
    if (powerSpectrumSlopeEffective <= -3.0d0) then
       reed2007Differential=0.0d0
    else
       !! Reed et al. (2007; equation 12). Note that the published version has some typos. Specifically, in the exponential term
       !! the "w" (scaled peak height) parameter should be squared, but is not. The expression below has been validated against
       !! Darren Reed's "genmf" code.
       reed2007Differential=+self%cosmologyParameters_%OmegaMatter    ()         &
            &               *self%cosmologyParameters_%densityCritical()         &
            &               /mass**2                                             &
            &               *(-alpha)                                            &
            &               *normalization                                       &
            &               *peakHeightScaled                                    &
            &               *sqrt(                                               &
            &                     +2.0d0                                         &
            &                     /Pi                                            &
            &                    )                                               &
            &               *(                                                   &
            &                 +1.0d0                                             &
            &                 +1.02d0                                            &
            &                 /peakHeightScaled**(2.0d0*p)                       &
            &                 +0.6d0*G1                                          &
            &                 +0.4d0*G2                                          &
            &                )                                                   &
            &               *exp(                                                &
            &                    -0.5000d0*peakHeightScaled          ** 2        &
            &                    -0.0325d0*peakHeightScaled          **(2.0d0*p) &
            &                    /(3.0d0+powerSpectrumSlopeEffective)** 2        &
            &                   )
    end if
    return
  end function reed2007Differential
