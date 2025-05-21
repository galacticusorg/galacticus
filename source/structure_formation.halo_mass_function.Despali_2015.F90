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
Implements a \cite{despali_universality_2015} dark matter halo mass function class.
!!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass, virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt

  !![
  <haloMassFunction name="haloMassFunctionDespali2015">
   <description>
    A dark matter halo mass function class using the function given by \cite{despali_universality_2015}. This uses the
    functional form proposed by \cite{sheth_ellipsoidal_2001} but with parameters $a$, $p$, and $A$ set using eqn.~(12) of
    \cite{despali_universality_2015}.
   </description>
   <deepCopy>
    <functionClass variables="referenceDensityContrast"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="referenceDensityContrast"/>
   </stateStorable>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionShethTormen) :: haloMassFunctionDespali2015
     !!{
     A halo mass function class using the fitting function of \cite{despali_universality_2015}.
     !!}
     private
     class(cosmologyFunctionsClass                                       ), pointer :: cosmologyFunctions_      => null()
     class(virialDensityContrastClass                                    ), pointer :: virialDensityContrast_   => null()
     type (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: referenceDensityContrast => null()
    contains
     !![
     <methods>
       <method description="Return the parameter $x$ in the \cite{despali_universality_2015} halo mass function fit." method="x" />
     </methods>
     !!]
     final     ::                  despali2015Destructor
     procedure :: x             => despali2015X
     procedure :: a             => despali2015A
     procedure :: p             => despali2015P
     procedure :: normalization => despali2015Normalization
  end type haloMassFunctionDespali2015

  interface haloMassFunctionDespali2015
     !!{
     Constructors for the \refClass{haloMassFunctionDespali2015} halo mass function class.
     !!}
     module procedure despali2015ConstructorParameters
     module procedure despali2015ConstructorInternal
  end interface haloMassFunctionDespali2015

contains

  function despali2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionDespali2015} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (haloMassFunctionDespali2015  )                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class(virialDensityContrastClass   ), pointer       :: virialDensityContrast_
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class(criticalOverdensityClass     ), pointer       :: criticalOverdensity_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="virialDensityContrast"    name="virialDensityContrast_"    source="parameters"/>
    !!]
    self=haloMassFunctionDespali2015(cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="virialDensityContrast_"   />
    !!]
    return
  end function despali2015ConstructorParameters

  function despali2015ConstructorInternal(cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionDespali2015} halo mass function class.
    !!}
    implicit none
    type (haloMassFunctionDespali2015  )                        :: self
    class(cosmologyParametersClass     ), target, intent(in   ) :: cosmologyParameters_
    class(cosmologyFunctionsClass      ), target, intent(in   ) :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass), target, intent(in   ) :: cosmologicalMassVariance_
    class(criticalOverdensityClass     ), target, intent(in   ) :: criticalOverdensity_
    class(virialDensityContrastClass   ), target, intent(in   ) :: virialDensityContrast_
    !![
    <constructorAssign variables="*cosmologyParameters_,*cosmologyFunctions_,*cosmologicalMassVariance_,*criticalOverdensity_,*virialDensityContrast_"/>
    !!]

    allocate(self%referenceDensityContrast)
    !![
    <referenceConstruct isResult="yes" owner="self" object="referenceDensityContrast" constructor="virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(.true.,cosmologyFunctions_)"/>
    !!]
    return
  end function despali2015ConstructorInternal

  subroutine despali2015Destructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionDespali2015} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionDespali2015), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"     />
    <objectDestructor name="self%virialDensityContrast_"  />
    <objectDestructor name="self%referenceDensityContrast"/>
    !!]
    return
  end subroutine despali2015Destructor

  double precision function despali2015A(self,time,mass)
    !!{
    Return the parameter $a$ in the {\normalfont \ttfamily despali2015} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionDespali2015), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass
    double precision                                             :: x

    x=self%x(time,mass)
    despali2015A=+0.4332d0*x**2 &
         &       +0.2263d0*x    &
         &       +0.7665d0
    return
  end function despali2015A

  double precision function despali2015P(self,time,mass)
    !!{
    Return the parameter $p$ in the {\normalfont \ttfamily despali2015} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionDespali2015), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass
    double precision                                             :: x

    x=self%x(time,mass)
    despali2015P=-0.1151d0*x**2 &
         &       +0.2554d0*x    &
         &       +0.2488d0
    return
  end function despali2015P

  double precision function despali2015Normalization(self,time,mass)
    !!{
    Return the normalization, $A$, in the {\normalfont \ttfamily despali2015} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionDespali2015), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass
    double precision                                             :: x

    x=self%x(time,mass)
    despali2015Normalization=-0.1362d0*x    &
         &                   +0.3292d0
    return
  end function despali2015Normalization

  double precision function despali2015X(self,time,mass)
    !!{
    Return the parameter $x$ in the {\normalfont \ttfamily despali2015} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionDespali2015), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass

    despali2015X=log10(                                                                    &
         &             +self%virialDensityContrast_  %densityContrast(mass=mass,time=time) &
         &             /self%referenceDensityContrast%densityContrast(mass=mass,time=time) &
         &            )
    return
  end function despali2015X

