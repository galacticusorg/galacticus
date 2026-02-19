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
  Implements a \cite{bhattacharya_mass_2011} dark matter halo mass function class.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <haloMassFunction name="haloMassFunctionBhattacharya2011">
   <description>The halo mass function is computed from the function given by \cite{bhattacharya_mass_2011}.</description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionBhattacharya2011
     !!{
     A halo mass function class using the fitting function of \cite{bhattacharya_mass_2011}.
     !!}
     private
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     double precision                                         :: a_                                 , p_, &
          &                                                      normalization_                     , q_
   contains
     !![
     <methods>
       <method description="Return the parameter $\bar{a}$ in the \cite{bhattacharya_mass_2011} halo mass function fit." method="a" />
       <method description="Return the parameter $\bar{p}$ in the \cite{bhattacharya_mass_2011} halo mass function fit." method="p" />
       <method description="Return the parameter $\bar{q}$ in the \cite{bhattacharya_mass_2011} halo mass function fit." method="q" />
       <method description="Return the parameter $\bar{A}$ in the \cite{bhattacharya_mass_2011} halo mass function fit." method="normalization" />
     </methods>
     !!]
     final     ::                  bhattacharya2011Destructor
     procedure :: differential  => bhattacharya2011Differential
     procedure :: a             => bhattacharya2011A
     procedure :: p             => bhattacharya2011P
     procedure :: q             => bhattacharya2011Q
     procedure :: normalization => bhattacharya2011Normalization
  end type haloMassFunctionBhattacharya2011

  interface haloMassFunctionBhattacharya2011
     !!{
     Constructors for the \refClass{haloMassFunctionBhattacharya2011} halo mass function class.
     !!}
     module procedure bhattacharya2011ConstructorParameters
     module procedure bhattacharya2011ConstructorInternal
  end interface haloMassFunctionBhattacharya2011

contains

  function bhattacharya2011ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionBhattacharya2011} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionBhattacharya2011)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass   ), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass        ), pointer       :: criticalOverdensity_
    double precision                                                  :: a                        , p, &
         &                                                               normalization            , q

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <inputParameter>
      <name>a</name>
      <source>parameters</source>
      <defaultValue>0.788d0</defaultValue>
      <defaultSource>\citep{comparat_accurate_2017}</defaultSource>
      <description>The parameter $\bar{a}$ in the \cite{bhattacharya_mass_2011} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>p</name>
      <source>parameters</source>
      <defaultValue>0.807d0</defaultValue>
      <defaultSource>\citep{comparat_accurate_2017}</defaultSource>
      <description>The parameter $\bar{p}$ in the \cite{bhattacharya_mass_2011} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>q</name>
      <source>parameters</source>
      <defaultValue>1.795d0</defaultValue>
      <defaultSource>\citep{comparat_accurate_2017}</defaultSource>
      <description>The parameter $\bar{q}$ in the \cite{bhattacharya_mass_2011} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <defaultValue>0.333d0</defaultValue>
      <defaultSource>\citep{comparat_accurate_2017}</defaultSource>
      <description>The normalization parameter $\bar{A}$ in the \cite{bhattacharya_mass_2011} halo mass function fit.</description>
    </inputParameter>
    !!]
    self=haloMassFunctionBhattacharya2011(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,a,p,q,normalization)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    !!]
    return
  end function bhattacharya2011ConstructorParameters

  function bhattacharya2011ConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,a,p,q,normalization) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionBhattacharya2011} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionBhattacharya2011)                        :: self
    class           (cosmologyParametersClass        ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass   ), target, intent(in   ) :: cosmologicalMassVariance_
    class           (criticalOverdensityClass        ), target, intent(in   ) :: criticalOverdensity_
    double precision                                          , intent(in   ) :: a                        , p, &
         &                                                                       normalization            , q
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologicalMassVariance_, *criticalOverdensity_"/>
    !!]

    self%            a_=a
    self%            p_=p
    self%            q_=q
    self%normalization_=normalization
    return
  end function bhattacharya2011ConstructorInternal

  subroutine bhattacharya2011Destructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionBhattacharya2011} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionBhattacharya2011), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"      />
    <objectDestructor name="self%cosmologicalMassVariance_" />
    <objectDestructor name="self%criticalOverdensity_"      />
    !!]
    return
  end subroutine bhattacharya2011Destructor

  double precision function bhattacharya2011Differential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (haloMassFunctionBhattacharya2011), intent(inout), target   :: self
    double precision                                  , intent(in   )           :: time   , mass
    type            (treeNode                        ), intent(inout), optional :: node
    double precision                                                            :: alpha  , nu          , &
         &                                                                         nuPrime, massVariance

    ! Set a default value.
    bhattacharya2011Differential=0.0d0
    ! Determine the mass variance. If zero, return zero mass function.
    massVariance=self%cosmologicalMassVariance_%rootVariance(mass,time)
    if (massVariance <=    0.0d0) return
    ! Compute the mass function.
    nu                     =+(                                                                &
         &                    +self%criticalOverdensity_%value(time=time,mass=mass,node=node) &
         &                    /massVariance                                                   &
         &                   )**2
    if (nu           <=    0.0d0) return
    nuPrime                =+self%a(time,mass)                                                &
         &                  *nu
    if (nuPrime      >  1500.0d0) return ! Exponential term will be zero beyond this point.
    alpha                  =+abs(self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time))
    bhattacharya2011Differential=+self%cosmologyParameters_%OmegaMatter    () &
         &                       *self%cosmologyParameters_%densityCritical() &
         &                       /mass**2                                     &
         &                       *alpha                                       &
         &                       *sqrt(                                       &
         &                             +2.0d0                                 &
         &                             *nuPrime**self%q(time,mass)            &
         &                             /Pi                                    &
         &                            )                                       &
         &                       *self%normalization(time,mass)               &
         &                       *(                                           &
         &                         +1.0d0                                     &
         &                         +1.0d0                                     &
         &                         /nuPrime**self%p(time,mass)                &
         &                       )                                            &
         &                       *exp(                                        &
         &                            -0.5d0                                  &
         &                            *nuPrime                                &
         &                       )
    return
  end function bhattacharya2011Differential

  double precision function bhattacharya2011A(self,time,mass)
    !!{
    Return the parameter $\bar{a}$ in the {\normalfont \ttfamily bhattacharya2011} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionBhattacharya2011), intent(inout) :: self
    double precision                                  , intent(in   ) :: time , mass
    !$GLC attributes unused :: time, mass

    bhattacharya2011A=self%a_
    return
  end function bhattacharya2011A

  double precision function bhattacharya2011P(self,time,mass)
    !!{
    Return the parameter $\bar{p}$ in the {\normalfont \ttfamily bhattacharya2011} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionBhattacharya2011), intent(inout) :: self
    double precision                                  , intent(in   ) :: time , mass
    !$GLC attributes unused :: time, mass

    bhattacharya2011P=self%p_
    return
  end function bhattacharya2011P

  double precision function bhattacharya2011Q(self,time,mass)
    !!{
    Return the parameter $\bar{q}$ in the {\normalfont \ttfamily bhattacharya2011} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionBhattacharya2011), intent(inout) :: self
    double precision                                  , intent(in   ) :: time , mass
    !$GLC attributes unused :: time, mass

    bhattacharya2011Q=self%q_
    return
  end function bhattacharya2011Q

  double precision function bhattacharya2011Normalization(self,time,mass)
    !!{
    Return the normalization, $\bar{A}$, in the {\normalfont \ttfamily bhattacharya2011} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionBhattacharya2011), intent(inout) :: self
    double precision                                  , intent(in   ) :: time , mass
    !$GLC attributes unused :: time, mass

    bhattacharya2011Normalization=self%normalization_
    return
  end function bhattacharya2011Normalization

