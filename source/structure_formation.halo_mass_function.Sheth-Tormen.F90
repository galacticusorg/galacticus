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
Implements a \cite{sheth_ellipsoidal_2001} dark matter halo mass function class.
!!}
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <haloMassFunction name="haloMassFunctionShethTormen">
   <description>
    A dark matter halo mass function class using the function given by \cite{sheth_ellipsoidal_2001}.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionShethTormen
     !!{
     A halo mass function class using the fitting function of \cite{sheth_ellipsoidal_2001}.
     !!}
     private
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     double precision                                         :: a_                                 , p_, &
          &                                                      normalization_
   contains
     !![
     <methods>
       <method description="Return the parameter $a$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit." method="a" />
       <method description="Return the parameter $p$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit." method="p" />
       <method description="Return the parameter $A$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit." method="normalization" />
     </methods>
     !!]
     final     ::                  shethTormenDestructor
     procedure :: differential  => shethTormenDifferential
     procedure :: a             => shethTormenA
     procedure :: p             => shethTormenP
     procedure :: normalization => shethTormenNormalization
  end type haloMassFunctionShethTormen

  interface haloMassFunctionShethTormen
     !!{
     Constructors for the \refClass{haloMassFunctionShethTormen} halo mass function class.
     !!}
     module procedure shethTormenConstructorParameters
     module procedure shethTormenConstructorInternal
  end interface haloMassFunctionShethTormen

contains

  function shethTormenConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionShethTormen} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionShethTormen  )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    double precision                                               :: a                        , p, &
         &                                                            normalization

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <inputParameter>
      <name>a</name>
      <source>parameters</source>
      <defaultValue>0.707d0</defaultValue>
      <description>The parameter $a$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>p</name>
      <source>parameters</source>
      <defaultValue>0.3d0</defaultValue>
      <description>The parameter $p$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <defaultValue>0.3221836349d0</defaultValue>
      <description>The normalization parameter $A$ in the \cite{sheth_ellipsoidal_2001} halo mass function fit.</description>
    </inputParameter>
    !!]
    self=haloMassFunctionShethTormen(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,a,p,normalization)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    !!]
    return
  end function shethTormenConstructorParameters

  function shethTormenConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,criticalOverdensity_,a,p,normalization) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionShethTormen} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionShethTormen  )                        :: self
    class           (cosmologyParametersClass     ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass), target, intent(in   ) :: cosmologicalMassVariance_
    class           (criticalOverdensityClass     ), target, intent(in   ) :: criticalOverdensity_
    double precision                                       , intent(in   ) :: a                        , p, &
         &                                                                    normalization
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologicalMassVariance_, *criticalOverdensity_"/>
    !!]

    self%a_            =a
    self%p_            =p
    self%normalization_=normalization
   return
  end function shethTormenConstructorInternal

  subroutine shethTormenDestructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionShethTormen} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionShethTormen), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"      />
    <objectDestructor name="self%cosmologicalMassVariance_" />
    <objectDestructor name="self%criticalOverdensity_"      />
    !!]
    return
  end subroutine shethTormenDestructor

  double precision function shethTormenDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (haloMassFunctionShethTormen), intent(inout), target   :: self
    double precision                             , intent(in   )           :: time   , mass
    type            (treeNode                   ), intent(inout), optional :: node
    double precision                                                       :: alpha  , nu          , &
         &                                                                    nuPrime, massVariance

    ! Determine the mass variance. If zero, return zero mass function.
    massVariance=self%cosmologicalMassVariance_%rootVariance(mass,time)
    if (massVariance <= 0.0d0) then
       shethTormenDifferential=0.0d0
       return
    end if
    ! Compute the mass function.
    nu                     =+(                                                                &
         &                    +self%criticalOverdensity_%value(time=time,mass=mass,node=node) &
         &                    /massVariance                                                   &
         &                   )**2
    nuPrime                =+self%a(time,mass)                                                &
         &                  *nu
    alpha                  =+abs(self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time))
    shethTormenDifferential=+self%cosmologyParameters_%OmegaMatter    () &
         &                  *self%cosmologyParameters_%densityCritical() &
         &                  /mass**2                                     &
         &                  *alpha                                       &
         &                  *sqrt(                                       &
         &                        +2.0d0                                 &
         &                        *nuPrime                               &
         &                        /Pi                                    &
         &                       )                                       &
         &                  *self%normalization(time,mass)               &
         &                  *(                                           &
         &                    +1.0d0                                     &
         &                    +1.0d0                                     &
         &                    /nuPrime**self%p(time,mass)                &
         &                  )                                            &
         &                  *exp(                                        &
         &                       -0.5d0                                  &
         &                       *nuPrime                                &
         &                  )
    return
  end function shethTormenDifferential

  double precision function shethTormenA(self,time,mass)
    !!{
    Return the parameter $a$ in the {\normalfont \ttfamily shethTormen} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormen), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass
    !$GLC attributes unused :: time, mass

    shethTormenA=self%a_
    return
  end function shethTormenA

  double precision function shethTormenP(self,time,mass)
    !!{
    Return the parameter $p$ in the {\normalfont \ttfamily shethTormen} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormen), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass
    !$GLC attributes unused :: time, mass

    shethTormenP=self%p_
    return
  end function shethTormenP

  double precision function shethTormenNormalization(self,time,mass)
    !!{
    Return the normalization, $A$, in the {\normalfont \ttfamily shethTormen} halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionShethTormen), intent(inout) :: self
    double precision                             , intent(in   ) :: time , mass
    !$GLC attributes unused :: time, mass

    shethTormenNormalization=self%normalization_
    return
  end function shethTormenNormalization

