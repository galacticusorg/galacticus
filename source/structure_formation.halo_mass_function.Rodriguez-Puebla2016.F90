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
  Implements a \cite{tinker_towardhalo_2008} dark matter halo mass function class.
  !!}

  !![
  <haloMassFunction name="haloMassFunctionRodriguezPuebla2016">
   <description>The halo mass function is computed from the function given by \cite{tinker_towardhalo_2008}, and using the fits of \cite{rodriguez-puebla_halo_2016} for the parameter values.</description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionTinker2008Form) :: haloMassFunctionRodriguezPuebla2016
     !!{
     A halo mass function class using the fitting function of \cite{tinker_towardhalo_2008}, and using the fits of \cite{rodriguez-puebla_halo_2016} for the parameter values.
     !!}
     private
   contains
     final     ::                  rodriguezPuebla2016Destructor
     procedure :: normalization => rodriguezPuebla2016Normalization
     procedure :: a             => rodriguezPuebla2016A
     procedure :: b             => rodriguezPuebla2016B
     procedure :: c             => rodriguezPuebla2016C
  end type haloMassFunctionRodriguezPuebla2016

  interface haloMassFunctionRodriguezPuebla2016
     !!{
     Constructors for the \refClass{haloMassFunctionRodriguezPuebla2016} halo mass function class.
     !!}
     module procedure rodriguezPuebla2016ConstructorParameters
     module procedure rodriguezPuebla2016ConstructorInternal
  end interface haloMassFunctionRodriguezPuebla2016

contains

  function rodriguezPuebla2016ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionRodriguezPuebla2016} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (haloMassFunctionRodriguezPuebla2016   )                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class(linearGrowthClass            ), pointer       :: linearGrowth_
    class(cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !!]
    self=haloMassFunctionRodriguezPuebla2016(                           &
         &                                   cosmologyParameters_     , &
         &                                   cosmologicalMassVariance_, &
         &                                   linearGrowth_            , &
         &                                   cosmologyFunctions_        &
         &                                  )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologyFunctions_"      />
    !!]
    return
  end function rodriguezPuebla2016ConstructorParameters

  function rodriguezPuebla2016ConstructorInternal(cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionRodriguezPuebla2016} halo mass function class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (haloMassFunctionRodriguezPuebla2016)                             :: self
    class           (cosmologyParametersClass           ), target     , intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass      ), target     , intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass                  ), target     , intent(in   ) :: linearGrowth_
    class           (cosmologyFunctionsClass            ), target     , intent(in   ) :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologicalMassVariance_, *linearGrowth_, *cosmologyFunctions_"/>
    !!]

    self%time=-1.0d0
    self%mass=-1.0d0
    return
  end function rodriguezPuebla2016ConstructorInternal

  subroutine rodriguezPuebla2016Destructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionRodriguezPuebla2016} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionRodriguezPuebla2016), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%cosmologyParameters_"     />
    !!]
    return
  end subroutine rodriguezPuebla2016Destructor

  double precision function rodriguezPuebla2016Normalization(self,time,mass)
    !!{
    Return the normalization for the {\normalfont \ttfamily rodriguezPuebla2016} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionRodriguezPuebla2016), intent(inout) :: self
    double precision                                     , intent(in   ) :: time             , mass
    double precision                                                     :: redshift
    double precision                                     , parameter     :: chi0    =+0.144d0,      &
         &                                                                  chi1    =-0.011d0,      &
         &                                                                  chi2    =+0.003d0
    !$GLC attributes unused :: mass

    redshift                        =self%cosmologyFunctions_%redshiftFromExpansionFactor(      &
         &                           self%cosmologyFunctions_%expansionFactor             (     &
         &                                                                                 time &
         &                                                                                )     &
         &                                                                               )
    rodriguezPuebla2016Normalization=+chi0             &
         &                           +chi1*redshift    &
         &                           +chi2*redshift**2
    return
  end function rodriguezPuebla2016Normalization

  double precision function rodriguezPuebla2016A(self,time,mass)
    !!{
    Return the normalization for the {\normalfont \ttfamily rodriguezPuebla2016} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionRodriguezPuebla2016), intent(inout) :: self
    double precision                                     , intent(in   ) :: time             , mass
    double precision                                                     :: redshift
    double precision                                     , parameter     :: chi0    =+1.351d0,      &
         &                                                                  chi1    =+0.068d0,      &
         &                                                                  chi2    =+0.006d0
    !$GLC attributes unused :: mass

    redshift                        =self%cosmologyFunctions_%redshiftFromExpansionFactor(      &
         &                           self%cosmologyFunctions_%expansionFactor             (     &
         &                                                                                 time &
         &                                                                                )     &
         &                                                                               )
    rodriguezPuebla2016A            =+chi0             &
         &                           +chi1*redshift    &
         &                           +chi2*redshift**2

    return
  end function rodriguezPuebla2016A

  double precision function rodriguezPuebla2016B(self,time,mass)
    !!{
    Return the normalization for the {\normalfont \ttfamily rodriguezPuebla2016} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionRodriguezPuebla2016), intent(inout) :: self
    double precision                                     , intent(in   ) :: time             , mass
    double precision                                                     :: redshift
    double precision                                     , parameter     :: chi0    =+3.113d0,      &
         &                                                                  chi1    =-0.077d0,      &
         &                                                                  chi2    =-0.013d0
    !$GLC attributes unused :: mass

    redshift                        =self%cosmologyFunctions_%redshiftFromExpansionFactor(      &
         &                           self%cosmologyFunctions_%expansionFactor             (     &
         &                                                                                 time &
         &                                                                                )     &
         &                                                                               )
    rodriguezPuebla2016B            =+chi0             &
         &                           +chi1*redshift    &
         &                           +chi2*redshift**2

    return
  end function rodriguezPuebla2016B

  double precision function rodriguezPuebla2016C(self,time,mass)
    !!{
    Return the normalization for the {\normalfont \ttfamily rodriguezPuebla2016} halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionRodriguezPuebla2016), intent(inout) :: self
    double precision                                     , intent(in   ) :: time             , mass
    double precision                                                     :: redshift
    double precision                                     , parameter     :: chi0    =+1.187d0,      &
         &                                                                  chi1    =+0.000d0,      &
         &                                                                  chi2    =+0.000d0
    !$GLC attributes unused :: mass

    redshift                        =self%cosmologyFunctions_%redshiftFromExpansionFactor(      &
         &                           self%cosmologyFunctions_%expansionFactor             (     &
         &                                                                                 time &
         &                                                                                )     &
         &                                                                               )
    rodriguezPuebla2016C            =+chi0             &
         &                           +chi1*redshift    &
         &                           +chi2*redshift**2

    return
  end function rodriguezPuebla2016C

