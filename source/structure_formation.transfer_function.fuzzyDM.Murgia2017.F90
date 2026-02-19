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
  Implements a transfer function class for fuzzy dark matter using the fitting function of
  \cite{murgia_non-cold_2017}.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <transferFunction name="transferFunctionFuzzyDMMurgia2017">
   <description>A transfer function class for fuzzy dark matter using the fitting function of \cite{murgia_non-cold_2017}.</description>
  </transferFunction>
  !!]
  type, extends(transferFunctionMurgia2017) :: transferFunctionFuzzyDMMurgia2017
     !!{
     A transfer function class for fuzzy dark matter using the fitting function of \cite{murgia_non-cold_2017}.
     !!}
     private
     double precision                                   :: m22
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
  end type transferFunctionFuzzyDMMurgia2017
   
  interface transferFunctionFuzzyDMMurgia2017
     !!{
     Constructors for the \refClass{transferFunctionFuzzyDMMurgia2017} transfer function class.
     !!}
     module procedure fuzzyDMMurgia2017ConstructorParameters
     module procedure fuzzyDMMurgia2017ConstructorInternal
  end interface transferFunctionFuzzyDMMurgia2017

contains

  function fuzzyDMMurgia2017ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{transferFunctionFuzzyDMMurgia2017} transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionFuzzyDMMurgia2017)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (transferFunctionClass            ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (darkMatterParticleClass          ), pointer       :: darkMatterParticle_
    double precision                                                   :: beta                , gamma, &
         &                                                                redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunction')) call Error_Report("an explicit 'transferFunction' must be given"//{introspection:location})
    ! Read parameters.

    !![
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>5.475d0</defaultValue>
      <defaultSource>\citep[][average of values in Table~4]{murgia_non-cold_2017}</defaultSource>
      <description>The parameter $\beta$, which controls the shape of the cut-off, appearing in the transfer function \citep{murgia_non-cold_2017}.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <defaultValue>-2.0d0</defaultValue>
      <defaultSource>\citep[][average of values in Table~4]{murgia_non-cold_2017}</defaultSource>
      <description>The parameter $\gamma$, which controls the shape of the cut-off, appearing in the transfer function \citep{murgia_non-cold_2017}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%equalityEpochMatterRadiation(requestTypeExpansionFactor))</defaultValue>
      <description>The redshift of the epoch at which the transfer function is defined.</description>
    </inputParameter>
    !!]
    self=transferFunctionFuzzyDMMurgia2017(transferFunctionCDM,beta,gamma,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunctionCDM" />
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function fuzzyDMMurgia2017ConstructorParameters

  function fuzzyDMMurgia2017ConstructorInternal(transferFunctionCDM,beta,gamma,time,cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the \refClass{transferFunctionFuzzyDMMurgia2017} transfer function class.
    !!}
    use :: Dark_Matter_Particles       , only : darkMatterParticleFuzzyDarkMatter
    use :: Error                       , only : Error_Report
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    type            (transferFunctionFuzzyDMMurgia2017)                         :: self
    class           (transferFunctionClass            ), target, intent(in   ) :: transferFunctionCDM
    double precision                                           , intent(in   ) :: beta                , gamma, &
         &                                                                        time
    class           (cosmologyParametersClass         ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterParticleClass          ), target, intent(in   ) :: darkMatterParticle_
    double precision                                                           :: alpha
    double precision                                                           :: wavenumberHalfMode
    !![
    <constructorAssign variables="*darkMatterParticle_"/>
    !!]

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       if (darkMatterParticle__%densityFraction() == 1.0d0) then
          self%m22=+darkMatterParticle__%mass() &
               &   *kilo                        &
               &   /1.0d-22
       else
          call Error_Report('transfer function is not implemented for a mixed CDM and fuzzy dark matter model'//{introspection:location})
       end if
    class default
       call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
    end select
    wavenumberHalfMode=+4.5d0                   &
         &             *self%m22**(4.0d0/9.0d0)
    alpha             =+(                       &
         &               -1.0d0                 &
         &               +2.0d0**(-0.5d0/gamma) &
         &              )**(1.0d0/beta)         &
         &             /wavenumberHalfMode
    self%transferFunctionMurgia2017=transferFunctionMurgia2017(transferFunctionCDM,alpha,beta,gamma,time,cosmologyParameters_,cosmologyFunctions_)
    return
  end function fuzzyDMMurgia2017ConstructorInternal
