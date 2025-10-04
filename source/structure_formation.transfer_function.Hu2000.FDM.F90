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
  Implements a transfer function class based on the fuzzy dark matter modifier of \cite{hu_fuzzy_2000}.
  !!}

  use :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <transferFunction name="transferFunctionHu2000FDM">
   <description>Provides a transfer function based on the fuzzy dark matter modifier of \cite{hu_fuzzy_2000}.</description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionHu2000FDM
     !!{
     A transfer function class which modifies another transfer function using the fuzzy dark matter modifier of
     \cite{hu_fuzzy_2000}.
     !!}
     private
     double precision                                    :: jeansWavenumberEq             , m22     , &
          &                                                 time                          , redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
   contains
     final     ::                          hu2000FDMDestructor
     procedure :: value                 => hu2000FDMValue
     procedure :: logarithmicDerivative => hu2000FDMLogarithmicDerivative
     procedure :: halfModeMass          => hu2000FDMHalfModeMass
     procedure :: quarterModeMass       => hu2000FDMQuarterModeMass
     procedure :: epochTime             => hu2000FDMEpochTime
  end type transferFunctionHu2000FDM

  interface transferFunctionHu2000FDM
     !!{
     Constructors for the \refClass{transferFunctionHu2000FDM} transfer function class.
     !!}
     module procedure hu2000FDMConstructorParameters
     module procedure hu2000FDMConstructorInternal
  end interface transferFunctionHu2000FDM

contains

  function hu2000FDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{transferFunctionHu2000FDM} transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Error                         , only : Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionHu2000FDM)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    class           (transferFunctionClass    ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_
    class           (darkMatterParticleClass  ), pointer       :: darkMatterParticle_
    double precision                                           :: redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunction')) call Error_Report("an explicit 'transferFunction' must be given"//{introspection:location})
    ! Read parameters.
    !![
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
    self=transferFunctionHu2000FDM(transferFunctionCDM,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunctionCDM" />
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function hu2000FDMConstructorParameters
  
  function hu2000FDMConstructorInternal(transferFunctionCDM,time,cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the \refClass{transferFunctionHu2000FDM} transfer function class.
    !!}
    use :: Cosmology_Parameters        , only : hubbleUnitsLittleH
    use :: Error                       , only : Error_Report
    use :: Dark_Matter_Particles       , only : darkMatterParticleFuzzyDarkMatter
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    type            (transferFunctionHu2000FDM)                        :: self
    class           (transferFunctionClass    ), target, intent(in   ) :: transferFunctionCDM
    double precision                                   , intent(in   ) :: time
    class           (cosmologyParametersClass ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass  ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterParticleClass  ), target, intent(in   ) :: darkMatterParticle_
    !![
    <constructorAssign variables="*transferFunctionCDM, time, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterParticle_"/>
    !!]

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       if (darkMatterParticle__%densityFraction() == 1.0d0) then
          self%m22              =+darkMatterParticle__%mass() &
               &                 *kilo                        &
               &                 /1.0d-22
          self%jeansWavenumberEq=9.00d0*sqrt(self%m22)
       else
          call Error_Report('transfer function is not implemented for a mixed CDM and fuzzy dark matter model'//{introspection:location})
       end if
    class default
       call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
    end select
    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    return
  end function hu2000FDMConstructorInternal

  subroutine hu2000FDMDestructor(self)
    !!{
    Destructor for the \refClass{transferFunctionHu2000FDM} transfer function class.
    !!}
    implicit none
    type(transferFunctionHu2000FDM), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%transferFunctionCDM" />
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine hu2000FDMDestructor

  double precision function hu2000FDMValue(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionHu2000FDM), intent(inout) :: self
    double precision                           , intent(in   ) :: wavenumber
    double precision                                           :: x

    hu2000FDMValue   =self%transferFunctionCDM%value(wavenumber)
    x                =+1.61d0                   &
         &            *self%m22**(1.0d0/18.0d0) &
         &            *(                        &
         &              +          wavenumber   &
         &              /self%jeansWavenumberEq &
         &             )
    hu2000FDMValue   =+hu2000FDMValue   &
         &            *(                &
         &              +cos(x**3)      &
         &              /(              &
         &                +1.0d0        &
         &                +x**8         &
         &               )              &
         &             )
    return
  end function hu2000FDMValue

  double precision function hu2000FDMLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionHu2000FDM), intent(inout) :: self
    double precision                           , intent(in   ) :: wavenumber
    double precision                                           :: x

    hu2000FDMLogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    x                             =+1.61d0                   &
         &                         *self%m22**(1.0d0/18.0d0) &
         &                         *(                        &
         &                           +          wavenumber   &
         &                           /self%jeansWavenumberEq &
         &                          )
    hu2000FDMLogarithmicDerivative=+hu2000FDMLogarithmicDerivative &
         &                         +(                              &
         &                           -8.0d0                        &
         &                           *       x**8                  &
         &                           /(1.0d0+x**8)                 &
         &                           -3.0d0                        &
         &                           *       x**3                  &
         &                           *tan(x**3)                    &
         &                          )
    return
  end function hu2000FDMLogarithmicDerivative

  double precision function hu2000FDMHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionHu2000FDM), intent(inout), target   :: self
    integer                                    , intent(  out), optional :: status
    double precision                                                     :: matterDensity, wavenumberHalfMode

    matterDensity        =+self%cosmologyParameters_%OmegaMatter    () &
         &                *self%cosmologyParameters_%densityCritical()
    wavenumberHalfMode   =+1.108d0                 &
         &                *4.5d0                   &
         &                *self%m22**(4.0d0/9.0d0)
    ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
    ! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].
    hu2000FDMHalfModeMass=+4.0d0                &
         &                *Pi                   &
         &                /3.0d0                &
         &                *matterDensity        &
         &                *(                    &
         &                  +Pi                 &
         &                  /wavenumberHalfMode &
         &                 )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function hu2000FDMHalfModeMass

  double precision function hu2000FDMQuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionHu2000FDM), intent(inout), target   :: self
    integer                                    , intent(  out), optional :: status
    double precision                                                     :: matterDensity, wavenumberQuarterMode

    matterDensity           =+self%cosmologyParameters_%OmegaMatter    () &
         &                   *self%cosmologyParameters_%densityCritical()
    wavenumberQuarterMode   =+1.230d0                 &
         &                   *4.5d0                   &
         &                   *self%m22**(4.0d0/9.0d0)
    ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
    ! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].
    hu2000FDMQuarterModeMass=+4.0d0                   &
         &                   *Pi                      &
         &                   /3.0d0                   &
         &                   *matterDensity           &
         &                   *(                       &
         &                     +Pi                    &
         &                     /wavenumberQuarterMode &
         &                    )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function hu2000FDMQuarterModeMass

  double precision function hu2000FDMEpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionHu2000FDM), intent(inout) :: self

    hu2000FDMEpochTime=self%time
    return
  end function hu2000FDMEpochTime
