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
Implements a transfer function class based on the \gls{wdm} modifier of \cite{bardeen_statistics_1986}.
!!}

  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <transferFunction name="transferFunctionBBKSWDM">
   <description>
    A transfer function class providing the \gls{wdm} transfer function fitting function of \cite{bardeen_statistics_1986}. The
    free-streaming length is computed from the dark matter particle (which must be of the {\normalfont \ttfamily
    darkMatterParticleWDMThermal} class) properties, and the resulting modifier of \cite{bardeen_statistics_1986} is applied to
    the provided \gls{cdm} transfer function.
   </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionBBKSWDM
     !!{
     A transfer function class which modifies another transfer function using the \gls{wdm} modifier of \cite{bardeen_statistics_1986}.
     !!}
     private
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
     double precision                                    :: lengthFreeStreaming           , time
   contains
     final     ::                          bbksWDMDestructor
     procedure :: value                 => bbksWDMValue
     procedure :: logarithmicDerivative => bbksWDMLogarithmicDerivative
     procedure :: halfModeMass          => bbksWDMHalfModeMass
     procedure :: epochTime             => bbksWDMEpochTime
  end type transferFunctionBBKSWDM

  interface transferFunctionBBKSWDM
     !!{
     Constructors for the ``{\normalfont \ttfamily bbksWDM}'' transfer function class.
     !!}
     module procedure bbksWDMConstructorParameters
     module procedure bbksWDMConstructorInternal
  end interface transferFunctionBBKSWDM

contains

  function bbksWDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``{\normalfont \ttfamily bbksWDM}'' transfer function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (transferFunctionBBKSWDM )                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(transferFunctionClass   ), pointer       :: transferFunctionCDM
    class(cosmologyParametersClass), pointer       :: cosmologyParameters_
    class(darkMatterParticleClass ), pointer       :: darkMatterParticle_
    class(cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=transferFunctionBBKSWDM(transferFunctionCDM,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="transferFunctionCDM" />
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function bbksWDMConstructorParameters

  function bbksWDMConstructorInternal(transferFunctionCDM,cosmologyParameters_,darkMatterParticle_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the ``{\normalfont \ttfamily bbksWDM}'' transfer function class.
    !!}
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleWDMThermal
    use :: Error                , only : Error_Report
    implicit none
    type            (transferFunctionBBKSWDM )                        :: self
    class           (transferFunctionClass   ), target, intent(in   ) :: transferFunctionCDM
    class           (cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_
    class           (darkMatterParticleClass ), target, intent(in   ) :: darkMatterParticle_
    double precision                                                  :: degreesOfFreedomEffectiveDecoupling
    !![
    <constructorAssign variables="*transferFunctionCDM, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterParticle_"/>
    !!]

    ! Get degrees of freedom at the time at which the dark matter particle decoupled.
    select type (particle => self%darkMatterParticle_)
    class is (darkMatterParticleWDMThermal)
       degreesOfFreedomEffectiveDecoupling=particle%degreesOfFreedomEffectiveDecoupling()
    class default
       degreesOfFreedomEffectiveDecoupling=0.0d0
       call Error_Report('transfer function expects a thermal warm dark matter particle'//{introspection:location})
    end select
    ! Compute the epoch - the transfer function is assumed to be for z=0.
    self%time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(0.0d0))
    ! Compute the free-streaming length-like parameter (equation G6 of BBKS).
    self%lengthFreeStreaming=+0.2d0                                                                                     &
         &                                         /(                                                                   &
         &                                           +degreesOfFreedomEffectiveDecoupling                               &
         &                                           /100.d0                                                            &
         &                                          )**(4.0d0/3.0d0)                                                    &
         &                                         /(                                                                   &
         &                                           +(                                                                 &
         &                                             +self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                             -self%cosmologyParameters_%OmegaBaryon   (                  )    &
         &                                            )                                                                 &
         &                                           /  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                          )
    return
  end function bbksWDMConstructorInternal

  subroutine bbksWDMDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily bbksWDM} transfer function class.
    !!}
    implicit none
    type(transferFunctionBBKSWDM), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterParticle_" />
    <objectDestructor name="self%transferFunctionCDM" />
    !!]
    return
  end subroutine bbksWDMDestructor

  double precision function bbksWDMValue(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionBBKSWDM), intent(inout) :: self
    double precision                         , intent(in   ) :: wavenumber
    double precision                                         :: wavenumberScaleFree

    wavenumberScaleFree=+wavenumber               &
         &              *self%lengthFreeStreaming
    bbksWDMValue       =+self%transferFunctionCDM%value(wavenumber) &
         &              *exp(                                       &
         &                   -0.5d0                                 &
         &                   *(                                     &
         &                     +  wavenumberScaleFree               &
         &                     *(                                   &
         &                       +1.0d0                             &
         &                       +wavenumberScaleFree               &
         &                      )                                   &
         &                    )                                     &
         &                  )
    return
  end function bbksWDMValue

  double precision function bbksWDMLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionBBKSWDM), intent(inout) :: self
    double precision                         , intent(in   ) :: wavenumber
    double precision                                         :: wavenumberScaleFree

    wavenumberScaleFree=+wavenumber               &
         &              *self%lengthFreeStreaming
    bbksWDMLogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber) &
         &                       -  wavenumberScaleFree                                      &
         &                       *(                                                          &
         &                         +0.5d0                                                    &
         &                         +wavenumberScaleFree                                      &
         &                        )
    return
  end function bbksWDMLogarithmicDerivative

  double precision function bbksWDMHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionBBKSWDM), intent(inout), target   :: self
    integer                                  , intent(  out), optional :: status
    double precision                         , parameter               :: wavenumberHalfModeScaleFree=sqrt(0.25d0+2.0d0*log(2.0d0))-0.5d0
    double precision                                                   :: matterDensity                                                  , wavenumberHalfMode

    wavenumberHalfMode =+wavenumberHalfModeScaleFree                 &
         &              /self%lengthFreeStreaming
    matterDensity      =+self%cosmologyParameters_%OmegaMatter    () &
         &              *self%cosmologyParameters_%densityCritical()
    bbksWDMHalfModeMass=+4.0d0                &
         &              *Pi                   &
         &              /3.0d0                &
         &              *matterDensity        &
         &              *(                    &
         &                +Pi                 &
         &                /wavenumberHalfMode &
         &              )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function bbksWDMHalfModeMass

  double precision function bbksWDMEpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionBBKSWDM), intent(inout) :: self

    bbksWDMEpochTime=self%time
    return
  end function bbksWDMEpochTime
