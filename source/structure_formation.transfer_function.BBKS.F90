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
  Implements the transfer function fitting function of \cite{bardeen_statistics_1986}.
  !!}

  use :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters , only : cosmologyParametersClass
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <transferFunction name="transferFunctionBBKS">
   <description>
    A transfer function class implementing the \cite{bardeen_statistics_1986} fitting function to compute the \gls{cdm}
    transfer function.
   </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionBBKS
     !!{
     A bbks transfer function class.
     !!}
     private
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_ => null()
     double precision                                    :: Gamma                        , time
   contains
     final     ::                          destructor
     procedure :: value                 => bbksValue
     procedure :: logarithmicDerivative => bbksLogarithmicDerivative
     procedure :: halfModeMass          => bbksHalfModeMass
     procedure :: quarterModeMass       => bbksQuarterModeMass
     procedure :: epochTime             => bbksEpochTime
  end type transferFunctionBBKS

  interface transferFunctionBBKS
     !!{
     Constructors for the \refClass{transferFunctionBBKS} transfer function class.
     !!}
     module procedure constructorParameters
     module procedure constructorInternal
  end interface transferFunctionBBKS

  ! Fitting function parameters.
  double precision, parameter :: p=2.34d0, a=3.89d0, b=16.10d0, c=5.46d0, d=6.71d0

contains

  function constructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{transferFunctionBBKS} transfer function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (transferFunctionBBKS    )                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(cosmologyParametersClass), pointer       :: cosmologyParameters_
    class(darkMatterParticleClass ), pointer       :: darkMatterParticle_
    class(cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=constructorInternal(darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function constructorParameters

  function constructorInternal(darkMatterParticle_,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{transferFunctionBBKS} transfer function class.
    !!}
    use :: Cosmology_Parameters , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles, only : darkMatterParticleCDM
    use :: Error                , only : Error_Report
    implicit none
    type (transferFunctionBBKS    )                        :: self
    class(darkMatterParticleClass ), intent(in   ), target :: darkMatterParticle_
    class(cosmologyParametersClass), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*darkMatterParticle_, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Error_Report('transfer function expects a cold dark matter particle'//{introspection:location})
    end select
    ! Compute the epoch - the transfer function is assumed to be for z=0.
    self%time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(0.0d0))
    ! Compute the Gamma parameter.
    self%Gamma=+             self%cosmologyParameters_%OmegaMatter   (                  ) &
         &     *             self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH) &
         &     *exp(                                                                      &
         &          -        self%cosmologyParameters_%OmegaBaryon   (                  ) &
         &          *(                                                                    &
         &            +1.0d0                                                              &
         &            +sqrt(                                                              &
         &                  +2.0d0                                                        &
         &                  *self%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH) &
         &                 )                                                              &
         &            /      self%cosmologyParameters_%OmegaMatter   (                  ) &
         &           )                                                                    &
         &         )                                                                      &
         &     /(                                                                         &
         &       +           self%cosmologyParameters_%temperatureCMB(                  ) &
         &       /2.7d0                                                                   &
         &      )**2
    return
  end function constructorInternal

  subroutine destructor(self)
    !!{
    Destructor for the \refClass{transferFunctionBBKS} transfer function class.
    !!}
    implicit none
    type(transferFunctionBBKS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine destructor

  double precision function bbksValue(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    implicit none
    class           (transferFunctionBBKS), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber
    double precision                                      :: wavenumberHUnits, q

    ! Get wavenumber in "little-h" units.
    wavenumberHUnits=+wavenumber                                                   &
         &           /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    q               =+wavenumberHUnits &
         &           /self%Gamma
    bbksValue       =+(             &
         &             +log(        &
         &                  +1.00d0 &
         &                  +p      &
         &                  *q      &
         &                 )        &
         &             /p           &
         &             /q           &
         &            )             &
         &           /(             &
         &             +1.0d0       &
         &             +(a*q)       &
         &             +(b*q)**2    &
         &             +(c*q)**3    &
         &             +(d*q)**4    &
         &            )**0.25d0
    return
  end function bbksValue

  double precision function bbksLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    implicit none
    class           (transferFunctionBBKS), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber
    double precision                                      :: wavenumberHUnits, q

    ! Get wavenumber in "little-h" units.
    wavenumberHUnits=+wavenumber                                                   &
         &           /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    q               =+wavenumberHUnits &
         &           /self%Gamma
    bbksLogarithmicDerivative=-2.0d0            &
         &                    +(                &
         &                      +4.0d0          &
         &                      +3.0d0*(a*q)    &
         &                      +2.0d0*(b*q)**2 &
         &                      +      (c*q)**3 &
         &                    )                 &
         &                    /(                &
         &                      +4.0d0          &
         &                      *(              &
         &                        +1.0d0        &
         &                        +      q      &
         &                        *(            &
         &                          +    a      &
         &                          +    q      &
         &                          *(          &
         &                            +  b**2   &
         &                            +  q      &
         &                            *(        &
         &                              +c**3   &
         &                              +d**4   &
         &                              *q      &
         &                             )        &
         &                           )          &
         &                         )            &
         &                       )              &
         &                     )                &
         &                    +(                &
         &                      +p              &
         &                      *q              &
         &                     )                &
         &                    /(                &
         &                      +(              &
         &                        +1.0d0        &
         &                        +p            &
         &                        *q            &
         &                       )              &
         &                      *log(           &
         &                           +1.0d0     &
         &                           +p         &
         &                           *q         &
         &                          )           &
         &                     )
    return
  end function bbksLogarithmicDerivative

  double precision function bbksHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function. Not supported in this implementation.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionBBKS), intent(inout), target   :: self
    integer                      , intent(  out), optional :: status
    !$GLC attributes unused :: self

    bbksHalfModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function bbksHalfModeMass

  double precision function bbksQuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative
    to a \gls{cdm} transfer function. Not supported in this implementation.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionBBKS), intent(inout), target   :: self
    integer                      , intent(  out), optional :: status
    !$GLC attributes unused :: self

    bbksQuarterModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function bbksQuarterModeMass

  double precision function bbksEpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionBBKS), intent(inout) :: self

    bbksEpochTime=self%time
    return
  end function bbksEpochTime
