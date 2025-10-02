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
  Implements a transfer function class based on the non-cold dark matter fitting function of
  \cite{murgia_non-cold_2017}.
  !!}

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <transferFunction name="transferFunctionMurgia2017">
   <description>Provides a transfer function based on the non-cold dark matter fitting function of \cite{murgia_non-cold_2017}.</description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionMurgia2017
     !!{
     A transfer function class based on the non-cold dark matter fitting function of \cite{murgia_non-cold_2017}.
     !!}
     private
     double precision                                    :: alpha                         , beta, &
          &                                                 gamma                         , time, &
          &                                                 redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
   contains
     final     ::                          murgia2017Destructor
     procedure :: value                 => murgia2017Value
     procedure :: logarithmicDerivative => murgia2017LogarithmicDerivative
     procedure :: halfModeMass          => murgia2017HalfModeMass
     procedure :: quarterModeMass       => murgia2017QuarterModeMass
     procedure :: epochTime             => murgia2017EpochTime
  end type transferFunctionMurgia2017

  interface transferFunctionMurgia2017
     !!{
     Constructors for the \refClass{transferFunctionMurgia2017} transfer function class.
     !!}
     module procedure murgia2017ConstructorParameters
     module procedure murgia2017ConstructorInternal
  end interface transferFunctionMurgia2017

contains

function murgia2017ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{transferFunctionMurgia2017} transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionMurgia2017)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (transferFunctionClass     ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    double precision                                            :: alpha               , beta    , &
         &                                                         gamma               , redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunction')) call Error_Report("an explicit 'transferFunction' must be given"//{introspection:location})
    ! Read parameters.
    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <defaultValue>0.0075d0</defaultValue>
      <description>The parameter $\alpha$, which sets the cut-off scale length, appearing in the warm dark matter transfer function \citep{murgia_non-cold_2017}.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>1.5d0</defaultValue>
      <description>The parameter $\beta$, which controls the shape of the cut-off, appearing in the warm dark matter transfer function \citep{murgia_non-cold_2017}.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <defaultValue>-10.0d0</defaultValue>
      <description>The parameter $\gamma$, which controls the shape of the cut-off, appearing in the warm dark matter transfer function \citep{murgia_non-cold_2017}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%equalityEpochMatterRadiation(requestTypeExpansionFactor))</defaultValue>
      <description>The redshift of the epoch at which the transfer function is defined.</description>
    </inputParameter>
    !!]
    self=transferFunctionMurgia2017(transferFunctionCDM,alpha,beta,gamma,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunctionCDM" />
    !!]
    return
  end function murgia2017ConstructorParameters

  function murgia2017ConstructorInternal(transferFunctionCDM,alpha,beta,gamma,time,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{transferFunctionMurgia2017} transfer function class.
    !!}
    implicit none
    type            (transferFunctionMurgia2017)                        :: self
    class           (transferFunctionClass     ), target, intent(in   ) :: transferFunctionCDM
    double precision                                    , intent(in   ) :: alpha               , beta, &
         &                                                                 gamma               , time
    class           (cosmologyParametersClass  ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), target, intent(in   ) :: cosmologyFunctions_
    !![
    <constructorAssign variables="*transferFunctionCDM, alpha, beta, gamma, time, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    return
  end function murgia2017ConstructorInternal

  subroutine murgia2017Destructor(self)
    !!{
    Destructor for the \refClass{transferFunctionMurgia2017} transfer function class.
    !!}
    implicit none
    type(transferFunctionMurgia2017), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%transferFunctionCDM" />
    !!]
    return
  end subroutine murgia2017Destructor

  double precision function murgia2017Value(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionMurgia2017), intent(inout) :: self
    double precision                            , intent(in   ) :: wavenumber

    murgia2017Value=+self%transferFunctionCDM%value(wavenumber)
    if (self%alpha > 0.0d0)                 &
         & murgia2017Value=+murgia2017Value &
         &                 *(               &
         &                   +1.0d0         &
         &                   +(             &
         &                     +self%alpha  &
         &                     *wavenumber  &
         &                    )**self%beta  &
         &                  )  **self%gamma
    return
  end function murgia2017Value

  double precision function murgia2017LogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionMurgia2017), intent(inout) :: self
    double precision                            , intent(in   ) :: wavenumber

    murgia2017LogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    if (self%alpha > 0.0d0)                                                 &
         & murgia2017LogarithmicDerivative=+murgia2017LogarithmicDerivative &
         &                                 +self%beta                       &
         &                                 *self%gamma                      &
         &                                 *(                               &
         &                                   +self%alpha                    &
         &                                   *wavenumber                    &
         &                                  )**self%beta                    &
         &                                 /(                               &
         &                                   +(                             &
         &                                     +1.0d0                       &
         &                                     +(                           &
         &                                       +self%alpha                &
         &                                       *wavenumber                &
         &                                      )**self%beta                &
         &                                    )                             &
         &                                   *wavenumber                    &
         &                                  )
    return
  end function murgia2017LogarithmicDerivative

  double precision function murgia2017HalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionMurgia2017), intent(inout), target   :: self
    integer                                     , intent(  out), optional :: status
    double precision                                                      :: matterDensity, wavenumberHalfMode

    matterDensity         =+self%cosmologyParameters_%OmegaMatter    () &
         &                 *self%cosmologyParameters_%densityCritical()
    wavenumberHalfMode    =+(                         &
         &                   +(                       &
         &                     +1.0d0                 &
         &                     /self%alpha            &
         &                    )                       &
         &                   *(                       &
         &                     +(                     &
         &                       +1.0d0               &
         &                       /2.0d0               &
         &                      )**(1.0d0/self%gamma) &
         &                     -1.0d0                 &
         &                    )**(1.0d0/self%beta)    &
         &                  )
    ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
    ! R = λ/2 = π/k [see Eq.(9) of \cite{schneider_non-linear_2012}].
    murgia2017HalfModeMass=+4.0d0                &
         &                 *Pi                   &
         &                 /3.0d0                &
         &                 *matterDensity        &
         &                 *(                    &
         &                   +Pi                 &
         &                   /wavenumberHalfMode &
         &                 )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function murgia2017HalfModeMass

  double precision function murgia2017QuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionMurgia2017), intent(inout), target   :: self
    integer                                     , intent(  out), optional :: status
    double precision                                                      :: matterDensity, wavenumberQuarterMode

    matterDensity            =+self%cosmologyParameters_%OmegaMatter    () &
         &                    *self%cosmologyParameters_%densityCritical()
    wavenumberQuarterMode    =+(                         &
         &                      +(                       &
         &                        +1.0d0                 &
         &                        /self%alpha            &
         &                       )                       &
         &                      *(                       &
         &                        +(                     &
         &                          +1.0d0               &
         &                          /4.0d0               &
         &                         )**(1.0d0/self%gamma) &
         &                        -1.0d0                 &
         &                       )**(1.0d0/self%beta)    &
         &                     )
    ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
    ! R = λ/2 = π/k [see Eq.(9) of \cite{schneider_non-linear_2012}].
    murgia2017QuarterModeMass=+4.0d0                   &
         &                    *Pi                      &
         &                    /3.0d0                   &
         &                    *matterDensity           &
         &                    *(                       &
         &                      +Pi                    &
         &                      /wavenumberQuarterMode &
         &                    )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function murgia2017QuarterModeMass

  double precision function murgia2017EpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionMurgia2017), intent(inout) :: self

    murgia2017EpochTime=self%time
    return
  end function murgia2017EpochTime
