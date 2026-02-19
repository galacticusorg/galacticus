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
Implements the ETHOS \citep{cyr-racine_ethoseffective_2016} transfer function, using the specific form given by
\cite[][eqn.~3]{bohr_halo_2021}.
!!}

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <transferFunction name="transferFunctionETHOSDM">
    <description>
      Implements the ETHOS \citep{cyr-racine_ethoseffective_2016} transfer function, using the specific form given by
      \cite[][eqn.~3]{bohr_halo_2021}.
    </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionETHOSDM
     !!{
     Implements the ETHOS \citep{cyr-racine_ethoseffective_2016} transfer function, using the specific form given by
     \cite[][eqn.~3]{bohr_halo_2021}.
     !!}
     private
     double precision                                    :: alpha                         , beta    , &
          &                                                 gamma                         , sigma   , &
          &                                                 tau                           , kPeak   , &
          &                                                 hPeak                         , h2      , &
          &                                                 time                          , redshift
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM  => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
   contains
     final     ::                          ETHOSDMDestructor
     procedure :: value                 => ETHOSDMValue
     procedure :: logarithmicDerivative => ETHOSDMLogarithmicDerivative
     procedure :: halfModeMass          => ETHOSDMHalfModeMass
     procedure :: quarterModeMass       => ETHOSDMQuarterModeMass
     procedure :: fractionModeMass      => ETHOSDMFractionModeMass
     procedure :: epochTime             => ETHOSDMEpochTime
  end type transferFunctionETHOSDM

  interface transferFunctionETHOSDM
     !!{
     Constructors for the \refClass{transferFunctionETHOSDM} transfer function class.
     !!}
     module procedure ETHOSDMConstructorParameters
     module procedure ETHOSDMConstructorInternal
  end interface transferFunctionETHOSDM

  ! Submodule-scope variables used in root finding.
  class           (transferFunctionETHOSDM), pointer :: self_
  double precision                                   :: modeTarget
  !$omp threadprivate(self_,modeTarget)

contains
  
  function ETHOSDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{transferFunctionETHOSDM} transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Error                         , only : Error_Report
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionETHOSDM )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    double precision                                          :: alpha               , beta , &
         &                                                       gamma               , sigma, & 
         &                                                       tau                 , kPeak, &
         &                                                       hPeak               , h2   , &
         &                                                       redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunction')) call Error_Report("an explicit 'transferFunction' must be given"//{introspection:location})
    ! Read parameters.
    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <defaultValue>40.0d0</defaultValue>
      <description>The parameter $\alpha$ appearing in the ETHOS transfer function \citep{cyr-racine_ethoseffective_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>1.5d0</defaultValue>
      <description>The parameter $\beta$ appearing in the ETHOS transfer function \citep{cyr-racine_ethoseffective_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <defaultValue>-10.0d0</defaultValue>
      <description>The parameter $\gamma$ appearing in the ETHOS transfer function \citep{cyr-racine_ethoseffective_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>sigma</name>
      <source>parameters</source>
      <defaultValue>-10.0d0</defaultValue>
      <description>The parameter $\sigma$ appearing in the ETHOS transfer function \citep{cyr-racine_ethoseffective_2016}, determines width of first peak in transfer function.</description>
    </inputParameter>
    <inputParameter>
      <name>tau</name>
      <source>parameters</source>
      <defaultValue>-10.0d0</defaultValue>
      <description>The parameter $\tau$ appearing in the ETHOS transfer function \citep{cyr-racine_ethoseffective_2016}, determines damping of DAO.</description>
    </inputParameter>
    <inputParameter>
      <name>kPeak</name>
      <source>parameters</source>
      <defaultValue>-10.0d0</defaultValue>
      <description>The parameter $k_\mathrm{peak}$, the wavenumber of first peak in ETHOS transfer function, appearing in the ETHOS transfer function \citep{cyr-racine_ethoseffective_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>hPeak</name>
      <source>parameters</source>
      <defaultValue>-10.0d0</defaultValue>
      <description>The parameter $h_\mathrm{peak}$, the amplitude of the first peak, appearing in the ETHOS transfer function \citep{cyr-racine_ethoseffective_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>h2</name>
      <source>parameters</source>
      <defaultValue>-10.0d0</defaultValue>
      <description>The parameter $h_2$, the amplitude of the second peak, appearing in the ETHOS transfer function \citep{cyr-racine_ethoseffective_2016}.</description>
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
    self=transferFunctionETHOSDM(transferFunctionCDM,alpha,beta,gamma,sigma,tau,kPeak,hPeak,h2,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunctionCDM" />
    !!]
    return
  end function ETHOSDMConstructorParameters
  
  function ETHOSDMConstructorInternal(transferFunctionCDM,alpha,beta,gamma,sigma,tau,kPeak,hPeak,h2,time,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{transferFunctionETHOSDM} transfer function class.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Error               , only : Error_Report
    implicit none
    type            (transferFunctionETHOSDM )                        :: self
    class           (transferFunctionClass   ), target, intent(in   ) :: transferFunctionCDM
    double precision                                  , intent(in   ) :: alpha              , beta , &
         &                                                               gamma              , sigma, &
         &                                                               tau                , kPeak, &
         &                                                               hPeak              , h2   , &
         &                                                               time
    class           (cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_
    !![
    <constructorAssign variables="*transferFunctionCDM, alpha, beta, gamma,sigma, tau, kPeak, hPeak, h2, time, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]
    
    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    return
  end function ETHOSDMConstructorInternal

  subroutine ETHOSDMDestructor(self)
    !!{
    Destructor for the \refClass{transferFunctionETHOSDM} transfer function class.
    !!}
    implicit none
    type(transferFunctionETHOSDM), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%transferFunctionCDM" />
    !!]
    return
  end subroutine ETHOSDMDestructor

  double precision function ETHOSDMValue(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber, using the specific form given by \cite[][eqn.~3]{bohr_halo_2021}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionETHOSDM), intent(inout) :: self
    double precision                         , intent(in   ) :: wavenumber

    ETHOSDMValue=+self%transferFunctionCDM%value(wavenumber)
    if (self%alpha > 0.0d0)                    &
         & ETHOSDMValue=+ETHOSDMValue          &
         &              *(                     &
         &                +(                   &
         &                  +1.0d0             &
         &                  +(                 &
         &                    +self%alpha      &
         &                    *wavenumber      &
         &                   )**self%beta      &
         &                 )  **self%gamma     &
         &                -sqrt(self%hPeak)    &
         &                *exp(                &
         &                     -0.5d0          &
         &                     *(              &
         &                       +(            &
         &                         +wavenumber &
         &                         -self%kPeak &
         &                        )            &
         &                       /(            &
         &                         +self%sigma &
         &                         *self%kPeak &
         &                        )            &
         &                      )**2           &
         &                    )                &
         &                +0.25d0              &
         &                *sqrt(self%h2)       &
         &                *erfc(               &
         &                      +(             &
         &                        +wavenumber  &
         &                        -1.805d0     &
         &                        *self%kPeak  &
         &                       )             &
         &                      /(             &
         &                        +self%tau    &
         &                        *self%kPeak  &
         &                       )             &
         &                      -2.0d0         &
         &                     )               &
         &                *erf(                &
         &                     -(              &
         &                       +wavenumber   &
         &                       -1.805d0      &
         &                       *self%kPeak   &
         &                      )              &
         &                     /(              &
         &                       +self%sigma   &
         &                       *self%kPeak   &
         &                      )              &
         &                     -2.0d0          &
         &                    )                &                    
         &                *cos(                &
         &                     +1.1083d0       &
         &                     *Pi             &
         &                     *wavenumber     &
         &                     /self%kPeak     &
         &                    )                &
         &               )
    return
  end function ETHOSDMValue

  double precision function ETHOSDMLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber, using the specific form given by \cite[][eqn.~3]{bohr_halo_2021}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (transferFunctionETHOSDM), intent(inout) :: self
    double precision                         , intent(in   ) :: wavenumber

    ETHOSDMLogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)
    if (self%alpha > 0.0d0)                                                                   &
         & ETHOSDMLogarithmicDerivative=+ETHOSDMLogarithmicDerivative                         &
         &                              +(                                                    &
         &                                +self%alpha                                         &
         &                                *(                                                  &
         &                                  +wavenumber                                       &
         &                                  *self%alpha                                       &
         &                                 )**(-1.0d0+self%beta)                              &
         &                                *(                                                  &
         &                                  +1.0d0                                            &
         &                                  +(                                                &
         &                                    +wavenumber                                     &
         &                                    *self%alpha                                     &
         &                                   )**self%beta                                     &
         &                                 )**(-1.0d0+self%gamma)                             &
         &                                *self%beta                                          &
         &                                *self%gamma                                         &
         &                                +(                                                  &
         &                                  +sqrt(self%hPeak)                                 &
         &                                  *(                                                &
         &                                    +wavenumber                                     &
         &                                    -self%kPeak                                     &
         &                                   )                                                &
         &                                 )                                                  &
         &                                /(                                                  &
         &                                  +exp(                                             &
         &                                       +(                                           &
         &                                         +wavenumber                                &
         &                                         -self%kPeak                                &
         &                                        )**2                                        &
         &                                       /(                                           &
         &                                         +2.0d0                                     &
         &                                         *self%kPeak**2                             &
         &                                         *self%sigma**2                             &
         &                                        )                                           &
         &                                      )                                             &
         &                                  *self%kPeak**2                                    &
         &                                  *self%sigma**2                                    &
         &                                 )                                                  &
         &                                -(                                                  &
         &                                  +sqrt(self%h2)                                    &
         &                                  *cos (                                            &
         &                                        +(                                          &
         &                                          +3.4818271379735677d0                     &
         &                                          *wavenumber                               &
         &                                         )                                          &
         &                                        /self%kPeak                                 &
         &                                       )                                            &
         &                                  *erfc(                                            &
         &                                        -2.0d0                                      &
         &                                        +(                                          &
         &                                          -wavenumber                               &
         &                                          +1.805d0                                  &
         &                                          *self%kPeak                               &
         &                                         )                                          &
         &                                        /(                                          &
         &                                          +self%kPeak                               &
         &                                          *self%sigma                               &
         &                                         )                                          &
         &                                       )                                            &
         &                                 )                                                  &
         &                                /(                                                  &
         &                                  +2.0d0                                            &
         &                                  *exp(                                             &
         &                                       -2.0d0                                       &
         &                                       +(                                           &
         &                                         +wavenumber                                &
         &                                         -1.805d0                                   &
         &                                         *self%kPeak                                &
         &                                        )                                           &
         &                                       /(                                           &
         &                                         +self%kPeak                                &
         &                                         *self%tau                                  &
         &                                        )                                           &
         &                                      )**2                                          &
         &                                  *self%kPeak                                       &
         &                                  *self%tau                                         &
         &                                  *sqrt(Pi)                                         &
         &                                 )                                                  &
         &                                +(                                                  &
         &                                  +sqrt(self%h2)                                    &
         &                                  *cos (                                            &
         &                                        +(                                          &
         &                                          +3.4818271379735677d0                     &
         &                                          *wavenumber                               &
         &                                         )                                          &
         &                                        /self%kPeak                                 &
         &                                       )                                            &
         &                                  *erfc(                                            &
         &                                        -2.0d0                                      &
         &                                        +(                                          &
         &                                          +wavenumber                               &
         &                                          -1.805d0                                  &
         &                                          *self%kPeak                               &
         &                                         )                                          &
         &                                        /(                                          &
         &                                          +self%kPeak                               &
         &                                          *self%tau                                 &
         &                                         )                                          &
         &                                       )                                            &
         &                                 )                                                  &
         &                                /(                                                  &
         &                                  +2.0d0                                            &
         &                                  *exp(                                             &
         &                                       -2.0d0                                       &
         &                                       +(                                           &
         &                                         -wavenumber                                &
         &                                         +1.805d0                                   &
         &                                         *self%kPeak                                &
         &                                        )                                           &
         &                                       /(                                           &
         &                                         +self%kPeak                                &
         &                                         *self%sigma                                &
         &                                        )                                           &
         &                                      )**2                                          &
         &                                  *self%kPeak                                       &
         &                                  *self%sigma                                       &
         &                                  *sqrt(Pi)                                         &
         &                                 )                                                  &
         &                                -(                                                  &
         &                                  +0.8704567844933919d0                             &
         &                                  *sqrt(self%h2)                                    &
         &                                  *erfc(                                            &
         &                                        -2.0d0                                      &
         &                                        +(                                          &
         &                                          -wavenumber                               &
         &                                          +1.805d0                                  &
         &                                          *self%kPeak                               &
         &                                         )                                          &
         &                                        /(                                          &
         &                                          +self%kPeak                               &
         &                                          *self%sigma                               &
         &                                         )                                          &
         &                                       )                                            &
         &                                  *erfc(                                            &
         &                                        -2.0d0                                      &
         &                                        +(                                          &
         &                                          +wavenumber                               &
         &                                          -1.805d0                                  &
         &                                          *self%kPeak                               &
         &                                         )                                          &
         &                                        /(                                          &
         &                                          +self%kPeak                               &
         &                                          *self%tau                                 &
         &                                         )                                          &
         &                                       )                                            &
         &                                  *sin (                                            &
         &                                        +(                                          &
         &                                          +3.4818271379735677d0                     &
         &                                          *wavenumber                               &
         &                                         )                                          &
         &                                        /self%kPeak                                 &
         &                                       )                                            &
         &                                 )                                                  &
         &                                /self%kPeak                                         &
         &                               )                                                    &
         &                              /(                                                    &
         &                                -(                                                  &
         &                                  +sqrt(self%hPeak)                                 &
         &                                  /exp (                                            &
         &                                        +(                                          &
         &                                          +wavenumber                               &
         &                                          -self%kPeak                               &
         &                                         )**2                                       &
         &                                        /(                                          &
         &                                          +2.0d0                                    &
         &                                          *self%kPeak**2                            &
         &                                          *self%sigma**2                            &
         &                                         )                                          &
         &                                       )                                            &
         &                                 )                                                  &
         &                                +(                                                  &
         &                                  +1.0d0                                            &
         &                                  +(                                                &
         &                                    +wavenumber                                     &
         &                                    *self%alpha                                     &
         &                                   )**self%beta                                     &
         &                                 )**self%gamma                                      &
         &                                +(                                                  &
         &                                  +sqrt(self%h2)                                    &
         &                                  *cos (                                            &
         &                                        +(                                          &
         &                                          +3.4818271379735677                       &
         &                                          *wavenumber                               &
         &                                         )                                          &
         &                                        /self%kPeak                                 &
         &                                       )                                            &
         &                                  *erfc(                                            &
         &                                        -2.0d0                                      &
         &                                        +(                                          &
         &                                          -wavenumber                               &
         &                                          +1.805d0                                  &
         &                                          *self%kPeak                               &
         &                                         )                                          &
         &                                        /(                                          &
         &                                          +self%kPeak                               &
         &                                          *self%sigma                               &
         &                                         )                                          &
         &                                       )                                            &
         &                                  *erfc(                                            &
         &                                        -2.0d0                                      &
         &                                        +(                                          &
         &                                          +wavenumber                               &
         &                                          -1.805d0                                  &
         &                                          *self%kPeak                               &
         &                                         )                                          &
         &                                        /(                                          &
         &                                          +self%kPeak                               &
         &                                          *self%tau                                 &
         &                                         )                                          &
         &                                       )                                            &
         &                                 )                                                  &
         &                                /4.0d0                                              &
         &                               )
    return
  end function ETHOSDMLogarithmicDerivative

  double precision function ETHOSDMHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    implicit none
    class  (transferFunctionETHOSDM), intent(inout), target   :: self
    integer                         , intent(  out), optional :: status

    ETHOSDMHalfModeMass=self%fractionModeMass(0.50d0,status)
    return
  end function ETHOSDMHalfModeMass

  double precision function ETHOSDMQuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative
    to a \gls{cdm} transfer function.
    !!}
    implicit none
    class  (transferFunctionETHOSDM), intent(inout), target   :: self
    integer                         , intent(  out), optional :: status

    ETHOSDMQuarterModeMass=self%fractionModeMass(0.25d0,status)
    return
  end function ETHOSDMQuarterModeMass

  double precision function ETHOSDMFractionModeMass(self,fraction,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (transferFunctionETHOSDM), intent(inout), target   :: self
    double precision                         , intent(in   )           :: fraction
    integer                                  , intent(  out), optional :: status
    double precision                                                   :: matterDensity, wavenumberFractionMode
    type            (rootFinder             )                          :: finder

    ! There is no analytic solution for the fraction-mode mass so we resort to numerical root finding. This is complicated by the fact
    ! that the transfer function oscillates. Our approach is to start at a wavenumber much smaller than the cut-off scale, 1/α,
    ! and slowly increase the wavenumber until the root is bracketed.
    finder                  =   rootFinder(                                                             &
         &                                 rootFunction                 =modeSolver                   , &
         &                                 toleranceRelative            =1.000d-3                     , &
         &                                 rangeExpandUpward            =1.001d+0                     , &
         &                                 rangeExpandDownward          =0.500d+0                     , &
         &                                 rangeExpandType              =rangeExpandMultiplicative    , &
         &                                 rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                                 rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
         &                                )
    self_                   =>  self
    modeTarget              =  +                         fraction        &
         &                     *self                    %value   (0.0d0) &
         &                     /self%transferFunctionCDM%value   (0.0d0)
    wavenumberFractionMode  =   finder%find(rootGuess=1.0d-2/self%alpha)
    matterDensity           =  +self%cosmologyParameters_%OmegaMatter    () &
         &                     *self%cosmologyParameters_%densityCritical()
    ETHOSDMFractionModeMass =  +4.0d0                    &
         &                     *Pi                       &
         &                     /3.0d0                    &
         &                     *matterDensity            &
         &                     *(                        &
         &                       +Pi                     &
         &                       /wavenumberFractionMode &
         &                      )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function ETHOSDMFractionModeMass
  
  double precision function modeSolver(wavenumber)
    !!{
    Function used in solving for half- and quarter-mode masses in the ETHOS transfer function.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber
  
    modeSolver=+self_                    %value     (wavenumber) &
         &     /self_%transferFunctionCDM%value     (wavenumber) &
         &     -                           modeTarget
    return
  end function modeSolver
  
  double precision function ETHOSDMEpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionETHOSDM), intent(inout) :: self

    ETHOSDMEpochTime=self%time
    return
  end function ETHOSDMEpochTime
