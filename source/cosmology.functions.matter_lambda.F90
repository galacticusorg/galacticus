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
  An implementation of the cosmological functions class for cosmologies consisting of collisionless
  matter plus a cosmological constant.
  !!}

  use :: Cosmology_Parameters   , only : cosmologyParameters, cosmologyParametersClass
  use :: Numerical_Interpolation, only : interpolator

  integer         , parameter :: ageTableNPointsPerDecade     =300
  double precision, parameter :: ageTableNPointsPerOctave     =dble(ageTableNPointsPerDecade)*log(2.0d0)/log(10.0d0)
  double precision, parameter :: ageTableIncrementFactor      =exp(int(ageTableNPointsPerOctave+1.0d0)*log(10.0d0)/dble(ageTableNPointsPerDecade))
  integer         , parameter :: distanceTableNPointsPerDecade=100

  ! Factor by which one component of Universe must dominate others such that we can ignore the others.
  double precision, parameter :: matterLambdaDominateFactor               =100.0d0

  !![
  <cosmologyFunctions name="cosmologyFunctionsMatterLambda">
   <description>
    Cosmological relations are computed assuming a universe that contains only collisionless matter and a cosmological
    constant.
   </description>
  </cosmologyFunctions>
  !!]
  type, extends(cosmologyFunctionsClass) :: cosmologyFunctionsMatterLambda
     !!{
     A cosmological functions class for cosmologies consisting of matter plus a cosmological constant.
     !!}
     private
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_                            => null()
     logical                                                               :: collapsingUniverse                              =.false.
     integer                                                               :: iTableTurnaround
     double precision                                                      :: expansionFactorMaximum                                   , expansionFactorPrevious                             =-1.0d0, &
          &                                                                   timeMaximum                                              , timePrevious                                        =-1.0d0, &
          &                                                                   timeTurnaround
     double precision                                       , dimension(2) :: expansionRatePrevious                                    , expansionRateExpansionFactorPrevious
     logical                                                               :: ageTableInitialized                             =.false.
     integer                                                               :: ageTableNumberPoints
     double precision                                                      :: ageTableTimeMaximum                             =20.0d0  , ageTableTimeMinimum                                 =1.0d-4
     double precision                                                      :: ageTableTimeLogarithmicMinimum                           , ageTableInverseDeltaLogTime
     double precision                          , allocatable, dimension(:) :: ageTableExpansionFactor                                  , ageTableTime
     type            (interpolator            ), allocatable               :: interpolatorTime                                         , interpolatorDistance                                       , &
          &                                                                   interpolatorDistanceInverse                              , interpolatorLuminosityDistance                             , &
          &                                                                   interpolatorLuminosityDistanceKCorrected
     logical                                                               :: distanceTableInitialized                        =.false.
     integer                                                               :: distanceTableNumberPoints
     double precision                                                      :: distanceTableTimeMaximum                                 , distanceTableTimeMinimum                            =1.0d+0
     double precision                          , allocatable, dimension(:) :: distanceTableComovingDistance                            , distanceTableComovingDistanceNegated                       , &
          &                                                                   distanceTableLuminosityDistanceNegated                   , distanceTableTime                                          , &
          &                                                                   distanceTableLuminosityDistanceKCorrectedNegated
     logical                                                               :: enableRangeChecks
   contains
     !![
     <methods>
       <method description="Tabulate comoving distance as a function of cosmic time." method="distanceTabulate" />
       <method description="Tabulate expansion factor as a function of cosmic time." method="expansionFactorTabulate" />
     </methods>
     !!]
     final     ::                                    matterLambdaDestructor
     procedure :: epochValidate                   => matterLambdaEpochValidate
     procedure :: cosmicTime                      => matterLambdaCosmicTime
     procedure :: timeBigCrunch                   => matterLambdaTimeBigCrunch
     procedure :: expansionFactor                 => matterLambdaExpansionFactor
     procedure :: expansionRate                   => matterLambdaExpansionRate
     procedure :: hubbleParameterEpochal          => matterLambdaHubbleParameterEpochal
     procedure :: hubbleParameterRateOfChange     => matterLambdaHubbleParameterRateOfChange
     procedure :: densityScalingEarlyTime         => matterLambdaDensityScalingEarlyTime
     procedure :: omegaMatterEpochal              => matterLambdaOmegaMatterEpochal
     procedure :: omegaMatterRateOfChange         => matterLambdaOmegaMatterRateOfChange
     procedure :: omegaDarkEnergyEpochal          => matterLambdaOmegaDarkEnergyEpochal
     procedure :: equationOfStateDarkEnergy       => matterLambdaEquationOfStateDarkEnergy
     procedure :: exponentDarkEnergy              => matterLambdaExponentDarkEnergy
     procedure :: equalityEpochMatterDarkEnergy   => matterLambdaEqualityEpochMatterDarkEnergy
     procedure :: equalityEpochMatterCurvature    => matterLambdaEqualityEpochMatterCurvature
     procedure :: equalityEpochMatterRadiation    => matterLambdaEqualityEpochMatterRadiation
     procedure :: dominationEpochMatter           => matterLambdaDominationEpochMatter
     procedure :: temperatureCMBEpochal           => matterLambdaTemperatureCMBEpochal
     procedure :: distanceComoving                => matterLambdaDistanceComoving
     procedure :: distanceLuminosity              => matterLambdaDistanceLuminosity
     procedure :: distanceAngular                 => matterLambdaDistanceAngular
     procedure :: timeAtDistanceComoving          => matterLambdaTimeAtDistanceComoving
     procedure :: distanceComovingConvert         => matterLambdaDistanceComovingConvert
     procedure :: expansionFactorTabulate         => matterLambdaMakeExpansionFactorTable
     procedure :: distanceTabulate                => matterLambdaMakeDistanceTable
     procedure :: matterDensityEpochal            => matterLambdaMatterDensityEpochal
     procedure :: distanceParticleHorizonComoving => matterLambdaDistanceParticleHorizonComoving
  end type cosmologyFunctionsMatterLambda

  ! Module scope pointer to the current object.
  class(cosmologyFunctionsMatterLambda), pointer :: self_
  !$omp threadprivate(self_)
  !$GLC ignore outlive :: self_

  interface cosmologyFunctionsMatterLambda
     !!{
     Constructors for the matter plus cosmological constant cosmological functions class.
     !!}
     module procedure matterLambdaConstructorParameters
     module procedure matterLambdaConstructorInternal
  end interface cosmologyFunctionsMatterLambda

contains

  function matterLambdaConstructorParameters(parameters) result(self)
    !!{
    Parameter-based constructor for the matter plus cosmological constant cosmological functions class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (cosmologyFunctionsMatterLambda)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(cosmologyParametersClass      ), pointer       :: cosmologyParameters_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=cosmologyFunctionsMatterLambda(cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function matterLambdaConstructorParameters

  function matterLambdaConstructorInternal(cosmologyParameters_) result(self)
    !!{
    Constructor for the matter plus cosmological constant cosmological functions class.
    !!}
    use :: Cosmology_Parameters , only : hubbleUnitsTime
    use :: Numerical_Comparison , only : Values_Agree
    use :: Numerical_ODE_Solvers, only : odeSolver
    implicit none
    type            (cosmologyFunctionsMatterLambda)               , target :: self
    class           (cosmologyParametersClass      ), intent(in   ), target :: cosmologyParameters_
    double precision                                , parameter             :: odeToleranceAbsolute  =1.0d-09, odeToleranceRelative   =1.0d-09
    double precision                                , parameter             :: omegaTolerance        =1.0d-09
    double complex                                  , dimension(3)          :: expansionFactorMaximum
    double precision                                , dimension(1)          :: timeMaximum
    double precision                                , parameter             :: toleranceRelative     =1.0d-10
    double precision                                                        :: OmegaDominant                 , expansionFactorDominant        , &
         &                                                                     densityPower
    type            (odeSolver                     )                        :: solver
    integer                                                                 :: i
    double complex                                                          :: omegaMatter                   , omegaDarkEnergy                , &
         &                                                                     omegaCurvature                , rootTerm
    !![
    <constructorAssign variables="*cosmologyParameters_"/>
    !!]

    ! Determine if this universe will collapse. We take the Friedmann equation, which gives H²(a) as a function of expansion
    ! factor, a, and solve for where H²(a)=0. If this has a real solution, then we have a collapsing universe.
    self%collapsingUniverse                  =.false.
    self%enableRangeChecks                   =.true.
    self%expansionFactorMaximum              =0.0d0
    self%timeTurnaround                      =0.0d0
    self%timeMaximum                         =0.0d0
    self%expansionRatePrevious               =-1.0d0
    self%expansionRateExpansionFactorPrevious=-1.0d0
    if    (Values_Agree(self%cosmologyParameters_%OmegaCurvature (),0.0d0,absTol=omegaTolerance)) then
       if (Values_Agree(self%cosmologyParameters_%OmegaDarkEnergy(),0.0d0,absTol=omegaTolerance)) then
          ! Einstein-de Sitter case. Always expands to infinity.
          self%collapsingUniverse=.false.
       else
          ! Flat Universe with cosmological constant.
          if (self%cosmologyParameters_%OmegaDarkEnergy() > 0.0d0) then
             ! Never collapses.
             self%collapsingUniverse=.false.
          else
             self%collapsingUniverse=.true.
             self%expansionFactorMaximum                               &
                  & =-(                                                &
                  &     self%cosmologyParameters_%OmegaMatter    ()    &
                  &    *self%cosmologyParameters_%OmegaDarkEnergy()**2 &
                  &   )**(1.0d0/3.0d0)                                 &
                  &   / self%cosmologyParameters_%OmegaDarkEnergy()
          end if
       end if
    else if (Values_Agree(self%cosmologyParameters_%OmegaDarkEnergy(),0.0d0,absTol=omegaTolerance)) then
       ! Non-flat universe with zero cosmological constant.
       if (self%cosmologyParameters_%OmegaMatter() > 1.0d0) then
          ! Universe is closed.
          self%collapsingUniverse    =.true.
          self%expansionFactorMaximum=-self%cosmologyParameters_%OmegaMatter   () &
               &                      /self%cosmologyParameters_%OmegaCurvature()
       else
          ! Universe is open.
          self%collapsingUniverse    =.false.
       end if
    else
       ! Non-flat universe. Seeking roots of the Friedmann equation - this is a cubic equation was we have three roots. We seek
       ! real, positive roots. If multiple such roots exist we take the one with smallest maximum expansion factor.
       omegaCurvature           = dcmplx(self%cosmologyParameters_%OmegaCurvature (),0.0d0)
       omegaMatter              = dcmplx(self%cosmologyParameters_%OmegaMatter    (),0.0d0)
       omegaDarkEnergy          = dcmplx(self%cosmologyParameters_%OmegaDarkEnergy(),0.0d0)
       rootTerm                 =+sqrt(                    &
            &                          + 4.0d0             &
            &                          *omegaCurvature **3 &
            &                          *omegaDarkEnergy**3 &
            &                          +27.0d0             &
            &                          *omegaMatter    **2 &
            &                          *omegaDarkEnergy**4 &
            &                         )
       expansionFactorMaximum(1)=-(2.0d0/3.0d0)    **(1.0d0/3.0d0) &
            &                    *  omegaCurvature                 &
            &                    /(                                &
            &                      -9.0d0                          &
            &                      *omegaMatter                    &
            &                      *omegaDarkEnergy** 2            &
            &                      +sqrt(3.0d0)                    &
            &                      *rootTerm                       &
            &                     )                **(1.0d0/3.0d0) &
            &                    +(                                &
            &                      -9.0d0                          &
            &                      *omegaMatter                    &
            &                      *omegaDarkEnergy** 2            &
            &                      +sqrt(3.0d0)                    &
            &                      *rootTerm                       &
            &                     )**(1.0d0/3.0d0)                 &
            &                    /2.0d0            **(1.0d0/3.0d0) &
            &                    /3.0d0            **(2.0d0/3.0d0) &
            &                    /  omegaDarkEnergy
       expansionFactorMaximum(2)=+dcmplx(1.0d0,sqrt(3.0d0))        &
            &                    *omegaCurvature                   &
            &                    /2.0d0            **(2.0d0/3.0d0) &
            &                    /3.0d0            **(1.0d0/3.0d0) &
            &                    /(                                &
            &                      -9.0d0                          &
            &                      *omegaMatter                    &
            &                      *omegaDarkEnergy** 2            &
            &                      +sqrt(3.0d0)                    &
            &                      *rootTerm                       &
            &                     )                **(1.0d0/3.0d0) &
            &                    -dcmplx(1.0d0,-sqrt(3.0d0))       &
            &                    *(                                &
            &                      -9.0d0                          &
            &                      *omegaMatter                    &
            &                      *omegaDarkEnergy** 2            &
            &                      +sqrt(3.0d0)                    &
            &                      *rootTerm                       &
            &                     )                **(1.0d0/3.0d0) &
            &                    /2.0d0            **(4.0d0/3.0d0) &
            &                    /3.0d0            **(2.0d0/3.0d0) &
            &                    /omegaDarkEnergy
       expansionFactorMaximum(3)=+dcmplx(1.0d0,-sqrt(3.0d0))       &
            &                    *omegaCurvature                   &
            &                    /2.0d0            **(2.0d0/3.0d0) &
            &                    /3.0d0            **(1.0d0/3.0d0) &
            &                    /(                                &
            &                      -9.0d0                          &
            &                      *omegaMatter                    &
            &                      *omegaDarkEnergy** 2            &
            &                      +sqrt(3.0d0)                    &
            &                      *rootTerm                       &
            &                     )                **(1.0d0/3.0d0) &
            &                    -dcmplx(1.0d0,+sqrt(3.0d0))       &
            &                    *(                                &
            &                      -9.0d0                          &
            &                      *omegaMatter                    &
            &                      *omegaDarkEnergy** 2            &
            &                      +sqrt(3.0d0)                    &
            &                      *rootTerm                       &
            &                     )                **(1.0d0/3.0d0) &
            &                    /2.0d0            **(4.0d0/3.0d0) &
            &                    /3.0d0            **(2.0d0/3.0d0) &
            &                    /omegaDarkEnergy
       do i=1,3
          if (real(expansionFactorMaximum(i)) > 0.0d0 .and. abs(imag(expansionFactorMaximum(i))) < toleranceRelative*real(expansionFactorMaximum(i))) then
             if (self%collapsingUniverse) then
                self%expansionFactorMaximum=min(self%expansionFactorMaximum,real(expansionFactorMaximum(i)))
             else
                self%collapsingUniverse    =.true.
                self%expansionFactorMaximum=                                real(expansionFactorMaximum(i))
             end if
          end if
       end do
    end if
    ! If we have a collapsing Universe, find time of turnaround, and maximum time.
    if (self%collapsingUniverse) then
       ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
       call self%densityScalingEarlyTime(matterLambdaDominateFactor,densityPower,expansionFactorDominant,OmegaDominant)
       ! Find the corresponding time. Note that we use the absolute value of the Hubble parameter here - in cases where the
       ! universe is collapsing at the present epoch we need to know the expansion rate (i.e. Hubble parameter) at the equivalent
       ! expansion factor during the expansion phase.
       timeMaximum(1)=1.0d0/abs(self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime))/sqrt(OmegaDominant)/expansionFactorDominant**(0.5d0*densityPower)
       ! Solve Friedmann equation to get time at turnaround.
       self_  => self
       solver =  odeSolver(1_c_size_t,matterLambdaCollapseODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
       call solver%solve(expansionFactorDominant,self%expansionFactorMaximum*(1.0d0-1.0d-4),timeMaximum)
       ! Extract turnaround time from ODE variables and set maximum time to twice turnaround time.
       self%timeTurnaround=timeMaximum(1)
       self%timeMaximum   =2.0d0*self%timeTurnaround
    end if
    return
  end function matterLambdaConstructorInternal

  subroutine matterLambdaDestructor(self)
    !!{
    Default constructor for the matter plus cosmological constant cosmological functions class.
    !!}
    implicit none
    type(cosmologyFunctionsMatterLambda), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine matterLambdaDestructor

  subroutine matterLambdaEpochValidate(self,timeIn,expansionFactorIn,collapsingIn,timeOut,expansionFactorOut,collapsingOut)
    !!{
    Validate a cosmic epoch, specified either by time or expansion factor, and optionally return time, expansion factor, and
    collapsing status.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactorIn , timeIn
    logical                                         , intent(in   ), optional :: collapsingIn
    double precision                                , intent(  out), optional :: expansionFactorOut, timeOut
    logical                                         , intent(  out), optional :: collapsingOut
    logical                                                                   :: collapsingActual

    ! Check that we have a uniquely specified epoch.
    if (      present(timeIn).and.present(expansionFactorIn) ) &
         & call Error_Report('either "time" or "expansionFactor" should be specified, not both'               //{introspection:location})
    if (.not.(present(timeIn).or. present(expansionFactorIn))) &
         & call Error_Report('one of "time" or "expansionFactor" should be specified'                         //{introspection:location})
    if (      present(timeIn).and.present(collapsingIn     ) ) &
         & call Error_Report('collapsing status of universe cannot be specified when epoch is defined by time'//{introspection:location})
    ! If we have a time, check that it is a valid, and compute outputs as required.
    if (present(timeIn)) then
       ! Validate.
       if (self%enableRangeChecks) then
          if (                              timeIn < 0.0d0           )                      &
               & call Error_Report('time preceeds the Big Bang' //{introspection:location})
          if (self%collapsingUniverse .and. timeIn > self%timeMaximum)                      &
               & call Error_Report('time exceeds the Big Crunch'//{introspection:location})
       end if
       ! Set outputs.
       if (present(timeOut           )) timeOut           =                           timeIn
       if (present(expansionFactorOut)) expansionFactorOut= self%expansionFactor     (timeIn)
       if (present(collapsingOut     )) collapsingOut     = self%collapsingUniverse           &
            &                                              .and.                              &
            &                                               self%timeTurnaround     < timeIn
       return
    end if
    ! If we have an expansion factor, check that it is valid, and compute outputs as required.
    if (present(expansionFactorIn)) then
       ! Validate.
       if (self%enableRangeChecks) then
          if (                              expansionFactorIn <                       0.0d0) &
               & call Error_Report('expansion factor preceeds the Big Bang'    //{introspection:location})
          if (self%collapsingUniverse .and. expansionFactorIn > self%expansionFactorMaximum) &
               & call Error_Report('expansion factor exceeds maximum expansion'//{introspection:location})
       end if
       ! Determine collapse status.
       collapsingActual=.false.
       if (present(collapsingIn)) collapsingActual=collapsingIn
       ! Validate collapse status.
       if (collapsingActual .and. .not.self%collapsingUniverse) &
            & call Error_Report('epoch during collapsing phase specified, but universe does not collapse'//{introspection:location})
       ! Set outputs.
       if (present(timeOut           )) timeOut           =self%cosmicTime(expansionFactorIn,collapsingActual)
       if (present(expansionFactorOut)) expansionFactorOut=expansionFactorIn
       if (present(collapsingOut     )) collapsingOut     =collapsingActual
       return
    end if
    return
  end subroutine matterLambdaEpochValidate

  double precision function matterLambdaCosmicTime(self,expansionFactor,collapsingPhase)
    !!{
    Return the cosmological matter density in units of the critical density at the present day.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   )           :: expansionFactor
    logical                                         , intent(in   ), optional :: collapsingPhase
    logical                                                                   :: collapsingPhaseActual

    ! Validate the input.
    call self%epochValidate(                                         &
         &                  expansionFactorIn=expansionFactor      , &
         &                  collapsingIn     =collapsingPhase      , &
         &                  collapsingOut    =collapsingPhaseActual  &
         &                 )
    ! Ensure tabulation is initialized.
    if (.not.self%ageTableInitialized) call self%expansionFactorTabulate(self%ageTableTimeMinimum)
    ! Ensure that the tabulation spans a sufficient range of expansion factors.
    if (collapsingPhaseActual) then
       ! In collapsing phase just ensure that a sufficiently large expansion factor has been reached.
       do while (self%ageTableExpansionFactor(self%ageTableNumberPoints) < expansionFactor)
          self%ageTableTimeMaximum=min(self%ageTableTimeMaximum*ageTableIncrementFactor,self%timeTurnaround)
          call self%expansionFactorTabulate()
       end do
    else
       ! In expanding phase ensure that sufficiently small and large expansion factors have been reached.
       do while (self%ageTableExpansionFactor(                    1) > expansionFactor)
          self%ageTableTimeMinimum=    self%ageTableTimeMinimum/ageTableIncrementFactor
          call self%expansionFactorTabulate()
       end do
       do while (self%ageTableExpansionFactor(self%iTableTurnaround) < expansionFactor)
          self%ageTableTimeMaximum=max(self%ageTableTimeMaximum*ageTableIncrementFactor,self%timeTurnaround)
          call self%expansionFactorTabulate()
       end do
    end if
    ! Interpolate to get cosmic time.
    matterLambdaCosmicTime=self%interpolatorTime%interpolate(expansionFactor)
    ! Adjust for collapsing phase.
    if (collapsingPhaseActual) matterLambdaCosmicTime=self%timeMaximum-matterLambdaCosmicTime
    return
  end function matterLambdaCosmicTime

  double precision function matterLambdaTimeBigCrunch(self)
    !!{
    Return the time of the Big Crunch (or a negative value if no Big Crunch occurs).
    !!}
    implicit none
    class(cosmologyFunctionsMatterLambda), intent(inout) :: self

    if (self%collapsingUniverse) then
       matterLambdaTimeBigCrunch=self%timeMaximum
    else
       matterLambdaTimeBigCrunch=-1.0d0
    end if
    return
  end function matterLambdaTimeBigCrunch

  integer function matterLambdaCollapseODEs(a,t,dtda)
    !!{
    System of differential equations to solve for age vs. expansion factor.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision              , intent(in   ) :: a
    double precision, dimension(:), intent(in   ) :: t
    double precision, dimension(:), intent(  out) :: dtda
    !$GLC attributes unused :: t

    ! Compare the rate of change of time with expansion factor. For this ODE system we are always interested in the expanding
    ! phase of the Universe, so we use the absolute value of the expansion rate in case the universe is defined during a
    ! collapsing phase.
    dtda(1)=1.0d0/a/abs(self_%expansionRate(a))
    matterLambdaCollapseODEs=GSL_Success
    return
  end function matterLambdaCollapseODEs

  double precision function matterLambdaExpansionFactor(self,time)
    !!{
    Returns the expansion factor at cosmological time {\normalfont \ttfamily time}.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : varying_string
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    double precision                                                :: timeEffective, h
    logical                                                         :: remakeTable
    integer                                                         :: i
    type            (varying_string                ), save          :: message
    !$omp threadprivate(message)
    character       (len=13                        )                :: label

    ! Check if the time differs from the previous time.
    if (time /= self%timePrevious) then
       ! Quit on invalid input.
       if (time < 0.0d0) call Error_Report('cosmological time must be positive'//{introspection:location})
       ! Check if we need to recompute our table.
       if (self%ageTableInitialized) then
          remakeTable=(time < self%ageTableTime(1).or.time > self%ageTableTime(self%ageTableNumberPoints))
       else
          remakeTable=.true.
       end if
       if (remakeTable) call self%expansionFactorTabulate(time)
       ! Quit on invalid input.
       if (self%collapsingUniverse.and.time > self%timeMaximum) then
          write (label,'(e12.6)')      time
          message="cosmological time ["//trim(adjustl(label))//" Gyr] exceeds that at the Big Crunch ["
          write (label,'(e12.6)') self%timeMaximum
          message=message//trim(adjustl(label))//" Gyr]"
          call Error_Report(message//{introspection:location})
       end if
       ! Find the effective time to which to interpolate.
       if (self%collapsingUniverse) then
          if (time <= self%timeTurnaround) then
             timeEffective=                 time
          else
             timeEffective=self%timeMaximum-time
          end if
       else
          timeEffective   =                 time
       end if
       ! Perform the interpolation. We use a custom interpolator here. The expansion factor vs. time table is distributed almost
       ! uniformly in log(time) - it's not perfectly uniform in log time as we have to preserve numerical tolerance errors arising
       ! from expansion of the table. Therefore, we attempt to identify the index of the entry in the table by directly computing
       ! it from the logarithm of the effective time, but allow for the possibility that we might have to adjust that initial
       ! guess to find the correct index. After that, a standard linear interpolation is used.
       ! Initial guess at the index for interpolation.
       i=int((log(timeEffective)-self%ageTableTimeLogarithmicMinimum)*self%ageTableInverseDeltaLogTime)+1
       ! Check that we've found the correct index, adjust as necessary.
       do while (timeEffective < self%ageTableTime(i  ))
          i=i-1
       end do
       do while (timeEffective > self%ageTableTime(i+1))
          i=i+1
       end do
       ! Compute interpolating factor.
       h=     +(     timeEffective     -self%ageTableTime(i)) &
            & /(self%ageTableTime (i+1)-self%ageTableTime(i))
       ! Evaluate the interpolation.
       self%expansionFactorPrevious=+(1.0d0-h)*self%ageTableExpansionFactor(i  ) &
            &                       +       h *self%ageTableExpansionFactor(i+1)
       self%timePrevious           = time
    end if
    ! Return the stored expansion factor.
    matterLambdaExpansionFactor=self%expansionFactorPrevious
    return
  end function matterLambdaExpansionFactor

  double precision function matterLambdaExpansionRate(self,expansionFactor)
    !!{
    Returns the cosmological expansion rate, $\dot{a}/a$ at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsStandard, hubbleUnitsTime
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: expansionFactor

    if      (expansionFactor == self%expansionRateExpansionFactorPrevious(1)) then
       matterLambdaExpansionRate=self%expansionRatePrevious(1)
    else if (expansionFactor == self%expansionRateExpansionFactorPrevious(2)) then
       matterLambdaExpansionRate=self%expansionRatePrevious(2)
    else
       ! Required value is simply the Hubble parameter but expressed in units of inverse Gyr.
       self%expansionRateExpansionFactorPrevious(1)=self%expansionRateExpansionFactorPrevious(2)
       self%expansionRatePrevious               (1)=self%expansionRatePrevious               (2)
       self%expansionRateExpansionFactorPrevious(2)=+                                                                     expansionFactor
       self%expansionRatePrevious               (2)=+self                     %hubbleParameterEpochal(expansionFactor    =expansionFactor) &
            &                                       *self%cosmologyParameters_%HubbleConstant        (hubbleUnitsTime                    ) &
            &                                       /self%cosmologyParameters_%HubbleConstant        (hubbleUnitsStandard                )
       matterLambdaExpansionRate                  =+self%expansionRatePrevious(2)
    end if
    return
  end function matterLambdaExpansionRate

  double precision function matterLambdaHubbleParameterEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Returns the Hubble parameter at the request cosmological time, {\normalfont \ttfamily time}, or expansion factor, {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsStandard
    use :: Error               , only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual, sqrtArgument

    ! Validate the epoch.
    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    ! Compute the Hubble parameter at the specified expansion factor.
    sqrtArgument                      =max(                                              &
         &                                  self%cosmologyParameters_%OmegaMatter    ()  &
         &                                 /expansionFactorActual**3                     &
         &                                 +self%cosmologyParameters_%OmegaDarkEnergy()  &
         &                                 +self%cosmologyParameters_%OmegaCurvature ()  &
         &                                 /expansionFactorActual**2                   , &
         &                                 0.0d0                                         &
         &                                )
    matterLambdaHubbleParameterEpochal= self%cosmologyParameters_%HubbleConstant(hubbleUnitsStandard) &
         &                             *sqrt(sqrtArgument)
    ! Make the Hubble parameter negative if we are in the collapsing phase of the Universe.
    if (self%collapsingUniverse) then
       if    (present(time           )) then
          if    (time>self%timeTurnaround) then
             matterLambdaHubbleParameterEpochal=-abs(matterLambdaHubbleParameterEpochal)
          else
             matterLambdaHubbleParameterEpochal=+abs(matterLambdaHubbleParameterEpochal)
          end if
       else
          if (present(collapsingPhase)) then
             if (collapsingPhase         ) then
                matterLambdaHubbleParameterEpochal=-abs(matterLambdaHubbleParameterEpochal)
             else
                matterLambdaHubbleParameterEpochal=+abs(matterLambdaHubbleParameterEpochal)
             end if
          end if
       end if
    end if
    return
  end function matterLambdaHubbleParameterEpochal

  double precision function matterLambdaHubbleParameterRateOfChange(self,time,expansionFactor,collapsingPhase)
    !!{
    Returns the rate of change of the Hubble parameter at the request cosmological time, {\normalfont \ttfamily time}, or expansion factor, {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual

    ! Validate the epoch.
    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    matterLambdaHubbleParameterRateOfChange                                                                            &
         & = +0.5d0                                                                                                    &
         & *        self%hubbleParameterEpochal(expansionFactor=expansionFactorActual,collapsingPhase=collapsingPhase) &
         & *        self%expansionRate         (                expansionFactorActual                                ) &
         & *(                                                                                                          &
         &   -3.0d0*self%cosmologyParameters_%OmegaMatter    ()/expansionFactorActual**3                               &
         &   -2.0d0*self%cosmologyParameters_%OmegaCurvature ()/expansionFactorActual**2                               &
         &  )                                                                                                          &
         & /(                                                                                                          &
         &   +      self%cosmologyParameters_%OmegaMatter    ()/expansionFactorActual**3                               &
         &   +      self%cosmologyParameters_%OmegaDarkEnergy()                                                        &
         &   +      self%cosmologyParameters_%OmegaCurvature ()/expansionFactorActual**2                               &
         &  )
    return
  end function matterLambdaHubbleParameterRateOfChange

  double precision function matterLambdaOmegaMatterEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the matter density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsStandard
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual

    ! Validate the epoch.
    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    matterLambdaOmegaMatterEpochal                                                                         &
         & =   self%cosmologyParameters_%OmegaMatter           (                                         ) &
         &  *(                                                                                             &
         &     self%cosmologyParameters_%HubbleConstant        (hubbleUnitsStandard                      ) &
         &    /self                     %HubbleParameterEpochal(expansionFactor    =expansionFactorActual) &
         &   )**2                                                                                          &
         &  /expansionFactorActual**3
    return
  end function matterLambdaOmegaMatterEpochal

  double precision function matterLambdaMatterDensityEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the matter density at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual
    !$GLC attributes unused :: collapsingPhase

    ! Determine the actual expansion factor to use.
    if (present(time)) then
       if (present(expansionFactor)) then
          expansionFactorActual=-1.0d0
          call Error_Report('only one of time or expansion factor can be specified'//{introspection:location})
       else
          expansionFactorActual=self%expansionFactor(time)
       end if
    else
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else
          expansionFactorActual=-1.0d0
          call Error_Report('either a time or expansion factor must be specified'//{introspection:location})
       end if
    end if
    matterLambdaMatterDensityEpochal=self%cosmologyParameters_%omegaMatter()*self%cosmologyParameters_%densityCritical()/expansionFactorActual**3
    return
  end function matterLambdaMatterDensityEpochal

  double precision function matterLambdaOmegaMatterRateOfChange(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the rate of change of the matter density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual

    ! Validate the epoch.
    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    matterLambdaOmegaMatterRateOfChange                                                          &
         & =self%omegaMatterEpochal(time,expansionFactor,collapsingPhase)                        &
         & *(                                                                                    &
         &   -3.0d0*self%expansionRate              (     expansionFactorActual                ) &
         &   -2.0d0*self%hubbleParameterRateOfChange(time,expansionFactor      ,collapsingPhase) &
         &   /      self%hubbleParameterEpochal     (time,expansionFactor      ,collapsingPhase) &
         &  )
    return
  end function matterLambdaOmegaMatterRateOfChange

  double precision function matterLambdaOmegaDarkEnergyEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the dark energy density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsStandard
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual

    ! Validate the epoch.
    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    matterLambdaOmegaDarkEnergyEpochal                                                                     &
         & =   self%cosmologyParameters_%OmegaDarkEnergy       (                                         ) &
         &  *(                                                                                             &
         &     self%cosmologyParameters_%HubbleConstant        (hubbleUnitsStandard                      ) &
         &    /self%                     HubbleParameterEpochal(expansionFactor    =expansionFactorActual) &
         &   )**2
    return
  end function matterLambdaOmegaDarkEnergyEpochal

  double precision function matterLambdaTemperatureCMBEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the temperature of the CMB at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual

    ! Validate the epoch.
    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    matterLambdaTemperatureCMBEpochal=self%cosmologyParameters_%temperatureCMB()/expansionFactorActual
    return
  end function matterLambdaTemperatureCMBEpochal

  subroutine matterLambdaDensityScalingEarlyTime(self,dominateFactor,densityPower,expansionFactorDominant,OmegaDominant)
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   )           :: dominateFactor
    double precision                                , intent(  out)           :: densityPower  , expansionFactorDominant
    double precision                                , intent(  out), optional :: OmegaDominant

    ! For matter and cosmological constant, matter always dominates at early times.
    densityPower=-3.0d0 ! Power-law scaling of matter density with expansion factor.
    ! Choose present day as default - will be used if no other densities present (i.e. Einstein-de Sitter).
    expansionFactorDominant=self%dominationEpochMatter(dominateFactor)
    ! Return the density parameter in the dominant species if required.
    if (present(OmegaDominant)) OmegaDominant=self%cosmologyParameters_%OmegaMatter()
    return
  end subroutine matterLambdaDensityScalingEarlyTime

  double precision function matterLambdaDominationEpochMatter(self,dominateFactor)
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: dominateFactor
    double precision                                                :: aMatterEquality                  , expansionFactorDominantCurvature, &
         &                                                             expansionFactorDominantDarkEnergy

    ! Choose present day as default - will be used if no other densities present (i.e. Einstein-de Sitter).
    matterLambdaDominationEpochMatter=1.0d0
    if (self%cosmologyParameters_%OmegaDarkEnergy() /= 0.0d0) then
       ! Find the expansion factor of matter-dark energy equality.
       aMatterEquality=self%equalityEpochMatterDarkEnergy(requestTypeExpansionFactor)
       ! Find the earlier expansion factor at which matter dominates by the specified amount (ratio of matter
       ! to dark energy density scales as the cube of expansion factor).
       expansionFactorDominantDarkEnergy=aMatterEquality/dominateFactor**(1.0d0/3.0d0)
       ! Choose earliest expansion factor.
       matterLambdaDominationEpochMatter=min(matterLambdaDominationEpochMatter,expansionFactorDominantDarkEnergy)
    end if
    if (self%cosmologyParameters_%OmegaCurvature() /= 0.0d0) then
       ! Find the expansion factor of matter-curvature equality.
       aMatterEquality=self%equalityEpochMatterCurvature(requestTypeExpansionFactor)
       ! Find the earlier expansion factor at which matter dominates by the specified amount (ratio of matter
       ! to curvature density scales as the expansion factor).
       expansionFactorDominantCurvature=aMatterEquality/dominateFactor
       ! Choose earliest expansion factor.
       matterLambdaDominationEpochMatter=min(matterLambdaDominationEpochMatter,expansionFactorDominantCurvature)
    end if
    return
  end function matterLambdaDominationEpochMatter

  double precision function matterLambdaEqualityEpochMatterDarkEnergy(self,requestType)
    !!{
    Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    !!}
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor, requestTypeTime
    implicit none
    class  (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    integer                                , intent(in   ), optional :: requestType
    integer                                                          :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if
    matterLambdaEqualityEpochMatterDarkEnergy       &
         & =(                                       &
         &        self%cosmologyParameters_%OmegaMatter    ()  &
         &   /abs(self%cosmologyParameters_%OmegaDarkEnergy()) &
         &  )**(1.0d0/3.0d0)
    if (requestTypeActual == requestTypeTime)                          &
         &                   matterLambdaEqualityEpochMatterDarkEnergy &
         &  =self%cosmicTime(matterLambdaEqualityEpochMatterDarkEnergy)
    return
  end function matterLambdaEqualityEpochMatterDarkEnergy

  double precision function matterLambdaEqualityEpochMatterCurvature(self,requestType)
    !!{
    Return the epoch of matter-curvature magnitude equality (either expansion factor or cosmic time).
    !!}
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor, requestTypeTime
    implicit none
    class  (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    integer                                , intent(in   ), optional :: requestType
    integer                                                          :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if
    matterLambdaEqualityEpochMatterCurvature=self%cosmologyParameters_%OmegaMatter()/abs(self%cosmologyParameters_%OmegaCurvature())
    if (requestTypeActual == requestTypeTime)                         &
         &                  matterLambdaEqualityEpochMatterCurvature  &
         & =self%cosmicTime(matterLambdaEqualityEpochMatterCurvature)
    return
  end function matterLambdaEqualityEpochMatterCurvature

  double precision function matterLambdaEqualityEpochMatterRadiation(self,requestType)
    !!{
    Return the epoch of matter-radiation magnitude equality (either expansion factor or cosmic time).
    !!}
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor, requestTypeTime
    implicit none
    class  (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    integer                                , intent(in   ), optional :: requestType
    integer                                                          :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if
    matterLambdaEqualityEpochMatterRadiation=self%cosmologyParameters_%OmegaRadiation()/self%cosmologyParameters_%OmegaMatter()
    if (requestTypeActual == requestTypeTime)                         &
         &                  matterLambdaEqualityEpochMatterRadiation  &
         & =self%cosmicTime(matterLambdaEqualityEpochMatterRadiation)
    return
  end function matterLambdaEqualityEpochMatterRadiation

  subroutine matterLambdaMakeExpansionFactorTable(self,time)
    !!{
    Builds a table of expansion factor vs. time.
    !!}
    use :: Cosmology_Parameters , only : hubbleUnitsTime
    use :: Numerical_Ranges     , only : Make_Range     , rangeTypeLogarithmic
    use :: Numerical_ODE_Solvers, only : odeSolver
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout), target       :: self
    double precision                                , intent(in   ), optional     :: time
    double precision                                , parameter                   :: odeToleranceAbsolute            =1.0d-9, odeToleranceRelative =1.0d-9
    double precision                                , allocatable  , dimension(:) :: ageTableExpansionFactorTemporary       , ageTableTimeTemporary
    double precision                                               , dimension(1) :: expansionFactor
    integer                                                                       :: iTime                                  , prefixPointCount
    double precision                                                              :: OmegaDominant                          , densityPower                , &
         &                                                                           tDominant                              , timeActual                  , &
         &                                                                           expansionFactorDominant
    type            (odeSolver                     )                              :: solver

    ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
    call self%densityScalingEarlyTime(matterLambdaDominateFactor,densityPower,expansionFactorDominant,OmegaDominant)
    ! Find the corresponding time. Note that we use the absolute value of the Hubble parameter here - in cases where the universe
    ! is collapsing at the present epoch we need to know the expansion rate (i.e. Hubble parameter) at the equivalent expansion
    ! factor during the expansion phase.
    tDominant=-2.0d0/densityPower/abs(self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime))/sqrt(OmegaDominant)/expansionFactorDominant**(0.5d0*densityPower)
    ! Find minimum and maximum times to tabulate.
    if (present(time)) then
       timeActual=time
       do while (self%ageTableTimeMinimum > min(timeActual,tDominant)/2.0d0)
          self%ageTableTimeMinimum=self%ageTableTimeMinimum/ageTableIncrementFactor
       end do
       do while (self%ageTableTimeMaximum < max(timeActual,tDominant)*2.0d0)
          self%ageTableTimeMaximum=self%ageTableTimeMaximum*ageTableIncrementFactor
       end do
    else
       do while (self%ageTableTimeMinimum > tDominant/2.0d0)
          self%ageTableTimeMinimum=self%ageTableTimeMinimum/ageTableIncrementFactor
       end do
       do while (self%ageTableTimeMaximum < tDominant*2.0d0)
          self%ageTableTimeMaximum=self%ageTableTimeMaximum*ageTableIncrementFactor
       end do
    end if
    if (self%collapsingUniverse) self%ageTableTimeMaximum=min(self%ageTableTimeMaximum,self%timeTurnaround)

    ! Determine number of points to tabulate.
    self%ageTableNumberPoints=int(log10(self%ageTableTimeMaximum/self%ageTableTimeMinimum)*dble(ageTableNPointsPerDecade))+1
    self%ageTableTimeMaximum =self%ageTableTimeMinimum*10.0d0**(dble(self%ageTableNumberPoints)/dble(ageTableNPointsPerDecade))
    if (self%collapsingUniverse) self%ageTableTimeMaximum=min(self%ageTableTimeMaximum,self%timeTurnaround)

    ! Deallocate arrays if currently allocated.
    if (allocated(self%ageTableTime)) then
       ! Determine number of points that are being added at the start of the array.
       prefixPointCount=int(log10(self%ageTableTime(1)/self%ageTableTimeMinimum)*dble(ageTableNPointsPerDecade)+0.5d0)
       call Move_Alloc(self%ageTableTime           ,ageTableTimeTemporary           )
       call Move_Alloc(self%ageTableExpansionFactor,ageTableExpansionFactorTemporary)
       ! Allocate the arrays to current required size.
       allocate(self%ageTableTime           (self%ageTableNumberPoints))
       allocate(self%ageTableExpansionFactor(self%ageTableNumberPoints))
       ! Create set of grid points in time variable.
       self%ageTableTime=Make_Range(self%ageTableTimeMinimum,self%ageTableTimeMaximum,self%ageTableNumberPoints,rangeTypeLogarithmic)
       ! Set the expansion factors to a negative value to indicate they are not yet computed.
       self%ageTableExpansionFactor=-1.0d0
       ! Paste in the previously computed regions.
       self%ageTableTime           (prefixPointCount+1:prefixPointCount+size(ageTableTimeTemporary))=ageTableTimeTemporary
       self%ageTableExpansionFactor(prefixPointCount+1:prefixPointCount+size(ageTableTimeTemporary))=ageTableExpansionFactorTemporary
       ! Deallocate the temporary arrays.
       deallocate(ageTableTimeTemporary           )
       deallocate(ageTableExpansionFactorTemporary)
    else
       ! Allocate the arrays to current required size.
       allocate(self%ageTableTime           (self%ageTableNumberPoints))
       allocate(self%ageTableExpansionFactor(self%ageTableNumberPoints))
       ! Create set of grid points in time variable.
       self%ageTableTime=Make_Range(self%ageTableTimeMinimum,self%ageTableTimeMaximum,self%ageTableNumberPoints,rangeTypeLogarithmic)
       ! Set the expansion factors to a negative value to indicate they are not yet computed.
       self%ageTableExpansionFactor=-1.0d0
    end if
    ! Compute quantities required for table interpolation.
    self%ageTableTimeLogarithmicMinimum=log(self%ageTableTimeMinimum)
    self%ageTableInverseDeltaLogTime   =dble(self%ageTableNumberPoints-1)/log(self%ageTableTimeMaximum/self%ageTableTimeMinimum)
    ! For the initial time, we approximate that we are at sufficiently early times that a single component dominates the Universe
    ! and use the appropriate analytic solution. Note that we use the absolute value of the Hubble parameter here - in cases where
    ! the universe is collapsing at the present epoch we need to know the expansion rate (i.e. Hubble parameter) at the equivalent
    ! expansion factor during the expansion phase.
    if (self%ageTableExpansionFactor(1) < 0.0d0)                             &
         &    self%ageTableExpansionFactor                (               1) &
         & =(                                                                &
         &   -0.5d0                                                          &
         &   *densityPower                                                   &
         &   *self%ageTableTime                            (              1) &
         &   *abs(self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime)) &
         &   *sqrt(OmegaDominant)                                            &
         &  )**(-2.0d0/densityPower)
    ! Solve ODE to get corresponding expansion factors.
    self%iTableTurnaround =  self%ageTableNumberPoints
    self_                 => self
    solver                =  odeSolver(1_c_size_t,matterLambdaAgeTableODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
    do iTime=2,self%ageTableNumberPoints
       ! Find the position in the table corresponding to turn around if we have a collapsing Universe.
       if     (                                                   &
            &   self%collapsingUniverse                           &
            &  .and.                                              &
            &   self%ageTableTime(iTime-1) <  self%timeTurnaround &
            &  .and.                                              &
            &   self%ageTableTime(iTime  ) >= self%timeTurnaround &
            & ) self%iTableTurnaround=iTime
       ! Compute the expansion factor if it is not already computed.
       if (self%ageTableExpansionFactor(iTime) < 0.0d0) then
          timeActual        =self%ageTableTime           (iTime-1)
          expansionFactor(1)=self%ageTableExpansionFactor(iTime-1)
          call solver%solve(timeActual,self%ageTableTime(iTime),expansionFactor)
          self%ageTableExpansionFactor(iTime)=expansionFactor(1)
       end if
    end do
    ! Build the interpolator.
    if (allocated(self%interpolatorTime)) deallocate(self%interpolatorTime)
    allocate(self%interpolatorTime)
    self%interpolatorTime=interpolator(self%ageTableExpansionFactor,self%ageTableTime)
    ! Flag that the table is now initialized.
    self%ageTableInitialized=.true.
    return
  end subroutine matterLambdaMakeExpansionFactorTable

  integer function matterLambdaAgeTableODEs(t,a,dadt)
    !!{
    System of differential equations to solve for expansion factor vs. age.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision              , intent(in   ) :: t
    double precision, dimension(:), intent(in   ) :: a
    double precision, dimension(:), intent(  out) :: dadt
    !$GLC attributes unused :: t

    ! For this ODE system we are always interested in the expanding phase of the Universe, so we use the absolute value of the
    ! expansion rate in case the universe is defined during a collapsing phase.
    dadt(1)=a(1)*abs(self_%expansionRate(a(1)))
    matterLambdaAgeTableODEs=GSL_Success
  end function matterLambdaAgeTableODEs

  double precision function matterLambdaTimeAtDistanceComoving(self,comovingDistance)
    !!{
    Returns the cosmological time corresponding to given {\normalfont \ttfamily comovingDistance}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: comovingDistance
    double precision                                                :: time
    logical                                                         :: remakeTable

    ! Quit on invalid input.
    if (comovingDistance < 0.0d0) call Error_Report('comoving distance must be positive'//{introspection:location})
    ! Check if we need to recompute our table.
    remakeTable=.true.
    do while (remakeTable)
       if (self%distanceTableInitialized) then
          remakeTable=self%distanceTableComovingDistance(1) < comovingDistance
          time       =0.5d0*self%distanceTableTime(1)
       else
          remakeTable=.true.
          time       =self%distanceTableTimeMinimum
       end if
       ! Remake table if necessary.
       if (remakeTable) call self%distanceTabulate(time)
    end do
    ! Interpolate to get the comoving distance.
    matterLambdaTimeAtDistanceComoving=self%interpolatorDistanceInverse%interpolate(-comovingDistance)
    return
  end function matterLambdaTimeAtDistanceComoving

  double precision function matterLambdaDistanceComoving(self,time)
    !!{
    Returns the comoving distance to cosmological time {\normalfont \ttfamily time}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    double precision                                , parameter     :: toleranceRelative=1.0d-6
    logical                                                         :: remakeTable
    character       (len=12                        )                :: timeNowLabel            , timeLabel

    ! Quit on invalid input.
    call self%epochValidate(timeIn=time)
    if (time > self%cosmicTime(1.0d0)*(1.0d0+toleranceRelative)) then
       write (timeLabel   ,'(e12.6)') time
       write (timeNowLabel,'(e12.6)') self%cosmicTime(1.0d0)
       call Error_Report('cosmological time ['//trim(timeLabel)//'] must be in the past [≤'//trim(timeNowLabel)//']'//{introspection:location})
    end if
    ! Check if we need to recompute our table.
    if (self%distanceTableInitialized) then
       remakeTable=(time < self%distanceTableTime(1).or.time > self%distanceTableTime(self%distanceTableNumberPoints))
    else
       remakeTable=.true.
    end if
    if (remakeTable) call self%distanceTabulate(time)
    ! Quit on invalid input.
    if (self%collapsingUniverse.and.time>self%timeMaximum) &
         & call Error_Report('cosmological time exceeds that at the Big Crunch'//{introspection:location})
    ! Interpolate to get the comoving distance.
    matterLambdaDistanceComoving=self%interpolatorDistance%interpolate(time)
    return
  end function matterLambdaDistanceComoving

  double precision function matterLambdaDistanceLuminosity(self,time)
    !!{
    Returns the luminosity distance to cosmological time {\normalfont \ttfamily time}.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    ! Compute the luminosity distance.
    matterLambdaDistanceLuminosity=self%distanceComoving(time)/self%expansionFactor(time)
    return
  end function matterLambdaDistanceLuminosity

  double precision function matterLambdaDistanceAngular(self,time,timeOrigin) result(distance)
    !!{
    Returns the angular diameter distance to cosmological time {\normalfont \ttfamily time}.
    !!}
    use :: Cosmology_Parameters            , only : hubbleUnitsTime
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : megaParsec     , gigaYear
    use :: Numerical_Constants_Physical    , only : speedLight
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   )           :: time
    double precision                                , intent(in   ), optional :: timeOrigin
    double precision                                                          :: distanceComoving, distanceComovingOrigin, &
         &                                                                       distanceHubble  , OmegaCurvature

    if (present(timeOrigin)) then
       ! Case with the origin of the distance at some z>0.
       if (timeOrigin < time) call Error_Report('expected timeOrigin ≥ time'//{introspection:location})
       ! Use the solution from (Peebles, 1993, "Principles of Physical Cosmology", pp 336–33; see also Hogg 1999; equation 19;
       ! arXiv:astro-ph/9905116; https://ui.adsabs.harvard.edu/abs/1999astro.ph..5116H).
       distanceComoving      =+self%distanceComoving(time      )
       distanceComovingOrigin=+self%distanceComoving(timeOrigin)
       distanceHubble        =+speedLight                                                &
            &                 /self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime) &
            &                 *gigaYear                                                  &
            &                 /megaParsec
       OmegaCurvature        =+self%cosmologyParameters_%OmegaCurvature()
       if (OmegaCurvature >= 0.0d0) then
          distance=+(                                &
               &     +distanceComoving               &
               &     *sqrt(                          &
               &           +1.0d0                    &
               &           +OmegaCurvature           &
               &           *(                        &
               &             +distanceComovingOrigin &
               &             /distanceHubble         &
               &            )**2                     &
               &          )                          &
               &     -distanceComovingOrigin         &
               &     *sqrt(                          &
               &           +1.0d0                    &
               &           +OmegaCurvature           &
               &           *(                        &
               &             +distanceComoving       &
               &             /distanceHubble         &
               &            )**2                     &
               &          )                          &
               &    )                                &
               &   *self%expansionFactor(time)
       else
          ! For negative curvature cosmologies the above equation is not valid - see, for example, Hogg (1999;
          ! arXiv:astro-ph/9905116; https://ui.adsabs.harvard.edu/abs/1999astro.ph..5116H).
          distance=+0.0d0
          call Error_Report('angular diameter distance between two times not implemented for Ωₖ < 0'//{introspection:location})
       end if
    else
       ! Standard case with the origin of the distance at z=0.
       distance=+self%distanceComoving(time) &
            &   *self%expansionFactor (time)
    end if
    return
  end function matterLambdaDistanceAngular

  double precision function matterLambdaDistanceComovingConvert(self,output,distanceLuminosity,distanceModulus,distanceModulusKCorrected,redshift)
    !!{
    Convert between different measures of distance.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    integer                                         , intent(in   )           :: output
    double precision                                , intent(in   ), optional :: distanceModulus    , distanceModulusKCorrected, &
         &                                                                       redshift           , distanceLuminosity
    logical                                                                   :: gotComovingDistance
    double precision                                                          :: comovingDistance   , luminosityDistance

    ! Check if we need to recompute our table.
    if (.not.self%distanceTableInitialized) call self%distanceTabulate(self%cosmicTime(1.0d0))
    ! Convert to comoving distance from whatever was supplied.
    gotComovingDistance=.false.
    comovingDistance   =-1.0d0
    if (present(distanceLuminosity).or.present(distanceModulus)) then
       if (present(distanceLuminosity)) then
          luminosityDistance=distanceLuminosity
       else
          luminosityDistance=10.0d0**((distanceModulus-25.0d0)/5.0d0)
       end if
       do while (luminosityDistance > -self%distanceTableLuminosityDistanceNegated(1))
          call self%distanceTabulate(0.5d0*self%distanceTableTimeMinimum)
       end do
       comovingDistance   =-self%interpolatorLuminosityDistance%interpolate(-luminosityDistance)
       gotComovingDistance=.true.
    else if (present(distanceModulusKCorrected)) then
       luminosityDistance=10.0d0**((distanceModulusKCorrected-25.0d0)/5.0d0)
       do while (luminosityDistance > -self%distanceTableLuminosityDistanceKCorrectedNegated(1))
          call self%distanceTabulate(0.5d0*self%distanceTableTimeMinimum)
       end do
       comovingDistance   =-self%interpolatorLuminosityDistanceKCorrected%interpolate(-luminosityDistance)
       gotComovingDistance=.true.
    end if
    if (present(redshift)) then
       comovingDistance   =self%distanceComoving(self%cosmicTime(self%expansionFactorFromRedshift(redshift)))
       gotComovingDistance=.true.
    end if
    if (.not.gotComovingDistance) &
         & call Error_Report('no distance measure provided'//{introspection:location})
    ! Convert to required distance measure.
    select case (output)
    case (distanceTypeComoving)
       matterLambdaDistanceComovingConvert=comovingDistance
    case default
       matterLambdaDistanceComovingConvert=-1.0d0
       call Error_Report('unrecognized output option'//{introspection:location})
    end select
    return
  end function matterLambdaDistanceComovingConvert

  double precision function matterLambdaDistanceParticleHorizonComoving(self,time)
    !!{
    Returns the comoving distance to the particle horizon at cosmological time {\normalfont \ttfamily time}.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    type            (integrator                    )                :: integrator_

    integrator_                                =integrator           (integrandParticleHorizon,toleranceRelative=1.0d-6)
    matterLambdaDistanceParticleHorizonComoving=integrator_%integrate(0.0d0                   ,                  time  )
    return

  contains

    double precision function integrandParticleHorizon(time)
      !!{
      Integrand used to compute the distance to the comoving particle horizon.
      !!}
      use :: Numerical_Constants_Physical    , only : speedLight
      use :: Numerical_Constants_Astronomical, only : megaParsec, gigaYear
      implicit none
      double precision, intent(in   ) :: time

      if (time > 0.0d0) then
         integrandParticleHorizon=+     speedLight            &
              &                   *     gigaYear              &
              &                   /     megaParsec            &
              &                   /self%expansionFactor(time)
      else
         integrandParticleHorizon=0.0d0
      end if
      return
    end function integrandParticleHorizon
    
  end function matterLambdaDistanceParticleHorizonComoving
  
  subroutine matterLambdaMakeDistanceTable(self,time)
    !!{
    Builds a table of distance vs. time.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Make_Range   , rangeTypeLogarithmic
    use :: Numerical_Constants_Astronomical, only : gigaYear  , megaParsec
    use :: Numerical_Constants_Physical    , only : speedLight
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout), target :: self
    double precision                                , intent(in   )         :: time
    double precision                                , parameter             :: toleranceAbsolute=1.0d-5, toleranceRelative=1.0d-5
    integer                                                                 :: iTime
    type            (integrator                    )                        :: integrator_

    ! Find minimum and maximum times to tabulate.
    self%distanceTableTimeMinimum=min(self%distanceTableTimeMinimum,0.5d0*time)
    self%distanceTableTimeMaximum=    self%cosmicTime(1.0d0)
    ! Determine number of points to tabulate.
    self%distanceTableNumberPoints=int(log10(self%distanceTableTimeMaximum/self%distanceTableTimeMinimum)*dble(distanceTableNPointsPerDecade))+1
    ! Deallocate arrays if currently allocated.
    if (allocated(self%distanceTableTime                               )) deallocate(self%distanceTableTime                               )
    if (allocated(self%distanceTableComovingDistance                   )) deallocate(self%distanceTableComovingDistance                   )
    if (allocated(self%distanceTableComovingDistanceNegated            )) deallocate(self%distanceTableComovingDistanceNegated            )
    if (allocated(self%distanceTableLuminosityDistanceNegated          )) deallocate(self%distanceTableLuminosityDistanceNegated          )
    if (allocated(self%distanceTableLuminosityDistanceKCorrectedNegated)) deallocate(self%distanceTableLuminosityDistanceKCorrectedNegated)
    ! Allocate the arrays to current required size.
    allocate(self%distanceTableTime                               (self%distanceTableNumberPoints))
    allocate(self%distanceTableComovingDistance                   (self%distanceTableNumberPoints))
    allocate(self%distanceTableComovingDistanceNegated            (self%distanceTableNumberPoints))
    allocate(self%distanceTableLuminosityDistanceNegated          (self%distanceTableNumberPoints))
    allocate(self%distanceTableLuminosityDistanceKCorrectedNegated(self%distanceTableNumberPoints))
    ! Create the range of times.
    self% distanceTableTime=Make_Range(self%distanceTableTimeMinimum,self%distanceTableTimeMaximum,self%distanceTableNumberPoints,rangeTypeLogarithmic)
    ! Integrate to get the comoving distance.
    integrator_ =  integrator(                                                         &
         &                                      matterLambdaComovingDistanceIntegrand, &
         &                    toleranceAbsolute=toleranceAbsolute                    , &
         &                    toleranceRelative=toleranceRelative                      &
         &                   )
    self_       => self
    do iTime=1,self%distanceTableNumberPoints
       self%distanceTableComovingDistance(iTime)=+integrator_%integrate(                                                                              &
            &                                                           self%expansionFactor(self%distanceTableTime(iTime                         )), &
            &                                                           self%expansionFactor(self%distanceTableTime(self%distanceTableNumberPoints))  &
            &                                                          )                                                                              &
            &                                    *speedLight                                                                                          &
            &                                    *gigaYear                                                                                            &
            &                                    /megaParsec       
       self              %distanceTableLuminosityDistanceNegated              (iTime)   &
            & =      self%distanceTableComovingDistance                       (iTime)   &
            &       /self%expansionFactor              (self%distanceTableTime(iTime))
       self              %distanceTableLuminosityDistanceKCorrectedNegated    (iTime)   &
            & =      self%distanceTableComovingDistance                       (iTime)   &
            &  /sqrt(self%expansionFactor              (self%distanceTableTime(iTime)))
    end do
    ! Make a negated copy of the distances so that we have an increasing array for use in interpolation routines.
    self%distanceTableComovingDistanceNegated            =-self%distanceTableComovingDistance
    self%distanceTableLuminosityDistanceNegated          =-self%distanceTableLuminosityDistanceNegated
    self%distanceTableLuminosityDistanceKCorrectedNegated=-self%distanceTableLuminosityDistanceKCorrectedNegated
    ! Build interpolators.
    if (allocated(self%interpolatorDistance                    )) deallocate(self%interpolatorDistance                    )
    if (allocated(self%interpolatorDistanceInverse             )) deallocate(self%interpolatorDistanceInverse             )
    if (allocated(self%interpolatorLuminosityDistance          )) deallocate(self%interpolatorLuminosityDistance          )
    if (allocated(self%interpolatorLuminosityDistanceKCorrected)) deallocate(self%interpolatorLuminosityDistanceKCorrected)
    allocate(self%interpolatorDistance                    )
    allocate(self%interpolatorDistanceInverse             )
    allocate(self%interpolatorLuminosityDistance          )
    allocate(self%interpolatorLuminosityDistanceKCorrected)
    self%interpolatorDistance                    =interpolator(self%distanceTableTime                               ,self%distanceTableComovingDistance       )
    self%interpolatorDistanceInverse             =interpolator(self%distanceTableComovingDistanceNegated            ,self%distanceTableTime                   )
    self%interpolatorLuminosityDistance          =interpolator(self%distanceTableLuminosityDistanceNegated          ,self%distanceTableComovingDistanceNegated)
    self%interpolatorLuminosityDistanceKCorrected=interpolator(self%distanceTableLuminosityDistanceKCorrectedNegated,self%distanceTableComovingDistanceNegated)
    ! Flag that the table is now initialized.
    self%distanceTableInitialized=.true.
    return
  end subroutine matterLambdaMakeDistanceTable

  double precision function matterLambdaComovingDistanceIntegrand(expansionFactor)
    !!{
    Integrand function used in computing the comoving distance.
    !!}
    implicit none
    double precision, intent(in   ) :: expansionFactor

    matterLambdaComovingDistanceIntegrand=1.0d0/expansionFactor**2/self_%expansionRate(expansionFactor)
    return
  end function matterLambdaComovingDistanceIntegrand

  double precision function matterLambdaEquationOfStateDarkEnergy(self,time,expansionFactor)
    !!{
    Return the dark energy equation of state.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor, time
    !$GLC attributes unused :: self, time, expansionFactor

    matterLambdaEquationOfStateDarkEnergy=-1.0d0
    return
  end function matterLambdaEquationOfStateDarkEnergy

  double precision function matterLambdaExponentDarkEnergy(self,time,expansionFactor)
    !!{
    Return the dark energy equation of state.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor, time
    !$GLC attributes unused :: self, time, expansionFactor

    matterLambdaExponentDarkEnergy=0.0d0
    return
  end function matterLambdaExponentDarkEnergy
