  !! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of the cosmological functions class for cosmologies consisting of collisionless
  !% matter plus a cosmological constant.

  !# <cosmologyFunctions name="cosmologyFunctionsMatterLambda" />
  use FGSL
  use Cosmology_Parameters

  integer         , parameter :: matterLambdaAgeTableNPointsPerDecade     =300
  double precision, parameter :: matterLambdaAgeTableNPointsPerOctave     =dble(matterLambdaAgeTableNPointsPerDecade)*log(2.0d0)/log(10.0d0)
  double precision, parameter :: matterLambdaAgeTableIncrementFactor      =exp(int(matterLambdaAgeTableNPointsPerOctave+1.0d0)*log(10.0d0)/dble(matterLambdaAgeTableNPointsPerDecade))
  integer         , parameter :: matterLambdaDistanceTableNPointsPerDecade=100

  ! Factor by which one component of Universe must dominate others such that we can ignore the others.
  double precision, parameter :: matterLambdaDominateFactor               =100.0d0

  type, extends(cosmologyFunctionsClass) :: cosmologyFunctionsMatterLambda
     !% A cosmological functions class for cosmologies consisting of matter plus a cosmological constant.
     private
     class           (cosmologyParametersClass), pointer                   :: cosmology
     logical                                                               :: collapsingUniverse                        =.false.
     integer                                                               :: iTableTurnaround
     double precision                                                      :: expansionFactorMaximum                            , expansionFactorPrevious                =-1.0d0, &
          &                                                                   timeMaximum                                       , timePrevious                           =-1.0d0, &
          &                                                                   timeTurnaround
     logical                                                               :: ageTableInitialized                       =.false.
     integer                                                               :: ageTableNumberPoints
     double precision                                                      :: ageTableTimeMaximum                       =20.0d0 , ageTableTimeMinimum                    =1.0d-4
     double precision                          , allocatable, dimension(:) :: ageTableExpansionFactor                           , ageTableTime
     type            (fgsl_interp             )                            :: interpolationObject                               , interpolationObjectInverse
     type            (fgsl_interp_accel       )                            :: interpolationAccelerator                          , interpolationAcceleratorInverse
     logical                                                               :: resetInterpolation                        =.true.
     logical                                                               :: resetInterpolationInverse                 =.true.
     logical                                                               :: distanceTableInitialized                  =.false.
     integer                                                               :: distanceTableNumberPoints
     double precision                                                      :: distanceTableTimeMaximum                          , distanceTableTimeMinimum               =1.0d-4
     double precision                          , allocatable, dimension(:) :: distanceTableComovingDistance                     , distanceTableComovingDistanceNegated          , &
          &                                                                   distanceTableLuminosityDistanceNegated            , distanceTableTime
     type            (fgsl_interp             )                            :: interpolationObjectDistance                       , interpolationObjectDistanceInverse            , &
          &                                                                   interpolationObjectLuminosityDistance
     type            (fgsl_interp_accel       )                            :: interpolationAcceleratorDistance                  , interpolationAcceleratorDistanceInverse       , &
          &                                                                   interpolationAcceleratorLuminosityDistance
     logical                                                               :: resetInterpolationDistance                =.true. , resetInterpolationDistanceInverse      =.true., &
          &                                                                   resetInterpolationLuminosityDistance      =.true.
   contains
     !@ <objectMethods>
     !@   <object>cosmologyFunctionsMatterLambda</object>
     !@   <objectMethod>
     !@     <method>distancetabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate comoving distance as a function of cosmic time.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>expansionFactorTabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate expansion factor as a function of cosmic time.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     !# <workaround type="gfortran" PR="58356" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58356">
     !# final     :: matterLambdaDestructor
     !# </workaround>
     procedure :: stateStore                   =>matterLambdaStateStore
     procedure :: stateRestore                 =>matterLambdaStateRestore
     procedure :: expansionFactorIsValid       =>matterLambdaExpansionFactorIsValid
     procedure :: cosmicTimeIsValid            =>matterLambdaCosmicTimeIsValid
     procedure :: cosmicTime                   =>matterLambdaCosmicTime
     procedure :: expansionFactor              =>matterLambdaExpansionFactor
     procedure :: expansionRate                =>matterLambdaExpansionRate
     procedure :: hubbleParameterEpochal       =>matterLambdaHubbleParameterEpochal
     procedure :: densityScalingEarlyTime      =>matterLambdaDensityScalingEarlyTime
     procedure :: omegaMatterEpochal           =>matterLambdaOmegaMatterEpochal
     procedure :: omegaDarkEnergyEpochal       =>matterLambdaOmegaDarkEnergyEpochal
     procedure :: equationOfStateDarkEnergy    =>matterLambdaEquationOfStateDarkEnergy
     procedure :: exponentDarkEnergy           =>matterLambdaExponentDarkEnergy
     procedure :: equalityEpochMatterDarkEnergy=>matterLambdaEqualityEpochMatterDarkEnergy
     procedure :: equalityEpochMatterCurvature =>matterLambdaEqualityEpochMatterCurvature
     procedure :: dominationEpochMatter        =>matterLambdaDominationEpochMatter
     procedure :: temperatureCMBEpochal        =>matterLambdaTemperatureCMBEpochal
     procedure :: distanceComoving             =>matterLambdaDistanceComoving
     procedure :: timeAtDistanceComoving       =>matterLambdaTimeAtDistanceComoving
     procedure :: distanceComovingConvert      =>matterLambdaDistanceComovingConvert
     procedure :: expansionFactorTabulate      =>matterLambdaMakeExpansionFactorTable
     procedure :: distanceTabulate             =>matterLambdaMakeDistanceTable
  end type cosmologyFunctionsMatterLambda

  ! Module scope pointer to the current object.
  class(cosmologyFunctionsMatterLambda), pointer :: matterLambdaSelfGlobal
  !$omp threadprivate(matterLambdaSelfGlobal)

  interface cosmologyFunctionsMatterLambda
     !% Constructors for the matter plus cosmological constant cosmological functions class.
     module procedure matterLambdaDefaultConstructor
     module procedure matterLambdaConstructor
  end interface cosmologyFunctionsMatterLambda

contains

  function matterLambdaDefaultConstructor()
    !% Default constructor for the matter plus cosmological constant cosmological functions class.
    use Cosmology_Parameters
    implicit none
    type (cosmologyFunctionsMatterLambda), target  :: matterLambdaDefaultConstructor
    class(cosmologyParametersClass      ), pointer :: thisCosmologyParameters

    ! Get the default cosmological parameters.
    thisCosmologyParameters        => cosmologyParameters()
    ! Use it to construct a matter plus cosmological constant cosmological functions class.
    matterLambdaDefaultConstructor = matterLambdaConstructor(thisCosmologyParameters)
    return
  end function matterLambdaDefaultConstructor

  function matterLambdaConstructor(thisCosmologyParameters)
    !% Constructor for the matter plus cosmological constant cosmological functions class.
    use Numerical_Comparison
    use Cosmology_Parameters
    use ISO_Varying_String
    use ODE_Solver
    use, intrinsic :: ISO_C_Binding
    implicit none
    type            (cosmologyFunctionsMatterLambda)               , target :: matterLambdaConstructor
    class           (cosmologyParametersClass      ), intent(in   ), target :: thisCosmologyParameters
    double precision                                , parameter             :: odeToleranceAbsolute   =1.0d-9, odeToleranceRelative   =1.0d-9
    double precision                                , parameter             :: omegaTolerance         =1.0d-9
    double precision                                                        :: OmegaDominant                 , aMaximum                      , &
         &                                                                     cubicTerm1                    , cubicTerm21                   , &
         &                                                                     cubicTerm21Squared            , cubicTerm25                   , &
         &                                                                     cubicTerm25Cubed              , cubicTerm5                    , &
         &                                                                     cubicTerm9                    , densityPower                  , &
         &                                                                     expansionFactorDominant       , timeMaximumimum     (1)
    type            (c_ptr                         )                        :: parameterPointer
    type            (fgsl_odeiv_step               )                        :: odeStepper
    type            (fgsl_odeiv_control            )                        :: odeController
    type            (fgsl_odeiv_evolve             )                        :: odeEvolver
    type            (fgsl_odeiv_system             )                        :: odeSystem
    logical                                                                 :: odeReset               =.true.

    ! Store a pointer to the cosmological parameters object.
    matterLambdaConstructor%cosmology => thisCosmologyParameters
    ! Determine if this universe will collapse. We take the Friedmann equation, which gives H^2 as a function of expansion
    ! factor, a, and solve for where H^2=0. If this has a real solution, then we have a collapsing universe.
    matterLambdaConstructor%collapsingUniverse    =.false.
    matterLambdaConstructor%expansionFactorMaximum=0.0d0
    matterLambdaConstructor%timeTurnaround        =0.0d0
    matterLambdaConstructor%timeMaximum           =0.0d0
    if    (Values_Agree(matterLambdaConstructor%cosmology%OmegaCurvature (),0.0d0,absTol=omegaTolerance)) then
       if (Values_Agree(matterLambdaConstructor%cosmology%OmegaDarkEnergy(),0.0d0,absTol=omegaTolerance)) then
          ! Einstein-de Sitter case. Always expands to infinity.
          matterLambdaConstructor%collapsingUniverse=.false.
       else
          ! Flat Universe with cosmological constant.
          if (matterLambdaConstructor%cosmology%OmegaDarkEnergy() > 0.0d0) then
             ! Never collapses.
             matterLambdaConstructor%collapsingUniverse=.false.
          else
             matterLambdaConstructor%collapsingUniverse=.true.
             matterLambdaConstructor%expansionFactorMaximum                    &
                  & =-(                                                        &
                  &     matterLambdaConstructor%cosmology%OmegaMatter    ()    &
                  &    *matterLambdaConstructor%cosmology%OmegaDarkEnergy()**2 &
                  &   )**(1.0d0/3.0d0)                                         &
                  &   / matterLambdaConstructor%cosmology%OmegaDarkEnergy()
          end if
       end if
    else
       if (Values_Agree(matterLambdaConstructor%cosmology%OmegaDarkEnergy(),0.0d0,absTol=omegaTolerance)) then
          ! Simple case for a matter-only universe.
          matterLambdaConstructor%collapsingUniverse=matterLambdaConstructor%cosmology%OmegaMatter() > 1.0d0
          if     ( matterLambdaConstructor%collapsingUniverse           ) &
               &   matterLambdaConstructor%expansionFactorMaximum         &
               & = matterLambdaConstructor%cosmology%OmegaMatter()        &
               & /(matterLambdaConstructor%cosmology%OmegaMatter()-1.0d0)
       else
          ! Case of matter plus dark energy.
          cubicTerm1        =1.0d0/matterLambdaConstructor%cosmology%OmegaDarkEnergy()
          cubicTerm5        =      matterLambdaConstructor%cosmology%OmegaMatter    ()**2
          cubicTerm9        =      matterLambdaConstructor%cosmology%OmegaDarkEnergy()**2
          cubicTerm21Squared                                                                                                                                 &
               & =-(                                                                                                                                         &
               &    -12.0d0                                                                                                                                  &
               &    +36.0d0           *matterLambdaConstructor%cosmology%OmegaMatter()                                                                       &
               &    +36.0d0                                                                             *matterLambdaConstructor%cosmology%OmegaDarkEnergy() &
               &    -36.0d0*cubicTerm5                                                                                                                       &
               &    -72.0d0           *matterLambdaConstructor%cosmology%OmegaMatter()*matterLambdaConstructor%cosmology%OmegaDarkEnergy()                   &
               &    -36.0d0*cubicTerm9                                                                                                                       &
               &    +12.0d0*cubicTerm5*matterLambdaConstructor%cosmology%OmegaMatter()                                                                       &
               &    -45.0d0*cubicTerm5                                                                  *matterLambdaConstructor%cosmology%OmegaDarkEnergy() &
               &    +36.0d0*cubicTerm9*matterLambdaConstructor%cosmology%OmegaMatter()                                                                       &
               &    +12.0d0*cubicTerm9                                                                  *matterLambdaConstructor%cosmology%OmegaDarkEnergy() &
               &   )&
               &  *cubicTerm1
          if (cubicTerm21Squared > 0.0d0) then
             cubicTerm21     =sqrt(cubicTerm21Squared)
             cubicTerm25Cubed=(                                                         &
                  &            -108.0d0*matterLambdaConstructor%cosmology%OmegaMatter() &
                  &            + 12.0d0*cubicTerm21                                     &
                  &           )                                                         &
                  &           *cubicTerm9
             if (cubicTerm25Cubed >= 0.0d0) then
                cubicTerm25=     cubicTerm25Cubed **(1.0d0/3.0d0)
             else
                cubicTerm25=-abs(cubicTerm25Cubed)**(1.0d0/3.0d0)
             end if
             aMaximum                                                        &
                  &  = cubicTerm1                                            &
                  &   *cubicTerm25                                           &
                  &   /6.0d0                                                 &
                  &   +2.0d0                                                 &
                  &   *(                                                     &
                  &     -1.0d0                                               &
                  &     +matterLambdaConstructor%cosmology%OmegaMatter    () &
                  &     +matterLambdaConstructor%cosmology%OmegaDarkEnergy() &
                  &   )                                                      &
                  &   /cubicTerm25
             matterLambdaConstructor%collapsingUniverse=aMaximum > 0.0d0
             if (matterLambdaConstructor%collapsingUniverse) &
                  & matterLambdaConstructor%expansionFactorMaximum=aMaximum
          end if
       end if
    end if
    ! If we have a collapsing Universe, find time of turnaround, and maximum time.
    if (matterLambdaConstructor%collapsingUniverse) then
       ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
       call matterLambdaConstructor%densityScalingEarlyTime(matterLambdaDominateFactor,densityPower,expansionFactorDominant,OmegaDominant)
       ! Find the corresponding time.
       timeMaximumimum(1)=1.0d0/matterLambdaConstructor%cosmology%HubbleConstant(unitsTime)/sqrt(OmegaDominant)/expansionFactorDominant**(0.5d0*densityPower)
       ! Solve Friedmann equation to get time at turnaround.
       matterLambdaSelfGlobal => matterLambdaConstructor
       odeReset=.true.
       call ODE_Solve(                                                               &
            &         odeStepper                                                   , &
            &         odeController                                                , &
            &         odeEvolver                                                   , &
            &         odeSystem                                                    , &
            &         expansionFactorDominant                                      , &
            &         matterLambdaConstructor%expansionFactorMaximum*(1.0d0-1.0d-4), &
            &         1                                                            , &
            &         timeMaximumimum                                              , &
            &         matterLambdaCollapseODEs                                     , &
            &         parameterPointer                                             , &
            &         odeToleranceAbsolute                                         , &
            &         odeToleranceRelative                                         , &
            &         reset=odeReset                                                 &
            &        )
       call ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
       odeReset=.true.
       ! Extract turnaround time from ODE variables and set maximum time to twice turnaround time.
       matterLambdaConstructor%timeTurnaround=timeMaximumimum(1)
       matterLambdaConstructor%timeMaximum       =2.0d0*matterLambdaConstructor%timeTurnaround
    end if
    return
  end function matterLambdaConstructor

  subroutine matterLambdaDestructor(self)
    !% Default constructor for the matter plus cosmological constant cosmological functions class.
    use Numerical_Interpolation
    implicit none
    type(cosmologyFunctionsMatterLambda), intent(inout) :: self

    if (allocated(self%ageTableExpansionFactor               )) deallocate(self%ageTableExpansionFactor               )
    if (allocated(self%ageTableTime                          )) deallocate(self%ageTableTime                          )
    if (allocated(self%distanceTableComovingDistance         )) deallocate(self%distanceTableComovingDistance         )
    if (allocated(self%distanceTableComovingDistanceNegated  )) deallocate(self%distanceTableComovingDistanceNegated  )
    if (allocated(self%distanceTableLuminosityDistanceNegated)) deallocate(self%distanceTableLuminosityDistanceNegated)
    if (allocated(self%distanceTableTime                     )) deallocate(self%distanceTableTime                     )
    call Interpolate_Done(self%interpolationObject       ,self%interpolationAccelerator       ,self%resetInterpolation       )
    call Interpolate_Done(self%interpolationObjectInverse,self%interpolationAcceleratorInverse,self%resetInterpolationInverse)
    return
  end subroutine matterLambdaDestructor

  double precision function matterLambdaCosmicTime(self,expansionFactor,collapsingPhase)
    !% Return the cosmological matter density in units of the critical density at the present day.
    use Galacticus_Error
    use Numerical_Interpolation
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   )           :: expansionFactor
    logical                                         , intent(in   ), optional :: collapsingPhase
    logical                                                                   :: collapsingPhaseActual

    ! Validate the input.
    if (.not.self%expansionFactorIsValid(expansionFactor)) &
         & call Galacticus_Error_Report('matterLambdaCosmicTime','expansion factor is invalid')
    ! Determine if we are in the expanding or collapsing phase for this universe.
    if (present(collapsingPhase)) then
       collapsingPhaseActual=collapsingPhase
    else
       collapsingPhaseActual=.false. ! Assume expanding phase by default.
    end if
    ! Ensure tabulation is initialized.
    if (.not.self%ageTableInitialized) call self%expansionFactorTabulate(self%ageTableTimeMinimum)
    ! Ensure that the tabulation spans a sufficient range of expansion factors.
    if (collapsingPhaseActual) then
       ! In collapsing phase just ensure that a sufficiently large expansion factor has been reached.
       do while (self%ageTableExpansionFactor(self%ageTableNumberPoints) < expansionFactor)
          self%ageTableTimeMaximum=min(self%ageTableTimeMaximum*matterLambdaAgeTableIncrementFactor,self%timeTurnaround)
          call self%expansionFactorTabulate()
       end do
    else
       ! In expanding phase ensure that sufficiently small and large expansion factors have been reached.
       do while (self%ageTableExpansionFactor(                    1) > expansionFactor)
          self%ageTableTimeMinimum=    self%ageTableTimeMinimum/matterLambdaAgeTableIncrementFactor
          call self%expansionFactorTabulate()
       end do
       do while (self%ageTableExpansionFactor(self%iTableTurnaround) < expansionFactor)
          self%ageTableTimeMaximum=max(self%ageTableTimeMaximum*matterLambdaAgeTableIncrementFactor,self%timeTurnaround)
          call self%expansionFactorTabulate()
       end do
    end if
    ! Interpolate to get cosmic time.
    matterLambdaCosmicTime                             &
         & =Interpolate(                               &
         &              self%ageTableNumberPoints    , &
         &              self%ageTableExpansionFactor , &
         &              self%ageTableTime            , &
         &              self%interpolationObject     , &
         &              self%interpolationAccelerator, &
         &              expansionFactor              , &
         &              reset=self%resetInterpolation  &
         &             )
    if (collapsingPhaseActual) matterLambdaCosmicTime=self%timeMaximum-matterLambdaCosmicTime
    return
  end function matterLambdaCosmicTime

  function matterLambdaCollapseODEs(a,t,dtda,parameterPointer) bind(c)
    !% System of differential equations to solve for age vs. expansion factor.
    use, intrinsic :: ISO_C_Binding
    implicit none
    integer(kind=c_int   )                              :: matterLambdaCollapseODEs
    real   (kind=c_double)              , value         :: a
    real   (kind=c_double), dimension(1), intent(in   ) :: t
    real   (kind=c_double), dimension(1)                :: dtda
    type   (c_ptr        )              , value         :: parameterPointer

    dtda(1)=1.0d0/a/matterLambdaSelfGlobal%expansionRate(a)
    matterLambdaCollapseODEs=FGSL_Success
    return
  end function matterLambdaCollapseODEs

  logical function matterLambdaExpansionFactorIsValid(self,expansionFactor)
    !% Checks that the expansion factor falls within allowed ranges.
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: expansionFactor

    matterLambdaExpansionFactorIsValid=expansionFactor > 0.0d0 .and. (expansionFactor < self%expansionFactorMaximum .or. .not.self%collapsingUniverse)
    return
  end function matterLambdaExpansionFactorIsValid

  logical function matterLambdaCosmicTimeIsValid(self,time)
    !% Checks that the time falls within allowed ranges.
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: time

    matterLambdaCosmicTimeIsValid=time > 0.0d0 .and. (time < self%timeMaximum .or. .not.self%collapsingUniverse)
    return
  end function matterLambdaCosmicTimeIsValid

  double precision function matterLambdaExpansionFactor(self,time)
    !% Returns the expansion factor at cosmological time {\tt time}.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    double precision                                                :: timeEffective
    logical                                                         :: remakeTable

    ! Check if the time differs from the previous time.
    if (time /= self%timePrevious) then
       ! Quit on invalid input.
       if (time < 0.0d0) call Galacticus_Error_Report('matterLambdaExpansionFactor','cosmological time must be positive')
       ! Check if we need to recompute our table.
       if (self%ageTableInitialized) then
          remakeTable=(time < self%ageTableTime(1).or.time > self%ageTableTime(self%ageTableNumberPoints))
       else
          remakeTable=.true.
       end if
       if (remakeTable)           call self%expansionFactorTabulate(time)
       ! Quit on invalid input.
       if (self%collapsingUniverse.and.time > self%timeMaximum) &
            & call Galacticus_Error_Report('matterLambdaExpansionFactor','cosmological time exceeds that at the Big Crunch')
       ! Interpolate to get the expansion factor.
       if (self%collapsingUniverse) then
          if (time <= self%timeTurnaround) then
             timeEffective=                      time
          else
             timeEffective=self%timeMaximum-time
          end if
       else
          timeEffective   =                      time
       end if
       self%expansionFactorPrevious                              &
            & =Interpolate(                                      &
            &              self%ageTableNumberPoints           , &
            &              self%ageTableTime                   , &
            &              self%ageTableExpansionFactor        , &
            &              self%interpolationObjectInverse     , &
            &              self%interpolationAcceleratorInverse, &
            &              timeEffective                       , &
            &              reset=self%resetInterpolationInverse  &
            &             )
       self%timePrevious=time
    end if
    ! Return the stored expansion factor.
    matterLambdaExpansionFactor=self%expansionFactorPrevious
    return
  end function matterLambdaExpansionFactor

  double precision function matterLambdaExpansionRate(self,expansionFactor)
    !% Returns the cosmological expansion rate, $\dot{a}/a$ at expansion factor {\tt expansionFactor}.
    use Cosmology_Parameters
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: expansionFactor

    ! Required value is simply the Hubble parameter but expressed in units of inverse Gyr.
    matterLambdaExpansionRate&
         & = self          %hubbleParameterEpochal(expansionFactor=expansionFactor) &
         &  *self%cosmology%HubbleConstant        (unitsTime                      ) &
         &  /self%cosmology%HubbleConstant        (unitsStandard                  )
    return
  end function matterLambdaExpansionRate

  double precision function matterLambdaHubbleParameterEpochal(self,time,expansionFactor,collapsingPhase)
    !% Returns the Hubble parameter at the request cosmological time, {\tt time}, or expansion factor, {\tt expansionFactor}.
    use Galacticus_Error
    use Cosmology_Parameters
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual, sqrtArgument

    ! Determine the actual expansion factor to use.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('matterLambdaHubbleParameterEpochal','only one of time or expansion factor can be specified')
       else
          expansionFactorActual=matterLambdaSelfGlobal%expansionFactor(time)
       end if
    else
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else
          call Galacticus_Error_Report('matterLambdaHubbleParameterEpochal','either a time or expansion factor must be specified')
       end if
    end if
    ! Compute the Hubble parameter at the specified expansion factor.
    sqrtArgument                      =max(                                   &
         &                                  self%cosmology%OmegaMatter    ()  &
         &                                 /expansionFactorActual**3          &
         &                                 +self%cosmology%OmegaDarkEnergy()  &
         &                                 +self%cosmology%OmegaCurvature ()  &
         &                                 /expansionFactorActual**2        , &
         &                                 0.0d0                              &
         &                                )
    matterLambdaHubbleParameterEpochal= self%cosmology%HubbleConstant(unitsStandard)&
         &                             *sqrt(sqrtArgument)
    ! Make the Hubble parameter negative if we are in the collapsing phase of the Universe.
    if (self%collapsingUniverse) then
       if    (present(time           )) then
          if    (time>self%timeTurnaround) matterLambdaHubbleParameterEpochal=-matterLambdaHubbleParameterEpochal
       else
          if (present(collapsingPhase)) then
             if (collapsingPhase         ) matterLambdaHubbleParameterEpochal=-matterLambdaHubbleParameterEpochal
          end if
       end if
    end if
    return
  end function matterLambdaHubbleParameterEpochal

  double precision function matterLambdaOmegaMatterEpochal(self,time,expansionFactor,collapsingPhase)
    !% Return the matter density parameter at expansion factor {\tt expansionFactor}.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual

    ! Determine the actual expansion factor to use.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('matterLambdaOmegaMatterEpochal','only one of time or expansion factor can be specified')
       else
          expansionFactorActual=self%expansionFactor(time)
       end if
    else
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else
          call Galacticus_Error_Report('matterLambdaOmegaMatterEpochal','either a time or expansion factor must be specified')
       end if
    end if
    matterLambdaOmegaMatterEpochal                                                          &
         & =   self%cosmology%OmegaMatter           (                                     ) &
         &  *(                                                                              &
         &     self%cosmology%HubbleConstant        (unitsStandard                        ) &
         &    /self          %HubbleParameterEpochal(expansionFactor=expansionFactorActual) &
         &   )**2                                                                           &
         &  /expansionFactorActual**3
    return
  end function matterLambdaOmegaMatterEpochal

  double precision function matterLambdaOmegaDarkEnergyEpochal(self,time,expansionFactor,collapsingPhase)
    !% Return the dark energy density parameter at expansion factor {\tt expansionFactor}.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual

    ! Determine the actual expansion factor to use.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('matterLambdaOmegaDarkEnergyEpochal','only one of time or expansion factor can be specified')
       else
          expansionFactorActual=self%expansionFactor(time)
       end if
    else
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else
          call Galacticus_Error_Report('matterLambdaOmegaDarkEnergyEpochal','either a time or expansion factor must be specified')
       end if
    end if
    matterLambdaOmegaDarkEnergyEpochal                                                      &
         & =   self%cosmology%OmegaDarkEnergy       (                                     ) &
         &  *(                                                                              &
         &     self%cosmology%HubbleConstant        (unitsStandard                        ) &
         &    /self%          HubbleParameterEpochal(expansionFactor=expansionFactorActual) &
         &   )**2
    return
  end function matterLambdaOmegaDarkEnergyEpochal

  double precision function matterLambdaTemperatureCMBEpochal(self,time,expansionFactor,collapsingPhase)
    !% Return the temperature of the CMB at expansion factor {\tt expansionFactor}.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor      , time
    logical                                         , intent(in   ), optional :: collapsingPhase
    double precision                                                          :: expansionFactorActual

    ! Determine the actual expansion factor to use.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('matterLambdaTemperatureCMBEpochal','only one of time or expansion factor can be specified')
       else
          expansionFactorActual=self%expansionFactor(time)
       end if
    else
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else
          call Galacticus_Error_Report('matterLambdaTemperatureCMBEpochal','either a time or expansion factor must be specified')
       end if
    end if
    matterLambdaTemperatureCMBEpochal=self%cosmology%temperatureCMB()/expansionFactorActual
    return
  end function matterLambdaTemperatureCMBEpochal

  subroutine matterLambdaDensityScalingEarlyTime(self,dominateFactor,densityPower,expansionFactorDominant,OmegaDominant)
    use Cosmology_Parameters
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   )           :: dominateFactor
    double precision                                , intent(  out)           :: densityPower  , expansionFactorDominant
    double precision                                , intent(  out), optional :: OmegaDominant

    ! For matter and cosmological constant, matter always dominates at early times.
    densityPower=-3.0d0 ! Power-law scaling of matter density with expansion factor.
    ! Choose present day as default - will be used if no other densities present (i.e. Einsetin-de Sitter).
    expansionFactorDominant=self%dominationEpochMatter(dominateFactor)
    ! Return the density parameter in the dominant species if required.
    if (present(OmegaDominant)) OmegaDominant=self%cosmology%OmegaMatter()
    return
  end subroutine matterLambdaDensityScalingEarlyTime

  double precision function matterLambdaDominationEpochMatter(self,dominateFactor)
    use Cosmology_Functions_Parameters
    use Cosmology_Parameters
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: dominateFactor
    double precision                                                :: aMatterEquality                  , expansionFactorDominantCurvature, &
         &                                                             expansionFactorDominantDarkEnergy

    ! Choose present day as default - will be used if no other densities present (i.e. Einsetin-de Sitter).
    matterLambdaDominationEpochMatter=1.0d0
    if (self%cosmology%OmegaDarkEnergy()/=0.0d0) then
       ! Find the expansion factor of matter-dark energy equality.
       aMatterEquality=self%equalityEpochMatterDarkEnergy(requestTypeExpansionFactor)
       ! Find the earlier expansion factor at which matter dominates by the specified amount (ratio of matter
       ! to dark energy density scales as the cube of expansion factor).
       expansionFactorDominantDarkEnergy=aMatterEquality/dominateFactor**(1.0d0/3.0d0)
       ! Choose earliest expansion factor.
       matterLambdaDominationEpochMatter=min(matterLambdaDominationEpochMatter,expansionFactorDominantDarkEnergy)
    end if
    if (self%cosmology%OmegaCurvature() /= 0.0d0) then
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
    !% Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    use Cosmology_Functions_Parameters
    use Cosmology_Parameters
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
         &        self%cosmology%OmegaMatter    ()  &
         &   /abs(self%cosmology%OmegaDarkEnergy()) &
         &  )**(1.0d0/3.0d0)
    if (requestType == requestTypeTime)                                &
         &                   matterLambdaEqualityEpochMatterDarkEnergy &
         &  =self%cosmicTime(matterLambdaEqualityEpochMatterDarkEnergy)
    return
  end function matterLambdaEqualityEpochMatterDarkEnergy

  double precision function matterLambdaEqualityEpochMatterCurvature(self,requestType)
    !% Return the epoch of matter-curvature magnitude equality (either expansion factor or cosmic time).
    use Cosmology_Functions_Parameters
    use Cosmology_Parameters
    implicit none
    class  (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    integer                                , intent(in   ), optional :: requestType
    integer                                                          :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if

    matterLambdaEqualityEpochMatterCurvature=self%cosmology%OmegaMatter()/abs(self%cosmology%OmegaCurvature())
    if (requestType == requestTypeTime)                               &
         &                  matterLambdaEqualityEpochMatterCurvature  &
         & =self%cosmicTime(matterLambdaEqualityEpochMatterCurvature)
    return
  end function matterLambdaEqualityEpochMatterCurvature

  subroutine matterLambdaMakeExpansionFactorTable(self,time)
    !% Builds a table of expansion factor vs. time.
    use Numerical_Interpolation
    use Numerical_Ranges
    use ODE_Solver
    use Memory_Management
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Parameters
    use FGSL
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout), target       :: self
    double precision                                , intent(in   ), optional     :: time
    double precision                                , parameter                   :: odeToleranceAbsolute               =1.0d-9, odeToleranceRelative   =1.0d-9
    double precision                                , allocatable  , dimension(:) :: ageTableExpansionFactorTemporary          , ageTableTimeTemporary
    integer                                                                       :: iTime                                     , prefixPointCount
    double precision                                                              :: OmegaDominant                             , densityPower                  , &
         &                                                                           expansionFactor                 (1)       , expansionFactorDominant       , &
         &                                                                           tDominant                                 , timeActual
    type            (c_ptr                         )                              :: parameterPointer
    type            (fgsl_odeiv_step               )                              :: odeStepper
    type            (fgsl_odeiv_control            )                              :: odeController
    type            (fgsl_odeiv_evolve             )                              :: odeEvolver
    type            (fgsl_odeiv_system             )                              :: odeSystem
    logical                                                                       :: odeReset                           =.true.

    ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
    call self%densityScalingEarlyTime(matterLambdaDominateFactor,densityPower,expansionFactorDominant,OmegaDominant)
    ! Find the corresponding time.
    tDominant=-2.0d0/densityPower/self%cosmology%HubbleConstant(unitsTime)/sqrt(OmegaDominant)/expansionFactorDominant**(0.5d0*densityPower)
    ! Find minimum and maximum times to tabulate.
    if (present(time)) then
       timeActual=time
       do while (self%ageTableTimeMinimum > min(timeActual,tDominant)/2.0d0)
          self%ageTableTimeMinimum=self%ageTableTimeMinimum/matterLambdaAgeTableIncrementFactor
       end do
       do while (self%ageTableTimeMaximum < max(timeActual,tDominant)*2.0d0)
          self%ageTableTimeMaximum=self%ageTableTimeMaximum*matterLambdaAgeTableIncrementFactor
       end do
    else
       do while (self%ageTableTimeMinimum > tDominant/2.0d0)
          self%ageTableTimeMinimum=self%ageTableTimeMinimum/matterLambdaAgeTableIncrementFactor
       end do
       do while (self%ageTableTimeMaximum < tDominant*2.0d0)
          self%ageTableTimeMaximum=self%ageTableTimeMaximum*matterLambdaAgeTableIncrementFactor
       end do
    end if
    if (self%collapsingUniverse) self%ageTableTimeMaximum=min(self%ageTableTimeMaximum,self%timeTurnaround)

    ! Determine number of points to tabulate.
    self%ageTableNumberPoints=int(log10(self%ageTableTimeMaximum/self%ageTableTimeMinimum)*dble(matterLambdaAgeTableNPointsPerDecade))+1
    self%ageTableTimeMaximum =self%ageTableTimeMinimum*10.0d0**(dble(self%ageTableNumberPoints)/dble(matterLambdaAgeTableNPointsPerDecade))
    if (self%collapsingUniverse) self%ageTableTimeMaximum=min(self%ageTableTimeMaximum,self%timeTurnaround)

    ! Deallocate arrays if currently allocated.
    if (allocated(self%ageTableTime)) then
       ! Determine number of points that are being added at the start of the array.
       prefixPointCount=int(log10(self%ageTableTime(1)/self%ageTableTimeMinimum)*dble(matterLambdaAgeTableNPointsPerDecade)+0.5d0)
       call Move_Alloc(self%ageTableTime           ,ageTableTimeTemporary           )
       call Move_Alloc(self%ageTableExpansionFactor,ageTableExpansionFactorTemporary)
       ! Allocate the arrays to current required size.
       call Alloc_Array(self%ageTableTime,           [self%ageTableNumberPoints])
       call Alloc_Array(self%ageTableExpansionFactor,[self%ageTableNumberPoints])
       ! Create set of grid points in time variable.
       self%ageTableTime=Make_Range(self%ageTableTimeMinimum,self%ageTableTimeMaximum,self%ageTableNumberPoints,rangeTypeLogarithmic)
       ! Set the expansion factors to a negative value to indicate they are not yet computed.
       self%ageTableExpansionFactor=-1.0d0
       ! Paste in the previously computed regions.
       self%ageTableTime           (prefixPointCount+1:prefixPointCount+size(ageTableTimeTemporary))=ageTableTimeTemporary
       self%ageTableExpansionFactor(prefixPointCount+1:prefixPointCount+size(ageTableTimeTemporary))=ageTableExpansionFactorTemporary
       ! Deallocate the temporary arrays.
       call Dealloc_Array(ageTableTimeTemporary           )
       call Dealloc_Array(ageTableExpansionFactorTemporary)
    else
       ! Allocate the arrays to current required size.
       call Alloc_Array(self%ageTableTime,           [self%ageTableNumberPoints])
       call Alloc_Array(self%ageTableExpansionFactor,[self%ageTableNumberPoints])
       ! Create set of grid points in time variable.
       self%ageTableTime=Make_Range(self%ageTableTimeMinimum,self%ageTableTimeMaximum,self%ageTableNumberPoints,rangeTypeLogarithmic)
       ! Set the expansion factors to a negative value to indicate they are not yet computed.
       self%ageTableExpansionFactor=-1.0d0
    end if
    ! For the initial time, we approximate that we are at sufficiently early times that a single component dominates the
    ! Universe and use the appropriate analytic solution.
    if (self%ageTableExpansionFactor(1) < 0.0d0)       &
         &    self%ageTableExpansionFactor(         1) &
         & =(                                          &
         &   -0.5d0                                    &
         &   *densityPower                             &
         &   *self%ageTableTime            (        1) &
         &   *self%cosmology%HubbleConstant(unitsTime) &
         &   *sqrt(OmegaDominant)                      &
         &  )**(-2.0d0/densityPower)
    ! Solve ODE to get corresponding expansion factors.
    self%iTableTurnaround  =  self%ageTableNumberPoints
    matterLambdaSelfGlobal => self
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
          call ODE_Solve(                          &
               &         odeStepper              , &
               &         odeController           , &
               &         odeEvolver              , &
               &         odeSystem               , &
               &         timeActual              , &
               &         self%ageTableTime(iTime), &
               &         1                       , &
               &         expansionFactor         , &
               &         matterLambdaAgeTableODEs, &
               &         parameterPointer        , &
               &         odeToleranceAbsolute    , &
               &         odeToleranceRelative    , &
               &         reset=odeReset            &
               &        )
          self%ageTableExpansionFactor(iTime)=expansionFactor(1)
       end if
    end do
    call ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
    call Interpolate_Done(self%interpolationObject       ,self%interpolationAccelerator       ,self%resetInterpolation       )
    call Interpolate_Done(self%interpolationObjectInverse,self%interpolationAcceleratorInverse,self%resetInterpolationInverse)
    self%resetInterpolation       =.true.
    self%resetInterpolationInverse=.true.
    ! Flag that the table is now initialized.
    self%ageTableInitialized=.true.
    return
  end subroutine matterLambdaMakeExpansionFactorTable

  function matterLambdaAgeTableODEs(t,a,dadt,parameterPointer) bind(c)
    !% System of differential equations to solve for expansion factor vs. age.
    use, intrinsic :: ISO_C_Binding
    integer(kind=c_int   )                              :: matterLambdaAgeTableODEs
    real   (kind=c_double)              , value         :: t
    real   (kind=c_double), dimension(1), intent(in   ) :: a
    real   (kind=c_double), dimension(1)                :: dadt
    type   (c_ptr        )              , value         :: parameterPointer

    dadt(1)=a(1)*matterLambdaSelfGlobal%expansionRate(a(1))
    matterLambdaAgeTableODEs=FGSL_Success
  end function matterLambdaAgeTableODEs

  double precision function matterLambdaTimeAtDistanceComoving(self,comovingDistance)
    !% Returns the cosmological time corresponding to given {\tt comovingDistance}.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: comovingDistance
    double precision                                                :: time
    logical                                                         :: remakeTable

    ! Quit on invalid input.
    if (comovingDistance < 0.0d0) call Galacticus_Error_Report('matterLambdaTimeFromComovingDistance','comoving distance must be positive')
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
    matterLambdaTimeAtDistanceComoving                                &
         & =Interpolate(                                              &
         &              self%distanceTableNumberPoints              , &
         &              self%distanceTableComovingDistanceNegated   , &
         &              self%distanceTableTime                      , &
         &              self%interpolationObjectDistanceInverse     , &
         &              self%interpolationAcceleratorDistanceInverse, &
         &              -comovingDistance                           , &
         &              reset=self%resetInterpolationDistanceInverse  &
         &             )
    return
  end function matterLambdaTimeAtDistanceComoving

  double precision function matterLambdaDistanceComoving(self,time)
    !% Returns the comoving distance to cosmological time {\tt time}.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    logical                                                         :: remakeTable

    ! Quit on invalid input.
    if (time <                 0.0d0 ) call Galacticus_Error_Report('matterLambdaComovingDistance','cosmological time must be positive'   )
    if (time > self%cosmicTime(1.0d0)) call Galacticus_Error_Report('matterLambdaComovingDistance','cosmological time must be in the past')
    ! Check if we need to recompute our table.
    if (self%distanceTableInitialized) then
       remakeTable=(time < self%distanceTableTime(1).or.time > self%distanceTableTime(self%distanceTableNumberPoints))
    else
       remakeTable=.true.
    end if
    if (remakeTable) call self%distanceTabulate(time)
    ! Quit on invalid input.
    if (self%collapsingUniverse.and.time>self%timeMaximum) &
         & call Galacticus_Error_Report('matterLambdaComovingDistance','cosmological time exceeds that at the Big Crunch')
    ! Interpolate to get the comoving distance.
    matterLambdaDistanceComoving                               &
         & =Interpolate(                                       &
         &              self%distanceTableNumberPoints       , &
         &              self%distanceTableTime               , &
         &              self%distanceTableComovingDistance   , &
         &              self%interpolationObjectDistance     , &
         &              self%interpolationAcceleratorDistance, &
         &              time                                 , &
         &              reset=self%resetInterpolationDistance  &
         &             )
    return
   end function matterLambdaDistanceComoving

  double precision function matterLambdaDistanceComovingConvert(self,output,distanceModulus,redshift)
    !% Convert bewteen different measures of distance.
    use Numerical_Interpolation
    use Galacticus_Error
    use Cosmology_Functions_Options
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    integer                                         , intent(in   )           :: output
    double precision                                , intent(in   ), optional :: distanceModulus    , redshift
    logical                                                                   :: gotComovingDistance
    double precision                                                          :: comovingDistance   , luminosityDistance

    ! Check if we need to recompute our table.
    if (.not.self%distanceTableInitialized) call self%distanceTabulate(self%cosmicTime(1.0d0))
    ! Convert to comoving distance from whatever was supplied.
    gotComovingDistance=.false.
    if (present(distanceModulus)) then
       luminosityDistance=10.0d0**((distanceModulus-25.0d0)/5.0d0)
       do while (luminosityDistance > -self%distanceTableLuminosityDistanceNegated(1))
          call self%distanceTabulate(0.5d0*self%distanceTableTimeMinimum)
       end do
       comovingDistance                                                      &
            & =-Interpolate(                                                 &
            &               self%distanceTableNumberPoints                 , &
            &               self%distanceTableLuminosityDistanceNegated    , &
            &               self%distanceTableComovingDistanceNegated      , &
            &               self%interpolationObjectLuminosityDistance     , &
            &               self%interpolationAcceleratorLuminosityDistance, &
            &               -luminosityDistance                            , &
            &               reset=self%resetInterpolationLuminosityDistance  &
            &              )
       gotComovingDistance=.true.
    end if
    if (present(redshift)) then
       comovingDistance   =self%distanceComoving(self%cosmicTime(self%expansionFactorFromRedshift(redshift)))
       gotComovingDistance=.true.
    end if
    if (.not.gotComovingDistance) &
         & call Galacticus_Error_Report('matterLambdaComovingDistanceConvert','no distance measure provided')
    ! Convert to required distance measure.
    select case (output)
    case (distanceTypeComoving)
       matterLambdaDistanceComovingConvert=comovingDistance
    case default
       call Galacticus_Error_Report('matterLambdaComovingDistanceConvert','unrecognized output option')
    end select
    return
  end function matterLambdaDistanceComovingConvert

  subroutine matterLambdaMakeDistanceTable(self,time)
    !% Builds a table of distance vs. time.
    use Numerical_Interpolation
    use Numerical_Ranges
    use Numerical_Integration
    use Memory_Management
    use, intrinsic :: ISO_C_Binding
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout), target :: self
    double precision                                , intent(in   )         :: time
    double precision                                , parameter             :: toleranceAbsolute   =1.0d-5, toleranceRelative=1.0d-5
    integer                                                                 :: iTime
    logical                                                                 :: resetIntegration
    type            (c_ptr                         )                        :: parameterPointer
    type            (fgsl_function                 )                        :: integrandFunction
    type            (fgsl_integration_workspace    )                        :: integrationWorkspace

    ! Find minimum and maximum times to tabulate.
    self%distanceTableTimeMinimum=min(self%distanceTableTimeMinimum,0.5d0*time)
    self%distanceTableTimeMaximum=    self%cosmicTime(1.0d0)
    ! Determine number of points to tabulate.
    self%distanceTableNumberPoints=int(log10(self%distanceTableTimeMaximum/self%distanceTableTimeMinimum)*dble(matterLambdaDistanceTableNPointsPerDecade))+1
    ! Deallocate arrays if currently allocated.
    if (allocated(self%distanceTableTime                     )) call Dealloc_Array(self%distanceTableTime                     )
    if (allocated(self%distanceTableComovingDistance         )) call Dealloc_Array(self%distanceTableComovingDistance         )
    if (allocated(self%distanceTableComovingDistanceNegated  )) call Dealloc_Array(self%distanceTableComovingDistanceNegated  )
    if (allocated(self%distanceTableLuminosityDistanceNegated)) call Dealloc_Array(self%distanceTableLuminosityDistanceNegated)
    ! Allocate the arrays to current required size.
    call Alloc_Array(self%distanceTableTime                     ,[self%distanceTableNumberPoints])
    call Alloc_Array(self%distanceTableComovingDistance         ,[self%distanceTableNumberPoints])
    call Alloc_Array(self%distanceTableComovingDistanceNegated  ,[self%distanceTableNumberPoints])
    call Alloc_Array(self%distanceTableLuminosityDistanceNegated,[self%distanceTableNumberPoints])
    ! Create the range of times.
    self% distanceTableTime=Make_Range(self%distanceTableTimeMinimum,self%distanceTableTimeMaximum,self%distanceTableNumberPoints,rangeTypeLogarithmic)
    ! Integrate to get the comoving distance.
    resetIntegration       =  .true.
    matterLambdaSelfGlobal => self
    do iTime=1,self%distanceTableNumberPoints
       self%distanceTableComovingDistance(iTime)                                 &
            & =Integrate(                                                        &
            &            self%distanceTableTime(iTime)                         , &
            &            self%distanceTableTime(self%distanceTableNumberPoints), &
            &            matterLambdaComovingDistanceIntegrand                 , &
            &            parameterPointer                                      , &
            &            integrandFunction                                     , &
            &            integrationWorkspace                                  , &
            &            toleranceAbsolute=toleranceAbsolute                   , &
            &            toleranceRelative=toleranceRelative                   , &
            &            reset=resetIntegration                                  &
            &           )
       self         %distanceTableLuminosityDistanceNegated              (iTime)  &
            & = self%distanceTableComovingDistance                       (iTime)  &
            &  /self%expansionFactor              (self%distanceTableTime(iTime))
    end do
    ! Make a negated copy of the distances so that we have an increasing array for use in interpolation routines.
    self%distanceTableComovingDistanceNegated  =-self%distanceTableComovingDistance
    self%distanceTableLuminosityDistanceNegated=-self%distanceTableLuminosityDistanceNegated
    ! Reset interpolators.
    call Interpolate_Done(self%interpolationObjectDistance          ,self%interpolationAcceleratorDistance          ,self%resetInterpolationDistance          )
    call Interpolate_Done(self%interpolationObjectDistanceInverse   ,self%interpolationAcceleratorDistanceInverse   ,self%resetInterpolationDistanceInverse   )
    call Interpolate_Done(self%interpolationObjectLuminosityDistance,self%interpolationAcceleratorLuminosityDistance,self%resetInterpolationLuminosityDistance)
    self%resetInterpolationDistance          =.true.
    self%resetInterpolationDistanceInverse   =.true.
    self%resetInterpolationLuminosityDistance=.true.
    ! Flag that the table is now initialized.
    self%distanceTableInitialized=.true.
    return
  end subroutine matterLambdaMakeDistanceTable

  function matterLambdaComovingDistanceIntegrand(time,parameterPointer) bind(c)
    !% Integrand function used in computing the comoving distance.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(kind=c_double)        :: matterLambdaComovingDistanceIntegrand
    real(kind=c_double), value :: time
    type(c_ptr        ), value :: parameterPointer

    matterLambdaComovingDistanceIntegrand=speedLight*gigaYear/megaParsec/matterLambdaSelfGlobal%expansionFactor(time)
    return
  end function matterLambdaComovingDistanceIntegrand

  double precision function matterLambdaEquationOfStateDarkEnergy(self,time,expansionFactor)
    !% Return the dark energy equation of state.
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor, time

    matterLambdaEquationOfStateDarkEnergy=-1.0d0
    return
  end function matterLambdaEquationOfStateDarkEnergy

  double precision function matterLambdaExponentDarkEnergy(self,time,expansionFactor)
    !% Return the dark energy equation of state.
    implicit none
    class           (cosmologyFunctionsMatterLambda), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: expansionFactor, time

    matterLambdaExponentDarkEnergy=0.0d0
    return
  end function matterLambdaExponentDarkEnergy

  subroutine matterLambdaStateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    class  (cosmologyFunctionsMatterLambda), intent(inout) :: self
    integer                                , intent(in   ) :: stateFile
    type   (fgsl_file                     ), intent(in   ) :: fgslStateFile

    ! Store the full tables, as they are hysteretic and cannot be reconstructed precisely without
    ! knowing the path by which they were originally constructed.
    write (stateFile) self%ageTableNumberPoints,self%ageTableTimeMinimum,self%ageTableTimeMaximum
    write (stateFile) self%ageTableTime,self%ageTableExpansionFactor
    write (stateFile) self%distanceTableNumberPoints,self%distanceTableTimeMinimum,self%distanceTableTimeMaximum
    write (stateFile) self%distanceTableTime,self%distanceTableComovingDistance,self%distanceTableComovingDistanceNegated
    return
  end subroutine matterLambdaStateStore

  subroutine matterLambdaStateRestore(self,stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use Memory_Management
    implicit none
    class  (cosmologyFunctionsMatterLambda), intent(inout) :: self
    integer                                , intent(in   ) :: stateFile
    type   (fgsl_file                     ), intent(in   ) :: fgslStateFile

    ! Read the tabulations.
    read (stateFile) self%ageTableNumberPoints,self%ageTableTimeMinimum,self%ageTableTimeMaximum
    if (allocated(self%ageTableTime           )) call Dealloc_Array(self%ageTableTime           )
    if (allocated(self%ageTableExpansionFactor)) call Dealloc_Array(self%ageTableExpansionFactor)
    call Alloc_Array(self%ageTableTime           ,[self%ageTableNumberPoints])
    call Alloc_Array(self%ageTableExpansionFactor,[self%ageTableNumberPoints])
    read (stateFile) self%ageTableTime,self%ageTableExpansionFactor
    read (stateFile) self%distanceTableNumberPoints,self%distanceTableTimeMinimum,self%distanceTableTimeMaximum
    if (allocated(self%distanceTableTime                   )) call Dealloc_Array(self%distanceTableTime                   )
    if (allocated(self%distanceTableComovingDistance       )) call Dealloc_Array(self%distanceTableComovingDistance       )
    if (allocated(self%distanceTableComovingDistanceNegated)) call Dealloc_Array(self%distanceTableComovingDistanceNegated)
    call Alloc_Array(self%distanceTableTime                   ,[self%distanceTableNumberPoints])
    call Alloc_Array(self%distanceTableComovingDistance       ,[self%distanceTableNumberPoints])
    call Alloc_Array(self%distanceTableComovingDistanceNegated,[self%distanceTableNumberPoints])
    read (stateFile) self%distanceTableTime,self%distanceTableComovingDistance,self%distanceTableComovingDistanceNegated
    ! Ensure that interpolation objects will get reset.
    self%resetInterpolation               =.true.
    self%resetInterpolationInverse        =.true.
    self%resetInterpolationDistance       =.true.
    self%resetInterpolationDistanceInverse=.true.
    return
  end subroutine matterLambdaStateRestore
