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

  !!{RST
  An implementation of the cosmological functions class for cosmologies consisting of collisionless matter and dark energy with an equation of state of the form: :math:`P=\rho^w` with :math:`w(a)=w_0+w_1 a (1-a)`.
  !!}

  use :: Root_Finder, only : rootFinder

  integer         , parameter :: ageTableNPointsPerDecade     =300
  ! Factor by which the earliest epoch of the expansion factor tabulation is placed below the epoch of single component
  ! domination - see the equivalent tabulation in the parent class for why that bound is seeded rather than left to follow
  ! the request.
  double precision, parameter :: ageTableTimeEarliestFactor   =100.0d0

  ! Factor by which one component of Universe must dominate others such that we can ignore the others.
  double precision, parameter :: factorDominate               =100.0d0

  ! Variables used in root finding.
  double precision            :: factorDominateCurrent
  !$omp threadprivate(factorDominateCurrent)

  !![
  <cosmologyFunctions name="cosmologyFunctionsMatterDarkEnergy" docformat="rst">
   <description>
   Cosmological relations are computed assuming a universe that contains only collisionless matter and dark energy with an equation of state :math:`w(a)=w_0+w_1a(1-a)` :cite:p:`jassal_wmap_2005`, with :math:`w_0=`\ ``[darkEnergyEquationOfStateW0]``, and :math:`w_1=`\ ``[darkEnergyEquationOfStateW1]``.
   </description>
  </cosmologyFunctions>
  !!]
  type, extends(cosmologyFunctionsMatterLambda) :: cosmologyFunctionsMatterDarkEnergy
     !!{RST
     A cosmological functions class for cosmologies consisting of matter plus dark energy with equation of state :math:`w(a)=w_0+a(1-a)w_1`.
     !!}
     private
     double precision             :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1
     type            (rootFinder) :: finderDomination           , finderEquality
   contains
     !![
     <methods docformat="rst">
       <method description="Set a module-scope pointer to ``self``."                                                        method="targetSelf"                  />
       <method description="Return the derivative of the dark energy exponent with respect to expansion factor."            method="exponentDarkEnergyDerivative"/>
       <method description="Initialize the state of the object, including the root finders and the expansion factor table." method="initialize"                  />
     </methods>
     !!]
     procedure :: cosmicTime                    => matterDarkEnergyCosmicTime
     procedure :: omegaDarkEnergyEpochal        => matterDarkEnergyOmegaDarkEnergyEpochal
     procedure :: hubbleParameterEpochal        => matterDarkEnergyHubbleParameterEpochal
     procedure :: hubbleParameterRateOfChange   => matterDarkEnergyHubbleParameterRateOfChange
     procedure :: equationOfStateDarkEnergy     => matterDarkEnergyEquationOfStateDarkEnergy
     procedure :: exponentDarkEnergy            => matterDarkEnergyExponentDarkEnergy
     procedure :: exponentDarkEnergyDerivative  => matterDarkEnergyExponentDarkEnergyDerivative
     procedure :: equalityEpochMatterDarkEnergy => matterDarkEnergyEqualityEpochMatterDarkEnergy
     procedure :: dominationEpochMatter         => matterDarkEnergyDominationEpochMatter
     procedure :: distanceComoving              => matterDarkEnergyDistanceComoving
     procedure :: timeAtDistanceComoving        => matterDarkEnergyTimeAtDistanceComoving
     procedure :: distanceComovingConvert       => matterDarkEnergyDistanceComovingConvert
     procedure :: expansionFactorTabulate       => matterDarkEnergyMakeExpansionFactorTable
     procedure :: targetSelf                    => matterDarkEnergyTargetSelf
     procedure :: initialize                    => matterDarkEnergyInitialize
  end type cosmologyFunctionsMatterDarkEnergy

  ! Module scope pointer to the current object.
  class(cosmologyFunctionsMatterDarkEnergy), pointer :: self_ => null()
  !$omp threadprivate(self_)

  interface cosmologyFunctionsMatterDarkEnergy
     !!{RST
     Constructors for the matter plus dark energy cosmological functions class.
     !!}
     module procedure matterDarkEnergyConstructorParameters
     module procedure matterDarkEnergyConstructorInternal
  end interface cosmologyFunctionsMatterDarkEnergy

contains

  function matterDarkEnergyConstructorParameters(parameters) result(self)
    !!{RST
    Default constructor for the matter plus dark energy cosmological functions class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (cosmologyFunctionsMatterDarkEnergy)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    double precision                                                    :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1

    !![
    <inputParameter docformat="rst">
      <name>darkEnergyEquationOfStateW0</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>
      The equation of state parameter for dark energy, :math:`w_0`, defined such that :math:`P=\rho^w` with :math:`w(a)=w_0+w_1 a (1-a)`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>darkEnergyEquationOfStateW1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The equation of state parameter for dark energy, :math:`w_1`, defined such that :math:`P=\rho^w` with :math:`w(a)=w_0+w_1 a (1-a)`.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    ! Use it to construct a matter plus dark energy cosmological functions class.
    self=cosmologyFunctionsMatterDarkEnergy(                             &
         &                                  cosmologyParameters_       , &
         &                                  darkEnergyEquationOfStateW0, &
         &                                  darkEnergyEquationOfStateW1  &
         &                                 )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function matterDarkEnergyConstructorParameters

  function matterDarkEnergyConstructorInternal(cosmologyParameters_,darkEnergyEquationOfStateW0,darkEnergyEquationOfStateW1) result(self)
    !!{RST
    Constructor for the matter plus dark energy cosmological functions class.
    !!}
    use :: Cosmology_Parameters, only : cosmologyParametersClass
    implicit none
    type            (cosmologyFunctionsMatterDarkEnergy)                        :: self
    class           (cosmologyParametersClass          ), intent(in   ), target :: cosmologyParameters_
    double precision                                    , intent(in   )         :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1
    !![
    <constructorAssign variables="*cosmologyParameters_, darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1"/>
    !!]

    call self%initialize()
    return
  end function matterDarkEnergyConstructorInternal

  subroutine matterDarkEnergyInitialize(self)
    !!{RST
    Initialize the state of a matter plus dark energy cosmological functions object. This is called from the constructors of
    this class and of any class which extends it. As ``self`` is declared to be polymorphic here, any dark energy equation of
    state functions overridden by an extending class are correctly used when initializing the root finders and the expansion
    factor table.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout) :: self
    double precision                                                    :: rangeExpandDownward      , rangeExpandUpward, &
         &                                                                 darkEnergyExponentCurrent

    self%collapsingUniverse                  =.false.
    self%enableRangeChecks                   =.true.
    self%expansionFactorMaximum              =0.0d0
    self%timeTurnaround                      =0.0d0
    self%timeMaximum                         =0.0d0
    self%expansionRatePrevious               =-1.0d0
    self%expansionRateExpansionFactorPrevious=-1.0d0
    ! Build root finder for epoch of matter domination.
    if (self%cosmologyParameters_%OmegaDarkEnergy() /= 0.0d0)                             &
         & self%finderDomination=rootFinder(                                              &
         &                                  rootFunction     =matterDarkEnergyDomination, &
         &                                  toleranceAbsolute=0.0d+0                    , &
         &                                  toleranceRelative=1.0d-6                      &
         &                                 )
    ! Build root finder for epoch of matter equality.
    darkEnergyExponentCurrent=self%exponentDarkEnergy(expansionFactor=1.0d0)
    if (darkEnergyExponentCurrent > -3.0d0) then
       ! Dark energy density is decaying less rapidly than matter.
       if (self%cosmologyParameters_%OmegaMatter() < self%cosmologyParameters_%OmegaDarkEnergy()) then
          ! Matter density is less than dark energy density. Search backward for epoch of domination.
          rangeExpandUpward  =1.0d0
          rangeExpandDownward=0.5d0
       else
          ! Matter density is greater than dark energy density. Search forward for epoch of domination.
          rangeExpandUpward  =2.0d0
          rangeExpandDownward=1.0d0
       end if
    else
       ! Dark energy density is decaying more rapidly than matter.
       if (self%cosmologyParameters_%OmegaMatter() < self%cosmologyParameters_%OmegaDarkEnergy()) then
          ! Matter density is less than dark energy density. Search forward for epoch of domination.
          rangeExpandUpward  =2.0d0
          rangeExpandDownward=1.0d0
       else
          ! Matter density is greater than dark energy density. Search backward for epoch of domination.
          rangeExpandUpward  =1.0d0
          rangeExpandDownward=0.50d0
       end if
    end if
    self%finderEquality=rootFinder(                                                &
         &                         rootFunction       =matterDarkEnergyDomination, &
         &                         toleranceAbsolute  =0.0d0                     , &
         &                         toleranceRelative  =1.0d-6                    , &
         &                         rangeExpandUpward  =rangeExpandUpward         , &
         &                         rangeExpandDownward=rangeExpandDownward       , &
         &                         rangeExpandType    =rangeExpandMultiplicative   &
         &                        )
    ! Force a build of the expansion factor table, which will determine if this Universe collapses.
    call self%expansionFactorTabulate()
    return
  end subroutine matterDarkEnergyInitialize

  double precision function matterDarkEnergyCosmicTime(self,expansionFactor,collapsingPhase)
    !!{RST
    Return the cosmological matter density in units of the critical density at the present day.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   )           :: expansionFactor
    logical                                             , intent(in   ), optional :: collapsingPhase
    logical                                                                       :: collapsingPhaseActual

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
       ! We assume that the universe does not collapse.
       call Error_Report('non-collapsing universe assumed'//{introspection:location})
    else
       ! In expanding phase ensure that sufficiently small and large expansion factors have been reached.
       ! Request a wider range rather than adjusting the bounds directly: the tabulation derives its range from the request
       ! and from the lattice it already holds, so a bound assigned here would simply be overwritten and the loop never end.
       do while (self%ageTableExpansionFactor(                        1) > expansionFactor)
          call self%expansionFactorTabulate(self%ageTableTimeMinimum/2.0d0)
       end do
       do while (self%ageTableExpansionFactor(self%ageTableNumberPoints) < expansionFactor)
          call self%expansionFactorTabulate(self%ageTableTimeMaximum*2.0d0)
       end do
    end if
    ! Interpolate to get cosmic time.
    matterDarkEnergyCosmicTime=self%interpolatorTime%interpolate(expansionFactor)
    return
  end function matterDarkEnergyCosmicTime

  double precision function matterDarkEnergyOmegaDarkEnergyEpochal(self,time,expansionFactor,collapsingPhase)
    !!{RST
    Return the dark energy density parameter at expansion factor ``expansionFactor``.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsStandard
    use :: Error               , only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    logical                                             , intent(in   ), optional :: collapsingPhase
    double precision                                                              :: expansionFactorActual
    !$GLC attributes unused :: collapsingPhase

    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    matterDarkEnergyOmegaDarkEnergyEpochal                                                                                      &
         & =                        self%cosmologyParameters_%OmegaDarkEnergy       (                                         ) &
         &  *expansionFactorActual**self                     %exponentDarkEnergy    (expansionFactor    =expansionFactorActual) &
         &  *(                                                                                                                  &
         &                          self%cosmologyParameters_%HubbleConstant        (hubbleUnitsStandard                      ) &
         &    /                     self                     %HubbleParameterEpochal(expansionFactor    =expansionFactorActual) &
         &   )**2
    return
  end function matterDarkEnergyOmegaDarkEnergyEpochal

  double precision function matterDarkEnergyDominationEpochMatter(self,dominateFactor)
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Root_Finder                   , only : rangeExpandMultiplicative
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout) :: self
    double precision                                    , intent(in   ) :: dominateFactor
    double precision                                                    :: aDominantCurvature , aDominantDarkEnergy      , &
         &                                                                 aMatterEquality    , darkEnergyExponentCurrent, &
         &                                                                 rangeExpandDownward, rangeExpandUpward

    ! Choose present day as default - will be used if no other densities present (i.e. Einstein-de Sitter).
    matterDarkEnergyDominationEpochMatter=1.0d0
    ! Case where dark energy is present.
    if (self%cosmologyParameters_%OmegaDarkEnergy() /= 0.0d0) then
       darkEnergyExponentCurrent=self%exponentDarkEnergy(expansionFactor=1.0d0)
       if (darkEnergyExponentCurrent > -3.0d0) then
          ! Dark energy density is decaying less rapidly than matter.
          if (self%cosmologyParameters_%OmegaMatter() < dominateFactor*self%cosmologyParameters_%OmegaDarkEnergy()) then
             ! Matter density is less than dominated dark energy density. Search backward for epoch of domination.
             rangeExpandUpward  =1.0d0
             rangeExpandDownward=0.5d0
          else
             ! Matter density is greater than dominated dark energy density. Search forward for epoch of domination.
             rangeExpandUpward  =2.0d0
             rangeExpandDownward=1.0d0
          end if
       else
          ! Dark energy density is decaying more rapidly than matter.
          if (self%cosmologyParameters_%OmegaMatter() < dominateFactor*self%cosmologyParameters_%OmegaDarkEnergy()) then
             ! Matter density is less than dominated dark energy density. Search forward for epoch of domination.
             rangeExpandUpward  =2.0d0
             rangeExpandDownward=1.0d0
          else
             ! Matter density is greater than dominated dark energy density. Search backward for epoch of domination.
             rangeExpandUpward  =1.0d0
             rangeExpandDownward=0.50d0
          end if
       end if
       call self%finderDomination%rangeExpand(                                               &
            &                                 rangeExpandUpward  =rangeExpandUpward        , &
            &                                 rangeExpandDownward=rangeExpandDownward      , &
            &                                 rangeExpandType    =rangeExpandMultiplicative  &
            &                                )
       factorDominateCurrent = dominateFactor       
       call self%targetSelf()
       aDominantDarkEnergy  =self%finderDomination%find(rootGuess=1.0d0)
       ! Choose earliest expansion factor.
       matterDarkEnergyDominationEpochMatter=min(matterDarkEnergyDominationEpochMatter,aDominantDarkEnergy)
    end if
    if (self%cosmologyParameters_%OmegaCurvature() /= 0.0d0) then
       ! Find the expansion factor of matter-curvature equality.
       aMatterEquality=self%equalityEpochMatterCurvature(requestTypeExpansionFactor)
       ! Find the earlier expansion factor at which matter dominates by the specified amount (ratio of matter
       ! to curvature density scales as the expansion factor).
       aDominantCurvature=aMatterEquality/dominateFactor
       ! Choose earliest expansion factor.
       matterDarkEnergyDominationEpochMatter=min(matterDarkEnergyDominationEpochMatter,aDominantCurvature)
    end if
    return
  end function matterDarkEnergyDominationEpochMatter

  double precision function matterDarkEnergyHubbleParameterEpochal(self,time,expansionFactor,collapsingPhase)
    !!{RST
    Returns the Hubble parameter at the request cosmological time, ``time``, or expansion factor, ``expansionFactor``.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsStandard
    use :: Error               , only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    logical                                             , intent(in   ), optional :: collapsingPhase
    double precision                                                              :: expansionFactorActual, sqrtArgument
    
    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    ! Compute the Hubble parameter at the specified expansion factor.
    sqrtArgument=                                                                                            &
         &       max(                                                                                        &
         &            self%cosmologyParameters_%OmegaMatter    ()                                            &
         &           /expansionFactorActual**3                                                               &
         &           +self%cosmologyParameters_%OmegaDarkEnergy()                                            &
         &           *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)  &
         &           +self%cosmologyParameters_%OmegaCurvature ()                                            &
         &           /expansionFactorActual**2                                                             , &
         &           0.0d0                                                                                   &
         &          )
    matterDarkEnergyHubbleParameterEpochal=self%cosmologyParameters_%HubbleConstant(hubbleUnitsStandard)*sqrt(sqrtArgument)
    ! Make the Hubble parameter negative if we are in the collapsing phase of the Universe.
    if (self%collapsingUniverse) then
       if    (present(time           )) then
          if    (time>self%timeTurnaround) matterDarkEnergyHubbleParameterEpochal=-matterDarkEnergyHubbleParameterEpochal
       else
          if (present(collapsingPhase)) then
             if (collapsingPhase         ) matterDarkEnergyHubbleParameterEpochal=-matterDarkEnergyHubbleParameterEpochal
          end if
       end if
    end if
    return
  end function matterDarkEnergyHubbleParameterEpochal

  double precision function matterDarkEnergyHubbleParameterRateOfChange(self,time,expansionFactor,collapsingPhase)
    !!{RST
    Returns the rate of change of the Hubble parameter at the requested cosmological time, ``time``, or expansion factor, ``expansionFactor``.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    logical                                             , intent(in   ), optional :: collapsingPhase
    double precision                                                              :: expansionFactorActual

    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    ! Compute the rat of change of the Hubble parameter.
    matterDarkEnergyHubbleParameterRateOfChange                                                                &
         & =0.5d0                                                                                              &
         & *self%hubbleParameterEpochal(expansionFactor=expansionFactorActual,collapsingPhase=collapsingPhase) &
         & *self%expansionRate         (                expansionFactorActual                                ) &
         & /(                                                                                                  &
         &   +self%cosmologyParameters_%OmegaMatter    ()                                                      &
         &   /expansionFactorActual**3                                                                         &
         &   +self%cosmologyParameters_%OmegaDarkEnergy()                                                      &
         &   *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)            &
         &   +self%cosmologyParameters_%OmegaCurvature ()                                                      &
         &   /expansionFactorActual**2                                                                         &
         &  )                                                                                                  &
         & *(                                                                                                  &
         &   -3.0d0*self%cosmologyParameters_%OmegaMatter()                                                    &
         &   /expansionFactorActual**3                                                                         &
         &   +self%cosmologyParameters_%OmegaDarkEnergy  ()                                                    &
         &   *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)            &
         &   *(                                                                                                &
         &     +self%exponentDarkEnergy(expansionFactor=expansionFactorActual)                                 &
         &     +expansionFactorActual                                                                          &
         &     *log(expansionFactorActual)                                                                     &
         &     *self%exponentDarkEnergyDerivative(expansionFactor=expansionFactorActual)                       &
         &    )                                                                                                &
         &   -2.0d0*self%cosmologyParameters_%OmegaCurvature()                                                 &
         &   /expansionFactorActual**2                                                                         &
         & )
    return
  end function matterDarkEnergyHubbleParameterRateOfChange

  double precision function matterDarkEnergyEqualityEpochMatterDarkEnergy(self,requestType)
    !!{RST
    Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    !!}
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor, requestTypeTime
    implicit none
    class  (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    integer                                    , intent(in   ), optional :: requestType
    integer                                                              :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if
    factorDominateCurrent =  1.0d0
    call self%targetSelf()
    matterDarkEnergyEqualityEpochMatterDarkEnergy=self%finderEquality%find(rootGuess=1.0d0)
    if (present(requestType)) then
       if (requestType == requestTypeTime)                                    &
            &                  matterDarkEnergyEqualityEpochMatterDarkEnergy  &
            & =self%cosmicTime(matterDarkEnergyEqualityEpochMatterDarkEnergy)
    end if
    return
  end function matterDarkEnergyEqualityEpochMatterDarkEnergy

  double precision function matterDarkEnergyDomination(expansionFactor)
    !!{RST
    Function used in root finding when seeking epoch at which matter dominates over dark energy.
    !!}
    implicit none
    double precision, intent(in   ) :: expansionFactor

    matterDarkEnergyDomination=+self_%cosmologyParameters_%OmegaMatter    ()                               &
         &                     /expansionFactor**3                                                         &
         &                     -factorDominateCurrent                                                      &
         &                     *self_%cosmologyParameters_%OmegaDarkEnergy()                               &
         &                     *expansionFactor**self_%exponentDarkEnergy(expansionFactor=expansionFactor)
    return
  end function matterDarkEnergyDomination

  subroutine matterDarkEnergyTargetSelf(self)
    !!{RST
    Set a module-scope pointer to the current dark energy cosmology functions object.
    !!}
    implicit none
    class(cosmologyFunctionsMatterDarkEnergy), intent(in   ), target :: self

    self_ => self
    return
  end subroutine matterDarkEnergyTargetSelf

  subroutine matterDarkEnergyMakeExpansionFactorTable(self,time)
    !!{RST
    Builds a table of expansion factor vs. time for dark energy universes.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsTime
    use :: Numerical_Ranges    , only : Range_Pinned   , Range_Lattice_Extend, rangeLattice, gridSchemePerDecade
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy)             , intent(inout), target   :: self
    double precision                                                 , intent(in   ), optional :: time
    double precision                                    , parameter                            :: turnaroundTimeTolerance         =1.0d-12
    double precision                                    , parameter                            :: odeToleranceAbsolute            =1.0d-9 , odeToleranceRelative    =1.0d-9
    ! The bounds are pinned to single lattice steps rather than to whole decades: at three hundred points per decade, with an
    ! ordinary differential equation integrated across each of them, rounding a bound out to a whole decade would add three
    ! hundred epochs for the sake of reaching one.
    integer                                             , parameter                            :: ageTableAnchorEvery             =1
    double precision                                    , allocatable, dimension(:)            :: ageTableExpansionFactorTemporary        , ageTableTimeTemporary
    logical                                             , allocatable, dimension(:)            :: isComputed
    type            (rangeLattice                      )                                       :: latticeTime
    integer                                                                                    :: iTime
    logical                                                                                    :: carryOver
    double precision                                                                           :: OmegaDominant                           , deltaTime                      , &
         &                                                                                        densityPower                            , timeDominant                   , &
         &                                                                                        expansionFactorDominant                 , timeActual
    logical                                                                                    :: solutionFound                           , timeExceeded
    double precision                                                 , dimension(3)            :: timesTarget

    ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
    call self%densityScalingEarlyTime(factorDominate,densityPower,expansionFactorDominant,OmegaDominant)
    ! Find the corresponding time.
    timeDominant=-2.0d0/densityPower/self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime)/sqrt(OmegaDominant)/expansionFactorDominant**(0.5d0*densityPower)
    ! Assume initially that this Universe does not collapse; the tabulation below detects otherwise. Once found, the epoch of
    ! turnaround is a property of the cosmology rather than of the tabulation, so it is determined on the first build and
    ! retained thereafter. Redetermining it on every call, as was formerly attempted here, cannot work once previously computed
    ! values are carried over, since the detection below examines only those points which are computed afresh.
    if (.not.self%ageTableInitialized) then
       self%collapsingUniverse    =.false.
       self%expansionFactorMaximum=0.0d0
       self%timeTurnaround        =0.0d0
       self%timeMaximum           =0.0d0
    end if
    ! Find the range of epochs to tabulate, pinning it to an absolute lattice so that the epochs tabulated - and therefore every
    ! expansion factor interpolated between them - depend only on which lattice points are spanned, and not on the sequence of
    ! epochs which happened to be requested. The bounds formerly grew by repeated multiplication by a fixed factor from an origin
    ! set by the first request, and the abscissae were then redistributed over the result by `Make_Range`, so that the previously
    ! computed points pasted back into the new table carried the *old* spacing. The table thereby became a patchwork of spacings
    ! which no single `ageTableInverseDeltaLogTime` could describe, and the index computed from it in `expansionFactor` drifted
    ! further from the truth with every extension.
    if (present(time)) then
       timeActual=time
    else
       timeActual=timeDominant
    end if
    timesTarget=[timeActual,timeDominant,timeDominant/ageTableTimeEarliestFactor]
    if (self%collapsingUniverse) then
       latticeTime=Range_Pinned(                                          &
            &                                   timesTarget             , &
            &                                   ageTableNPointsPerDecade, &
            &                                   gridSchemePerDecade     , &
            &                   marginFactor  = 2.0d0                   , &
            &                   anchorEvery   = ageTableAnchorEvery     , &
            &                   limitMaximum  = self%timeTurnaround     , &
            &                   latticeCurrent= self%latticeAgeTable      &
            &                  )
    else
       latticeTime=Range_Pinned(                                          &
            &                                   timesTarget             , &
            &                                   ageTableNPointsPerDecade, &
            &                                   gridSchemePerDecade     , &
            &                   marginFactor  = 2.0d0                   , &
            &                   anchorEvery   = ageTableAnchorEvery     , &
            &                   latticeCurrent= self%latticeAgeTable      &
            &                  )
    end if
    ! The expansion factor is integrated forward from an initial condition imposed at the earliest tabulated epoch, so every
    ! value depends on where that integration began. The tabulation can therefore be carried over only while the earliest epoch
    ! is unchanged; where it moves, everything is discarded by presenting an undefined lattice to the extension.
    carryOver=self%ageTableInitialized .and. self%latticeAgeTable%isDefined()
    if (carryOver) carryOver=latticeTime%indexMinimum == self%latticeAgeTable%indexMinimum
    if (.not.carryOver) then
       if (allocated(self%ageTableExpansionFactor)) deallocate(self%ageTableExpansionFactor)
       self%latticeAgeTable=rangeLattice()
    end if
    call Range_Lattice_Extend(self%latticeAgeTable,latticeTime,self%ageTableExpansionFactor,isComputed)
    self%latticeAgeTable     =latticeTime
    self%ageTableNumberPoints=latticeTime%count
    if (allocated(self%ageTableTime)) deallocate(self%ageTableTime)
    self%ageTableTime        =latticeTime%values()
    self%ageTableTimeMinimum =self%ageTableTime(                        1)
    self%ageTableTimeMaximum =self%ageTableTime(self%ageTableNumberPoints)
    ! Compute quantities required for table interpolation. These are evaluated again after the tabulation below, so that they
    ! describe the table actually held even where it is truncated at turnaround.
    self%ageTableTimeLogarithmicMinimum=log(self%ageTableTimeMinimum)
    self%ageTableInverseDeltaLogTime   =dble(self%ageTableNumberPoints-1)/log(self%ageTableTimeMaximum/self%ageTableTimeMinimum)
    ! For the initial time, we approximate that we are at sufficiently early times that a single component dominates the
    ! Universe and use the appropriate analytic solution.
    if (.not.isComputed(1))                                             &
         &    self%ageTableExpansionFactor           (               1) &
         & =(                                                           &
         &   -0.5d0                                                     &
         &   *densityPower                                              &
         &   *self%ageTableTime                       (              1) &
         &   *self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime) &
         &   *sqrt(OmegaDominant)                                       &
         &  )**(-2.0d0/densityPower)
    ! Solve ODE to get corresponding expansion factors.
    self%iTableTurnaround  =  self%ageTableNumberPoints
    call self%targetSelf()
    do iTime=2,self%ageTableNumberPoints
       ! Compute the expansion factor if it is not already computed.
       if (.not.isComputed(iTime)) then
          self%ageTableExpansionFactor(iTime)=matterDarkEnergyExpansionFactorChange(                                       &
               &                                                                    self%ageTableTime           (iTime-1), &
               &                                                                    self%ageTableTime           (iTime  ), &
               &                                                                    self%ageTableExpansionFactor(iTime-1)  &
               &                                                                   )
          ! Check for a universe which is no longer expanding (i.e. has reached its maximum expansion).
          if (self%ageTableExpansionFactor(iTime) == self%ageTableExpansionFactor(iTime-1)) then
             ! Record that we have a collapsing Universe.
             self%collapsingUniverse=.true.
             ! Record the maximum expansion factor.
             self%expansionFactorMaximum=self%ageTableExpansionFactor(iTime-1)
             ! Find the time of maximum expansion by bisection. Disable checks of epoch ranges while we search for the maximum
             ! expansion factor as we may exceed the maximum while searching.
             self%timeTurnaround   =(self%ageTableTime(iTime-1)+self%ageTableTime(iTime-2))/2.0d0
             deltaTime             =(self%ageTableTime(iTime-1)-self%ageTableTime(iTime-2))/2.0d0
             solutionFound         =.false.
             self%enableRangeChecks=.false.
             do while (.not.solutionFound)
                timeExceeded=(                                                                &
                     &         matterDarkEnergyExpansionFactorChange(                         &
                     &                                 self%ageTableTime           (iTime-2), &
                     &                                 self%timeTurnaround                  , &
                     &                                 self%ageTableExpansionFactor(iTime-2)  &
                     &                                )                                       &
                     &        >=                                                              &
                     &         self%ageTableExpansionFactor(iTime-1)                          &
                     &       )
                solutionFound=timeExceeded .and. deltaTime < turnaroundTimeTolerance*self%timeTurnaround
                if (.not.solutionFound) then
                   deltaTime=0.5d0*deltaTime
                   if (timeExceeded) then
                      self%timeTurnaround=self%timeTurnaround-deltaTime
                   else
                      self%timeTurnaround=self%timeTurnaround+deltaTime
                   end if
                end if
             end do
             self%enableRangeChecks=.true.
             self%timeMaximum      =2.0d0*self%timeTurnaround
             ! Limit the tables to the expanding part of the evolution. The lattice held is shortened to match: its points are
             ! unchanged, so the shortened table remains a set of points of the same absolute lattice and can be extended - or
             ! carried over - exactly as before. Now that the epoch of turnaround is known it is imposed as a hard limit on the
             ! range above, so this truncation happens only on the build which discovers it.
             self%iTableTurnaround     =iTime-2
             self%ageTableNumberPoints =self%iTableTurnaround
             self%latticeAgeTable%count=self%ageTableNumberPoints
             call Move_Alloc(self%ageTableTime           ,ageTableTimeTemporary           )
             call Move_Alloc(self%ageTableExpansionFactor,ageTableExpansionFactorTemporary)
             allocate(self%ageTableTime           (self%ageTableNumberPoints))
             allocate(self%ageTableExpansionFactor(self%ageTableNumberPoints))
             self%ageTableTime           =ageTableTimeTemporary           (1:self%ageTableNumberPoints)
             self%ageTableExpansionFactor=ageTableExpansionFactorTemporary(1:self%ageTableNumberPoints)
             self%ageTableTimeMaximum    =self%ageTableTime               (  self%ageTableNumberPoints)
             exit
          end if
       end if
    end do
    ! Recompute the quantities required for table interpolation, so that they describe the table actually held even where it was
    ! truncated at turnaround above.
    self%ageTableInverseDeltaLogTime=dble(self%ageTableNumberPoints-1)/log(self%ageTableTimeMaximum/self%ageTableTimeMinimum)
    if (allocated(self%interpolatorTime)) deallocate(self%interpolatorTime)
    allocate(self%interpolatorTime)
    self%interpolatorTime=interpolator(self%ageTableExpansionFactor,self%ageTableTime)
    ! Flag that the table is now initialized.
    self%ageTableInitialized=.true.
    return
  end subroutine matterDarkEnergyMakeExpansionFactorTable

  double precision function matterDarkEnergyExpansionFactorChange(timeStart,timeEnd,expansionFactorStart)
    !!{RST
    Compute the expansion factor at time ``timeEnd`` given an initial value ``expansionFactorStart`` at time ``timeStart``.
    !!}
    use :: Numerical_ODE_Solvers, only : odeSolver
    implicit none
    double precision           , intent(in   ) :: expansionFactorStart       , timeEnd                     , &
         &                                        timeStart
    double precision           , dimension(1)  :: y
    double precision           , parameter     :: odeToleranceAbsolute=1.0d-9, odeToleranceRelative=1.0d-12
    double precision                           :: time
    type            (odeSolver)                :: solver

    time     =timeStart
    y     (1)=expansionFactorStart
    solver   =odeSolver(1_c_size_t,matterDarkEnergyAgeTableODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
    call solver%solve(time,timeEnd,y)
    matterDarkEnergyExpansionFactorChange=y(1)
    return
  end function matterDarkEnergyExpansionFactorChange

  integer function matterDarkEnergyAgeTableODEs(t,a,dadt)
    !!{RST
    System of differential equations to solve for expansion factor vs. age.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision              , intent(in   ) :: t
    double precision, dimension(:), intent(in   ) :: a
    double precision, dimension(:), intent(  out) :: dadt
    !$GLC attributes unused :: t

    if (a(1) <= 0.0d0) then
       dadt(1)=0.0d0
    else
       dadt(1)=a(1)*self_%expansionRate(a(1))
    end if
    matterDarkEnergyAgeTableODEs=GSL_Success
  end function matterDarkEnergyAgeTableODEs

  double precision function matterDarkEnergyTimeAtDistanceComoving(self,comovingDistance)
    !!{RST
    Returns the cosmological time corresponding to given ``comovingDistance``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout) :: self
    double precision                                    , intent(in   ) :: comovingDistance
    !$GLC attributes unused :: self, comovingDistance

    matterDarkEnergyTimeAtDistanceComoving=0.0d0
    call Error_Report('functionality not implemented'//{introspection:location})
    return
  end function matterDarkEnergyTimeAtDistanceComoving

  double precision function matterDarkEnergyDistanceComoving(self,time)
    !!{RST
    Returns the comoving distance to cosmological time ``time``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout) :: self
    double precision                                    , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    matterDarkEnergyDistanceComoving=0.0d0
    call Error_Report('functionality not implemented'//{introspection:location})
    return
   end function matterDarkEnergyDistanceComoving

  double precision function matterDarkEnergyDistanceComovingConvert(self,output,distanceLuminosity,distanceModulus,distanceModulusKCorrected,redshift)
    !!{RST
    Convert between different measures of distance.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    integer                                             , intent(in   )           :: output
    double precision                                    , intent(in   ), optional :: distanceModulus, distanceModulusKCorrected, &
         &                                                                           redshift       , distanceLuminosity
    !$GLC attributes unused :: self, output, distanceModulus, distanceModulusKCorrected, redshift, distanceLuminosity

    matterDarkEnergyDistanceComovingConvert=0.0d0
    call Error_Report('functionality not implemented'//{introspection:location})
    return
  end function matterDarkEnergyDistanceComovingConvert

  double precision function matterDarkEnergyEquationOfStateDarkEnergy(self,time,expansionFactor)
    !!{RST
    Return the dark energy equation of state.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    double precision                                                              :: expansionFactorActual

    matterDarkEnergyEquationOfStateDarkEnergy=self%darkEnergyEquationOfStateW0
    if (self%darkEnergyEquationOfStateW1 /= 0.0d0) then
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else if (present(time)) then
          expansionFactorActual=self%expansionFactor(time)
       else
          call Error_Report('equation of state is time dependent, but no time given'//{introspection:location})
          expansionFactorActual=1.0d0
       end if
       matterDarkEnergyEquationOfStateDarkEnergy=+matterDarkEnergyEquationOfStateDarkEnergy &
            &                                    +self%darkEnergyEquationOfStateW1          &
            &                                    *       expansionFactorActual              &
            &                                    *(1.0d0-expansionFactorActual)
    end if
   return
  end function matterDarkEnergyEquationOfStateDarkEnergy

  double precision function matterDarkEnergyExponentDarkEnergy(self,time,expansionFactor)
    !!{RST
    Return the dark energy exponent.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    double precision                                                              :: expansionFactorActual

    matterDarkEnergyExponentDarkEnergy=-3.0d0*(1.0d0+self%darkEnergyEquationOfStateW0)
    if (self%darkEnergyEquationOfStateW1 /= 0.0d0) then
       if      (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else if (present(time           )) then
          expansionFactorActual=self%expansionFactor(time)
       else
          call Error_Report('equation of state is time dependent, but no time given'//{introspection:location})
          expansionFactorActual=1.0d0
       end if
       if (expansionFactorActual /= 1.0d0) &
            & matterDarkEnergyExponentDarkEnergy=matterDarkEnergyExponentDarkEnergy+3.0d0*self%darkEnergyEquationOfStateW1*(1.0d0-expansionFactorActual)**2/2.0d0/log(expansionFactorActual)
    end if
    return
  end function matterDarkEnergyExponentDarkEnergy

  double precision function matterDarkEnergyExponentDarkEnergyDerivative(self,time,expansionFactor)
    !!{RST
    Return the derivative of the dark energy exponent with respect to expansion factor.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    double precision                                                              :: expansionFactorActual

    if      (present(expansionFactor)) then
       expansionFactorActual=expansionFactor
    else if (present(time           )) then
       expansionFactorActual=self%expansionFactor(time)
    else
       if (self%darkEnergyEquationOfStateW1 /= 0.0d0) call Error_Report('equation of state is time dependent, but no time given'//{introspection:location})
       expansionFactorActual=1.0d0
    end if
    if (expansionFactorActual == 1.0d0) then
       matterDarkEnergyExponentDarkEnergyDerivative= &
            & -1.5d0&
            & *self%darkEnergyEquationOfStateW1
    else
       matterDarkEnergyExponentDarkEnergyDerivative= &
            & -1.5d0                                 &
            & *self%darkEnergyEquationOfStateW1      &
            & *(                                     &
            &       (1.0d0-expansionFactorActual)**2 &
            &   /          expansionFactorActual     &
            &   /log(      expansionFactorActual)**2 &
            &   +2.0d0                               &
            &   *   (1.0d0-expansionFactorActual)    &
            &   /log(      expansionFactorActual)    &
            &  )
    end if
    return
  end function matterDarkEnergyExponentDarkEnergyDerivative
