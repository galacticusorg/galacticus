!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  !% matter and dark energy with an equation of state of the form: $P=\rho^w$ with $w(a)=w_0+w_1 a (1-a)$.

  !# <cosmologyFunctions name="cosmologyFunctionsMatterDarkEnergy">
  !#  <description>Cosmological relations are computed assuming a universe that contains only matter and dark energy with an equation of state $w(a)=w_0+w_1a(1-a)$.</description>
  !# </cosmologyFunctions>
  use FGSL
  use Cosmology_Parameters

  integer         , parameter :: matterDarkEnergyAgeTableNPointsPerDecade     =300
  double precision, parameter :: matterDarkEnergyAgeTableNPointsPerOctave     =dble(matterDarkEnergyAgeTableNPointsPerDecade)*log(2.0d0)/log(10.0d0)
  double precision, parameter :: matterDarkEnergyAgeTableIncrementFactor      =exp(int(matterDarkEnergyAgeTableNPointsPerOctave+1.0d0)*log(10.0d0)/dble(matterDarkEnergyAgeTableNPointsPerDecade))
  integer         , parameter :: matterDarkEnergyDistanceTableNPointsPerDecade=100

  ! Factor by which one component of Universe must dominate others such that we can ignore the others.
  double precision, parameter :: matterDarkEnergyDominateFactor               =100.0d0

  ! Variables used in root finding.
  double precision            :: matterDarkEnergyDominateFactorCurrent
  !$omp threadprivate(matterDarkEnergyDominateFactorCurrent)

  ! Default parameters of the dark energy equation of state
  logical                     :: matterDarkEnergyDefaultsRead=.false.
  double precision            :: matterDarkEnergyEquationOfStateW0Default, matterDarkEnergyEquationOfStateW1Default

  type, extends(cosmologyFunctionsMatterLambda) :: cosmologyFunctionsMatterDarkEnergy
     !% A cosmological functions class for cosmologies consisting of matter plus dark energy with equation of state $w(a)=w_0+a(1-a)w_1$.
     private
     double precision :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1
   contains
     !@ <objectMethods>
     !@   <object>cosmologyFunctionsMatterDarkEnergy</object>
     !@   <objectMethod>
     !@     <method>targetSelf</method>
     !@     <type>void</type>
     !@     <arguments></arguments>
     !@     <description>Set a module-scope pointer to {\tt self}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>exponentDarkEnergyDerivative</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ [time]\argin, \doublezero\ [expansionFactor]\argin</arguments>
     !@     <description>Return the derivative of the dark energy exponent with respect to expansion factor.</description>
     !@   </objectMethod>
     !@ </objectMethods>
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
  end type cosmologyFunctionsMatterDarkEnergy

  ! Module scope pointer to the current object.
  class(cosmologyFunctionsMatterDarkEnergy), pointer :: matterDarkEnergySelfGlobal
  !$omp threadprivate(matterDarkEnergySelfGlobal)

  interface cosmologyFunctionsMatterDarkEnergy
     !% Constructors for the matter plus dark energy cosmological functions class.
     module procedure matterDarkEnergyDefaultConstructor
     module procedure matterDarkEnergyConstructor
  end interface cosmologyFunctionsMatterDarkEnergy

contains

  function matterDarkEnergyDefaultConstructor()
    !% Default constructor for the matter plus dark energy cosmological functions class.
    use Input_Parameters
    use Cosmology_Parameters
    implicit none
    type (cosmologyFunctionsMatterDarkEnergy), target  :: matterDarkEnergyDefaultConstructor
    class(cosmologyParametersClass          ), pointer :: thisCosmologyParameters

    if (.not.matterDarkEnergyDefaultsRead) then
       !$omp critical(matterDarkEnergyDefaultConstructorParameters)
       if (.not.matterDarkEnergyDefaultsRead) then
          ! Read the dark energy equation of state.
          !@ <inputParameter>
          !@   <name>darkEnergyEquationOfStateW0</name>
          !@   <defaultValue>-1 (cosmological constant)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The equation of state parameter for dark energy, $w_0$, defined such that $P=\rho^w$ with $w(a)=w_0+w_1 a (1-a)$.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>cosmology</group>
          !@ </inputParameter>
          call Get_Input_Parameter('darkEnergyEquationOfStateW0',matterDarkEnergyEquationOfStateW0Default,defaultValue=-1.0d0)
          !@ <inputParameter>
          !@   <name>darkEnergyEquationOfStateW1</name>
          !@   <defaultValue>0 (constant equation of state)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The equation of state parameter for dark energy, $w_1$, defined such that $P=\rho^w$ with $w(a)=w_0+w_1 a (1-a)$.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>cosmology</group>
          !@ </inputParameter>
          call Get_Input_Parameter('darkEnergyEquationOfStateW1',matterDarkEnergyEquationOfStateW1Default,defaultValue=0.0d0)
          matterDarkEnergyDefaultsRead=.true.
       end if
       !$omp end critical(matterDarkEnergyDefaultConstructorParameters)
    end if
    ! Get the default cosmological parameters.
    thisCosmologyParameters => cosmologyParameters()
    ! Use it to construct a matter plus dark energy cosmological functions class.
    matterDarkEnergyDefaultConstructor                                            &
         & =matterDarkEnergyConstructor(                                          &
         &                              thisCosmologyParameters                 , &
         &                              matterDarkEnergyEquationOfStateW0Default, &
         &                              matterDarkEnergyEquationOfStateW1Default  &
         &                             )
    return
  end function matterDarkEnergyDefaultConstructor

  function matterDarkEnergyConstructor(thisCosmologyParameters,darkEnergyEquationOfStateW0,darkEnergyEquationOfStateW1)
    !% Constructor for the matter plus dark energy cosmological functions class.
    use Numerical_Comparison
    use Cosmology_Parameters
    use ISO_Varying_String
    use ODE_Solver
    use, intrinsic :: ISO_C_Binding
    implicit none
    type            (cosmologyFunctionsMatterDarkEnergy)                        :: matterDarkEnergyConstructor
    class           (cosmologyParametersClass          ), intent(in   ), target :: thisCosmologyParameters
    double precision                                    , intent(in   )         :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1

    ! Store a pointer to the cosmological parameters object.
    matterDarkEnergyConstructor%cosmology => thisCosmologyParameters
    ! Store equation of state.
    matterDarkEnergyConstructor%darkEnergyEquationOfStateW0=darkEnergyEquationOfStateW0
    matterDarkEnergyConstructor%darkEnergyEquationOfStateW1=darkEnergyEquationOfStateW1
    ! Force a build of the expansion factor table, which will determine if this Universe collapses.
    call matterDarkEnergyConstructor%expansionFactorTabulate()
    return
  end function matterDarkEnergyConstructor

  double precision function matterDarkEnergyCosmicTime(self,expansionFactor,collapsingPhase)
    !% Return the cosmological matter density in units of the critical density at the present day.
    use Galacticus_Error
    use Numerical_Interpolation
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   )           :: expansionFactor
    logical                                             , intent(in   ), optional :: collapsingPhase
    logical                                                                       :: collapsingPhaseActual

    ! Validate the input.
    if (.not.self%expansionFactorIsValid(expansionFactor)) &
         & call Galacticus_Error_Report('matterDarkEnergyCosmicTime','expansion factor is invalid')
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
       ! We assume that the universe does not collapse.
       call Galacticus_Error_Report('matterDarkEnergyCosmicTime','non-collapsing universe assumed')
    else
       ! In expanding phase ensure that sufficiently small and large expansion factors have been reached.
       do while (self%ageTableExpansionFactor(                        1) > expansionFactor)
          self%ageTableTimeMinimum=    self%ageTableTimeMinimum/matterDarkEnergyAgeTableIncrementFactor
          call self%expansionFactorTabulate()
       end do
       do while (self%ageTableExpansionFactor(self%ageTableNumberPoints) < expansionFactor)
          self%ageTableTimeMaximum=max(self%ageTableTimeMaximum*matterDarkEnergyAgeTableIncrementFactor,self%timeTurnaround)
          call self%expansionFactorTabulate()
       end do
    end if
    ! Interpolate to get cosmic time.
    matterDarkEnergyCosmicTime                             &
         & =Interpolate(                               &
         &              self%ageTableNumberPoints    , &
         &              self%ageTableExpansionFactor , &
         &              self%ageTableTime            , &
         &              self%interpolationObject     , &
         &              self%interpolationAccelerator, &
         &              expansionFactor              , &
         &              reset=self%resetInterpolation  &
         &             )
    return
  end function matterDarkEnergyCosmicTime

  double precision function matterDarkEnergyOmegaDarkEnergyEpochal(self,time,expansionFactor,collapsingPhase)
    !% Return the dark energy density parameter at expansion factor {\tt expansionFactor}.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    logical                                             , intent(in   ), optional :: collapsingPhase
    double precision                                                              :: expansionFactorActual

   ! Determine the actual expansion factor to use.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('matterDarkEnergyOmegaDarkEnergyEpochal','only one of time or expansion factor can be specified')
       else
          expansionFactorActual=self%expansionFactor(time)
       end if
    else
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else
          call Galacticus_Error_Report('matterDarkEnergyOmegaDarkEnergyEpochal','either a time or expansion factor must be specified')
       end if
    end if
    matterDarkEnergyOmegaDarkEnergyEpochal                                                                       &
         & =                        self%cosmology%OmegaDarkEnergy       (                                     ) &
         &  *expansionFactorActual**self          %exponentDarkEnergy    (expansionFactor=expansionFactorActual) &
         &  *(                                                                                                   &
         &                          self%cosmology%HubbleConstant        (unitsStandard                        ) &
         &    /                     self%          HubbleParameterEpochal(expansionFactor=expansionFactorActual) &
         &   )**2
    return
  end function matterDarkEnergyOmegaDarkEnergyEpochal

  double precision function matterDarkEnergyDominationEpochMatter(self,dominateFactor)
    use Cosmology_Functions_Parameters
    use Root_Finder
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout) :: self
    type            (rootFinder                        ), save          :: finder
    !$omp threadprivate(finder)
    double precision                                    , intent(in   ) :: dominateFactor
    double precision                                                    :: aDominantCurvature              , aDominantDarkEnergy              , &
         &                                                                 aMatterEquality                 , darkEnergyExponentCurrent        , &
         &                                                                 expansionFactorDominantCurvature, expansionFactorDominantDarkEnergy, &
         &                                                                 rangeExpandDownward             , rangeExpandUpward

    ! Choose present day as default - will be used if no other densities present (i.e. Einsetin-de Sitter).
    matterDarkEnergyDominationEpochMatter=1.0d0
    ! Case where dark energy is present.
    if (self%cosmology%OmegaDarkEnergy() /= 0.0d0) then
       if (.not.finder%isInitialized()) then
          darkEnergyExponentCurrent=self%exponentDarkEnergy(expansionFactor=1.0d0)
          if (darkEnergyExponentCurrent > -3.0d0) then
             ! Dark energy density is decaying less rapidly than matter.
             if (self%cosmology%OmegaMatter() < dominateFactor*self%cosmology%OmegaDarkEnergy()) then
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
             if (self%cosmology%OmegaMatter() < dominateFactor*self%cosmology%OmegaDarkEnergy()) then
                ! Matter density is less than dominated dark energy density. Search forward for epoch of domination.
                rangeExpandUpward  =2.0d0
                rangeExpandDownward=1.0d0
             else
                ! Matter density is greater than dominated dark energy density. Search backward for epoch of domination.
                rangeExpandUpward  =1.0d0
                rangeExpandDownward=0.50d0
             end if
          end if
          call finder%rootFunction(                                               &
               &                   matterDarkEnergyDomination                     &
               &                  )
          call finder%tolerance   (                                               &
               &                   toleranceAbsolute  =0.0d0                    , &
               &                   toleranceRelative  =1.0d-6                     &
               &                  )
          call finder%rangeExpand (                                               &
               &                   rangeExpandUpward  =rangeExpandUpward        , &
               &                   rangeExpandDownward=rangeExpandDownward      , &
               &                   rangeExpandType    =rangeExpandMultiplicative  &
               &                  )
       end if
       matterDarkEnergyDominateFactorCurrent = dominateFactor
       call self%targetSelf()
       aDominantDarkEnergy  =finder%find(rootGuess=1.0d0)
       ! Choose earliest expansion factor.
       matterDarkEnergyDominationEpochMatter=min(matterDarkEnergyDominationEpochMatter,aDominantDarkEnergy)
    end if
    if (self%cosmology%OmegaCurvature() /= 0.0d0) then
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
    !% Returns the Hubble parameter at the request cosmological time, {\tt time}, or expansion factor, {\tt expansionFactor}.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    logical                                             , intent(in   ), optional :: collapsingPhase
    double precision                                                              :: expansionFactorActual, sqrtArgument

    ! Determine the actual expansion factor to use.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('matterDarkEnergyHubbleParameterEpochal','only one of time or expansion factor can be specified')
       else
          expansionFactorActual=self%expansionFactor(time)
       end if
    else
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else
          call Galacticus_Error_Report('matterDarkEnergyHubbleParameterEpochal','either a time or expansion factor must be specified')
       end if
    end if
    ! Compute the Hubble parameter at the specified expansion factor.
    sqrtArgument=                                                                                            &
         &       max(                                                                                        &
         &            self%cosmology%OmegaMatter    ()                                                       &
         &           /expansionFactorActual**3                                                               &
         &           +self%cosmology%OmegaDarkEnergy()                                                       &
         &           *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)  &
         &           +self%cosmology%OmegaCurvature ()                                                       &
         &           /expansionFactorActual**2                                                             , &
         &           0.0d0                                                                                   &
         &          )
    matterDarkEnergyHubbleParameterEpochal=self%cosmology%HubbleConstant(unitsStandard)*sqrt(sqrtArgument)
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
    !% Returns the rate of change of the Hubble parameter at the requested cosmological time, {\tt time}, or expansion factor, {\tt expansionFactor}.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    logical                                             , intent(in   ), optional :: collapsingPhase
    double precision                                                              :: expansionFactorActual, sqrtArgument

    ! Determine the actual expansion factor to use.
    if (present(time)) then
       if (present(expansionFactor)) then
          call Galacticus_Error_Report('matterDarkEnergyHubbleParameterRateOfChange','only one of time or expansion factor can be specified')
       else
          expansionFactorActual=self%expansionFactor(time)
       end if
    else
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else
          call Galacticus_Error_Report('matterDarkEnergyHubbleParameterRateOfChange','either a time or expansion factor must be specified')
       end if
    end if
    ! Compute the rat of change of the Hubble parameter.
    matterDarkEnergyHubbleParameterRateOfChange                                                                &
         & =0.5d0                                                                                              &
         & *self%hubbleParameterEpochal(expansionFactor=expansionFactorActual,collapsingPhase=collapsingPhase) &
         & *self%expansionRate         (                expansionFactorActual                                ) &
         & /(                                                                                                  &
         &   +self%cosmology%OmegaMatter    ()                                                                 &
         &   /expansionFactorActual**3                                                                         &
         &   +self%cosmology%OmegaDarkEnergy()                                                                 &
         &   *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)            &
         &   +self%cosmology%OmegaCurvature ()                                                                 &
         &   /expansionFactorActual**2                                                                         &
         &  )                                                                                                  &
         & *(                                                                                                  &
         &   -3.0d0*self%cosmology%OmegaMatter()                                                               &
         &   /expansionFactorActual**3                                                                         &
         &   +self%cosmology%OmegaDarkEnergy  ()                                                               &
         &   *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)            &
         &   *(                                                                                                &
         &     +self%exponentDarkEnergy(expansionFactor=expansionFactorActual)                                 &
         &     +expansionFactorActual                                                                          &
         &     *log(expansionFactorActual)                                                                     &
         &     *self%exponentDarkEnergyDerivative(expansionFactor=expansionFactorActual)                       &
         &    )                                                                                                &
         &   -2.0d0*self%cosmology%OmegaCurvature()                                                            &
         &   /expansionFactorActual**2                                                                         &
         & )
    return
  end function matterDarkEnergyHubbleParameterRateOfChange

  double precision function matterDarkEnergyEqualityEpochMatterDarkEnergy(self,requestType)
    !% Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    use Cosmology_Functions_Parameters
    use Cosmology_Parameters
    use Root_Finder
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    integer                                             , intent(in   ), optional :: requestType
    type            (rootFinder                        ), save                    :: finder
    !$omp threadprivate(finder)
    integer                                                                       :: requestTypeActual
    double precision                                                              :: darkEnergyExponentCurrent, rangeExpandDownward, &
         &                                                                           rangeExpandUpward

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if
    if (.not.finder%isInitialized()) then
       darkEnergyExponentCurrent=self%exponentDarkEnergy(expansionFactor=1.0d0)
       if (darkEnergyExponentCurrent > -3.0d0) then
          ! Dark energy density is decaying less rapidly than matter.
          if (self%cosmology%OmegaMatter() < self%cosmology%OmegaDarkEnergy()) then
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
          if (self%cosmology%OmegaMatter() < self%cosmology%OmegaDarkEnergy()) then
             ! Matter density is less than dark energy density. Search forward for epoch of domination.
             rangeExpandUpward  =2.0d0
             rangeExpandDownward=1.0d0
          else
             ! Matter density is greater than dark energy density. Search backward for epoch of domination.
             rangeExpandUpward  =1.0d0
             rangeExpandDownward=0.50d0
          end if
       end if
       call finder%rootFunction(                                               &
            &                   matterDarkEnergyDomination                     &
            &                  )
       call finder%tolerance   (                                               &
            &                   toleranceAbsolute  =0.0d0                    , &
            &                   toleranceRelative  =1.0d-6                     &
            &                  )
       call finder%rangeExpand (                                               &
            &                   rangeExpandUpward  =rangeExpandUpward        , &
            &                   rangeExpandDownward=rangeExpandDownward      , &
            &                   rangeExpandType    =rangeExpandMultiplicative  &
            &                  )
    end if
    matterDarkEnergyDominateFactorCurrent =  1.0d0
    call self%targetSelf()
    matterDarkEnergyEqualityEpochMatterDarkEnergy=finder%find(rootGuess=1.0d0)
    if (present(requestType)) then
       if (requestType == requestTypeTime)                                    &
            &                  matterDarkEnergyEqualityEpochMatterDarkEnergy  &
            & =self%cosmicTime(matterDarkEnergyEqualityEpochMatterDarkEnergy)
    end if
    return
  end function matterDarkEnergyEqualityEpochMatterDarkEnergy

  double precision function matterDarkEnergyDomination(expansionFactor)
    !% Function used in root finding when seeking epoch at which matter dominates over dark energy.
    implicit none
    double precision, intent(in   ) :: expansionFactor

    matterDarkEnergyDomination=                                                                                                 &
         &                      matterDarkEnergySelfGlobal%cosmology%OmegaMatter    ()                                          &
         &                     /expansionFactor**3                                                                              &
         &                     -matterDarkEnergyDominateFactorCurrent                                                           &
         &                     *matterDarkEnergySelfGlobal%cosmology%OmegaDarkEnergy()                                          &
         &                     *expansionFactor**matterDarkEnergySelfGlobal%exponentDarkEnergy(expansionFactor=expansionFactor)
    return
  end function matterDarkEnergyDomination

  subroutine matterDarkEnergyTargetSelf(self)
    !% Set a module-scope pointer to the current dark energy cosmology functions object.
    implicit none
    class(cosmologyFunctionsMatterDarkEnergy), intent(in   ), target :: self

    matterDarkEnergySelfGlobal => self
    return
  end subroutine matterDarkEnergyTargetSelf

  subroutine matterDarkEnergyMakeExpansionFactorTable(self,time)
    !% Builds a table of expansion factor vs. time for dark energy universes.
    use Numerical_Interpolation
    use Numerical_Ranges
    use Memory_Management
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Parameters
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy)             , intent(inout), target   :: self
    double precision                                                 , intent(in   ), optional :: time
    double precision                                    , parameter                            :: turnaroundTimeTolerance         =1.0d-12
    double precision                                    , parameter                            :: odeToleranceAbsolute            =1.0d-9 , odeToleranceRelative    =1.0d-9
    double precision                                    , allocatable, dimension(:)            :: ageTableExpansionFactorTemporary        , ageTableTimeTemporary
    integer                                                                                    :: iTime                                   , prefixPointCount
    double precision                                                                           :: OmegaDominant                           , deltaTime                      , &
         &                                                                                        densityPower                            , expansionFactor      (1)       , &
         &                                                                                        expansionFactorDominant                 , timeActual                     , &
         &                                                                                        timeDominant
    logical                                                                                    :: solutionFound                           , timeExceeded

    ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
    call self%densityScalingEarlyTime(matterDarkEnergyDominateFactor,densityPower,expansionFactorDominant,OmegaDominant)
    ! Find the corresponding time.
    timeDominant=-2.0d0/densityPower/self%cosmology%HubbleConstant(unitsTime)/sqrt(OmegaDominant)/expansionFactorDominant**(0.5d0*densityPower)
    ! Find minimum and maximum times to tabulate.
    if (present(time)) then
       timeActual=time
       do while (self%ageTableTimeMinimum > min(timeActual,timeDominant)/2.0d0)
          self%ageTableTimeMinimum=self%ageTableTimeMinimum/matterDarkEnergyAgeTableIncrementFactor
       end do
       do while (self%ageTableTimeMaximum < max(timeActual,timeDominant)*2.0d0)
          self%ageTableTimeMaximum=self%ageTableTimeMaximum*matterDarkEnergyAgeTableIncrementFactor
       end do
    else
       do while (self%ageTableTimeMinimum > timeDominant/2.0d0)
          self%ageTableTimeMinimum=self%ageTableTimeMinimum/matterDarkEnergyAgeTableIncrementFactor
       end do
       do while (self%ageTableTimeMaximum < timeDominant*2.0d0)
          self%ageTableTimeMaximum=self%ageTableTimeMaximum*matterDarkEnergyAgeTableIncrementFactor
       end do
    end if
    ! Determine number of points to tabulate.
    self%ageTableNumberPoints=int(log10(self%ageTableTimeMaximum/self%ageTableTimeMinimum)     *dble(matterDarkEnergyAgeTableNPointsPerDecade))+1
    self%ageTableTimeMaximum =self%ageTableTimeMinimum*10.0d0**(dble(self%ageTableNumberPoints)/dble(matterDarkEnergyAgeTableNPointsPerDecade))
    ! Assume this Universe does not collapse initially.
    self%collapsingUniverse    =.false.
    self%expansionFactorMaximum=0.0d0
    self%timeTurnaround        =0.0d0
    self%timeMaximum           =0.0d0
    ! Deallocate arrays if currently allocated.
    if (allocated(self%ageTableTime)) then
       ! Determine number of points that are being added at the start of the array.
       prefixPointCount=int(log10(self%ageTableTime(1)/self%ageTableTimeMinimum)*dble(matterDarkEnergyAgeTableNPointsPerDecade)+0.5d0)
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
    call self%targetSelf()
    do iTime=2,self%ageTableNumberPoints
       ! Compute the expansion factor if it is not already computed.
       if (self%ageTableExpansionFactor(iTime) < 0.0d0) then
          self%ageTableExpansionFactor(iTime)=matterDarkEnergyExpansionFactorChange(                                  &
               &                                                 self%ageTableTime           (iTime-1), &
               &                                                 self%ageTableTime           (iTime  ), &
               &                                                 self%ageTableExpansionFactor(iTime-1)  &
               &                                                )
          ! Check for a universe which is no longer expanding (i.e. has reached its maximum expansion).
          if (self%ageTableExpansionFactor(iTime) == self%ageTableExpansionFactor(iTime-1)) then
             ! Record that we have a collapsing Universe.
             self%collapsingUniverse=.true.
             ! Record the maximum expansion factor.
             self%expansionFactorMaximum=self%ageTableExpansionFactor(iTime-1)
             ! Find the time of maximum expansion by bisection.
             self%timeTurnaround=(self%ageTableTime(iTime-1)+self%ageTableTime(iTime-2))/2.0d0
             deltaTime          =(self%ageTableTime(iTime-1)-self%ageTableTime(iTime-2))/2.0d0
             solutionFound      =.false.
             do while (.not.solutionFound)
                timeExceeded=(                                                                &
                &              matterDarkEnergyExpansionFactorChange(                         &
                &                                      self%ageTableTime           (iTime-2), &
                &                                      self%timeTurnaround                  , &
                &                                      self%ageTableExpansionFactor(iTime-2)  &
                &                                     )                                       &
                &             >=                                                              &
                &              self%ageTableExpansionFactor(iTime-1)                          &
                &            )
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
             self%timeMaximum=2.0d0*self%timeTurnaround
             ! Limit the tables to the expanding part of the evolution.
             self%iTableTurnaround=iTime-2
             call Move_Alloc(self%ageTableTime           ,ageTableTimeTemporary           )
             call Move_Alloc(self%ageTableExpansionFactor,ageTableExpansionFactorTemporary)
             self%ageTableNumberPoints=self%iTableTurnaround
             call Alloc_Array(self%ageTableTime,           [self%ageTableNumberPoints])
             call Alloc_Array(self%ageTableExpansionFactor,[self%ageTableNumberPoints])
             self%ageTableTime           =ageTableTimeTemporary           (1:self%ageTableNumberPoints)
             self%ageTableExpansionFactor=ageTableExpansionFactorTemporary(1:self%ageTableNumberPoints)
             exit
          end if
       end if
    end do
    call Interpolate_Done(self%interpolationObject       ,self%interpolationAccelerator       ,self%resetInterpolation       )
    call Interpolate_Done(self%interpolationObjectInverse,self%interpolationAcceleratorInverse,self%resetInterpolationInverse)
    self%resetInterpolation       =.true.
    self%resetInterpolationInverse=.true.
    ! Flag that the table is now initialized.
    self%ageTableInitialized      =.true.
    return
  end subroutine matterDarkEnergyMakeExpansionFactorTable

  double precision function matterDarkEnergyExpansionFactorChange(timeStart,timeEnd,expansionFactorStart)
    !% Compute the expansion factor at time {\tt timeEnd} given an initial value {\tt expansionFactorStart} at time {\tt
    !% timeStart}.
    use ODE_Solver
    use, intrinsic :: ISO_C_Binding
    implicit none
    double precision                    , intent(in   ) :: expansionFactorStart       , timeEnd                     , &
         &                                                 timeStart
    double precision                    , dimension(1)  :: y
    double precision                    , parameter     :: odeToleranceAbsolute=1.0d-9, odeToleranceRelative=1.0d-12
    double precision                                    :: time
    type            (c_ptr             )                :: parameterPointer
    type            (fgsl_odeiv_step   )                :: odeStepper
    type            (fgsl_odeiv_control)                :: odeController
    type            (fgsl_odeiv_evolve )                :: odeEvolver
    type            (fgsl_odeiv_system )                :: odeSystem
    logical                                             :: odeReset            =.true.

    time=timeStart
    y(1)=expansionFactorStart
    call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,time,timeEnd,1,y,matterDarkEnergyAgeTableODEs,parameterPointer&
         &,odeToleranceAbsolute,odeToleranceRelative,reset=odeReset)
    call ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
    matterDarkEnergyExpansionFactorChange=y(1)
    return
  end function matterDarkEnergyExpansionFactorChange

  function matterDarkEnergyAgeTableODEs(t,a,dadt,parameterPointer) bind(c)
    !% System of differential equations to solve for expansion factor vs. age.
    use, intrinsic :: ISO_C_Binding
    integer(kind=c_int   )                              :: matterDarkEnergyAgeTableODEs
    real   (kind=c_double)              , value         :: t
    real   (kind=c_double), dimension(1), intent(in   ) :: a
    real   (kind=c_double), dimension(1)                :: dadt
    type   (c_ptr        )              , value         :: parameterPointer

    if (a(1) <= 0.0d0) then
       dadt(1)=0.0d0
    else
       dadt(1)=a(1)*matterDarkEnergySelfGlobal%expansionRate(a(1))
    end if
    matterDarkEnergyAgeTableODEs=FGSL_Success
  end function matterDarkEnergyAgeTableODEs

  double precision function matterDarkEnergyTimeAtDistanceComoving(self,comovingDistance)
    !% Returns the cosmological time corresponding to given {\tt comovingDistance}.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout) :: self
    double precision                                    , intent(in   ) :: comovingDistance

    call Galacticus_Error_Report('matterDarkEnergyTimeAtDistanceComoving','functionality not implemented')
    return
  end function matterDarkEnergyTimeAtDistanceComoving

  double precision function matterDarkEnergyDistanceComoving(self,time)
    !% Returns the comoving distance to cosmological time {\tt time}.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout) :: self
    double precision                                    , intent(in   ) :: time

    call Galacticus_Error_Report('matterDarkEnergyDistanceComoving','functionality not implemented')
    return
   end function matterDarkEnergyDistanceComoving

  double precision function matterDarkEnergyDistanceComovingConvert(self,output,distanceModulus,redshift)
    !% Convert bewteen different measures of distance.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    integer                                             , intent(in   )           :: output
    double precision                                    , intent(in   ), optional :: distanceModulus, redshift

    call Galacticus_Error_Report('matterDarkEnergyDistanceComovingConvert','functionality not implemented')
    return
  end function matterDarkEnergyDistanceComovingConvert

  double precision function matterDarkEnergyEquationOfStateDarkEnergy(self,time,expansionFactor)
    !% Return the dark energy equation of state.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    double precision                                                              :: expansionFactorActual

    if (present(expansionFactor)) then
       expansionFactorActual=expansionFactor
    else if (present(time)) then
       expansionFactorActual=self%expansionFactor(time)
    else
       if (self%darkEnergyEquationOfStateW1 /= 0.0d0) call Galacticus_Error_Report('matterDarkEnergyEquationOfStateDarkEnergy','equation of state is time dependent, but no time given')
       expansionFactorActual=1.0d0
    end if
    matterDarkEnergyEquationOfStateDarkEnergy= &
         &  self%darkEnergyEquationOfStateW0   &
         & +self%darkEnergyEquationOfStateW1   &
         & *       expansionFactorActual       &
         & *(1.0d0-expansionFactorActual)
   return
  end function matterDarkEnergyEquationOfStateDarkEnergy

  double precision function matterDarkEnergyExponentDarkEnergy(self,time,expansionFactor)
    !% Return the dark energy exponent.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    double precision                                                              :: expansionFactorActual

    if      (present(expansionFactor)) then
       expansionFactorActual=expansionFactor
    else if (present(time           )) then
       expansionFactorActual=self%expansionFactor(time)
    else
       if (self%darkEnergyEquationOfStateW1 /= 0.0d0) call Galacticus_Error_Report('matterDarkEnergyExponentDarkEnergy','equation of state is time dependent, but no time given')
       expansionFactorActual=1.0d0
    end if
    if (expansionFactorActual == 1.0d0) then
       matterDarkEnergyExponentDarkEnergy=-3.0d0*(1.0d0+self%darkEnergyEquationOfStateW0)
    else
       matterDarkEnergyExponentDarkEnergy=-3.0d0*(1.0d0+self%darkEnergyEquationOfStateW0)+3.0d0*self%darkEnergyEquationOfStateW1*(1.0d0-expansionFactorActual)**2/2.0d0/log(expansionFactorActual)
    end if
    return
  end function matterDarkEnergyExponentDarkEnergy

  double precision function matterDarkEnergyExponentDarkEnergyDerivative(self,time,expansionFactor)
    !% Return the derivative of the dark energy exponent with respect to expansion factor.
    use Galacticus_Error
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergy), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    double precision                                                              :: expansionFactorActual

    if      (present(expansionFactor)) then
       expansionFactorActual=expansionFactor
    else if (present(time           )) then
       expansionFactorActual=self%expansionFactor(time)
    else
       if (self%darkEnergyEquationOfStateW1 /= 0.0d0) call Galacticus_Error_Report('matterDarkEnergyExponentDarkEnergyDerivative','equation of state is time dependent, but no time given')
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
