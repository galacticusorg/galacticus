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
  An implementation of dark matter halo profile scale radii in which radii are computed from the concentration.
  !!}

  use :: Cosmology_Functions               , only : cosmologyFunctions            , cosmologyFunctionsClass
  use :: Cosmology_Parameters              , only : cosmologyParameters           , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScale           , darkMatterHaloScaleClass           , darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_Concentration, only : darkMatterProfileConcentration, darkMatterProfileConcentrationClass
  use :: Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMO          , darkMatterProfileDMOClass
  use :: Galacticus_Nodes                  , only : nodeComponentBasic            , nodeComponentDarkMatterProfile     , treeNode
  use :: Virial_Density_Contrast           , only : virialDensityContrast         , virialDensityContrastClass

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusConcentration">
    <description>Dark matter halo scale radii are computed from the concentration.</description>
    <deepCopy>
      <functionClass variables="darkMatterHaloScaleDefinition"/>
    </deepCopy>
    <stateStorable>
      <functionClass variables="darkMatterHaloScaleDefinition"/>
    </stateStorable>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusConcentration
     !!{
     A dark matter halo profile scale radius class which computes radii from the concentration.
     !!}
     private
     class           (cosmologyParametersClass                          ), pointer :: cosmologyParameters_              => null()
     class           (cosmologyFunctionsClass                           ), pointer :: cosmologyFunctions_               => null()
     class           (darkMatterHaloScaleClass                          ), pointer :: darkMatterHaloScale_              => null()
     class           (darkMatterProfileDMOClass                         ), pointer :: darkMatterProfileDMO_             => null(), darkMatterProfileDMODefinition  => null()
     class           (virialDensityContrastClass                        ), pointer :: virialDensityContrast_            => null(), virialDensityContrastDefinition => null()
     class           (darkMatterProfileConcentrationClass               ), pointer :: darkMatterProfileConcentration_   => null()
     type            (darkMatterHaloScaleVirialDensityContrastDefinition), pointer :: darkMatterHaloScaleDefinition     => null()
     logical                                                                       :: correctForConcentrationDefinition          , useMeanConcentration
     double precision                                                              :: massRatioPrevious
   contains
     final     ::           concentrationDestructor
     procedure :: radius => concentrationRadius
  end type darkMatterProfileScaleRadiusConcentration

  interface darkMatterProfileScaleRadiusConcentration
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusConcentration} dark matter halo profile scale radius class.
     !!}
     module procedure concentrationConstructorParameters
     module procedure concentrationConstructorInternal
  end interface darkMatterProfileScaleRadiusConcentration

  ! Container type used to maintain a stack of state.
  type :: concentrationState
     class           (darkMatterProfileScaleRadiusConcentration), pointer :: self                  => null()
     type            (treeNode                                 ), pointer :: nodeWork              => null(), node => null()
     class           (nodeComponentBasic                       ), pointer :: basic                 => null()
     class           (nodeComponentDarkMatterProfile           ), pointer :: darkMatterProfile     => null()
     double precision                                                     :: concentrationOriginal          , mass
  end type concentrationState

  ! State stack.
  type   (concentrationState), allocatable, dimension(:) :: state_
  integer                                                :: stateCount=0
  !$omp threadprivate(stateCount,state_)

contains

  function concentrationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusConcentration} dark matter halo profile scale radius class which takes a
    parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (darkMatterProfileScaleRadiusConcentration)                :: self
    type   (inputParameters                          ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class  (cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    class  (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    class  (darkMatterProfileDMOClass                ), pointer       :: darkMatterProfileDMO_
    class  (virialDensityContrastClass               ), pointer       :: virialDensityContrast_
    class  (darkMatterProfileConcentrationClass      ), pointer       :: darkMatterProfileConcentration_
    logical                                                           :: correctForConcentrationDefinition, useMeanConcentration

    !![
    <inputParameter>
      <name>correctForConcentrationDefinition</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, then when computing dark matter profile scale radii using concentrations, any difference between the current definition of halo scales
        (i.e. typically virial density contrast definitions) and density profiles and those assumed in measuring the concentrations will be taken into account.
        If false, the concentration is applied blindly.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>useMeanConcentration</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, then when computing dark matter profile scale radii using concentrations do not account for any possible scatter in the concentration-mass relation.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"            name="cosmologyParameters_"            source="parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"           name="darkMatterProfileDMO_"           source="parameters"/>
    <objectBuilder class="virialDensityContrast"          name="virialDensityContrast_"          source="parameters"/>
    <objectBuilder class="darkMatterProfileConcentration" name="darkMatterProfileConcentration_" source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusConcentration(correctForConcentrationDefinition,useMeanConcentration,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,virialDensityContrast_,darkMatterProfileConcentration_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"           />
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="darkMatterHaloScale_"           />
    <objectDestructor name="darkMatterProfileDMO_"          />
    <objectDestructor name="virialDensityContrast_"         />
    <objectDestructor name="darkMatterProfileConcentration_"/>
    !!]
    return
  end function concentrationConstructorParameters

  function concentrationConstructorInternal(correctForConcentrationDefinition,useMeanConcentration,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,virialDensityContrast_,darkMatterProfileConcentration_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileScaleRadiusConcentration} dark matter halo profile scale radius class.
    !!}
    implicit none
    type   (darkMatterProfileScaleRadiusConcentration)                        :: self
    class  (cosmologyParametersClass                 ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    class  (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    class  (darkMatterProfileDMOClass                ), intent(in   ), target :: darkMatterProfileDMO_
    class  (virialDensityContrastClass               ), intent(in   ), target :: virialDensityContrast_
    class  (darkMatterProfileConcentrationClass      ), intent(in   ), target :: darkMatterProfileConcentration_
    logical                                           , intent(in   )         :: correctForConcentrationDefinition, useMeanConcentration
    !![
    <constructorAssign variables="correctForConcentrationDefinition, useMeanConcentration, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *virialDensityContrast_, *darkMatterProfileConcentration_"/>
    !!]

    ! Get definitions of virial density contrast and dark matter profile as used by the concentration definition.
    allocate(self%darkMatterHaloScaleDefinition)
    !![
    <referenceAcquire   isResult="yes" owner="self" target="virialDensityContrastDefinition" source="self%darkMatterProfileConcentration_%     densityContrastDefinition()"/>
    <referenceAcquire   isResult="yes" owner="self" target="darkMatterProfileDMODefinition"  source="self%darkMatterProfileConcentration_%darkMatterProfileDMODefinition()"/>
    <referenceConstruct isResult="yes" owner="self" object="darkMatterHaloScaleDefinition"   constructor="darkMatterHaloScaleVirialDensityContrastDefinition(self%cosmologyParameters_,self%cosmologyFunctions_,self%virialDensityContrastDefinition)"/>
    !!]
    self%massRatioPrevious=2.0d0
    return
  end function concentrationConstructorInternal

  subroutine concentrationDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileScaleRadiusConcentration} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusConcentration), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScaleDefinition"  />
    <objectDestructor name="self%cosmologyParameters_"           />
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%darkMatterHaloScale_"           />
    <objectDestructor name="self%darkMatterProfileDMO_"          />
    <objectDestructor name="self%virialDensityContrast_"         />
    <objectDestructor name="self%darkMatterProfileConcentration_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition" />
    <objectDestructor name="self%virialDensityContrastDefinition"/>
    !!]
    return
  end subroutine concentrationDestructor

  double precision function concentrationRadius(self,node)
    !!{
    Compute the scale radius of the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Calculations_Resets , only : Calculations_Reset
    use :: Numerical_Comparison, only : Values_Differ
    use :: Root_Finder         , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (darkMatterProfileScaleRadiusConcentration), intent(inout), target        :: self
    type            (treeNode                                 ), intent(inout), target        :: node
    class           (nodeComponentBasic                       ), pointer                      :: basic
    double precision                                           , parameter                    :: massRatioBuffer      =1.1d0, massRatioShrink=0.99d0
    type            (concentrationState                       ), allocatable   , dimension(:) :: concentrationStateTmp
    type            (rootFinder                               )                               :: finder
    double precision                                                                          :: concentration              , massDefinition        , &
         &                                                                                       massRatio
    integer                                                                                   :: i

    ! Increment the state stack.
    if (.not.allocated(state_)) then
       allocate(state_(1))
    else if (stateCount == size(state_)) then
       call move_alloc(state_,concentrationStateTmp)
       allocate(state_(size(concentrationStateTmp)+1))
       state_(1:size(concentrationStateTmp))=concentrationStateTmp
       do i=1,size(concentrationStateTmp)
          nullify(concentrationStateTmp(i)%self             )
          nullify(concentrationStateTmp(i)%nodeWork         )
          nullify(concentrationStateTmp(i)%node             )
          nullify(concentrationStateTmp(i)%basic            )
          nullify(concentrationStateTmp(i)%darkMatterProfile)
       end do
    end if
    stateCount=stateCount+1
    ! Find the original concentration.
    if (self%useMeanConcentration) then
       state_(stateCount)%concentrationOriginal=self%darkMatterProfileConcentration_%concentrationMean(node)
    else
       state_(stateCount)%concentrationOriginal=self%darkMatterProfileConcentration_%concentration    (node)
    end if
    ! Determine if concentration must be corrected.
    if (self%correctForConcentrationDefinition) then
       ! Get the basic component of the supplied node and extract its mass.
       basic                   => node %basic()
       state_(stateCount)%mass =  basic%mass ()
       ! If there is no difference between the alt and non-alt virial density contrasts, then no correction need be made.
       if     (                                                                                                      &
            &  Values_Differ(                                                                                        &
            &                       self%virialDensityContrast_         %densityContrast(basic%mass(),basic%time()), &
            &                       self%virialDensityContrastDefinition%densityContrast(basic%mass(),basic%time()), &
            &                relTol=1.0d-6                                                                           &
            &               )                                                                                        &
            & ) then
          ! Create a node and set the mass and time.
          state_(stateCount)%self                       => self
          state_(stateCount)%node                       => node
          state_(stateCount)%nodeWork                   => treeNode                                     (                 )
          state_(stateCount)%nodeWork         %hostTree => node%hostTree
          state_(stateCount)%basic                      => state_(stateCount)%nodeWork%basic            (autoCreate=.true.)
          state_(stateCount)%darkMatterProfile          => state_(stateCount)%nodeWork%darkMatterProfile(autoCreate=.true.)
          call state_(stateCount)%basic%timeSet            (basic%time())
          call state_(stateCount)%basic%timeLastIsolatedSet(basic%time())
          ! The finder is initialized each time as it is allocated on the stack - this allows this function to be called recursively.
          finder=rootFinder(                                                             &
               &            rootFunction                 =concentrationMassRoot        , &
               &            toleranceRelative            =1.0d-3                       , &
               &            rangeExpandUpward            =1.0d0*self%massRatioPrevious , &
               &            rangeExpandDownward          =1.0d0/self%massRatioPrevious , &
               &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &            rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &            rangeExpandType              =rangeExpandMultiplicative      &
               &           )
          massDefinition=finder%find(rootGuess=state_(stateCount)%mass)
          ! Find the ratio of the recovered mass under the given definition to the input mass, defined to be always greater than
          ! unity. This will be used as the basis of the range expansion for the next solution.
          if (massDefinition > state_(stateCount)%mass) then
             massRatio=1.0d0*(massDefinition/state_(stateCount)%mass)
          else
             massRatio=1.0d0/(massDefinition/state_(stateCount)%mass)
          end if
          ! Increase this mass ratio by a small factor.
          massRatio=massRatioBuffer*massRatio
          ! If this new mass ratio exceeds our previous mass ratio, update the previous mass ratio for use in the next
          ! solution. Otherwise, shrink the previous mass ratio by a small amount.
          if (massRatio > self%massRatioPrevious) then
             self%massRatioPrevious=massRatio
          else
             self%massRatioPrevious=massRatioShrink*self%massRatioPrevious
          end if
          ! Update the work node properties and computed concentration.
          call state_(stateCount)%basic%massSet(massDefinition)
          call Calculations_Reset(state_(stateCount)%nodeWork)
          ! Find the concentration.
          if (self%useMeanConcentration) then
             ! We are simply using the mean concentration-mass relation here.
             concentration=+self%darkMatterProfileConcentration_%concentrationMean(state_(stateCount)%nodeWork)
          else
             ! In this case we need to allow for possible scatter in the concentration mass relation. Therefore, we take the original
             ! concentration (which may include some scatter away from the mean relation) and scale it by the ratio of the mean
             ! concentrations for the corrected and original nodes.
             concentration=+                                                       state_(stateCount)%concentrationOriginal  &
                  &        *self%darkMatterProfileConcentration_%concentrationMean(state_(stateCount)%nodeWork             ) &
                  &        /self%darkMatterProfileConcentration_%concentrationMean(                   node                 )
          end if
          concentrationRadius=+self%darkMatterHaloScaleDefinition%radiusVirial (state_(stateCount)%nodeWork) &
               &              /                                   concentration
          call state_(stateCount)%nodeWork%destroy()
          deallocate(state_(stateCount)%nodeWork)
       else
          concentrationRadius=+self%darkMatterHaloScale_%radiusVirial(node) &
               &              /state_(stateCount)%concentrationOriginal
       end if
    else
       concentrationRadius=+self%darkMatterHaloScale_%radiusVirial(node) &
            &              /state_(stateCount)%concentrationOriginal
    end if
    ! Release stack.
    nullify(state_(stateCount)%self             )
    nullify(state_(stateCount)%nodeWork         )
    nullify(state_(stateCount)%node             )
    nullify(state_(stateCount)%basic            )
    nullify(state_(stateCount)%darkMatterProfile)
    stateCount=stateCount-1
    return
  end function concentrationRadius

  double precision function concentrationMassRoot(massDefinitionTrial)
    !!{
    Root function used to find the mass of a halo corresponding to the definition used for a particular concentration class.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Mass_Distributions , only : massDistributionClass
    implicit none
    double precision                       , intent(in   ) :: massDefinitionTrial
    class           (massDistributionClass), pointer       :: massDistribution_
    double precision                                       :: radiusOuterDefinition, concentrationDefinition, &
         &                                                    radiusCore           , massOuter              , &
         &                                                    radiusOuter          , densityOuter

    ! Set the mass of the worker node.
    call state_(stateCount)%basic%massSet(massDefinitionTrial)
    call Calculations_Reset(state_(stateCount)%nodeWork)
    ! Get outer radius for this trial definition mass.
    radiusOuterDefinition=state_(stateCount)%self%darkMatterHaloScaleDefinition%radiusVirial(state_(stateCount)%nodeWork)
    ! Get concentration for this a trial definition mass.
    if (state_(stateCount)%self%useMeanConcentration) then
       ! We are simply using the mean concentration-mass relation here.
       concentrationDefinition=state_(stateCount)%self%darkMatterProfileConcentration_%concentrationMean(state_(stateCount)%nodeWork)
    else
       ! In this case we need to allow for possible scatter in the concentration mass relation. Therefore, we take the original
       ! concentration (which may include some scatter away from the mean relation) and scale it by the ratio of the mean
       ! concentrations for the corrected and original nodes.
       concentrationDefinition=+state_(stateCount)                                     %concentrationOriginal                              &
            &                  *state_(stateCount)%self%darkMatterProfileConcentration_%concentrationMean    (state_(stateCount)%nodeWork) &
            &                  /state_(stateCount)%self%darkMatterProfileConcentration_%concentrationMean    (state_(stateCount)%node    )
    end if
    ! Get core radius.
    radiusCore=radiusOuterDefinition/concentrationDefinition
    call state_(stateCount)%darkMatterProfile%scaleSet(radiusCore)
    call Calculations_Reset(state_(stateCount)%nodeWork)
    ! Find the non-alt density.
    densityOuter=+state_(stateCount)%self%cosmologyFunctions_   %matterDensityEpochal(                                state_(stateCount)%basic%time()) &
         &       *state_(stateCount)%self%virialDensityContrast_%densityContrast     (state_(stateCount)%basic%mass(),state_(stateCount)%basic%time())
    ! Get the current mass distribution.
    massDistribution_ => state_(stateCount)%self%darkMatterProfileDMODefinition%get(state_(stateCount)%nodeWork)
    ! Solve for radius which encloses required non-alt density.
    radiusOuter=massDistribution_%radiusEnclosingDensity(densityOuter)
    ! Get the mass within this radius.
    massOuter  =massDistribution_%massEnclosedBySphere  ( radiusOuter)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Return root function.
    concentrationMassRoot=massOuter-state_(stateCount)%mass
    return
  end function concentrationMassRoot

