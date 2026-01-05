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
  Implementation of a cooling radius class for isothermal halos, assuming collisional ionization equilibrium such that cooling
  time scales as inverse density.
  !!}

  use :: Cooling_Times          , only : coolingTimeClass
  use :: Cooling_Times_Available, only : coolingTimeAvailableClass
  use :: Cosmology_Functions    , only : cosmologyFunctions       , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Kind_Numbers           , only : kind_int8

  !![
  <coolingRadius name="coolingRadiusIsothermal">
   <description>
  !!]

  !![
    A cooling radius class that computes the cooling radius by assuming an isothermal density profile, and a cooling rate
    proportional to density squared. This implies a cooling time:
    \begin{equation}
     t_\mathrm{cool} \equiv {E\over\dot{E}} \propto \rho(r)^{-1}.
    \end{equation}
    The cooling radius is then derived using
    \begin{equation}
     \rho(r_\mathrm{cool}) \propto t_\mathrm{available}^{-1}
    \end{equation}
    which implies
    \begin{equation}
     r_\mathrm{cool} = r_\mathrm{virial} \left( { t_\mathrm{available} \over t_\mathrm{cool,virial}} \right)^{1/2},
    \end{equation}
    where $t_\mathrm{cool,virial}$ is the cooling time at the virial radius.
   </description>
   <deepCopy>
    <functionClass variables="radiation"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="radiation"/>
   </stateStorable>
  </coolingRadius>
  !!]
  type, extends(coolingRadiusClass) :: coolingRadiusIsothermal
     !!{
     Implementation of cooling radius class in which the cooling radius is defined as that radius at which the time available
     for cooling equals the cooling time.
     !!}
     private
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_        => null()
     class           (darkMatterHaloScaleClass               ), pointer :: darkMatterHaloScale_       => null()
     class           (coolingTimeAvailableClass              ), pointer :: coolingTimeAvailable_      => null()
     class           (coolingTimeClass                       ), pointer :: coolingTime_               => null()
     type            (radiationFieldCosmicMicrowaveBackground), pointer :: radiation                  => null()
     integer         (kind=kind_int8                         )          :: lastUniqueID               =  -1
     integer                                                            :: abundancesCount                     , chemicalsCount
     ! Stored values of cooling radius.
     logical                                                            :: radiusComputed                      , radiusGrowthRateComputed
     double precision                                                   :: radiusGrowthRateStored              , radiusStored
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                     isothermalDestructor
     procedure :: autoHook         => isothermalAutoHook
     procedure :: radius           => isothermalRadius
     procedure :: radiusGrowthRate => isothermalRadiusGrowthRate
     procedure :: calculationReset => isothermalCalculationReset
  end type coolingRadiusIsothermal

  interface coolingRadiusIsothermal
     !!{
     Constructors for the isothermal cooling radius class.
     !!}
     module procedure isothermalConstructorParameters
     module procedure isothermalConstructorInternal
  end interface coolingRadiusIsothermal

contains

  function isothermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the isothermal cooling radius class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingRadiusIsothermal  )                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(coolingTimeAvailableClass), pointer       :: coolingTimeAvailable_
    class(coolingTimeClass         ), pointer       :: coolingTime_
    class(darkMatterHaloScaleClass ), pointer       :: darkMatterHaloScale_
    class(cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="coolingTimeAvailable" name="coolingTimeAvailable_" source="parameters"/>
    <objectBuilder class="coolingTime"          name="coolingTime_"          source="parameters"/>
    !!]
    self=coolingRadiusIsothermal(cosmologyFunctions_,darkMatterHaloScale_,coolingTimeAvailable_,coolingTime_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="coolingTimeAvailable_"/>
    <objectDestructor name="coolingTime_"         />
    !!]
    return
  end function isothermalConstructorParameters

  function isothermalConstructorInternal(cosmologyFunctions_,darkMatterHaloScale_,coolingTimeAvailable_,coolingTime_) result(self)
    !!{
    Internal constructor for the isothermal cooling radius class.
    !!}
    use :: Abundances_Structure         , only : Abundances_Property_Count, abundances
    use :: Array_Utilities              , only : operator(.intersection.)
    use :: Chemical_Abundances_Structure, only : Chemicals_Property_Count
    use :: Error                        , only : Component_List           , Error_Report
    use :: Galacticus_Nodes             , only : defaultHotHaloComponent
    implicit none
    type (coolingRadiusIsothermal  )                        :: self
    class(cosmologyFunctionsClass  ), intent(in   ), target :: cosmologyFunctions_
    class(darkMatterHaloScaleClass ), intent(in   ), target :: darkMatterHaloScale_
    class(coolingTimeAvailableClass), intent(in   ), target :: coolingTimeAvailable_
    class(coolingTimeClass         ), intent(in   ), target :: coolingTime_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *darkMatterHaloScale_, *coolingTimeAvailable_, *coolingTime_"/>
    !!]

    ! Initial state of stored solutions.
    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    ! Get a count of the number of abundances and chemicals properties.
    self%abundancesCount=Abundances_Property_Count()
    self%chemicalsCount =Chemicals_Property_Count ()
    ! Initialize radiation field.
    allocate(self%radiation)
    !![
    <referenceConstruct isResult="yes" owner="self" object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
    !!]
    ! Check that required components are gettable.
    if     (                                                                                                             &
         &  .not.(                                                                                                       &
         &         defaultHotHaloComponent%       massIsGettable() .and.                                                 &
         &         defaultHotHaloComponent% abundancesIsGettable() .and.                                                 &
         &         defaultHotHaloComponent%outerRadiusIsGettable() .and.                                                 &
         &        (defaultHotHaloComponent%  chemicalsIsGettable() .or.  self%chemicalsCount == 0)                       &
         &       )                                                                                                       &
         & ) call Error_Report                                                                                           &
         & (                                                                                                             &
         &  'This method requires that the "mass", "abundances", "outerRadius", and "chemicals" '//                      &
         &  '(if any chemicals are being used) properties of the hot halo are gettable.'         //                      &
         &  Component_List(                                                                                              &
         &                 'hotHalo'                                                                                  ,  &
         &                  defaultHotHaloComponent%massAttributeMatch       (requireGettable=.true.                 )   &
         &                 .intersection.                                                                                &
         &                  defaultHotHaloComponent%abundancesAttributeMatch (requireGettable=.true.                 )   &
         &                 .intersection.                                                                                &
         &                  defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.                 )   &
         &                 .intersection.                                                                                &
         &                  defaultHotHaloComponent%chemicalsAttributeMatch  (requireGettable=self%chemicalsCount > 0)   &
         &                )                                                                                           // &
         &  {introspection:location}                                                                                     &
         & )
    return
  end function isothermalConstructorInternal

  subroutine isothermalAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(coolingRadiusIsothermal), intent(inout) :: self

    call calculationResetEvent%attach(self,isothermalCalculationReset,openMPThreadBindingAllLevels,label='coolingRadiusIsothermal')
    return
  end subroutine isothermalAutoHook

  subroutine isothermalDestructor(self)
    !!{
    Destructor for the isothermal cooling radius class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(coolingRadiusIsothermal), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%coolingTimeAvailable_"/>
    <objectDestructor name="self%coolingTime_"         />
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%radiation"            />
    !!]
    if (calculationResetEvent%isAttached(self,isothermalCalculationReset)) call calculationResetEvent%detach(self,isothermalCalculationReset)
    return
  end subroutine isothermalDestructor

  subroutine isothermalCalculationReset(self,node,uniqueID)
    !!{
    Reset the cooling radius calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (coolingRadiusIsothermal), intent(inout) :: self
    type   (treeNode               ), intent(inout) :: node
    integer(kind_int8              ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    self%lastUniqueID            =uniqueID
    return
  end subroutine isothermalCalculationReset

  double precision function isothermalRadiusGrowthRate(self,node)
    !!{
    Returns the cooling radius growth rate (in Mpc/Gyr) in the hot atmosphere.
    !!}
    implicit none
    class           (coolingRadiusIsothermal), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    double precision                                         :: radiusCooling, radiusVirial

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())

    ! Check if cooling radius growth rate is already computed.
    if (.not.self%radiusGrowthRateComputed) then
      ! Get the cooling radius.
       radiusCooling=self                     %      radius(node)
       ! Get the virial radius.
       radiusVirial =self%darkMatterHaloScale_%radiusVirial(node)
       ! Check if cooling radius has reached virial radius.
       if (radiusCooling >= radiusVirial) then
          self%radiusGrowthRateStored=0.0d0
       else
          ! Compute the growth rate of the cooling radius.
          self%radiusGrowthRateStored=+0.5d0                                                      &
               &                      *radiusCooling                                              &
               &                      *self%coolingTimeAvailable_%timeAvailableIncreaseRate(node) &
               &                      /self%coolingTimeAvailable_%timeAvailable            (node)
       end if
       ! Flag that cooling radius is now computed.
       self%radiusGrowthRateComputed=.true.
    end if
    ! Return the stored value.
    isothermalRadiusGrowthRate=self%radiusGrowthRateStored
    return
  end function isothermalRadiusGrowthRate

  double precision function isothermalRadius(self,node)
    !!{
    Return the cooling radius in the isothermal model.
    !!}
    use :: Abundances_Structure             , only : abundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Fraction_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                   , nodeComponentHotHalo       , treeNode
    use :: Mass_Distributions               , only : massDistributionClass                , kinematicsDistributionClass
    use :: Coordinates                      , only : coordinateSpherical                  , assignment(=)
    use :: Galactic_Structure_Options       , only : componentTypeHotHalo                 , massTypeGaseous
    implicit none
    class           (coolingRadiusIsothermal    ), intent(inout), target  :: self
    type            (treeNode                   ), intent(inout), target  :: node
    class           (nodeComponentBasic         )               , pointer :: basic
    class           (nodeComponentHotHalo       )               , pointer :: hotHalo
    class           (massDistributionClass      )               , pointer :: massDistribution_
    class           (kinematicsDistributionClass)               , pointer :: kinematicsDistribution_
    type            (coordinateSpherical        )                         :: coordinates
    double precision                                                      :: coolingTime            , timeAvailable          , &
         &                                                                   density                , massToDensityConversion, &
         &                                                                   temperature            , radiusVirial
    type            (abundances                 )                         :: hotAbundances
    type            (chemicalAbundances         )                         :: chemicalFractions      , chemicalMasses

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if cooling radius is already computed.
    if (.not.self%radiusComputed) then
       ! Get the time available for cooling in node.
       timeAvailable                   =  self%coolingTimeAvailable_%timeAvailable(node)
       ! Get the abundances for this node.
       hotHalo                         => node    %hotHalo  ()
       hotAbundances                   =  hotHalo%abundances()
       call hotAbundances%massToMassFraction(hotHalo%mass())
       ! Get the chemicals for this node.
       if (self%chemicalsCount > 0) then
          chemicalMasses=hotHalo%chemicals()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(chemicalFractions)
          if (hotHalo%mass() > 0.0d0) then
             massToDensityConversion=Chemicals_Mass_To_Fraction_Conversion(hotHalo%mass())
          else
             massToDensityConversion=0.0d0
          end if          
          ! Convert to number density per unit total mass density.
          chemicalFractions=chemicalFractions*massToDensityConversion
       end if
       ! Set epoch for radiation field.
       basic => node%basic()
       call self%radiation%timeSet(basic%time())
       ! Get the virial radius.
       radiusVirial=self%darkMatterHaloScale_%radiusVirial(node)
       ! Get the mass distribution.
       massDistribution_       => node             %massDistribution      (componentTypeHotHalo,massTypeGaseous)
       kinematicsDistribution_ => massDistribution_%kinematicsDistribution(                                    )
       ! Compute density, temperature and abundances.
       coordinates=[radiusVirial,0.0d0,0.0d0]
       density    =massDistribution_                   %density    (coordinates)
       temperature=kinematicsDistribution_             %temperature(coordinates)
       !![
       <objectDestructor name="massDistribution_"      />
       <objectDestructor name="kinematicsDistribution_"/>
       !!]          
       ! Compute the cooling time at the virial radius.
       coolingTime =self                  %coolingTime_%time       (node,temperature,density,hotAbundances,chemicalFractions*density,self%radiation)
       if (coolingTime < timeAvailable) then
          ! Cooling time available exceeds cooling time at virial radius, return virial radius.
          self%radiusStored=radiusVirial
       else
          ! Cooling radius is between zero and virial radii.
          self%radiusStored=radiusVirial*sqrt(timeAvailable/coolingTime)
       end if
       ! Flag that cooling radius is now computed.
       self%radiusComputed=.true.
    end if
    isothermalRadius=self%radiusStored
    return
  end function isothermalRadius
