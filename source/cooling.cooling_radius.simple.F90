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
  Implementation of a simple cooling radius class.
  !!}

  use :: Abundances_Structure         , only : abundances
  use :: Chemical_Abundances_Structure, only : chemicalAbundances
  use :: Cooling_Times                , only : coolingTimeClass
  use :: Cooling_Times_Available      , only : coolingTimeAvailableClass
  use :: Cosmology_Functions          , only : cosmologyFunctions                     , cosmologyFunctionsClass
  use :: Kind_Numbers                 , only : kind_int8
  use :: Radiation_Fields             , only : radiationFieldCosmicMicrowaveBackground
  use :: Root_Finder                  , only : rootFinder

  !![
  <coolingRadius name="coolingRadiusSimple">
   <description>
    A cooling radius class that computes the cooling radius by seeking the radius at which the time available for cooling (see
    \refPhysics{coolingTimeAvailable}) equals the cooling time (see \refPhysics{coolingTime}). The growth rate is determined
    consistently based on the slope of the density profile, the density dependence of the cooling function and the rate at
    which the time available for cooling is increasing. This method assumes that the cooling time is a monotonic function of
    radius.
   </description>
   <deepCopy>
    <functionClass variables="radiation"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="radiation"/>
   </stateStorable>
  </coolingRadius>
  !!]
  type, extends(coolingRadiusClass) :: coolingRadiusSimple
     !!{
     Implementation of cooling radius class in which the cooling radius is defined as that radius at which the time available
     for cooling equals the cooling time.
     !!}
     private
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_        => null()
     class           (coolingTimeClass                       ), pointer :: coolingTime_               => null()
     class           (coolingTimeAvailableClass              ), pointer :: coolingTimeAvailable_      => null()
     type            (radiationFieldCosmicMicrowaveBackground), pointer :: radiation                  => null()
     type            (rootFinder                             )          :: finder
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
     final     ::                     simpleDestructor
     procedure :: autoHook         => simpleAutoHook
     procedure :: radius           => simpleRadius
     procedure :: radiusGrowthRate => simpleRadiusGrowthRate
     procedure :: calculationReset => simpleCalculationReset
  end type coolingRadiusSimple

  interface coolingRadiusSimple
     !!{
     Constructors for the simple cooling radius class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface coolingRadiusSimple

  ! Module scope variables used in root finding.
  class           (coolingRadiusSimple), pointer :: self_
  type            (treeNode           ), pointer :: node_
  double precision                               :: coolingTimeAvailable_
  type            (abundances         )          :: abundancesGas_
  type            (chemicalAbundances )          :: fractionsChemical_
  !$omp threadprivate(self_,node_,coolingTimeAvailable_,abundancesGas_,fractionsChemical_)

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the simple cooling radius class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingRadiusSimple      )                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(coolingTimeAvailableClass), pointer       :: coolingTimeAvailable_
    class(coolingTimeClass         ), pointer       :: coolingTime_
    class(cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="coolingTimeAvailable" name="coolingTimeAvailable_"      source="parameters"/>
    <objectBuilder class="coolingTime"          name="coolingTime_"               source="parameters"/>
    !!]
    self=coolingRadiusSimple(cosmologyFunctions_,coolingTimeAvailable_,coolingTime_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="coolingTimeAvailable_"/>
    <objectDestructor name="coolingTime_"         />
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(cosmologyFunctions_,coolingTimeAvailable_,coolingTime_) result(self)
    !!{
    Internal constructor for the simple cooling radius class.
    !!}
    use :: Abundances_Structure         , only : Abundances_Property_Count, abundances
    use :: Array_Utilities              , only : operator(.intersection.)
    use :: Chemical_Abundances_Structure, only : Chemicals_Property_Count
    use :: Error                        , only : Component_List           , Error_Report
    use :: Galacticus_Nodes             , only : defaultHotHaloComponent
    implicit none
    type            (coolingRadiusSimple      )                        :: self
    class           (cosmologyFunctionsClass  ), intent(in   ), target :: cosmologyFunctions_
    class           (coolingTimeAvailableClass), intent(in   ), target :: coolingTimeAvailable_
    class           (coolingTimeClass         ), intent(in   ), target :: coolingTime_
    double precision                           , parameter             :: toleranceAbsolute         =0.0d0, toleranceRelative=1.0d-6
    !![
    <constructorAssign variables="*cosmologyFunctions_, *coolingTimeAvailable_, *coolingTime_"/>
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
    ! Initialize a root finder.
    self%finder=rootFinder(                                     &
         &                 rootFunction     =coolingRadiusRoot, &
         &                 toleranceAbsolute=toleranceAbsolute, &
         &                 toleranceRelative=toleranceRelative  &
         &                )
    return
  end function simpleConstructorInternal

  subroutine simpleAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(coolingRadiusSimple), intent(inout) :: self

    call calculationResetEvent%attach(self,simpleCalculationReset,openMPThreadBindingAllLevels,label='coolingRadiusSimple')
    return
  end subroutine simpleAutoHook

  subroutine simpleDestructor(self)
    !!{
    Destructor for the simple cooling radius class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(coolingRadiusSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingTimeAvailable_"/>
    <objectDestructor name="self%coolingTime_"         />
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%radiation"            />
    !!]
    if (calculationResetEvent%isAttached(self,simpleCalculationReset)) call calculationResetEvent%detach(self,simpleCalculationReset)
    return
  end subroutine simpleDestructor

  subroutine simpleCalculationReset(self,node,uniqueID)
    !!{
    Reset the cooling radius calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (coolingRadiusSimple), intent(inout) :: self
    type   (treeNode           ), intent(inout) :: node
    integer(kind_int8          ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    self%lastUniqueID            =uniqueID
    return
  end subroutine simpleCalculationReset

  double precision function simpleRadiusGrowthRate(self,node)
    !!{
    Returns the cooling radius growth rate (in Mpc/Gyr) in the hot atmosphere.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic   , nodeComponentHotHalo       , treeNode
    use :: Mass_Distributions        , only : massDistributionClass, kinematicsDistributionClass
    use :: Coordinates               , only : coordinateSpherical  , assignment(=)
    use :: Galactic_Structure_Options, only : componentTypeHotHalo , massTypeGaseous
    implicit none
    class           (coolingRadiusSimple        ), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    class           (nodeComponentBasic         ), pointer       :: basic
    class           (nodeComponentHotHalo       ), pointer       :: hotHalo
    class           (massDistributionClass      ), pointer       :: massDistribution_
    class           (kinematicsDistributionClass), pointer       :: kinematicsDistribution_
    type            (coordinateSpherical        )                :: coordinates
    double precision                                             :: coolingRadius                   , coolingTimeAvailable      , &
         &                                                          coolingTimeAvailableIncreaseRate, coolingTimeDensityLogSlope, &
         &                                                          coolingTimeTemperatureLogSlope  , density                   , &
         &                                                          densityLogSlope                 , outerRadius               , &
         &                                                          temperature                     , temperatureLogSlope       , &
         &                                                          slope

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if cooling radius growth rate is already computed.
    if (.not.self%radiusGrowthRateComputed) then
       ! Flag that cooling radius is now computed.
       self%radiusGrowthRateComputed=.true.
       ! Get node components.
       hotHalo => node%hotHalo()
       ! Get the outer radius.
       outerRadius=hotHalo%outerRadius()
       ! Get the cooling radius.
       coolingRadius=self%radius(node)
       ! Check if cooling radius has reached the outer radius.
       if (coolingRadius >= outerRadius) then
          self%radiusGrowthRateStored=0.0d0
       else
          ! Set epoch for radiation field.
          basic => node%basic()
          call self%radiation%timeSet(basic%time())
          ! Get density and temperature  at the cooling radius, plus their gradients.
          coordinates             =  [coolingRadius,0.0d0,0.0d0]
          massDistribution_       => node                   %massDistribution              (componentTypeHotHalo,massTypeGaseous)
          kinematicsDistribution_ => massDistribution_      %kinematicsDistribution        (                                    )
          density                 =  massDistribution_      %density                       (coordinates                         )
          temperature             =  kinematicsDistribution_%temperature                   (coordinates                         )
          densityLogSlope         =  massDistribution_      %densityGradientRadial         (coordinates,logarithmic=.true.      )
          temperatureLogSlope     =  kinematicsDistribution_%temperatureGradientLogarithmic(coordinates                         )
          !![
	  <objectDestructor name="massDistribution_"      />
	  <objectDestructor name="kinematicsDistribution_"/>
          !!]          
          ! Get the time available for cooling in node and its rate of increase.
          coolingTimeAvailable            =self%coolingTimeAvailable_%timeAvailable            (node)
          coolingTimeAvailableIncreaseRate=self%coolingTimeAvailable_%timeAvailableIncreaseRate(node)
          ! Get gradients of cooling time with density and temperature.
          coolingTimeDensityLogSlope    =self%coolingTime_%gradientDensityLogarithmic    (node,temperature,density,abundancesGas_,fractionsChemical_*density,self%radiation)
          coolingTimeTemperatureLogSlope=self%coolingTime_%gradientTemperatureLogarithmic(node,temperature,density,abundancesGas_,fractionsChemical_*density,self%radiation)
          ! Compute rate at which cooling radius grows.
          if (coolingRadius > 0.0d0) then
             slope                      =+    densityLogSlope*coolingTimeDensityLogSlope     &
                  &                      +temperatureLogSlope*coolingTimeTemperatureLogSlope
             if (slope /= 0.0d0) then
                self%radiusGrowthRateStored=+coolingRadius                                      &
                     &                      /coolingTimeAvailable                               &
                     &                      *coolingTimeAvailableIncreaseRate                   &
                     &                      /slope
             else
                ! The slope of the cooling radius vs. time available curve can be zero in cases where there are discontinuities in
                ! the cooling function. Catch such behavior here.
                self%radiusGrowthRateStored=0.0d0
             end if
          else
             self%radiusGrowthRateStored=0.0d0
          end if
       end if
    end if
    ! Return the stored value.
    simpleRadiusGrowthRate=self%radiusGrowthRateStored
    return
  end function simpleRadiusGrowthRate

  double precision function simpleRadius(self,node)
    !!{
    Return the cooling radius in the simple model.
    !!}
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Fraction_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                   , nodeComponentHotHalo, treeNode
    implicit none
    class           (coolingRadiusSimple ), intent(inout), target :: self
    type            (treeNode            ), intent(inout), target :: node
    class           (nodeComponentBasic  ), pointer               :: basic
    class           (nodeComponentHotHalo), pointer               :: hotHalo
    double precision                      , parameter             :: zeroRadius    =0.0d0
    type            (chemicalAbundances  )                        :: chemicalMasses
    double precision                                              :: outerRadius         , massToDensityConversion, &
         &                                                           rootZero            , rootOuter

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if cooling radius is already computed.
    if (.not.self%radiusComputed) then
       ! Flag that cooling radius is now computed.
       self%radiusComputed=.true.
       ! Get the time available for cooling in node.
       coolingTimeAvailable_=self%coolingTimeAvailable_%timeAvailable(node)
       ! Get node components.
       hotHalo => node%hotHalo()
       ! Get the abundances for this node.
       abundancesGas_=hotHalo%abundances()
       call abundancesGas_%massToMassFraction(hotHalo%mass())
       ! Get the chemicals for this node.
       if (self%chemicalsCount > 0) then
          chemicalMasses=hotHalo%chemicals()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(fractionsChemical_)
          ! Compute factor converting mass of chemicals in (M☉) to number density per unit total mass density (in cm⁻³ / M☉
          ! Mpc⁻³).
          if (hotHalo%mass() > 0.0d0) then
             massToDensityConversion=Chemicals_Mass_To_Fraction_Conversion(hotHalo%mass())
          else
             massToDensityConversion=0.0d0
          end if          
          ! Convert to number density per unit total mass density.
          call fractionsChemical_%scale(massToDensityConversion)
       end if
       ! Set epoch for radiation field.
       basic => node%basic()
       call self%radiation%timeSet(basic%time())
       ! Set module-scope pointers.
       self_ => self
       node_ => node
       ! Check if cooling time at hot halo outer radius is reached.
       outerRadius=hotHalo%outerRadius()
       rootOuter=coolingRadiusRoot(outerRadius)
       if (rootOuter < 0.0d0) then
          ! Cooling time available exceeds cooling time at outer radius radius, return outer radius.
          self%radiusStored=outerRadius
          simpleRadius     =self%radiusStored
          return
       end if
       ! Check if cooling time at halo center is reached.
       rootZero=coolingRadiusRoot(zeroRadius)
       if (rootZero > 0.0d0) then
          ! Cooling time at halo center exceeds the time available, return zero radius.
          self%radiusStored=zeroRadius
          simpleRadius     =self%radiusStored
          return
       end if
       ! Cooling radius is between zero and outer radii. Search for the cooling radius.
       self%radiusStored=self%finder      %find(rootRange=[zeroRadius,outerRadius],rootRangeValues=[rootZero,rootOuter])
       simpleRadius     =self%radiusStored
    else
       simpleRadius     =self%radiusStored
    end if
    return
  end function simpleRadius

  double precision function coolingRadiusRoot(radius)
    !!{
    Root function which evaluates the difference between the cooling time at {\normalfont \ttfamily radius} and the time available for cooling.
    !!}
    use :: Mass_Distributions        , only : massDistributionClass, kinematicsDistributionClass
    use :: Coordinates               , only : coordinateSpherical  , assignment(=)
    use :: Galactic_Structure_Options, only : componentTypeHotHalo , massTypeGaseous
    implicit none
    double precision                             , intent(in   ) :: radius
    double precision                                             :: coolingTime            , density, &
         &                                                          temperature
    class           (massDistributionClass      ), pointer       :: massDistribution_
    class           (kinematicsDistributionClass), pointer       :: kinematicsDistribution_
    type            (chemicalAbundances         ), save          :: densityChemicals
    !$omp threadprivate(densityChemicals)
    type            (coordinateSpherical        )                :: coordinates
   
    ! Compute density, temperature and abundances.
    coordinates             =   [radius,0.0d0,0.0d0]
    massDistribution_       =>  node_                  %massDistribution      (componentTypeHotHalo,massTypeGaseous)
    kinematicsDistribution_ =>  massDistribution_      %kinematicsDistribution(                                    )
    density                 =   massDistribution_      %density               (coordinates                         )
    if (associated(kinematicsDistribution_)) then
       temperature          =   kinematicsDistribution_%temperature           (coordinates                         )
    else
       temperature          =   0.0d0
    end if
    !![
    <objectDestructor name="massDistribution_"      />
    <objectDestructor name="kinematicsDistribution_"/>
    !!]          
    densityChemicals=fractionsChemical_
    call densityChemicals%scale(density)
    ! Compute the cooling time at the specified radius.
    coolingTime=self_%coolingTime_%time(node_,temperature,density,abundancesGas_,densityChemicals,self_%radiation)
    ! Return the difference between cooling time and time available.
    coolingRadiusRoot=coolingTime-coolingTimeAvailable_
    return
  end function coolingRadiusRoot
