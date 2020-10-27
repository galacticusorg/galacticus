!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of a simple cooling radius class.

  use :: Abundances_Structure         , only : abundances
  use :: Chemical_Abundances_Structure, only : chemicalAbundances
  use :: Cooling_Times                , only : coolingTimeClass
  use :: Cooling_Times_Available      , only : coolingTimeAvailableClass
  use :: Cosmology_Functions          , only : cosmologyFunctions                     , cosmologyFunctionsClass
  use :: Hot_Halo_Mass_Distributions  , only : hotHaloMassDistributionClass
  use :: Hot_Halo_Temperature_Profiles, only : hotHaloTemperatureProfileClass
  use :: Kind_Numbers                 , only : kind_int8
  use :: Radiation_Fields             , only : radiationFieldCosmicMicrowaveBackground

  !# <coolingRadius name="coolingRadiusSimple">
  !#  <description>
  !#   A cooling radius class that computes the cooling radius by seeking the radius at which the time available for cooling (see
  !#   \S\ref{sec:TimeAvailableCooling}) equals the cooling time (see \S\ref{sec:CoolingTime}). The growth rate is determined
  !#   consistently based on the slope of the density profile, the density dependence of the cooling function and the rate at
  !#   which the time available for cooling is increasing. This method assumes that the cooling time is a monotonic function of
  !#   radius.
  !#  </description>
  !#  <deepCopy>
  !#   <functionClass variables="radiation"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="radiation"/>
  !#  </stateStorable>
  !# </coolingRadius>
  type, extends(coolingRadiusClass) :: coolingRadiusSimple
     !% Implementation of cooling radius class in which the cooling radius is defined as that radius at which the time available
     !% for cooling equals the cooling time.
     private
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_        => null()
     class           (coolingTimeClass                       ), pointer :: coolingTime_               => null()
     class           (coolingTimeAvailableClass              ), pointer :: coolingTimeAvailable_      => null()
     class           (hotHaloMassDistributionClass           ), pointer :: hotHaloMassDistribution_   => null()
     class           (hotHaloTemperatureProfileClass         ), pointer :: hotHaloTemperatureProfile_ => null()
     type            (radiationFieldCosmicMicrowaveBackground), pointer :: radiation                  => null()
     integer         (kind=kind_int8                         )          :: lastUniqueID               =  -1
     integer                                                            :: abundancesCount                     , chemicalsCount
     ! Stored values of cooling radius.
     logical                                                            :: radiusComputed                      , radiusGrowthRateComputed
     double precision                                                   :: radiusGrowthRateStored              , radiusStored
   contains
     !@ <objectMethods>
     !@   <object>coolingRadiusSimple</object>
     !@   <objectMethod>
     !@     <method>calculationReset</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(table)\textgreater} node\arginout</arguments>
     !@     <description>Reset memoized calculations.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                     simpleDestructor
     procedure :: autoHook         => simpleAutoHook
     procedure :: radius           => simpleRadius
     procedure :: radiusGrowthRate => simpleRadiusGrowthRate
     procedure :: calculationReset => simpleCalculationReset
  end type coolingRadiusSimple

  interface coolingRadiusSimple
     !% Constructors for the simple cooling radius class.
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface coolingRadiusSimple

  ! Module scope variables used in root finding.
  class           (coolingRadiusSimple), pointer :: simpleSelf_
  type            (treeNode           ), pointer :: simpleNode_
  double precision                               :: simpleCoolingTimeAvailable_
  type            (abundances         )          :: simpleGasAbundances_
  type            (chemicalAbundances )          :: simpleChemicalDensities_
  !$omp threadprivate(simpleSelf_,simpleNode_,simpleCoolingTimeAvailable_,simpleGasAbundances_,simpleChemicalDensities_)

contains

  function simpleConstructorParameters(parameters) result(self)
    !% Constructor for the simple cooling radius class which builds the object from a parameter set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingRadiusSimple           )                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(coolingTimeAvailableClass     ), pointer       :: coolingTimeAvailable_
    class(coolingTimeClass              ), pointer       :: coolingTime_
    class(hotHaloTemperatureProfileClass), pointer       :: hotHaloTemperatureProfile_
    class(hotHaloMassDistributionClass  ), pointer       :: hotHaloMassDistribution_
    class(cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_

    !# <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    !# <objectBuilder class="coolingTimeAvailable"      name="coolingTimeAvailable_"      source="parameters"/>
    !# <objectBuilder class="coolingTime"               name="coolingTime_"               source="parameters"/>
    !# <objectBuilder class="hotHaloTemperatureProfile" name="hotHaloTemperatureProfile_" source="parameters"/>
    !# <objectBuilder class="hotHaloMassDistribution"   name="hotHaloMassDistribution_"   source="parameters"/>
    self=coolingRadiusSimple(cosmologyFunctions_,coolingTimeAvailable_,coolingTime_,hotHaloTemperatureProfile_,hotHaloMassDistribution_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"       />
    !# <objectDestructor name="coolingTimeAvailable_"     />
    !# <objectDestructor name="coolingTime_"              />
    !# <objectDestructor name="hotHaloTemperatureProfile_"/>
    !# <objectDestructor name="hotHaloMassDistribution_"  />
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(cosmologyFunctions_,coolingTimeAvailable_,coolingTime_,hotHaloTemperatureProfile_,hotHaloMassDistribution_) result(self)
    !% Internal constructor for the simple cooling radius class.
    use :: Abundances_Structure         , only : Abundances_Property_Count, abundances
    use :: Array_Utilities              , only : operator(.intersection.)
    use :: Chemical_Abundances_Structure, only : Chemicals_Property_Count
    use :: Galacticus_Error             , only : Galacticus_Component_List, Galacticus_Error_Report
    use :: Galacticus_Nodes             , only : defaultHotHaloComponent
    implicit none
    type (coolingRadiusSimple           )                        :: self
    class(cosmologyFunctionsClass       ), intent(in   ), target :: cosmologyFunctions_
    class(coolingTimeAvailableClass     ), intent(in   ), target :: coolingTimeAvailable_
    class(coolingTimeClass              ), intent(in   ), target :: coolingTime_
    class(hotHaloTemperatureProfileClass), intent(in   ), target :: hotHaloTemperatureProfile_
    class(hotHaloMassDistributionClass  ), intent(in   ), target :: hotHaloMassDistribution_
    !# <constructorAssign variables="*cosmologyFunctions_, *coolingTimeAvailable_, *coolingTime_, *hotHaloTemperatureProfile_, *hotHaloMassDistribution_"/>

    ! Initial state of stored solutions.
    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    ! Get a count of the number of abundances and chemicals properties.
    self%abundancesCount=Abundances_Property_Count()
    self%chemicalsCount =Chemicals_Property_Count ()
    ! Initialize radiation field.
    allocate(self%radiation)
    !# <referenceConstruct isResult="yes" owner="self" object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
    ! Check that required components are gettable.
    if     (                                                                                                                        &
         &  .not.(                                                                                                                  &
         &         defaultHotHaloComponent%       massIsGettable() .and.                                                            &
         &         defaultHotHaloComponent% abundancesIsGettable() .and.                                                            &
         &         defaultHotHaloComponent%outerRadiusIsGettable() .and.                                                            &
         &        (defaultHotHaloComponent%  chemicalsIsGettable() .or.  self%chemicalsCount == 0)                                  &
         &       )                                                                                                                  &
         & ) call Galacticus_Error_Report                                                                                           &
         & (                                                                                                                        &
         &  'This method requires that the "mass", "abundances", "outerRadius", and "chemicals" '//                                 &
         &  '(if any chemicals are being used) properties of the hot halo are gettable.'         //                                 &
         &  Galacticus_Component_List(                                                                                              &
         &                            'hotHalo'                                                                                  ,  &
         &                             defaultHotHaloComponent%massAttributeMatch       (requireGettable=.true.                 )   &
         &                            .intersection.                                                                                &
         &                             defaultHotHaloComponent%abundancesAttributeMatch (requireGettable=.true.                 )   &
         &                            .intersection.                                                                                &
         &                             defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.                 )   &
         &                            .intersection.                                                                                &
         &                             defaultHotHaloComponent%chemicalsAttributeMatch  (requireGettable=self%chemicalsCount > 0)   &
         &                           )                                                                                           // &
         &  {introspection:location}                                                                                                &
         & )
    return
  end function simpleConstructorInternal

  subroutine simpleAutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(coolingRadiusSimple), intent(inout) :: self

    call calculationResetEvent%attach(self,simpleCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine simpleAutoHook

  subroutine simpleDestructor(self)
    !% Destructor for the simple cooling radius class.
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(coolingRadiusSimple), intent(inout) :: self

    !# <objectDestructor name="self%coolingTimeAvailable_"     />
    !# <objectDestructor name="self%coolingTime_"              />
    !# <objectDestructor name="self%hotHaloTemperatureProfile_"/>
    !# <objectDestructor name="self%hotHaloMassDistribution_"  />
    !# <objectDestructor name="self%cosmologyFunctions_"       />
    !# <objectDestructor name="self%radiation"                 />
    call calculationResetEvent%detach(self,simpleCalculationReset)
    return
  end subroutine simpleDestructor

  subroutine simpleCalculationReset(self,node)
    !% Reset the cooling radius calculation.
    implicit none
    class(coolingRadiusSimple), intent(inout) :: self
    type (treeNode           ), intent(inout) :: node

    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    self%lastUniqueID            =node%uniqueID()
    return
  end subroutine simpleCalculationReset

  double precision function simpleRadiusGrowthRate(self,node)
    !% Returns the cooling radius growth rate (in Mpc/Gyr) in the hot atmosphere.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentHotHalo, treeNode
    implicit none
    class           (coolingRadiusSimple ), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node
    class           (nodeComponentBasic  ), pointer       :: basic
    class           (nodeComponentHotHalo), pointer       :: hotHalo
    double precision                                      :: coolingRadius                   , coolingTimeAvailable      , &
         &                                                   coolingTimeAvailableIncreaseRate, coolingTimeDensityLogSlope, &
         &                                                   coolingTimeTemperatureLogSlope  , density                   , &
         &                                                   densityLogSlope                 , outerRadius               , &
         &                                                   temperature                     , temperatureLogSlope

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
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
          density            =self%hotHaloMassDistribution_  %density            (node,coolingRadius)
          temperature        =self%hotHaloTemperatureProfile_%temperature        (node,coolingRadius)
          densityLogSlope    =self%hotHaloMassDistribution_  %densityLogSlope    (node,coolingRadius)
          temperatureLogSlope=self%hotHaloTemperatureProfile_%temperatureLogSlope(node,coolingRadius)
          ! Get the time available for cooling in node and its rate of increase.
          coolingTimeAvailable            =self%coolingTimeAvailable_%timeAvailable            (node)
          coolingTimeAvailableIncreaseRate=self%coolingTimeAvailable_%timeAvailableIncreaseRate(node)
          ! Get gradients of cooling time with density and temperature.
          coolingTimeDensityLogSlope    =self%coolingTime_%gradientDensityLogarithmic    (temperature,density,simpleGasAbundances_,simpleChemicalDensities_,self%radiation)
          coolingTimeTemperatureLogSlope=self%coolingTime_%gradientTemperatureLogarithmic(temperature,density,simpleGasAbundances_,simpleChemicalDensities_,self%radiation)
          ! Compute rate at which cooling radius grows.
          if (coolingRadius > 0.0d0) then
             self%radiusGrowthRateStored=+coolingRadius                                        &
                  &                      /coolingTimeAvailable                                 &
                  &                      *coolingTimeAvailableIncreaseRate                     &
                  &                      /(                                                    &
                  &                        +    densityLogSlope*coolingTimeDensityLogSlope     &
                  &                        +temperatureLogSlope*coolingTimeTemperatureLogSlope &
                  &                       )
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
    !% Return the cooling radius in the simple model.
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                  , nodeComponentHotHalo, treeNode
    use :: Root_Finder                      , only : rootFinder
    implicit none
    class           (coolingRadiusSimple ), intent(inout), target :: self
    type            (treeNode            ), intent(inout), target :: node
    class           (nodeComponentBasic  ), pointer               :: basic
    class           (nodeComponentHotHalo), pointer               :: hotHalo
    double precision                      , parameter             :: zeroRadius       =0.0d0
    double precision                      , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative      =1.0d-6
    type            (rootFinder          ), save                  :: finder
    !$omp threadprivate(finder)
    type            (chemicalAbundances  )                        :: chemicalMasses
    double precision                                              :: outerRadius            , massToDensityConversion

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if cooling radius is already computed.
    if (.not.self%radiusComputed) then
       ! Flag that cooling radius is now computed.
       self%radiusComputed=.true.
       ! Get the time available for cooling in node.
       simpleCoolingTimeAvailable_=self%coolingTimeAvailable_%timeAvailable(node)
       ! Get node components.
       hotHalo => node%hotHalo()
       ! Get the abundances for this node.
       simpleGasAbundances_=hotHalo%abundances()
       call simpleGasAbundances_%massToMassFraction(hotHalo%mass())
       ! Get the chemicals for this node.
       if (self%chemicalsCount > 0) then
          chemicalMasses=hotHalo%chemicals()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(simpleChemicalDensities_)
          ! Compute factor converting mass of chemicals in (M_Solar/M_Atomic) to number density in cm^-3.
          if (hotHalo%outerRadius() > 0.0d0) then
             massToDensityConversion=Chemicals_Mass_To_Density_Conversion(hotHalo%outerRadius())
          else
             massToDensityConversion=0.0d0
          end if
          ! Convert to number density.
          simpleChemicalDensities_=simpleChemicalDensities_*massToDensityConversion
       end if
       ! Set epoch for radiation field.
       basic => node%basic()
       call self%radiation%timeSet(basic%time())
       ! Set module-scope pointers.
       simpleSelf_ => self
       simpleNode_ => node
       ! Check if cooling time at halo center is reached.
       if (coolingRadiusRoot(zeroRadius) > 0.0d0) then
          ! Cooling time at halo center exceeds the time available, return zero radius.
          self%radiusStored=zeroRadius
          simpleRadius     =self%radiusStored
          return
       end if
       ! Check if cooling time at hot halo outer radius is reached.
       outerRadius=hotHalo%outerRadius()
       if (coolingRadiusRoot(outerRadius) < 0.0d0) then
          ! Cooling time available exceeds cooling time at outer radius radius, return outer radius.
          self%radiusStored=outerRadius
          simpleRadius     =self%radiusStored
          return
       end if
       ! Cooling radius is between zero and outer radii. Search for the cooling radius.
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(coolingRadiusRoot                  )
          call finder%tolerance   (toleranceAbsolute,toleranceRelative)
       end if
       self%radiusStored =  finder%find(rootRange=[zeroRadius,outerRadius])
       simpleRadius      =  self%radiusStored
    else
       simpleRadius      =  self%radiusStored
    end if
    return
  end function simpleRadius

  double precision function coolingRadiusRoot(radius)
    !% Root function which evaluates the difference between the cooling time at {\normalfont \ttfamily radius} and the time available for cooling.
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: coolingTime, density, temperature

    ! Compute density, temperature and abundances.
    density    =simpleSelf_%hotHaloMassDistribution_  %density    (simpleNode_,radius)
    temperature=simpleSelf_%hotHaloTemperatureProfile_%temperature(simpleNode_,radius)
    ! Compute the cooling time at the specified radius.
    coolingTime=simpleSelf_%coolingTime_              %time       (temperature,density,simpleGasAbundances_,simpleChemicalDensities_,simpleSelf_%radiation)
    ! Return the difference between cooling time and time available.
    coolingRadiusRoot=coolingTime-simpleCoolingTimeAvailable_
    return
  end function coolingRadiusRoot
