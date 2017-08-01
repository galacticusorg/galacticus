!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  
  use Kind_Numbers
  use Cooling_Times_Available
  use Hot_Halo_Mass_Distributions
  use Hot_Halo_Temperature_Profiles
  use Radiation_Structure
  use Abundances_Structure
  use Chemical_Abundances_Structure

  !# <coolingRadius name="coolingRadiusSimple" defaultThreadPrivate="yes">
  !#  <description>
  !#   A cooling radius class Computes the cooling radius by seeking the radius at which the time available for cooling equals the
  !#   cooling time. The growth rate is determined consistently based on the slope of the density profile, the density dependence
  !#   of the cooling function and the rate at which the time available for cooling is increasing. This method assumes that the
  !#   cooling time is a monotonic function of radius.
  !#  </description>
  !# </coolingRadius>
  type, extends(coolingRadiusClass) :: coolingRadiusSimple
     !% Implementation of cooling radius class in which the cooling radius is defined as that radius at which the time available
     !% for cooling equals the cooling time.
     private
     class           (coolingTimeAvailableClass     ), pointer :: coolingTimeAvailable_
     class           (hotHaloMassDistributionClass  ), pointer :: hotHaloMassDistribution_
     class           (hotHaloTemperatureProfileClass), pointer :: hotHaloTemperatureProfile_
     integer         (kind=kind_int8                )          :: lastUniqueID                  =-1
     integer                                                   :: abundancesCount                  , chemicalsCount
     ! Stored values of cooling radius.
     logical                                                   :: radiusComputed                   , radiusGrowthRateComputed
     double precision                                          :: radiusGrowthRateStored           , radiusStored
   contains
     final     ::                     simpleDestructor
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
  type            (radiationStructure )          :: simpleRadiation_
  type            (chemicalAbundances )          :: simpleChemicalDensities_
  !$omp threadprivate(simpleSelf_,simpleNode_,simpleCoolingTimeAvailable_,simpleGasAbundances_,simpleRadiation_,simpleChemicalDensities_)

contains

  function simpleConstructorParameters(parameters) result(self)
    !% Constructor for the simple cooling radius class which builds the object from a parameter set.
    use Input_Parameters2
    implicit none
    type (coolingRadiusSimple      )                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(coolingTimeAvailableClass), pointer       :: coolingTimeAvailable_

    !# <objectBuilder class="coolingTimeAvailable" name="coolingTimeAvailable_" source="parameters"/>
    self=coolingRadiusSimple(coolingTimeAvailable_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(coolingTimeAvailable_) result(self)
    !% Internal constructor for the simple cooling radius class.
    use ISO_Varying_String
    use Galacticus_Error
    use Array_Utilities
    use String_Handling
    use Abundances_Structure
    use Chemical_Abundances_Structure
    implicit none
    type (coolingRadiusSimple      )                        :: self
    class(coolingTimeAvailableClass), intent(in   ), target :: coolingTimeAvailable_
    !# <constructorAssign variables="*coolingTimeAvailable_"/>
    
    ! Initial state of stored solutions.
    self%radiusComputed          =.false.
    self%radiusGrowthRateComputed=.false.
    ! Get a count of the number of abundances and chemicals properties.
    self%abundancesCount=Abundances_Property_Count()
    self%chemicalsCount =Chemicals_Property_Count ()
    ! Check that required components are gettable.
    if     (                                                                                                                       &
         &  .not.(                                                                                                                 &
         &         defaultHotHaloComponent%       massIsGettable() .and.                                                           &
         &         defaultHotHaloComponent% abundancesIsGettable() .and.                                                           &
         &         defaultHotHaloComponent%outerRadiusIsGettable() .and.                                                           &
         &        (defaultHotHaloComponent%  chemicalsIsGettable() .or.  self%chemicalsCount == 0)                                 &
         &       )                                                                                                                 &
         & ) call Galacticus_Error_Report                                                                                          &
         & (                                                                                                                       &
         &  'simpleConstructorInternal'                                                                                          , &
         &  'This method requires that the "mass", "abundances", "outerRadius", and "chemicals" '//                                &
         &  '(if any chemicals are being used) properties of the hot halo are gettable.'         //                                &
         &  Galacticus_Component_List(                                                                                             &
         &                            'hotHalo'                                                                                  , &
         &                             defaultHotHaloComponent%massAttributeMatch       (requireGettable=.true.                 )  &
         &                            .intersection.                                                                               &
         &                             defaultHotHaloComponent%abundancesAttributeMatch (requireGettable=.true.                 )  &
         &                            .intersection.                                                                               &
         &                             defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.                 )  &
         &                            .intersection.                                                                               &
         &                             defaultHotHaloComponent%chemicalsAttributeMatch  (requireGettable=self%chemicalsCount > 0)  &
         &                           )                                                                                             &
         & )
    return
  end function simpleConstructorInternal
  
  subroutine simpleDestructor(self)
    !% Destructor for the simple cooling radius class.
    implicit none
    type(coolingRadiusSimple), intent(inout) :: self

    !# <objectDestructor name="self%coolingTimeAvailable_"/>
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
    use Radiation_Structure
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Cooling_Times
    implicit none
    class           (coolingRadiusSimple ), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node
    class           (nodeComponentHotHalo), pointer       :: hotHalo
    double precision                                      :: coolingRadius                   , coolingTimeAvailable      , &
         &                                                   coolingTimeAvailableIncreaseRate, coolingTimeDensityLogSlope, &
         &                                                   coolingTimeTemperatureLogSlope  , density                   , &
         &                                                   densityLogSlope                 , outerRadius               , &
         &                                                   temperature                     , temperatureLogSlope
    type            (abundances          )                :: gasAbundances
    type            (chemicalAbundances  )                :: chemicalDensities
    type            (radiationStructure  )                :: radiation

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
          ! Get required objects.
          self%hotHaloTemperatureProfile_ => hotHaloTemperatureProfile()
          self%hotHaloMassDistribution_   => hotHaloMassDistribution  ()
          ! Define radiation field.
          call radiation%define([radiationTypeCMB])
          ! Get density and temperature  at the cooling radius, plus their gradients.
          density            =self%hotHaloMassDistribution_  %density            (node,coolingRadius)
          temperature        =self%hotHaloTemperatureProfile_%temperature        (node,coolingRadius)
          densityLogSlope    =self%hotHaloMassDistribution_  %densityLogSlope    (node,coolingRadius)
          temperatureLogSlope=self%hotHaloTemperatureProfile_%temperatureLogSlope(node,coolingRadius)
          ! Get the time available for cooling in node and its rate of increase.
          coolingTimeAvailable            =self%coolingTimeAvailable_%timeAvailable            (node)
          coolingTimeAvailableIncreaseRate=self%coolingTimeAvailable_%timeAvailableIncreaseRate(node)
          ! Get gradients of cooling time with density and temperature.
          coolingTimeDensityLogSlope    =Cooling_Time_Density_Log_Slope    (temperature,density,gasAbundances,chemicalDensities,radiation)
          coolingTimeTemperatureLogSlope=Cooling_Time_Temperature_Log_Slope(temperature,density,gasAbundances,chemicalDensities,radiation)
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
   use Chemical_Reaction_Rates_Utilities
    use Root_Finder
    implicit none
    class           (coolingRadiusSimple ), intent(inout), target :: self
    type            (treeNode            ), intent(inout), target :: node
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
       ! Get required objects.
       self%hotHaloMassDistribution_   => hotHaloMassDistribution  ()
       self%hotHaloTemperatureProfile_ => hotHaloTemperatureProfile()
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
       ! Set the radiation field.
       call simpleRadiation_%define([radiationTypeCMB])
       call simpleRadiation_%set(node)
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
    use Cooling_Times
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: coolingTime, density, temperature

    ! Compute density, temperature and abundances.
    density    =simpleSelf_%hotHaloMassDistribution_  %density    (simpleNode_,radius)
    temperature=simpleSelf_%hotHaloTemperatureProfile_%temperature(simpleNode_,radius)
    ! Compute the cooling time at the specified radius.
    coolingTime=Cooling_Time(temperature,density,simpleGasAbundances_,simpleChemicalDensities_,simpleRadiation_)
    ! Return the difference between cooling time and time available.
    coolingRadiusRoot=coolingTime-simpleCoolingTimeAvailable_
    return
  end function coolingRadiusRoot
