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

!% Contains a module which implements a simple cooling radius calculation (finds the radius at which the time available for
!% cooling equals the cooling time).

module Cooling_Radii_Simple
  !% Implements a simple cooling radius calculation (finds the radius at which the time available for cooling equals the cooling
  !% time).
  use Galacticus_Nodes
  use Kind_Numbers
  use Radiation_Structure
  use Abundances_Structure
  use Chemical_Abundances_Structure
  implicit none
  private
  public :: Cooling_Radius_Simple_Initialize, Cooling_Radius_Simple_Reset

  ! Module global variable that stores the time available for cooling.
  double precision                              :: coolingTimeAvailable
  !$omp threadprivate(coolingTimeAvailable)
  ! Module global pointer to the active node.
  type            (treeNode          ), pointer :: activeNode
  !$omp threadprivate(activeNode)
  ! Internal record of the number of abundance and chemical properties.
  integer                                       :: abundancesCount                      , chemicalsCount

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8    )          :: lastUniqueID                 =-1
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not cooling radius has already been computed for this node.
  logical                                       :: coolingRadiusComputed        =.false., coolingRadiusGrowthRateComputed=.false.
  !$omp threadprivate(coolingRadiusComputed,coolingRadiusGrowthRateComputed)
  ! Stored values of cooling radius.
  double precision                              :: coolingRadiusGrowthRateStored        , coolingRadiusStored
  !$omp threadprivate(coolingRadiusStored,coolingRadiusGrowthRateStored)
  ! Abundances and chemical objects used in cooling calculations.
  type            (abundances        )          :: gasAbundances
  !$omp threadprivate(gasAbundances)
  type            (chemicalAbundances)          :: chemicalDensities                    , chemicalMasses
  !$omp threadprivate(chemicalMasses,chemicalDensities)
  ! Radiation structure used in cooling calculations.
  type            (radiationStructure)          :: radiation
  !$omp threadprivate(radiation)
contains

  !# <coolingRadiusMethod>
  !#  <unitName>Cooling_Radius_Simple_Initialize</unitName>
  !# </coolingRadiusMethod>
  subroutine Cooling_Radius_Simple_Initialize(coolingRadiusMethod,Cooling_Radius_Get,Cooling_Radius_Growth_Rate_Get)
    !% Initializes the ``simple'' cooling radius module.
    use ISO_Varying_String
    use Galacticus_Error
    use Array_Utilities
    use String_Handling
    implicit none
    type     (varying_string                   ), intent(in   )               :: coolingRadiusMethod
    procedure(Cooling_Radius_Simple            ), intent(inout), pointer      :: Cooling_Radius_Get
    procedure(Cooling_Radius_Growth_Rate_Simple), intent(inout), pointer      :: Cooling_Radius_Growth_Rate_Get

    if (coolingRadiusMethod == 'simple') then
       Cooling_Radius_Get             => Cooling_Radius_Simple
       Cooling_Radius_Growth_Rate_Get => Cooling_Radius_Growth_Rate_Simple
       ! Get a count of the number of abundances and chemicals properties.
       abundancesCount=Abundances_Property_Count()
       chemicalsCount =Chemicals_Property_Count ()
       ! Check that required components are gettable.
       if     (                                                                                                                  &
            &  .not.(                                                                                                            &
            &         defaultHotHaloComponent%       massIsGettable() .and.                                                      &
            &         defaultHotHaloComponent% abundancesIsGettable() .and.                                                      &
            &         defaultHotHaloComponent%outerRadiusIsGettable() .and.                                                      &
            &        (defaultHotHaloComponent%  chemicalsIsGettable() .or.  chemicalsCount == 0)                                 &
            &       )                                                                                                            &
            & ) call Galacticus_Error_Report                                                                                     &
            & (                                                                                                                  &
            &  'Cooling_Radius_Simple_Initialize'                                                                              , &
            &  'This method requires that the "mass", "abundances", "outerRadius", and "chemicals" '//                           &
            &  '(if any chemicals are being used) properties of the hot halo are gettable.'         //                           &
            &  Galacticus_Component_List(                                                                                        &
            &                            'hotHalo'                                                                             , &
            &                             defaultHotHaloComponent%massAttributeMatch       (requireGettable=.true.            )  &
            &                            .intersection.                                                                          &
            &                             defaultHotHaloComponent%abundancesAttributeMatch (requireGettable=.true.            )  &
            &                            .intersection.                                                                          &
            &                             defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.            )  &
            &                            .intersection.                                                                          &
            &                             defaultHotHaloComponent%chemicalsAttributeMatch  (requireGettable=chemicalsCount > 0)  &
            &                           )                                                                                        &
            & )
    end if
    return
  end subroutine Cooling_Radius_Simple_Initialize

  !# <calculationResetTask>
  !# <unitName>Cooling_Radius_Simple_Reset</unitName>
  !# </calculationResetTask>
  subroutine Cooling_Radius_Simple_Reset(thisNode)
    !% Reset the cooling radius calculation.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    coolingRadiusComputed          =.false.
    coolingRadiusGrowthRateComputed=.false.
    lastUniqueID                   =thisNode%uniqueID()

    ! Ensure the radiation structure is defined to use the cosmic microwave background.
    if (.not.radiation%isDefined()) call radiation%define([radiationTypeCMB])
    return
  end subroutine Cooling_Radius_Simple_Reset

  double precision function Cooling_Radius_Growth_Rate_Simple(thisNode)
    !% Return the growth rate of the cooling radius in the ``simple'' model in Mpc/Gyr.
    use Hot_Halo_Temperature_Profile
    use Hot_Halo_Density_Profile
    use Cooling_Times
    use Cooling_Times_Available
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    double precision                                               :: coolingRadius                   , coolingTimeAvailable      , &
         &                                                            coolingTimeAvailableIncreaseRate, coolingTimeDensityLogSlope, &
         &                                                            coolingTimeTemperatureLogSlope  , density                   , &
         &                                                            densityLogSlope                 , outerRadius               , &
         &                                                            temperature                     , temperatureLogSlope

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Cooling_Radius_Simple_Reset(thisNode)

    ! Check if cooling radius growth rate is already computed.
    if (.not.coolingRadiusGrowthRateComputed) then
       ! Flag that cooling radius is now computed.
       coolingRadiusGrowthRateComputed=.true.

       ! Get node components.
       thisHotHaloComponent => thisNode%hotHalo()

       ! Get the outer radius.
       outerRadius=thisHotHaloComponent%outerRadius()

       ! Get the cooling radius.
       coolingRadius=Cooling_Radius_Simple(thisNode)

       ! Check if cooling radius has reached the outer radius.
       if (coolingRadius >= outerRadius) then
          coolingRadiusGrowthRateStored=0.0d0
       else
          ! Get the time available for cooling in thisNode.
          coolingTimeAvailable=Cooling_Time_Available(thisNode)

          ! Get the rate of increase of the time available for cooling.
          coolingTimeAvailableIncreaseRate=Cooling_Time_Available_Increase_Rate(thisNode)

          ! Logarithmic slope of density profile.
          densityLogSlope=Hot_Halo_Density_Log_Slope(thisNode,coolingRadius)

          ! Logarithmic slope of density profile.
          temperatureLogSlope=Hot_Halo_Temperature_Logarithmic_Slope(thisNode,coolingRadius)

          ! Get cooling density, temperature and metallicity.
          density=Hot_Halo_Density(activeNode,coolingRadius)
          temperature=Hot_Halo_Temperature(activeNode,coolingRadius)

          ! Logarithmic slope of the cooling time-density relation.
          coolingTimeDensityLogSlope=Cooling_Time_Density_Log_Slope(temperature,density,gasAbundances,chemicalDensities,radiation)

          ! Logarithmic slope of the cooling time-temperature relation.
          coolingTimeTemperatureLogSlope=Cooling_Time_Temperature_Log_Slope(temperature,density,gasAbundances,chemicalDensities,radiation)

          ! Compute rate at which cooling radius grows.
          if (coolingRadius > 0.0d0) then
             coolingRadiusGrowthRateStored=(coolingRadius/coolingTimeAvailable)*coolingTimeAvailableIncreaseRate&
                  &/(densityLogSlope*coolingTimeDensityLogSlope+temperatureLogSlope*coolingTimeTemperatureLogSlope)
          else
             coolingRadiusGrowthRateStored=0.0d0
          end if
       end if
    end if
    ! Return the stored value.
    Cooling_Radius_Growth_Rate_Simple=coolingRadiusGrowthRateStored
    return
  end function Cooling_Radius_Growth_Rate_Simple

  double precision function Cooling_Radius_Simple(thisNode)
    !% Return the cooling radius in the simple model.
    use Cooling_Times_Available
    use Root_Finder
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    double precision                      , parameter              :: zeroRadius          =0.0d0
    double precision                      , parameter              :: toleranceAbsolute   =0.0d0, toleranceRelative=1.0d-6
    type            (rootFinder          ), save                   :: finder
    !$omp threadprivate(finder)
    double precision                                               :: outerRadius

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Cooling_Radius_Simple_Reset(thisNode)

    ! Check if cooling radius is already computed.
    if (.not.coolingRadiusComputed) then
       ! Flag that cooling radius is now computed.
       coolingRadiusComputed=.true.

       ! Get the time available for cooling in thisNode.
       coolingTimeAvailable=Cooling_Time_Available(thisNode)

       ! Make the node available to the root finding routine.
       activeNode => thisNode

       ! Initialize quantities needed by the solver.
       call Cooling_Radius_Solver_Initialize(thisNode)

       ! Check if cooling time at halo center is reached.
       if (Cooling_Radius_Root(zeroRadius) > 0.0d0) then
          ! Cooling time at halo center exceeds the time available, return zero radius.
          coolingRadiusStored  =zeroRadius
          Cooling_Radius_Simple=coolingRadiusStored
         return
       end if

       ! Get node components.
       thisHotHaloComponent => thisNode%hotHalo()

       ! Check if cooling time at hot halo outer radius is reached.
       outerRadius=thisHotHaloComponent%outerRadius()
       if (Cooling_Radius_Root(outerRadius) < 0.0d0) then
          ! Cooling time available exceeds cooling time at outer radius radius, return outer radius.
          coolingRadiusStored=outerRadius
          Cooling_Radius_Simple=coolingRadiusStored
          return
       end if

       ! Cooling radius is between zero and outer radii. Search for the cooling radius.
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(Cooling_Radius_Root                )
          call finder%tolerance   (toleranceAbsolute,toleranceRelative)
       end if
       coolingRadiusStored=finder%find(rootRange=[zeroRadius,outerRadius])
       Cooling_Radius_Simple=coolingRadiusStored
    else
       Cooling_Radius_Simple=coolingRadiusStored
    end if
    return
  end function Cooling_Radius_Simple

  subroutine Cooling_Radius_Solver_Initialize(thisNode)
    !% Initialize the abundances, chemical properties and radiation field for {\tt thisNode} for use in cooling radius
    !% calculations.
    use Chemical_Reaction_Rates_Utilities
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    double precision                                               :: massToDensityConversion

    ! Get node components.
    thisHotHaloComponent => thisNode%hotHalo()

    ! Get the abundances for this node.
    gasAbundances=thisHotHaloComponent%abundances()
    call gasAbundances%massToMassFraction(thisHotHaloComponent%mass())

    ! Get the chemicals for this node.
    if (chemicalsCount > 0) then
       chemicalMasses=thisHotHaloComponent%chemicals()
       ! Scale all chemical masses by their mass in atomic mass units to get a number density.
       call chemicalMasses%massToNumber(chemicalDensities)
       ! Compute factor converting mass of chemicals in (M_Solar/M_Atomic) to number density in cm^-3.
       if (thisHotHaloComponent%outerRadius() > 0.0d0) then
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(thisHotHaloComponent%outerRadius())
       else
          massToDensityConversion=0.0d0
       end if
       ! Convert to number density.
       chemicalDensities=chemicalDensities*massToDensityConversion
    end if

    ! Set the radiation field.
    call radiation%set(thisNode)
    return
  end subroutine Cooling_Radius_Solver_Initialize

  double precision function Cooling_Radius_Root(radius)
    !% Root function which evaluates the difference between the cooling time at {\tt radius} and the time available for cooling.
    use Cooling_Times
    use Hot_Halo_Density_Profile
    use Hot_Halo_Temperature_Profile
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: coolingTime, density, temperature

    ! Compute density, temperature and abundances.
    density    =Hot_Halo_Density    (activeNode,radius)
    temperature=Hot_Halo_Temperature(activeNode,radius)
    ! Compute the cooling time at the specified radius.
    coolingTime=Cooling_Time(temperature,density,gasAbundances,chemicalDensities,radiation)
    ! Return the difference between cooling time and time available.
    Cooling_Radius_Root=coolingTime-coolingTimeAvailable
    return
  end function Cooling_Radius_Root

end module Cooling_Radii_Simple
