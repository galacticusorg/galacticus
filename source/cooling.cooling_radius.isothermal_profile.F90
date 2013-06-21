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

!% Contains a module which implements a calculation of cooling radius appropriate for an
!% isothermal halo, assuming collisional ionization equilibrium such that cooling time
!% scales as inverse density.

module Cooling_Radii_Isothermal
  !% Implements calculation of cooling radius appropriate for an isothermal halo, assuming
  !% collisional ionization equilibrium such that cooling time scales as inverse density.
  use Galacticus_Nodes
  use Kind_Numbers
  use Radiation_Structure
  use Abundances_Structure
  use Chemical_Abundances_Structure
  implicit none
  private
  public :: Cooling_Radius_Isothermal_Initialize, Cooling_Radius_Isothermal_Reset

  ! Internal record of the number of abundance and chemical properties.
  integer                              :: abundancesCount                      , chemicalsCount                          
  
  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8    ) :: lastUniqueID                 =-1                                               
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not cooling radius has already been computed for this node.
  logical                              :: coolingRadiusComputed        =.false., coolingRadiusGrowthRateComputed=.false. 
  !$omp threadprivate(coolingRadiusComputed,coolingRadiusGrowthRateComputed)
  ! Stored values of cooling radius.
  double precision                     :: coolingRadiusGrowthRateStored        , coolingRadiusStored                     
  !$omp threadprivate(coolingRadiusStored,coolingRadiusGrowthRateStored)
  ! Abundances and chemical objects used in cooling calculations.
  type            (abundances        ) :: hotAbundances                                                                  
  !$omp threadprivate(hotAbundances)
  type            (chemicalAbundances) :: chemicalDensities                    , chemicalMasses                          
  !$omp threadprivate(chemicalMasses,chemicalDensities)
  ! Radiation structure used in cooling calculations.
  type            (radiationStructure) :: radiation                                                                      
  !$omp threadprivate(radiation)
contains

  !# <coolingRadiusMethod>
  !#  <unitName>Cooling_Radius_Isothermal_Initialize</unitName>
  !# </coolingRadiusMethod>
  subroutine Cooling_Radius_Isothermal_Initialize(coolingRadiusMethod,Cooling_Radius_Get,Cooling_Radius_Growth_Rate_Get)
    !% Initializes the ``isothermal'' cooling radius module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                       ), intent(in   )          :: coolingRadiusMethod            
    procedure(Cooling_Radius_Isothermal            ), intent(inout), pointer :: Cooling_Radius_Get             
    procedure(Cooling_Radius_Growth_Rate_Isothermal), intent(inout), pointer :: Cooling_Radius_Growth_Rate_Get 
    
    if (coolingRadiusMethod == 'isothermal') then
       Cooling_Radius_Get             => Cooling_Radius_Isothermal
       Cooling_Radius_Growth_Rate_Get => Cooling_Radius_Growth_Rate_Isothermal
       ! Get a count of the number of abundances and chemicals properties.
       abundancesCount=Abundances_Property_Count()
       chemicalsCount =Chemicals_Property_Count ()
    end if
    return
  end subroutine Cooling_Radius_Isothermal_Initialize

  !# <calculationResetTask>
  !# <unitName>Cooling_Radius_Isothermal_Reset</unitName>
  !# </calculationResetTask>
  subroutine Cooling_Radius_Isothermal_Reset(thisNode)
    !% Reset the cooling radius calculation.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode 
    
    coolingRadiusComputed          =.false.
    coolingRadiusGrowthRateComputed=.false.
    lastUniqueID                   =thisNode%uniqueID()

    ! Ensure the radiation structure is defined to be null.
    if (.not.radiation%isDefined()) call radiation%define([radiationTypeNull])
    return
  end subroutine Cooling_Radius_Isothermal_Reset

  double precision function Cooling_Radius_Growth_Rate_Isothermal(thisNode)
    !% Return the growth rate of the cooling radius in the ``isothermal'' model in Mpc/Gyr.
    use Dark_Matter_Halo_Scales
    use Cooling_Times_Available
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode                                                  
    double precision                                   :: coolingRadius                   , coolingTimeAvailable, & 
         &                                                coolingTimeAvailableIncreaseRate, virialRadius            
    
    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Cooling_Radius_Isothermal_Reset(thisNode)

    ! Check if cooling radius growth rate is already computed.
    if (.not.coolingRadiusGrowthRateComputed) then
       ! Flag that cooling radius is now computed.
       coolingRadiusGrowthRateComputed=.true.

       ! Get the cooling radius.
       coolingRadius=Cooling_Radius_Isothermal(thisNode)
       
       ! Get the virial radius.
       virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
       
       ! Check if cooling radius has reached virial radius.
       if (coolingRadius >= virialRadius) then
          coolingRadiusGrowthRateStored=0.0d0
       else
          ! Get the time available for cooling in thisNode.
          coolingTimeAvailable=Cooling_Time_Available(thisNode)
          
          ! Get the rate of increase of the time available for cooling.
          coolingTimeAvailableIncreaseRate=Cooling_Time_Available_Increase_Rate(thisNode)

          ! Compute the growth rate of the cooling radius.
          coolingRadiusGrowthRateStored=0.5d0*coolingRadius*coolingTimeAvailableIncreaseRate/coolingTimeAvailable

       end if
    end if
    ! Return the stored value.
    Cooling_Radius_Growth_Rate_Isothermal=coolingRadiusGrowthRateStored       
    return
  end function Cooling_Radius_Growth_Rate_Isothermal

  double precision function Cooling_Radius_Isothermal(thisNode)
    !% Return the cooling radius in the isothermal model.
    use Dark_Matter_Halo_Scales
    use Cooling_Times_Available
    use Chemical_Reaction_Rates_Utilities
    use Cooling_Times
    use Hot_Halo_Density_Profile
    use Hot_Halo_Temperature_Profile
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode                                         
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent                             
    double precision                                               :: coolingTime         , coolingTimeAvailable   , & 
         &                                                            density             , massToDensityConversion, & 
         &                                                            temperature         , virialRadius               
    
    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Cooling_Radius_Isothermal_Reset(thisNode)

    ! Check if cooling radius is already computed.
    if (.not.coolingRadiusComputed) then
       ! Flag that cooling radius is now computed.
       coolingRadiusComputed=.true.

       ! Get the time available for cooling in thisNode.
       coolingTimeAvailable=Cooling_Time_Available(thisNode)

       ! Get the hot halo component.
       thisHotHaloComponent => thisNode%hotHalo()

       ! Get the abundances for this node.
       hotAbundances=thisHotHaloComponent%abundances()
       call hotAbundances%massToMassFraction(thisHotHaloComponent%mass())
       
       ! Get the chemicals for this node.
       if (chemicalsCount > 0) then
          chemicalMasses=thisHotHaloComponent%chemicals()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(chemicalDensities)
          ! Compute factor converting mass of chemicals in (M_Solar/M_Atomic) to number density in cm^-3.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(Dark_Matter_Halo_Virial_Radius(thisNode))
          ! Convert to number density.
          chemicalDensities=chemicalDensities*massToDensityConversion
       end if
       
       ! Set the radiation field.
       call radiation%set(thisNode)

       ! Get the virial radius.
       virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode             )
       
       ! Compute density, temperature and abundances.
       density     =Hot_Halo_Density              (thisNode,virialRadius)
       temperature =Hot_Halo_Temperature          (thisNode,virialRadius)
       
       ! Compute the cooling time at the virial radius.
       coolingTime=Cooling_Time(temperature,density,hotAbundances,chemicalDensities,radiation)

       if (coolingTime < coolingTimeAvailable) then
          ! Cooling time available exceeds cooling time at virial radius, return virial radius.
          coolingRadiusStored=virialRadius
       else       
          ! Cooling radius is between zero and virial radii.
          coolingRadiusStored=virialRadius*sqrt(coolingTimeAvailable/coolingTime)
       end if
    end if
    Cooling_Radius_Isothermal=coolingRadiusStored
    return
  end function Cooling_Radius_Isothermal

end module Cooling_Radii_Isothermal
