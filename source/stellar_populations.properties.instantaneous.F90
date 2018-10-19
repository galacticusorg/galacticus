!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements stellar population properties in the instantaneous recycling approximation.

module Stellar_Population_Properties_Instantaneous
  !% Implements stellar population properties in the instantaneous recycling approximation.
  implicit none
  private
  public :: Stellar_Population_Properties_Instantaneous_Initialize

  ! Index of abundance pattern to use for elemental abundances.
  integer :: abundanceIndex

contains

  !# <stellarPopulationPropertiesMethod>
  !#  <unitName>Stellar_Population_Properties_Instantaneous_Initialize</unitName>
  !# </stellarPopulationPropertiesMethod>
  subroutine Stellar_Population_Properties_Instantaneous_Initialize(stellarPopulationPropertiesMethod &
       &,Stellar_Population_Properties_Rates_Get,Stellar_Population_Properties_Scales_Get&
       &,Stellar_Population_Properties_History_Count_Get ,Stellar_Population_Properties_History_Create_Do)
    !% Initializes the instantaneous recycling approximation stellar population properties module.
    use ISO_Varying_String
    use Atomic_Data
    implicit none
    type     (varying_string                                            ), intent(in   )          :: stellarPopulationPropertiesMethod
    procedure(Stellar_Population_Properties_History_Count_Instantaneous ), intent(inout), pointer :: Stellar_Population_Properties_History_Count_Get
    procedure(Stellar_Population_Properties_History_Create_Instantaneous), intent(inout), pointer :: Stellar_Population_Properties_History_Create_Do
    procedure(Stellar_Population_Properties_Rates_Instantaneous         ), intent(inout), pointer :: Stellar_Population_Properties_Rates_Get
    procedure(Stellar_Population_Properties_Scales_Instantaneous        ), intent(inout), pointer :: Stellar_Population_Properties_Scales_Get

    if (stellarPopulationPropertiesMethod == 'instantaneous') then
       Stellar_Population_Properties_Rates_Get         => Stellar_Population_Properties_Rates_Instantaneous
       Stellar_Population_Properties_Scales_Get        => Stellar_Population_Properties_Scales_Instantaneous
       Stellar_Population_Properties_History_Count_Get => Stellar_Population_Properties_History_Count_Instantaneous
       Stellar_Population_Properties_History_Create_Do => Stellar_Population_Properties_History_Create_Instantaneous

       ! Get index of abundance pattern to use.
       abundanceIndex=Abundance_Pattern_Lookup(abundanceName="solar")
    end if
    return
  end subroutine Stellar_Population_Properties_Instantaneous_Initialize

  integer function Stellar_Population_Properties_History_Count_Instantaneous()
    !% Returns the number of histories required by the instantaneous stellar populations properties module.
    implicit none

    ! We require no histories.
    Stellar_Population_Properties_History_Count_Instantaneous=0
    return
  end function Stellar_Population_Properties_History_Count_Instantaneous

  subroutine Stellar_Population_Properties_Rates_Instantaneous(starFormationRate,fuelAbundances,component_,node,propertiesHistory,stellarMassRate&
       &,stellarAbundancesRates,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate,stellarLuminositiesRatesCompute)
    !% Return an array of stellar population property rates of change given a star formation rate and fuel abundances.
    use, intrinsic :: ISO_C_Binding
    use            :: Galactic_Structure_Options
    use            :: Galacticus_Nodes
    use            :: Abundances_Structure
    use            :: Histories
    use            :: Stellar_Luminosities_Structure
    use            :: Stellar_Feedback
    use            :: Stellar_Populations
    use            :: Stellar_Population_Selectors
    implicit none
    double precision                                , intent(  out)      :: energyInputRate                              , fuelMassRate                         , &
         &                                                                  stellarMassRate
    type            (abundances                    ), intent(inout)      :: fuelAbundancesRates                          , stellarAbundancesRates
    type            (stellarLuminosities           ), intent(inout)      :: stellarLuminositiesRates
    double precision                                , intent(in   )      :: starFormationRate
    type            (abundances                    ), intent(in   )      :: fuelAbundances
    class           (nodeComponent                 ), intent(in   )      :: component_
    type            (treeNode                      ), intent(inout)      :: node
    type            (history                       ), intent(inout)      :: propertiesHistory
    logical                                         , intent(in   )      :: stellarLuminositiesRatesCompute
    class           (nodeComponentBasic            ), pointer            :: basic
    class           (stellarPopulationClass        ), pointer            :: stellarPopulation_
    class           (stellarPopulationSelectorClass), pointer            :: stellarPopulationSelector_
    type            (stellarLuminosities           ), dimension(2), save :: stellarLuminositiesRatesPrevious
    double precision                                , dimension(2), save :: starFormationRatePrevious       =-huge(0.0d0), fuelMetallicityPrevious=-huge(0.0d0), &
         &                                                                  timePrevious                    =-huge(0.0d0)
    integer         (c_size_t                      ), dimension(2), save :: populationIDPrevious            =-1_c_size_t
    !$omp threadprivate(stellarLuminositiesRatesPrevious,starFormationRatePrevious,fuelMetallicityPrevious,timePrevious,populationIDPrevious)
    integer         (c_size_t                      )                     :: populationID
    integer                                                              :: componentIndex
    double precision                                                     :: fuelMetallicity                              , fuelMetalsRateOfChange              , &
         &                                                                  recycledFractionInstantaneous                , stellarMetalsRateOfChange           , &
         &                                                                  time                                         , yieldInstantaneous
    !GCC$ attributes unused :: propertiesHistory
    
    ! Get the instantaneous recycled fraction and yields for the selected stellar population.
    stellarPopulationSelector_    => stellarPopulationSelector                               (                                           )
    stellarPopulation_            => stellarPopulationSelector_%select                       (starFormationRate,fuelAbundances,component_)
    populationID                  =  stellarPopulation_        %uniqueID                     (                                           )
    recycledFractionInstantaneous =  stellarPopulation_        %recycledFractionInstantaneous(                                           )
    yieldInstantaneous            =  stellarPopulation_        %yieldInstantaneous           (                                           )
    ! Get the metallicity of the fuel supply.
    fuelMetallicity=Abundances_Get_Metallicity(fuelAbundances)
    ! Set the stellar and fuel mass rates of change.
    stellarMassRate= (1.0d0-recycledFractionInstantaneous)*starFormationRate
    fuelMassRate   =-stellarMassRate
    ! Set energy input rate to canonical value assuming that all energy is injected instantaneously.
    energyInputRate=starFormationRate*feedbackEnergyInputAtInfinityCanonical
    ! Set the rates of change of the stellar and fuel metallicities.
    stellarMetalsRateOfChange=stellarMassRate*fuelMetallicity
    fuelMetalsRateOfChange   =-stellarMetalsRateOfChange+yieldInstantaneous*starFormationRate
    call stellarAbundancesRates%metallicitySet(stellarMetalsRateOfChange,adjustElements=adjustElementsReset,abundanceIndex=abundanceIndex)
    call fuelAbundancesRates   %metallicitySet(fuelMetalsRateOfChange   ,adjustElements=adjustElementsReset,abundanceIndex=abundanceIndex)
    ! Get the current cosmological time for this node.
    basic => node %basic()
    time  =  basic%time ()
    ! Set luminosity rates of change.
    if (stellarLuminositiesRatesCompute) then
       select type (component_)
       class is (nodeComponentDisk    )
          componentIndex=componentTypeDisk
       class is (nodeComponentSpheroid)
          componentIndex=componentTypeSpheroid
       class default
          componentIndex=-1
       end select
       select case (componentIndex)
       case (componentTypeDisk,componentTypeSpheroid)
          if     (                                                                &
               &   populationID       /=     populationIDPrevious(componentIndex) &
               &  .or.                                                            &
               &   starFormationRate /= starFormationRatePrevious(componentIndex) &
               &  .or.                                                            &
               &   time              /=              timePrevious(componentIndex) &
               &  .or.                                                            &
               &   fuelMetallicity   /=   fuelMetallicityPrevious(componentIndex) &
               & ) then
             call stellarLuminositiesRatesPrevious(componentIndex)%setLuminosities(starFormationRate,stellarPopulation_,time,fuelAbundances)
             populationIDPrevious     (componentIndex)=populationID
             starFormationRatePrevious(componentIndex)=starFormationRate
             timePrevious             (componentIndex)=time
             fuelMetallicityPrevious  (componentIndex)=fuelMetallicity
          end if
          stellarLuminositiesRates=stellarLuminositiesRatesPrevious(componentIndex)
       case default
          call stellarLuminositiesRates%setLuminosities(starFormationRate,stellarPopulation_,time,fuelAbundances)
       end select
    end if
    return
  end subroutine Stellar_Population_Properties_Rates_Instantaneous

  subroutine Stellar_Population_Properties_Scales_Instantaneous(propertiesHistory,stellarMass,stellarAbundances)
    !% Set the scalings for error control on the absolute values of stellar population properties. The instantaneous method
    !% requires none, so just return.
    use Histories
    use Abundances_Structure
    implicit none
    double precision            , intent(in   ) :: stellarMass
    type            (abundances), intent(in   ) :: stellarAbundances
    type            (history   ), intent(inout) :: propertiesHistory
    !GCC$ attributes unused :: propertiesHistory, stellarMass, stellarAbundances

    ! No history is used in this case, so simply return.
    return
  end subroutine Stellar_Population_Properties_Scales_Instantaneous

  subroutine Stellar_Population_Properties_History_Create_Instantaneous(node,propertiesHistory)
    !% Create any history required for storing stellar population properties. The instantaneous method requires none, so don't
    !% create one.
    use Galacticus_Nodes
    use Histories
    implicit none
    type(treeNode), intent(inout) :: node
    type(history ), intent(inout) :: propertiesHistory
    !GCC$ attributes unused :: node, propertiesHistory
    
    return
  end subroutine Stellar_Population_Properties_History_Create_Instantaneous

end module Stellar_Population_Properties_Instantaneous
