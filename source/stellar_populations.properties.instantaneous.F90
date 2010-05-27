!% Contains a module which implements stellar population properties in the instantaneous recycling approximation.

module Stellar_Population_Properties_Instantaneous
  !% Implements stellar population properties in the instantaneous recycling approximation.
  private
  public :: Stellar_Population_Properties_Instantaneous_Initialize

contains

  !# <stellarPopulationPropertiesMethod>
  !#  <unitName>Stellar_Population_Properties_Instantaneous_Initialize</unitName>
  !# </stellarPopulationPropertiesMethod>
  subroutine Stellar_Population_Properties_Instantaneous_Initialize(stellarPopulationPropertiesMethod &
       &,Stellar_Population_Properties_Rates_Get,Stellar_Population_Properties_History_Count_Get&
       &,Stellar_Population_Properties_History_Create_Do)
    !% Initializes the instantaneous recycling approximation stellar population properties module.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: stellarPopulationPropertiesMethod
    procedure(),          pointer, intent(inout) :: Stellar_Population_Properties_Rates_Get&
         &,Stellar_Population_Properties_History_Count_Get,Stellar_Population_Properties_History_Create_Do
    
    if (stellarPopulationPropertiesMethod == 'instantaneous') then
       Stellar_Population_Properties_Rates_Get         =>  Stellar_Population_Properties_Rates_Instantaneous    
       Stellar_Population_Properties_History_Count_Get =>  Stellar_Population_Properties_History_Count_Instantaneous    
       Stellar_Population_Properties_History_Create_Do => Stellar_Population_Properties_History_Create_Instantaneous
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
  
  subroutine Stellar_Population_Properties_Rates_Instantaneous(starFormationRate,fuelAbundances,thisNode,thisHistory,stellarMassRate&
       &,stellarAbundancesRates,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
    !% Return an array of stellar population property rates of change given a star formation rate and fuel abundances.
    use Tree_Nodes
    use Tree_Node_Methods
    use Abundances_Structure
    use Histories
    use Star_Formation_IMF
    use Stellar_Population_Properties_Luminosities
    use Stellar_Feedback
    implicit none
    double precision,          intent(out)                 :: stellarMassRate,fuelMassRate,energyInputRate
    type(abundancesStructure), intent(out)                 :: stellarAbundancesRates,fuelAbundancesRates
    double precision,          intent(out),   dimension(:) :: stellarLuminositiesRates
    double precision,          intent(in)                  :: starFormationRate
    type(abundancesStructure), intent(in)                  :: fuelAbundances
    type(treeNode),            intent(inout), pointer      :: thisNode
    type(history),             intent(inout)               :: thisHistory
    integer                                                :: imfSelected
    double precision                                       :: recycledFractionInstantaneous,yieldInstantaneous,fuelMetallicity&
         &,stellarMetalsRateOfChange,fuelMetalsRateOfChange,time
 
    ! Get the instantaneous recycling rate for the IMF.
    recycledFractionInstantaneous=IMF_Recycled_Fraction_Instantaneous(starFormationRate,fuelAbundances)

    ! Get the yield for this IMF.
    yieldInstantaneous=IMF_Yield_Instantaneous(starFormationRate,fuelAbundances)

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
    call stellarAbundancesRates%metallicitySet(stellarMetalsRateOfChange)
    call fuelAbundancesRates   %metallicitySet(fuelMetalsRateOfChange   )

    ! Get the IMF.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances)

    ! Get the current cosmological time for this node.
    time=Tree_Node_Time(thisNode)

    ! Set luminosity rates of change.
    if (size(stellarLuminositiesRates) > 0) stellarLuminositiesRates=starFormationRate&
         &*Stellar_Population_Luminosities_Get(imfSelected,time,fuelAbundances)

    return
  end subroutine Stellar_Population_Properties_Rates_Instantaneous

  subroutine Stellar_Population_Properties_History_Create_Instantaneous(thisNode,thisHistory)
    !% Create any history required for storing stellar population properties. The instantaneous method requires none, so just return.
    use Tree_Nodes
    use Histories
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(history),  intent(inout) :: thisHistory
  
    return
  end subroutine Stellar_Population_Properties_History_Create_Instantaneous

end module Stellar_Population_Properties_Instantaneous
