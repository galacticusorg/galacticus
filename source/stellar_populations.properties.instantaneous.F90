!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
    type(varying_string),          intent(in)    :: stellarPopulationPropertiesMethod
    procedure(integer),   pointer, intent(inout) :: Stellar_Population_Properties_History_Count_Get
    procedure(),          pointer, intent(inout) :: Stellar_Population_Properties_History_Create_Do&
         &,Stellar_Population_Properties_Rates_Get,Stellar_Population_Properties_Scales_Get

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
  
  subroutine Stellar_Population_Properties_Rates_Instantaneous(starFormationRate,fuelAbundances,component,thisNode,thisHistory,stellarMassRate&
       &,stellarAbundancesRates,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
    !% Return an array of stellar population property rates of change given a star formation rate and fuel abundances.
    use Tree_Nodes
    use Abundances_Structure
    use Histories
    use Star_Formation_IMF
    use Stellar_Population_Properties_Luminosities
    use Stellar_Feedback
    implicit none
    double precision,          intent(out)                 :: stellarMassRate,fuelMassRate,energyInputRate
    type(abundancesStructure), intent(inout)               :: stellarAbundancesRates,fuelAbundancesRates
    double precision,          intent(out),   dimension(:) :: stellarLuminositiesRates
    double precision,          intent(in)                  :: starFormationRate
    type(abundancesStructure), intent(in)                  :: fuelAbundances
    integer,                   intent(in)                  :: component
    type(treeNode),            intent(inout), pointer      :: thisNode
    type(history),             intent(inout)               :: thisHistory
    integer                                                :: imfSelected
    double precision                                       :: recycledFractionInstantaneous,yieldInstantaneous,fuelMetallicity&
         &,stellarMetalsRateOfChange,fuelMetalsRateOfChange,time
 
    ! Get the instantaneous recycling rate for the IMF.
    recycledFractionInstantaneous=IMF_Recycled_Fraction_Instantaneous(starFormationRate,fuelAbundances,component)

    ! Get the yield for this IMF.
    yieldInstantaneous=IMF_Yield_Instantaneous(starFormationRate,fuelAbundances,component)

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

    ! Get the IMF.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances,component)

    ! Get the current cosmological time for this node.
    time=Tree_Node_Time(thisNode)

    ! Set luminosity rates of change.
    if (size(stellarLuminositiesRates) > 0) stellarLuminositiesRates=starFormationRate&
         &*Stellar_Population_Luminosities_Get(imfSelected,time,fuelAbundances)

    return
  end subroutine Stellar_Population_Properties_Rates_Instantaneous

  subroutine Stellar_Population_Properties_Scales_Instantaneous(thisHistory,stellarMass,stellarAbundances)
    !% Set the scalings for error control on the absolute values of stellar population properties. The instantaneous method
    !% requires none, so just return.
    use Histories
    use Abundances_Structure
    implicit none
    double precision,          intent(in)    :: stellarMass
    type(abundancesStructure), intent(in)    :: stellarAbundances
    type(history),             intent(inout) :: thisHistory

    ! No history is used in this case, so simply return.

    return
  end subroutine Stellar_Population_Properties_Scales_Instantaneous

  subroutine Stellar_Population_Properties_History_Create_Instantaneous(thisNode,thisHistory)
    !% Create any history required for storing stellar population properties. The instantaneous method requires none, so just
    !% return.
    use Tree_Nodes
    use Histories
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(history),  intent(inout)          :: thisHistory
  
    return
  end subroutine Stellar_Population_Properties_History_Create_Instantaneous

end module Stellar_Population_Properties_Instantaneous
