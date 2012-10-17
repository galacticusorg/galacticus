!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements stellar population properties with noninstantaneous recycling.

module Stellar_Population_Properties_Noninstantaneous
  !% Implements stellar population properties with noninstantaneous recycling.
  implicit none
  private
  public :: Stellar_Population_Properties_Noninstantaneous_Initialize

  ! Count of number of elements (plus total metals) that are to be tracked.
  integer            :: elementsCount
  ! Count of the number of histories required by this implementation.
  integer            :: historyCount
  ! Indices for histories.
  integer, parameter :: recycledRateIndex          =1
  integer, parameter :: energyInputRateIndex       =2
  integer, parameter :: returnedMetalRateBeginIndex=3
  integer            :: returnedMetalRateEndIndex,metalYieldRateBeginIndex,metalYieldRateEndIndex

  ! Number of times to store in histories.
  integer            :: noninstantHistoryTimesCount

contains

  !# <stellarPopulationPropertiesMethod>
  !#  <unitName>Stellar_Population_Properties_Noninstantaneous_Initialize</unitName>
  !# </stellarPopulationPropertiesMethod>
  subroutine Stellar_Population_Properties_Noninstantaneous_Initialize(stellarPopulationPropertiesMethod &
       &,Stellar_Population_Properties_Rates_Get,Stellar_Population_Properties_Scales_Get&
       &,Stellar_Population_Properties_History_Count_Get ,Stellar_Population_Properties_History_Create_Do)
    !% Initializes the noninstantaneous recycling stellar population properties module.
    use ISO_Varying_String
    use Input_Parameters
    use Abundances_Structure
    implicit none
    type(varying_string),          intent(in)    :: stellarPopulationPropertiesMethod
    procedure(integer),   pointer, intent(inout) :: Stellar_Population_Properties_History_Count_Get
    procedure(),          pointer, intent(inout) :: Stellar_Population_Properties_Rates_Get&
         &,Stellar_Population_Properties_History_Create_Do,Stellar_Population_Properties_Scales_Get
    
    if (stellarPopulationPropertiesMethod == 'noninstantaneous') then
       Stellar_Population_Properties_Rates_Get         => Stellar_Population_Properties_Rates_Noninstantaneous    
       Stellar_Population_Properties_Scales_Get        => Stellar_Population_Properties_Scales_Noninstantaneous    
       Stellar_Population_Properties_History_Count_Get => Stellar_Population_Properties_History_Count_Noninstantaneous    
       Stellar_Population_Properties_History_Create_Do => Stellar_Population_Properties_History_Create_Noninstantaneous
       !@ <inputParameter>
       !@   <name>noninstantHistoryTimesCount</name>
       !@   <defaultValue>10</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The number of times at which a galaxy's stellar properties history is stored.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('noninstantHistoryTimesCount',noninstantHistoryTimesCount,defaultValue=10)

       ! Get a count of the number of elements (plus total metals) that will be tracked.
       elementsCount=Abundances_Property_Count()

       ! Determine the number of histories that we must store. We require one for recycled mass, one for energy input and two
       ! (return rate and yield rate) for each element.
       historyCount=2+2*elementsCount

       ! Establish indices in the array of histories where returned metal rates and yield rates will be stored.
       returnedMetalRateEndIndex=returnedMetalRateBeginIndex+elementsCount-1
       metalYieldRateBeginIndex =returnedMetalRateEndIndex                +1
       metalYieldRateEndIndex   =metalYieldRateBeginIndex   +elementsCount-1
      end if
    return
  end subroutine Stellar_Population_Properties_Noninstantaneous_Initialize

  integer function Stellar_Population_Properties_History_Count_Noninstantaneous()
    !% Returns the number of histories required by the noninstantaneous stellar populations properties module.
    implicit none
  
    ! Return number of histories required.
    Stellar_Population_Properties_History_Count_Noninstantaneous=historyCount
    return
  end function Stellar_Population_Properties_History_Count_Noninstantaneous
  
  subroutine Stellar_Population_Properties_Rates_Noninstantaneous(starFormationRate,fuelAbundances,component,thisNode,thisHistory&
       &,stellarMassRate ,stellarAbundancesRates,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
    !% Return an array of stellar population property rates of change given a star formation rate and fuel abundances.
    use Tree_Nodes
    use Abundances_Structure
    use Histories
    use Star_Formation_IMF
    use Stellar_Population_Properties_Luminosities
    use Numerical_Interpolation
    use FGSL
    implicit none
    double precision,          intent(out)                             :: stellarMassRate,fuelMassRate,energyInputRate
    type(abundancesStructure), intent(inout)                           :: stellarAbundancesRates,fuelAbundancesRates
    double precision,          intent(out),   dimension(:)             :: stellarLuminositiesRates
    double precision,          intent(in)                              :: starFormationRate
    type(abundancesStructure), intent(in)                              :: fuelAbundances
    integer,                   intent(in)                              :: component
    type(treeNode),            intent(inout), pointer                  :: thisNode
    type(history),             intent(inout)                           :: thisHistory
    double precision,                         dimension(elementsCount) :: fuelMetallicity,stellarMetalsRateOfChange&
         &,fuelMetalsRateOfChange,metalReturnRate,metalYieldRate
    integer                                                            :: imfSelected,iHistory,iElement
    double precision                                                   :: ageMinimum,ageMaximum,currentTime,recyclingRate
    type(fgsl_interp_accel)                                            :: interpolationAccelerator
    logical                                                            :: interpolationReset

    ! Get the current time.
    currentTime=Tree_Node_Time(thisNode)

    ! Get interpolating factors in stellar population history.
    interpolationReset=.true.
    iHistory          =Interpolate_Locate(size(thisHistory%time),thisHistory%time,interpolationAccelerator,currentTime,interpolationReset)
    call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)

    ! Get recycling, energy input, metal recycling and metal yield rates.
    recyclingRate  =thisHistory%data(iHistory,          recycledRateIndex                                            )
    energyInputRate=thisHistory%data(iHistory,       energyInputRateIndex                                            )
    metalReturnRate=thisHistory%data(iHistory,returnedMetalRateBeginIndex:returnedMetalRateBeginIndex+elementsCount-1)
    metalYieldRate =thisHistory%data(iHistory,   metalYieldRateBeginIndex:   metalYieldRateBeginIndex+elementsCount-1)

    ! Get the metallicity of the fuel supply.
    call fuelAbundances%unpack(fuelMetallicity)

    ! Set the stellar and fuel mass rates of change.
    stellarMassRate=starFormationRate-recyclingRate
    fuelMassRate   =-stellarMassRate

    ! Set the rates of change of the stellar and fuel metallicities.
    stellarMetalsRateOfChange=starFormationRate*fuelMetallicity-metalReturnRate
    fuelMetalsRateOfChange   =-stellarMetalsRateOfChange+metalYieldRate
    call stellarAbundancesRates%pack(stellarMetalsRateOfChange)
    call fuelAbundancesRates   %pack(fuelMetalsRateOfChange   )

    ! Get the IMF.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances,component)

    ! Set luminosity rates of change.
    if (size(stellarLuminositiesRates) > 0) stellarLuminositiesRates=starFormationRate&
         &*Stellar_Population_Luminosities_Get(imfSelected,currentTime,fuelAbundances)

    ! Set rates of change in the stellar populations properties future history.
    do iHistory=1,size(thisHistory%time)-1
       ! Find the age of the forming stellar population at the future time. We average over the time between successive timesteps
       ! to ensure that the rates will integrate to the correct values.
       ageMinimum=max(thisHistory%time(iHistory  )-currentTime,0.0d0)
       ageMaximum=max(thisHistory%time(iHistory+1)-currentTime,0.0d0)
       ! Check that it really is in the future and that the timestep over which the contribution to be made is non-zero.
       if (ageMaximum >= 0.0d0 .and. ageMaximum > ageMinimum) then
          ! Get the recycling rate.
          recyclingRate=IMF_Recycling_Rate_NonInstantaneous(starFormationRate,fuelAbundances,component,ageMinimum,ageMaximum)&
               &*starFormationRate
          ! Accumulate the mass recycling rate from this population at the future time.
          thisHistory%rates(iHistory,recycledRateIndex     )=thisHistory%rates(iHistory,recycledRateIndex     ) +recyclingRate
          ! Get the (normalized) energy input rate.
          thisHistory%rates(iHistory,energyInputRateIndex  )=thisHistory%rates(iHistory,energyInputRateIndex  ) &
               &+IMF_Energy_Input_Rate_NonInstantaneous(starFormationRate,fuelAbundances,component,ageMinimum,ageMaximum)*starFormationRate
          ! Accumulate the metal return rate from this population at the future time.
          thisHistory%rates(iHistory,returnedMetalRateBeginIndex:returnedMetalRateEndIndex)=thisHistory%rates(iHistory&
               &,returnedMetalRateBeginIndex:returnedMetalRateEndIndex)+recyclingRate*fuelMetallicity
          ! Loop over all elements (and total metallicity).
          do iElement=1,elementsCount
             ! Get the metal yield rate.
             thisHistory%rates(iHistory,metalYieldRateBeginIndex+iElement-1)=thisHistory%rates(iHistory,metalYieldRateBeginIndex &
                  &+iElement-1)+IMF_Metal_Yield_Rate_NonInstantaneous(starFormationRate,fuelAbundances,component,ageMinimum,ageMaximum&
                  &,iElement)*starFormationRate
          end do
       end if
    end do
    return
  end subroutine Stellar_Population_Properties_Rates_Noninstantaneous

  subroutine Stellar_Population_Properties_Scales_Noninstantaneous(thisHistory,stellarMass,stellarAbundances)
    !% Set the scalings for error control on the absolute values of stellar population properties.
    use Histories
    use Stellar_Feedback
    use Abundances_Structure
    use Memory_Management
    implicit none
    double precision,          intent(in)                            :: stellarMass
    type(abundancesStructure), intent(in)                            :: stellarAbundances
    type(history),             intent(inout)                         :: thisHistory
    double precision,          parameter                             :: stellarMassMinimum      =1.0d0
    double precision,          parameter                             :: stellarAbundancesMinimum=1.0d0
    double precision,          dimension(elementsCount)              :: abundances
    double precision,          dimension(:            ), allocatable :: timeSteps
    integer                                                          :: scaleIndex

    ! Get timesteps.
    call thisHistory%timeSteps(timeSteps)

    ! Get abundances.
    call stellarAbundances%unpack(abundances)

    ! Set scaling factors for recycled mass.
    thisHistory   %scales(:,recycledRateIndex                       )=max(stellarMass           ,stellarMassMinimum      )                                       /timeSteps

    ! Set scaling factors for metal recycling rates.
    forall(scaleIndex=1:elementsCount)
       thisHistory%scales(:,returnedMetalRateBeginIndex-1+scaleIndex)=max(abundances(scaleIndex),stellarAbundancesMinimum)                                       /timeSteps
    end forall

    ! Set scaling factors for metal yield rates.
    forall(scaleIndex=1:elementsCount)
       thisHistory%scales(:,metalYieldRateBeginIndex   -1+scaleIndex)=max(abundances(scaleIndex),stellarAbundancesMinimum)                                       /timeSteps
    end forall
    
    ! Set scaling factors for energy input rates.
    thisHistory   %scales(:,energyInputRateIndex                    )=max(stellarMass           ,stellarMassMinimum      )*feedbackEnergyInputAtInfinityCanonical/timeSteps

    ! Destroy temporary array.
    call Dealloc_Array(timeSteps)

    return
  end subroutine Stellar_Population_Properties_Scales_Noninstantaneous

  subroutine Stellar_Population_Properties_History_Create_Noninstantaneous(thisNode,thisHistory)
    !% Create any history required for storing stellar population properties.
    use Histories
    use Numerical_Ranges
    use Tree_Nodes
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode
    type(history),   intent(inout)          :: thisHistory
    double precision                        :: timeBegin,timeEnd  

    ! Decide on start and end times for the history.
    timeBegin=Tree_Node_Time(thisNode)
    timeEnd  =historyStorageLatestTime

    ! Create the history.
    call thisHistory%create(historyCount,noninstantHistoryTimesCount,timeBegin,timeEnd,rangeTypeLogarithmic)

    return
  end subroutine Stellar_Population_Properties_History_Create_Noninstantaneous

end module Stellar_Population_Properties_Noninstantaneous
