!% Contains a module which implements stellar population properties with noninstantaneous recycling.

module Stellar_Population_Properties_Noninstantaneous
  !% Implements stellar population properties with noninstantaneous recycling.
  private
  public :: Stellar_Population_Properties_Noninstantaneous_Initialize

  ! Count of the number of histories required by this implementation.
  integer, parameter :: historyCount          =4
  ! Indices for histories.
  integer, parameter :: recycledRateIndex     =1
  integer, parameter :: returnedMetalRateIndex=2
  integer, parameter :: metalYieldRateIndex   =3
  integer, parameter :: energyInputRateIndex  =4

  ! Number of times to store in histories.
  integer            :: noninstantHistoryTimesCount

contains

  !# <stellarPopulationPropertiesMethod>
  !#  <unitName>Stellar_Population_Properties_Noninstantaneous_Initialize</unitName>
  !# </stellarPopulationPropertiesMethod>
  subroutine Stellar_Population_Properties_Noninstantaneous_Initialize(stellarPopulationPropertiesMethod &
       &,Stellar_Population_Properties_Rates_Get,Stellar_Population_Properties_History_Count_Get&
       &,Stellar_Population_Properties_History_Create_Do)
    !% Initializes the noninstantaneous recycling stellar population properties module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: stellarPopulationPropertiesMethod
    procedure(),          pointer, intent(inout) :: Stellar_Population_Properties_Rates_Get&
         &,Stellar_Population_Properties_History_Count_Get,Stellar_Population_Properties_History_Create_Do
    
    if (stellarPopulationPropertiesMethod == 'noninstantaneous') then
       Stellar_Population_Properties_Rates_Get         => Stellar_Population_Properties_Rates_Noninstantaneous    
       Stellar_Population_Properties_History_Count_Get => Stellar_Population_Properties_History_Count_Noninstantaneous    
       Stellar_Population_Properties_History_Create_Do => Stellar_Population_Properties_History_Create_Noninstantaneous
       !@ <inputParameter>
       !@   <name>noninstantHistoryTimesCount</name>
       !@   <defaultValue>10</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The number of times at which a galaxy's stellar properties history is stored.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('noninstantHistoryTimesCount',noninstantHistoryTimesCount,defaultValue=10)
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
  
  subroutine Stellar_Population_Properties_Rates_Noninstantaneous(starFormationRate,fuelAbundances,thisNode,thisHistory,stellarMassRate&
       &,stellarAbundancesRates,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
    !% Return an array of stellar population property rates of change given a star formation rate and fuel abundances.
    use Tree_Nodes
    use Tree_Node_Methods
    use Abundances_Structure
    use Histories
    use Star_Formation_IMF
    use Stellar_Population_Properties_Luminosities
    use Numerical_Interpolation
    use FGSL
    implicit none
    double precision,          intent(out)                 :: stellarMassRate,fuelMassRate,energyInputRate
    type(abundancesStructure), intent(out)                 :: stellarAbundancesRates,fuelAbundancesRates
    double precision,          intent(out),   dimension(:) :: stellarLuminositiesRates
    double precision,          intent(in)                  :: starFormationRate
    type(abundancesStructure), intent(in)                  :: fuelAbundances
    type(treeNode),            intent(inout), pointer      :: thisNode
    type(history),             intent(inout)               :: thisHistory
    integer                                                :: imfSelected,iHistory
    double precision                                       :: fuelMetallicity ,stellarMetalsRateOfChange ,fuelMetalsRateOfChange&
         &,ageMinimum,ageMaximum,currentTime,recyclingRate,metalReturnRate,historyFactors(2),metalYieldRate
    type(fgsl_interp_accel)                                :: interpolationAccelerator
    logical                                                :: interpolationReset

    ! Get the current time.
    currentTime=Tree_Node_Time(thisNode)

    ! Get interpolating factors in stellar population history.
    interpolationReset=.true.
    iHistory      =Interpolate_Locate                 (size(thisHistory%time),thisHistory%time,interpolationAccelerator,currentTime,interpolationReset)
    historyFactors=Interpolate_Linear_Generate_Factors(size(thisHistory%time),thisHistory%time,iHistory                ,currentTime                   )

    ! Interpolate to get recycling rate.
    recyclingRate  =Interpolate_Linear_Do(size(thisHistory%time),thisHistory%data(:,recycledRateIndex     ),iHistory,historyFactors)

    ! Interpolate to get metal return rate.
    metalReturnRate=Interpolate_Linear_Do(size(thisHistory%time),thisHistory%data(:,returnedMetalRateIndex),iHistory,historyFactors)

    ! Interpolate to get metal yield rate.
    metalYieldRate =Interpolate_Linear_Do(size(thisHistory%time),thisHistory%data(:,metalYieldRateIndex   ),iHistory,historyFactors)

    ! Interpolate to get metal yield rate.
    energyInputRate =Interpolate_Linear_Do(size(thisHistory%time),thisHistory%data(:,energyInputRateIndex ),iHistory,historyFactors)

    ! Get the metallicity of the fuel supply.
    fuelMetallicity=Abundances_Get_Metallicity(fuelAbundances)

    ! Set the stellar and fuel mass rates of change.
    stellarMassRate=starFormationRate-recyclingRate
    fuelMassRate   =-stellarMassRate

    ! Set the rates of change of the stellar and fuel metallicities.
    stellarMetalsRateOfChange=starFormationRate*fuelMetallicity-metalReturnRate
    fuelMetalsRateOfChange   =-stellarMetalsRateOfChange+metalYieldRate
    call stellarAbundancesRates%metallicitySet(stellarMetalsRateOfChange)
    call fuelAbundancesRates   %metallicitySet(fuelMetalsRateOfChange   )

    ! Get the IMF.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances)

    ! Set luminosity rates of change.
    if (size(stellarLuminositiesRates) > 0) stellarLuminositiesRates=starFormationRate&
         &*Stellar_Population_Luminosities_Get(imfSelected,currentTime,fuelAbundances)

    ! Set rates of change in the stellar populations properties future history.
    do iHistory=1,size(thisHistory%time)
       ! Find the age of the forming stellar population at the future time. We average over the time between successive timesteps
       ! to ensure that the rates will integrate to the correct values.
       !
       ! Minimum age is zero in the first timestep, and geometric mean of current and previous timestep otherwise.
       if (iHistory == 1) then
          ageMinimum=0.0d0
       else
          ageMinimum=max(dsqrt(thisHistory%time(iHistory)*thisHistory%time(iHistory-1))-currentTime,0.0d0)
       end if
       ! Maximum timestep is final time in final timestep, and geometric mean of current and next timestep otherwise.
       if (iHistory == size(thisHistory%time)) then
          ageMaximum=thisHistory%time(iHistory)-currentTime
       else
          ageMaximum=dsqrt(thisHistory%time(iHistory)*thisHistory%time(iHistory+1))-currentTime
       end if
       ! Check that it really is in the future.
       if (ageMaximum >= 0.0d0) then
          ! Get the recycling rate.
          recyclingRate=IMF_Recycling_Rate_NonInstantaneous(starFormationRate,fuelAbundances,ageMinimum,ageMaximum)&
               &*starFormationRate
          ! Accumulate the mass recycling rate from this population at the future time.
          thisHistory%rates(iHistory,recycledRateIndex     )=thisHistory%rates(iHistory,recycledRateIndex     ) +recyclingRate
          ! Accumulate the metal return rate from this population at the future time.
          thisHistory%rates(iHistory,returnedMetalRateIndex)=thisHistory%rates(iHistory,returnedMetalRateIndex) +recyclingRate&
               &*fuelMetallicity
          ! Get the metal yield rate.
          thisHistory%rates(iHistory,metalYieldRateIndex   )=thisHistory%rates(iHistory,metalYieldRateIndex   ) &
               &+IMF_Metal_Yield_Rate_NonInstantaneous(starFormationRate,fuelAbundances,ageMinimum,ageMaximum)*starFormationRate
          ! Get the (normalized) energy input rate.
          thisHistory%rates(iHistory,energyInputRateIndex  )=thisHistory%rates(iHistory,energyInputRateIndex  ) &
               &+IMF_Energy_Input_Rate_NonInstantaneous(starFormationRate,fuelAbundances,ageMinimum,ageMaximum)*starFormationRate
       end if
    end do

    return
  end subroutine Stellar_Population_Properties_Rates_Noninstantaneous

  subroutine Stellar_Population_Properties_History_Create_Noninstantaneous(thisNode,thisHistory)
    !% Create any history required for storing stellar population properties.
    use Histories
    use Numerical_Ranges
    use Tree_Nodes
    use Tree_Node_Methods
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
