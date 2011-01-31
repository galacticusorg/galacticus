!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module of hot halo tree node methods.

module Tree_Node_Methods_Hot_Halo
  !% Implement hot halo tree node methods.
  use Tree_Nodes
  use Components
  use Tree_Node_Methods_Hot_Halo_Data
  use Radiation_Structure
  private
  public :: Tree_Node_Methods_Hot_Halo_Initialize, Tree_Node_Methods_Hot_Halo_Thread_Initialize, Hot_Halo_Starve,&
       & Hot_Halo_Remove_Before_Satellite_Merging, Tree_Node_Hot_Halo_Promote, Hot_Halo_Subresolution_Initialize,&
       & Galacticus_Output_Tree_Hot_Halo_Standard, Galacticus_Output_Tree_Hot_Halo_Standard_Property_Count,&
       & Galacticus_Output_Tree_Hot_Halo_Standard_Names, Tree_Node_Hot_Halo_Reset_Standard,&
       & Tree_Node_Hot_Halo_Post_Evolve_Standard, Tree_Node_Methods_Hot_Halo_Standard_Dump, Hot_Halo_Scale_Set
  
  ! Internal count of abundances and molecules.
  integer                                     :: abundancesCount,moleculesCount
  double precision, allocatable, dimension(:) :: abundancesWork,abundancesParent,abundancesCoolingRate,abundancesReturnRate,abundancesHost
  !$omp threadprivate(abundancesWork,abundancesParent,abundancesCoolingRate,abundancesReturnRate,abundancesHost)
  double precision, allocatable, dimension(:) :: moleculesValue,moleculesCoolingRate,moleculesAccretionRate,moleculesChemicalRates
  !$omp threadprivate(moleculesValue,moleculesCoolingRate,moleculesAccretionRate,moleculesChemicalRates)

  ! Property indices.
  integer, parameter :: propertyCountBase=5, dataCount=0, historyCount=0
  integer            :: propertyCount
  integer, parameter :: massIndex                    =1
  integer, parameter :: angularMomentumIndex         =2
  integer, parameter :: outflowedMassIndex           =3
  integer, parameter :: outflowedAngularMomentumIndex=4
  integer, parameter :: unaccretedMassIndex          =5
  integer            :: hotAbundancesIndex      ,hotAbundancesIndexEnd
  integer            :: outflowedAbundancesIndex,outflowedAbundancesIndexEnd
  integer            :: hotMoleculesIndex       ,hotMoleculesIndexEnd

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Hot_Halo_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Hot_Halo_Unaccreted_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Hot_Halo_Angular_Momentum</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Hot_Halo_Abundances</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Hot_Halo_Molecules</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Hot_Halo_Outflowed_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Hot_Halo_Outflowed_Ang_Mom</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Hot_Halo_Outflowed_Abundances</methodName>
  !# </treeNodeMethodsPointer>

  ! Define pipes.
  !
  ! Pointer to procedure which will receive gas cooling out of the hot halo.
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Hot_Halo_Cooling_Mass_To</pipeName>
  !# </treeNodePipePointer>
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Hot_Halo_Cooling_Angular_Momentum_To</pipeName>
  !# </treeNodePipePointer>
  !# <treeNodePipePointer type="array">
  !#  <pipeName>Tree_Node_Hot_Halo_Cooling_Abundances_To</pipeName>
  !# </treeNodePipePointer>
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Hot_Halo_Outflow_Mass_To</pipeName>
  !# </treeNodePipePointer>
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Hot_Halo_Outflow_Angular_Momentum_To</pipeName>
  !# </treeNodePipePointer>
  !# <treeNodePipePointer type="array">
  !#  <pipeName>Tree_Node_Hot_Halo_Outflow_Abundances_To</pipeName>
  !# </treeNodePipePointer>
  ! Pointer to procedure which handles generic mass sinks in the hot halo.
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Hot_Halo_Hot_Gas_Sink</pipeName>
  !# </treeNodePipePointer>
  ! Pointer to procedure for heat input to the hot halo.
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Hot_Halo_Heat_Input</pipeName>
  !# </treeNodePipePointer>

  ! Configuration variables.
  logical          :: starveSatellites
  double precision :: hotHaloOutflowReturnRate

  ! Quantities stored to avoid repeated computation.
  logical          :: gotCoolingRate=.false., gotCoolingConversions=.false.
  double precision :: coolingRate,massHeatingRateRemaining,angularMomentumCoolingConversion
  !$omp threadprivate(gotCoolingRate,coolingRate,massHeatingRateRemaining,gotCoolingConversions,angularMomentumCoolingConversion)

  ! Radiation structure.
  type(radiationStructure) :: radiation
  !$omp threadprivate(radiation)

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Hot_Halo_Initialize</unitName>
  !#  <optionName default="standard">treeNodeMethodHotHalo</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Hot_Halo_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node hot halo methods module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    use Abundances_Structure
    use Molecular_Abundances_Structure
    use Memory_Management
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'standard') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Standard hot halo method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Get numbers of abundance and molecules properties.
       abundancesCount=Abundances_Property_Count()
       moleculesCount =Molecules_Property_Count ()

       ! Assign indices to properties.
       hotAbundancesIndex         =propertyCountBase                          +1
       hotAbundancesIndexEnd      =hotAbundancesIndex         +abundancesCount-1
       outflowedAbundancesIndex   =hotAbundancesIndexEnd                      +1
       outflowedAbundancesIndexEnd=outflowedAbundancesIndex   +abundancesCount-1
       hotMoleculesIndex          =outflowedAbundancesIndexEnd                +1
       hotMoleculesIndexEnd       =hotMoleculesIndex          +moleculesCount -1
       propertyCount              =hotMoleculesIndexEnd

       ! Set up procedure pointers.
       ! Unaccreted mass reservoir:
       Tree_Node_Hot_Halo_Unaccreted_Mass                   => Tree_Node_Hot_Halo_Unaccreted_Mass_Standard
       Tree_Node_Hot_Halo_Unaccreted_Mass_Set               => Tree_Node_Hot_Halo_Unaccreted_Mass_Set_Standard
       Tree_Node_Hot_Halo_Unaccreted_Mass_Rate_Adjust       => Tree_Node_Hot_Halo_Unaccreted_Mass_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Unaccreted_Mass_Rate_Compute      => Tree_Node_Rate_Rate_Compute_Dummy
       ! Hot gas reservoir:
       Tree_Node_Hot_Halo_Mass                              => Tree_Node_Hot_Halo_Mass_Standard
       Tree_Node_Hot_Halo_Mass_Set                          => Tree_Node_Hot_Halo_Mass_Set_Standard
       Tree_Node_Hot_Halo_Mass_Rate_Adjust                  => Tree_Node_Hot_Halo_Mass_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Mass_Rate_Compute                 => Tree_Node_Hot_Halo_Mass_Rate_Compute_Standard
       Tree_Node_Hot_Halo_Angular_Momentum                  => Tree_Node_Hot_Halo_Angular_Momentum_Standard
       Tree_Node_Hot_Halo_Angular_Momentum_Set              => Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard
       Tree_Node_Hot_Halo_Angular_Momentum_Rate_Adjust      => Tree_Node_Hot_Halo_Angular_Momentum_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Angular_Momentum_Rate_Compute     => Tree_Node_Hot_Halo_Angular_Momentum_Rate_Compute_Standard
       Tree_Node_Hot_Halo_Abundances                        => Tree_Node_Hot_Halo_Abundances_Standard
       Tree_Node_Hot_Halo_Abundances_Set                    => Tree_Node_Hot_Halo_Abundances_Set_Standard
       Tree_Node_Hot_Halo_Abundances_Rate_Adjust            => Tree_Node_Hot_Halo_Abundances_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Abundances_Rate_Compute           => Tree_Node_Hot_Halo_Abundances_Rate_Compute_Standard
       Tree_Node_Hot_Halo_Molecules                         => Tree_Node_Hot_Halo_Molecules_Standard
       Tree_Node_Hot_Halo_Molecules_Set                     => Tree_Node_Hot_Halo_Molecules_Set_Standard
       Tree_Node_Hot_Halo_Molecules_Rate_Adjust             => Tree_Node_Hot_Halo_Molecules_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Molecules_Rate_Compute            => Tree_Node_Hot_Halo_Molecules_Rate_Compute_Standard
       ! Outflowed gas reservoier:
       Tree_Node_Hot_Halo_Outflowed_Mass                    => Tree_Node_Hot_Halo_Outflowed_Mass_Standard
       Tree_Node_Hot_Halo_Outflowed_Mass_Set                => Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard
       Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Adjust        => Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Compute       => Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Compute_Standard
       Tree_Node_Hot_Halo_Outflowed_Ang_Mom                 => Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard
       Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set             => Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard
       Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Adjust     => Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Compute    => Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Compute_Standard
       Tree_Node_Hot_Halo_Outflowed_Abundances              => Tree_Node_Hot_Halo_Outflowed_Abundances_Standard
       Tree_Node_Hot_Halo_Outflowed_Abundances_Set          => Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard
       Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Adjust  => Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Compute => Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Compute_Standard
       
       ! Set externally connected pipes to dummy procedures if not already connected elsewhere.
       if (.not.associated(Tree_Node_Hot_Halo_Cooling_Mass_To            )) &
            & Tree_Node_Hot_Halo_Cooling_Mass_To             => Tree_Node_Rate_Adjust_Dummy
       if (.not.associated(Tree_Node_Hot_Halo_Cooling_Angular_Momentum_To)) &
            & Tree_Node_Hot_Halo_Cooling_Angular_Momentum_To => Tree_Node_Rate_Adjust_Dummy
       if (.not.associated(Tree_Node_Hot_Halo_Cooling_Abundances_To      )) &
            & Tree_Node_Hot_Halo_Cooling_Abundances_To       => Tree_Node_Rate_Adjust_Array_Dummy

       ! Set internally connected pipes to our procedures.
       Tree_Node_Hot_Halo_Outflow_Mass_To             => Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Outflow_Angular_Momentum_To => Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Outflow_Abundances_To       => Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Hot_Gas_Sink                => Tree_Node_Hot_Halo_Hot_Gas_Sink_Rate_Adjust_Standard
       Tree_Node_Hot_Halo_Heat_Input                  => Tree_Node_Hot_Halo_Heat_Input_Rate_Adjust_Standard

       ! Determine whether satellite nodes will be starved of gas.
       !@ <inputParameter>
       !@   <name>starveSatellites</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not the hot halo should be removed (``starved'') when a node becomes a satellite.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starveSatellites',starveSatellites,defaultValue=.true.)

       ! Get rate (in units of halo inverse dynamical time) at which outflowed gas returns to the hot gas reservoir.
       !@ <inputParameter>
       !@   <name>hotHaloOutflowReturnRate</name>
       !@   <defaultValue>1.26027</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies the rate at which reheated mass is returned to the hot phase in units of the inverse halo dynamical time.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloOutflowReturnRate',hotHaloOutflowReturnRate,defaultValue=1.26027d0)

    end if
    return
  end subroutine Tree_Node_Methods_Hot_Halo_Initialize
  
  !# <treeNodeCreateThreadInitialize>
  !#  <unitName>Tree_Node_Methods_Hot_Halo_Thread_Initialize</unitName>
  !# </treeNodeCreateThreadInitialize>
  subroutine Tree_Node_Methods_Hot_Halo_Thread_Initialize
    !% Initializes the tree node hot halo methods module.
    use Memory_Management
    implicit none
    
    ! Check if this implementation is selected.
    if (methodSelected.and..not.allocated(abundancesWork)) then

       ! Allocate work arrays for abundances.
       call Alloc_Array(abundancesWork       ,[abundancesCount])
       call Alloc_Array(abundancesHost       ,[abundancesCount])
       call Alloc_Array(abundancesParent     ,[abundancesCount])
       call Alloc_Array(abundancesCoolingRate,[abundancesCount])
       call Alloc_Array(abundancesReturnRate ,[abundancesCount])

       ! Allocate work arrays for molecules.
       call Alloc_Array(moleculesValue        ,[moleculesCount])
       call Alloc_Array(moleculesCoolingRate  ,[moleculesCount])
       call Alloc_Array(moleculesAccretionRate,[moleculesCount])
       call Alloc_Array(moleculesChemicalRates,[moleculesCount])

       ! Define the radiation component to include both the CMB and the intergalactic background.
       call radiation%define([radiationTypeCMB,radiationTypeIGB])

    end if
    return
  end subroutine Tree_Node_Methods_Hot_Halo_Thread_Initialize
  
  !# <calculationResetTask>
  !# <unitName>Tree_Node_Hot_Halo_Reset_Standard</unitName>
  !# </calculationResetTask>
  subroutine Tree_Node_Hot_Halo_Reset_Standard(thisNode)
    !% Remove memory of stored computed values as we're about to begin computing derivatives anew.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    gotCoolingRate       =.false.
    gotCoolingConversions=.false.
    return
  end subroutine Tree_Node_Hot_Halo_Reset_Standard

  !# <postEvolveTask>
  !#  <unitName>Tree_Node_Hot_Halo_Post_Evolve_Standard</unitName>
  !# </postEvolveTask>
  subroutine Tree_Node_Hot_Halo_Post_Evolve_Standard(thisNode)
    !% Do processing of the node required after evolution.
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: parentNode

    if (thisNode%isSatellite().and.starveSatellites) then
       ! Transfer any outflowed gas to the hot halo of the parent node.
       parentNode => thisNode%parentNode
       do while (parentNode%isSatellite())
          parentNode => parentNode%parentNode
       end do
       call Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard(parentNode,Tree_Node_Hot_Halo_Outflowed_Mass_Standard(parentNode) &
            &+Tree_Node_Hot_Halo_Outflowed_Mass_Standard(thisNode))
       call Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard(thisNode,0.0d0)
       call Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard(parentNode&
            &,Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(parentNode)+Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(thisNode))
       call Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard(thisNode,0.0d0)
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(parentNode,abundancesParent)
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode  ,abundancesWork  )
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard(parentNode,abundancesWork+abundancesParent)
       abundancesWork(:)=0.0d0
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard(thisNode,abundancesWork)
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Post_Evolve_Standard

  double precision function Tree_Node_Hot_Halo_Mass_Standard(thisNode)
    !% Return the node hot halo mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       Tree_Node_Hot_Halo_Mass_Standard=thisNode%components(thisIndex)%properties(massIndex,propertyValue)
    else
       Tree_Node_Hot_Halo_Mass_Standard=0.0d0
    end if
    return
  end function Tree_Node_Hot_Halo_Mass_Standard

  subroutine Tree_Node_Hot_Halo_Mass_Set_Standard(thisNode,mass)
    !% Set the node hot halo mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(massIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Hot_Halo_Mass_Set_Standard

  subroutine Tree_Node_Hot_Halo_Hot_Gas_Sink_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Account for a sink of gaseous material in the hot halo hot gas.
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    
    ! If no hot halo component currently exists and we have some sink from then there is a problem.
    ! spheroid. If sink has zero rate, just return instead.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) call Galacticus_Error_Report('Tree_Node_Hot_Halo_Hot_Gas_Sink_Rate_Adjust_Standard','mass sink from non-existant hot halo')
       return
    end if
    ! Trap cases where an attempt is made to add gas via this sink function.
    if (rateAdjustment > 0.0d0) call Galacticus_Error_Report('Tree_Node_Hot_Halo_Hot_Gas_Sink_Rate_Adjust_Standard','attempt to add mass via sink in hot halo')

    ! Proportionally adjust the rates of all components of the hot gas reservoir.
    call Tree_Node_Hot_Halo_Hot_Gas_All_Rate_Adjust_Standard(thisNode,rateAdjustment)

    return
  end subroutine Tree_Node_Hot_Halo_Hot_Gas_Sink_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Hot_Gas_All_Rate_Adjust_Standard(thisNode,gasMassRate)
    !% Adjusts the rates of all components of the hot gas reservoir under the assumption of uniformly distributed properties
    !% (e.g. fully-mixed metals).
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: gasMassRate
    integer                                  :: thisIndex
    double precision                         :: gasMass

    ! Exit immediately for zero rate.
    if (gasMassRate == 0.0d0) return

    ! Get the index of the component.
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)    

    ! Get the gas mass present.
    gasMass=thisNode%components(thisIndex)%properties(massIndex,propertyValue)

    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0) then
       ! Mass.
       thisNode         %components(thisIndex)%properties(massIndex                               ,propertyDerivative) &
            & = thisNode%components(thisIndex)%properties(massIndex                               ,propertyDerivative) &
            &  +gasMassRate
       ! Angular momentum.
       thisNode         %components(thisIndex)%properties(angularMomentumIndex                    ,propertyDerivative) &
            & = thisNode%components(thisIndex)%properties(angularMomentumIndex                    ,propertyDerivative) &
            &  +thisNode%components(thisIndex)%properties(angularMomentumIndex                    ,propertyValue     ) &
            &  *(gasMassRate/gasMass)
       ! Metal abundances.
       thisNode         %components(thisIndex)%properties(hotAbundancesIndex:hotAbundancesIndexEnd,propertyDerivative) &
            & = thisNode%components(thisIndex)%properties(hotAbundancesIndex:hotAbundancesIndexEnd,propertyDerivative) & 
            &  +thisNode%components(thisIndex)%properties(hotAbundancesIndex:hotAbundancesIndexEnd,propertyValue     ) &
            &  *(gasMassRate/gasMass)
       ! Molecular abundances.
       thisNode         %components(thisIndex)%properties(hotMoleculesIndex:hotMoleculesIndexEnd  ,propertyDerivative) &
            & = thisNode%components(thisIndex)%properties(hotMoleculesIndex:hotMoleculesIndexEnd  ,propertyDerivative) & 
            &  *thisNode%components(thisIndex)%properties(hotMoleculesIndex:hotMoleculesIndexEnd  ,propertyValue     ) &
            &  +(gasMassRate/gasMass)
    end if
  
    return
  end subroutine Tree_Node_Hot_Halo_Hot_Gas_All_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Heat_Input_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% An incoming pipe that for sources of heating to the hot halo.
    use Galacticus_Error
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    procedure(),      pointer                :: interruptProcedurePassed
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    double precision                         :: massHeatingRate
    
    ! Trap cases where an attempt is made to remove energy via this input function.
    if (rateAdjustment < 0.0d0) call Galacticus_Error_Report('Tree_Node_Hot_Halo_Heat_Input_Rate_Adjust_Standard','attempt to remove energy via heat input pipe in hot halo')

    ! Get a local copy of the interrupt procedure.
    interruptProcedurePassed => interruptProcedure

    ! Get the index of the component.
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)    

    ! Ensure that the cooling rate has been computed.
    call Get_Cooling_Rate(thisNode)

    ! Compute mass heating rate from energy heating rate, but don't allow it to exceed the remaining budget.
    massHeatingRate=min(rateAdjustment/Dark_Matter_Halo_Virial_Velocity(thisNode)**2,massHeatingRateRemaining)
    
    ! Update the remaining budget of allowed mass heating rate.
    massHeatingRateRemaining=max(massHeatingRateRemaining-massHeatingRate,0.0d0)

    ! Call routine to apply this mass heating rate to all hot halo cooling pipes.
    call Hot_Halo_Standard_Push_To_Cooling_Pipes(thisNode,interrupt,interruptProcedurePassed,-massHeatingRate)

    ! Return our local copy of the interrupt procedure.
    interruptProcedure => interruptProcedurePassed

    return
  end subroutine Tree_Node_Hot_Halo_Heat_Input_Rate_Adjust_Standard

  subroutine Hot_Halo_Standard_Push_To_Cooling_Pipes(thisNode,interrupt,interruptProcedure,massRate)
    !% Push mass through the cooling pipes (along with appropriate amounts of metals and angular momentum) at the given rate.
    use Dark_Matter_Halo_Spins
    use Cooling_Radii
    use Dark_Matter_Profiles
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: massRate
    procedure(),      pointer                :: interruptProcedurePassed
    double precision                         :: angularMomentumCoolingRate

    ! Ignore zero rates.
    if (massRate /= 0.0d0) then
       
       ! Get a local copy of the interrupt procedure.
       interruptProcedurePassed => interruptProcedure

       ! Remove mass from the hot component.
       call Tree_Node_Hot_Halo_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedurePassed,-massRate)
       ! Pipe the mass rate to whatever component claimed it.
       if (associated(Tree_Node_Hot_Halo_Cooling_Mass_To)) call Tree_Node_Hot_Halo_Cooling_Mass_To(thisNode,interrupt &
            &,interruptProcedurePassed,massRate)
       
       ! Get the corresponding rate of change of angular momentum.
       if (.not.gotCoolingConversions) then
          angularMomentumCoolingConversion=Cooling_Radius(thisNode)*Dark_Matter_Profile_Rotation_Normalization(thisNode)&
               &*Tree_Node_Hot_Halo_Angular_Momentum(thisNode)/Tree_Node_Hot_Halo_Mass(thisNode)
          
          ! Flag that cooling conversion factors have now been computed.
          gotCoolingConversions=.true.
       end if
       angularMomentumCoolingRate=massRate*angularMomentumCoolingConversion
       call Tree_Node_Hot_Halo_Angular_Momentum_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedurePassed, &
            &-angularMomentumCoolingRate)
       ! Pipe the cooling rate to which ever component claimed it.
       if (associated(Tree_Node_Hot_Halo_Cooling_Angular_Momentum_To)) call&
            & Tree_Node_Hot_Halo_Cooling_Angular_Momentum_To(thisNode,interrupt,interruptProcedurePassed&
            &,sign(angularMomentumCoolingRate,massRate))
       
       ! Get the rate of change of abundances.
       call Tree_Node_Hot_Halo_Abundances_Standard(thisNode,abundancesWork)
       abundancesCoolingRate=massRate*abundancesWork/Tree_Node_Hot_Halo_Mass(thisNode)
       call Tree_Node_Hot_Halo_Abundances_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedurePassed,-abundancesCoolingRate)
       ! Pipe the cooling rate to which ever component claimed it.
       if (associated(Tree_Node_Hot_Halo_Cooling_Abundances_To)) call Tree_Node_Hot_Halo_Cooling_Abundances_To(thisNode,interrupt&
            & ,interruptProcedurePassed,abundancesCoolingRate)
       
       ! Return our local copy of the interrupt procedure.
       interruptProcedure => interruptProcedurePassed

    end if
    return
  end subroutine Hot_Halo_Standard_Push_To_Cooling_Pipes

  subroutine Tree_Node_Hot_Halo_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node hot halo mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(massIndex,propertyDerivative)=thisNode%components(thisIndex)%properties(massIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Mass_Rate_Adjust_Standard

  double precision function Tree_Node_Hot_Halo_Unaccreted_Mass_Standard(thisNode)
    !% Return the node unaccreted hot halo mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       Tree_Node_Hot_Halo_Unaccreted_Mass_Standard=thisNode%components(thisIndex)%properties(unaccretedMassIndex,propertyValue)
    else
       Tree_Node_Hot_Halo_Unaccreted_Mass_Standard=0.0d0
    end if
    return
  end function Tree_Node_Hot_Halo_Unaccreted_Mass_Standard

  subroutine Tree_Node_Hot_Halo_Unaccreted_Mass_Set_Standard(thisNode,mass)
    !% Set the node unaccreted hot halo mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(unaccretedMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Hot_Halo_Unaccreted_Mass_Set_Standard

  subroutine Tree_Node_Hot_Halo_Unaccreted_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node hot halo unaccreted mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(unaccretedMassIndex,propertyDerivative)=thisNode%components(thisIndex)%properties(unaccretedMassIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Unaccreted_Mass_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Mass_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Accretion_Halos
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    procedure(),    pointer                :: interruptProcedurePassed
    double precision                       :: massAccretionRate,failedMassAccretionRate

    ! Find the rate of gas mass accretion onto the halo.
    massAccretionRate      =Halo_Baryonic_Accretion_Rate       (thisNode)
    failedMassAccretionRate=Halo_Baryonic_Failed_Accretion_Rate(thisNode)

    ! If no hot halo component currently exists and we have some accretion then interrupt and create a hot halo.
    if (.not.thisNode%componentExists(componentIndex)) then    
       if (massAccretionRate /= 0.0d0 .or. failedMassAccretionRate /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Hot_Halo_Create
       end if
       return
    end if
    
    call Tree_Node_Hot_Halo_Mass_Rate_Adjust_Standard           (thisNode,interrupt,interruptProcedurePassed,massAccretionRate      )
    call Tree_Node_Hot_Halo_Unaccreted_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedurePassed,failedMassAccretionRate)

    ! Next compute the cooling rate in this halo.
    call Get_Cooling_Rate(thisNode)
    ! Pipe the cooling rate to which ever component claimed it.
    call Hot_Halo_Standard_Push_To_Cooling_Pipes(thisNode,interrupt,interruptProcedurePassed,coolingRate)

    ! Return a copy of our local interrupt pointer.
    interruptProcedure => interruptProcedurePassed

    return
  end subroutine Tree_Node_Hot_Halo_Mass_Rate_Compute_Standard

  double precision function Tree_Node_Hot_Halo_Angular_Momentum_Standard(thisNode)
    !% Return the node hot halo angular momentum.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       Tree_Node_Hot_Halo_Angular_Momentum_Standard=thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyValue)
    else
       Tree_Node_Hot_Halo_Angular_Momentum_Standard=0.0d0
    end if
    return
  end function Tree_Node_Hot_Halo_Angular_Momentum_Standard

  subroutine Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard(thisNode,angularMomentum)
    !% Set the node hot halo angular momentum.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: angularMomentum
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyValue)=angularMomentum
    return
  end subroutine Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard

  subroutine Tree_Node_Hot_Halo_Angular_Momentum_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node hot halo mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(angularMomentumIndex ,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Angular_Momentum_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Angular_Momentum_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node angular momentum rate of change. Note that the rate of change due to cooling is not included here
    !% as it is handled elsewhere.
    use Dark_Matter_Halo_Spins
    use Cooling_Radii
    use Dark_Matter_Profiles
    use Accretion_Halos
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    procedure(),    pointer                :: interruptProcedurePassed
    double precision                       :: massAccretionRate,angularMomentumAccretionRate,angularMomentumCoolingRate

    ! Find the rate of gas mass accretion onto the halo.
    massAccretionRate=Halo_Baryonic_Accretion_Rate(thisNode)
    ! If no hot halo component currently exists and we have some accretion then interrupt and create a hot halo.
    if (.not.thisNode%componentExists(componentIndex)) then    
       if (massAccretionRate /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Hot_Halo_Create
       end if
       return
    end if

    if (massAccretionRate > 0.0d0) then
       angularMomentumAccretionRate=Dark_Matter_Halo_Angular_Momentum_Growth_Rate(thisNode)*(massAccretionRate &
            &/Tree_Node_Mass_Accretion_Rate(thisNode))
       call Tree_Node_Hot_Halo_Angular_Momentum_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure &
            &,angularMomentumAccretionRate)
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Angular_Momentum_Rate_Compute_Standard

  subroutine Tree_Node_Hot_Halo_Abundances_Standard(thisNode,abundances)
    !% Return the node hot halo abundances.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: abundances(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       abundances(:)=thisNode%components(thisIndex)%properties(hotAbundancesIndex:hotAbundancesIndexEnd,propertyValue)
    else
       abundances(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Abundances_Standard

  subroutine Tree_Node_Hot_Halo_Abundances_Set_Standard(thisNode,abundances)
    !% Set the node hot halo angular momentum.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: abundances(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(hotAbundancesIndex:hotAbundancesIndexEnd,propertyValue)=abundances
    return
  end subroutine Tree_Node_Hot_Halo_Abundances_Set_Standard

  subroutine Tree_Node_Hot_Halo_Abundances_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node hot halo mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment(:)
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(hotAbundancesIndex:hotAbundancesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(hotAbundancesIndex:hotAbundancesIndexEnd,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Abundances_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Abundances_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Dark_Matter_Halo_Spins
    use Cooling_Radii
    use Dark_Matter_Profiles
    use Abundances_Structure
    use Accretion_Halos
    implicit none
    type(treeNode),            pointer, intent(inout)     :: thisNode
    logical,                            intent(inout)     :: interrupt
    procedure(),               pointer, intent(inout)     :: interruptProcedure
    procedure(),               pointer                    :: interruptProcedurePassed
    double precision,          dimension(abundancesCount) :: abundances,abundancesCoolingRate,abundancesAccretionRate
    type(abundancesStructure), save                       :: accretionRateAbundances
    !$omp threadprivate(accretionRateAbundances)

    ! Get the rate at which abundances are accreted onto this halo.
    call Halo_Baryonic_Accretion_Rate_Abundances(thisNode,accretionRateAbundances)
    call accretionRateAbundances%unpack(abundancesAccretionRate)
    if (any(abundancesAccretionRate /= 0.0d0)) then
       if (.not.thisNode%componentExists(componentIndex)) then    
          interrupt=.true.
          interruptProcedure => Hot_Halo_Create
          return
       end if
       call Tree_Node_Hot_Halo_Abundances_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,abundancesAccretionRate)
    end if
    
    ! Get the cooling rate.
    call Get_Cooling_Rate(thisNode)
    ! Adjust abundance rates for cooling.
    if (coolingRate > 0.0d0) then
       call Tree_Node_Hot_Halo_Abundances_Standard(thisNode,abundances)
       abundancesCoolingRate=coolingRate*abundances/Tree_Node_Hot_Halo_Mass(thisNode)
       call Tree_Node_Hot_Halo_Abundances_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,-abundancesCoolingRate)
       ! Pipe the cooling rate to which ever component claimed it.
       if (associated(Tree_Node_Hot_Halo_Cooling_Abundances_To)) call Tree_Node_Hot_Halo_Cooling_Abundances_To(thisNode,interrupt&
            &,interruptProcedurePassed,abundancesCoolingRate)

       ! Point the interrupt procedure that is returned to our internal one.
       interruptProcedure => interruptProcedurePassed
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Abundances_Rate_Compute_Standard

  subroutine Tree_Node_Hot_Halo_Molecules_Standard(thisNode,molecules)
    !% Return the node hot halo molecules.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: molecules(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       molecules(:)=thisNode%components(thisIndex)%properties(hotMoleculesIndex:hotMoleculesIndexEnd,propertyValue)
    else
       molecules(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Molecules_Standard

  subroutine Tree_Node_Hot_Halo_Molecules_Set_Standard(thisNode,molecules)
    !% Set the node hot halo angular momentum.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: molecules(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(hotMoleculesIndex:hotMoleculesIndexEnd,propertyValue)=molecules
    return
  end subroutine Tree_Node_Hot_Halo_Molecules_Set_Standard

  subroutine Tree_Node_Hot_Halo_Molecules_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node hot halo mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment(:)
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(hotMoleculesIndex:hotMoleculesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(hotMoleculesIndex:hotMoleculesIndexEnd,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Molecules_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Molecules_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node molecules rate of change.
    use Dark_Matter_Halo_Scales
    use Molecular_Abundances_Structure
    use Accretion_Halos
    use Molecular_Reaction_Rates
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Atomic
    use Numerical_Constants_Math
    use Numerical_Constants_Prefixes
    implicit none
    type(treeNode),                     pointer, intent(inout) :: thisNode
    logical,                                     intent(inout) :: interrupt
    procedure(),                        pointer, intent(inout) :: interruptProcedure
    type(molecularAbundancesStructure), save                   :: accretionRateMolecules,molecularMasses,molecularDensities&
         &,molecularDensitiesRates,molecularMassesRates
    !$omp threadprivate(accretionRateMolecules,molecularMasses,molecularDensities,molecularDensitiesRates,molecularMassesRates)
    double precision                                           :: massToDensityConversion,temperature

    ! If no molecules are being tracked, simply return.
    if (moleculesCount == 0) return

    ! Get the rate at which molecules are accreted onto this halo.
    call Halo_Baryonic_Accretion_Rate_Molecules(thisNode,accretionRateMolecules)
    call accretionRateMolecules%unpack(moleculesAccretionRate)
    if (any(moleculesAccretionRate /= 0.0d0)) then
       if (.not.thisNode%componentExists(componentIndex)) then    
          interrupt=.true.
          interruptProcedure => Hot_Halo_Create
          return
       end if
       call Tree_Node_Hot_Halo_Molecules_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,moleculesAccretionRate)
    end if
    
    ! Ensure that the cooling rate has been computed.
    call Get_Cooling_Rate(thisNode)
    ! For non-zero cooling rate, adjust the rates of molecule masses.
    if (coolingRate > 0.0d0) then
       ! Get the masses of molecule components.
       call Tree_Node_Hot_Halo_Molecules_Standard(thisNode,moleculesValue)
       ! Compute the rate at which molecules are lost via cooling.
       moleculesCoolingRate=coolingRate*moleculesValue/Tree_Node_Hot_Halo_Mass(thisNode)
       ! Adjust the rates of molecule masses accordingly.
       call Tree_Node_Hot_Halo_Molecules_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,-moleculesCoolingRate)     
    end if

    ! Compute the rates of change of molecular masses due to chemical reactions.
    ! Get the temperature of the hot reservoir.
    temperature=Dark_Matter_Halo_Virial_Temperature(thisNode)
    ! Set the radiation background.
    call radiation%set(thisNode)
    ! Get the masses of molecules.
    call Tree_Node_Hot_Halo_Molecules_Standard(thisNode,moleculesValue)
    ! Truncate masses to zero to avoid unphysical behavior.
    where (moleculesValue < 0.0d0)
       moleculesValue=0.0d0
    end where
    call molecularMasses%pack(moleculesValue)
    ! Scale all molecular masses by their mass in atomic mass units to get a number density.
    call molecularMasses%massToNumber(molecularDensities)
    ! Compute factor converting mass of molecules in (M_Solar/M_Atomic) to number density in cm^-3.
    massToDensityConversion=3.0d0*massSolar/atomicMassUnit/4.0d0/Pi/(hecto*megaParsec*Dark_Matter_Halo_Virial_Radius(thisNode))**3
    ! Convert to number density.
    call molecularDensities%multiply(massToDensityConversion)
    ! Compute the molecular reaction rates.
    call Molecular_Reaction_Rate(molecularDensitiesRates,temperature,molecularDensities,radiation)
    ! Convert to mass change rates.
    call molecularDensitiesRates%numberToMass(molecularMassesRates)
    call molecularMassesRates%multiply(gigaYear/massToDensityConversion)
    ! Adjust rates appropriately.
    call molecularMassesRates%unpack(moleculesChemicalRates)
    call Tree_Node_Hot_Halo_Molecules_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,moleculesChemicalRates)
    return
  end subroutine Tree_Node_Hot_Halo_Molecules_Rate_Compute_Standard

  double precision function Tree_Node_Hot_Halo_Outflowed_Mass_Standard(thisNode)
    !% Return the node hot halo mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       Tree_Node_Hot_Halo_Outflowed_Mass_Standard=thisNode%components(thisIndex)%properties(outflowedMassIndex,propertyValue)
    else
       Tree_Node_Hot_Halo_Outflowed_Mass_Standard=0.0d0
    end if
    return
  end function Tree_Node_Hot_Halo_Outflowed_Mass_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard(thisNode,mass)
    !% Set the node hot halo mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(outflowedMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node hot halo mass rate of change.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(outflowedMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(outflowedMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Dark_Matter_Halo_Scales
    use Abundances_Structure
    use Molecular_Abundances_Structure
    use Numerical_Constants_Math
    use Numerical_Constants_Atomic
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Ionization_States
    implicit none
    type(treeNode),                     pointer, intent(inout) :: thisNode
    logical,                                     intent(inout) :: interrupt
    procedure(),                        pointer, intent(inout) :: interruptProcedure
    procedure(),                        pointer                :: interruptProcedurePassed
    type(abundancesStructure),          save                   :: outflowedAbundances
    type(molecularAbundancesStructure), save                   :: molecularDensities,molecularRates,molecularMasses
    !$omp threadprivate(outflowedAbundances,molecularDensities,molecularRates,molecularMasses)
    double precision                       :: massReturnRate,temperature,hydrogenByMass,massToDensityConversion&
         &,numberDensityHydrogen,outflowedMass

    ! If a hot component exists, compute rate of return of outflowed gas to the hot gas reservoir.
    if (thisNode%componentExists(componentIndex).and.(.not.starveSatellites.or..not.thisNode%isSatellite())) then    
       outflowedMass=Tree_Node_Hot_Halo_Outflowed_Mass_Standard(thisNode)
       massReturnRate=hotHaloOutflowReturnRate*outflowedMass/Dark_Matter_Halo_Dynamical_Timescale(thisNode)
       call Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,-massReturnRate)
       call Tree_Node_Hot_Halo_Mass_Rate_Adjust_Standard          (thisNode,interrupt,interruptProcedure, massReturnRate)

       ! If we have a non-zero return rate, compute associated molecular rates.
       if (moleculesCount > 0 .and. massReturnRate /= 0.0d0) then
          
          ! Compute coefficient in conversion of mass to density for this node.
          massToDensityConversion=massSolar/4.0d0/Pi/(hecto*megaParsec*Dark_Matter_Halo_Virial_Radius(thisNode))**3
          
          ! Get the abundances of the outflowed material.
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode,abundancesWork)
          ! Convert to mass fractions and pack.
          abundancesWork=abundancesWork/outflowedMass
          call outflowedAbundances%pack(abundancesWork)
          
          ! Get the hydrogen mass fraction in outflowed gas.
          hydrogenByMass=outflowedAbundances%hydrogenMassFraction()
          
          ! Compute the temperature and density of material in the hot halo.
          temperature          =Dark_Matter_Halo_Virial_Temperature(thisNode)
          numberDensityHydrogen=hydrogenByMass*outflowedMass*massToDensityConversion/atomicMassUnit/atomicMassHydrogen
          
          ! Set the radiation field.
          call radiation%set(thisNode)
          
          ! Get the molecule densities.
          call Molecular_Densities(molecularDensities,temperature,numberDensityHydrogen,outflowedAbundances,radiation)
          
          ! Convert from densities to masses.
          call molecularDensities%numberToMass(molecularMasses)
          molecularRates=molecularMasses
          call molecularRates%multiply(massReturnRate*hydrogenByMass/numberDensityHydrogen/atomicMassHydrogen)
          
          ! Compute the rate at which molecules are returned to the hot reservoir.
          call molecularRates%unpack(moleculesValue)
          call Tree_Node_Hot_Halo_Molecules_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,moleculesValue)
          
       end if

    end if
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Mass_Rate_Compute_Standard

  double precision function Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(thisNode)
    !% Return the node hot halo angular momentum.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard&
            &=thisNode%components(thisIndex)%properties(outflowedAngularMomentumIndex,propertyValue)
    else
       Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard=0.0d0
    end if
    return
  end function Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard(thisNode,angularMomentum)
    !% Set the node hot halo angular momentum.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: angularMomentum
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(outflowedAngularMomentumIndex,propertyValue)=angularMomentum
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node hot halo mass rate of change.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(outflowedAngularMomentumIndex,propertyDerivative) &
         &=thisNode%components(thisIndex)%properties(outflowedAngularMomentumIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    procedure(),    pointer                :: interruptProcedurePassed
    double precision                       :: angularMomentumReturnRate

    ! If a hot halo exists, compute rate of return of outflowed gas angular momentum to the hot gas reservoir.
    if (thisNode%componentExists(componentIndex).and.(.not.starveSatellites.or..not.thisNode%isSatellite())) then
       angularMomentumReturnRate=hotHaloOutflowReturnRate*Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(thisNode)&
            &/Dark_Matter_Halo_Dynamical_Timescale(thisNode)
       call Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure, &
            &-angularMomentumReturnRate)
       call Tree_Node_Hot_Halo_Angular_Momentum_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,&
            & angularMomentumReturnRate)
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Rate_Compute_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode,abundances)
    !% Return the node hot halo angular momentum.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: abundances(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       abundances(:)=thisNode%components(thisIndex)%properties(outflowedAbundancesIndex:outflowedAbundancesIndexEnd,propertyValue)
    else
       abundances(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Abundances_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard(thisNode,abundances)
    !% Set the node hot halo angular momentum.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: abundances(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(outflowedAbundancesIndex:outflowedAbundancesIndexEnd,propertyValue)=abundances
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node hot halo mass rate of change.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment(:)
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%properties(outflowedAbundancesIndex:outflowedAbundancesIndexEnd,propertyDerivative) &
         &=thisNode%components(thisIndex)%properties(outflowedAbundancesIndex:outflowedAbundancesIndexEnd,propertyDerivative)&
         &+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Adjust_Standard

  subroutine Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    procedure(),      pointer                :: interruptProcedurePassed

    ! If a hot halo exists, compute rate of return of outflowed gas angular momentum to the hot gas reservoir.
    if (thisNode%componentExists(componentIndex).and.(.not.starveSatellites.or..not.thisNode%isSatellite())) then
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode,abundancesWork)
       abundancesReturnRate=hotHaloOutflowReturnRate*abundancesWork/Dark_Matter_Halo_Dynamical_Timescale(thisNode)
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure, &
            &-abundancesReturnRate)
       call Tree_Node_Hot_Halo_Abundances_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,abundancesReturnRate)
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Outflowed_Abundances_Rate_Compute_Standard

  !# <scaleSetTask>
  !#  <unitName>Hot_Halo_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Hot_Halo_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: scaleMassRelative=1.0d-6
    integer                                  :: thisIndex
    double precision                         :: massVirial,radiusVirial,velocityVirial

    ! Determine if method is active and a hot halo component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       massVirial    =Tree_Node_Mass                  (thisNode)
       radiusVirial  =Dark_Matter_Halo_Virial_Radius  (thisNode)
       velocityVirial=Dark_Matter_Halo_Virial_Velocity(thisNode)
       thisNode%components(thisIndex)%properties(massIndex                                           ,propertyScale)=massVirial                            *scaleMassRelative
       thisNode%components(thisIndex)%properties(outflowedMassIndex                                  ,propertyScale)=massVirial                            *scaleMassRelative
       thisNode%components(thisIndex)%properties(unaccretedMassIndex                                 ,propertyScale)=massVirial                            *scaleMassRelative
       thisNode%components(thisIndex)%properties(hotAbundancesIndex      :hotAbundancesIndexEnd      ,propertyScale)=massVirial                            *scaleMassRelative
       thisNode%components(thisIndex)%properties(outflowedAbundancesIndex:outflowedAbundancesIndexEnd,propertyScale)=massVirial                            *scaleMassRelative
       thisNode%components(thisIndex)%properties(hotMoleculesIndex       :hotMoleculesIndexEnd       ,propertyScale)=massVirial                            *scaleMassRelative
       thisNode%components(thisIndex)%properties(angularMomentumIndex                                ,propertyScale)=massVirial*radiusVirial*velocityVirial*scaleMassRelative
       thisNode%components(thisIndex)%properties(outflowedAngularMomentumIndex                       ,propertyScale)=massVirial*radiusVirial*velocityVirial*scaleMassRelative
    end if
    return
  end subroutine Hot_Halo_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Hot_Halo_Subresolution_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Hot_Halo_Subresolution_Initialize(thisNode)
    !% Initialize the contents of the hot halo component for any sub-resolution accretion (i.e. the gas that would have been
    !% accreted if the merger tree had infinite resolution).
    use Accretion_Halos
    use Dark_Matter_Halo_Spins
    use Molecular_Abundances_Structure
    use Abundances_Structure
    implicit none
    type(treeNode),                     pointer, intent(inout)     :: thisNode
    type(abundancesStructure),          save                       :: accretedAbundances
    type(molecularAbundancesStructure), save                       :: accretedMolecules
    !$omp threadprivate(accretedAbundances,accretedMolecules)
    double precision,                   dimension(abundancesCount) :: abundances
    double precision,                   dimension(moleculesCount ) :: molecules
    double precision                                               :: hotHaloMass,failedHotHaloMass,angularMomentum

    ! If this method is selected and the node has no child then initialize it.
    if (methodSelected.and..not.associated(thisNode%childNode)) then
       ! Get the mass of hot gas accreted and the mass that failed to accrete.
       hotHaloMass      =Halo_Baryonic_Accreted_Mass       (thisNode)
       failedHotHaloMass=Halo_Baryonic_Failed_Accreted_Mass(thisNode)
       ! If either is non-zero, then create a hot halo component and add these masses to it.
       if (hotHaloMass > 0.0d0 .or. failedHotHaloMass > 0.0d0) then
          call Hot_Halo_Create                                 (thisNode                   )
          call Tree_Node_Hot_Halo_Mass_Set_Standard            (thisNode,hotHaloMass       )
          call Tree_Node_Hot_Halo_Unaccreted_Mass_Set_Standard (thisNode,failedHotHaloMass )
          ! Also add the appropriate angular momentum.
          angularMomentum=hotHaloMass*Dark_Matter_Halo_Angular_Momentum(thisNode)/Tree_Node_Mass(thisNode)
          call Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard(thisNode,angularMomentum   )
          ! Add the appropriate abundances.
          call Halo_Baryonic_Accreted_Abundances               (thisNode,accretedAbundances)
          call accretedAbundances%unpack(abundances)
          call Tree_Node_Hot_Halo_Abundances_Set_Standard      (thisNode,abundances        )
          ! Also add the appropriate molecular masses.
          call Halo_Baryonic_Accreted_Molecules                (thisNode,accretedMolecules )
          call accretedMolecules%unpack(molecules)
          call Tree_Node_Hot_Halo_Molecules_Set_Standard       (thisNode,molecules         )
       end if
    end if
    return
  end subroutine Hot_Halo_Subresolution_Initialize

  !# <nodeMergerTask>
  !#  <unitName>Hot_Halo_Starve</unitName>
  !# </nodeMergerTask>
  subroutine Hot_Halo_Starve(thisNode)
    !% Starve {\tt thisNode} by transferring its hot halo to its parent.
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: parentNode
    double precision, dimension(moleculesCount ) :: molecules ,moleculesParent

    ! Determine if method is active and a hot halo component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       
       ! Find the parent node.
       parentNode => thisNode%parentNode

       ! Any gas that failed to be accreted by this halo is always transferred to the parent.
       call Tree_Node_Hot_Halo_Unaccreted_Mass_Set_Standard(parentNode, Tree_Node_Hot_Halo_Unaccreted_Mass_Standard(parentNode) &
            &                                                          +Tree_Node_Hot_Halo_Unaccreted_Mass_Standard(thisNode  ))
       call Tree_Node_Hot_Halo_Unaccreted_Mass_Set_Standard(thisNode  , 0.0d0                                                  )
       
       ! Determine if starvation is to be applied.
       if (starveSatellites) then
          
          ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate to
          ! this hot halo (and will be moved to the parent at the end of the evolution timestep).
          call Tree_Node_Hot_Halo_Mass_Set_Standard(parentNode,Tree_Node_Hot_Halo_Mass_Standard(parentNode) &
               &+Tree_Node_Hot_Halo_Mass_Standard(thisNode))
          call Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard(parentNode,Tree_Node_Hot_Halo_Angular_Momentum_Standard(parentNode) &
               &+Tree_Node_Hot_Halo_Angular_Momentum_Standard(thisNode))
          call Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard(parentNode,Tree_Node_Hot_Halo_Outflowed_Mass_Standard(parentNode) &
               &+Tree_Node_Hot_Halo_Outflowed_Mass_Standard(thisNode))
          call Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard(parentNode &
               &,Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(parentNode) &
               &+Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(thisNode))
          call Tree_Node_Hot_Halo_Mass_Set_Standard             (thisNode,0.0d0)
          call Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard (thisNode,0.0d0)
          call Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard   (thisNode,0.0d0)
          call Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard(thisNode,0.0d0)
          call Tree_Node_Hot_Halo_Abundances_Standard(thisNode  ,abundancesWork  )
          call Tree_Node_Hot_Halo_Abundances_Standard(parentNode,abundancesParent)
          abundancesParent=abundancesParent+abundancesWork
          abundancesWork  =0.0d0
          call Tree_Node_Hot_Halo_Abundances_Set_Standard(thisNode  ,abundancesWork  )
          call Tree_Node_Hot_Halo_Abundances_Set_Standard(parentNode,abundancesParent)
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode  ,abundancesWork      )
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(parentNode,abundancesParent)
          abundancesParent=abundancesParent+abundancesWork
          abundancesWork  =0.0d0
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard(thisNode  ,abundancesWork  )
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard(parentNode,abundancesParent)
          ! Hot reservoir molcules.
          call Tree_Node_Hot_Halo_Molecules_Standard               (thisNode  , molecules                                                )
          call Tree_Node_Hot_Halo_Molecules_Standard               (parentNode, moleculesParent                                          )
          moleculesParent =moleculesParent +molecules
          molecules       =0.0d0
          call Tree_Node_Hot_Halo_Molecules_Set_Standard           (thisNode  , molecules                                                )
          call Tree_Node_Hot_Halo_Molecules_Set_Standard           (parentNode, moleculesParent                                          )
       end if
    end if
    return
  end subroutine Hot_Halo_Starve

  !# <satelliteMergerTask>
  !#  <unitName>Hot_Halo_Remove_Before_Satellite_Merging</unitName>
  !# </satelliteMergerTask>
  subroutine Hot_Halo_Remove_Before_Satellite_Merging(thisNode)
    !% Remove any hot halo associated with {\tt thisNode} before it merges with its host halo.
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: hostNode
    double precision, dimension(moleculesCount ) :: molecules ,moleculesHost

    ! Determine if starvation is to be applied.
    if (methodSelected.and..not.starveSatellites.and.thisNode%componentExists(componentIndex)) then
       
       ! Find the node to merge with.
       call thisNode%mergesWith(hostNode)
       
       ! Move the hot halo to the host.
       call Tree_Node_Hot_Halo_Mass_Set_Standard(hostNode,Tree_Node_Hot_Halo_Mass_Standard(hostNode) &
            &+Tree_Node_Hot_Halo_Mass_Standard(thisNode))
       call Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard(hostNode,Tree_Node_Hot_Halo_Angular_Momentum_Standard(hostNode) &
            &+Tree_Node_Hot_Halo_Angular_Momentum_Standard(thisNode))
       call Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard(hostNode,Tree_Node_Hot_Halo_Outflowed_Mass_Standard(hostNode) &
            &+Tree_Node_Hot_Halo_Outflowed_Mass_Standard(thisNode))
       call Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard(hostNode&
            &,Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(hostNode) &
            &+Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(thisNode))
       call Tree_Node_Hot_Halo_Mass_Set_Standard(thisNode,0.0d0)
       call Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard(thisNode,0.0d0)
       call Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard(thisNode,0.0d0)
       call Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard(thisNode,0.0d0)
       call Tree_Node_Hot_Halo_Abundances_Standard(thisNode,abundancesWork)
       call Tree_Node_Hot_Halo_Abundances_Standard(hostNode,abundancesHost)
       abundancesHost=abundancesHost+abundancesWork
       abundancesWork=0.0d0
       call Tree_Node_Hot_Halo_Abundances_Set_Standard(thisNode,abundancesWork)
       call Tree_Node_Hot_Halo_Abundances_Set_Standard(hostNode,abundancesHost)
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode,abundancesWork)
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(hostNode,abundancesHost)
       abundancesHost=abundancesHost+abundancesWork
       abundancesWork=0.0d0
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard(thisNode,abundancesWork)
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard(hostNode,abundancesHost)
       ! Hot reservoir molecules.
       call Tree_Node_Hot_Halo_Molecules_Standard               (thisNode, molecules                                              )
       call Tree_Node_Hot_Halo_Molecules_Standard               (hostNode, moleculesHost                                          )
       moleculesHost =moleculesHost +molecules
       molecules     =0.0d0
       call Tree_Node_Hot_Halo_Molecules_Set_Standard           (thisNode, molecules                                              )
       call Tree_Node_Hot_Halo_Molecules_Set_Standard           (hostNode, moleculesHost                                          )
    end if
    return
  end subroutine Hot_Halo_Remove_Before_Satellite_Merging

  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Hot_Halo_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Hot_Halo_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the hot halo mass of {\tt
    !% thisNode} to account for any hot halo already in the parent.
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: parentNode
    double precision, dimension(moleculesCount ) :: molecules ,moleculesParent
    double precision                             :: hotHaloMass,angularMomentum

    ! Check if this method is selected.
    if (methodSelected) then
       ! Get the parent node of this node.
       parentNode => thisNode%parentNode
       ! If the parent node has a hot halo component, then add it to that of this node.
       if (parentNode%componentExists(componentIndex)) then
          hotHaloMass=Tree_Node_Hot_Halo_Mass_Standard(thisNode)+Tree_Node_Hot_Halo_Mass_Standard(parentNode)
          call Tree_Node_Hot_Halo_Mass_Set_Standard(thisNode,hotHaloMass)
          angularMomentum=Tree_Node_Hot_Halo_Angular_Momentum_Standard(thisNode)+Tree_Node_Hot_Halo_Angular_Momentum_Standard(parentNode)
          call Tree_Node_Hot_Halo_Angular_Momentum_Set_Standard(thisNode,angularMomentum)
          hotHaloMass=Tree_Node_Hot_Halo_Outflowed_Mass_Standard(thisNode)+Tree_Node_Hot_Halo_Outflowed_Mass_Standard(parentNode)
          call Tree_Node_Hot_Halo_Outflowed_Mass_Set_Standard(thisNode,hotHaloMass)
          angularMomentum=Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(thisNode)&
               &+Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(parentNode)
          call Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Set_Standard(thisNode,angularMomentum)
          call Tree_Node_Hot_Halo_Abundances_Standard(thisNode  ,abundancesWork  )
          call Tree_Node_Hot_Halo_Abundances_Standard(parentNode,abundancesParent)
          abundancesWork=abundancesParent+abundancesWork
          call Tree_Node_Hot_Halo_Abundances_Set_Standard(thisNode,abundancesWork)
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode  ,abundancesWork  )
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(parentNode,abundancesParent)
          abundancesWork=abundancesParent+abundancesWork
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Set_Standard(thisNode,abundancesWork)
          ! Hot reservoir molecules.
          call Tree_Node_Hot_Halo_Molecules_Standard               (thisNode  ,molecules       )
          call Tree_Node_Hot_Halo_Molecules_Standard               (parentNode,moleculesParent )
          molecules      =moleculesParent +molecules
          call Tree_Node_Hot_Halo_Molecules_Set_Standard           (thisNode  ,molecules       )
       end if
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Promote

  subroutine Get_Cooling_Rate(thisNode)
    !% Get and store the cooling rate for {\tt thisNode}.
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    if (.not.gotCoolingRate) then
       if (Tree_Node_Hot_Halo_Mass_Standard(thisNode) > 0.0d0) then
          ! Get the cooling time.
          coolingRate=Cooling_Rate(thisNode)
       else
          coolingRate=0.0d0
       end if

       ! Store a copy of this cooling rate as the remaining mass heating rate budget. This is used to ensure that we never heat
       ! gas at a rate greater than it is cooling.
       massHeatingRateRemaining=coolingRate
       
       ! Flag that cooling rate has now been computed.
       gotCoolingRate=.true.
    end if
    return
  end subroutine Get_Cooling_Rate

  integer function Tree_Node_Hot_Halo_Index(thisNode)
    !% Ensure the hot halo component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Hot_Halo_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Hot_Halo_Index

  subroutine Hot_Halo_Create(thisNode)
    !% Creates a hot halo component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message

    ! Display a message.
    message='Creating hot halo component for node '
    message=message//thisNode%index()
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Create the component.
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    return
  end subroutine Hot_Halo_Create

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Hot_Halo_Standard_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hot_Halo_Standard</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Hot_Halo_Standard_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of hot halo properties to be written to the \glc\ output file.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Abundances_Structure
    use ISO_Varying_String
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    integer                                       :: iAbundance
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='hotHaloMass'
       doublePropertyComments(doubleProperty)='Mass of gas in the hot halo.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='hotHaloUnaccretedMass'
       doublePropertyComments(doubleProperty)='Mass of gas that failed to accrete into the hot halo.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='hotHaloAngularMomentum'
       doublePropertyComments(doubleProperty)='Angular momentum of gas in the hot halo.'
       doublePropertyUnitsSI (doubleProperty)=massSolar*megaParsec*kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='outflowedMass'
       doublePropertyComments(doubleProperty)='Mass of outflowed gas in the hot halo.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='outflowedAngularMomentum'
       doublePropertyComments(doubleProperty)='Angular momentum of outflowed gas in the hot halo.'
       doublePropertyUnitsSI (doubleProperty)=massSolar*megaParsec*kilo
       do iAbundance=1,abundancesCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='hotHalo'//Abundances_Names(iAbundance)
          doublePropertyComments(doubleProperty)='Hot halo abundance property.'
          doublePropertyUnitsSI (doubleProperty)=massSolar
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='outflowed'//Abundances_Names(iAbundance)
          doublePropertyComments(doubleProperty)='Outflowed gas abundance property.'
          doublePropertyUnitsSI (doubleProperty)=massSolar
       end do
     end if
    return
  end subroutine Galacticus_Output_Tree_Hot_Halo_Standard_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Hot_Halo_Standard_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hot_Halo_Standard</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Hot_Halo_Standard_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of hot halo properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+outflowedAbundancesIndexEnd
    return
  end subroutine Galacticus_Output_Tree_Hot_Halo_Standard_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Hot_Halo_Standard</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hot_Halo_Standard</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Hot_Halo_Standard(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store hot halo properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Kind_Numbers
    implicit none
    double precision,        intent(in)                 :: time
    type(treeNode),          intent(inout), pointer     :: thisNode
    integer,                 intent(inout)              :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)              :: integerBuffer(:,:)
    double precision,        intent(inout)              :: doubleBuffer(:,:)
    double precision,        dimension(abundancesCount) :: hotAbundanceMasses,outflowedAbundanceMasses
    integer                                             :: iAbundance

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Hot_Halo_Mass             (thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Hot_Halo_Unaccreted_Mass  (thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Hot_Halo_Angular_Momentum (thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Hot_Halo_Outflowed_Mass   (thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Hot_Halo_Outflowed_Ang_Mom(thisNode)
       call Tree_Node_Hot_Halo_Abundances_Standard          (thisNode,hotAbundanceMasses      )
       call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode,outflowedAbundanceMasses)
       do iAbundance=1,abundancesCount
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=hotAbundanceMasses      (iAbundance)
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=outflowedAbundanceMasses(iAbundance)
       end do
   end if
    return
  end subroutine Galacticus_Output_Tree_Hot_Halo_Standard

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Hot_Halo_Standard_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Hot_Halo_Standard_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    use Abundances_Structure
    use ISO_Varying_String
    implicit none
    type(treeNode),   intent(inout), pointer     :: thisNode
    double precision, dimension(abundancesCount) :: hotAbundanceMasses,outflowedAbundanceMasses
    integer                                      :: iAbundance

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'hot halo component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'hot halo mass:'                      ,Tree_Node_Hot_Halo_Mass_Standard             (thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'hot halo unaccreted mass:'           ,Tree_Node_Hot_Halo_Unaccreted_Mass_Standard  (thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'hot halo outflowed mass:'            ,Tree_Node_Hot_Halo_Outflowed_Mass_Standard   (thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'hot halo angular momentum:'          ,Tree_Node_Hot_Halo_Angular_Momentum_Standard (thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'hot halo outflowed angular momentum:',Tree_Node_Hot_Halo_Outflowed_Ang_Mom_Standard(thisNode)
          call Tree_Node_Hot_Halo_Abundances_Standard          (thisNode,hotAbundanceMasses      )
          call Tree_Node_Hot_Halo_Outflowed_Abundances_Standard(thisNode,outflowedAbundanceMasses)
          do iAbundance=1,abundancesCount
             write (0,'(2x,a50,1x,e12.6)') 'hot halo '          //char(Abundances_Names(iAbundance))//':',hotAbundanceMasses      (iAbundance)
             write (0,'(2x,a50,1x,e12.6)') 'hot halo outflowed '//char(Abundances_Names(iAbundance))//':',outflowedAbundanceMasses(iAbundance)
          end do
       else
          write (0,'(1x,a)'           ) 'hot halo component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Hot_Halo_Standard_Dump

end module Tree_Node_Methods_Hot_Halo
