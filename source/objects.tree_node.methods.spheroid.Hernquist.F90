!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module of Hernquist-profile spheroid tree node methods.

module Tree_Node_Methods_Hernquist_Spheroid
  !% Implement Hernquist spheroid tree node methods.
  use Tree_Nodes
  use Histories
  use Components
  use Stellar_Population_Properties
  private
  public :: Tree_Node_Methods_Hernquist_Spheroid_Initialize, Hernquist_Spheroid_Satellite_Merging,&
       & Galacticus_Output_Tree_Spheroid_Hernquist, Galacticus_Output_Tree_Spheroid_Hernquist_Property_Count,&
       & Galacticus_Output_Tree_Spheroid_Hernquist_Names, Hernquist_Spheroid_Radius_Solver, Hernquist_Spheroid_Enclosed_Mass,&
       & Hernquist_Spheroid_Density, Hernquist_Spheroid_Rotation_Curve, Tree_Node_Spheroid_Post_Evolve_Hernquist,&
       & Tree_Node_Methods_Hernquist_Spheroid_Dump, Hernquist_Spheroid_Radius_Solver_Plausibility
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Internal count of abundances.
  integer :: abundancesCount

  ! Internal count of luminosities.
  integer :: luminositiesCount

  ! Property indices.
  integer, parameter :: propertyCountBase=3, dataCountBase=2, historyCount=1
  integer            :: propertyCount      , dataCount
  integer, parameter :: angularMomentumIndex=1
  integer, parameter :: gasMassIndex        =2
  integer, parameter :: stellarMassIndex    =3
  integer            :: gasAbundancesIndex      ,gasAbundancesIndexEnd
  integer            :: stellarAbundancesIndex  ,stellarAbundancesIndexEnd
  integer            :: stellarLuminositiesIndex,stellarLuminositiesIndexEnd
  integer, parameter :: radiusIndex          =1
  integer, parameter :: velocityIndex        =2
  integer, parameter :: stellarHistoryIndex  =1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spheroid_Radius</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spheroid_Half_Mass_Radius</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spheroid_Velocity</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spheroid_Gas_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spheroid_Angular_Momentum</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spheroid_Stellar_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Spheroid_Gas_Abundances</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Spheroid_Stellar_Abundances</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Spheroid_Stellar_Luminosities</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spheroid_SFR</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="history">
  !#  <methodName>Tree_Node_Spheroid_Stellar_Properties_History</methodName>
  !# </treeNodeMethodsPointer>

  ! Define pipes.
  !
  ! Pointer to procedure which handles generic mass sinks in the spheroid.
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Spheroid_Gas_Sink</pipeName>
  !# </treeNodePipePointer>
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Spheroid_Gas_Energy_Input</pipeName>
  !# </treeNodePipePointer>

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

  ! Parameters controlling the physical implementation.
  double precision :: spheroidEnergeticOutflowMassRate
  double precision :: spheroidMassToleranceAbsolute

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Hernquist_Spheroid_Initialize</unitName>
  !#  <optionName default="Hernquist">treeNodeMethodSpheroid</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Hernquist_Spheroid_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node Hernquist spheroid methods module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    use Galacticus_Error
    use Stellar_Population_Properties_Luminosities
    use Abundances_Structure
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption.eq.'Hernquist') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Hernquist spheroid method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Get number of abundance properties.
       abundancesCount=Abundances_Property_Count()

       ! Get number of luminosity properties.
       luminositiesCount=Stellar_Population_Luminosities_Count()

       ! Determine number of properties needed, including those for stars etc.
       propertyCount=propertyCountBase+2*abundancesCount+luminositiesCount
       dataCount    =dataCountBase

       ! Assign indices to properties.
       gasAbundancesIndex         =propertyCountBase                          +1
       gasAbundancesIndexEnd      =gasAbundancesIndex       +abundancesCount  -1
       stellarAbundancesIndex     =gasAbundancesIndexEnd                      +1
       stellarAbundancesIndexEnd  =stellarAbundancesIndex   +abundancesCount  -1
       stellarLuminositiesIndex   =stellarAbundancesIndexEnd                  +1
       stellarLuminositiesIndexEnd=stellarLuminositiesIndex +luminositiesCount-1

       ! Set up procedure pointers.
       Tree_Node_Spheroid_Radius                                  => Hernquist_Spheroid_Radius
       Tree_Node_Spheroid_Radius_Set                              => null()
       Tree_Node_Spheroid_Radius_Rate_Adjust                      => null()
       Tree_Node_Spheroid_Radius_Rate_Compute                     => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Half_Mass_Radius                        => Hernquist_Spheroid_Half_Mass_Radius
       Tree_Node_Spheroid_Half_Mass_Radius_Set                    => null()
       Tree_Node_Spheroid_Half_Mass_Radius_Rate_Adjust            => null()
       Tree_Node_Spheroid_Half_Mass_Radius_Rate_Compute           => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Velocity                                => Hernquist_Spheroid_Velocity
       Tree_Node_Spheroid_Velocity_Set                            => null()
       Tree_Node_Spheroid_Velocity_Rate_Adjust                    => null()
       Tree_Node_Spheroid_Velocity_Rate_Compute                   => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Angular_Momentum                        => Tree_Node_Spheroid_Angular_Momentum_Hernquist
       Tree_Node_Spheroid_Angular_Momentum_Set                    => Tree_Node_Spheroid_Angular_Momentum_Set_Hernquist
       Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust            => Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Angular_Momentum_Rate_Compute           => Tree_Node_Spheroid_Angular_Momentum_Rate_Compute_Hernquist

       Tree_Node_Spheroid_Gas_Mass                                => Tree_Node_Spheroid_Gas_Mass_Hernquist
       Tree_Node_Spheroid_Gas_Mass_Set                            => Tree_Node_Spheroid_Gas_Mass_Set_Hernquist
       Tree_Node_Spheroid_Gas_Mass_Rate_Adjust                    => Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Gas_Mass_Rate_Compute                   => Tree_Node_Spheroid_Gas_Mass_Rate_Compute_Hernquist

       Tree_Node_Spheroid_Stellar_Mass                            => Tree_Node_Spheroid_Stellar_Mass_Hernquist
       Tree_Node_Spheroid_Stellar_Mass_Set                        => Tree_Node_Spheroid_Stellar_Mass_Set_Hernquist
       Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust                => Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Stellar_Mass_Rate_Compute               => Tree_Node_Spheroid_Stellar_Mass_Rate_Compute_Hernquist

       Tree_Node_Spheroid_Stellar_Abundances                      => Tree_Node_Spheroid_Stellar_Abundances_Hernquist
       Tree_Node_Spheroid_Stellar_Abundances_Set                  => Tree_Node_Spheroid_Stellar_Abundances_Set_Hernquist
       Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust          => Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Stellar_Abundances_Rate_Compute         => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Gas_Abundances                          => Tree_Node_Spheroid_Gas_Abundances_Hernquist
       Tree_Node_Spheroid_Gas_Abundances_Set                      => Tree_Node_Spheroid_Gas_Abundances_Set_Hernquist
       Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust              => Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Gas_Abundances_Rate_Compute             => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Stellar_Luminosities                    => Tree_Node_Spheroid_Stellar_Luminosities_Hernquist
       Tree_Node_Spheroid_Stellar_Luminosities_Set                => Tree_Node_Spheroid_Stellar_Luminosities_Set_Hernquist
       Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust        => Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Stellar_Luminosities_Rate_Compute       => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_SFR                                     => Hernquist_Spheroid_SFR
       Tree_Node_Spheroid_SFR_Set                                 => null()
       Tree_Node_Spheroid_SFR_Rate_Adjust                         => null()
       Tree_Node_Spheroid_SFR_Rate_Compute                        => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Stellar_Properties_History              => Tree_Node_Spheroid_Stellar_Properties_History_Hernquist
       Tree_Node_Spheroid_Stellar_Properties_History_Set          => Tree_Node_Spheroid_Stellar_Properties_History_Set_Hernquist
       Tree_Node_Spheroid_Stellar_Properties_History_Rate_Adjust  => null()
       Tree_Node_Spheroid_Stellar_Properties_History_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       ! Associate pipes with procedures.
       Tree_Node_Spheroid_Gas_Sink                                => Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Gas_Energy_Input                        => Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Hernquist

       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>spheroidEnergeticOutflowMassRate</name>
       !@   <defaultValue>1.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The proportionallity factor relating mass outflow rate from the spheroid to the energy input rate divided by $V_{\rm spheroid}^2$.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidEnergeticOutflowMassRate',spheroidEnergeticOutflowMassRate,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>spheroidMassToleranceAbsolute</name>
       !@   <defaultValue>$10^{-6} M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The mass tolerance used to judge whether the spheroid is physically plausible.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidMassToleranceAbsolute',spheroidMassToleranceAbsolute,defaultValue=1.0d-6)

    end if
    return
  end subroutine Tree_Node_Methods_Hernquist_Spheroid_Initialize
  
  
  !# <postEvolveTask>
  !# <unitName>Tree_Node_Spheroid_Post_Evolve_Hernquist</unitName>
  !# </postEvolveTask>
  subroutine Tree_Node_Spheroid_Post_Evolve_Hernquist(thisNode)
    !% Trim histories attached to the spheroid.
    use Histories
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    ! Trim the stellar populations properties future history.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       call thisNode%components(thisIndex)%histories(stellarHistoryIndex)%trim(Tree_Node_Time(thisNode))
    end if
    return
  end subroutine Tree_Node_Spheroid_Post_Evolve_Hernquist

  double precision function Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
    !% Return the node Hernquist spheroid gas mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Gas_Mass_Hernquist=thisNode%components(thisIndex)%properties(gasMassIndex,propertyValue)
    else
       Tree_Node_Spheroid_Gas_Mass_Hernquist=0.0d0
    end if
    return
  end function Tree_Node_Spheroid_Gas_Mass_Hernquist

  subroutine Tree_Node_Spheroid_Gas_Mass_Set_Hernquist(thisNode,mass)
    !% Set the node Hernquist spheroid gas mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(gasMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Spheroid_Gas_Mass_Set_Hernquist

  subroutine Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Account for a sink of gaseous material in the Hernquist spheroid.
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    double precision                         :: gasMass,stellarMass
    
    ! If no Hernquist spheroid component currently exists and we have some sink from then there is a problem.
    ! spheroid. If sink has zero rate, just return instead.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) call Galacticus_Error_Report('Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Hernquist','mass sink from non-existant Hernquist spheroid')
       return
    end if
    ! Trap cases where an attempt is made to add gas via this sink function.
    if (rateAdjustment > 0.0d0) call Galacticus_Error_Report('Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Hernquist','attempt to add mass via sink in Hernquist spheroid')
    ! Return if no adjustment is being made.
    if (rateAdjustment == 0.0d0) return
    ! Get the index of the component.
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)    
    ! Get the gas mass present.
    gasMass=thisNode%components(thisIndex)%properties(gasMassIndex,propertyValue)
    ! Get the stellar mass present.
    stellarMass=thisNode%components(thisIndex)%properties(stellarMassIndex,propertyValue)
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0) then
       thisNode%components(thisIndex)%properties(gasMassIndex,propertyDerivative)&
            &=thisNode%components(thisIndex)%properties(gasMassIndex,propertyDerivative)+rateAdjustment
       thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyDerivative)&
            &=thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyDerivative) &
            & +rateAdjustment/(gasMass+stellarMass) &
            & *thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyValue)
       thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
            &=thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative) & 
            & +(rateAdjustment/gasMass) &
            & *thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)
    end if
    return
  end subroutine Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Hernquist

  subroutine Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Handles input of energy into the spheroid gas from other components (e.g. black holes). The energy input rate should be in
    !% units of $M_\odot$ km$^2$ s$^{-2}$ Gyr$^{-1}$.
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer, intent(inout)              :: thisNode
    logical,                   intent(inout)              :: interrupt
    procedure(),      pointer, intent(inout)              :: interruptProcedure
    double precision,          intent(in)                 :: rateAdjustment
    double precision,          dimension(abundancesCount) :: abundancesOutflowRate
    integer                                               :: thisIndex
    double precision                                      :: gasMass,stellarMass,massOutflowRate,angularMomentumOutflowRate&
         &,spheroidVelocity

    ! If no Hernquist spheroid component currently exists then energy input to it has no effect.
    if (.not.thisNode%componentExists(componentIndex)) return

    ! Trap cases where an attempt is made to remove energy via this input function.
    if (rateAdjustment < 0.0d0) call Galacticus_Error_Report('Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Hernquist' &
         & ,'attempt to remove energy via input pipe in Hernquist spheroid')
    
    ! Return if no adjustment is being made.
    if (rateAdjustment == 0.0d0) return
    ! Get the index of the component.
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)    
    ! Get the gas mass present.
    gasMass    =thisNode%components(thisIndex)%properties(gasMassIndex    ,propertyValue)
    ! Get the stellar mass present.
    stellarMass=thisNode%components(thisIndex)%properties(stellarMassIndex,propertyValue)
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0) then
       ! Compute outflow rates of quantities and adjust rates in the spheroid appropriately.
       spheroidVelocity=Hernquist_Spheroid_Velocity(thisNode)
       if (spheroidVelocity > 0.0d0) then
          massOutflowRate=spheroidEnergeticOutflowMassRate*rateAdjustment/spheroidVelocity**2
          thisNode%components(thisIndex)%properties(gasMassIndex,propertyDerivative) &
               &=thisNode%components(thisIndex)%properties(gasMassIndex,propertyDerivative)-massOutflowRate
          angularMomentumOutflowRate=massOutflowRate/(gasMass+stellarMass) &
               & *thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyValue)
          thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyDerivative)&
               &=thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyDerivative) &
               & -angularMomentumOutflowRate
          abundancesOutflowRate=(massOutflowRate/gasMass) &
               & *thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)
          thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
               &=thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative) & 
               & -abundancesOutflowRate
          ! Add outflowing rates to the hot halo component.
          call Tree_Node_Hot_Halo_Outflow_Mass_To            (thisNode,interrupt,interruptProcedure,massOutflowRate           )
          call Tree_Node_Hot_Halo_Outflow_Angular_Momentum_To(thisNode,interrupt,interruptProcedure,angularMomentumOutflowRate)
          call Tree_Node_Hot_Halo_Outflow_Abundances_To      (thisNode,interrupt,interruptProcedure,abundancesOutflowRate     )
       end if
    end if
    return
  end subroutine Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Hernquist

  subroutine Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node Hernquist spheroid gas mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    ! If no Hernquist spheroid component currently exists and we have some cooling into it then interrupt and create an Hernquist
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Hernquist_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(gasMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(gasMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Hernquist

  subroutine Tree_Node_Spheroid_Gas_Mass_Rate_Compute_Hernquist(thisNode,interrupt,interruptProcedure)
    !% Compute the Hernquist spheroid node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure

    ! Do nothing here - work will be done elsewhere.
    return
  end subroutine Tree_Node_Spheroid_Gas_Mass_Rate_Compute_Hernquist

  double precision function Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
    !% Return the node Hernquist spheroid stellar mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Stellar_Mass_Hernquist=thisNode%components(thisIndex)%properties(stellarMassIndex,propertyValue)
    else
       Tree_Node_Spheroid_Stellar_Mass_Hernquist=0.0d0
    end if
    return
  end function Tree_Node_Spheroid_Stellar_Mass_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Mass_Set_Hernquist(thisNode,mass)
    !% Set the node Hernquist spheroid stellar mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(stellarMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Spheroid_Stellar_Mass_Set_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node Hernquist spheroid stellar mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    ! If no Hernquist spheroid component currently exists and we have some mass flow into it then interrupt and create a Hernquist
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Hernquist_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(stellarMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(stellarMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Mass_Rate_Compute_Hernquist(thisNode,interrupt,interruptProcedure)
    !% Compute the Hernquist spheroid node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    use Star_Formation_Feedback_Spheroids
    use Abundances_Structure
    implicit none
    type(treeNode),            pointer, intent(inout)       :: thisNode
    logical,                            intent(inout)       :: interrupt
    procedure(),               pointer, intent(inout)       :: interruptProcedure
    double precision,          dimension(luminositiesCount) :: stellarLuminositiesRates
    double precision,          dimension(abundancesCount)   :: abundanceMasses,abundancesOutflowRate
    double precision,          parameter                    :: starFormationRateMinimum=1.0d-10
    integer                                                 :: thisIndex
    double precision                                        :: starFormationRate,stellarMassRate ,fuelMassRate,fuelMass&
         &,massOutflowRate,spheroidMass,angularMomentumOutflowRate,energyInputRate,gasMass
    type(abundancesStructure), save                         :: fuelAbundances,stellarAbundancesRates,fuelAbundancesRates
    !$omp threadprivate(fuelAbundances,stellarAbundancesRates,fuelAbundancesRates)

    ! Compute star formation rate if this node has a spheroid.
    if (thisNode%componentExists(componentIndex)) then
       ! Find the star formation timescale.
       starFormationRate=Hernquist_Spheroid_SFR(thisNode)

       ! If rate is above minimum, then compute related rates. We ignore tiny rates (which can occur as gas is depleted towards
       ! zero) as they can result in divergences in metallicity etc.
       if (starFormationRate > starFormationRateMinimum) then
          ! Get the available fuel mass.
          fuelMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
          
          ! Find the metallicity of the fuel supply.
          call Tree_Node_Spheroid_Gas_Abundances_Hernquist(thisNode,abundanceMasses)
          abundanceMasses=abundanceMasses/fuelMass
          call fuelAbundances%pack(abundanceMasses)
          
          ! Get the component index.
          thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)

          ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
          call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,thisNode&
               &,thisNode%components(thisIndex)%histories(stellarHistoryIndex),stellarMassRate,stellarAbundancesRates &
               &,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
          
          ! Adjust rates.
          call Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Hernquist        (thisNode,interrupt,interruptProcedure&
               &,stellarMassRate)
          call Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Hernquist            (thisNode,interrupt,interruptProcedure,fuelMassRate &
               &  )
          call stellarAbundancesRates%unpack(abundanceMasses)
          call Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Hernquist  (thisNode,interrupt,interruptProcedure&
               &,abundanceMasses)
          call fuelAbundancesRates%unpack(abundanceMasses)
          call Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Hernquist      (thisNode,interrupt,interruptProcedure&
               &,abundanceMasses)
          call Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure&
               &,stellarLuminositiesRates)

          ! Find rate of outflow of material from the spheroid and pipe it to the outflowed reservoir.
          massOutflowRate=Star_Formation_Feedback_Spheroid_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
          if (massOutflowRate > 0.0d0) then
             gasMass     =Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
             spheroidMass=gasMass+Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
             angularMomentumOutflowRate=massOutflowRate*Tree_Node_Spheroid_Angular_Momentum_Hernquist(thisNode)/spheroidMass
             call Tree_Node_Spheroid_Gas_Abundances_Hernquist(thisNode,abundanceMasses)
             abundancesOutflowRate=massOutflowRate*abundanceMasses/gasMass
             call Tree_Node_Hot_Halo_Outflow_Mass_To                     (thisNode,interrupt,interruptProcedure,&
                  & massOutflowRate           )
             call Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Hernquist        (thisNode,interrupt,interruptProcedure,&
                  &-massOutflowRate           )
             call Tree_Node_Hot_Halo_Outflow_Angular_Momentum_To         (thisNode,interrupt,interruptProcedure,&
                  & angularMomentumOutflowRate)
             call Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,&
                  &-angularMomentumOutflowRate)
             call Tree_Node_Hot_Halo_Outflow_Abundances_To               (thisNode,interrupt,interruptProcedure,&
                  & abundancesOutflowRate)
             call Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Hernquist  (thisNode,interrupt,interruptProcedure,&
                  &-abundancesOutflowRate)
          end if
                
       end if
    end if
    return
  end subroutine Tree_Node_Spheroid_Stellar_Mass_Rate_Compute_Hernquist

  subroutine Tree_Node_Spheroid_Gas_Abundances_Hernquist(thisNode,abundanceMasses)
    !% Return the node Hernquist spheroid gas abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: abundanceMasses(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       abundanceMasses(:)=thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd&
            &,propertyValue)
    else
       abundanceMasses(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Spheroid_Gas_Abundances_Hernquist

  subroutine Tree_Node_Spheroid_Gas_Abundances_Set_Hernquist(thisNode,abundanceMasses)
    !% Set the node Hernquist spheroid gas abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: abundanceMasses(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)=abundanceMasses(:)
    return
  end subroutine Tree_Node_Spheroid_Gas_Abundances_Set_Hernquist

  subroutine Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node Hernquist spheroid gas abundance masses rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no Hernquist spheroid component currently exists and we have some cooling into it then interrupt and create an Hernquist
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Hernquist_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Abundances_Hernquist(thisNode,abundanceMasses)
    !% Return the node Hernquist spheroid stellar abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: abundanceMasses(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       abundanceMasses(:)=thisNode%components(thisIndex)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd&
            &,propertyValue)
    else
       abundanceMasses(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Spheroid_Stellar_Abundances_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Abundances_Set_Hernquist(thisNode,abundanceMasses)
    !% Set the node Hernquist spheroid stellar abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: abundanceMasses(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyValue)=abundanceMasses(:)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Abundances_Set_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node Hernquist spheroid stellar abundance masses rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no Hernquist spheroid component currently exists and we have some cooling into it then interrupt and create an Hernquist
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Hernquist_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Luminosities_Hernquist(thisNode,luminosities)
    !% Return the node Hernquist spheroid stellar luminosities.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: luminosities(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       luminosities(:)=thisNode%components(thisIndex)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd&
            &,propertyValue)
    else
       luminosities(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Spheroid_Stellar_Luminosities_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Luminosities_Set_Hernquist(thisNode,luminosities)
    !% Set the node Hernquist spheroid stellar luminosities.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: luminosities(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyValue)=luminosities(:)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Luminosities_Set_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node Hernquist spheroid stellar luminosity rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no Hernquist spheroid component currently exists and we have some cooling into it then interrupt and create an Hernquist
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Hernquist_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Hernquist

  double precision function Tree_Node_Spheroid_Angular_Momentum_Hernquist(thisNode)
    !% Return the node Hernquist spheroid gas mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Angular_Momentum_Hernquist=thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyValue)
    else
       Tree_Node_Spheroid_Angular_Momentum_Hernquist=0.0d0
    end if
    return
  end function Tree_Node_Spheroid_Angular_Momentum_Hernquist

  subroutine Tree_Node_Spheroid_Angular_Momentum_Set_Hernquist(thisNode,angularMomentum)
    !% Set the node Hernquist spheroid gas mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: angularMomentum
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyValue)=angularMomentum
    return
  end subroutine Tree_Node_Spheroid_Angular_Momentum_Set_Hernquist

  subroutine Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node Hernquist spheroid gas mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    ! If no Hernquist spheroid component currently exists and we have some mass flow into it then interrupt and create a Hernquist
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Hernquist_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyDerivative) &
         &=thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Hernquist

  subroutine Tree_Node_Spheroid_Angular_Momentum_Rate_Compute_Hernquist(thisNode,interrupt,interruptProcedure)
    !% Compute the Hernquist spheroid node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure

    ! No need to do anything here - the work will be done elsewhere.
    return
  end subroutine Tree_Node_Spheroid_Angular_Momentum_Rate_Compute_Hernquist

  function Tree_Node_Spheroid_Stellar_Properties_History_Hernquist(thisNode)
    !% Return spheroid stellar properties history.
    use Histories
    implicit none
    type(history)                          :: Tree_Node_Spheroid_Stellar_Properties_History_Hernquist
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Stellar_Properties_History_Hernquist=thisNode%components(thisIndex)%histories(stellarHistoryIndex)
    else
       Tree_Node_Spheroid_Stellar_Properties_History_Hernquist=nullHistory
    end if
    return
  end function Tree_Node_Spheroid_Stellar_Properties_History_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Properties_History_Set_Hernquist(thisNode,thisHistory)
    !% Set the node spheroid stellar properties history.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(history),           intent(in)    :: thisHistory
    integer                                :: thisIndex

    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    call thisNode%components(thisIndex)%histories(stellarHistoryIndex)%destroy()
    thisNode%components(thisIndex)%histories(stellarHistoryIndex)=thisHistory
    return
  end subroutine Tree_Node_Spheroid_Stellar_Properties_History_Set_Hernquist


  !# <satelliteMergerTask>
  !#  <unitName>Hernquist_Spheroid_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Hernquist_Spheroid_Satellite_Merging(thisNode)
    !% Transfer any Hernquist spheroid associated with {\tt thisNode} to its host halo.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    use Satellite_Merging_Remnant_Sizes_Properties
    implicit none
    type(treeNode),   pointer, intent(inout)       :: thisNode
    type(treeNode),   pointer                      :: hostNode
    double precision, dimension(abundancesCount)   :: thisAbundances,hostAbundances
    double precision, dimension(luminositiesCount) :: thisLuminosities,hostLuminosities
    double precision, parameter                    :: spheroidAngularMomentumRatio=0.5595552243d0 ! Ratio of specific angular
                                                                                                  ! momentum at half-mass radius to
                                                                                                  ! global mean.
    double precision                               :: spheroidSpecificAngularMomentum,diskSpecificAngularMomentum,angularMomentum&
         &,spheroidMass
    type(history)                                  :: historyDisk,historySpheroid,thisHistory

    ! Check that method is selected.
    if (methodSelected) then
       
       ! Find the host node.
       hostNode => thisNode%parentNode

       ! Get specific angular momentum of the host spheroid and disk material.
       if (Tree_Node_Spheroid_Gas_Mass_Hernquist(hostNode)+Tree_Node_Spheroid_Stellar_Mass_Hernquist(hostNode) > 0.0d0) then
          spheroidSpecificAngularMomentum=Tree_Node_Spheroid_Angular_Momentum_Hernquist(hostNode)&
               &/(Tree_Node_Spheroid_Gas_Mass_Hernquist(hostNode)+Tree_Node_Spheroid_Stellar_Mass_Hernquist(hostNode))
       else
          spheroidSpecificAngularMomentum=0.0d0
       end if

       if (Tree_Node_Disk_Gas_Mass(hostNode)+Tree_Node_Disk_Stellar_Mass(hostNode) > 0.0d0) then
          diskSpecificAngularMomentum=Tree_Node_Disk_Angular_Momentum(hostNode)/(Tree_Node_Disk_Gas_Mass(hostNode)+Tree_Node_Disk_Stellar_Mass(hostNode))
       else
          diskSpecificAngularMomentum=0.0d0
       end if

       ! Move gas material within the host if necessary.
       select case (thisHostGasMovesTo)
       case (movesToDisk)
          call Tree_Node_Disk_Gas_Mass_Set                    (hostNode, Tree_Node_Disk_Gas_Mass              (hostNode) &
               &                                                        +Tree_Node_Spheroid_Gas_Mass_Hernquist(hostNode))
          call Tree_Node_Disk_Gas_Abundances                  (hostNode, hostAbundances                                 )
          call Tree_Node_Spheroid_Gas_Abundances_Hernquist    (hostNode, thisAbundances                                 )
          call Tree_Node_Disk_Gas_Abundances_Set              (hostNode, hostAbundances+thisAbundances                  )
          call Tree_Node_Disk_Angular_Momentum_Set            (hostNode, Tree_Node_Disk_Angular_Momentum(hostNode)       &
               &+spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Gas_Mass_Hernquist(hostNode)                        )
          call Tree_Node_Spheroid_Angular_Momentum_Set        (hostNode, Tree_Node_Spheroid_Angular_Momentum(hostNode)   &
               &-spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Gas_Mass_Hernquist(hostNode)                        )
          call Tree_Node_Spheroid_Gas_Mass_Set_Hernquist      (hostNode, 0.0d0                                          )
          thisAbundances=0.0d0
          call Tree_Node_Spheroid_Gas_Abundances_Set_Hernquist(hostNode, thisAbundances                                 )
       case (movesToSpheroid)
          call Tree_Node_Spheroid_Gas_Mass_Set_Hernquist      (hostNode, Tree_Node_Spheroid_Gas_Mass_Hernquist(hostNode) &
               &                                                        +Tree_Node_Disk_Gas_Mass              (hostNode))
          call Tree_Node_Disk_Angular_Momentum_Set            (hostNode, Tree_Node_Disk_Angular_Momentum(hostNode)       &
               &-diskSpecificAngularMomentum*Tree_Node_Disk_Gas_Mass(hostNode)                                          )
          call Tree_Node_Spheroid_Gas_Abundances_Hernquist    (hostNode,hostAbundances                                  )
          call Tree_Node_Disk_Gas_Abundances                  (hostNode,thisAbundances                                  )
          call Tree_Node_Spheroid_Gas_Abundances_Set_Hernquist(hostNode,hostAbundances+thisAbundances                   )
          call Tree_Node_Disk_Gas_Mass_Set                    (hostNode,0.0d0                                           )
          thisAbundances=0.0d0
          call Tree_Node_Disk_Gas_Abundances_Set              (hostNode,thisAbundances                                  )
       case (doesNotMove)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('Hernquist_Spheroid_Satellite_Merging','unrecognized movesTo descriptor')
       end select

       ! Move stellar material within the host if necessary.
       select case (thisHostStarsMoveTo)
       case (movesToDisk)
          call Tree_Node_Disk_Stellar_Mass_Set                      (hostNode, Tree_Node_Disk_Stellar_Mass              (hostNode) &
               &                                                              +Tree_Node_Spheroid_Stellar_Mass_Hernquist(hostNode))
          call Tree_Node_Disk_Stellar_Abundances                    (hostNode, hostAbundances                                     )
          call Tree_Node_Spheroid_Stellar_Abundances_Hernquist      (hostNode, thisAbundances                                     )
          call Tree_Node_Disk_Stellar_Abundances_Set                (hostNode, hostAbundances+thisAbundances                      )
          call Tree_Node_Disk_Stellar_Luminosities                  (hostNode, hostLuminosities                                   )
          call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist    (hostNode, thisLuminosities                                   )
          call Tree_Node_Disk_Stellar_Luminosities_Set              (hostNode, hostLuminosities+thisLuminosities                  )
          call Tree_Node_Disk_Angular_Momentum_Set                  (hostNode, Tree_Node_Disk_Angular_Momentum(hostNode)           &
               &+spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Stellar_Mass_Hernquist(hostNode)                              )
          call Tree_Node_Spheroid_Angular_Momentum_Set              (hostNode, Tree_Node_Spheroid_Angular_Momentum(hostNode)        &
               &-spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Stellar_Mass_Hernquist(hostNode)                              )
          call Tree_Node_Spheroid_Stellar_Mass_Set_Hernquist        (hostNode, 0.0d0                                              )
          thisAbundances=0.0d0
          call Tree_Node_Spheroid_Stellar_Abundances_Set_Hernquist  (hostNode, thisAbundances                                     )
          thisLuminosities=0.0d0
          call Tree_Node_Spheroid_Stellar_Luminosities_Set_Hernquist(hostNode, thisLuminosities                                   )
          ! Also add stellar properties histories.
          historyDisk    =Tree_Node_Disk_Stellar_Properties_History              (hostNode                )
          historySpheroid=Tree_Node_Spheroid_Stellar_Properties_History_Hernquist(hostNode                )
          call historyDisk%add(historySpheroid)
          call Tree_Node_Disk_Stellar_Properties_History_Set                     (hostNode,historyDisk    )
          call historySpheroid%reset()
          call Tree_Node_Spheroid_Stellar_Properties_History_Set_Hernquist       (hostNode,historySpheroid)
       case (movesToSpheroid)
          call Tree_Node_Spheroid_Stellar_Mass_Set_Hernquist        (hostNode, Tree_Node_Spheroid_Stellar_Mass_Hernquist(hostNode) &
               &                                                              +Tree_Node_Disk_Stellar_Mass              (hostNode))
          call Tree_Node_Disk_Angular_Momentum_Set                  (hostNode, Tree_Node_Disk_Angular_Momentum(hostNode)           &
               &-diskSpecificAngularMomentum*Tree_Node_Disk_Stellar_Mass(hostNode)                                                )
          call Tree_Node_Spheroid_Stellar_Abundances_Hernquist      (hostNode, hostAbundances                                     )
          call Tree_Node_Disk_Stellar_Abundances                    (hostNode, thisAbundances                                     )
          call Tree_Node_Spheroid_Stellar_Abundances_Set_Hernquist  (hostNode, hostAbundances+thisAbundances                      )
          call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist    (hostNode, hostLuminosities                                   )
          call Tree_Node_Disk_Stellar_Luminosities                  (hostNode, thisLuminosities                                   )
          call Tree_Node_Spheroid_Stellar_Luminosities_Set_Hernquist(hostNode, hostLuminosities+thisLuminosities                  )
          call Tree_Node_Disk_Stellar_Mass_Set                      (hostNode, 0.0d0                                              )
          thisAbundances=0.0d0
          call Tree_Node_Disk_Stellar_Abundances_Set                (hostNode, thisAbundances                                     )
          thisLuminosities=0.0d0
          call Tree_Node_Disk_Stellar_Luminosities_Set              (hostNode, thisLuminosities                                   )
          ! Also add stellar properties histories.
          historyDisk    =Tree_Node_Disk_Stellar_Properties_History              (hostNode                )
          historySpheroid=Tree_Node_Spheroid_Stellar_Properties_History_Hernquist(hostNode                )
          call historySpheroid%add(historyDisk)
          call Tree_Node_Spheroid_Stellar_Properties_History_Set_Hernquist       (hostNode,historySpheroid)
          call historyDisk%reset()
          call Tree_Node_Disk_Stellar_Properties_History_Set                     (hostNode,historyDisk    )
       case (doesNotMove)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('Hernquist_Spheroid_Satellite_Merging','unrecognized movesTo descriptor')
       end select
   
       ! Get specific angular momentum of the spheroid material.
       spheroidMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
       if (spheroidMass > 0.0d0) then
          spheroidSpecificAngularMomentum=Tree_Node_Spheroid_Angular_Momentum_Hernquist(thisNode)/spheroidMass

          ! Move the gas component of the Hernquist spheroid to the host.
          select case (thisMergerGasMovesTo)
          case (movesToDisk)
             call Tree_Node_Disk_Gas_Mass_Set                (hostNode, Tree_Node_Disk_Gas_Mass              (hostNode) &
                  &                                                    +Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode))
             call Tree_Node_Disk_Gas_Abundances              (hostNode,hostAbundances)
             call Tree_Node_Spheroid_Gas_Abundances_Hernquist(thisNode,thisAbundances)
             call Tree_Node_Disk_Gas_Abundances_Set          (hostNode,hostAbundances+thisAbundances)
             call Tree_Node_Disk_Angular_Momentum_Set        (hostNode,Tree_Node_Disk_Angular_Momentum(hostNode) &
                  &+spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode))
          case (movesToSpheroid)
             call Tree_Node_Spheroid_Gas_Mass_Set_Hernquist  (hostNode, Tree_Node_Spheroid_Gas_Mass_Hernquist(hostNode) &
                  &                                                    +Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode))
             call Tree_Node_Spheroid_Gas_Abundances_Hernquist(hostNode,hostAbundances)
             call Tree_Node_Spheroid_Gas_Abundances_Hernquist(thisNode,thisAbundances)
             call Tree_Node_Spheroid_Gas_Abundances_Set_Hernquist(hostNode,hostAbundances+thisAbundances)
          case default
             call Galacticus_Error_Report('Hernquist_Spheroid_Satellite_Merging','unrecognized movesTo descriptor')
          end select
          call Tree_Node_Spheroid_Gas_Mass_Set_Hernquist      (thisNode,0.0d0         )
          thisAbundances=0.0d0
          call Tree_Node_Spheroid_Gas_Abundances_Set_Hernquist(thisNode,thisAbundances)
          
          ! Move the stellar component of the Hernquist spheroid to the host.
          select case (thisMergerStarsMoveTo)
          case (movesToDisk)
             call Tree_Node_Disk_Stellar_Mass_Set                  (hostNode, Tree_Node_Disk_Stellar_Mass(hostNode) &
                  &                                                          +Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode))
             call Tree_Node_Disk_Stellar_Abundances                (hostNode,hostAbundances)
             call Tree_Node_Spheroid_Stellar_Abundances_Hernquist  (thisNode,thisAbundances)
             call Tree_Node_Disk_Stellar_Abundances_Set            (hostNode,hostAbundances+thisAbundances)
             call Tree_Node_Disk_Stellar_Luminosities              (hostNode,hostLuminosities)
             call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist(thisNode,thisLuminosities)
             call Tree_Node_Disk_Stellar_Luminosities_Set          (hostNode,hostLuminosities+thisLuminosities)
             call Tree_Node_Disk_Angular_Momentum_Set              (hostNode,Tree_Node_Disk_Angular_Momentum(hostNode) &
                  &+spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode))
             ! Also add stellar properties histories.
             historySpheroid=Tree_Node_Spheroid_Stellar_Properties_History_Hernquist(thisNode                )
             thisHistory    =Tree_Node_Disk_Stellar_Properties_History              (hostNode                )
             call thisHistory%add(historySpheroid)
             call Tree_Node_Disk_Stellar_Properties_History_Set                     (hostNode,thisHistory    )
             call historySpheroid%reset()
             call Tree_Node_Spheroid_Stellar_Properties_History_Set_Hernquist       (thisNode,historySpheroid)
         case (movesToSpheroid)
             call Tree_Node_Spheroid_Stellar_Mass_Set_Hernquist        (hostNode, Tree_Node_Spheroid_Stellar_Mass_Hernquist(hostNode) &
                  &                                                              +Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode))
             call Tree_Node_Spheroid_Stellar_Abundances_Hernquist      (hostNode,hostAbundances)
             call Tree_Node_Spheroid_Stellar_Abundances_Hernquist      (thisNode,thisAbundances)
             call Tree_Node_Spheroid_Stellar_Abundances_Set_Hernquist  (hostNode,hostAbundances+thisAbundances)
             call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist    (hostNode,hostLuminosities)
             call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist    (thisNode,thisLuminosities)
             call Tree_Node_Spheroid_Stellar_Luminosities_Set_Hernquist(hostNode,hostLuminosities+thisLuminosities)
             ! Also add stellar properties histories.
             historySpheroid=Tree_Node_Spheroid_Stellar_Properties_History_Hernquist(thisNode                )
             thisHistory    =Tree_Node_Spheroid_Stellar_Properties_History_Hernquist(hostNode                )
             call thisHistory%add(historySpheroid)
             call Tree_Node_Spheroid_Stellar_Properties_History_Set_Hernquist       (hostNode,thisHistory    )
             call historySpheroid%reset()
             call Tree_Node_Spheroid_Stellar_Properties_History_Set_Hernquist       (thisNode,historySpheroid)
          case default
             call Galacticus_Error_Report('Hernquist_Spheroid_Satellite_Merging','unrecognized movesTo descriptor')
          end select
          call Tree_Node_Spheroid_Stellar_Mass_Set_Hernquist        (thisNode,0.0d0           )
          thisAbundances=0.0d0
          call Tree_Node_Spheroid_Stellar_Abundances_Set_Hernquist  (thisNode,thisAbundances  )
          thisLuminosities=0.0d0
          call Tree_Node_Spheroid_Stellar_Luminosities_Set_Hernquist(thisNode,thisLuminosities)       
          call Tree_Node_Spheroid_Angular_Momentum_Set_Hernquist    (thisNode,0.0d0           )
       end if
       
       ! Set the angular momentum of the spheroid.
       if (remnantSpecificAngularMomentum /= remnantNoChangeValue) then
          angularMomentum=(remnantSpecificAngularMomentum/spheroidAngularMomentumRatio)*(Tree_Node_Spheroid_Gas_Mass(hostNode)&
               &+Tree_Node_Spheroid_Stellar_Mass(hostNode))
          call Tree_Node_Spheroid_Angular_Momentum_Set_Hernquist(hostNode,angularMomentum)
       end if

    end if
    return
  end subroutine Hernquist_Spheroid_Satellite_Merging
  
  !# <enclosedMassTask>
  !#  <unitName>Hernquist_Spheroid_Enclosed_Mass</unitName>
  !# </enclosedMassTask>
  subroutine Hernquist_Spheroid_Enclosed_Mass(thisNode,radius,massType,componentType,componentMass)
    !% Computes the mass within a given radius for an Hernquist spheroid.
    use Galactic_Structure_Options
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass
    double precision                         :: fractionalRadius,spheroidRadius
    
    componentMass=0.0d0
    if (.not.methodSelected                           ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
    if (.not.thisNode%componentExists(componentIndex)) return
    select case (massType)
    case (massTypeAll,massTypeBaryonic,massTypeGalactic)
       componentMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
    case (massTypeGaseous)
       componentMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
    case (massTypeStellar)
       componentMass=Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
    end select
    ! Return if total mass was requested.  
    if (radius >= radiusLarge)                          return
    ! Return if mass is zero.
    if (componentMass <= 0.0d0)                         return
    ! Compute actual mass.
    spheroidRadius=Hernquist_Spheroid_Radius(thisNode)
    if (spheroidRadius > 0.0d0) then
       fractionalRadius=radius/spheroidRadius
       componentMass=componentMass*(fractionalRadius/(fractionalRadius+1.0d0))**2
    end if
    return
  end subroutine Hernquist_Spheroid_Enclosed_Mass

  !# <rotationCurveTask>
  !#  <unitName>Hernquist_Spheroid_Rotation_Curve</unitName>
  !# </rotationCurveTask>
  subroutine Hernquist_Spheroid_Rotation_Curve(thisNode,radius,massType,componentType,componentVelocity)
    !% Computes the rotation curve at a given radius for an Hernquist spheroid.
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentVelocity
    double precision                         :: componentMass
    
    componentVelocity=0.0d0
    if (radius > 0.0d0) then
       call Hernquist_Spheroid_Enclosed_Mass(thisNode,radius,massType,componentType,componentMass)
       if (componentMass >= 0.0d0) componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass)/dsqrt(radius)
    end if
    return
  end subroutine Hernquist_Spheroid_Rotation_Curve

  !# <densityTask>
  !#  <unitName>Hernquist_Spheroid_Density</unitName>
  !# </densityTask>
  subroutine Hernquist_Spheroid_Density(thisNode,positionSpherical,massType,componentType,componentDensity)
    !% Computes the density at a given position for an Hernquist spheroid.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: positionSpherical(3)
    double precision, intent(out)            :: componentDensity
    double precision                         :: fractionalRadius
    
    componentDensity=0.0d0
    if (.not.methodSelected                           ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
    if (.not.thisNode%componentExists(componentIndex)) return
    if (Hernquist_Spheroid_Radius(thisNode) <= 0.0d0) return
    select case (massType)
    case (massTypeAll,massTypeBaryonic,massTypeGalactic)
       componentDensity=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
    case (massTypeGaseous)
       componentDensity=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
    case (massTypeStellar)
       componentDensity=Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
    end select
    ! Return if density is zero.
    if (componentDensity <= 0.0d0) return
    ! Compute actual density.
    fractionalRadius=positionSpherical(1)/Hernquist_Spheroid_Radius(thisNode)
    componentDensity=componentDensity/2.0d0/Pi/fractionalRadius/(Hernquist_Spheroid_Radius(thisNode)*(1.0d0+fractionalRadius))**3
    return
  end subroutine Hernquist_Spheroid_Density

  !# <radiusSolverPlausibility>
  !#  <unitName>Hernquist_Spheroid_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Hernquist_Spheroid_Radius_Solver_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
    !% Determines whether the spheroid is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: galaxyIsPhysicallyPlausible

    ! Return immediately if our method is not selected.
    if (.not.methodSelected) return

    if (Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)+Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode) < -spheroidMassToleranceAbsolute) then
       galaxyIsPhysicallyPlausible=.false.
    else
       if (Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)+Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode) >&
            & spheroidMassToleranceAbsolute .and. Tree_Node_Spheroid_Angular_Momentum_Hernquist(thisNode) < 0.0d0)&
            & galaxyIsPhysicallyPlausible=.false.
    end if

    return
  end subroutine Hernquist_Spheroid_Radius_Solver_Plausibility

  !# <radiusSolverTask>
  !#  <unitName>Hernquist_Spheroid_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Hernquist_Spheroid_Radius_Solver(thisNode,componentActive,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !% Interface for the size solver algorithm.
    implicit none
    type(treeNode),          pointer, intent(inout) :: thisNode
    logical,                          intent(out)   :: componentActive
    double precision,                 intent(out)   :: specificAngularMomentum
    procedure(),             pointer, intent(out)   :: Radius_Get,Velocity_Get
    procedure(),             pointer, intent(out)   :: Radius_Set,Velocity_Set
    double precision,        parameter              :: spheroidAngularMomentumRatio=0.2546479089d0 ! Ratio of specific angular
                                                                                                   ! momentum at scale radius to
                                                                                                   ! global mean.
    double precision                                :: specificAngularMomentumMean,angularMomentum,spheroidMass

    ! Determine if thisNode has an active spheroid component supported by this module.    
    componentActive=methodSelected

    if (methodSelected) componentActive=thisNode%componentExists(componentIndex)    
    if (componentActive) then
       ! Get the angular momentum.
       angularMomentum=Tree_Node_Spheroid_Angular_Momentum_Hernquist(thisNode)
       if (angularMomentum >= 0.0d0) then
          ! Compute the specific angular momentum at the scale radius, assuming a flat rotation curve.
          spheroidMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
          if (spheroidMass > 0.0d0) then
             specificAngularMomentumMean=angularMomentum/spheroidMass
          else
             specificAngularMomentumMean=0.0d0
          end if
          specificAngularMomentum=spheroidAngularMomentumRatio*specificAngularMomentumMean
          ! Associate the pointers with the appropriate property routines.
          Radius_Get   => Hernquist_Spheroid_Radius
          Radius_Set   => Hernquist_Spheroid_Radius_Set
          Velocity_Get => Hernquist_Spheroid_Velocity
          Velocity_Set => Hernquist_Spheroid_Velocity_Set
      else
          call Hernquist_Spheroid_Radius_Set  (thisNode,0.0d0)
          call Hernquist_Spheroid_Velocity_Set(thisNode,0.0d0)
          componentActive=.false.
       end if
    end if
    return
  end subroutine Hernquist_Spheroid_Radius_Solver


  double precision function Hernquist_Spheroid_SFR(thisNode)
    !% Return the star formation rate of the Hernquist spheroid.
    use Star_Formation_Timescales_Spheroids
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex
    double precision                         :: starFormationTimescale,gasMass

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       
       ! Get the star formation timescale.
       starFormationTimescale=Star_Formation_Timescale_Spheroid(thisNode)
       
       ! Get the gas mass.
       gasMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)

       ! If timescale is finite and gas mass is positive, then compute star formation rate.
       if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0) then
          Hernquist_Spheroid_SFR=gasMass/starFormationTimescale
       else
          Hernquist_Spheroid_SFR=0.0d0
       end if
    else
       Hernquist_Spheroid_SFR=0.0d0
    end if
    return
  end function Hernquist_Spheroid_SFR

  double precision function Hernquist_Spheroid_Radius(thisNode)
    !% Return the scale radius of the Hernquist spheroid.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       Hernquist_Spheroid_Radius=thisNode%components(thisIndex)%data(radiusIndex)
    else
       Hernquist_Spheroid_Radius=0.0d0
    end if
    return
  end function Hernquist_Spheroid_Radius

  double precision function Hernquist_Spheroid_Half_Mass_Radius(thisNode)
    !% Return the scale radius of the Hernquist spheroid.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: halfMassRadiusToScaleRadius=1.0d0/(dsqrt(2.0d0)-1.0d0)

    Hernquist_Spheroid_Half_Mass_Radius=Hernquist_Spheroid_Radius(thisNode)*halfMassRadiusToScaleRadius
    return
  end function Hernquist_Spheroid_Half_Mass_Radius

  subroutine Hernquist_Spheroid_Radius_Set(thisNode,radius)
    !% Set the scale radius of the Hernquist spheroid.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: radius
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       thisNode%components(thisIndex)%data(radiusIndex)=max(radius,0.0d0)
    end if
    return
  end subroutine Hernquist_Spheroid_Radius_Set

  double precision function Hernquist_Spheroid_Velocity(thisNode)
    !% Return the circular velocity of the Hernquist spheroid.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       Hernquist_Spheroid_Velocity=thisNode%components(thisIndex)%data(velocityIndex)
    else
       Hernquist_Spheroid_Velocity=0.0d0
    end if
    return
  end function Hernquist_Spheroid_Velocity

  subroutine Hernquist_Spheroid_Velocity_Set(thisNode,velocity)
    !% Set the circular velocity of the Hernquist spheroid.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: velocity
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       thisNode%components(thisIndex)%data(velocityIndex)=velocity
    end if
    return
  end subroutine Hernquist_Spheroid_Velocity_Set

  integer function Tree_Node_Hernquist_Spheroid_Index(thisNode)
    !% Ensure the Hernquist spheroid component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (.not.thisNode%componentExists(componentIndex)) then
       ! Create the component.
       call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
       ! Get the index for this component.
       thisIndex=thisNode%componentIndex(componentIndex)
       ! Create the stellar properties history.
       call Stellar_Population_Properties_History_Create(thisNode,thisNode%components(thisIndex)%histories(stellarHistoryIndex))
    else
       ! Get the index for this component.
       thisIndex=thisNode%componentIndex(componentIndex)
    end if
    Tree_Node_Hernquist_Spheroid_Index=thisIndex
    return
  end function Tree_Node_Hernquist_Spheroid_Index

  subroutine Hernquist_Spheroid_Create(thisNode)
    !% Creates an Hernquist spheroid component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message
    integer                                     :: thisIndex

    ! Display a message.
    message='Creating Hernquist spheroid component for node '
    message=message//thisNode%index()
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Get the index of the component (which will also ensure that the component is created).
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    return
  end subroutine Hernquist_Spheroid_Create

  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Spheroid_Hernquist_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spheroid_Hernquist</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Spheroid_Hernquist_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of Hernquist spheroid properties to be written to the \glc\ output file.
    use Abundances_Structure
    use ISO_Varying_String
    use Stellar_Population_Properties_Luminosities
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    integer                                       :: iAbundance,iLuminosity

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidGasMass'
       doublePropertyComments(doubleProperty)='Mass of gas in the Hernquist spheroid.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidStellarMass'
       doublePropertyComments(doubleProperty)='Mass of stars in the Hernquist spheroid at scale length.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidAngularMomentum'
       doublePropertyComments(doubleProperty)='Angular momentum of the Hernquist spheroid.'
       doublePropertyUnitsSI (doubleProperty)=massSolar*megaParsec*kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidScaleLength'
       doublePropertyComments(doubleProperty)='Radial scale length in the Hernquist spheroid.'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidCircularVelocity'
       doublePropertyComments(doubleProperty)='Circular velocity of the Hernquist spheroid at scale length.'
       doublePropertyUnitsSI (doubleProperty)=kilo
       do iAbundance=1,abundancesCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='spheroidGas'//Abundances_Names(iAbundance)
          doublePropertyComments(doubleProperty)='Spheroid gas phase abundance property.'
          doublePropertyUnitsSI (doubleProperty)=massSolar
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='spheroidStellar'//Abundances_Names(iAbundance)
          doublePropertyComments(doubleProperty)='Spheroid stellar abundance property.'
          doublePropertyUnitsSI (doubleProperty)=massSolar
       end do
       do iLuminosity=1,luminositiesCount
          if (Stellar_Population_Luminosities_Output(iLuminosity,time)) then
             doubleProperty=doubleProperty+1
             doublePropertyNames   (doubleProperty)='spheroidStellar'//Stellar_Population_Luminosities_Name(iLuminosity)
             doublePropertyComments(doubleProperty)='Spheroid stellar luminosity property.'
             doublePropertyUnitsSI (doubleProperty)=luminosityZeroPointAB
          end if
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Spheroid_Hernquist_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Spheroid_Hernquist_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spheroid_Hernquist</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Spheroid_Hernquist_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of Hernquist spheroid properties to be written to the the \glc\ output file.
    use Stellar_Population_Properties_Luminosities
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+propertyCount+dataCount-luminositiesCount&
         &+Stellar_Population_Luminosities_Output_Count(time)

    return
  end subroutine Galacticus_Output_Tree_Spheroid_Hernquist_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Spheroid_Hernquist</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spheroid_Hernquist</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Spheroid_Hernquist(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store Hernquist spheroid properties in the \glc\ output file buffers.
    use Stellar_Population_Properties_Luminosities
    use Tree_Nodes
    implicit none
    double precision, intent(in)                   :: time
    type(treeNode),   intent(inout), pointer       :: thisNode
    integer,          intent(inout)                :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer,          intent(inout)                :: integerBuffer(:,:)
    double precision, intent(inout)                :: doubleBuffer(:,:)
    double precision, dimension(abundancesCount)   :: gasAbundanceMasses,stellarAbundanceMasses
    double precision, dimension(luminositiesCount) :: stellarLuminosities
    integer                                        :: iAbundance,iLuminosity

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spheroid_Angular_Momentum_Hernquist(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Hernquist_Spheroid_Radius(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Hernquist_Spheroid_Velocity(thisNode)
       call Tree_Node_Spheroid_Gas_Abundances_Hernquist    (thisNode,gasAbundanceMasses)
       call Tree_Node_Spheroid_Stellar_Abundances_Hernquist(thisNode,stellarAbundanceMasses)
       do iAbundance=1,abundancesCount
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=gasAbundanceMasses(iAbundance)
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=stellarAbundanceMasses(iAbundance)
       end do
       call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist(thisNode,stellarLuminosities)
       do iLuminosity=1,luminositiesCount
          if (Stellar_Population_Luminosities_Output(iLuminosity,time)) then
             doubleProperty=doubleProperty+1
             doubleBuffer(doubleBufferCount,doubleProperty)=stellarLuminosities(iLuminosity)
          end if
       end do
      end if
    return
  end subroutine Galacticus_Output_Tree_Spheroid_Hernquist

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Hernquist_Spheroid_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Hernquist_Spheroid_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    use Abundances_Structure
    use ISO_Varying_String
    use Stellar_Population_Properties_Luminosities
    implicit none
    type(treeNode),   intent(inout), pointer       :: thisNode
    double precision, dimension(abundancesCount)   :: gasAbundanceMasses,stellarAbundanceMasses
    double precision, dimension(luminositiesCount) :: stellarLuminosities
    integer                                        :: iAbundance,iLuminosity

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)') 'spheroid component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'spheroid gas mass:',Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'spheroid stellar mass:',Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'spheroid angular momentum:',Tree_Node_Spheroid_Angular_Momentum_Hernquist(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'spheroid scale length:',Hernquist_Spheroid_Radius(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'spheroid circular velocity:',Hernquist_Spheroid_Velocity(thisNode)
          call Tree_Node_Spheroid_Gas_Abundances_Hernquist    (thisNode,gasAbundanceMasses    )
          call Tree_Node_Spheroid_Stellar_Abundances_Hernquist(thisNode,stellarAbundanceMasses)
          do iAbundance=1,abundancesCount
             write (0,'(2x,a50,1x,e12.6)') 'spheroid gas '//char(Abundances_Names(iAbundance))//':',gasAbundanceMasses(iAbundance)
             write (0,'(2x,a50,1x,e12.6)') 'spheroid stellar '//char(Abundances_Names(iAbundance))//':',stellarAbundanceMasses(iAbundance)
          end do
          call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist(thisNode,stellarLuminosities)
          do iLuminosity=1,luminositiesCount
             write (0,'(2x,a50,1x,e12.6)') 'spheroid stellar '//char(Stellar_Population_Luminosities_Name(iLuminosity))//':',stellarLuminosities(iLuminosity)
          end do
       else
          write (0,'(1x,a)') 'spheroid component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Hernquist_Spheroid_Dump

end module Tree_Node_Methods_Hernquist_Spheroid
