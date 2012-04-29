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


!% Contains a module of black hole tree node methods.

module Tree_Node_Methods_Black_Hole
  !% Implement black hole tree node methods.
  use Tree_Nodes
  use Tree_Node_Methods_Black_Hole_Data
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Black_Hole_Initialize, Galacticus_Output_Tree_Black_Hole_Standard,&
       & Galacticus_Output_Tree_Black_Hole_Standard_Property_Count, Galacticus_Output_Tree_Black_Hole_Standard_Names,&
       & Tree_Node_Black_Hole_Reset_Standard, Black_Hole_Satellite_Merging, Tree_Node_Methods_Black_Hole_Standard_Dump,&
       & Black_Hole_Standard_Scale_Set, Black_Hole_Standard_Property_Identifiers_Decode,                              &
       & Galacticus_Output_Tree_Black_Hole_Properties,Galacticus_Output_Tree_Black_Hole_Merger,                       &
       & Black_Hole_Triple_Interaction_Black_Holes
  
  ! Property indices.
  integer, parameter :: propertyCount=3, dataCount=1, historyCount=0
  integer, parameter :: massIndex                 =1
  integer, parameter :: spinIndex                 =2
  integer, parameter :: radiusIndex               =3
  integer, parameter :: tripleInteractionTimeIndex=1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Black_Hole_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Black_Hole_Spin</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Black_Hole_Radial_Position</methodName>
  !# </treeNodeMethodsPointer>

  ! Seed mass for black holes.
  double precision :: blackHoleSeedMass

  ! Seed spin for black holes.
  double precision, parameter :: blackHoleSeedSpin=1.0d-3

  ! Accretion model parameters.
  ! Enhancement factors for the accretion rate.
  double precision :: bondiHoyleAccretionEnhancementSpheroid,bondiHoyleAccretionEnhancementHotHalo
  ! Temperature of accreting gas.
  double precision :: bondiHoyleAccretionTemperatureSpheroid
  ! Control for hot mode only accretion.
  logical          :: bondiHoyleAccretionHotModeOnly

  ! Feedback parameters.
  double precision :: blackHoleWindEfficiency
  logical          :: blackHoleHeatsHotHalo

  ! Quantities stored to avoid repeated computation.
  logical         , allocatable, dimension(:) :: gotAccretionRate
  double precision, allocatable, dimension(:) :: accretionRate,accretionRateHotHalo
  !$omp threadprivate(gotAccretionRate,accretionRate,accretionRateHotHalo)

  ! Output options.
  logical          :: blackHoleOutputAccretion
  logical          :: blackHoleOutputData
  logical          :: blackHoleOutputMergers
  
  ! Option specifying whether the triple black hole interaction should be used.
  logical          :: tripleBlackHoleInteraction

  ! Index of black hole instance about to merge.
  integer          :: mergingInstance
  !$omp threadprivate(mergingInstance)

  ! Index of black hole involved in three-body interactions
  integer          :: binaryInstance,tripleInstance
  !$omp threadprivate(binaryInstance)
  !$omp threadprivate(tripleInstance)


contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Black_Hole_Initialize</unitName>
  !#  <optionName default="standard">treeNodeMethodBlackHole</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Black_Hole_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node black hole methods module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
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
       message='Standard black hole method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Black_Hole_Mass                         => Tree_Node_Black_Hole_Mass_Standard
       Tree_Node_Black_Hole_Mass_Set                     => Tree_Node_Black_Hole_Mass_Set_Standard
       Tree_Node_Black_Hole_Mass_Rate_Adjust             => Tree_Node_Black_Hole_Mass_Rate_Adjust_Standard
       Tree_Node_Black_Hole_Mass_Rate_Compute            => Tree_Node_Black_Hole_Mass_Rate_Compute_Standard
       Tree_Node_Black_Hole_Spin                         => Tree_Node_Black_Hole_Spin_Standard
       Tree_Node_Black_Hole_Spin_Set                     => Tree_Node_Black_Hole_Spin_Set_Standard
       Tree_Node_Black_Hole_Spin_Rate_Adjust             => Tree_Node_Black_Hole_Spin_Rate_Adjust_Standard
       Tree_Node_Black_Hole_Spin_Rate_Compute            => Tree_Node_Black_Hole_Spin_Rate_Compute_Standard
       Tree_Node_Black_Hole_Radial_Position              => Tree_Node_Black_Hole_Radial_Position_Standard
       Tree_Node_Black_Hole_Radial_Position_Set          => Tree_Node_Black_Hole_Radial_Position_Set_Standard
       Tree_Node_Black_Hole_Radial_Position_Rate_Adjust  => Tree_Node_Black_Hole_Radial_Position_Rate_Adjust_Standard
       Tree_Node_Black_Hole_Radial_Position_Rate_Compute => Tree_Node_Black_Hole_Radial_Position_Rate_Compute_Standard

       ! Get the black hole seed mass.
       !@ <inputParameter>
       !@   <name>blackHoleSeedMass</name>
       !@   <defaultValue>100</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass of the seed black hole placed at the center of each newly formed galaxy.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleSeedMass",blackHoleSeedMass,defaultValue=100.0d0)

       ! Get accretion rate enhancement factors.
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionEnhancementSpheroid</name>
       !@   <defaultValue>3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which the Bondi-Hoyle accretion rate of spheroid gas onto black holes in enhanced.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionEnhancementSpheroid",bondiHoyleAccretionEnhancementSpheroid,defaultValue&
            &=3.0d0)
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionEnhancementHotHalo</name>
       !@   <defaultValue>11.4486</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which the Bondi-Hoyle accretion rate of hot halo gas onto black holes in enhanced.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionEnhancementHotHalo",bondiHoyleAccretionEnhancementHotHalo,defaultValue=11.4486d0)
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionHotModeOnly</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether accretion from the hot halo should only occur if the halo is in the hot accretion mode.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionHotModeOnly",bondiHoyleAccretionHotModeOnly,defaultValue=.false.)

       ! Get temperature of accreting gas.
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionTemperatureSpheroid</name>
       !@   <defaultValue>$10^2$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The assumed temperature (in Kelvin) of gas in the spheroid when computing Bondi-Hoyle accretion rates onto black holes.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionTemperatureSpheroid",bondiHoyleAccretionTemperatureSpheroid,defaultValue&
            &=1.0d2)

       ! Get temperature of accreting gas.
       !@ <inputParameter>
       !@   <name>blackHoleWindEfficiency</name>
       !@   <defaultValue>$2.2157\times 10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The efficiency of the black hole-driven wind: $L_{\rm wind} = \epsilon_{\rm wind} \dot{M}_\bullet \clight^2$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleWindEfficiency",blackHoleWindEfficiency,defaultValue=2.2157d-3)

       ! Options controlling AGN feedback.
       !@ <inputParameter>
       !@   <name>blackHoleHeatsHotHalo</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not the black hole launched jets should heat the hot halo.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleHeatsHotHalo",blackHoleHeatsHotHalo,defaultValue=.true.)

       ! Get options controlling output.
       !@ <inputParameter>
       !@   <name>blackHoleOutputAccretion</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not accretion rates and jet powers will be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleOutputAccretion",blackHoleOutputAccretion,defaultValue=.false.)

       ! Get options controlling output.
       !@ <inputParameter>
       !@   <name>blackHoleOutputData</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not properties for all black holes (rather than just the central black hole) will be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleOutputData",blackHoleOutputData,defaultValue=.false.)

       !@ <inputParameter>
       !@   <name>blackHoleOutputMergers</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not properties of black hole mergers will be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleOutputMergers",blackHoleOutputMergers,defaultValue=.false.)

       ! Get options controlling three body interactions.
       !@ <inputParameter>
       !@   <name>tripleBlackHoleInteraction</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not triple black hole interactions will be accounted for.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("tripleBlackHoleInteraction",tripleBlackHoleInteraction,defaultValue=.false.)

    end if
    return
  end subroutine Tree_Node_Methods_Black_Hole_Initialize
  
  !# <calculationResetTask>
  !# <unitName>Tree_Node_Black_Hole_Reset_Standard</unitName>
  !# </calculationResetTask>
  subroutine Tree_Node_Black_Hole_Reset_Standard(thisNode)
    !% Remove memory of stored computed values as we're about to begin computing derivatives anew.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    if (allocated(gotAccretionRate)) gotAccretionRate=.false.
    return
  end subroutine Tree_Node_Black_Hole_Reset_Standard

  double precision function Tree_Node_Black_Hole_Mass_Standard(thisNode,instance)
    !% Return the node black hole mass.
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(in),    optional :: instance
    integer                                 :: thisIndex,thisInstance

    if (thisNode%componentExists(componentIndex)) then
       thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
       thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
       Tree_Node_Black_Hole_Mass_Standard=thisNode%components(thisIndex)%instance(thisInstance)%properties(massIndex,propertyValue)
    else
       Tree_Node_Black_Hole_Mass_Standard=blackHoleSeedMass
    end if
    return
  end function Tree_Node_Black_Hole_Mass_Standard

  subroutine Tree_Node_Black_Hole_Mass_Set_Standard(thisNode,mass,instance)
    !% Set the node black hole mass.
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    double precision, intent(in)              :: mass
    integer,          intent(in),    optional :: instance
    integer                                   :: thisIndex,thisInstance
    
    thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
    thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
    thisNode%components(thisIndex)%instance(thisInstance)%properties(massIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Black_Hole_Mass_Set_Standard

  subroutine Tree_Node_Black_Hole_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node black hole mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    logical,          intent(inout)           :: interrupt
    procedure(),      intent(inout), pointer  :: interruptProcedure
    double precision, intent(in)              :: rateAdjustment
    integer,          intent(in),    optional :: instance
    integer                                   :: thisIndex,thisInstance
    
    thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
    thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
    thisNode%components(thisIndex)%instance(thisInstance)%properties(massIndex,propertyDerivative)=thisNode%components(thisIndex)%instance(thisInstance)%properties(massIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Black_Hole_Mass_Rate_Adjust_Standard

  subroutine Tree_Node_Black_Hole_Mass_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the black hole node mass rate of change.
    use Accretion_Disks
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    procedure(),      pointer                :: interruptProcedurePassed
    double precision, parameter              :: windVelocity=1.0d4   ! Velocity of disk wind.
    double precision, parameter              :: ismTemperature=1.0d4 ! Temperature of the ISM.
    double precision, parameter              :: criticalDensityNormalization=2.0d0*massHydrogenAtom*speedLight**2*megaParsec/3.0d0&
         &/Pi/boltzmannsConstant/gigaYear/ismTemperature/kilo/windVelocity
    integer                                  :: iInstance,instanceCount,thisIndex
    double precision                         :: restMassAccretionRate,massAccretionRate,radiativeEfficiency,energyInputRate &
         &,spheroidDensity,spheroidGasMass,spheroidRadius,criticalDensity,windFraction,spheroidDensityOverCriticalDensity&
         &,heatingRate,couplingEfficiency,jetEfficiency

    ! Get a local copy of the interrupt procedure.
    interruptProcedurePassed => interruptProcedure

    ! Loop over instances.
    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Black_Hole_Index(thisNode)
       instanceCount=size(thisNode%components(thisIndex)%instance)
    else
       instanceCount=0
    end if
    do iInstance=1,max(instanceCount,1)
       if (instanceCount > 0) call thisNode%components(thisIndex)%activeInstanceSet(iInstance)
       
       ! Find the rate of rest mass accretion onto the black hole.
       restMassAccretionRate=Mass_Accretion_Rate(thisNode)
       
       ! Finish if there is no accretion.
       if (restMassAccretionRate <= 0.0d0) cycle
    
       ! Find the radiative efficiency of the accretion.
       radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisNode,restMassAccretionRate)
       
       ! Find the jet efficiency.
       if (restMassAccretionRate > 0.0d0) then
          jetEfficiency=Accretion_Disk_Jet_Power(thisNode,restMassAccretionRate)/restMassAccretionRate/(speedLight/kilo)**2
       else
          jetEfficiency=0.0d0
       end if

       ! Find the rate of increase in mass of the black hole.
       massAccretionRate=restMassAccretionRate*(1.0d0-radiativeEfficiency-jetEfficiency)
       
       ! If no black hole component currently exists and we have some accretion then interrupt and create a black hole.
       if (.not.thisNode%componentExists(componentIndex)) then    
          if (massAccretionRate /= 0.0d0) then
             interrupt=.true.
             interruptProcedure => Black_Hole_Create
          end if
          return
       end if
       call Tree_Node_Black_Hole_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedurePassed, massAccretionRate,instance=iInstance)
    
       ! Remove the accreted mass from the spheroid component.
       call Tree_Node_Spheroid_Gas_Sink                   (thisNode,interrupt,interruptProcedurePassed,-accretionRate       (iInstance))
    
       ! Remove the accreted mass from the hot halo component.
       call Tree_Node_Hot_Halo_Hot_Gas_Sink               (thisNode,interrupt,interruptProcedurePassed,-accretionRateHotHalo(iInstance))
       
       ! Add heating to the hot halo component.
       if (blackHoleHeatsHotHalo) then

          ! Compute jet coupling efficiency based on whether halo is cooling quasistatically. Reduce this efficiency as the gas
          ! content in the halo drops below the cosmological mean.
          couplingEfficiency=Hot_Mode_Fraction(thisNode)*((Omega_Matter()/Omega_b())*Tree_Node_Hot_Halo_Mass(thisNode)/Tree_Node_Mass(thisNode))**2

          ! Get jet power.
          heatingRate=jetEfficiency*restMassAccretionRate*(speedLight/kilo)**2*couplingEfficiency

          ! Pipe this power to the hot halo.
          call Tree_Node_Hot_Halo_Heat_Input(thisNode,interrupt,interruptProcedurePassed,heatingRate)
         
       end if
       
       ! Add energy to the spheroid component.
       if (blackHoleWindEfficiency > 0.0d0) then
          spheroidGasMass=Tree_Node_Spheroid_Gas_Mass(thisNode)
          if (spheroidGasMass > 0.0d0) then
             spheroidRadius=Tree_Node_Spheroid_Radius(thisNode)
             if (spheroidRadius > 0.0d0) then
                spheroidDensity=3.0d0*spheroidGasMass/4.0d0/Pi/spheroidRadius**3
                criticalDensity=criticalDensityNormalization*blackHoleWindEfficiency*restMassAccretionRate/spheroidRadius**2
                
                ! Construct an interpolating factor such that the energy input from the wind drops to zero below half of the critical density.
                spheroidDensityOverCriticalDensity=spheroidDensity/criticalDensity-0.5d0
                if (spheroidDensityOverCriticalDensity <= 0.0d0) then
                   ! No energy input below half of critical density.
                   windFraction=0.0d0
                else if (spheroidDensityOverCriticalDensity >= 1.0d0) then
                   ! Full energy input above 1.5 times critical density.
                   windFraction=1.0d0
                else
                   ! Smooth polynomial interpolating function between these limits.
                   windFraction=3.0d0*spheroidDensityOverCriticalDensity**2-2.0d0*spheroidDensityOverCriticalDensity**3
                end if
                
                ! Compute the energy input and send it down the spheroid gas energy input pipe.
                energyInputRate=windFraction*blackHoleWindEfficiency*restMassAccretionRate*(speedLight/kilo)**2
                call Tree_Node_Spheroid_Gas_Energy_Input(thisNode,interrupt,interruptProcedurePassed,energyInputRate)

             end if
          end if
       end if
    end do
    if (instanceCount > 0) call thisNode%components(thisIndex)%activeInstanceNullify()
    
    ! Return our local copy of the interrupt procedure.
    interruptProcedure => interruptProcedurePassed
    
    return
  end subroutine Tree_Node_Black_Hole_Mass_Rate_Compute_Standard
  
  double precision function Tree_Node_Black_Hole_Spin_Standard(thisNode,instance)
    !% Return the node black hole angular momentum.
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode   
    integer,          intent(in),    optional :: instance
    double precision, parameter               :: blackHoleSpinMaximum=0.9999d0 ! Maximum allowed spin (useful to avoid infinities).
    integer                                   :: thisIndex,thisInstance

    if (thisNode%componentExists(componentIndex)) then
       thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
       thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
       Tree_Node_Black_Hole_Spin_Standard=thisNode%components(thisIndex)%instance(thisInstance)%properties(spinIndex,propertyValue)
       ! Ensure that the spin is within range.
       if (Tree_Node_Black_Hole_Spin_Standard < 0.0d0) then
          Tree_Node_Black_Hole_Spin_Standard=0.0d0
          call Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,0.0d0)
       end if
       if (Tree_Node_Black_Hole_Spin_Standard > blackHoleSpinMaximum) then
          Tree_Node_Black_Hole_Spin_Standard=blackHoleSpinMaximum
          call Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,blackHoleSpinMaximum)
       end if
    else
       Tree_Node_Black_Hole_Spin_Standard=blackHoleSeedSpin
    end if
    return
  end function Tree_Node_Black_Hole_Spin_Standard

  subroutine Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,spin,instance)
    !% Set the node black hole angular momentum.
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    double precision, intent(in)              :: spin
    integer,          intent(in),    optional :: instance
    integer                                   :: thisIndex,thisInstance

    thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
    thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
    thisNode%components(thisIndex)%instance(thisInstance)%properties(spinIndex,propertyValue)=spin
    return
  end subroutine Tree_Node_Black_Hole_Spin_Set_Standard

  subroutine Tree_Node_Black_Hole_Spin_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node black hole mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    logical,          intent(inout)           :: interrupt
    procedure(),      intent(inout), pointer  :: interruptProcedure
    double precision, intent(in)              :: rateAdjustment
    integer,          intent(in),    optional :: instance
    integer                                   :: thisIndex,thisInstance
    
    thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
    thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
    thisNode%components(thisIndex)%instance(thisInstance)%properties(spinIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(thisInstance)%properties(spinIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Black_Hole_Spin_Rate_Adjust_Standard

  subroutine Tree_Node_Black_Hole_Spin_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the black hole node mass rate of change.
    use Accretion_Disks
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    procedure(),    pointer                :: interruptProcedurePassed
    integer                                :: iInstance,instanceCount,thisIndex
    double precision                       :: massAccretionRate,spinUpRate

    ! Loop over instances.
    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Black_Hole_Index(thisNode)
       instanceCount=size(thisNode%components(thisIndex)%instance)
    else
       instanceCount=0
    end if
    do iInstance=1,max(instanceCount,1)
       if (instanceCount > 0) call thisNode%components(thisIndex)%activeInstanceSet(iInstance)
       
       ! Find the rate of gas mass accretion onto the halo.
       massAccretionRate=Mass_Accretion_Rate(thisNode)
       ! If no black hole component currently exists and we have some accretion then interrupt and create a black hole.
       if (.not.thisNode%componentExists(componentIndex)) then    
          if (massAccretionRate /= 0.0d0) then
             interrupt=.true.
             interruptProcedure => Black_Hole_Create
          end if
          return
       end if
       
       if (massAccretionRate > 0.0d0) then
          spinUpRate=Black_Hole_Spin_Up_Rate(thisNode,massAccretionRate)
          call Tree_Node_Black_Hole_Spin_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,spinUpRate,instance=iInstance)
       end if
    end do
     if (instanceCount > 0) call thisNode%components(thisIndex)%activeInstanceNullify()
    return
  end subroutine Tree_Node_Black_Hole_Spin_Rate_Compute_Standard

  double precision function Tree_Node_Black_Hole_Radial_Position_Standard(thisNode,instance)
    !% Return the node black hole radial position.
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(in),    optional :: instance
    integer                                 :: thisIndex,thisInstance

    if (thisNode%componentExists(componentIndex)) then
       thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
       thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
       Tree_Node_Black_Hole_Radial_Position_Standard=thisNode%components(thisIndex)%instance(thisInstance)%properties(radiusIndex,propertyValue)
    else
       Tree_Node_Black_Hole_Radial_Position_Standard=0.0d0
    end if
    return
  end function Tree_Node_Black_Hole_Radial_Position_Standard

  subroutine Tree_Node_Black_Hole_Radial_Position_Set_Standard(thisNode,radius,instance)
    !% Set the node black hole radial position.
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    double precision, intent(in)              :: radius
    integer,          intent(in),    optional :: instance
    integer                                   :: thisIndex,thisInstance

    thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
    thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
    thisNode%components(thisIndex)%instance(thisInstance)%properties(radiusIndex,propertyValue)=radius
    return
  end subroutine Tree_Node_Black_Hole_Radial_Position_Set_Standard

  subroutine Tree_Node_Black_Hole_Radial_Position_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Adjust the rate of evolution of the black hole radial position.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    logical,          intent(inout)           :: interrupt
    procedure(),      intent(inout), pointer  :: interruptProcedure
    double precision, intent(in)              :: rateAdjustment
    integer,          intent(in),    optional :: instance
    integer                                   :: thisIndex,thisInstance
    
    thisInstance=Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
    thisIndex   =Tree_Node_Black_Hole_Index       (thisNode)
    thisNode                                                    &
         &          %components(thisIndex)                      &
         &          %instance(thisInstance)                     &
         &          %properties(radiusIndex,propertyDerivative) &
         & = thisNode                                           &
         &          %components(thisIndex                     ) &
         &          %instance  (thisInstance                  ) &
         &          %properties(radiusIndex,propertyDerivative) &
         &  +rateAdjustment
    return
  end subroutine Tree_Node_Black_Hole_Radial_Position_Rate_Adjust_Standard

  subroutine Tree_Node_Black_Hole_Radial_Position_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the black hole node radial position rate of change.
    use Tree_Nodes
    use Black_Hole_Binary_Separations
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    use Merger_Tree_Active
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    integer                                  :: iInstance,thisIndex
    double precision                         :: radialMigrationRate,radiusHardBinary,binaryRadius
    logical                                  :: binaryRadiusFound

    ! If a black hole exists, compute the rate of change of its radial position.
    if (thisNode%componentExists(componentIndex)) then    
       thisIndex=Tree_Node_Black_Hole_Index(thisNode)
       do iInstance=1,size(thisNode%components(thisIndex)%instance)
          if (iInstance > 1) then
             call thisNode%components(thisIndex)%activeInstanceSet(iInstance)
             ! Compute the hard binary radius.
             radiusHardBinary= (                                                                 &
                  &              gravitationalConstantGalacticus                                 &
                  &             *(                                                               &
                  &                Tree_Node_Black_Hole_Mass       (thisNode,instance=1        ) &
                  &               +Tree_Node_Black_Hole_Mass       (thisNode,instance=iInstance) &
                  &              )                                                               &
                  &            )                                                                 &
                  &           /(                                                                 &
                  &              4.0d0                                                           &
                  &             *  Dark_Matter_Halo_Virial_Velocity(thisNode)**2                 &
                  &            )
             ! Places a new black hole in the center of the galaxy in case there is no central one.
             if (          Tree_Node_Black_Hole_Mass (thisNode,instance=1        ) == 0.0d0               .and. &
               & Tree_Node_Black_Hole_Radial_Position(thisNode,instance=iInstance) <= radiusHardBinary  ) then
                mergingInstance=iInstance
                interrupt=.true.
                interruptProcedure => Black_Hole_Standard_Merge_Black_Holes
                return
             end if
             ! Check for a black hole that is about to merge.
             if (Tree_Node_Black_Hole_Radial_Position(thisNode,instance=iInstance) <= 0.0d0) then    
                ! Record which instance is merging, then trigger an interrupt.
                mergingInstance=iInstance
                interrupt=.true.
                interruptProcedure => Black_Hole_Standard_Merge_Black_Holes
                return
             end if
             ! Set the rate of radial migration.
             radialMigrationRate=Black_Hole_Binary_Separation_Growth_Rate(thisNode)
             call Tree_Node_Black_Hole_Radial_Position_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,radialMigrationRate,instance=iInstance)
          end if
       end do
       call thisNode%components(thisIndex)%activeInstanceNullify()
    end if

    ! Loop over black holes, testing for triple black hole interactions. Find the three closest black holes then check if a three
    ! body interaction occurs using the radial condition derived in Hoffman and Loeb (2007).
    binaryRadiusFound=.false.
    if (thisNode%componentExists(componentIndex) .and. tripleBlackHoleInteraction) then
       thisIndex=Tree_Node_Black_Hole_Index(thisNode)
       if (size(thisNode%components(thisIndex)%instance) >= 3 .and. .not. Tree_Node_Black_Hole_Mass(thisNode,instance=1) <= 0.0d0 ) then
          do iInstance=2,size(thisNode%components(thisIndex)%instance)
             if     (                                                                                             &
                  &  (          Tree_Node_Black_Hole_Radial_Position(thisNode,instance=iInstance) <= binaryRadius &
                  &   .or. .not.binaryRadiusFound                                                                 &
                  &  )                                                                                            &
                  &  .and. .not.Tree_Node_Black_Hole_Mass           (thisNode,instance=iInstance) <= 0.0d0        &
                  &  .and. .not.Tree_Node_Black_Hole_Radial_Position(thisNode,instance=iInstance) <= 0.0d0        &
                  & ) then
                binaryRadius     =Tree_Node_Black_Hole_Radial_Position(thisNode,instance=iInstance)
                binaryInstance   =iInstance
                binaryRadiusFound=.true.
             end if
          end do
          if (binaryRadiusFound) then
             ! Compute the hard binary radius.
             radiusHardBinary= (                                                                      &
                  &              gravitationalConstantGalacticus                                      &
                  &             *(                                                                    &
                  &                Tree_Node_Black_Hole_Mass       (thisNode,instance=1             ) &
                  &               +Tree_Node_Black_Hole_Mass       (thisNode,instance=binaryInstance) &
                  &              )                                                                    &
                  &            )                                                                      &
                  &           /(                                                                      &
                  &              4.0d0                                                                &
                  &             *  Dark_Matter_Halo_Virial_Velocity(thisNode)**2                      &
                  &            )
             ! Search for a third black hole.  
             do iInstance=2,size(thisNode%components(thisIndex)%instance)
                if     (      .not.                                                        iInstance                    == binaryInstance   &
                     &  .and. .not. Tree_Node_Black_Hole_Mass           (thisNode,instance=iInstance)                   <= 0.0d0            &
                     &  .and. .not. Tree_Node_Black_Hole_Radial_Position(thisNode,instance=iInstance)                   <= 0.0d0            &
                     &  .and.       Tree_Node_Black_Hole_Radial_Position(thisNode,instance=iInstance)                   <= radiusHardBinary &
                     &  .and.       thisNode%components(thisIndex)%instance(iInstance)%data(tripleInteractionTimeIndex) == 0                &
                     & ) then
                   tripleInstance=iInstance
                   interrupt=.true.
                   interruptProcedure => Black_Hole_Triple_Interaction_Black_Holes
                   return
                end if
             end do
          end if
       end if
    end if
    return
  end subroutine Tree_Node_Black_Hole_Radial_Position_Rate_Compute_Standard

  integer function Tree_Node_Black_Hole_Get_Instance(thisNode,instance)
    !% Returns the index of the instance to use for get, set and rate adjust functions.
    use Tree_Nodes
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(in),    optional :: instance
    integer                                 :: thisIndex

    if (present(instance)) then
       ! A specific instance was requested, so simply return it.
       Tree_Node_Black_Hole_Get_Instance=instance
    else
       ! No specific instance was requested - find the active instance.
       if (thisNode%componentExists(componentIndex)) then
          thisIndex=Tree_Node_Black_Hole_Index(thisNode)
          Tree_Node_Black_Hole_Get_Instance=thisNode%components(thisIndex)%activeInstance()
          ! If there is no active instance, return the first instance.
          if (Tree_Node_Black_Hole_Get_Instance == instanceNull) Tree_Node_Black_Hole_Get_Instance=1
       else
          Tree_Node_Black_Hole_Get_Instance=1
      end if
    end if
    return
  end function Tree_Node_Black_Hole_Get_Instance

  !# <scaleSetTask>
  !#  <unitName>Black_Hole_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Black_Hole_Standard_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: scaleMassRelative=1.0d-4
    double precision, parameter              :: scaleSizeRelative=1.0d-4
    double precision, parameter              :: scaleSizeAbsolute=1.0d-6
    integer                                  :: thisIndex,iInstance
 
    ! Determine if method is active and a black hole component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Black_Hole_Index(thisNode)
       
       ! Loop over instances.
       do iInstance=1,size(thisNode%components(thisIndex)%instance)
          
          ! Set scale for mass.
          thisNode                                                                                                                         &
               & %components(thisIndex                )                                                                                    &
               & %instance  (iInstance                )                                                                                    &
               & %properties(massIndex  ,propertyScale)=max( scaleMassRelative*Tree_Node_Spheroid_Stellar_Mass(thisNode          )         &
               &                                            ,                  Tree_Node_Black_Hole_Mass      (thisNode,iInstance)         &
               &                                           )
          
          ! Set scale for spin.
          thisNode                                                                                                                         &
               & %components(thisIndex                )                                                                                    &
               & %instance  (iInstance                )                                                                                    &
               & %properties(spinIndex  ,propertyScale)=1.0d0
          
          ! Set scale for radius.
          thisNode&
               & %components(thisIndex                )                                                                                     &
               & %instance  (iInstance                )                                                                                     &
               & %properties(radiusIndex,propertyScale)=maxval(                                                                             &
               &                                               [                                                                            &
               &                                                 scaleSizeAbsolute                                                          &
               &                                                ,scaleSizeRelative*Tree_Node_Spheroid_Half_Mass_Radius (thisNode          ) &
               &                                                ,                  Tree_Node_Black_Hole_Radial_Position(thisNode,iInstance) &
               &                                               ]                                                                            &
               &                                              )
          
       end do       
    end if
    return
  end subroutine Black_Hole_Standard_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Black_Hole_Satellite_Merging</unitName>
  !# </satelliteMergerTask>
  subroutine Black_Hole_Satellite_Merging(thisNode)
    !% Merge (instantaneously) any black hole associated with {\tt thisNode} before it merges with its host halo.
    use Black_Hole_Binary_Mergers
    use Black_Hole_Binary_Initial_Radii
    use Black_Hole_Binary_Recoil_Velocities
    use Galactic_Structure_Potentials
    use Galactic_Structure_Options
    use Components
    implicit none
    type(treeNode),   pointer,     intent(inout) :: thisNode
    type(treeNode),   pointer                    :: hostNode
    type(component),  allocatable, dimension(:)  :: instancesTemporary
    integer                                      :: thisIndex,hostIndex,iInstance,nonNullBlackHoleCount,firstNonNullBlackHole
    double precision                             :: radiusInitial,blackHoleMassNew,blackHoleSpinNew
    double precision                             :: recoilVelocity,massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2

    if (methodSelected) then
       
       ! Find the node to merge with.
       call thisNode%mergesWith(hostNode)
       
       ! Find the initial radius of the satellite black hole in the remnant.
       radiusInitial=Black_Hole_Binary_Initial_Radius(thisNode,hostNode)
       
       ! If the separation is non-positive, assume that the black holes merge instantaneously.
       if (radiusInitial <= 0.0d0) then
          ! Loop over all black holes in the satellite galaxy.
          thisIndex=Tree_Node_Black_Hole_Index(thisNode)
          do iInstance=1,size(thisNode%components(thisIndex)%instance)
             call Black_Hole_Binary_Merger(Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=iInstance), &
                  &                        Tree_Node_Black_Hole_Mass_Standard(hostNode                   ), &
                  &                        Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=iInstance), &
                  &                        Tree_Node_Black_Hole_Spin_Standard(hostNode                   ), &
                  &                        blackHoleMassNew                                               , &
                  &                        blackHoleSpinNew                                                 &
                  &                       )     
             ! Merge the black holes instantaneously.
             ! Check which black hole is more massive in order to compute an appropriate recoil velocity
             if (Tree_Node_Black_Hole_Mass_Standard(hostNode) >= Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=iInstance)) then
                massBlackHole1=Tree_Node_Black_Hole_Mass_Standard(hostNode                   ) 
                massBlackHole2=Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=iInstance)
                spinBlackHole1=Tree_Node_Black_Hole_Spin_Standard(hostNode                   )
                spinBlackHole2=Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=iInstance)
             else
                massBlackHole2=Tree_Node_Black_Hole_Mass_Standard(hostNode                   ) 
                massBlackHole1=Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=iInstance)
                spinBlackHole2=Tree_Node_Black_Hole_Spin_Standard(hostNode                   )
                spinBlackHole1=Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=iInstance)
             end if
             ! Now calculate the recoil velocity of the binary black hole and check wether it escapes the galaxy. (Note that
             ! we subtract the black hole's own contribution to the potential here.)
             recoilVelocity=Black_Hole_Binary_Recoil_Velocity(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)
             if (recoilVelocity > 0.0d0) then
                if (0.5d0*recoilVelocity**2+Galactic_Structure_Potential(thisNode,0.0d0)-Galactic_Structure_Potential(thisNode&
                     &,0.0d0,componentType=componentTypeBlackHole) > 0.0d0) then
                   blackHoleMassNew=0.0d0
                   blackHoleSpinNew=0.0d0
                end if
             end if
             ! Move the black hole to the host.
             call Galacticus_Output_Tree_Black_Hole_Merger(thisNode,massBlackHole1,massBlackHole2)
             call Tree_Node_Black_Hole_Mass_Set_Standard(hostNode,blackHoleMassNew)
             call Tree_Node_Black_Hole_Spin_Set_Standard(hostNode,blackHoleSpinNew)
          end do
       else
          ! Move black holes from the satellite to the host, giving them the appropriate initial radius.
          if (hostNode%componentExists(componentIndex)) then
             hostIndex=Tree_Node_Black_Hole_Index(hostNode)
             thisIndex=Tree_Node_Black_Hole_Index(thisNode)
             ! Adjust the radii of the black holes in the satellite galaxy.
             forall(iInstance=1:size(thisNode%components(thisIndex)%instance))
                thisNode%components(thisIndex)%instance(iInstance)%properties(radiusIndex,propertyValue)=radiusInitial
                ! Declares them as not having interacted in a triple black hole interaction.
                thisNode%components(thisIndex)%instance(iInstance)%data(tripleInteractionTimeIndex)=0
             end forall
             ! Determine how many non-null black holes the satellite contains.
             nonNullBlackHoleCount=size(thisNode%components(thisIndex)%instance)
             firstNonNullBlackHole=1
             if (thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyValue) <= 0.0d0) then
                nonNullBlackHoleCount=nonNullBlackHoleCount-1
                firstNonNullBlackHole=2
             end if
             if (nonNullBlackHoleCount > 0) then
                ! The host already has some black holes, so append those from the satellite to this list.
                allocate(                                                                  &
                     &   instancesTemporary(                                               &
                     &                       nonNullBlackHoleCount                         &
                     &                      +size(hostNode%components(hostIndex)%instance) &
                     &                     )                                               &
                     &  )
                instancesTemporary(                                               1                     &
                     &             :size(hostNode%components(hostIndex)%instance)                       &
                     &            )=hostNode%components(hostIndex)%instance
                instancesTemporary( size(hostNode%components(hostIndex)%instance)+1                     &
                     &             :size(hostNode%components(hostIndex)%instance)+nonNullBlackHoleCount &
                     &            )=thisNode%components(thisIndex)%instance(firstNonNullBlackHole:firstNonNullBlackHole+nonNullBlackHoleCount-1)
                deallocate(thisNode%components(thisIndex)%instance)
                deallocate(hostNode%components(hostIndex)%instance)
                call Move_Alloc(instancesTemporary,hostNode%components(hostIndex)%instance)
             end if
          else
             ! The host has no black hole of its own. Simply move those from the satellite to the host.
             hostIndex=Tree_Node_Black_Hole_Index(hostNode)
             thisIndex=Tree_Node_Black_Hole_Index(thisNode)
             call Move_Alloc(thisNode%components(thisIndex)%instance,hostNode%components(hostIndex)%instance)
             forall(iInstance=1:size(hostNode%components(hostIndex)%instance))
                hostNode%components(hostIndex)%instance(iInstance)%properties(radiusIndex,propertyValue)=radiusInitial
             end forall
          end if
       end if
       ! Destroy the black hole component in the satellite.
       call thisNode%destroyComponent(componentIndex)
    end if
    return
  end subroutine Black_Hole_Satellite_Merging
  
  subroutine Black_Hole_Standard_Merge_Black_Holes(thisNode)
    !% Merge two black holes.
    use Black_Hole_Binary_Recoil_Velocities
    use Black_Hole_Binary_Mergers
    use Galactic_Structure_Options
    use Galactic_Structure_Potentials
    implicit none
    type(treeNode),   pointer,     intent(inout) :: thisNode
    type(component),  allocatable, dimension(:)  :: componentsTemporary
    integer                                      :: thisIndex,instanceCount
    double precision                             :: blackHoleMassNew,blackHoleSpinNew
    double precision                             :: recoilVelocity,massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2
    
    call Black_Hole_Binary_Merger(Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=mergingInstance), &
         &                        Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=1              ), &
         &                        Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=mergingInstance), &
         &                        Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=1              ), &
         &                        blackHoleMassNew                                                     , &
         &                        blackHoleSpinNew                                                       &
         &                       )
 
    ! Check which black hole is more massive in order to compute an appropriate recoil velocity.
    if (Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=1) >= Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=mergingInstance)) then
       massBlackHole1=Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=1              ) 
       massBlackHole2=Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=mergingInstance)
       spinBlackHole1=Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=1              )
       spinBlackHole2=Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=mergingInstance)
    else
       massBlackHole2=Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=1              ) 
       massBlackHole1=Tree_Node_Black_Hole_Mass_Standard(thisNode,instance=mergingInstance)
       spinBlackHole2=Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=1              )
       spinBlackHole1=Tree_Node_Black_Hole_Spin_Standard(thisNode,instance=mergingInstance)
    end if
    
    ! Calculate the recoil velocity of the binary black hole and check wether it escapes the galaxy
    recoilVelocity=Black_Hole_Binary_Recoil_Velocity(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)
    
    ! Compare the recoil velocity to the potential and determine wether the binary is ejected or stays in the galaxy.           
    if (0.5d0*recoilVelocity**2+Galactic_Structure_Potential(thisNode,0.0d0)-Galactic_Structure_Potential(thisNode&
         &,0.0d0,componentType=componentTypeBlackHole) > 0.0d0) then
       blackHoleMassNew=0.0d0
       blackHoleSpinNew=0.0d0
    end if
    
    ! Set the mass and spin of the central black hole.
    call Galacticus_Output_Tree_Black_Hole_Merger(thisNode,massBlackHole1,massBlackHole2)
    call Tree_Node_Black_Hole_Mass_Set_Standard(thisNode,blackHoleMassNew,instance=1)
    call Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,blackHoleSpinNew,instance=1)
    
    ! Remove the merging black hole from the list.
    thisIndex=Tree_Node_Black_Hole_Index(thisNode)
    instanceCount=size(thisNode%components(thisIndex)%instance)
    call Move_Alloc(thisNode%components(thisIndex)%instance,componentsTemporary)
    allocate(thisNode%components(thisIndex)%instance(instanceCount-1))
    thisNode%components(thisIndex)%instance(1:mergingInstance-1)=componentsTemporary(1:mergingInstance-1)
    if (instanceCount > 2) thisNode%components(thisIndex)%instance(mergingInstance:instanceCount-1)=componentsTemporary(mergingInstance+1:instanceCount)
    deallocate(componentsTemporary)
    return
  end subroutine Black_Hole_Standard_Merge_Black_Holes

  subroutine Black_Hole_Triple_Interaction_Black_Holes(thisNode,interrupt,interruptProcedure)
    !% Handles triple black holes interactions, using conditions similar to those of \cite{volonteri_assembly_2003}.
    use Black_Hole_Binary_Mergers
    use Tree_Nodes
    use Memory_Management
    use Numerical_Constants_Physical
    use Galactic_Structure_Options
    use Galactic_Structure_Potentials
    use Merger_Tree_Active
    use Components
    implicit none
    type(component),  allocatable, dimension(:)  :: componentsTemporary
    logical,                       intent(inout) :: interrupt
    procedure(),      pointer,     intent(inout) :: interruptProcedure
    type(treeNode),   pointer,     intent(inout) :: thisNode
    integer                                      :: thisIndex,ejectedInstance,instanceCount,newBinaryInstance
    double precision                             :: massEjected,massBinary,massRatioIntruder,velocityEjected,kineticEnergyChange,bindingEnergy,newRadius,velocityBinary

    ! We have to distinguish two cases, where a different black hole is ejected, the one with the lowest mass.
    massRatioIntruder=   Tree_Node_Black_Hole_Mass(thisNode,instance=tripleInstance)  &
         &            /( Tree_Node_Black_Hole_Mass(thisNode,instance=1             )  &
         &              +Tree_Node_Black_Hole_Mass(thisNode,instance=binaryInstance)  &
         &             )

    thisIndex=Tree_Node_Black_Hole_Index(thisNode)
    thisNode%components(thisIndex)%instance(tripleInstance)%data(tripleInteractionTimeIndex)=Tree_Node_Time(thisNode)
    if (massRatioIntruder <= 2.0d0) then
       if (Tree_Node_Black_Hole_Mass(thisNode,instance=tripleInstance) <= Tree_Node_Black_Hole_Mass(thisNode,instance=binaryInstance)) then
          newRadius           = Tree_Node_Black_Hole_Radial_Position(thisNode,instance=binaryInstance)/(1.0d0+0.4d0*massRatioIntruder)
          call Tree_Node_Black_Hole_Radial_Position_Set(thisNode,newRadius,instance=binaryInstance)
          bindingEnergy       = gravitationalConstantGalacticus*( Tree_Node_Black_Hole_Mass(thisNode,instance=tripleInstance)                     &
               &                                                 *Tree_Node_Black_Hole_Mass(thisNode,instance=1             )                     &
               &                                                )                                                                                 &
               &               / Tree_Node_Black_Hole_Radial_Position(thisNode,instance=tripleInstance)
          kineticEnergyChange=0.4d0*massRatioIntruder*bindingEnergy
          ejectedInstance    =tripleInstance
          newBinaryInstance  =binaryInstance
       else  
          newRadius           = Tree_Node_Black_Hole_Radial_Position(thisNode,instance=tripleInstance)/(1.0d0+0.4d0*massRatioIntruder )
          call Tree_Node_Black_Hole_Radial_Position_Set(thisNode,newRadius,instance=tripleInstance)
          bindingEnergy       = gravitationalConstantGalacticus*( Tree_Node_Black_Hole_Mass(thisNode,instance=binaryInstance)                     &
               &                                                 *Tree_Node_Black_Hole_Mass(thisNode,instance=1             )                     &
               &                                                )                                                                                 &
               &               / Tree_Node_Black_Hole_Radial_Position(thisNode,instance=binaryInstance)
          kineticEnergyChange=0.4d0*massRatioIntruder*bindingEnergy
          ejectedInstance    =binaryInstance
          newBinaryInstance  =tripleInstance
       end if
    else
       ! This latter case can be referred to as head-on collision. 
       newRadius             =0.53d0*Tree_Node_Black_Hole_Radial_Position(thisNode,instance=tripleInstance)
       call Tree_Node_Black_Hole_Radial_Position_Set(thisNode,newRadius,instance=tripleInstance)
       bindingEnergy         = gravitationalConstantGalacticus*(                                                                                 &
            &                                                    Tree_Node_Black_Hole_Mass(thisNode,instance=binaryInstance)                     &
            &                                                   *Tree_Node_Black_Hole_Mass(thisNode,instance=1             )                     &
            &                                                  )                                                                                 &
            &                 /Tree_Node_Black_Hole_Radial_Position(thisNode,instance=binaryInstance)
       kineticEnergyChange=0.9d0*massRatioIntruder*bindingEnergy
       ejectedInstance    =binaryInstance
       newBinaryInstance  =tripleInstance
    end if          
    ! First we find the lightest black hole and tag it as being ejected.
    massEjected= Tree_Node_Black_Hole_Mass(thisNode,instance=ejectedInstance  )
    massBinary = Tree_Node_Black_Hole_Mass(thisNode,instance=newBinaryInstance) &
         &      +Tree_Node_Black_Hole_Mass(thisNode,instance=1                )
    velocityEjected=dsqrt(kineticEnergyChange/(1.0d0+massEjected/massBinary )/massEjected*2.0d0)
    velocityBinary =dsqrt(kineticEnergyChange/(1.0d0+massBinary /massEjected)/massBinary *2.0d0)
    if (0.5d0*velocityEjected**2+Galactic_Structure_Potential(thisNode,Tree_Node_Black_Hole_Radial_Position(thisNode,instance=ejectedInstance)) > 0.0d0) then
       ! Remove the ejected black hole from the list.
       instanceCount=size(thisNode%components(thisIndex)%instance)
       call Move_Alloc(thisNode%components(thisIndex)%instance,componentsTemporary)
       allocate(thisNode%components(thisIndex)%instance(instanceCount-1))
       thisNode%components(thisIndex)%instance(1:ejectedInstance-1)=componentsTemporary(1:ejectedInstance-1)
       if (instanceCount > 2) thisNode%components(thisIndex)%instance(ejectedInstance:instanceCount-1)=componentsTemporary(ejectedInstance+1:instanceCount)
       deallocate(componentsTemporary)
    end if 
    if (0.5d0*velocityBinary**2+Galactic_Structure_Potential(thisNode,Tree_Node_Black_Hole_Radial_Position(thisNode,instance&
         &=newBinaryInstance))-Galactic_Structure_Potential(thisNode ,Tree_Node_Black_Hole_Radial_Position(thisNode,instance&
         &=newBinaryInstance),componentType=componentTypeBlackHole) > 0.0d0) then
       ! Remove the binary black hole from the list.
       instanceCount=size(thisNode%components(thisIndex)%instance)
       call Move_Alloc(thisNode%components(thisIndex)%instance,componentsTemporary)
       allocate(thisNode%components(thisIndex)%instance(instanceCount-1))
       thisNode%components(thisIndex)%instance(1:newBinaryInstance-1)=componentsTemporary(1:newBinaryInstance-1)
       if (instanceCount > 2) thisNode%components(thisIndex)%instance(newBinaryInstance:instanceCount-1)=componentsTemporary(newBinaryInstance+1:instanceCount)
       deallocate(componentsTemporary)
       ! And set the central black hole as a zero mass component.
       call Tree_Node_Black_Hole_Mass_Set_Standard(thisNode,0.0d0,instance=1)
       call Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,0.0d0,instance=1)
    end if   
    return
  end subroutine Black_Hole_Triple_Interaction_Black_Holes
  
  double precision function Mass_Accretion_Rate(thisNode)
    !% Returns the rate of mass accretion onto the black hole in {\tt thisNode}.
    use Cosmological_Parameters
    use Bondi_Hoyle_Lyttleton_Accretion
    use Galactic_Structure_Densities
    use Galactic_Structure_Options
    use Ideal_Gases_Thermodynamics
    use Black_Hole_Fundamentals
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use Accretion_Disks
    use Hot_Halo_Temperature_Profile
    use Memory_Management
    use Black_Hole_Binary_Separations
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: gasDensityMinimum=1.0d0 ! Lowest gas density to consider when computing accretion rates onto black hole (in units of M_Solar/Mpc^3).
    integer                                  :: iInstance
    double precision                         :: blackHoleMass,gasDensity,relativeVelocity,accretionRadius,jeansLength&
         &,radiativeEfficiency,position(3),hotHaloTemperature,hotModeFraction,accretionRateMaximum
    
    ! Get the active instance.
    iInstance=Tree_Node_Black_Hole_Get_Instance(thisNode)

    ! Ensure arrays are sufficiently large.
    if (allocated(gotAccretionRate)) then
       if (size(gotAccretionRate) < iInstance) then
          call Dealloc_Array(gotAccretionRate                )
          call Dealloc_Array(accretionRate                   )
          call Dealloc_Array(accretionRateHotHalo            )
          call Alloc_Array  (gotAccretionRate    ,[iInstance])
          call Alloc_Array  (accretionRate       ,[iInstance])
          call Alloc_Array  (accretionRateHotHalo,[iInstance])
          gotAccretionRate=.false.
       end if
    else
       call Alloc_Array(gotAccretionRate    ,[iInstance])
       call Alloc_Array(accretionRate       ,[iInstance])
       call Alloc_Array(accretionRateHotHalo,[iInstance])
       gotAccretionRate=.false.
    end if

    if (.not.gotAccretionRate(iInstance)) then
       ! Get black hole mass.
       blackHoleMass=Tree_Node_Black_Hole_Mass_Standard(thisNode)

       ! Check black hole mass is positive.
       if (blackHoleMass > 0.0d0) then

          ! Compute the relative velocity of black hole and gas. We assume that relative motion arises only from the radial
          ! migration of the black hole.
          relativeVelocity=Black_Hole_Binary_Separation_Growth_Rate(thisNode)*Mpc_per_km_per_s_To_Gyr

          ! Contribution from spheroid:

          ! Get the accretion radius. We take this to be the larger of the Bondi-Hoyle radius and the current radius position of
          ! the black hole.
          accretionRadius=max(                                                                                              &
               &               Bondi_Hoyle_Lyttleton_Accretion_Radius(blackHoleMass,bondiHoyleAccretionTemperatureSpheroid) &
               &              ,Tree_Node_Black_Hole_Radial_Position  (thisNode                                            ) &
               &             )

          ! Set the position.
          position=[accretionRadius,0.0d0,0.0d0]
          ! Get density of gas at the galactic center.
          gasDensity=Galactic_Structure_Density(thisNode,position,coordinateSystem=coordinateSystemSpherical,massType&
               &=massTypeGaseous,componentType=componentTypeSpheroid)
          
          ! Check if we have a non-negligible gas density.
          if (gasDensity > gasDensityMinimum) then

             ! Get the Jeans length scale.
             jeansLength=Ideal_Gas_Jeans_Length(bondiHoyleAccretionTemperatureSpheroid,gasDensity)
             ! Limit the smoothing scale to the scale of the spheroid.
             jeansLength=min(jeansLength,Tree_Node_Spheroid_Radius(thisNode))

             ! If the Jeans length exceeds the Bondi-Hoyle-Lyttleton accretion radius, then recompute gas density for a larger
             ! radius, as the gas should be smoothly distributed on scales below the Jeans length.
             if (jeansLength > accretionRadius) then
                ! Set the position.
                position=[jeansLength,0.0d0,0.0d0]
                ! Get density of gas at the galactic center.
                gasDensity=Galactic_Structure_Density(thisNode,position,coordinateSystem=coordinateSystemCylindrical,massType&
                     &=massTypeGaseous,componentType=componentTypeSpheroid)
             end if
             ! Compute the accretion rate.
             accretionRate(iInstance)=bondiHoyleAccretionEnhancementSpheroid*Bondi_Hoyle_Lyttleton_Accretion_Rate(blackHoleMass,gasDensity&
                  &,relativeVelocity,bondiHoyleAccretionTemperatureSpheroid)
             
             ! Get the radiative efficiency of the accretion.
             radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisNode,accretionRate(iInstance))
             
             ! Limit the accretion rate to the Eddington limit.
             if (radiativeEfficiency > 0.0d0) accretionRate(iInstance)=min(accretionRate(iInstance),Black_Hole_Eddington_Accretion_Rate(thisNode)&
                  &/radiativeEfficiency)

          else
             ! Gas density is negative - set zero accretion rate.
             accretionRate(iInstance)=0.0d0
          end if

          ! Contribution from hot halo:

          ! Get halo gas temperature.
          hotHaloTemperature=Hot_Halo_Temperature(thisNode,radius=0.0d0)

          ! Get the accretion radius.
          accretionRadius=Bondi_Hoyle_Lyttleton_Accretion_Radius(blackHoleMass,hotHaloTemperature)
          accretionRadius=min(accretionRadius,Tree_Node_Hot_Halo_Outer_Radius(thisNode))
          
          ! Set the position.
          position=[accretionRadius,0.0d0,0.0d0]

          ! Find the fraction of gas in the halo which is in the hot mode. Set this to unity if hot/cold mode is not to be considered.
          select case (bondiHoyleAccretionHotModeOnly)
          case (.true.)
             hotModeFraction=Hot_Mode_Fraction(thisNode)
          case (.false.)
             hotModeFraction=1.0d0
          end select
             
          ! Get density of gas at the galactic center - scaled by the fraction in the hot accretion mode.
          gasDensity=hotModeFraction*Galactic_Structure_Density(thisNode,position,coordinateSystem=coordinateSystemSpherical,massType&
               &=massTypeGaseous,componentType=componentTypeHotHalo)

          ! Check if we have a non-zero gas density.
          if (gasDensity > gasDensityMinimum) then

             ! Compute the accretion rate.
             accretionRateHotHalo(iInstance)=bondiHoyleAccretionEnhancementHotHalo*Bondi_Hoyle_Lyttleton_Accretion_Rate(blackHoleMass,gasDensity&
                  &,relativeVelocity,hotHaloTemperature,accretionRadius)

             ! Limit the accretion rate to the total mass of the hot halo, divided by the sound crossing time.
             accretionRateMaximum=Tree_Node_Hot_Halo_Mass(thisNode)/(Tree_Node_Hot_Halo_Outer_Radius(thisNode)/(kilo*gigaYear/megaParsec)/Ideal_Gas_Sound_Speed(hotHaloTemperature))
             accretionRateHotHalo(iInstance)=min(accretionRateHotHalo(iInstance),accretionRateMaximum)

             ! Get the radiative efficiency of the accretion.
             radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisNode,accretionRateHotHalo(iInstance))
             
             ! Limit the accretion rate to the Eddington limit.
             if (radiativeEfficiency > 0.0d0) accretionRateHotHalo(iInstance)=min(accretionRateHotHalo(iInstance)&
                  &,Black_Hole_Eddington_Accretion_Rate(thisNode)/radiativeEfficiency)

          else
             ! No gas density, so zero accretion rate.
             accretionRateHotHalo(iInstance)=0.0d0
          end if

       else
          accretionRate       (iInstance)=0.0d0
          accretionRateHotHalo(iInstance)=0.0d0
       end if
       gotAccretionRate(iInstance)=.true.
    end if
    Mass_Accretion_Rate=accretionRate(iInstance)+accretionRateHotHalo(iInstance)
    return
  end function Mass_Accretion_Rate

  integer function Tree_Node_Black_Hole_Index(thisNode)
    !% Ensure the black hole component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Black_Hole_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Black_Hole_Index

  subroutine Black_Hole_Create(thisNode)
    !% Creates a black hole component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message

    ! Display a message.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Creating black hole component for node '
       message=message//thisNode%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if
    ! Create the component.
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    ! Set to the seed mass.
    call Tree_Node_Black_Hole_Mass_Set           (thisNode,blackHoleSeedMass)
    call Tree_Node_Black_Hole_Spin_Set           (thisNode,blackHoleSeedSpin)
    call Tree_Node_Black_Hole_Radial_Position_Set(thisNode,0.0d0            )
    return
  end subroutine Black_Hole_Create

  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Standard_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Standard</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Black_Hole_Standard_Names(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI&
       &,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of black hole properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    use ISO_Varying_String
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    
    if (methodSelected) then
       !@ <outputPropertyGroup>
       !@   <name>blackHole</name>
       !@   <description>Black hole properities</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>blackHoleCount</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Number of super-massive black holes in the galaxy.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>blackHole</group>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='blackHoleCount'
       integerPropertyComments(integerProperty)='Number of super-massive black holes in the galaxy.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>blackHoleMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Mass of the black hole.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>blackHole</group>
       !@ </outputProperty>
       doublePropertyNames    (doubleProperty )='blackHoleMass'
       doublePropertyComments (doubleProperty )='Mass of the black hole.'
       doublePropertyUnitsSI  (doubleProperty )=massSolar
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>blackHoleSpin</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Spin of the black hole.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>blackHole</group>
       !@ </outputProperty>
       doublePropertyNames    (doubleProperty )='blackHoleSpin'
       doublePropertyComments (doubleProperty )='Spin of the black hole.'
       doublePropertyUnitsSI  (doubleProperty )=0.0d0
       if (blackHoleOutputAccretion) then
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>blackHoleAccretionRate</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>Rest-mass accretion rate onto the black hole.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@   <group>blackHole</group>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='blackHoleAccretionRate'
          doublePropertyComments(doubleProperty)='Rest-mass accretion rate onto the black hole.'
          doublePropertyUnitsSI (doubleProperty)=massSolar/gigaYear
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>blackHoleJetPower</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>Power of the black hole-driven jet.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@   <group>blackHole</group>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='blackHoleJetPower'
          doublePropertyComments(doubleProperty)='Power of the black hole-driven jet.'
          doublePropertyUnitsSI (doubleProperty)=massSolar*kilo**2/gigaYear
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>blackHoleRadiativeEfficiency</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>The radiative efficiency of the black hole accretion system.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@   <group>blackHole</group>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='blackHoleRadiativeEfficiency'
          doublePropertyComments(doubleProperty)='The radiative efficiency of the black hole accretion system.'
          doublePropertyUnitsSI (doubleProperty)=0.0d0
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Standard_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Standard_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Standard</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Black_Hole_Standard_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of black hole properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount
    integer,          parameter     :: extraPropertyCount=3

    if (methodSelected) then
       integerPropertyCount=integerPropertyCount+1
       doublePropertyCount =doublePropertyCount +2
       if (blackHoleOutputAccretion) doublePropertyCount=doublePropertyCount+extraPropertyCount
    end if
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Standard_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Standard</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Standard</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Black_Hole_Standard(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store black hole properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Kind_Numbers
    use Accretion_Disks
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)
    integer                                         :: thisIndex,blackHoleCount
    double precision                                :: restMassAccretionRate

    if (methodSelected) then
       ! Ensure the calculation is reset.
       call Tree_Node_Black_Hole_Reset_Standard(thisNode)
       ! Get the rest mass accretion rate.
       restMassAccretionRate=Mass_Accretion_Rate(thisNode)
       
       ! Store the properties.
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Black_Hole_Mass_Standard(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Black_Hole_Spin_Standard(thisNode)
       if (blackHoleOutputAccretion) then
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=restMassAccretionRate
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Accretion_Disk_Jet_Power           (thisNode,restMassAccretionRate)
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Accretion_Disk_Radiative_Efficiency(thisNode,restMassAccretionRate)
       end if
       ! Count number of black holes associated with this galaxy.
       blackHoleCount=0
       if (thisNode%componentExists(componentIndex)) then
          thisIndex=Tree_Node_Black_Hole_Index(thisNode)
          blackHoleCount=size(thisNode%components(thisIndex)%instance)
       end if
       integerProperty=integerProperty+1
       integerBuffer(integerBufferCount,integerProperty)=blackHoleCount
    end if
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Standard

  subroutine Galacticus_Output_Tree_Black_Hole_Merger(thisNode,massBlackHole1,massBlackHole2)
    !% Outputs properties of merging black holes.
    use Tree_Nodes
    use IO_HDF5
    use Galacticus_HDF5
    use Memory_Management
    use Kind_Numbers
    use ISO_Varying_String
    use String_Handling
    use Merger_Tree_Active
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode   
    double precision, intent(in)             :: massBlackHole1,massBlackHole2
    type(hdf5Object)                         :: mergersGroup

    ! Exit if merger data is not to be output.
    if (.not.blackHoleOutputMergers) return

    ! Ignore mergers with zero mass black holes.
    if (massBlackHole2 <= 0.0d0) return

    ! Open the group to which black hole mergers should be written.
    !$omp critical (Galacticus_Output_Tree_Black_Hole_Merger)
    mergersGroup=galacticusOutputFile%openGroup("blackHoleMergers","Black hole mergers data.")
    ! Append to the datasets.
    call mergersGroup%writeDataset([massBlackHole1          ],"massBlackHole1","Mass of the first merging black hole." ,appendTo=.true.)
    call mergersGroup%writeDataset([massBlackHole2          ],"massBlackHole2","Mass of the second merging black hole.",appendTo=.true.)
    call mergersGroup%writeDataset([Tree_Node_Time(thisNode)],"timeOfMerger"  ,"The time of the black hole merger."    ,appendTo=.true.)
    call mergersGroup%writeDataset([activeTreeWeight        ],"volumeWeight"  ,"The weight for the black hole merger." ,appendTo=.true.)
    ! Close the group.
    call mergersGroup%close()  
    !$omp end critical (Galacticus_Output_Tree_Black_Hole_Merger)
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Merger

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Properties</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Galacticus_Output_Tree_Black_Hole_Properties(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Output properties for all black holes in {\tt thisNode}.
    use Tree_Nodes
    use IO_HDF5
    use Galacticus_HDF5
    use Memory_Management
    use Kind_Numbers
    use ISO_Varying_String
    use String_Handling
    use Dark_Matter_Profiles
    use Cosmology_Functions
    use Black_Hole_Binary_Separations
    use Black_Hole_Binary_Separations_Standard
    use Accretion_Disks
    use Cooling_Radii
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),          intent(inout), pointer      :: thisNode
    integer(kind=kind_int8), intent(in)                  :: treeIndex
    integer,                 intent(in)                  :: iOutput
    logical,                 intent(in)                  :: nodePassesFilter
    integer(kind=kind_int8), allocatable,   dimension(:) :: nodeIndex, mergerTreeIndex
    double precision,        allocatable,   dimension(:) :: massAccretionRate, radiativeEfficiency, mass, spin, timescale, radius
    integer                                              :: iBlackHoleNumber,blackHoleCount,thisIndex
    type(hdf5Object)                                     :: blackHolesGroup, outputGroup
    type(varying_string)                                 :: groupName

    ! If black hole output was requested , output their properties.
    if (nodePassesFilter .and. blackHoleOutputData) then

       ! Get a count of the number of black holes present.
       blackHoleCount=1
       if (thisNode%componentExists(componentIndex)) then
          thisIndex=Tree_Node_Black_Hole_Index(thisNode)
          blackHoleCount=size(thisNode%components(thisIndex)%instance)
       end if

       ! Open the output group.
       !$omp critical (HDF5_Access)
       blackHolesGroup=galacticusOutputFile%openGroup("blackHole","Black hole data.")
       groupName="Output"
       groupName=groupName//iOutput
       outputGroup=blackHolesGroup%openGroup(char(groupName),"Properties of black holes for all trees at each output.")  
       !$omp end critical (HDF5_Access)

       ! Allocate array to store profile.
       call Alloc_Array(radius,             [blackHoleCount])
       call Alloc_Array(spin,               [blackHoleCount])
       call Alloc_Array(mass,               [blackHoleCount])
       call Alloc_Array(timescale,          [blackHoleCount])
       call Alloc_Array(massAccretionRate,  [blackHoleCount])
       call Alloc_Array(radiativeEfficiency,[blackHoleCount])
       call Alloc_Array(nodeIndex,          [blackHoleCount])
       call Alloc_Array(mergerTreeIndex,    [blackHoleCount])

       ! Construct arrays of black hole properties.
       if (thisNode%componentExists(componentIndex)) thisIndex=Tree_Node_Black_Hole_Index(thisNode)
       do iBlackHoleNumber=1,blackHoleCount
          if (blackHoleCount > 1) call thisNode%components(thisIndex)%activeInstanceSet(iBlackHoleNumber)
          mass               (iBlackHoleNumber)=Tree_Node_Black_Hole_Mass           (thisNode                              )
          spin               (iBlackHoleNumber)=Tree_Node_Black_Hole_Spin           (thisNode                              )
          radius             (iBlackHoleNumber)=Tree_Node_Black_Hole_Radial_Position(thisNode                              )
          radiativeEfficiency(iBlackHoleNumber)=Accretion_Disk_Radiative_Efficiency (thisNode,Mass_Accretion_Rate(thisNode))
          massAccretionRate  (iBlackHoleNumber)=Mass_Accretion_Rate                 (thisNode                              ) 
          nodeIndex          (iBlackHoleNumber)=thisNode%index()
          mergerTreeIndex    (iBlackHoleNumber)=treeIndex     
          if (iBlackHoleNumber > 1) then
             if (Black_Hole_Binary_Separation_Growth_Rate(thisNode) /= 0.0d0)then
                timescale(iBlackHoleNumber)=-Tree_Node_Black_Hole_Radial_Position    (thisNode) &
                     &                      /Black_Hole_Binary_Separation_Growth_Rate(thisNode) 
             else
                timescale(iBlackHoleNumber)=0.0d0
             end if
          else
             timescale   (iBlackHoleNumber)=0.0d0
          end if
       end do
       if (blackHoleCount > 1) call thisNode%components(thisIndex)%activeInstanceNullify()
       ! Write dataset to the group, first the arrays containing all data.
       !$omp critical (HDF5_Access)
       call outputGroup%writeDataset(mass,               "mass"               ,"The black hole masses.",                appendTo=.true.)
       call outputGroup%writeDataset(spin,               "spin"               ,"The black hole spins.",                 appendTo=.true.)
       call outputGroup%writeDataset(radius,             "radius"             ,"The black hole radial positions.",      appendTo=.true.)
       call outputGroup%writeDataset(timescale,          "timescale"          ,"The black hole timescales for merger.", appendTo=.true.)
       call outputGroup%writeDataset(radiativeEfficiency,"radiativeEfficiency","The black hole radiative efficiencies.",appendTo=.true.)
       call outputGroup%writeDataset(massAccretionRate,  "accretionRate"      ,"The black hole accretion rates.",       appendTo=.true.)
       call outputGroup%writeDataset(nodeIndex,          "nodeIndex"          ,"The black hole host galaxy inices.",    appendTo=.true.)
       call outputGroup%writeDataset(mergerTreeIndex,    "mergerTreeIndex"    ,"The black hole merger tree indices.",   appendTo=.true.)

       ! Deallocatate profile arrays.
       call Dealloc_Array(mass               )
       call Dealloc_Array(spin               )
       call Dealloc_Array(radius             )
       call Dealloc_Array(timescale          )
       call Dealloc_Array(radiativeEfficiency)
       call Dealloc_Array(massAccretionRate  )
       call Dealloc_Array(nodeIndex          )
       call Dealloc_Array(mergerTreeIndex    )
       call outputGroup    %close()
       call blackHolesGroup%close()   
       !$omp end critical (HDF5_Access)
    end if
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Properties

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Black_Hole_Standard_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Black_Hole_Standard_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    integer                                :: iInstance,thisIndex

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'black hole component -> properties:'
          thisIndex=Tree_Node_Black_Hole_Index(thisNode)
          do iInstance=1,size(thisNode%components(thisIndex)%instance)
             write (0,'(2x,a50,1x,e12.6,a2,i4,a2)') 'black hole mass:',Tree_Node_Black_Hole_Mass(thisNode),' [',iInstance,']'
             write (0,'(2x,a50,1x,e12.6,a2,i4,a2)') 'black hole spin:',Tree_Node_Black_Hole_Spin(thisNode),' [',iInstance,']'
          end do
       else
          write (0,'(1x,a)'           ) 'black hole component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Black_Hole_Standard_Dump

  !# <decodePropertyIdentifiersTask>
  !#  <unitName>Black_Hole_Standard_Property_Identifiers_Decode</unitName>
  !# </decodePropertyIdentifiersTask>
  subroutine Black_Hole_Standard_Property_Identifiers_Decode(propertyComponent,propertyObject,propertyIndex,matchedProperty,propertyName)
    !% Decodes property identifiers to property names for the standard black hole module.
    use ISO_Varying_String
    implicit none
    integer,              intent(in)    :: propertyComponent,propertyObject,propertyIndex
    logical,              intent(inout) :: matchedProperty
    type(varying_string), intent(inout) :: propertyName

    if (methodSelected.and..not.matchedProperty) then
       if (propertyComponent == componentIndex) then
          matchedProperty=.true.
          propertyName="blackHole:"
          select case (propertyObject)
          case (objectTypeProperty)
             select case (propertyIndex)
             case (massIndex  )
                propertyName=propertyName//":mass"
             case (spinIndex  )
                propertyName=propertyName//":spin"
             case (radiusIndex)
                propertyName=propertyName//":radius"
             end select
          end select
       end if
    end if

    return
  end subroutine Black_Hole_Standard_Property_Identifiers_Decode

  double precision function Hot_Mode_Fraction(thisNode)
    !% A simple interpolating function which is used as a measure of the fraction of a halo which is in the hot accretion mode.
    use Cooling_Radii
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: coolingRadiusFractionalTransitionMinimum=0.9d0
    double precision, parameter              :: coolingRadiusFractionalTransitionMaximum=1.0d0
    double precision                         :: x,coolingRadiusFractional

    coolingRadiusFractional=Cooling_Radius(thisNode)/Dark_Matter_Halo_Virial_Radius(thisNode)
    if      (coolingRadiusFractional < coolingRadiusFractionalTransitionMinimum) then
       Hot_Mode_Fraction=1.0d0
    else if (coolingRadiusFractional > coolingRadiusFractionalTransitionMaximum) then
       Hot_Mode_Fraction=0.0d0
    else
       x=      (coolingRadiusFractional                 -coolingRadiusFractionalTransitionMinimum) &
            & /(coolingRadiusFractionalTransitionMaximum-coolingRadiusFractionalTransitionMinimum)
       Hot_Mode_Fraction=x**2*(2.0d0*x-3.0d0)+1.0d0
    end if
    return
  end function Hot_Mode_Fraction
  
end module Tree_Node_Methods_Black_Hole
