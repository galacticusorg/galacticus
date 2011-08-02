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


!% Contains a module of Sérsic-profile spheroid tree node methods.

module Tree_Node_Methods_Sersic_Spheroid
  !% Implement Sérsic spheroid tree node methods.
  use Tree_Nodes
  use Histories
  use Components
  use Stellar_Population_Properties
  use FGSL
  implicit none
  private
  public :: Tree_Node_Methods_Sersic_Spheroid_Initialize, Sersic_Spheroid_Satellite_Merging,&
       & Galacticus_Output_Tree_Spheroid_Sersic, Galacticus_Output_Tree_Spheroid_Sersic_Property_Count,&
       & Galacticus_Output_Tree_Spheroid_Sersic_Names, Sersic_Spheroid_Radius_Solver, Sersic_Spheroid_Enclosed_Mass,&
       & Sersic_Spheroid_Density, Sersic_Spheroid_Rotation_Curve, Tree_Node_Spheroid_Post_Evolve_Sersic,&
       & Tree_Node_Methods_Sersic_Spheroid_Dump, Sersic_Spheroid_Radius_Solver_Plausibility, Sersic_Spheroid_Scale_Set&
       &,Sersic_Spheroid_Star_Formation_History_Output,Tree_Node_Methods_Sersic_Spheroid_Thread_Initialize,&
       & Sersic_Spheroid_Property_Identifiers_Decode, Sersic_Profile_Tabulate_State_Store, Sersic_Profile_Tabulate_State_Retrieve
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Internal count of abundances and work arrays.
  integer                                     :: abundancesCount
  double precision, allocatable, dimension(:) :: abundancesOutflowRate,abundancesValue,abundancesDisk,abundancesSpheroid
  !$omp threadprivate(abundancesOutflowRate,abundancesValue,abundancesDisk,abundancesSpheroid)

  ! Internal count of luminosities and work arrays.
  integer                                     :: luminositiesCount
  double precision, allocatable, dimension(:) :: stellarLuminositiesRates,luminositiesDisk,luminositiesSpheroid
  !$omp threadprivate(stellarLuminositiesRates,luminositiesDisk,luminositiesSpheroid)

  ! Property indices.
  integer, parameter :: propertyCountBase=3, dataCountBase=2, historyCount=2
  integer            :: propertyCount      , dataCount
  integer, parameter :: angularMomentumIndex=1
  integer, parameter :: gasMassIndex        =2
  integer, parameter :: stellarMassIndex    =3
  integer            :: gasAbundancesIndex      ,gasAbundancesIndexEnd
  integer            :: stellarAbundancesIndex  ,stellarAbundancesIndexEnd
  integer            :: stellarLuminositiesIndex,stellarLuminositiesIndexEnd
  integer, parameter :: radiusIndex              =1
  integer, parameter :: velocityIndex            =2
  integer, parameter :: stellarHistoryIndex      =1
  integer, parameter :: starFormationHistoryIndex=2

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
  !# <treeNodeMethodsPointer type="history">
  !#  <methodName>Tree_Node_Spheroid_Star_Formation_History</methodName>
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
  integer          :: spheroidSersicIndex
  double precision :: spheroidSersicCoefficient
  double precision :: spheroidEnergeticOutflowMassRate,spheroidOutflowTimescaleMinimum
  double precision :: spheroidMassToleranceAbsolute

  ! Storage for the star formation history time range, used whene extending this range.
  double precision, allocatable, dimension(:) :: starFormationHistoryTemplate
  !$omp threadprivate(starFormationHistoryTemplate)

  ! Options controlling output.
  logical          :: spheroidOutputStarFormationRate

  ! Global variables used in integrations.
  double precision :: radiusStart

  ! Arrays to store the tabulated Sérsic profile.
  logical                                     :: sersicTableInitialized=.false.
  integer                                     :: sersicTableCount
  double precision                            :: sersicTableRadiusMinimum   =1.0d-3, &
       &                                         sersicTableRadiusMaximum   =1.0d+3
  double precision                            :: sersicTable3dHalfMassRadius=1.0d+0
  double precision                            :: sersicTableAngularMomentumHalfMassRadiusToMeanRatio
  integer,          parameter                 :: sersicTablePointsPerDecade=100
  double precision, allocatable, dimension(:) :: sersicTableRadius,sersicTableDensity,sersicTableEnclosedMass
  type(fgsl_interp)                           :: sersicTableInterpolationObject
  type(fgsl_interp_accel)                     :: sersicTableInterpolationAccelerator
  logical                                     :: sersicTableInterpolationReset=.true.

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Sersic_Spheroid_Initialize</unitName>
  !#  <optionName>treeNodeMethodSpheroid</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Sersic_Spheroid_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node Sérsic spheroid methods module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    use Galacticus_Error
    use Stellar_Population_Properties_Luminosities
    use Abundances_Structure
    use Memory_Management
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'Sersic') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Sérsic spheroid method selected [component index '
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
       Tree_Node_Spheroid_Radius                                  => Sersic_Spheroid_Radius
       Tree_Node_Spheroid_Radius_Set                              => null()
       Tree_Node_Spheroid_Radius_Rate_Adjust                      => null()
       Tree_Node_Spheroid_Radius_Rate_Compute                     => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Half_Mass_Radius                        => Sersic_Spheroid_Half_Mass_Radius
       Tree_Node_Spheroid_Half_Mass_Radius_Set                    => null()
       Tree_Node_Spheroid_Half_Mass_Radius_Rate_Adjust            => null()
       Tree_Node_Spheroid_Half_Mass_Radius_Rate_Compute           => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Velocity                                => Sersic_Spheroid_Velocity
       Tree_Node_Spheroid_Velocity_Set                            => null()
       Tree_Node_Spheroid_Velocity_Rate_Adjust                    => null()
       Tree_Node_Spheroid_Velocity_Rate_Compute                   => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Angular_Momentum                        => Tree_Node_Spheroid_Angular_Momentum_Sersic
       Tree_Node_Spheroid_Angular_Momentum_Set                    => Tree_Node_Spheroid_Angular_Momentum_Set_Sersic
       Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust            => Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Angular_Momentum_Rate_Compute           => Tree_Node_Spheroid_Angular_Momentum_Rate_Compute_Sersic

       Tree_Node_Spheroid_Gas_Mass                                => Tree_Node_Spheroid_Gas_Mass_Sersic
       Tree_Node_Spheroid_Gas_Mass_Set                            => Tree_Node_Spheroid_Gas_Mass_Set_Sersic
       Tree_Node_Spheroid_Gas_Mass_Rate_Adjust                    => Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Gas_Mass_Rate_Compute                   => Tree_Node_Spheroid_Gas_Mass_Rate_Compute_Sersic

       Tree_Node_Spheroid_Stellar_Mass                            => Tree_Node_Spheroid_Stellar_Mass_Sersic
       Tree_Node_Spheroid_Stellar_Mass_Set                        => Tree_Node_Spheroid_Stellar_Mass_Set_Sersic
       Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust                => Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Stellar_Mass_Rate_Compute               => Tree_Node_Spheroid_Stellar_Mass_Rate_Compute_Sersic

       Tree_Node_Spheroid_Stellar_Abundances                      => Tree_Node_Spheroid_Stellar_Abundances_Sersic
       Tree_Node_Spheroid_Stellar_Abundances_Set                  => Tree_Node_Spheroid_Stellar_Abundances_Set_Sersic
       Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust          => Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Stellar_Abundances_Rate_Compute         => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Gas_Abundances                          => Tree_Node_Spheroid_Gas_Abundances_Sersic
       Tree_Node_Spheroid_Gas_Abundances_Set                      => Tree_Node_Spheroid_Gas_Abundances_Set_Sersic
       Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust              => Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Gas_Abundances_Rate_Compute             => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Stellar_Luminosities                    => Tree_Node_Spheroid_Stellar_Luminosities_Sersic
       Tree_Node_Spheroid_Stellar_Luminosities_Set                => Tree_Node_Spheroid_Stellar_Luminosities_Set_Sersic
       Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust        => Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Stellar_Luminosities_Rate_Compute       => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_SFR                                     => Sersic_Spheroid_SFR
       Tree_Node_Spheroid_SFR_Set                                 => null()
       Tree_Node_Spheroid_SFR_Rate_Adjust                         => null()
       Tree_Node_Spheroid_SFR_Rate_Compute                        => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Stellar_Properties_History              => Tree_Node_Spheroid_Stellar_Properties_History_Sersic
       Tree_Node_Spheroid_Stellar_Properties_History_Set          => Tree_Node_Spheroid_Stellar_Properties_History_Set_Sersic
       Tree_Node_Spheroid_Stellar_Properties_History_Rate_Adjust  => Tree_Node_Spheroid_Stellar_Prprts_History_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Stellar_Properties_History_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Star_Formation_History                  => Tree_Node_Spheroid_Star_Formation_History_Sersic
       Tree_Node_Spheroid_Star_Formation_History_Set              => Tree_Node_Spheroid_Star_Formation_History_Set_Sersic
       Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust      => Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Star_Formation_History_Rate_Compute     => Tree_Node_Rate_Rate_Compute_Dummy

       ! Associate pipes with procedures.
       Tree_Node_Spheroid_Gas_Sink                                => Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Sersic
       Tree_Node_Spheroid_Gas_Energy_Input                        => Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Sersic

       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>spheroidSersicIndex</name>
       !@   <defaultValue>4</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The index in the Sérsic profile describing the spheroid component.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidSersicIndex',spheroidSersicIndex,defaultValue=4)
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
       !@ <inputParameter>
       !@   <name>spheroidOutflowTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the spheroid dynamical time) on which outflows may deplete gas in the spheroid.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidOutflowTimescaleMinimum',spheroidOutflowTimescaleMinimum,defaultValue=1.0d-3)
       !@ <inputParameter>
       !@   <name>spheroidOutputStarFormationRate</name>
       !@   <defaultValue>false.</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Determines whether or not the star formation rate in the spheroid of each galaxy will be output.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidOutputStarFormationRate',spheroidOutputStarFormationRate,defaultValue=.false.)

    end if
    return
  end subroutine Tree_Node_Methods_Sersic_Spheroid_Initialize
  
  !# <treeNodeCreateThreadInitialize>
  !#  <unitName>Tree_Node_Methods_Sersic_Spheroid_Thread_Initialize</unitName>
  !# </treeNodeCreateThreadInitialize>
  subroutine Tree_Node_Methods_Sersic_Spheroid_Thread_Initialize
    !% Initializes each thread for the tree node Sersic spheroid methods module.
    use Memory_Management
    implicit none

    ! Check if this implementation is selected.
    if (methodSelected.and..not.allocated(abundancesOutflowRate)) then

       ! Allocate work arrays for abundances.
       call Alloc_Array(abundancesOutflowRate,[abundancesCount])
       call Alloc_Array(abundancesValue      ,[abundancesCount])
       call Alloc_Array(abundancesDisk       ,[abundancesCount])
       call Alloc_Array(abundancesSpheroid   ,[abundancesCount])

       ! Allocate work arrays for luminosities.
       call Alloc_Array(stellarLuminositiesRates,[luminositiesCount])
       call Alloc_Array(luminositiesDisk        ,[luminositiesCount])
       call Alloc_Array(luminositiesSpheroid    ,[luminositiesCount])

    end if
    return
  end subroutine Tree_Node_Methods_Sersic_Spheroid_Thread_Initialize
  
  !# <postEvolveTask>
  !# <unitName>Tree_Node_Spheroid_Post_Evolve_Sersic</unitName>
  !# </postEvolveTask>
  subroutine Tree_Node_Spheroid_Post_Evolve_Sersic(thisNode)
    !% Trim histories attached to the spheroid.
    use Galacticus_Display
    use String_Handling
    use ISO_Varying_String
    use Histories
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, save                   :: fractionalErrorMaximum=0.0d0
    integer                                  :: thisIndex
    double precision                         :: specificAngularMomentum,fractionalError,spheroidMass
    character(len=20)                        :: valueString
    type(varying_string)                     :: message

    if (methodSelected .and. thisNode%componentExists(componentIndex)) then
       ! Get the index of the spheroid component.
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)

       ! Trim the stellar populations properties future history.
       call thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)%trim(Tree_Node_Time(thisNode))

       ! Trap negative gas masses.
       if (Tree_Node_Spheroid_Gas_Mass(thisNode) < 0.0d0) then
          
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=dabs(thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue))&
               &/(thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)&
               &+dabs(thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)))
          !$omp critical (Sersic_Spheroid_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.          
             message='Warning: spheroid has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index            = '//thisNode%index() //char(10)
             write (valueString,'(e12.6)') thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex    ,propertyValue)
             message=message//'  Spheroid gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
             message=message//'  Spheroid stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure         = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Gas mass will be reset to zero (in future cases also).'//char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerically inaccuracy in the ODE solutions.'//char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call Galacticus_Display_Message(message,verbosityWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Sersic_Spheroid_Post_Evolve_Check)
          
          ! Get the specific angular momentum of the spheroid material
          spheroidMass= thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex    ,propertyValue) &
               &       +thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
          if (spheroidMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex                                    ,propertyValue)=0.0d0
             thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex  :stellarAbundancesIndexEnd  ,propertyValue)=0.0d0  
             thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyValue)=0.0d0
          else
             specificAngularMomentum=thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)&
                  &/spheroidMass
          end if
          
          ! Reset the gas, abundances and angular momentum of the spheroid.
          thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex                            ,propertyValue)=0.0d0
          thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)=0.0d0
          thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex                    ,propertyValue)= &
               & specificAngularMomentum*thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
       end if
       
    end if
    return
  end subroutine Tree_Node_Spheroid_Post_Evolve_Sersic

  double precision function Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode,instance)
    !% Return the node Sérsic spheroid gas mass.
    implicit none
    type(treeNode), pointer,  intent(inout) :: thisNode
    integer,        optional, intent(in)    :: instance
    integer                                 :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Gas_Mass_Sersic=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)
    else
       Tree_Node_Spheroid_Gas_Mass_Sersic=0.0d0
    end if
    return
  end function Tree_Node_Spheroid_Gas_Mass_Sersic

  subroutine Tree_Node_Spheroid_Gas_Mass_Set_Sersic(thisNode,mass,instance)
    !% Set the node Sérsic spheroid gas mass.
    implicit none
    type(treeNode),   pointer,  intent(inout) :: thisNode
    double precision,           intent(in)    :: mass
    integer,          optional, intent(in)    :: instance
    integer                                   :: thisIndex

    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Spheroid_Gas_Mass_Set_Sersic

  subroutine Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Account for a sink of gaseous material in the Sérsic spheroid.
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer,  intent(inout) :: thisNode
    logical,                    intent(inout) :: interrupt
    procedure(),      pointer,  intent(inout) :: interruptProcedure
    double precision,           intent(in)    :: rateAdjustment
    integer,          optional, intent(in)    :: instance
    integer                                   :: thisIndex
    double precision                          :: gasMass,stellarMass
    
    ! If no Sérsic spheroid component currently exists and we have some sink from then there is a problem.
    ! spheroid. If sink has zero rate, just return instead.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) call Galacticus_Error_Report('Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Sersic','mass sink from non-existant Sérsic spheroid')
       return
    end if
    ! Trap cases where an attempt is made to add gas via this sink function.
    if (rateAdjustment > 0.0d0) call Galacticus_Error_Report('Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Sersic','attempt to add mass via sink in Sérsic spheroid')
    ! Return if no adjustment is being made.
    if (rateAdjustment == 0.0d0) return
    ! Get the index of the component.
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)    
    ! Get the gas mass present.
    gasMass=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)
    ! Get the stellar mass present.
    stellarMass=thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0 .and. gasMass+stellarMass > 0.0d0) then
       thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)&
            &=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)+rateAdjustment
       thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyDerivative)&
            &=thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyDerivative) &
            & +rateAdjustment/(gasMass+stellarMass) &
            & *thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)
       thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
            &=thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative) & 
            & +(rateAdjustment/gasMass) &
            & *thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)
    end if
    return
  end subroutine Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Sersic

  subroutine Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Handles input of energy into the spheroid gas from other components (e.g. black holes). The energy input rate should be in
    !% units of $M_\odot$ km$^2$ s$^{-2}$ Gyr$^{-1}$.
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer,  intent(inout)    :: thisNode
    logical,                    intent(inout)    :: interrupt
    procedure(),      pointer,  intent(inout)    :: interruptProcedure
    double precision,           intent(in)       :: rateAdjustment
    integer,          optional, intent(in)       :: instance
    double precision, dimension(abundancesCount) :: abundancesOutflowRate
    integer                                      :: thisIndex
    double precision                             :: gasMass,stellarMass,massOutflowRate,angularMomentumOutflowRate&
         &,spheroidVelocity

    ! If no Sérsic spheroid component currently exists then energy input to it has no effect.
    if (.not.thisNode%componentExists(componentIndex)) return

    ! Trap cases where an attempt is made to remove energy via this input function.
    if (rateAdjustment < 0.0d0) call Galacticus_Error_Report('Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Sersic' &
         & ,'attempt to remove energy via input pipe in Sérsic spheroid')
    
    ! Return if no adjustment is being made.
    if (rateAdjustment == 0.0d0) return
    ! Get the index of the component.
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)    
    ! Get the gas mass present.
    gasMass    =thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex    ,propertyValue)
    ! Get the stellar mass present.
    stellarMass=thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0 .and. gasMass+stellarMass > 0.0d0) then
       ! Compute outflow rates of quantities and adjust rates in the spheroid appropriately.
       spheroidVelocity=Sersic_Spheroid_Velocity(thisNode)
       if (spheroidVelocity > 0.0d0) then
          massOutflowRate=spheroidEnergeticOutflowMassRate*rateAdjustment/spheroidVelocity**2
          thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative) &
               &=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)-massOutflowRate
          angularMomentumOutflowRate=massOutflowRate/(gasMass+stellarMass) &
               & *thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)
          thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyDerivative)&
               &=thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyDerivative) &
               & -angularMomentumOutflowRate
          abundancesOutflowRate=(massOutflowRate/gasMass) &
               & *thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)
          thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
               &=thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative) & 
               & -abundancesOutflowRate
          ! Add outflowing rates to the hot halo component.
          call Tree_Node_Hot_Halo_Outflow_Mass_To            (thisNode,interrupt,interruptProcedure,massOutflowRate           )
          call Tree_Node_Hot_Halo_Outflow_Angular_Momentum_To(thisNode,interrupt,interruptProcedure,angularMomentumOutflowRate)
          call Tree_Node_Hot_Halo_Outflow_Abundances_To      (thisNode,interrupt,interruptProcedure,abundancesOutflowRate     )
       end if
    end if
    return
  end subroutine Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Sersic

  subroutine Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node Sérsic spheroid gas mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer,  intent(inout) :: thisNode
    logical,                    intent(inout) :: interrupt
    procedure(),      pointer,  intent(inout) :: interruptProcedure
    double precision,           intent(in)    :: rateAdjustment
    integer,          optional, intent(in)    :: instance
    integer                                   :: thisIndex
    
    ! If no Sérsic spheroid component currently exists and we have some cooling into it then interrupt and create an Sérsic
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Sersic_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Sersic

  subroutine Tree_Node_Spheroid_Gas_Mass_Rate_Compute_Sersic(thisNode,interrupt,interruptProcedure)
    !% Compute the Sérsic spheroid node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure

    ! Do nothing here - work will be done elsewhere.
    return
  end subroutine Tree_Node_Spheroid_Gas_Mass_Rate_Compute_Sersic

  double precision function Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode,instance)
    !% Return the node Sérsic spheroid stellar mass.
    implicit none
    type(treeNode), pointer,  intent(inout) :: thisNode
    integer,        optional, intent(in)    :: instance
    integer                                 :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Stellar_Mass_Sersic=thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
    else
       Tree_Node_Spheroid_Stellar_Mass_Sersic=0.0d0
    end if
    return
  end function Tree_Node_Spheroid_Stellar_Mass_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Mass_Set_Sersic(thisNode,mass,instance)
    !% Set the node Sérsic spheroid stellar mass.
    implicit none
    type(treeNode),   pointer,  intent(inout) :: thisNode
    double precision,           intent(in)    :: mass
    integer,          optional, intent(in)    :: instance
    integer                                   :: thisIndex

    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Spheroid_Stellar_Mass_Set_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node Sérsic spheroid stellar mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer,  intent(inout) :: thisNode
    logical,                    intent(inout) :: interrupt
    procedure(),      pointer,  intent(inout) :: interruptProcedure
    double precision,           intent(in)    :: rateAdjustment
    integer,          optional, intent(in)    :: instance
    integer                                   :: thisIndex
    
    ! If no Sérsic spheroid component currently exists and we have some mass flow into it then interrupt and create a Sérsic
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Sersic_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Mass_Rate_Compute_Sersic(thisNode,interrupt,interruptProcedure)
    !% Compute the Sérsic spheroid node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    use Star_Formation_Feedback_Spheroids
    use Star_Formation_Feedback_Expulsion_Spheroids
    use Abundances_Structure
    use Galactic_Structure_Options
    use Galacticus_Output_Star_Formation_Histories
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),            pointer, intent(inout) :: thisNode
    logical,                            intent(inout) :: interrupt
    procedure(),               pointer, intent(inout) :: interruptProcedure
    integer                                           :: thisIndex
    double precision                                  :: starFormationRate,stellarMassRate ,fuelMassRate,fuelMass &
         &,massOutflowRate,spheroidMass,angularMomentumOutflowRate,energyInputRate,gasMass,spheroidDynamicalTime&
         &,massOutflowRateToHotHalo,massOutflowRateFromHalo,outflowToHotHaloFraction
    type(abundancesStructure), save                   :: fuelAbundances,stellarAbundancesRates,fuelAbundancesRates
    !$omp threadprivate(fuelAbundances,stellarAbundancesRates,fuelAbundancesRates)

    ! Compute star formation rate if this node has a spheroid.
    if (thisNode%componentExists(componentIndex)) then

       ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
       if     (    Tree_Node_Spheroid_Angular_Momentum(thisNode) < 0.0d0 &
            & .or. Tree_Node_Spheroid_Radius          (thisNode) < 0.0d0 &
            & .or. Tree_Node_Spheroid_Gas_Mass        (thisNode) < 0.0d0 &
            & ) return

       ! Find the star formation timescale.
       starFormationRate=Sersic_Spheroid_SFR(thisNode)

       ! Get the available fuel mass.
       fuelMass=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)
       
       ! Find the metallicity of the fuel supply.
       call Tree_Node_Spheroid_Gas_Abundances_Sersic(thisNode,abundancesValue)
       call fuelAbundances%pack(abundancesValue)
       call fuelAbundances%massToMassFraction(fuelMass)

       ! Get the component index.
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       
       ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
       call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,componentTypeSpheroid,thisNode&
            &,thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex),stellarMassRate,stellarAbundancesRates &
            &,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
       
       ! Adjust rates.
       call Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Sersic        (thisNode,interrupt,interruptProcedure,stellarMassRate         )
       call Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Sersic            (thisNode,interrupt,interruptProcedure,fuelMassRate            )
       call stellarAbundancesRates%unpack(abundancesValue)
       call Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Sersic  (thisNode,interrupt,interruptProcedure,abundancesValue         )
       call fuelAbundancesRates%unpack(abundancesValue)
       call Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Sersic      (thisNode,interrupt,interruptProcedure,abundancesValue         )
       call Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,stellarLuminositiesRates)

       ! Record the star formation history.
       call Star_Formation_History_Record(thisNode,thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)&
            &,fuelAbundances,starFormationRate)
       
       ! Find rate of outflow of material from the spheroid and pipe it to the outflowed reservoir.
       massOutflowRateToHotHalo=Star_Formation_Feedback_Spheroid_Outflow_Rate          (thisNode,starFormationRate,energyInputRate)
       massOutflowRateFromHalo =Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
       massOutflowRate         =massOutflowRateToHotHalo+massOutflowRateFromHalo
       if (massOutflowRate > 0.0d0) then
          ! Find the fraction of material which outflows to the hot halo.
          outflowToHotHaloFraction=massOutflowRateToHotHalo/massOutflowRate
          
          ! Get the masses of the spheroid.
          gasMass=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)
          spheroidMass=gasMass+Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
          
          ! Limit the outflow rate timescale to a multiple of the dynamical time.
          spheroidDynamicalTime=Mpc_per_km_per_s_To_Gyr*Tree_Node_Spheroid_Radius(thisNode)/Tree_Node_Spheroid_Velocity(thisNode)
          ! Limit the mass outflow rate.
          massOutflowRate=min(massOutflowRate,gasMass/spheroidOutflowTimescaleMinimum/spheroidDynamicalTime)
          call Tree_Node_Hot_Halo_Outflow_Mass_To                       (thisNode,interrupt,interruptProcedure,&
               & massOutflowRate*outflowToHotHaloFraction)
          call Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Sersic        (thisNode,interrupt,interruptProcedure,&
               &-massOutflowRate                         )
          
          ! Compute the angular momentum outflow rate.
          if (spheroidMass > 0.0d0) then
             angularMomentumOutflowRate=(massOutflowRate/spheroidMass)*Tree_Node_Spheroid_Angular_Momentum_Sersic(thisNode)
             call Tree_Node_Hot_Halo_Outflow_Angular_Momentum_To           (thisNode,interrupt,interruptProcedure,&
                  & angularMomentumOutflowRate*outflowToHotHaloFraction)
             call Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,&
                  &-angularMomentumOutflowRate                         )
          end if
          
          ! Compute the abundances outflow rate.
          call Tree_Node_Spheroid_Gas_Abundances_Sersic(thisNode,abundancesValue)
          call Abundances_Mass_To_Mass_Fraction(abundancesValue,gasMass)
          abundancesOutflowRate=massOutflowRate*abundancesValue
          call Tree_Node_Hot_Halo_Outflow_Abundances_To                 (thisNode,interrupt,interruptProcedure,&
               & abundancesOutflowRate*outflowToHotHaloFraction)
          call Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Sersic  (thisNode,interrupt,interruptProcedure,&
               &-abundancesOutflowRate                         )
       end if
    end if
    return
  end subroutine Tree_Node_Spheroid_Stellar_Mass_Rate_Compute_Sersic

  subroutine Tree_Node_Spheroid_Gas_Abundances_Sersic(thisNode,abundanceMasses)
    !% Return the node Sérsic spheroid gas abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: abundanceMasses(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       abundanceMasses(:)=thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd&
            &,propertyValue)
    else
       abundanceMasses(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Spheroid_Gas_Abundances_Sersic

  subroutine Tree_Node_Spheroid_Gas_Abundances_Set_Sersic(thisNode,abundanceMasses)
    !% Set the node Sérsic spheroid gas abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: abundanceMasses(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)=abundanceMasses(:)
    return
  end subroutine Tree_Node_Spheroid_Gas_Abundances_Set_Sersic

  subroutine Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node Sérsic spheroid gas abundance masses rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no Sérsic spheroid component currently exists and we have some cooling into it then interrupt and create an Sérsic
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Sersic_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Abundances_Sersic(thisNode,abundanceMasses)
    !% Return the node Sérsic spheroid stellar abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: abundanceMasses(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       abundanceMasses(:)=thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd&
            &,propertyValue)
    else
       abundanceMasses(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Spheroid_Stellar_Abundances_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Abundances_Set_Sersic(thisNode,abundanceMasses)
    !% Set the node Sérsic spheroid stellar abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: abundanceMasses(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyValue)=abundanceMasses(:)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Abundances_Set_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node Sérsic spheroid stellar abundance masses rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no Sérsic spheroid component currently exists and we have some cooling into it then interrupt and create an Sérsic
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Sersic_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Luminosities_Sersic(thisNode,luminosities)
    !% Return the node Sérsic spheroid stellar luminosities.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: luminosities(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       luminosities(:)=thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd&
            &,propertyValue)
    else
       luminosities(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Spheroid_Stellar_Luminosities_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Luminosities_Set_Sersic(thisNode,luminosities)
    !% Set the node Sérsic spheroid stellar luminosities.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: luminosities(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyValue)=luminosities(:)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Luminosities_Set_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node Sérsic spheroid stellar luminosity rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no Sérsic spheroid component currently exists and we have some cooling into it then interrupt and create an Sérsic
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Sersic_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Sersic

  double precision function Tree_Node_Spheroid_Angular_Momentum_Sersic(thisNode,instance)
    !% Return the node Sérsic spheroid gas mass.
    implicit none
    type(treeNode), pointer,  intent(inout) :: thisNode
    integer,        optional, intent(in)    :: instance
    integer                                 :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Angular_Momentum_Sersic=thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)
    else
       Tree_Node_Spheroid_Angular_Momentum_Sersic=0.0d0
    end if
    return
  end function Tree_Node_Spheroid_Angular_Momentum_Sersic

  subroutine Tree_Node_Spheroid_Angular_Momentum_Set_Sersic(thisNode,angularMomentum,instance)
    !% Set the node Sérsic spheroid gas mass.
    implicit none
    type(treeNode),   pointer,  intent(inout) :: thisNode
    double precision,           intent(in)    :: angularMomentum
    integer,          optional, intent(in)    :: instance
    integer                                   :: thisIndex

    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)=angularMomentum
    return
  end subroutine Tree_Node_Spheroid_Angular_Momentum_Set_Sersic

  subroutine Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node Sérsic spheroid gas mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer,  intent(inout) :: thisNode
    logical,                    intent(inout) :: interrupt
    procedure(),      pointer,  intent(inout) :: interruptProcedure
    double precision,           intent(in)    :: rateAdjustment
    integer,          optional, intent(in)    :: instance
    integer                                   :: thisIndex
    
    ! If no Sérsic spheroid component currently exists and we have some mass flow into it then interrupt and create a Sérsic
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Sersic_Spheroid_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyDerivative) &
         &=thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Sersic

  subroutine Tree_Node_Spheroid_Angular_Momentum_Rate_Compute_Sersic(thisNode,interrupt,interruptProcedure)
    !% Compute the Sérsic spheroid node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure

    ! No need to do anything here - the work will be done elsewhere.
    return
  end subroutine Tree_Node_Spheroid_Angular_Momentum_Rate_Compute_Sersic

  function Tree_Node_Spheroid_Stellar_Properties_History_Sersic(thisNode)
    !% Return spheroid stellar properties history.
    use Histories
    implicit none
    type(history)                          :: Tree_Node_Spheroid_Stellar_Properties_History_Sersic
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Stellar_Properties_History_Sersic=thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)
    else
       Tree_Node_Spheroid_Stellar_Properties_History_Sersic=nullHistory
    end if
    return
  end function Tree_Node_Spheroid_Stellar_Properties_History_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Properties_History_Set_Sersic(thisNode,thisHistory)
    !% Set the node spheroid stellar properties history.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(history),           intent(in)    :: thisHistory
    integer                                :: thisIndex

    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    call thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)%destroy()
    thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)=thisHistory
    return
  end subroutine Tree_Node_Spheroid_Stellar_Properties_History_Set_Sersic

  function Tree_Node_Spheroid_Star_Formation_History_Sersic(thisNode)
    !% Return the spheroid star formation history.
    use Histories
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(history)                          :: Tree_Node_Spheroid_Star_Formation_History_Sersic
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Star_Formation_History_Sersic=thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)
    else
       Tree_Node_Spheroid_Star_Formation_History_Sersic=nullHistory
    end if
    return
  end function Tree_Node_Spheroid_Star_Formation_History_Sersic

  subroutine Tree_Node_Spheroid_Star_Formation_History_Set_Sersic(thisNode,thisHistory)
    !% Set the node spheroid star formation history.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(history),           intent(in)    :: thisHistory
    integer                                :: thisIndex

    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    call thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%destroy()
    thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)=thisHistory
    return
  end subroutine Tree_Node_Spheroid_Star_Formation_History_Set_Sersic

  subroutine Tree_Node_Spheroid_Stellar_Prprts_History_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the rates for the stellar properties history.
    use Histories
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    logical,                  intent(inout) :: interrupt
    procedure(),     pointer, intent(inout) :: interruptProcedure
    type(history),            intent(in)    :: rateAdjustments
    integer                                 :: thisIndex

    ! If no Sérsic spheroid component currently exists and we have some non-zero rate into it then interrupt and create an Sérsic
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments%data /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Sersic_Spheroid_Create
       end if
       return
    end if
    ! Get the index for this component.
    thisIndex=thisNode%componentIndex(componentIndex)
    ! Adjust the rate.
    call thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)%add(rateAdjustments,addTo=historyRates)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Prprts_History_Rate_Adjust_Sersic

  subroutine Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Sersic(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the rates for the star formation history.
    use Histories
    use Galacticus_Error
    use Memory_Management
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    logical,                  intent(inout) :: interrupt
    procedure(),     pointer, intent(inout) :: interruptProcedure
    type(history),            intent(in)    :: rateAdjustments
    integer                                 :: thisIndex

    ! If no Sérsic spheroid component currently exists and we have some non-zero rate into it then interrupt and create an Sérsic
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments%data /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Sersic_Spheroid_Create
       end if
       return
    end if
    ! Get the index for this component.
    thisIndex=thisNode%componentIndex(componentIndex)
    ! Ensure that a history already exists in the spheroid.
    if (.not.thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%exists())                &
         & call Galacticus_Error_Report('Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Sersic' &
         & ,'no star formation history has been created in spheroid')
    ! Check if the star formation history in the spheroid spans a sufficient range to accept the input rates.
    if (rateAdjustments%time(1) < thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%time(1) .or.                      &
         & rateAdjustments%time(size(rateAdjustments%time)) > thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%time( &
         & size(thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%time))) then
       ! It does not, so interrupt evolution and extend the history.
       if (allocated(starFormationHistoryTemplate)) call Dealloc_Array(starFormationHistoryTemplate)
       call Alloc_Array(starFormationHistoryTemplate,shape(rateAdjustments%time))
       starFormationHistoryTemplate=rateAdjustments%time
       interrupt=.true.
       interruptProcedure => Sersic_Spheroid_Star_Formation_History_Extend
       return
    end if

    ! Adjust the rate.
    call thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%combine(rateAdjustments,addTo=historyRates)
    return
  end subroutine Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Sersic

  !# <scaleSetTask>
  !#  <unitName>Sersic_Spheroid_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Sersic_Spheroid_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}. Note that gas masses get an additional scaling down since they can approach
    !% zero and we'd like to prevent them from becoming negative.
    use Abundances_Structure
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(treeNode),            pointer, intent(inout) :: thisNode
    double precision,          parameter              :: massMinimum           =1.0d0
    double precision,          parameter              :: angularMomentumMinimum=0.1d0
    double precision,          parameter              :: luminosityMinimum     =1.0d0
    double precision,          parameter              :: gasMassScaling        =0.1d0
    type(abundancesStructure), save                   :: stellarAbundances
    !$omp threadprivate(stellarAbundances)
    integer                                           :: thisIndex
    double precision                                  :: mass,angularMomentum

    ! Determine if method is active and a spheroid component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)

       ! Set scale for angular momentum.
       angularMomentum=Tree_Node_Spheroid_Angular_Momentum(thisNode)+Tree_Node_Disk_Angular_Momentum(thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyScale)=max(angularMomentum,angularMomentumMinimum)

       ! Set scale for gas mass.
       mass=Tree_Node_Spheroid_Gas_Mass(thisNode)+Tree_Node_Disk_Gas_Mass(thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyScale)=gasMassScaling*max(mass,massMinimum)

       ! Set scale for stellar mass.
       mass=Tree_Node_Spheroid_Stellar_Mass(thisNode)+Tree_Node_Disk_Stellar_Mass(thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyScale)=max(mass,massMinimum)

       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for gas abundances.
          call Tree_Node_Spheroid_Gas_Abundances_Sersic(thisNode,abundancesSpheroid)
          call Tree_Node_Disk_Gas_Abundances              (thisNode,abundancesDisk    )
          thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyScale)=gasMassScaling*max(abundancesDisk+abundancesSpheroid,massMinimum)
          
          ! Set scale for stellar abundances.
          call Tree_Node_Spheroid_Stellar_Abundances_Sersic(thisNode,abundancesSpheroid)
          call Tree_Node_Disk_Stellar_Abundances              (thisNode,abundancesDisk    )
          thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyScale)=max(abundancesDisk+abundancesSpheroid,massMinimum)
       end if

       ! Set scales for stellar luminosities if necessary.
       if (luminositiesCount > 0) then        
          ! Set scale for stellar luminosities.
          call Tree_Node_Spheroid_Stellar_Luminosities_Sersic(thisNode,luminositiesSpheroid)
          call Tree_Node_Disk_Stellar_Luminosities              (thisNode,luminositiesDisk    )
          thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyScale)=max(luminositiesDisk+luminositiesSpheroid,luminosityMinimum)
       end if

       ! Set scales for stellar population properties history.
       call Tree_Node_Spheroid_Stellar_Abundances_Sersic(thisNode,abundancesSpheroid)
       call stellarAbundances%pack(abundancesSpheroid)
       call Stellar_Population_Properties_Scales(thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex      ),Tree_Node_Spheroid_Stellar_Mass(thisNode),stellarAbundances)
       call Star_Formation_History_Scales       (thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex),Tree_Node_Spheroid_Stellar_Mass(thisNode),stellarAbundances)
    end if
    return
  end subroutine Sersic_Spheroid_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Sersic_Spheroid_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Sersic_Spheroid_Satellite_Merging(thisNode)
    !% Transfer any Sérsic spheroid associated with {\tt thisNode} to its host halo.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    use Satellite_Merging_Remnant_Sizes_Properties
    implicit none
    type(treeNode),   pointer, intent(inout)       :: thisNode
    type(treeNode),   pointer                      :: hostNode
    double precision, dimension(abundancesCount)   :: thisAbundances,hostAbundances
    double precision, dimension(luminositiesCount) :: thisLuminosities,hostLuminosities
    type(history)                                  :: historyDisk,historySpheroid,thisHistory
    double precision                               :: spheroidSpecificAngularMomentum,diskSpecificAngularMomentum,angularMomentum&
         &,spheroidMass

    ! Check that method is selected.
    if (methodSelected) then
       
       ! Find the node to merge with.
       call thisNode%mergesWith(hostNode)

       ! Get specific angular momentum of the host spheroid and disk material.
       if (Tree_Node_Spheroid_Gas_Mass_Sersic(hostNode)+Tree_Node_Spheroid_Stellar_Mass_Sersic(hostNode) > 0.0d0) then
          spheroidSpecificAngularMomentum=Tree_Node_Spheroid_Angular_Momentum_Sersic(hostNode)&
               &/(Tree_Node_Spheroid_Gas_Mass_Sersic(hostNode)+Tree_Node_Spheroid_Stellar_Mass_Sersic(hostNode))
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
               &                                                        +Tree_Node_Spheroid_Gas_Mass_Sersic(hostNode))
          call Tree_Node_Disk_Gas_Abundances                  (hostNode, hostAbundances                                 )
          call Tree_Node_Spheroid_Gas_Abundances_Sersic    (hostNode, thisAbundances                                 )
          call Tree_Node_Disk_Gas_Abundances_Set              (hostNode, hostAbundances+thisAbundances                  )
          call Tree_Node_Disk_Angular_Momentum_Set            (hostNode, Tree_Node_Disk_Angular_Momentum(hostNode)       &
               &+spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Gas_Mass_Sersic(hostNode)                        )
          call Tree_Node_Spheroid_Angular_Momentum_Set        (hostNode, Tree_Node_Spheroid_Angular_Momentum(hostNode)   &
               &-spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Gas_Mass_Sersic(hostNode)                        )
          call Tree_Node_Spheroid_Gas_Mass_Set_Sersic      (hostNode, 0.0d0                                          )
          thisAbundances=0.0d0
          call Tree_Node_Spheroid_Gas_Abundances_Set_Sersic(hostNode, thisAbundances                                 )
       case (movesToSpheroid)
          call Tree_Node_Spheroid_Gas_Mass_Set_Sersic      (hostNode, Tree_Node_Spheroid_Gas_Mass_Sersic(hostNode) &
               &                                                        +Tree_Node_Disk_Gas_Mass              (hostNode))
          call Tree_Node_Disk_Angular_Momentum_Set            (hostNode, Tree_Node_Disk_Angular_Momentum(hostNode)       &
               &-diskSpecificAngularMomentum*Tree_Node_Disk_Gas_Mass(hostNode)                                          )
          call Tree_Node_Spheroid_Gas_Abundances_Sersic    (hostNode,hostAbundances                                  )
          call Tree_Node_Disk_Gas_Abundances                  (hostNode,thisAbundances                                  )
          call Tree_Node_Spheroid_Gas_Abundances_Set_Sersic(hostNode,hostAbundances+thisAbundances                   )
          call Tree_Node_Disk_Gas_Mass_Set                    (hostNode,0.0d0                                           )
          thisAbundances=0.0d0
          call Tree_Node_Disk_Gas_Abundances_Set              (hostNode,thisAbundances                                  )
       case (doesNotMove)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('Sersic_Spheroid_Satellite_Merging','unrecognized movesTo descriptor')
       end select

       ! Move stellar material within the host if necessary.
       select case (thisHostStarsMoveTo)
       case (movesToDisk)
          call Tree_Node_Disk_Stellar_Mass_Set                      (hostNode, Tree_Node_Disk_Stellar_Mass              (hostNode) &
               &                                                              +Tree_Node_Spheroid_Stellar_Mass_Sersic(hostNode))
          call Tree_Node_Disk_Stellar_Abundances                    (hostNode, hostAbundances                                     )
          call Tree_Node_Spheroid_Stellar_Abundances_Sersic      (hostNode, thisAbundances                                     )
          call Tree_Node_Disk_Stellar_Abundances_Set                (hostNode, hostAbundances+thisAbundances                      )
          call Tree_Node_Disk_Stellar_Luminosities                  (hostNode, hostLuminosities                                   )
          call Tree_Node_Spheroid_Stellar_Luminosities_Sersic    (hostNode, thisLuminosities                                   )
          call Tree_Node_Disk_Stellar_Luminosities_Set              (hostNode, hostLuminosities+thisLuminosities                  )
          call Tree_Node_Disk_Angular_Momentum_Set                  (hostNode, Tree_Node_Disk_Angular_Momentum(hostNode)           &
               &+spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Stellar_Mass_Sersic(hostNode)                              )
          call Tree_Node_Spheroid_Angular_Momentum_Set              (hostNode, Tree_Node_Spheroid_Angular_Momentum(hostNode)        &
               &-spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Stellar_Mass_Sersic(hostNode)                              )
          call Tree_Node_Spheroid_Stellar_Mass_Set_Sersic        (hostNode, 0.0d0                                              )
          thisAbundances=0.0d0
          call Tree_Node_Spheroid_Stellar_Abundances_Set_Sersic  (hostNode, thisAbundances                                     )
          thisLuminosities=0.0d0
          call Tree_Node_Spheroid_Stellar_Luminosities_Set_Sersic(hostNode, thisLuminosities                                   )
          ! Also add stellar properties histories.
          historyDisk    =Tree_Node_Disk_Stellar_Properties_History              (hostNode                )
          historySpheroid=Tree_Node_Spheroid_Stellar_Properties_History_Sersic(hostNode                )
          call historyDisk%add(historySpheroid)
          call Tree_Node_Disk_Stellar_Properties_History_Set                     (hostNode,historyDisk    )
          call historySpheroid%reset()
          call Tree_Node_Spheroid_Stellar_Properties_History_Set_Sersic       (hostNode,historySpheroid)
          ! Also add stellar properties histories.
          historyDisk    =Tree_Node_Disk_Star_formation_History                  (hostNode                )
          historySpheroid=Tree_Node_Spheroid_Star_formation_History_Sersic    (hostNode                )
          call historyDisk%combine(historySpheroid)
          call Tree_Node_Disk_Star_formation_History_Set                         (hostNode,historyDisk    )
          call historySpheroid%reset()
          call Tree_Node_Spheroid_Star_formation_History_Set_Sersic           (hostNode,historySpheroid)
          call historyDisk    %destroy(recordMemory=.false.)
          call historySpheroid%destroy(recordMemory=.false.)
       case (movesToSpheroid)
          call Tree_Node_Spheroid_Stellar_Mass_Set_Sersic        (hostNode, Tree_Node_Spheroid_Stellar_Mass_Sersic(hostNode) &
               &                                                              +Tree_Node_Disk_Stellar_Mass              (hostNode))
          call Tree_Node_Disk_Angular_Momentum_Set                  (hostNode, Tree_Node_Disk_Angular_Momentum(hostNode)           &
               &-diskSpecificAngularMomentum*Tree_Node_Disk_Stellar_Mass(hostNode)                                                )
          call Tree_Node_Spheroid_Stellar_Abundances_Sersic      (hostNode, hostAbundances                                     )
          call Tree_Node_Disk_Stellar_Abundances                    (hostNode, thisAbundances                                     )
          call Tree_Node_Spheroid_Stellar_Abundances_Set_Sersic  (hostNode, hostAbundances+thisAbundances                      )
          call Tree_Node_Spheroid_Stellar_Luminosities_Sersic    (hostNode, hostLuminosities                                   )
          call Tree_Node_Disk_Stellar_Luminosities                  (hostNode, thisLuminosities                                   )
          call Tree_Node_Spheroid_Stellar_Luminosities_Set_Sersic(hostNode, hostLuminosities+thisLuminosities                  )
          call Tree_Node_Disk_Stellar_Mass_Set                      (hostNode, 0.0d0                                              )
          thisAbundances=0.0d0
          call Tree_Node_Disk_Stellar_Abundances_Set                (hostNode, thisAbundances                                     )
          thisLuminosities=0.0d0
          call Tree_Node_Disk_Stellar_Luminosities_Set              (hostNode, thisLuminosities                                   )
          ! Also add stellar properties histories.
          historyDisk    =Tree_Node_Disk_Stellar_Properties_History              (hostNode                )
          historySpheroid=Tree_Node_Spheroid_Stellar_Properties_History_Sersic(hostNode                )
          call historySpheroid%add(historyDisk)
          call Tree_Node_Spheroid_Stellar_Properties_History_Set_Sersic       (hostNode,historySpheroid)
          call historyDisk%reset()
          call Tree_Node_Disk_Stellar_Properties_History_Set                     (hostNode,historyDisk    )
          ! Also add stellar properties histories.
          historyDisk    =Tree_Node_Disk_Star_formation_History                  (hostNode                )
          historySpheroid=Tree_Node_Spheroid_Star_formation_History_Sersic    (hostNode                )
          call historySpheroid%combine(historyDisk)
          call Tree_Node_Spheroid_Star_formation_History_Set_Sersic           (hostNode,historySpheroid)
          call historyDisk%reset()
          call Tree_Node_Disk_Star_formation_History_Set                         (hostNode,historyDisk    )
          call historyDisk    %destroy(recordMemory=.false.)
          call historySpheroid%destroy(recordMemory=.false.)
       case (doesNotMove)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('Sersic_Spheroid_Satellite_Merging','unrecognized movesTo descriptor')
       end select
   
       ! Get specific angular momentum of the spheroid material.
       spheroidMass=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
       if (spheroidMass > 0.0d0) then
          spheroidSpecificAngularMomentum=Tree_Node_Spheroid_Angular_Momentum_Sersic(thisNode)/spheroidMass

          ! Move the gas component of the Sersic spheroid to the host.
          select case (thisMergerGasMovesTo)
          case (movesToDisk)
             call Tree_Node_Disk_Gas_Mass_Set                (hostNode, Tree_Node_Disk_Gas_Mass              (hostNode) &
                  &                                                    +Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode))
             call Tree_Node_Disk_Gas_Abundances              (hostNode,hostAbundances)
             call Tree_Node_Spheroid_Gas_Abundances_Sersic(thisNode,thisAbundances)
             call Tree_Node_Disk_Gas_Abundances_Set          (hostNode,hostAbundances+thisAbundances)
             call Tree_Node_Disk_Angular_Momentum_Set        (hostNode,Tree_Node_Disk_Angular_Momentum(hostNode) &
                  &+spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode))
          case (movesToSpheroid)
             call Tree_Node_Spheroid_Gas_Mass_Set_Sersic  (hostNode, Tree_Node_Spheroid_Gas_Mass_Sersic(hostNode) &
                  &                                                    +Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode))
             call Tree_Node_Spheroid_Gas_Abundances_Sersic(hostNode,hostAbundances)
             call Tree_Node_Spheroid_Gas_Abundances_Sersic(thisNode,thisAbundances)
             call Tree_Node_Spheroid_Gas_Abundances_Set_Sersic(hostNode,hostAbundances+thisAbundances)
          case default
             call Galacticus_Error_Report('Sersic_Spheroid_Satellite_Merging','unrecognized movesTo descriptor')
          end select
          call Tree_Node_Spheroid_Gas_Mass_Set_Sersic      (thisNode,0.0d0         )
          thisAbundances=0.0d0
          call Tree_Node_Spheroid_Gas_Abundances_Set_Sersic(thisNode,thisAbundances)
          
          ! Move the stellar component of the Sérsic spheroid to the host.
          select case (thisMergerStarsMoveTo)
          case (movesToDisk)
             call Tree_Node_Disk_Stellar_Mass_Set                  (hostNode, Tree_Node_Disk_Stellar_Mass(hostNode) &
                  &                                                          +Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode))
             call Tree_Node_Disk_Stellar_Abundances                (hostNode,hostAbundances)
             call Tree_Node_Spheroid_Stellar_Abundances_Sersic  (thisNode,thisAbundances)
             call Tree_Node_Disk_Stellar_Abundances_Set            (hostNode,hostAbundances+thisAbundances)
             call Tree_Node_Disk_Stellar_Luminosities              (hostNode,hostLuminosities)
             call Tree_Node_Spheroid_Stellar_Luminosities_Sersic(thisNode,thisLuminosities)
             call Tree_Node_Disk_Stellar_Luminosities_Set          (hostNode,hostLuminosities+thisLuminosities)
             call Tree_Node_Disk_Angular_Momentum_Set              (hostNode,Tree_Node_Disk_Angular_Momentum(hostNode) &
                  &+spheroidSpecificAngularMomentum*Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode))
             ! Also add stellar properties histories.
             historySpheroid=Tree_Node_Spheroid_Stellar_Properties_History_Sersic(thisNode                )
             thisHistory    =Tree_Node_Disk_Stellar_Properties_History              (hostNode                )
             call thisHistory%add(historySpheroid)
             call Tree_Node_Disk_Stellar_Properties_History_Set                     (hostNode,thisHistory    )
             call historySpheroid%reset()
             call Tree_Node_Spheroid_Stellar_Properties_History_Set_Sersic       (thisNode,historySpheroid)
             ! Also add star formation histories.
             historySpheroid=Tree_Node_Spheroid_Star_formation_History_Sersic    (thisNode                )
             thisHistory    =Tree_Node_Disk_Star_formation_History                  (hostNode                )
             call thisHistory%combine(historySpheroid)
             call Tree_Node_Disk_Star_formation_History_Set                         (hostNode,thisHistory    )
             call historySpheroid%reset()
             call Tree_Node_Spheroid_Star_formation_History_Set_Sersic           (thisNode,historySpheroid)
             call thisHistory    %destroy(recordMemory=.false.)
             call historySpheroid%destroy(recordMemory=.false.)
         case (movesToSpheroid)
             call Tree_Node_Spheroid_Stellar_Mass_Set_Sersic        (hostNode, Tree_Node_Spheroid_Stellar_Mass_Sersic(hostNode) &
                  &                                                              +Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode))
             call Tree_Node_Spheroid_Stellar_Abundances_Sersic      (hostNode,hostAbundances)
             call Tree_Node_Spheroid_Stellar_Abundances_Sersic      (thisNode,thisAbundances)
             call Tree_Node_Spheroid_Stellar_Abundances_Set_Sersic  (hostNode,hostAbundances+thisAbundances)
             call Tree_Node_Spheroid_Stellar_Luminosities_Sersic    (hostNode,hostLuminosities)
             call Tree_Node_Spheroid_Stellar_Luminosities_Sersic    (thisNode,thisLuminosities)
             call Tree_Node_Spheroid_Stellar_Luminosities_Set_Sersic(hostNode,hostLuminosities+thisLuminosities)
             ! Also add stellar properties histories.
             historySpheroid=Tree_Node_Spheroid_Stellar_Properties_History_Sersic(thisNode                )
             thisHistory    =Tree_Node_Spheroid_Stellar_Properties_History_Sersic(hostNode                )
             call thisHistory%add(historySpheroid)
             call Tree_Node_Spheroid_Stellar_Properties_History_Set_Sersic       (hostNode,thisHistory    )
             call historySpheroid%reset()
             call Tree_Node_Spheroid_Stellar_Properties_History_Set_Sersic       (thisNode,historySpheroid)
             ! Also add star formation histories.
             historySpheroid=Tree_Node_Spheroid_Star_formation_History_Sersic    (thisNode                )
             thisHistory    =Tree_Node_Spheroid_Star_formation_History_Sersic    (hostNode                )
             call thisHistory%combine(historySpheroid)
             call Tree_Node_Spheroid_Star_formation_History_Set_Sersic           (hostNode,thisHistory    )
             call historySpheroid%reset()
             call Tree_Node_Spheroid_Star_formation_History_Set_Sersic           (thisNode,historySpheroid)
             call thisHistory    %destroy(recordMemory=.false.)
             call historySpheroid%destroy(recordMemory=.false.)
          case default
             call Galacticus_Error_Report('Sersic_Spheroid_Satellite_Merging','unrecognized movesTo descriptor')
          end select
          call Tree_Node_Spheroid_Stellar_Mass_Set_Sersic        (thisNode,0.0d0           )
          thisAbundances=0.0d0
          call Tree_Node_Spheroid_Stellar_Abundances_Set_Sersic  (thisNode,thisAbundances  )
          thisLuminosities=0.0d0
          call Tree_Node_Spheroid_Stellar_Luminosities_Set_Sersic(thisNode,thisLuminosities)       
          call Tree_Node_Spheroid_Angular_Momentum_Set_Sersic    (thisNode,0.0d0           )
       end if
       
       ! Set the angular momentum of the spheroid.
       if (remnantSpecificAngularMomentum /= remnantNoChangeValue) then
          ! Note that the remnant specific angular momentum computed by the merger remnant modules automatically gives the mean
          ! specific angular momentum of the component by virtue of the fact that it computes the ratio of the actual angular
          ! momentum to the contribution from the component's own rotation curve at its scale radius.
          angularMomentum=remnantSpecificAngularMomentum*(Tree_Node_Spheroid_Gas_Mass(hostNode) &
               &+Tree_Node_Spheroid_Stellar_Mass(hostNode))
          call Tree_Node_Spheroid_Angular_Momentum_Set_Sersic(hostNode,angularMomentum)
       end if

    end if
    return
  end subroutine Sersic_Spheroid_Satellite_Merging
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Sersic_Profile_Tabulate_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Sersic_Profile_Tabulate_State_Store(stateFile,fgslStateFile)
    !% Write the history state to file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) sersicTableRadiusMinimum,sersicTableRadiusMaximum,sersicTable3dHalfMassRadius
    return
  end subroutine Sersic_Profile_Tabulate_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Sersic_Profile_Tabulate_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Sersic_Profile_Tabulate_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the history state from the file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) sersicTableRadiusMinimum,sersicTableRadiusMaximum,sersicTable3dHalfMassRadius
    sersicTableInitialized=.false.
    return
  end subroutine Sersic_Profile_Tabulate_State_Retrieve

  subroutine Sersic_Profile_Tabulate(radius)
    !% Tabulate the density and enclosed mass in a dimensionless Sérsic profile.
    use, intrinsic :: ISO_C_Binding                             
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Integration
    use Numerical_Interpolation
    use Numerical_Constants_Math         
    use FGSL
    implicit none
    double precision,                intent(in), optional :: radius
    double precision,                parameter            :: radiusMaximumToTabulate=1000.0d0
    logical                                               :: rebuildTable,tableHasSufficientExtent
    integer                                               :: iRadius
    double precision                                      :: deltaRadius,massPrevious,radiusInfinity,integrand,previousIntegrand &
         &,angularMomentum,integrandAngularMomentum,previousIntegrandAngularMomentum,angularMomentumMean,radiusActual
    type(c_ptr)                                           :: parameterPointer
    type(fgsl_function)                                   :: integrandFunction
    type(fgsl_integration_workspace)                      :: integrationWorkspace

    ! Check if a radius was specified. Use it if so, otherwise use a midpoint radius.
    if (present(radius)) then
       radiusActual=min(radius,radiusMaximumToTabulate)
    else
       radiusActual=dsqrt(sersicTableRadiusMinimum*sersicTableRadiusMaximum)
    end if

    !$omp critical (Sersic_Profile_Tabulation)
    ! Determine if the table must be rebuilt.
    if (sersicTableInitialized) then
       rebuildTable=(radiusActual*sersicTable3dHalfMassRadius < sersicTableRadiusMinimum) &
            &        .or.                               &
            &       (radiusActual*sersicTable3dHalfMassRadius > sersicTableRadiusMaximum)
    else
       rebuildTable=.true.
    end if

    ! Rebuild the table if necessary.
    if (rebuildTable) then
       ! Try building the table until it has sufficient extent to encompass the requested radius.
       tableHasSufficientExtent=.false.
       do while (.not.tableHasSufficientExtent)
          ! Find suitable radius limits.
          sersicTableRadiusMinimum=min(sersicTableRadiusMinimum,0.5d0*radiusActual*sersicTable3dHalfMassRadius)
          sersicTableRadiusMaximum=max(sersicTableRadiusMaximum,2.0d0*radiusActual*sersicTable3dHalfMassRadius)
          ! Determine the number of points at which to tabulate the profile.
          sersicTableCount=int(dlog10(sersicTableRadiusMaximum/sersicTableRadiusMinimum)*dble(sersicTablePointsPerDecade))+1
          ! Allocate arrays for storing the tables.
          if (allocated(sersicTableRadius)) then
             call Dealloc_Array(sersicTableRadius      )
             call Dealloc_Array(sersicTableDensity     )
             call Dealloc_Array(sersicTableEnclosedMass)
          end if
          call Alloc_Array(sersicTableRadius      ,[sersicTableCount])
          call Alloc_Array(sersicTableDensity     ,[sersicTableCount])
          call Alloc_Array(sersicTableEnclosedMass,[sersicTableCount])
          ! Create an array of logarithmically distributed radii.
          sersicTableRadius=Make_Range(sersicTableRadiusMinimum,sersicTableRadiusMaximum,sersicTableCount,rangeType=rangeTypeLogarithmic)
          ! Compute the coefficient appearing in the Sérsic profile.
          spheroidSersicCoefficient=2.303d0*(0.8689d0*dble(spheroidSersicIndex)-0.1447d0) ! Wadadekar et al. (1999; AJ; 117; 1219)
          ! Compute a suitably large approximation to infinite radius for use in integration.
          radiusInfinity=10.0d0*sersicTableRadiusMaximum
          ! Initialize the total angular momentum of the profile to zero.
          angularMomentum=0.0d0
          ! Loop over radii and compute the inverse Abel integral required to get the 3D Sérsic profile.
          do iRadius=1,sersicTableCount
             radiusStart=sersicTableRadius(iRadius)
             sersicTableDensity(iRadius)=Integrate( radiusStart,radiusInfinity                       &
                  &                                ,Sersic_Abel_Integrand,parameterPointer           &
                  &                                ,integrandFunction,integrationWorkspace           &
                  &                                ,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3 &
                  &                               )
             call Integrate_Done(integrandFunction,integrationWorkspace)       
             ! Accumulate the enclosed mass using a simple trapezoidal integration.
             if (iRadius == 1) then
                deltaRadius                     =sersicTableRadius(iRadius)
                massPrevious                    =0.0d0
                previousIntegrand               =0.0d0
                previousIntegrandAngularMomentum=0.0d0
             else
                deltaRadius                     =sersicTableRadius(iRadius)-sersicTableRadius      (iRadius-1)
                massPrevious                    =                           sersicTableEnclosedMass(iRadius-1)
                previousIntegrand               =4.0d0*Pi*(sersicTableRadius(iRadius-1)**2)*sersicTableDensity(iRadius-1)
                previousIntegrandAngularMomentum=4.0d0*Pi*(sersicTableRadius(iRadius-1)**3)*sersicTableDensity(iRadius-1)
             end if
             integrand               =4.0d0*Pi*(sersicTableRadius(iRadius)**2)*sersicTableDensity(iRadius)
             integrandAngularMomentum=4.0d0*Pi*(sersicTableRadius(iRadius)**3)*sersicTableDensity(iRadius)
             sersicTableEnclosedMass(iRadius)=0.5d0*(previousIntegrand               +integrand               )*deltaRadius+massPrevious
             angularMomentum                 =0.5d0*(previousIntegrandAngularMomentum+integrandAngularMomentum)*deltaRadius+angularMomentum
          end do
          ! Compute the mean angular momentum of the profile.
          angularMomentumMean=angularMomentum/sersicTableEnclosedMass(sersicTableCount)
          ! Normalize the mass and density to unit total mass.
          sersicTableDensity     =sersicTableDensity     /sersicTableEnclosedMass(sersicTableCount)
          sersicTableEnclosedMass=sersicTableEnclosedMass/sersicTableEnclosedMass(sersicTableCount)
          ! Find the half mass radius.
          call Interpolate_Done(sersicTableInterpolationObject,sersicTableInterpolationAccelerator,sersicTableInterpolationReset)
          sersicTableInterpolationReset=.true.
          sersicTable3dHalfMassRadius=Interpolate(                                                                  &
               &                                                            maxloc(sersicTableEnclosedMass,dim=1),  &
               &                                  sersicTableEnclosedMass(1:maxloc(sersicTableEnclosedMass,dim=1)), &
               &                                  sersicTableRadius      (1:maxloc(sersicTableEnclosedMass,dim=1)), &
               &                                  sersicTableInterpolationObject,                                   &
               &                                  sersicTableInterpolationAccelerator,                              &
               &                                  0.5d0,                                                            &
               &                                  reset=sersicTableInterpolationReset                               &
               &                                 )
          call Interpolate_Done(sersicTableInterpolationObject,sersicTableInterpolationAccelerator,sersicTableInterpolationReset)
          sersicTableInterpolationReset=.true.
          ! Compute the ratio of angular momentum at the half-mass radius to the mean value.
          sersicTableAngularMomentumHalfMassRadiusToMeanRatio=sersicTable3dHalfMassRadius/angularMomentumMean
          ! Scale radii and densities to be in units of the 3D half mass radius.
          sersicTableRadius =sersicTableRadius /sersicTable3dHalfMassRadius
          sersicTableDensity=sersicTableDensity*sersicTable3dHalfMassRadius**3
          ! Test that the table has sufficient extent for the requested radius.
          tableHasSufficientExtent=(radiusActual >= sersicTableRadiusMinimum) & 
               &                    .and. &
               &                   (radiusActual <= sersicTableRadiusMaximum)
       end do
       ! Flag that the table is initialized.
       sersicTableInitialized=.true.
    end if
    !$omp end critical (Sersic_Profile_Tabulation)
   return
  end subroutine Sersic_Profile_Tabulate

  function Sersic_Abel_Integrand(radius,parameterPointer) bind(c)
    !% The integrand in the Abel integral used to invert the Sérsic profile to get the corresponding 3-D profile.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    implicit none
    real(c_double)        :: Sersic_Abel_Integrand
    real(c_double), value :: radius
    type(c_ptr)           :: parameterPointer

    if (radius > radiusStart) then
       Sersic_Abel_Integrand= (1.0d0/Pi)                                                                         &
            &                *      spheroidSersicCoefficient*(radius**(1.0d0/dble(spheroidSersicIndex) -1.0d0)) &
            &                *dexp(-spheroidSersicCoefficient*(radius**(1.0d0/dble(spheroidSersicIndex))-1.0d0)) &
            &                /dble(spheroidSersicIndex)                                                          &
            &                /dsqrt(radius**2-radiusStart**2)
    else
       Sersic_Abel_Integrand=0.0d0
    end if
    return
  end function Sersic_Abel_Integrand

  !# <enclosedMassTask>
  !#  <unitName>Sersic_Spheroid_Enclosed_Mass</unitName>
  !# </enclosedMassTask>
  subroutine Sersic_Spheroid_Enclosed_Mass(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for an Sérsic spheroid.
    use Galactic_Structure_Options
    use Numerical_Interpolation
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType,weightBy,weightIndex
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass
    double precision                         :: fractionalRadius,spheroidRadius
    
    componentMass=0.0d0
    if (radius <= 0.0d0                               ) return
    if (.not.methodSelected                           ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
    if (.not.thisNode%componentExists(componentIndex)) return
    select case (weightBy)
    case (weightByMass      )
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic)
          componentMass=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
       case (massTypeGaseous)
          componentMass=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)
       case (massTypeStellar)
          componentMass=Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
       end select
    case (weightByLuminosity)
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
          call Tree_Node_Spheroid_Stellar_Luminosities_Sersic(thisNode,luminositiesSpheroid)
          componentMass=luminositiesSpheroid(weightIndex)
       end select
    end select
    ! Return if total mass was requested.  
    if (radius >= radiusLarge)                          return
    ! Return if mass is zero.
    if (componentMass <= 0.0d0)                         return
    ! Compute actual mass.
    spheroidRadius=Sersic_Spheroid_Radius(thisNode)
    if (spheroidRadius > 0.0d0) then
       fractionalRadius=radius/spheroidRadius
       call Sersic_Profile_Tabulate(fractionalRadius)
       ! Scale by the fractional enclosed mass if the radius is less than the largest radius tabulated.
       if (fractionalRadius < sersicTableRadius(sersicTableCount)) then
          componentMass=componentMass*Interpolate(sersicTableCount,sersicTableRadius,sersicTableEnclosedMass &
               &,sersicTableInterpolationObject,sersicTableInterpolationAccelerator,fractionalRadius,reset&
               &=sersicTableInterpolationReset)
       end if
    end if
    return
  end subroutine Sersic_Spheroid_Enclosed_Mass

  !# <rotationCurveTask>
  !#  <unitName>Sersic_Spheroid_Rotation_Curve</unitName>
  !# </rotationCurveTask>
  subroutine Sersic_Spheroid_Rotation_Curve(thisNode,radius,massType,componentType,componentVelocity)
    !% Computes the rotation curve at a given radius for an Sérsic spheroid.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentVelocity
    double precision                         :: componentMass

    ! Set to zero by default.
    componentVelocity=0.0d0

    ! Compute if a spheroid is present.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then       
       if (radius > 0.0d0) then
          call Sersic_Spheroid_Enclosed_Mass(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
          if (componentMass > 0.0d0) componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass)/dsqrt(radius)
       end if
    end if
    return
  end subroutine Sersic_Spheroid_Rotation_Curve

  !# <densityTask>
  !#  <unitName>Sersic_Spheroid_Density</unitName>
  !# </densityTask>
  subroutine Sersic_Spheroid_Density(thisNode,positionSpherical,massType,componentType,componentDensity)
    !% Computes the density at a given position for an Sérsic spheroid.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    use Numerical_Interpolation
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
    if (Sersic_Spheroid_Radius(thisNode) <= 0.0d0) return
    select case (massType)
    case (massTypeAll,massTypeBaryonic,massTypeGalactic)
       componentDensity=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
    case (massTypeGaseous)
       componentDensity=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)
    case (massTypeStellar)
       componentDensity=Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
    end select
    ! Return if density is zero.
    if (componentDensity <= 0.0d0) then
       componentDensity=0.0d0
       return
    end if
    ! Compute actual density.
    fractionalRadius=positionSpherical(1)/Sersic_Spheroid_Radius(thisNode)
    call Sersic_Profile_Tabulate(fractionalRadius)
    if (fractionalRadius < sersicTableRadius(sersicTableCount)) then
       componentDensity=componentDensity*Interpolate(sersicTableCount,sersicTableRadius,sersicTableDensity &
            &,sersicTableInterpolationObject,sersicTableInterpolationAccelerator,fractionalRadius,reset&
            &=sersicTableInterpolationReset)/Sersic_Spheroid_Radius(thisNode)**3
    else
       ! Assume zero density beyond the largest radius tabulated.
       componentDensity=0.0d0
    end if
    return
  end subroutine Sersic_Spheroid_Density

  !# <radiusSolverPlausibility>
  !#  <unitName>Sersic_Spheroid_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Sersic_Spheroid_Radius_Solver_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
    !% Determines whether the spheroid is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: galaxyIsPhysicallyPlausible

    ! Return immediately if our method is not selected.
    if (.not.methodSelected) return

    if (Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)+Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode) < -spheroidMassToleranceAbsolute) then
       galaxyIsPhysicallyPlausible=.false.
    else
       if (Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)+Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode) >&
            & spheroidMassToleranceAbsolute .and. Tree_Node_Spheroid_Angular_Momentum_Sersic(thisNode) < 0.0d0)&
            & galaxyIsPhysicallyPlausible=.false.
    end if

    return
  end subroutine Sersic_Spheroid_Radius_Solver_Plausibility

  !# <radiusSolverTask>
  !#  <unitName>Sersic_Spheroid_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Sersic_Spheroid_Radius_Solver(thisNode,componentActive,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !% Interface for the size solver algorithm.
    implicit none
    type(treeNode),              pointer, intent(inout) :: thisNode
    logical,                              intent(out)   :: componentActive
    double precision,                     intent(out)   :: specificAngularMomentum
    procedure(double precision), pointer, intent(out)   :: Radius_Get,Velocity_Get
    procedure(),                 pointer, intent(out)   :: Radius_Set,Velocity_Set
    double precision                                    :: specificAngularMomentumMean,angularMomentum,spheroidMass

    ! Determine if thisNode has an active spheroid component supported by this module.    
    componentActive=methodSelected

    if (methodSelected) componentActive=thisNode%componentExists(componentIndex)    
    if (componentActive) then
       ! Get the angular momentum.
       angularMomentum=Tree_Node_Spheroid_Angular_Momentum_Sersic(thisNode)
       if (angularMomentum >= 0.0d0) then
          ! Compute the specific angular momentum at the scale radius, assuming a flat rotation curve.
          spheroidMass=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
          if (spheroidMass > 0.0d0) then
             specificAngularMomentumMean=angularMomentum/spheroidMass
          else
             specificAngularMomentumMean=0.0d01
          end if
          call Sersic_Profile_Tabulate()
          specificAngularMomentum=sersicTableAngularMomentumHalfMassRadiusToMeanRatio*specificAngularMomentumMean
          ! Associate the pointers with the appropriate property routines.
          Radius_Get   => Sersic_Spheroid_Radius
          Radius_Set   => Sersic_Spheroid_Radius_Set
          Velocity_Get => Sersic_Spheroid_Velocity
          Velocity_Set => Sersic_Spheroid_Velocity_Set
      else
          call Sersic_Spheroid_Radius_Set  (thisNode,0.0d0)
          call Sersic_Spheroid_Velocity_Set(thisNode,0.0d0)
          componentActive=.false.
       end if
    end if
    return
  end subroutine Sersic_Spheroid_Radius_Solver

  double precision function Sersic_Spheroid_SFR(thisNode,instance)
    !% Return the star formation rate of the S\'ersic spheroid.
    use Star_Formation_Timescales_Spheroids
    implicit none
    type(treeNode),   pointer,  intent(inout) :: thisNode
    integer,          optional, intent(in)    :: instance
    integer                                   :: thisIndex
    double precision                          :: starFormationTimescale,gasMass

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       
       ! Get the star formation timescale.
       starFormationTimescale=Star_Formation_Timescale_Spheroid(thisNode)
       
       ! Get the gas mass.
       gasMass=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)

       ! If timescale is finite and gas mass is positive, then compute star formation rate.
       if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0) then
          Sersic_Spheroid_SFR=gasMass/starFormationTimescale
       else
          Sersic_Spheroid_SFR=0.0d0
       end if
    else
       Sersic_Spheroid_SFR=0.0d0
    end if
    return
  end function Sersic_Spheroid_SFR

  double precision function Sersic_Spheroid_Radius(thisNode,instance)
    !% Return the scale radius of the Sérsic spheroid.
    implicit none
    type(treeNode), pointer,  intent(inout) :: thisNode
    integer,        optional, intent(in)    :: instance
    integer                                 :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       Sersic_Spheroid_Radius=thisNode%components(thisIndex)%instance(1)%data(radiusIndex)
    else
       Sersic_Spheroid_Radius=0.0d0
    end if
    return
  end function Sersic_Spheroid_Radius

  double precision function Sersic_Spheroid_Half_Mass_Radius(thisNode,instance)
    !% Return the half-mass radius of the Sérsic spheroid. This is just the same as the scale radius,
    implicit none
    type(treeNode), pointer,  intent(inout) :: thisNode
    integer,        optional, intent(in)    :: instance

    Sersic_Spheroid_Half_Mass_Radius=Sersic_Spheroid_Radius(thisNode)
    return
  end function Sersic_Spheroid_Half_Mass_Radius

  subroutine Sersic_Spheroid_Radius_Set(thisNode,radius)
    !% Set the scale radius of the Sérsic spheroid.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: radius
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       thisNode%components(thisIndex)%instance(1)%data(radiusIndex)=max(radius,0.0d0)
    end if
    return
  end subroutine Sersic_Spheroid_Radius_Set

  double precision function Sersic_Spheroid_Velocity(thisNode,instance)
    !% Return the circular velocity of the Sérsic spheroid.
    implicit none
    type(treeNode), pointer,  intent(inout) :: thisNode
    integer,        optional, intent(in)    :: instance
    integer                                 :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       Sersic_Spheroid_Velocity=thisNode%components(thisIndex)%instance(1)%data(velocityIndex)
    else
       Sersic_Spheroid_Velocity=0.0d0
    end if
    return
  end function Sersic_Spheroid_Velocity

  subroutine Sersic_Spheroid_Velocity_Set(thisNode,velocity)
    !% Set the circular velocity of the Sérsic spheroid.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: velocity
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       thisNode%components(thisIndex)%instance(1)%data(velocityIndex)=velocity
    end if
    return
  end subroutine Sersic_Spheroid_Velocity_Set

  integer function Tree_Node_Sersic_Spheroid_Index(thisNode)
    !% Ensure the Sérsic spheroid component exists and return its position in the components array.
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (.not.thisNode%componentExists(componentIndex)) then
       ! Create the component.
       call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
       ! Get the index for this component.
       thisIndex=thisNode%componentIndex(componentIndex)
       ! Create the stellar properties history.
       call Stellar_Population_Properties_History_Create(thisNode,thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex))
       ! Create the star formation history.
       call Star_Formation_History_Create(thisNode,thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex))
    else
       ! Get the index for this component.
       thisIndex=thisNode%componentIndex(componentIndex)
    end if
    Tree_Node_Sersic_Spheroid_Index=thisIndex
    return
  end function Tree_Node_Sersic_Spheroid_Index

  subroutine Sersic_Spheroid_Create(thisNode)
    !% Creates an Sérsic spheroid component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message
    integer                                     :: thisIndex

    ! Display a message.
    message='Creating Sérsic spheroid component for node '
    message=message//thisNode%index()
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Get the index of the component (which will also ensure that the component is created).
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    return
  end subroutine Sersic_Spheroid_Create

  subroutine Sersic_Spheroid_Star_Formation_History_Extend(thisNode)
    !% Extend the range of a star formation history in a Sersic spheroid component for {\tt thisNode}.
    use Histories
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    integer                                     :: thisIndex

    ! Get the index of the component.
    thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
    ! Extend the range as necessary.
    call thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%extend(times=starFormationHistoryTemplate)
    return
  end subroutine Sersic_Spheroid_Star_Formation_History_Extend

  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Spheroid_Sersic_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spheroid_Sersic</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Spheroid_Sersic_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of Sérsic spheroid properties to be written to the \glc\ output file.
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
       doublePropertyComments(doubleProperty)='Mass of gas in the Sérsic spheroid.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidStellarMass'
       doublePropertyComments(doubleProperty)='Mass of stars in the Sérsic spheroid at scale length.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidAngularMomentum'
       doublePropertyComments(doubleProperty)='Angular momentum of the Sérsic spheroid.'
       doublePropertyUnitsSI (doubleProperty)=massSolar*megaParsec*kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidScaleLength'
       doublePropertyComments(doubleProperty)='Radial scale length in the Sérsic spheroid.'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='spheroidCircularVelocity'
       doublePropertyComments(doubleProperty)='Circular velocity of the Sérsic spheroid at scale length.'
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
       if (spheroidOutputStarFormationRate) then
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='spheroidStarFormationRate'
          doublePropertyComments(doubleProperty)='Spheroid star formation rate.'
          doublePropertyUnitsSI (doubleProperty)=massSolar/gigaYear
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Spheroid_Sersic_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Spheroid_Sersic_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spheroid_Sersic</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Spheroid_Sersic_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of Sérsic spheroid properties to be written to the the \glc\ output file.
    use Stellar_Population_Properties_Luminosities
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) then
       doublePropertyCount=doublePropertyCount+propertyCount+dataCount-luminositiesCount&
            &+Stellar_Population_Luminosities_Output_Count(time)
       if (spheroidOutputStarFormationRate) doublePropertyCount=doublePropertyCount+1
    end if
    return
  end subroutine Galacticus_Output_Tree_Spheroid_Sersic_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Spheroid_Sersic</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spheroid_Sersic</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Spheroid_Sersic(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store Sérsic spheroid properties in the \glc\ output file buffers.
    use Stellar_Population_Properties_Luminosities
    use Tree_Nodes
    use Kind_Numbers
    implicit none
    double precision,        intent(in)                   :: time
    type(treeNode),          intent(inout), pointer       :: thisNode
    integer,                 intent(inout)                :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)                :: integerBuffer(:,:)
    double precision,        intent(inout)                :: doubleBuffer(:,:)
    double precision,        dimension(abundancesCount)   :: gasAbundanceMasses,stellarAbundanceMasses
    double precision,        dimension(luminositiesCount) :: stellarLuminosities
    integer                                               :: iAbundance,iLuminosity

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spheroid_Angular_Momentum_Sersic(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Sersic_Spheroid_Radius(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Sersic_Spheroid_Velocity(thisNode)
       call Tree_Node_Spheroid_Gas_Abundances_Sersic    (thisNode,gasAbundanceMasses)
       call Tree_Node_Spheroid_Stellar_Abundances_Sersic(thisNode,stellarAbundanceMasses)
       do iAbundance=1,abundancesCount
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=gasAbundanceMasses(iAbundance)
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=stellarAbundanceMasses(iAbundance)
       end do
       call Tree_Node_Spheroid_Stellar_Luminosities_Sersic(thisNode,stellarLuminosities)
       do iLuminosity=1,luminositiesCount
          if (Stellar_Population_Luminosities_Output(iLuminosity,time)) then
             doubleProperty=doubleProperty+1
             doubleBuffer(doubleBufferCount,doubleProperty)=stellarLuminosities(iLuminosity)
          end if
       end do
       if (spheroidOutputStarFormationRate) then
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spheroid_SFR(thisNode)
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Spheroid_Sersic

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Sersic_Spheroid_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Sersic_Spheroid_Dump(thisNode)
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
          write (0,'(2x,a50,1x,e12.6)') 'spheroid gas mass:',Tree_Node_Spheroid_Gas_Mass_Sersic(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'spheroid stellar mass:',Tree_Node_Spheroid_Stellar_Mass_Sersic(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'spheroid angular momentum:',Tree_Node_Spheroid_Angular_Momentum_Sersic(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'spheroid scale length:',Sersic_Spheroid_Radius(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'spheroid circular velocity:',Sersic_Spheroid_Velocity(thisNode)
          call Tree_Node_Spheroid_Gas_Abundances_Sersic    (thisNode,gasAbundanceMasses    )
          call Tree_Node_Spheroid_Stellar_Abundances_Sersic(thisNode,stellarAbundanceMasses)
          do iAbundance=1,abundancesCount
             write (0,'(2x,a50,1x,e12.6)') 'spheroid gas '//char(Abundances_Names(iAbundance))//':',gasAbundanceMasses(iAbundance)
             write (0,'(2x,a50,1x,e12.6)') 'spheroid stellar '//char(Abundances_Names(iAbundance))//':',stellarAbundanceMasses(iAbundance)
          end do
          call Tree_Node_Spheroid_Stellar_Luminosities_Sersic(thisNode,stellarLuminosities)
          do iLuminosity=1,luminositiesCount
             write (0,'(2x,a50,1x,e12.6)') 'spheroid stellar '//char(Stellar_Population_Luminosities_Name(iLuminosity))//':',stellarLuminosities(iLuminosity)
          end do
       else
          write (0,'(1x,a)') 'spheroid component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Sersic_Spheroid_Dump

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Sersic_Spheroid_Star_Formation_History_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Sersic_Spheroid_Star_Formation_History_Output(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Store the star formation history in the output file.
    use Kind_Numbers
    use Tree_Nodes
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(treeNode),          intent(inout), pointer      :: thisNode
    integer,                 intent(in)                  :: iOutput
    integer(kind=kind_int8), intent(in)                  :: treeIndex
    logical,                 intent(in)                  :: nodePassesFilter
    integer                                              :: thisIndex

    ! Output the star formation history if a spheroid exists for this component.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Sersic_Spheroid_Index(thisNode)
       call Star_Formation_History_Output(thisNode,nodePassesFilter,thisNode%components(thisIndex)%instance(1)&
            &%histories(starFormationHistoryIndex),iOutput ,treeIndex,'spheroid')
    end if
    return
  end subroutine Sersic_Spheroid_Star_Formation_History_Output

  !# <decodePropertyIdentifiersTask>
  !#  <unitName>Sersic_Spheroid_Property_Identifiers_Decode</unitName>
  !# </decodePropertyIdentifiersTask>
  subroutine Sersic_Spheroid_Property_Identifiers_Decode(propertyComponent,propertyObject,propertyIndex,matchedProperty,propertyName)
    !% Decodes property identifiers to property names for the Sérsic spheroid module.
    use ISO_Varying_String
    implicit none
    integer,              intent(in)    :: propertyComponent,propertyObject,propertyIndex
    logical,              intent(inout) :: matchedProperty
    type(varying_string), intent(inout) :: propertyName

    if (methodSelected.and..not.matchedProperty) then
       if (propertyComponent == componentIndex) then
          matchedProperty=.true.
          propertyName="sersicSpheroid:"
          select case (propertyObject)
          case (objectTypeProperty)
             if      (propertyIndex == angularMomentumIndex                                                       ) then
                propertyName=propertyName//":angularMomentum"
             else if (propertyIndex == gasMassIndex                                                               ) then
                propertyName=propertyName//":gasMass"
             else if (propertyIndex == stellarMassIndex                                                           ) then
                propertyName=propertyName//":stellarMass"
             else if (propertyIndex >= gasAbundancesIndex       .and. propertyIndex <= gasAbundancesIndexEnd      ) then
                propertyName=propertyName//":gasAbundances"
             else if (propertyIndex >= stellarAbundancesIndex   .and. propertyIndex <= stellarAbundancesIndexEnd  ) then
                propertyName=propertyName//":stellarAbundances"
             else if (propertyIndex >= stellarLuminositiesIndex .and. propertyIndex <= stellarLuminositiesIndexEnd) then
                propertyName=propertyName//":stellarLuminosities"
             end if
          case (objectTypeHistory)
             select case (propertyIndex)
             case (stellarHistoryIndex)
                propertyName=propertyName//":stellarHistory"
             case (starFormationHistoryIndex)
                propertyName=propertyName//":starFormationHistory"
             end select
          end select
       end if
    end if

    return
  end subroutine Sersic_Spheroid_Property_Identifiers_Decode

end module Tree_Node_Methods_Sersic_Spheroid
