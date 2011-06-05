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


!% Contains a module of Hernquist-profile spheroid tree node methods.

module Tree_Node_Methods_Hernquist_Spheroid
  !% Implement Hernquist spheroid tree node methods.
  use Tree_Nodes
  use Histories
  use Components
  use Stellar_Population_Properties
  private
  public :: Tree_Node_Methods_Hernquist_Spheroid_Initialize, Tree_Node_Methods_Hernquist_Spheroid_Thread_Initialize,&
       & Hernquist_Spheroid_Satellite_Merging, Galacticus_Output_Tree_Spheroid_Hernquist,&
       & Galacticus_Output_Tree_Spheroid_Hernquist_Property_Count, Galacticus_Output_Tree_Spheroid_Hernquist_Names,&
       & Hernquist_Spheroid_Radius_Solver, Hernquist_Spheroid_Enclosed_Mass, Hernquist_Spheroid_Density,&
       & Hernquist_Spheroid_Rotation_Curve, Tree_Node_Spheroid_Post_Evolve_Hernquist, Tree_Node_Methods_Hernquist_Spheroid_Dump,&
       & Hernquist_Spheroid_Radius_Solver_Plausibility, Hernquist_Spheroid_Scale_Set, Hernquist_Spheroid_Post_Evolve,&
       & Hernquist_Spheroid_Star_Formation_History_Output
  
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
  double precision :: spheroidEnergeticOutflowMassRate,spheroidOutflowTimescaleMinimum
  double precision :: spheroidMassToleranceAbsolute,spheroidAngularMomentumAtScaleRadius

  ! Storage for the star formation history time range, used whene extending this range.
  double precision, allocatable, dimension(:) :: starFormationHistoryTemplate
  !$omp threadprivate(starFormationHistoryTemplate)

  ! Options controlling output.
  logical          :: spheroidOutputStarFormationRate

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
    use Memory_Management
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
       Tree_Node_Spheroid_Stellar_Properties_History_Rate_Adjust  => Tree_Node_Spheroid_Stellar_Prprts_History_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Stellar_Properties_History_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Spheroid_Star_Formation_History                  => Tree_Node_Spheroid_Star_Formation_History_Hernquist
       Tree_Node_Spheroid_Star_Formation_History_Set              => Tree_Node_Spheroid_Star_Formation_History_Set_Hernquist
       Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust      => Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Star_Formation_History_Rate_Compute     => Tree_Node_Rate_Rate_Compute_Dummy

       ! Associate pipes with procedures.
       Tree_Node_Spheroid_Gas_Sink                                => Tree_Node_Spheroid_Gas_Sink_Rate_Adjust_Hernquist
       Tree_Node_Spheroid_Gas_Energy_Input                        => Tree_Node_Spheroid_Gas_Energy_Input_Rate_Adjust_Hernquist

       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>spheroidAngularMomentumAtScaleRadius</name>
       !@   <defaultValue>0.2546479089 [value for a self-gravitating Hernquist spheroid]</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The assumed ratio of the specific angular momentum at the scale radius to the mean specific angular momentum of a Hernquist spheroid .
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidAngularMomentumAtScaleRadius',spheroidAngularMomentumAtScaleRadius,defaultValue=0.2546479089d0)
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
  end subroutine Tree_Node_Methods_Hernquist_Spheroid_Initialize
  
  !# <treeNodeCreateThreadInitialize>
  !#  <unitName>Tree_Node_Methods_Hernquist_Spheroid_Thread_Initialize</unitName>
  !# </treeNodeCreateThreadInitialize>
  subroutine Tree_Node_Methods_Hernquist_Spheroid_Thread_Initialize
    !% Initializes each thread for the tree node Hernquist spheroid methods module.
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
  end subroutine Tree_Node_Methods_Hernquist_Spheroid_Thread_Initialize
  
  !# <postEvolveTask>
  !# <unitName>Tree_Node_Spheroid_Post_Evolve_Hernquist</unitName>
  !# </postEvolveTask>
  subroutine Tree_Node_Spheroid_Post_Evolve_Hernquist(thisNode)
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
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)

       ! Trim the stellar populations properties future history.
       call thisNode%components(thisIndex)%histories(stellarHistoryIndex)%trim(Tree_Node_Time(thisNode))

       ! Trap negative gas masses.
       if (Tree_Node_Spheroid_Gas_Mass(thisNode) < 0.0d0) then
          
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=dabs(thisNode%components(thisIndex)%properties(gasMassIndex,propertyValue))&
               &/(thisNode%components(thisIndex)%properties(stellarMassIndex,propertyValue)&
               &+dabs(thisNode%components(thisIndex)%properties(gasMassIndex,propertyValue)))
          !$omp critical (Hernquist_Spheroid_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.          
             message='Warning: spheroid has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index            = '//thisNode%index() //char(10)
             write (valueString,'(e12.6)') thisNode%components(thisIndex)%properties(gasMassIndex    ,propertyValue)
             message=message//'  Spheroid gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') thisNode%components(thisIndex)%properties(stellarMassIndex,propertyValue)
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
          !$omp end critical (Hernquist_Spheroid_Post_Evolve_Check)
          
          ! Get the specific angular momentum of the spheroid material
          spheroidMass= thisNode%components(thisIndex)%properties(gasMassIndex    ,propertyValue) &
               &       +thisNode%components(thisIndex)%properties(stellarMassIndex,propertyValue)
          if (spheroidMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             thisNode%components(thisIndex)%properties(stellarMassIndex                                    ,propertyValue)=0.0d0
             thisNode%components(thisIndex)%properties(stellarAbundancesIndex  :stellarAbundancesIndexEnd  ,propertyValue)=0.0d0  
             thisNode%components(thisIndex)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyValue)=0.0d0
          else
             specificAngularMomentum=thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyValue)&
                  &/spheroidMass
          end if
          
          ! Reset the gas, abundances and angular momentum of the spheroid.
          thisNode%components(thisIndex)%properties(gasMassIndex                            ,propertyValue)=0.0d0
          thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)=0.0d0
          thisNode%components(thisIndex)%properties(angularMomentumIndex                    ,propertyValue)= &
               & specificAngularMomentum*thisNode%components(thisIndex)%properties(stellarMassIndex,propertyValue)
       end if
       
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
    if (gasMass > 0.0d0 .and. gasMass+stellarMass > 0.0d0) then
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
    if (gasMass > 0.0d0 .and. gasMass+stellarMass > 0.0d0) then
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
    procedure(), pointer, intent(inout) :: interruptProcedure

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
    procedure(), pointer, intent(inout) :: interruptProcedure
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
    use Galactic_Structure_Options
    use Galacticus_Output_Star_Formation_Histories
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),            pointer, intent(inout) :: thisNode
    logical,                            intent(inout) :: interrupt
    procedure(),               pointer, intent(inout) :: interruptProcedure
    integer                                           :: thisIndex
    double precision                                  :: starFormationRate,stellarMassRate ,fuelMassRate,fuelMass&
         &,massOutflowRate,spheroidMass,angularMomentumOutflowRate,energyInputRate,gasMass,spheroidDynamicalTime
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
       starFormationRate=Hernquist_Spheroid_SFR(thisNode)

       ! Get the available fuel mass.
       fuelMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
       
       ! Find the metallicity of the fuel supply.
       call Tree_Node_Spheroid_Gas_Abundances_Hernquist(thisNode,abundancesValue)
       call fuelAbundances%pack(abundancesValue)
       call fuelAbundances%massToMassFraction(fuelMass)

       ! Get the component index.
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       
       ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
       call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,componentTypeSpheroid,thisNode&
            &,thisNode%components(thisIndex)%histories(stellarHistoryIndex),stellarMassRate,stellarAbundancesRates &
            &,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
       
       ! Adjust rates.
       call Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust_Hernquist        (thisNode,interrupt,interruptProcedure,stellarMassRate         )
       call Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Hernquist            (thisNode,interrupt,interruptProcedure,fuelMassRate            )
       call stellarAbundancesRates%unpack(abundancesValue)
       call Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust_Hernquist  (thisNode,interrupt,interruptProcedure,abundancesValue         )
       call fuelAbundancesRates%unpack(abundancesValue)
       call Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Hernquist      (thisNode,interrupt,interruptProcedure,abundancesValue         )
       call Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,stellarLuminositiesRates)

       ! Record the star formation history.
       call Star_Formation_History_Record(thisNode,thisNode%components(thisIndex)%histories(starFormationHistoryIndex)&
            &,fuelAbundances,starFormationRate)
       
       ! Find rate of outflow of material from the spheroid and pipe it to the outflowed reservoir.
       massOutflowRate=Star_Formation_Feedback_Spheroid_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
       if (massOutflowRate > 0.0d0) then
          ! Get the masses of the spheroid.
          gasMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
          spheroidMass=gasMass+Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
          
          ! Limit the outflow rate timescale to a multiple of the dynamical time.
          spheroidDynamicalTime=Mpc_per_km_per_s_To_Gyr*Tree_Node_Spheroid_Radius(thisNode)/Tree_Node_Spheroid_Velocity(thisNode)
          ! Limit the mass outflow rate.
          massOutflowRate=min(massOutflowRate,gasMass/spheroidOutflowTimescaleMinimum/spheroidDynamicalTime)
          call Tree_Node_Hot_Halo_Outflow_Mass_To                       (thisNode,interrupt,interruptProcedure,&
               & massOutflowRate           )
          call Tree_Node_Spheroid_Gas_Mass_Rate_Adjust_Hernquist        (thisNode,interrupt,interruptProcedure,&
               &-massOutflowRate           )
          
          ! Compute the angular momentum outflow rate.
          if (spheroidMass > 0.0d0) then
             angularMomentumOutflowRate=(massOutflowRate/spheroidMass)*Tree_Node_Spheroid_Angular_Momentum_Hernquist(thisNode)
             call Tree_Node_Hot_Halo_Outflow_Angular_Momentum_To           (thisNode,interrupt,interruptProcedure,&
                  & angularMomentumOutflowRate)
             call Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,&
                  &-angularMomentumOutflowRate)
          end if
          
          ! Compute the abundances outflow rate.
          call Tree_Node_Spheroid_Gas_Abundances_Hernquist(thisNode,abundancesValue)
          call Abundances_Mass_To_Mass_Fraction(abundancesValue,gasMass)
          abundancesOutflowRate=massOutflowRate*abundancesValue
          call Tree_Node_Hot_Halo_Outflow_Abundances_To                 (thisNode,interrupt,interruptProcedure,&
               & abundancesOutflowRate)
          call Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust_Hernquist  (thisNode,interrupt,interruptProcedure,&
               &-abundancesOutflowRate)
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

  function Tree_Node_Spheroid_Star_Formation_History_Hernquist(thisNode)
    !% Return the spheroid star formation history.
    use Histories
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(history)                          :: Tree_Node_Spheroid_Star_Formation_History_Hernquist
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       Tree_Node_Spheroid_Star_Formation_History_Hernquist=thisNode%components(thisIndex)%histories(starFormationHistoryIndex)
    else
       Tree_Node_Spheroid_Star_Formation_History_Hernquist=nullHistory
    end if
    return
  end function Tree_Node_Spheroid_Star_Formation_History_Hernquist

  subroutine Tree_Node_Spheroid_Star_Formation_History_Set_Hernquist(thisNode,thisHistory)
    !% Set the node spheroid star formation history.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(history),           intent(in)    :: thisHistory
    integer                                :: thisIndex

    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    call thisNode%components(thisIndex)%histories(starFormationHistoryIndex)%destroy()
    thisNode%components(thisIndex)%histories(starFormationHistoryIndex)=thisHistory
    return
  end subroutine Tree_Node_Spheroid_Star_Formation_History_Set_Hernquist

  subroutine Tree_Node_Spheroid_Stellar_Prprts_History_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the rates for the stellar properties history.
    use Histories
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    logical,                  intent(inout) :: interrupt
    procedure(),     pointer, intent(inout) :: interruptProcedure
    type(history),            intent(in)    :: rateAdjustments
    integer                                 :: thisIndex

    ! If no Hernquist spheroid component currently exists and we have some non-zero rate into it then interrupt and create an Hernquist
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments%data /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Hernquist_Spheroid_Create
       end if
       return
    end if
    ! Get the index for this component.
    thisIndex=thisNode%componentIndex(componentIndex)
    ! Adjust the rate.
    call thisNode%components(thisIndex)%histories(stellarHistoryIndex)%add(rateAdjustments,addTo=historyRates)
    return
  end subroutine Tree_Node_Spheroid_Stellar_Prprts_History_Rate_Adjust_Hernquist

  subroutine Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Hernquist(thisNode,interrupt,interruptProcedure,rateAdjustments)
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

    ! If no Hernquist spheroid component currently exists and we have some non-zero rate into it then interrupt and create an Hernquist
    ! spheroid.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments%data /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Hernquist_Spheroid_Create
       end if
       return
    end if
    ! Get the index for this component.
    thisIndex=thisNode%componentIndex(componentIndex)
    ! Ensure that a history already exists in the spheroid.
    if (.not.thisNode%components(thisIndex)%histories(starFormationHistoryIndex)%exists())                &
         & call Galacticus_Error_Report('Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Hernquist' &
         & ,'no star formation history has been created in spheroid')
    ! Check if the star formation history in the spheroid spans a sufficient range to accept the input rates.
    if (rateAdjustments%time(1) < thisNode%components(thisIndex)%histories(starFormationHistoryIndex)%time(1) .or.                      &
         & rateAdjustments%time(size(rateAdjustments%time)) > thisNode%components(thisIndex)%histories(starFormationHistoryIndex)%time( &
         & size(thisNode%components(thisIndex)%histories(starFormationHistoryIndex)%time))) then
       ! It does not, so interrupt evolution and extend the history.
       if (allocated(starFormationHistoryTemplate)) call Dealloc_Array(starFormationHistoryTemplate)
       call Alloc_Array(starFormationHistoryTemplate,shape(rateAdjustments%time))
       starFormationHistoryTemplate=rateAdjustments%time
       interrupt=.true.
       interruptProcedure => Hernquist_Spheroid_Star_Formation_History_Extend
       return
    end if

    ! Adjust the rate.
    call thisNode%components(thisIndex)%histories(starFormationHistoryIndex)%combine(rateAdjustments,addTo=historyRates)
    return
  end subroutine Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Hernquist

  !# <scaleSetTask>
  !#  <unitName>Hernquist_Spheroid_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Hernquist_Spheroid_Scale_Set(thisNode)
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
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)

       ! Set scale for angular momentum.
       angularMomentum=Tree_Node_Spheroid_Angular_Momentum(thisNode)+Tree_Node_Disk_Angular_Momentum(thisNode)
       thisNode%components(thisIndex)%properties(angularMomentumIndex,propertyScale)=max(angularMomentum,angularMomentumMinimum)

       ! Set scale for gas mass.
       mass=Tree_Node_Spheroid_Gas_Mass(thisNode)+Tree_Node_Disk_Gas_Mass(thisNode)
       thisNode%components(thisIndex)%properties(gasMassIndex,propertyScale)=gasMassScaling*max(mass,massMinimum)

       ! Set scale for stellar mass.
       mass=Tree_Node_Spheroid_Stellar_Mass(thisNode)+Tree_Node_Disk_Stellar_Mass(thisNode)
       thisNode%components(thisIndex)%properties(stellarMassIndex,propertyScale)=max(mass,massMinimum)

       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for gas abundances.
          call Tree_Node_Spheroid_Gas_Abundances_Hernquist(thisNode,abundancesSpheroid)
          call Tree_Node_Disk_Gas_Abundances              (thisNode,abundancesDisk    )
          thisNode%components(thisIndex)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyScale)=gasMassScaling*max(abundancesDisk+abundancesSpheroid,massMinimum)
          
          ! Set scale for stellar abundances.
          call Tree_Node_Spheroid_Stellar_Abundances_Hernquist(thisNode,abundancesSpheroid)
          call Tree_Node_Disk_Stellar_Abundances              (thisNode,abundancesDisk    )
          thisNode%components(thisIndex)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyScale)=max(abundancesDisk+abundancesSpheroid,massMinimum)
       end if

       ! Set scales for stellar luminosities if necessary.
       if (luminositiesCount > 0) then        
          ! Set scale for stellar luminosities.
          call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist(thisNode,luminositiesSpheroid)
          call Tree_Node_Disk_Stellar_Luminosities              (thisNode,luminositiesDisk    )
          thisNode%components(thisIndex)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyScale)=max(luminositiesDisk+luminositiesSpheroid,luminosityMinimum)
       end if

       ! Set scales for stellar population properties history.
       call Tree_Node_Spheroid_Stellar_Abundances_Hernquist(thisNode,abundancesSpheroid)
       call stellarAbundances%pack(abundancesSpheroid)
       call Stellar_Population_Properties_Scales(thisNode%components(thisIndex)%histories(stellarHistoryIndex      ),Tree_Node_Spheroid_Stellar_Mass(thisNode),stellarAbundances)
       call Star_Formation_History_Scales       (thisNode%components(thisIndex)%histories(starFormationHistoryIndex),Tree_Node_Spheroid_Stellar_Mass(thisNode),stellarAbundances)
    end if
    return
  end subroutine Hernquist_Spheroid_Scale_Set

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
    type(history)                                  :: historyDisk,historySpheroid,thisHistory
    double precision                               :: spheroidSpecificAngularMomentum,diskSpecificAngularMomentum,angularMomentum&
         &,spheroidMass

    ! Check that method is selected.
    if (methodSelected) then
       
       ! Find the node to merge with.
       call thisNode%mergesWith(hostNode)

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
          ! Also add stellar properties histories.
          historyDisk    =Tree_Node_Disk_Star_formation_History                  (hostNode                )
          historySpheroid=Tree_Node_Spheroid_Star_formation_History_Hernquist    (hostNode                )
          call historyDisk%combine(historySpheroid)
          call Tree_Node_Disk_Star_formation_History_Set                         (hostNode,historyDisk    )
          call historySpheroid%reset()
          call Tree_Node_Spheroid_Star_formation_History_Set_Hernquist           (hostNode,historySpheroid)
          call historyDisk%destroy()
          call historySpheroid%destroy()
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
          ! Also add stellar properties histories.
          historyDisk    =Tree_Node_Disk_Star_formation_History                  (hostNode                )
          historySpheroid=Tree_Node_Spheroid_Star_formation_History_Hernquist    (hostNode                )
          call historySpheroid%combine(historyDisk)
          call Tree_Node_Spheroid_Star_formation_History_Set_Hernquist           (hostNode,historySpheroid)
          call historyDisk%reset()
          call Tree_Node_Disk_Star_formation_History_Set                         (hostNode,historyDisk    )
          call historyDisk%destroy()
          call historySpheroid%destroy()
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
             ! Also add star formation histories.
             historySpheroid=Tree_Node_Spheroid_Star_formation_History_Hernquist    (thisNode                )
             thisHistory    =Tree_Node_Disk_Star_formation_History                  (hostNode                )
             call thisHistory%combine(historySpheroid)
             call Tree_Node_Disk_Star_formation_History_Set                         (hostNode,thisHistory    )
             call historySpheroid%reset()
             call Tree_Node_Spheroid_Star_formation_History_Set_Hernquist           (thisNode,historySpheroid)
             call thisHistory%destroy()
             call historySpheroid%destroy()
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
             ! Also add star formation histories.
             historySpheroid=Tree_Node_Spheroid_Star_formation_History_Hernquist    (thisNode                )
             thisHistory    =Tree_Node_Spheroid_Star_formation_History_Hernquist    (hostNode                )
             call thisHistory%combine(historySpheroid)
             call Tree_Node_Spheroid_Star_formation_History_Set_Hernquist           (hostNode,thisHistory    )
             call historySpheroid%reset()
             call Tree_Node_Spheroid_Star_formation_History_Set_Hernquist           (thisNode,historySpheroid)
             call thisHistory%destroy()
             call historySpheroid%destroy()
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
          ! Note that the remnant specific angular momentum computed by the merger remnant modules automatically gives the mean
          ! specific angular momentum of the component by virtue of the fact that it computes the ratio of the actual angular
          ! momentum to the contribution from the component's own rotation curve at its scale radius.
          angularMomentum=remnantSpecificAngularMomentum*(Tree_Node_Spheroid_Gas_Mass(hostNode) &
               &+Tree_Node_Spheroid_Stellar_Mass(hostNode))
          call Tree_Node_Spheroid_Angular_Momentum_Set_Hernquist(hostNode,angularMomentum)
       end if

    end if
    return
  end subroutine Hernquist_Spheroid_Satellite_Merging
  
  !# <enclosedMassTask>
  !#  <unitName>Hernquist_Spheroid_Enclosed_Mass</unitName>
  !# </enclosedMassTask>
  subroutine Hernquist_Spheroid_Enclosed_Mass(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for an Hernquist spheroid.
    use Galactic_Structure_Options
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType,weightBy,weightIndex
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass
    double precision                         :: fractionalRadius,spheroidRadius
    
    componentMass=0.0d0
    if (.not.methodSelected                           ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
    if (.not.thisNode%componentExists(componentIndex)) return
    select case (weightBy)
    case (weightByMass      )
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic)
          componentMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)+Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
       case (massTypeGaseous)
          componentMass=Tree_Node_Spheroid_Gas_Mass_Hernquist(thisNode)
       case (massTypeStellar)
          componentMass=Tree_Node_Spheroid_Stellar_Mass_Hernquist(thisNode)
       end select
    case (weightByLuminosity)
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
          call Tree_Node_Spheroid_Stellar_Luminosities_Hernquist(thisNode,luminositiesSpheroid)
          componentMass=luminositiesSpheroid(weightIndex)
       end select
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
          call Hernquist_Spheroid_Enclosed_Mass(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
          if (componentMass > 0.0d0) componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass)/dsqrt(radius)
       end if
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
    if (componentDensity <= 0.0d0) then
       componentDensity=0.0d0
       return
    end if
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
    type(treeNode),              pointer, intent(inout) :: thisNode
    logical,                              intent(out)   :: componentActive
    double precision,                     intent(out)   :: specificAngularMomentum
    procedure(),                 pointer, intent(out)   :: Radius_Set,Velocity_Set
    procedure(double precision), pointer, intent(out)   :: Radius_Get,Velocity_Get
    double precision                                    :: specificAngularMomentumMean,angularMomentum,spheroidMass

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
          specificAngularMomentum=spheroidAngularMomentumAtScaleRadius*specificAngularMomentumMean
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
       call Stellar_Population_Properties_History_Create(thisNode,thisNode%components(thisIndex)%histories(stellarHistoryIndex))
       ! Create the star formation history.
       call Star_Formation_History_Create(thisNode,thisNode%components(thisIndex)%histories(starFormationHistoryIndex))
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

  subroutine Hernquist_Spheroid_Star_Formation_History_Extend(thisNode)
    !% Extend the range of a star formation history in a Hernquist spheroid component for {\tt thisNode}.
    use Histories
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    integer                                     :: thisIndex

    ! Get the index of the component.
    thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
    ! Extend the range as necessary.
    call thisNode%components(thisIndex)%histories(starFormationHistoryIndex)%extend(times=starFormationHistoryTemplate)
    return
  end subroutine Hernquist_Spheroid_Star_Formation_History_Extend

  
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
       if (spheroidOutputStarFormationRate) then
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='spheroidStarFormationRate'
          doublePropertyComments(doubleProperty)='Spheroid star formation rate.'
          doublePropertyUnitsSI (doubleProperty)=massSolar/gigaYear
       end if
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

    if (methodSelected) then
       doublePropertyCount=doublePropertyCount+propertyCount+dataCount-luminositiesCount&
            &+Stellar_Population_Luminosities_Output_Count(time)
       if (spheroidOutputStarFormationRate) doublePropertyCount=doublePropertyCount+1
    end if
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
       if (spheroidOutputStarFormationRate) then
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spheroid_SFR(thisNode)
       end if
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

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Hernquist_Spheroid_Star_Formation_History_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Hernquist_Spheroid_Star_Formation_History_Output(thisNode,iOutput,treeIndex,nodePassesFilter)
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
       thisIndex=Tree_Node_Hernquist_Spheroid_Index(thisNode)
       call Star_Formation_History_Output(thisNode,nodePassesFilter,thisNode%components(thisIndex)&
            &%histories(starFormationHistoryIndex),iOutput ,treeIndex,'spheroid')
    end if
    return
  end subroutine Hernquist_Spheroid_Star_Formation_History_Output

end module Tree_Node_Methods_Hernquist_Spheroid
