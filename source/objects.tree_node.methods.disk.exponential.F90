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


!% Contains a module of exponential disk tree node methods.

module Tree_Node_Methods_Exponential_Disk
  !% Implement exponential disk tree node methods.
  use Tree_Node_Methods_Disk_Exponential_Data
  use Tree_Nodes
  use Histories
  use Components
  use Stellar_Population_Properties
  private
  public :: Tree_Node_Methods_Exponential_Disk_Initialize, Tree_Node_Methods_Exponential_Disk_Thread_Initialize,&
       & Exponential_Disk_Satellite_Merging, Galacticus_Output_Tree_Disk_Exponential,&
       & Galacticus_Output_Tree_Disk_Exponential_Property_Count, Galacticus_Output_Tree_Disk_Exponential_Names,&
       & Exponential_Disk_Radius_Solver, Exponential_Disk_Enclosed_Mass, Exponential_Disk_Rotation_Curve,&
       & Tree_Node_Disk_Post_Evolve_Exponential, Tree_Node_Methods_Exponential_Disk_Dump,&
       & Exponential_Disk_Radius_Solver_Plausibility, Tree_Node_Methods_Exponential_Disk_State_Store,&
       & Tree_Node_Methods_Exponential_Disk_State_Retrieve, Exponential_Disk_Scale_Set,&
       & Exponential_Disk_Star_Formation_History_Output, Exponential_Disk_Property_Identifiers_Decode
  
  ! Internal count of abundances and work arrays.
  integer                                     :: abundancesCount
  double precision, allocatable, dimension(:) :: abundancesTransferRate,abundancesOutflowRate,abundancesValue,abundancesDisk,abundancesSpheroid
  !$omp threadprivate(abundancesTransferRate,abundancesOutflowRate,abundancesValue,abundancesDisk,abundancesSpheroid)

  ! Internal count of luminosities and work arrays.
  integer                                     :: luminositiesCount
  double precision, allocatable, dimension(:) :: stellarLuminositiesRates,luminositiesTransferRate,luminositiesValue&
       &,luminositiesDisk,luminositiesSpheroid
  !$omp threadprivate(stellarLuminositiesRates,luminositiesTransferRate,luminositiesValue,luminositiesDisk,luminositiesSpheroid)

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
  !#  <methodName>Tree_Node_Disk_Radius</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Disk_Half_Mass_Radius</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Disk_Velocity</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Disk_Gas_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Disk_Angular_Momentum</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Disk_Stellar_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Disk_Gas_Abundances</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Disk_Stellar_Abundances</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Disk_Stellar_Luminosities</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Disk_SFR</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="history">
  !#  <methodName>Tree_Node_Disk_Stellar_Properties_History</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="history">
  !#  <methodName>Tree_Node_Disk_Star_Formation_History</methodName>
  !# </treeNodeMethodsPointer>

  ! Parameters controlling the physical implementation.
  double precision :: diskMassToleranceAbsolute,diskOutflowTimescaleMinimum

  ! Tabulation of the exponential disk rotation curve.
  integer,          parameter                 :: rotationCurvePointsPerDecade=10
  integer                                     :: rotationCurvePointsCount
  logical                                     :: rotationCurveInitialized=.false.
  double precision, parameter                 :: rotationCurveHalfRadiusMinimumDefault=1.0d-6,rotationCurveHalfRadiusMaximumDefault=10.0d0
  double precision                            :: rotationCurveHalfRadiusMinimum=rotationCurveHalfRadiusMinimumDefault&
       &,rotationCurveHalfRadiusMaximum=rotationCurveHalfRadiusMaximumDefault
  double precision, allocatable, dimension(:) :: rotationCurveHalfRadius,rotationCurveBesselFactors
  double precision                            :: scaleLengthFactor
  logical                                     :: scaleLengthFactorSet=.false.

  ! Options controlling output.
  logical                                     :: diskOutputStarFormationRate

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Exponential_Disk_Initialize</unitName>
  !#  <optionName default="exponential">treeNodeMethodDisk</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Exponential_Disk_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node exponential disk methods module.
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
    if (componentOption.eq.'exponential') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Exponential disk method selected [component index '
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
       Tree_Node_Disk_Radius                                  => Exponential_Disk_Radius
       Tree_Node_Disk_Radius_Set                              => null()
       Tree_Node_Disk_Radius_Rate_Adjust                      => null()
       Tree_Node_Disk_Radius_Rate_Compute                     => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Disk_Half_Mass_Radius                        => Exponential_Disk_Half_Mass_Radius
       Tree_Node_Disk_Half_Mass_Radius_Set                    => null()
       Tree_Node_Disk_Half_Mass_Radius_Rate_Adjust            => null()
       Tree_Node_Disk_Half_Mass_Radius_Rate_Compute           => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Disk_Velocity                                => Exponential_Disk_Velocity
       Tree_Node_Disk_Velocity_Set                            => null()
       Tree_Node_Disk_Velocity_Rate_Adjust                    => null()
       Tree_Node_Disk_Velocity_Rate_Compute                   => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Disk_Angular_Momentum                        => Tree_Node_Disk_Angular_Momentum_Exponential
       Tree_Node_Disk_Angular_Momentum_Set                    => Tree_Node_Disk_Angular_Momentum_Set_Exponential
       Tree_Node_Disk_Angular_Momentum_Rate_Adjust            => Tree_Node_Disk_Angular_Momentum_Rate_Adjust_Exponential
       Tree_Node_Disk_Angular_Momentum_Rate_Compute           => Tree_Node_Disk_Angular_Momentum_Rate_Compute_Exponential
      
       Tree_Node_Disk_Gas_Mass                                => Tree_Node_Disk_Gas_Mass_Exponential
       Tree_Node_Disk_Gas_Mass_Set                            => Tree_Node_Disk_Gas_Mass_Set_Exponential
       Tree_Node_Disk_Gas_Mass_Rate_Adjust                    => Tree_Node_Disk_Gas_Mass_Rate_Adjust_Exponential
       Tree_Node_Disk_Gas_Mass_Rate_Compute                   => Tree_Node_Disk_Gas_Mass_Rate_Compute_Exponential

       Tree_Node_Disk_Stellar_Mass                            => Tree_Node_Disk_Stellar_Mass_Exponential
       Tree_Node_Disk_Stellar_Mass_Set                        => Tree_Node_Disk_Stellar_Mass_Set_Exponential
       Tree_Node_Disk_Stellar_Mass_Rate_Adjust                => Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Exponential
       Tree_Node_Disk_Stellar_Mass_Rate_Compute               => Tree_Node_Disk_Stellar_Mass_Rate_Compute_Exponential

       Tree_Node_Disk_Stellar_Abundances                      => Tree_Node_Disk_Stellar_Abundances_Exponential
       Tree_Node_Disk_Stellar_Abundances_Set                  => Tree_Node_Disk_Stellar_Abundances_Set_Exponential
       Tree_Node_Disk_Stellar_Abundances_Rate_Adjust          => Tree_Node_Disk_Stellar_Abundances_Rate_Adjust_Exponential
       Tree_Node_Disk_Stellar_Abundances_Rate_Compute         => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Disk_Gas_Abundances                          => Tree_Node_Disk_Gas_Abundances_Exponential
       Tree_Node_Disk_Gas_Abundances_Set                      => Tree_Node_Disk_Gas_Abundances_Set_Exponential
       Tree_Node_Disk_Gas_Abundances_Rate_Adjust              => Tree_Node_Disk_Gas_Abundances_Rate_Adjust_Exponential
       Tree_Node_Disk_Gas_Abundances_Rate_Compute             => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Disk_Stellar_Luminosities                    => Tree_Node_Disk_Stellar_Luminosities_Exponential
       Tree_Node_Disk_Stellar_Luminosities_Set                => Tree_Node_Disk_Stellar_Luminosities_Set_Exponential
       Tree_Node_Disk_Stellar_Luminosities_Rate_Adjust        => Tree_Node_Disk_Stellar_Luminosities_Rate_Adjust_Exponential
       Tree_Node_Disk_Stellar_Luminosities_Rate_Compute       => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Disk_SFR                                     => Exponential_Disk_SFR
       Tree_Node_Disk_SFR_Set                                 => null()
       Tree_Node_Disk_SFR_Rate_Adjust                         => null()
       Tree_Node_Disk_SFR_Rate_Compute                        => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Disk_Stellar_Properties_History              => Tree_Node_Disk_Stellar_Properties_History_Exponential
       Tree_Node_Disk_Stellar_Properties_History_Set          => Tree_Node_Disk_Stellar_Properties_History_Set_Exponential
       Tree_Node_Disk_Stellar_Properties_History_Rate_Adjust  => null()
       Tree_Node_Disk_Stellar_Properties_History_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Disk_Star_Formation_History                  => Tree_Node_Disk_Star_Formation_History_Exponential
       Tree_Node_Disk_Star_Formation_History_Set              => Tree_Node_Disk_Star_Formation_History_Set_Exponential
       Tree_Node_Disk_Star_Formation_History_Rate_Adjust      => null()
       Tree_Node_Disk_Star_Formation_History_Rate_Compute     => Tree_Node_Rate_Rate_Compute_Dummy

       ! Grab the cooling mass/angular momentum pipes from the hot halo component.
       if (.not.associated(Tree_Node_Hot_Halo_Cooling_Mass_To)) then
          Tree_Node_Hot_Halo_Cooling_Mass_To => Tree_Node_Disk_Gas_Mass_Rate_Adjust_Exponential
       else
          call Galacticus_Error_Report('Tree_Node_Methods_Exponential_Disk_Initialize','expected to find unclaimed hot halo cooling pipe')
       end if
       if (.not.associated(Tree_Node_Hot_Halo_Cooling_Angular_Momentum_To)) then
          Tree_Node_Hot_Halo_Cooling_Angular_Momentum_To => Tree_Node_Disk_Angular_Momentum_Rate_Adjust_Exponential
       else
          call Galacticus_Error_Report('Tree_Node_Methods_Exponential_Disk_Initialize','expected to find unclaimed hot halo cooling angular momentum pipe')
       end if
       if (.not.associated(Tree_Node_Hot_Halo_Cooling_Abundances_To)) then
          Tree_Node_Hot_Halo_Cooling_Abundances_To => Tree_Node_Disk_Gas_Abundances_Rate_Adjust_Exponential
       else
          call Galacticus_Error_Report('Tree_Node_Methods_Exponential_Disk_Initialize','expected to find unclaimed hot halo cooling abundances pipe')
       end if

       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>diskMassToleranceAbsolute</name>
       !@   <defaultValue>$10^{-6} M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The mass tolerance used to judge whether the disk is physically plausible.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('diskMassToleranceAbsolute',diskMassToleranceAbsolute,defaultValue=1.0d-6)
       !@ <inputParameter>
       !@   <name>diskOutflowTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the disk dynamical time) on which outflows may deplete gas in the disk.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowTimescaleMinimum',diskOutflowTimescaleMinimum,defaultValue=1.0d-3)
       !@ <inputParameter>
       !@   <name>diskOutputStarFormationRate</name>
       !@   <defaultValue>false.</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Determines whether or not the star formation rate in the disk of each galaxy will be output.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutputStarFormationRate',diskOutputStarFormationRate,defaultValue=.false.)

    end if
    return
  end subroutine Tree_Node_Methods_Exponential_Disk_Initialize
  
  
  !# <treeNodeCreateThreadInitialize>
  !#  <unitName>Tree_Node_Methods_Exponential_Disk_Thread_Initialize</unitName>
  !# </treeNodeCreateThreadInitialize>
  subroutine Tree_Node_Methods_Exponential_Disk_Thread_Initialize
    !% Initializes each thread for the tree node exponential disk methods module.
    use Memory_Management
    implicit none

    ! Check if this implementation is selected.
    if (methodSelected.and..not.allocated(abundancesTransferRate)) then

       ! Allocate work arrays for abundances.
       call Alloc_Array(abundancesTransferRate,[abundancesCount])
       call Alloc_Array(abundancesOutflowRate ,[abundancesCount])
       call Alloc_Array(abundancesValue       ,[abundancesCount])
       call Alloc_Array(abundancesDisk        ,[abundancesCount])
       call Alloc_Array(abundancesSpheroid    ,[abundancesCount])

       ! Allocate work arrays for luminosities.
       call Alloc_Array(stellarLuminositiesRates,[luminositiesCount])
       call Alloc_Array(luminositiesTransferRate,[luminositiesCount])
       call Alloc_Array(luminositiesValue       ,[luminositiesCount])
       call Alloc_Array(luminositiesDisk        ,[luminositiesCount])
       call Alloc_Array(luminositiesSpheroid    ,[luminositiesCount])

    end if
    return
  end subroutine Tree_Node_Methods_Exponential_Disk_Thread_Initialize
    
  !# <postEvolveTask>
  !# <unitName>Tree_Node_Disk_Post_Evolve_Exponential</unitName>
  !# </postEvolveTask>
  subroutine Tree_Node_Disk_Post_Evolve_Exponential(thisNode)
    !% Trim histories attached to the disk.
    use Galacticus_Display
    use String_Handling
    use ISO_Varying_String
    use Histories
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, save                   :: fractionalErrorMaximum=0.0d0
    integer                                  :: thisIndex
    double precision                         :: specificAngularMomentum,fractionalError,diskMass
    character(len=20)                        :: valueString
    type(varying_string)                     :: message

    ! Check if an exponential disk component exists.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then
       ! Get the index of the component.
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)

       ! Trim the stellar populations properties future history.
       call thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)%trim(Tree_Node_Time(thisNode))

       ! Trap negative gas masses.
       if (Tree_Node_Disk_Gas_Mass(thisNode) < 0.0d0) then
          
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=dabs(thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue))&
               &/(thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)&
               &+dabs(thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)))
          !$omp critical (Exponential_Disk_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.          
             message='Warning: disk has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//thisNode%index() //char(10)
             write (valueString,'(e12.6)') thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex    ,propertyValue)
             message=message//'  Disk gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
             message=message//'  Disk stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure     = '//trim(valueString)//char(10)
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
          !$omp end critical (Exponential_Disk_Post_Evolve_Check)
          
          ! Get the specific angular momentum of the disk material
          diskMass= thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex    ,propertyValue) &
               &   +thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
          if (diskMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex                                    ,propertyValue)=0.0d0
             thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex  :stellarAbundancesIndexEnd  ,propertyValue)=0.0d0  
             thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyValue)=0.0d0
          else
             specificAngularMomentum=thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)&
                  &/diskMass
          end if

          ! Reset the gas, abundances and angular momentum of the disk.
          thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex                            ,propertyValue)=0.0d0
          thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)=0.0d0
          thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex                    ,propertyValue)= &
               & specificAngularMomentum*thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)


       end if
       
    end if
    return
  end subroutine Tree_Node_Disk_Post_Evolve_Exponential

  double precision function Tree_Node_Disk_Gas_Mass_Exponential(thisNode,instance)
    !% Return the node exponential disk gas mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       Tree_Node_Disk_Gas_Mass_Exponential=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)
    else
       Tree_Node_Disk_Gas_Mass_Exponential=0.0d0
    end if
    return
  end function Tree_Node_Disk_Gas_Mass_Exponential

  subroutine Tree_Node_Disk_Gas_Mass_Set_Exponential(thisNode,mass,instance)
    !% Set the node exponential disk gas mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Disk_Gas_Mass_Set_Exponential

  subroutine Tree_Node_Disk_Gas_Mass_Rate_Adjust_Exponential(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node exponential disk gas mass rate of change.
    use Cosmological_Parameters
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    ! If no exponential disk component currently exists and we have some cooling into it then interrupt and create an exponential
    ! disk.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Exponential_Disk_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Disk_Gas_Mass_Rate_Adjust_Exponential

  subroutine Tree_Node_Disk_Gas_Mass_Rate_Compute_Exponential(thisNode,interrupt,interruptProcedure)
    !% Compute the exponential disk node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure

    ! No need to do anything here - the hot halo will pipe a cooling rate to us.
    return
  end subroutine Tree_Node_Disk_Gas_Mass_Rate_Compute_Exponential

  double precision function Tree_Node_Disk_Stellar_Mass_Exponential(thisNode,instance)
    !% Return the node exponential disk stellar mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       Tree_Node_Disk_Stellar_Mass_Exponential=thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
    else
       Tree_Node_Disk_Stellar_Mass_Exponential=0.0d0
    end if
    return
  end function Tree_Node_Disk_Stellar_Mass_Exponential

  subroutine Tree_Node_Disk_Stellar_Mass_Set_Exponential(thisNode,mass,instance)
    !% Set the node exponential disk stellar mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Disk_Stellar_Mass_Set_Exponential

  subroutine Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Exponential(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node exponential disk stellar mass rate of change.
    use Cosmological_Parameters
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    ! If no exponential disk component currently exists and we have some cooling into it then interrupt and create an exponential
    ! disk.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Exponential_Disk_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Exponential

  subroutine Tree_Node_Disk_Stellar_Mass_Rate_Compute_Exponential(thisNode,interrupt,interruptProcedureReturn)
    !% Compute the exponential disk node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    use Star_Formation_Feedback_Disks
    use Star_Formation_Feedback_Expulsion_Disks
    use Abundances_Structure
    use Galactic_Structure_Options
    use Galactic_Dynamics_Bar_Instabilities
    use Galacticus_Output_Star_Formation_Histories 
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),            pointer, intent(inout) :: thisNode
    logical,                            intent(inout) :: interrupt
    procedure(),               pointer, intent(inout) :: interruptProcedureReturn
    procedure(),               pointer                :: interruptProcedure
    integer                                           :: thisIndex,iHistory
    double precision                                  :: starFormationRate,stellarMassRate,fuelMassRate,fuelMass &
         &,massOutflowRate,diskMass,angularMomentumOutflowRate,transferRate,barInstabilityTimescale,gasMass,energyInputRate&
         &,diskDynamicalTime,massOutflowRateToHotHalo,massOutflowRateFromHalo,outflowToHotHaloFraction
    type(abundancesStructure), save                   :: fuelAbundances,stellarAbundancesRates,fuelAbundancesRates
    !$omp threadprivate(fuelAbundances,stellarAbundancesRates,fuelAbundancesRates)
    type(history)                                     :: historyTransferRate

    ! Get a local copy of the interrupt procedure.
    interruptProcedure => interruptProcedureReturn

    ! Compute star formation rate if this node has a disk.
    if (thisNode%componentExists(componentIndex)) then

       ! Check for a realistic disk, return immediately if disk is unphysical.
       if     (    Tree_Node_Disk_Angular_Momentum    (thisNode) < 0.0d0 &
            & .or. Tree_Node_Disk_Radius              (thisNode) < 0.0d0 &
            & .or. Tree_Node_Disk_Gas_Mass            (thisNode) < 0.0d0 &
            & ) return

       ! Compute the star formation rate.
       starFormationRate=Exponential_Disk_SFR(thisNode)
       
       ! Get the available fuel mass.
       fuelMass=Tree_Node_Disk_Gas_Mass_Exponential(thisNode)
       
       ! Find the metallicity of the fuel supply.
       call Tree_Node_Disk_Gas_Abundances_Exponential(thisNode,abundancesValue)
       call fuelAbundances%pack(abundancesValue)
       call fuelAbundances%massToMassFraction(fuelMass)

       ! Get the component index.
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       
       ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
       call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,componentTypeDisk,thisNode&
            &,thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex),stellarMassRate,stellarAbundancesRates&
            &,stellarLuminositiesRates,fuelMassRate,fuelAbundancesRates,energyInputRate)

       ! Adjust rates.
       call Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Exponential        (thisNode,interrupt,interruptProcedure,stellarMassRate         )
       call Tree_Node_Disk_Gas_Mass_Rate_Adjust_Exponential            (thisNode,interrupt,interruptProcedure,fuelMassRate            )
       call stellarAbundancesRates%unpack(abundancesValue)
       call Tree_Node_Disk_Stellar_Abundances_Rate_Adjust_Exponential  (thisNode,interrupt,interruptProcedure,abundancesValue         )
       call fuelAbundancesRates%unpack(abundancesValue)
       call Tree_Node_Disk_Gas_Abundances_Rate_Adjust_Exponential      (thisNode,interrupt,interruptProcedure,abundancesValue         )
       call Tree_Node_Disk_Stellar_Luminosities_Rate_Adjust_Exponential(thisNode,interrupt,interruptProcedure,stellarLuminositiesRates)

       ! Record the star formation history.
       call Star_Formation_History_Record(thisNode,thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)&
            &,fuelAbundances,starFormationRate)
       
       ! Find rate of outflow of material from the disk and pipe it to the outflowed reservoir.
       massOutflowRateToHotHalo=Star_Formation_Feedback_Disk_Outflow_Rate          (thisNode,starFormationRate,energyInputRate)
       massOutflowRateFromHalo =Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
       massOutflowRate         =massOutflowRateToHotHalo+massOutflowRateFromHalo
       if (massOutflowRate > 0.0d0) then
          ! Find the fraction of material which outflows to the hot halo.
          outflowToHotHaloFraction=massOutflowRateToHotHalo/massOutflowRate

          ! Get the masses of the disk.
          gasMass =        Tree_Node_Disk_Gas_Mass_Exponential    (thisNode)
          diskMass=gasMass+Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)
          
          ! Limit the outflow rate timescale to a multiple of the dynamical time.
          diskDynamicalTime=Mpc_per_km_per_s_To_Gyr*Tree_Node_Disk_Radius(thisNode)/Tree_Node_Disk_Velocity(thisNode)   
          massOutflowRate=min(massOutflowRate,gasMass/diskOutflowTimescaleMinimum/diskDynamicalTime)
          
          ! Compute the angular momentum outflow rate.
          if (diskMass > 0.0d0) then
             angularMomentumOutflowRate=massOutflowRate*Tree_Node_Disk_Angular_Momentum_Exponential(thisNode)/diskMass
          else
             angularMomentumOutflowRate=0.0d0
          end if
          if (gasMass  > 0.0d0) then
             call Tree_Node_Disk_Gas_Abundances_Exponential(thisNode,abundancesValue)
             call Abundances_Mass_To_Mass_Fraction(abundancesValue,gasMass)
             abundancesOutflowRate=massOutflowRate*abundancesValue
          else
             abundancesOutflowRate=0.0d0
          end if
          call Tree_Node_Hot_Halo_Outflow_Mass_To                     (thisNode,interrupt,interruptProcedure, &
               &  massOutflowRate           *outflowToHotHaloFraction)
          call Tree_Node_Disk_Gas_Mass_Rate_Adjust_Exponential        (thisNode,interrupt,interruptProcedure, &
               & -massOutflowRate                                    )
          call Tree_Node_Hot_Halo_Outflow_Angular_Momentum_To         (thisNode,interrupt,interruptProcedure, &
               &  angularMomentumOutflowRate*outflowToHotHaloFraction)
          call Tree_Node_Disk_Angular_Momentum_Rate_Adjust_Exponential(thisNode,interrupt,interruptProcedure, &
               & -angularMomentumOutflowRate                         )
          call Tree_Node_Hot_Halo_Outflow_Abundances_To               (thisNode,interrupt,interruptProcedure, &
               &  abundancesOutflowRate     *outflowToHotHaloFraction)
          call Tree_Node_Disk_Gas_Abundances_Rate_Adjust_Exponential  (thisNode,interrupt,interruptProcedure, &
               & -abundancesOutflowRate                              )
       end if

       ! Determine if the disk is bar unstable and, if so, the rate at which material is moved to the pseudo-bulge.
       barInstabilityTimescale=Bar_Instability_Timescale(thisNode)

       ! Negative timescale indicates no bar instability.
       if (barInstabilityTimescale >= 0.0d0) then
          ! Disk is unstable, so compute rates at which material is transferred to the spheroid.
          ! Gas mass.
          transferRate=max(0.0d0,Tree_Node_Disk_Gas_Mass_Exponential        (thisNode))/barInstabilityTimescale
          call Tree_Node_Disk_Gas_Mass_Rate_Adjust_Exponential              (thisNode,interrupt,interruptProcedure,-transferRate            )
          call Tree_Node_Spheroid_Gas_Mass_Rate_Adjust                      (thisNode,interrupt,interruptProcedure, transferRate            )
          ! Stellar mass.
          transferRate=max(0.0d0,Tree_Node_Disk_Stellar_Mass_Exponential    (thisNode))/barInstabilityTimescale
          call Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Exponential          (thisNode,interrupt,interruptProcedure,-transferRate            )
          call Tree_Node_Spheroid_Stellar_Mass_Rate_Adjust                  (thisNode,interrupt,interruptProcedure, transferRate            )
          ! Angular momentum.
          transferRate=max(0.0d0,Tree_Node_Disk_Angular_Momentum_Exponential(thisNode))/barInstabilityTimescale
          call Tree_Node_Disk_Angular_Momentum_Rate_Adjust_Exponential      (thisNode,interrupt,interruptProcedure,-transferRate            )
          call Tree_Node_Spheroid_Angular_Momentum_Rate_Adjust              (thisNode,interrupt,interruptProcedure, transferRate            )
          ! Gas abundances.
          call Tree_Node_Disk_Gas_Abundances_Exponential                    (thisNode,abundancesValue)
          abundancesTransferRate=max(0.0d0,abundancesValue)/barInstabilityTimescale
          call Tree_Node_Disk_Gas_Abundances_Rate_Adjust_Exponential        (thisNode,interrupt,interruptProcedure,-abundancesTransferRate  )
          call Tree_Node_Spheroid_Gas_Abundances_Rate_Adjust                (thisNode,interrupt,interruptProcedure, abundancesTransferRate  )
          ! Stellar abundances.
          call Tree_Node_Disk_Stellar_Abundances_Exponential                (thisNode,abundancesValue)
          abundancesTransferRate=max(0.0d0,abundancesValue)/barInstabilityTimescale
          call Tree_Node_Disk_Stellar_Abundances_Rate_Adjust_Exponential    (thisNode,interrupt,interruptProcedure,-abundancesTransferRate  )
          call Tree_Node_Spheroid_Stellar_Abundances_Rate_Adjust            (thisNode,interrupt,interruptProcedure, abundancesTransferRate  )
          ! Stellar luminosities.
          call Tree_Node_Disk_Stellar_Luminosities_Exponential              (thisNode,luminositiesValue)
          luminositiesTransferRate=max(0.0d0,luminositiesValue)/barInstabilityTimescale
          call Tree_Node_Disk_Stellar_Luminosities_Rate_Adjust_Exponential  (thisNode,interrupt,interruptProcedure,-luminositiesTransferRate)
          call Tree_Node_Spheroid_Stellar_Luminosities_Rate_Adjust          (thisNode,interrupt,interruptProcedure, luminositiesTransferRate)
          ! Stellar properties history.
          if (thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)%exists()) then
             thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
             historyTransferRate=thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)/barInstabilityTimescale
             thisNode                    %components(thisIndex)%instance(1)%histories(stellarHistoryIndex)%rates                          &
                  &             =thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)%rates-historyTransferRate%data
             call Tree_Node_Spheroid_Stellar_Properties_History_Rate_Adjust (thisNode,interrupt,interruptProcedure, historyTransferRate     )
             call historyTransferRate%destroy()
          end if
          ! Star formation history.
          if (thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%exists()) then
             thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
             historyTransferRate=thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)/barInstabilityTimescale
             thisNode                    %components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%rates                          &
                  &             =thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%rates-historyTransferRate%data
             call Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust     (thisNode,interrupt,interruptProcedure, historyTransferRate     )
          end if
          
       end if

    end if

    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure
    
    return
  end subroutine Tree_Node_Disk_Stellar_Mass_Rate_Compute_Exponential

  subroutine Tree_Node_Disk_Gas_Abundances_Exponential(thisNode,abundanceMasses)
    !% Return the node exponential disk gas abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: abundanceMasses(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       abundanceMasses(:)=thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd&
            &,propertyValue)
    else
       abundanceMasses(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Disk_Gas_Abundances_Exponential

  subroutine Tree_Node_Disk_Gas_Abundances_Set_Exponential(thisNode,abundanceMasses)
    !% Set the node exponential disk gas abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: abundanceMasses(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyValue)=abundanceMasses(:)
    return
  end subroutine Tree_Node_Disk_Gas_Abundances_Set_Exponential

  subroutine Tree_Node_Disk_Gas_Abundances_Rate_Adjust_Exponential(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node exponential disk gas abundance masses rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no exponential disk component currently exists and we have some cooling into it then interrupt and create an exponential
    ! disk.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Exponential_Disk_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Disk_Gas_Abundances_Rate_Adjust_Exponential

  subroutine Tree_Node_Disk_Stellar_Abundances_Exponential(thisNode,abundanceMasses)
    !% Return the node exponential disk stellar abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: abundanceMasses(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       abundanceMasses(:)=thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd&
            &,propertyValue)
    else
       abundanceMasses(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Disk_Stellar_Abundances_Exponential

  subroutine Tree_Node_Disk_Stellar_Abundances_Set_Exponential(thisNode,abundanceMasses)
    !% Set the node exponential disk stellar abundance masses.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: abundanceMasses(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyValue)=abundanceMasses(:)
    return
  end subroutine Tree_Node_Disk_Stellar_Abundances_Set_Exponential

  subroutine Tree_Node_Disk_Stellar_Abundances_Rate_Adjust_Exponential(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node exponential disk stellar abundance masses rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no exponential disk component currently exists and we have some cooling into it then interrupt and create an exponential
    ! disk.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Exponential_Disk_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Disk_Stellar_Abundances_Rate_Adjust_Exponential

  subroutine Tree_Node_Disk_Stellar_Luminosities_Exponential(thisNode,luminosities)
    !% Return the node exponential disk stellar luminosities.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: luminosities(:)
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       luminosities(:)=thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd&
            &,propertyValue)
    else
       luminosities(:)=0.0d0
    end if
    return
  end subroutine Tree_Node_Disk_Stellar_Luminosities_Exponential

  subroutine Tree_Node_Disk_Stellar_Luminosities_Set_Exponential(thisNode,luminosities)
    !% Set the node exponential disk stellar luminosities.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: luminosities(:)
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyValue)=luminosities(:)
    return
  end subroutine Tree_Node_Disk_Stellar_Luminosities_Set_Exponential

  subroutine Tree_Node_Disk_Stellar_Luminosities_Rate_Adjust_Exponential(thisNode,interrupt,interruptProcedure,rateAdjustments)
    !% Adjust the node exponential disk stellar luminosity rates of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustments(:)
    integer                                  :: thisIndex
    
    ! If no exponential disk component currently exists and we have some cooling into it then interrupt and create an exponential
    ! disk.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (any(rateAdjustments /= 0.0d0)) then
          interrupt=.true.
          interruptProcedure => Exponential_Disk_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyDerivative)&
         &+rateAdjustments(:)
    return
  end subroutine Tree_Node_Disk_Stellar_Luminosities_Rate_Adjust_Exponential

  double precision function Tree_Node_Disk_Angular_Momentum_Exponential(thisNode,instance)
    !% Return the node exponential disk gas mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       ! Force angular momentum to be positive. Can become negative due to rounding errors in ODE solver.
       if (thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue) < 0.0d0)&
            & thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)=0.0d0
       Tree_Node_Disk_Angular_Momentum_Exponential=thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)
    else
       Tree_Node_Disk_Angular_Momentum_Exponential=0.0d0
    end if
    return
  end function Tree_Node_Disk_Angular_Momentum_Exponential

  subroutine Tree_Node_Disk_Angular_Momentum_Set_Exponential(thisNode,angularMomentum,instance)
    !% Set the node exponential disk gas mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: angularMomentum
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyValue)=angularMomentum
    return
  end subroutine Tree_Node_Disk_Angular_Momentum_Set_Exponential

  subroutine Tree_Node_Disk_Angular_Momentum_Rate_Adjust_Exponential(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node exponential disk gas mass rate of change.
    use Cosmological_Parameters
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    ! If no exponential disk component currently exists and we have some cooling into it then interrupt and create an exponential
    ! disk.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Exponential_Disk_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyDerivative) &
         &=thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Disk_Angular_Momentum_Rate_Adjust_Exponential

  subroutine Tree_Node_Disk_Angular_Momentum_Rate_Compute_Exponential(thisNode,interrupt,interruptProcedure)
    !% Compute the exponential disk node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure

    ! No need to do anything here - the hot halo will pipe a cooling rate to us.
    return
  end subroutine Tree_Node_Disk_Angular_Momentum_Rate_Compute_Exponential

  function Tree_Node_Disk_Stellar_Properties_History_Exponential(thisNode)
    !% Return the disk stellar properties history.
    use Histories
    implicit none
    type(history)                          :: Tree_Node_Disk_Stellar_Properties_History_Exponential
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       Tree_Node_Disk_Stellar_Properties_History_Exponential=thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)
    else
       Tree_Node_Disk_Stellar_Properties_History_Exponential=nullHistory
    end if
    return
  end function Tree_Node_Disk_Stellar_Properties_History_Exponential

  subroutine Tree_Node_Disk_Stellar_Properties_History_Set_Exponential(thisNode,thisHistory)
    !% Set the node disk stellar properties history.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(history),           intent(in)    :: thisHistory
    integer                                :: thisIndex

    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    call thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)%destroy()
    thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex)=thisHistory
    return
  end subroutine Tree_Node_Disk_Stellar_Properties_History_Set_Exponential

  function Tree_Node_Disk_Star_Formation_History_Exponential(thisNode)
    !% Return the disk star formation history.
    use Histories
    implicit none
    type(history)                          :: Tree_Node_Disk_Star_Formation_History_Exponential
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       Tree_Node_Disk_Star_Formation_History_Exponential=thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)
    else
       Tree_Node_Disk_Star_Formation_History_Exponential=nullHistory
    end if
    return
  end function Tree_Node_Disk_Star_Formation_History_Exponential

  subroutine Tree_Node_Disk_Star_Formation_History_Set_Exponential(thisNode,thisHistory)
    !% Set the node disk star formation history.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(history),           intent(in)    :: thisHistory
    integer                                :: thisIndex

    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    call thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)%destroy()
    thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex)=thisHistory
    return
  end subroutine Tree_Node_Disk_Star_Formation_History_Set_Exponential

  !# <scaleSetTask>
  !#  <unitName>Exponential_Disk_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Exponential_Disk_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    use Abundances_Structure
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(treeNode),            pointer, intent(inout) :: thisNode
    double precision,          parameter              :: massMinimum           =1.0d0
    double precision,          parameter              :: angularMomentumMinimum=0.1d0
    double precision,          parameter              :: luminosityMinimum     =1.0d0
    type(abundancesStructure), save                   :: stellarAbundances
    !$omp threadprivate(stellarAbundances)
    integer                                           :: thisIndex
    double precision                                  :: mass,angularMomentum

    ! Determine if method is active and a disk component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)

       ! Set scale for angular momentum.
       angularMomentum=Tree_Node_Spheroid_Angular_Momentum(thisNode)+Tree_Node_Disk_Angular_Momentum(thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(angularMomentumIndex,propertyScale)=max(angularMomentum,angularMomentumMinimum)

       ! Set scale for gas mass.
       mass=Tree_Node_Spheroid_Gas_Mass(thisNode)+Tree_Node_Disk_Gas_Mass(thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyScale)=max(mass,massMinimum)

       ! Set scale for stellar mass.
       mass=Tree_Node_Spheroid_Stellar_Mass(thisNode)+Tree_Node_Disk_Stellar_Mass(thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyScale)=max(mass,massMinimum)

       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for gas abundances.
          call Tree_Node_Spheroid_Gas_Abundances(thisNode,abundancesSpheroid)
          call Tree_Node_Disk_Gas_Abundances    (thisNode,abundancesDisk    )
          thisNode%components(thisIndex)%instance(1)%properties(gasAbundancesIndex:gasAbundancesIndexEnd,propertyScale)=max(abundancesDisk+abundancesSpheroid,massMinimum)
          
          ! Set scale for stellar abundances.
          call Tree_Node_Spheroid_Stellar_Abundances(thisNode,abundancesSpheroid)
          call Tree_Node_Disk_Stellar_Abundances    (thisNode,abundancesDisk    )
          thisNode%components(thisIndex)%instance(1)%properties(stellarAbundancesIndex:stellarAbundancesIndexEnd,propertyScale)=max(abundancesDisk+abundancesSpheroid,massMinimum)
       end if

       ! Set scales for stellar luminosities if necessary.
       if (luminositiesCount > 0) then        
          ! Set scale for stellar luminosities.
          call Tree_Node_Spheroid_Stellar_Luminosities(thisNode,luminositiesSpheroid)
          call Tree_Node_Disk_Stellar_Luminosities    (thisNode,luminositiesDisk    )
          thisNode%components(thisIndex)%instance(1)%properties(stellarLuminositiesIndex:stellarLuminositiesIndexEnd,propertyScale)=max(luminositiesDisk+luminositiesSpheroid,luminosityMinimum)
       end if

       ! Set scales for stellar population properties and star formation histories.
       call Tree_Node_Disk_Stellar_Abundances_Exponential(thisNode,abundancesDisk)
       call stellarAbundances%pack(abundancesDisk)
       call Stellar_Population_Properties_Scales(thisNode%components(thisIndex)%instance(1)%histories(stellarHistoryIndex      ),Tree_Node_Disk_Stellar_Mass(thisNode),stellarAbundances)
       call Star_Formation_History_Scales       (thisNode%components(thisIndex)%instance(1)%histories(starFormationHistoryIndex),Tree_Node_Disk_Stellar_Mass(thisNode),stellarAbundances)
    end if
    return
  end subroutine Exponential_Disk_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Exponential_Disk_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Exponential_Disk_Satellite_Merging(thisNode)
    !% Transfer any exponential disk associated with {\tt thisNode} to its host halo.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer, intent(inout)       :: thisNode
    type(treeNode),   pointer                      :: hostNode
    double precision, dimension(abundancesCount)   :: thisAbundances,hostAbundances
    double precision, dimension(luminositiesCount) :: thisLuminosities,hostLuminosities
    type(history)                                  :: thisHistory,hostHistory
    double precision                               :: specificAngularMomentum

    ! Check that method is selected.
    if (methodSelected) then
       
       ! Find the node to merge with.
       call thisNode%mergesWith(hostNode)
       
       ! Get specific angular momentum of the disk material.

       if (Tree_Node_Disk_Gas_Mass(thisNode)+Tree_Node_Disk_Stellar_Mass(thisNode) > 0.0d0) then
          specificAngularMomentum=Tree_Node_Disk_Angular_Momentum_Exponential(thisNode)&
               &/(Tree_Node_Disk_Gas_Mass_Exponential(thisNode)+Tree_Node_Disk_Stellar_Mass_Exponential(thisNode))
       else
          specificAngularMomentum=0.0d0
       end if

       ! Move the gas component of the exponential disk to the host.
       select case (thisMergerGasMovesTo)
       case (movesToDisk)
          call Tree_Node_Disk_Gas_Mass_Set_Exponential        (hostNode, Tree_Node_Disk_Gas_Mass_Exponential(hostNode)         &
               &                                                        +Tree_Node_Disk_Gas_Mass_Exponential(thisNode)        )
          call Tree_Node_Disk_Gas_Abundances_Exponential      (hostNode, hostAbundances                                       )
          call Tree_Node_Disk_Gas_Abundances_Exponential      (thisNode, thisAbundances                                       )
          call Tree_Node_Disk_Gas_Abundances_Set_Exponential  (hostNode, hostAbundances+thisAbundances                        )
          call Tree_Node_Disk_Angular_Momentum_Set_Exponential(hostNode, Tree_Node_Disk_Angular_Momentum_Exponential(hostNode) &
               &+specificAngularMomentum*Tree_Node_Disk_Gas_Mass_Exponential(thisNode)                                        )
       case (movesToSpheroid)
          call Tree_Node_Spheroid_Gas_Mass_Set                (hostNode, Tree_Node_Spheroid_Gas_Mass        (hostNode) &
               &                                                        +Tree_Node_Disk_Gas_Mass_Exponential(thisNode)        )
          call Tree_Node_Spheroid_Gas_Abundances              (hostNode, hostAbundances                                       )
          call Tree_Node_Disk_Gas_Abundances_Exponential      (thisNode, thisAbundances                                       )
          call Tree_Node_Spheroid_Gas_Abundances_Set          (hostNode, hostAbundances+thisAbundances                        )
       case default
          call Galacticus_Error_Report('Exponential_Disk_Satellite_Merging','unrecognized movesTo descriptor')
       end select
       call Tree_Node_Disk_Gas_Mass_Set_Exponential           (thisNode, 0.0d0                                                )
       thisAbundances=0.0d0
       call Tree_Node_Disk_Gas_Abundances_Set_Exponential     (thisNode, thisAbundances                                       )

       ! Move the stellar component of the exponential disk to the host.
       select case (thisMergerStarsMoveTo)
       case (movesToDisk)
          call Tree_Node_Disk_Stellar_Mass_Set_Exponential        (hostNode, Tree_Node_Disk_Stellar_Mass_Exponential(hostNode)     &
               &                                                            +Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)    )
          call Tree_Node_Disk_Stellar_Abundances_Exponential      (hostNode, hostAbundances                                       )
          call Tree_Node_Disk_Stellar_Abundances_Exponential      (thisNode, thisAbundances                                       )
          call Tree_Node_Disk_Stellar_Abundances_Set_Exponential  (hostNode, hostAbundances+thisAbundances                        )
          call Tree_Node_Disk_Stellar_Luminosities_Exponential    (hostNode, hostLuminosities                                     )
          call Tree_Node_Disk_Stellar_Luminosities_Exponential    (thisNode, thisLuminosities                                     )
          call Tree_Node_Disk_Stellar_Luminosities_Set_Exponential(hostNode, hostLuminosities+thisLuminosities                    )
          call Tree_Node_Disk_Angular_Momentum_Set_Exponential    (hostNode, Tree_Node_Disk_Angular_Momentum_Exponential(hostNode) &
               &+specificAngularMomentum*Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)                                        )
          ! Also add stellar properties histories.
          thisHistory=Tree_Node_Disk_Stellar_Properties_History(thisNode)
          hostHistory=Tree_Node_Disk_Stellar_Properties_History(hostNode)
          call hostHistory%add(thisHistory)
          call Tree_Node_Disk_Stellar_Properties_History_Set(hostNode,hostHistory)
          call thisHistory%reset()
          call Tree_Node_Disk_Stellar_Properties_History_Set(thisNode,thisHistory)
          ! Also add star formation histories.
          thisHistory=Tree_Node_Disk_Star_Formation_History(thisNode)
          hostHistory=Tree_Node_Disk_Star_Formation_History(hostNode)
          call hostHistory%combine(thisHistory)
          call Tree_Node_Disk_Star_Formation_History_Set(hostNode,hostHistory)
          call thisHistory%reset()
          call Tree_Node_Disk_Star_Formation_History_Set(thisNode,thisHistory)
          call thisHistory%destroy()
          call hostHistory%destroy()
       case (movesToSpheroid)
          call Tree_Node_Spheroid_Stellar_Mass_Set                (hostNode, Tree_Node_Spheroid_Stellar_Mass        (hostNode)     &
               &                                                            +Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)    )
          call Tree_Node_Spheroid_Stellar_Abundances              (hostNode, hostAbundances                                       )
          call Tree_Node_Disk_Stellar_Abundances_Exponential      (thisNode, thisAbundances                                       )
          call Tree_Node_Spheroid_Stellar_Abundances_Set          (hostNode, hostAbundances+thisAbundances                        )
          call Tree_Node_Spheroid_Stellar_Luminosities            (hostNode, hostLuminosities                                     )
          call Tree_Node_Disk_Stellar_Luminosities_Exponential    (thisNode, thisLuminosities                                     )
          call Tree_Node_Spheroid_Stellar_Luminosities_Set        (hostNode, hostLuminosities+thisLuminosities                    )
          ! Also add stellar properties histories.
          thisHistory=Tree_Node_Disk_Stellar_Properties_History(thisNode)
          hostHistory=Tree_Node_Spheroid_Stellar_Properties_History(hostNode)
          call hostHistory%add(thisHistory)
          call Tree_Node_Spheroid_Stellar_Properties_History_Set(hostNode,hostHistory)
          call thisHistory%reset()
          call Tree_Node_Disk_Stellar_Properties_History_Set    (thisNode,thisHistory)
          ! Also add star formation histories.
          thisHistory=Tree_Node_Disk_Star_Formation_History(thisNode)
          hostHistory=Tree_Node_Spheroid_Star_Formation_History(hostNode)
          call hostHistory%combine(thisHistory)
          call Tree_Node_Spheroid_Star_Formation_History_Set(hostNode,hostHistory)
          call thisHistory%reset()
          call Tree_Node_Disk_Star_Formation_History_Set    (thisNode,thisHistory)
          call thisHistory%destroy()
          call hostHistory%destroy()
       case default
          call Galacticus_Error_Report('Exponential_Disk_Satellite_Merging','unrecognized movesTo descriptor')
       end select
       call Tree_Node_Disk_Stellar_Mass_Set_Exponential           (thisNode, 0.0d0                                                )
       thisAbundances=0.0d0
       call Tree_Node_Disk_Stellar_Abundances_Set_Exponential     (thisNode, thisAbundances                                       )
       thisLuminosities=0.0d0
       call Tree_Node_Disk_Stellar_Luminosities_Set_Exponential   (thisNode, thisLuminosities                                     )
       call Tree_Node_Disk_Angular_Momentum_Set_Exponential       (thisNode, 0.0d0                                                )
    end if
    return
  end subroutine Exponential_Disk_Satellite_Merging
  
  !# <enclosedMassTask>
  !#  <unitName>Exponential_Disk_Enclosed_Mass</unitName>
  !# </enclosedMassTask>
  subroutine Exponential_Disk_Enclosed_Mass(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for an exponential disk.
    use Galactic_Structure_Options
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType,weightBy,weightIndex
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass
    double precision                         :: fractionalRadius,diskRadius
    
    componentMass=0.0d0
    if (.not.methodSelected                           ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
    if (.not.thisNode%componentExists(componentIndex)) return
    select case (weightBy)
    case (weightByMass      )
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic)
          componentMass=Tree_Node_Disk_Gas_Mass_Exponential(thisNode)+Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)
       case (massTypeGaseous)
          componentMass=Tree_Node_Disk_Gas_Mass_Exponential(thisNode)
       case (massTypeStellar)
          componentMass=Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)
       end select
    case (weightByLuminosity)
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
          call Tree_Node_Disk_Stellar_Luminosities_Exponential(thisNode,luminositiesDisk)
          componentMass=luminositiesDisk(weightIndex)
       end select
    end select
    ! Return if no mass.
    if (componentMass <= 0.0d0)                         return
    ! Return if the total mass was requested.
    if (radius >= radiusLarge)                          return
    ! Compute the actual mass.
    diskRadius=Exponential_Disk_Radius(thisNode)
    if (diskRadius > 0.0d0) then
       fractionalRadius=radius/diskRadius
       componentMass=componentMass*(1.0d0-(1.0d0+fractionalRadius)*dexp(-fractionalRadius))
    end if
    return
  end subroutine Exponential_Disk_Enclosed_Mass

  !# <rotationCurveTask>
  !#  <unitName>Exponential_Disk_Rotation_Curve</unitName>
  !# </rotationCurveTask>
  subroutine Exponential_Disk_Rotation_Curve(thisNode,radius,massType,componentType,componentVelocity)
    !% Computes the rotation curve at a given radius for an exponential disk.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentVelocity
    double precision, parameter              :: fractionalRadiusMaximum=30.0d0
    double precision                         :: fractionalRadius,fractionalRadiusFactor,diskRadius,componentMass,halfRadius

    ! Set to zero by default.
    componentVelocity=0.0d0

    ! Compute if a disk is present.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then
       
       ! Get the mass of the disk.
       call Exponential_Disk_Enclosed_Mass(thisNode,radiusLarge,massType,componentType,weightByMass,weightIndexNull,componentMass)
       if (componentMass <= 0.0d0) return
       
       ! Compute the actual velocity.
       diskRadius=Exponential_Disk_Radius(thisNode)
       if (diskRadius > 0.0d0) then
          fractionalRadius=radius/diskRadius
          if (fractionalRadius > fractionalRadiusMaximum) then
             ! Beyond some maximum radius, approximate the disk as a spherical distribution to avoid evaluating Bessel functions for
             ! very large arguments.
             componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass/radius)
          else
             ! We are often called at precisely one scale length. Use pre-computed factors in that case.
             if (fractionalRadius == 1.0d0) then
                !$omp critical (ExponentialDiskFactorCompute)
                if (.not.scaleLengthFactorSet) then
                   halfRadius          =0.5d0
                   scaleLengthFactor   =Exponential_Disk_Rotation_Curve_Bessel_Factors(halfRadius)
                   scaleLengthFactorSet=.true.
                end if
                !$omp end critical (ExponentialDiskFactorCompute)
                fractionalRadiusFactor=scaleLengthFactor
             else
                halfRadius=0.5d0*fractionalRadius
                fractionalRadiusFactor=Exponential_Disk_Rotation_Curve_Bessel_Factors(halfRadius)
             end if
             componentVelocity=dsqrt(2.0d0*(gravitationalConstantGalacticus*componentMass/diskRadius)*fractionalRadiusFactor)
          end if
       end if
    end if
    return
  end subroutine Exponential_Disk_Rotation_Curve

  double precision function Exponential_Disk_Rotation_Curve_Bessel_Factors(halfRadius)
    !% Compute Bessel function factors appearing in the expression for an razor-thin exponential disk rotation curve.
    use Bessel_Functions
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Interpolation
    use FGSL
    implicit none
    double precision,        intent(in) :: halfRadius
    integer                             :: iPoint
    type(fgsl_interp),       save       :: interpolationObject
    type(fgsl_interp_accel), save       :: interpolationAccelerator
    logical,                 save       :: interpolationReset=.true.

    !$omp critical(Exponential_Disk_Rotation_Curve_Tabulate)
    if (.not.rotationCurveInitialized .or. halfRadius <  rotationCurveHalfRadiusMinimum .or. halfRadius > rotationCurveHalfRadiusMaximum) then
       ! Find the minimum and maximum half-radii to tabulate.
       rotationCurveHalfRadiusMinimum=min(rotationCurveHalfRadiusMinimum,0.5d0*halfRadius)
       rotationCurveHalfRadiusMaximum=max(rotationCurveHalfRadiusMaximum,2.0d0*halfRadius)

       ! Determine how many points to tabulate.
       rotationCurvePointsCount=int(dlog10(rotationCurveHalfRadiusMaximum/rotationCurveHalfRadiusMinimum)*dble(rotationCurvePointsPerDecade))+1

       ! Allocate table arrays.
       if (allocated(rotationCurveHalfRadius   )) call Dealloc_Array(rotationCurveHalfRadius   )
       if (allocated(rotationCurveBesselFactors)) call Dealloc_Array(rotationCurveBesselFactors)
       call Alloc_Array(rotationCurveHalfRadius   ,[rotationCurvePointsCount])
       call Alloc_Array(rotationCurveBesselFactors,[rotationCurvePointsCount])

       ! Create range of half-radii.
       rotationCurveHalfRadius=Make_Range(rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,rotationCurvePointsCount,rangeType=rangeTypeLogarithmic)

       ! Compute Bessel factors.
       do iPoint=1,rotationCurvePointsCount
          rotationCurveBesselFactors(iPoint)=rotationCurveHalfRadius(iPoint)**2&
               &*(Bessel_Function_I0(rotationCurveHalfRadius(iPoint))*Bessel_Function_K0(rotationCurveHalfRadius(iPoint)) &
               &-Bessel_Function_I1(rotationCurveHalfRadius(iPoint))*Bessel_Function_K1(rotationCurveHalfRadius(iPoint)))
       end do

       ! Reset the interpolations.
       call Interpolate_Done(interpolationObject,interpolationAccelerator,interpolationReset)
       interpolationReset=.true.

       ! Flag that the rotation curve is now initialized.
       rotationCurveInitialized=.true.
    end if
    !$omp end critical(Exponential_Disk_Rotation_Curve_Tabulate)

    ! Interpolate in the tabulated function.
    !$omp critical(Exponential_Disk_Rotation_Curve_Interpolate)
    Exponential_Disk_Rotation_Curve_Bessel_Factors=Interpolate(rotationCurvePointsCount,rotationCurveHalfRadius&
         &,rotationCurveBesselFactors,interpolationObject,interpolationAccelerator,halfRadius,reset=interpolationReset)
    !$omp end critical(Exponential_Disk_Rotation_Curve_Interpolate)

    return
  end function Exponential_Disk_Rotation_Curve_Bessel_Factors

  !# <radiusSolverPlausibility>
  !#  <unitName>Exponential_Disk_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Exponential_Disk_Radius_Solver_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
    !% Determines whether the disk is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: galaxyIsPhysicallyPlausible

    ! Return immediately if our method is not selected.
    if (.not.methodSelected) return
    
    if (Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)+Tree_Node_Disk_Gas_Mass_Exponential(thisNode) < -diskMassToleranceAbsolute) then
       galaxyIsPhysicallyPlausible=.false.
    else
       if (Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)+Tree_Node_Disk_Gas_Mass_Exponential(thisNode) > 0.0d0 .and.&
            & Tree_Node_Disk_Angular_Momentum_Exponential(thisNode) < 0.0d0) galaxyIsPhysicallyPlausible=.false.
    end if
    return
  end subroutine Exponential_Disk_Radius_Solver_Plausibility
  
  !# <radiusSolverTask>
  !#  <unitName>Exponential_Disk_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Exponential_Disk_Radius_Solver(thisNode,componentActive,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !% Interface for the size solver algorithm.
    implicit none
    type(treeNode),              pointer, intent(inout) :: thisNode
    logical,                              intent(out)   :: componentActive
    double precision,                     intent(out)   :: specificAngularMomentum
    procedure(double precision), pointer, intent(out)   :: Radius_Get,Velocity_Get
    procedure(),                 pointer, intent(out)   :: Radius_Set,Velocity_Set
    double precision                                    :: specificAngularMomentumMean,angularMomentum,diskMass

    ! Determine if thisNode has an active disk component supported by this module.    
    componentActive=methodSelected
    if (methodSelected) componentActive=thisNode%componentExists(componentIndex)    
    if (componentActive) then
       ! Get the angular momentum.
       angularMomentum=Tree_Node_Disk_Angular_Momentum_Exponential(thisNode)
       if (angularMomentum >= 0.0d0) then
          ! Compute the specific angular momentum at the scale radius, assuming a flat rotation curve.
          diskMass=Tree_Node_Disk_Gas_Mass_Exponential(thisNode)+Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)
          if (diskMass > 0.0d0) then
             specificAngularMomentumMean=angularMomentum/diskMass
          else
             specificAngularMomentumMean=0.0d0
          end if
          specificAngularMomentum=0.5d0*specificAngularMomentumMean
          ! Associate the pointers with the appropriate property routines.
          Radius_Get   => Exponential_Disk_Radius
          Radius_Set   => Exponential_Disk_Radius_Set
          Velocity_Get => Exponential_Disk_Velocity
          Velocity_Set => Exponential_Disk_Velocity_Set
       else
          call Exponential_Disk_Radius_Set  (thisNode,0.0d0)
          call Exponential_Disk_Velocity_Set(thisNode,0.0d0)
          componentActive=.false.
       end if
    end if
    return
  end subroutine Exponential_Disk_Radius_Solver


  double precision function Exponential_Disk_Radius(thisNode,instance)
    !% Return the scale radius of the exponential disk.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer,          intent(in), optional   :: instance
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       Exponential_Disk_Radius=thisNode%components(thisIndex)%instance(1)%data(radiusIndex)
    else
       Exponential_Disk_Radius=0.0d0
    end if
    return
  end function Exponential_Disk_Radius

  double precision function Exponential_Disk_Half_Mass_Radius(thisNode,instance)
    !% Return the half-mass radius of the exponential disk.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer,          intent(in), optional   :: instance
    double precision, parameter              :: halfMassRadiusToScaleRadius=1.678346990d0

    Exponential_Disk_Half_Mass_Radius=Exponential_Disk_Radius(thisNode)*halfMassRadiusToScaleRadius
    return
  end function Exponential_Disk_Half_Mass_Radius

  double precision function Exponential_Disk_SFR(thisNode,instance)
    !% Return the star formation rate of the exponential disk.
    use Star_Formation_Timescales_Disks
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer,          intent(in), optional   :: instance
    integer                                  :: thisIndex
    double precision                         :: starFormationTimescale,gasMass

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       
       ! Get the star formation timescale.
       starFormationTimescale=Star_Formation_Timescale_Disk(thisNode)
       
       ! Get the gas mass.
       gasMass=Tree_Node_Disk_Gas_Mass_Exponential(thisNode)
       
       ! If timescale is finite and gas mass is positive, then compute star formation rate.
       if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0) then
          Exponential_Disk_SFR=gasMass/starFormationTimescale
       else
          Exponential_Disk_SFR=0.0d0
       end if
    else
       Exponential_Disk_SFR=0.0d0
    end if
    return
  end function Exponential_Disk_SFR

  subroutine Exponential_Disk_Radius_Set(thisNode,radius,instance)
    !% Set the scale radius of the exponential disk.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: radius
    integer,          intent(in), optional   :: instance
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       thisNode%components(thisIndex)%instance(1)%data(radiusIndex)=max(radius,0.0d0)
    end if
    return
  end subroutine Exponential_Disk_Radius_Set

  double precision function Exponential_Disk_Velocity(thisNode,instance)
    !% Return the circular velocity of the exponential disk.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer,          intent(in), optional   :: instance
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       Exponential_Disk_Velocity=thisNode%components(thisIndex)%instance(1)%data(velocityIndex)
    else
       Exponential_Disk_Velocity=0.0d0
    end if
    return
  end function Exponential_Disk_Velocity

  subroutine Exponential_Disk_Velocity_Set(thisNode,velocity,instance)
    !% Set the circular velocity of the exponential disk.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: velocity
    integer,          intent(in), optional   :: instance
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       thisNode%components(thisIndex)%instance(1)%data(velocityIndex)=velocity
    end if
    return
  end subroutine Exponential_Disk_Velocity_Set

  integer function Tree_Node_Exponential_Disk_Index(thisNode)
    !% Ensure the exponential disk component exists and return its position in the components array.
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex

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
    Tree_Node_Exponential_Disk_Index=thisIndex
    return
  end function Tree_Node_Exponential_Disk_Index

  subroutine Exponential_Disk_Create(thisNode)
    !% Creates an exponential disk component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message
    integer                                     :: thisIndex

    ! Display a message.
    message='Creating exponential disk component for node '
    message=message//thisNode%index()
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Get the index of the component (which will also ensure that the component is created).
    thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
    return
  end subroutine Exponential_Disk_Create

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Disk_Exponential_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Disk_Exponential</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Disk_Exponential_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of exponential disk properties to be written to the \glc\ output file.
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
       doublePropertyNames   (doubleProperty)='diskGasMass'
       doublePropertyComments(doubleProperty)='Mass of gas in the exponential disk.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='diskStellarMass'
       doublePropertyComments(doubleProperty)='Mass of stars in the exponential disk at scale length.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='diskAngularMomentum'
       doublePropertyComments(doubleProperty)='Angular momentum of the exponential disk.'
       doublePropertyUnitsSI (doubleProperty)=massSolar*megaParsec*kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='diskScaleLength'
       doublePropertyComments(doubleProperty)='Radial scale length in the exponential disk.'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='diskCircularVelocity'
       doublePropertyComments(doubleProperty)='Circular velocity of the exponential disk at scale length.'
       doublePropertyUnitsSI (doubleProperty)=kilo
       do iAbundance=1,abundancesCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='diskGas'//Abundances_Names(iAbundance)
          doublePropertyComments(doubleProperty)='Disk gas phase abundance property.'
          doublePropertyUnitsSI (doubleProperty)=massSolar
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='diskStellar'//Abundances_Names(iAbundance)
          doublePropertyComments(doubleProperty)='Disk stellar abundance property.'
          doublePropertyUnitsSI (doubleProperty)=massSolar
       end do
       do iLuminosity=1,luminositiesCount
          if (Stellar_Population_Luminosities_Output(iLuminosity,time)) then
             doubleProperty=doubleProperty+1
             doublePropertyNames   (doubleProperty)='diskStellar'//Stellar_Population_Luminosities_Name(iLuminosity)
             doublePropertyComments(doubleProperty)='Disk stellar luminosity property.'
             doublePropertyUnitsSI (doubleProperty)=luminosityZeroPointAB
          end if
       end do
       if (diskOutputStarFormationRate) then
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='diskStarFormationRate'
          doublePropertyComments(doubleProperty)='Disk star formation rate.'
          doublePropertyUnitsSI (doubleProperty)=massSolar/gigaYear
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Disk_Exponential_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Disk_Exponential_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Disk_Exponential</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Disk_Exponential_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of exponential disk properties to be written to the the \glc\ output file.
    use Stellar_Population_Properties_Luminosities
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) then
       doublePropertyCount=doublePropertyCount+propertyCount+dataCount-luminositiesCount &
            & +Stellar_Population_Luminosities_Output_Count(time)
       if (diskOutputStarFormationRate) doublePropertyCount=doublePropertyCount+1
    end if

    return
  end subroutine Galacticus_Output_Tree_Disk_Exponential_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Disk_Exponential</unitName>
  !#  <sortName>Galacticus_Output_Tree_Disk_Exponential</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Disk_Exponential(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store exponential disk properties in the \glc\ output file buffers.
    use Stellar_Population_Properties_Luminosities
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
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Disk_Gas_Mass_Exponential(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Disk_Angular_Momentum_Exponential(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Exponential_Disk_Radius(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Exponential_Disk_Velocity(thisNode)
       call Tree_Node_Disk_Gas_Abundances_Exponential    (thisNode,gasAbundanceMasses)
       call Tree_Node_Disk_Stellar_Abundances_Exponential(thisNode,stellarAbundanceMasses)
       do iAbundance=1,abundancesCount
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=gasAbundanceMasses(iAbundance)
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=stellarAbundanceMasses(iAbundance)
       end do
       call Tree_Node_Disk_Stellar_Luminosities_Exponential(thisNode,stellarLuminosities)
       do iLuminosity=1,luminositiesCount
          if (Stellar_Population_Luminosities_Output(iLuminosity,time)) then
             doubleProperty=doubleProperty+1
             doubleBuffer(doubleBufferCount,doubleProperty)=stellarLuminosities(iLuminosity)
          end if
       end do
       if (diskOutputStarFormationRate) then
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Disk_SFR(thisNode)
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Disk_Exponential

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Exponential_Disk_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Exponential_Disk_Dump(thisNode)
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
          write (0,'(1x,a)') 'disk component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'disk gas mass:',Tree_Node_Disk_Gas_Mass_Exponential(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'disk stellar mass:',Tree_Node_Disk_Stellar_Mass_Exponential(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'disk angular momentum:',Tree_Node_Disk_Angular_Momentum_Exponential(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'disk scale length:',Exponential_Disk_Radius(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'disk circular velocity:',Exponential_Disk_Velocity(thisNode)
          call Tree_Node_Disk_Gas_Abundances_Exponential    (thisNode,gasAbundanceMasses    )
          call Tree_Node_Disk_Stellar_Abundances_Exponential(thisNode,stellarAbundanceMasses)
          do iAbundance=1,abundancesCount
             write (0,'(2x,a50,1x,e12.6)') 'disk gas '//char(Abundances_Names(iAbundance))//':',gasAbundanceMasses(iAbundance)
             write (0,'(2x,a50,1x,e12.6)') 'disk stellar '//char(Abundances_Names(iAbundance))//':',stellarAbundanceMasses(iAbundance)
          end do
          call Tree_Node_Disk_Stellar_Luminosities_Exponential(thisNode,stellarLuminosities)
          do iLuminosity=1,luminositiesCount
             write (0,'(2x,a50,1x,e12.6)') 'disk stellar '//char(Stellar_Population_Luminosities_Name(iLuminosity))//':',stellarLuminosities(iLuminosity)
          end do
       else
          write (0,'(1x,a)') 'disk component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Exponential_Disk_Dump

  !# <galacticusStateStoreTask>
  !#  <unitName>Tree_Node_Methods_Exponential_Disk_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Tree_Node_Methods_Exponential_Disk_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,scaleLengthFactor,scaleLengthFactorSet
    return
  end subroutine Tree_Node_Methods_Exponential_Disk_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Tree_Node_Methods_Exponential_Disk_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Tree_Node_Methods_Exponential_Disk_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,scaleLengthFactor,scaleLengthFactorSet
    ! Flag that the table is now uninitialized.
    rotationCurveInitialized=.false.
    return
  end subroutine Tree_Node_Methods_Exponential_Disk_State_Retrieve
  
  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Exponential_Disk_Star_Formation_History_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Exponential_Disk_Star_Formation_History_Output(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Store the star formation history in the output file.
    use Kind_Numbers
    use Tree_Nodes
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(treeNode),          intent(inout), pointer  :: thisNode
    integer,                 intent(in)              :: iOutput
    integer(kind=kind_int8), intent(in)              :: treeIndex
    logical,                 intent(in)              :: nodePassesFilter
    integer                                          :: thisIndex

    ! Output the star formation history if a disk exists for this component.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Exponential_Disk_Index(thisNode)
       call Star_Formation_History_Output(thisNode,nodePassesFilter,thisNode%components(thisIndex)%instance(1)&
            &%histories(starFormationHistoryIndex),iOutput ,treeIndex,'disk')
    end if
    return
  end subroutine Exponential_Disk_Star_Formation_History_Output

  !# <decodePropertyIdentifiersTask>
  !#  <unitName>Exponential_Disk_Property_Identifiers_Decode</unitName>
  !# </decodePropertyIdentifiersTask>
  subroutine Exponential_Disk_Property_Identifiers_Decode(propertyComponent,propertyObject,propertyIndex,matchedProperty,propertyName)
    !% Decodes property identifiers to property names for the exponential disk module.
    use ISO_Varying_String
    implicit none
    integer,              intent(in)    :: propertyComponent,propertyObject,propertyIndex
    logical,              intent(inout) :: matchedProperty
    type(varying_string), intent(inout) :: propertyName

    if (methodSelected.and..not.matchedProperty) then
       if (propertyComponent == componentIndex) then
          matchedProperty=.true.
          propertyName="exponentialDisk:"
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
  end subroutine Exponential_Disk_Property_Identifiers_Decode
  
end module Tree_Node_Methods_Exponential_Disk
