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


!% Contains a module of very simple disk tree node methods.

module Tree_Node_Methods_Very_Simple_Disk
  !% Implement very simple disk tree node methods.
  use ISO_Varying_String
  use Tree_Nodes
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Very_Simple_Disk_Initialize, Tree_Node_Disk_Post_Evolve_Very_Simple, Very_Simple_Disk_Scale_Set,&
       & Very_Simple_Disk_Satellite_Merging, Very_Simple_Disk_Enclosed_Mass, Galacticus_Output_Tree_Disk_Very_Simple_Names,&
       & Galacticus_Output_Tree_Disk_Very_Simple_Property_Count, Galacticus_Output_Tree_Disk_Very_Simple,&
       & Tree_Node_Methods_Very_Simple_Disk_Dump, Very_Simple_Disk_Property_Identifiers_Decode

  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

  ! Property indices.
  integer, parameter :: propertyCount=2, dataCount=0, historyCount=0
  integer, parameter :: gasMassIndex    =1
  integer, parameter :: stellarMassIndex=2

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Disk_Gas_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Disk_Stellar_Mass</methodName>
  !# </treeNodeMethodsPointer>

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Very_Simple_Disk_Initialize</unitName>
  !#  <optionName>treeNodeMethodDisk</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Very_Simple_Disk_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node very simple disk methods module.
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    use Galacticus_Error
    use Stellar_Population_Properties_Luminosities
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'verySimple') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Very simple disk method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Disk_Gas_Mass                  => Tree_Node_Disk_Gas_Mass_Very_Simple
       Tree_Node_Disk_Gas_Mass_Set              => Tree_Node_Disk_Gas_Mass_Set_Very_Simple
       Tree_Node_Disk_Gas_Mass_Rate_Adjust      => Tree_Node_Disk_Gas_Mass_Rate_Adjust_Very_Simple
       Tree_Node_Disk_Gas_Mass_Rate_Compute     => Tree_Node_Disk_Gas_Mass_Rate_Compute_Very_Simple

       Tree_Node_Disk_Stellar_Mass              => Tree_Node_Disk_Stellar_Mass_Very_Simple
       Tree_Node_Disk_Stellar_Mass_Set          => Tree_Node_Disk_Stellar_Mass_Set_Very_Simple
       Tree_Node_Disk_Stellar_Mass_Rate_Adjust  => Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Very_Simple
       Tree_Node_Disk_Stellar_Mass_Rate_Compute => Tree_Node_Disk_Stellar_Mass_Rate_Compute_Very_Simple

       ! Grab the cooling mass pipes from the hot halo component.
       if (.not.associated(Tree_Node_Hot_Halo_Cooling_Mass_To)) then
          Tree_Node_Hot_Halo_Cooling_Mass_To => Tree_Node_Disk_Gas_Mass_Rate_Adjust_Very_Simple
       else
          call Galacticus_Error_Report('Tree_Node_Methods_Very_Simple_Disk_Initialize','expected to find unclaimed hot halo cooling pipe')
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Very_Simple_Disk_Initialize
  
  !# <postEvolveTask>
  !# <unitName>Tree_Node_Disk_Post_Evolve_Very_Simple</unitName>
  !# </postEvolveTask>
  subroutine Tree_Node_Disk_Post_Evolve_Very_Simple(thisNode)
    !% Catch rounding errors in the very simple disk gas evolution.
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, save                   :: fractionalErrorMaximum=0.0d0
    integer                                  :: thisIndex
    double precision                         :: specificAngularMomentum,fractionalError,diskMass
    character(len=20)                        :: valueString
    type(varying_string)                     :: message

    ! Check if an very simple disk component exists.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then
       ! Get the index of the component.
       thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)

       ! Trap negative gas masses.
       if (Tree_Node_Disk_Gas_Mass(thisNode) < 0.0d0) then
          
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=dabs(thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue))&
               &/(thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)&
               &+dabs(thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)))
          !$omp critical (Very_Simple_Disk_Post_Evolve_Check)
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
          !$omp end critical (Very_Simple_Disk_Post_Evolve_Check)
          
          ! Get the specific angular momentum of the disk material
          diskMass= thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex    ,propertyValue) &
               &   +thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
          if (diskMass == 0.0d0) then
             thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)=0.0d0
          end if

          ! Reset the gas mass of the disk.
          thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)=0.0d0
       
       end if
       
    end if
    return
  end subroutine Tree_Node_Disk_Post_Evolve_Very_Simple

  double precision function Tree_Node_Disk_Gas_Mass_Very_Simple(thisNode,instance)
    !% Return the node very simple disk gas mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
       Tree_Node_Disk_Gas_Mass_Very_Simple=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)
    else
       Tree_Node_Disk_Gas_Mass_Very_Simple=0.0d0
    end if
    return
  end function Tree_Node_Disk_Gas_Mass_Very_Simple

  subroutine Tree_Node_Disk_Gas_Mass_Set_Very_Simple(thisNode,mass,instance)
    !% Set the node very simple disk gas mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Disk_Gas_Mass_Set_Very_Simple

  subroutine Tree_Node_Disk_Gas_Mass_Rate_Adjust_Very_Simple(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node very simple disk gas mass rate of change.
    use Cosmological_Parameters
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    ! If no very simple disk component currently exists and we have some cooling into it then interrupt and create an very simple
    ! disk.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Very_Simple_Disk_Create
       end if
       return
    end if

    thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Disk_Gas_Mass_Rate_Adjust_Very_Simple

  subroutine Tree_Node_Disk_Gas_Mass_Rate_Compute_Very_Simple(thisNode,interrupt,interruptProcedure)
    !% Compute the very simple disk node mass rate of change.
    use Cosmological_Parameters
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure

    ! No need to do anything here - the hot halo will pipe a cooling rate to us.
    return
  end subroutine Tree_Node_Disk_Gas_Mass_Rate_Compute_Very_Simple

  double precision function Tree_Node_Disk_Stellar_Mass_Very_Simple(thisNode,instance)
    !% Return the node very simple disk stellar mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
       Tree_Node_Disk_Stellar_Mass_Very_Simple=thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)
    else
       Tree_Node_Disk_Stellar_Mass_Very_Simple=0.0d0
    end if
    return
  end function Tree_Node_Disk_Stellar_Mass_Very_Simple

  subroutine Tree_Node_Disk_Stellar_Mass_Set_Very_Simple(thisNode,mass,instance)
    !% Set the node very simple disk stellar mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Disk_Stellar_Mass_Set_Very_Simple

  subroutine Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Very_Simple(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node very simple disk stellar mass rate of change.
    use Cosmological_Parameters
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    ! If no very simple disk component currently exists and we have some cooling into it then interrupt and create an very simple
    ! disk.
    if (.not.thisNode%componentExists(componentIndex)) then
       if (rateAdjustment /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Very_Simple_Disk_Create
       end if
       return
    end if
    thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Very_Simple

  subroutine Tree_Node_Disk_Stellar_Mass_Rate_Compute_Very_Simple(thisNode,interrupt,interruptProcedureReturn)
    !% Compute the very simple disk node mass rate of change.
    use Star_Formation_Feedback_Disks
    use Stellar_Feedback
    implicit none
    type(treeNode),            pointer, intent(inout) :: thisNode
    logical,                            intent(inout) :: interrupt
    procedure(),               pointer, intent(inout) :: interruptProcedureReturn
    procedure(),               pointer                :: interruptProcedure
    integer                                           :: thisIndex
    double precision                                  :: starFormationRate,fuelMass &
         &,massOutflowRate,gasMass,diskMass,energyInputRate
  
    ! Get a local copy of the interrupt procedure.
    interruptProcedure => interruptProcedureReturn

    ! Compute star formation rate if this node has a disk.
    if (thisNode%componentExists(componentIndex)) then

       ! Check for a realistic disk, return immediately if disk is unphysical.
       if (Tree_Node_Disk_Gas_Mass(thisNode) < 0.0d0) return

       ! Compute the star formation rate.
       starFormationRate=Very_Simple_Disk_SFR(thisNode)
       
       ! Get the available fuel mass.
       fuelMass=Tree_Node_Disk_Gas_Mass_Very_Simple(thisNode)
       
       ! Get the component index.
       thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
       
       ! Adjust rates.
       call Tree_Node_Disk_Stellar_Mass_Rate_Adjust_Very_Simple(thisNode,interrupt,interruptProcedure,+starFormationRate)
       call Tree_Node_Disk_Gas_Mass_Rate_Adjust_Very_Simple    (thisNode,interrupt,interruptProcedure,-starFormationRate)

       ! Find rate of outflow of material from the disk and pipe it to the outflowed reservoir.
       energyInputRate=feedbackEnergyInputAtInfinityCanonical*starFormationRate
       massOutflowRate=Star_Formation_Feedback_Disk_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
     
       if (massOutflowRate > 0.0d0) then
        
          ! Get the masses of the disk.
          gasMass =        Tree_Node_Disk_Gas_Mass_Very_Simple    (thisNode)
          diskMass=gasMass+Tree_Node_Disk_Stellar_Mass_Very_Simple(thisNode)
          
          call Tree_Node_Hot_Halo_Outflow_Mass_To             (thisNode,interrupt,interruptProcedure,+massOutflowRate)
          call Tree_Node_Disk_Gas_Mass_Rate_Adjust_Very_Simple(thisNode,interrupt,interruptProcedure,-massOutflowRate)
       end if

    end if

    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure
    
    return
  end subroutine Tree_Node_Disk_Stellar_Mass_Rate_Compute_Very_Simple

  !# <scaleSetTask>
  !#  <unitName>Very_Simple_Disk_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Very_Simple_Disk_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type(treeNode),            pointer, intent(inout) :: thisNode
    double precision,          parameter              :: massMinimum=1.0d0
    integer                                           :: thisIndex
    double precision                                  :: mass

    ! Determine if method is active and a disk component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)

       ! Set scale for gas mass.
       mass=Tree_Node_Disk_Gas_Mass(thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(gasMassIndex,propertyScale)=max(mass,massMinimum)

       ! Set scale for stellar mass.
       mass=Tree_Node_Disk_Stellar_Mass(thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(stellarMassIndex,propertyScale)=max(mass,massMinimum)

    end if
    return
  end subroutine Very_Simple_Disk_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Very_Simple_Disk_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Very_Simple_Disk_Satellite_Merging(thisNode)
    !% Transfer any very simple disk associated with {\tt thisNode} to its host halo.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer, intent(inout)       :: thisNode
    type(treeNode),   pointer                      :: hostNode

    ! Check that method is selected.
    if (methodSelected) then
       
       ! Find the node to merge with.
       call thisNode%mergesWith(hostNode)
       
       ! Move the gas component of the very simple disk to the host.
       select case (thisMergerGasMovesTo)
       case (movesToDisk)
          call Tree_Node_Disk_Gas_Mass_Set_Very_Simple        (hostNode, Tree_Node_Disk_Gas_Mass_Very_Simple(hostNode)         &
               &                                                        +Tree_Node_Disk_Gas_Mass_Very_Simple(thisNode)        )
       case (movesToSpheroid)
          call Galacticus_Error_Report('Very_Simple_Disk_Satellite_Merging','this component does not work with spheroids')
       case default
          call Galacticus_Error_Report('Very_Simple_Disk_Satellite_Merging','unrecognized movesTo descriptor')
       end select
       call Tree_Node_Disk_Gas_Mass_Set_Very_Simple           (thisNode, 0.0d0                                                )

       ! Move the stellar component of the very simple disk to the host.
       select case (thisMergerStarsMoveTo)
       case (movesToDisk)
          call Tree_Node_Disk_Stellar_Mass_Set_Very_Simple        (hostNode, Tree_Node_Disk_Stellar_Mass_Very_Simple(hostNode)     &
               &                                                            +Tree_Node_Disk_Stellar_Mass_Very_Simple(thisNode)    )
       case (movesToSpheroid)
          call Galacticus_Error_Report('Very_Simple_Disk_Satellite_Merging','this component does not work with spheroids')
       case default
          call Galacticus_Error_Report('Very_Simple_Disk_Satellite_Merging','unrecognized movesTo descriptor')
       end select
       call Tree_Node_Disk_Stellar_Mass_Set_Very_Simple           (thisNode, 0.0d0                                                )
    end if
    return
  end subroutine Very_Simple_Disk_Satellite_Merging
  
  !# <enclosedMassTask>
  !#  <unitName>Very_Simple_Disk_Enclosed_Mass</unitName>
  !# </enclosedMassTask>
  subroutine Very_Simple_Disk_Enclosed_Mass(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for an very simple disk.
    use Galactic_Structure_Options
    use Galacticus_Error
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType,weightBy,weightIndex
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass
    
    componentMass=0.0d0
    if (.not.methodSelected                           ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
    if (.not.thisNode%componentExists(componentIndex)) return
    select case (weightBy)
    case (weightByMass      )
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic)
          componentMass=Tree_Node_Disk_Gas_Mass_Very_Simple(thisNode)+Tree_Node_Disk_Stellar_Mass_Very_Simple(thisNode)
       case (massTypeGaseous)
          componentMass=Tree_Node_Disk_Gas_Mass_Very_Simple(thisNode)
       case (massTypeStellar)
          componentMass=Tree_Node_Disk_Stellar_Mass_Very_Simple(thisNode)
       end select
    case default
       call Galacticus_Error_Report('Very_Simple_Disk_Enclosed_Mass','this component does not track luminosity')
    end select
    ! Return if no mass.
    if (componentMass <= 0.0d0)                         return
    ! Return if the total mass was requested.
    if (radius >= radiusLarge)                          return
    ! Otherwise we have an error.
    call Galacticus_Error_Report('Very_Simple_Disk_Enclosed_Mass','this component does not specify a mass profile')
    return
  end subroutine Very_Simple_Disk_Enclosed_Mass

  double precision function Very_Simple_Disk_SFR(thisNode,instance)
    !% Return the star formation rate of the very simple disk.
    use Star_Formation_Timescales_Disks
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer,          intent(in), optional   :: instance
    integer                                  :: thisIndex
    double precision                         :: starFormationTimescale,gasMass

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
       
       ! Get the star formation timescale.
       starFormationTimescale=Star_Formation_Timescale_Disk(thisNode)
       
       ! Get the gas mass.
       gasMass=Tree_Node_Disk_Gas_Mass_Very_Simple(thisNode)
       
       ! If timescale is finite and gas mass is positive, then compute star formation rate.
       if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0) then
          Very_Simple_Disk_SFR=gasMass/starFormationTimescale
       else
          Very_Simple_Disk_SFR=0.0d0
       end if
    else
       Very_Simple_Disk_SFR=0.0d0
    end if
    return
  end function Very_Simple_Disk_SFR

  integer function Tree_Node_Very_Simple_Disk_Index(thisNode)
    !% Ensure the very simple disk component exists and return its position in the components array.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex

    if (.not.thisNode%componentExists(componentIndex)) then
       ! Create the component.
       call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
       ! Get the index for this component.
       thisIndex=thisNode%componentIndex(componentIndex)
    else
       ! Get the index for this component.
       thisIndex=thisNode%componentIndex(componentIndex)
    end if
    Tree_Node_Very_Simple_Disk_Index=thisIndex
    return
  end function Tree_Node_Very_Simple_Disk_Index

  subroutine Very_Simple_Disk_Create(thisNode)
    !% Creates an very simple disk component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message
    integer                                     :: thisIndex

    ! Display a message.
    message='Creating very simple disk component for node '
    message=message//thisNode%index()
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Get the index of the component (which will also ensure that the component is created).
    thisIndex=Tree_Node_Very_Simple_Disk_Index(thisNode)
    return
  end subroutine Very_Simple_Disk_Create

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Disk_Very_Simple_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Disk_Very_Simple</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Disk_Very_Simple_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of very simple disk properties to be written to the \glc\ output file.
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

    if (methodSelected) then
       !@ <outputPropertyGroup>
       !@   <name>disk</name>
       !@   <description>Disk properites</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>diskGasMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Mass of gas in the very simple disk.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>disk</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='diskGasMass'
       doublePropertyComments(doubleProperty)='Mass of gas in the very simple disk.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>diskStellarMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Mass of stars in the very simple disk at scale length.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>disk</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='diskStellarMass'
       doublePropertyComments(doubleProperty)='Mass of stars in the very simple disk.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
    end if
    return
  end subroutine Galacticus_Output_Tree_Disk_Very_Simple_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Disk_Very_Simple_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Disk_Very_Simple</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Disk_Very_Simple_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of very simple disk properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Galacticus_Output_Tree_Disk_Very_Simple_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Disk_Very_Simple</unitName>
  !#  <sortName>Galacticus_Output_Tree_Disk_Very_Simple</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Disk_Very_Simple(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store very simple disk properties in the \glc\ output file buffers.
    use Kind_Numbers
    implicit none
    double precision,        intent(in)                   :: time
    type(treeNode),          intent(inout), pointer       :: thisNode
    integer,                 intent(inout)                :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)                :: integerBuffer(:,:)
    double precision,        intent(inout)                :: doubleBuffer(:,:)
 
    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Disk_Gas_Mass_Very_Simple(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Disk_Stellar_Mass_Very_Simple(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Disk_Very_Simple

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Very_Simple_Disk_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Very_Simple_Disk_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    use ISO_Varying_String
    implicit none
    type(treeNode),   intent(inout), pointer       :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)') 'disk component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'disk gas mass:',Tree_Node_Disk_Gas_Mass_Very_Simple(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'disk stellar mass:',Tree_Node_Disk_Stellar_Mass_Very_Simple(thisNode)
        else
          write (0,'(1x,a)') 'disk component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Very_Simple_Disk_Dump

  !# <decodePropertyIdentifiersTask>
  !#  <unitName>Very_Simple_Disk_Property_Identifiers_Decode</unitName>
  !# </decodePropertyIdentifiersTask>
  subroutine Very_Simple_Disk_Property_Identifiers_Decode(propertyComponent,propertyObject,propertyIndex,matchedProperty,propertyName)
    !% Decodes property identifiers to property names for the very simple disk module.
    use ISO_Varying_String
    implicit none
    integer,              intent(in)    :: propertyComponent,propertyObject,propertyIndex
    logical,              intent(inout) :: matchedProperty
    type(varying_string), intent(inout) :: propertyName

    if (methodSelected.and..not.matchedProperty) then
       if (propertyComponent == componentIndex) then
          matchedProperty=.true.
          propertyName="verySimpleDisk:"
          select case (propertyObject)
          case (objectTypeProperty)
             if      (propertyIndex == gasMassIndex    ) then
                propertyName=propertyName//":gasMass"
             else if (propertyIndex == stellarMassIndex) then
                propertyName=propertyName//":stellarMass"
             end if
          end select
       end if
    end if

    return
  end subroutine Very_Simple_Disk_Property_Identifiers_Decode
  
end module Tree_Node_Methods_Very_Simple_Disk
