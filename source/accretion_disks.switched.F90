!% Contains a module which implements calculations of properties of accretion disks which switch between thin and ADAF depening on
!% the accretion rate.

module Accretion_Disks_Switched
  !% Implements calculations of properties of accretion disks which switch between thin and ADAF depening on
  !% the accretion rate.
  use Accretion_Disks_ADAF
  use Accretion_Disks_Shakura_Sunyaev
  private
  public :: Accretion_Disks_Switched_Initialize

  ! Values used to indicate which type of accretion disk is being used.
  integer, parameter :: accretionDiskThin=0
  integer, parameter :: accretionDiskADAF=1

  ! Parameters controlling the range of accretion rates over which the accretion disk will be an ADAF.
  double precision :: accretionRateAdafMinimum,accretionRateAdafMaximum
  logical          :: accretionRateAdafMinimumExists,accretionRateAdafMaximumExists

contains

  !# <accretionDisksMethod>
  !#  <unitName>Accretion_Disks_Switched_Initialize</unitName>
  !# </accretionDisksMethod>
  subroutine Accretion_Disks_Switched_Initialize(accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get&
       &,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: accretionDisksMethod
    procedure(),          pointer, intent(inout) :: Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get&
         &,Accretion_Disk_Jet_Power_Get
    character(len=30)                            :: accretionRateAdaf
    
    if (accretionDisksMethod == 'switched') then
       Accretion_Disk_Radiative_Efficiency_Get => Accretion_Disk_Radiative_Efficiency_Switched
       Black_Hole_Spin_Up_Rate_Get             => Black_Hole_Spin_Up_Rate_Switched
       Accretion_Disk_Jet_Power_Get            => Accretion_Disk_Jet_Power_Switched
       !@ <inputParameter>
       !@   <name>accretionRateAdafMinimum</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The accretion rate (in Eddington units) below which a switched accretion disk stop being an ADAF.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionRateAdafMinimum",accretionRateAdaf,defaultValue='0.01d0')
       if (trim(accretionRateAdaf) == "none") then
          accretionRateAdafMinimumExists=.false.
       else
          accretionRateAdafMinimumExists=.true.
          read (accretionRateAdaf,*) accretionRateAdafMinimum
       end if
       !@ <inputParameter>
       !@   <name>accretionRateAdafMaximum</name>
       !@   <defaultValue>0.3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The accretion rate (in Eddington units) above which a switched accretion disk stop being an ADAF.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionRateAdafMaximum",accretionRateAdaf,defaultValue="0.30d0")
       if (trim(accretionRateAdaf) == "none") then
          accretionRateAdafMaximumExists=.false.
       else
          accretionRateAdafMaximumExists=.true.
          read (accretionRateAdaf,*) accretionRateAdafMaximum
       end if
    end if
    return
  end subroutine Accretion_Disks_Switched_Initialize

  double precision function Accretion_Disk_Radiative_Efficiency_Switched(thisNode,massAccretionRate)
    !% Computes the radiative efficiency for a switching accretion disk.
    use Black_Hole_Fundamentals
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    select case (Accretion_Disk_Switched_Type(thisNode,massAccretionRate))
    case (accretionDiskThin)
       Accretion_Disk_Radiative_Efficiency_Switched=Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev(thisNode,massAccretionRate)
    case (accretionDiskADAF)
       Accretion_Disk_Radiative_Efficiency_Switched=Accretion_Disk_Radiative_Efficiency_ADAF(thisNode,massAccretionRate)
    end select

    return
  end function Accretion_Disk_Radiative_Efficiency_Switched

  double precision function Accretion_Disk_Jet_Power_Switched(thisNode,massAccretionRate)
    !% Computes the jet power for a switching accretion disk.
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    select case (Accretion_Disk_Switched_Type(thisNode,massAccretionRate))
    case (accretionDiskThin)
       Accretion_Disk_Jet_Power_Switched=Accretion_Disk_Jet_Power_Shakura_Sunyaev(thisNode,massAccretionRate)
    case (accretionDiskADAF)
       Accretion_Disk_Jet_Power_Switched=Accretion_Disk_Jet_Power_ADAF(thisNode,massAccretionRate)
    end select

    return
  end function Accretion_Disk_Jet_Power_Switched

  double precision function Black_Hole_Spin_Up_Rate_Switched(thisNode,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisNode} due to accretion from a switching accretion disk.
    !% disk.
    use Tree_Nodes
    use Tree_Node_Methods
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    select case (Accretion_Disk_Switched_Type(thisNode,massAccretionRate))
    case (accretionDiskThin)
       Black_Hole_Spin_Up_Rate_Switched=Black_Hole_Spin_Up_Rate_Shakura_Sunyaev(thisNode,massAccretionRate)
    case (accretionDiskADAF)
       Black_Hole_Spin_Up_Rate_Switched=Black_Hole_Spin_Up_Rate_ADAF(thisNode,massAccretionRate)
    end select

    return
  end function Black_Hole_Spin_Up_Rate_Switched

  integer function Accretion_Disk_Switched_Type(thisNode,massAccretionRate)
    !% Decide which type of accretion disk to use.
    use Tree_Nodes
    use Tree_Node_Methods
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision                         :: eddingtonAccretionRate,massAccretionRateDimensionless

    ! Get the Eddington accretion rate.
    eddingtonAccretionRate=Black_Hole_Eddington_Accretion_Rate(thisNode)

    ! Check that a black hole is present.
    if (eddingtonAccretionRate > 0.0d0) then

       ! Compute the accretion rate in Eddington units.
       massAccretionRateDimensionless=massAccretionRate/eddingtonAccretionRate
       
       ! Decide which type of accretion disk to use.
       if (        (accretionRateAdafMinimumExists .and. massAccretionRateDimensionless < accretionRateAdafMinimum) &
            & .or. (accretionRateAdafMaximumExists .and. massAccretionRateDimensionless > accretionRateAdafMaximum)) then
          Accretion_Disk_Switched_Type=accretionDiskThin
       else
          Accretion_Disk_Switched_Type=accretionDiskADAF
       end if

    else
       
       ! No black hole present: assume a thin disk.
       Accretion_Disk_Switched_Type=accretionDiskThin

    end if

    return
  end function Accretion_Disk_Switched_Type

end module Accretion_Disks_Switched
