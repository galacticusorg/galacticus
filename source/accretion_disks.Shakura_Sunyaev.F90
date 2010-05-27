!% Contains a module which implements calculations of properties of thin Shakura-Sunyaev accretion disks.

module Accretion_Disks_Shakura_Sunyaev
  !% Implements calculations of properties of thin Shakura-Sunyaev accretion disks.
  private
  public :: Accretion_Disks_Shakura_Sunyaev_Initialize, Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev,&
       & Black_Hole_Spin_Up_Rate_Shakura_Sunyaev, Accretion_Disk_Jet_Power_Shakura_Sunyaev

contains

  !# <accretionDisksMethod>
  !#  <unitName>Accretion_Disks_Shakura_Sunyaev_Initialize</unitName>
  !# </accretionDisksMethod>
  subroutine Accretion_Disks_Shakura_Sunyaev_Initialize(accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get&
       &,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: accretionDisksMethod
    procedure(),          pointer, intent(inout) :: Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get&
         &,Accretion_Disk_Jet_Power_Get
    
    if (accretionDisksMethod == 'Shakura-Sunyaev') then
       Accretion_Disk_Radiative_Efficiency_Get => Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev
       Black_Hole_Spin_Up_Rate_Get             => Black_Hole_Spin_Up_Rate_Shakura_Sunyaev
       Accretion_Disk_Jet_Power_Get            => Accretion_Disk_Jet_Power_Shakura_Sunyaev
    end if
    return
  end subroutine Accretion_Disks_Shakura_Sunyaev_Initialize

  double precision function Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev(thisNode,massAccretionRate)
    !% Computes the radiative efficiency for a Shakura-Sunyaev (thin) accretion disk.
    use Black_Hole_Fundamentals
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev=1.0d0-Black_Hole_ISCO_Specific_Energy(thisNode,units=unitsGravitational,orbit=orbitPrograde)
    return
  end function Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev

  double precision function Accretion_Disk_Jet_Power_Shakura_Sunyaev(thisNode,massAccretionRate)
    !% Computes the jet power for a Shakura-Sunyaev (thin) accretion disk.
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    ! Assume no jet is produced.
    Accretion_Disk_Jet_Power_Shakura_Sunyaev=0.0d0
    return
  end function Accretion_Disk_Jet_Power_Shakura_Sunyaev

  double precision function Black_Hole_Spin_Up_Rate_Shakura_Sunyaev(thisNode,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisNode} due to accretion from a Shakura-Sunyaev (thin) accretion
    !% disk.
    use Tree_Nodes
    use Tree_Node_Methods
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision                         :: accretionEfficiency,spinToMassRateOfChangeRatio

    spinToMassRateOfChangeRatio=Black_Hole_ISCO_Specific_Angular_Momentum(thisNode,units=unitsGravitational,orbit=orbitPrograde)&
         &-2.0d0*Tree_Node_Black_Hole_Spin(thisNode)*Black_Hole_ISCO_Specific_Energy(thisNode,units=unitsGravitational,orbit&
         &=orbitPrograde)
    Black_Hole_Spin_Up_Rate_Shakura_Sunyaev=spinToMassRateOfChangeRatio*massAccretionRate/Tree_Node_Black_Hole_Mass(thisNode)
    return
  end function Black_Hole_Spin_Up_Rate_Shakura_Sunyaev

end module Accretion_Disks_Shakura_Sunyaev
