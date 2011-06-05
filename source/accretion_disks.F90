!% Contains a module which implements calculations related to accretion disks.

module Accretion_Disks
  !% Implements calculations related to accretion disks.
  use ISO_Varying_String
  private
  public :: Accretion_Disk_Radiative_Efficiency,Black_Hole_Spin_Up_Rate,Accretion_Disk_Jet_Power

  ! Flag to indicate if this module has been initialized.  
  logical                                      :: accretionDisksInitialized=.false.

  ! Name of mass movement method used.
  type(varying_string)                         :: accretionDisksMethod

  ! Pointer to the subroutine that returns descriptors for mass movement.
  procedure(Accretion_Disk_Radiative_Efficiency), pointer :: Accretion_Disk_Radiative_Efficiency_Get => null()
  procedure(Black_Hole_Spin_Up_Rate),             pointer :: Black_Hole_Spin_Up_Rate_Get             => null()
  procedure(Accretion_Disk_Jet_Power),            pointer :: Accretion_Disk_Jet_Power_Get            => null()

contains

  subroutine Accretion_Disks_Initialize
    !% Initalize the accretion disk module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="accretionDisksMethod" type="moduleUse">
    include 'accretion_disks.modules.inc'
    !# </include>
    implicit none

    !$omp critical(accretionDisksInitialize)
    if (.not.accretionDisksInitialized) then
       ! Do the binary black hole merger method parameter.
       !@ <inputParameter>
       !@   <name>accretionDisksMethod</name>
       !@   <defaultValue>switched</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Selects which accretion disk method should be used.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('accretionDisksMethod',accretionDisksMethod,defaultValue='switched')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="accretionDisksMethod" type="code" action="subroutine">
       !#  <subroutineArgs>accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get</subroutineArgs>
       include 'accretion_disks.inc'
       !# </include>
       if (.not.(associated(Accretion_Disk_Radiative_Efficiency_Get).and.associated(Black_Hole_Spin_Up_Rate_Get) &
            & .and.associated(Accretion_Disk_Jet_Power_Get))) call&
            & Galacticus_Error_Report('Accretion_Disks_Initialize','method ' //char(accretionDisksMethod)//' is unrecognized')
       ! Flag that the module is now initialized.
       accretionDisksInitialized=.true.
    end if
    !$omp end critical(accretionDisksInitialized)

    return
  end subroutine Accretion_Disks_Initialize
  
  double precision function Accretion_Disk_Radiative_Efficiency(thisNode,massAccretionRate)
    !% Computes the radiative efficiency for an accretion disk.
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    ! Ensure the module is initalized.
    call Accretion_Disks_Initialize

    ! Get the radiative efficiency.
    Accretion_Disk_Radiative_Efficiency=Accretion_Disk_Radiative_Efficiency_Get(thisNode,massAccretionRate)

    return
  end function Accretion_Disk_Radiative_Efficiency

  double precision function Accretion_Disk_Jet_Power(thisNode,massAccretionRate)
    !% Computes the jet power for an accretion disk in units of $M_\odot$ (km/s)$^2$ Gyr$^{-1}$.
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    ! Ensure the module is initalized.
    call Accretion_Disks_Initialize

    ! Get the radiative efficiency.
    Accretion_Disk_Jet_Power=Accretion_Disk_Jet_Power_Get(thisNode,massAccretionRate)

    return
  end function Accretion_Disk_Jet_Power

  double precision function Black_Hole_Spin_Up_Rate(thisNode,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisNode} due to accretion from an accretion
    !% disk.
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    ! Ensure the module is initalized.
    call Accretion_Disks_Initialize

    ! Get the spin up rate.
    Black_Hole_Spin_Up_Rate=Black_Hole_Spin_Up_Rate_Get(thisNode,massAccretionRate)

    return
  end function Black_Hole_Spin_Up_Rate

end module Accretion_Disks
