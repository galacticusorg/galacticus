!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations related to accretion disks.

module Accretion_Disks
  !% Implements calculations related to accretion disks.
  use ISO_Varying_String
  implicit none
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

    if (.not.accretionDisksInitialized) then
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
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('accretionDisksMethod',accretionDisksMethod,defaultValue='switched')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="accretionDisksMethod" type="functionCall" functionType="void">
          !#  <functionArgs>accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get</functionArgs>
          include 'accretion_disks.inc'
          !# </include>
          if (.not.(associated(Accretion_Disk_Radiative_Efficiency_Get).and.associated(Black_Hole_Spin_Up_Rate_Get) &
               & .and.associated(Accretion_Disk_Jet_Power_Get))) call&
               & Galacticus_Error_Report('Accretion_Disks_Initialize','method ' //char(accretionDisksMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          accretionDisksInitialized=.true.
       end if
       !$omp end critical(accretionDisksInitialize)
    end if
    return
  end subroutine Accretion_Disks_Initialize
  
  double precision function Accretion_Disk_Radiative_Efficiency(thisBlackHole,massAccretionRate)
    !% Computes the radiative efficiency for an accretion disk.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate

    ! Ensure the module is initalized.
    call Accretion_Disks_Initialize()

    ! Get the radiative efficiency.
    Accretion_Disk_Radiative_Efficiency=Accretion_Disk_Radiative_Efficiency_Get(thisBlackHole,massAccretionRate)

    return
  end function Accretion_Disk_Radiative_Efficiency

  double precision function Accretion_Disk_Jet_Power(thisBlackHole,massAccretionRate)
    !% Computes the jet power for an accretion disk in units of $M_\odot$ (km/s)$^2$ Gyr$^{-1}$.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate

    ! Ensure the module is initalized.
    call Accretion_Disks_Initialize

    ! Get the radiative efficiency.
    Accretion_Disk_Jet_Power=Accretion_Disk_Jet_Power_Get(thisBlackHole,massAccretionRate)

    return
  end function Accretion_Disk_Jet_Power

  double precision function Black_Hole_Spin_Up_Rate(thisBlackHole,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisNode} due to accretion from an accretion
    !% disk.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate

    ! Ensure the module is initalized.
    call Accretion_Disks_Initialize

    ! Get the spin up rate.
    Black_Hole_Spin_Up_Rate=Black_Hole_Spin_Up_Rate_Get(thisBlackHole,massAccretionRate)

    return
  end function Black_Hole_Spin_Up_Rate

end module Accretion_Disks
