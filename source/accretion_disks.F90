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
       !@   <defaultValue>ADAF</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Selects which accretion disk method should be used.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('accretionDisksMethod',accretionDisksMethod,defaultValue='ADAF')
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
    !$omp end critical(accretionDisksInitialize)

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
