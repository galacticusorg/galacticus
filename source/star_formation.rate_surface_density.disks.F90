!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of star formation rate surface densities for galactic disks.

module Star_Formation_Rate_Surface_Density_Disks
  !% Implements calculations of star formation rate surface densities for galactic disks.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Star_Formation_Rate_Surface_Density_Disk
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: moduleInitialized=.false.

  ! Name of method to use.
  type(varying_string) :: starFormationRateSurfaceDensityDisksMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Rate_Surface_Density_Disk), pointer :: Star_Formation_Rate_Surface_Density_Disk_Get => null()

contains

  subroutine Star_Formation_Rate_Surface_Density_Disks_Initialize
    !% Initialize the disk star formation rate surface density module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationRateSurfaceDensityDisksMethod" type="moduleUse">
    include 'star_formation.rate_surface_density.disks.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Star_Formation_Rate_Surface_Density_Disks_Initialization) 
       if (.not.moduleInitialized) then
          ! Get the disk star formation timescale method parameter.
          !@ <inputParameter>
          !@   <name>starFormationRateSurfaceDensityDisksMethod</name>
          !@   <defaultValue>KMT09</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing star formation timescales in disks.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('starFormationRateSurfaceDensityDisksMethod',starFormationRateSurfaceDensityDisksMethod,defaultValue='KMT09')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="starFormationRateSurfaceDensityDisksMethod" type="functionCall" functionType="void">
          !#  <functionArgs>starFormationRateSurfaceDensityDisksMethod,Star_Formation_Rate_Surface_Density_Disk_Get</functionArgs>
          include 'star_formation.rate_surface_density.disks.inc'
          !# </include>
          if (.not.associated(Star_Formation_Rate_Surface_Density_Disk_Get)) call&
               & Galacticus_Error_Report('Star_Formation_Rate_Surface_Density_Disks' ,'method '&
               &//char(starFormationRateSurfaceDensityDisksMethod)//' is unrecognized')
          moduleInitialized=.true.
       end if
       !$omp end critical(Star_Formation_Rate_Surface_Density_Disks_Initialization) 
    end if
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_Initialize

  double precision function Star_Formation_Rate_Surface_Density_Disk(thisNode,radius)
    !% Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) in the disk component of
    !% {\tt thisNode} at the given {\tt radius}.
    implicit none
    type            (treeNode), pointer, intent(inout) :: thisNode
    double precision          ,          intent(in   ) :: radius

    ! Initialize the module.
    call Star_Formation_Rate_Surface_Density_Disks_Initialize()
    ! Get the star formation rate surface density.
    Star_Formation_Rate_Surface_Density_Disk=Star_Formation_Rate_Surface_Density_Disk_Get(thisNode,radius)
    return
  end function Star_Formation_Rate_Surface_Density_Disk
  
end module Star_Formation_Rate_Surface_Density_Disks
