!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
!!    Andrew Benson <abenson@carnegiescience.edu>
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
  public :: Star_Formation_Rate_Surface_Density_Disk          , Star_Formation_Rate_Surface_Density_Disk_Intervals, &
       &    Star_Formation_Rate_Surface_Density_Disk_Unchanged

  ! Flag to indicate if this module has been initialized.
  logical                                                                :: moduleInitialized                                     =.false.

  ! Name of method to use.
  type     (varying_string                                    )          :: starFormationRateSurfaceDensityDisksMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Rate_Surface_Density_Disk          ), pointer :: Star_Formation_Rate_Surface_Density_Disk_Get          =>null()
  procedure(Star_Formation_Rate_Surface_Density_Disk_Intervals), pointer :: Star_Formation_Rate_Surface_Density_Disk_Intervals_Get=>null()
  procedure(Star_Formation_Rate_Surface_Density_Disk_Unchanged), pointer :: Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get=>null()

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
          !# <inputParameter>
          !#   <name>starFormationRateSurfaceDensityDisksMethod</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>var_str('KMT09')</defaultValue>
          !#   <description>The name of the method to be used for computing star formation timescales in disks.</description>
          !#   <group>starFormation</group>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !# </inputParameter>
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="starFormationRateSurfaceDensityDisksMethod" type="functionCall" functionType="void">
          !#  <functionArgs>starFormationRateSurfaceDensityDisksMethod,Star_Formation_Rate_Surface_Density_Disk_Get,Star_Formation_Rate_Surface_Density_Disk_Intervals_Get,Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get</functionArgs>
          include 'star_formation.rate_surface_density.disks.inc'
          !# </include>
          if (.not.(associated(Star_Formation_Rate_Surface_Density_Disk_Get).and.associated(Star_Formation_Rate_Surface_Density_Disk_Intervals_Get).and.associated(Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get))) &
               & call Galacticus_Error_Report('method '//char(starFormationRateSurfaceDensityDisksMethod)//' is unrecognized'//{introspection:location})
          moduleInitialized=.true.
       end if
       !$omp end critical(Star_Formation_Rate_Surface_Density_Disks_Initialization)
    end if
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_Initialize

  function Star_Formation_Rate_Surface_Density_Disk_Intervals(thisNode,radiusInner,radiusOuter)
    !% Return a set of integration intervals to use when integrating over the surface density of star formation rate.
    implicit none
    double precision          , allocatable  , dimension(:,:) :: Star_Formation_Rate_Surface_Density_Disk_Intervals
    type            (treeNode), intent(inout), target         :: thisNode
    double precision          , intent(in   )                 :: radiusInner, radiusOuter

    ! Initialize the module.
    call Star_Formation_Rate_Surface_Density_Disks_Initialize()
    ! Get the star formation rate surface density.
    Star_Formation_Rate_Surface_Density_Disk_Intervals=Star_Formation_Rate_Surface_Density_Disk_Intervals_Get(thisNode,radiusInner,radiusOuter)
    return
  end function Star_Formation_Rate_Surface_Density_Disk_Intervals

  logical function Star_Formation_Rate_Surface_Density_Disk_Unchanged(thisNode)
    !% Return true if the surface density rate of star formation is unchanged since the previous evaluation.
    implicit none
    type(treeNode), intent(inout) :: thisNode

    ! Initialize the module.
    call Star_Formation_Rate_Surface_Density_Disks_Initialize()
    ! Get the star formation rate surface density.
    Star_Formation_Rate_Surface_Density_Disk_Unchanged=Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get(thisNode)
    return
  end function Star_Formation_Rate_Surface_Density_Disk_Unchanged

  double precision function Star_Formation_Rate_Surface_Density_Disk(thisNode,radius)
    !% Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) in the disk component of
    !% {\normalfont \ttfamily thisNode} at the given {\normalfont \ttfamily radius}.
    implicit none
    type            (treeNode), intent(inout) :: thisNode
    double precision          , intent(in   ) :: radius

    ! Initialize the module.
    call Star_Formation_Rate_Surface_Density_Disks_Initialize()
    ! Get the star formation rate surface density.
    Star_Formation_Rate_Surface_Density_Disk=Star_Formation_Rate_Surface_Density_Disk_Get(thisNode,radius)
    return
  end function Star_Formation_Rate_Surface_Density_Disk

end module Star_Formation_Rate_Surface_Density_Disks
