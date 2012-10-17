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

!% Contains a module which implements structure tasks for exponential disk components.

module Tree_Node_Methods_Exponential_Disk_Structure_Tasks
  !% Implement structure tasks exponential disk tree node methods.
  use Tree_Node_Methods_Disk_Exponential_Data
  use Kind_Numbers
  implicit none
  private
  public :: Exponential_Disk_Density, Exponential_Disk_Surface_Density, Exponential_Disk_Structure_Tasks_Reset
  
  ! Record of unique ID of node which we last computed results for.
  integer(kind=kind_int8) :: lastUniqueID=-1
  !$omp threadprivate(lastUniqueID)

  ! Records of previously computed and stored quantities.
  logical                 :: surfaceDensityCentralGasComputed,surfaceDensityCentralStellarComputed,surfaceDensityCentralTotalComputed
  !$omp threadprivate(surfaceDensityCentralGasComputed,surfaceDensityCentralStellarComputed,surfaceDensityCentralTotalComputed)
  double precision        :: surfaceDensityCentralGas,surfaceDensityCentralStellar,surfaceDensityCentralTotal
  !$omp threadprivate(surfaceDensityCentralGas,surfaceDensityCentralStellar,surfaceDensityCentralTotal)
  logical                 :: radiusScaleDiskComputed
  !$omp threadprivate(radiusScaleDiskComputed)
  double precision        :: radiusScaleDisk
  !$omp threadprivate(radiusScaleDisk)

contains

  !# <calculationResetTask>
  !#   <unitName>Exponential_Disk_Structure_Tasks_Reset</unitName>
  !# </calculationResetTask>
  subroutine Exponential_Disk_Structure_Tasks_Reset(thisNode)
    !% Reset exponential disk structure calculations.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    radiusScaleDiskComputed             =.false.
    surfaceDensityCentralGasComputed    =.false.
    surfaceDensityCentralStellarComputed=.false.
    surfaceDensityCentralTotalComputed  =.false.
    lastUniqueID                        =thisNode%uniqueID()
    return
  end subroutine Exponential_Disk_Structure_Tasks_Reset

  !# <densityTask>
  !#  <unitName>Exponential_Disk_Density</unitName>
  !# </densityTask>
  subroutine Exponential_Disk_Density(thisNode,positionSpherical,massType,componentType,componentDensity)
    !% Computes the density at a given position for an exponential disk.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    use Coordinate_Systems
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: positionSpherical(3)
    double precision, intent(out)            :: componentDensity
    double precision, parameter              :: diskHeightToRadiusRatio=0.1d0
    double precision                         :: fractionalRadius,fractionalHeight,positionCylindrical(3)
    
    componentDensity=0.0d0
    if (.not.methodSelected                           ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
    if (.not.thisNode%componentExists(componentIndex)) return

    ! Determine mass type.
    select case (massType)
    case (massTypeAll,massTypeBaryonic,massTypeGalactic)
       componentDensity=Tree_Node_Disk_Gas_Mass(thisNode)+Tree_Node_Disk_Stellar_Mass(thisNode)
    case (massTypeGaseous)
       componentDensity=Tree_Node_Disk_Gas_Mass(thisNode)
    case (massTypeStellar)
       componentDensity=Tree_Node_Disk_Stellar_Mass(thisNode)
    end select
    ! Return if no density.
    if (componentDensity                <= 0.0d0) return
    ! Return if zero size disk.
    if (Tree_Node_Disk_Radius(thisNode) <= 0.0d0) return
    ! Compute the actual density.
    positionCylindrical=Coordinates_Spherical_To_Cylindrical(positionSpherical)
    fractionalRadius=positionCylindrical(1)/Tree_Node_Disk_Radius(thisNode)
    fractionalHeight=positionCylindrical(2)/(diskHeightToRadiusRatio*Tree_Node_Disk_Radius(thisNode))
    componentDensity=componentDensity*dexp(-fractionalRadius)/cosh(0.5d0*fractionalHeight)**2/4.0d0/Pi&
         &/Tree_Node_Disk_Radius(thisNode)**3/diskHeightToRadiusRatio
    return
  end subroutine Exponential_Disk_Density

  !# <surfaceDensityTask>
  !#  <unitName>Exponential_Disk_Surface_Density</unitName>
  !# </surfaceDensityTask>
  subroutine Exponential_Disk_Surface_Density(thisNode,positionCylindrical,massType,componentType,componentDensity)
    !% Computes the surface density at a given position for an exponential disk.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: positionCylindrical(3)
    double precision, intent(out)            :: componentDensity
    double precision                         :: fractionalRadius
    
    componentDensity=0.0d0
    if (.not.methodSelected                                                             ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
    if (.not.thisNode%componentExists(componentIndex)                                   ) return

    ! Check whether this is a new node.
    if (thisNode%uniqueID() /= lastUniqueID) call Exponential_Disk_Structure_Tasks_Reset(thisNode)

    ! Determine disk radius.
    if (.not.radiusScaleDiskComputed) then
       radiusScaleDisk        =Tree_Node_Disk_Radius(thisNode)
       radiusScaleDiskComputed=.true.
    end if

    ! Determine mass type.
    select case (massType)
    case (massTypeAll,massTypeBaryonic,massTypeGalactic)
       if (.not.surfaceDensityCentralTotalComputed  ) then
          surfaceDensityCentralTotal          =(Tree_Node_Disk_Gas_Mass(thisNode)+Tree_Node_Disk_Stellar_Mass(thisNode))/2.0d0/Pi&
               & /radiusScaleDisk**2
          surfaceDensityCentralTotalComputed  =.true.
       end if
       componentDensity=surfaceDensityCentralTotal
    case (massTypeGaseous)
       if (.not.surfaceDensityCentralGasComputed    ) then
          surfaceDensityCentralGas            = Tree_Node_Disk_Gas_Mass(thisNode)                                       /2.0d0/Pi&
               & /radiusScaleDisk**2
          surfaceDensityCentralGasComputed    =.true.
       end if
       componentDensity=surfaceDensityCentralGas
    case (massTypeStellar)
       if (.not.surfaceDensityCentralStellarComputed) then
          surfaceDensityCentralStellar        =                                   Tree_Node_Disk_Stellar_Mass(thisNode)/2.0d0/Pi&
               & /radiusScaleDisk**2
          surfaceDensityCentralStellarComputed=.true.
       end if
       componentDensity=surfaceDensityCentralStellar
    end select

    ! Return if no density.
    if (componentDensity <= 0.0d0) return

    ! Compute the actual density.
    fractionalRadius=positionCylindrical(1)/radiusScaleDisk
    componentDensity=componentDensity*dexp(-fractionalRadius)
    return
  end subroutine Exponential_Disk_Surface_Density

end module Tree_Node_Methods_Exponential_Disk_Structure_Tasks
