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

!% Contains a module which implements structure tasks for the exponential disk node component.

module Node_Component_Disk_Exponential_Structure_Tasks
  use Kind_Numbers
  implicit none
  private
  public :: Node_Component_Disk_Exponential_Calculation_Reset, Node_Component_Disk_Exponential_Surface_Density

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
  !#   <unitName>Node_Component_Disk_Exponential_Calculation_Reset</unitName>
  !# </calculationResetTask>
  subroutine Node_Component_Disk_Exponential_Calculation_Reset(thisNode)
    !% Reset exponential disk structure calculations.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    radiusScaleDiskComputed             =.false.
    surfaceDensityCentralGasComputed    =.false.
    surfaceDensityCentralStellarComputed=.false.
    surfaceDensityCentralTotalComputed  =.false.
    lastUniqueID                        =thisNode%uniqueID()
    return
  end subroutine Node_Component_Disk_Exponential_Calculation_Reset

  !# <surfaceDensityTask>
  !#  <unitName>Node_Component_Disk_Exponential_Surface_Density</unitName>
  !# </surfaceDensityTask>
  double precision function Node_Component_Disk_Exponential_Surface_Density(thisNode,positionCylindrical,massType,componentType,haloLoaded)
    !% Computes the surface density at a given position for an exponential disk.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    use Galacticus_Nodes
    implicit none
    type (treeNode         ), intent(inout), pointer  :: thisNode
    class(nodeComponentDisk),                pointer  :: thisDiskComponent
    integer                 , intent(in   )           :: massType,componentType
    double precision        , intent(in   )           :: positionCylindrical(3)
    logical                 , intent(in   ), optional :: haloLoaded
    double precision                                  :: fractionalRadius
    
    ! Return immediately if disk component is not requested.    
    Node_Component_Disk_Exponential_Surface_Density=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return

    ! Get the disk component and check that it is of the exponential class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)
    
       ! Check whether this is a new node.
       if (thisNode%uniqueID() /= lastUniqueID) call Node_Component_Disk_Exponential_Calculation_Reset(thisNode)

       ! Determine disk radius.
       if (.not.radiusScaleDiskComputed) then
          radiusScaleDisk        =thisDiskComponent%radius()
          radiusScaleDiskComputed=.true.
       end if

       ! Determine mass type.
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic)
          if (.not.surfaceDensityCentralTotalComputed  ) then
             surfaceDensityCentralTotal          =(thisDiskComponent%massGas()+thisDiskComponent%massStellar())/2.0d0/Pi/radiusScaleDisk**2
             surfaceDensityCentralTotalComputed  =.true.
          end if
          Node_Component_Disk_Exponential_Surface_Density=surfaceDensityCentralTotal
       case (massTypeGaseous)
          if (.not.surfaceDensityCentralGasComputed    ) then
             surfaceDensityCentralGas            = thisDiskComponent%massGas()                                 /2.0d0/Pi/radiusScaleDisk**2
             surfaceDensityCentralGasComputed    =.true.
          end if
          Node_Component_Disk_Exponential_Surface_Density=surfaceDensityCentralGas
       case (massTypeStellar)
          if (.not.surfaceDensityCentralStellarComputed) then
             surfaceDensityCentralStellar        =                             thisDiskComponent%massStellar() /2.0d0/Pi/radiusScaleDisk**2
             surfaceDensityCentralStellarComputed=.true.
          end if
          Node_Component_Disk_Exponential_Surface_Density=surfaceDensityCentralStellar
       end select

       ! Return if no density.
       if (Node_Component_Disk_Exponential_Surface_Density <= 0.0d0) return
       
       ! Compute the actual density.
       fractionalRadius=positionCylindrical(1)/radiusScaleDisk
       Node_Component_Disk_Exponential_Surface_Density=Node_Component_Disk_Exponential_Surface_Density*exp(-fractionalRadius)
    end select
    return
  end function Node_Component_Disk_Exponential_Surface_Density
  
end module Node_Component_Disk_Exponential_Structure_Tasks
