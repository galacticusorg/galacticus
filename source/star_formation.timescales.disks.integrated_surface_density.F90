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

!% Contains a module which implements a global star formation timescale
!% for galactic disks by integrating over the surface density of star
!% formation rate.

module Star_Formation_Timescale_Disks_Integrated_SD
  !% Implements the a global star formation timescale for galactic disks by integrating over the surface density of star
  !% formation rate.
  use Galacticus_Nodes
  implicit none
  private
  public :: Star_Formation_Timescale_Disks_Integrated_SD_Initialize

  ! Pointer to active node used in integral functions, plus variables needed by integral function.
  type(treeNode),   pointer :: activeNode
  !$omp threadprivate(activeNode)

contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Integrated_SD_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Integrated_SD_Initialize(starFormationTimescaleDisksMethod&
       &,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``integrated surface density'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: starFormationTimescaleDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Timescale_Disk_Get
    
    if (starFormationTimescaleDisksMethod == 'integratedSurfaceDensity') Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Integrated_SD
    return
  end subroutine Star_Formation_Timescale_Disks_Integrated_SD_Initialize

  double precision function Star_Formation_Timescale_Disk_Integrated_SD(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the
    !% galactic disk of {\tt thisNode}, by integrating over the surface
    !% density of star formation rate.
    use Galacticus_Nodes
    use Numerical_Constants_Math
    use FGSL
    use Numerical_Integration
    use, intrinsic :: ISO_C_Binding
    implicit none
    type (treeNode                  ), intent(inout), pointer :: thisNode
    class(nodeComponentDisk         ),                pointer :: thisDiskComponent
    double precision                 , parameter              :: radiusInnerDimensionless=0.0d0,radiusOuterDimensionless=10.0d0
    double precision                                          :: gasMass,diskScaleRadius,starFormationRate,radiusInner,radiusOuter
    type (c_ptr                     )                         :: parameterPointer
    type (fgsl_function             )                         :: integrandFunction
    type (fgsl_integration_workspace)                         :: integrationWorkspace

    ! Get the disk properties.
    thisDiskComponent => thisNode%disk()
    gasMass             =thisDiskComponent%massGas()
    diskScaleRadius     =thisDiskComponent%radius ()
    ! Check if the disk is physical.
    if (gasMass <= 0.0d0 .or. diskScaleRadius <= 0.0d0) then
       ! It is not, so return zero timescale.
       Star_Formation_Timescale_Disk_Integrated_SD=0.0d0
    else
       ! Set a pointer to the node that is accessible by integral function.
       activeNode => thisNode
       ! Compute suitable limits for the integration.
       radiusInner=diskScaleRadius*radiusInnerDimensionless
       radiusOuter=diskScaleRadius*radiusOuterDimensionless
       ! Compute the star formation rate. A low order integration rule (FGSL_Integ_Gauss15) works well here.
       starFormationRate=2.0d0*Pi*Integrate(radiusInner,radiusOuter,Star_Formation_Rate_Integrand_Surface_Density&
            &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3&
            &,integrationRule=FGSL_Integ_Gauss15)
       call Integrate_Done(integrandFunction,integrationWorkspace)
       ! Compute the star formation timescale.
       if (starFormationRate > 0.0d0) then
          Star_Formation_Timescale_Disk_Integrated_SD=gasMass/starFormationRate
       else
          Star_Formation_Timescale_Disk_Integrated_SD=0.0d0
       end if
    end if
    return
  end function Star_Formation_Timescale_Disk_Integrated_SD
  
  function Star_Formation_Rate_Integrand_Surface_Density(radius,parameterPointer) bind(c)
    !% Integrand function for the ``integrated surface density'' star formation rate calculation.
    use, intrinsic :: ISO_C_Binding
    use Star_Formation_Rate_Surface_Density_Disks
    implicit none
    real(c_double)        :: Star_Formation_Rate_Integrand_Surface_Density
    real(c_double), value :: radius
    type(c_ptr   ), value :: parameterPointer

    ! Compute the star formation rate integrand.
    Star_Formation_Rate_Integrand_Surface_Density=radius*Star_Formation_Rate_Surface_Density_Disk(activeNode,radius)
    return
  end function Star_Formation_Rate_Integrand_Surface_Density

end module Star_Formation_Timescale_Disks_Integrated_SD
