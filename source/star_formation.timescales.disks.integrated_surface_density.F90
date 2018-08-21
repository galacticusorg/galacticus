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
  type            (treeNode), pointer :: activeNode
  !$omp threadprivate(activeNode)

  ! Stored previous rate.
  double precision                    :: starFormationRatePrevious                              =-1.0d0
  !$omp threadprivate(starFormationRatePrevious)
  
  ! Tolerance parameter.
  double precision                    :: starFormationTimescaleIntegratedSurfaceDensityTolerance
  
contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Integrated_SD_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Integrated_SD_Initialize(starFormationTimescaleDisksMethod&
       &,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``integrated surface density'' disk star formation timescale module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type     (varying_string                             ), intent(in   )          :: starFormationTimescaleDisksMethod
    procedure(Star_Formation_Timescale_Disk_Integrated_SD), intent(inout), pointer :: Star_Formation_Timescale_Disk_Get

    if (starFormationTimescaleDisksMethod == 'integratedSurfaceDensity') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Integrated_SD
       !# <inputParameter>
       !#   <name>starFormationTimescaleIntegratedSurfaceDensityTolerance</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d-3</defaultValue>
       !#   <description>Relative tolerance to use when integrating star formation rate surface densities over the disk.</description>
       !#   <source>globalParameters</source>
       !#   <type>float</type>
       !# </inputParameter>
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Integrated_SD_Initialize

  double precision function Star_Formation_Timescale_Disk_Integrated_SD(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the
    !% galactic disk of {\normalfont \ttfamily thisNode}, by integrating over the surface
    !% density of star formation rate.
    use Numerical_Constants_Math
    use Numerical_Integration
    use Star_Formation_Rate_Surface_Density_Disks
    implicit none
    type            (treeNode                                 ), intent(inout), target         :: thisNode
    double precision                                           , allocatable  , dimension(:,:) :: intervals
    class           (nodeComponentDisk                        ), pointer                       :: thisDiskComponent
    class           (starFormationRateSurfaceDensityDisksClass), pointer                       :: starFormationRateSurfaceDensityDisks_
    double precision                                           , parameter                     :: radiusInnerDimensionless             =0.0d0, radiusOuterDimensionless=10.0d0
    double precision                                                                           :: diskScaleRadius                            , gasMass                        , &
         &                                                                                        radiusInner                                , radiusOuter                    , &
         &                                                                                        starFormationRate
    type            (fgsl_function                            )                                :: integrandFunction
    type            (fgsl_integration_workspace               )                                :: integrationWorkspace
    logical                                                                                    :: integrationReset
    integer                                                                                    :: i

    ! Get the disk properties.
    starFormationRateSurfaceDensityDisks_ => starFormationRateSurfaceDensityDisks()
    thisDiskComponent => thisNode%disk()
    gasMass             =thisDiskComponent%massGas()
    diskScaleRadius     =thisDiskComponent%radius ()
    ! Test whether the star formation rate surface density function changed. If it did not we can re-use the previous integral.
    if (starFormationRateSurfaceDensityDisks_%unchanged(thisNode)) then
       starFormationRate=starFormationRatePrevious
    else
       ! Check if the disk is physical.
       if (gasMass <= 0.0d0 .or. diskScaleRadius <= 0.0d0) then
          ! It is not, so return zero timescale.
          starFormationRate=0.0d0
       else
          ! Set a pointer to the node that is accessible by integral function.
          activeNode => thisNode
          ! Compute suitable limits for the integration.
          radiusInner=diskScaleRadius*radiusInnerDimensionless
          radiusOuter=diskScaleRadius*radiusOuterDimensionless
          ! Get a set of intervals into which this integral should be broken.
          intervals=starFormationRateSurfaceDensityDisks_%intervals(thisNode,radiusInner,radiusOuter)
          ! Compute the star formation rate. A low order integration rule (FGSL_Integ_Gauss15) works well here.       
          starFormationRate=0.0d0
          do i=1,size(intervals,dim=2)       
             integrationReset=.true.
             starFormationRate=+starFormationRate                                                                                                &
                  &            +Integrate(                                                                                                       &
                  &                       intervals(1,i)                                                                                       , &
                  &                       intervals(2,i)                                                                                       , &
                  &                       Star_Formation_Rate_Integrand_Surface_Density                                                        , &
                  &                       integrandFunction                                                                                    , &
                  &                       integrationWorkspace                                                                                 , &
                  &                       reset                                        =integrationReset                                       , &
                  &                       toleranceAbsolute                            =0.0d+0                                                 , &
                  &                       toleranceRelative                            =starFormationTimescaleIntegratedSurfaceDensityTolerance, &
                  &                       integrationRule                              =FGSL_Integ_Gauss15                                       &
                  &                      )
             call Integrate_Done(integrandFunction,integrationWorkspace)
          end do
          starFormationRate=2.0d0*Pi*starFormationRate
       end if
       starFormationRatePrevious=starFormationRate
    end if
    ! Compute the star formation timescale.
    if (starFormationRate > 0.0d0) then
       Star_Formation_Timescale_Disk_Integrated_SD=gasMass/starFormationRate
    else
       Star_Formation_Timescale_Disk_Integrated_SD=0.0d0
    end if
    return
  end function Star_Formation_Timescale_Disk_Integrated_SD
  
  double precision function Star_Formation_Rate_Integrand_Surface_Density(radius)
    !% Integrand function for the ``integrated surface density'' star formation rate calculation.
    use Star_Formation_Rate_Surface_Density_Disks
    implicit none
    double precision                                           , intent(in   ) :: radius
    class           (starFormationRateSurfaceDensityDisksClass), pointer       :: starFormationRateSurfaceDensityDisks_

    ! Compute the star formation rate integrand.
    starFormationRateSurfaceDensityDisks_ => starFormationRateSurfaceDensityDisks()
    Star_Formation_Rate_Integrand_Surface_Density=radius*starFormationRateSurfaceDensityDisks_%rate(activeNode,radius)
    return
  end function Star_Formation_Rate_Integrand_Surface_Density
  
end module Star_Formation_Timescale_Disks_Integrated_SD
