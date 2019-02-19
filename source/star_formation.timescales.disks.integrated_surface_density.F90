!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implementation of a timescale for star formation in galactic disks which computes the timescale by integrating a
  !% star formation rate over the disk.

  use Star_Formation_Rate_Surface_Density_Disks

  !# <starFormationTimescaleDisks name="starFormationTimescaleDisksIntgrtdSurfaceDensity">
  !#  <description>A timescale for star formation in galactic disks which computes the timescale by integrating a star formation rate over the disk.</description>
  !# </starFormationTimescaleDisks>
  type, extends(starFormationTimescaleDisksClass) :: starFormationTimescaleDisksIntgrtdSurfaceDensity
     !% Implementation of a timescale for star formation in galactic disks which computes the timescale by integrating a
     !% star formation rate over the disk.
     private
     class           (starFormationRateSurfaceDensityDisksClass), pointer :: starFormationRateSurfaceDensityDisks_ => null()
     double precision                                                     :: tolerance                            , starFormationRatePrevious
   contains
     final     ::                     intgrtdSurfaceDensityDestructor
     procedure :: timescale        => intgrtdSurfaceDensityTimescale
  end type starFormationTimescaleDisksIntgrtdSurfaceDensity

  interface starFormationTimescaleDisksIntgrtdSurfaceDensity
     !% Constructors for the {\normalfont \ttfamily intgrtdSurfaceDensity} timescale for star formation in disks class.
     module procedure intgrtdSurfaceDensityConstructorParameters
     module procedure intgrtdSurfaceDensityConstructorInternal
  end interface starFormationTimescaleDisksIntgrtdSurfaceDensity

  ! Pointer to active node used in integral functions, plus variables needed by integral function.
  class(starFormationTimescaleDisksIntgrtdSurfaceDensity), pointer :: intgrtdSurfaceDensitySelf
  type (treeNode                                        ), pointer :: intgrtdSurfaceDensityNode
  !$omp threadprivate(intgrtdSurfaceDensitySelf,intgrtdSurfaceDensityNode)

contains

  function intgrtdSurfaceDensityConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily intgrtdSurfaceDensity} timescale for star formation in disks class which takes a
    !% parameter set as input.
    use Input_Parameters
    implicit none
    type            (starFormationTimescaleDisksIntgrtdSurfaceDensity)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (starFormationRateSurfaceDensityDisksClass       ), pointer       :: starFormationRateSurfaceDensityDisks_
    double precision                                                                  :: tolerance

    !# <inputParameter>
    !#   <name>tolerance</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d-3</defaultValue>
    !#   <description>Relative tolerance to use when integrating star formation rate surface densities over the disk.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <objectBuilder class="starFormationRateSurfaceDensityDisks" name="starFormationRateSurfaceDensityDisks_" source="parameters"/>
    self=starFormationTimescaleDisksIntgrtdSurfaceDensity(tolerance,starFormationRateSurfaceDensityDisks_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="starFormationRateSurfaceDensityDisks_"/>
    return
  end function intgrtdSurfaceDensityConstructorParameters

  function intgrtdSurfaceDensityConstructorInternal(tolerance,starFormationRateSurfaceDensityDisks_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily intgrtdSurfaceDensity} timescale for star formation in disks class.
    implicit none
    type            (starFormationTimescaleDisksIntgrtdSurfaceDensity)                        :: self
    double precision                                                  , intent(in   )         :: tolerance
    class           (starFormationRateSurfaceDensityDisksClass       ), intent(in   ), target :: starFormationRateSurfaceDensityDisks_
    !# <constructorAssign variables="tolerance, *starFormationRateSurfaceDensityDisks_"/>
    
    return
  end function intgrtdSurfaceDensityConstructorInternal

  subroutine intgrtdSurfaceDensityDestructor(self)
    !% Destructor for the {\normalfont \ttfamily intgrtdSurfaceDensity} timescale for star formation in disks class.
    implicit none
    type(starFormationTimescaleDisksIntgrtdSurfaceDensity), intent(inout) :: self

    !# <objectDestructor name="self%starFormationRateSurfaceDensityDisks_" />
    return
  end subroutine intgrtdSurfaceDensityDestructor

  double precision function intgrtdSurfaceDensityTimescale(self,node)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\normalfont \ttfamily node}, by integrating
    !% over the surface density of star formation rate.
    use FGSL                    , only : FGSL_Integ_Gauss61, fgsl_function, fgsl_integration_workspace
    use Galacticus_Nodes        , only : nodeComponentDisk
    use Numerical_Constants_Math
    use Numerical_Integration
    implicit none
    class           (starFormationTimescaleDisksIntgrtdSurfaceDensity), intent(inout), target         :: self
    type            (treeNode                                        ), intent(inout), target         :: node
    double precision                                                  , allocatable  , dimension(:,:) :: intervals
    class           (nodeComponentDisk                               ), pointer                       :: disk
    double precision                                                  , parameter                     :: radiusInnerDimensionless=0.0d0, radiusOuterDimensionless=10.0d0
    double precision                                                                                  :: radiusDisk                    , massGas                        , &
         &                                                                                               radiusInner                   , radiusOuter                    , &
         &                                                                                               starFormationRate
    type            (fgsl_function                                   )                                :: integrandFunction
    type            (fgsl_integration_workspace                      )                                :: integrationWorkspace
    logical                                                                                           :: integrationReset
    integer                                                                                           :: i

    ! Get the disk properties.
    disk       => node%disk   ()
    massGas    =  disk%massGas()
    radiusDisk =  disk%radius ()
    ! Test whether the star formation rate surface density function changed. If it did not we can re-use the previous integral.
    if (self%starFormationRateSurfaceDensityDisks_%unchanged(node)) then
       starFormationRate=self%starFormationRatePrevious
    else
       ! Check if the disk is physical.
       if (massGas <= 0.0d0 .or. radiusDisk <= 0.0d0) then
          ! It is not, so return zero timescale.
          starFormationRate=0.0d0
       else
          ! Set a pointer to self and to the node that is accessible by integral function.
          intgrtdSurfaceDensitySelf => self
          intgrtdSurfaceDensityNode => node
          ! Compute suitable limits for the integration.
          radiusInner=radiusDisk*radiusInnerDimensionless
          radiusOuter=radiusDisk*radiusOuterDimensionless
          ! Get a set of intervals into which this integral should be broken.
          intervals=self%starFormationRateSurfaceDensityDisks_%intervals(node,radiusInner,radiusOuter)
          ! Compute the star formation rate. A low order integration rule (FGSL_Integ_Gauss15) works well here.       
          starFormationRate=0.0d0
          do i=1,size(intervals,dim=2)       
             integrationReset=.true.
             starFormationRate=+starFormationRate                                               &
                  &            +Integrate(                                                      &
                  &                       intervals(1,i)                                      , &
                  &                       intervals(2,i)                                      , &
                  &                       intgrtdSurfaceDensityIntegrand                      , &
                  &                       integrandFunction                                   , &
                  &                       integrationWorkspace                                , &
                  &                       reset                            =integrationReset  , &
                  &                       toleranceAbsolute                =0.0d+0            , &
                  &                       toleranceRelative                =self%tolerance    , &
                  &                       integrationRule                  =FGSL_Integ_Gauss61  &
                  &                      )
             call Integrate_Done(integrandFunction,integrationWorkspace)
          end do
          starFormationRate=2.0d0*Pi*starFormationRate
       end if
       self%starFormationRatePrevious=starFormationRate
    end if
    ! Compute the star formation timescale.
    if (starFormationRate > 0.0d0) then
       intgrtdSurfaceDensityTimescale=massGas/starFormationRate
    else
       intgrtdSurfaceDensityTimescale=0.0d0
    end if
    return
  end function intgrtdSurfaceDensityTimescale
  
  double precision function intgrtdSurfaceDensityIntegrand(radius)
    !% Integrand function for the ``integrated surface density'' star formation rate calculation.
    implicit none
    double precision, intent(in   ) :: radius

    intgrtdSurfaceDensityIntegrand=radius*intgrtdSurfaceDensitySelf%starFormationRateSurfaceDensityDisks_%rate(intgrtdSurfaceDensityNode,radius)
    return
  end function intgrtdSurfaceDensityIntegrand
