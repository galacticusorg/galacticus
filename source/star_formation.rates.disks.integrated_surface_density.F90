!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  Implementation of a star formation rate in galactic disks which computes the rate by integrating a star formation rate over
  the disk.
  !!}

  use :: Numerical_Integration                    , only : integrator
  use :: Star_Formation_Rate_Surface_Density_Disks, only : starFormationRateSurfaceDensityDisksClass

  !![
  <starFormationRateDisks name="starFormationRateDisksIntgrtdSurfaceDensity">
   <description>
    A star formation rate in galactic disks which computes the rate by integrating a star formation rate over the
    disk. Specifically, the star formation rate is given by
    \begin{equation}
     \dot{M}_\star = \int_0^\infty 2 \pi r \dot{\Sigma}_\star(r) \mathrm{d}r,
    \end{equation}
    where $\dot{\Sigma}_\star(r)$ is the surface density of star formation rate.
   </description>
  </starFormationRateDisks>
  !!]
  type, extends(starFormationRateDisksClass) :: starFormationRateDisksIntgrtdSurfaceDensity
     !!{
     Implementation of a rate for star formation in galactic disks which computes the rate by integrating a star formation rate
     over the disk.
     !!}
     private
     class           (starFormationRateSurfaceDensityDisksClass), pointer :: starFormationRateSurfaceDensityDisks_ => null()
     type            (integrator                               )          :: integrator_
     double precision                                                     :: tolerance                                      , starFormationRatePrevious
   contains
     final     ::         intgrtdSurfaceDensityDestructor
     procedure :: rate => intgrtdSurfaceDensityRate
  end type starFormationRateDisksIntgrtdSurfaceDensity

  interface starFormationRateDisksIntgrtdSurfaceDensity
     !!{
     Constructors for the \refClass{starFormationRateDisksIntgrtdSurfaceDensity} star formation rate in disks class.
     !!}
     module procedure intgrtdSurfaceDensityConstructorParameters
     module procedure intgrtdSurfaceDensityConstructorInternal
  end interface starFormationRateDisksIntgrtdSurfaceDensity

  ! Pointer to active node used in integral functions, plus variables needed by integral function.
  class(starFormationRateDisksIntgrtdSurfaceDensity), pointer :: self_
  type (treeNode                                   ), pointer :: node_
  !$omp threadprivate(self_,node_)
    
contains

  function intgrtdSurfaceDensityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationRateDisksIntgrtdSurfaceDensity} star formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationRateDisksIntgrtdSurfaceDensity)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (starFormationRateSurfaceDensityDisksClass  ), pointer       :: starFormationRateSurfaceDensityDisks_
    double precision                                                             :: tolerance

    !![
    <inputParameter>
      <name>tolerance</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>Relative tolerance to use when integrating star formation rate surface densities over the disk.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateSurfaceDensityDisks" name="starFormationRateSurfaceDensityDisks_" source="parameters"/>
    !!]
    self=starFormationRateDisksIntgrtdSurfaceDensity(tolerance,starFormationRateSurfaceDensityDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateSurfaceDensityDisks_"/>
    !!]
    return
  end function intgrtdSurfaceDensityConstructorParameters

  function intgrtdSurfaceDensityConstructorInternal(tolerance,starFormationRateSurfaceDensityDisks_) result(self)
    !!{
    Internal constructor for the \refClass{starFormationRateDisksIntgrtdSurfaceDensity} star formation rate in disks class.
    !!}
    use :: Numerical_Integration, only : GSL_Integ_Gauss15
    implicit none
    type            (starFormationRateDisksIntgrtdSurfaceDensity)                        :: self
    double precision                                             , intent(in   )         :: tolerance
    class           (starFormationRateSurfaceDensityDisksClass  ), intent(in   ), target :: starFormationRateSurfaceDensityDisks_
    ! Set an absolute tolerance for star formation rate integration. The value chosen is such that below this tolerance much less
    ! than a single star would form over the age of the Universe.
    double precision                                             , parameter             :: toleranceAbsolute                    =1.0d-12
    !![
    <constructorAssign variables="tolerance, *starFormationRateSurfaceDensityDisks_"/>
    !!]

    self%integrator_=integrator(intgrtdSurfaceDensityIntegrand,toleranceAbsolute=toleranceAbsolute,toleranceRelative=self%tolerance,integrationRule=GSL_Integ_Gauss15)
    return
  end function intgrtdSurfaceDensityConstructorInternal

  subroutine intgrtdSurfaceDensityDestructor(self)
    !!{
    Destructor for the \refClass{starFormationRateDisksIntgrtdSurfaceDensity} star formation rate in disks class.
    !!}
    implicit none
    type(starFormationRateDisksIntgrtdSurfaceDensity), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateSurfaceDensityDisks_" />
    !!]
    return
  end subroutine intgrtdSurfaceDensityDestructor

  double precision function intgrtdSurfaceDensityRate(self,node)
    !!{
    Returns the star formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily node}, by
    integrating over the surface density of star formation rate.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentDisk, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (starFormationRateDisksIntgrtdSurfaceDensity), intent(inout), target         :: self
    type            (treeNode                                   ), intent(inout), target         :: node
    double precision                                             , allocatable  , dimension(  :) :: integralsAnalytic
    logical                                                      , allocatable  , dimension(  :) :: intervalIsAnalytic
    double precision                                             , allocatable  , dimension(:,:) :: intervals
    class           (nodeComponentDisk                          ), pointer                       :: disk
    double precision                                             , parameter                     :: radiusInnerDimensionless=0.0d+00, radiusOuterDimensionless=10.0d0
    double precision                                                                             :: radiusDisk                      , massGas                        , &
         &                                                                                          radiusInner                     , radiusOuter
    integer                                                                                      :: i

    ! Test whether the star formation rate surface density function changed. If it did not we can re-use the previous integral.
    if (self%starFormationRateSurfaceDensityDisks_%unchanged(node)) then
       intgrtdSurfaceDensityRate=self%starFormationRatePrevious
    else
       ! Get the disk properties.
       disk       => node%disk   ()
       massGas    =  disk%massGas()
       radiusDisk =  disk%radius ()
       ! Check if the disk is physical.
       if (massGas <= 0.0d0 .or. radiusDisk <= 0.0d0) then
          ! It is not, so return zero rate.
          intgrtdSurfaceDensityRate=0.0d0
       else
          ! Set a pointer to self and to the node that is accessible by the integrand function.
          self_ => self
          node_ => node
          ! Compute suitable limits for the integration.
          radiusInner=radiusDisk*radiusInnerDimensionless
          radiusOuter=radiusDisk*radiusOuterDimensionless
          ! Get a set of intervals into which this integral should be broken.
          intervals=self%starFormationRateSurfaceDensityDisks_%intervals(node,radiusInner,radiusOuter,intervalIsAnalytic,integralsAnalytic)
          ! Compute the star formation rate. A low order integration rule (GSL_Integ_Gauss15) works well here.
          intgrtdSurfaceDensityRate=0.0d0
          do i=1,size(intervals,dim=2)
             if (intervalIsAnalytic(i)) then
                intgrtdSurfaceDensityRate=+                 intgrtdSurfaceDensityRate                                &
                     &                    +                 integralsAnalytic        (                           i )
             else
                intgrtdSurfaceDensityRate=+                 intgrtdSurfaceDensityRate                                &
                     &                    +self%integrator_%integrate                (intervals(1,i),intervals(2,i))
             end if
          end do
          intgrtdSurfaceDensityRate=+2.0d0                     &
               &                    *Pi                        &
               &                    *intgrtdSurfaceDensityRate
          if (allocated(intervalIsAnalytic)) deallocate(intervalIsAnalytic)
          if (allocated(integralsAnalytic )) deallocate(integralsAnalytic )
       end if
       self%starFormationRatePrevious=intgrtdSurfaceDensityRate
    end if
    return
  end function intgrtdSurfaceDensityRate

  double precision function intgrtdSurfaceDensityIntegrand(radius)
    !!{
    Integrand function for the ``integrated surface density'' star formation rate calculation.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    intgrtdSurfaceDensityIntegrand=+                                                       radius  &
         &                         *self_%starFormationRateSurfaceDensityDisks_%rate(node_,radius)
    return
  end function intgrtdSurfaceDensityIntegrand
  
