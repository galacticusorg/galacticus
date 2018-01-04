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

!% Contains a module which implements the Kennicutt-Schmidt star formation rate surface density for galactic disks.

module Star_Formation_Rate_Surface_Density_Disks_KS
  !% Implements the Kennicutt-Schmidt star formation rate surface density for galactic disks.
  use Galacticus_Nodes
  use Kind_Numbers
  implicit none
  private
  public :: Star_Formation_Rate_Surface_Density_Disks_KS_Reset,&
       & Star_Formation_Rate_Surface_Density_Disks_KS_Initialize

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8) :: lastUniqueID                              =-1
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not factors have been precomputed.
  logical                          :: factorsComputed                           =.false.
  !$omp threadprivate(factorsComputed)
  ! Precomputed factors.
  double precision                 :: criticalDensityFactor                             , hydrogenMassFraction
  !$omp threadprivate(criticalDensityFactor,hydrogenMassFraction)
  ! Parameters of the model.
  double precision                 :: starFormationKennicuttSchmidtExponent             , starFormationKennicuttSchmidtExponentTruncated, &
       &                              starFormationKennicuttSchmidtNormalization        , toomreParameterCritical                       , &
       &                              velocityDispersionDiskGas
  logical                          :: starFormationKennicuttSchmidtTruncate

contains

  !# <calculationResetTask>
  !# <unitName>Star_Formation_Rate_Surface_Density_Disks_KS_Reset</unitName>
  !# </calculationResetTask>
  subroutine Star_Formation_Rate_Surface_Density_Disks_KS_Reset(thisNode)
    !% Reset the Kennicutt-Schmidt relation calculation.
    implicit none
    type(treeNode), intent(inout) :: thisNode

    factorsComputed=.false.
    lastUniqueID   =thisNode%uniqueID()
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_KS_Reset

  !# <starFormationRateSurfaceDensityDisksMethod>
  !#  <unitName>Star_Formation_Rate_Surface_Density_Disks_KS_Initialize</unitName>
  !# </starFormationRateSurfaceDensityDisksMethod>
  subroutine Star_Formation_Rate_Surface_Density_Disks_KS_Initialize(starFormationRateSurfaceDensityDisksMethod&
       &,Star_Formation_Rate_Surface_Density_Disk_Get,Star_Formation_Rate_Surface_Density_Disk_Intervals_Get)
    !% Initializes the ``Kennicutt-Schmidt'' disk star formation rate surface density.
    use ISO_Varying_String
    use Input_Parameters
    use Numerical_Constants_Prefixes
    implicit none
    type     (varying_string                                       ), intent(in   )          :: starFormationRateSurfaceDensityDisksMethod
    procedure(Star_Formation_Rate_Surface_Density_Disk_KS          ), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Get
    procedure(Star_Formation_Rate_Surface_Density_Disk_Intervals_KS), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Intervals_Get

    if (starFormationRateSurfaceDensityDisksMethod == 'Kennicutt-Schmidt') then
       Star_Formation_Rate_Surface_Density_Disk_Get           => Star_Formation_Rate_Surface_Density_Disk_KS
       Star_Formation_Rate_Surface_Density_Disk_Intervals_Get => Star_Formation_Rate_Surface_Density_Disk_Intervals_KS
       ! Get parameters of for the timescale calculation.
       !# <inputParameter>
       !#   <name>starFormationKennicuttSchmidtNormalization</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{kennicutt_global_1998}</defaultSource>
       !#   <defaultValue>0.147d0</defaultValue>
       !#   <description>The normalization of the Kennicutt-Schmidt star formation law [$M_\odot$ Gyr$^{-1}$pc$^{-2}$].</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>starFormationKennicuttSchmidtExponent</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{kennicutt_global_1998}</defaultSource>
       !#   <defaultValue>1.400d0</defaultValue>
       !#   <description>The exponent in the Kennicutt-Schmidt star formation law.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>starFormationKennicuttSchmidtTruncate</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>.true.</defaultValue>
       !#   <description>Specifies whether or not to truncate star formation below a critical surface density in disks.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>boolean</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>starFormationKennicuttSchmidtExponentTruncated</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>6.0d0</defaultValue>
       !#   <description>The exponent of the $\Sigma_{\mathrm gas}/\Sigma_{\mathrm crit}$ term used in truncating the Kennicutt-Schmidt star formation law.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>velocityDispersionDiskGas</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{leroy_star_2008}</defaultSource>
       !#   <defaultValue>10.0d0</defaultValue>
       !#   <description>The velocity dispersion of gas in disks.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>toomreParameterCritical</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{kennicutt_star_1989}</defaultSource>
       !#   <defaultValue>0.4d0</defaultValue>
       !#   <description>The critical Toomre parameter for star formation in disks.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       ! Renormalize the Kennicutt-Schmidt relation to our internal units.
       starFormationKennicuttSchmidtNormalization=                       &
            &  starFormationKennicuttSchmidtNormalization                &
            & *mega**(2.0d0-2.0d0*starFormationKennicuttSchmidtExponent)
    end if
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_KS_Initialize

  double precision function Star_Formation_Rate_Surface_Density_Disk_KS(thisNode,radius)
    !% Returns the star formation rate surface density  (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation in the galactic disk of {\normalfont \ttfamily thisNode}. The disk is assumed to obey the Kennicutt-Schmidt law:
    !% \begin{equation}
    !% \Sigma_\star = A \left(x_{\mathrm H} {\Sigma_{\mathrm gas}\over M_\odot \hbox{pc}^{-2}}\right)^N,
    !% \end{equation}
    !% where $A=${\normalfont \ttfamily [starFormationKennicuttSchmidtNormalization]} and $N=${\tt
    !% [starFormationKennicuttSchmidtExponent]}. Optionally, star formation is truncated for gas surface densities below a critical density of:
    !% \begin{equation}
    !% \Sigma_{\mathrm crit} = {q_{\mathrm crit} \kappa \sigma_{\mathrm gas} \over \pi \G},
    !% \end{equation}
    !% where $\kappa$ is the epicyclic frequency in the disk, $\sigma_{\mathrm gas}$ is the velocity dispersion of gas in the disk and
    !% $q_{\mathrm crit}=${\normalfont \ttfamily [toomreParameterCritical]} is a dimensionless constant of order unity which controls where the critical
    !% density occurs. $\sigma_{\mathrm gas}$ is assumed to be a constant equal to {\normalfont \ttfamily [velocityDispersionDiskGas]} and the disk is
    !% assumed to have a flat rotation curve such that $\kappa = \sqrt{2} V/R$.
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    use Abundances_Structure
    use Galactic_Structure_Surface_Densities
    use Galactic_Structure_Options
    implicit none
    type            (treeNode         ), intent(inout) :: thisNode
    double precision                   , intent(in   ) :: radius
    class           (nodeComponentDisk), pointer       :: thisDiskComponent
    type            (abundances       ), save          :: fuelAbundances
    !$omp threadprivate(fuelAbundances)
    double precision                                   :: criticalDensity  , gasMass, surfaceDensityGas

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Star_Formation_Rate_Surface_Density_Disks_KS_Reset(thisNode)
    ! Check if factors have been precomputed.
    if (.not.factorsComputed) then
       ! Get the disk properties.
       thisDiskComponent => thisNode         %disk   ()
       gasMass             =thisDiskComponent%massGas()
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       fuelAbundances=thisDiskComponent%abundancesGas()
       call fuelAbundances%massToMassFraction(gasMass)
       hydrogenMassFraction=fuelAbundances%hydrogenMassFraction()
       ! Compute the constant factor appearing in the critical density.
       criticalDensityFactor=toomreParameterCritical*sqrt(2.0d0)&
            &*velocityDispersionDiskGas*thisDiskComponent%velocity()/Pi/gravitationalConstantGalacticus
       ! Record that factors have now been computed.
       factorsComputed=.true.
    end if
    ! Get gas surface density.
    surfaceDensityGas=Galactic_Structure_Surface_Density(thisNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeGaseous)
    ! Compute the star formation rate surface density.
    Star_Formation_Rate_Surface_Density_Disk_KS=                             &
         &  starFormationKennicuttSchmidtNormalization                                      &
         & *(hydrogenMassFraction*surfaceDensityGas)**starFormationKennicuttSchmidtExponent
    ! Check if we are applying a truncation radius.
    if (starFormationKennicuttSchmidtTruncate) then
       ! Always return zero star formation rate at zero radius, as critical density will be infinite.
       if (radius <= 0.0d0) then
          Star_Formation_Rate_Surface_Density_Disk_KS=0.0d0
          return
       end if
       ! Compute the critical density for star formation.
       criticalDensity=criticalDensityFactor/radius
       ! Check if gas is above the critical density. Return zero star formation rate if it is not.
       if (surfaceDensityGas < criticalDensity) Star_Formation_Rate_Surface_Density_Disk_KS&
            &=Star_Formation_Rate_Surface_Density_Disk_KS*(surfaceDensityGas/criticalDensity)&
            &**starFormationKennicuttSchmidtExponentTruncated
    end if
    return
  end function Star_Formation_Rate_Surface_Density_Disk_KS

  function Star_Formation_Rate_Surface_Density_Disk_Intervals_KS(thisNode,radiusInner,radiusOuter)
    !% Returns intervals to use for integrating the Kennicutt-Schmidt star formation rate over a galactic disk.
    implicit none
    double precision          , allocatable  , dimension(:,:) :: Star_Formation_Rate_Surface_Density_Disk_Intervals_KS
    type            (treeNode), intent(inout), target         :: thisNode
    double precision          , intent(in   )                 :: radiusInner, radiusOuter
    !GCC$ attributes unused :: thisNode

    allocate(Star_Formation_Rate_Surface_Density_Disk_Intervals_KS(2,1))
    Star_Formation_Rate_Surface_Density_Disk_Intervals_KS=reshape([radiusInner,radiusOuter],[2,1])
    return
  end function Star_Formation_Rate_Surface_Density_Disk_Intervals_KS

end module Star_Formation_Rate_Surface_Density_Disks_KS
