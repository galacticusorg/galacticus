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

!+    Contributions to this file made by: Arya Farahi, Andrew Benson.

!% Contains a module which implements the extended Schmidt star formation rate surface density law of \cite{shi_extended_2011} for galactic disks.

module Star_Formation_Rate_Surface_Density_Disks_ExSchmidt
  !% Implements the extended Schmidt star formation rate surface density law of \cite{shi_extended_2011} for galactic disks.
  use Galacticus_Nodes
  use Kind_Numbers
  implicit none
  private
  public :: Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Reset,&
       & Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Initialize

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8) :: lastUniqueID                            =-1
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not factors have been precomputed.
  logical                          :: factorsComputed                         =.false.
  !$omp threadprivate(factorsComputed)
  ! Precomputed factors.
  double precision                 :: hydrogenMassFraction
  !$omp threadprivate(hydrogenMassFraction)
  ! Parameters of the model.
  double precision                 :: starFormationExtendedSchmidtGasExponent         , starFormationExtendedSchmidtNormalization, &
       &                              starFormationExtendedSchmidtStarExponent

contains

  !# <calculationResetTask>
  !# <unitName>Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Reset</unitName>
  !# </calculationResetTask>
  subroutine Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Reset(thisNode)
    !% Reset the extended Schmidt relation calculation.
    implicit none
    type(treeNode), intent(inout) :: thisNode

    factorsComputed=.false.
    lastUniqueID   =thisNode%uniqueID()
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Reset

  !# <starFormationRateSurfaceDensityDisksMethod>
  !#  <unitName>Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Initialize</unitName>
  !# </starFormationRateSurfaceDensityDisksMethod>
  subroutine Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Initialize(starFormationRateSurfaceDensityDisksMethod&
       &,Star_Formation_Rate_Surface_Density_Disk_Get,Star_Formation_Rate_Surface_Density_Disk_Intervals_Get,Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get)
    !% Initializes the ``extended Schmidt'' disk star formation rate surface density.
    use ISO_Varying_String
    use Input_Parameters
    use Numerical_Constants_Prefixes
    implicit none
    type     (varying_string                                              ), intent(in   )          :: starFormationRateSurfaceDensityDisksMethod
    procedure(Star_Formation_Rate_Surface_Density_Disk_ExSchmidt          ), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Get
    procedure(Star_Formation_Rate_Surface_Density_Disk_Intervals_ExSchmidt), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Intervals_Get
    procedure(Star_Formation_Rate_Surface_Density_Disk_Unchanged_ExSchmidt), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get

    if (starFormationRateSurfaceDensityDisksMethod == 'extendedSchmidt') then
       Star_Formation_Rate_Surface_Density_Disk_Get           => Star_Formation_Rate_Surface_Density_Disk_ExSchmidt
       Star_Formation_Rate_Surface_Density_Disk_Intervals_Get => Star_Formation_Rate_Surface_Density_Disk_Intervals_ExSchmidt
       Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get => Star_Formation_Rate_Surface_Density_Disk_Unchanged_ExSchmidt
       ! Get parameters of for the timescale calculation.
       !# <inputParameter>
       !#   <name>starFormationExtendedSchmidtNormalization</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{shi_extended_2011}</defaultSource>
       !#   <defaultValue>0.5248d-10</defaultValue>
       !#   <description>The normalization of the extended Schmidt star formation law [$M_\odot$ yr$^{-1}$pc$^{-2}$].</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>starFormationExtendedSchmidtGasExponent</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{shi_extended_2011}</defaultSource>
       !#   <defaultValue>1.0000d+0</defaultValue>
       !#   <description>The exponent of gas surface density in the extended Schmidt star formation law.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>starFormationExtendedSchmidtStarExponent</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{shi_extended_2011}</defaultSource>
       !#   <defaultValue>0.4800d+0</defaultValue>
       !#   <description>The exponent of stellar surface density in the extended Schmidt star formation law.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       ! Renormalize the relation to internal units.
       starFormationExtendedSchmidtNormalization=                          &
            &  starFormationExtendedSchmidtNormalization                   &
            & *(mega**2)*giga                                              & ! Convert from Msun/pc^2/yr to Msun/Mpc^2/Gyr
            & *((1.0d0/mega**2)**starFormationExtendedSchmidtStarExponent) & ! Unit conversion for stars.
            & *((1.0d0/mega**2)**starFormationExtendedSchmidtGasExponent )   ! Hydrogen fraction and unit conversion for gas.
    end if
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Initialize

  double precision function Star_Formation_Rate_Surface_Density_Disk_ExSchmidt(thisNode,radius)
    !% Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation
    !% in the galactic disk of {\normalfont \ttfamily thisNode}. The disk is assumed to obey the extended Schmidt law of \cite{shi_extended_2011}:
    !% \begin{equation}
    !% \dot{\Sigma}_\star = A \left(x_{\mathrm H} {\Sigma_{\mathrm gas}\over M_\odot \hbox{pc}^{-2}}\right)
    !% ^{N_1} \left({\Sigma_{\star}\over M_\odot \hbox{pc}^{-2}}\right)^{N_2},
    !% \end{equation}
    !% where $A=${\normalfont \ttfamily [starFormationExtendedSchmidtNormalization]} and $N_1=${\tt
    !% [starFormationExtendedSchmidtGasExponent]}. $N_2=${\normalfont \ttfamily [starFormationExtendedSchmidtStarExponent]}.
    use Abundances_Structure
    use Galactic_Structure_Surface_Densities
    use Galactic_Structure_Options
    implicit none
    type            (treeNode         ), intent(inout) :: thisNode
    double precision                   , intent(in   ) :: radius
    class           (nodeComponentDisk), pointer       :: thisDiskComponent
    type            (abundances       ), save          :: fuelAbundances
    !$omp threadprivate(fuelAbundances)
    double precision                                   :: gasMass          , surfaceDensityGas, surfaceDensityStar

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Reset(thisNode)
    ! Check if factors have been precomputed.
    if (.not.factorsComputed) then
       ! Get the disk properties.
       thisDiskComponent => thisNode         %disk   ()
       gasMass             =thisDiskComponent%massGas()
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       fuelAbundances=thisDiskComponent%abundancesGas()
       call fuelAbundances%massToMassFraction(gasMass)
       hydrogenMassFraction=fuelAbundances%hydrogenMassFraction()
       ! Record that factors have now been computed.
       factorsComputed=.true.
    end if
    ! Return zero rate for non-positive radius.
    if (radius <= 0.0d0) then
       Star_Formation_Rate_Surface_Density_Disk_ExSchmidt=0.0d0
       return
    end if
    ! Get gas surface density.
    surfaceDensityGas=Galactic_Structure_Surface_Density(thisNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeGaseous)
    ! Get stellar surface density.
    surfaceDensityStar=Galactic_Structure_Surface_Density(thisNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeStellar)
    ! Compute the star formation rate surface density.
    Star_Formation_Rate_Surface_Density_Disk_ExSchmidt=                                           &
         &  starFormationExtendedSchmidtNormalization                                             & ! Normalization of the star formation rate.
         & *((hydrogenMassFraction*surfaceDensityGas )**starFormationExtendedSchmidtGasExponent ) &
         & *(                      surfaceDensityStar **starFormationExtendedSchmidtStarExponent)
    return
  end function Star_Formation_Rate_Surface_Density_Disk_ExSchmidt

  function Star_Formation_Rate_Surface_Density_Disk_Intervals_ExSchmidt(thisNode,radiusInner,radiusOuter)
    !% Returns intervals to use for integrating the extended Kennicutt-Schmidt star formation rate over a galactic disk.
    implicit none
    double precision          , allocatable  , dimension(:,:) :: Star_Formation_Rate_Surface_Density_Disk_Intervals_ExSchmidt
    type            (treeNode), intent(inout), target         :: thisNode
    double precision          , intent(in   )                 :: radiusInner, radiusOuter
    !GCC$ attributes unused :: thisNode

    allocate(Star_Formation_Rate_Surface_Density_Disk_Intervals_ExSchmidt(2,1))
    Star_Formation_Rate_Surface_Density_Disk_Intervals_ExSchmidt=reshape([radiusInner,radiusOuter],[2,1])
    return
  end function Star_Formation_Rate_Surface_Density_Disk_Intervals_ExSchmidt
  
  logical function Star_Formation_Rate_Surface_Density_Disk_Unchanged_ExSchmidt(thisNode)
    !% Claim that the surface rate density of star formation is unchanged so that it is always re-evaluated.
    implicit none
    type(treeNode), intent(inout) :: thisNode
    !GCC$ attributes unused :: thisNode

    Star_Formation_Rate_Surface_Density_Disk_Unchanged_ExSchmidt=.false.
    return
  end function Star_Formation_Rate_Surface_Density_Disk_Unchanged_ExSchmidt

end module Star_Formation_Rate_Surface_Density_Disks_ExSchmidt
