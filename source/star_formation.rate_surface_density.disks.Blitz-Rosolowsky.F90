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

!% Contains a module which implements the \cite{blitz_role_2006} star formation rate surface density law for galactic disks.

module Star_Formation_Rate_Surface_Density_Disks_BR
  !% Implements the \cite{blitz_role_2006} star formation rate surface density law for galactic disks.
  use Galacticus_Nodes
  use Kind_Numbers
  implicit none
  private
  public :: Star_Formation_Rate_Surface_Density_Disks_BR_Reset,&
       & Star_Formation_Rate_Surface_Density_Disks_BR_Initialize

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8) :: lastUniqueID                          =-1
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not factors have been precomputed.
  logical                          :: factorsComputed                       =.false.
  !$omp threadprivate(factorsComputed)
  ! Precomputed factors.
  double precision                 :: diskScaleRadius                               , gasMass                                           , &
       &                              hydrogenMassFraction                          , stellarMass
  !$omp threadprivate(hydrogenMassFraction,diskScaleRadius,gasMass,stellarMass)
  ! Parameters of the model.
  double precision                 :: heightToRadialScaleDiskBlitzRosolowsky        , pressureCharacteristicBlitzRosolowsky             , &
       &                              pressureExponentBlitzRosolowsky               , starFormationFrequencyNormalizationBlitzRosolowsky, &
       &                              surfaceDensityCriticalBlitzRosolowsky         , surfaceDensityExponentBlitzRosolowsky             , &
       &                              velocityDispersionDiskGas

contains

  !# <calculationResetTask>
  !# <unitName>Star_Formation_Rate_Surface_Density_Disks_BR_Reset</unitName>
  !# </calculationResetTask>
  subroutine Star_Formation_Rate_Surface_Density_Disks_BR_Reset(thisNode)
    !% Reset the extended Schmidt relation calculation.
    implicit none
    type(treeNode), intent(inout) :: thisNode

    factorsComputed=.false.
    lastUniqueID   =thisNode%uniqueID()
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_BR_Reset

  !# <starFormationRateSurfaceDensityDisksMethod>
  !#  <unitName>Star_Formation_Rate_Surface_Density_Disks_BR_Initialize</unitName>
  !# </starFormationRateSurfaceDensityDisksMethod>
  subroutine Star_Formation_Rate_Surface_Density_Disks_BR_Initialize(starFormationRateSurfaceDensityDisksMethod&
       &,Star_Formation_Rate_Surface_Density_Disk_Get,Star_Formation_Rate_Surface_Density_Disk_Intervals_Get,Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get)
    !% Initializes the ``extended Schmidt'' disk star formation rate surface density.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    implicit none
    type     (varying_string                                       ), intent(in   )          :: starFormationRateSurfaceDensityDisksMethod
    procedure(Star_Formation_Rate_Surface_Density_Disk_BR          ), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Get
    procedure(Star_Formation_Rate_Surface_Density_Disk_Intervals_BR), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Intervals_Get
    procedure(Star_Formation_Rate_Surface_Density_Disk_Unchanged_BR), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get

    if (starFormationRateSurfaceDensityDisksMethod == 'Blitz-Rosolowsky2006') then
       Star_Formation_Rate_Surface_Density_Disk_Get           => Star_Formation_Rate_Surface_Density_Disk_BR
       Star_Formation_Rate_Surface_Density_Disk_Intervals_Get => Star_Formation_Rate_Surface_Density_Disk_Intervals_BR
       Star_Formation_Rate_Surface_Density_Disk_Unchanged_Get => Star_Formation_Rate_Surface_Density_Disk_Unchanged_BR
       ! Get parameters of for the timescale calculation.
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
       !#   <name>heightToRadialScaleDiskBlitzRosolowsky</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{kregel_flattening_2002}</defaultSource>
       !#   <defaultValue>0.137d0</defaultValue>
       !#   <description>The ratio of scale height to scale radius for disks in the ``Blitz-Rosolowsky'' star formation timescale calculation.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>surfaceDensityCriticalBlitzRosolowsky</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{bigiel_star_2008}</defaultSource>
       !#   <defaultValue>200.0d0</defaultValue>
       !#   <description>The surface density (in units of $M_\odot$ pc$^{-2}$) in the ``Blitz-Rosolowsky'' star formation timescale calculation at which low-density truncation begins.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>surfaceDensityExponentBlitzRosolowsky</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{bigiel_star_2008}</defaultSource>
       !#   <defaultValue>0.4d0</defaultValue>
       !#   <description>The exponent for surface density in the ``Blitz-Rosolowsky'' star formation timescale calculation at in the high density regime.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>starFormationFrequencyNormalizationBlitzRosolowsky</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{leroy_star_2008}</defaultSource>
       !#   <defaultValue>5.25d-10</defaultValue>
       !#   <description>The star formation frequency (in the low-density limit and in units of yr$^{-1}$) in the ``Blitz-Rosolowsky'' star formation timescale calculation.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>pressureCharacteristicBlitzRosolowsky</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{blitz_role_2006}</defaultSource>
       !#   <defaultValue>4.54d0</defaultValue>
       !#   <description>The characteristic pressure (given as $P_0/k_\mathrm{B}$ in units of K cm$^{-3}$) in the scaling relation of molecular hydrogen fraction with disk pressure in the ``Blitz-Rosolowsky'' star formation timescale calculation.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>pressureExponentBlitzRosolowsky</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{blitz_role_2006}</defaultSource>
       !#   <defaultValue>0.92d0</defaultValue>
       !#   <description>The exponent in the scaling relation of molecular hydrogen fraction with disk pressure in the ``Blitz-Rosolowsky'' star formation timescale calculation.</description>
       !#   <group>starFormation</group>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       if (pressureExponentBlitzRosolowsky < 0.0d0) call Galacticus_Error_Report('pressureExponentBlitzRosolowsky < 0 violates assumptions'//{introspection:location})
       ! Convert parameters to internal units.
       surfaceDensityCriticalBlitzRosolowsky             =surfaceDensityCriticalBlitzRosolowsky*(mega**2)                                                    ! Convert to M_Solar/Mpc^2.
       starFormationFrequencyNormalizationBlitzRosolowsky=starFormationFrequencyNormalizationBlitzRosolowsky*giga                                            ! Convert to Gyr^-1.
       pressureCharacteristicBlitzRosolowsky             =pressureCharacteristicBlitzRosolowsky*boltzmannsConstant*((hecto*megaParsec)**3)/massSolar/kilo**2 ! Convert to M_Solar (km/s)^2 / Mpc.
    end if
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_BR_Initialize

  double precision function Star_Formation_Rate_Surface_Density_Disk_BR(thisNode,radius)
    !% Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation
    !% in the galactic disk of {\normalfont \ttfamily thisNode}. The disk is assumed to obey the
    !% \cite{blitz_role_2006} star formation rule.
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
    double precision                                   :: molecularFraction , pressureRatio, surfaceDensityGas, &
         &                                                surfaceDensityStar

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Star_Formation_Rate_Surface_Density_Disks_BR_Reset(thisNode)
    ! Check if factors have been precomputed.
    if (.not.factorsComputed) then
       ! Get the disk properties.
       thisDiskComponent => thisNode       %disk       ()
       gasMass           =thisDiskComponent%massGas    ()
       stellarMass       =thisDiskComponent%massStellar()
       diskScaleRadius   =thisDiskComponent%radius     ()
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       fuelAbundances=thisDiskComponent%abundancesGas()
       call fuelAbundances%massToMassFraction(gasMass)
       hydrogenMassFraction=fuelAbundances%hydrogenMassFraction()
       ! Record that factors have now been computed.
       factorsComputed=.true.
    end if
    ! Return zero rate for non-positive radius or mass.
    if (gasMass <= 0.0d0 .or. stellarMass < 0.0d0 .or. diskScaleRadius <= 0.0d0) then
       Star_Formation_Rate_Surface_Density_Disk_BR=0.0d0
       return
    end if
    ! Get gas surface density.
    surfaceDensityGas=Galactic_Structure_Surface_Density(thisNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeGaseous)
    ! Get stellar surface density.
    surfaceDensityStar=Galactic_Structure_Surface_Density(thisNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeStellar)
    ! Compute the pressure ratio that Blitz & Rosolowsky use to compute the molecular fraction.
    pressureRatio=0.5d0*Pi*gravitationalConstantGalacticus*surfaceDensityGas*(surfaceDensityGas+velocityDispersionDiskGas &
         &*sqrt(surfaceDensityStar/Pi/gravitationalConstantGalacticus/heightToRadialScaleDiskBlitzRosolowsky&
         &/diskScaleRadius))/pressureCharacteristicBlitzRosolowsky
    ! Compute the molecular fraction, limited to 100% molecular.
    if (pressureRatio >= 1.0d0) then
       molecularFraction=1.0d0
    else
       molecularFraction=min(pressureRatio**pressureExponentBlitzRosolowsky,1.0d0)
    end if
    ! Compute the star formation rate surface density.
    Star_Formation_Rate_Surface_Density_Disk_BR=surfaceDensityGas*hydrogenMassFraction*molecularFraction&
         &*starFormationFrequencyNormalizationBlitzRosolowsky*(1.0d0+(hydrogenMassFraction*surfaceDensityGas&
         &/surfaceDensityCriticalBlitzRosolowsky)**surfaceDensityExponentBlitzRosolowsky)
    return
  end function Star_Formation_Rate_Surface_Density_Disk_BR

  function Star_Formation_Rate_Surface_Density_Disk_Intervals_BR(thisNode,radiusInner,radiusOuter)
    !% Returns intervals to use for integrating the Blitz-Rosolowsky star formation rate over a galactic disk.
    implicit none
    double precision          , allocatable  , dimension(:,:) :: Star_Formation_Rate_Surface_Density_Disk_Intervals_BR
    type            (treeNode), intent(inout), target         :: thisNode
    double precision          , intent(in   )                 :: radiusInner, radiusOuter
    !GCC$ attributes unused :: thisNode
    
    allocate(Star_Formation_Rate_Surface_Density_Disk_Intervals_BR(2,1))
    Star_Formation_Rate_Surface_Density_Disk_Intervals_BR=reshape([radiusInner,radiusOuter],[2,1])
    return
  end function Star_Formation_Rate_Surface_Density_Disk_Intervals_BR
  
  logical function Star_Formation_Rate_Surface_Density_Disk_Unchanged_BR(thisNode)
    !% Claim that the surface rate density of star formation is unchanged so that it is always re-evaluated.
    implicit none
    type(treeNode), intent(inout) :: thisNode
    !GCC$ attributes unused :: thisNode

    Star_Formation_Rate_Surface_Density_Disk_Unchanged_BR=.false.
    return
  end function Star_Formation_Rate_Surface_Density_Disk_Unchanged_BR

end module Star_Formation_Rate_Surface_Density_Disks_BR
