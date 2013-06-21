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

!% Contains a module which implements the \cite{krumholz_star_2009} star formation rate surface density law for galactic disks.

module Star_Formation_Rate_Surface_Density_Disks_KMT09
  !% Implements the \cite{krumholz_star_2009} star formation rate surface density law for galactic disks.
  use Galacticus_Nodes
  use Kind_Numbers
  implicit none
  private
  public :: Star_Formation_Rate_Surface_Density_Disks_KMT09_Reset,&
       & Star_Formation_Rate_Surface_Density_Disks_KMT09_Initialize

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8  )            :: lastUniqueID                       =-1                                      
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not factors have been precomputed.
  logical                                       :: factorsComputed                    =.false.                                 
  !$omp threadprivate(factorsComputed)
  ! Precomputed factors.
  double precision                              :: chi                                        , diskScaleRadius            , & 
       &                                           gasMass                                    , hydrogenMassFraction       , & 
       &                                           metallicityRelativeToSolar                 , sNormalization             , & 
       &                                           sigmaMolecularComplexNormalization                                          
  !$omp threadprivate(hydrogenMassFraction,gasMass,diskScaleRadius,metallicityRelativeToSolar)  !$omp threadprivate(chi,sigmaMolecularComplexNormalization,sNormalization)
  ! Parameters of the model.
  double precision                              :: molecularComplexClumpingFactorKMT09        , starFormationFrequencyKMT09    
  
  ! Pointer the the molecular fraction function that is to be used.
  procedure       (double precision), pointer   :: KMT09_Molecular_Fraction           =>null()                                 
  
  ! Minimum fraction of molecular hydrogen allowed.
  double precision                  , parameter :: molecularFractionMinimum           =1.0d-4                                  
  
contains

  !# <calculationResetTask>
  !# <unitName>Star_Formation_Rate_Surface_Density_Disks_KMT09_Reset</unitName>
  !# </calculationResetTask>
  subroutine Star_Formation_Rate_Surface_Density_Disks_KMT09_Reset(thisNode)
    !% Reset the extended Schmidt relation calculation.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode 
    
    factorsComputed=.false.
    lastUniqueID   =thisNode%uniqueID()
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_KMT09_Reset

  !# <starFormationRateSurfaceDensityDisksMethod>
  !#  <unitName>Star_Formation_Rate_Surface_Density_Disks_KMT09_Initialize</unitName>
  !# </starFormationRateSurfaceDensityDisksMethod>
  subroutine Star_Formation_Rate_Surface_Density_Disks_KMT09_Initialize(starFormationRateSurfaceDensityDisksMethod&
       &,Star_Formation_Rate_Surface_Density_Disk_Get)
    !% Initializes the ``KMT09'' disk star formation rate surface density.
    use ISO_Varying_String
    use Input_Parameters
    use Abundances_Structure
    use Numerical_Constants_Prefixes
    implicit none
    type     (varying_string                                ), intent(in   )          :: starFormationRateSurfaceDensityDisksMethod   
    procedure(Star_Formation_Rate_Surface_Density_Disk_KMT09), intent(inout), pointer :: Star_Formation_Rate_Surface_Density_Disk_Get 
    logical                                                                           :: molecularFractionFastKMT09                   
    
    if (starFormationRateSurfaceDensityDisksMethod == 'KMT09') then
       Star_Formation_Rate_Surface_Density_Disk_Get => Star_Formation_Rate_Surface_Density_Disk_KMT09
       ! Get parameters of our model.
       !@ <inputParameter>
       !@   <name>starFormationFrequencyKMT09</name>
       !@   <defaultValue>$0.385$ \citep{krumholz_star_2009}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The star formation frequency (in units of Gyr$^{-1}$) in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationFrequencyKMT09',starFormationFrequencyKMT09,defaultValue=0.385d0)
       !@ <inputParameter>
       !@   <name>molecularComplexClumpingFactorKMT09</name>
       !@   <defaultValue>5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The density enhancement (relative to mean disk density) for molecular complexes in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('molecularComplexClumpingFactorKMT09',molecularComplexClumpingFactorKMT09,defaultValue=5.0d0)
       !@ <inputParameter>
       !@   <name>molecularFractionFastKMT09</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@      Selects whether the fast (but less accurate) fitting formula for molecular hydrogen should be used in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('molecularFractionFastKMT09',molecularFractionFastKMT09,defaultValue=.false.)
       ! Set a pointer to the molecular hydrogen fraction fitting function to be used.
       select case (molecularFractionFastKMT09)
       case (.true.)
          KMT09_Molecular_Fraction => KMT09_Molecular_Fraction_Fast
       case(.false.)
          KMT09_Molecular_Fraction => KMT09_Molecular_Fraction_Slow
       end select
    end if
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_KMT09_Initialize

  double precision function Star_Formation_Rate_Surface_Density_Disk_KMT09(thisNode,radius)
    !% Returns the star formation rate surface density (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation
    !% in the galactic disk of {\tt thisNode}. The disk is assumed to obey the
    !% \cite{krumholz_star_2009} star formation rule.
    use Galacticus_Nodes
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    use Abundances_Structure
    use Galactic_Structure_Surface_Densities
    use Galactic_Structure_Options
    use Numerical_Constants_Prefixes
    implicit none
    type            (treeNode         ), intent(inout), pointer :: thisNode                                                                              
    double precision                   , intent(in   )          :: radius                                                                                
    class           (nodeComponentDisk)               , pointer :: thisDiskComponent                                                                     
    type            (abundances       ), save                   :: fuelAbundances                                                                        
    !$omp threadprivate(fuelAbundances)
    double precision                   , parameter              :: surfaceDensityTransition=85.0d12                                 !   M_Solar/Mpc^2    
    double precision                                            :: cloudFactor                     , molecularFraction                               , & 
         &                                                         s                               , sigmaMolecularComplex                           , & 
         &                                                         surfaceDensityGas               , surfaceDensityGasDimensionless                      
    
    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Star_Formation_Rate_Surface_Density_Disks_KMT09_Reset(thisNode)
    ! Check if factors have been precomputed.
    if (.not.factorsComputed) then
       ! Get the disk properties.
       thisDiskComponent => thisNode       %disk   ()
       gasMass           =thisDiskComponent%massGas()
       diskScaleRadius   =thisDiskComponent%radius ()
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       fuelAbundances=thisDiskComponent%abundancesGas()
       call fuelAbundances%massToMassFraction(gasMass)
       hydrogenMassFraction=fuelAbundances%hydrogenMassFraction()
       ! Get the metallicity in Solar units, and related quantities.
       metallicityRelativeToSolar=fuelAbundances%metallicity(linearByMassSolar)
       if (metallicityRelativeToSolar > 0.0d0) then
          chi                               =0.77d0*(1.0d0+3.1d0*metallicityRelativeToSolar**0.365d0)
          sigmaMolecularComplexNormalization=hydrogenMassFraction*molecularComplexClumpingFactorKMT09/mega**2
          sNormalization                    =log(1.0d0+0.6d0*chi+0.01d0*chi**2)/(0.04d0*metallicityRelativeToSolar)
       end if
       ! Record that factors have now been computed.
       factorsComputed=.true.
    end if
    ! Check if the disk is physical.
    if (gasMass <= 0.0d0 .or. diskScaleRadius <= 0.0d0) then
       ! It is not, so return zero rate.
       Star_Formation_Rate_Surface_Density_Disk_KMT09=0.0d0
    else
       ! Get gas surface density.
       surfaceDensityGas=Galactic_Structure_Surface_Density(thisNode,[radius,0.0d0,0.0d0],coordinateSystem&
            &=coordinateSystemCylindrical,componentType=componentTypeDisk,massType=massTypeGaseous)
       ! Compute the molecular fraction.
       if (metallicityRelativeToSolar > 0.0d0) then
          sigmaMolecularComplex=sigmaMolecularComplexNormalization*surfaceDensityGas
          s                    =sNormalization/sigmaMolecularComplex
          molecularFraction    =KMT09_Molecular_Fraction(s)
       else
          molecularFraction    =molecularFractionMinimum
       end if
       ! Compute the cloud density factor.
       surfaceDensityGasDimensionless=hydrogenMassFraction*surfaceDensityGas/surfaceDensityTransition
       if (surfaceDensityGasDimensionless < 1.0d0) then
          cloudFactor=surfaceDensityGasDimensionless**(-0.33d0)
       else
          cloudFactor=surfaceDensityGasDimensionless**(+0.33d0)
       end if
       ! Compute the star formation rate surface density.
       Star_Formation_Rate_Surface_Density_Disk_KMT09=                             &
            &                                          starFormationFrequencyKMT09 &
            &                                         *surfaceDensityGas           &
            &                                         *cloudFactor                 &
            &                                         *molecularFraction
    end if
    return
  end function Star_Formation_Rate_Surface_Density_Disk_KMT09

  double precision function KMT09_Molecular_Fraction_Slow(s)
    !% Slow (but more accurate at low molecular fraction) fitting function from \cite{krumholz_star_2009} for the molecular
    !% hydrogen fraction.
    implicit none
    double precision, intent(in   ) :: s               
    double precision, parameter     :: sMinimum=1.0d-6 
    double precision, parameter     :: sMaximum=8.0d0  
    double precision                :: delta           
    
    ! Check if s is below maximum. If not, simply truncate to the minimum fraction that we allow. Also use a simple series    ! expansion for cases of very small s.
    if      (s <  sMinimum) then
       KMT09_Molecular_Fraction_Slow=1.0d0-0.75d0*s
    else if (s >= sMaximum) then
       KMT09_Molecular_Fraction_Slow=                                                               molecularFractionMinimum
    else
       delta                        =0.0712d0/((0.1d0/s+0.675d0)**2.8d0)
       KMT09_Molecular_Fraction_Slow=max(1.0d0-1.0d0/((1.0d0+(((1.0d0+delta)/0.75d0/s)**5))**0.2d0),molecularFractionMinimum)
    end if
    return
  end function KMT09_Molecular_Fraction_Slow

  double precision function KMT09_Molecular_Fraction_Fast(s)
    !% Fast (but less accurate at low molecular fraction) fitting function from \cite{mckee_atomic--molecular_2010} for the
    !% molecular hydrogen fraction.
    implicit none
    double precision, intent(in   ) :: s 
    
    ! Check that s is below 2 - if it is, compute the molecular fraction, otherwise truncate to the minimum.
    if (s < 2.0d0) then
       KMT09_Molecular_Fraction_Fast=max(1.0d0-0.75d0*s/(1.0d0+0.25d0*s),molecularFractionMinimum)
    else
       KMT09_Molecular_Fraction_Fast=                                    molecularFractionMinimum
    end if
    return
  end function KMT09_Molecular_Fraction_Fast

end module Star_Formation_Rate_Surface_Density_Disks_KMT09
