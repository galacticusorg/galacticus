!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements the \cite{krumholz_star_2009} star formation timescale for galactic disks.

module Star_Formation_Timescale_Disks_KMT09
  !% Implements the \cite{krumholz_star_2009} star formation timescale for galactic disks.
  use Tree_Nodes
  private
  public :: Star_Formation_Timescale_Disks_KMT09_Initialize

  ! Internal copy of the number of abundances properties.
  integer          :: abundancesCount

  ! Parameters of the model.
  double precision :: starFormationFrequencyKMT09,molecularComplexClumpingFactorKMT09

  ! Module global variables use in integrations.
  double precision :: hydrogenMassFraction,chi,sNormalization,sigmaMolecularComplexNormalization,metallicityRelativeToSolar
  !$omp threadprivate(hydrogenMassFraction,chi,sNormalization,sigmaMolecularComplexNormalization,metallicityRelativeToSolar)

  ! Pointer the the molecular fraction function that is to be used.
  procedure(double precision), pointer   :: KMT09_Molecular_Fraction => null()

  ! Minimum fraction of molecular hydrogen allowed.
  double precision,            parameter :: molecularFractionMinimum=1.0d-4

  ! Pointer to active node used in integral functions, plus variables needed by integral function.
  type(treeNode),              pointer   :: activeNode
  !$omp threadprivate(activeNode)

contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_KMT09_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_KMT09_Initialize(starFormationTimescaleDisksMethod&
       &,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``Krumholz-McKee-Tumlinson'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    use Abundances_Structure
    implicit none
    type(varying_string),                 intent(in)    :: starFormationTimescaleDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Timescale_Disk_Get
    logical                                             :: molecularFractionFastKMT09
    
    if (starFormationTimescaleDisksMethod == 'KMT09') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_KMT09

       !@ <inputParameter>
       !@   <name>starFormationFrequencyKMT09</name>
       !@   <defaultValue>$0.385$ \citep{krumholz_star_2009}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The star formation frequency (in units of Gyr$^{-1}$) in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationFrequencyKMT09',starFormationFrequencyKMT09,defaultValue=0.385d0)
       !@ <inputParameter>
       !@   <name>molecularComplexClumpingFactorKMT09</name>
       !@   <defaultValue>5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The density enhancement (relative to mean disk density) for molecular complexes in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('molecularComplexClumpingFactorKMT09',molecularComplexClumpingFactorKMT09,defaultValue=5.0d0)
       !@ <inputParameter>
       !@   <name>molecularFractionFastKMT09</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@      Selects whether the fast (but less accurate) fitting formula for molecular hydrogen should be used in the ``Krumholz-McKee-Tumlinson'' star formation timescale calculation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('molecularFractionFastKMT09',molecularFractionFastKMT09,defaultValue=.false.)

       ! Set a pointer to the molecular hydrogen fraction fitting function to be used.
       select case (molecularFractionFastKMT09)
       case (.true.)
          KMT09_Molecular_Fraction => KMT09_Molecular_Fraction_Fast
       case(.false.)
          KMT09_Molecular_Fraction => KMT09_Molecular_Fraction_Slow
       end select

       ! Get the number of abundance properties.
       abundancesCount=Abundances_Property_Count()
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_KMT09_Initialize

  double precision function Star_Formation_Timescale_Disk_KMT09(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt thisNode}. The disk is assumed to obey the
    !% \cite{krumholz_star_2009} star formation rule.
    use Tree_Nodes
    use Numerical_Constants_Math
    use Numerical_Constants_Prefixes
    use Abundances_Structure
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use FGSL
    implicit none
    type(treeNode),            intent(inout), pointer     :: thisNode
    type(abundancesStructure), save                       :: fuelAbundances
    !$omp threadprivate(fuelAbundances)
    double precision,          dimension(abundancesCount) :: abundanceMasses
    ! Inner and outer radii (in units of disk scale length) between which the star formation surface density rate is
    ! integrated. The outer radius should be large enough that the exponential decline in surface density implies little star
    ! formation is missed at larger radius.
    double precision,          parameter                  :: radiusInnerDimensionless=0.0d0,radiusOuterDimensionless=10.0d0
    double precision                                      :: gasMass,diskScaleRadius,starFormationRate,radiusInner,radiusOuter
    type(c_ptr)                                           :: parameterPointer
    type(fgsl_function)                                   :: integrandFunction
    type(fgsl_integration_workspace)                      :: integrationWorkspace
  
    ! Get the disk properties.
    gasMass        =Tree_Node_Disk_Gas_Mass(thisNode)
    diskScaleRadius=Tree_Node_Disk_Radius  (thisNode)

    ! Check if the disk is physical.
    if (gasMass <= 0.0d0 .or. diskScaleRadius <= 0.0d0) then
       ! It is not, so return zero timescale.
       Star_Formation_Timescale_Disk_KMT09=0.0d0
    else
       ! Find the hydrogen fraction in the disk gas.
       call Tree_Node_Disk_Gas_Abundances(thisNode,abundanceMasses)
       call fuelAbundances%pack(abundanceMasses)
       call fuelAbundances%massToMassFraction(gasMass)
       hydrogenMassFraction      =fuelAbundances%hydrogenMassFraction()
       metallicityRelativeToSolar=fuelAbundances%metallicity(linearByMassSolar)

       ! Set a pointer to the node that is accessible by integral function.
       activeNode => thisNode

       ! Compute suitable limits for the integration.
       radiusInner=diskScaleRadius*radiusInnerDimensionless
       radiusOuter=diskScaleRadius*radiusOuterDimensionless

       ! Precompute constant factors appearing in the integrand.
       if (metallicityRelativeToSolar > 0.0d0) then
          chi                               =0.77d0*(1.0d0+3.1d0*metallicityRelativeToSolar**0.365d0)
          sigmaMolecularComplexNormalization=hydrogenMassFraction*molecularComplexClumpingFactorKMT09/mega**2
          sNormalization                    =dlog(1.0d0+0.6d0*chi+0.01d0*chi**2)/(0.04d0*metallicityRelativeToSolar)
       end if

       ! Compute the star formation rate. A low order integration rule works best here as the integrand can be discontinuous.
       starFormationRate=2.0d0*Pi*hydrogenMassFraction*starFormationFrequencyKMT09*Integrate(radiusInner,radiusOuter &
            &,Star_Formation_Rate_Integrand_KMT09,parameterPointer ,integrandFunction,integrationWorkspace,toleranceAbsolute &
            &=0.0d0,toleranceRelative=1.0d-3,integrationRule=FGSL_Integ_Gauss15)
       call Integrate_Done(integrandFunction,integrationWorkspace)

       ! Infer the star formation timescale.
       if (starFormationRate > 0.0d0) then
          Star_Formation_Timescale_Disk_KMT09=gasMass/starFormationRate
       else
          Star_Formation_Timescale_Disk_KMT09=0.0d0
       end if
    end if
    return
  end function Star_Formation_Timescale_Disk_KMT09

  function Star_Formation_Rate_Integrand_KMT09(radius,parameterPointer) bind(c)
    !% Integrand function for the ``Krumholz-McKee-Tumlinson'' star formation rate calculation.
    use, intrinsic :: ISO_C_Binding
    use Galactic_Structure_Surface_Densities
    use Galactic_Structure_Options
    use Numerical_Constants_Prefixes
    implicit none
    real(c_double)              :: Star_Formation_Rate_Integrand_KMT09
    real(c_double),   value     :: radius
    type(c_ptr),      value     :: parameterPointer
    double precision, parameter :: surfaceDensityTransition=85.0d12 ! M_Solar/Mpc^2
    double precision            :: molecularFraction,cloudFactor,surfaceDensityGasDimensionless,exponentialFactor&
         &,sigmaMolecularComplex,s,surfaceDensityGas

    ! Get gas surface density.
    surfaceDensityGas=Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,massType=massTypeGaseous,componentType=componentTypeDisk)

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

    ! Compute the star formation rate integrand.
    Star_Formation_Rate_Integrand_KMT09=radius*surfaceDensityGas*cloudFactor*molecularFraction
    return
  end function Star_Formation_Rate_Integrand_KMT09

  double precision function KMT09_Molecular_Fraction_Slow(s)
    !% Slow (but more accurate at low molecular fraction) fitting function from \cite{krumholz_star_2009} for the molecular hydrogen fraction.
    implicit none
    double precision, intent(in) :: s
    double precision, parameter  :: sMaximum=8.0d0
    double precision             :: delta
    
    ! Check if s is below maximum. If not, simply truncate to the minimum fraction that we allow.
    if (s < sMaximum) then
       delta                        =0.0712d0/((0.1d0/s+0.675d0)**2.8d0)
       KMT09_Molecular_Fraction_Slow=max(1.0d0-1.0d0/((1.0d0+(((1.0d0+delta)/0.75d0/s)**5))**0.2d0),molecularFractionMinimum)
    else
       KMT09_Molecular_Fraction_Slow=                                                               molecularFractionMinimum
    end if
    return
  end function KMT09_Molecular_Fraction_Slow

  double precision function KMT09_Molecular_Fraction_Fast(s)
    !% Slow (but less accurate at low molecular fraction) fitting function from \cite{mckee_atomic--molecular_2010} for the molecular hydrogen fraction.
    implicit none
    double precision, intent(in) :: s

    ! Check that s is below 2 - if it is, compute the molecular fraction, otherwise truncate to the minimum.
    if (s < 2.0d0) then
       KMT09_Molecular_Fraction_Fast=max(1.0d0-0.75d0*s/(1.0d0+0.25d0*s),molecularFractionMinimum)
    else
       KMT09_Molecular_Fraction_Fast=                                    molecularFractionMinimum
    end if
    return
  end function KMT09_Molecular_Fraction_Fast

end module Star_Formation_Timescale_Disks_KMT09
