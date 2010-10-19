!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements the \cite{blitz_role_2006} star formation timescale for galactic disks.

module Star_Formation_Timescale_Disks_Blitz_Rosolowsky
  !% Implements the \cite{blitz_role_2006} star formation timescale for galactic disks.
  use Tree_Nodes
  private
  public :: Star_Formation_Timescale_Disks_Blitz_Rosolowsky_Initialize

  ! Internal copy of the number of abundances properties.
  integer          :: abundancesCount

  ! Parameters of the model.
  double precision :: velocityDispersionDiskGas,heightToRadialScaleDisk,surfaceDensityCriticalBlitzRosolowsky&
       &,surfaceDensityExponentBlitzRosolowsky,starFormationFrequencyNormalizationBlitzRosolowsky&
       &,pressureCharacteristicBlitzRosolowsky,pressureExponentBlitzRosolowsky

  ! Pointer to active node used in integral functions, plus variables needed by integral function.
  type(treeNode),   pointer :: activeNode
  double precision          :: hydrogenMassFraction,diskScaleRadius
  !$omp threadprivate(activeNode,hydrogenMassFraction,diskScaleRadius)

contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Blitz_Rosolowsky_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Blitz_Rosolowsky_Initialize(starFormationTimescaleDisksMethod&
       &,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``Blitz-Rosolowsky'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use Abundances_Structure
    use Galacticus_Error
    implicit none
    type(varying_string),                 intent(in)    :: starFormationTimescaleDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Timescale_Disk_Get
    
    if (starFormationTimescaleDisksMethod == 'Blitz-Rosolowsky') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Blitz_Rosolowsky

       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>velocityDispersionDiskGas</name>
       !@   <defaultValue>10 \citep{leroy_star_2008}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity dispersion of gas in disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('velocityDispersionDiskGas',velocityDispersionDiskGas,defaultValue=10.0d0)
       !@ <inputParameter>
       !@   <name>heightToRadialScaleDisk</name>
       !@   <defaultValue>0.137 \citep{kregel_flattening_2002}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The ratio of scale height to scale radius for disks in the ``Blitz-Rosolowsky'' star formation timescale calculation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('heightToRadialScaleDisk',heightToRadialScaleDisk,defaultValue=0.137d0)
       !@ <inputParameter>
       !@   <name>surfaceDensityCriticalBlitzRosolowsky</name>
       !@   <defaultValue>200 \citep{bigiel_star_2008}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The surface density (in units of $M_\odot$ pc$^{-2}$) in the ``Blitz-Rosolowsky'' star formation timescale calculation at which low-density truncation begins.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('surfaceDensityCriticalBlitzRosolowsky',surfaceDensityCriticalBlitzRosolowsky,defaultValue=200.0d0)
       !@ <inputParameter>
       !@   <name>surfaceDensityExponentBlitzRosolowsky</name>
       !@   <defaultValue>0.4 \citep{bigiel_star_2008}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent for surface density in the ``Blitz-Rosolowsky'' star formation timescale calculation at in the high density regime.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('surfaceDensityExponentBlitzRosolowsky',surfaceDensityExponentBlitzRosolowsky,defaultValue=0.4d0)
       !@ <inputParameter>
       !@   <name>starFormationFrequencyNormalizationBlitzRosolowsky</name>
       !@   <defaultValue>$5.25\times 10^{-10}$ \citep{leroy_star_2008}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The star formation frequency (in the low-density limit and in units of yr$^{-1}$) in the ``Blitz-Rosolowsky'' star formation timescale calculation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationFrequencyNormalizationBlitzRosolowsky',starFormationFrequencyNormalizationBlitzRosolowsky,defaultValue=5.25d-10)
       !@ <inputParameter>
       !@   <name>pressureCharacteristicBlitzRosolowsky</name>
       !@   <defaultValue>4.54 \citep{blitz_role_2006}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The characteristic pressure (given as $P_0/k_{\rm B}$ in units of K cm$^{-3}$) in the scaling relation of molecular hydrogen fraction with disk pressure in the ``Blitz-Rosolowsky'' star formation timescale calculation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('pressureCharacteristicBlitzRosolowsky',pressureCharacteristicBlitzRosolowsky,defaultValue=4.54d0)
       !@ <inputParameter>
       !@   <name>pressureExponentBlitzRosolowsky</name>
       !@   <defaultValue>0.92 \citep{blitz_role_2006}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent in the scaling relation of molecular hydrogen fraction with disk pressure in the ``Blitz-Rosolowsky'' star formation timescale calculation.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('pressureExponentBlitzRosolowsky',pressureExponentBlitzRosolowsky,defaultValue=0.92d0)
       if (pressureExponentBlitzRosolowsky < 0.0d0) call Galacticus_Error_Report('Star_Formation_Timescale_Disks_Blitz_Rosolowsky_Initialize','pressureExponentBlitzRosolowsky < 0 violates assumptions')

       ! Convert parameters to internal units.
       surfaceDensityCriticalBlitzRosolowsky             =surfaceDensityCriticalBlitzRosolowsky*(mega**2)                                                    ! Convert to M_Solar/Mpc^2.
       starFormationFrequencyNormalizationBlitzRosolowsky=starFormationFrequencyNormalizationBlitzRosolowsky*giga                                            ! Convert to Gyr^-1.
       pressureCharacteristicBlitzRosolowsky             =pressureCharacteristicBlitzRosolowsky*boltzmannsConstant*((hecto*megaParsec)**3)/massSolar/kilo**2 ! Convert to M_Solar (km/s)^2 / Mpc.

       ! Get the number of abundance properties.
       abundancesCount=Abundances_Property_Count()
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Blitz_Rosolowsky_Initialize

  double precision function Star_Formation_Timescale_Disk_Blitz_Rosolowsky(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt thisNode}. The disk is assumed to obey the
    !% \cite{blitz_role_2006} star formation rule.
    use Tree_Nodes
    use Numerical_Constants_Math
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Physical
    use Abundances_Structure
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    implicit none
    type(treeNode),            intent(inout), pointer     :: thisNode
    type(abundancesStructure), save                       :: fuelAbundances
    !$omp threadprivate(fuelAbundances)
    double precision,          dimension(abundancesCount) :: abundanceMasses
    ! Inner and outer radii (in units of disk scale length) between which the star formation surface density rate is
    ! integrated. The outer radius should be large enough that the decline in surface density implies little star
    ! formation is missed at larger radius.
    double precision,          parameter                  :: radiusInnerDimensionless=0.0d0,radiusOuterDimensionless=10.0d0
    double precision                                      :: gasMass,stellarMass,starFormationRate,radiusInner,radiusOuter
    type(c_ptr)                                           :: parameterPointer
    type(fgsl_function)                                   :: integrandFunction
    type(fgsl_integration_workspace)                      :: integrationWorkspace
  
    ! Get the disk properties.
    gasMass        =Tree_Node_Disk_Gas_Mass    (thisNode)
    stellarMass    =Tree_Node_Disk_Stellar_Mass(thisNode)
    diskScaleRadius=Tree_Node_Disk_Radius      (thisNode)

    ! Check if the disk is physical.
    if (gasMass <= 0.0d0 .or. stellarMass < 0.0d0 .or. diskScaleRadius <= 0.0d0) then
       ! It is not, so return zero timescale.
       Star_Formation_Timescale_Disk_Blitz_Rosolowsky=0.0d0
    else
       ! Find the hydrogen fraction in the disk gas.
       call Tree_Node_Disk_Gas_Abundances(thisNode,abundanceMasses)
       call fuelAbundances%pack(abundanceMasses)
       call fuelAbundances%massToMassFraction(gasMass)
       hydrogenMassFraction=fuelAbundances%hydrogenMassFraction()

       ! Set a pointer to the node that is accessible by integral function.
       activeNode => thisNode

       ! Compute suitable limits for the integration.
       radiusInner=diskScaleRadius*radiusInnerDimensionless
       radiusOuter=diskScaleRadius*radiusOuterDimensionless

       ! Compute the star formation rate. A low order integration rule works best here as the integrand can be discontinuous.
       starFormationRate=Integrate(radiusInner,radiusOuter,Star_Formation_Rate_Integrand_BR,parameterPointer ,integrandFunction&
            &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3,integrationRule=FGSL_Integ_Gauss15)
       call Integrate_Done(integrandFunction,integrationWorkspace)

       ! Infer the star formation timescale.
       Star_Formation_Timescale_Disk_Blitz_Rosolowsky=gasMass/starFormationRate

    end if
    return
  end function Star_Formation_Timescale_Disk_Blitz_Rosolowsky

  function Star_Formation_Rate_Integrand_BR(radius,parameterPointer) bind(c)
    !% Integrand function for the ``Blitz-Rosolowsky'' star formation rate calculation.
    use Galactic_Structure_Surface_Densities
    use Galactic_Structure_Options
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    implicit none
    real(c_double)          :: Star_Formation_Rate_Integrand_BR
    real(c_double),   value :: radius
    type(c_ptr),      value :: parameterPointer
    double precision        :: surfaceDensityGas,surfaceDensityStellar,pressureRatio,molecularFraction

    ! Get gas and stellar surface densities.
    surfaceDensityGas    =Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,massType=massTypeGaseous,componentType=componentTypeDisk)
    surfaceDensityStellar=Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,massType=massTypeStellar,componentType=componentTypeDisk)

    ! Compute the pressure ratio that Blitz & Rosolowsky use to compute the molecular fraction.
    pressureRatio=0.5d0*Pi*gravitationalConstantGalacticus*surfaceDensityGas*(surfaceDensityGas+velocityDispersionDiskGas&
         &*dsqrt(surfaceDensityStellar/Pi/gravitationalConstantGalacticus/heightToRadialScaleDisk/diskScaleRadius))&
         &/pressureCharacteristicBlitzRosolowsky

    ! Compute the molecular fraction, limited to 100% molecular.
    if (pressureRatio >= 1.0d0) then
       molecularFraction=1.0d0
    else
       molecularFraction=min(pressureRatio**pressureExponentBlitzRosolowsky,1.0d0)
    end if

    ! Compute the star formation rate integrand.
    Star_Formation_Rate_Integrand_BR=2.0d0*Pi*radius*surfaceDensityGas*hydrogenMassFraction*molecularFraction&
         &*starFormationFrequencyNormalizationBlitzRosolowsky*(1.0d0+(hydrogenMassFraction*surfaceDensityGas&
         &/surfaceDensityCriticalBlitzRosolowsky)**surfaceDensityExponentBlitzRosolowsky)

    return
  end function Star_Formation_Rate_Integrand_BR

end module Star_Formation_Timescale_Disks_Blitz_Rosolowsky
