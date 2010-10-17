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


!% Contains a module which implements the Kennicutt-Schmidt star formation timescale for galactic disks.

module Star_Formation_Timescale_Disks_Kennicutt_Schmidt
  !% Implements the Kennicutt-Schmidt star formation timescale for galactic disks.
  use Tree_Nodes
  private
  public :: Star_Formation_Timescale_Disks_Kennicutt_Schmidt_Initialize

  ! Internal copy of the number of abundances properties.
  integer          :: abundancesCount

  ! Parameters of the model.
  double precision :: starFormationKennicuttSchmidtNormalization,starFormationKennicuttSchmidtExponent,velocityDispersionDiskGas&
       &,toomreParameterCritical,starFormationKennicuttSchmidtExponentTruncated
  logical          :: starFormationKennicuttSchmidtTruncate

  ! Module global parameter used in root finding.
  double precision :: criticalDensityFactor
  !$omp threadprivate(criticalDensityFactor)

  ! Pointer to active node used in integral functions, plus variables needed by integral function.
  type(treeNode),   pointer :: activeNode
  !$omp threadprivate(activeNode)

contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Kennicutt_Schmidt_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Kennicutt_Schmidt_Initialize(starFormationTimescaleDisksMethod&
       &,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``Kennicutt-Schmidt'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    use Abundances_Structure
    implicit none
    type(varying_string),          intent(in)    :: starFormationTimescaleDisksMethod
    procedure(),          pointer, intent(inout) :: Star_Formation_Timescale_Disk_Get
    
    if (starFormationTimescaleDisksMethod == 'Kennicutt-Schmidt') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Kennicutt_Schmidt
       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationKennicuttSchmidtNormalization</name>
       !@   <defaultValue>$0.147$ \citep{kennicutt_global_1998}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The normalization of the Kennicutt-Schmidt star formation law [$M_\odot$ Gyr$^{-1}$pc$^{-2}$].
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationKennicuttSchmidtNormalization'    ,starFormationKennicuttSchmidtNormalization    ,defaultValue=0.147d0)
       !@ <inputParameter>
       !@   <name>starFormationKennicuttSchmidtExponent</name>
       !@   <defaultValue>$1.4$ \citep{kennicutt_global_1998}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent in the Kennicutt-Schmidt star formation law.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationKennicuttSchmidtExponent'         ,starFormationKennicuttSchmidtExponent         ,defaultValue=1.400d0)
       !@ <inputParameter>
       !@   <name>starFormationKennicuttSchmidtTruncate</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to truncate star formation below a critical surface density in disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationKennicuttSchmidtTruncate'         ,starFormationKennicuttSchmidtTruncate         ,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>starFormationKennicuttSchmidtExponentTruncated</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of the $\Sigma_{\rm gas}/\Sigma_{\rm crit}$ term used in truncating the Kennicutt-Schmidt star formation law.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationKennicuttSchmidtExponentTruncated',starFormationKennicuttSchmidtExponentTruncated,defaultValue= 6.0d0)
       !@ <inputParameter>
       !@   <name>velocityDispersionDiskGas</name>
       !@   <defaultValue>10 \citep{leroy_star_2008}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity dispersion of gas in disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('velocityDispersionDiskGas'                     ,velocityDispersionDiskGas                     ,defaultValue=10.0d0)
       !@ <inputParameter>
       !@   <name>toomreParameterCritical</name>
       !@   <defaultValue>0.4 \citep{kennicutt_star_1989}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The critical Toomre parameter for star formation in disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('toomreParameterCritical'                       ,toomreParameterCritical                       ,defaultValue= 0.4d0)

       ! Get the number of abundance properties.
       abundancesCount=Abundances_Property_Count()
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Kennicutt_Schmidt_Initialize

  double precision function Star_Formation_Timescale_Disk_Kennicutt_Schmidt(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt thisNode}. The disk is assumed to obey the Kennicutt-Schmidt law:
    !% \begin{equation}
    !% \Sigma_\star = A \left(x_{\rm H} {\Sigma_{\rm gas}\over M_\odot \hbox{pc}^{-2}}\right)^N,
    !% \end{equation}
    !% where $A=${\tt [starFormationKennicuttSchmidtNormalization]} and $N=${\tt
    !% [starFormationKennicuttSchmidtExponent]}. Optionally, star formation is truncated for gas surface densities below a critical density of:
    !% \begin{equation}
    !% \Sigma_{\rm crit} = {q_{\rm crit} \kappa \sigma_{\rm gas} \over \pi \G},
    !% \end{equation}
    !% where $\kappa$ is the epicyclic frequency in the disk, $\sigma_{\rm gas}$ is the velocity dispersion of gas in the disk and
    !% $q_{\rm crit}=${\tt [toomreParameterCritical]} is a dimensionless constant of order unity which controls where the critical
    !% density occurs. $\sigma_{\rm gas}$ is assumed to be a constant equal to {\tt [velocityDispersionDiskGas]} and the disk is
    !% assumed to have a flat rotation curve such that $\kappa = \sqrt{2} V/R$.
    use Tree_Nodes
    use Numerical_Constants_Math
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Physical
    use Abundances_Structure
    use FGSL
    use Numerical_Integration
    use, intrinsic :: ISO_C_Binding
    implicit none
    type(treeNode),            intent(inout), pointer     :: thisNode
    type(abundancesStructure), save                       :: fuelAbundances
    !$omp threadprivate(fuelAbundances)
    double precision,          dimension(abundancesCount) :: abundanceMasses
    double precision,          parameter                  :: radiusInnerDimensionless=0.0d0,radiusOuterDimensionless=10.0d0
    double precision                                      :: gasMass,diskScaleRadius,starFormationRate,radiusInner,radiusOuter&
         &,hydrogenMassFraction
    type(c_ptr)                                           :: parameterPointer
    type(fgsl_function)                                   :: integrandFunction
    type(fgsl_integration_workspace)                      :: integrationWorkspace

    ! Get the disk properties.
    gasMass             =Tree_Node_Disk_Gas_Mass(thisNode)
    diskScaleRadius     =Tree_Node_Disk_Radius  (thisNode)

    ! Check if the disk is physical.
    if (gasMass <= 0.0d0 .or. diskScaleRadius <= 0.0d0) then
       ! It is not, so return zero timescale.
       Star_Formation_Timescale_Disk_Kennicutt_Schmidt=0.0d0
    else
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       call Tree_Node_Disk_Gas_Abundances(thisNode,abundanceMasses)
       call fuelAbundances%pack(abundanceMasses)
       call fuelAbundances%massToMassFraction(gasMass)
       hydrogenMassFraction=fuelAbundances%hydrogenMassFraction()

       ! Set a pointer to the node that is accessible by integral function.
       activeNode => thisNode

       ! Compute critical surface density factor if necessary.
       if (starFormationKennicuttSchmidtTruncate) criticalDensityFactor=toomreParameterCritical*dsqrt(2.0d0)&
            &*velocityDispersionDiskGas*Tree_Node_Disk_Velocity(activeNode)/Pi/gravitationalConstantGalacticus

       ! Compute suitable limits for the integration.
       radiusInner=diskScaleRadius*radiusInnerDimensionless
       radiusOuter=diskScaleRadius*radiusOuterDimensionless

       ! Compute the star formation rate. A low order integration rule (FGSL_Integ_Gauss15) works well here, particularly when a
       ! truncation surface density is used (since that truncation introduces a discontinuity in the integrand).
       starFormationRate=2.0d0*Pi*starFormationKennicuttSchmidtNormalization*(mega**2)*((hydrogenMassFraction/mega**2)&
            &**starFormationKennicuttSchmidtExponent)*Integrate(radiusInner,radiusOuter,Star_Formation_Rate_Integrand_KS&
            &,parameterPointer ,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3&
            &,integrationRule=FGSL_Integ_Gauss15)
       call Integrate_Done(integrandFunction,integrationWorkspace)

       ! Compute the star formation timescale.
       if (starFormationRate > 0.0d0) then
          Star_Formation_Timescale_Disk_Kennicutt_Schmidt=gasMass/starFormationRate
       else
          Star_Formation_Timescale_Disk_Kennicutt_Schmidt=0.0d0
       end if

    end if
    return
  end function Star_Formation_Timescale_Disk_Kennicutt_Schmidt
  
  function Star_Formation_Rate_Integrand_KS(radius,parameterPointer) bind(c)
    !% Integrand function for the ``Kennicutt-Schmidt'' star formation rate calculation.
    use Galactic_Structure_Surface_Densities
    use Galactic_Structure_Options
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    implicit none
    real(c_double)          :: Star_Formation_Rate_Integrand_KS
    real(c_double),   value :: radius
    type(c_ptr),      value :: parameterPointer
    double precision        :: surfaceDensityGas,criticalDensity

    ! Get gas surface density.
    surfaceDensityGas=Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,massType=massTypeGaseous,componentType=componentTypeDisk)

    ! Compute the star formation rate integrand.
    Star_Formation_Rate_Integrand_KS=radius*(surfaceDensityGas**starFormationKennicuttSchmidtExponent)

    ! Check if we are applying a truncation radius.
    if (starFormationKennicuttSchmidtTruncate) then

       ! Always return zero star formation rate at zero radius, as critical density will be infinite.
       if (radius <= 0.0d0) then
          Star_Formation_Rate_Integrand_KS=0.0d0
          return
       end if

       ! Compute the critical density for star formation.
       criticalDensity=criticalDensityFactor/radius
       
       ! Check if gas is above the critical density. Return zero star formation rate if it is not.
       if (surfaceDensityGas < criticalDensity) Star_Formation_Rate_Integrand_KS=Star_Formation_Rate_Integrand_KS&
            &*(surfaceDensityGas/criticalDensity)**starFormationKennicuttSchmidtExponentTruncated

    end if

    return
  end function Star_Formation_Rate_Integrand_KS

end module Star_Formation_Timescale_Disks_Kennicutt_Schmidt
