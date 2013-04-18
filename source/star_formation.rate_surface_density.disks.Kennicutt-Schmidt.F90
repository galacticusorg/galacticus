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
  integer         (kind=kind_int8) :: lastUniqueID=-1
  !$omp threadprivate(lastUniqueID)

  ! Record of whether or not factors have been precomputed.
  logical                          :: factorsComputed=.false.
  !$omp threadprivate(factorsComputed)

  ! Precomputed factors.
  double precision                 :: criticalDensityFactor,hydrogenMassFraction
  !$omp threadprivate(criticalDensityFactor,hydrogenMassFraction)
  
  ! Parameters of the model.
  double precision                 :: starFormationKennicuttSchmidtNormalization,starFormationKennicuttSchmidtExponent&
       &,velocityDispersionDiskGas ,toomreParameterCritical,starFormationKennicuttSchmidtExponentTruncated
  logical                          :: starFormationKennicuttSchmidtTruncate

contains

  !# <calculationResetTask>
  !# <unitName>Star_Formation_Rate_Surface_Density_Disks_KS_Reset</unitName>
  !# </calculationResetTask>
  subroutine Star_Formation_Rate_Surface_Density_Disks_KS_Reset(thisNode)
    !% Reset the Kennicutt-Schmidt relation calculation.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    factorsComputed=.false.
    lastUniqueID   =thisNode%uniqueID()
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_KS_Reset

  !# <starFormationRateSurfaceDensityDisksMethod>
  !#  <unitName>Star_Formation_Rate_Surface_Density_Disks_KS_Initialize</unitName>
  !# </starFormationRateSurfaceDensityDisksMethod>
  subroutine Star_Formation_Rate_Surface_Density_Disks_KS_Initialize(starFormationRateSurfaceDensityDisksMethod&
       &,Star_Formation_Rate_Surface_Density_Disk_Get)
    !% Initializes the ``Kennicutt-Schmidt'' disk star formation rate surface density.
    use ISO_Varying_String
    use Input_Parameters
    use Abundances_Structure
    use Numerical_Constants_Prefixes
    implicit none
    type     (varying_string  ),          intent(in   ) :: starFormationRateSurfaceDensityDisksMethod
    procedure(Star_Formation_Rate_Surface_Density_Disk_KS), pointer, intent(inout) :: Star_Formation_Rate_Surface_Density_Disk_Get
    
    if (starFormationRateSurfaceDensityDisksMethod == 'Kennicutt-Schmidt') then
       Star_Formation_Rate_Surface_Density_Disk_Get => Star_Formation_Rate_Surface_Density_Disk_KS
       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationKennicuttSchmidtNormalization</name>
       !@   <defaultValue>$0.147$ \citep{kennicutt_global_1998}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The normalization of the Kennicutt-Schmidt star formation law [$M_\odot$ Gyr$^{-1}$pc$^{-2}$].
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationKennicuttSchmidtNormalization'    ,starFormationKennicuttSchmidtNormalization    ,defaultValue=0.147d0)
       !@ <inputParameter>
       !@   <name>starFormationKennicuttSchmidtExponent</name>
       !@   <defaultValue>$1.4$ \citep{kennicutt_global_1998}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent in the Kennicutt-Schmidt star formation law.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationKennicuttSchmidtExponent'         ,starFormationKennicuttSchmidtExponent         ,defaultValue=1.400d0)
       !@ <inputParameter>
       !@   <name>starFormationKennicuttSchmidtTruncate</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to truncate star formation below a critical surface density in disks.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationKennicuttSchmidtTruncate'         ,starFormationKennicuttSchmidtTruncate         ,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>starFormationKennicuttSchmidtExponentTruncated</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of the $\Sigma_{\rm gas}/\Sigma_{\rm crit}$ term used in truncating the Kennicutt-Schmidt star formation law.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationKennicuttSchmidtExponentTruncated',starFormationKennicuttSchmidtExponentTruncated,defaultValue= 6.0d0)
       !@ <inputParameter>
       !@   <name>velocityDispersionDiskGas</name>
       !@   <defaultValue>10 \citep{leroy_star_2008}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity dispersion of gas in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('velocityDispersionDiskGas'                     ,velocityDispersionDiskGas                     ,defaultValue=10.0d0)
       !@ <inputParameter>
       !@   <name>toomreParameterCritical</name>
       !@   <defaultValue>0.4 \citep{kennicutt_star_1989}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The critical Toomre parameter for star formation in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('toomreParameterCritical'                       ,toomreParameterCritical                       ,defaultValue= 0.4d0)
       ! Renormalize the Kennicutt-Schmidt relation to our internal units.
       starFormationKennicuttSchmidtNormalization=                       &
            &  starFormationKennicuttSchmidtNormalization                &
            & *mega**(2.0d0-2.0d0*starFormationKennicuttSchmidtExponent)
    end if
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_KS_Initialize

  double precision function Star_Formation_Rate_Surface_Density_Disk_KS(thisNode,radius)
    !% Returns the star formation rate surface density  (in $M_\odot$ Gyr$^{-1}$ Mpc$^{-2}$) for star formation in the galactic disk of {\tt thisNode}. The disk is assumed to obey the Kennicutt-Schmidt law:
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
    class           (nodeComponentDisk),                pointer :: thisDiskComponent
    type            (abundances       ), save                   :: fuelAbundances
    !$omp threadprivate(fuelAbundances)
    double precision                                            :: gasMass,criticalDensity,surfaceDensityGas

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
       criticalDensityFactor=toomreParameterCritical*dsqrt(2.0d0)&
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

end module Star_Formation_Rate_Surface_Density_Disks_KS
