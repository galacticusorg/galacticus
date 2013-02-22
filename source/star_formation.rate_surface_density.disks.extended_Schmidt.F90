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
  integer         (kind=kind_int8) :: lastUniqueID=-1
  !$omp threadprivate(lastUniqueID)

  ! Record of whether or not factors have been precomputed.
  logical                          :: factorsComputed=.false.
  !$omp threadprivate(factorsComputed)

  ! Precomputed factors.
  double precision                 :: hydrogenMassFraction
  !$omp threadprivate(hydrogenMassFraction)
  
  ! Parameters of the model.
  double precision                 :: starFormationExtendedSchmidtNormalization,starFormationExtendedSchmidtGasExponent&
       &,starFormationExtendedSchmidtStarExponent

contains

  !# <calculationResetTask>
  !# <unitName>Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Reset</unitName>
  !# </calculationResetTask>
  subroutine Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Reset(thisNode)
    !% Reset the extended Schmidt relation calculation.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    factorsComputed=.false.
    lastUniqueID   =thisNode%uniqueID()
    return
  end subroutine Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Reset

  !# <starFormationRateSurfaceDensityDisksMethod>
  !#  <unitName>Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Initialize</unitName>
  !# </starFormationRateSurfaceDensityDisksMethod>
  subroutine Star_Formation_Rate_Surface_Density_Disks_ExSchmidt_Initialize(starFormationRateSurfaceDensityDisksMethod&
       &,Star_Formation_Rate_Surface_Density_Disk_Get)
    !% Initializes the ``extended Schmidt'' disk star formation rate surface density.
    use ISO_Varying_String
    use Input_Parameters
    use Abundances_Structure
    use Numerical_Constants_Prefixes
    implicit none
    type     (varying_string  ),          intent(in   ) :: starFormationRateSurfaceDensityDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Rate_Surface_Density_Disk_Get
    
    if (starFormationRateSurfaceDensityDisksMethod == 'extendedSchmidt') then
       Star_Formation_Rate_Surface_Density_Disk_Get => Star_Formation_Rate_Surface_Density_Disk_ExSchmidt
    ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationExtendedSchmidtNormalization</name>
       !@   <defaultValue>$10^{-10.28}$ \citep{shi_extended_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The normalization of the extended Schmidt star formation law [$M_\odot$ yr$^{-1}$pc$^{-2}$].
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationExtendedSchmidtNormalization',starFormationExtendedSchmidtNormalization,defaultValue=0.5248d-10)
       !@ <inputParameter>
       !@   <name>starFormationExtendedSchmidtGasExponent</name>
       !@   <defaultValue>$1.0$ \citep{shi_extended_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of gas surface density in the extended Schmidt star formation law.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationExtendedSchmidtGasExponent'  ,starFormationExtendedSchmidtGasExponent  ,defaultValue=1.0000d+0 )
       !@ <inputParameter>
       !@   <name>starFormationExtendedSchmidtStarExponent</name>
       !@   <defaultValue>$0.48$ \citep{shi_extended_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of stellar surface density in the extended Schmidt star formation law.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationExtendedSchmidtStarExponent' ,starFormationExtendedSchmidtStarExponent ,defaultValue=0.4800d+0 )
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
    !% in the galactic disk of {\tt thisNode}. The disk is assumed to obey the extended Schmidt law of \cite{shi_extended_2011}:
    !% \begin{equation}
    !% \dot{\Sigma}_\star = A \left(x_{\rm H} {\Sigma_{\rm gas}\over M_\odot \hbox{pc}^{-2}}\right)
    !% ^{N_1} \left({\Sigma_{\star}\over M_\odot \hbox{pc}^{-2}}\right)^{N_2},
    !% \end{equation}
    !% where $A=${\tt [starFormationExtendedSchmidtNormalization]} and $N_1=${\tt
    !% [starFormationExtendedSchmidtGasExponent]}. $N_2=${\tt [starFormationExtendedSchmidtStarExponent]}.
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
    double precision                                            :: gasMass,diskScaleRadius,surfaceDensityGas,surfaceDensityStar

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

end module Star_Formation_Rate_Surface_Density_Disks_ExSchmidt
