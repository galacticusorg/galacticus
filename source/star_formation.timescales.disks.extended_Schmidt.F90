!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the extended Schmidt star formation timescale of \cite{shi_extended_2011} for galactic disks.

module Star_Formation_Timescale_Disks_Extended_Schmidt
  !% Implements the extended Schmidt star formation timescale of \cite{shi_extended_2011} for galactic disks.
  use Galacticus_Nodes
  implicit none
  private
  public :: Star_Formation_Timescale_Disks_Extended_Schmidt_Initialize

  ! Internal copy of the number of abundances properties.
  integer          :: abundancesCount

  ! Parameters of the model.
  double precision :: starFormationExtendedSchmidtNormalization,starFormationExtendedSchmidtGasExponent,starFormationExtendedSchmidtStarExponent

  ! Pointer to active node used in integral functions, plus variables needed by integral function.
  type(treeNode),   pointer :: activeNode
  !$omp threadprivate(activeNode)

contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Extended_Schmidt_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Extended_Schmidt_Initialize(starFormationTimescaleDisksMethod&
       &,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``extended Schmidt'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    use Abundances_Structure
    implicit none
    type(varying_string),                 intent(in)    :: starFormationTimescaleDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Timescale_Disk_Get

    if (starFormationTimescaleDisksMethod == 'extendedSchmidt') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Extended_Schmidt
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

       ! Get the number of abundance properties.
       abundancesCount=Abundances_Property_Count()

    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Extended_Schmidt_Initialize

  double precision function Star_Formation_Timescale_Disk_Extended_Schmidt(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt thisNode}. The disk is assumed to obey the extended Schmidt law of  \cite{shi_extended_2011}:
    !% \begin{equation}
    !% \dot{\Sigma}_\star = A \left(x_{\rm H} {\Sigma_{\rm gas}\over M_\odot \hbox{pc}^{-2}}\right)
    !% ^{N_1} \left({\Sigma_{\star}\over M_\odot \hbox{pc}^{-2}}\right)^{N_2},
    !% \end{equation}
    !% where $A=${\tt [starFormationExtendedSchmidtNormalization]} and $N_1=${\tt
    !% [starFormationExtendedSchmidtGasExponent]}. $N_2=${\tt [starFormationExtendedSchmidtStarExponent]}.
    use Galacticus_Nodes
    use Numerical_Constants_Math
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Physical
    use Abundances_Structure
    use FGSL
    use Numerical_Integration
    use, intrinsic :: ISO_C_Binding
    implicit none
    type (treeNode                  ), intent(inout), pointer :: thisNode
    class(nodeComponentDisk         ),                pointer :: thisDiskComponent
    type (abundances                ), save                   :: fuelAbundances
    !$omp threadprivate(fuelAbundances)
    double precision                 , parameter              :: radiusInnerDimensionless=0.0d0,radiusOuterDimensionless=10.0d0
    double precision                                          :: gasMass,diskScaleRadius,starFormationRate,radiusInner,radiusOuter&
         &,hydrogenMassFraction
    type (c_ptr                     )                         :: parameterPointer
    type (fgsl_function             )                         :: integrandFunction
    type (fgsl_integration_workspace)                         :: integrationWorkspace

    ! Get the disk properties.
    thisDiskComponent => thisNode%disk()
    gasMass           =  thisDiskComponent%massGas()
    diskScaleRadius   =  thisDiskComponent%radius ()

    ! Check if the disk is physical.
    if (gasMass <= 0.0d0 .or. diskScaleRadius <= 0.0d0) then
       ! It is not, so return zero timescale.
       Star_Formation_Timescale_Disk_Extended_Schmidt=0.0d0
    else
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       fuelAbundances=thisDiskComponent%abundancesGas()
       call fuelAbundances%massToMassFraction(gasMass)
       hydrogenMassFraction=fuelAbundances%hydrogenMassFraction()

       ! Set a pointer to the node that is accessible by integral function.
       activeNode => thisNode

       ! Compute suitable limits for the integration.
       radiusInner=diskScaleRadius*radiusInnerDimensionless
       radiusOuter=diskScaleRadius*radiusOuterDimensionless

       ! Compute the star formation rate.
       starFormationRate= 2.0d0*Pi                                                                   & ! Geometric factor.
            &            *(mega**2)*giga                                                             & ! Convert from Msun/pc^2/yr to Msun/Mpc^2/Gyr
            &            *starFormationExtendedSchmidtNormalization                                  & ! Normalization of the star formation rate.
            &            *((hydrogenMassFraction/mega**2)**starFormationExtendedSchmidtGasExponent ) & ! Hydrogen fraction and unit conversion for gas.
            &            *((1.0d0               /mega**2)**starFormationExtendedSchmidtStarExponent) & ! Unit conversion for stars.
            &            *Integrate(                                                                 & ! Integral over the disk.
            &                       radiusInner                                                     ,&
            &                       radiusOuter                                                     ,&
            &                       Star_Formation_Rate_Integrand_ES                                ,&
            &                       parameterPointer                                                ,&
            &                       integrandFunction                                               ,&
            &                       integrationWorkspace                                            ,&
            &                       toleranceAbsolute=0.0d0                                         ,&
            &                       toleranceRelative=1.0d-3                                        ,&
            &                       integrationRule  =FGSL_Integ_Gauss15                             &
            &                      )
       call Integrate_Done(integrandFunction,integrationWorkspace)

       ! Compute the star formation timescale.
       if (starFormationRate > 0.0d0) then
          Star_Formation_Timescale_Disk_Extended_Schmidt=gasMass/starFormationRate
       else
          Star_Formation_Timescale_Disk_Extended_Schmidt=0.0d0
       end if

    end if
    return
  end function Star_Formation_Timescale_Disk_Extended_Schmidt

  function Star_Formation_Rate_Integrand_ES(radius,parameterPointer) bind(c)
    !% Integrand function for the ``extended Schmidt'' star formation rate calculation.
    use Galactic_Structure_Surface_Densities
    use Galactic_Structure_Options
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    implicit none
    real(c_double)          :: Star_Formation_Rate_Integrand_ES
    real(c_double),   value :: radius
    type(c_ptr),      value :: parameterPointer
    double precision        :: surfaceDensityGas,surfaceDensityStar

    ! Get gas surface density.
    surfaceDensityGas=Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,massType=massTypeGaseous,componentType=componentTypeDisk)

    ! Get stellar surface density.
    surfaceDensityStar=Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystem&
         &=coordinateSystemCylindrical,massType=massTypeStellar,componentType=componentTypeDisk)

    ! Compute the star formation rate integrand.
    Star_Formation_Rate_Integrand_ES= radius&
         &*(surfaceDensityGas**starFormationExtendedSchmidtGasExponent)&
         &*(surfaceDensityStar**starFormationExtendedSchmidtStarExponent)

    if (radius <= 0.0d0) then
       Star_Formation_Rate_Integrand_ES=0.0d0
       return
    end if

    return
  end function Star_Formation_Rate_Integrand_ES

end module Star_Formation_Timescale_Disks_Extended_Schmidt
