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


!% Contains a module which implements a calculation of the core radius in the hot halo density profile that is a fraction of the
!% dark matter profile scale radius, but which grows as gas is depleted to keep the density at the virial radius constant (similar
!% to the algorithm of \citep{cole_hierarchical_2000}).

module Hot_Halo_Density_Cored_Isothermal_Core_Radii_Growing_Core
  !% Implements a calculation of the core radius in the hot halo density profile that is a fraction of the
  !% dark matter profile scale radius, but which grows as gas is depleted to keep the density at the virial radius constant (similar
  !% to the algorithm of \citep{cole_hierarchical_2000}).
  use Tree_Nodes
  use FGSL
  private
  public :: Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize

  ! Parameters of the model.
  double precision                            :: isothermalCoreRadiusOverScaleRadius,isothermalCoreRadiusOverVirialRadiusMaximum
  
  ! Minimum and maximum radii (in units of the virial radius) allowed for cores.
  double precision                            :: coreRadiusMinimum,coreRadiusMaximum
  
  ! Tabulation of core radius vs. virial density relation.
  integer,          parameter                 :: coreRadiusTablePointsPerDecade=100
  integer                                     :: coreRadiusTableCount
  double precision, allocatable, dimension(:) :: coreRadiusTableCoreRadius,coreRadiusTableDensityFactor
  logical                                     :: coreRadiusTableInitialized=.false.

  ! Interpolator variables.
  type(fgsl_interp)                           :: interpolationObject
  type(fgsl_interp_accel)                     :: interpolationAccelerator
  logical                                     :: interpolationReset=.true.

contains

  !# <hotHaloCoredIsothermalCoreRadiiMethod>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize</unitName>
  !# </hotHaloCoredIsothermalCoreRadiiMethod>
  subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize(hotHaloCoredIsothermalCoreRadiiMethod&
       &,Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get)
    !% Initializes the ``growing core'' cored isothermal hot halo profile core radius module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: hotHaloCoredIsothermalCoreRadiiMethod
    procedure(),          pointer, intent(inout) :: Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get
    
    if (hotHaloCoredIsothermalCoreRadiiMethod == 'growing core') then
       Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get => Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core
       !@ <inputParameter>
       !@   <name>isothermalCoreRadiusOverScaleRadius</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The core radius in the ``cored isothermal'' hot halo density profile in units of the dark matter profile scale radius.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('isothermalCoreRadiusOverScaleRadius',isothermalCoreRadiusOverScaleRadius,defaultValue=0.1d0)
       !@ <inputParameter>
       !@   <name>isothermalCoreRadiusOverVirialRadiusMaximum</name>
       !@   <defaultValue>3.2</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum core radius in the ``cored isothermal'' hot halo density profile in units of the virial radius.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('isothermalCoreRadiusOverVirialRadiusMaximum',isothermalCoreRadiusOverVirialRadiusMaximum,defaultValue=3.2d0)

       ! Ensure that the dark matter profile supports the scale property.
       if (.not.associated(Tree_Node_Dark_Matter_Profile_Scale)) call Galacticus_Error_Report( &
            & 'Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize',&
            &'method requires a component that supports the Tree_Node_Dark_Matter_Profile_Scale property')
    end if
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radii_GC_Initialize

  double precision function Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core(thisNode)
    !% Returns the radius (in Mpc) of the core radius in a cored isothermal hot halo density profile. Assumes that the radius is
    !% a fraction of the dark matter profile scale radius, but which grows as gas is depleted to keep the density at the virial
    !% radius constant (similar to the algorithm of \citep{cole_hierarchical_2000}).
    use Cosmological_Parameters
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    use Memory_Management
    use Numerical_Ranges
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: isothermalCoreRadiusOverVirialRadius,isothermalCoreRadiusOverVirialRadiusInitial,hotGasFraction,targetValue
    logical                                  :: makeTable

    ! Find the fraction of gas in the hot halo relative to that expected from the universal baryon fraction.
    hotGasFraction=(Tree_Node_Hot_Halo_Mass(thisNode)/Tree_Node_Mass(thisNode))*(Omega_0()/Omega_b())

    ! Return an arbitrary value for empty halos.
    if (hotGasFraction <= 0.0d0) then
       Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core=isothermalCoreRadiusOverScaleRadius
       return
    end if

    ! Comptue the desired core radius (in units of the virial radius) for a fully populated halo.
    isothermalCoreRadiusOverVirialRadiusInitial=isothermalCoreRadiusOverScaleRadius*Tree_Node_Dark_Matter_Profile_Scale(thisNode)&
         &/Dark_Matter_Halo_Virial_Radius(thisNode)

    ! Create a tabulation of core radius vs. virial density factor if necessary.
    !$omp critical (Hot_Halo_Density_Growing_Core_Interpolation)
    makeTable=(.not.coreRadiusTableInitialized) 
    if (.not.makeTable) makeTable=(isothermalCoreRadiusOverVirialRadiusInitial < coreRadiusTableCoreRadius(coreRadiusTableCount))
    if (makeTable) then
       coreRadiusMinimum   =min(isothermalCoreRadiusOverScaleRadius,isothermalCoreRadiusOverVirialRadiusInitial)
       coreRadiusMaximum   =isothermalCoreRadiusOverVirialRadiusMaximum
       coreRadiusTableCount=int(dlog10(coreRadiusMaximum/coreRadiusMinimum)*dble(coreRadiusTablePointsPerDecade))+1
       if (allocated(coreRadiusTableCoreRadius   )) call Dealloc_Array(coreRadiusTableCoreRadius   )
       if (allocated(coreRadiusTableDensityFactor)) call Dealloc_Array(coreRadiusTableDensityFactor)
       call Alloc_Array(coreRadiusTableCoreRadius   ,[coreRadiusTableCount])
       call Alloc_Array(coreRadiusTableDensityFactor,[coreRadiusTableCount])
       coreRadiusTableCoreRadius   =Make_Range(coreRadiusMaximum,coreRadiusMinimum,coreRadiusTableCount,rangeType=rangeTypeLogarithmic)
       coreRadiusTableDensityFactor=Growing_Core_Virial_Density_Function(coreRadiusTableCoreRadius)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,interpolationReset)
       interpolationReset=.true.
       coreRadiusTableInitialized=.true.
    end if
    !$omp end critical (Hot_Halo_Density_Growing_Core_Interpolation)

    ! Compute the target value of the function giving the density at the virial radius per unit gas mass.
    targetValue=Growing_Core_Virial_Density_Function(isothermalCoreRadiusOverVirialRadiusInitial)*hotGasFraction

    ! Interpolate to get the required core radius.
    if      (hotGasFraction >= 1.0d0                                             ) then
       isothermalCoreRadiusOverVirialRadius=isothermalCoreRadiusOverVirialRadiusInitial
    else if (targetValue    >= coreRadiusTableDensityFactor(coreRadiusTableCount)) then
       isothermalCoreRadiusOverVirialRadius=coreRadiusTableCoreRadius(coreRadiusTableCount)
    else
       !$omp critical (Hot_Halo_Density_Growing_Core_Interpolation)
       isothermalCoreRadiusOverVirialRadius=Interpolate(coreRadiusTableCount,coreRadiusTableDensityFactor,coreRadiusTableCoreRadius&
            &,interpolationObject,interpolationAccelerator,targetValue,reset=interpolationReset)
       !$omp end critical (Hot_Halo_Density_Growing_Core_Interpolation)
    end if
    
    ! Compute the resulting core radius.
    Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core=isothermalCoreRadiusOverVirialRadius*Dark_Matter_Halo_Virial_Radius(thisNode)

    return
  end function Hot_Halo_Density_Cored_Isothermal_Core_Radius_Growing_Core
  
  elemental double precision function Growing_Core_Virial_Density_Function(radiusOverVirialRadius)
    !% Returns the function $(1+r_{\rm c}^2)[1-r_{\rm c} \tan^{-1}(1/r_{\rm c}]$ which is proportional to the density at the
    !% virial radius of a cored isothermal profile with core radius $r_{\rm c}$ (in units of the virial radius) per unit mass.
    implicit none
    double precision, intent(in) :: radiusOverVirialRadius
    
    Growing_Core_Virial_Density_Function=(1.0d0+radiusOverVirialRadius**2)*(1.0d0-radiusOverVirialRadius*datan(1.0d0/radiusOverVirialRadius))
    return
  end function Growing_Core_Virial_Density_Function
  
end module Hot_Halo_Density_Cored_Isothermal_Core_Radii_Growing_Core
