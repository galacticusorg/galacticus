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


!% Contains a module which implements calculations of various scales for dark matter halos.

module Dark_Matter_Halo_Scales
  !% Implements calculations of various scales for dark matter halos.
  use Tree_Nodes
  use Kind_Numbers
  private
  public :: Dark_Matter_Halo_Dynamical_Timescale, Dark_Matter_Halo_Virial_Velocity, Dark_Matter_Halo_Virial_Velocity_Growth_Rate,&
       & Dark_Matter_Halo_Virial_Radius, Dark_Matter_Halo_Virial_Radius_Growth_Rate, Dark_Matter_Halo_Mean_Density,&
       & Dark_Matter_Halo_Mean_Density_Growth_Rate, Dark_Matter_Halo_Virial_Temperature, Dark_Matter_Halo_Scales_Reset

  ! Record of unique ID of node which we last computed results for.
  integer(kind=kind_int8) :: lastUniqueID=-1
  !$omp threadprivate(lastUniqueID)

  ! Record of whether or not halo scales have already been computed for this node.
  logical :: virialRadiusComputed=.false.,virialTemperatureComputed=.false.,virialVelocityComputed=.false.,dynamicalTimescaleComputed=.false.
  !$omp threadprivate(virialRadiusComputed,virialTemperatureComputed,virialVelocityComputed,dynamicalTimescaleComputed)

  ! Stored values of halo scales.
  double precision :: virialRadiusStored,virialTemperatureStored,virialVelocityStored,dynamicalTimescaleStored
  !$omp threadprivate(virialRadiusStored,virialTemperatureStored,virialVelocityStored,dynamicalTimescaleStored)

contains

  !# <calculationResetTask>
  !# <unitName>Dark_Matter_Halo_Scales_Reset</unitName>
  !# </calculationResetTask>
  subroutine Dark_Matter_Halo_Scales_Reset(thisNode)
    !% Reset the cooling radius calculation.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    virialRadiusComputed      =.false.
    virialTemperatureComputed =.false.
    virialVelocityComputed    =.false.
    dynamicalTimescaleComputed=.false.
    lastUniqueID              =thisNode%uniqueID()
    return
  end subroutine Dark_Matter_Halo_Scales_Reset

  double precision function Dark_Matter_Halo_Dynamical_Timescale(thisNode)
    !% Returns the dynamical timescale for {\tt thisNode}.
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Dark_Matter_Halo_Scales_Reset(thisNode)

    ! Check if halo dynamical timescale is already computed. Compute and store if not.
    if (.not.dynamicalTimescaleComputed) then
       dynamicalTimescaleComputed=.true.
       dynamicalTimescaleStored=Dark_Matter_Halo_Virial_Radius(thisNode)*(megaParsec/kilo/gigaYear) &
            &/Dark_Matter_Halo_Virial_Velocity(thisNode)
    end if
    
    ! Return the stored timescale.
    Dark_Matter_Halo_Dynamical_Timescale=dynamicalTimescaleStored
    return
  end function Dark_Matter_Halo_Dynamical_Timescale

  double precision function Dark_Matter_Halo_Virial_Velocity(thisNode)
    !% Returns the virial velocity scale for {\tt thisNode}.
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Dark_Matter_Halo_Scales_Reset(thisNode)

    ! Check if virial velocity is already computed. Compute and store if not.
    if (.not.virialVelocityComputed) then
       virialVelocityComputed=.true.
       virialVelocityStored=dsqrt(gravitationalConstantGalacticus*Tree_Node_Mass(thisNode)&
            &/Dark_Matter_Halo_Virial_Radius(thisNode))
    end if
    
    ! Return the stored virial velocity.
    Dark_Matter_Halo_Virial_Velocity=virialVelocityStored
    return
  end function Dark_Matter_Halo_Virial_Velocity

  double precision function Dark_Matter_Halo_Virial_Velocity_Growth_Rate(thisNode)
    !% Returns the growth rate of the virial velocity scale for {\tt thisNode}.
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode

    Dark_Matter_Halo_Virial_Velocity_Growth_Rate=0.5d0*Dark_Matter_Halo_Virial_Velocity(thisNode)&
         &*(Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)-Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)&
         &/Dark_Matter_Halo_Virial_Radius(thisNode))
    return
  end function Dark_Matter_Halo_Virial_Velocity_Growth_Rate

  double precision function Dark_Matter_Halo_Virial_Temperature(thisNode)
    !% Returns the virial temperature (in Kelvin) for {\tt thisNode}.
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Dark_Matter_Halo_Scales_Reset(thisNode)

    ! Check if virial temperature is already computed. Compute and store if not.
    if (.not.virialTemperatureComputed) then
       virialTemperatureComputed=.true.
       virialTemperatureStored=0.5d0*atomicMassUnit*meanAtomicMassPrimordial*((kilo&
            &*Dark_Matter_Halo_Virial_Velocity(thisNode))**2)/boltzmannsConstant
    end if

    ! Return the stored temperature.
    Dark_Matter_Halo_Virial_Temperature=virialTemperatureStored
    return
  end function Dark_Matter_Halo_Virial_Temperature

  double precision function Dark_Matter_Halo_Virial_Radius(thisNode)
    !% Returns the virial radius scale for {\tt thisNode}.
    use Numerical_Constants_Math
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Dark_Matter_Halo_Scales_Reset(thisNode)

    ! Check if virial radius is already computed. Compute and store if not.
    if (.not.virialRadiusComputed) then
       virialRadiusComputed=.true.
       virialRadiusStored=(3.0d0*Tree_Node_Mass(thisNode)/4.0d0/Pi/Dark_Matter_Halo_Mean_Density(thisNode))**(1.0d0/3.0d0)
    end if

    ! Return the stored value.
    Dark_Matter_Halo_Virial_Radius=virialRadiusStored
    return
  end function Dark_Matter_Halo_Virial_Radius

  double precision function Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)
    !% Returns the growth rate of the virial radius scale for {\tt thisNode}.
    use Numerical_Constants_Math
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode

    Dark_Matter_Halo_Virial_Radius_Growth_Rate=(1.0d0/3.0d0)*Dark_Matter_Halo_Virial_Radius(thisNode)&
         &*(Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)-Dark_Matter_Halo_Mean_Density_Growth_Rate(thisNode)&
         &/Dark_Matter_Halo_Mean_Density(thisNode))
    return
  end function Dark_Matter_Halo_Virial_Radius_Growth_Rate
  
  double precision function Dark_Matter_Halo_Mean_Density(thisNode)
    !% Returns the mean density for {\tt thisNode}.
    use Cosmological_Parameters
    use Cosmology_Functions
    use Virial_Density_Contrast
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: time
    double precision, save                   :: timePrevious=-1.0d0,densityPrevious
    !$omp threadprivate(timePrevious,densityPrevious)

    ! Get the time at which this halo was last an isolated halo.
    time=Tree_Node_Time_Last_Isolated(thisNode)
    ! If time is not the same as the one previously used then compute its mean density based on mean cosmological density and
    ! overdensity of a collapsing halo, and store it.
    if (time /= timePrevious) then
       timePrevious=time
       densityPrevious=Halo_Virial_Density_Contrast(time)*Omega_0()*Critical_Density()/Expansion_Factor(time)**3
    end if
    ! Return the stored value.
    Dark_Matter_Halo_Mean_Density=densityPrevious
    return
  end function Dark_Matter_Halo_Mean_Density

  double precision function Dark_Matter_Halo_Mean_Density_Growth_Rate(thisNode)
    !% Returns the growth rate of the mean density for {\tt thisNode}.
    use Cosmological_Parameters
    use Cosmology_Functions
    use Virial_Density_Contrast
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: time,aExpansion
    double precision, save                   :: timePrevious=-1.0d0,densityGrowthRatePrevious
    !$omp threadprivate(timePrevious,densityGrowthRatePrevious)

    if (thisNode%isSatellite()) then
       ! Satellite halo is not growing, return zero rate.
       Dark_Matter_Halo_Mean_Density_Growth_Rate=0.0d0
    else
       ! Get the time at which this halo was last an isolated halo.
       time=Tree_Node_Time_Last_Isolated(thisNode)
       ! Check if the time is different from that one previously used.
       if (time /= timePrevious) then
          ! It is not, so recompute the density growth rate.
          timePrevious=time
          ! Get the expansion factor at this time.
          aExpansion=Expansion_Factor(time)
          ! Compute growth rate of its mean density based on mean cosmological density and overdensity of a collapsing halo.
          Dark_Matter_Halo_Mean_Density_Growth_Rate=Dark_Matter_Halo_Mean_Density(thisNode)&
               &*(Halo_Virial_Density_Contrast_Rate_of_Change(time)/Halo_Virial_Density_Contrast(time)-3.0d0&
               &*Expansion_Rate(aExpansion))
       end if
       ! Return the stored value.
       Dark_Matter_Halo_Mean_Density_Growth_Rate=densityGrowthRatePrevious
    end if
    return
  end function Dark_Matter_Halo_Mean_Density_Growth_Rate

end module Dark_Matter_Halo_Scales
