!% Contains a module which implements calculations of various scales for dark matter halos.

module Dark_Matter_Halo_Scales
  !% Implements calculations of various scales for dark matter halos.
  use Tree_Nodes
  use Tree_Node_Methods
  private
  public :: Dark_Matter_Halo_Dynamical_Timescale, Dark_Matter_Halo_Virial_Velocity, Dark_Matter_Halo_Virial_Velocity_Growth_Rate,&
       & Dark_Matter_Halo_Virial_Radius, Dark_Matter_Halo_Virial_Radius_Growth_Rate, Dark_Matter_Halo_Mean_Density,&
       & Dark_Matter_Halo_Mean_Density_Growth_Rate, Dark_Matter_Halo_Virial_Temperature

contains

  double precision function Dark_Matter_Halo_Dynamical_Timescale(thisNode)
    !% Returns the dynamical timescale for {\tt thisNode}.
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Dark_Matter_Halo_Dynamical_Timescale=Dark_Matter_Halo_Virial_Radius(thisNode)*(megaParsec/kilo/gigaYear)&
         &/Dark_Matter_Halo_Virial_Velocity(thisNode)
    return
  end function Dark_Matter_Halo_Dynamical_Timescale

  double precision function Dark_Matter_Halo_Virial_Velocity(thisNode)
    !% Returns the virial velocity scale for {\tt thisNode}.
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode

    Dark_Matter_Halo_Virial_Velocity=dsqrt(gravitationalConstantGalacticus*Tree_Node_Mass(thisNode)&
         &/Dark_Matter_Halo_Virial_Radius(thisNode))
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

    Dark_Matter_Halo_Virial_Temperature=0.5d0*atomicMassUnit*meanAtomicMassPrimordial*((kilo&
         &*Dark_Matter_Halo_Virial_Velocity(thisNode))**2)/boltzmannsConstant
    return
  end function Dark_Matter_Halo_Virial_Temperature

  double precision function Dark_Matter_Halo_Virial_Radius(thisNode)
    !% Returns the virial radius scale for {\tt thisNode}.
    use Numerical_Constants_Math
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode

    Dark_Matter_Halo_Virial_Radius=(3.0d0*Tree_Node_Mass(thisNode)/4.0d0/Pi/Dark_Matter_Halo_Mean_Density(thisNode))**(1.0d0&
         &/3.0d0)
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
