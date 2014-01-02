!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of various scales for dark matter halos.

module Dark_Matter_Halo_Scales
  !% Implements calculations of various scales for dark matter halos.
  use Galacticus_Nodes
  use Kind_Numbers
  use Tables
  implicit none
  private
  public :: Dark_Matter_Halo_Dynamical_Timescale, Dark_Matter_Halo_Virial_Velocity,&
       & Dark_Matter_Halo_Virial_Velocity_Growth_Rate, Dark_Matter_Halo_Virial_Radius,&
       & Dark_Matter_Halo_Virial_Radius_Growth_Rate, Dark_Matter_Halo_Mean_Density,&
       & Dark_Matter_Halo_Mean_Density_Growth_Rate,&
       & Dark_Matter_Halo_Virial_Temperature, Dark_Matter_Halo_Scales_Reset, &
       & Dark_Matter_Halo_Scales_State_Store,Dark_Matter_Halo_Scales_State_Retrieve

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8          )            :: lastUniqueID                   =-1
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not halo scales have already been computed for this node.
  logical                                               :: dynamicalTimescaleComputed     =.false., virialRadiusComputed  =.false., &
       &                                                   virialTemperatureComputed      =.false., virialVelocityComputed=.false.
  !$omp threadprivate(virialRadiusComputed,virialTemperatureComputed,virialVelocityComputed,dynamicalTimescaleComputed)
  ! Stored values of halo scales.
  double precision                                      :: dynamicalTimescaleStored               , virialRadiusStored            , &
       &                                                   virialTemperatureStored                , virialVelocityStored
  !$omp threadprivate(virialRadiusStored,virialTemperatureStored,virialVelocityStored,dynamicalTimescaleStored)
  ! Table for fast lookup of the mean density of halos.
  double precision                                      :: meanDensityTimeMaximum         =-1.0d0 , meanDensityTimeMinimum=-1.0d0
  integer                                   , parameter :: meanDensityTablePointsPerDecade=100
  type            (table1DLogarithmicLinear)            :: meanDensityTable
  logical                                               :: resetMeanDensityTable
  !$omp threadprivate(meanDensityTable,meanDensityTimeMinimum,meanDensityTimeMaximum,resetMeanDensityTable)
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
    type (treeNode          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic)               , pointer :: thisBasicComponent

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Dark_Matter_Halo_Scales_Reset(thisNode)

    ! Check if virial velocity is already computed. Compute and store if not.
    if (.not.virialVelocityComputed) then
       ! Get the basic component.
       thisBasicComponent => thisNode%basic()
       ! Compute the virial velocity.
       virialVelocityStored=sqrt(gravitationalConstantGalacticus*thisBasicComponent%mass() &
            &/Dark_Matter_Halo_Virial_Radius(thisNode))
       ! Record that virial velocity has now been computed.
       virialVelocityComputed=.true.
    end if

    ! Return the stored virial velocity.
    Dark_Matter_Halo_Virial_Velocity=virialVelocityStored
    return
  end function Dark_Matter_Halo_Virial_Velocity

  double precision function Dark_Matter_Halo_Virial_Velocity_Growth_Rate(thisNode)
    !% Returns the growth rate of the virial velocity scale for {\tt thisNode}.
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic)               , pointer :: thisBasicComponent

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Halo_Virial_Velocity_Growth_Rate=0.5d0*Dark_Matter_Halo_Virial_Velocity(thisNode) &
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()-Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode) &
         &/Dark_Matter_Halo_Virial_Radius(thisNode))
    return
  end function Dark_Matter_Halo_Virial_Velocity_Growth_Rate

  double precision function Dark_Matter_Halo_Virial_Temperature(thisNode)
    !% Returns the virial temperature (in Kelvin) for {\tt thisNode}.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

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
    type (treeNode          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic)               , pointer :: thisBasicComponent

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Dark_Matter_Halo_Scales_Reset(thisNode)

    ! Check if virial radius is already computed. Compute and store if not.
    if (.not.virialRadiusComputed) then
       ! Get the basic component.
       thisBasicComponent => thisNode%basic()
       ! Compute the virial radius.
       virialRadiusStored=(3.0d0*thisBasicComponent%mass()/4.0d0/Pi/Dark_Matter_Halo_Mean_Density(thisNode))**(1.0d0/3.0d0)
       ! Record that the virial radius has been computed.
       virialRadiusComputed=.true.
    end if

    ! Return the stored value.
    Dark_Matter_Halo_Virial_Radius=virialRadiusStored
    return
  end function Dark_Matter_Halo_Virial_Radius

  double precision function Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)
    !% Returns the growth rate of the virial radius scale for {\tt thisNode}.
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic)               , pointer :: thisBasicComponent

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Halo_Virial_Radius_Growth_Rate=(1.0d0/3.0d0)*Dark_Matter_Halo_Virial_Radius(thisNode)&
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()-Dark_Matter_Halo_Mean_Density_Growth_Rate(thisNode)&
         &/Dark_Matter_Halo_Mean_Density(thisNode))
    return
  end function Dark_Matter_Halo_Virial_Radius_Growth_Rate

  double precision function Dark_Matter_Halo_Mean_Density(thisNode)
    !% Returns the mean density for {\tt thisNode}.
    use Cosmology_Parameters
    use Cosmology_Functions
    use Virial_Density_Contrast
    implicit none
    type            (treeNode                ), intent(inout), pointer :: thisNode
    class           (nodeComponentBasic      )               , pointer :: thisBasicComponent
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    class           (cosmologyFunctionsClass )               , pointer :: cosmologyFunctionsDefault
    integer                                                            :: i                        , meanDensityTablePoints
    double precision                                                   :: time

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()
    ! Get the time at which this halo was last an isolated halo.
    time=thisBasicComponent%timeLastIsolated()
    if (time <= 0.0d0) time=thisBasicComponent%time()
    ! Retabulate the mean density vs. time if necessary.
    if (resetMeanDensityTable .or. time < meanDensityTimeMinimum .or. time > meanDensityTimeMaximum) then
       resetMeanDensityTable=.false.
       if (meanDensityTimeMinimum <= 0.0d0) then
          meanDensityTimeMinimum=                           time/10.0d0
          meanDensityTimeMaximum=                           time* 2.0d0
       else
          meanDensityTimeMinimum=min(meanDensityTimeMinimum,time/10.0d0)
          meanDensityTimeMaximum=max(meanDensityTimeMaximum,time* 2.0d0)
       end if
       meanDensityTablePoints=int(log10(meanDensityTimeMaximum/meanDensityTimeMinimum)*dble(meanDensityTablePointsPerDecade))+1
       call meanDensityTable%destroy()
       call meanDensityTable%create(meanDensityTimeMinimum,meanDensityTimeMaximum,meanDensityTablePoints)
       ! Get the default cosmology.
       thisCosmologyParameters => cosmologyParameters()
       ! Get the default cosmology functions object.
       cosmologyFunctionsDefault => cosmologyFunctions()
       do i=1,meanDensityTablePoints
          call meanDensityTable%populate(Halo_Virial_Density_Contrast(meanDensityTable%x(i))*thisCosmologyParameters%OmegaMatter()*thisCosmologyParameters%densityCritical()/cosmologyFunctionsDefault%expansionFactor(meanDensityTable%x(i))**3,i)
       end do
    end if
    ! Return the stored value.
    Dark_Matter_Halo_Mean_Density=meanDensityTable%interpolate(time)
    return
  end function Dark_Matter_Halo_Mean_Density

  double precision function Dark_Matter_Halo_Mean_Density_Growth_Rate(thisNode)
    !% Returns the growth rate of the mean density for {\tt thisNode}.
    use Cosmology_Functions
    use Virial_Density_Contrast
    implicit none
    type            (treeNode               ), intent(inout), pointer :: thisNode
    class           (nodeComponentBasic     )               , pointer :: thisBasicComponent
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctionsDefault
    double precision                                                  :: aExpansion               , time
    double precision                         , save                   :: densityGrowthRatePrevious, timePrevious=-1.0d0
    !$omp threadprivate(timePrevious,densityGrowthRatePrevious)
    if (thisNode%isSatellite()) then
       ! Satellite halo is not growing, return zero rate.
       Dark_Matter_Halo_Mean_Density_Growth_Rate=0.0d0
    else
       ! Get the basic component.
       thisBasicComponent => thisNode%basic()
       ! Get the time at which this halo was last an isolated halo.
       time=thisBasicComponent%timeLastIsolated()
       ! Check if the time is different from that one previously used.
       if (time /= timePrevious) then
          ! It is not, so recompute the density growth rate.
          timePrevious=time
          ! Get the default cosmology functions object.
          cosmologyFunctionsDefault => cosmologyFunctions()
          ! Get the expansion factor at this time.
          aExpansion=cosmologyFunctionsDefault%expansionFactor(time)
          ! Compute growth rate of its mean density based on mean cosmological density and overdensity of a collapsing halo.
          densityGrowthRatePrevious=Dark_Matter_Halo_Mean_Density(thisNode)&
               &*(Halo_Virial_Density_Contrast_Rate_of_Change(time)/Halo_Virial_Density_Contrast(time)-3.0d0&
               &*cosmologyFunctionsDefault%expansionRate(aExpansion))
       end if
       ! Return the stored value.
       Dark_Matter_Halo_Mean_Density_Growth_Rate=densityGrowthRatePrevious
    end if
    return
  end function Dark_Matter_Halo_Mean_Density_Growth_Rate

  !# <galacticusStateStoreTask>
  !#  <unitName>Dark_Matter_Halo_Scales_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Dark_Matter_Halo_Scales_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    write (stateFile) meanDensityTimeMinimum,meanDensityTimeMaximum
    return
  end subroutine Dark_Matter_Halo_Scales_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Dark_Matter_Halo_Scales_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Dark_Matter_Halo_Scales_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    read (stateFile) meanDensityTimeMinimum,meanDensityTimeMaximum
    ! Ensure that interpolation objects will get reset.
    resetMeanDensityTable=.true.
    return
  end subroutine Dark_Matter_Halo_Scales_State_Retrieve

end module Dark_Matter_Halo_Scales
