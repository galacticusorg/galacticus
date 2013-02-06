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

!% Contains a module which implements a cored isothermal profile for hot gas halos.

module Hot_Halo_Density_Profile_Cored_Isothermal
  !% Implements a cored isothermal profile for hot gas halos.
  implicit none
  private
  public :: Hot_Halo_Density_Cored_Isothermal

  ! Record of whether an active component can supply the hot halo mass property.
  logical :: hotHaloActive

contains

  !# <hotHaloDensityMethod>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal</unitName>
  !# </hotHaloDensityMethod>
  subroutine Hot_Halo_Density_Cored_Isothermal(hotHaloDensityMethod,Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get&
       &,Hot_Halo_Enclosed_Mass_Get,Hot_Halo_Profile_Rotation_Normalization_Get,Hot_Halo_Profile_Radial_Moment_Get)
    !% Initialize the cored isothermal hot halo density profile module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Nodes
    implicit none
    type(varying_string),                 intent(in)    :: hotHaloDensityMethod
    procedure(double precision), pointer, intent(inout) :: Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get&
         &,Hot_Halo_Enclosed_Mass_Get,Hot_Halo_Profile_Rotation_Normalization_Get,Hot_Halo_Profile_Radial_Moment_Get
    
    if (hotHaloDensityMethod == 'coredIsothermal') then
       Hot_Halo_Density_Get                        => Hot_Halo_Density_Cored_Isothermal_Get
       Hot_Halo_Density_Log_Slope_Get              => Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get
       Hot_Halo_Enclosed_Mass_Get                  => Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get
       Hot_Halo_Profile_Rotation_Normalization_Get => Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get
       Hot_Halo_Profile_Radial_Moment_Get          => Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get
       ! Detect whether there is an active component which can provide a hot gas mass.
       hotHaloActive=defaultHotHaloComponent%massIsGettable()
    end if
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal

  double precision function Density_Normalization_Factor(coreRadius,outerRadius)
    !% Computes the density profile normalization factor for a given core radius and outer radius.
    implicit none
    double precision, intent(in) :: coreRadius,outerRadius
    double precision             :: outerRadiusOverCoreRadius    
    double precision, save       :: outerRadiusOverCoreRadiusPrevious=-1.0d0,densityNormalizationPrevious
    !$omp threadprivate(outerRadiusOverCoreRadiusPrevious,densityNormalizationPrevious)
    double precision, parameter  :: outerRadiusOverCoreRadiusSmall=1.0d-6

    outerRadiusOverCoreRadius=outerRadius/coreRadius
    if (outerRadiusOverCoreRadius /= outerRadiusOverCoreRadiusPrevious) then
       outerRadiusOverCoreRadiusPrevious=outerRadiusOverCoreRadius
       if      (outerRadiusOverCoreRadius <= 0.0d0                         ) then
          ! For zero or negative outer radius, return a zero normalization.
          densityNormalizationPrevious=0.0d0
       else if (outerRadiusOverCoreRadius <  outerRadiusOverCoreRadiusSmall) then 
          ! For small outer radii, use a series approximation to the exact solution.
          densityNormalizationPrevious=3.0d0/outerRadiusOverCoreRadius**3+9.0d0/5.0d0/outerRadiusOverCoreRadius
       else
          ! For larger outer radii, use the exact solution.
          densityNormalizationPrevious=1.0d0/(outerRadiusOverCoreRadius-datan(outerRadiusOverCoreRadius))
       end if
    end if
    Density_Normalization_Factor=densityNormalizationPrevious
    return
  end function Density_Normalization_Factor

  double precision function Hot_Halo_Density_Cored_Isothermal_Get(thisNode,radius)
    !% Compute the density at radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Galacticus_Nodes
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    use Numerical_Constants_Math
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    double precision           , intent(in   )          :: radius
    class(nodeComponentHotHalo),                pointer :: thisHotHaloComponent
    double precision                                    :: hotGasMass,outerRadius,coreRadius,densityNormalization

    thisHotHaloComponent => thisNode%hotHalo()
    hotGasMass           =  thisHotHaloComponent%mass       ()
    outerRadius          =  thisHotHaloComponent%outerRadius()
    coreRadius           =  Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)
    densityNormalization =  Density_Normalization_Factor(coreRadius,outerRadius)*hotGasMass/4.0d0/Pi/(coreRadius**3)
    Hot_Halo_Density_Cored_Isothermal_Get=densityNormalization/(1.0d0+(radius/coreRadius)**2)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Get
  
  double precision function Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get(thisNode,radius)
    !% Compute the density at radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Galacticus_Nodes
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: coreRadius,radiusInCoreUnitsSquared

    coreRadius              =Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)
    radiusInCoreUnitsSquared=(radius/coreRadius)**2
    Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get=-2.0d0*radiusInCoreUnitsSquared/(1.0d0+radiusInCoreUnitsSquared)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get
  
  double precision function Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get(thisNode,radius)
    !% Compute the mass enclosed within radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Galacticus_Nodes
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    use Galactic_Structure_Options
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    double precision,            intent(in   )          :: radius
    class(nodeComponentHotHalo),                pointer :: thisHotHaloComponent
    double precision                                    :: hotGasMass,outerRadius,coreRadius

    ! Return immediately with zero mass if no active component can supply a hot halo mass.
    if (.not.hotHaloActive) then
       Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get=0.0d0
       return
    end if

    thisHotHaloComponent => thisNode%hotHalo()
    hotGasMass           =  thisHotHaloComponent%mass       ()

    if (hotGasMass <= 0.0d0) then
       Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get=0.0d0
       return
    end if
    outerRadius          =  thisHotHaloComponent%outerRadius()
    ! Truncate the profile at the virial radius.
    if (radius >= outerRadius) then
       Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get=hotGasMass
       return
    end if
    coreRadius =Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)
    Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get=hotGasMass*Density_Normalization_Factor(coreRadius,outerRadius)*(radius&
         &/coreRadius-datan(radius/coreRadius))
    return
  end function Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get

  double precision function Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get(thisNode)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentHotHalo),                pointer :: thisHotHaloComponent
    double precision                                    :: radiusOuter,radiusCoreOverRadiusOuter

    ! Get outer radius and ratio of core radius to virial radius.
    thisHotHaloComponent => thisNode            %hotHalo    ()
    radiusOuter          =  thisHotHaloComponent%outerRadius()
    if (radiusOuter <= 0.0d0) then
       Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get=0.0d0
       return
    end if
    radiusCoreOverRadiusOuter=Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)/radiusOuter

    ! Compute the normalization.
    Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get=                                                  &
         &                                                       (                                                 &
         &                                                        1.0d0                                            &
         &                                                       - radiusCoreOverRadiusOuter                       &
         &                                                        *datan(1.0d0/radiusCoreOverRadiusOuter)          &
         &                                                      )                                                  &
         &                                                     /(                                                  &
         &                                                        0.5d0                                            &
         &                                                       + radiusCoreOverRadiusOuter**2                    &
         &                                                        *dlog(                                           &
         &                                                               radiusCoreOverRadiusOuter                 &
         &                                                              /dsqrt(1.0d0+radiusCoreOverRadiusOuter**2) &
         &                                                             )                                           &
         &                                                      )                                                  &
         &                                                     /radiusOuter
    return
  end function Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get
  
  double precision function Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get(thisNode,moment,radius)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    use Galacticus_Error
    use Numerical_Comparison
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo),                pointer :: thisHotHaloComponent
    double precision                      , intent(in   )          :: moment,radius
    double precision                                               :: radiusOuter,radiusCore,radiusOuterOverRadiusCore,hotGasMass&
         & ,densityNormalization

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()

    ! Get outer radius and ratio of core radius to outer radius.
    radiusOuter=max(radius,thisHotHaloComponent%outerRadius())
    if (radiusOuter <= 0.0d0) then
       Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get=0.0d0
       return
    end if
    radiusCore               =Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)
    radiusOuterOverRadiusCore=radiusOuter/radiusCore

    ! Get the density normalization
    hotGasMass=thisHotHaloComponent%mass()
    if (hotGasMass <= 0.0d0) then
       Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get=0.0d0
       return
    end if
    densityNormalization=Density_Normalization_Factor(radiusCore,radiusOuter)*hotGasMass/4.0d0/Pi/(radiusCore**3)

    ! Compute the moment.
    if      (Values_Agree(moment,2.0d0,absTol=1.0d-3)) then
       Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get=(radiusOuterOverRadiusCore-atan(radiusOuterOverRadiusCore))
    else if (Values_Agree(moment,3.0d0,absTol=1.0d-3)) then
       Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get=0.5d0*(radiusOuterOverRadiusCore**2-log(1.0d0+radiusOuterOverRadiusCore**2))
    else
       ! Abort for unsupported moments.
       call Galacticus_Error_Report('Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get','only 2nd and 3rd moments are supported')
    end if
    ! Make the result dimensionful.
    Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get=                                                 &
         &   Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get*densityNormalization*radiusCore**moment
    return
  end function Hot_Halo_Profile_Radial_Moment_Cored_Isothermal_Get
  
end module Hot_Halo_Density_Profile_Cored_Isothermal
