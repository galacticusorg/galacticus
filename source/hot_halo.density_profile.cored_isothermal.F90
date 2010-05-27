!% Contains a module which implements a cored isothermal profile for hot gas halos.

module Hot_Halo_Density_Profile_Cored_Isothermal
  !% Implements a cored isothermal profile for hot gas halos.
  private
  public :: Hot_Halo_Density_Cored_Isothermal

  ! Parameter which specifies the ratio of the cored isothermal profile core radius to the virial radius.
  double precision :: isothermalCoreRadiusOverVirialRadius

  ! Pre-computed factor that appears in the density normalization.
  double precision :: denistyNormalizationFactor

contains

  !# <hotHaloDensityMethod>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal</unitName>
  !# </hotHaloDensityMethod>
  subroutine Hot_Halo_Density_Cored_Isothermal(hotHaloDensityMethod,Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get)
    !% Initialize the cored isothermal hot halo density profile module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: hotHaloDensityMethod
    procedure(),          pointer, intent(inout) :: Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get
    
    if (hotHaloDensityMethod == 'cored isothermal') then
       Hot_Halo_Density_Get => Hot_Halo_Density_Cored_Isothermal_Get
       Hot_Halo_Density_Log_Slope_Get => Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get
       !@ <inputParameter>
       !@   <name>isothermalCoreRadiusOverVirialRadius</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The core radius in the ``cored isothermal'' hot halo density profile in units of the virial radius.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('isothermalCoreRadiusOverVirialRadius',isothermalCoreRadiusOverVirialRadius,defaultValue=0.1d0)
       denistyNormalizationFactor=1.0d0/(1.0d0/isothermalCoreRadiusOverVirialRadius-datan(1.0d0&
            &/isothermalCoreRadiusOverVirialRadius))
    end if
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal

  double precision function Hot_Halo_Density_Cored_Isothermal_Get(thisNode,radius)
    !% Compute the density at radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Tree_Node_Methods
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: hotGasMass,virialRadius,coreRadius,densityNormalization

    hotGasMass=Tree_Node_Hot_Halo_Mass(thisNode)
    virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
    coreRadius=isothermalCoreRadiusOverVirialRadius*virialRadius
    densityNormalization=denistyNormalizationFactor*hotGasMass/4.0d0/Pi/(coreRadius**3)
    Hot_Halo_Density_Cored_Isothermal_Get=densityNormalization/(1.0d0+(radius/coreRadius)**2)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Get
  
  double precision function Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get(thisNode,radius)
    !% Compute the density at radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Tree_Node_Methods
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: virialRadius,coreRadius,radiusInCoreUnitsSquared

    virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
    coreRadius=isothermalCoreRadiusOverVirialRadius*virialRadius
    radiusInCoreUnitsSquared=(radius/coreRadius)**2
    Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get=-2.0d0*radiusInCoreUnitsSquared/(1.0d0+radiusInCoreUnitsSquared)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get
  
end module Hot_Halo_Density_Profile_Cored_Isothermal
