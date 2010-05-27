!% Contains a module which implements calculations of baryonic accretion onto halos using a simple truncation to mimic
!% reionization.

module Accretion_Halos_Simple
  !% Implements calculations of baryonic accretion onto halos using a simple truncation to mimic reionization.
  private
  public :: Accretion_Halos_Simple_Initialize

  ! Parameters controlling when accretion is suppressed.
  double precision :: reionizationSuppressionRedshift,reionizationSuppressionTime,reionizationSuppressionVelocity

contains

  !# <accretionHalosMethod>
  !#  <unitName>Accretion_Halos_Simple_Initialize</unitName>
  !# </accretionHalosMethod>
  subroutine Accretion_Halos_Simple_Initialize(accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get &
       &,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Cosmology_Functions
    implicit none
    type(varying_string),          intent(in)    :: accretionHalosMethod
    procedure(),          pointer, intent(inout) :: Halo_Baryonic_Accretion_Rate_Get,Halo_Baryonic_Accreted_Mass_Get&
         &,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get
    
    if (accretionHalosMethod == 'simple') then
       Halo_Baryonic_Accretion_Rate_Get        => Halo_Baryonic_Accretion_Rate_Simple_Get
       Halo_Baryonic_Accreted_Mass_Get         => Halo_Baryonic_Accreted_Mass_Simple_Get
       Halo_Baryonic_Failed_Accretion_Rate_Get => Halo_Baryonic_Failed_Accretion_Rate_Simple_Get
       Halo_Baryonic_Failed_Accreted_Mass_Get  => Halo_Baryonic_Failed_Accreted_Mass_Simple_Get
       !@ <inputParameter>
       !@   <name>reionizationSuppressionRedshift</name>
       !@   <defaultValue>8.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The redshift below which baryonic accretion is suppressed.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("reionizationSuppressionRedshift",reionizationSuppressionRedshift,defaultValue= 8.0d0)
       reionizationSuppressionTime=Cosmology_Age(Expansion_Factor_from_Redshift(reionizationSuppressionRedshift))
       !@ <inputParameter>
       !@   <name>reionizationSuppressionVelocity</name>
       !@   <defaultValue>30.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The velocity scale below which baryonic accretion is suppressed.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("reionizationSuppressionVelocity",reionizationSuppressionVelocity,defaultValue=30.0d0)
    end if
    return
  end subroutine Accretion_Halos_Simple_Initialize

  double precision function Halo_Baryonic_Accretion_Rate_Simple_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: growthRate,unaccretedMass

    if (thisNode%isSatellite()) then
       Halo_Baryonic_Accretion_Rate_Simple_Get=0.0d0
    else
       if (Tree_Node_Time(thisNode) > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
            & reionizationSuppressionVelocity) then
          Halo_Baryonic_Accretion_Rate_Simple_Get=0.0d0
       else
          Halo_Baryonic_Accretion_Rate_Simple_Get=(Omega_b()/Omega_0())*Tree_Node_Mass_Accretion_Rate(thisNode)
          unaccretedMass=Tree_Node_Hot_Halo_Unaccreted_Mass(thisNode)
          if (unaccretedMass > 0.0d0) then
             growthRate=Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)
             Halo_Baryonic_Accretion_Rate_Simple_Get=Halo_Baryonic_Accretion_Rate_Simple_Get+unaccretedMass*growthRate
          end if
       end if
    end if
    return
  end function Halo_Baryonic_Accretion_Rate_Simple_Get

  double precision function Halo_Baryonic_Accreted_Mass_Simple_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (Tree_Node_Time(thisNode) > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
         & reionizationSuppressionVelocity) then
       Halo_Baryonic_Accreted_Mass_Simple_Get=0.0d0
    else
       Halo_Baryonic_Accreted_Mass_Simple_Get=(Omega_b()/Omega_0())*Tree_Node_Mass(thisNode)
    end if
    return
  end function Halo_Baryonic_Accreted_Mass_Simple_Get
  
  double precision function Halo_Baryonic_Failed_Accretion_Rate_Simple_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: growthRate,unaccretedMass

    if (thisNode%isSatellite()) then
       Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=0.0d0
    else
       if (Tree_Node_Time(thisNode) > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
            & reionizationSuppressionVelocity) then
          Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=(Omega_b()/Omega_0())*Tree_Node_Mass_Accretion_Rate(thisNode)
       else
          unaccretedMass=Tree_Node_Hot_Halo_Unaccreted_Mass(thisNode)
          if (unaccretedMass > 0.0d0) then
             growthRate=Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)
             Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=-unaccretedMass*growthRate
          else
             Halo_Baryonic_Failed_Accretion_Rate_Simple_Get=0.0d0
          end if
       end if
    end if
    return
  end function Halo_Baryonic_Failed_Accretion_Rate_Simple_Get

  double precision function Halo_Baryonic_Failed_Accreted_Mass_Simple_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (Tree_Node_Time(thisNode) > reionizationSuppressionTime .and. Dark_Matter_Halo_Virial_Velocity(thisNode) <&
         & reionizationSuppressionVelocity) then
       Halo_Baryonic_Failed_Accreted_Mass_Simple_Get=(Omega_b()/Omega_0())*Tree_Node_Mass(thisNode)
    else
       Halo_Baryonic_Failed_Accreted_Mass_Simple_Get=0.0d0
    end if
    return
  end function Halo_Baryonic_Failed_Accreted_Mass_Simple_Get
  
end module Accretion_Halos_Simple
