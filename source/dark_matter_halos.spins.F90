!% Contains a module which implements calculations of dark matter halo angular momentum.

module Dark_Matter_Halo_Spins
  !% Implements calculations of dark matter halo angular momentum.
  private
  public :: Dark_Matter_Halo_Angular_Momentum, Dark_Matter_Halo_Angular_Momentum_Growth_Rate

contains

  double precision function Dark_Matter_Halo_Angular_Momentum(thisNode)
    !% Returns the total anuglar momentum of {\tt thisNode} based on its mass, energy and spin parameter.
    use Tree_Nodes
    use Tree_Node_Methods
    use Numerical_Constants_Physical
    use Dark_Matter_Profiles
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    Dark_Matter_Halo_Angular_Momentum=Tree_Node_Spin(thisNode)*gravitationalConstantGalacticus*Tree_Node_Mass(thisNode)**2.5d0 &
         &/dsqrt(dabs(Dark_Matter_Profile_Energy(thisNode)))
    return
  end function Dark_Matter_Halo_Angular_Momentum

  double precision function Dark_Matter_Halo_Angular_Momentum_Growth_Rate(thisNode)
    !% Returns the rate of change of the total anuglar momentum of {\tt thisNode} based on its mass, energy and spin parameter.
    use Tree_Nodes
    use Tree_Node_Methods
    use Numerical_Constants_Physical
    use Dark_Matter_Profiles
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    Dark_Matter_Halo_Angular_Momentum_Growth_Rate=Dark_Matter_Halo_Angular_Momentum(thisNode)&
         &*(Tree_Node_Spin_Growth_Rate(thisNode)/Tree_Node_Spin(thisNode)+2.5d0*Tree_Node_Mass_Accretion_Rate(thisNode)&
         &/Tree_Node_Mass(thisNode)-0.5d0*Dark_Matter_Profile_Energy_Growth_Rate(thisNode)/Dark_Matter_Profile_Energy(thisNode))

    return
  end function Dark_Matter_Halo_Angular_Momentum_Growth_Rate

end module Dark_Matter_Halo_Spins
