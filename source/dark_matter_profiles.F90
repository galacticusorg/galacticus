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






!% Contains a module which implements calculations related to the dark matter halo density profile.

module Dark_Matter_Profiles
  !% Implements calculations related to the dark matter halo density profile.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Dark_Matter_Profile_Rotation_Normalization, Dark_Matter_Profile_Energy, Dark_Matter_Profile_Energy_Growth_Rate,&
       & Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum,Dark_Matter_Profile_Circular_Velocity&
       &,Dark_Matter_Profile_Potential,Dark_Matter_Profile_Enclosed_Mass


  ! Flag to indicate if this module has been initialized.  
  logical              :: darkMatterProfileInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: darkMatterProfileMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Dark_Matter_Profile_Template), pointer :: Dark_Matter_Profile_Rotation_Normalization_Get => null()
  procedure(Dark_Matter_Profile_Template), pointer :: Dark_Matter_Profile_Energy_Get => null()
  procedure(Dark_Matter_Profile_Template), pointer :: Dark_Matter_Profile_Energy_Growth_Rate_Get => null()
  abstract interface
     double precision function Dark_Matter_Profile_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Dark_Matter_Profile_Template
  end interface
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get =>&
       & null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Circular_Velocity_Get => null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Potential_Get         => null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Enclosed_Mass_Get     => null()
  abstract interface
     double precision function Dark_Matter_Profile_Parameter_Template(thisNode,inputParameter)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: inputParameter
     end function Dark_Matter_Profile_Parameter_Template
  end interface

contains

  subroutine Dark_Matter_Profile_Initialize
    !% Initialize the dark matter profile module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="darkMatterProfileMethod" type="moduleUse">
    include 'dark_matter_profiles.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Dark_Matter_Profile_Initialization) 
    ! Initialize if necessary.
    if (.not.darkMatterProfileInitialized) then
       ! Get the halo spin distribution method parameter.
       !@ <inputParameter>
       !@   <name>darkMatterProfileMethod</name>
       !@   <defaultValue>isothermal</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of dark matter halo density profiles.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterProfileMethod',darkMatterProfileMethod,defaultValue='isothermal')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="darkMatterProfileMethod" type="code" action="subroutine">
       !#  <subroutineArgs>darkMatterProfileMethod,Dark_Matter_Profile_Energy_Get,Dark_Matter_Profile_Energy_Growth_Rate_Get,Dark_Matter_Profile_Rotation_Normalization_Get,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get,Dark_Matter_Profile_Circular_Velocity_Get,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get</subroutineArgs>
       include 'dark_matter_profiles.inc'
       !# </include>
       if (.not.(     associated(Dark_Matter_Profile_Energy_Get                               )   &
            &    .and.associated(Dark_Matter_Profile_Energy_Growth_Rate_Get                   )   & 
            &    .and.associated(Dark_Matter_Profile_Rotation_Normalization_Get               )   &
            &    .and.associated(Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get)   &
            &    .and.associated(Dark_Matter_Profile_Circular_Velocity_Get                    )   &
            &    .and.associated(Dark_Matter_Profile_Potential_Get                            )   &
            &    .and.associated(Dark_Matter_Profile_Enclosed_Mass_Get                        ))) &
            & call Galacticus_Error_Report('Dark_Matter_Profile','method ' //char(darkMatterProfileMethod)//' is unrecognized')
       darkMatterProfileInitialized=.true.
    end if
    !$omp end critical(Dark_Matter_Profile_Initialization) 

    return
  end subroutine Dark_Matter_Profile_Initialize

  double precision function Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum(thisNode,specificAngularMomentum)
    !% Returns the radius (in Mpc) in the dark matter profile of {\tt thisNode} at which the specific angular momentum of a
    !% circular orbit equals {\tt specificAngularMomentum} (specified in units of km s$^{-1}$ Mpc.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: specificAngularMomentum

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the energy using the selected method.
    Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum&
         &=Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get(thisNode,specificAngularMomentum)

    return
  end function Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum

  double precision function Dark_Matter_Profile_Enclosed_Mass(thisNode,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: radius

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the enclosed mass using the selected method.
    Dark_Matter_Profile_Enclosed_Mass=Dark_Matter_Profile_Enclosed_Mass_Get(thisNode,radius)

    return
  end function Dark_Matter_Profile_Enclosed_Mass

  double precision function Dark_Matter_Profile_Circular_Velocity(thisNode,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: radius

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the circular velocity using the selected method.
    Dark_Matter_Profile_Circular_Velocity=Dark_Matter_Profile_Circular_Velocity_Get(thisNode,radius)

    return
  end function Dark_Matter_Profile_Circular_Velocity

  double precision function Dark_Matter_Profile_Potential(thisNode,radius)
    !% Returns the gravitational potential (in (km/s)$^2$) in the dark matter profile of {\tt thisNode} at the given {\tt radius}
    !% (given in units of Mpc).
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: radius

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the circular velocity using the selected method.
    Dark_Matter_Profile_Potential=Dark_Matter_Profile_Potential_Get(thisNode,radius)

    return
  end function Dark_Matter_Profile_Potential

  double precision function Dark_Matter_Profile_Rotation_Normalization(thisNode)
    !% Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in
    !% radius) for {\tt thisNode}. Specifically, the normalization, $A$, returned is such that $V_{\rm rot} = A J/M$.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the energy using the selected method.
    Dark_Matter_Profile_Rotation_Normalization=Dark_Matter_Profile_Rotation_Normalization_Get(thisNode)

    return
  end function Dark_Matter_Profile_Rotation_Normalization

  double precision function Dark_Matter_Profile_Energy(thisNode)
    !% Returns the total energy of {\tt thisNode} in units of $M_\odot$ km$^2$ s$^{-1}$.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the energy using the selected method.
    Dark_Matter_Profile_Energy=Dark_Matter_Profile_Energy_Get(thisNode)

    return
  end function Dark_Matter_Profile_Energy

  double precision function Dark_Matter_Profile_Energy_Growth_Rate(thisNode)
    !% Returns the rate of chance of the total energy of {\tt thisNode} in units of $M_\odot$ km$^2$ s$^{-1}$ Gyr$^{-1}$.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the energy using the selected method.
    Dark_Matter_Profile_Energy_Growth_Rate=Dark_Matter_Profile_Energy_Growth_Rate_Get(thisNode)

    return
  end function Dark_Matter_Profile_Energy_Growth_Rate

end module Dark_Matter_Profiles
