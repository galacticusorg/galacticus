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

!% Contains a module which implements calculations related to the dark matter halo density profile.

module Dark_Matter_Profiles
  !% Implements calculations related to the dark matter halo density profile.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Dark_Matter_Profile_Rotation_Normalization, Dark_Matter_Profile_Energy, Dark_Matter_Profile_Energy_Growth_Rate,&
       & Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum,Dark_Matter_Profile_Circular_Velocity &
       &,Dark_Matter_Profile_Potential,Dark_Matter_Profile_Enclosed_Mass,Dark_Matter_Profile_kSpace,&
       & Dark_Matter_Profile_Density_Task, Dark_Matter_Profile_Density,Dark_Matter_Profile_Rotation_Curve_Task&
       &,Dark_Matter_Profile_Potential_Task, Dark_Matter_Profile_Enclosed_Mass_Task,Dark_Matter_Profile_Freefall_Radius &
       &,Dark_Matter_Profile_Freefall_Radius_Increase_Rate,Dark_Matter_Profile_Rotation_Curve_Gradient_Task

  ! Flag to indicate if this module has been initialized.  
  logical              :: darkMatterProfileInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: darkMatterProfileMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Dark_Matter_Profile_Template), pointer :: Dark_Matter_Profile_Rotation_Normalization_Get => null()
  procedure(Dark_Matter_Profile_Template), pointer :: Dark_Matter_Profile_Energy_Get                 => null()
  procedure(Dark_Matter_Profile_Template), pointer :: Dark_Matter_Profile_Energy_Growth_Rate_Get     => null()
  abstract interface
     double precision function Dark_Matter_Profile_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Dark_Matter_Profile_Template
  end interface
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get =>&
       & null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Circular_Velocity_Get             => null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Potential_Get                     => null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Enclosed_Mass_Get                 => null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_kSpace_Get                        => null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Density_Get                       => null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Freefall_Radius_Get               => null()
  procedure(Dark_Matter_Profile_Parameter_Template), pointer :: Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get => null()
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
    
    ! Initialize if necessary.
    if (.not.darkMatterProfileInitialized) then
       !$omp critical(Dark_Matter_Profile_Initialization) 
       if (.not.darkMatterProfileInitialized) then
          ! Get the halo spin distribution method parameter.
          !@ <inputParameter>
          !@   <name>darkMatterProfileMethod</name>
          !@   <defaultValue>NFW</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of dark matter halo density profiles.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('darkMatterProfileMethod',darkMatterProfileMethod,defaultValue='NFW')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="darkMatterProfileMethod" type="code" action="subroutine">
          !#  <subroutineArgs>darkMatterProfileMethod,Dark_Matter_Profile_Density_Get,Dark_Matter_Profile_Energy_Get,Dark_Matter_Profile_Energy_Growth_Rate_Get,Dark_Matter_Profile_Rotation_Normalization_Get,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get,Dark_Matter_Profile_Circular_Velocity_Get,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get,Dark_Matter_Profile_kSpace_Get,Dark_Matter_Profile_Freefall_Radius_Get,Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get</subroutineArgs>
          include 'dark_matter_profiles.inc'
          !# </include>
          if (.not.(     associated(Dark_Matter_Profile_Density_Get                              )   &
               &    .and.associated(Dark_Matter_Profile_Energy_Get                               )   & 
               &    .and.associated(Dark_Matter_Profile_Energy_Growth_Rate_Get                   )   & 
               &    .and.associated(Dark_Matter_Profile_Rotation_Normalization_Get               )   &
               &    .and.associated(Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get)   &
               &    .and.associated(Dark_Matter_Profile_Circular_Velocity_Get                    )   &
               &    .and.associated(Dark_Matter_Profile_Potential_Get                            )   &
               &    .and.associated(Dark_Matter_Profile_Enclosed_Mass_Get                        )   &
               &    .and.associated(Dark_Matter_Profile_kSpace_Get                               )   &
               &    .and.associated(Dark_Matter_Profile_Freefall_Radius_Get                      )   &
               &    .and.associated(Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get        ))) &
               & call Galacticus_Error_Report('Dark_Matter_Profile','method ' //char(darkMatterProfileMethod)//' is unrecognized')
          darkMatterProfileInitialized=.true.
       end if
       !$omp end critical(Dark_Matter_Profile_Initialization) 
    end if
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
         & =Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get(thisNode,specificAngularMomentum)

    return
  end function Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum

  double precision function Dark_Matter_Profile_Density(thisNode,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given
    !% in units of Mpc).
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: radius

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the enclosed mass using the selected method.
    Dark_Matter_Profile_Density=Dark_Matter_Profile_Density_Get(thisNode,radius)

    return
  end function Dark_Matter_Profile_Density

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

  double precision function Dark_Matter_Profile_kSpace(thisNode,waveNumber)
    !% Returns the normalized Fourier space density profile of the dark matter profile of {\tt thisNode} at the given {\tt waveNumber}
    !% (given in units of Mpc$^{-1}$).
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: waveNumber

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the circular velocity using the selected method.
    Dark_Matter_Profile_kSpace=Dark_Matter_Profile_kSpace_Get(thisNode,waveNumber)

    return
  end function Dark_Matter_Profile_kSpace

  double precision function Dark_Matter_Profile_Freefall_Radius(thisNode,time)
    !% Returns the freefall radius (in Mpc) corresponding to the given {\tt time} (in Gyr) in {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: time

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the freefall radius using the selected method.
    Dark_Matter_Profile_Freefall_Radius=Dark_Matter_Profile_Freefall_Radius_Get(thisNode,time)

    return
  end function Dark_Matter_Profile_Freefall_Radius

  double precision function Dark_Matter_Profile_Freefall_Radius_Increase_Rate(thisNode,time)
    !% Returns the rate of increase of the freefall radius (in Mpc/Gyr) corresponding to the given {\tt time} (in Gyr) in {\tt
    !% thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: time

    ! Initialize the module.
    call Dark_Matter_Profile_Initialize

    ! Get the increase rate using the selected method.
    Dark_Matter_Profile_Freefall_Radius_Increase_Rate=Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get(thisNode,time)

    return
  end function Dark_Matter_Profile_Freefall_Radius_Increase_Rate

  !# <enclosedMassTask>
  !#  <unitName>Dark_Matter_Profile_Enclosed_Mass_Task</unitName>
  !# </enclosedMassTask>
  subroutine Dark_Matter_Profile_Enclosed_Mass_Task(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Cosmological_Parameters
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType,weightBy,weightIndex
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass
    
    componentMass=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    if (.not.(weightBy      == weightByMass                                                )) return

    if (radius >= radiusLarge) then
       ! Return the total mass of the halo in this case.
       componentMass=Tree_Node_Mass(thisNode)
    else
       ! Return the mass within the specified radius.
       componentMass=Dark_Matter_Profile_Enclosed_Mass(thisNode,radius)
    end if
    ! Scale to account for just the dark component.
    componentMass=componentMass*(Omega_Matter()-Omega_b())/Omega_Matter()
    return
  end subroutine Dark_Matter_Profile_Enclosed_Mass_Task

  !# <rotationCurveTask>
  !#  <unitName>Dark_Matter_Profile_Rotation_Curve_Task</unitName>
  !# </rotationCurveTask>
  subroutine Dark_Matter_Profile_Rotation_Curve_Task(thisNode,radius,massType,componentType,componentVelocity)
    !% Computes the rotation curve at a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentVelocity
    double precision                         :: componentMass

    ! Set to zero by default.
    componentVelocity=0.0d0

    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       call Dark_Matter_Profile_Enclosed_Mass_Task(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
       if (componentMass > 0.0d0) componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass)/dsqrt(radius)
    end if
    return
  end subroutine Dark_Matter_Profile_Rotation_Curve_Task

  !# <densityTask>
  !#  <unitName>Dark_Matter_Profile_Density_Task</unitName>
  !# </densityTask>
  subroutine Dark_Matter_Profile_Density_Task(thisNode,positionSpherical,massType,componentType,componentDensity)
    !% Computes the density at a given position for a dark matter profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: positionSpherical(3)
    double precision, intent(out)            :: componentDensity
    
    componentDensity=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return

    componentDensity=Dark_Matter_Profile_Density(thisNode,positionSpherical(1))
    return
  end subroutine Dark_Matter_Profile_Density_Task

  !# <rotationCurveGradientTask>
  !#  <unitName>Dark_Matter_Profile_Rotation_Curve_Gradient_Task</unitName>
  !# </rotationCurveGradientTask>
  subroutine Dark_Matter_Profile_Rotation_Curve_Gradient_Task(thisNode,radius,massType,componentType,componentRotationCurveGradient)
    !% Computes the rotation curve gradient for the dark matter.
    use Tree_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentRotationCurveGradient
    double precision                         :: positionSpherical(3),componentMass,componentDensity

    ! Set to zero by default.
    componentRotationCurveGradient=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(     massType == massTypeAll      .or.      massType == massTypeDark         )) return
    if (radius <= 0.0d0) return
    positionSpherical=[radius,0.0d0,0.0d0]
    call Dark_Matter_Profile_Enclosed_Mass_Task(thisNode,radius           ,massType,componentType,weightByMass,weightIndexNull,componentMass   )
    call Dark_Matter_Profile_Density_Task      (thisNode,positionSpherical,massType,componentType                             ,componentDensity)
    if (componentMass ==0.0d0 .or. componentDensity == 0.0d0) return
    componentRotationCurveGradient = gravitationalConstantGalacticus    &
                 &                  *(-componentMass/radius**2          &
                 &                    +4.0d0*Pi*radius*componentDensity &
                 &                   )  
    return
  end subroutine Dark_Matter_Profile_Rotation_Curve_Gradient_Task

  !# <potentialTask>
  !#  <unitName>Dark_Matter_Profile_Potential_Task</unitName>
  !# </potentialTask>
  subroutine Dark_Matter_Profile_Potential_Task(thisNode,radius,componentType,massType,componentPotential)
    !% Return the potential due to dark matter.
    use Galactic_Structure_Options
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: componentType,massType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentPotential

    componentPotential=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return

    componentPotential=Dark_Matter_Profile_Potential(thisNode,radius)
    return
  end subroutine Dark_Matter_Profile_Potential_Task

end module Dark_Matter_Profiles
