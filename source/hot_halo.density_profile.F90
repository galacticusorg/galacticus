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

!% Contains a module that implements calculations of the hot halo gas density profile.

module Hot_Halo_Density_Profile
  !% Implements calculations of the hot halo gas density profile.
  use ISO_Varying_String
  use Galacticus_Nodes
  !# <include directive="hotHaloDensityMethod" type="moduleUse">
  include 'hot_halo.density_profile.modules.inc'
  !# </include>
  implicit none
  private
  public :: Hot_Halo_Density, Hot_Halo_Density_Log_Slope, Hot_Halo_Enclosed_Mass, Hot_Halo_Profile_Density_Task,&
       & Hot_Halo_Profile_Rotation_Curve_Task, Hot_Halo_Profile_Enclosed_Mass_Task, Hot_Halo_Profile_Rotation_Normalization,&
       & Hot_Halo_Profile_Rotation_Curve_Gradient_Task,Hot_Halo_Profile_Radial_Moment
  ! Flag to indicate if this module has been initialized.
  logical                                           :: hotHaloDensityInitialized     =.false.

  ! Name of cooling time available method used.
  type     (varying_string               )          :: hotHaloDensityMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Hot_Halo_Density_Get_Template), pointer :: Hot_Halo_Density_Get          =>null()
  procedure(Hot_Halo_Density_Get_Template), pointer :: Hot_Halo_Density_Log_Slope_Get=>null()
  procedure(Hot_Halo_Density_Get_Template), pointer :: Hot_Halo_Enclosed_Mass_Get    =>null()
  abstract interface
     double precision function Hot_Halo_Density_Get_Template(thisNode,radius)
       import treeNode
       type            (treeNode), intent(inout), pointer :: thisNode
       double precision          , intent(in   )          :: radius
     end function Hot_Halo_Density_Get_Template
  end interface
  procedure(Hot_Halo_Profile_Rotation_Normalization_Template), pointer :: Hot_Halo_Profile_Rotation_Normalization_Get=>null()
  abstract interface
     double precision function Hot_Halo_Profile_Rotation_Normalization_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Hot_Halo_Profile_Rotation_Normalization_Template
  end interface
  procedure(Hot_Halo_Profile_Radial_Moment), pointer :: Hot_Halo_Profile_Radial_Moment_Get

contains

  subroutine Hot_Halo_Density_Initialize
    !% Initialize the hot halo density profile module.
    use Input_Parameters
    use Galacticus_Error
    implicit none

    ! Initialize if necessary.
    if (.not.hotHaloDensityInitialized) then
       !$omp critical(Hot_Halo_Density_Initialization)
       if (.not.hotHaloDensityInitialized) then
          ! Get the cooling time available method parameter.
          !@ <inputParameter>
          !@   <name>hotHaloDensityMethod</name>
          !@   <defaultValue>coredIsothermal</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of the hot halo density profile.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloDensityMethod',hotHaloDensityMethod,defaultValue='coredIsothermal')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="hotHaloDensityMethod" type="functionCall" functionType="void">
          !#  <functionArgs>hotHaloDensityMethod,Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get,Hot_Halo_Enclosed_Mass_Get,Hot_Halo_Profile_Rotation_Normalization_Get,Hot_Halo_Profile_Radial_Moment_Get</functionArgs>
          include 'hot_halo.density_profile.inc'
          !# </include>
          if     (.not.(     associated(Hot_Halo_Density_Get                       )                    &
               &        .and.associated(Hot_Halo_Density_Log_Slope_Get             )                    &
               &        .and.associated(Hot_Halo_Enclosed_Mass_Get                 )                    &
               &        .and.associated(Hot_Halo_Profile_Rotation_Normalization_Get)                    &
               &        .and.associated(Hot_Halo_Profile_Radial_Moment_Get         )                    &
               &       )                                                                                &
               & )                                                                                      &
               & call Galacticus_Error_Report(                                                          &
               &                               'Hot_Halo_Density_Initialize'                            &
               &                              ,'method '//char(hotHaloDensityMethod)//' is unrecognized'&
               &                             )
          hotHaloDensityInitialized=.true.
       end if
       !$omp end critical(Hot_Halo_Density_Initialization)
    end if
    return
  end subroutine Hot_Halo_Density_Initialize

  double precision function Hot_Halo_Density(thisNode,radius)
    !% Return the density of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Density=Hot_Halo_Density_Get(thisNode,radius)

    return
  end function Hot_Halo_Density

  double precision function Hot_Halo_Density_Log_Slope(thisNode,radius)
    !% Return the density of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Density_Log_Slope=Hot_Halo_Density_Log_Slope_Get(thisNode,radius)

    return
  end function Hot_Halo_Density_Log_Slope

  double precision function Hot_Halo_Enclosed_Mass(thisNode,radius)
    !% Return the enclosed mass of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Enclosed_Mass=Hot_Halo_Enclosed_Mass_Get(thisNode,radius)

    return
  end function Hot_Halo_Enclosed_Mass

  !# <enclosedMassTask>
  !#  <unitName>Hot_Halo_Profile_Enclosed_Mass_Task</unitName>
  !# </enclosedMassTask>
  double precision function Hot_Halo_Profile_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the mass within a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Cosmological_Parameters
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType, massType, weightBy, weightIndex
    double precision          , intent(in   )           :: radius
    logical                   , intent(in   ), optional :: haloLoaded

    ! Return zero mass if the requested mass type or component is not matched.
    Hot_Halo_Profile_Enclosed_Mass_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return

    ! Return the enclosed mass.
    Hot_Halo_Profile_Enclosed_Mass_Task=max(Hot_Halo_Enclosed_Mass(thisNode,radius),0.0d0)
    return
  end function Hot_Halo_Profile_Enclosed_Mass_Task

  !# <rotationCurveTask>
  !#  <unitName>Hot_Halo_Profile_Rotation_Curve_Task</unitName>
  !# </rotationCurveTask>
  double precision function Hot_Halo_Profile_Rotation_Curve_Task(thisNode,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve at a given radius for the hot halo density profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType, massType
    double precision          , intent(in   )           :: radius
    logical                   , intent(in   ), optional :: haloLoaded
    double precision                                    :: componentMass

    ! Set to zero by default.
    Hot_Halo_Profile_Rotation_Curve_Task=0.0d0

    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Hot_Halo_Profile_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
       if (componentMass > 0.0d0) Hot_Halo_Profile_Rotation_Curve_Task=sqrt(gravitationalConstantGalacticus*componentMass)/sqrt(radius)
    end if
    return
  end function Hot_Halo_Profile_Rotation_Curve_Task

  !# <rotationCurveGradientTask>
  !#  <unitName>Hot_Halo_Profile_Rotation_Curve_Gradient_Task</unitName>
  !# </rotationCurveGradientTask>
  double precision function Hot_Halo_Profile_Rotation_Curve_Gradient_Task(thisNode,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve gradient at a given radius for the hot halo density profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType   , massType
    double precision          , intent(in   )           :: radius
    logical                   , intent(in   ), optional :: haloLoaded
    double precision                                    :: componentDensity, componentMass

    ! Set to zero by default.
    Hot_Halo_Profile_Rotation_Curve_Gradient_Task=0.0d0

    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Hot_Halo_Profile_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
       if (componentMass > 0.0d0) then
          componentDensity=Hot_Halo_Density(thisNode,radius)
          Hot_Halo_Profile_Rotation_Curve_Gradient_Task=gravitationalConstantGalacticus*(-componentMass/radius**2+4.0d0*Pi*radius&
               &*componentDensity)
       end if
    end if
    return
  end function Hot_Halo_Profile_Rotation_Curve_Gradient_Task

  !# <densityTask>
  !#  <unitName>Hot_Halo_Profile_Density_Task</unitName>
  !# </densityTask>
  double precision function Hot_Halo_Profile_Density_Task(thisNode,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the density at a given position for a dark matter profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType       , massType, weightBy, &
         &                                                 weightIndex
    double precision          , intent(in   )           :: positionSpherical(3)
    logical                   , intent(in   ), optional :: haloLoaded

    Hot_Halo_Profile_Density_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return

    Hot_Halo_Profile_Density_Task=max(Hot_Halo_Density(thisNode,positionSpherical(1)),0.0d0)
    return
  end function Hot_Halo_Profile_Density_Task

  double precision function Hot_Halo_Profile_Rotation_Normalization(thisNode)
    !% Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in
    !% radius) for {\tt thisNode}. Specifically, the normalization, $A$, returned is such that $V_{\rm rot} = A J/M$.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the energy using the selected method.
    Hot_Halo_Profile_Rotation_Normalization=Hot_Halo_Profile_Rotation_Normalization_Get(thisNode)

    return
  end function Hot_Halo_Profile_Rotation_Normalization

  double precision function Hot_Halo_Profile_Radial_Moment(thisNode,moment,radius)
    !% Returns a radial {\tt moment} of the hot gas profile in {\tt thisNode} to the specified {\tt radius}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: moment  , radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize()
    ! Get the energy using the selected method.
    Hot_Halo_Profile_Radial_Moment=Hot_Halo_Profile_Radial_Moment_Get(thisNode,moment,radius)
    return
  end function Hot_Halo_Profile_Radial_Moment

end module Hot_Halo_Density_Profile
