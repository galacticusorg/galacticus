!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module that implements calculations of the hot halo gas density profile.

module Hot_Halo_Density_Profile
  !% Implements calculations of the hot halo gas density profile.
  use ISO_Varying_String
  use Tree_Nodes
  !# <include directive="hotHaloDensityMethod" type="moduleUse">
  include 'hot_halo.density_profile.modules.inc'
  !# </include>
  private
  public :: Hot_Halo_Density, Hot_Halo_Density_Log_Slope, Hot_Halo_Enclosed_Mass, Hot_Halo_Profile_Density_Task,&
       & Hot_Halo_Profile_Rotation_Curve_Task, Hot_Halo_Profile_Enclosed_Mass_Task, Hot_Halo_Profile_Rotation_Normalization

  ! Flag to indicate if this module has been initialized.  
  logical              :: hotHaloDensityInitialized=.false.

  ! Name of cooling time available method used.
  type(varying_string) :: hotHaloDensityMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Hot_Halo_Density_Get_Template), pointer :: Hot_Halo_Density_Get                        => null()
  procedure(Hot_Halo_Density_Get_Template), pointer :: Hot_Halo_Density_Log_Slope_Get              => null()
  procedure(Hot_Halo_Density_Get_Template), pointer :: Hot_Halo_Enclosed_Mass_Get                  => null()
  abstract interface
     double precision function Hot_Halo_Density_Get_Template(thisNode,radius)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: radius
     end function Hot_Halo_Density_Get_Template
  end interface
  procedure(Hot_Halo_Profile_Rotation_Normalization_Template), pointer :: Hot_Halo_Profile_Rotation_Normalization_Get => null()
  abstract interface
     double precision function Hot_Halo_Profile_Rotation_Normalization_Template(thisNode)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
     end function Hot_Halo_Profile_Rotation_Normalization_Template
  end interface

contains

  subroutine Hot_Halo_Density_Initialize
    !% Initialize the hot halo density profile module.
    use Input_Parameters
    use Galacticus_Error
    implicit none

    !$omp critical(Hot_Halo_Density_Initialization) 
    ! Initialize if necessary.
    if (.not.hotHaloDensityInitialized) then
       ! Get the cooling time available method parameter.
       !@ <inputParameter>
       !@   <name>hotHaloDensityMethod</name>
       !@   <defaultValue>cored isothermal</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of the hot halo density profile.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloDensityMethod',hotHaloDensityMethod,defaultValue='cored isothermal')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="hotHaloDensityMethod" type="code" action="subroutine">
       !#  <subroutineArgs>hotHaloDensityMethod,Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get,Hot_Halo_Enclosed_Mass_Get,Hot_Halo_Profile_Rotation_Normalization_Get</subroutineArgs>
       include 'hot_halo.density_profile.inc'
       !# </include>
       if     (.not.(     associated(Hot_Halo_Density_Get                       )                    &
            &        .and.associated(Hot_Halo_Density_Log_Slope_Get             )                    &
            &        .and.associated(Hot_Halo_Enclosed_Mass_Get                 )                    &
            &        .and.associated(Hot_Halo_Profile_Rotation_Normalization_Get)                    &
            &       )                                                                                &
            & )                                                                                      &
            & call Galacticus_Error_Report(                                                          &
            &                               'Hot_Halo_Density_Initialize'                            &
            &                              ,'method '//char(hotHaloDensityMethod)//' is unrecognized'&
            &                             )
       hotHaloDensityInitialized=.true.
    end if
    !$omp end critical(Hot_Halo_Density_Initialization) 
    return
  end subroutine Hot_Halo_Density_Initialize

  double precision function Hot_Halo_Density(thisNode,radius)
    !% Return the density of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Density=Hot_Halo_Density_Get(thisNode,radius)

    return
  end function Hot_Halo_Density

  double precision function Hot_Halo_Density_Log_Slope(thisNode,radius)
    !% Return the density of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Density_Log_Slope=Hot_Halo_Density_Log_Slope_Get(thisNode,radius)

    return
  end function Hot_Halo_Density_Log_Slope

  double precision function Hot_Halo_Enclosed_Mass(thisNode,radius)
    !% Return the enclosed mass of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Enclosed_Mass=Hot_Halo_Enclosed_Mass_Get(thisNode,radius)
    
    return
  end function Hot_Halo_Enclosed_Mass
  
  !# <enclosedMassTask>
  !#  <unitName>Hot_Halo_Profile_Enclosed_Mass_Task</unitName>
  !# </enclosedMassTask>
  subroutine Hot_Halo_Profile_Enclosed_Mass_Task(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Cosmological_Parameters
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType,weightBy,weightIndex
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass

    ! Return zero mass if the requested mass type or component is not matched.
    componentMass=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return

    ! Return the enclosed mass.
    componentMass=Hot_Halo_Enclosed_Mass(thisNode,radius)
    return
  end subroutine Hot_Halo_Profile_Enclosed_Mass_Task

  !# <rotationCurveTask>
  !#  <unitName>Hot_Halo_Profile_Rotation_Curve_Task</unitName>
  !# </rotationCurveTask>
  subroutine Hot_Halo_Profile_Rotation_Curve_Task(thisNode,radius,massType,componentType,componentVelocity)
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
       call Hot_Halo_Profile_Enclosed_Mass_Task(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
       if (componentMass > 0.0d0) componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass)/dsqrt(radius)
    end if
    return
  end subroutine Hot_Halo_Profile_Rotation_Curve_Task

  !# <densityTask>
  !#  <unitName>Hot_Halo_Profile_Density_Task</unitName>
  !# </densityTask>
  subroutine Hot_Halo_Profile_Density_Task(thisNode,positionSpherical,massType,componentType,componentDensity)
    !% Computes the density at a given position for a dark matter profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: positionSpherical(3)
    double precision, intent(out)            :: componentDensity
    
    componentDensity=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return

    componentDensity=Hot_Halo_Density(thisNode,positionSpherical(1))
    return
  end subroutine Hot_Halo_Profile_Density_Task

  double precision function Hot_Halo_Profile_Rotation_Normalization(thisNode)
    !% Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in
    !% radius) for {\tt thisNode}. Specifically, the normalization, $A$, returned is such that $V_{\rm rot} = A J/M$.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the energy using the selected method.
    Hot_Halo_Profile_Rotation_Normalization=Hot_Halo_Profile_Rotation_Normalization_Get(thisNode)

    return
  end function Hot_Halo_Profile_Rotation_Normalization

end module Hot_Halo_Density_Profile
