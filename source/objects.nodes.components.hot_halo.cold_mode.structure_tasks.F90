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

!% Contains a module which implements structure tasks for the cold mode hot halo component.

module Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks
  !% Implements structure tasks for the cold mode hot halo component.
  use Mass_Distributions
  implicit none
  private
  public :: Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task          , Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task, &
       &    Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task, Node_Component_Hot_Halo_Cold_Mode_Density_Task

  ! The mass distribution object.
  class  (massDistribution), pointer, public :: coldModeMassDistribution
  !$omp threadprivate(coldModeMassDistribution)

contains

  !# <enclosedMassTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task</unitName>
  !# </enclosedMassTask>
  double precision function Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the mass within a given radius for the cold mode hot halo component.
    use Galactic_Structure_Options
    use Galacticus_Error
    use Galacticus_Nodes
    use Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radii
    implicit none
    type            (treeNode            ), intent(inout), pointer  :: thisNode
    integer                               , intent(in   )           :: componentType, massType   , &
         &                                                             weightBy     , weightIndex
    double precision                      , intent(in   )           :: radius
    logical                               , intent(in   ), optional :: haloLoaded
    class           (nodeComponentHotHalo)               , pointer  :: thisHotHalo
    double precision                                                :: radiusOuter  , radiusCore

    ! Return zero mass if the requested mass type or component is not matched.
    Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task=0.0d0
    if (.not.defaultHotHaloComponent%coldModeIsActive()                                                                     ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeColdHalo                                )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    ! Get the hot halo component.
    thisHotHalo => thisNode   %hotHalo    ()
    ! Get the outer radius.
    radiusOuter =  thisHotHalo%outerRadius()
    if (radiusOuter <= 0.0d0) return
    ! Compute the enclosed mass.
    select type (coldModeMassDistribution)
    type is (massDistributionBetaProfile)
       ! Find the scale length of the cold mode halo.
       radiusCore  =  Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radius(thisNode)
       ! Initialize the mass profile
       call coldModeMassDistribution%initialize(beta=2.0d0/3.0d0,coreRadius=radiusCore,mass=thisHotHalo%massCold(),outerRadius=thisHotHalo%outerRadius())
       ! Compute the enclosed mass.
       Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task=coldModeMassDistribution%massEnclosedBySphere(radius)
    class default
       call Galacticus_Error_Report('Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task','unsupported mass distribution')
    end select
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task

  !# <rotationCurveTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task</unitName>
  !# </rotationCurveTask>
  double precision function Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task(thisNode,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve at a given radius for the hot halo density profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer  :: thisNode
    integer                   , intent(in   )           :: componentType, massType
    double precision          , intent(in   )           :: radius
    logical                   , intent(in   ), optional :: haloLoaded
    double precision                                    :: componentMass

    ! Set to zero by default.
    Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task=0.0d0
    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
       if (componentMass > 0.0d0) Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task=sqrt(gravitationalConstantGalacticus*componentMass)/sqrt(radius)
    end if
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task

  !# <rotationCurveGradientTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task</unitName>
  !# </rotationCurveGradientTask>
  double precision function Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task(thisNode,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve gradient at a given radius for the hot halo density profile.
    use Galacticus_Nodes
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
    Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task=0.0d0
    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
       if (componentMass > 0.0d0) then
          componentDensity=Node_Component_Hot_Halo_Cold_Mode_Density_Task(thisNode,[radius,0.0d0,0.0d0],componentType,massType,weightByMass,weightIndexNull,haloLoaded)
          Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task=gravitationalConstantGalacticus*(-componentMass/radius**2+4.0d0*Pi*radius&
               &*componentDensity)
       end if
    end if
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task

  !# <densityTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Density_Task</unitName>
  !# </densityTask>
  double precision function Node_Component_Hot_Halo_Cold_Mode_Density_Task(thisNode,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the density at a given position for a dark matter profile.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Coordinates
    use Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radii
    use Galacticus_Error
    implicit none
    type            (treeNode            ), intent(inout), pointer  :: thisNode
    integer                               , intent(in   )           :: componentType       , massType   , &
         &                                                             weightBy            , weightIndex
    double precision                      , intent(in   )           :: positionSpherical(3)
    logical                               , intent(in   ), optional :: haloLoaded
    class           (nodeComponentHotHalo)               , pointer  :: thisHotHalo
    type            (coordinateSpherical )                          :: position
    double precision                                                :: radiusOuter         , radiusCore

    Node_Component_Hot_Halo_Cold_Mode_Density_Task=0.0d0
    if (.not.defaultHotHaloComponent%coldModeIsActive()                                                                     ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeColdHalo                                )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    ! Get the hot halo component.
    thisHotHalo => thisNode   %hotHalo    ()
    ! Get the outer radius.
    radiusOuter =  thisHotHalo%outerRadius()
    if (radiusOuter <= 0.0d0) return
    ! Compute the enclosed mass.
    select type (coldModeMassDistribution)
    type is (massDistributionBetaProfile)
       ! Find the scale length of the cold mode halo.
       radiusCore  =  Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radius(thisNode)
       ! Initialize the mass profile
       call coldModeMassDistribution%initialize(beta=2.0d0/3.0d0,coreRadius=radiusCore,mass=thisHotHalo%massCold(),outerRadius=thisHotHalo%outerRadius())
       ! Compute the density.
       position=[positionSpherical(1)/radiusCore,0.0d0,0.0d0]
       Node_Component_Hot_Halo_Cold_Mode_Density_Task=coldModeMassDistribution%density(position)
    class default
       call Galacticus_Error_Report('Node_Component_Hot_Halo_Cold_Mode_Density_Task','unsupported mass distribution')
    end select
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Density_Task

end module Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks
