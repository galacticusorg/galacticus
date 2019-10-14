!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module which implements structure tasks related to the dark matter halo density profile.

module Dark_Matter_Profile_Structure_Tasks
  !% Implements structure tasks related to the dark matter halo density profile.
  use Dark_Matter_Profiles
  private
  public :: Dark_Matter_Profile_Enclosed_Mass_Task,Dark_Matter_Profile_Density_Task&
       &,Dark_Matter_Profile_Rotation_Curve_Task ,Dark_Matter_Profile_Potential_Task&
       &,Dark_Matter_Profile_Rotation_Curve_Gradient_Task

contains

  !# <enclosedMassTask>
  !#  <unitName>Dark_Matter_Profile_Enclosed_Mass_Task</unitName>
  !# </enclosedMassTask>
  double precision function Dark_Matter_Profile_Enclosed_Mass_Task(node,radius,componentType,massType,weightBy,weightIndex)
    !% Computes the mass within a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Galacticus_Nodes          , only : treeNode, nodeComponentBasic
    implicit none
    type            (treeNode              ), intent(inout)           :: node
    integer                                 , intent(in   )           :: componentType     , massType   , &
         &                                                               weightBy          , weightIndex
    double precision                        , intent(in   )           :: radius
    class           (nodeComponentBasic    )               , pointer  :: basic
    class           (darkMatterProfileClass)               , pointer  :: darkMatterProfile_
    !GCC$ attributes unused :: weightIndex

    Dark_Matter_Profile_Enclosed_Mass_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    if (.not.(weightBy      == weightByMass                                                )) return

    ! Get required objects.
    darkMatterProfile_ => darkMatterProfile()
    ! Test radius.
    if (radius >= radiusLarge) then
       ! Return the total mass of the halo in this case.
       basic => node%basic()
       Dark_Matter_Profile_Enclosed_Mass_Task=basic%mass()
    else if (radius <= 0.0d0) then
       ! Zero radius. Return zero mass.
       Dark_Matter_Profile_Enclosed_Mass_Task=0.0d0
    else       
       ! Return the mass within the radius.
       Dark_Matter_Profile_Enclosed_Mass_Task=darkMatterProfile_%enclosedMass(node,radius)
    end if
    return
  end function Dark_Matter_Profile_Enclosed_Mass_Task

  !# <rotationCurveTask>
  !#  <unitName>Dark_Matter_Profile_Rotation_Curve_Task</unitName>
  !# </rotationCurveTask>
  double precision function Dark_Matter_Profile_Rotation_Curve_Task(node,radius,componentType,massType)
    !% Computes the rotation curve at a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Galacticus_Nodes            , only : treeNode
    implicit none
    type            (treeNode), intent(inout)           :: node
    integer                   , intent(in   )           :: componentType, massType
    double precision          , intent(in   )           :: radius
    double precision                                    :: componentMass

    ! Set to zero by default.
    Dark_Matter_Profile_Rotation_Curve_Task=0.0d0

    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Dark_Matter_Profile_Enclosed_Mass_Task(node,radius,componentType,massType,weightByMass,weightIndexNull)
       if (componentMass > 0.0d0) Dark_Matter_Profile_Rotation_Curve_Task=sqrt(gravitationalConstantGalacticus*componentMass)&
            &/sqrt(radius)
    end if
    return
  end function Dark_Matter_Profile_Rotation_Curve_Task

  !# <densityTask>
  !#  <unitName>Dark_Matter_Profile_Density_Task</unitName>
  !# </densityTask>
  double precision function Dark_Matter_Profile_Density_Task(node,positionSpherical,componentType,massType,weightBy,weightIndex)
    !% Computes the density at a given position for a dark matter profile.
    use Galactic_Structure_Options
    use Galacticus_Nodes          , only : treeNode
    implicit none
    type            (treeNode              ), intent(inout)           :: node
    integer                                 , intent(in   )           :: componentType          , massType   , &
         &                                                               weightBy               , weightIndex
    double precision                        , intent(in   )           :: positionSpherical   (3)
    class           (darkMatterProfileClass)               , pointer  :: darkMatterProfile_
    !GCC$ attributes unused :: weightIndex
    
    ! Return zero if the component and mass type is not matched.
    Dark_Matter_Profile_Density_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    if (.not.(weightBy      == weightByMass                                                )) return
    ! Get required objects.
    darkMatterProfile_ => darkMatterProfile()
    ! Compute the density
    Dark_Matter_Profile_Density_Task=darkMatterProfile_%density(node,positionSpherical(1))
    return
  end function Dark_Matter_Profile_Density_Task

  !# <rotationCurveGradientTask>
  !#  <unitName>Dark_Matter_Profile_Rotation_Curve_Gradient_Task</unitName>
  !# </rotationCurveGradientTask>
  double precision function Dark_Matter_Profile_Rotation_Curve_Gradient_Task(node,radius,componentType,massType)
    !% Computes the rotation curve gradient for the dark matter.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    use Galacticus_Nodes            , only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    integer                   , intent(in   ) :: componentType   , massType
    double precision          , intent(in   ) :: radius
    double precision                          :: componentDensity, componentMass, positionSpherical(3)

    ! Set to zero by default.
    Dark_Matter_Profile_Rotation_Curve_Gradient_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(     massType == massTypeAll      .or.      massType == massTypeDark         )) return
    if (radius <= 0.0d0) return

    positionSpherical=[radius,0.0d0,0.0d0]
    componentMass   =Dark_Matter_Profile_Enclosed_Mass_Task(node,radius           ,componentType,massType,weightByMass,weightIndexNull)
    componentDensity=Dark_Matter_Profile_Density_Task      (node,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
    if (componentMass ==0.0d0 .or. componentDensity == 0.0d0) return
    Dark_Matter_Profile_Rotation_Curve_Gradient_Task=           &
         &                   gravitationalConstantGalacticus    &
         &                  *(-componentMass/radius**2          &
         &                    +4.0d0*Pi*radius*componentDensity &
         &                   )
    return
  end function Dark_Matter_Profile_Rotation_Curve_Gradient_Task

  !# <potentialTask>
  !#  <unitName>Dark_Matter_Profile_Potential_Task</unitName>
  !# </potentialTask>
  double precision function Dark_Matter_Profile_Potential_Task(node,radius,componentType,massType,status)
    !% Return the potential due to dark matter.
    use Galactic_Structure_Options
    use Galacticus_Error, only : Galacticus_Error_Report
    use Galacticus_Nodes          , only : treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer  :: node
    integer                                 , intent(in   )           :: componentType, massType
    double precision                        , intent(in   )           :: radius
    integer                                 , intent(inout), optional :: status
    class           (darkMatterProfileClass)               , pointer  :: darkMatterProfile_
    integer                                                           :: statusLocal

    Dark_Matter_Profile_Potential_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    darkMatterProfile_                 => darkMatterProfile           (                       )
    Dark_Matter_Profile_Potential_Task =  darkMatterProfile_%potential(node,radius,statusLocal)
    if (present(status).and.statusLocal /= structureErrorCodeSuccess) status=structureErrorCodeSuccess
    return
  end function Dark_Matter_Profile_Potential_Task

end module Dark_Matter_Profile_Structure_Tasks
