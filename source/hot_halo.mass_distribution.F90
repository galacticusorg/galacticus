!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which provides an object that implements hot halo mass distributions.

module Hot_Halo_Mass_Distributions
  !% Provides an object that implements hot halo density profiles.
  use ISO_Varying_String
  use Mass_Distributions
  use Galacticus_Nodes
  private
  public :: Hot_Halo_Mass_Distribution_Density_Task      , Hot_Halo_Mass_Distribution_Rotation_Curve_Task         , &
       &    Hot_Halo_Mass_Distribution_Enclosed_Mass_Task, Hot_Halo_Mass_Distribution_Rotation_Curve_Gradient_Task

  !# <include directive="hotHaloMassDistribution" type="function" >
  !#  <descriptiveName>Hot Halo Mass Distributions</descriptiveName>
  !#  <description>Object implementing hot halo mass distributions.</description>
  !#  <default>betaProfile</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="density" >
  !#   <description>Return the density of the hot halo at the given {\tt radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: radius</argument>
  !#  </method>
  !#  <method name="densityLogSlope" >
  !#   <description>Return the logarithmic slope of the density of the hot halo at the given {\tt radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: radius</argument>
  !#  </method>
  !#  <method name="enclosedMass" >
  !#   <description>Return the mass enclosed in the hot halo at the given {\tt radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: radius</argument>
  !#  </method>
  !#  <method name="radialMoment" >
  !#   <description>Return the density of the hot halo at the given {\tt radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), pointer :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: moment, radius</argument>
  !#  </method>
  !#  <method name="rotationNormalization" >
  !#   <description>Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in radius) for {\tt node}. Specifically, the normalization, $A$, returned is such that $V_{\rm rot} = A J/M$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout), pointer :: node</argument>
  !#  </method>
  include 'hotHaloMassDistribution.type.inc'
  !# </include>

  !# <enclosedMassTask>
  !#  <unitName>Hot_Halo_Mass_Distribution_Enclosed_Mass_Task</unitName>
  !# </enclosedMassTask>
  double precision function Hot_Halo_Mass_Distribution_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the mass within a given radius for a dark matter profile.
    use Galactic_Structure_Options
    implicit none
    type            (treeNode                    ), intent(inout), pointer  :: thisNode
    integer                                       , intent(in   )           :: componentType, massType, weightBy, weightIndex
    double precision                              , intent(in   )           :: radius
    logical                                       , intent(in   ), optional :: haloLoaded
    class           (hotHaloMassDistributionClass)               , pointer  :: hotHalo

    ! Return zero mass if the requested mass type or component is not matched.
    Hot_Halo_Mass_Distribution_Enclosed_Mass_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    ! Return the enclosed mass.
    hotHalo => hotHaloMassDistribution()
    Hot_Halo_Mass_Distribution_Enclosed_Mass_Task=max(hotHalo%enclosedMass(thisNode,radius),0.0d0)
    return
  end function Hot_Halo_Mass_Distribution_Enclosed_Mass_Task

  !# <rotationCurveTask>
  !#  <unitName>Hot_Halo_Mass_Distribution_Rotation_Curve_Task</unitName>
  !# </rotationCurveTask>
  double precision function Hot_Halo_Mass_Distribution_Rotation_Curve_Task(thisNode,radius,componentType,massType,haloLoaded)
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
    Hot_Halo_Mass_Distribution_Rotation_Curve_Task=0.0d0

    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Hot_Halo_Mass_Distribution_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
       if (componentMass > 0.0d0) Hot_Halo_Mass_Distribution_Rotation_Curve_Task=sqrt(gravitationalConstantGalacticus*componentMass)/sqrt(radius)
    end if
    return
  end function Hot_Halo_Mass_Distribution_Rotation_Curve_Task

  !# <rotationCurveGradientTask>
  !#  <unitName>Hot_Halo_Mass_Distribution_Rotation_Curve_Gradient_Task</unitName>
  !# </rotationCurveGradientTask>
  double precision function Hot_Halo_Mass_Distribution_Rotation_Curve_Gradient_Task(thisNode,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve gradient at a given radius for the hot halo density profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    type            (treeNode                    ), intent(inout), pointer  :: thisNode
    integer                                       , intent(in   )           :: componentType   , massType
    double precision                              , intent(in   )           :: radius
    logical                                       , intent(in   ), optional :: haloLoaded
    class           (hotHaloMassDistributionClass)               , pointer  :: hotHalo
    double precision                                                        :: componentDensity, componentMass

    ! Set to zero by default.
    Hot_Halo_Mass_Distribution_Rotation_Curve_Gradient_Task=0.0d0

    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Hot_Halo_Mass_Distribution_Enclosed_Mass_Task(thisNode,radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
       if (componentMass > 0.0d0) then
          hotHalo => hotHaloMassDistribution()
          componentDensity=hotHalo%density(thisNode,radius)
          Hot_Halo_Mass_Distribution_Rotation_Curve_Gradient_Task=gravitationalConstantGalacticus*(-componentMass/radius**2+4.0d0*Pi*radius&
               &*componentDensity)
       end if
    end if
    return
  end function Hot_Halo_Mass_Distribution_Rotation_Curve_Gradient_Task

  !# <densityTask>
  !#  <unitName>Hot_Halo_Mass_Distribution_Density_Task</unitName>
  !# </densityTask>
  double precision function Hot_Halo_Mass_Distribution_Density_Task(thisNode,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the density at a given position for a dark matter profile.
    use Galactic_Structure_Options
    implicit none
    type            (treeNode                    ), intent(inout), pointer  :: thisNode
    integer                                       , intent(in   )           :: componentType       , massType, weightBy, &
         &                                                                     weightIndex
    double precision                              , intent(in   )           :: positionSpherical(3)
    class           (hotHaloMassDistributionClass)               , pointer  :: hotHalo
    logical                                       , intent(in   ), optional :: haloLoaded

    Hot_Halo_Mass_Distribution_Density_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    
    hotHalo => hotHaloMassDistribution()
    Hot_Halo_Mass_Distribution_Density_Task=max(hotHalo%density(thisNode,positionSpherical(1)),0.0d0)
    return
  end function Hot_Halo_Mass_Distribution_Density_Task

end module Hot_Halo_Mass_Distributions
