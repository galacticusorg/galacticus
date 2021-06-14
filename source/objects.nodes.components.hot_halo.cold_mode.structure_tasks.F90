!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements structure tasks for the cold mode hot halo component.

module Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks
  !% Implements structure tasks for the cold mode hot halo component.
  use :: Hot_Halo_Cold_Mode_Density_Core_Radii, only : hotHaloColdModeCoreRadiiClass
  use :: Mass_Distributions                   , only : massDistributionBetaProfile
  use :: Dark_Matter_Halo_Scales              , only : darkMatterHaloScaleClass
  implicit none
  private
  public :: Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task          , Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task, &
       &    Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task, Node_Component_Hot_Halo_Cold_Mode_Density_Task       , &
       &    Node_Component_Hot_Halo_Cold_Mode_Acceleration_Task           , Node_Component_Hot_Halo_Cold_Mode_Tidal_Tensor_Task  , &
       &    Node_Component_Hot_Halo_Cold_Mode_Chandrasekhar_Integral

  type (massDistributionBetaProfile  ), public          :: coldModeMassDistribution
  class(hotHaloColdModeCoreRadiiClass), public, pointer :: hotHaloColdModeCoreRadii_
  class(darkMatterHaloScaleClass     ), public, pointer :: darkMatterHaloScale_
  !$omp threadprivate(coldModeMassDistribution,hotHaloColdModeCoreRadii_,darkMatterHaloScale_)

contains

  !# <enclosedMassTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task</unitName>
  !# </enclosedMassTask>
  double precision function Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task(node,radius,componentType,massType,weightBy,weightIndex)
    !% Computes the mass within a given radius for the cold mode hot halo component.
    use :: Galactic_Structure_Options, only : componentTypeAll       , componentTypeColdHalo, massTypeAll , massTypeBaryonic, &
          &                                   massTypeGaseous        , radiusLarge          , weightByMass
    use :: Galacticus_Nodes          , only : defaultHotHaloComponent, nodeComponentHotHalo , treeNode
    implicit none
    type            (treeNode            ), intent(inout)           :: node
    integer                               , intent(in   )           :: componentType, massType   , &
         &                                                             weightBy     , weightIndex
    double precision                      , intent(in   )           :: radius
    class           (nodeComponentHotHalo)               , pointer  :: hotHalo
    double precision                                                :: radiusOuter  , radiusCore
    !$GLC attributes unused :: weightIndex

    ! Return zero mass if the requested mass type or component is not matched.
    Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task=0.0d0
    if (.not.defaultHotHaloComponent%coldModeIsActive()                                                                     ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeColdHalo                                )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    ! Get the hot halo component.
    hotHalo => node   %hotHalo    ()
    ! Check for total mass request.
    if (radius >= radiusLarge) then
       Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task=hotHalo%massCold()
       return
    end if
    ! Get the outer radius.
    radiusOuter =  hotHalo%outerRadius()
    if (radiusOuter <= 0.0d0) return
    ! Compute the enclosed mass.
    ! Find the scale length of the cold mode halo.
    radiusCore=hotHaloColdModeCoreRadii_%radius(node)
    ! Initialize the mass profile
    coldModeMassDistribution=massDistributionBetaProfile(beta=2.0d0/3.0d0,coreRadius=radiusCore,mass=hotHalo%massCold(),outerRadius=hotHalo%outerRadius())
    ! Compute the enclosed mass.
    Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task=coldModeMassDistribution%massEnclosedBySphere(radius)
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task

  !# <accelerationTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Acceleration_Task</unitName>
  !# </accelerationTask>
  function Node_Component_Hot_Halo_Cold_Mode_Acceleration_Task(node,positionCartesian,componentType,massType)
    !% Computes the acceleration due to a dark matter profile.
    use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes                , only : treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear                       , megaParsec
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                         , dimension(3) :: Node_Component_Hot_Halo_Cold_Mode_Acceleration_Task
    type            (treeNode), intent(inout)               :: node
    integer                   , intent(in   )               :: componentType                                      , massType
    double precision          , intent(in   ), dimension(3) :: positionCartesian
    double precision                                        :: radius

    radius                                             =+sqrt(sum(positionCartesian**2))
    Node_Component_Hot_Halo_Cold_Mode_Acceleration_Task=-kilo                                                                                                                 &
         &                                              *gigaYear                                                                                                               &
         &                                              /megaParsec                                                                                                            &
         &                                              *gravitationalConstantGalacticus                                                                                       &
         &                                              *Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task(node,radius,componentType,massType,weightByMass,weightIndexNull) &
         &                                              *positionCartesian                                                                                                     &
         &                                              /radius**3
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Acceleration_Task

  !# <tidalTensorTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Tidal_Tensor_Task</unitName>
  !# </tidalTensorTask>
  function Node_Component_Hot_Halo_Cold_Mode_Tidal_Tensor_Task(node,positionCartesian,componentType,massType)
    !% Computes the tidalTensor due to the cold mode halo.
    use :: Galactic_Structure_Options  , only : weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes            , only : treeNode
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Tensors                     , only : tensorRank2Dimension3Symmetric , tensorIdentityR2D3Sym, assignment(=), operator(*)
    use :: Vectors                     , only : Vector_Outer_Product
    implicit none
    type            (tensorRank2Dimension3Symmetric)                              :: Node_Component_Hot_Halo_Cold_Mode_Tidal_Tensor_Task
    type            (treeNode                      ), intent(inout)               :: node
    integer                                         , intent(in   )               :: componentType                        , massType
    double precision                                , intent(in   ), dimension(3) :: positionCartesian
    double precision                                               , dimension(3) :: positionSpherical
    double precision                                                              :: radius                               , massEnclosed, &
         &                                                                           density
    type            (tensorRank2Dimension3Symmetric)                              :: positionTensor
    
    radius                                             =sqrt(sum(positionCartesian**2))
    positionSpherical                                  =[radius,0.0d0,0.0d0]
    massEnclosed                                       =Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task(node,radius           ,componentType,massType,weightByMass,weightIndexNull)
    density                                            =Node_Component_Hot_Halo_Cold_Mode_Density_Task      (node,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
    positionTensor                                     =Vector_Outer_Product                                (     positionCartesian,symmetrize=.true.                                  )
    Node_Component_Hot_Halo_Cold_Mode_Tidal_Tensor_Task=+gravitationalConstantGalacticus                           &
         &                                              *(                                                         &
         &                                                -(massEnclosed         /radius**3)*tensorIdentityR2D3Sym &
         &                                                +(massEnclosed*3.0d0   /radius**5)*positionTensor        &
         &                                                -(density     *4.0d0*Pi/radius**2)*positionTensor        &
         &                                               )
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Tidal_Tensor_Task

  !# <rotationCurveTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task</unitName>
  !# </rotationCurveTask>
  double precision function Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task(node,radius,componentType,massType)
    !% Computes the rotation curve at a given radius for the hot halo density profile.
    use :: Galactic_Structure_Options  , only : weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes            , only : treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode), intent(inout)           :: node
    integer                   , intent(in   )           :: componentType, massType
    double precision          , intent(in   )           :: radius
    double precision                                    :: componentMass

    ! Set to zero by default.
    Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task=0.0d0
    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task(node,radius,componentType,massType,weightByMass,weightIndexNull)
       if (componentMass > 0.0d0) Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task=sqrt(gravitationalConstantGalacticus*componentMass)/sqrt(radius)
    end if
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Task

  !# <rotationCurveGradientTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task</unitName>
  !# </rotationCurveGradientTask>
  double precision function Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task(node,radius,componentType,massType)
    !% Computes the rotation curve gradient at a given radius for the hot halo density profile.
    use :: Galactic_Structure_Options  , only : weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes            , only : treeNode
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode), intent(inout) :: node
    integer                   , intent(in   ) :: componentType   , massType
    double precision          , intent(in   ) :: radius
    double precision                          :: componentDensity, componentMass

    ! Set to zero by default.
    Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task=0.0d0
    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=Node_Component_Hot_Halo_Cold_Mode_Enclosed_Mass_Task(node,radius,componentType,massType,weightByMass,weightIndexNull)
       if (componentMass > 0.0d0) then
          componentDensity=Node_Component_Hot_Halo_Cold_Mode_Density_Task(node,[radius,0.0d0,0.0d0],componentType,massType,weightByMass,weightIndexNull)
          Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task=gravitationalConstantGalacticus*(-componentMass/radius**2+4.0d0*Pi*radius&
               &*componentDensity)
       end if
    end if
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Rotation_Curve_Gradient_Task

  !# <densityTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Density_Task</unitName>
  !# </densityTask>
  double precision function Node_Component_Hot_Halo_Cold_Mode_Density_Task(node,positionSpherical,componentType,massType,weightBy,weightIndex)
    !% Computes the density at a given position for a dark matter profile.
    use :: Coordinates               , only : assignment(=)          , coordinateSpherical
    use :: Galactic_Structure_Options, only : componentTypeAll       , componentTypeColdHalo, massTypeAll, massTypeBaryonic, &
          &                                   massTypeGaseous        , weightByMass
    use :: Galacticus_Nodes          , only : defaultHotHaloComponent, nodeComponentHotHalo , treeNode
    implicit none
    type            (treeNode            ), intent(inout)           :: node
    integer                               , intent(in   )           :: componentType       , massType   , &
         &                                                             weightBy            , weightIndex
    double precision                      , intent(in   )           :: positionSpherical(3)
    class           (nodeComponentHotHalo)               , pointer  :: hotHalo
    type            (coordinateSpherical )                          :: position
    double precision                                                :: radiusOuter         , radiusCore
    !$GLC attributes unused :: weightIndex

    Node_Component_Hot_Halo_Cold_Mode_Density_Task=0.0d0
    if (.not.defaultHotHaloComponent%coldModeIsActive()                                                                     ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeColdHalo                                )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Get the outer radius.
    radiusOuter =  hotHalo%outerRadius()
    if (radiusOuter <= 0.0d0) return
    ! Compute the enclosed mass.
    ! Find the scale length of the cold mode halo.
    radiusCore=hotHaloColdModeCoreRadii_%radius(node)
    ! Initialize the mass profile
    coldModeMassDistribution=massDistributionBetaProfile(beta=2.0d0/3.0d0,coreRadius=radiusCore,mass=hotHalo%massCold(),outerRadius=hotHalo%outerRadius())
    ! Compute the density.
    position=[positionSpherical(1)/radiusCore,0.0d0,0.0d0]
    Node_Component_Hot_Halo_Cold_Mode_Density_Task=coldModeMassDistribution%density(position)
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Density_Task

  !# <chandrasekharIntegralTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Chandrasekhar_Integral</unitName>
  !# </chandrasekharIntegralTask>
  function Node_Component_Hot_Halo_Cold_Mode_Chandrasekhar_Integral(node,positionCartesian,velocityCartesian,componentType,massType)
    !% Computes the Chandrasekhar integral due to a dark matter profile.
    use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale, darkMatterHaloScaleClass
    use :: Galactic_Structure_Options, only : weightByMass       , weightIndexNull
    use :: Galacticus_Nodes          , only : treeNode
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    double precision                         , dimension(3) :: Node_Component_Hot_Halo_Cold_Mode_Chandrasekhar_Integral
    type            (treeNode), intent(inout)               :: node
    integer                   , intent(in   )               :: componentType                                                  , massType
    double precision          , intent(in   ), dimension(3) :: positionCartesian                                              , velocityCartesian
    double precision                         , dimension(3) :: positionSpherical
    double precision          , parameter                   :: XvMaximum                                               =10.0d0
    double precision                                        :: radius                                                         , velocity         , &
         &                                                     density                                                        , xV
    
    radius                                                   =  sqrt(sum(positionCartesian**2))
    velocity                                                 =  sqrt(sum(velocityCartesian**2))
    positionSpherical                                        =  [radius,0.0d0,0.0d0]
    density                                                  =  Node_Component_Hot_Halo_Cold_Mode_Density_Task(node,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
    xV                                                       = +                           velocity       &
         &                                                     /darkMatterHaloScale_%virialVelocity(node) &
         &                                                     /sqrt(2.0d0)
    Node_Component_Hot_Halo_Cold_Mode_Chandrasekhar_Integral = -density              &
         &                                                     *velocityCartesian    &
         &                                                     /velocity         **3
    if (Xv <= XvMaximum)                                                                                                      &
         & Node_Component_Hot_Halo_Cold_Mode_Chandrasekhar_Integral=+Node_Component_Hot_Halo_Cold_Mode_Chandrasekhar_Integral &
         &                                                          *(                                                        &
         &                                                            +erf ( xV   )                                           &
         &                                                            -2.0d0                                                  &
         &                                                            *      xV                                               &
         &                                                            *exp (-xV**2)                                           &
         &                                                            /sqrt( Pi   )                                           &
         &                                                          )
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Chandrasekhar_Integral

end module Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks
