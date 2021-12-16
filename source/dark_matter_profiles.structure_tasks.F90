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

!!{
Contains a module which implements structure tasks related to the dark matter halo density profile.
!!}

module Dark_Matter_Profile_Structure_Tasks
  !!{
  Implements structure tasks related to the dark matter halo density profile.
  !!}
  use :: Dark_Matter_Profiles, only : darkMatterProfileClass
  private
  public :: Dark_Matter_Profile_Enclosed_Mass_Task               , Dark_Matter_Profile_Density_Task                       , Dark_Matter_Profile_Rotation_Curve_Task        , Dark_Matter_Profile_Potential_Task               , &
       &    Dark_Matter_Profile_Rotation_Curve_Gradient_Task     , Dark_Matter_Profile_Acceleration_Task                  , Dark_Matter_Profile_Tidal_Tensor_Task          , Dark_Matter_Profile_Chandrasekhar_Integral_Task  , &
       &    Dark_Matter_Profile_Structure_Tasks_Thread_Initialize, Dark_Matter_Profile_Structure_Tasks_Thread_Uninitialize, Dark_Matter_Profile_Structure_Tasks_State_Store, Dark_Matter_Profile_Structure_Tasks_State_Restore

  class(darkMatterProfileClass), pointer  :: darkMatterProfile_
  !$omp threadprivate(darkMatterProfile_)

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Dark_Matter_Profile_Structure_Tasks_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Dark_Matter_Profile_Structure_Tasks_Thread_Initialize(parameters_)
    !!{
    Initializes the dark matter profile structure tasks module.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    !![
    <objectBuilder class="darkMatterProfile" name="darkMatterProfile_" source="parameters_"/>
    !!]
    return
  end subroutine Dark_Matter_Profile_Structure_Tasks_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Dark_Matter_Profile_Structure_Tasks_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Dark_Matter_Profile_Structure_Tasks_Thread_Uninitialize()
    !!{
    Uninitializes the dark matter profile structure tasks module.
    !!}
    implicit none

    !![
    <objectDestructor name="darkMatterProfile_"/>
    !!]
    return
  end subroutine Dark_Matter_Profile_Structure_Tasks_Thread_Uninitialize

  !![
  <enclosedMassTask>
   <unitName>Dark_Matter_Profile_Enclosed_Mass_Task</unitName>
  </enclosedMassTask>
  !!]
  double precision function Dark_Matter_Profile_Enclosed_Mass_Task(node,radius,componentType,massType,weightBy,weightIndex)
    !!{
    Computes the mass within a given radius for a dark matter profile.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll  , componentTypeDarkHalo , massTypeAll, massTypeDark, &
          &                                   radiusLarge       , weightByMass
    use :: Galacticus_Nodes          , only : nodeComponentBasic, treeNode
    implicit none
    type            (treeNode          ), intent(inout)           :: node
    integer                             , intent(in   )           :: componentType, massType   , &
         &                                                           weightBy     , weightIndex
    double precision                    , intent(in   )           :: radius
    class           (nodeComponentBasic)               , pointer  :: basic
    !$GLC attributes unused :: weightIndex

    Dark_Matter_Profile_Enclosed_Mass_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    if (.not.(weightBy      == weightByMass                                                )) return

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

  !![
  <accelerationTask>
   <unitName>Dark_Matter_Profile_Acceleration_Task</unitName>
  </accelerationTask>
  !!]
  function Dark_Matter_Profile_Acceleration_Task(node,positionCartesian,componentType,massType)
    !!{
    Computes the acceleration due to a dark matter profile.
    !!}
    use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes                , only : treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear                       , megaParsec
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                         , dimension(3) :: Dark_Matter_Profile_Acceleration_Task
    type            (treeNode), intent(inout)               :: node
    integer                   , intent(in   )               :: componentType                        , massType
    double precision          , intent(in   ), dimension(3) :: positionCartesian
    double precision                                        :: radius

    radius=sqrt(sum(positionCartesian**2))
    Dark_Matter_Profile_Acceleration_Task=-kilo                                                                                                    &
         &                                *gigaYear                                                                                                &
         &                                /megaParsec                                                                                              &
         &                                *gravitationalConstantGalacticus                                                                         &
         &                                *Dark_Matter_Profile_Enclosed_Mass_Task(node,radius,componentType,massType,weightByMass,weightIndexNull) &
         &                                *positionCartesian                                                                                       &
         &                                /radius**3
    return
  end function Dark_Matter_Profile_Acceleration_Task

  !![
  <chandrasekharIntegralTask>
   <unitName>Dark_Matter_Profile_Chandrasekhar_Integral_Task</unitName>
  </chandrasekharIntegralTask>
  !!]
  function Dark_Matter_Profile_Chandrasekhar_Integral_Task(node,positionCartesian,velocityCartesian,componentType,massType)
    !!{
    Computes the Chandrasekhar integral due to a dark matter profile.
    !!}
    use :: Galactic_Structure_Options, only : weightByMass, weightIndexNull
    use :: Galacticus_Nodes          , only : treeNode
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    double precision                         , dimension(3) :: Dark_Matter_Profile_Chandrasekhar_Integral_Task
    type            (treeNode), intent(inout)               :: node
    integer                   , intent(in   )               :: componentType                                         , massType
    double precision          , intent(in   ), dimension(3) :: positionCartesian                                     , velocityCartesian
    double precision                         , dimension(3) :: positionSpherical
    double precision          , parameter                   :: XvMaximum                                      =10.0d0
    double precision                                        :: radius                                                , velocity         , &
         &                                                     density                                               , xV

    radius                                          =  sqrt(sum(positionCartesian**2))
    velocity                                        =  sqrt(sum(velocityCartesian**2))
    positionSpherical                               =  [radius,0.0d0,0.0d0]
    density                                         =  Dark_Matter_Profile_Density_Task(node,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
    xV                                              = +                         velocity                        &
         &                                            /darkMatterProfile_%radialVelocityDispersion(node,radius) &
         &                                            /sqrt(2.0d0)
    Dark_Matter_Profile_Chandrasekhar_Integral_Task = -density              &
         &                                            *velocityCartesian    &
         &                                            /velocity         **3
    if (Xv <= XvMaximum)                                                                                    &
         & Dark_Matter_Profile_Chandrasekhar_Integral_Task=+Dark_Matter_Profile_Chandrasekhar_Integral_Task &
         &                                                 *(                                               &
         &                                                   +erf ( xV   )                                  &
         &                                                   -2.0d0                                         &
         &                                                   *      xV                                      &
         &                                                   *exp (-xV**2)                                  &
         &                                                   /sqrt( Pi   )                                  &
         &                                                 )
    return
  end function Dark_Matter_Profile_Chandrasekhar_Integral_Task

  !![
  <tidalTensorTask>
   <unitName>Dark_Matter_Profile_Tidal_Tensor_Task</unitName>
  </tidalTensorTask>
  !!]
  function Dark_Matter_Profile_Tidal_Tensor_Task(node,positionCartesian,componentType,massType)
    !!{
    Computes the tidalTensor due to a dark matter profile.
    !!}
    use :: Galactic_Structure_Options  , only : weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes            , only : treeNode
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Tensors                     , only : tensorRank2Dimension3Symmetric , tensorIdentityR2D3Sym, assignment(=), operator(*)
    use :: Vectors                     , only : Vector_Outer_Product
    implicit none
    type            (tensorRank2Dimension3Symmetric)                              :: Dark_Matter_Profile_Tidal_Tensor_Task
    type            (treeNode                      ), intent(inout)               :: node
    integer                                         , intent(in   )               :: componentType                        , massType
    double precision                                , intent(in   ), dimension(3) :: positionCartesian
    double precision                                               , dimension(3) :: positionSpherical
    double precision                                                              :: radius                               , massEnclosed, &
         &                                                                           density
    type            (tensorRank2Dimension3Symmetric)                              :: positionTensor
    
    radius           =sqrt(sum(positionCartesian**2))
    positionSpherical=[radius,0.0d0,0.0d0]
    massEnclosed     =Dark_Matter_Profile_Enclosed_Mass_Task(node,radius           ,componentType,massType,weightByMass,weightIndexNull)
    density          =Dark_Matter_Profile_Density_Task      (node,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
    positionTensor   =Vector_Outer_Product                  (     positionCartesian,symmetrize=.true.                                  )
    Dark_Matter_Profile_Tidal_Tensor_Task=+gravitationalConstantGalacticus                           &
         &                                *(                                                         &
         &                                  -(massEnclosed         /radius**3)*tensorIdentityR2D3Sym &
         &                                  +(massEnclosed*3.0d0   /radius**5)*positionTensor        &
         &                                  -(density     *4.0d0*Pi/radius**2)*positionTensor        &
         &                                )
    return
  end function Dark_Matter_Profile_Tidal_Tensor_Task

  !![
  <rotationCurveTask>
   <unitName>Dark_Matter_Profile_Rotation_Curve_Task</unitName>
  </rotationCurveTask>
  !!]
  double precision function Dark_Matter_Profile_Rotation_Curve_Task(node,radius,componentType,massType)
    !!{
    Computes the rotation curve at a given radius for a dark matter profile.
    !!}
    use :: Galactic_Structure_Options  , only : weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes            , only : treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
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

  !![
  <densityTask>
   <unitName>Dark_Matter_Profile_Density_Task</unitName>
  </densityTask>
  !!]
  double precision function Dark_Matter_Profile_Density_Task(node,positionSpherical,componentType,massType,weightBy,weightIndex)
    !!{
    Computes the density at a given position for a dark matter profile.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll, componentTypeDarkHalo , massTypeAll, massTypeDark, &
          &                                   weightByMass
    use :: Galacticus_Nodes          , only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    integer                   , intent(in   ) :: componentType       , massType   , &
         &                                       weightBy            , weightIndex
    double precision          , intent(in   ) :: positionSpherical(3)
    !$GLC attributes unused :: weightIndex

    ! Return zero if the component and mass type is not matched.
    Dark_Matter_Profile_Density_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    if (.not.(weightBy      == weightByMass                                                )) return
    ! Compute the density
    Dark_Matter_Profile_Density_Task=darkMatterProfile_%density(node,positionSpherical(1))
    return
  end function Dark_Matter_Profile_Density_Task

  !![
  <rotationCurveGradientTask>
   <unitName>Dark_Matter_Profile_Rotation_Curve_Gradient_Task</unitName>
  </rotationCurveGradientTask>
  !!]
  double precision function Dark_Matter_Profile_Rotation_Curve_Gradient_Task(node,radius,componentType,massType)
    !!{
    Computes the rotation curve gradient for the dark matter.
    !!}
    use :: Galactic_Structure_Options  , only : componentTypeAll               , componentTypeDarkHalo, massTypeAll, massTypeDark, &
          &                                     weightByMass                   , weightIndexNull
    use :: Galacticus_Nodes            , only : treeNode
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
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

  !![
  <potentialTask>
   <unitName>Dark_Matter_Profile_Potential_Task</unitName>
  </potentialTask>
  !!]
  double precision function Dark_Matter_Profile_Potential_Task(node,radius,componentType,massType,status)
    !!{
    Return the potential due to dark matter.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll         , componentTypeDarkHalo , massTypeAll, massTypeDark, &
          &                                   structureErrorCodeSuccess
    use :: Galacticus_Error          , only : Galacticus_Error_Report
    use :: Galacticus_Nodes          , only : treeNode
    implicit none
    type            (treeNode), intent(inout)           :: node
    integer                   , intent(in   )           :: componentType, massType
    double precision          , intent(in   )           :: radius
    integer                   , intent(inout), optional :: status
    integer                                             :: statusLocal

    Dark_Matter_Profile_Potential_Task=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDarkHalo)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeDark         )) return
    Dark_Matter_Profile_Potential_Task=darkMatterProfile_%potential(node,radius,statusLocal)
    if (present(status).and.statusLocal /= structureErrorCodeSuccess) status=structureErrorCodeSuccess
    return
  end function Dark_Matter_Profile_Potential_Task

  !![
  <galacticusStateStoreTask>
   <unitName>Dark_Matter_Profile_Structure_Tasks_State_Store</unitName>
  </galacticusStateStoreTask>
  !!]
  subroutine Dark_Matter_Profile_Structure_Tasks_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDisk -> verySimple',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterProfile_"/>
    !!]
    return
  end subroutine Dark_Matter_Profile_Structure_Tasks_State_Store

  !![
  <galacticusStateRetrieveTask>
   <unitName>Dark_Matter_Profile_Structure_Tasks_State_Restore</unitName>
  </galacticusStateRetrieveTask>
  !!]
  subroutine Dark_Matter_Profile_Structure_Tasks_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDisk -> verySimple',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterProfile_"/>
    !!]
    return
  end subroutine Dark_Matter_Profile_Structure_Tasks_State_Restore

end module Dark_Matter_Profile_Structure_Tasks
