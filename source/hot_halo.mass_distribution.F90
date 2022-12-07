!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which provides a hot halo mass distribution class.
!!}

module Hot_Halo_Mass_Distributions
  !!{
  Provides an object which provides a hot halo mass distribution class.
  !!}
  use :: Galacticus_Nodes             , only : treeNode
  use :: Hot_Halo_Temperature_Profiles, only : hotHaloTemperatureProfileClass
  private
  public :: hotHaloMassDistributionDensity              , hotHaloMassDistributionRotationCurve          , &
       &    hotHaloMassDistributionEnclosedMass         , hotHaloMassDistributionRotationCurveGradient  , &
       &    hotHaloMassDistributionAcceleration         , hotHaloMassDistributionAccelerationTidalTensor, &
       &    hotHaloMassDistributionChandrasekharIntegral, hotHaloMassDistributionThreadInitialize       , &
       &    hotHaloMassDistributionThreadUninitialize   , hotHaloMassDistributionDefaultStateStore      , &
       &    hotHaloMassDistributionDefaultStateRestore  , hotHaloMassDistributionDensitySphericalAverage

  !![
  <functionClass>
   <name>hotHaloMassDistribution</name>
   <descriptiveName>Hot Halo Mass Distributions</descriptiveName>
   <description>
    Object implementing hot halo mass distributions.
   </description>
   <default>betaProfile</default>
   <method name="density" >
    <description>Return the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: radius</argument>
   </method>
   <method name="densityLogSlope" >
    <description>Return the logarithmic slope of the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   )          :: radius</argument>
   </method>
   <method name="enclosedMass" >
    <description>Return the mass enclosed in the hot halo at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: radius</argument>
   </method>
   <method name="radialMoment" >
    <description>Return the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: moment, radius</argument>
   </method>
   <method name="rotationNormalization" >
    <description>Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in radius) for {\normalfont \ttfamily node}. Specifically, the normalization, $A$, returned is such that $V_\mathrm{rot} = A J/M$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

  class(hotHaloMassDistributionClass  ), pointer :: hotHaloMassDistribution_
  class(hotHaloTemperatureProfileClass), pointer :: hotHaloTemperatureProfile_
  !$omp threadprivate(hotHaloMassDistribution_,hotHaloTemperatureProfile_)
  
  contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>hotHaloMassDistributionThreadInitialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine hotHaloMassDistributionThreadInitialize(parameters_)
    !!{
    Initializes the hot halo profile structure tasks module.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    !![
    <objectBuilder class="hotHaloMassDistribution"   name="hotHaloMassDistribution_"   source="parameters_"/>
    <objectBuilder class="hotHaloTemperatureProfile" name="hotHaloTemperatureProfile_" source="parameters_"/>
    !!]
    return
  end subroutine hotHaloMassDistributionThreadInitialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>hotHaloMassDistributionThreadUninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine hotHaloMassDistributionThreadUninitialize()
    !!{
    Uninitializes the hot halo profile structure tasks module.
    !!}
    implicit none

    !![
    <objectDestructor name="hotHaloMassDistribution_"  />
    <objectDestructor name="hotHaloTemperatureProfile_"/>
    !!]
    return
  end subroutine hotHaloMassDistributionThreadUninitialize

  !![
  <enclosedMassTask>
   <unitName>hotHaloMassDistributionEnclosedMass</unitName>
  </enclosedMassTask>
  !!]
  double precision function hotHaloMassDistributionEnclosedMass(node,radius,componentType,massType,weightBy,weightIndex)
    !!{
    Computes the mass within a given radius for a hot halo profile.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll       , componentTypeHotHalo   , massTypeAll , massTypeBaryonic            , &
         &                                    massTypeGaseous        , radiusLarge            , weightByMass, enumerationComponentTypeType, &
         &                                    enumerationMassTypeType, enumerationWeightByType
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo   , treeNode
    implicit none
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationComponentTypeType), intent(in   ) :: componentType
    type            (enumerationMassTypeType     ), intent(in   ) :: massType
    type            (enumerationWeightByType     ), intent(in   ) :: weightBy
    integer                                       , intent(in   ) :: weightIndex
    double precision                              , intent(in   ) :: radius
    class           (nodeComponentHotHalo        ), pointer       :: hotHalo
    !$GLC attributes unused :: weightIndex

    ! Return zero mass if the requested mass type or component is not matched.
    hotHaloMassDistributionEnclosedMass=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    ! Return the enclosed mass.
    if (radius >= radiusLarge) then
       hotHalo                             => node   %hotHalo()
       hotHaloMassDistributionEnclosedMass =  hotHalo%mass   ()
    else
       hotHaloMassDistributionEnclosedMass =  max(hotHaloMassDistribution_%enclosedMass(node,radius),0.0d0)
    end if
    return
  end function hotHaloMassDistributionEnclosedMass

  !![
  <accelerationTask>
   <unitName>hotHaloMassDistributionAcceleration</unitName>
  </accelerationTask>
  !!]
  function hotHaloMassDistributionAcceleration(node,positionCartesian,componentType,massType)
    !!{
    Computes the acceleration due to a hot halo profile.
    !!}
    use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull, enumerationComponentTypeType, enumerationMassTypeType
    use :: Galacticus_Nodes                , only : treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear                       , megaParsec
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                             , dimension(3) :: hotHaloMassDistributionAcceleration
    type            (treeNode                    ), intent(inout)               :: node
    type            (enumerationComponentTypeType), intent(in   )               :: componentType
    type            (enumerationMassTypeType     ), intent(in   )               :: massType
    double precision                              , intent(in   ), dimension(3) :: positionCartesian
    double precision                                                            :: radius

    radius                             =+sqrt(sum(positionCartesian**2))
    hotHaloMassDistributionAcceleration=-kilo                                                                                                 &
         &                              *gigaYear                                                                                             &
         &                              /megaParsec                                                                                           &
         &                              *gravitationalConstantGalacticus                                                                      &
         &                              *hotHaloMassDistributionEnclosedMass(node,radius,componentType,massType,weightByMass,weightIndexNull) &
         &                              *positionCartesian                                                                                    &
         &                              /radius**3
    return
  end function hotHaloMassDistributionAcceleration

  !![
  <tidalTensorTask>
   <unitName>hotHaloMassDistributionAccelerationTidalTensor</unitName>
  </tidalTensorTask>
  !!]
  function hotHaloMassDistributionAccelerationTidalTensor(node,positionCartesian,componentType,massType)
    !!{
    Computes the tidalTensor due to the cold mode halo.
    !!}
    use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull      , enumerationComponentTypeType, enumerationMassTypeType
    use :: Galacticus_Nodes                , only : treeNode
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Tensors                         , only : tensorRank2Dimension3Symmetric , tensorIdentityR2D3Sym, assignment(=)               , operator(*)
    use :: Vectors                         , only : Vector_Outer_Product
    implicit none
    type            (tensorRank2Dimension3Symmetric)                              :: hotHaloMassDistributionAccelerationTidalTensor
    type            (treeNode                      ), intent(inout)               :: node
    type            (enumerationComponentTypeType  ), intent(in   )               :: componentType
    type            (enumerationMassTypeType       ), intent(in   )               :: massType
    double precision                                , intent(in   ), dimension(3) :: positionCartesian
    double precision                                               , dimension(3) :: positionSpherical
    double precision                                                              :: radius                                        , massEnclosed, &
         &                                                                           density
    type            (tensorRank2Dimension3Symmetric)                              :: positionTensor
    
    radius           =sqrt(sum(positionCartesian**2))
    positionSpherical=[radius,0.0d0,0.0d0]
    massEnclosed     =hotHaloMassDistributionEnclosedMass(node,radius           ,componentType,massType,weightByMass,weightIndexNull)
    density          =hotHaloMassDistributionDensity     (node,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
    positionTensor   =Vector_Outer_Product               (     positionCartesian,symmetrize=.true.                                  )
    hotHaloMassDistributionAccelerationTidalTensor=+gravitationalConstantGalacticus                           &
         &                                         *(                                                         &
         &                                           -(massEnclosed         /radius**3)*tensorIdentityR2D3Sym &
         &                                           +(massEnclosed*3.0d0   /radius**5)*positionTensor        &
         &                                           -(density     *4.0d0*Pi/radius**2)*positionTensor        &
         &                                          )
    return
  end function hotHaloMassDistributionAccelerationTidalTensor

  !![
  <chandrasekharIntegralTask>
   <unitName>hotHaloMassDistributionChandrasekharIntegral</unitName>
  </chandrasekharIntegralTask>
  !!]
  function hotHaloMassDistributionChandrasekharIntegral(node,nodeSatellite,positionCartesian,velocityCartesian,componentType,massType)
    !!{
    Computes the Chandrasekhar integral due to the hot halo.
    !!}
    use :: Galactic_Structure_Options, only : weightByMass         , weightIndexNull, enumerationComponentTypeType, enumerationMassTypeType
    use :: Galacticus_Nodes          , only : treeNode
    use :: Ideal_Gases_Thermodynamics, only : Ideal_Gas_Sound_Speed
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    double precision                                             , dimension(3) :: hotHaloMassDistributionChandrasekharIntegral
    type            (treeNode                    ), intent(inout)               :: node                                                , nodeSatellite
    type            (enumerationComponentTypeType), intent(in   )               :: componentType
    type            (enumerationMassTypeType     ), intent(in   )               :: massType
    double precision                              , intent(in   ), dimension(3) :: positionCartesian                                   , velocityCartesian
    double precision                                             , dimension(3) :: positionSpherical
    double precision                              , parameter                   :: XvMaximum                                    =10.0d0
    double precision                                                            :: radius                                              , velocity         , &
         &                                                                         density                                             , xV
    !$GLC attributes unused :: radiusHalfMass, nodeSatellite
    
    radius                                      = sqrt(sum(positionCartesian**2))
    velocity                                    = sqrt(sum(velocityCartesian**2))
    positionSpherical                           = [radius,0.0d0,0.0d0]
    density                                     = hotHaloMassDistributionDensity(node,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
    xV                                          =+velocity                                                                       &
         &                                       /Ideal_Gas_Sound_Speed(hotHaloTemperatureProfile_%temperature(node,radius)) &
         &                                       /sqrt(2.0d0)
    hotHaloMassDistributionChandrasekharIntegral=-density              &
         &                                       *velocityCartesian    &
         &                                       /velocity         **3
    if (Xv <= XvMaximum)                                                                              &
         & hotHaloMassDistributionChandrasekharIntegral=+hotHaloMassDistributionChandrasekharIntegral &
         &                                                 *(                                         &
         &                                                   +erf ( xV   )                            &
         &                                                   -2.0d0                                   &
         &                                                   *      xV                                &
         &                                                   *exp (-xV**2)                            &
         &                                                   /sqrt( Pi   )                            &
         &                                                 )
    return
  end function hotHaloMassDistributionChandrasekharIntegral
  
  !![
  <rotationCurveTask>
   <unitName>hotHaloMassDistributionRotationCurve</unitName>
  </rotationCurveTask>
  !!]
  double precision function hotHaloMassDistributionRotationCurve(node,radius,componentType,massType)
    !!{
    Computes the rotation curve at a given radius for the hot halo density profile.
    !!}
    use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull, enumerationComponentTypeType, enumerationMassTypeType
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationComponentTypeType), intent(in   ) :: componentType
    type            (enumerationMassTypeType     ), intent(in   ) :: massType
    double precision                              , intent(in   ) :: radius
    double precision                                              :: componentMass

    ! Compute rotation curve if radius is non-zero.
    hotHaloMassDistributionRotationCurve=0.0d0
    if (radius > 0.0d0) then
       componentMass=hotHaloMassDistributionEnclosedMass(node,radius,componentType,massType,weightByMass,weightIndexNull)
       if (componentMass > 0.0d0)                                                         &
            & hotHaloMassDistributionRotationCurve=+sqrt(                                 &
            &                                            +gravitationalConstantGalacticus &
            &                                            *componentMass                   &
            &                                            /radius                          &
            &                                           )
    end if
    return
  end function hotHaloMassDistributionRotationCurve

  !![
  <rotationCurveGradientTask>
   <unitName>hotHaloMassDistributionRotationCurveGradient</unitName>
  </rotationCurveGradientTask>
  !!]
  double precision function hotHaloMassDistributionRotationCurveGradient(node,radius,componentType,massType)
    !!{
    Computes the rotation curve gradient at a given radius for the hot halo density profile.
    !!}
    use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull, enumerationComponentTypeType, enumerationMassTypeType
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationComponentTypeType), intent(in   ) :: componentType
    type            (enumerationMassTypeType     ), intent(in   ) :: massType
    double precision                              , intent(in   ) :: radius
    double precision                                              :: componentDensity, componentMass

    ! Set to zero by default.
    hotHaloMassDistributionRotationCurveGradient=0.0d0
    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=hotHaloMassDistributionEnclosedMass(node,radius,componentType,massType,weightByMass,weightIndexNull)
       if (componentMass > 0.0d0) then
          componentDensity                            = hotHaloMassDistribution_%density(node,radius)
          hotHaloMassDistributionRotationCurveGradient=+gravitationalConstantGalacticus    &
               &                                       *(                                  &
               &                                         -componentMass                    &
               &                                         /radius                       **2 &
               &                                         +4.0d0                            &
               &                                         *Pi                               &
               &                                         *radius                           &
               &                                         *componentDensity                 &
               &                                        )
       end if
    end if
    return
  end function hotHaloMassDistributionRotationCurveGradient

  !![
  <densityTask>
   <unitName>hotHaloMassDistributionDensity</unitName>
  </densityTask>
  !!]
  double precision function hotHaloMassDistributionDensity(node,positionSpherical,componentType,massType,weightBy,weightIndex)
    !!{
    Computes the density at a given position for a hot halo profile.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll       , componentTypeHotHalo, massTypeAll                 , massTypeBaryonic       , &
         &                                    massTypeGaseous        , weightByMass        , enumerationComponentTypeType, enumerationMassTypeType, &
         &                                    enumerationWeightByType
    implicit none
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationComponentTypeType), intent(in   ) :: componentType
    type            (enumerationMassTypeType     ), intent(in   ) :: massType
    type            (enumerationWeightByType     ), intent(in   ) :: weightBy
    integer                                       , intent(in   ) :: weightIndex
    double precision                              , intent(in   ) :: positionSpherical(3)
    !$GLC attributes unused :: weightIndex

    hotHaloMassDistributionDensity=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return

    hotHaloMassDistributionDensity=max(hotHaloMassDistribution_%density(node,positionSpherical(1)),0.0d0)
    return
  end function hotHaloMassDistributionDensity

  !![
  <densitySphericalAverageTask>
   <unitName>hotHaloMassDistributionDensitySphericalAverage</unitName>
  </densitySphericalAverageTask>
  !!]
  double precision function hotHaloMassDistributionDensitySphericalAverage(node,radius,componentType,massType,weightBy,weightIndex)
    !!{
    Computes the sphreically-averaged density at a given radius for a hot halo profile.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll       , componentTypeHotHalo, massTypeAll                 , massTypeBaryonic       , &
         &                                    massTypeGaseous        , weightByMass        , enumerationComponentTypeType, enumerationMassTypeType, &
         &                                    enumerationWeightByType
    implicit none
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationComponentTypeType), intent(in   ) :: componentType
    type            (enumerationMassTypeType     ), intent(in   ) :: massType
    type            (enumerationWeightByType     ), intent(in   ) :: weightBy
    integer                                       , intent(in   ) :: weightIndex
    double precision                              , intent(in   ) :: radius
    !$GLC attributes unused :: weightIndex

    hotHaloMassDistributionDensitySphericalAverage=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return

    hotHaloMassDistributionDensitySphericalAverage=max(hotHaloMassDistribution_%density(node,radius),0.0d0)
    return
  end function hotHaloMassDistributionDensitySphericalAverage

  !![
  <stateStoreTask>
   <unitName>hotHaloMassDistributionDefaultStateStore</unitName>
  </stateStoreTask>
  !!]
  subroutine hotHaloMassDistributionDefaultStateStore(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: hot halo mass distribution',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="hotHaloMassDistribution_ hotHaloTemperatureProfile_"/>
    !!]
    return
  end subroutine hotHaloMassDistributionDefaultStateStore

  !![
  <stateRetrieveTask>
   <unitName>hotHaloMassDistributionDefaultStateRestore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine hotHaloMassDistributionDefaultStateRestore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: hot halo mass distribution',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="hotHaloMassDistribution_ hotHaloTemperatureProfile_"/>
    !!]
    return
  end subroutine hotHaloMassDistributionDefaultStateRestore

end module Hot_Halo_Mass_Distributions
