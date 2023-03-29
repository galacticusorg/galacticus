!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a program to galactic structure functions.
!!}

program Test_Galactic_Structure
  !!{
  Tests galactic structure functions.
  !!}
  use :: Coordinates                     , only : assignment(=)                                                 , coordinateSpherical             , coordinateCartesian
  use :: Display                         , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Functions_Global_Utilities      , only : Functions_Global_Set
  use :: Galacticus_Nodes                , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize     , treeNode                           , nodeComponentDisk           , &
       &                                          nodeComponentBasic                                            , nodeComponentBlackHole
  use :: Input_Parameters                , only : inputParameters
  use :: Node_Components                 , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Galactic_Structure_Options      , only : componentTypeSpheroid                                         , componentTypeDisk                , massTypeAll                        , componentTypeDisk           , &
       &                                          coordinateSystemCartesian                                     , massTypeStellar                  , massTypeGaseous                    , enumerationMassTypeType     , &
       &                                          componentTypeBlackHole                                        , massTypeBlackHole
  use :: Unit_Tests                      , only : Assert                                                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  use :: Mass_Distributions              , only : massDistributionClass
  use :: Cosmology_Parameters            , only : cosmologyParametersSimple
  use :: Cosmology_Functions             , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMONFW
  use :: Dark_Matter_Profiles            , only : darkMatterProfileAdiabaticGnedin2004
  use :: Dark_Matter_Profiles_Generic    , only : nonAnalyticSolversNumerical
  use :: Bessel_Functions                , only : Bessel_Function_I0                                            , Bessel_Function_I1               , Bessel_Function_K0                 , Bessel_Function_K1
  use :: Galactic_Structure              , only : galacticStructureStandard
  use :: Tensors                         , only : assignment(=)
  use :: Virial_Density_Contrast         , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Numerical_Constants_Math        , only : Pi
  use :: Numerical_Constants_Prefixes    , only : kilo
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  use :: Numerical_Constants_Astronomical, only : gigaYear                                                      , gravitationalConstantGalacticus  , megaParsec
  implicit none 
  type            (inputParameters                                               )                 :: parameters
  type            (treeNode                                                      ), pointer        :: node_
  class           (nodeComponentBasic                                            ), pointer        :: basic_
  class           (nodeComponentDisk                                             ), pointer        :: disk_
  class           (nodeComponentBlackHole                                        ), pointer        :: blackHole_
  class           (massDistributionClass                                         ), pointer        :: massDistribution_
  type            (darkMatterProfileDMONFW                                       ), pointer        :: darkMatterProfileDMO_
  type            (darkMatterProfileAdiabaticGnedin2004                          ), pointer        :: darkMatterProfile_
  type            (galacticStructureStandard                                     ), pointer        :: galacticStructure_
  type            (cosmologyParametersSimple                                     ), pointer        :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                ), pointer        :: cosmologyFunctions_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer        :: darkMatterHaloScale_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer        :: virialDensityContrast_
  type            (coordinateCartesian                                           )                 :: positionCartesian                    , positionCartesianOuter
  type            (enumerationMassTypeType                                       )                 :: massType
  double precision                                                                , dimension(3  ) :: coordinatesCartesian                 , coordinatesCartesianOuter              , &
       &                                                                                              accelerationDirect                   , accelerationIndirect                   , &
       &                                                                                              accelerationTarget
  double precision                                                                , dimension(3,3) :: tidalTensorComponentsDirect          , tidalTensorSphericalComponents         , &
       &                                                                                              tidalTensorComponentsIndirect
  double precision                                                                , parameter      :: massDiskStellar              =1.0d+10, massDiskGas                    =1.000d9, &
       &                                                                                              radiusDisk                   =3.5d-03, scaleHeightFractionalDisk      =0.137d0, &
       &                                                                                              massBlackHole                =1.0d+05
  double precision                                                                                 :: densityTarget                        , surfaceDensityTarget                   , &
       &                                                                                              massTarget                           , radius                                 , &
       &                                                                                              densityDirect                        , densityIndirect                        , &
       &                                                                                              surfaceDensityDirect                 , surfaceDensityIndirect                 , &
       &                                                                                              massEnclosedDirect                   , massEnclosedIndirect                   , &
       &                                                                                              massEnclosedTarget                   , densitySphericalAverageTarget          , &
       &                                                                                              densitySphericalAverageDirect        , densitySphericalAverageIndirect        , &
       &                                                                                              rotationCurveDirect                  , rotationCurveIndirect                  , &
       &                                                                                              rotationCurveTarget                  , radiusHalf                             , &
       &                                                                                              rotationCurveGradientDirect          , rotationCurveGradientIndirect          , &
       &                                                                                              rotationCurveGradientTarget          , potentialTarget                        , &
       &                                                                                              potentialDirect                      , potentialIndirect                      , &
       &                                                                                              radiusOuter
  integer                                                                                          :: i                                    , j
  character       (len=15                                                        )                 :: label                                , labelRadius                            , &
       &                                                                                              labelRadiusOuter
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize.
  parameters=inputParameters('testSuite/parameters/galacticStructure.xml')
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build required objects.
  allocate(cosmologyParameters_                 )
  allocate(cosmologyFunctions_                  )
  allocate(virialDensityContrast_               )
  allocate(galacticStructure_           )
  allocate(darkMatterHaloScale_                 )
  allocate(darkMatterProfileDMO_             )
  allocate(darkMatterProfile_)
  !![
  <referenceConstruct object="cosmologyParameters_"  >
   <constructor>
    cosmologyParametersSimple                                     (                                                                  &amp;
     &amp;                                                         OmegaMatter                         = 0.3000d0                  , &amp;
     &amp;                                                         OmegaBaryon                         = 0.0450d0                  , &amp;
     &amp;                                                         OmegaDarkEnergy                     = 0.7000d0                  , &amp;
     &amp;                                                         temperatureCMB                      = 2.7800d0                  , &amp;
     &amp;                                                         HubbleConstant                      =69.3000d0                    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"   >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                                                  &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_         &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_">
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                  &amp;
     &amp;                                                         tableStore                          =.true.                     , &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_          &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"  >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                                                  &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_       , &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_        , &amp;
     &amp;                                                         virialDensityContrast_              =virialDensityContrast_       &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="galacticStructure_"    >
   <constructor>
    galacticStructureStandard                                     (                                                                  &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_        , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_       , &amp;
     &amp;                                                         darkMatterProfile_                  =darkMatterProfile_           &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMO_" >
   <constructor>
    darkMatterProfileDMONFW                                       (                                                                  &amp;
     &amp;                                                         velocityDispersionUseSeriesExpansion=.false.                    , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_         &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfile_"    >
   <constructor>
     darkMatterProfileAdiabaticGnedin2004                         (                                                                  &amp;
     &amp;                                                         A                                   =0.85d0                     , &amp;
     &amp;                                                         omega                               =0.80d0                     , &amp;
     &amp;                                                         radiusFractionalPivot               =1.00d0                     , &amp;
     &amp;                                                         toleranceRelative                   =1.0d-2                     , &amp;
     &amp;                                                         nonAnalyticSolver                   =nonAnalyticSolversNumerical, &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_       , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_       , &amp;
     &amp;                                                         darkMatterProfileDMO_               =darkMatterProfileDMO_      , &amp;
     &amp;                                                         galacticStructure_                  =galacticStructure_           &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  ! Build a node.
  node_      =>  treeNode          (                 )
  basic_     =>  node_   %basic    (autoCreate=.true.)
  disk_      =>  node_   %disk     (autoCreate=.true.)
  blackHole_ =>  node_   %blackHole(autoCreate=.true.)
  ! Set properties of the basic component.
  call basic_    %massSet            (1.00d12        )
  call basic_    %timeSet            (1.38d01        )
  call basic_    %timeLastIsolatedSet(1.38d01        )
  ! Set properties of the disk component.
  call disk_     %          radiusSet(radiusDisk     )
  call disk_     %     massStellarSet(massDiskStellar)
  call disk_     %         massGasSet(massDiskGas    )
  ! Set properties of the black hole component.
  call blackHole_%            massSet(massBlackHole  )
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Galactic structure functions")
  ! Get the mass distribution.
  massDistribution_ => node_%massDistribution()
  ! Test disk mass distribution.
  call Unit_Tests_Begin_Group("Disk")
  do j=1,3
     select case (j)
     case (1)
        radius= 5.0d-3
     case (2)
        radius= 8.0d-3
     case (3)
        radius=15.0d-3
     end select
     radiusHalf=+0.5d0      &
          &     *radius     &
          &     /radiusDisk
     write (labelRadius,'(f4.1)') kilo*radius
     do i=1,3
        select case (i)
        case (1)
           massTarget=+massDiskStellar
           massType  = massTypeStellar
           label     ='★    '
        case (2)
           massTarget=                +massDiskGas
           massType  = massTypeGaseous
           label     ='☁    '
        case (3)
           massTarget=+massDiskStellar+massDiskGas
           massType  = massTypeAll
           label     ='ₜₒₜₐₗ'
        end select
        coordinatesCartesian=[radius,0.0d0,0.0d0]
        positionCartesian   =coordinatesCartesian
        ! Surface density.
        surfaceDensityTarget  =+massTarget                      &
             &                 /2.0d0                           &
             &                 /Pi                              &
             &                 /radiusDisk**2                   &
             &                 *exp (                           &
             &                       -coordinatesCartesian(1)   &
             &                       /radiusDisk                &
             &                      )
        surfaceDensityDirect  =+massDistribution_ %surfaceDensity(      positionCartesian   ,componentType=componentTypeDisk,massType=massType                                           )
        surfaceDensityIndirect=+galacticStructure_%surfaceDensity(node_,coordinatesCartesian,componentType=componentTypeDisk,massType=massType,coordinateSystem=coordinateSystemCartesian)
        call Assert("Σ"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",surfaceDensityDirect  ,surfaceDensityTarget,relTol=1.0d-6)
        call Assert("Σ"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",surfaceDensityIndirect,surfaceDensityTarget,relTol=1.0d-6)
        ! Density.
        densityTarget  =+surfaceDensityTarget            &
             &          /cosh(                           &
             &                +coordinatesCartesian(3)   &
             &                /radiusDisk                &
             &                /scaleHeightFractionalDisk &
             &               )**2                        &
             &          /2.0d0                           &
             &          /radiusDisk                      &
             &          /scaleHeightFractionalDisk
        densityDirect  =+massDistribution_ %density(      positionCartesian   ,componentType=componentTypeDisk,massType=massType                                           )
        densityIndirect=+galacticStructure_%density(node_,coordinatesCartesian,componentType=componentTypeDisk,massType=massType,coordinateSystem=coordinateSystemCartesian)
        call Assert("ρ"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",densityDirect  ,densityTarget,relTol=1.0d-6)
        call Assert("ρ"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",densityIndirect,densityTarget,relTol=1.0d-6)
        ! Enclosed mass.
        massEnclosedTarget  =+massTarget               &
             &               *(                        &
             &                 +  1.0d0                &
             &                 -(                      &
             &                   +1.0d0                &
             &                   +   radius/radiusDisk &
             &                  )                      &
             &                 *exp(                   &
             &                      -radius/radiusDisk &
             &                     )                   &
             &                )
        massEnclosedDirect  =+massDistribution_ %massEnclosedBySphere(      radius,componentType=componentTypeDisk,massType=massType)
        massEnclosedIndirect=+galacticStructure_%massEnclosed        (node_,radius,componentType=componentTypeDisk,massType=massType)
        call Assert("M"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",massEnclosedDirect  ,massEnclosedTarget,relTol=1.0d-6)
        call Assert("M"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",massEnclosedIndirect,massEnclosedTarget,relTol=1.0d-6)
        ! Density spherical average.
        densitySphericalAverageTarget=+massTarget             &
             &                        *     radius/radiusDisk &
             &                        *exp(                   &
             &                             -radius/radiusDisk &
             &                            )                   &
             &                        /4.0d0                  &
             &                        /Pi                     &
             &                        /radiusDisk**3
        densitySphericalAverageDirect  =+massDistribution_ %densitySphericalAverage(      radius,componentType=componentTypeDisk,massType=massType)
        densitySphericalAverageIndirect=+galacticStructure_%densitySphericalAverage(node_,radius,componentType=componentTypeDisk,massType=massType)
        call Assert("⟨ρ̂⟩"//char(8)//char(8)//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",massEnclosedDirect  ,massEnclosedTarget,relTol=1.0d-6)
        call Assert("⟨ρ̂⟩"//char(8)//char(8)//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",massEnclosedIndirect,massEnclosedTarget,relTol=1.0d-6)
        ! Rotation curve.
        rotationCurveTarget  =+sqrt(                                                                           &
             &                      +2.0d0                                                                     &
             &                      *gravitationalConstantGalacticus                                           &
             &                      *massTarget                                                                &
             &                      /radiusDisk                                                                &
             &                      *radiusHalf**2                                                             &
             &                      *(                                                                         &
             &                                  +Bessel_Function_I0(radiusHalf)*Bessel_Function_K0(radiusHalf) &
             &                                  -Bessel_Function_I1(radiusHalf)*Bessel_Function_K1(radiusHalf) &
             &                       )                                                                         &
             &                     )
        rotationCurveDirect  =+massDistribution_ %rotationCurve   (      radius,componentType=componentTypeDisk,massType=massType)
        rotationCurveIndirect=+galacticStructure_%velocityRotation(node_,radius,componentType=componentTypeDisk,massType=massType)
        call Assert("V"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",rotationCurveDirect  ,rotationCurveTarget,relTol=1.0d-2)
        call Assert("V"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",rotationCurveIndirect,rotationCurveTarget,relTol=1.0d-2)
        ! Rotation curve gradient.
        rotationCurveGradientTarget  =+gravitationalConstantGalacticus                                                                                &
             &                        *massTarget                                                                                                     &
             &                        /radiusDisk**2                                                                                                  &
             &                        *(                                                                                                              &
             &                          +2.0d0                                                                                                        &
             &                          *radiusHalf                                                                                                   &
             &                          *(                                                                                                            &
             &                            +  Bessel_Function_I0(radiusHalf)                                           *Bessel_Function_K0(radiusHalf) &
             &                            -  Bessel_Function_I1(radiusHalf)                                           *Bessel_Function_K1(radiusHalf) &
             &                           )                                                                                                            &
             &                          +radiusHalf**2                                                                                                &
             &                          *(                                                                                                            &
             &                            +  Bessel_Function_I1(radiusHalf)                                           *Bessel_Function_K0(radiusHalf) &
             &                            -  Bessel_Function_K1(radiusHalf)                                           *Bessel_Function_I0(radiusHalf) &
             &                            -(+Bessel_Function_I0(radiusHalf)-Bessel_Function_I1(radiusHalf)/radiusHalf)*Bessel_Function_K1(radiusHalf) &
             &                            -(-Bessel_Function_K0(radiusHalf)-Bessel_Function_K1(radiusHalf)/radiusHalf)*Bessel_Function_I1(radiusHalf) &
             &                           )                                                                                                            &
             &                         )
        rotationCurveGradientDirect  =+massDistribution_ %rotationCurveGradient   (      radius,componentType=componentTypeDisk,massType=massType)
        rotationCurveGradientIndirect=+galacticStructure_%velocityRotationGradient(node_,radius,componentType=componentTypeDisk,massType=massType) &
             &                        *2.0d0                                                                                                       &
             &                        *rotationCurveIndirect
        call Assert("dV"//trim(label)//"/dr(r="//trim(labelRadius)//"kpc) [  direct]",rotationCurveGradientDirect  ,rotationCurveGradientTarget,relTol=1.0d-2)
        call Assert("dV"//trim(label)//"/dr(r="//trim(labelRadius)//"kpc) [indirect]",rotationCurveGradientIndirect,rotationCurveGradientTarget,relTol=1.0d-2)
     end do
     ! Evaluate non-analytic properties at large radii and compare the the result for a spherical mass distribution.
     radius              =50.0d0*radiusDisk
     coordinatesCartesian=[radius,0.0d0,0.0d0]
     positionCartesian   =coordinatesCartesian
     write (labelRadius,'(f5.1)') kilo*radius
     ! Potential.
     radiusOuter              =75.0d0*radiusDisk
     coordinatesCartesianOuter=[radiusOuter,0.0d0,0.0d0]
     positionCartesianOuter   =coordinatesCartesianOuter
     write (labelRadiusOuter,'(f5.1)') kilo*radiusOuter
     potentialTarget  =-gravitationalConstantGalacticus*massTarget*(1.0d0/radius-1.0d0/radiusOuter)
     potentialDirect  =+massDistribution_ %potential(      positionCartesian     ,componentType=componentTypeDisk,massType=massType) &
          &            -massDistribution_ %potential(      positionCartesianOuter,componentType=componentTypeDisk,massType=massType)
     potentialIndirect=+galacticStructure_%potential(node_,radius                ,componentType=componentTypeDisk,massType=massType) &
          &            -galacticStructure_%potential(node_,radiusOuter           ,componentType=componentTypeDisk,massType=massType)
     call Assert("Φ"//trim(label)//"(r="//trim(labelRadius)//"kpc)-Φ"//trim(label)//"(r="//trim(labelRadiusOuter)//"kpc) [  direct]",potentialDirect  ,potentialTarget,relTol=2.0d-2)
     call Assert("Φ"//trim(label)//"(r="//trim(labelRadius)//"kpc)-Φ"//trim(label)//"(r="//trim(labelRadiusOuter)//"kpc) [indirect]",potentialIndirect,potentialTarget,relTol=2.0d-2)
     ! Acceleration.
     accelerationTarget  =-gravitationalConstantGalacticus &
          &               *massTarget                      &
          &               /radius**3                       &
          &               *coordinatesCartesian            &
          &               *kilo                            &
          &               *gigaYear                        &
          &               /megaParsec     
     accelerationDirect  =+massDistribution_ %acceleration(      positionCartesian   ,componentType=componentTypeDisk,massType=massType)
     accelerationIndirect=+galacticStructure_%acceleration(node_,coordinatesCartesian,componentType=componentTypeDisk,massType=massType)
     call Assert("a⃗"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",accelerationDirect  ,accelerationTarget,relTol=2.0d-2)
     call Assert("a⃗"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",accelerationIndirect,accelerationTarget,relTol=2.0d-2)
     ! Tidal tensor.
     tidalTensorSphericalComponents=reshape([+2.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0]*gravitationalConstantGalacticus*massTarget/radius**3,[3,3])
     tidalTensorComponentsDirect   =massDistribution_ %tidalTensor(      positionCartesian   ,componentType=componentTypeDisk,massType=massType)
     tidalTensorComponentsIndirect =galacticStructure_%tidalTensor(node_,coordinatesCartesian,componentType=componentTypeDisk,massType=massType)
     call Assert("g⃡"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",tidalTensorComponentsDirect  ,tidalTensorSphericalComponents,absTol=1.0d-6,relTol=3.0d-3)
     call Assert("g⃡"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",tidalTensorComponentsIndirect,tidalTensorSphericalComponents,absTol=1.0d-6,relTol=3.0d-3)
  end do
  call Unit_Tests_End_Group()
  ! Test black hole distribution.
  call Unit_Tests_Begin_Group("Black hole")
  do j=1,3
     select case (j)
     case (1)
        radius= 5.0d-3
     case (2)
        radius= 8.0d-3
     case (3)
        radius=15.0d-3
     end select
     radiusHalf=+0.5d0      &
          &     *radius     &
          &     /radiusDisk
     write (labelRadius,'(f4.1)') kilo*radius
     massTarget=+massBlackHole
     massType  = massTypeBlackHole
     label     ='∙    '
     coordinatesCartesian=[radius,0.0d0,0.0d0]
     positionCartesian   =coordinatesCartesian
     ! Enclosed mass.
     massEnclosedTarget  =+massTarget
     massEnclosedDirect  =+massDistribution_ %massEnclosedBySphere(      radius,componentType=componentTypeBlackHole,massType=massType)
     massEnclosedIndirect=+galacticStructure_%massEnclosed        (node_,radius,componentType=componentTypeBlackHole,massType=massType)
     call Assert("M"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",massEnclosedDirect  ,massEnclosedTarget,relTol=1.0d-6)
     call Assert("M"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",massEnclosedIndirect,massEnclosedTarget,relTol=1.0d-6)
     ! Rotation curve.
     rotationCurveTarget  =+sqrt(                                 &
          &                      +gravitationalConstantGalacticus &
          &                      *massTarget                      &
          &                      /radius                          &
          &                     )
     rotationCurveDirect  =+massDistribution_ %rotationCurve   (      radius,componentType=componentTypeBlackHole,massType=massType)
     rotationCurveIndirect=+galacticStructure_%velocityRotation(node_,radius,componentType=componentTypeBlackHole,massType=massType)
     call Assert("V"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",rotationCurveDirect  ,rotationCurveTarget,relTol=1.0d-2)
     call Assert("V"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",rotationCurveIndirect,rotationCurveTarget,relTol=1.0d-2)
     ! Rotation curve gradient.
     rotationCurveGradientTarget  =-gravitationalConstantGalacticus &
          &                        *massTarget                      &
          &                        /radius**2
     rotationCurveGradientDirect  =+massDistribution_ %rotationCurveGradient   (      radius,componentType=componentTypeBlackHole,massType=massType)
     rotationCurveGradientIndirect=+galacticStructure_%velocityRotationGradient(node_,radius,componentType=componentTypeBlackHole,massType=massType) &
          &                        *2.0d0                                                                                                            &
          &                        *rotationCurveIndirect
     call Assert("dV"//trim(label)//"/dr(r="//trim(labelRadius)//"kpc) [  direct]",rotationCurveGradientDirect  ,rotationCurveGradientTarget,relTol=1.0d-2)
     call Assert("dV"//trim(label)//"/dr(r="//trim(labelRadius)//"kpc) [indirect]",rotationCurveGradientIndirect,rotationCurveGradientTarget,relTol=1.0d-2)
     ! Acceleration.
     accelerationTarget  =-gravitationalConstantGalacticus &
          &               *massTarget                      &
          &               /radius**3                       &
          &               *coordinatesCartesian            &
          &               *kilo                            &
          &               *gigaYear                        &
          &               /megaParsec     
     accelerationDirect  =+massDistribution_ %acceleration(      positionCartesian   ,componentType=componentTypeBlackHole,massType=massType)
     accelerationIndirect=+galacticStructure_%acceleration(node_,coordinatesCartesian,componentType=componentTypeBlackHole,massType=massType)
     call Assert("a⃗"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",accelerationDirect  ,accelerationTarget,relTol=2.0d-2)
     call Assert("a⃗"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",accelerationIndirect,accelerationTarget,relTol=2.0d-2)
     ! Tidal tensor.
     tidalTensorSphericalComponents=reshape([+2.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0]*gravitationalConstantGalacticus*massTarget/radius**3,[3,3])
     tidalTensorComponentsDirect   =massDistribution_ %tidalTensor(      positionCartesian   ,componentType=componentTypeBlackHole,massType=massType)
     tidalTensorComponentsIndirect =galacticStructure_%tidalTensor(node_,coordinatesCartesian,componentType=componentTypeBlackHole,massType=massType)
     call Assert("g⃡"//trim(label)//"(r="//trim(labelRadius)//"kpc) [  direct]",tidalTensorComponentsDirect  ,tidalTensorSphericalComponents,absTol=1.0d-6,relTol=3.0d-3)
     call Assert("g⃡"//trim(label)//"(r="//trim(labelRadius)//"kpc) [indirect]",tidalTensorComponentsIndirect,tidalTensorSphericalComponents,absTol=1.0d-6,relTol=3.0d-3)
  end do
  call Unit_Tests_End_Group()
  ! Clean up objects.
  !![
  <objectDestructor name="massDistribution_"     />
  <objectDestructor name="cosmologyParameters_"  />
  <objectDestructor name="cosmologyFunctions_"   />
  <objectDestructor name="virialDensityContrast_"/>
  <objectDestructor name="darkMatterHaloScale_"  />
  <objectDestructor name="darkMatterProfileDMO_" />
  <objectDestructor name="galacticStructure_"    />
  !!]
  call node_%destroy()
  deallocate(node_)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Galactic_Structure
