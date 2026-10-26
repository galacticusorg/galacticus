!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program to test mass distributions.
!!}

program Test_Mass_Distributions
  !!{RST
  Tests mass distributions.
  !!}
  use :: Coordinates                   , only : assignment(=)                        , coordinateSpherical                   , coordinateCartesian
  use :: Display                       , only : displayVerbositySet                  , verbosityLevelStandard                , displayMessage
  use :: Error_Functions               , only : Error_Function
  use :: Error                         , only : Error_Report
  use :: Events_Hooks                  , only : eventsHooksInitialize
  use :: Galactic_Component_Truncations, only : galacticComponentTruncationClass     , galacticComponentTruncationExponential
  use :: Galactic_Structure_Options    , only : componentTypeSpheroid                , componentTypeDisk                     , radiusLarge
  use :: IO_HDF5                       , only : ioHDF5AccessInitialize
  use :: Linear_Algebra                , only : assignment(=)                        , vector
  use :: Mass_Distributions            , only : massDistributionBetaProfile          , massDistributionClass                 , massDistributionExponentialDisk        , massDistributionGaussianEllipsoid           , &
       &                                        massDistributionHernquist            , massDistributionSersic                , massDistributionSpherical              , massDistributionComposite                   , &
       &                                        massDistributionList                 , massDistributionSymmetryCylindrical   , enumerationMassDistributionSymmetryType, massDistributionSphericalScaler             , &
       &                                        massDistributionCylindricalScaler    , massDistributionCylindrical           , massDistributionPatejLoeb2015          , massDistributionNFW                         , &
       &                                        massDistributionSIDMParametricProfile, massDistributionBlackHole             , massDistributionSoliton                , massDistributionSolitonNFW                  , &
       &                                        massDistributionIsothermal           , kinematicsDistributionSphericalScaler , kinematicsDistributionNFW              , kinematicsDistributionCollisionlessTabulated, &
       &                                        kinematicsDistributionClass          , kinematicsDistributionLocal           , massDistributionSphericalHeated        , massDistributionHeatingTidal                , &
       &                                        massDistributionHeatingClass         , massDistributionSphericalHeatedTracer , nonAnalyticSolversNumerical            , massDistributionSphericalAccelerator
  use :: Numerical_Constants_Math      , only : Pi                                   , e
  use :: Tensors                       , only : assignment(=)
  use :: Unit_Tests                    , only : Assert                               , Unit_Tests_Begin_Group                , Unit_Tests_End_Group                   , Unit_Tests_Finish                           , &
       &                                        compareEquals                        , compareNotEqual                       , compareLessThan
  implicit none
  class           (massDistributionClass                  )                                    , allocatable :: massDistribution_                                                                                                       , massDistributionRotated         , &
       &                                                                                                        massDistributionDisk                                                                                                    , massDistributionSpheroid        , &
       &                                                                                                        massDistributionDMO                                                                                                     , massDistributionWrappedSpherical, &
       &                                                                                                        massDistributionWrappedCylindrical
  ! Dedicated objects for the spherical scaler memoization tests, so that those tests do not disturb (or
  ! depend on the allocation state of) the shared objects above.
  class           (massDistributionClass                  )                                    , pointer     :: massDistributionMomentInner
  class           (massDistributionClass                  )                                    , allocatable :: massDistributionScalerTest                                                                                              , massDistributionScalerInner
  class           (massDistributionClass                  )                                    , pointer     :: massDistributionDisk_                                                                                                   , massDistributionSpheroid_       , &
       &                                                                                                        massDistributionAll_                                                                                                    , massDistributionDiskAgain_
  class           (kinematicsDistributionClass            )                                    , pointer     :: kinematicsDistribution_
  ! Objects for the spherical scaler re-initialization (pool re-use) test.
  class           (massDistributionClass                  )                                    , pointer     :: massDistributionScalerInnerA                                                                                            , massDistributionScalerInnerB
  class           (massDistributionClass                  )                                    , allocatable :: massDistributionScalerFresh
  class           (kinematicsDistributionClass            )                                    , pointer     :: kinematicsDistributionInnerA                                                                                            , kinematicsDistributionInnerB    , &
       &                                                                                                        kinematicsDistributionScaler_
  double precision                                         , dimension(5)                                    :: scalerReinitializeFresh                                                                                                 , scalerReinitializeReused
  ! Objects for the galactic component truncation tests.
  class           (galacticComponentTruncationClass       )                                    , allocatable :: galacticComponentTruncation_
  class           (massDistributionSpherical              )                                    , pointer     :: massDistributionUntruncated_                                                                                            , massDistributionTruncated_
  ! Reference values for exponential truncation of a dimensionless Hernquist profile (with r_decay/r_t=0.1), computed
  ! independently by root-finding the truncation radius and direct numerical quadrature (scipy) of the truncated profile.
  integer                                                  , parameter                                       :: truncationTableCount                 =4
  double precision                                         , dimension(truncationTableCount  )              :: truncationTableFraction              =[0.90d0           ,0.50d0           ,0.10d0           ,0.01d0           ], &
       &                                                                                                        truncationTableRadiusHalfMassTarget  =[2.0378552605d+00,1.0000000000d+00,2.8800715553d-01,7.6092056823d-02]
  type            (massDistributionList                   )                                    , pointer     :: massDistributions
  integer                                                  , parameter                                       :: sidmParametricTableCount              =7
  double precision                                         , dimension(sidmParametricTableCount)             :: sidmParametricTableRadius             =[+0.125000000000000d+0,+0.250000000000000d+0,+0.500000000000000d+0,+1.000000000000000d+0,+2.000000000000000d+0,+4.000000000000000d+0,+8.000000000000000d+0], &
       &                                                                                                        sidmParametricTableDensityTarget      =[+1.578707460302960d+0,+1.260746394893958d+0,+7.474634802255240d-1,+2.462395302527262d-1,+5.550143415127595d-2,+9.999389741552710d-3,+1.543203989720854d-3], &
       &                                                                                                        sidmParametricTableGradientTarget     =[-2.261132728058798d-1,-4.588235294117647d-1,-1.166666666666667d+0,-1.941176470588235d+0,-2.329442282749676d+0,-2.599755918965097d+0,-2.777762519221542d+0], &
       &                                                                                                        sidmParametricTableMassTarget         =[+1.369278576289049d-2,+9.273816539469000d-2,+5.224473772681247d-1,+1.957581650325344d+0,+4.944657616168121d+0,+9.686918108280710d+0,+1.595612873885020d+1], &
       &                                                                                                        sidmParametricTableRadialMoment0Target=[+2.221821944239244d-1,+3.989488972124282d-1,+6.471785700600648d-1,+8.687852495825560d-1,+9.890703047718140d-1,+1.038038927921240d+0,+1.054510130769266d+0], &
       &                                                                                                        sidmParametricTableRadialMoment1Target=[+1.333970313751376d-2,+4.607079153523620d-2,+1.364760414161601d-1,+2.924871265742027d-1,+4.582720079481481d-1,+5.915584441874174d-1,+6.804453157808887d-1], &
       &                                                                                                        sidmParametricTableRadialMoment2Target=[+1.089637269431175d-3,+7.379868717919331d-3,+4.157504129880918d-2,+1.557793980776343d-1,+3.934833507550720d-1,+7.708604501296330d-1,+1.269748380699330d+0], &
       &                                                                                                        sidmParametricTableRadialMoment3Target=[+1.009434956307133d-4,+1.351333811789922d-3,+1.474631367295585d-2,+1.015426667994900d-1,+4.558765935959807d-1,+1.567266971756937d+0,+4.480533036621432d+0], &
       &                                                                                                        sidmParametricTableFourierTarget      =[+1.223827808259200d-3,+1.809381194427616d-2,+8.926657506109750d-2,+2.404779846784714d-1,+5.298489986770916d-1,+8.402142769464090d-1,+9.566992274102220d-1]
  double precision                                         , dimension(sidmParametricTableCount)             :: sidmParametricTableDensity                                                                                                                                                                        , &
       &                                                                                                        sidmParametricTableGradient                                                                                                                                                                       , &
       &                                                                                                        sidmParametricTableMass                                                                                                                                                                           , &
       &                                                                                                        sidmParametricTableRadialMoment0                                                                                                                                                                  , &
       &                                                                                                        sidmParametricTableRadialMoment1                                                                                                                                                                  , &
       &                                                                                                        sidmParametricTableRadialMoment2                                                                                                                                                                  , &
       &                                                                                                        sidmParametricTableRadialMoment3                                                                                                                                                                  , &
       &                                                                                                        sidmParametricTableFourier                                                                                                                                                                        , &
       &                                                                                                        sidmParametricTableRadiusMass                                                                                                                                                                     , &
       &                                                                                                        sidmParametricTableRadiusDensity
  integer                                                  , parameter                                       :: sersicTableCount                     =8
  double precision                                         , dimension(sersicTableCount       )              :: sersicTableRadius                    =[1.0000d-06,1.0000d-5,1.0000d-4,1.0000d-3,1.0000d-2,1.0000d-1,1.0000d+0,1.0000d+1]
  ! Mass targets for Sersic profile from Mazure & Capelato (2001).
  double precision                                         , dimension(sersicTableCount       )              :: sersicTableMassTarget                =[1.4730d-11,2.1130d-9,2.5959d-7,2.4545d-5,1.4961d-3,4.4102d-2,4.1536d-1,9.4308d-1]
  ! Density targets for Sersic profile from Mazure & Capelato (2001).
  double precision                                         , dimension(sersicTableCount       )              :: sersicTableDensityTarget             =[2.5553d+06,3.5797d+5,4.2189d+4,3.7044d+3,1.9679d+2,4.4047d+0,2.1943d-2,7.8166d-6]
  ! Potential targets for Sersic profile from Young (1976).
  double precision                                         , dimension(sersicTableCount       )              :: sersicTablePotentialTarget           =[1.0000d+00,9.9993d-1,9.9908d-1,9.9027d-1,9.2671d-1,6.7129d-1,2.4945d-1,3.7383d-2]
  double precision                                         , dimension(sersicTableCount       )              :: sersicTableDensity                                                                                                      , sersicTableMass                                  , &
       &                                                                                                        sersicTablePotential
  type            (coordinateSpherical                    )                                                  :: position                                                                                                                , positionZero                                     , &
       &                                                                                                        positionReference
  type            (coordinateCartesian                    )                                                  :: positionCartesian
  type            (enumerationMassDistributionSymmetryType)                                                  :: symmetry_
  integer                                                                                                    :: i                                                                                                                       , referenceCountDiscard                            , &
       &                                                                                                        indexDisk
  double precision                                                                                           :: radiusInProjection                                                                                                      , radius                                           , &
       &                                                                                                        rotationCurveGradientAnalytic                                                                                           , rotationCurveGradientNumerical                   , &
       &                                                                                                        massFraction                                                                                                            , time                                             , &
       &                                                                                                        radiusCylindrical                                                                                                       , height                                           , &
       &                                                                                                        scaleHeight                                                                                                             , potentialGradientRadial                          , &
       &                                                                                                        potentialGradientVertical
  ! Reference values for the exponential disk profile.
  double precision                                         , dimension(3,4)                                  :: accelerationReferences
  double precision                                         , dimension(  6)                                  :: potentialReferences
  double precision                                         , dimension(  4)                                  :: rotationCurveReferences
  ! Parameters and reference values for the black hole re-initialization (pool re-use) test.
  double precision                                         , parameter                                       :: massBlackHoleA                       =1.0d+6                                                                            , massBlackHoleB                             =4.0d8, &
       &                                                                                                        radiusBlackHoleTest                  =1.0d-3
  double precision                                                                                           :: enclosedReferenceA                                                                                                      , enclosedReferenceB                               , &
       &                                                                                                        rotationReferenceA                                                                                                      , rotationReferenceB
  double precision                                         , parameter                                       :: epsilonFiniteDifference              =0.01d0
  double precision                                         , parameter                                       :: epsilonFiniteDifferencePotential     =0.10d0
  character       (len=4                                  )                                                  :: label
  character       (len=5                                  )                                                  :: labelHeight                                                                                                             , labelRadius
  double precision                                         , dimension(3,3)                                  :: tidalTensorComponents                                                                                                   , tidalTensorSphericalComponents
  double precision                                         , dimension(3,3)                                  :: tidalTensorScalerReference                                                                                              , tidalTensorScalerDisplaced
  double precision                                         , dimension(3  )                                  :: acceleration                                                                                                            , accelerationReference
  double precision                                         , dimension(4  )                                  :: massPatejLoeb                                                                                                           , densityPatejLoeb                                 , &
       &                                                                                                        densitySlopePatejLoeb                                                                                                   , densityMomentPatejLoeb                           , &
       &                                                                                                        potentialPatejLoeb
  double precision                                         , dimension(4  )                                  :: massIsothermal                                                                                                          , densityIsothermal                                , &
       &                                                                                                        densitySlopeIsothermal                                                                                                  , densityMomentIsothermal                          , &
       &                                                                                                        potentialIsothermal                                                                                                     , fourierTransformIsothermal                       , &
       &                                                                                                        radiusFreeFallIsothermal                                                                                                , radiusFreeFallGrowthRateIsothermal               , &
       &                                                                                                        massIsothermalNumerical                                                                                                 , potentialIsothermalDifferenceNumerical           , &
       &                                                                                                        fourierTransformIsothermalNumerical                                                                                     , radiusFreeFallIsothermalNumerical                , &
       &                                                                                                        radiusFreeFallGrowthRateIsothermalNumerical                                                                             , densitySlopeIsothermalNumerical                  , &
       &                                                                                                        radiusEnclosingMassIsothermal                                                                                           , radiusEnclosingMassIsothermalNumerical           , &
       &                                                                                                        radiusEnclosingDensityIsothermal                                                                                        , radiusEnclosingDensityIsothermalNumerical        , &
       &                                                                                                        radiiIsothermal                                                                                                         , radiusFromSpecificAngularMomentumIsothermal      , &
       &                                                                                                        radiusFromSpecificAngularMomentumIsothermalNumerical                                                                    , velocityCircularIsothermal
  type            (vector                                 ), dimension(:  )                    , allocatable :: axes
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Initialize HDF5 lock.
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Mass distributions")

  ! SIDM parametric profile.
  call Unit_Tests_Begin_Group("SIDM parametric profile")
  allocate(massDistributionSIDMParametricProfile :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionSIDMParametricProfile)
     allocate(kinematicsDistributionCollisionlessTabulated :: kinematicsDistribution_)
     select type (kinematicsDistribution_)
     type is (kinematicsDistributionCollisionlessTabulated)
        kinematicsDistribution_=kinematicsDistributionCollisionlessTabulated(toleranceRelativeVelocityDispersion=1.0d-3,toleranceRelativeVelocityDispersionMaximum=1.0d-3)
     end select
     massDistribution_=massDistributionSIDMParametricProfile(beta=4.0d0,densityNormalization=1.0d0,radiusScale=1.0d0,radiusCore=0.5d0)
     call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,sidmParametricTableCount
        position                           =[sidmParametricTableRadius(i),0.0d0,0.0d0]
        sidmParametricTableDensity      (i)=massDistribution_%density               (                   coordinates  =      position                                                                 )
        sidmParametricTableGradient     (i)=massDistribution_%densityGradientRadial (logarithmic=.true.,coordinates  =      position                                                                 )
        sidmParametricTableMass         (i)=massDistribution_%massEnclosedBySphere  (                   radius       =      sidmParametricTableRadius    (i)                                         )
        sidmParametricTableRadialMoment0(i)=massDistribution_%densityRadialMoment   (moment     =0.0d0 ,radiusMaximum=      sidmParametricTableRadius    (i)                                         )
        sidmParametricTableRadialMoment1(i)=massDistribution_%densityRadialMoment   (moment     =1.0d0 ,radiusMaximum=      sidmParametricTableRadius    (i)                                         )
        sidmParametricTableRadialMoment2(i)=massDistribution_%densityRadialMoment   (moment     =2.0d0 ,radiusMaximum=      sidmParametricTableRadius    (i)                                         )
        sidmParametricTableRadialMoment3(i)=massDistribution_%densityRadialMoment   (moment     =3.0d0 ,radiusMaximum=      sidmParametricTableRadius    (i)                                         )
        sidmParametricTableFourier      (i)=massDistribution_%fourierTransform      (radiusOuter=8.0d0 ,wavenumber   =1.0d0/sidmParametricTableRadius    (i)                                         )
        sidmParametricTableRadiusMass   (i)=massDistribution_%radiusEnclosingMass   (                   mass         =      sidmParametricTableMassTarget(i)                                         )
        sidmParametricTableRadiusDensity(i)=massDistribution_%radiusEnclosingDensity(                   density      =3.0d0*sidmParametricTableMassTarget(i)/4.0d0/Pi/sidmParametricTableRadius(i)**3)
     end do
     call Assert("Density"                 ,sidmParametricTableDensity      ,sidmParametricTableDensityTarget      ,relTol=1.0d-6)
     call Assert("Density gradient"        ,sidmParametricTableGradient     ,sidmParametricTableGradientTarget     ,relTol=1.0d-6)
     call Assert("Mass enclosed"           ,sidmParametricTableMass         ,sidmParametricTableMassTarget         ,relTol=1.0d-6)
     call Assert("Radial moment (m=0) "    ,sidmParametricTableRadialMoment0,sidmParametricTableRadialMoment0Target,relTol=1.0d-6)
     call Assert("Radial moment (m=1) "    ,sidmParametricTableRadialMoment1,sidmParametricTableRadialMoment1Target,relTol=1.0d-6)
     call Assert("Radial moment (m=2) "    ,sidmParametricTableRadialMoment2,sidmParametricTableRadialMoment2Target,relTol=1.0d-6)
     call Assert("Radial moment (m=3) "    ,sidmParametricTableRadialMoment3,sidmParametricTableRadialMoment3Target,relTol=1.0d-6)
     call Assert("Fourier transform"       ,sidmParametricTableFourier      ,sidmParametricTableFourierTarget      ,relTol=1.0d-6)
     call Assert("Radius enclosing mass"   ,sidmParametricTableRadiusMass   ,sidmParametricTableRadius             ,relTol=1.0d-6)
     call Assert("Radius enclosing density",sidmParametricTableRadiusDensity,sidmParametricTableRadius             ,relTol=1.0d-2)
     ! Also test the central density is respected when seeking the radius enclosing a density that exceeds the central density.
     call Assert("Radius enclosing density (central)",massDistribution_%radiusEnclosingDensity(1.0d2),0.0d0,absTol=1.0d-6)
  end select
  nullify   (kinematicsDistribution_)
  deallocate(massDistribution_      )
  call Unit_Tests_End_Group()

  ! Hernquist profile.
  call Unit_Tests_Begin_Group("Hernquist profile")
  allocate(massDistributionHernquist :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionHernquist)
     massDistribution_=massDistributionHernquist(dimensionless=.true.)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     position=[1.0d0,0.0d0,0.0d0]
     call Assert("Half mass radius          (        dimensionless)",massDistribution_%radiusHalfMass      (        ), 1.0d0+sqrt(2.0d0)        ,absTol=1.0d-6)
     call Assert("Mass within scale radius  (        dimensionless)",massDistribution_%massEnclosedBySphere(1.0d0   ), 0.25d0                   ,absTol=1.0d-6)
     call Assert("Potential at scale radius (        dimensionless)",massDistribution_%potential           (position),-0.50d0                   ,absTol=1.0d-6)
  end select
  deallocate(massDistribution_)
  allocate(massDistributionHernquist :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionHernquist)
     massDistribution_=massDistributionHernquist(dimensionless=.false.,mass=2.0d0,scaleLength=2.0d0)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     call Assert("Half mass radius          (        dimensionful )",massDistribution_%radiusHalfMass         (        ),2.0d0     *(+1.0d0+sqrt(2.0d0))                              ,absTol=1.0d-6)
     call Assert("Mass within scale radius  (        dimensionful )",massDistribution_%massEnclosedBySphere   (2.0d0   ),2.0d0     *(+0.25d0           )                              ,absTol=1.0d-6)
     call Assert("Mass within cylinder      (r=0.5rₛ; dimensionful )",massDistribution_%massEnclosedByCylinder(1.0d0   ),2.0d0/3.0d0*(-1.0d0+2.0d0*log(2.0d0+sqrt(3.0d0))/sqrt(3.0d0)),absTol=1.0d-6)
     call Assert("Mass within cylinder      (r=1.0rₛ; dimensionful )",massDistribution_%massEnclosedByCylinder(2.0d0   ),2.0d0/3.0d0                                                  ,absTol=1.0d-6)
     call Assert("Mass within cylinder      (r=2.0rₛ; dimensionful )",massDistribution_%massEnclosedByCylinder(4.0d0   ),8.0d0/3.0d0*(+1.0d0-Pi/3.0d0                    /sqrt(3.0d0)),absTol=1.0d-6)
  end select
  deallocate(massDistribution_)
  call Unit_Tests_End_Group()

  ! Exponential truncation of galactic components.
  call Unit_Tests_Begin_Group("Galactic component truncation (exponential)")
  allocate(massDistributionHernquist :: massDistributionUntruncated_)
  select type (massDistributionUntruncated_)
  type is (massDistributionHernquist)
     massDistributionUntruncated_=massDistributionHernquist(dimensionless=.true.)
  end select
  ! Hold a reference to the untruncated distribution, so that it is not destroyed along with the truncated distributions which
  ! wrap it.
  !![
  <referenceCountIncrement object="massDistributionUntruncated_"/>
  !!]
  allocate(galacticComponentTruncationExponential :: galacticComponentTruncation_)
  select type (galacticComponentTruncation_)
  type is (galacticComponentTruncationExponential)
     galacticComponentTruncation_=galacticComponentTruncationExponential(radiusDecayOverRadiusTruncate=0.1d0)
  end select
  do i=1,truncationTableCount
     write (label,'(f4.2)') truncationTableFraction(i)
     massDistributionTruncated_ => galacticComponentTruncation_%massDistribution(massDistributionUntruncated_,truncationTableFraction(i))
     call Assert("Mass retained           (f="//label//")",massDistributionTruncated_  %massTotal            (                                                         ),truncationTableFraction            (i),relTol=1.0d-3)
     call Assert("Half-mass radius        (f="//label//")",massDistributionTruncated_  %radiusHalfMass       (                                                         ),truncationTableRadiusHalfMassTarget(i),relTol=1.0d-3)
     !![
     <objectDestructor name="massDistributionTruncated_"/>
     !!]
  end do
  ! A negligible truncation should leave the distribution unchanged.
  massDistributionTruncated_ => galacticComponentTruncation_%massDistribution(massDistributionUntruncated_,1.0d0-1.0d-6)
  call Assert("Negligible truncation (mass retained)"          ,massDistributionTruncated_  %massTotal            (                                               ),1.0d0,absTol=1.0d-12)
  !![
  <objectDestructor name="massDistributionTruncated_"  />
  <objectDestructor name="massDistributionUntruncated_"/>
  !!]
  deallocate(galacticComponentTruncation_)
  call Unit_Tests_End_Group()

  ! Heating of a tracer population along with a heated reference (dark matter) distribution.
  call Unit_Tests_Begin_Group("Heated tracer")
  block
    class           (massDistributionSpherical   ), pointer      :: darkMatter_                                    , darkMatterHeated_          , &
         &                                                          tracerNFW_                                     , tracerHernquist_           , &
         &                                                          tracerHeated_
    class           (massDistributionHeatingClass), pointer      :: heating_
    double precision                              , parameter    :: heatingNormalized=1.0d+7                       , massTracer          =1.0d10, &
         &                                                          fractionStep     =1.0d-2
    double precision                              , dimension(4) :: radii            =[3.0d-3,1.0d-2,3.0d-2,1.0d-1]
    double precision                                             :: radiusHalfMass                                 , radiusHalfMassHeated       , &
         &                                                          massDerivative
    type            (coordinateSpherical         )               :: coordinates_
    integer                                                      :: j
    character       (len=7                       )               :: labelRadius_

    ! Build an NFW dark matter distribution, and a heated version of it.
    allocate(massDistributionNFW :: darkMatter_)
    select type (darkMatter_)
    type is (massDistributionNFW)
       !![
       <referenceConstruct object="darkMatter_" constructor="massDistributionNFW(scaleLength=30.0d-3,radiusVirial=300.0d-3,mass=1.0d12,dimensionless=.false.)"/>
       !!]
    end select
    allocate(massDistributionHeatingTidal :: heating_)
    select type (heating_)
    type is (massDistributionHeatingTidal)
       !![
       <referenceConstruct object="heating_" constructor="massDistributionHeatingTidal(heatingNormalized,0.0d0,0.0d0,0.0d0,0.0d0)"/>
       !!]
    end select
    allocate(massDistributionSphericalHeated :: darkMatterHeated_)
    select type (darkMatterHeated_)
    type is (massDistributionSphericalHeated)
       !![
       <referenceConstruct object="darkMatterHeated_" constructor="massDistributionSphericalHeated(nonAnalyticSolversNumerical,.false.,.false.,.false.,1.0d-12,1.0d-3,darkMatter_,heating_)"/>
       !!]
    end select
    ! A tracer with the same profile as the dark matter must be heated identically to the dark matter.
    allocate(massDistributionNFW :: tracerNFW_)
    select type (tracerNFW_)
    type is (massDistributionNFW)
       !![
       <referenceConstruct object="tracerNFW_" constructor="massDistributionNFW(scaleLength=30.0d-3,radiusVirial=300.0d-3,mass=massTracer,dimensionless=.false.)"/>
       !!]
    end select
    allocate(massDistributionSphericalHeatedTracer :: tracerHeated_)
    select type (tracerHeated_)
    type is (massDistributionSphericalHeatedTracer)
       !![
       <referenceConstruct object="tracerHeated_" constructor="massDistributionSphericalHeatedTracer(tracerNFW_,darkMatterHeated_)"/>
       !!]
    end select
    do j=1,size(radii)
       write (labelRadius_,'(es7.1)') radii(j)
       coordinates_=[radii(j),0.0d0,0.0d0]
       call Assert("Same profile as dark matter: mass enclosed (r="//labelRadius_//")",tracerHeated_%massEnclosedBySphere(radii(j)    ),massTracer/1.0d12*darkMatterHeated_%massEnclosedBySphere(radii(j)    ),relTol=1.0d-6)
       call Assert("Same profile as dark matter: density       (r="//labelRadius_//")",tracerHeated_%density             (coordinates_),massTracer/1.0d12*darkMatterHeated_%density             (coordinates_),relTol=1.0d-6)
    end do
    !![
    <objectDestructor name="tracerHeated_"/>
    !!]
    ! A compact Hernquist tracer.
    allocate(massDistributionHernquist :: tracerHernquist_)
    select type (tracerHernquist_)
    type is (massDistributionHernquist)
       !![
       <referenceConstruct object="tracerHernquist_" constructor="massDistributionHernquist(dimensionless=.false.,mass=massTracer,scaleLength=5.0d-3)"/>
       !!]
    end select
    allocate(massDistributionSphericalHeatedTracer :: tracerHeated_)
    select type (tracerHeated_)
    type is (massDistributionSphericalHeatedTracer)
       !![
       <referenceConstruct object="tracerHeated_" constructor="massDistributionSphericalHeatedTracer(tracerHernquist_,darkMatterHeated_)"/>
       !!]
    end select
    radiusHalfMass      =tracerHernquist_%radiusHalfMass()
    radiusHalfMassHeated=tracerHeated_   %radiusHalfMass()
    call Assert("Hernquist tracer: mass is conserved"                      ,tracerHeated_%massTotal           (                    ),massTracer      ,relTol=1.0d-12)
    call Assert("Hernquist tracer: half-mass radius increases"             ,radiusHalfMassHeated > radiusHalfMass                   ,.true.                         )
    call Assert("Hernquist tracer: half-mass radius encloses half the mass",tracerHeated_%massEnclosedBySphere(radiusHalfMassHeated),0.5d0*massTracer,relTol=1.0d-4 )
    ! The density must be consistent with the derivative of the enclosed mass.
    do j=1,size(radii)
       write (labelRadius_,'(es7.1)') radii(j)
       coordinates_  =[radii(j),0.0d0,0.0d0]
       massDerivative=+(                                                                    &
            &           +tracerHeated_%massEnclosedBySphere(radii(j)*(1.0d0+fractionStep)) &
            &           -tracerHeated_%massEnclosedBySphere(radii(j)*(1.0d0-fractionStep)) &
            &          )                                                                    &
            &         /2.0d0                                                                &
            &         /fractionStep                                                         &
            &         /radii(j)
       call Assert("Hernquist tracer: density is consistent with enclosed mass (r="//labelRadius_//")",massDerivative,4.0d0*Pi*radii(j)**2*tracerHeated_%density(coordinates_),relTol=3.0d-3)
    end do
    !![
    <objectDestructor name="tracerHeated_"    />
    <objectDestructor name="darkMatterHeated_"/>
    <objectDestructor name="heating_"         />
    !!]
    ! With strong heating the outer shells of the dark matter (and so of the tracer) become unbound, and are moved to
    ! infinity. All of the tracer mass must still be enclosed by the largest radius, and the radius enclosing all but a little
    ! of the mass must be found without error.
    allocate(massDistributionHeatingTidal :: heating_)
    select type (heating_)
    type is (massDistributionHeatingTidal)
       !![
       <referenceConstruct object="heating_" constructor="massDistributionHeatingTidal(1.0d3*heatingNormalized,0.0d0,0.0d0,0.0d0,0.0d0)"/>
       !!]
    end select
    allocate(massDistributionSphericalHeated :: darkMatterHeated_)
    select type (darkMatterHeated_)
    type is (massDistributionSphericalHeated)
       !![
       <referenceConstruct object="darkMatterHeated_" constructor="massDistributionSphericalHeated(nonAnalyticSolversNumerical,.false.,.false.,.false.,1.0d-12,1.0d-3,darkMatter_,heating_)"/>
       !!]
    end select
    allocate(massDistributionSphericalHeatedTracer :: tracerHeated_)
    select type (tracerHeated_)
    type is (massDistributionSphericalHeatedTracer)
       !![
       <referenceConstruct object="tracerHeated_" constructor="massDistributionSphericalHeatedTracer(tracerHernquist_,darkMatterHeated_)"/>
       !!]
    end select
    call Assert("Strong heating: the outermost shells are unbound"            ,tracerHeated_%massEnclosedBySphere        (1.0d3                ) < massTracer,.true.                   )
    call Assert("Strong heating: all mass is enclosed by large radii"         ,tracerHeated_%massEnclosedBySphere        (2.0d0*radiusLarge    )             ,massTracer,relTol=1.0d-12)
    call Assert("Strong heating: radius enclosing 99% of the mass"            ,tracerHeated_%radiusEnclosingMass         (massFractional=0.99d0) > 0.0d0     ,.true.                   )
    call Assert("Strong heating: radius enclosing 99% of the mass (numerical)",tracerHeated_%radiusEnclosingMassNumerical(massFractional=0.99d0) > 0.0d0     ,.true.                   )
    !![
    <objectDestructor name="tracerHeated_"    />
    <objectDestructor name="darkMatterHeated_"/>
    <objectDestructor name="heating_"         />
    !!]
    ! With zero heating the tracer must be unchanged.
    allocate(massDistributionHeatingTidal :: heating_)
    select type (heating_)
    type is (massDistributionHeatingTidal)
       !![
       <referenceConstruct object="heating_" constructor="massDistributionHeatingTidal(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0)"/>
       !!]
    end select
    allocate(massDistributionSphericalHeated :: darkMatterHeated_)
    select type (darkMatterHeated_)
    type is (massDistributionSphericalHeated)
       !![
       <referenceConstruct object="darkMatterHeated_" constructor="massDistributionSphericalHeated(nonAnalyticSolversNumerical,.false.,.false.,.false.,1.0d-12,1.0d-3,darkMatter_,heating_)"/>
       !!]
    end select
    allocate(massDistributionSphericalHeatedTracer :: tracerHeated_)
    select type (tracerHeated_)
    type is (massDistributionSphericalHeatedTracer)
       !![
       <referenceConstruct object="tracerHeated_" constructor="massDistributionSphericalHeatedTracer(tracerHernquist_,darkMatterHeated_)"/>
       !!]
    end select
    call Assert("Zero heating: mass enclosed is unchanged",tracerHeated_%massEnclosedBySphere(radii(2)),tracerHernquist_%massEnclosedBySphere(radii(2)),relTol=1.0d-12)
    !![
    <objectDestructor name="tracerHeated_"    />
    <objectDestructor name="darkMatterHeated_"/>
    <objectDestructor name="heating_"         />
    <objectDestructor name="tracerHernquist_" />
    <objectDestructor name="tracerNFW_"       />
    <objectDestructor name="darkMatter_"      />
    !!]
  end block
  call Unit_Tests_End_Group()

  ! Radial moments. Reference values were computed by direct numerical quadrature (scipy.integrate.quad) of ρ(r) rᵐ (or
  ! Σ(R) Rᵐ) for the dimensionful profiles. Over the full range, ∫ 4π ρ r² dr = M and ∫ 2π Σ R dR = M.
  call Unit_Tests_Begin_Group("Radial moments")
  allocate(massDistributionHernquist :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionHernquist)
     massDistribution_=massDistributionHernquist(dimensionless=.false.,mass=2.0d0,scaleLength=2.0d0)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     call Assert("Hernquist: m=2, r=[0,∞)    ",massDistribution_%densityRadialMoment(2.0d0                                        ),2.0d0/4.0d0/Pi                                                   ,relTol=1.0d-6)
     call Assert("Hernquist: m=1, r=[0.5,3]  ",massDistribution_%densityRadialMoment(1.0d0,radiusMinimum=0.5d0,radiusMaximum=3.0d0),0.03819718634205488d0                                            ,relTol=1.0d-6)
     call Assert("Hernquist: m=2, r=[0.5,3]  ",massDistribution_%densityRadialMoment(2.0d0,radiusMinimum=0.5d0,radiusMaximum=3.0d0),0.050929581789406514d0                                           ,relTol=1.0d-6)
     call Assert("Hernquist: m=3, r=[0.5,3]  ",massDistribution_%densityRadialMoment(3.0d0,radiusMinimum=0.5d0,radiusMaximum=3.0d0),0.08476412777945765d0                                            ,relTol=1.0d-6)
     ! Small radii, where analytic solutions written as differences of terms of order unity would cancel catastrophically.
     call Assert("Hernquist: m=2, r=[0,1e-6rₛ]",massDistribution_%densityRadialMoment(2.0d0,radiusMinimum=0.0d0,radiusMaximum=2.0d-6),4.999990000015000d-13/Pi                                       ,relTol=1.0d-9)
     call Assert("Hernquist: m=3, r=[0,1e-6rₛ]",massDistribution_%densityRadialMoment(3.0d0,radiusMinimum=0.0d0,radiusMaximum=2.0d-6),2.0d0*3.3333258333453333d-19/Pi                                 ,relTol=1.0d-9)
     call Assert("Hernquist: m=1, r=[1,∞)    ",massDistribution_%densityRadialMoment(1.0d0,radiusMinimum=1.0d0                    ),0.0353677651315323d0                                             ,relTol=1.0d-6)
     call Assert("Hernquist: m=2.5, r=[0.5,3]",massDistribution_%densityRadialMoment(2.5d0,radiusMinimum=0.5d0,radiusMaximum=3.0d0),massDistribution_%densityRadialMomentNumerical(2.5d0,0.5d0,3.0d0),relTol=1.0d-3)
  end select
  allocate(massDistributionHernquist :: massDistributionMomentInner)
  select type (massDistributionMomentInner)
  type is (massDistributionHernquist)
     massDistributionMomentInner=massDistributionHernquist(dimensionless=.true.)
  end select
  allocate(massDistributionSphericalScaler :: massDistributionScalerTest)
  select type (massDistributionScalerTest)
  type is (massDistributionSphericalScaler)
     select type (massDistributionMomentInner)
     class is (massDistributionSpherical)
        massDistributionScalerTest=massDistributionSphericalScaler(factorScalingLength=2.0d0,factorScalingMass=2.0d0,massDistribution_=massDistributionMomentInner)
     end select
  end select
  select type (massDistributionScalerTest)
  class is (massDistributionSpherical)
     call Assert("Scaled Hernquist: m=2, r=[0,∞)  ",massDistributionScalerTest%densityRadialMoment(2.0d0                                        ),2.0d0/4.0d0/Pi       ,relTol=1.0d-6)
     call Assert("Scaled Hernquist: m=3, r=[0.5,3]",massDistributionScalerTest%densityRadialMoment(3.0d0,radiusMinimum=0.5d0,radiusMaximum=3.0d0),0.08476412777945765d0,relTol=1.0d-6)
  end select
  ! The scaler holds a counted reference to `massDistributionMomentInner` and destroys it in its own destructor, so only the
  ! scaler itself is deallocated here.
  deallocate(massDistribution_          )
  deallocate(massDistributionScalerTest )
  allocate(massDistributionExponentialDisk :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionExponentialDisk)
     massDistribution_=massDistributionExponentialDisk(mass=2.0d0,radiusScale=2.0d0,scaleHeight=0.01d0)
  end select
  call Assert("Exponential disk: m=1, R=[0,∞)",massDistribution_%surfaceDensityRadialMoment(1.0d0                                        ),2.0d0/2.0d0/Pi      ,relTol=1.0d-6)
  call Assert("Exponential disk: m=2, R=[1,3]",massDistribution_%surfaceDensityRadialMoment(2.0d0,radiusMinimum=1.0d0,radiusMaximum=3.0d0),0.2250648139159444d0,relTol=1.0d-6)
  allocate(massDistributionExponentialDisk :: massDistributionMomentInner)
  select type (massDistributionMomentInner)
  type is (massDistributionExponentialDisk)
     massDistributionMomentInner=massDistributionExponentialDisk(scaleHeight=0.005d0,dimensionless=.true.)
  end select
  allocate(massDistributionCylindricalScaler :: massDistributionScalerTest)
  select type (massDistributionScalerTest)
  type is (massDistributionCylindricalScaler)
     select type (massDistributionMomentInner)
     class is (massDistributionCylindrical)
        massDistributionScalerTest=massDistributionCylindricalScaler(factorScalingLength=2.0d0,factorScalingMass=2.0d0,massDistribution_=massDistributionMomentInner)
     end select
  end select
  call Assert("Scaled exponential disk: m=1, R=[0,∞)",massDistributionScalerTest%surfaceDensityRadialMoment(1.0d0                                        ),2.0d0/2.0d0/Pi      ,relTol=1.0d-6)
  call Assert("Scaled exponential disk: m=2, R=[1,3]",massDistributionScalerTest%surfaceDensityRadialMoment(2.0d0,radiusMinimum=1.0d0,radiusMaximum=3.0d0),0.2250648139159444d0,relTol=1.0d-6)
  ! The scaler holds a counted reference to `massDistributionMomentInner` and destroys it in its own destructor, so only the
  ! scaler itself is deallocated here.
  deallocate(massDistribution_          )
  deallocate(massDistributionScalerTest )
  call Unit_Tests_End_Group()

  ! Sersic profile.
  call Unit_Tests_Begin_Group("Sersic profile")
  allocate(massDistributionSersic :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionSersic)
     massDistribution_=massDistributionSersic(index=4.0d0,dimensionless=.true.)
     radiusInProjection=massDistribution_%radiusHalfMassProjected()
  class default
     radiusInProjection=0.0d0
     call Error_Report('unknown mass distribution'//{introspection:location})
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     position=[1.0d0,0.0d0,0.0d0]
     call Assert("Half mass radius",massDistribution_%radiusHalfMass(),1.0d0,absTol=1.0d-6)
     ! Tabulate the mass and potential.
     do i=1,sersicTableCount
        position=[sersicTableRadius(i)*radiusInProjection,0.0d0,0.0d0]
        sersicTableMass     (i)=massDistribution_%massEnclosedBySphere(sersicTableRadius(i)*radiusInProjection)
        sersicTableDensity  (i)=massDistribution_%density             (position                               )
        sersicTablePotential(i)=massDistribution_%potential           (position                               )
     end do
     ! Convert to reduced potential used by Young (1976).
     sersicTablePotential=sersicTablePotential/sersicTablePotential(1)
     ! Convert density to effective radius units.
     sersicTableDensity=sersicTableDensity*radiusInProjection**3
     call Assert("Density vs radius"  ,sersicTableDensity  ,sersicTableDensityTarget  ,relTol=1.0d-3)
     call Assert("Mass vs radius"     ,sersicTableMass     ,sersicTableMassTarget     ,relTol=2.1d-2)
     call Assert("Potential vs radius",sersicTablePotential,sersicTablePotentialTarget,relTol=1.0d-2)
  end select
  deallocate(massDistribution_)
  call Unit_Tests_End_Group()

  ! Beta profile.
  call Unit_Tests_Begin_Group("Beta profile")
  call Unit_Tests_Begin_Group("β = 2/3")
  allocate(massDistributionBetaProfile :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionBetaProfile)
     massDistribution_=massDistributionBetaProfile(beta=2.0d0/3.0d0,dimensionless=.true.)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     position    =[1.0d0,0.0d0,0.0d0]
     positionZero=[0.0d0,0.0d0,0.0d0]
     call Assert(                                                              &
          &      "Mass within scale radius"                                  , &
          &      +massDistribution_%massEnclosedBySphere (1.0d0             ), &
          &      +(4.0d0-Pi)*Pi                                              , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Potential at scale radius"                                 , &
          &      +massDistribution_%potential            (position          )  &
          &      -massDistribution_%potential            (positionZero      ), &
          &      +Pi*(Pi+log(4.0d0)-4.0d0)                                   , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Radial moment, m=1, from 0 to 1"                           , &
          &      +massDistribution_%densityRadialMoment  (1.0d0,0.0d0,1.0d0) , &
          &      +0.3465735902799727d0                                       , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Radial moment, m=2, from 0 to 1"                           , &
          &      +massDistribution_%densityRadialMoment  (2.0d0,0.0d0,1.0d0) , &
          &      +0.2146018366025517d0                                       , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Radial moment, m=3, from 0 to 1"                           , &
          &      +massDistribution_%densityRadialMoment  (3.0d0,0.0d0,1.0d0) , &
          &      +0.1534264097200273d0                                       , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Square integral"                                           , &
          &      +massDistribution_%densitySquareIntegral(      0.0d0,1.0d0) , &
          &      +1.793209546954886d0                                        , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Density gradient (logarithmic)"                            , &
          &      +massDistribution_%densityGradientRadial(position,.true.  ) , &
          &      -1.000000000000000d0                                        , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Density gradient (non-logarithmic)"                        , &
          &      +massDistribution_%densityGradientRadial(position,.false. ) , &
          &      +massDistribution_%densityGradientRadial(position,.true.  )   &
          &      *massDistribution_%density              (position         ) , &
          &      absTol=1.0d-6                                                 &
          &     )
  end select
  call Unit_Tests_End_Group()
  deallocate(massDistribution_)
  call Unit_Tests_Begin_Group("β = 2")
  allocate(massDistributionBetaProfile :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionBetaProfile)
     massDistribution_=massDistributionBetaProfile(beta=2.0d0,dimensionless=.true.)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     call Assert(                                                              &
          &      "Radial moment, m=1, unbounded"                             , &
          &      +massDistribution_%densityRadialMoment  (1.0d0            ) , &
          &      +1.0d0/4.0d0                                                , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Radial moment, m=2, unbounded"                             , &
          &      +massDistribution_%densityRadialMoment  (2.0d0            ) , &
          &      +Pi/16.0d0                                                  , &
          &      absTol=1.0d-6                                                 &
          &     )
     call Assert(                                                              &
          &      "Radial moment, m=3, unbounded"                             , &
          &      +massDistribution_%densityRadialMoment  (3.0d0            ) , &
          &      +1.0d0/4.0d0                                                , &
          &      absTol=1.0d-6                                                 &
          &     )
  end select
  call Unit_Tests_End_Group()
  ! Ensure that a dimensionful profile produces the correct mass inside of its outer radius.
  deallocate(massDistribution_)
  allocate(massDistributionBetaProfile :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionBetaProfile)
     massDistribution_=massDistributionBetaProfile(beta=2.0d0/3.0d0,radiusCore=3.8492316686261747d-002,mass=182582297.19568533d0,radiusOuter=0.12829569846196026d0)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     call Assert("Mass within outer radius",massDistribution_%massEnclosedBySphere(0.12829569846196026d0),182582297.19568533d0,relTol=1.0d-6)
  end select
  deallocate(massDistribution_)
  ! Ensure that the potential of a dimensionful profile truncated at its outer radius matches (Φ(r)-Φ(0))/G=∫₀ʳ M(r′)/r′² dr′,
  ! found by direct numerical integration of the density (using mpmath), at radii well inside the core, inside the outer radius,
  ! and beyond the outer radius (where the potential must be that of a point mass).
  block
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    double precision          , parameter                  :: radiusCore     =3.8492316686261747d-002, radiusOuter=0.12829569846196026d0, &
         &                                                    mass           =182582297.19568533d0
    double precision          , dimension(2  ), parameter  :: betas          =[2.0d0/3.0d0,0.9d0]
    character       (len=4   ), dimension(2  ), parameter  :: labelsBeta     =['2/3 ','0.9 ']
    double precision          , dimension(4  ), parameter  :: radii          =[5.0d-4*radiusCore,0.5d0*radiusCore,radiusOuter,2.0d0*radiusOuter]
    character       (len=11  ), dimension(4  ), parameter  :: labelsRadius   =['5e-4 r_core','0.5 r_core ','r_outer    ','2 r_outer  ']
    double precision          , dimension(4,2), parameter  :: potentialsTarget=reshape(                                                                                  &
         &                                                                             [                                                                                &
         &                                                                              96.235390151914037d0,89769133.903370879d0,1456940246.0073782d0,2168508519.0848309d0, &
         &                                                                             165.60492933421483d0,150858194.57747654d0,1882837097.8388988d0,2594405370.9163515d0  &
         &                                                                             ]                                                                              , &
         &                                                                             [4,2]                                                                            &
         &                                                                            )
    integer                                                :: j
    do i=1,2
       allocate(massDistributionBetaProfile :: massDistribution_)
       select type (massDistribution_)
       type is (massDistributionBetaProfile)
          massDistribution_=massDistributionBetaProfile(beta=betas(i),radiusCore=radiusCore,mass=mass,radiusOuter=radiusOuter,truncateAtOuterRadius=.true.)
       end select
       positionZero=[0.0d0,0.0d0,0.0d0]
       do j=1,4
          position=[radii(j),0.0d0,0.0d0]
          call Assert(                                                                                                        &
               &      "Truncated, dimensionful profile, β="//trim(labelsBeta(i))//", potential at r="//trim(labelsRadius(j)), &
               &      +(                                                                                                      &
               &        +massDistribution_%potential(position    )                                                            &
               &        -massDistribution_%potential(positionZero)                                                            &
               &       )                                                                                                      &
               &      /gravitationalConstant_internal                                                                       , &
               &      potentialsTarget(j,i)                                                                                 , &
               &      relTol=1.0d-6                                                                                           &
               &     )
       end do
       deallocate(massDistribution_)
    end do
  end block
  call Unit_Tests_End_Group()

  ! Exponential disk profile. Tests are carried out for a thin disk, and for a thicker disk with a scale height typical of those
  ! used in models. Reference values were computed by the `exponentialDiskAccelerations.py` script in the galacticusDevTools
  ! repository (physicsReferences/), with the `--scaleHeight` option set to the scale height of each disk.
  call Unit_Tests_Begin_Group("Exponential disk profile")
  do indexDisk=1,2
     select case (indexDisk)
     case (1)
        scaleHeight            =0.010d0
        accelerationReferences =reshape(                                                &
             &                          [                                               &
             &                           -2.73856896d-01,0.0d0,-2.77960050d-01,         &
             &                           -2.43542857d-01,0.0d0,-3.45809193d-01,         &
             &                           -2.55934175d-01,0.0d0,-5.83951896d-01,         &
             &                           -1.55488612d-01,0.0d0,-1.30061903d-01          &
             &                          ]                                             , &
             &                          [3,4]                                           &
             &                         )
        potentialReferences    =[-7.57472032d-01,-7.25924201d-01,-8.85012212d-01,-4.84181094d-01,-7.59056749d-01,-9.92098271d-01]
        rotationCurveReferences=[+1.11813634d-01,+3.66925192d-01,+5.24812965d-01,+6.20666151d-01]
     case (2)
        scaleHeight            =0.137d0
        accelerationReferences =reshape(                                                &
             &                          [                                               &
             &                           -2.46408607d-01,0.0d0,-2.46220802d-02,         &
             &                           -2.36349877d-01,0.0d0,-2.07432588d-01,         &
             &                           -2.24496978d-01,0.0d0,-1.12739450d-01,         &
             &                           -1.55813682d-01,0.0d0,-1.26417795d-01          &
             &                          ]                                             , &
             &                          [3,4]                                           &
             &                         )
        potentialReferences    =[-7.28232193d-01,-7.17065216d-01,-8.48335481d-01,-4.84137402d-01,-7.28355422d-01,-9.23329623d-01]
        rotationCurveReferences=[+9.18654861d-02,+3.36088741d-01,+4.96508338d-01,+6.03801051d-01]
     end select
     write (labelHeight,'(f5.3)') scaleHeight
     call Unit_Tests_Begin_Group("Scale height h="//labelHeight)
     allocate(massDistributionExponentialDisk :: massDistribution_)
     select type (massDistribution_)
     type is (massDistributionExponentialDisk)
        massDistribution_=massDistributionExponentialDisk(scaleHeight=scaleHeight,dimensionless=.true.)
        ! Test that the mass within a cylindrical radius matches analytical results.
        call Unit_Tests_Begin_Group("Mass within cylindrical radius matches analytical results")
        call Assert("Mass enclosed by cylinder at r=½",massDistribution_%massEnclosedByCylinder(0.5d0),(2.0d0-3.0d0/sqrt(e))/2.0d0,relTol=1.0d-6)
        call Assert("Mass enclosed by cylinder at r=1",massDistribution_%massEnclosedByCylinder(1.0d0),(e-2.0d0)/e                ,relTol=1.0d-6)
        call Assert("Mass enclosed by cylinder at r=2",massDistribution_%massEnclosedByCylinder(2.0d0),1.0d0-3.0d0/e**2           ,relTol=1.0d-6)
        call Unit_Tests_End_Group()
        ! Test that the rotation curve gradient matches a finite-difference estimate.
        call Unit_Tests_Begin_Group("Rotation curve gradient matches finite-difference estimate")
        do i=1,3
           select case (i)
           case (1)
              radius=0.5d0
           case (2)
              radius=1.0d0
           case (3)
              radius=4.0d0
           end select
           write (label,'(f4.1)') radius
           rotationCurveGradientAnalytic =+  massDistribution_%rotationCurveGradient(radius=radius                                )
           rotationCurveGradientNumerical=+(                                                                                           &
                &                           +massDistribution_%rotationCurve        (radius=radius*(1.0d0+epsilonFiniteDifference))**2 &
                &                           -massDistribution_%rotationCurve        (radius=radius                                )**2 &
                &                           )                                                                                          &
                &                           /epsilonFiniteDifference                                                                   &
                &                           /radius
           call Assert("Rotation curve gradient at r="//trim(adjustl(label)),rotationCurveGradientAnalytic,rotationCurveGradientNumerical,relTol=5.0d0*epsilonFiniteDifference)
        end do
        call Unit_Tests_End_Group()
        ! Test that gravitational acceleration matches expectations for a point mass distribution at large radii.
        call Unit_Tests_Begin_Group("Acceleration approaches point mass solution at large radii")
        do i=1,2
           select case (i)
           case (1)
              radius  =45.0d0
           case (2)
              radius  =55.0d0
           end select
           write (label,'(f4.1)') radius
           ! Vertically above the disk.
           position=[radius,0.0d0,0.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",0  , 0  ]",accelerationReference,[0.0d0,0.0d0,-1.0d0/radius**2],absTol=1.0d-6,relTol=1.0d-3)
           ! Vertically below the disk.
           position=[radius,Pi,0.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π  , 0  ]",accelerationReference,[0.0d0,0.0d0,+1.0d0/radius**2],absTol=1.0d-6,relTol=1.0d-3)
           ! Disk plane, along +x-axis.
           position=[radius,Pi/2.0d0,0.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",accelerationReference,[-1.0d0/radius**2,0.0d0,0.0d0],absTol=1.0d-6,relTol=1.0d-2)
           ! Disk plane, along -y-axis.
           position=[radius,Pi/2.0d0,-Pi/2.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",accelerationReference,[0.0d0,+1.0d0/radius**2,0.0d0],absTol=1.0d-6,relTol=1.0d-2)
        end do
        call Unit_Tests_End_Group()
        ! Test that gravitational tidal tensor matches expectations for a point mass distribution at large radii.
        call Unit_Tests_Begin_Group("Tidal tensor approaches point mass solution at large radii")
        do i=1,2
           select case (i)
           case (1)
              radius  =45.0d0
           case (2)
              radius  =55.0d0
           end select
           write (label,'(f4.1)') radius
           ! Vertically above the disk.
           position                      =[radius,0.0d0,0.0d0]
           tidalTensorSphericalComponents=reshape([-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,+2.0d0]/radius**3,[3,3])
           tidalTensorComponents         =massDistribution_%tidalTensor         (position)
           call Assert("[r,θ,φ] = ["//trim(label)//",0  , 0  ]",tidalTensorComponents,tidalTensorSphericalComponents,absTol=1.0d-6,relTol=3.0d-3)
           ! Vertically below the disk.
           position=[radius,Pi,0.0d0]
           tidalTensorSphericalComponents=reshape([-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,+2.0d0]/radius**3,[3,3])
           tidalTensorComponents         =massDistribution_%tidalTensor         (position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π  , 0  ]",tidalTensorComponents,tidalTensorSphericalComponents,absTol=1.0d-6,relTol=3.0d-3)
           ! Disk plane, along +x-axis.
           position                      =[radius,Pi/2.0d0,0.0d0]
           tidalTensorSphericalComponents=reshape([+2.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0]/radius**3,[3,3])
           tidalTensorComponents         =massDistribution_%tidalTensor         (position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",tidalTensorComponents,tidalTensorSphericalComponents,absTol=1.0d-6,relTol=2.0d-2)
           ! Disk plane, along -y-axis.
           position                      =[radius,Pi/2.0d0,-Pi/2.0d0]
           tidalTensorSphericalComponents=reshape([-1.0d0,0.0d0,0.0d0,0.0d0,+2.0d0,0.0d0,0.0d0,0.0d0,-1.0d0]/radius**3,[3,3])
           tidalTensorComponents         =massDistribution_%tidalTensor         (position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",tidalTensorComponents,tidalTensorSphericalComponents,absTol=1.0d-6,relTol=1.0d-2)
        end do
        call Unit_Tests_End_Group()
        ! Test that the gravitational acceleration in the mid-plane is consistent with the rotation curve. These are tabulated
        ! separately (the acceleration on a grid in radius and height, the rotation curve in radius in the mid-plane).
        call Unit_Tests_Begin_Group("Mid-plane acceleration is consistent with rotation curve")
        do i=1,2
           select case (i)
           case (1)
              radius=0.1d0
           case (2)
              radius=1.0d0
           end select
           write (label,'(f4.1)') radius
           ! Disk plane, along +x-axis.
           position=[radius,Pi/2.0d0,0.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",accelerationReference,[-1.0d0,0.0d0,0.0d0]*massDistribution_%rotationCurve(radius)**2/radius,absTol=1.0d-6,relTol=1.0d-3)
           ! Disk plane, along -y-axis.
           position=[radius,Pi/2.0d0,-Pi/2.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",accelerationReference,[0.0d0,+1.0d0,0.0d0]*massDistribution_%rotationCurve(radius)**2/radius,absTol=1.0d-6,relTol=1.0d-3)
        end do
        call Unit_Tests_End_Group()
        ! Test that the tabulated tidal tensor satisfies the Poisson equation. The tidal tensor is T = -∇∇Φ, so its trace must equal
        ! -4πGρ (with G=1 for this dimensionless disk). This is tested close to the disk plane, where the vertical structure of the
        ! disk dominates the trace.
        call Unit_Tests_Begin_Group("Tidal tensor satisfies the Poisson equation")
        do i=1,4
           select case (i)
           case (1)
              radiusCylindrical=0.5d0
              height           =0.010d0
           case (2)
              radiusCylindrical=1.0d0
              height           =0.010d0
           case (3)
              radiusCylindrical=1.0d0
              height           =0.020d0
           case (4)
              radiusCylindrical=2.0d0
              height           =0.015d0
           end select
           write (label      ,'(f3.1)') radiusCylindrical
           write (labelHeight,'(f5.3)') height
           position             =[sqrt(radiusCylindrical**2+height**2),atan2(radiusCylindrical,height),0.0d0]
           tidalTensorComponents=massDistribution_%tidalTensor(position)
           call Assert(                                               &
                &      "R="//trim(label)//", z="//labelHeight       , &
                &      +tidalTensorComponents(1,1)                    &
                &      +tidalTensorComponents(2,2)                    &
                &      +tidalTensorComponents(3,3)                  , &
                &      -4.0d0*Pi*massDistribution_%density(position), &
                &      relTol=5.0d-2                                  &
                &     )
        end do
        call Unit_Tests_End_Group()
        ! Test the acceleration against reference values. These were computed as Hankel transforms of the surface density, with the
        ! vertical kernel evaluated by direct quadrature of the convolution of the sech² vertical density profile with exp(-k|z-z'|) - i.e.
        ! independently of the form of the kernel used in the tabulation. The tolerance allows for the bilinear interpolation in the
        ! tabulated solution.
        call Unit_Tests_Begin_Group("Acceleration matches reference values")
        do i=1,4
           select case (i)
           case (1)
              radiusCylindrical=1.0d0
              height           =0.01d0
           case (2)
              radiusCylindrical=1.0d0
              height           =0.10d0
           case (3)
              radiusCylindrical=0.5d0
              height           =0.03d0
           case (4)
              radiusCylindrical=2.0d0
              height           =0.30d0
           end select
           write (label      ,'(f3.1)') radiusCylindrical
           write (labelHeight,'(f5.3)') height
           position    =[sqrt(radiusCylindrical**2+height**2),atan2(radiusCylindrical,height),0.0d0]
           acceleration=massDistribution_%acceleration(position)
           call Assert("R="//trim(label)//", z="//labelHeight,acceleration,accelerationReferences(:,i),absTol=1.0d-6,relTol=2.0d-3)
        end do
        call Unit_Tests_End_Group()
        ! Test the potential against reference values. The
        ! final two points lie below the smallest height (and, for the last, radius) tabulated, and so test the extension of the
        ! potential to those regions. The tolerance allows for the bilinear interpolation in the tabulated solution.
        call Unit_Tests_Begin_Group("Potential matches reference values")
        do i=1,6
           select case (i)
           case (1)
              radiusCylindrical=1.000d0
              height           =0.010d0
           case (2)
              radiusCylindrical=1.000d0
              height           =0.100d0
           case (3)
              radiusCylindrical=0.500d0
              height           =0.030d0
           case (4)
              radiusCylindrical=2.000d0
              height           =0.300d0
           case (5)
              radiusCylindrical=1.000d0
              height           =0.000d0
           case (6)
              radiusCylindrical=0.005d0
              height           =0.005d0
           end select
           write (labelRadius,'(f5.3)') radiusCylindrical
           write (labelHeight,'(f5.3)') height
           position=[sqrt(radiusCylindrical**2+height**2),atan2(radiusCylindrical,height),0.0d0]
           call Assert("R="//labelRadius//", z="//labelHeight,massDistribution_%potential(position),potentialReferences(i),relTol=1.0d-3)
        end do
        call Unit_Tests_End_Group()
        ! Test that the gradient of the potential, estimated by finite differences, is consistent with the acceleration. The potential
        ! is interpolated bilinearly (in logarithmic radius and height) in its tabulation, so that its gradient within a single cell
        ! of the table is only first-order accurate (in error by up to 2%). A step spanning several
        ! cells is therefore used for the finite difference estimate.
        call Unit_Tests_Begin_Group("Potential gradient is consistent with acceleration")
        do i=1,3
           select case (i)
           case (1)
              radiusCylindrical=1.0d0
              height           =0.10d0
           case (2)
              radiusCylindrical=0.5d0
              height           =0.05d0
           case (3)
              radiusCylindrical=2.0d0
              height           =0.30d0
           end select
           write (labelRadius,'(f5.3)') radiusCylindrical
           write (labelHeight,'(f5.3)') height
           position                 =[sqrt(radiusCylindrical**2+height**2),atan2(radiusCylindrical,height),0.0d0]
           acceleration             =massDistribution_%acceleration(position)
           position                 =[sqrt((radiusCylindrical*(1.0d0+epsilonFiniteDifferencePotential))**2+height**2),atan2(radiusCylindrical*(1.0d0+epsilonFiniteDifferencePotential),height),0.0d0]
           potentialGradientRadial  =+massDistribution_%potential(position)
           position                 =[sqrt((radiusCylindrical*(1.0d0-epsilonFiniteDifferencePotential))**2+height**2),atan2(radiusCylindrical*(1.0d0-epsilonFiniteDifferencePotential),height),0.0d0]
           potentialGradientRadial  =+(                                       &
                &                      +potentialGradientRadial               &
                &                      -massDistribution_%potential(position) &
                &                     )                                       &
                &                    /2.0d0                                   &
                &                    /epsilonFiniteDifferencePotential        &
                &                    /radiusCylindrical
           position                 =[sqrt(radiusCylindrical**2+(height*(1.0d0+epsilonFiniteDifferencePotential))**2),atan2(radiusCylindrical,height*(1.0d0+epsilonFiniteDifferencePotential)),0.0d0]
           potentialGradientVertical=+massDistribution_%potential(position)
           position                 =[sqrt(radiusCylindrical**2+(height*(1.0d0-epsilonFiniteDifferencePotential))**2),atan2(radiusCylindrical,height*(1.0d0-epsilonFiniteDifferencePotential)),0.0d0]
           potentialGradientVertical=+(                                       &
                &                      +potentialGradientVertical             &
                &                      -massDistribution_%potential(position) &
                &                     )                                       &
                &                    /2.0d0                                   &
                &                    /epsilonFiniteDifferencePotential        &
                &                    /height
           call Assert("R="//labelRadius//", z="//labelHeight,[-potentialGradientRadial,-potentialGradientVertical],[acceleration(1),acceleration(3)],relTol=1.0d-2)
        end do
        call Unit_Tests_End_Group()
        ! Test the mid-plane rotation curve against reference values.
        call Unit_Tests_Begin_Group("Rotation curve matches reference values")
        do i=1,4
           select case (i)
           case (1)
              radius=0.1d0
           case (2)
              radius=0.5d0
           case (3)
              radius=1.0d0
           case (4)
              radius=2.2d0
           end select
           write (label,'(f4.1)') radius
           call Assert("R="//trim(adjustl(label)),massDistribution_%rotationCurve(radius),rotationCurveReferences(i),relTol=1.0d-5)
        end do
        call Unit_Tests_End_Group()
     class default
        call Error_Report('unknown mass distribution'//{introspection:location})
     end select
     deallocate(massDistribution_)
     call Unit_Tests_End_Group()
  end do
  call Unit_Tests_End_Group()

  ! Gaussian ellipsoid profile.
  call Unit_Tests_Begin_Group("Gaussian ellipsoid profile")
  allocate(massDistributionGaussianEllipsoid :: massDistribution_      )
  allocate(massDistributionGaussianEllipsoid :: massDistributionRotated)
  select type (massDistribution_)
  type is (massDistributionGaussianEllipsoid)
     select type (massDistributionRotated)
     type is (massDistributionGaussianEllipsoid)
        ! Mass distribution aligned with the principle Cartesian axes.
        allocate(axes(3))
        axes(1)=vector([+1.0d0,+0.0d0,+0.0d0])
        axes(2)=vector([+0.0d0,+1.0d0,+0.0d0])
        axes(3)=vector([+0.0d0,+0.0d0,+1.0d0])
        massDistribution_      =massDistributionGaussianEllipsoid(scaleLength=[1.0d0,0.5d0,1.0d0],axes=axes,dimensionless=.true.)
        deallocate(axes)
        ! An identical mass distribution rotated by +π/2 around the z-axis.
        allocate(axes(3))
        axes(1)=vector([+0.0d0,-1.0d0,+0.0d0])
        axes(2)=vector([+1.0d0,+0.0d0,+0.0d0])
        axes(3)=vector([+0.0d0,+0.0d0,+1.0d0])
        massDistributionRotated=massDistributionGaussianEllipsoid(scaleLength=[1.0d0,0.5d0,1.0d0],axes=axes,dimensionless=.true.)
        deallocate(axes)
        ! Test that gravitational acceleration matches expectations for a point mass distribution at large radii.
        call Unit_Tests_Begin_Group("Acceleration approaches point mass solution at large radii")
        do i=1,2
           select case (i)
           case (1)
              radius  =45.0d0
           case (2)
              radius  =55.0d0
           end select
           write (label,'(f4.1)') radius
           ! Along the +z axis.
           position=[radius,0.0d0,0.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",0  , 0  ]",accelerationReference,[0.0d0,0.0d0,-1.0d0/radius**2],absTol=1.0d-6,relTol=3.0d-2)
           ! Along the -z axis.
           position=[radius,Pi   ,0.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π  , 0  ]",accelerationReference,[0.0d0,0.0d0,+1.0d0/radius**2],absTol=1.0d-6,relTol=3.0d-2)
           ! x-y plane, along +x-axis.
           position                      =[radius,Pi/2.0d0,0.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",accelerationReference,[-1.0d0/radius**2,0.0d0,0.0d0],absTol=1.0d-6,relTol=3.0d-2)
           ! x-y plane, along -y-axis.
           position                      =[radius,Pi/2.0d0,-Pi/2.0d0]
           accelerationReference=massDistribution_%acceleration(position)
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",accelerationReference,[0.0d0,+1.0d0/radius**2,0.0d0],absTol=1.0d-6,relTol=3.0d-2)
        end do
        call Unit_Tests_End_Group()
        ! Test that acceleration due to rotated ellipsoid matches the rotated acceleration of an unrotated ellipsoid.
        call Unit_Tests_Begin_Group("Acceleration of rotated ellipsoid")
        position    =[1.0d0,Pi/2.0d0,3.0d0*Pi/4.0d0]
        acceleration=massDistributionRotated%acceleration(position)     ! Evaluate the acceleration in the rotated distribution at the rotated position.
        acceleration=[acceleration(2),-acceleration(1),acceleration(3)] ! De-rotate the acceleration due to the rotated distribution.
        position   =[1.0d0,Pi/2.0d0,Pi/4.0d0]                           ! Evaluate the acceleration in the un-rotated distribution at the un-rotated position.
        accelerationReference=massDistribution_%acceleration(position)
        call Assert("[r,θ,φ] = [1,π/2,π/4]",accelerationReference,acceleration,absTol=1.0d-6,relTol=1.0d-6)
        call Unit_Tests_End_Group()
        class default
        call Error_Report('unknown mass distribution'//{introspection:location})
     end select
     class default
     call Error_Report('unknown mass distribution'//{introspection:location})
  end select
  deallocate(massDistribution_      )
  deallocate(massDistributionRotated)
  allocate(massDistributionGaussianEllipsoid :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionGaussianEllipsoid)
     allocate(axes(3))
     axes(1)=vector([1.0d0,0.0d0,0.0d0])
     axes(2)=vector([0.0d0,1.0d0,0.0d0])
     axes(3)=vector([0.0d0,0.0d0,1.0d0])
     massDistribution_=massDistributionGaussianEllipsoid(scaleLength=[1.0d0,1.0d0,1.0d0],axes=axes,dimensionless=.true.)
     deallocate(axes)
     ! Test that gravitational acceleration matches expectations for a Gaussian spheroid.
     call Unit_Tests_Begin_Group("Acceleration matches expectation for Gaussian spheroid")
     do i=1,2
        select case (i)
        case (1)
           radius  =1.0d0
        case (2)
           radius  =2.0d0
        end select
        write (label,'(f4.1)') radius
        ! Compute the enclosed mass fraction for a Gaussian spheroid.
        massFraction=+Error_Function(radius/sqrt(2.0d0)) &
             &       -sqrt(2.0d0/Pi)                     &
             &       *radius*exp(-0.5d0*radius**2)
        ! Along the +z axis.
        position=[radius,0.0d0,0.0d0]
        accelerationReference=massDistribution_%acceleration(position)
        call Assert("[r,θ,φ] = ["//trim(label)//",0  , 0  ]",accelerationReference,[0.0d0,0.0d0,-massFraction/radius**2],absTol=1.0d-6,relTol=2.0d-2)
        ! Along the -z axis.
        position=[radius,Pi   ,0.0d0]
        accelerationReference=massDistribution_%acceleration(position)
        call Assert("[r,θ,φ] = ["//trim(label)//",π  , 0  ]",accelerationReference,[0.0d0,0.0d0,+massFraction/radius**2],absTol=1.0d-6,relTol=2.0d-2)
        ! x-y plane, along +x-axis.
        position                      =[radius,Pi/2.0d0,0.0d0]
        accelerationReference=massDistribution_%acceleration(position)
        call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",accelerationReference,[-massFraction/radius**2,0.0d0,0.0d0],absTol=1.0d-6,relTol=2.0d-2)
        ! x-y plane, along -y-axis.
        position                      =[radius,Pi/2.0d0,-Pi/2.0d0]
        accelerationReference=massDistribution_%acceleration(position)
        call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",accelerationReference,[0.0d0,+massFraction/radius**2,0.0d0],absTol=1.0d-6,relTol=2.0d-2)
     end do
     call Unit_Tests_End_Group()
  class default
     call Error_Report('unknown mass distribution'//{introspection:location})
  end select
  deallocate(massDistribution_)
  call Unit_Tests_End_Group()
  
  ! Composite profile.
  call Unit_Tests_Begin_Group("Composite profile (Hernquist + Exponential disk)")
  allocate(                                   massDistributions                       )
  allocate(                                   massDistributions%next                  )
  allocate(massDistributionHernquist       :: massDistributions     %massDistribution_)
  allocate(massDistributionExponentialDisk :: massDistributions%next%massDistribution_)
  select type (massDistribution_ => massDistributions     %massDistribution_)
  type is (massDistributionHernquist      )
     massDistribution_=massDistributionHernquist      (mass=9.0d09,scaleLength=0.7d-3                   ,componentType=componentTypeSpheroid)
  end select
  select type (massDistribution_ => massDistributions%next%massDistribution_)
  type is (massDistributionExponentialDisk)
     massDistribution_=massDistributionExponentialDisk(mass=5.7d10,radiusScale=2.9d-3,scaleHeight=0.3d-3,componentType=componentTypeDisk    )
  end select
  allocate(massDistributionComposite :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionComposite)
     massDistribution_=massDistributionComposite(massDistributions)
  end select
  symmetry_=massDistribution_%symmetry()
  call Assert("Maximal symmetry [cylindrical]",symmetry_%ID,massDistributionSymmetryCylindrical%ID)
  massDistributionDisk_     => massDistribution_%subset(componentType=componentTypeDisk    )
  massDistributionSpheroid_ => massDistribution_%subset(componentType=componentTypeSpheroid)
  positionCartesian         =  [1.0d-3,0.0d0,0.0d0]
  call Assert("Density at (x,y,z)=(1,0,0) kpc"         ,massDistribution_        %density(positionCartesian),1.47756308872d18,relTol=1.0d-6)
  call Assert("Spheroid density at (x,y,z)=(1,0,0) kpc",massDistributionSpheroid_%density(positionCartesian),2.04086330446d17,relTol=1.0d-6)
  call Assert("Disk density at (x,y,z)=(1,0,0) kpc"    ,massDistributionDisk_    %density(positionCartesian),1.27347675828d18,relTol=1.0d-6)
  ! Mass within cylindrical radius.
  call Assert("Mass within cylinder of radius 1.3 kpc" ,massDistribution_%massEnclosedByCylinder(1.3d-3),4.557847074d9+4.271864472d9,relTol=1.0d-6)
  ! Mass within cylindrical radius.
  call Assert("Cylindrical half-mass radius"           ,massDistribution_%massEnclosedByCylinder(massDistribution_%radiusCylindricalEnclosingMass(massFractional=0.50d0)),0.5d0*(9.0d9+5.7d10),relTol=1.0d-6)
  nullify   (massDistributions)
  deallocate(massDistribution_)
  !![
  <objectDestructor name="massDistributionDisk_"    />
  <objectDestructor name="massDistributionSpheroid_"/>
  !!]
  call Unit_Tests_End_Group()
  
  ! Composite profile with scaled components.
  call Unit_Tests_Begin_Group("Composite profile (Hernquist + Exponential disk; scaled)")
  allocate(massDistributionExponentialDisk   :: massDistributionDisk                    )
  allocate(massDistributionHernquist         :: massDistributionSpheroid                )
  allocate(                                     massDistributions                       )
  allocate(                                     massDistributions%next                  )
  allocate(massDistributionSphericalScaler   :: massDistributions     %massDistribution_)
  allocate(massDistributionCylindricalScaler :: massDistributions%next%massDistribution_)
  select type (massDistribution_ => massDistributionSpheroid                )
  type is (massDistributionHernquist        )
     massDistribution_   =massDistributionHernquist        (                                                           dimensionless    =.true.                  )
  end select
  select type (massDistribution_ => massDistributionDisk                    )
  type is (massDistributionExponentialDisk  )
     massDistribution_   =massDistributionExponentialDisk  (scaleHeight        =0.3d-3/2.9d-3                         ,dimensionless    =.true.                  )
  end select
  select type (massDistribution_ => massDistributions     %massDistribution_)
  type is (massDistributionSphericalScaler                  )
     select type (massDistributionSpheroid)
     class is (massDistributionSpherical  )
        massDistribution_=massDistributionSphericalScaler  (factorScalingLength=0.7d-3       ,factorScalingMass=9.0d09,massDistribution_=massDistributionSpheroid)
     end select
  end select
  select type (massDistribution_ => massDistributions%next%massDistribution_)
  type is (massDistributionCylindricalScaler)
     select type (massDistributionDisk    )
     class is (massDistributionCylindrical)
        massDistribution_=massDistributionCylindricalScaler(factorScalingLength=2.9d-3       ,factorScalingMass=5.7d10,massDistribution_=massDistributionDisk    )
     end select
  end select
  allocate(massDistributionComposite :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionComposite)
     massDistribution_=massDistributionComposite(massDistributions)
  end select
  symmetry_=massDistribution_%symmetry()
  call Assert("Maximal symmetry [cylindrical]",symmetry_        %ID                        ,massDistributionSymmetryCylindrical%ID              )
  positionCartesian=[1.0d-3,0.0d0,0.0d0]
  call Assert("Density at (x,y,z)=(1,0,0) kpc",massDistribution_%density(positionCartesian),1.47756308872d18                      ,relTol=1.0d-6)
  nullify   (massDistributions)
  deallocate(massDistribution_)
  call Unit_Tests_End_Group()

  ! Spherical scaler tidal tensor memoization.
  ! `sphericalScalerTidalTensor` memoizes its result, keyed on the Cartesian components of the coordinates at
  ! which it was last evaluated. These tests pin the invariants of that cache: a repeated evaluation at one
  ! position must reproduce the cached tensor exactly; evaluating at a second position and then returning to
  ! the first must not serve a stale result; and, because the key is Cartesian, a spherical coordinate object
  ! must give the same tensor as the Cartesian object describing the same point.
  call Unit_Tests_Begin_Group("Spherical scaler tidal tensor (memoization)")
  allocate(massDistributionHernquist       :: massDistributionScalerInner)
  select type (massDistributionScalerInner)
  type is (massDistributionHernquist      )
     massDistributionScalerInner=massDistributionHernquist      (mass=9.0d09,scaleLength=0.7d-3,componentType=componentTypeSpheroid)
  end select
  allocate(massDistributionSphericalScaler :: massDistributionScalerTest )
  select type (massDistributionScalerTest)
  type is (massDistributionSphericalScaler)
     select type (massDistributionScalerInner)
     class is (massDistributionSpherical  )
        massDistributionScalerTest=massDistributionSphericalScaler(factorScalingLength=2.0d0,factorScalingMass=3.0d0,massDistribution_=massDistributionScalerInner)
     end select
  end select
  ! Evaluate at a first position, then again at the same position: the memoized result must be returned
  ! unchanged.
  position                      =[1.0d-3,Pi/3.0d0,Pi/5.0d0]
  tidalTensorComponents         =massDistributionScalerTest%tidalTensor(position         )
  tidalTensorScalerReference    =tidalTensorComponents
  tidalTensorComponents         =massDistributionScalerTest%tidalTensor(position         )
  call Assert("repeat evaluation returns cached tensor"      ,tidalTensorComponents,tidalTensorScalerReference,absTol=0.0d0                )
  ! Evaluate at a second, distinct position. This must miss the cache and give a different tensor.
  positionReference             =[4.0d-3,Pi/3.0d0,Pi/5.0d0]
  tidalTensorScalerDisplaced    =massDistributionScalerTest%tidalTensor(positionReference)
  call Assert("distinct position gives distinct tensor"      ,maxval(abs(tidalTensorScalerDisplaced-tidalTensorScalerReference)) > 0.0d0,.true.)
  ! Return to the first position: the tensor must match the original value, not the displaced one.
  tidalTensorComponents         =massDistributionScalerTest%tidalTensor(position         )
  call Assert("returning to first position is not stale"     ,tidalTensorComponents,tidalTensorScalerReference,absTol=0.0d0                )
  ! The same point expressed in Cartesian coordinates must give the same tensor.
  positionCartesian             =position
  tidalTensorComponents         =massDistributionScalerTest%tidalTensor(positionCartesian)
  call Assert("Cartesian input matches spherical input"      ,tidalTensorComponents,tidalTensorScalerReference,absTol=0.0d0                )
  ! Note that the scaler holds a counted reference to `massDistributionScalerInner` and destroys it in its own
  ! destructor, so only the scaler itself is deallocated here.
  deallocate(massDistributionScalerTest )
  call Unit_Tests_End_Group()

  ! Spherical scaler re-initialization (pool re-use).
  ! `darkMatterProfileDarkMatterOnly` re-uses pooled scalers (wrapping, for example, NFW distributions) by re-pointing them (and
  ! their attached scaled kinematics distributions) to the distribution for a new node via their "reinitialize" methods. Verify
  ! that a scaler which has been used for one distribution, and then re-initialized to a second, returns results identical to
  ! those of a scaler freshly constructed for the second. The tidal tensor is evaluated at the same position for both
  ! distributions, so that a stale memoized tensor would be detected.
  call Unit_Tests_Begin_Group("Spherical scaler (re-initialization for pool re-use)")
  allocate(massDistributionNFW       :: massDistributionScalerInnerA)
  allocate(massDistributionNFW       :: massDistributionScalerInnerB)
  allocate(kinematicsDistributionNFW :: kinematicsDistributionInnerA)
  allocate(kinematicsDistributionNFW :: kinematicsDistributionInnerB)
  select type (massDistributionScalerInnerA)
  type is (massDistributionNFW)
     massDistributionScalerInnerA=massDistributionNFW(mass=1.0d11,radiusVirial=0.1d0,scaleLength=0.01d0)
  end select
  select type (massDistributionScalerInnerB)
  type is (massDistributionNFW)
     massDistributionScalerInnerB=massDistributionNFW(mass=1.0d12,radiusVirial=0.2d0,scaleLength=0.03d0)
  end select
  select type (kinematicsDistributionInnerA)
  type is (kinematicsDistributionNFW)
     kinematicsDistributionInnerA=kinematicsDistributionNFW(useSeriesApproximation=.false.)
  end select
  select type (kinematicsDistributionInnerB)
  type is (kinematicsDistributionNFW)
     kinematicsDistributionInnerB=kinematicsDistributionNFW(useSeriesApproximation=.false.)
  end select
  call massDistributionScalerInnerA%setKinematicsDistribution(kinematicsDistributionInnerA)
  call massDistributionScalerInnerB%setKinematicsDistribution(kinematicsDistributionInnerB)
  position=[2.0d-2,Pi/3.0d0,Pi/5.0d0]
  ! Reference results from a scaler freshly constructed around distribution B.
  allocate(massDistributionSphericalScaler :: massDistributionScalerFresh)
  select type (massDistributionScalerFresh)
  type is (massDistributionSphericalScaler)
     select type (massDistributionScalerInnerB)
     class is (massDistributionSpherical)
        massDistributionScalerFresh=massDistributionSphericalScaler(factorScalingLength=2.0d0,factorScalingMass=3.0d0,massDistribution_=massDistributionScalerInnerB)
     end select
  end select
  allocate(kinematicsDistributionSphericalScaler :: kinematicsDistributionScaler_)
  select type (kinematicsDistributionScaler_)
  type is (kinematicsDistributionSphericalScaler)
     kinematicsDistributionScaler_=kinematicsDistributionSphericalScaler(factorScalingLength=2.0d0,factorScalingMass=3.0d0,kinematicsDistribution_=kinematicsDistributionInnerB)
  end select
  call massDistributionScalerFresh%setKinematicsDistribution(kinematicsDistributionScaler_)
  call scalerReinitializeEvaluate(massDistributionScalerFresh,scalerReinitializeFresh,tidalTensorScalerReference)
  ! A scaler constructed around distribution A, used, and then re-initialized to distribution B.
  allocate(massDistributionSphericalScaler :: massDistributionScalerTest)
  select type (massDistributionScalerTest)
  type is (massDistributionSphericalScaler)
     select type (massDistributionScalerInnerA)
     class is (massDistributionSpherical)
        massDistributionScalerTest=massDistributionSphericalScaler(factorScalingLength=2.0d0,factorScalingMass=3.0d0,massDistribution_=massDistributionScalerInnerA)
     end select
  end select
  allocate(kinematicsDistributionSphericalScaler :: kinematicsDistributionScaler_)
  select type (kinematicsDistributionScaler_)
  type is (kinematicsDistributionSphericalScaler)
     kinematicsDistributionScaler_=kinematicsDistributionSphericalScaler(factorScalingLength=2.0d0,factorScalingMass=3.0d0,kinematicsDistribution_=kinematicsDistributionInnerA)
  end select
  call massDistributionScalerTest%setKinematicsDistribution(kinematicsDistributionScaler_)
  call scalerReinitializeEvaluate(massDistributionScalerTest,scalerReinitializeReused,tidalTensorScalerDisplaced)
  call Assert("distribution A differs from distribution B",any(scalerReinitializeReused /= scalerReinitializeFresh),.true.)
  select type (massDistributionScalerTest)
  type is (massDistributionSphericalScaler)
     select type (massDistributionScalerInnerB)
     class is (massDistributionSpherical)
        call massDistributionScalerTest%reinitialize(massDistributionScalerInnerB)
     end select
  end select
  select type (kinematicsDistributionScaler_)
  type is (kinematicsDistributionSphericalScaler)
     call kinematicsDistributionScaler_%reinitialize(kinematicsDistributionInnerB)
  end select
  call scalerReinitializeEvaluate(massDistributionScalerTest,scalerReinitializeReused,tidalTensorComponents)
  call Assert("re-use A->B: mass enclosed"                   ,scalerReinitializeReused(1),scalerReinitializeFresh(1),absTol=0.0d0)
  call Assert("re-use A->B: density"                         ,scalerReinitializeReused(2),scalerReinitializeFresh(2),absTol=0.0d0)
  call Assert("re-use A->B: potential"                       ,scalerReinitializeReused(3),scalerReinitializeFresh(3),absTol=0.0d0)
  call Assert("re-use A->B: radius enclosing mass"           ,scalerReinitializeReused(4),scalerReinitializeFresh(4),absTol=0.0d0)
  call Assert("re-use A->B: velocity dispersion"             ,scalerReinitializeReused(5),scalerReinitializeFresh(5),absTol=0.0d0)
  call Assert("re-use A->B: tidal tensor"                    ,tidalTensorComponents      ,tidalTensorScalerReference,absTol=0.0d0)
  ! The scalers hold counted references to their inner distributions, and to their attached kinematics distributions, and
  ! destroy them in their own destructors, so only the scalers themselves are deallocated here.
  deallocate(massDistributionScalerFresh)
  deallocate(massDistributionScalerTest )
  call Unit_Tests_End_Group()

  ! Composite subset: all-match fast path & list-node free-list re-use.
  ! Exercises the optimizations in compositeSubset: (1) when every member matches the requested
  ! filter the whole composite is returned directly ("self") instead of rebuilding a duplicate, and
  ! (2) the thread-private massDistributionList free-list recycles member-list nodes. The composite
  ! is rebuilt and destroyed on each iteration so the free-list acquire/release path is exercised
  ! across object lifetimes; densities must stay correct and stable. Reference values match the
  ! "Composite profile" test above (same component parameters).
  call Unit_Tests_Begin_Group("Composite subset (fast path & node re-use)")
  do i=1,4
     allocate(                                   massDistributions                       )
     allocate(                                   massDistributions%next                  )
     allocate(massDistributionHernquist       :: massDistributions     %massDistribution_)
     allocate(massDistributionExponentialDisk :: massDistributions%next%massDistribution_)
     select type (massDistribution_ => massDistributions     %massDistribution_)
     type is (massDistributionHernquist      )
        massDistribution_=massDistributionHernquist      (mass=9.0d09,scaleLength=0.7d-3                   ,componentType=componentTypeSpheroid)
     end select
     select type (massDistribution_ => massDistributions%next%massDistribution_)
     type is (massDistributionExponentialDisk)
        massDistribution_=massDistributionExponentialDisk(mass=5.7d10,radiusScale=2.9d-3,scaleHeight=0.3d-3,componentType=componentTypeDisk    )
     end select
     allocate(massDistributionComposite :: massDistribution_)
     select type (massDistribution_)
     type is (massDistributionComposite)
        massDistribution_=massDistributionComposite(massDistributions)
     end select
     positionCartesian=[1.0d-3,0.0d0,0.0d0]
     ! All-match subset (default componentTypeAll, massTypeAll): every member matches, so the fast
     ! path returns a distribution equivalent to the full composite.
     massDistributionAll_       => massDistribution_%subset(                                   )
     ! Partial subsets must return only the matching component (i.e. NOT the whole composite).
     massDistributionDisk_      => massDistribution_%subset(componentType=componentTypeDisk    )
     massDistributionSpheroid_  => massDistribution_%subset(componentType=componentTypeSpheroid)
     ! Re-subsetting the same composite must give an identical result (free-list re-use path).
     massDistributionDiskAgain_ => massDistribution_%subset(componentType=componentTypeDisk    )
     call Assert("all-match subset = full composite"            ,loc(massDistributionAll_      ),loc(massDistribution_    ),compare=compareEquals  )
     call Assert("disk subset      ≠ full composite"            ,loc(massDistributionDisk_     ),loc(massDistribution_    ),compare=compareNotEqual)
     call Assert("spheroid subset  ≠ full composite"            ,loc(massDistributionSpheroid_ ),loc(massDistribution_    ),compare=compareNotEqual)
     call Assert("disk re-subset   = disk subset"               ,loc(massDistributionDiskAgain_),loc(massDistributionDisk_),compare=compareEquals  )
     call Assert("all-match subset density = full composite"    ,massDistributionAll_      %density(positionCartesian),1.47756308872d18,relTol=1.0d-6)
     call Assert("disk subset density      = disk component"    ,massDistributionDisk_     %density(positionCartesian),1.27347675828d18,relTol=1.0d-6)
     call Assert("spheroid subset density  = spheroid component",massDistributionSpheroid_ %density(positionCartesian),2.04086330446d17,relTol=1.0d-6)
     call Assert("re-subset disk density consistent"            ,massDistributionDiskAgain_%density(positionCartesian),1.27347675828d18,relTol=1.0d-6)
     ! Release references. The all-match result is "self" (the composite), so balance its reference
     ! count without destroying it; the partial subsets are counted references to members.
     referenceCountDiscard=massDistributionAll_%referenceCountDecrement()
     !![
     <objectDestructor name="massDistributionDisk_"     />
     <objectDestructor name="massDistributionSpheroid_" />
     <objectDestructor name="massDistributionDiskAgain_"/>
     !!]
     nullify   (massDistributionAll_)
     nullify   (massDistributions   )
     deallocate(massDistribution_   )
  end do
  call Unit_Tests_End_Group()

  ! Patej & Loeb (2015) profile.
  call Unit_Tests_Begin_Group("Patej-Loeb (2015) profile")
  allocate(massDistributionNFW :: massDistributionDMO)
  select type (massDistributionDMO)
  type is (massDistributionNFW)
     massDistributionDMO=massDistributionNFW(scaleLength=30.0d-3,radiusVirial=300.0d-3,mass=1.0d12,dimensionless=.false.)
  end select
  allocate(massDistributionPatejLoeb2015 :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionPatejLoeb2015)
     massDistribution_=massDistributionPatejLoeb2015(gamma=1.15d0,massDistribution_=massDistributionDMO,mass=1.0d11,radiusOuter=450.0d-3,radiusShock=450.0d-3)
  end select
  do i=1,4
     radius                  =450.0d-3/2.0d0**(4-i)
     position                =[radius,0.0d0,0.0d0]
     massPatejLoeb         (i)=massDistribution_%massEnclosedBySphere (radius                     )
     densityPatejLoeb      (i)=massDistribution_%density              (position                   )
     densitySlopePatejLoeb (i)=massDistribution_%densityGradientRadial(position,logarithmic=.true.)
     potentialPatejLoeb    (i)=massDistribution_%potential            (position                   )
     densityMomentPatejLoeb(i)=massDistribution_%densityRadialMoment  (-dble(i-1),450.0d-3/8.0d0,450.0d-3/1.0d0)
  end do
  call Assert(                                                                                            &
       &      "M(r) at r=[⅛,¼,½,1]rₛ"                                                                   , &
       &      massPatejLoeb                                                                             , &
       &      [+1.555565950166512d10,+3.514143268178815d10,+6.418101338664305d10,+1.0000000000000000d11], &
       &      relTol=1.0d-6                                                                               &
       &     )
  call Assert(                                                                                            &
       &      "ρ(r) at r=[⅛,¼,½,1]rₛ"                                                                   , &
       &      densityPatejLoeb                                                                          , &
       &      [+9.377714853155970d12,+1.985078164671267d12,+3.3223318757716390d11,+4.809898607847662d10], &
       &      relTol=1.0d-6                                                                               &
       &     )
  call Assert(                                                                                            &
       &      "α(r) at r=[⅛,¼,½,1]rₛ"                                                                   , &
       &      densitySlopePatejLoeb                                                                     , &
       &      [-2.030591309692208d00,-2.431529802045667d00,-2.7035845062827400d00,-2.856250000000000d00], &
       &      relTol=1.0d-6                                                                               &
       &     )
  call Assert(                                                                                            &
       &      "φ(r) at r=[⅛,¼,½,1]rₛ"                                                                   , &
       &      potentialPatejLoeb                                                                        , &
       &      [-4.079034939868978d03,-3.184928656647174d03,-2.280711782065997d03,-1.5198632412001940d03], &
       &      relTol=1.0d-3                                                                               &
       &     )
  call Assert(                                                                                            &
       &      "ℛᵢ(⅛rₛ,rₛ)"                                                                              , &
       &      densityMomentPatejLoeb                                                                    , &
       &      [+3.765562121061061d11,+4.129007833522932d12,+5.196937885662541d13,+7.1052241045332330d14], &
       &      relTol=1.0d-6                                                                               &
       &     )
  deallocate(massDistribution_)
  call Unit_Tests_End_Group()

  ! Isothermal profile.
  call Unit_Tests_Begin_Group("Isothermal profile")
  allocate(massDistributionIsothermal  :: massDistribution_      )
  allocate(kinematicsDistributionLocal :: kinematicsDistribution_)
  select type (kinematicsDistribution_)
  type is (kinematicsDistributionLocal)
     kinematicsDistribution_=kinematicsDistributionLocal(alpha=1.0d0/sqrt(2.0d0))
  end select
  select type (massDistribution_)
  type is (massDistributionIsothermal)
     massDistribution_=massDistributionIsothermal(mass=1.0d12,lengthReference=300.0d-3)
     call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
     do i=1,4
        radius                                                 =300.0d-3/2.0d0**(4-i)
        position                                               =[radius,0.0d0,0.0d0]
        positionReference                                      =[300.0d-3,0.0d0,0.0d0]
        time                                                   =         2.0d0**(i-1)
        radiiIsothermal                                     (i)=radius
        massIsothermal                                      (i)=massDistribution_%massEnclosedBySphere                      (radius                                    )
        massIsothermalNumerical                             (i)=massDistribution_%massEnclosedBySphereNumerical             (radius                                    )
        densityIsothermal                                   (i)=massDistribution_%density                                   (position                                  )
        velocityCircularIsothermal                          (i)=massDistribution_%rotationCurve                             (radius                                    )
        radiusEnclosingMassIsothermal                       (i)=massDistribution_%radiusEnclosingMass                       (      massIsothermal(i)                   )
        radiusEnclosingMassIsothermalNumerical              (i)=massDistribution_%radiusEnclosingMassNumerical              (      massIsothermal(i)                   )
        radiusEnclosingDensityIsothermal                    (i)=massDistribution_%radiusEnclosingDensity                    (3.0d0*massIsothermal(i)/4.0d0/Pi/radius**3)
        radiusEnclosingDensityIsothermalNumerical           (i)=massDistribution_%radiusEnclosingDensityNumerical           (3.0d0*massIsothermal(i)/4.0d0/Pi/radius**3)
        radiusFromSpecificAngularMomentumIsothermal         (i)=massDistribution_%radiusFromSpecificAngularMomentum         (velocityCircularIsothermal(i)*radius      )
        radiusFromSpecificAngularMomentumIsothermalNumerical(i)=massDistribution_%radiusFromSpecificAngularMomentumNumerical(velocityCircularIsothermal(i)*radius      )
        densitySlopeIsothermal                              (i)=massDistribution_%densityGradientRadial                     (position,logarithmic=.true.               )
        densitySlopeIsothermalNumerical                     (i)=massDistribution_%densityGradientRadialNumerical            (position,logarithmic=.true.               )
        potentialIsothermal                                 (i)=massDistribution_%potential                                 (position                                  )
        potentialIsothermalDifferenceNumerical              (i)=massDistribution_%potentialDifferenceNumerical              (position,positionReference                )
        densityMomentIsothermal                             (i)=massDistribution_%densityRadialMoment                       (-dble(i-1),300.0d-3/8.0d0,300.0d-3/1.0d0  )
        fourierTransformIsothermal                          (i)=massDistribution_%fourierTransform                          (300.0d-3,1.0d0/radius                     )
        fourierTransformIsothermalNumerical                 (i)=massDistribution_%fourierTransformNumerical                 (300.0d-3,1.0d0/radius                     )
        radiusFreefallIsothermal                            (i)=massDistribution_%radiusFreefall                            (time                                      )
        radiusFreefallIsothermalNumerical                   (i)=massDistribution_%radiusFreefallNumerical                   (time                                      )
        radiusFreefallGrowthRateIsothermal                  (i)=massDistribution_%radiusFreefallIncreaseRate                (time                                      )
        radiusFreefallGrowthRateIsothermalNumerical         (i)=massDistribution_%radiusFreefallIncreaseRateNumerical       (time                                      )
     end do
  end select
  call Assert(                                                                                             &
       &      "M(r) at r=[⅛,¼,½,1]rₗ"                                                                    , &
       &      massIsothermal                                                                             , &
       &      [+1.250000000000000d+11,+2.50000000000000d+11,+5.00000000000000d+11,+1.000000000000000d+12], &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "M(r) at r=[⅛,¼,½,1]rₗ numerical"                                                          , &
       &      massIsothermal                                                                             , &
       &      massIsothermalNumerical                                                                    , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "ρ(r) at r=[⅛,¼,½,1]rₗ"                                                                    , &
       &      densityIsothermal                                                                          , &
       &      [+1.886280807015056d+14,+4.71570201753764d+13,+1.17892550438441d+13,+2.947313760961025d+12], &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "V(r) at r=[⅛,¼,½,1]rₗ"                                                                    , &
       &      velocityCircularIsothermal                                                                 , &
       &      [+1.19735820315671d+02,+1.19735820315671d+02,+1.197358203156711d+02,+1.197358203156711d+02], &
       &      relTol=1.0d-4                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "r(M) at r=[⅛,¼,½,1]rₗ"                                                                    , &
       &      radiusEnclosingMassIsothermal                                                              , &
       &      radiiIsothermal                                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "r(M) at r=[⅛,¼,½,1]rₗ numerical"                                                          , &
       &      radiusEnclosingMassIsothermalNumerical                                                     , &
       &      radiiIsothermal                                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "r(ρ) at r=[⅛,¼,½,1]rₗ"                                                                    , &
       &      radiusEnclosingDensityIsothermal                                                           , &
       &      radiiIsothermal                                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "r(ρ) at r=[⅛,¼,½,1]rₗ numerical"                                                          , &
       &      radiusEnclosingDensityIsothermalNumerical                                                  , &
       &      radiiIsothermal                                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "r(j) at r=[⅛,¼,½,1]rₗ"                                                                     , &
       &      radiusFromSpecificAngularMomentumIsothermal                                                , &
       &      radiiIsothermal                                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "r(j) at r=[⅛,¼,½,1]rₗ numerical"                                                           , &
       &      radiusFromSpecificAngularMomentumIsothermalNumerical                                       , &
       &      radiiIsothermal                                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "α(r) at r=[⅛,¼,½,1]rₗ"                                                                     , &
       &      densitySlopeIsothermal                                                                     , &
       &      [-2.000000000000000d+00,-2.00000000000000d+00,-2.00000000000000d+00,-2.000000000000000d+00], &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "α(r) at r=[⅛,¼,½,1]rₗ numerical"                                                           , &
       &      densitySlopeIsothermal                                                                     , &
       &      densitySlopeIsothermalNumerical                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "φ(r) at r=[⅛,¼,½,1]rₗ"                                                                     , &
       &      potentialIsothermal                                                                        , &
       &      [-2.981226023588325d+04,-1.98748401572555d+04,-9.93742007862775d+03,+0.000000000000000d+00], &       
       &      relTol=1.0d-3                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "φ(r)-φ(rₗ) at r=[⅛,¼,½,1]rₗ numerical"                                                     , &
       &      potentialIsothermal                   -potentialIsothermal(4)                              , &
       &      potentialIsothermalDifferenceNumerical                                                     , &
       &      relTol=1.0d-3                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "ℛᵢ(⅛rₗ,rₗ)"                                                                                , &
       &      densityMomentIsothermal                                                                    , &
       &      [+6.189358898018186d+12,+9.28403834702773d+13,+1.67341925761232d+15,+3.352569403093177d+16], &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "ℱ(rₗ,k) at k=[8,4,2,1]/rₗ"                                                                 , &
       &      fourierTransformIsothermal                                                                 , &
       &      fourierTransformIsothermalNumerical                                                        , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "ℱ(rₗ,k) at k=[8,4,2,1]/rₗ numerical"                                                      , &
       &      fourierTransformIsothermalNumerical                                                        , &
       &      [+1.967733527133679d-01,+4.395507847372633d-01,+8.02706488401347d-01,+9.46083070367183d-01], &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "rᵩ(t) at t=[1,2,4,8]Gyr"                                                                  , &
       &      radiusFreefallIsothermal                                                                   , &
       &      [+9.769496930104140d-02,+1.953899386020827d-01,+3.907798772041654d-01,+7.8155975440833d-01], &
       &      relTol=1.0d-3                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "rᵩ(t) at t=[1,2,4,8]Gyr numerical"                                                        , &
       &      radiusFreefallIsothermal                                                                   , &
       &      radiusFreefallIsothermalNumerical                                                          , &
       &      relTol=1.0d-3                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "drᵩ/dt(t) at t=[1,2,4,8]Gyr numerical"                                                    , &
       &      radiusFreefallGrowthRateIsothermal                                                         , &
       &      radiusFreefallGrowthRateIsothermalNumerical                                                , &
       &      relTol=1.0d-3                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "E"                                                                                        , &
       &      massDistribution_%energy(300.0d-3,massDistribution_)                                       , &
       &      -3.58417d15                                                                                , &
       &      relTol=1.0d-3                                                                                &
       &     )
  call Unit_Tests_End_Group()

  ! Black hole mass distribution: verify that re-initialization (as used by the object pool to re-use
  ! a pooled object for a new node) fully resets the object's state, leaving no stale residue. This
  ! directly exercises the "initialize" method on which the black hole component's object pool relies.
  call Unit_Tests_Begin_Group("Black hole mass distribution (re-initialization for pool re-use)")
  if (allocated(massDistribution_)) deallocate(massDistribution_)
  allocate(massDistributionBlackHole :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionBlackHole)
     ! Reference outputs from freshly-constructed objects of each mass.
     massDistribution_ =massDistributionBlackHole             (mass=massBlackHoleA)
     enclosedReferenceA=massDistribution_%massEnclosedBySphere(radiusBlackHoleTest)
     rotationReferenceA=massDistribution_%rotationCurve       (radiusBlackHoleTest)
     massDistribution_ =massDistributionBlackHole             (mass=massBlackHoleB)
     enclosedReferenceB=massDistribution_%massEnclosedBySphere(radiusBlackHoleTest)
     rotationReferenceB=massDistribution_%rotationCurve       (radiusBlackHoleTest)
     ! The object now holds state for massBlackHoleB. Re-initialize it to massBlackHoleA (as the pool
     ! does on re-use) and confirm every output matches a freshly-constructed massBlackHoleA object.
     call massDistribution_%initialize(mass=massBlackHoleA)
     call Assert("re-use B->A: total mass"    ,massDistribution_%massTotal           (                   ),massBlackHoleA    ,relTol=1.0d-9)
     call Assert("re-use B->A: mass enclosed" ,massDistribution_%massEnclosedBySphere(radiusBlackHoleTest),enclosedReferenceA,relTol=1.0d-9)
     call Assert("re-use B->A: rotation curve",massDistribution_%rotationCurve       (radiusBlackHoleTest),rotationReferenceA,relTol=1.0d-9)
     ! Re-initialize back to massBlackHoleB and confirm no stale state survived from massBlackHoleA.
     call massDistribution_%initialize(mass=massBlackHoleB)
     call Assert("re-use A->B: total mass"    ,massDistribution_%massTotal           (                   ),massBlackHoleB    ,relTol=1.0d-9)
     call Assert("re-use A->B: mass enclosed" ,massDistribution_%massEnclosedBySphere(radiusBlackHoleTest),enclosedReferenceB,relTol=1.0d-9)
     call Assert("re-use A->B: rotation curve",massDistribution_%rotationCurve       (radiusBlackHoleTest),rotationReferenceB,relTol=1.0d-9)
  end select
  deallocate(massDistribution_)
  call Unit_Tests_End_Group()

  ! Optional `logarithmic` argument in radial density gradients.
  ! The `<optionalArgument>` directive generates a guarded local copy (`logarithmic_`) of the optional
  ! `logarithmic` dummy argument. Several classes read the dummy argument itself instead, which is
  ! invalid when the argument is not supplied by the caller. Check that omitting the argument gives
  ! the same answer as passing `.false.` explicitly (its documented default).
  call Unit_Tests_Begin_Group("Optional `logarithmic` argument in radial density gradients")
  position=[1.0d0,0.0d0,0.0d0]
  allocate(massDistributionSoliton :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionSoliton)
     massDistribution_=massDistributionSoliton(radiusCore=1.0d0,densitySolitonCentral=1.0d0,dimensionless=.true.)
  end select
  call Assert(                                                                        &
       &      "soliton"                                                             , &
       &      +massDistribution_%densityGradientRadial(position                    ), &
       &      +massDistribution_%densityGradientRadial(position,logarithmic=.false.), &
       &      relTol=1.0d-12                                                          &
       &     )
  deallocate(massDistribution_)
  allocate(massDistributionSolitonNFW :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionSolitonNFW)
     massDistribution_=massDistributionSolitonNFW(radiusScale=1.0d0,radiusCore=0.1d0,radiusSoliton=0.3d0,densitySolitonCentral=1.0d0,densityNormalizationNFW=1.0d0,dimensionless=.true.)
  end select
  ! Test both sides of the soliton-to-NFW transition.
  position=[0.2d0,0.0d0,0.0d0]
  call Assert(                                                                        &
       &      "solitonNFW (soliton regime)"                                         , &
       &      +massDistribution_%densityGradientRadial(position                    ), &
       &      +massDistribution_%densityGradientRadial(position,logarithmic=.false.), &
       &      relTol=1.0d-12                                                          &
       &     )
  position=[1.0d0,0.0d0,0.0d0]
  call Assert(                                                                        &
       &      "solitonNFW (NFW regime)"                                             , &
       &      +massDistribution_%densityGradientRadial(position                    ), &
       &      +massDistribution_%densityGradientRadial(position,logarithmic=.false.), &
       &      relTol=1.0d-12                                                          &
       &     )
  deallocate(massDistribution_)
  ! Scaled distributions forward the argument to the distribution that they wrap. Note that the
  ! scaler takes a reference to the wrapped distribution and releases it when the scaler itself is
  ! destroyed - so the wrapped distributions must not be deallocated here.
  allocate(massDistributionHernquist       :: massDistributionWrappedSpherical)
  select type (massDistributionWrappedSpherical)
  type is (massDistributionHernquist      )
     massDistributionWrappedSpherical=massDistributionHernquist      (dimensionless=.true.)
  end select
  allocate(massDistributionSphericalScaler :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionSphericalScaler)
     select type (massDistributionWrappedSpherical)
     class is (massDistributionSpherical  )
        massDistribution_=massDistributionSphericalScaler(factorScalingLength=0.7d-3,factorScalingMass=9.0d09,massDistribution_=massDistributionWrappedSpherical)
     end select
  end select
  position=[1.0d-3,0.0d0,0.0d0]
  call Assert(                                                                        &
       &      "sphericalScaler"                                                     , &
       &      +massDistribution_%densityGradientRadial(position                    ), &
       &      +massDistribution_%densityGradientRadial(position,logarithmic=.false.), &
       &      relTol=1.0d-12                                                          &
       &     )
  deallocate(massDistribution_)
  allocate(massDistributionExponentialDisk   :: massDistributionWrappedCylindrical)
  select type (massDistributionWrappedCylindrical)
  type is (massDistributionExponentialDisk  )
     massDistributionWrappedCylindrical=massDistributionExponentialDisk  (scaleHeight=0.3d-3/2.9d-3,dimensionless=.true.)
  end select
  allocate(massDistributionCylindricalScaler :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionCylindricalScaler)
     select type (massDistributionWrappedCylindrical)
     class is (massDistributionCylindrical)
        massDistribution_=massDistributionCylindricalScaler(factorScalingLength=2.9d-3,factorScalingMass=5.7d10,massDistribution_=massDistributionWrappedCylindrical)
     end select
  end select
  positionCartesian=[1.0d-3,0.0d0,0.0d0]
  call Assert(                                                                                 &
       &      "cylindricalScaler"                                                            , &
       &      +massDistribution_%densityGradientRadial(positionCartesian                    ), &
       &      +massDistribution_%densityGradientRadial(positionCartesian,logarithmic=.false.), &
       &      relTol=1.0d-12                                                                   &
       &     )
  deallocate(massDistribution_)
  call Unit_Tests_End_Group()

  ! The potential tabulation rounds its bounds outward to whole octaves of an absolute lattice. A radius which already lies on
  ! one is the case where those bounds can coincide, which is rejected as a range of fewer than two points - so evaluate there.
  call Unit_Tests_Begin_Group("Potential at a radius on a whole octave")
  block
    type            (massDistributionHernquist) :: hernquist_
    type            (coordinateSpherical      ) :: coordinatesOctave
    double precision                            :: potentialOctave
    hernquist_       =massDistributionHernquist(mass=1.0d0,scaleLength=1.0d0)
    coordinatesOctave=[1.0d0,0.0d0,0.0d0]
    potentialOctave  =hernquist_%potential(coordinatesOctave)
    call Assert("the potential is finite at a radius lying on a whole octave",potentialOctave < 0.0d0,.true.)
  end block
  call Unit_Tests_End_Group()

  ! Test that the tabulations which are pinned to absolute lattices behave as intended. These properties are what pinning is
  ! for, and none of the tests above constrain them: they exercise the profiles at fixed radii, never driving a tabulation to
  ! extend and so never revealing whether extension moves values already computed.
  call Unit_Tests_Begin_Group("Pinned tabulations")
  block
    type            (massDistributionExponentialDisk) :: diskExtended  , diskDirect
    type            (massDistributionSersic         ) :: sersicA       , sersicB
    type            (coordinateSpherical            ) :: coordinates_
    double precision                                  :: velocityProbe , velocityAfterExtension, &
         &                                               velocityDirect, densityA              , &
         &                                               densityB      , velocityFar
    ! The exponential disk rotation curve. The radii are chosen to avoid the three paths which bypass the tabulation: the
    ! small-argument series expansion below a half-radius of 1e-3, the cached factor used at exactly one scale radius, and the
    ! point-mass approximation beyond `radiusMaximum`=30 scale radii. r=4 gives a half-radius of 2, inside the seed range;
    ! r=28 gives 14, outside it, and so forces a genuine extension.
    diskExtended         =massDistributionExponentialDisk(radiusScale=1.0d0,mass=1.0d0)
    velocityProbe        =diskExtended%rotationCurve( 4.0d0)
    velocityFar          =diskExtended%rotationCurve(28.0d0)
    velocityAfterExtension=diskExtended%rotationCurve( 4.0d0)
    call Assert("exponential disk: extension preserves a tabulated rotation curve",velocityAfterExtension == velocityProbe,.true.)
    ! A disk driven straight to the wider range must agree bit-for-bit with one which reached it by extension.
    diskDirect           =massDistributionExponentialDisk(radiusScale=1.0d0,mass=1.0d0)
    velocityDirect       =diskDirect  %rotationCurve(28.0d0)
    call Assert("exponential disk: rotation curve is independent of the order of requests",velocityDirect == velocityFar,.true.)
    call Assert("exponential disk: probe value is independent of the order of requests",diskDirect%rotationCurve(4.0d0) == velocityProbe,.true.)
    ! The Sersic profile. Nothing computed for it can be carried over when its range moves - the densities are Abel integrals
    ! taken out to a bound-derived radius, and the enclosed masses accumulate from the innermost point - so what is tested is
    ! reproducibility: two profiles driven to nearby radii must build identical tables. The radii are large enough to push the
    ! bounds beyond the seed range, where under the previous scheme they would have differed by the ratio of the two radii.
    sersicA     =massDistributionSersic(index=4.0d0,radiusHalfMass=1.0d0,mass=1.0d0)
    sersicB     =massDistributionSersic(index=4.0d0,radiusHalfMass=1.0d0,mass=1.0d0)
    coordinates_=[900.0d0,0.0d0,0.0d0]
    densityA    =sersicA%density(coordinates_)
    coordinates_=[901.0d0,0.0d0,0.0d0]
    densityB    =sersicB%density(coordinates_)
    coordinates_=[  1.0d0,0.0d0,0.0d0]
    densityA    =sersicA%density(coordinates_)
    densityB    =sersicB%density(coordinates_)
    call Assert("Sersic: tabulation is reproducible across nearby requested radii",densityA == densityB,.true.)
    call Assert("Sersic: enclosed mass is reproducible across nearby requested radii",sersicA%massEnclosedBySphere(1.0d0) == sersicB%massEnclosedBySphere(1.0d0),.true.)
  end block
  call Unit_Tests_End_Group()

  ! Test that enclosed masses interpolated by the spherical accelerator agree with the exact masses to within (about) its
  ! tolerance, a relative tolerance on the mass. The cache is populated on a coarse grid of radii, and then queried between those
  ! radii. The accelerator accepts an interpolated mass if two interpolations agree to within the tolerance, which estimates, but
  ! does not bound, the error of each, so the actual errors are required to be within twice the tolerance.
  call Unit_Tests_Begin_Group("Spherical accelerator")
  block
    class           (massDistributionSpherical), pointer :: nfwAccelerated_, accelerator_
    double precision                                     :: radius_        , errorMaximum
    integer                                              :: i_
    character       (len=12                   )          :: labelError
    allocate(massDistributionNFW                  :: nfwAccelerated_)
    allocate(massDistributionSphericalAccelerator :: accelerator_   )
    select type (nfwAccelerated_)
    type is (massDistributionNFW)
       !![
       <referenceConstruct object="nfwAccelerated_" constructor="massDistributionNFW(scaleLength=30.0d-3,radiusVirial=300.0d-3,mass=1.0d12,dimensionless=.false.)"/>
       !!]
    end select
    select type (accelerator_)
    type is (massDistributionSphericalAccelerator)
       !![
       <referenceConstruct object="accelerator_" constructor="massDistributionSphericalAccelerator(toleranceRelative=1.0d-2,factorRadiusMaximum=20.0d0,nonAnalyticSolver=nonAnalyticSolversNumerical,massDistribution_=nfwAccelerated_)"/>
       !!]
    end select
    do i_=0,8
       radius_     =1.0d-3*3.0d0**i_
       errorMaximum=accelerator_%massEnclosedBySphere(radius_)
    end do
    errorMaximum=0.0d0
    do i_=0,399
       radius_     =1.0d-3*10.0d0**((dble(i_)+0.37d0)/200.0d0)
       errorMaximum=max(errorMaximum,abs(accelerator_%massEnclosedBySphere(radius_)/nfwAccelerated_%massEnclosedBySphere(radius_)-1.0d0))
    end do
    write (labelError,'(e12.6)') errorMaximum
    call displayMessage('maximum fractional error in accelerated NFW enclosed mass: '//trim(labelError))
    call Assert("NFW enclosed mass within tolerance",errorMaximum,2.0d-2,compare=compareLessThan)
    !![
    <objectDestructor name="accelerator_"   />
    <objectDestructor name="nfwAccelerated_"/>
    !!]
  end block
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

contains

  subroutine scalerReinitializeEvaluate(massDistributionScaler,results,tidalTensor)
    !!{RST
    Evaluate the quantities compared in the spherical scaler re-initialization test.
    !!}
    implicit none
    class           (massDistributionClass      ), intent(inout), target        :: massDistributionScaler
    double precision                             , intent(  out), dimension(5  ) :: results
    double precision                             , intent(  out), dimension(3,3) :: tidalTensor
    class           (kinematicsDistributionClass), pointer                       :: kinematicsDistributionScaler

    kinematicsDistributionScaler => massDistributionScaler%kinematicsDistribution()
    results(1) =massDistributionScaler      %massEnclosedBySphere         (radius=2.0d-2                                         )
    results(2) =massDistributionScaler      %density                      (coordinates=position                                  )
    results(3) =massDistributionScaler      %potential                    (coordinates=position                                  )
    results(4) =massDistributionScaler      %radiusEnclosingMass          (mass=1.0d11                                           )
    results(5) =kinematicsDistributionScaler%velocityDispersion1D         (position,massDistributionScaler,massDistributionScaler)
    tidalTensor=massDistributionScaler      %tidalTensor                  (position                                              )
    !![
    <objectDestructor name="kinematicsDistributionScaler"/>
    !!]
    return
  end subroutine scalerReinitializeEvaluate

end program Test_Mass_Distributions
