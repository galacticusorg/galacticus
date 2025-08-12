!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a program to test mass distributions.
!!}

program Test_Mass_Distributions
  !!{
  Tests mass distributions.
  !!}
  use :: Coordinates               , only : assignment(=)                    , coordinateSpherical                , coordinateCartesian
  use :: Display                   , only : displayVerbositySet              , verbosityLevelStandard
  use :: Error_Functions           , only : Error_Function
  use :: Error                     , only : Error_Report
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Galactic_Structure_Options, only : componentTypeSpheroid            , componentTypeDisk
  use :: IO_HDF5                   , only : ioHDF5AccessInitialize
  use :: Linear_Algebra            , only : assignment(=)                    , vector
  use :: Mass_Distributions        , only : massDistributionBetaProfile      , massDistributionClass              , massDistributionExponentialDisk        , massDistributionGaussianEllipsoid   , &
       &                                    massDistributionHernquist        , massDistributionSersic             , massDistributionSpherical              , massDistributionComposite           , &
       &                                    massDistributionList             , massDistributionSymmetryCylindrical, enumerationMassDistributionSymmetryType, massDistributionSphericalScaler     , &
       &                                    massDistributionCylindricalScaler, massDistributionCylindrical        , massDistributionPatejLoeb2015          , massDistributionNFW                 , &
       &                                    massDistributionIsothermal       , kinematicsDistributionClass        , kinematicsDistributionLocal
  use :: Numerical_Constants_Math  , only : Pi
  use :: Tensors                   , only : assignment(=)
  use :: Unit_Tests                , only : Assert                           , Unit_Tests_Begin_Group             , Unit_Tests_End_Group                   , Unit_Tests_Finish
  implicit none
  class           (massDistributionClass                  )                             , allocatable :: massDistribution_                                                                                               , massDistributionRotated                   , &
       &                                                                                                 massDistributionDisk                                                                                            , massDistributionSpheroid                  , &
       &                                                                                                 massDistributionDMO
  class           (massDistributionClass                  )                             , pointer     :: massDistributionDisk_                                                                                           , massDistributionSpheroid_ 
  class           (kinematicsDistributionClass            )                             , allocatable :: kinematicsDistribution_
  type            (massDistributionList                   )                             , pointer     :: massDistributions
  integer                                                  , parameter                                :: sersicTableCount             =8
  double precision                                         , dimension(sersicTableCount)              :: sersicTableRadius            =[1.0000d-06,1.0000d-5,1.0000d-4,1.0000d-3,1.0000d-2,1.0000d-1,1.0000d+0,1.0000d+1]
  ! Mass targets for Sersic profile from Mazure & Capelato (2001).
  double precision                                         , dimension(sersicTableCount)              :: sersicTableMassTarget        =[1.4730d-11,2.1130d-9,2.5959d-7,2.4545d-5,1.4961d-3,4.4102d-2,4.1536d-1,9.4308d-1]
  ! Density targets for Sersic profile from Mazure & Capelato (2001).
  double precision                                         , dimension(sersicTableCount)              :: sersicTableDensityTarget     =[2.5553d+06,3.5797d+5,4.2189d+4,3.7044d+3,1.9679d+2,4.4047d+0,2.1943d-2,7.8166d-6]
  ! Potential targets for Sersic profile from Young (1976).
  double precision                                         , dimension(sersicTableCount)              :: sersicTablePotentialTarget   =[1.0000d+00,9.9993d-1,9.9908d-1,9.9027d-1,9.2671d-1,6.7129d-1,2.4945d-1,3.7383d-2]
  double precision                                         , dimension(sersicTableCount)              :: sersicTableDensity                                                                                              , sersicTableMass                            , &
       &                                                                                                 sersicTablePotential
  type            (coordinateSpherical                    )                                           :: position                                                                                                        , positionZero                               , &
       &                                                                                                 positionReference
  type            (coordinateCartesian                    )                                           :: positionCartesian
  type            (enumerationMassDistributionSymmetryType)                                           :: symmetry_
  integer                                                                                             :: i
  double precision                                                                                    :: radiusInProjection                                                                                              , radius                                     , &
       &                                                                                                 rotationCurveGradientAnalytic                                                                                   , rotationCurveGradientNumerical             , &
       &                                                                                                 massFraction                                                                                                    , time
  double precision                                         , parameter                                :: epsilonFiniteDifference      =0.01d0
  character       (len=4                                  )                                           :: label
  double precision                                         , dimension(3,3)                           :: tidalTensorComponents                                                                                           , tidalTensorSphericalComponents
  double precision                                         , dimension(3  )                           :: acceleration
  double precision                                         , dimension(4  )                           :: massPatejLoeb                                                                                                   , densityPatejLoeb                           , &
       &                                                                                                 densitySlopePatejLoeb                                                                                           , densityMomentPatejLoeb                     , &
       &                                                                                                 potentialPatejLoeb
  double precision                                         , dimension(4  )                           :: massIsothermal                                                                                                  , densityIsothermal                          , &
       &                                                                                                 densitySlopeIsothermal                                                                                          , densityMomentIsothermal                    , &
       &                                                                                                 potentialIsothermal                                                                                             , fourierTransformIsothermal                 , &
       &                                                                                                 radiusFreeFallIsothermal                                                                                        , radiusFreeFallGrowthRateIsothermal         , &
       &                                                                                                 massIsothermalNumerical                                                                                         , potentialIsothermalDifferenceNumerical     , &
       &                                                                                                 fourierTransformIsothermalNumerical                                                                             , radiusFreeFallIsothermalNumerical          , &
       &                                                                                                 radiusFreeFallGrowthRateIsothermalNumerical                                                                     , densitySlopeIsothermalNumerical            , &
       &                                                                                                 radiusEnclosingMassIsothermal                                                                                   , radiusEnclosingMassIsothermalNumerical     , &
       &                                                                                                 radiusEnclosingDensityIsothermal                                                                                , radiusEnclosingDensityIsothermalNumerical  , &
       &                                                                                                 radiiIsothermal                                                                                                 , radiusFromSpecificAngularMomentumIsothermal, &
       &                                                                                                 radiusFromSpecificAngularMomentumIsothermalNumerical                                                            , velocityCircularIsothermal
  type            (vector                                 ), dimension(:  )             , allocatable :: axes
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Initialize HDF5 lock.
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Mass distributions")

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
     call Assert("Half mass radius          (dimensionless)",massDistribution_%radiusHalfMass      (        ), 1.0d0+sqrt(2.0d0)        ,absTol=1.0d-6)
     call Assert("Mass within scale radius  (dimensionless)",massDistribution_%massEnclosedBySphere(1.0d0   ), 0.25d0                   ,absTol=1.0d-6)
     call Assert("Potential at scale radius (dimensionless)",massDistribution_%potential           (position),-0.50d0                   ,absTol=1.0d-6)
  end select
  deallocate(massDistribution_)
  allocate(massDistributionHernquist :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionHernquist)
     massDistribution_=massDistributionHernquist(dimensionless=.false.,mass=2.0d0,scaleLength=2.0d0)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     position=[2.0d0,0.0d0,0.0d0]
     call Assert("Half mass radius          (dimensionful )",massDistribution_%radiusHalfMass      (        ),2.0d0*( 1.0d0+sqrt(2.0d0)),absTol=1.0d-6)
     call Assert("Mass within scale radius  (dimensionful )",massDistribution_%massEnclosedBySphere(2.0d0   ),2.0d0*( 0.25d0           ),absTol=1.0d-6)
  end select
  deallocate(massDistribution_)
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
     massDistribution_=massDistributionBetaProfile(beta=2.0d0/3.0d0,coreRadius=3.8492316686261747d-002,mass=182582297.19568533d0,outerRadius=0.12829569846196026d0)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     call Assert("Mass within outer radius",massDistribution_%massEnclosedBySphere(0.12829569846196026d0),182582297.19568533d0,relTol=1.0d-6)
  end select
  deallocate(massDistribution_)
  call Unit_Tests_End_Group()

  ! Exponential disk profile.
  call Unit_Tests_Begin_Group("Exponential disk profile")
  allocate(massDistributionExponentialDisk :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionExponentialDisk)
     massDistribution_=massDistributionExponentialDisk(scaleHeight=0.01d0,dimensionless=.true.)
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
        call Assert("[r,θ,φ] = ["//trim(label)//",0  , 0  ]",massDistribution_%acceleration(position),[0.0d0,0.0d0,-1.0d0/radius**2],absTol=1.0d-6,relTol=1.0d-3)
        ! Vertically below the disk.
        position=[radius,Pi,0.0d0]
        call Assert("[r,θ,φ] = ["//trim(label)//",π  , 0  ]",massDistribution_%acceleration(position),[0.0d0,0.0d0,+1.0d0/radius**2],absTol=1.0d-6,relTol=1.0d-3)
        ! Disk plane, along +x-axis.
        position=[radius,Pi/2.0d0,0.0d0]
        call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",massDistribution_%acceleration(position),[-1.0d0/radius**2,0.0d0,0.0d0],absTol=1.0d-6,relTol=1.0d-2)
        ! Disk plane, along -y-axis.
        position=[radius,Pi/2.0d0,-Pi/2.0d0]
        call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",massDistribution_%acceleration(position),[0.0d0,+1.0d0/radius**2,0.0d0],absTol=1.0d-6,relTol=1.0d-2)
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
     ! Test that gravitational acceleration matches expectations for a thin disk in the midplane.
     call Unit_Tests_Begin_Group("Acceleration approaches razor-thin disk in midplane")
     radius  =0.1d0
     write (label,'(f4.1)') radius
     ! Disk plane, along +x-axis.
     position=[radius,Pi/2.0d0,0.0d0]
     call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",massDistribution_%acceleration(position),[-1.0d0,0.0d0,0.0d0]*massDistribution_%rotationCurve(radius)**2/radius,absTol=1.0d-6,relTol=1.0d-1)
     ! Disk plane, along -y-axis.
     position=[radius,Pi/2.0d0,-Pi/2.0d0]
     call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",massDistribution_%acceleration(position),[0.0d0,+1.0d0,0.0d0]*massDistribution_%rotationCurve(radius)**2/radius,absTol=1.0d-6,relTol=1.0d-1)
     call Unit_Tests_End_Group()
  class default
     call Error_Report('unknown mass distribution'//{introspection:location})
  end select
  deallocate(massDistribution_)
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
           call Assert("[r,θ,φ] = ["//trim(label)//",0  , 0  ]",massDistribution_%acceleration(position),[0.0d0,0.0d0,-1.0d0/radius**2],absTol=1.0d-6,relTol=3.0d-2)
           ! Along the -z axis.
           position=[radius,Pi   ,0.0d0]
           call Assert("[r,θ,φ] = ["//trim(label)//",π  , 0  ]",massDistribution_%acceleration(position),[0.0d0,0.0d0,+1.0d0/radius**2],absTol=1.0d-6,relTol=3.0d-2)
           ! x-y plane, along +x-axis.
           position                      =[radius,Pi/2.0d0,0.0d0]
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",massDistribution_%acceleration(position),[-1.0d0/radius**2,0.0d0,0.0d0],absTol=1.0d-6,relTol=3.0d-2)
           ! x-y plane, along -y-axis.
           position                      =[radius,Pi/2.0d0,-Pi/2.0d0]
           call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",massDistribution_%acceleration(position),[0.0d0,+1.0d0/radius**2,0.0d0],absTol=1.0d-6,relTol=3.0d-2)
        end do
        call Unit_Tests_End_Group()
        ! Test that acceleration due to rotated ellipsoid matches the rotated acceleration of an unrotated ellipsoid.
        call Unit_Tests_Begin_Group("Acceleration of rotated ellipsoid")
        position    =[1.0d0,Pi/2.0d0,3.0d0*Pi/4.0d0]
        acceleration=massDistributionRotated%acceleration(position)     ! Evaluate the acceleration in the rotated distribution at the rotated position.
        acceleration=[acceleration(2),-acceleration(1),acceleration(3)] ! De-rotate the acceleration due to the rotated distribution.
        position   =[1.0d0,Pi/2.0d0,Pi/4.0d0]                           ! Evaluate the acceleration in the un-rotated distribution at the un-rotated position.
        call Assert("[r,θ,φ] = [1,π/2,π/4]",massDistribution_%acceleration(position),acceleration,absTol=1.0d-6,relTol=1.0d-6)
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
        call Assert("[r,θ,φ] = ["//trim(label)//",0  , 0  ]",massDistribution_%acceleration(position),[0.0d0,0.0d0,-massFraction/radius**2],absTol=1.0d-6,relTol=2.0d-2)
        ! Along the -z axis.
        position=[radius,Pi   ,0.0d0]
        call Assert("[r,θ,φ] = ["//trim(label)//",π  , 0  ]",massDistribution_%acceleration(position),[0.0d0,0.0d0,+massFraction/radius**2],absTol=1.0d-6,relTol=2.0d-2)
        ! x-y plane, along +x-axis.
        position                      =[radius,Pi/2.0d0,0.0d0]
        call Assert("[r,θ,φ] = ["//trim(label)//",π/2, 0  ]",massDistribution_%acceleration(position),[-massFraction/radius**2,0.0d0,0.0d0],absTol=1.0d-6,relTol=2.0d-2)
        ! x-y plane, along -y-axis.
        position                      =[radius,Pi/2.0d0,-Pi/2.0d0]
        call Assert("[r,θ,φ] = ["//trim(label)//",π/2,-π/2]",massDistribution_%acceleration(position),[0.0d0,+massFraction/radius**2,0.0d0],absTol=1.0d-6,relTol=2.0d-2)
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
     massDistribution_=massDistributionExponentialDisk(mass=5.7d10,scaleRadius=2.9d-3,scaleHeight=0.3d-3,componentType=componentTypeDisk    )
  end select
  allocate(massDistributionComposite :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionComposite)
     massDistribution_=massDistributionComposite(massDistributions)
  end select
  symmetry_=massDistribution_%symmetry()
  call Assert("Maximal symmetry [cylindrical]"         ,symmetry_        %ID                                ,massDistributionSymmetryCylindrical%ID              )
  massDistributionDisk_     => massDistribution_%subset(componentType=componentTypeDisk    )
  massDistributionSpheroid_ => massDistribution_%subset(componentType=componentTypeSpheroid)
  positionCartesian         =  [1.0d-3,0.0d0,0.0d0]
  call Assert("Density at (x,y,z)=(1,0,0) kpc"         ,massDistribution_        %density(positionCartesian),1.47756308872d18                      ,relTol=1.0d-6)
  call Assert("Spheroid density at (x,y,z)=(1,0,0) kpc",massDistributionSpheroid_%density(positionCartesian),2.04086330446d17                      ,relTol=1.0d-6)
  call Assert("Disk density at (x,y,z)=(1,0,0) kpc"    ,massDistributionDisk_    %density(positionCartesian),1.27347675828d18                      ,relTol=1.0d-6)
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

  ! Patej & Loeb (2015) profile.
  call Unit_Tests_Begin_Group("Patej-Loeb (2015) profile")
  allocate(massDistributionNFW :: massDistributionDMO)
  select type (massDistributionDMO)
  type is (massDistributionNFW)
     massDistributionDMO=massDistributionNFW(scaleLength=30.0d-3,virialRadius=300.0d-3,mass=1.0d12,dimensionless=.false.)
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
       &      "r(j) at r=[⅛,¼,½,1]rₗ"                                                                    , &
       &      radiusFromSpecificAngularMomentumIsothermal                                                , &
       &      radiiIsothermal                                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "r(j) at r=[⅛,¼,½,1]rₗ numerical"                                                          , &
       &      radiusFromSpecificAngularMomentumIsothermalNumerical                                       , &
       &      radiiIsothermal                                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "α(r) at r=[⅛,¼,½,1]rₗ"                                                                    , &
       &      densitySlopeIsothermal                                                                     , &
       &      [-2.000000000000000d+00,-2.00000000000000d+00,-2.00000000000000d+00,-2.000000000000000d+00], &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "α(r) at r=[⅛,¼,½,1]rₗ numerical"                                                          , &
       &      densitySlopeIsothermal                                                                     , &
       &      densitySlopeIsothermalNumerical                                                            , &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "φ(r) at r=[⅛,¼,½,1]rₗ"                                                                    , &
       &      potentialIsothermal                                                                        , &
       &      [-2.981226023588325d+04,-1.98748401572555d+04,-9.93742007862775d+03,+0.000000000000000d+00], &       
       &      relTol=1.0d-3                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "φ(r)-φ(rₗ) at r=[⅛,¼,½,1]rₗ numerical"                                                    , &
       &      potentialIsothermal                   -potentialIsothermal(4)                              , &
       &      potentialIsothermalDifferenceNumerical                                                     , &
       &      relTol=1.0d-3                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "ℛᵢ(⅛rₗ,rₗ)"                                                                               , &
       &      densityMomentIsothermal                                                                    , &
       &      [+6.189358898018186d+12,+9.28403834702773d+13,+1.67341925761232d+15,+3.352569403093177d+16], &
       &      relTol=1.0d-6                                                                                &
       &     )
  call Assert(                                                                                             &
       &      "ℱ(rₗ,k) at k=[8,4,2,1]/rₗ"                                                                , &
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

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Test_Mass_Distributions
