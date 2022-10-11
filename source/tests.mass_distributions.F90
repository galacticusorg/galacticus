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
Contains a program to test mass distributions.
!!}

program Test_Mass_Distributions
  !!{
  Tests mass distributions.
  !!}
  use :: Coordinates             , only : assignment(=)              , coordinateSpherical
  use :: Display                 , only : displayVerbositySet        , verbosityLevelStandard
  use :: Error_Functions         , only : Error_Function
  use :: Error                   , only : Error_Report
  use :: Events_Hooks            , only : eventsHooksInitialize
  use :: Linear_Algebra          , only : assignment(=)              , vector
  use :: Mass_Distributions      , only : massDistributionBetaProfile, massDistributionClass , massDistributionExponentialDisk, massDistributionGaussianEllipsoid, &
          &                               massDistributionHernquist  , massDistributionSersic, massDistributionSpherical
  use :: Numerical_Constants_Math, only : Pi
  use :: Tensors                 , only : assignment(=)
  use :: Unit_Tests              , only : Assert                     , Unit_Tests_Begin_Group, Unit_Tests_End_Group           , Unit_Tests_Finish
  implicit none
  class           (massDistributionClass)                             , allocatable :: massDistribution_                                                                                               , massDistributionRotated
  integer                                , parameter                                :: sersicTableCount             =8
  double precision                       , dimension(sersicTableCount)              :: sersicTableRadius            =[1.0000d-06,1.0000d-5,1.0000d-4,1.0000d-3,1.0000d-2,1.0000d-1,1.0000d+0,1.0000d+1]
  ! Mass targets for Sersic profile from Mazure & Capelato (2001).
  double precision                       , dimension(sersicTableCount)              :: sersicTableMassTarget        =[1.4730d-11,2.1130d-9,2.5959d-7,2.4545d-5,1.4961d-3,4.4102d-2,4.1536d-1,9.4308d-1]
  ! Density targets for Sersic profile from Mazure & Capelato (2001).
  double precision                       , dimension(sersicTableCount)              :: sersicTableDensityTarget     =[2.5553d+06,3.5797d+5,4.2189d+4,3.7044d+3,1.9679d+2,4.4047d+0,2.1943d-2,7.8166d-6]
  ! Potential targets for Sersic profile from Young (1976).
  double precision                       , dimension(sersicTableCount)              :: sersicTablePotentialTarget   =[1.0000d+00,9.9993d-1,9.9908d-1,9.9027d-1,9.2671d-1,6.7129d-1,2.4945d-1,3.7383d-2]
  double precision                       , dimension(sersicTableCount)              :: sersicTableDensity                                                                                              , sersicTableMass               , &
       &                                                                               sersicTablePotential
  type            (coordinateSpherical  )                                           :: position                                                                                                        , positionZero
  integer                                                                           :: i
  double precision                                                                  :: radiusInProjection                                                                                              , radius                        , &
       &                                                                               rotationCurveGradientAnalytic                                                                                   , rotationCurveGradientNumerical, &
       &                                                                               massFraction
  double precision                       , parameter                                :: epsilonFiniteDifference      =0.01d0
  character       (len=4                )                                           :: label
  double precision                       , dimension(3,3)                           :: tidalTensorComponents                                                                                           , tidalTensorSphericalComponents
  double precision                       , dimension(3  )                           :: acceleration
  type            (vector               ), dimension(:  )             , allocatable :: axes
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
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
     call Assert("Half mass radius"         ,massDistribution_%radiusHalfMass      (        ), 1.0d0+sqrt(2.0d0),absTol=1.0d-6)
     call Assert("Mass within scale radius" ,massDistribution_%massEnclosedBySphere(1.0d0   ), 0.25d0           ,absTol=1.0d-6)
     call Assert("Potential at scale radius",massDistribution_%potential           (position),-0.50d0           ,absTol=1.0d-6)
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
  allocate(massDistributionBetaProfile :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionBetaProfile)
     massDistribution_=massDistributionBetaProfile(beta=2.0d0/3.0d0,dimensionless=.true.)
  end select
  select type (massDistribution_)
  class is (massDistributionSpherical)
     position    =[1.0d0,0.0d0,0.0d0]
     positionZero=[0.0d0,0.0d0,0.0d0]
     call Assert(                                                        &
          &      "Mass within scale radius"                            , &
          &      +massDistribution_%massEnclosedBySphere(1.0d0       ), &
          &      +(4.0d0-Pi)*Pi                                        , &
          &      absTol=1.0d-6                                           &
          &     )
     call Assert(                                                        &
          &      "Potential at scale radius"                           , &
          &      +massDistribution_%potential           (position    )  &
          &      -massDistribution_%potential           (positionZero), &
          &      +Pi*(Pi+log(4.0d0)-4.0d0)                             , &
          &      absTol=1.0d-6                                           &
          &     )
  end select
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
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Mass_Distributions
