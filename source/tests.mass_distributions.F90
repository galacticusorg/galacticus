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

!% Contains a program to test mass distributions.

program Test_Mass_Distributions
  !% Tests mass distributions.
  use Galacticus_Error
  use Unit_Tests
  use Memory_Management
  use Mass_Distributions
  use Coordinates
  use Numerical_Constants_Math
  use Galacticus_Display
  implicit none
  class           (massDistributionClass), allocatable                 :: massDistribution_
  integer                                , parameter                   :: sersicTableCount          =8
  double precision                       , dimension(sersicTableCount) :: sersicTableRadius         =[1.0000d-06,1.0000d-5,1.0000d-4,1.0000d-3,1.0000d-2,1.0000d-1,1.0000d+0,1.0000d+1]
  ! Mass targets for Sersic profile from Mazure & Capelato (2001).
  double precision                       , dimension(sersicTableCount) :: sersicTableMassTarget     =[1.4730d-11,2.1130d-9,2.5959d-7,2.4545d-5,1.4961d-3,4.4102d-2,4.1536d-1,9.4308d-1]
  ! Density targets for Sersic profile from Mazure & Capelato (2001).
  double precision                       , dimension(sersicTableCount) :: sersicTableDensityTarget  =[2.5553d+06,3.5797d+5,4.2189d+4,3.7044d+3,1.9679d+2,4.4047d+0,2.1943d-2,7.8166d-6]
  ! Potential targets for Sersic profile from Young (1976).
  double precision                       , dimension(sersicTableCount) :: sersicTablePotentialTarget=[1.0000d+00,9.9993d-1,9.9908d-1,9.9027d-1,9.2671d-1,6.7129d-1,2.4945d-1,3.7383d-2]
  double precision                       , dimension(sersicTableCount) :: sersicTableDensity                                                                                           , sersicTableMass, &
       &                                                                  sersicTablePotential
  type            (coordinateSpherical  )                              :: position                                                                                                     , positionZero
  integer                                                              :: i
  double precision                                                     :: radiusInProjection

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.mass_distributions.size')

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
     call Galacticus_Error_Report('unknown mass distribution'//{introspection:location})
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

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Mass_Distributions
