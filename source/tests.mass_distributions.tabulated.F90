!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
Contains a program to test tabulated mass distributions.
!!}

program Test_Mass_Distributions_Tabulated
  !!{
  Tests mass distributions.
  !!}
  use :: Coordinates                     , only : coordinateSpherical    , assignment(=)
  use :: Display                         , only : displayVerbositySet    , verbosityLevelWorking
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Mass_Distributions              , only : massDistributionClass  , massDistributionSpherical      , massDistributionCoredNFW   , kinematicsDistributionCollisionlessTabulated, &
       &                                          massDistributionCuspNFW, kinematicsDistributionCuspNFW  , kinematicsDistributionClass
  use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr      , gravitationalConstant_internal
  use :: Numerical_Constants_Math        , only : Pi
  use :: Unit_Tests                      , only : Assert                 , Unit_Tests_Begin_Group         , Unit_Tests_End_Group       , Unit_Tests_Finish                           , &
       &                                          Skip
  use :: IO_HDF5                         , only : ioHDF5AccessInitialize
  implicit none
  class           (massDistributionClass      ), allocatable             :: massDistribution_
  class           (kinematicsDistributionClass), pointer                 :: kinematicsDistribution_
  double precision                             , parameter               :: radiusVirial     =0.25d0, radiusScale               =0.025d0, &
       &                                                                    massVirial       =1.0d12, radiusCore                =0.010d0, &
       &                                                                    yCusp            =0.2d00
  double precision                                        , dimension(7) :: mass                    , densityMoment0                    , &
       &                                                                    densityMoment1          , densityMoment2                    , &
       &                                                                    densityMoment3          , energy                            , &
       &                                                                    radiusFreefall          , radiusFreefallIncreaseRate        , &
       &                                                                    potential               , fourierTransform                  , &
       &                                                                    velocityDispersion      , massTarget                        , &
       &                                                                    densityMoment0Target    , densityMoment1Target              , &
       &                                                                    densityMoment2Target    , densityMoment3Target
  type            (coordinateSpherical        )                          :: coordinates             , coordinatesReference
  double precision                                                       :: timeScale
  integer                                                                :: i                       , iProfile
  double precision                             , parameter, dimension(7) :: radiiScaleFree      =[                                           &
       &                                                                                          0.010000000000000d00,0.030000000000000d00, &
       &                                                                                          0.100000000000000d00,0.300000000000000d00, &
       &                                                                                          1.000000000000000d00,3.000000000000000d00, &
       &                                                                                          1.000000000000000d01                       &
       &                                                                                         ]
  double precision                                        , dimension(7) :: energyTarget
  double precision                                        , dimension(7) :: radiusFreefallTarget
  double precision                                        , dimension(7) :: radiusFreefallIncreaseRateTarget
  double precision                                        , dimension(7) :: potentialTarget
  double precision                                        , dimension(7) :: fourierTransformTarget
  double precision                                        , dimension(7) :: velocityDispersionTarget
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelWorking)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Initialize HDF5 lock.
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Tabulated mass distributions")

  ! Iterate over profiles.
  do iProfile=1,2
     select case (iProfile)
     case (1)
        ! Cored NFW profile.
        call Unit_Tests_Begin_Group("Cored NFW profile")
        allocate(massDistributionCoredNFW :: massDistribution_)
        massTarget          =[                                           &
             &                6.370508886341162d05,1.611292088903051d07, &
             &                4.820816988421535d08,7.794652959722677d09, &
             &                8.579226849001460d10,3.725560473368585d11, &
             &                1.000000000000000d12                       &
             &               ]
        densityMoment0Target=[                                           &
             &                2.460730293686941d12,7.069365704837000d12, &
             &                2.048945657552057d13,4.439368231418265d13, &
             &                7.258203857541593d13,8.492694861727370d13, &
             &                8.799525923615780d13                       &
             &               ]
        densityMoment1Target=[                                           &
             &                3.053052459464275d08,2.592947607342040d09, &
             &                2.384847303281488d10,1.367185830935901d11, &
             &                5.322663280946257d11,1.037860585072224d12, &
             &                1.407477797390162d12                       &
             &               ]
        densityMoment2Target=[                                           &
             &                5.069489879599125d04,1.282225499839590d06, &
             &                3.836284256845901d07,6.202787726771978d08, &
             &                6.827131793393528d09,2.964706822499870d10, &
             &                7.957747149596004d10                       &
             &               ]
        densityMoment3Target=[                                           &
             &                9.484085428575420d00,7.165022587728307d02, &
             &                7.043283943903859d04,3.296870652850856d06, &
             &                1.105431237182676d08,1.240614494009899d09, &
             &                8.818052733199510d09                       &
             &               ]
     case (2)
        ! Cusp-NFW profile.
        call Unit_Tests_Begin_Group("Cusp-NFW profile")
        allocate(massDistributionCuspNFW  :: massDistribution_)
        massTarget          =[                                           &
             &                9.374216590052220d07,5.326799595417917d08, &
             &                3.881856746823816d09,2.369796337315842d10, &
             &                1.344565548522473d11,4.319437432755695d11, &
             &                1.000000000000000d12                       &
             &               ]
        densityMoment0Target=[                                           &
             &                0.000000000000000d00,0.000000000000000d00, &
             &                0.000000000000000d00,0.000000000000000d00, &
             &                0.000000000000000d00,0.000000000000000d00, &
             &                0.000000000000000d00                       &
             &               ]
        densityMoment1Target=[                                           &
             &                8.732721039381150d10,1.596308717474747d11, &
             &                3.292682327043573d11,6.564239689274449d11, &
             &                1.246216942232268d12,1.781027285787356d12, &
             &                2.118611334134790d12                       &
             &               ]
        densityMoment2Target=[                                           &
             &                7.459764539604310d06,4.238932432353350d07, &
             &                3.089083448158171d08,1.885824006024425d09, &
             &                1.069971266792086d10,3.437299093996175d10, &
             &                7.957747150213412d10                       &
             &               ]
        densityMoment3Target=[                                           &
             &                1.130431624151061d03,1.952464146474319d04, &
             &                4.819043953615453d05,8.763585073867530d06, &
             &                1.551564332097499d08,1.307501276772168d09, &
             &                8.112997987118110d09                       &
             &               ]
     end select
     select type (massDistribution_)
     type is (massDistributionCoredNFW)
        allocate(kinematicsDistributionCollisionlessTabulated :: kinematicsDistribution_)
        select type (kinematicsDistribution_)
        type is (kinematicsDistributionCollisionlessTabulated)
           kinematicsDistribution_=kinematicsDistributionCollisionlessTabulated(toleranceRelativeVelocityDispersion=1.0d-3,toleranceRelativeVelocityDispersionMaximum=1.0d-3)
        end select
        massDistribution_=massDistributionCoredNFW(mass=massVirial,radiusVirial=radiusVirial,radiusScale=radiusScale,radiusCore=radiusCore,dimensionless=.false.,toleranceRelativePotential=1.0d-3)
        call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
     type is (massDistributionCuspNFW )
        allocate(kinematicsDistributionCuspNFW :: kinematicsDistribution_)
        select type (kinematicsDistribution_)
        type is (kinematicsDistributionCuspNFW)
           kinematicsDistribution_=kinematicsDistributionCuspNFW               (toleranceRelativeVelocityDispersion=1.0d-3,toleranceRelativeVelocityDispersionMaximum=1.0d-3)
        end select
        massDistribution_=massDistributionCuspNFW (mass=massVirial,radiusVirial=radiusVirial,radiusScale=radiusScale,y         =yCusp     ,dimensionless=.false.,toleranceRelativePotential=1.0d-3)
        call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
     end select
     timeScale=+1.0d0                                &
          &    /sqrt(                                &
          &          +gravitationalConstant_internal &
          &          *3.0d0                          &
          &          /4.0d0                          &
          &          /Pi                             &
          &          *massVirial                     &
          &          /radiusVirial**3                &
          &         )                                &
          &    *MpcPerKmPerSToGyr
     select type (massDistribution_)
     class is (massDistributionSpherical)
        coordinatesReference=[radiusScale*radiiScaleFree(1),0.0d0,0.0d0]
       do i=1,7
           coordinates                        =[radiusScale*radiiScaleFree(i),0.0d0,0.0d0]
           mass                            (i)=massDistribution_      %massEnclosedBySphere               (                         radiiScaleFree(i)*radiusScale                                       )
           potential                       (i)=massDistribution_      %potentialDifference                (                                           coordinates                  ,coordinatesReference)
           potentialTarget                 (i)=massDistribution_      %potentialDifferenceNumerical       (                                           coordinates                  ,coordinatesReference)
           energy                          (i)=massDistribution_      %energy                             (                         radiiScaleFree(i)*radiusScale,massDistribution_                     )
           energyTarget                    (i)=massDistribution_      %energyNumerical                    (                         radiiScaleFree(i)*radiusScale,massDistribution_                     )
           velocityDispersion              (i)=kinematicsDistribution_%velocityDispersion1D               (                                           coordinates,massDistribution_,massDistribution_   )
           velocityDispersionTarget        (i)=kinematicsDistribution_%velocityDispersion1DNumerical      (                                           coordinates,massDistribution_,massDistribution_   )
           if (all(densityMoment0Target > 0.0d0)) &
                & densityMoment0           (i)=massDistribution_      %densityRadialMoment                (0.0d0,0.0d0,             radiiScaleFree(i)*radiusScale                                       )
           densityMoment1                  (i)=massDistribution_      %densityRadialMoment                (1.0d0,0.0d0,             radiiScaleFree(i)*radiusScale                                       )
           densityMoment2                  (i)=massDistribution_      %densityRadialMoment                (2.0d0,0.0d0,             radiiScaleFree(i)*radiusScale                                       )
           densityMoment3                  (i)=massDistribution_      %densityRadialMoment                (3.0d0,0.0d0,             radiiScaleFree(i)*radiusScale                                       )
           radiusFreefall                  (i)=massDistribution_      %radiusFreefall                     (                         radiiScaleFree(i)*  timeScale                                       )
           radiusFreefallTarget            (i)=massDistribution_      %radiusFreefallNumerical            (                         radiiScaleFree(i)*  timeScale                                       )
           radiusFreefallIncreaseRate      (i)=massDistribution_      %radiusFreefallIncreaseRate         (                         radiiScaleFree(i)*  timeScale                                       )
           radiusFreefallIncreaseRateTarget(i)=massDistribution_      %radiusFreefallIncreaseRateNumerical(                         radiiScaleFree(i)*  timeScale                                       )
           fourierTransform                (i)=massDistribution_      %fourierTransform                   (            radiusVirial,radiiScaleFree(i)/radiusScale                                       )
           fourierTransformTarget          (i)=massDistribution_      %fourierTransformNumerical          (            radiusVirial,radiiScaleFree(i)/radiusScale                                       )
        end do
        ! Only potential differences are relevant.
        potential      =potential      -potential      (1)
        potentialTarget=potentialTarget-potentialTarget(1)
        ! Test assertions.
        call    Assert("Mass within radius"           ,mass                      ,massTarget                      ,relTol=3.0d-3                                     )
        call    Assert("Energy within radius"         ,energy                    ,energyTarget                    ,relTol=3.0d-2                                     )
        call    Assert("Potential"                    ,potential                 ,potentialTarget                 ,relTol=1.0d-2                                     )
        call    Assert("Velocity dispersion"          ,velocityDispersion        ,velocityDispersionTarget        ,relTol=4.0d-2                                     )
        if (all(densityMoment0Target > 0.0d0)) then
           call Assert("Radial density moment (m=0)"  ,densityMoment0            ,densityMoment0Target            ,relTol=1.0d-3                                     )
        else
           call Skip  ("Radial density moment (m=0)"  ,"divergent"                                                                                                   )
        end if
        call    Assert("Radial density moment (m=1)"  ,densityMoment1            ,densityMoment1Target            ,relTol=1.0d-3                                     )
        call    Assert("Radial density moment (m=2)"  ,densityMoment2            ,densityMoment2Target            ,relTol=1.0d-3                                     )
        call    Assert("Radial density moment (m=3)"  ,densityMoment3            ,densityMoment3Target            ,relTol=2.0d-3                                     )
        call    Assert("Fourier transform"            ,fourierTransform          ,fourierTransformTarget          ,relTol=1.2d-2                                     )
        call    Assert("Freefall radius"              ,radiusFreefall            ,radiusFreefallTarget            ,relTol=1.0d-1,absTol=1.0d-1*radiusCore            )
        call    Assert("Freefall radius increase rate",radiusFreefallIncreaseRate,radiusFreefallIncreaseRateTarget,relTol=1.0d-2,absTol=2.0d+0*radiusVirial/timeScale)
     end select
     deallocate(massDistribution_      )
     nullify   (kinematicsDistribution_)
     call Unit_Tests_End_Group()
  end do
  
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Mass_Distributions_Tabulated
