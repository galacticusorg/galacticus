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
Contains a program to test decaying dark matter calculations.
!!}

program Test_Decaying_Dark_Matter
  !!{
  Tests of decaying dark matter calculations.
  !!}
  use :: Display             , only : displayVerbositySet               , verbosityLevelStandard
  use :: Unit_Tests          , only : Assert                            , Unit_Tests_Begin_Group          , Unit_Tests_End_Group                         , Unit_Tests_Finish                          , &
       &                              compareLessThanOrEqual
  use :: HDF5_Access         , only : hdf5Access
  use :: IO_HDF5             , only : hdf5Object
  use :: Input_Paths         , only : inputPath                         , pathTypeExec
  use :: ISO_Varying_String  , only : char                              , operator(//)
  use :: Events_Hooks        , only : eventsHooksInitialize
  use :: Decaying_Dark_Matter, only : decayingDarkMatterFractionRetained, decayingDarkMatterEnergyRetained, decayingDarkMatterFractionRetainedDerivatives, decayingDarkMatterEnergyRetainedDerivatives
  implicit none
  type            (hdf5Object)                            :: file
  double precision            , dimension(:), allocatable :: fractionRetainedTarget                      , fractionRetainedTargetUncertainty             , &
       &                                                     energyRetainedTarget                        , energyRetainedTargetUncertainty               , &
       &                                                     fractionRetained                            , energyRetained                                , &
       &                                                     energyDerivativeVelocityKick                , fractionDerivativeVelocityKick                , &
       &                                                     energyDerivativeVelocityKickFiniteDifference, fractionDerivativeVelocityKickFiniteDifference, &
       &                                                     velocityKicks
  double precision                                        :: velocityDispersion                          , velocityEscape                                , &
       &                                                     fractionDerivativeVelocityEscapeScaleFree   , energyDerivativeVelocityEscapeScaleFree
  integer                                                 :: i
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize.
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Decaying dark matter")
  ! Retention fractions.
  call Unit_Tests_Begin_Group("Retention")
  !$ call hdf5Access%set()
  call file%openFile(char(inputPath(pathTypeExec)//"testSuite/data/decayingDarkMatterRetention.hdf5"),readOnly=.true.)
  call file%readDataset  ('velocityKick'               ,velocityKicks                    )
  call file%readDataset  ('fractionRetained'           ,fractionRetainedTarget           )
  call file%readDataset  ('fractionRetainedUncertainty',fractionRetainedTargetUncertainty)
  call file%readDataset  ('energyRetained'             ,energyRetainedTarget             )
  call file%readDataset  ('energyRetainedUncertainty'  ,energyRetainedTargetUncertainty  )
  call file%readAttribute('velocityDispersion'         ,velocityDispersion               )
  call file%readAttribute('velocityEscape'             ,velocityEscape                   )
  call file%close        (                                                               )
  !$ call hdf5Access%unset()
  allocate(energyRetained                                (size(velocityKicks)))
  allocate(fractionRetained                              (size(velocityKicks)))
  allocate(energyDerivativeVelocityKick                  (size(velocityKicks)))
  allocate(energyDerivativeVelocityKickFiniteDifference  (size(velocityKicks)))
  allocate(fractionDerivativeVelocityKick                (size(velocityKicks)))
  allocate(fractionDerivativeVelocityKickFiniteDifference(size(velocityKicks)))
  do i=1,size(velocityKicks)
     fractionRetained(i)=decayingDarkMatterFractionRetained(velocityDispersion,velocityEscape,velocityKicks(i))
     energyRetained  (i)=decayingDarkMatterEnergyRetained  (velocityDispersion,velocityEscape,velocityKicks(i))
     call decayingDarkMatterFractionRetainedDerivatives(velocityDispersion,velocityEscape,velocityKicks(i),fractionDerivativeVelocityEscapeScaleFree,fractionDerivativeVelocityKick(i))
     call decayingDarkMatterEnergyRetainedDerivatives  (velocityDispersion,velocityEscape,velocityKicks(i),energyDerivativeVelocityEscapeScaleFree  ,energyDerivativeVelocityKick  (i))
     if (i == 1) then
        fractionDerivativeVelocityKickFiniteDifference(i)=fractionDerivativeVelocityKick(i)
        energyDerivativeVelocityKickFiniteDifference  (i)=energyDerivativeVelocityKick  (i)
     else
        fractionDerivativeVelocityKickFiniteDifference(i)=(fractionRetained(i)-fractionRetained(i-1))/(velocityKicks(i)-velocityKicks(i-1))*velocityDispersion
        energyDerivativeVelocityKickFiniteDifference  (i)=(energyRetained  (i)-energyRetained  (i-1))/(velocityKicks(i)-velocityKicks(i-1))*velocityDispersion
     end if
  end do
  energyRetained                              =+energyRetained                              /(0.5d0*velocityKicks**2)
  energyDerivativeVelocityKick                =+energyDerivativeVelocityKick                /(0.5d0*velocityKicks**2)
  energyDerivativeVelocityKickFiniteDifference=+energyDerivativeVelocityKickFiniteDifference/(0.5d0*velocityKicks**2)
  call Assert('fraction; f     ',abs(fractionRetained-fractionRetainedTarget),max(1.0d-5,5.0d0*fractionRetainedTargetUncertainty),compareLessThanOrEqual)
  call Assert('  energy; ε     ',abs(energyRetained  -energyRetainedTarget  ),max(1.0d-5,5.0d0*energyRetainedTargetUncertainty  ),compareLessThanOrEqual)
  call Assert('          ∂f/∂vₖ',fractionDerivativeVelocityKick,fractionDerivativeVelocityKickFiniteDifference,relTol=1.5d-1,absTol=1.0d-3)
  call Assert('          ∂f/∂vₑ',energyDerivativeVelocityKick  ,energyDerivativeVelocityKickFiniteDifference  ,relTol=4.0d-1,absTol=2.0d-2)
  deallocate(fractionRetainedTarget                        )
  deallocate(fractionRetainedTargetUncertainty             )
  deallocate(energyRetainedTarget                          )
  deallocate(energyRetainedTargetUncertainty               )
  deallocate(energyDerivativeVelocityKick                  )
  deallocate(energyDerivativeVelocityKickFiniteDifference  )
  deallocate(fractionDerivativeVelocityKick                )
  deallocate(fractionDerivativeVelocityKickFiniteDifference)
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish    ()

end program Test_Decaying_Dark_Matter
