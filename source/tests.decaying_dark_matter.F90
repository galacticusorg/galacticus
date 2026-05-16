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
  use :: IO_HDF5             , only : hdf5Object                        , ioHDF5AccessInitialize
  use :: Input_Paths         , only : inputPath                         , pathTypeExec
  use :: ISO_Varying_String  , only : char                              , operator(//)                    , varying_string
  use :: Events_Hooks        , only : eventsHooksInitialize
  use :: Decaying_Dark_Matter, only : decayingDarkMatterFractionRetained, decayingDarkMatterEnergyRetained, decayingDarkMatterFractionRetainedDerivatives, decayingDarkMatterEnergyRetainedDerivatives
  implicit none
  double precision               , dimension(:  ), allocatable :: velocityKicks                                 , velocityEscapes
  double precision               , dimension(:,:), allocatable :: fractionRetainedTarget                        , fractionRetainedTargetUncertainty               , &
       &                                                          energyRetainedTarget                          , energyRetainedTargetUncertainty                 , &
       &                                                          fractionRetained                              , energyRetained                                  , &
       &                                                          energyDerivativeVelocityKick                  , fractionDerivativeVelocityKick                  , &
       &                                                          energyDerivativeVelocityEscape                , fractionDerivativeVelocityEscape                , &
       &                                                          energyDerivativeVelocityKickFiniteDifference  , fractionDerivativeVelocityKickFiniteDifference  , &
       &                                                          energyDerivativeVelocityEscapeFiniteDifference, fractionDerivativeVelocityEscapeFiniteDifference
  double precision                                             :: velocityDispersion
  integer                                                      :: i                                             , j
  type            (varying_string)                             :: fileName
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize.
  call eventsHooksInitialize()
  ! Initialize HDF5 lock.
  call ioHDF5AccessInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Decaying dark matter")
  ! Retention fractions.
  call Unit_Tests_Begin_Group("Retention")
  block
    type(hdf5Object) :: file
    !$ call hdf5Access%set()
    fileName=inputPath(pathTypeExec)//"testSuite/data/decayingDarkMatterRetention.hdf5"
    file    =hdf5Object(fileName,readOnly=.true.)
    call file%readDataset  ('velocityKick'               ,velocityKicks                    )
    call file%readDataset  ('velocityEscape'             ,velocityEscapes                  )
    call file%readDataset  ('fractionRetained'           ,fractionRetainedTarget           )
    call file%readDataset  ('fractionRetainedUncertainty',fractionRetainedTargetUncertainty)
    call file%readDataset  ('energyRetained'             ,energyRetainedTarget             )
    call file%readDataset  ('energyRetainedUncertainty'  ,energyRetainedTargetUncertainty  )
    call file%readAttribute('velocityDispersion'         ,velocityDispersion               )
    !$ call hdf5Access%unset()
  end block
  allocate(energyRetained                                  (size(velocityKicks),size(velocityEscapes)))
  allocate(fractionRetained                                (size(velocityKicks),size(velocityEscapes)))
  allocate(energyDerivativeVelocityKick                    (size(velocityKicks),size(velocityEscapes)))
  allocate(energyDerivativeVelocityKickFiniteDifference    (size(velocityKicks),size(velocityEscapes)))
  allocate(fractionDerivativeVelocityKick                  (size(velocityKicks),size(velocityEscapes)))
  allocate(fractionDerivativeVelocityKickFiniteDifference  (size(velocityKicks),size(velocityEscapes)))
  allocate(energyDerivativeVelocityEscape                  (size(velocityKicks),size(velocityEscapes)))
  allocate(energyDerivativeVelocityEscapeFiniteDifference  (size(velocityKicks),size(velocityEscapes)))
  allocate(fractionDerivativeVelocityEscape                (size(velocityKicks),size(velocityEscapes)))
  allocate(fractionDerivativeVelocityEscapeFiniteDifference(size(velocityKicks),size(velocityEscapes)))
  do i=1,size(velocityKicks)
     do j=1,size(velocityEscapes)
        fractionRetained(i,j)=decayingDarkMatterFractionRetained(velocityDispersion,velocityEscapes(j),velocityKicks(i))
        energyRetained  (i,j)=decayingDarkMatterEnergyRetained  (velocityDispersion,velocityEscapes(j),velocityKicks(i))
        call decayingDarkMatterFractionRetainedDerivatives(velocityDispersion,velocityEscapes(j),velocityKicks(i),fractionDerivativeVelocityEscape(i,j),fractionDerivativeVelocityKick(i,j))
        call decayingDarkMatterEnergyRetainedDerivatives  (velocityDispersion,velocityEscapes(j),velocityKicks(i),energyDerivativeVelocityEscape  (i,j),energyDerivativeVelocityKick  (i,j))
        if (i == 1) then
           fractionDerivativeVelocityKickFiniteDifference  (i,j)=fractionDerivativeVelocityKick  (i,j)
           energyDerivativeVelocityKickFiniteDifference    (i,j)=energyDerivativeVelocityKick    (i,j)
        else
           fractionDerivativeVelocityKickFiniteDifference  (i,j)=(fractionRetained(i,j)-fractionRetained(i-1,j))/(velocityKicks  (i)-velocityKicks  (i-1))*velocityDispersion
           energyDerivativeVelocityKickFiniteDifference    (i,j)=(energyRetained  (i,j)-energyRetained  (i-1,j))/(velocityKicks  (i)-velocityKicks  (i-1))*velocityDispersion
        end if
        if (j == 1) then
           fractionDerivativeVelocityEscapeFiniteDifference(i,j)=fractionDerivativeVelocityEscape(i,j)
           energyDerivativeVelocityEscapeFiniteDifference  (i,j)=energyDerivativeVelocityEscape  (i,j)
        else
           fractionDerivativeVelocityEscapeFiniteDifference(i,j)=(fractionRetained(i,j)-fractionRetained(i,j-1))/(velocityEscapes(j)-velocityEscapes(j-1))*velocityDispersion
           energyDerivativeVelocityEscapeFiniteDifference  (i,j)=(energyRetained  (i,j)-energyRetained  (i,j-1))/(velocityEscapes(j)-velocityEscapes(j-1))*velocityDispersion
        end if
     end do
  end do
  do j=1,size(velocityEscapes)
     energyRetained                              (:,j)=+energyRetained                              (:,j)/(0.5d0*velocityKicks**2)
     energyDerivativeVelocityKick                (:,j)=+energyDerivativeVelocityKick                (:,j)/(0.5d0*velocityKicks**2)
     energyDerivativeVelocityKickFiniteDifference(:,j)=+energyDerivativeVelocityKickFiniteDifference(:,j)/(0.5d0*velocityKicks**2)
  end do
  call Assert('fraction; f     ',abs(fractionRetained-fractionRetainedTarget),max(1.0d-5,5.0d0*fractionRetainedTargetUncertainty),compareLessThanOrEqual)
  call Assert('  energy; ε     ',abs(energyRetained  -energyRetainedTarget  ),max(1.0d-5,5.0d0*energyRetainedTargetUncertainty  ),compareLessThanOrEqual)
  call Assert('          ∂f/∂vₖ',fractionDerivativeVelocityKick  ,fractionDerivativeVelocityKickFiniteDifference  ,relTol=1.5d-1,absTol=2.0d-3)
  call Assert('          ∂ε/∂vₖ',energyDerivativeVelocityKick    ,energyDerivativeVelocityKickFiniteDifference    ,relTol=4.0d-1,absTol=2.0d-2)
  call Assert('          ∂f/∂vₑ',fractionDerivativeVelocityEscape,fractionDerivativeVelocityEscapeFiniteDifference,relTol=1.5d-1,absTol=1.0d-3)
  call Assert('          ∂ε/∂vₑ',energyDerivativeVelocityEscape  ,energyDerivativeVelocityEscapeFiniteDifference  ,relTol=6.0d-1,absTol=5.0d-2)
  deallocate(fractionRetainedTarget                        )
  deallocate(fractionRetainedTargetUncertainty             )
  deallocate(energyRetainedTarget                          )
  deallocate(energyRetainedTargetUncertainty               )
  deallocate(energyDerivativeVelocityKick                  )
  deallocate(energyDerivativeVelocityKickFiniteDifference  )
  deallocate(fractionDerivativeVelocityKick                )
  deallocate(fractionDerivativeVelocityKickFiniteDifference)
  call Unit_Tests_End_Group()
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Decaying_Dark_Matter
