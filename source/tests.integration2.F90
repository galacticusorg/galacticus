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
Contains a program to test integration routines.
!!}

program Test_Integration2
  !!{
  Tests that numerical integration routines work.
  !!}
  use :: Display                    , only : displayIndent                   , displayMessage                             , displayUnindent                                 , displayVerbositySet                            , &
          &                                  verbosityLevelStandard
  use :: Error                      , only : Error_Report
  use :: ISO_Varying_String         , only : assignment(=)                   , char                                       , len
  use :: Kind_Numbers               , only : kind_int8
  use :: Numerical_Integration2     , only : integrator1D                    , integrator2                                , integratorAdaptiveCompositeTrapezoidal1D        , integratorCompositeGaussKronrod1D              , &
          &                                  integratorCompositeTrapezoidal1D, integratorMultiVectorized1D                , integratorMultiVectorizedCompositeGaussKronrod1D, integratorMultiVectorizedCompositeTrapezoidal1D, &
          &                                  integratorVectorized1D          , integratorVectorizedCompositeGaussKronrod1D, integratorVectorizedCompositeTrapezoidal1D
  use :: Test_Integration2_Functions, only : testFunctions                   , testFunctionsInitialize                    , testFunctionsMulti                              , testIntegrator                                 , &
          &                                  testIntegratorMulti
  use :: Numerical_Integration      , only : integrator1                      => integrator, GSL_Integ_Gauss15                          , GSL_Integ_Gauss61
  use :: Display                    , only : displayIndent                   , displayMessage                             , displayUnindent                                 , displayVerbositySet                            , &
          &                                  verbosityLevelStandard
  use :: Error                      , only : Error_Report
  use :: ISO_Varying_String         , only : assignment(=)                   , char                                       , len
  use :: Kind_Numbers               , only : kind_int8
  use :: Numerical_Integration2     , only : integrator1D                    , integrator2                                , integratorAdaptiveCompositeTrapezoidal1D        , integratorCompositeGaussKronrod1D              , &
          &                                  integratorCompositeTrapezoidal1D, integratorMultiVectorized1D                , integratorMultiVectorizedCompositeGaussKronrod1D, integratorMultiVectorizedCompositeTrapezoidal1D, &
          &                                  integratorVectorized1D          , integratorVectorizedCompositeGaussKronrod1D, integratorVectorizedCompositeTrapezoidal1D
  use :: Test_Integration2_Functions, only : testFunctions                   , testFunctionsInitialize                    , testFunctionsMulti                              , testIntegrator                                 , &
          &                                  testIntegratorMulti
  implicit none
  double precision                                                                             :: timeMean                       , error
  integer         (kind=kind_int8            )                                                 :: countStart                     , countEnd                , &
       &                                                                                          countRate
  character       (len =   3                 )                                                 :: units
  character       (len = 128                 )                                                 :: label                          , formatter               , &
       &                                                                                          status
  character       (len =1024                 )                                                 :: message
  double precision                            , parameter                                      :: toleranceRelative       =1.0d-6, toleranceAbsolute=1.0d-6
  integer                                                                                      :: i                              , trial                   , &
       &                                                                                          integrationRule                , iFunction               , &
       &                                                                                          descriptionLengthMaximum
  type            (integrator1               )                                   , allocatable :: integrator1_
  integer                                     , parameter                                      :: integratorCount         =  7
  integer                                     , parameter                                      :: integratorMultiCount    =  2
  type            (testIntegrator            ), dimension(integratorCount       )              :: integrators
  type            (testIntegratorMulti       ), dimension(integratorMultiCount  )              :: integratorsMulti
  double precision                            , dimension(integratorCount       )              :: integral
  double precision                            , dimension(integratorMultiCount,2)              :: integralMulti
  double precision                            , dimension(                     2)              :: integralGSL
  double precision                            , dimension(                     :), allocatable :: tolerancesAbsolute, tolerancesRelative
  integer         (kind=kind_int8            ), dimension(integratorCount       )              :: time
  integer         (kind=kind_int8            ), dimension(integratorMultiCount  )              :: timeMulti
  integer         (kind=kind_int8            ), dimension(                     2)              :: timeGSL
  integer                                     , parameter                                      :: trialCount              =10

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Determine units of system clock.
  call System_Clock(countStart,countRate)
  select case (countRate)
  case (      1000)
     units="ms"
  case (   1000000)
     units="μs"
  case (1000000000)
     units="ns"
  end select
  ! Iterate over functions.
  call testFunctionsInitialize()
  do iFunction=1,size(testFunctions)
     write (message,'(a,a,a,f5.2,a,f5.2)') "Test: ∫ "                                , &
          &                                trim(testFunctions(iFunction)%description), &
          &                                " dx; from x="                            , &
          &                                testFunctions(iFunction)%rangeLow         , &
          &                                " to "                                    , &
          &                                testFunctions(iFunction)%rangeHigh
     call displayIndent(message)
     ! Allocate integrators.
     allocate(integratorCompositeTrapezoidal1D                 :: integrators     (1)%integrator_)
     integrators     (1)%description="Scalar composite trapezoidal"
     allocate(integratorAdaptiveCompositeTrapezoidal1D         :: integrators     (2)%integrator_)
     integrators     (2)%description="Scalar adaptive composite trapezoidal"
     allocate(integratorVectorizedCompositeTrapezoidal1D       :: integrators     (3)%integrator_)
     integrators     (3)%description="Vector adaptive composite trapezoidal"
     allocate(integratorCompositeGaussKronrod1D                :: integrators     (4)%integrator_)
     integrators     (4)%description="Scalar adaptive composite Gauss-Kronrod 15-point"
     integrators     (4)%order   =15
     allocate(integratorCompositeGaussKronrod1D                :: integrators     (5)%integrator_)
     integrators     (5)%description="Scalar adaptive composite Gauss-Kronrod 61-point"
     integrators     (5)%order   =61
     allocate(integratorVectorizedCompositeGaussKronrod1D      :: integrators     (6)%integrator_)
     integrators     (6)%description="Vector adaptive composite Gauss-Kronrod 15-point"
     integrators     (6)%order   =15
     allocate(integratorVectorizedCompositeGaussKronrod1D      :: integrators     (7)%integrator_)
     integrators     (7)%description="Vector adaptive composite Gauss-Kronrod 61-point"
     integrators     (7)%order   =61
     ! Find the longest description.
     descriptionLengthMaximum=26
     do i=1,integratorCount
        descriptionLengthMaximum=max(descriptionLengthMaximum,len(integrators(i)%description))
     end do
     write (formatter,'(a,i3.3,a)') '(a',descriptionLengthMaximum,',": ⏲=",f16.3,1x,a,1x,"∫=",f14.10,1x,"ℯ=",e14.8," [",a,"]")'
     ! Initialize integrators.
     do i=1,integratorCount
        select type (integrator_ => integrators(i)%integrator_)
        class is (integratorCompositeTrapezoidal1D           )
           call integrator_%initialize  (24                                    )
           call integrator_%toleranceSet(toleranceAbsolute,toleranceRelative   )
        class is (integratorAdaptiveCompositeTrapezoidal1D   )
           call integrator_%initialize  (24                                    )
           call integrator_%toleranceSet(toleranceAbsolute,toleranceRelative   )
        class is (integratorVectorizedCompositeTrapezoidal1D )
           call integrator_%initialize  (24                                    )
           call integrator_%toleranceSet(toleranceAbsolute,toleranceRelative   )
        class is (integratorCompositeGaussKronrod1D          )
           call integrator_%initialize  (24               ,integrators(i)%order)
           call integrator_%toleranceSet(toleranceAbsolute,toleranceRelative   )
        class is (integratorVectorizedCompositeGaussKronrod1D)
           call integrator_%initialize  (24               ,integrators(i)%order)
           call integrator_%toleranceSet(toleranceAbsolute,toleranceRelative   )
        class default
           call Error_Report('unknown integrator class [1]'//{introspection:location})
        end select
     end do
     ! Assign functions.
     do i=1,integratorCount
        select type (integrator_ => integrators(i)%integrator_)
        class is (integrator1D)
           call    integrator_%integrandSet(testFunctions(iFunction)%scalar)
        class is (integratorVectorized1D)
           call integrator_%integrandSet(testFunctions(iFunction)%vector)
        class default
           call Error_Report('unknown integrator class [2]'//{introspection:location})
        end select
     end do
     ! Initialize times to zero.
     time   =0
     timeGSL=0
     ! Iterate over trials.
     do trial=1,trialCount
        ! Evaluate internal integrators.
        do i=1,integratorCount
           select type (integrator_ => integrators(i)%integrator_)
           class is (integrator1D)
              call System_Clock(countStart,countRate)
              integral(i)=integrator_%evaluate(testFunctions(iFunction)%rangeLow,testFunctions(iFunction)%rangeHigh)
              call System_Clock(countEnd  ,countRate)
              time(i)=time(i)+(countEnd-countStart)
           class is (integratorVectorized1D)
              call System_Clock(countStart,countRate)
              integral(i)=integrator_%evaluate(testFunctions(iFunction)%rangeLow,testFunctions(iFunction)%rangeHigh)
              call System_Clock(countEnd  ,countRate)
              time(i)=time(i)+(countEnd-countStart)
           class default
              call Error_Report('unknown integrator class [3]'//{introspection:location})
           end select          
        end do
        ! Evaluate GSL integrators.
        do i=1,2
           allocate(integrator1_)
           select case (i)
           case (1)
              integrationRule=GSL_Integ_Gauss15
           case (2)
              integrationRule=GSL_Integ_Gauss61
           end select
           integrator1_=integrator1(                                                       &
                &                                     testFunctions    (iFunction)%scalar, &
                &                   toleranceAbsolute=toleranceAbsolute                  , &
                &                   toleranceRelative=toleranceRelative                  , &
                &                   integrationRule  =integrationRule                      &
                &                  )
           call System_Clock(countStart,countRate)
           integralGSL(i)=integrator1_%integrate(                                    &
                &                                testFunctions(iFunction)%rangeLow , &
                &                                testFunctions(iFunction)%rangeHigh  &
                &                               )
           deallocate(integrator1_)
           call System_Clock  (countEnd         ,countRate           )
           timeGSL(i)=timeGSL(i)+(countEnd-countStart)
        end do
     end do
     ! Report.
     do i=1,integratorCount
       timeMean=dble(time(i))/dble(trialCount)
        error= +abs(integral(i)-testFunctions(iFunction)%solution) &
             & /abs(            testFunctions(iFunction)%solution)
        if (error > toleranceAbsolute/abs(testFunctions(iFunction)%solution) .and. error > toleranceRelative) then
           status="FAIL"
        else
           status="pass"
        end if
        write (message,formatter) char(integrators(i)%description), &
             &                    timeMean                        , &
             &                    units                           , &
             &                    integral(i)                     , &
             &                    error                           , &
             &                    trim(status)
        call displayMessage(message)
     end do
     do i=1,2
        timeMean=dble(timeGSL(i))/dble(trialCount)
        select case (i)
        case (1)
           label="GSL Gauss-Kronrod 15-point"
        case (2)
           label="GSL Gauss-Kronrod 61-point"
        end select
         error= +abs(integralGSL(i)-testFunctions(iFunction)%solution) &
              & /abs(               testFunctions(iFunction)%solution)
        if (error > toleranceAbsolute/abs(testFunctions(iFunction)%solution) .and. error > toleranceRelative) then
           status="FAIL"
        else
           status="pass"
        end if
        write (message,formatter) trim(label)   , &
             &                    timeMean      , &
             &                    units         , &
             &                    integralGSL(i), &
             &                    error         , &
             &                    trim(status)
        call displayMessage(message)
     end do
     ! Destroy integrators.
     do i=1,integratorCount
        deallocate(integrators(i)%integrator_)
     end do
     call displayUnindent("done")
  end do
  ! Iterate over multi-integrand functions.
  do iFunction=1,size(testFunctionsMulti)
     write (message,'(a,a,a,f5.2,a,f5.2)') "Test: ∫ "                                     , &
          &                                trim(testFunctionsMulti(iFunction)%description), &
          &                                " dx; from x="                                 , &
          &                                testFunctionsMulti(iFunction)%rangeLow         , &
          &                                " to "                                         , &
          &                                testFunctionsMulti(iFunction)%rangeHigh
     call displayIndent(message)
     ! Allocate integrators.
     allocate(integratorMultiVectorizedCompositeTrapezoidal1D  :: integratorsMulti( 1)%integrator_)
     integratorsMulti( 1)%description="Multi-integrand vector adaptive composite trapezoidal"
     allocate(integratorMultiVectorizedCompositeGaussKronrod1D :: integratorsMulti( 2)%integrator_)
     integratorsMulti( 2)%description="Multi-integrand vector adaptive composite Gauss-Kronrod 61-point"
     integratorsMulti( 2)%order   =61
     ! Find the longest description.
     descriptionLengthMaximum=26
     do i=1,integratorMultiCount
        descriptionLengthMaximum=max(descriptionLengthMaximum,len(integratorsMulti(i)%description))
     end do
     write (formatter,'(a,i3.3,a)') '(a',descriptionLengthMaximum,',": ⏲=",f16.3,1x,a,1x,"∫=",f14.10,",",f14.10,1x,"ℯ=",e14.8," [",a,"]")'
     ! Initialize integrators.
     do i=1,integratorMultiCount
        select type (integrator_ => integratorsMulti(i)%integrator_)
        class is (integratorMultiVectorizedCompositeTrapezoidal1D)
           call integrator_%initialize   (1000000                                    )
           allocate(tolerancesAbsolute(size(testFunctionsMulti(iFunction)%solution)))
           allocate(tolerancesRelative(size(testFunctionsMulti(iFunction)%solution)))
           tolerancesAbsolute=toleranceAbsolute
           tolerancesRelative=toleranceRelative
           call integrator_%tolerancesSet(tolerancesAbsolute,tolerancesRelative)
           deallocate(tolerancesAbsolute)
           deallocate(tolerancesRelative)
        class is (integratorMultiVectorizedCompositeGaussKronrod1D)
           call integrator_%initialize   (1000000          ,integratorsMulti(i)%order)
           allocate(tolerancesAbsolute(size(testFunctionsMulti(iFunction)%solution)))
           allocate(tolerancesRelative(size(testFunctionsMulti(iFunction)%solution)))
           tolerancesAbsolute=toleranceAbsolute
           tolerancesRelative=toleranceRelative
           call integrator_%tolerancesSet(tolerancesAbsolute,tolerancesRelative)
           deallocate(tolerancesAbsolute)
           deallocate(tolerancesRelative)
        class default
           call Error_Report('unknown integrator class [1.m]'//{introspection:location})
        end select
     end do
      ! Assign functions.
      do i=1,integratorMultiCount
         select type (integrator_ => integratorsMulti(i)%integrator_)
         class is (integratorMultiVectorized1D)
            call integrator_%integrandSet(size(testFunctionsMulti(iFunction)%solution),testFunctionsMulti(iFunction)%vector)
         class default
            call Error_Report('unknown integrator class [2.m]'//{introspection:location})
         end select
      end do
      ! Initialize times to zero.
      timeMulti=0
      ! Iterate over trials.
      do trial=1,trialCount
         ! Evaluate internal integrators.
         do i=1,integratorMultiCount
            select type (integrator_ => integratorsMulti(i)%integrator_)
            class is (integratorMultiVectorized1D)
               call System_Clock(countStart,countRate)
               integralMulti(i,:)=integrator_%evaluate(testFunctionsMulti(iFunction)%rangeLow,testFunctionsMulti(iFunction)%rangeHigh)
               call System_Clock(countEnd  ,countRate)
               timeMulti(i)=timeMulti(i)+(countEnd-countStart)
               class default
               call Error_Report('unknown integrator class [3.m]'//{introspection:location})
            end select
         end do
      end do
      ! Report.
      do i=1,integratorMultiCount
         timeMean=+dble(timeMulti(i))/dble(trialCount)
         error   =+maxval(                                                                &
              &           +abs(integralMulti(i,:)-testFunctionsMulti(iFunction)%solution) &
              &           /abs(                  +testFunctionsMulti(iFunction)%solution) &
             &          )
        if (error > toleranceAbsolute/abs(testFunctions(iFunction)%solution) .and. error > toleranceRelative) then
           status="FAIL"
        else
           status="pass"
        end if
        write (message,formatter) char(integratorsMulti(i)%description), &
             &                    timeMean                             , &
             &                    units                                , &
             &                    integralMulti(i,:)                   , &
             &                    error                                , &
             &                    trim(status)
        call displayMessage(message)
     end do
     ! Destroy integrators.
     do i=1,integratorMultiCount
        deallocate(integratorsMulti(i)%integrator_)
     end do
     call displayUnindent("done")
  end do
end program Test_Integration2
