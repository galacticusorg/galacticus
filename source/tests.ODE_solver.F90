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
Contains a program to test ODE-IV2 solver routines.
!!}

program Test_ODE_Solver
  !!{
  Tests that ODE solver routines work.
  !!}
  use            :: Display                  , only : displayVerbositySet                             , verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding            , only : C_Null_FunPtr                                   , c_size_t
  use            :: Numerical_Integration2   , only : integratorMultiVectorizedCompositeGaussKronrod1D
  use            :: Numerical_ODE_Solvers    , only : gsl_odeiv2_step_msbdf                           , gsl_odeiv2_step_msbdfactive, odeSolver
  use            :: Test_ODE_Solver_Functions, only : Integrands_Set_2                                , Jacobian_Set_1             , Jacobian_Set_2      , ODE_Set_1        , &
          &                                           ODE_Set_2
  use            :: Unit_Tests               , only : Assert                                          , Unit_Tests_Begin_Group     , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                                                  , dimension(10  ) :: xEnd
  double precision                                                  , dimension(   2) :: y            , z
  ! Solutions of the latent variables integrals. Computed in Mathematica.
  double precision                                                  , dimension(10,2) :: latentSolution=reshape(                                                                                                                                                                            &
       &                                                                                                        [                                                                                                                                                                           &
       &                                                                                                         0.765949925574d0, 1.72843953777d0, 2.52176907333d0, 3.26649435359d0, 4.21693267110d0, 5.04253217602d0, 5.77282469398d0, 6.70065044178d0, 7.56127579166d0, 8.28419598523d0, &
       &                                                                                                         0.896393789463d0, 1.61919857807d0, 2.48093194847d0, 3.40790822044d0, 4.13786263655d0, 4.96452886147d0, 5.91440877877d0, 6.65857801213d0, 7.45285084203d0, 8.41514219895d0  &
       &                                                                                                        ]                                                                                                                                                                         , &
       &                                                                                                        [10,2]                                                                                                                                                                      &
       &                                                                                                       )
  type            (odeSolver                                       ), allocatable     :: solver
  integer                                                                             :: i
  double precision                                                                    :: xStart
  character       (len=32                                          )                  :: message
  type            (integratorMultiVectorizedCompositeGaussKronrod1D)                  :: integrator_

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("ODE-IV2 solver")

  ! Sinusoid.
  call Unit_Tests_Begin_Group("y′=sin(x)")
  allocate(solver)
  solver=odeSolver(1_c_size_t,ODE_Set_1,Jacobian_Set_1,toleranceAbsolute=1.0d-9,toleranceRelative=1.0d-9,scale=[1.0d0],stepperType=gsl_odeiv2_step_msbdf)
  xEnd=[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0]
  do i=1,size(xEnd)
     y(1:1)=[0.0d0]
     xStart=0.0d0
     call solver%solve(xStart,xEnd(i),y(1:1))
     write (message,'(a,f4.1)') "x=0 to ",xEnd(i)
     call Assert(trim(message),y(1:1),1.0d0-cos(xEnd(i:i)),absTol=1.0d-6,relTol=5.0d-6)
  end do
  deallocate(solver)
  call Unit_Tests_End_Group()

  ! Sinusoid - reverse.
  call Unit_Tests_Begin_Group("y′=sin(x) reversed")
  allocate(solver)
  solver=odeSolver(1_c_size_t,ODE_Set_1,toleranceAbsolute=1.0d-9,toleranceRelative=1.0d-9,scale=[1.0d0])
  xEnd=[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0]
  do i=1,size(xEnd)
     y(1:1)=[1.0d0-cos(xEnd(i:i))]
     xStart=xEnd(i)
     call solver%solve(xStart,0.0d0,y(1:1))
     write (message,'(a,f4.1,a)') "x=",xEnd(i)," to 0"
     call Assert(trim(message),y(1:1),[0.0d0],absTol=1.0d-6,relTol=5.0d-6)
  end do
  deallocate(solver)
  call Unit_Tests_End_Group()

  ! Harmonic oscillator.
  call Unit_Tests_Begin_Group("y″=-y; z′=1/√{1+y²}")
  xEnd=[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0]
  call integrator_%initialize   (24,61)
  call integrator_%tolerancesSet([1.0d-9,1.0d-9],[1.0d-9,1.0d-9])
  allocate(solver)
  solver=odeSolver(2_c_size_t,ODE_Set_2,Jacobian_Set_2,toleranceAbsolute=1.0d-9,toleranceRelative=1.0d-9,scale=[1.0d0,1.0d0],stepperType=gsl_odeiv2_step_msbdfactive,integrator=integrator_,integrands=Integrands_Set_2)
  do i=1,size(xEnd)
     y(1:2)=[1.0d0,0.0d0]
     z(1:2)=[0.0d0,0.0d0]
     xStart= 0.0d0
     call solver%solve(xStart,xEnd(i),y=y(1:2),z=z(1:2))
     write (message,'(a,f4.1)') "active: x=0 to ",xEnd(i)
     call Assert(trim(message),y(1:2),[cos(xEnd(i)),-sin(xEnd(i))],relTol=1.0d-6)
     write (message,'(a,f4.1)') "latent: x=0 to ",xEnd(i)
     call Assert(trim(message),z(1:2),latentSolution(i,:),relTol=1.0d-6)
  end do
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_ODE_Solver
