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
Contains a program which tests nonlinear collapse solution in an Einstein-de Sitter cosmology.
!!}

program Tests_Spherical_Collapse_NonLinear
  !!{
  Tests nonlinear collapse solution in an Einstein-de Sitter cosmology.
  !!}
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Display                   , only : displayVerbositySet                             , verbosityLevelStandard
  use :: Linear_Growth             , only : linearGrowthCollisionlessMatter
  use :: Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt
  use :: Tables                    , only : table2DLinLinLin
  use :: Unit_Tests                , only : Assert                                          , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  ! Target nonlinear overdensities are computed using the analytic solutions for expanding/collapsing perturbations in an
  ! Einstein-de Sitter universe.
  double precision                                                  , dimension(5), parameter :: overdensityLinear         =[-1.500000d0,-1.000000d0,+0.000000d0,+1.00000d0,+1.50000d0], &
       &                                                                                         overdensityNonLinearTarget=[-0.653016d0,-0.540387d0,+0.000000d0,+3.69477d0,+4.82042d1]
  double precision                                                  , dimension(5)            :: overdensityNonLinear
  type            (table2DLinLinLin                                )                          :: linearToNonLinear
  type            (cosmologyFunctionsMatterLambda                  )                          :: cosmologyFunctions_
  type            (cosmologyParametersSimple                       )                          :: cosmologyParameters_
  type            (linearGrowthCollisionlessMatter                 )                          :: linearGrowth_
  type            (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)                          :: sphericalCollapseSolver_
  integer                                                                                     :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize cosmology.
  cosmologyParameters_    =cosmologyParametersSimple                       (                                      &
       &                                                                    OmegaMatter    = 1.00d0             , &
       &                                                                    OmegaBaryon    = 0.00d0             , &
       &                                                                    OmegaDarkEnergy= 0.00d0             , &
       &                                                                    temperatureCMB = 2.73d0             , &
       &                                                                    HubbleConstant =70.00d0               &
       &                                                                   )
  cosmologyFunctions_     =cosmologyFunctionsMatterLambda                  (                                      &
       &                                                                                    cosmologyParameters_  &
       &                                                                   )
  linearGrowth_           =linearGrowthCollisionlessMatter                 (                                      &
       &                                                                                    cosmologyParameters_, &
       &                                                                                    cosmologyFunctions_   &
       &                                                                   )
  sphericalCollapseSolver_=sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt(                                      &
       &                                                                                    cosmologyFunctions_ , &
       &                                                                                    linearGrowth_         &
       &                                                                   )
  ! Get a table of linear to nonlinear overdensity mapping.
  call sphericalCollapseSolver_%linearNonlinearMap(cosmologyFunctions_%cosmicTime(1.0d0),linearToNonLinear)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: nonlinear solutions")
  do i=1,size(overdensityLinear)
     overdensityNonLinear(i)=linearToNonLinear%interpolate(overdensityLinear(i),cosmologyFunctions_%cosmicTime(1.0d0))
  end do
  call Assert('Linear → nonlinear spherical collapse mapping',overdensityNonLinear,overdensityNonLinearTarget,relTol=1.0d-2,absTol=1.0d-3)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Spherical_Collapse_NonLinear
