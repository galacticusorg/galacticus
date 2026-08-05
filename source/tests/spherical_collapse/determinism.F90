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

!!{RST
Contains a program which tests that spherical collapse tabulations are deterministic.
!!}

program Tests_Spherical_Collapse_Determinism
  !!{RST
  Tests that the tabulated spherical collapse solutions are pinned to an absolute lattice, and so are independent of the history
  of requests made of them. For each of the three spherical collapse solvers a table is built by requesting a single time, and
  is then extended by requesting a time far outside of the tabulated range. The test requires that:

  #. the tabulated range is pinned to whole decades (or octaves), so that it takes one of a small set of discrete values
     irrespective of the requested time;
  #. extending the table changes neither the abscissae nor the values of any previously computed point, bit-for-bit---it is the
     absence of this property which necessitated the warm-up hack in the ``Cole2000`` merger tree builder and the two-attempt
     workaround in the ``Font2008`` ram pressure stripping class;
  #. a table subsequently obtained over the wider range is identical to the table built by extension.
  !!}
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda                  , cosmologyFunctionsMatterDarkEnergy
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Display                   , only : displayVerbositySet                             , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Numerical_Ranges          , only : rangeLattice                                    , enumerationGridSchemeType                    , gridSchemePerDecade                               , &
       &                                    gridSchemePerOctave
  use :: Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt, sphericalCollapseSolverCllsnlssMttrDarkEnergy, sphericalCollapseSolverBaryonsDarkMatterDarkEnergy, &
       &                                    cllsnlssMttrDarkEnergyFixedAtTurnaround
  use :: Tables                    , only : table1D
  use :: Unit_Tests                , only : Assert                                          , Unit_Tests_Begin_Group                       , Unit_Tests_End_Group                              , &
       &                                    Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple                         )           :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                    )           :: cosmologyFunctionsLambda_
  type            (cosmologyFunctionsMatterDarkEnergy                )           :: cosmologyFunctionsDarkEnergy_
  type            (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt  )           :: solverCosmologicalConstant_
  type            (sphericalCollapseSolverCllsnlssMttrDarkEnergy     )           :: solverDarkEnergy_
  type            (sphericalCollapseSolverBaryonsDarkMatterDarkEnergy)           :: solverBaryons_
  ! Times, in Gyr, at which the tabulations are triggered. The second lies far outside of the range which will have been
  ! tabulated in response to the first, and so forces the table to be extended.
  double precision                                                   , parameter :: timeFirst                    =10.00d+0,                     &
       &                                                                            timeSecond                   = 5.00d-2
  integer                                                            , parameter :: pointsPerOctave              = 4
  integer                                                            , parameter :: solverCosmologicalConstant   = 1      , solverDarkEnergy=2, &
       &                                                                            solverBaryons                = 3

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: determinism of tabulations")
  ! Construct the cosmology.
  cosmologyParameters_        =cosmologyParametersSimple                       (                                                                       &
       &                                                                        OmegaMatter                = 0.2750d0                                , &
       &                                                                        OmegaBaryon                = 0.0458d0                                , &
       &                                                                        OmegaDarkEnergy            = 0.7250d0                                , &
       &                                                                        temperatureCMB             = 2.7800d0                                , &
       &                                                                        HubbleConstant             =70.2000d0                                  &
       &                                                                       )
  cosmologyFunctionsLambda_   =cosmologyFunctionsMatterLambda                  (                                                                       &
       &                                                                        cosmologyParameters_       =cosmologyParameters_                       &
       &                                                                       )
  cosmologyFunctionsDarkEnergy_=cosmologyFunctionsMatterDarkEnergy             (                                                                       &
       &                                                                        cosmologyParameters_       =cosmologyParameters_                     , &
       &                                                                        darkEnergyEquationOfStateW0=-1.0d0                                   , &
       &                                                                        darkEnergyEquationOfStateW1=+0.0d0                                     &
       &                                                                       )
  ! Construct the solvers. Table storage is disabled throughout so that the test neither depends upon, nor perturbs, any
  ! tabulations cached in the dynamic data path by other runs.
  solverCosmologicalConstant_ =sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt  (                                                                     &
       &                                                                          cosmologyFunctions_      =cosmologyFunctionsLambda_                  &
       &                                                                         )
  solverDarkEnergy_           =sphericalCollapseSolverCllsnlssMttrDarkEnergy     (                                                                     &
       &                                                                          energyFixedAt             =cllsnlssMttrDarkEnergyFixedAtTurnaround , &
       &                                                                          cosmologyFunctions_       =cosmologyFunctionsDarkEnergy_             &
       &                                                                         )
  solverBaryons_              =sphericalCollapseSolverBaryonsDarkMatterDarkEnergy(                                                                     &
       &                                                                          baryonsCluster            =.false.                                 , &
       &                                                                          tablePointsPerOctave      =pointsPerOctave                         , &
       &                                                                          energyFixedAt             =cllsnlssMttrDarkEnergyFixedAtTurnaround , &
       &                                                                          perturbationSmall         =1.0d-3                                  , &
       &                                                                          cosmologyParameters_      =cosmologyParameters_                    , &
       &                                                                          cosmologyFunctions_       =cosmologyFunctionsDarkEnergy_             &
       &                                                                         )
  ! Test each solver in turn. The final argument is the anchor interval, in lattice steps, which each solver pins its bounds to -
  ! the per-decade solvers anchor to half decades (500 of their 1000 points per decade), while the per-octave solver anchors to
  ! whole octaves.
  call testSolver(solverCosmologicalConstant,'collisionless matter + Λ' ,gridSchemePerDecade,500            )
  call testSolver(solverDarkEnergy          ,'collisionless matter + DE',gridSchemePerDecade,500            )
  call testSolver(solverBaryons             ,'baryons + DM + DE'        ,gridSchemePerOctave,pointsPerOctave)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  subroutine testSolver(solver,label,scheme,anchorEvery)
    !!{RST
    Test that the tabulation built by the given solver is pinned to an absolute lattice, and that extending it preserves both
    the abscissae and the values of every previously computed point.
    !!}
    implicit none
    integer                                    , intent(in   )               :: solver
    character       (len=*                    ), intent(in   )               :: label
    type            (enumerationGridSchemeType), intent(in   )               :: scheme
    integer                                    , intent(in   )               :: anchorEvery
    class           (table1D                  ), allocatable                 :: table_       , tableRepeat_
    type            (rangeLattice             )                              :: latticeFirst , latticeSecond
    double precision                           , allocatable, dimension(:  ) :: xValuesFirst , xValuesSecond, &
         &                                                                      xValuesRepeat
    double precision                           , allocatable, dimension(:,:) :: yValuesFirst , yValuesSecond, &
         &                                                                      yValuesRepeat
    integer                                                                  :: offset

    ! Build the table by requesting a single time.
    call getTable(solver,timeFirst ,table_)
    latticeFirst =table_%lattice
    xValuesFirst =table_%xs     ()
    yValuesFirst =table_%ys     ()
    call Assert(label//': table is built on an absolute lattice'  ,latticeFirst%isDefined()     ,.true.)
    call Assert(label//': table uses the expected gridding scheme',latticeFirst%scheme == scheme,.true.)
    ! Since the bounds are pinned, the lattice index of each must be an exact multiple of the anchor interval.
    call Assert(label//': lower bound is pinned to the anchor interval',mod(latticeFirst%indexMinimum  ,anchorEvery),0)
    call Assert(label//': upper bound is pinned to the anchor interval',mod(latticeFirst%indexMaximum(),anchorEvery),0)
    ! Request a time far outside of the tabulated range, forcing the table to be extended.
    call getTable(solver,timeSecond,table_)
    latticeSecond=table_%lattice
    xValuesSecond=table_%xs     ()
    yValuesSecond=table_%ys     ()
    call Assert(label//': extended table lies on the same lattice',latticeSecond%covers(latticeFirst)      ,.true.)
    call Assert(label//': extended table covers a larger range'   ,latticeSecond%count > latticeFirst%count,.true.)
    ! Extension must not have moved any abscissa, nor changed any previously computed value.
    offset=latticeFirst%indexMinimum-latticeSecond%indexMinimum
    call Assert(                                                                               &
         &      label//': extension preserves previously tabulated abscissae bit-for-bit'    , &
         &      all(xValuesSecond(offset+1:offset+latticeFirst%count  ) == xValuesFirst     ), &
         &      .true.                                                                         &
         &     )
    call Assert(                                                                               &
         &      label//': extension preserves previously tabulated values bit-for-bit'       , &
         &      all(yValuesSecond(offset+1:offset+latticeFirst%count,1) == yValuesFirst(:,1)), &
         &      .true.                                                                         &
         &     )
    ! A table subsequently obtained over the wider range must be identical to the one built by extension.
    call getTable(solver,timeSecond,tableRepeat_)
    xValuesRepeat=tableRepeat_%xs()
    yValuesRepeat=tableRepeat_%ys()
    call Assert(label//': a subsequently obtained table has identical abscissae',all(xValuesRepeat      == xValuesSecond     ),.true.)
    call Assert(label//': a subsequently obtained table has identical values'   ,all(yValuesRepeat(:,1) == yValuesSecond(:,1)),.true.)
    return
  end subroutine testSolver

  subroutine getTable(solver,time,table_)
    !!{RST
    Return the virial density contrast tabulation of the given solver, building or extending it as necessary.
    !!}
    use :: Error, only : Error_Report
    implicit none
    integer                  , intent(in   )               :: solver
    double precision         , intent(in   )               :: time
    class           (table1D), intent(inout), allocatable  :: table_

    select case (solver)
    case (solverCosmologicalConstant)
       call solverCosmologicalConstant_%virialDensityContrast(time,.false.,table_)
    case (solverDarkEnergy           )
       call solverDarkEnergy_          %virialDensityContrast(time,.false.,table_)
    case (solverBaryons              )
       call solverBaryons_             %virialDensityContrast(time,.false.,table_)
    case default
       call Error_Report('unknown solver'//{introspection:location})
    end select
    return
  end subroutine getTable

end program Tests_Spherical_Collapse_Determinism
