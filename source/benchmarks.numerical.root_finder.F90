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
Contains a program to benchmark the :galacticus-class:`rootFinder` class.

Each scenario constructs a ``rootFinder`` configured for a particular solver type and entry point. A precomputed deterministic array of per-call shifts perturbs the underlying function on every call so the optimizer cannot constant-fold the solver chain and the wrapper's endpoint-value cache sees realistic turnover. An inner loop of many ``find()`` calls per timed trial amortizes ``System_Clock`` resolution.

The scenarios are intentionally aligned with code paths in ``numerical.root_finder.F90``:

``*_bracket``
   hits ``rootFinderFindRange`` (two endpoint evaluations).

``*_bracket_values``
   hits ``rootFinderFindRangeValues``: the caller precomputes both endpoint function values (so ``find`` performs no endpoint evaluations of its own; instead the wrapper's cache serves them on the first two GSL callbacks). The two caller-side ``fQuad`` evaluations sit inside the timed region, just as they do for the ``*_bracket`` scenario (where the same two evaluations happen inside ``find``), so the difference between the two scenarios is the cost of wrapper-cache hits versus the cost of evaluating those same endpoints through ``rootFinderFindRange``'s procedure-pointer call path.

``*_bracket_valueLow``
   hits ``rootFinderFindRangeValueLow`` (one endpoint evaluation).

``*_guess``
   hits ``rootFinderFindGuess`` (the degenerate bracket case where two redundant endpoint evaluations happen) plus the bracket-expansion ``do while`` loop.

``newton``, ``steffenson``
   exercise the derivative-based callbacks (no wrapper cache).

A running checksum of the returned roots is accumulated across every call and printed at the end so the compiler cannot dead-code-eliminate the work.
!!}

program Benchmark_Numerical_Root_Finder
  use, intrinsic :: ISO_Fortran_Env                , only : output_unit
  use            :: Benchmark_Utilities            , only : Benchmark_Report             , Benchmark_Seed_RNG           , &
       &                                                    Benchmark_Sink_Add           , Benchmark_Sink_Print
  use            :: Benchmark_Root_Finder_Functions, only : bmShift                      , fQuad                        , &
       &                                                    fQuadDerivative              , fQuadBoth                    , &
       &                                                    fFar
  use            :: Display                        , only : displayVerbositySet          , verbosityLevelStandard
  use            :: Kind_Numbers                   , only : kind_int8
  use            :: Root_Finder                    , only : rootFinder                                                  , &
       &                                                    gsl_root_fsolver_brent       , gsl_root_fsolver_bisection   , &
       &                                                    gsl_root_fsolver_falsepos    , gsl_root_fdfsolver_newton    , &
       &                                                    gsl_root_fdfsolver_steffenson                               , &
       &                                                    rangeExpandMultiplicative    , rangeExpandAdditive          , &
       &                                                    rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, &
       &                                                    enumerationRangeExpandType   , stoppingCriterionDelta
  implicit none

  ! Benchmark control.
  integer         , parameter                 :: trialCount   =  100      ! Outer trials over which we report mean and standard error.
  integer         , parameter                 :: warmupTrials =    2      ! Trials dropped to absorb cache / branch-predictor / lazy-init warmup.
  integer         , parameter                 :: innerCount   = 1000      ! find() calls per trial — heavier than a single interp() so we use fewer.
  integer         , parameter                 :: shiftCount   = 4096      ! Size of the precomputed shift array (power of 2 for cheap masking).
  integer         , parameter                 :: shiftMask    = shiftCount-1

  double precision, dimension(0:shiftCount-1) :: shifts_

  call displayVerbositySet(verbosityLevelStandard)

  write (output_unit,'(a,i0,a,i0,a,i0,a)') '# Benchmark configuration: trialCount=',trialCount, &
       &                                   ' innerCount=',innerCount,' warmupTrials=',warmupTrials,' (dropped)'

  call buildShifts()

  ! Bracketed entry points — root well inside the bracket, no expansion.
  call benchmarkBracket        (gsl_root_fsolver_brent       ,id="rootFinder_brent_bracket"         ,description="Brent solver, bracket given"                          )
  call benchmarkBracketValues  (gsl_root_fsolver_brent       ,id="rootFinder_brent_bracket_values"  ,description="Brent solver, bracket + both f-values given"          )
  call benchmarkBracketValueLow(gsl_root_fsolver_brent       ,id="rootFinder_brent_bracket_valueLow",description="Brent solver, bracket + f(lower) given"               )
  call benchmarkBracket        (gsl_root_fsolver_bisection   ,id="rootFinder_bisection_bracket"     ,description="Bisection solver, bracket given"                      )
  call benchmarkBracket        (gsl_root_fsolver_falsepos    ,id="rootFinder_falsepos_bracket"      ,description="False-position solver, bracket given"                 )

  ! Guess entry point with bracket expansion (root well outside the initial guess).
  call benchmarkGuessExpand    (rangeExpandMultiplicative    ,id="rootFinder_brent_guess_multExpand",description="Brent solver, single guess + multiplicative expansion")
  call benchmarkGuessExpand    (rangeExpandAdditive          ,id="rootFinder_brent_guess_addExpand" ,description="Brent solver, single guess + additive expansion"      )

  ! Derivative-based solvers (no wrapper cache; exercises gsl_root_fdfsolver path).
  call benchmarkDerivative     (gsl_root_fdfsolver_newton    ,id="rootFinder_newton_guess"          ,description="Newton solver, single guess"                          )
  call benchmarkDerivative     (gsl_root_fdfsolver_steffenson,id="rootFinder_steffenson_guess"      ,description="Steffenson solver, single guess"                      )

  ! Make the accumulated checksum externally observable so the optimizer
  ! cannot elide the benchmark.
  call Benchmark_Sink_Print()

contains

  subroutine buildShifts()
    !!{RST
    Build a deterministic array of per-call shifts in :math:`[-0.4,0.4]` used to perturb the root location on every call.
    !!}
    integer :: i

    call Benchmark_Seed_RNG()
    call random_number(shifts_)
    do i=0,shiftCount-1
       shifts_(i)=0.8d0*shifts_(i)-0.4d0
    end do
    return
  end subroutine buildShifts

  subroutine benchmarkBracket(solverType,id,description)
    !!{RST
    Time ``find(rootRange)`` on a bracket that already contains the root. Each call evaluates the user function at both bracket endpoints (``rootFinderFindRange``) and then iterates the GSL solver to convergence.
    !!}
    integer                         , intent(in   )                      :: solverType
    character       (len=*         ), intent(in   )                      :: id        , description
    type            (rootFinder    )                                     :: finder_
    integer         (kind=kind_int8)                                     :: countStart, countEnd
    integer         (kind=kind_int8)             , dimension(trialCount) :: trialTime
    integer                                                              :: trial     , i
    double precision                                                     :: acc       , root
    double precision                             , dimension(2)          :: bracket

    finder_=rootFinder(                              &
         &             rootFunction     =fQuad     , &
         &             solverType       =solverType, &
         &             toleranceAbsolute=0.0d0     , &
         &             toleranceRelative=1.0d-7      &
         &            )
    bracket=[0.0d0,3.0d0]

    ! Warm: trigger lazy GSL allocation / first-touch caches.
    bmShift=shifts_(0)
    root   =finder_%find(rootRange=bracket)
    call Benchmark_Sink_Add(root)

    do trial=1,trialCount
       acc=0.0d0
       call System_Clock(countStart)
       do i=0,innerCount-1
          bmShift=shifts_(iand(i,shiftMask))
          root   =finder_%find(rootRange=bracket)
          acc    =acc+root
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       call Benchmark_Sink_Add(acc)
    end do
    call Benchmark_Report('numericalRootFinder',id,description,trialTime, &
         &                warmupTrials=warmupTrials,innerCount=innerCount,unitsLabel='ns/call')
    return
  end subroutine benchmarkBracket

  subroutine benchmarkBracketValues(solverType,id,description)
    !!{RST
    Time ``find(rootRange,rootRangeValues)`` where both endpoint function values are precomputed by the caller. Hits the ``rootFinderFindRangeValues`` entry point and avoids any duplicate evaluation of the user function at the bracket ends inside ``find``. The wrapper's cache of endpoint values should also be used by the first two GSL calls.

    The two caller-side ``fQuad`` evaluations are deliberately inside the timed region: they represent the real cost a caller pays to use this entry point. The matching ``*_bracket`` scenario pays the same two evaluations inside ``rootFinderFindRange``, so the comparison between the two scenarios isolates "wrapper cache hit" vs.\ "wrapper cache miss plus procedure-pointer call to ``finderFunction``" rather than introducing or removing user-function calls.
    !!}
    integer                         , intent(in   )                      :: solverType
    character       (len=*         ), intent(in   )                      :: id        , description
    type            (rootFinder    )                                     :: finder_
    integer         (kind=kind_int8)                                     :: countStart, countEnd
    integer         (kind=kind_int8)             , dimension(trialCount) :: trialTime
    integer                                                              :: trial     , i
    double precision                                                     :: acc       , root
    double precision                             , dimension(2)          :: bracket   , bracketValues

    finder_=rootFinder(                              &
         &             rootFunction     =fQuad     , &
         &             solverType       =solverType, &
         &             toleranceAbsolute=0.0d0     , &
         &             toleranceRelative=1.0d-7      &
         &            )
    bracket=[0.0d0,3.0d0]

    bmShift         =shifts_(0)
    bracketValues(1)=fQuad(bracket(1))
    bracketValues(2)=fQuad(bracket(2))
    root            =finder_%find(rootRange=bracket,rootRangeValues=bracketValues)
    call Benchmark_Sink_Add(root)

    do trial=1,trialCount
       acc=0.0d0
       call System_Clock(countStart)
       do i=0,innerCount-1
          bmShift         =shifts_(iand(i,shiftMask))
          bracketValues(1)=fQuad(bracket(1))
          bracketValues(2)=fQuad(bracket(2))
          root            =finder_%find(rootRange=bracket,rootRangeValues=bracketValues)
          acc             =acc+root
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       call Benchmark_Sink_Add(acc)
    end do
    call Benchmark_Report('numericalRootFinder',id,description,trialTime, &
         &                warmupTrials=warmupTrials,innerCount=innerCount,unitsLabel='ns/call')
    return
  end subroutine benchmarkBracketValues

  subroutine benchmarkBracketValueLow(solverType,id,description)
    !!{RST
    Time ``findWithFLower(rootRange,fLower)``. Hits ``rootFinderFindRangeValueLow``, which evaluates the user function at only the upper bracket endpoint.
    !!}
    integer                         , intent(in   )                      :: solverType
    character       (len=*         ), intent(in   )                      :: id        , description
    type            (rootFinder    )                                     :: finder_
    integer         (kind=kind_int8)                                     :: countStart, countEnd
    integer         (kind=kind_int8)             , dimension(trialCount) :: trialTime
    integer                                                              :: trial     , i
    double precision                                                     :: acc       , root       , &
         &                                                                  fLower
    double precision                             , dimension(2)          :: bracket

    finder_=rootFinder(                              &
         &             rootFunction     =fQuad     , &
         &             solverType       =solverType, &
         &             toleranceAbsolute=0.0d0     , &
         &             toleranceRelative=1.0d-7      &
         &            )
    bracket=[0.0d0,3.0d0]

    bmShift=shifts_(0)
    fLower =fQuad(bracket(1))
    root   =finder_%findWithFLower(rootRange=bracket,rootRangeValueLow=fLower)
    call Benchmark_Sink_Add(root)

    do trial=1,trialCount
       acc=0.0d0
       call System_Clock(countStart)
       do i=0,innerCount-1
          bmShift=shifts_(iand(i,shiftMask))
          fLower =fQuad(bracket(1))
          root   =finder_%findWithFLower(rootRange=bracket,rootRangeValueLow=fLower)
          acc    =acc+root
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       call Benchmark_Sink_Add(acc)
    end do
    call Benchmark_Report('numericalRootFinder',id,description,trialTime, &
         &                warmupTrials=warmupTrials,innerCount=innerCount,unitsLabel='ns/call')
    return
  end subroutine benchmarkBracketValueLow

  subroutine benchmarkGuessExpand(expandType,id,description)
    !!{RST
    Time ``find(rootGuess)`` on a function whose root sits well above the initial guess, forcing the bracket-expansion loop in ``rootFinderFind`` to run. With ``rangeExpandMultiplicative`` and ``rangeExpandUpward``\ =2 the upper bracket doubles from 0.5 until it crosses the root near 10 (:math:`\sim`\ 5 expansion steps). With ``rangeExpandAdditive`` and ``rangeExpandUpward``\ =2 the upper bracket grows by 2 each step.
    !!}
    type            (enumerationRangeExpandType), intent(in   )         :: expandType
    character       (len=*                     ), intent(in   )         :: id        , description
    type            (rootFinder                )                        :: finder_
    integer         (kind=kind_int8            )                        :: countStart, countEnd
    integer         (kind=kind_int8            ), dimension(trialCount) :: trialTime
    integer                                                             :: trial     , i
    double precision                                                    :: acc       , root

    finder_=rootFinder(                                                             &
         &             rootFunction                 =fFar                         , &
         &             solverType                   =gsl_root_fsolver_brent       , &
         &             toleranceAbsolute            =0.0d0                        , &
         &             toleranceRelative            =1.0d-7                       , &
         &             rangeExpandType              =expandType                   , &
         &             rangeExpandUpward            =2.0d0                        , &
         &             rangeExpandDownward          =1.0d0                        , &
         &             rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &             rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &            )

    bmShift=shifts_(0)
    root   =finder_%find(rootGuess=0.5d0)
    call Benchmark_Sink_Add(root)

    do trial=1,trialCount
       acc=0.0d0
       call System_Clock(countStart)
       do i=0,innerCount-1
          bmShift=shifts_(iand(i,shiftMask))
          root   =finder_%find(rootGuess=0.5d0)
          acc    =acc+root
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       call Benchmark_Sink_Add(acc)
    end do
    call Benchmark_Report('numericalRootFinder',id,description,trialTime, &
         &                warmupTrials=warmupTrials,innerCount=innerCount,unitsLabel='ns/call')
    return
  end subroutine benchmarkGuessExpand

  subroutine benchmarkDerivative(solverType,id,description)
    !!{RST
    Time ``find(rootGuess)`` with a derivative-based GSL solver. Both the value-only and value+derivative callbacks are exercised; the wrapper does *not* cache endpoint values in this path.
    !!}
    integer                         , intent(in   )                      :: solverType
    character       (len=*         ), intent(in   )                      :: id        , description
    type            (rootFinder    )                                     :: finder_
    integer         (kind=kind_int8)                                     :: countStart, countEnd
    integer         (kind=kind_int8)             , dimension(trialCount) :: trialTime
    integer                                                              :: trial     , i
    double precision                                                     :: acc       , root

    finder_=rootFinder(                                               &
         &             rootFunction          =fQuad                 , &
         &             rootFunctionDerivative=fQuadDerivative       , &
         &             rootFunctionBoth      =fQuadBoth             , &
         &             solverType            =solverType            , &
         &             toleranceAbsolute     =0.0d0                 , &
         &             toleranceRelative     =1.0d-7                , &
         &             stoppingCriterion     =stoppingCriterionDelta  &
         &            )

    bmShift=shifts_(0)
    root   =finder_%find(rootGuess=2.0d0)
    call Benchmark_Sink_Add(root)

    do trial=1,trialCount
       acc=0.0d0
       call System_Clock(countStart)
       do i=0,innerCount-1
          bmShift=shifts_(iand(i,shiftMask))
          root   =finder_%find(rootGuess=2.0d0)
          acc    =acc+root
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       call Benchmark_Sink_Add(acc)
    end do
    call Benchmark_Report('numericalRootFinder',id,description,trialTime, &
         &                warmupTrials=warmupTrials,innerCount=innerCount,unitsLabel='ns/call')
    return
  end subroutine benchmarkDerivative

end program Benchmark_Numerical_Root_Finder
