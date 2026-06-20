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
Contains a module providing shared utilities for Galacticus benchmark programs: statistics over per-trial timings, machine-parseable output formatting, deterministic random-number seeding, and a checksum sink for defeating dead-code elimination in microbenchmarks.
!!}

module Benchmark_Utilities
  !!{RST
  Shared utilities for benchmark programs.

  All output lines have the form

  .. code-block:: none

     BENCHMARK <suite> <id> "<description>" <mean> <stderr> "<units>"

  where ``<mean>`` and ``<stderr>`` are floating-point and the trailing ``<units>`` is the display unit. The format is whitespace-tokenized, grep/awk-friendly, and consistent across all benchmark programs.
  !!}
  use, intrinsic :: ISO_Fortran_Env, only : output_unit
  use            :: Kind_Numbers   , only : kind_int8
  implicit none
  private
  public :: Benchmark_Report, Benchmark_Seed_RNG, Benchmark_Sink_Add, Benchmark_Sink_Print

  ! Module-level checksum sink. Accumulated via Benchmark_Sink_Add and
  ! observed via Benchmark_Sink_Print, which must be called before the
  ! program exits to prevent the optimizer from eliding the additions.
  double precision :: sink_=0.0d0

contains

  subroutine Benchmark_Report(suite,id,description,trialTime,warmupTrials,innerCount,unitsLabel)
    !!{RST
    Compute the mean and standard error of the per-trial timings in ``trialTime`` (after discarding the first ``warmupTrials`` entries, default ``2``) and write a machine-parseable ``BENCHMARK`` line to standard output.

    Two reporting modes are supported:

    Auto units (``unitsLabel`` absent):
       the displayed value is the raw mean tick count per trial. The unit label is auto-picked from the ``System_Clock`` ``count_rate``: ``ms``, :math:`\mu`\ ``s`` or ``ns`` at :math:`10^3`, :math:`10^6`, or :math:`10^9` ticks/s respectively, falling back to ``ticks`` otherwise. Matches the original convention used by legacy benchmark programs.

    Explicit units (``unitsLabel`` present):
       ticks are converted to nanoseconds using ``count_rate`` and then divided by ``innerCount`` (default ``1``) to yield a per-iteration time. The supplied ``unitsLabel`` (e.g.\ ``"ns/call"``) is used verbatim.

    ``innerCount`` also divides the result in the auto-units branch, in case a caller wants per-iteration timings without an explicit label.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : operator(//), var_str
    use :: String_Handling   , only : operator(//)
    character       (len=*         )                , intent(in   )           :: suite       , id         , description
    integer         (kind=kind_int8), dimension(:)  , intent(in   )           :: trialTime
    integer                                         , intent(in   ), optional :: warmupTrials, innerCount
    character       (len=*         )                , intent(in   ), optional :: unitsLabel
    integer                                                                   :: warmup_     , innerCount_, kept       , &
         &                                                                       totalTrials
    integer         (kind=kind_int8)                                          :: countRate
    double precision                                                          :: meanTicks   , varTicks   , stdErrTicks, &
         &                                                                       mean        , stdErr     , scale
    character       (len=16        )                                          :: units_

    warmup_    =2
    innerCount_=1
    if (present(warmupTrials)) warmup_    =warmupTrials
    if (present(innerCount  )) innerCount_=innerCount

    totalTrials=size(trialTime)
    kept       =totalTrials-warmup_
    if (kept < 2) call Error_Report(var_str('Benchmark_Report needs >= 2 trials after warmup (have ')//totalTrials//', warmup='//warmup_//')'//{introspection:location})

    meanTicks  =      dble(sum(trialTime(warmup_+1:totalTrials)   ))/dble(kept)
    varTicks   =max(  0.0d0,                                                                  &
         &          ( dble(sum(trialTime(warmup_+1:totalTrials)**2))/dble(kept)-meanTicks**2) &
         &           *dble(kept)/dble(kept-1)                                                 &
         &         )
    stdErrTicks=sqrt(varTicks/dble(kept))

    call System_Clock(count_rate=countRate)
    if (countRate <= 0_kind_int8) then
       write (output_unit,'(a)') 'ERROR: System_Clock returned a non-positive tick rate; cannot benchmark.'
       stop 1
    end if

    if (present(unitsLabel)) then
       ! Convert ticks -> ns explicitly, then to per-iteration via innerCount.
       scale =1.0d9/dble(countRate)/dble(innerCount_)
       units_=unitsLabel
    else
       ! Leave the displayed value in raw ticks (per iteration); pick the
       ! human unit from count_rate. This matches legacy behavior.
       scale=1.0d0/dble(innerCount_)
       select case (countRate)
       case (        1000_kind_int8)
          units_='ms'
       case (     1000000_kind_int8)
          units_='μs'
       case (  1000000000_kind_int8)
          units_='ns'
       case default
          units_='ticks'
       end select
    end if

    mean  =meanTicks  *scale
    stdErr=stdErrTicks*scale

    write (output_unit,'(a,1x,a,1x,a,1x,a,1x,f14.3,1x,f14.3,1x,a1,a,a1)') &
         & 'BENCHMARK',trim(suite),trim(id),'"'//trim(description)//'"',  &
         & mean,stdErr,'"',trim(units_),'"'
    return
  end subroutine Benchmark_Report

  subroutine Benchmark_Seed_RNG(seedValue)
    !!{RST
    Seed the intrinsic random number generator deterministically so that benchmarks that draw random inputs are reproducible run-to-run. The optional ``seedValue`` lets callers shift the seed across scenarios if they want decorrelated input streams.
    !!}
    integer, intent(in   ), optional   :: seedValue
    integer                            :: seedSize, seedValue_, i
    integer, allocatable, dimension(:) :: seed

    seedValue_=20260513
    if (present(seedValue)) seedValue_=seedValue
    call random_seed(size=seedSize)
    allocate(seed(seedSize))
    do i=1,seedSize
       seed(i)=seedValue_+i
    end do
    call random_seed(put=seed)
    deallocate(seed)
    return
  end subroutine Benchmark_Seed_RNG

  subroutine Benchmark_Sink_Add(value)
    !!{RST
    Add a value to the module-level checksum sink. Microbenchmarks should pass evaluation results to this routine so that the optimizer cannot dead-code-eliminate the work. Call ``Benchmark_Sink_Print`` once before the program exits to make the accumulated value externally observable.
    !!}
    double precision, intent(in   ) :: value

    sink_=sink_+value
    return
  end subroutine Benchmark_Sink_Add

  subroutine Benchmark_Sink_Print()
    !!{RST
    Print the accumulated checksum sink. Must be called once before the program exits, or the optimizer is free to elide every ``Benchmark_Sink_Add`` call.
    !!}
    
    write (output_unit,'(a,es16.8)') '# Benchmark checksum (do not optimize away): ',sink_
    return
  end subroutine Benchmark_Sink_Print

end module Benchmark_Utilities
