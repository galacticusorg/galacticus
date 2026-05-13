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
Contains a program to benchmark the numerical interpolation primitives provided by the
\refClass{interpolator} and \refClass{interpolator2D} types.

Each scenario constructs an interpolator, builds a precomputed array of query points
(so we measure interpolation, not RNG), and runs an inner loop of many calls inside
each timed trial. We report time-per-call rather than time-per-trial: at typical clock
resolutions a single \mono{\%interpolate} call is only a handful of ticks, so timing
one call at a time is dominated by noise.

A running checksum is accumulated across every call and printed at the end so the
compiler cannot dead-code-eliminate the interpolation results.
!!}

program Benchmark_Numerical_Interpolation
  use, intrinsic :: ISO_Fortran_Env        , only : output_unit
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Display                , only : displayVerbositySet, verbosityLevelStandard
  use            :: Kind_Numbers           , only : kind_int8
  use            :: Numerical_Interpolation, only : interpolator       , interpolator2D        , &
       &                                            gsl_interp_linear  , gsl_interp_cspline
  implicit none

  ! Benchmark control.
  integer        , parameter :: trialCount    =100      ! Outer trials over which we report mean and standard error.
  integer        , parameter :: warmupTrials  =  2      ! Trials at the start dropped to absorb cache/branch-predictor warmup.
  integer        , parameter :: innerCount    =100000   ! Inner loop iterations per trial (amortizes clock resolution).
  integer        , parameter :: queryCount    =4096     ! Size of the precomputed query-point array (power of 2 for cheap masking).
  integer        , parameter :: queryMask     =queryCount-1

  ! Clock rate (ticks per second). All per-call times are converted explicitly to ns.
  integer        (kind=kind_int8) :: countRate

  ! Accumulated sink (printed at the end so calls cannot be elided).
  double precision :: sink=0.0d0

  call displayVerbositySet(verbosityLevelStandard)

  call System_Clock(count_rate=countRate)
  if (countRate <= 0_kind_int8) then
     write (output_unit,'(a)') 'ERROR: System_Clock returned a non-positive tick rate; cannot benchmark.'
     stop 1
  end if

  write (output_unit,'(a,i0,a,i0,a)') '# Benchmark configuration: trialCount=',trialCount, &
       &                              ' innerCount=',innerCount,' (first 2 trials dropped)'
  write (output_unit,'(a,i0,a)'     ) '# System_Clock rate: ',countRate,' ticks/s (all times reported in ns/call)'

  ! 1D linear interpolations across a range of array sizes and access patterns.
  call benchmarkInterpolate1D(           10,gsl_interp_linear ,sequential=.true. ,id="interp1D_linear_N10_seq"     ,description="1D linear interp, N=10, sequential access"     )
  call benchmarkInterpolate1D(           10,gsl_interp_linear ,sequential=.false.,id="interp1D_linear_N10_rand"    ,description="1D linear interp, N=10, random access"         )
  call benchmarkInterpolate1D(          100,gsl_interp_linear ,sequential=.true. ,id="interp1D_linear_N100_seq"    ,description="1D linear interp, N=100, sequential access"    )
  call benchmarkInterpolate1D(          100,gsl_interp_linear ,sequential=.false.,id="interp1D_linear_N100_rand"   ,description="1D linear interp, N=100, random access"        )
  call benchmarkInterpolate1D(        10000,gsl_interp_linear ,sequential=.true. ,id="interp1D_linear_N10000_seq"  ,description="1D linear interp, N=10000, sequential access"  )
  call benchmarkInterpolate1D(        10000,gsl_interp_linear ,sequential=.false.,id="interp1D_linear_N10000_rand" ,description="1D linear interp, N=10000, random access"      )

  ! 1D cspline path (smaller iteration count: spline eval is heavier per call).
  call benchmarkInterpolate1D(          100,gsl_interp_cspline,sequential=.false.,id="interp1D_cspline_N100_rand"  ,description="1D cspline interp, N=100, random access"       )

  ! Sub-primitives.
  call benchmarkLinearFactors (          100                                     ,id="linearFactors_N100_rand"     ,description="1D linearFactors, N=100, random access"        )
  call benchmarkLocate        (          100                                     ,id="locate_N100_rand"            ,description="1D locate, N=100, random access"               )

  ! 2D bilinear.
  call benchmarkInterpolate2D (           50                                     ,id="interp2D_N50x50_rand"        ,description="2D bilinear interp, 50x50, random access"      )

  ! Print the sink. This MUST be observable to the outside world or the optimizer
  ! is free to elide the entire benchmark.
  write (output_unit,'(a,es16.8)') '# Checksum (do not optimize away): ',sink

contains

  subroutine buildArray(n,x,y)
    !!{
    Build a smooth, monotonically increasing $x$-array and a corresponding $y$-array.
    The arrays span $x\in[0,1]$ with a deliberate non-uniform spacing so that the
    interpolator cannot exploit uniform-grid shortcuts (and so that linear and
    cubic splines have distinguishable behavior).
    !!}
    integer                                              , intent(in   ) :: n
    double precision                , allocatable, dimension(:), intent(  out) :: x, y
    integer                                                                   :: i
    double precision                                                          :: t

    allocate(x(n),y(n))
    do i=1,n
       t   =dble(i-1)/dble(n-1)
       ! Non-uniform but strictly increasing parametrization.
       x(i)=t+0.1d0*t*(1.0d0-t)
       y(i)=sin(6.0d0*t)+0.5d0*cos(2.5d0*t)
    end do
    return
  end subroutine buildArray

  subroutine buildQueries(n,xMin,xMax,sequential,queries)
    !!{
    Build a precomputed array of query points within $[\mathrm{xMin},\mathrm{xMax}]$.
    If \mono{sequential} is true the points walk monotonically through the range,
    cycling back to the start (this is the case where the GSL accelerator should
    hit its cached index). Otherwise the points are uniformly random, which forces
    the accelerator to fall back to a binary search on every call.
    !!}
    integer        , intent(in   )                  :: n
    double precision, intent(in   )                  :: xMin, xMax
    logical        , intent(in   )                  :: sequential
    double precision, intent(  out), dimension(0:n-1) :: queries
    integer                                          :: i
    double precision                                  :: t
    double precision, dimension(:), allocatable      :: r

    if (sequential) then
       do i=0,n-1
          t          =dble(i)/dble(n)
          ! Step monotonically through (xMin,xMax), shrunk slightly so we stay strictly in-range.
          queries(i)=xMin+(xMax-xMin)*(0.001d0+0.998d0*t)
       end do
    else
       allocate(r(n))
       call random_number(r)
       do i=0,n-1
          queries(i)=xMin+(xMax-xMin)*(0.001d0+0.998d0*r(i+1))
       end do
       deallocate(r)
    end if
    return
  end subroutine buildQueries

  subroutine seedRandom()
    !!{
    Seed the intrinsic RNG deterministically so the random-access scenarios are
    reproducible from run to run.
    !!}
    integer              :: seedSize
    integer, allocatable :: seed(:)
    integer              :: i

    call random_seed(size=seedSize)
    allocate(seed(seedSize))
    do i=1,seedSize
       seed(i)=20260513+i
    end do
    call random_seed(put=seed)
    deallocate(seed)
    return
  end subroutine seedRandom

  subroutine reportTiming(id,description,trialTime)
    !!{
    Compute and report the per-call mean time (in ns) and its standard error,
    following the line format used by other Galacticus benchmark programs.
    !!}
    character      (len=*       )  , intent(in   )                       :: id, description
    integer        (kind=kind_int8), intent(in   ), dimension(trialCount) :: trialTime
    integer                                                              :: kept
    double precision                                                     :: meanTicks    , varTicks   , stdErrTicks, &
         &                                                                  ticksToNs    , meanPerCall, stdErrPerCall

    kept         =trialCount-warmupTrials
    meanTicks    =       dble(sum(trialTime(warmupTrials+1:trialCount)   ))/dble(kept)
    varTicks     =max(0.0d0,                                                                                          &
         &            (dble(sum(trialTime(warmupTrials+1:trialCount)**2))/dble(kept)-meanTicks**2)*                    &
         &            dble(kept)/dble(max(kept-1,1))                                                                   &
         &           )
    stdErrTicks  =sqrt(varTicks/dble(kept))
    ! Convert (ticks per trial) -> (ns per single call). countRate is ticks/s, so 1 tick = 1e9/countRate ns;
    ! divide by innerCount to get per-call.
    ticksToNs    =1.0d9/dble(countRate)
    meanPerCall  =meanTicks  *ticksToNs/dble(innerCount)
    stdErrPerCall=stdErrTicks*ticksToNs/dble(innerCount)
    write (output_unit,'(a,1x,a,1x,a,1x,a,1x,f14.3,1x,f14.3,1x,a)') &
         & 'BENCHMARK','numericalInterpolation',trim(id),'"'//trim(description)//'"', &
         & meanPerCall,stdErrPerCall,'"ns/call"'
    return
  end subroutine reportTiming

  subroutine benchmarkInterpolate1D(n,interpType,sequential,id,description)
    integer        , intent(in   )                :: n, interpType
    logical        , intent(in   )                :: sequential
    character      (len=*), intent(in   )         :: id, description
    type            (interpolator   )              :: interp_
    double precision, allocatable, dimension(:)   :: x, y, queries
    integer        (kind=kind_int8)              :: countStart, countEnd
    integer        (kind=kind_int8), dimension(trialCount) :: trialTime
    integer                                       :: trial, i
    double precision                              :: acc, v

    call seedRandom()
    call buildArray  (n,x,y)
    allocate(queries(0:queryCount-1))
    call buildQueries(queryCount,x(1),x(n),sequential,queries)
    interp_=interpolator(x,y,interpolationType=interpType)

    ! Warm the interpolator (in particular trigger lazy initialization for
    ! linear-only deferred-init paths) and prime the accelerator.
    v=interp_%interpolate(queries(0))
    sink=sink+v

    do trial=1,trialCount
       acc=0.0d0
       call System_Clock(countStart)
       do i=0,innerCount-1
          acc=acc+interp_%interpolate(queries(iand(i,queryMask)))
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       sink           =sink+acc                ! Force the optimizer to keep 'acc'.
    end do

    call reportTiming(id,description,trialTime)
    deallocate(x,y,queries)
    return
  end subroutine benchmarkInterpolate1D

  subroutine benchmarkLinearFactors(n,id,description)
    integer        , intent(in   )                :: n
    character      (len=*), intent(in   )         :: id, description
    type            (interpolator   )              :: interp_
    double precision, allocatable, dimension(:)   :: x, y, queries
    integer        (kind=kind_int8)              :: countStart, countEnd
    integer        (kind=kind_int8), dimension(trialCount) :: trialTime
    integer                                       :: trial, i
    integer        (c_size_t      )              :: idx
    double precision, dimension(0:1)              :: h
    double precision                              :: acc

    call seedRandom()
    call buildArray  (n,x,y)
    allocate(queries(0:queryCount-1))
    call buildQueries(queryCount,x(1),x(n),.false.,queries)
    interp_=interpolator(x,y)

    call interp_%linearFactors(queries(0),idx,h)
    sink=sink+h(0)+h(1)+dble(idx)

    do trial=1,trialCount
       acc=0.0d0
       call System_Clock(countStart)
       do i=0,innerCount-1
          call interp_%linearFactors(queries(iand(i,queryMask)),idx,h)
          acc=acc+h(0)+dble(idx)              ! Use the results so the call cannot be elided.
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       sink           =sink+acc
    end do

    call reportTiming(id,description,trialTime)
    deallocate(x,y,queries)
    return
  end subroutine benchmarkLinearFactors

  subroutine benchmarkLocate(n,id,description)
    integer        , intent(in   )                :: n
    character      (len=*), intent(in   )         :: id, description
    type            (interpolator   )              :: interp_
    double precision, allocatable, dimension(:)   :: x, y, queries
    integer        (kind=kind_int8)              :: countStart, countEnd
    integer        (kind=kind_int8), dimension(trialCount) :: trialTime
    integer                                       :: trial, i
    integer        (c_size_t      )              :: idx, accInt

    call seedRandom()
    call buildArray  (n,x,y)
    allocate(queries(0:queryCount-1))
    call buildQueries(queryCount,x(1),x(n),.false.,queries)
    interp_=interpolator(x,y)

    idx=interp_%locate(queries(0))
    sink=sink+dble(idx)

    do trial=1,trialCount
       accInt=0_c_size_t
       call System_Clock(countStart)
       do i=0,innerCount-1
          idx   =interp_%locate(queries(iand(i,queryMask)))
          accInt=accInt+idx
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       sink           =sink+dble(accInt)
    end do

    call reportTiming(id,description,trialTime)
    deallocate(x,y,queries)
    return
  end subroutine benchmarkLocate

  subroutine benchmarkInterpolate2D(n,id,description)
    integer        , intent(in   )                :: n
    character      (len=*), intent(in   )         :: id, description
    type            (interpolator2D )              :: interp_
    double precision, allocatable, dimension(:)   :: x, y, yIgnored, qx, qy
    double precision, allocatable, dimension(:,:) :: z
    integer        (kind=kind_int8)              :: countStart, countEnd
    integer        (kind=kind_int8), dimension(trialCount) :: trialTime
    integer                                       :: trial, i, j
    double precision                              :: acc, v

    call seedRandom()
    ! Build the x and y axes. buildArray writes both an x and a y; we only need its
    ! monotonically increasing first array, so the second is discarded after each call.
    call buildArray(n,x,yIgnored); deallocate(yIgnored)
    call buildArray(n,y,yIgnored); deallocate(yIgnored)
    allocate(z(n,n))
    do j=1,n
       do i=1,n
          z(i,j)=sin(5.0d0*x(i))*cos(3.0d0*y(j))
       end do
    end do
    allocate(qx(0:queryCount-1),qy(0:queryCount-1))
    call buildQueries(queryCount,x(1),x(n),.false.,qx)
    call buildQueries(queryCount,y(1),y(n),.false.,qy)
    interp_=interpolator2D(x,y,z)

    v=interp_%interpolate(qx(0),qy(0))
    sink=sink+v

    do trial=1,trialCount
       acc=0.0d0
       call System_Clock(countStart)
       do i=0,innerCount-1
          acc=acc+interp_%interpolate(qx(iand(i,queryMask)),qy(iand(i,queryMask)))
       end do
       call System_Clock(countEnd)
       trialTime(trial)=countEnd-countStart
       sink           =sink+acc
    end do

    call reportTiming(id,description,trialTime)
    deallocate(x,y,z,qx,qy)
    return
  end subroutine benchmarkInterpolate2D

end program Benchmark_Numerical_Interpolation
