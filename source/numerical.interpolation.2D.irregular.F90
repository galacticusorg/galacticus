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
Contains a module which implements bivariate interpolation on irregularly distributed points using the Akima algorithm
(Algorithm 526, ACM Transactions on Mathematical Software, 1978).
!!}

module Numerical_Interpolation_2D_Irregular
  !!{
  Implements bivariate interpolation on irregularly distributed points (Akima 1978).  All state is encapsulated in
  \reftype{interp2dIrregularObject} so that separate instances are fully independent and thread-safe without any critical
  sections.
  !!}
  implicit none
  private
  public :: Interpolate_2D_Irregular, interp2dIrregularObject

  type :: interp2dIrregularObject
     !!{
     Type encapsulating all state for a 2-D irregular-grid interpolation.  Each instance is self-contained; operations on
     different instances have no shared mutable state and are therefore safe to call from concurrent OpenMP threads.
     !!}
     private
     ! Input data points
     integer                                        :: nData       = 0
     integer                                        :: nNeighbors  = 5
     double precision, allocatable, dimension(:)    :: xData, yData, zData
     ! Triangulation: ipt(3*nTriangles) stores vertex indices CCW; ipl(3*nBorder) stores border segment endpoints + triangle
     integer                                        :: nTriangles  = 0
     integer                                        :: nBorder     = 0
     integer         , allocatable, dimension(:)    :: ipt, ipl
     ! Closest-neighbour indices: ipc(nNeighbors*nData)
     integer         , allocatable, dimension(:)    :: ipc
     ! Partial derivatives: pd(5*nData) stores ZX,ZY,ZXX,ZXY,ZYY per point
     double precision, allocatable, dimension(:)    :: pd
     ! 9-section lookup grid (built once from the triangulation)
     logical                                        :: gridReady   = .false.
     double precision                               :: xs1, xs2, ys1, ys2
     integer         , allocatable, dimension(:)    :: ntsc          ! (9)  triangles per section
     integer         , allocatable, dimension(:)    :: sectionData   ! (9*nTriangles) packed triangle numbers per section
     double precision, allocatable, dimension(:)    :: triBounds     ! (4*nTriangles) xmn,xmx,ymn,ymx per triangle
     ! Cache: last triangle/border-segment code located (replaces module-level itpv/itipv)
     integer                                        :: lastTriangle = 0
     ! Initialisation flag
     logical                                        :: initialized = .false.
  end type interp2dIrregularObject

  interface Interpolate_2D_Irregular
     module procedure Interpolate_2D_Irregular_Array
     module procedure Interpolate_2D_Irregular_Scalar
  end interface Interpolate_2D_Irregular

contains

  function Interpolate_2D_Irregular_Array(dataX,dataY,dataZ,interpolateX,interpolateY,workspace,numberComputePoints,reset) &
       result(zi)
    !!{
    Perform interpolation on a set of points irregularly spaced on a 2D surface.
    !!}
    implicit none
    type            (interp2dIrregularObject)                               , intent(inout)           :: workspace
    double precision                         , dimension(:)                 , intent(in   )           :: dataX        , dataY       , &
         &                                                                                               dataZ        , interpolateX, &
         &                                                                                               interpolateY
    integer                                                                 , intent(in   ), optional :: numberComputePoints
    logical                                                                 , intent(inout), optional :: reset
    double precision                         , dimension(size(interpolateX))                          :: zi
    logical                                                                                           :: doReset
    integer                                                                                           :: i

    ! Determine whether a full reset/reinitialisation is needed.
    if (present(reset)) then
       doReset = reset
       reset   = .false.
    else
       doReset = .true.
    end if

    ! Set number of neighbours for derivative estimation.
    if (present(numberComputePoints)) then
       workspace%nNeighbors = numberComputePoints
    else
       workspace%nNeighbors = 5
    end if

    ! On reset (or first call), rebuild triangulation, neighbours, derivatives, and section grid.
    if (doReset .or. .not.workspace%initialized) then
       call initializeWorkspace(workspace, dataX, dataY, dataZ)
    end if

    ! Locate and interpolate each requested point.
    do i = 1, size(interpolateX)
       zi(i) = interpolateOne(workspace, interpolateX(i), interpolateY(i))
    end do

  end function Interpolate_2D_Irregular_Array

  double precision function Interpolate_2D_Irregular_Scalar(dataX,dataY,dataZ,interpolateX,interpolateY,workspace, &
       &                                                     numberComputePoints,reset)
    !!{
    Scalar wrapper: interpolate at a single point.
    !!}
    implicit none
    type            (interp2dIrregularObject)              , intent(inout)           :: workspace
    double precision                         , dimension(:), intent(in   )           :: dataX, dataY, dataZ
    double precision                                       , intent(in   )           :: interpolateX, interpolateY
    integer                                                , intent(in   ), optional :: numberComputePoints
    logical                                                , intent(inout), optional :: reset
    double precision                         , dimension(1)                          :: xArr, yArr, zArr

    xArr(1) = interpolateX
    yArr(1) = interpolateY
    zArr    = Interpolate_2D_Irregular_Array(dataX,dataY,dataZ,xArr,yArr,workspace,numberComputePoints,reset)
    Interpolate_2D_Irregular_Scalar = zArr(1)
  end function Interpolate_2D_Irregular_Scalar

  ! ── internal: top-level orchestration ────────────────────────────────────────

  subroutine initializeWorkspace(ws, xd, yd, zd)
    !!{
    Build triangulation, find closest neighbours, estimate derivatives, and build the 9-section lookup grid.
    All results are stored in {\normalfont \ttfamily ws}.
    !!}
    implicit none
    type            (interp2dIrregularObject), intent(inout) :: ws
    double precision, dimension(:)           , intent(in   ) :: xd, yd, zd
    integer                                                  :: ndp, ncp

    ndp = size(xd)
    ncp = ws%nNeighbors

    ! Store data.
    ws%nData = ndp
    if (allocated(ws%xData)) deallocate(ws%xData, ws%yData, ws%zData)
    allocate(ws%xData(ndp), ws%yData(ndp), ws%zData(ndp))
    ws%xData = xd;  ws%yData = yd;  ws%zData = zd

    ! Triangulate.
    call buildTriangulation(ws)

    ! Find closest neighbours.
    if (allocated(ws%ipc)) deallocate(ws%ipc)
    allocate(ws%ipc(ncp*ndp))
    call findClosestNeighbors(ndp, xd, yd, ncp, ws%ipc)

    ! Estimate partial derivatives.
    if (allocated(ws%pd)) deallocate(ws%pd)
    allocate(ws%pd(5*ndp))
    call estimateDerivatives(ndp, xd, yd, zd, ncp, ws%ipc, ws%pd)

    ! Build 9-section lookup grid.
    call buildSectionGrid(ws)

    ws%lastTriangle = 0
    ws%initialized  = .true.
  end subroutine initializeWorkspace

  double precision function interpolateOne(ws, xii, yii)
    !!{
    Locate the triangle containing {\normalfont \ttfamily (xii,yii)} and interpolate.
    !!}
    implicit none
    type            (interp2dIrregularObject), intent(inout) :: ws
    double precision                         , intent(in   ) :: xii, yii
    integer                                                  :: iti

    iti             = locatePoint(ws, xii, yii)
    ws%lastTriangle = iti
    interpolateOne  = interpolatePoint(ws, xii, yii, iti)
  end function interpolateOne

  ! ── Phase-2 placeholder: buildTriangulation ──────────────────────────────────
  subroutine buildTriangulation(ws)
    implicit none
    type(interp2dIrregularObject), intent(inout) :: ws
    ! TODO Phase 2
    call Error_Report('buildTriangulation not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end subroutine buildTriangulation

  ! ── Phase-3 placeholder: findClosestNeighbors ────────────────────────────────
  subroutine findClosestNeighbors(ndp, xd, yd, ncp, ipc)
    implicit none
    integer                      , intent(in )               :: ndp, ncp
    double precision, dimension(:), intent(in )               :: xd, yd
    integer         , dimension(:), intent(out)               :: ipc
    ! TODO Phase 3
    call Error_Report('findClosestNeighbors not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end subroutine findClosestNeighbors

  ! ── Phase-3 placeholder: estimateDerivatives ─────────────────────────────────
  subroutine estimateDerivatives(ndp, xd, yd, zd, ncp, ipc, pd)
    implicit none
    integer                      , intent(in )               :: ndp, ncp
    double precision, dimension(:), intent(in )               :: xd, yd, zd
    integer         , dimension(:), intent(in )               :: ipc
    double precision, dimension(:), intent(out)               :: pd
    ! TODO Phase 3
    call Error_Report('estimateDerivatives not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end subroutine estimateDerivatives

  ! ── Phase-4 placeholder: buildSectionGrid ────────────────────────────────────
  subroutine buildSectionGrid(ws)
    implicit none
    type(interp2dIrregularObject), intent(inout) :: ws
    ! TODO Phase 4
    call Error_Report('buildSectionGrid not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end subroutine buildSectionGrid

  ! ── Phase-4 placeholder: locatePoint ─────────────────────────────────────────
  integer function locatePoint(ws, xii, yii)
    implicit none
    type            (interp2dIrregularObject), intent(inout) :: ws
    double precision                         , intent(in   ) :: xii, yii
    locatePoint = 0
    ! TODO Phase 4
    call Error_Report('locatePoint not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end function locatePoint

  ! ── Phase-5 placeholder: interpolatePoint ────────────────────────────────────
  double precision function interpolatePoint(ws, xii, yii, iti)
    implicit none
    type            (interp2dIrregularObject), intent(in) :: ws
    double precision                         , intent(in) :: xii, yii
    integer                                  , intent(in) :: iti
    interpolatePoint = 0.0d0
    ! TODO Phase 5
    call Error_Report('interpolatePoint not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end function interpolatePoint

end module Numerical_Interpolation_2D_Irregular
