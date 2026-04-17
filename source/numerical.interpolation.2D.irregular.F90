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
     ! ── legacy workspace (used by idbvip; removed in Phase 6) ─────────────────
     integer         , allocatable, dimension(:) :: integerWork
     double precision, allocatable, dimension(:) :: realWork
     ! ── input data points ─────────────────────────────────────────────────────
     integer                                     :: nData      = 0
     integer                                     :: nNeighbors = 5
     double precision, allocatable, dimension(:) :: xData, yData, zData
     ! ── triangulation ─────────────────────────────────────────────────────────
     ! ipt(3*nTriangles): vertex indices listed CCW
     ! ipl(3*nBorder):    border segment endpoints + owning triangle number
     integer                                     :: nTriangles = 0
     integer                                     :: nBorder    = 0
     integer         , allocatable, dimension(:) :: ipt, ipl
     ! ── closest-neighbour indices: ipc(nNeighbors*nData) ─────────────────────
     integer         , allocatable, dimension(:) :: ipc
     ! ── partial derivatives: pd(5*nData) = ZX,ZY,ZXX,ZXY,ZYY per point ──────
     double precision, allocatable, dimension(:) :: pd
     ! ── 9-section lookup grid ─────────────────────────────────────────────────
     logical                                     :: gridReady  = .false.
     double precision                            :: xs1 = 0.0d0, xs2 = 0.0d0
     double precision                            :: ys1 = 0.0d0, ys2 = 0.0d0
     integer         , allocatable, dimension(:) :: ntsc        ! (9)
     integer         , allocatable, dimension(:) :: sectionData ! (9*nTriangles) packed per section
     double precision, allocatable, dimension(:) :: triBounds   ! (4*nTriangles) xmn,xmx,ymn,ymx
     ! ── inter-call cache (replaces module-level itpv/itipv) ──────────────────
     integer                                     :: lastTriangle = 0
     ! ── initialisation flag ───────────────────────────────────────────────────
     logical                                     :: initialized = .false.
  end type interp2dIrregularObject

  interface Interpolate_2D_Irregular
     module procedure Interpolate_2D_Irregular_Array
     module procedure Interpolate_2D_Irregular_Scalar
  end interface Interpolate_2D_Irregular

contains

  ! ════════════════════════════════════════════════════════════════════════════
  ! Public API  (Phase 6 will replace the idbvip call with the new path)
  ! ════════════════════════════════════════════════════════════════════════════

  function Interpolate_2D_Irregular_Array(dataX,dataY,dataZ,interpolateX,interpolateY,workspace,numberComputePoints,reset) &
       result(zi)
    !!{
    Perform interpolation on a set of points irregularly spaced on a 2D surface.
    !!}
    use :: Bivar, only : idbvip
    implicit none
    type            (interp2dIrregularObject)                               , intent(inout)           :: workspace
    double precision                         , dimension(:)                 , intent(in   )           :: dataX        , dataY       , &
         &                                                                                               dataZ        , interpolateX, &
         &                                                                                               interpolateY
    integer                                                                 , intent(in   ), optional :: numberComputePoints
    logical                                                                 , intent(inout), optional :: reset
    double precision                         , dimension(size(interpolateX))                          :: zi
    integer                                                                                           :: dataPointCount        , integerWorkspaceSize, &
         &                                                                                               interpolatedPointCount , numberComputePointsActual, &
         &                                                                                               realWorkspaceSize      , resetFlag
    logical                                                                                           :: resetActual

    ! Determine reset status.
    if (present(reset)) then
       resetActual = reset
       reset       = .false.
    else
       resetActual = .true.
    end if
    if (resetActual) then
       resetFlag = 1
    else
       resetFlag = 2
    end if

    ! Decide how many points to use for computing partial derivatives.
    if (present(numberComputePoints)) then
       numberComputePointsActual = numberComputePoints
    else
       numberComputePointsActual = 5
    end if

    dataPointCount         = size(dataX)
    interpolatedPointCount = size(interpolateX)

    ! Ensure legacy workspace is large enough.
    integerWorkspaceSize = max(31,27+numberComputePointsActual)*dataPointCount+interpolatedPointCount
    realWorkspaceSize    = 8*dataPointCount
    if (allocated(workspace%integerWork)) then
       if (size(workspace%integerWork) < integerWorkspaceSize) then
          deallocate(workspace%integerWork)
          allocate  (workspace%integerWork(integerWorkspaceSize))
          workspace%integerWork = 0
       end if
       if (size(workspace%realWork) < realWorkspaceSize) then
          deallocate(workspace%realWork)
          allocate  (workspace%realWork(realWorkspaceSize))
          workspace%realWork = 0.0d0
       end if
    else
       allocate(workspace%integerWork(integerWorkspaceSize))
       allocate(workspace%realWork   (realWorkspaceSize   ))
       workspace%integerWork = 0
       workspace%realWork    = 0.0d0
    end if

    ! TODO Phase 6: replace the idbvip call below with the new path (initializeWorkspace +
    ! interpolateOne loop) and remove the use :: Bivar above and the legacy workspace fields.
    !$omp critical (idbvip)
    call idbvip(resetFlag,numberComputePointsActual,dataPointCount,dataX,dataY,dataZ, &
         &      interpolatedPointCount,interpolateX,interpolateY,zi,                   &
         &      workspace%integerWork,workspace%realWork)
    !$omp end critical (idbvip)

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

  ! ════════════════════════════════════════════════════════════════════════════
  ! Internal orchestration (called from Phase 6 onwards)
  ! ════════════════════════════════════════════════════════════════════════════

  subroutine initializeWorkspace(ws, xd, yd, zd)
    !!{
    Build triangulation, find closest neighbours, estimate derivatives, and build the 9-section lookup grid.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (interp2dIrregularObject), intent(inout) :: ws
    double precision, dimension(:)           , intent(in   ) :: xd, yd, zd
    integer                                                  :: ndp, ncp

    ndp = size(xd)
    ncp = ws%nNeighbors

    ws%nData = ndp
    if (allocated(ws%xData)) deallocate(ws%xData, ws%yData, ws%zData)
    allocate(ws%xData(ndp), ws%yData(ndp), ws%zData(ndp))
    ws%xData = xd;  ws%yData = yd;  ws%zData = zd

    call buildTriangulation(ws)

    if (allocated(ws%ipc)) deallocate(ws%ipc)
    allocate(ws%ipc(ncp*ndp))
    call findClosestNeighbors(ndp, xd, yd, ncp, ws%ipc)

    if (allocated(ws%pd)) deallocate(ws%pd)
    allocate(ws%pd(5*ndp))
    call estimateDerivatives(ndp, xd, yd, zd, ncp, ws%ipc, ws%pd)

    call buildSectionGrid(ws)

    ws%lastTriangle = 0
    ws%initialized  = .true.
  end subroutine initializeWorkspace

  double precision function interpolateOne(ws, xii, yii)
    !!{
    Locate the triangle containing {\normalfont \ttfamily (xii,yii)} and interpolate.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (interp2dIrregularObject), intent(inout) :: ws
    double precision                         , intent(in   ) :: xii, yii
    integer                                                  :: iti

    iti             = locatePoint(ws, xii, yii)
    ws%lastTriangle = iti
    interpolateOne  = interpolatePoint(ws, xii, yii, iti)
  end function interpolateOne

  ! ════════════════════════════════════════════════════════════════════════════
  ! Phase 2 stub: triangulation
  ! ════════════════════════════════════════════════════════════════════════════

  subroutine buildTriangulation(ws)
    use :: Error, only : Error_Report
    implicit none
    type(interp2dIrregularObject), intent(inout) :: ws
    ! TODO Phase 2: port idtang + idxchg
    call Error_Report('buildTriangulation not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end subroutine buildTriangulation

  ! ════════════════════════════════════════════════════════════════════════════
  ! Phase 3 stubs: closest neighbours + partial derivatives
  ! ════════════════════════════════════════════════════════════════════════════

  subroutine findClosestNeighbors(ndp, xd, yd, ncp, ipc)
    use :: Error, only : Error_Report
    implicit none
    integer                       , intent(in ) :: ndp, ncp
    double precision, dimension(:), intent(in ) :: xd, yd
    integer         , dimension(:), intent(out) :: ipc
    ipc = 0
    ! TODO Phase 3: port idcldp
    call Error_Report('findClosestNeighbors not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end subroutine findClosestNeighbors

  subroutine estimateDerivatives(ndp, xd, yd, zd, ncp, ipc, pd)
    use :: Error, only : Error_Report
    implicit none
    integer                       , intent(in ) :: ndp, ncp
    double precision, dimension(:), intent(in ) :: xd, yd, zd
    integer         , dimension(:), intent(in ) :: ipc
    double precision, dimension(:), intent(out) :: pd
    pd = 0.0d0
    ! TODO Phase 3: port idpdrv
    call Error_Report('estimateDerivatives not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end subroutine estimateDerivatives

  ! ════════════════════════════════════════════════════════════════════════════
  ! Phase 4 stubs: 9-section grid + point location
  ! ════════════════════════════════════════════════════════════════════════════

  subroutine buildSectionGrid(ws)
    use :: Error, only : Error_Report
    implicit none
    type(interp2dIrregularObject), intent(inout) :: ws
    ! TODO Phase 4: extracted from idlctn initialisation block
    call Error_Report('buildSectionGrid not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end subroutine buildSectionGrid

  integer function locatePoint(ws, xii, yii)
    use :: Error, only : Error_Report
    implicit none
    type            (interp2dIrregularObject), intent(inout) :: ws
    double precision                         , intent(in   ) :: xii, yii
    locatePoint = 0
    ! TODO Phase 4: port idlctn lookup
    call Error_Report('locatePoint not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end function locatePoint

  ! ════════════════════════════════════════════════════════════════════════════
  ! Phase 5 stub: interpolation at a located point
  ! ════════════════════════════════════════════════════════════════════════════

  double precision function interpolatePoint(ws, xii, yii, iti)
    use :: Error, only : Error_Report
    implicit none
    type            (interp2dIrregularObject), intent(in) :: ws
    double precision                         , intent(in) :: xii, yii
    integer                                  , intent(in) :: iti
    interpolatePoint = 0.0d0
    ! TODO Phase 5: port idptip
    call Error_Report('interpolatePoint not yet implemented'//char(0),'numerical.interpolation.2D.irregular')
  end function interpolatePoint

end module Numerical_Interpolation_2D_Irregular
