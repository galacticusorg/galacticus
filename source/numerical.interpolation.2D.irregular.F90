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
  ! Phase 2: triangulation
  ! ════════════════════════════════════════════════════════════════════════════

  subroutine buildTriangulation(ws)
    !!{
    Triangulate the data points stored in {\normalfont \ttfamily ws}.  Produces a Delaunay-like triangulation using
    the max-min-angle criterion (Lawson).  On return, {\normalfont \ttfamily ws\%ipt}, {\normalfont \ttfamily ws\%ipl},
    {\normalfont \ttfamily ws\%nTriangles} and {\normalfont \ttfamily ws\%nBorder} are set.
    Port of {\normalfont \ttfamily idtang} from the original BIVAR package (Akima 1978).
    !!}
    implicit none
    type(interp2dIrregularObject), intent(inout) :: ws
    integer                                      :: ndp, ndpm1
    integer                                      :: ip, ip1, ip1p1, ip2, ip3
    integer                                      :: ipmn1, ipmn2
    integer                                      :: ipl1, ipl2, iplj1, iplj2
    integer                                      :: ipti, ipti1, ipti2, ipt1, ipt2, ipt3
    integer                                      :: it, it1t3, it2t3, itt3, itt3r
    integer                                      :: itf(2), ntf
    integer                                      :: irep, ilf, ilft2, jlt3
    integer                                      :: jp, jp1, jp2, jp2t3, jp3t3
    integer                                      :: jpc, jpmn, jpmx
    integer                                      :: jwl, jwl1, jwl1mn, nlft2
    integer                                      :: nl0, nln, nlf, nlfc, nlt3, nlnt3
    integer                                      :: nt0, ntt3, ntt3p3
    integer                                      :: nsh, nsht3
    integer         , allocatable, dimension(:)  :: iwl, iwp
    double precision, allocatable, dimension(:)  :: wk
    double precision                             :: ar, armn, armx
    double precision                             :: dsq12, dsqi, dsqmn, dsqmx
    double precision                             :: dx, dx21, dxmn, dxmx
    double precision                             :: dy, dy21, dymn, dymx
    double precision                             :: x1, xdmp, y1, ydmp
    integer         , parameter                  :: nrep  = 100
    double precision, parameter                  :: ratio = 1.0d-6

    ndp   = ws%nData
    ndpm1 = ndp - 1
    if (ndp < 4) then
       write (*,'(a,i0)') '  buildTriangulation: fewer than 4 data points, ndp = ', ndp
       ws%nTriangles = 0
       return
    end if

    ! Allocate temporary workspace and output arrays.
    allocate(iwl(18*ndp), iwp(ndp), wk(ndp))
    if (allocated(ws%ipt)) deallocate(ws%ipt)
    if (allocated(ws%ipl)) deallocate(ws%ipl)
    allocate(ws%ipt(max(6*ndp-15,9)), ws%ipl(6*ndp))

    ! ── Find closest pair of data points ──────────────────────────────────────
    dsqmn = dsqf2(ws%xData(1),ws%yData(1),ws%xData(2),ws%yData(2))
    ipmn1 = 1
    ipmn2 = 2
    do ip1 = 1, ndpm1
       x1    = ws%xData(ip1)
       y1    = ws%yData(ip1)
       ip1p1 = ip1 + 1
       do ip2 = ip1p1, ndp
          dsqi = dsqf2(x1,y1,ws%xData(ip2),ws%yData(ip2))
          if (dsqi == 0.0d0) then
             write (*,'(a,2(1x,i0))') '  buildTriangulation: identical data points', ip1, ip2
             ws%nTriangles = 0
             deallocate(iwl, iwp, wk)
             return
          end if
          if (dsqi < dsqmn) then
             dsqmn = dsqi;  ipmn1 = ip1;  ipmn2 = ip2
          end if
       end do
    end do
    dsq12 = dsqmn
    xdmp  = (ws%xData(ipmn1) + ws%xData(ipmn2)) / 2.0d0
    ydmp  = (ws%yData(ipmn1) + ws%yData(ipmn2)) / 2.0d0

    ! ── Sort remaining points by distance from midpoint ───────────────────────
    iwp(1) = ipmn1;  iwp(2) = ipmn2
    jp1    = 2
    do ip1 = 1, ndp
       if (ip1 == ipmn1 .or. ip1 == ipmn2) cycle
       jp1      = jp1 + 1
       iwp(jp1) = ip1
       wk(jp1)  = dsqf2(xdmp,ydmp,ws%xData(ip1),ws%yData(ip1))
    end do
    do jp1 = 3, ndpm1
       dsqmn = wk(jp1);  jpmn = jp1
       do jp2 = jp1+1, ndp
          if (wk(jp2) < dsqmn) then
             dsqmn = wk(jp2);  jpmn = jp2
          end if
       end do
       if (jpmn /= jp1) then
          it       = iwp(jp1);  iwp(jp1) = iwp(jpmn);  iwp(jpmn) = it
          wk(jpmn) = wk(jp1)
       end if
    end do

    ! ── Ensure first 3 points are not collinear ───────────────────────────────
    ar   = dsq12 * ratio
    x1   = ws%xData(ipmn1);  y1 = ws%yData(ipmn1)
    dx21 = ws%xData(ipmn2) - x1;  dy21 = ws%yData(ipmn2) - y1
    jp   = 3
    do while (jp <= ndp)
       ip = iwp(jp)
       if (abs((ws%yData(ip)-y1)*dx21 - (ws%xData(ip)-x1)*dy21) > ar) exit
       jp = jp + 1
    end do
    if (jp > ndp) then
       write (*,'(a,i0)') '  buildTriangulation: all collinear data points, ndp = ', ndp
       ws%nTriangles = 0
       deallocate(iwl, iwp, wk)
       return
    end if
    if (jp /= 3) then
       ip = iwp(jp)
       do jpc = jp, 4, -1
          iwp(jpc) = iwp(jpc-1)
       end do
       iwp(3) = ip
    end if

    ! ── Form the first triangle ───────────────────────────────────────────────
    ip1 = ipmn1;  ip2 = ipmn2;  ip3 = iwp(3)
    if (side2(ws%xData(ip1),ws%yData(ip1),ws%xData(ip2),ws%yData(ip2), &
              ws%xData(ip3),ws%yData(ip3)) < 0.0d0) then
       ip1 = ipmn2;  ip2 = ipmn1
    end if
    nt0  = 1;  ntt3 = 3
    ws%ipt(1) = ip1;  ws%ipt(2) = ip2;  ws%ipt(3) = ip3
    nl0  = 3;  nlt3 = 9
    ws%ipl(1) = ip1;  ws%ipl(2) = ip2;  ws%ipl(3) = 1
    ws%ipl(4) = ip2;  ws%ipl(5) = ip3;  ws%ipl(6) = 1
    ws%ipl(7) = ip3;  ws%ipl(8) = ip1;  ws%ipl(9) = 1

    ! ── Add remaining (ndp-3) points one at a time ───────────────────────────
    do jp1 = 4, ndp
       ip1 = iwp(jp1)
       x1  = ws%xData(ip1);  y1 = ws%yData(ip1)

       ! Find visible border line segments (rightmost jpmn, leftmost jpmx).
       ip2  = ws%ipl(1)
       jpmn = 1
       dxmn = ws%xData(ip2) - x1;  dymn = ws%yData(ip2) - y1
       dsqmn = dxmn**2 + dymn**2;  armn = dsqmn * ratio
       jpmx  = 1
       dxmx  = dxmn;  dymx = dymn;  dsqmx = dsqmn;  armx = armn
       do jp2 = 2, nl0
          ip2  = ws%ipl(3*jp2-2)
          dx   = ws%xData(ip2) - x1;  dy = ws%yData(ip2) - y1
          dsqi = dx**2 + dy**2
          ar   = dy*dxmn - dx*dymn
          if (ar <= armn) then
             if (ar < -armn .or. dsqi < dsqmn) then
                jpmn = jp2;  dxmn = dx;  dymn = dy
                dsqmn = dsqi;  armn = dsqmn * ratio
             end if
          end if
          ar = dy*dxmx - dx*dymx
          if (ar >= -armx) then
             if (ar > armx .or. dsqi < dsqmx) then
                jpmx = jp2;  dxmx = dx;  dymx = dy
                dsqmx = dsqi;  armx = dsqmx * ratio
             end if
          end if
       end do
       if (jpmx < jpmn) jpmx = jpmx + nl0

       ! Rotate ipl so invisible segments (1..jpmn-1) come first.
       nsh = jpmn - 1
       if (nsh > 0) then
          nsht3 = nsh * 3
          do jp2t3 = 3, nsht3, 3
             jp3t3 = jp2t3 + nlt3
             ws%ipl(jp3t3-2) = ws%ipl(jp2t3-2)
             ws%ipl(jp3t3-1) = ws%ipl(jp2t3-1)
             ws%ipl(jp3t3)   = ws%ipl(jp2t3)
          end do
          do jp2t3 = 3, nlt3, 3
             jp3t3 = jp2t3 + nsht3
             ws%ipl(jp2t3-2) = ws%ipl(jp3t3-2)
             ws%ipl(jp2t3-1) = ws%ipl(jp3t3-1)
             ws%ipl(jp2t3)   = ws%ipl(jp3t3)
          end do
          jpmx = jpmx - nsh
       end if

       ! Add triangles for each visible segment; update border; flag edges.
       jwl   = 0
       nln   = nl0    ! will be updated when jp2 == nl0
       nlnt3 = nlt3
       do jp2 = jpmx, nl0
          jp2t3 = jp2 * 3
          ipl1  = ws%ipl(jp2t3-2);  ipl2 = ws%ipl(jp2t3-1);  it = ws%ipl(jp2t3)
          nt0   = nt0 + 1;  ntt3 = ntt3 + 3
          ws%ipt(ntt3-2) = ipl2;  ws%ipt(ntt3-1) = ipl1;  ws%ipt(ntt3) = ip1
          if (jp2 == jpmx) then
             ws%ipl(jp2t3-1) = ip1;  ws%ipl(jp2t3) = nt0
          end if
          if (jp2 == nl0) then
             nln   = jpmx + 1;  nlnt3 = nln * 3
             ws%ipl(nlnt3-2) = ip1;  ws%ipl(nlnt3-1) = ws%ipl(1);  ws%ipl(nlnt3) = nt0
          end if
          ! Find the vertex of triangle 'it' not on segment (ipl1,ipl2).
          itt3 = it * 3
          ipti = ws%ipt(itt3-2)
          if (ipti == ipl1 .or. ipti == ipl2) then
             ipti = ws%ipt(itt3-1)
             if (ipti == ipl1 .or. ipti == ipl2) ipti = ws%ipt(itt3)
          end if
          if (triangleSwapCheck(ws%xData, ws%yData, ip1, ipti, ipl1, ipl2) /= 0) then
             ws%ipt(itt3-2) = ipti;  ws%ipt(itt3-1) = ipl1;  ws%ipt(itt3)   = ip1
             ws%ipt(ntt3-1) = ipti
             if (jp2 == jpmx) ws%ipl(jp2t3) = it
             if (jp2 == nl0 .and. ws%ipl(3) == it) ws%ipl(3) = nt0
             jwl = jwl + 4
             iwl(jwl-3) = ipl1;  iwl(jwl-2) = ipti
             iwl(jwl-1) = ipti;  iwl(jwl)   = ipl2
          end if
       end do
       nl0  = nln;  nlt3 = nlnt3

       ! Improve triangulation: swap flagged edges up to nrep times.
       nlf = jwl / 2
       if (nlf == 0) cycle
       ntt3p3 = ntt3 + 3
       do irep = 1, nrep
          do ilf = 1, nlf
             ilft2 = ilf * 2
             ipl1  = iwl(ilft2-1);  ipl2 = iwl(ilft2)
             ntf   = 0
             do itt3r = 3, ntt3, 3
                itt3 = ntt3p3 - itt3r
                ipt1 = ws%ipt(itt3-2);  ipt2 = ws%ipt(itt3-1);  ipt3 = ws%ipt(itt3)
                if (ipl1/=ipt1 .and. ipl1/=ipt2 .and. ipl1/=ipt3) cycle
                if (ipl2/=ipt1 .and. ipl2/=ipt2 .and. ipl2/=ipt3) cycle
                ntf      = ntf + 1;  itf(ntf) = itt3 / 3
                if (ntf == 2) exit
             end do
             if (ntf < 2) cycle
             it1t3 = itf(1) * 3
             ipti1 = ws%ipt(it1t3-2)
             if (ipti1==ipl1 .or. ipti1==ipl2) then
                ipti1 = ws%ipt(it1t3-1)
                if (ipti1==ipl1 .or. ipti1==ipl2) ipti1 = ws%ipt(it1t3)
             end if
             it2t3 = itf(2) * 3
             ipti2 = ws%ipt(it2t3-2)
             if (ipti2==ipl1 .or. ipti2==ipl2) then
                ipti2 = ws%ipt(it2t3-1)
                if (ipti2==ipl1 .or. ipti2==ipl2) ipti2 = ws%ipt(it2t3)
             end if
             if (triangleSwapCheck(ws%xData, ws%yData, ipti1, ipti2, ipl1, ipl2) == 0) cycle
             ws%ipt(it1t3-2) = ipti1;  ws%ipt(it1t3-1) = ipti2;  ws%ipt(it1t3)   = ipl1
             ws%ipt(it2t3-2) = ipti2;  ws%ipt(it2t3-1) = ipti1;  ws%ipt(it2t3)   = ipl2
             jwl = jwl + 8
             iwl(jwl-7) = ipl1;  iwl(jwl-6) = ipti1
             iwl(jwl-5) = ipti1; iwl(jwl-4) = ipl2
             iwl(jwl-3) = ipl2;  iwl(jwl-2) = ipti2
             iwl(jwl-1) = ipti2; iwl(jwl)   = ipl1
             do jlt3 = 3, nlt3, 3
                iplj1 = ws%ipl(jlt3-2);  iplj2 = ws%ipl(jlt3-1)
                if ((iplj1==ipl1.and.iplj2==ipti2).or.(iplj2==ipl1.and.iplj1==ipti2)) &
                     ws%ipl(jlt3) = itf(1)
                if ((iplj1==ipl2.and.iplj2==ipti1).or.(iplj2==ipl2.and.iplj1==ipti1)) &
                     ws%ipl(jlt3) = itf(2)
             end do
          end do
          nlfc = nlf;  nlf = jwl / 2
          if (nlf == nlfc) exit
          ! Compact iwl: discard already-processed edges, keep new ones.
          jwl    = 0;  jwl1mn = (nlfc+1)*2;  nlft2 = nlf*2
          do jwl1 = jwl1mn, nlft2, 2
             jwl = jwl + 2;  iwl(jwl-1) = iwl(jwl1-1);  iwl(jwl) = iwl(jwl1)
          end do
          nlf = jwl / 2
       end do
    end do   ! jp1 loop over remaining points

    ! ── Rearrange triangles so each is listed counter-clockwise ───────────────
    do itt3 = 3, ntt3, 3
       ip1 = ws%ipt(itt3-2);  ip2 = ws%ipt(itt3-1);  ip3 = ws%ipt(itt3)
       if (side2(ws%xData(ip1),ws%yData(ip1),ws%xData(ip2),ws%yData(ip2), &
                 ws%xData(ip3),ws%yData(ip3)) < 0.0d0) then
          ws%ipt(itt3-2) = ip2;  ws%ipt(itt3-1) = ip1
       end if
    end do

    ws%nTriangles = nt0
    ws%nBorder    = nl0
    deallocate(iwl, iwp, wk)

  contains

    double precision function dsqf2(u1,v1,u2,v2)
      double precision, intent(in) :: u1, v1, u2, v2
      dsqf2 = (u2-u1)**2 + (v2-v1)**2
    end function dsqf2

    double precision function side2(u1,v1,u2,v2,u3,v3)
      double precision, intent(in) :: u1, v1, u2, v2, u3, v3
      side2 = (v3-v1)*(u2-u1) - (u3-u1)*(v2-v1)
    end function side2

  end subroutine buildTriangulation

  integer function triangleSwapCheck(x, y, i1, i2, i3, i4)
    !!{
    Determine whether swapping the shared diagonal of the quadrilateral formed by points {\normalfont \ttfamily i1}--{\normalfont
    \ttfamily i4} improves the minimum angle (Lawson max-min-angle criterion).  Returns 1 if a swap is recommended, 0 otherwise.
    Port of {\normalfont \ttfamily idxchg} from the original BIVAR package (Akima 1978).
    !!}
    implicit none
    double precision, dimension(:), intent(in) :: x, y
    integer                       , intent(in) :: i1, i2, i3, i4
    double precision :: x1, y1, x2, y2, x3, y3, x4, y4
    double precision :: u1, u2, u3, u4
    double precision :: a1sq, b1sq, c1sq, a2sq, b2sq, c2sq
    double precision :: a3sq, b3sq, c3sq, a4sq, b4sq, c4sq
    double precision :: s1sq, s2sq, s3sq, s4sq

    x1 = x(i1);  y1 = y(i1);  x2 = x(i2);  y2 = y(i2)
    x3 = x(i3);  y3 = y(i3);  x4 = x(i4);  y4 = y(i4)

    triangleSwapCheck = 0
    u3 = (y2-y3)*(x1-x3) - (x2-x3)*(y1-y3)
    u4 = (y1-y4)*(x2-x4) - (x1-x4)*(y2-y4)
    if (.not. (u3*u4 > 0.0d0)) return

    u1   = (y3-y1)*(x4-x1) - (x3-x1)*(y4-y1)
    u2   = (y4-y2)*(x3-x2) - (x4-x2)*(y3-y2)
    a1sq = (x1-x3)**2+(y1-y3)**2;  b3sq = a1sq
    b1sq = (x4-x1)**2+(y4-y1)**2;  a4sq = b1sq
    c1sq = (x3-x4)**2+(y3-y4)**2;  c2sq = c1sq
    a2sq = (x2-x4)**2+(y2-y4)**2;  b4sq = a2sq
    b2sq = (x3-x2)**2+(y3-y2)**2;  a3sq = b2sq
    c3sq = (x2-x1)**2+(y2-y1)**2;  c4sq = c3sq
    s1sq = u1*u1 / (c1sq*max(a1sq,b1sq))
    s2sq = u2*u2 / (c2sq*max(a2sq,b2sq))
    s3sq = u3*u3 / (c3sq*max(a3sq,b3sq))
    s4sq = u4*u4 / (c4sq*max(a4sq,b4sq))
    if (min(s1sq,s2sq) < min(s3sq,s4sq)) triangleSwapCheck = 1
  end function triangleSwapCheck

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
