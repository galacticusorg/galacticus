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
Contains a program to test the numerical range making code.
!!}

program Test_Make_Ranges
  !!{RST
  Tests that numerical range making code works correctly.
  !!}
  use :: Array_Utilities , only : Array_Reverse
  use :: Display         , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Ranges, only : Make_Range         , rangeTypeLinear       , rangeTypeLogarithmic, Range_Pinned        , &
       &                          rangeLattice       , gridSchemePerDecade   , gridSchemePerOctave , gridSchemePerUnit   , &
       &                          Range_Lattice_Offset, Range_Lattice_Extend
  use :: Unit_Tests      , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                               , dimension(1:11) :: range1
  double precision                               , dimension(0:10) :: range2
  type            (rangeLattice)                                   :: lattice     , latticeNarrow, &
       &                                                              latticeWide , latticeUnion
  double precision              , allocatable    , dimension(:   ) :: values      , valuesNarrow , &
       &                                                              valuesWide
  integer                                                          :: offset
  double precision              , allocatable    , dimension(:   ) :: valuesRaw
  logical                       , allocatable    , dimension(:   ) :: isComputedRaw

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical ranges")

  ! Create a linear range.
  range1=Make_Range(1.0d0,3.0d0,size(range1),rangeType=rangeTypeLinear)
  call Assert("linear range creation",range1,[1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0])

  ! Create a linear range in an offset array.
  range2=Make_Range(1.0d0,3.0d0,size(range2),rangeType=rangeTypeLinear)
  call Assert("linear range creation (offset array)",range2,[1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0])

  ! Create a logarithmic range.
  range1=Make_Range(10.0d0,1000.0d0,size(range1),rangeType=rangeTypeLogarithmic)
  call Assert("logarithmic range creation",range1,10.0d0**([1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0]),relTol=1.0d-6)

  ! Create a reverse order range.
  range1=Make_Range(3.0d0,1.0d0,size(range1),rangeType=rangeTypeLinear)
  call Assert("reverse linear range creation",range1,Array_Reverse([1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0]))

  call Unit_Tests_End_Group()

  ! Pinned ranges.
  call Unit_Tests_Begin_Group("Pinned ranges")

  ! Margin-then-pin, points per decade: a request at x=3 with the default factor-of-2 margin spans [1.5,6], which pins outward to
  ! [10⁰,10¹].
  lattice=Range_Pinned(3.0d0,10,gridSchemePerDecade)
  call Assert("per-decade: margin then pin to whole decades"   ,[lattice%minimum(),lattice%maximum()],[1.0d0,10.0d0])
  call Assert("per-decade: exact point count"                  , lattice%count                       ,11            )
  call Assert("per-decade: lattice index of the lower bound"   , lattice%indexMinimum                , 0            )

  ! Pinning must be applied *after* the margin, so that the requested value is never left sitting at the table edge. A request at
  ! x=9.99 must not pin to [10⁰,10¹].
  lattice=Range_Pinned(9.99d0,10,gridSchemePerDecade)
  call Assert("per-decade: request near a decade boundary is not left at the table edge",[lattice%minimum(),lattice%maximum()],[1.0d0,100.0d0])

  ! Exact powers of ten must not be rounded down by a whole decade by round-off in the logarithm (log10(1000) evaluates to
  ! 2.9999999999999996 in double precision). Here a unit margin is used so that the bounds are exactly 10³ and 10⁵.
  lattice=Range_Pinned([1.0d3,1.0d5],10,gridSchemePerDecade,marginFactor=1.0d0)
  call Assert("per-decade: exact powers of ten are not shifted by round-off",[lattice%minimum(),lattice%maximum()],[1.0d3,1.0d5])
  call Assert("per-decade: exact point count across two decades"           , lattice%count                       ,21            )

  ! The anchor interval is independent of the grid density. Anchoring every lattice step pins the bounds to the lattice points
  ! themselves; a request at x=3 with the default margin spans [1.5,6], so the bounds pin to 10^0.1 and 10^0.8.
  lattice=Range_Pinned(3.0d0,10,gridSchemePerDecade,anchorEvery=1)
  call Assert("per-decade: anchored to every lattice point",[lattice%indexMinimum,lattice%indexMaximum(),lattice%count],[1,8,8])

  ! Anchoring every half decade. A request at x=15 with the default margin spans [7.5,30], which pins outward to [10^0.5,10^1.5]
  ! - whereas anchoring to whole decades gives the far wider [10⁰,10²].
  lattice    =Range_Pinned(15.0d0,10,gridSchemePerDecade,anchorEvery=5)
  latticeWide=Range_Pinned(15.0d0,10,gridSchemePerDecade              )
  call Assert("per-decade: anchored to half decades" ,[lattice    %indexMinimum,lattice    %indexMaximum(),lattice    %count],[5, 15,11])
  call Assert("per-decade: anchored to whole decades",[latticeWide%indexMinimum,latticeWide%indexMaximum(),latticeWide%count],[0, 20,21])
  call Assert("per-decade: half-decade bounds"       ,[lattice%minimum(),lattice%maximum()],[10.0d0**0.5d0,10.0d0**1.5d0],relTol=1.0d-12)

  ! Explicitly anchoring every `pointsPer` steps must reproduce the default behavior exactly.
  latticeUnion=Range_Pinned(15.0d0,10,gridSchemePerDecade,anchorEvery=10)
  call Assert("per-decade: an anchor of `pointsPer` reproduces the default",[latticeUnion%indexMinimum,latticeUnion%count],[latticeWide%indexMinimum,latticeWide%count])

  ! Whatever the anchor, the points remain those of the same absolute lattice, so ranges built with different anchors are still
  ! bit-identical wherever they overlap. Here the half-decade range lies within the whole-decade one.
  values    =lattice    %values()
  valuesWide=latticeWide%values()
  offset    =lattice%indexMinimum-latticeWide%indexMinimum
  call Assert("per-decade: abscissae are independent of the anchor interval" , &
       &      all(valuesWide(offset+1:offset+lattice%count) == values)       , &
       &      .true.                                                           &
       &     )

  ! Points per octave. A request at x=3 with the default margin spans [1.5,6], pinning outward to [2⁰,2³].
  lattice=Range_Pinned(3.0d0,4,gridSchemePerOctave)
  call Assert("per-octave: margin then pin to whole octaves",[lattice%minimum(),lattice%maximum()],[1.0d0,8.0d0])
  call Assert("per-octave: exact point count"               , lattice%count                       ,13           )

  ! Points per unit interval. A request at x=2.5 with a unit margin spans [1.5,3.5], pinning outward to [1,4].
  lattice=Range_Pinned(2.5d0,2,gridSchemePerUnit)
  call Assert("per-unit: margin then pin to whole units",[lattice%minimum(),lattice%maximum()],[1.0d0,4.0d0])
  call Assert("per-unit: exact point count"             , lattice%count                       ,7            )
  values =lattice%values()
  call Assert("per-unit: lattice values"                , values                              ,[1.0d0,1.5d0,2.0d0,2.5d0,3.0d0,3.5d0,4.0d0])

  ! Hard limits must clip the pinned bounds, with the clipped edge snapped *inward* so that it remains a lattice point. Without
  ! limits the range here would be [10⁰,10¹]; the limits clip it to the lattice points just inside 1.5 and 5.
  lattice=Range_Pinned(3.0d0,10,gridSchemePerDecade,limitMinimum=1.5d0,limitMaximum=5.0d0)
  call Assert("per-decade: hard limits snap inward to lattice points",[lattice%indexMinimum,lattice%indexMaximum()],[2,6])
  call Assert("per-decade: clipped lower bound lies inside the limit",lattice%minimum() > 1.5d0,.true.)
  call Assert("per-decade: clipped upper bound lies inside the limit",lattice%maximum() < 5.0d0,.true.)

  ! Two ranges built on the same lattice must have bit-identical abscissae wherever they overlap - this is the property which
  ! makes value reuse on extension, and the merging of stored tables, exact.
  latticeNarrow=Range_Pinned(3.0d+0,10,gridSchemePerDecade)
  latticeWide  =Range_Pinned(3.0d+3,10,gridSchemePerDecade,latticeCurrent=latticeNarrow)
  valuesNarrow =latticeNarrow%values()
  valuesWide   =latticeWide  %values()
  offset       =latticeNarrow%indexMinimum-latticeWide%indexMinimum
  call Assert("per-decade: union with an existing lattice contains it"      ,latticeWide%covers(latticeNarrow),.true.)
  call Assert("per-decade: abscissae are bit-identical on the overlap"                                               , &
       &      all(valuesWide(offset+1:offset+latticeNarrow%count) == valuesNarrow)                                   , &
       &      .true.                                                                                                   &
       &     )

  ! The union of lattices is exact, and independent of the order in which the ranges were requested.
  latticeUnion=Range_Pinned(3.0d+0,10,gridSchemePerDecade,latticeCurrent=latticeWide)
  call Assert("per-decade: union is independent of request order",[latticeUnion%indexMinimum,latticeUnion%count],[latticeWide%indexMinimum,latticeWide%count])

  ! Extension of a tabulation held in a raw array, rather than in a table object.
  latticeNarrow=Range_Pinned(15.0d0,10,gridSchemePerDecade,anchorEvery=5)
  call Range_Lattice_Extend(rangeLattice(),latticeNarrow,valuesRaw,isComputedRaw)
  call Assert("raw array: extension of an unallocated array computes every point",count(isComputedRaw),0                  )
  call Assert("raw array: extension of an unallocated array sizes to the lattice",size (valuesRaw    ),latticeNarrow%count)
  valuesRaw=latticeNarrow%values()
  ! Extend onto a wider lattice, and check that the previously computed values are preserved at the right offset.
  latticeWide=Range_Pinned(1.5d3,10,gridSchemePerDecade,anchorEvery=5,latticeCurrent=latticeNarrow)
  offset     =Range_Lattice_Offset(latticeNarrow,latticeWide)
  call Range_Lattice_Extend(latticeNarrow,latticeWide,valuesRaw,isComputedRaw)
  call Assert("raw array: extension preserves precisely the previously computed points",count(isComputedRaw)                                 ,latticeNarrow%count)
  call Assert("raw array: extension marks the correct window as computed"              ,all(isComputedRaw(offset+1:offset+latticeNarrow%count)),.true.           )
  call Assert("raw array: extension preserves the values bit-for-bit"                                       , &
       &      all(valuesRaw(offset+1:offset+latticeNarrow%count) == latticeNarrow%values())                 , &
       &      .true.                                                                                          &
       &     )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Make_Ranges
