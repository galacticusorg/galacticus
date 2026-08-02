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
       &                          rangeLattice       , gridSchemePerDecade   , gridSchemePerOctave , gridSchemePerUnit
  use :: Unit_Tests      , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                               , dimension(1:11) :: range1
  double precision                               , dimension(0:10) :: range2
  type            (rangeLattice)                                   :: lattice     , latticeNarrow, &
       &                                                              latticeWide , latticeUnion
  double precision              , allocatable    , dimension(:   ) :: values      , valuesNarrow , &
       &                                                              valuesWide
  integer                                                          :: offset

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

  ! Anchoring to the grid lattice itself, rather than to whole decades.
  lattice=Range_Pinned(3.0d0,10,gridSchemePerDecade,anchorToGrid=.true.)
  call Assert("per-decade: anchored to the grid",[lattice%indexMinimum,lattice%indexMaximum(),lattice%count],[1,8,8])

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

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Make_Ranges
