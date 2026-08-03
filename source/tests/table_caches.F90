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
Contains a program to test persistent, mergeable caches of tabulated functions.
!!}

program Test_Table_Caches
  !!{RST
  Tests that persistent caches of tabulated functions round-trip exactly, and that a cached tabulation whose range only
  partially covers that held in memory is *merged* with it rather than discarded---the property which allows the files under
  ``pathTypeDataDynamic`` to grow monotonically towards the union of every range ever requested.
  !!}
  use :: Cosmology_Functions                  , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                 , only : cosmologyParametersSimple
  use :: Display                              , only : displayVerbositySet                       , verbosityLevelStandard
  use :: Error                                , only : errorStatusSuccess
  use :: Events_Hooks                         , only : eventsHooksInitialize
  use :: Intergalactic_Medium_Filtering_Masses, only : intergalacticMediumFilteringMassGnedin2000
  use :: Intergalactic_Medium_State           , only : intergalacticMediumStateSimple
  use :: Linear_Growth                        , only : linearGrowthCollisionlessMatter
  use :: File_Utilities                       , only : Directory_Make                            , File_Exists             , File_Remove
  use :: IO_HDF5                              , only : ioHDF5AccessInitialize
  use :: ISO_Varying_String                   , only : varying_string                            , assignment(=)
  use :: Numerical_Ranges                     , only : Range_Pinned                              , rangeLattice            , gridSchemePerDecade
  use :: ISO_Varying_String                   , only : char
  use :: Table_Caches                         , only : Table_Cache_Restore                       , Table_Cache_Store       , Table_Cache_File_Name
  use :: Tables                               , only : table1D                                   , table1DLogarithmicLinear, table2DLogLogLin
  use :: Unit_Tests                           , only : Assert                                    , Unit_Tests_Begin_Group  , Unit_Tests_End_Group , Unit_Tests_Finish
  implicit none
  class           (table1D                                   ), allocatable                 :: tableStored              , tableRestored         , &
       &                                                                                       tableMerged              , tableFine
  type            (rangeLattice                              )                              :: latticeNarrow            , latticeWide           , &
       &                                                                                       latticePartial
  type            (varying_string                            )                              :: fileName                 , fileNameFiltering     , &
       &                                                                                       fileName2D
  type            (table2DLogLogLin                          )                              :: table2DStored            , table2DRestored
  type            (rangeLattice                              )                              :: latticeX2D               , latticeY2D
  logical                                                     , allocatable, dimension(:,:) :: isComputed2D
  double precision                                            , allocatable, dimension(:,:) :: z2DStored                , z2DRestored
  logical                                                     , allocatable, dimension(:  ) :: isComputed
  double precision                                            , allocatable, dimension(:  ) :: xStored                  , xRestored
  double precision                                            , allocatable, dimension(:,:) :: yStored                  , yRestored
  integer                                                                                   :: status                   , i                     , &
       &                                                                                       j
  type            (cosmologyParametersSimple                 )                              :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda            )                              :: cosmologyFunctions_
  type            (intergalacticMediumStateSimple            )                              :: intergalacticMediumState_
  type            (linearGrowthCollisionlessMatter           )                              :: linearGrowth_
  type            (intergalacticMediumFilteringMassGnedin2000)                              :: filteringMass_
  double precision                                                                          :: massFilteringLate        , massFilteringLateAgain, &
       &                                                                                       massFilteringEarly

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize the lock used to serialize access to the HDF5 library.
  call ioHDF5AccessInitialize()
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Table caches")
  ! Work in a scratch file, removing any left over from a previous run so that the test is hermetic.
  call Directory_Make('testSuite/outputs')
  fileName='testSuite/outputs/tableCache.hdf5'
  if (File_Exists(fileName)) call File_Remove(fileName)

  ! Build a tabulation of y=x² over a narrow range and store it.
  latticeNarrow=Range_Pinned(15.0d0,10,gridSchemePerDecade,anchorEvery=5)
  allocate(table1DLogarithmicLinear :: tableStored)
  call tableStored%extend(latticeNarrow,isComputed)
  call populate(tableStored)
  xStored=tableStored%xs()
  yStored=tableStored%ys()
  call Table_Cache_Store(tableStored,fileName)
  call Assert('a cache file is written',File_Exists(fileName),.true.)

  ! Restore into a fresh table: the abscissae, the values, and the lattice must all be recovered exactly.
  call Table_Cache_Restore(tableRestored,fileName,status)
  xRestored=tableRestored%xs()
  yRestored=tableRestored%ys()
  call Assert('a cached tabulation is restored'               ,status                                                          ,errorStatusSuccess                              )
  call Assert('the restored tabulation is on the same lattice',[tableRestored%lattice%indexMinimum,tableRestored%lattice%count],[latticeNarrow%indexMinimum,latticeNarrow%count])
  call Assert('the restored abscissae are bit-identical'      ,all(xRestored      == xStored      )                            ,.true.                                          )
  call Assert('the restored values are bit-identical'         ,all(yRestored(:,1) == yStored(:,1) )                            ,.true.                                          )

  ! Extend the restored table downwards, compute the newly added points, and write it back. The cache must then hold the union.
  latticeWide=Range_Pinned(0.05d0,10,gridSchemePerDecade,anchorEvery=5,latticeCurrent=tableRestored%lattice)
  call tableRestored%extend(latticeWide,isComputed)
  call Assert('extension of a restored table preserves the cached points',count(isComputed),latticeNarrow%count)
  call populate(tableRestored)
  call Table_Cache_Store(tableRestored,fileName)
  deallocate(tableRestored)
  call Table_Cache_Restore(tableRestored,fileName,status)
  call Assert('the cache grows to the union of the ranges written to it',[tableRestored%lattice%indexMinimum,tableRestored%lattice%count],[latticeWide%indexMinimum,latticeWide%count])
  call Assert('every value in the grown cache is correct'               ,agrees(tableRestored)                                           ,.true.                                      )

  ! Reset the cache to the narrow range, then store a tabulation which only partially overlaps it. Storing must merge the two,
  ! leaving both the caller's table and the file holding the union - with the values contributed by the file spliced in exactly.
  call File_Remove(fileName)
  call Table_Cache_Store(tableStored,fileName)
  latticePartial=rangeLattice(gridSchemePerDecade,10,latticeNarrow%indexMinimum+5,latticeNarrow%count+5)
  allocate(table1DLogarithmicLinear :: tableMerged)
  call tableMerged%extend(latticePartial,isComputed)
  call populate(tableMerged)
  call Table_Cache_Store(tableMerged,fileName)
  call Assert('storing a partially overlapping tabulation merges it with the cache', &
       &      [tableMerged%lattice%indexMinimum,tableMerged%lattice%indexMaximum()], &
       &      [latticeNarrow%indexMinimum      ,latticePartial%indexMaximum()     ]  &
       &     )
  call Assert('the values spliced in from the cache are correct',agrees(tableMerged),.true.)

  ! A cached tabulation on an incommensurate lattice must be ignored rather than merged.
  allocate(table1DLogarithmicLinear :: tableFine)
  call tableFine%extend(rangeLattice(gridSchemePerDecade,20,latticeNarrow%indexMinimum,latticeNarrow%count),isComputed)
  call populate(tableFine)
  call Table_Cache_Restore(tableFine,fileName,status)
  call Assert('a cached tabulation on an incommensurate lattice is ignored',[tableFine%lattice%pointsPer,tableFine%lattice%count],[20,latticeNarrow%count])

  ! The same round trip for a two-dimensional table, whose axes are pinned independently.
  fileName2D='testSuite/outputs/tableCache2D.hdf5'
  if (File_Exists(fileName2D)) call File_Remove(fileName2D)
  latticeX2D=Range_Pinned(15.0d0,4,gridSchemePerDecade,anchorEvery=2)
  latticeY2D=Range_Pinned( 3.0d0,4,gridSchemePerDecade,anchorEvery=2)
  call table2DStored%extend(latticeX2D,latticeY2D,isComputed2D)
  do i=1,latticeX2D%count
     do j=1,latticeY2D%count
        call table2DStored%populate(table2DStored%x(i)*table2DStored%y(j),i,j)
     end do
  end do
  z2DStored=table2DStored%zs()
  call Table_Cache_Store  (table2DStored  ,fileName2D       )
  call Table_Cache_Restore(table2DRestored,fileName2D,status)
  z2DRestored=table2DRestored%zs()
  call Assert('a cached 2D tabulation is restored'                ,status,errorStatusSuccess)
  call Assert('the restored 2D tabulation is on the same lattices',                                                                                        &
       &      [table2DRestored%latticeX%indexMinimum,table2DRestored%latticeX%count,table2DRestored%latticeY%indexMinimum,table2DRestored%latticeY%count], &
       &      [latticeX2D     %indexMinimum         ,latticeX2D     %count         ,latticeY2D     %indexMinimum         ,latticeY2D     %count         ]  &
       &     )
  call Assert('the restored 2D values are bit-identical'          ,all(z2DRestored == z2DStored),.true.)

  ! Exercise a real, persisted tabulation which uses the cache: the Gnedin2000 filtering mass. The key property is that
  ! requesting an earlier epoch, which forces the table to be extended downwards, must not change the value already tabulated at
  ! a later epoch - the previously computed points are reused, not recomputed on a shifted grid.
  cosmologyParameters_     =cosmologyParametersSimple                (                                                                 &
       &                                                              OmegaMatter               = 0.2750d0                           , &
       &                                                              OmegaBaryon               = 0.0458d0                           , &
       &                                                              OmegaDarkEnergy           = 0.7250d0                           , &
       &                                                              temperatureCMB            = 2.7800d0                           , &
       &                                                              HubbleConstant            =70.2000d0                             &
       &                                                             )
  cosmologyFunctions_      =cosmologyFunctionsMatterLambda           (                                                                 &
       &                                                              cosmologyParameters_      =cosmologyParameters_                  &
       &                                                             )
  intergalacticMediumState_=intergalacticMediumStateSimple           (                                                                 &
       &                                                              reionizationRedshift      = 8.00d0                             , &
       &                                                              reionizationTemperature   = 1.00d4                             , &
       &                                                              preReionizationTemperature= 1.00d4                             , &
       &                                                              cosmologyFunctions_       =cosmologyFunctions_                 , &
       &                                                              cosmologyParameters_      =cosmologyParameters_                  &
       &                                                             )
  linearGrowth_            =linearGrowthCollisionlessMatter          (                                                                 &
       &                                                              cosmologyParameters_      =cosmologyParameters_                , &
       &                                                              cosmologyFunctions_       =cosmologyFunctions_                   &
       &                                                             )
  filteringMass_           =intergalacticMediumFilteringMassGnedin2000(                                                                &
       &                                                              timeTooEarlyIsFatal       =.true.                              , &
       &                                                              cosmologyParameters_      =cosmologyParameters_                , &
       &                                                              cosmologyFunctions_       =cosmologyFunctions_                 , &
       &                                                              linearGrowth_             =linearGrowth_                       , &
       &                                                              intergalacticMediumState_ =intergalacticMediumState_             &
       &                                                             )
  ! Discard any tabulation cached by a previous run, so that this test always exercises the cold-cache path.
  fileNameFiltering=Table_Cache_File_Name(                                                                                             &
       &                                  subDirectory    ='intergalacticMedium'                                                     , &
       &                                  objectType      =char(filteringMass_%objectType      (                                   )), &
       &                                  hashedDescriptor=char(filteringMass_%hashedDescriptor(includeSourceDigest         =.true.  , &
       &                                                                                        includeFileModificationTimes=.true.))  &
       &                                 )
  if (File_Exists(fileNameFiltering)) call File_Remove(fileNameFiltering)
  massFilteringLate     =filteringMass_%massFiltering(1.0d0)
  massFilteringEarly    =filteringMass_%massFiltering(0.2d0)
  massFilteringLateAgain=filteringMass_%massFiltering(1.0d0)
  call Assert('the filtering mass is positive'                                     ,massFilteringLate      >  0.0d0            ,.true.)
  call Assert('extending the filtering mass table downwards reaches earlier epochs',massFilteringEarly     <  massFilteringLate,.true.)
  call Assert('extending the filtering mass table leaves earlier results unchanged',massFilteringLateAgain == massFilteringLate,.true.)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  subroutine populate(table_)
    !!{RST
    Populate any as-yet uncomputed points of ``table_`` with :math:`y=x^2`.
    !!}
    implicit none
    class  (table1D), intent(inout), allocatable :: table_
    integer                                      :: j

    select type (table_)
    type is (table1DLogarithmicLinear)
       do j=1,table_%size()
          call table_%populate(table_%x(j)**2,j)
       end do
    end select
    return
  end subroutine populate

  logical function agrees(table_)
    !!{RST
    Return true if every value in ``table_`` is that expected for :math:`y=x^2`.
    !!}
    implicit none
    class  (table1D), intent(inout), allocatable :: table_
    integer                                      :: j

    agrees=.true.
    do j=1,table_%size()
       if (table_%y(j) /= table_%x(j)**2) agrees=.false.
    end do
    return
  end function agrees

end program Test_Table_Caches
