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
Contains a program to test the binary search tree code.
!!}

program Test_Binary_Search_Trees
  !!{
  Tests that binary search tree code works correctly.
  !!}
  use :: Display            , only : displayVerbositySet      , verbosityLevelStandard
  use :: Binary_Search_Trees, only : binaryTree               , binaryTreeNode
  use :: Unit_Tests         , only : Assert                   , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
          &                          compareGreaterThanOrEqual, compareLessThanOrEqual, compareEquals
  implicit none
  type(binaryTree    ), allocatable :: tree
  type(binaryTreeNode), pointer     :: left     , right      , &
       &                               successor, predecessor
 
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Binary search trees")

  ! Construct a simple tree.
  allocate(tree)
  tree=binaryTree()
  call tree%insert(3.0d0,1.0d0/3.0d0)
  call tree%insert(2.0d0,1.0d0/2.0d0)
  call tree%insert(8.0d0,1.0d0/8.0d0)
  call tree%insert(6.0d0,1.0d0/6.0d0)
  call tree%insert(1.0d0,1.0d0/1.0d0)

  call tree%bracket( 4.0d0,left,right)
  predecessor => left %predecessor()
  successor   => right%successor  ()
  call Assert("bracket exists         ( 4.0)",[associated(left),associated(right)],[.true.,.true.])
  call Assert("bracket left           ( 4.0)",4.0d0  ,           left       %key ,compare=compareGreaterThanOrEqual)
  call Assert("bracket right          ( 4.0)",4.0d0  ,           right      %key ,compare=compareLessThanOrEqual   )
  call Assert("successor              ( 4.0)",8.0d0  ,           successor  %key ,compare=compareEquals            )
  call Assert("predecessor            ( 4.0)",2.0d0  ,           predecessor%key ,compare=compareEquals            )
  
  call tree%bracket( 1.0d0,left,right)
  predecessor => left %predecessor()
  successor   => right%successor  ()
  call Assert("bracket exists         ( 1.0)",[associated(left),associated(right)],[.true.,.true.])
  call Assert("bracket left           ( 1.0)",1.0d0  ,           left       %key ,compare=compareGreaterThanOrEqual)
  call Assert("bracket right          ( 1.0)",1.0d0  ,           right      %key ,compare=compareLessThanOrEqual   )
  call Assert("successor              ( 1.0)",2.0d0  ,           successor  %key ,compare=compareEquals            )
  call Assert("predecessor            ( 1.0)",.false.,associated(predecessor    )                                  )

  call tree%bracket( 8.0d0,left,right)
  predecessor => left %predecessor()
  successor   => right%successor  ()
  call Assert("bracket exists         ( 8.0)",[associated(left),associated(right)],[.true.,.true.])
  call Assert("bracket left           ( 8.0)",8.0d0  ,           left       %key ,compare=compareGreaterThanOrEqual)
  call Assert("bracket right          ( 8.0)",8.0d0  ,           right      %key ,compare=compareLessThanOrEqual   )
  call Assert("successor              ( 8.0)",.false.,associated(successor      )                                  )
  call Assert("predecessor            ( 8.0)",6.0d0  ,           predecessor%key ,compare=compareEquals            )

  call tree%bracket(12.0d0,left,right)
  call Assert("bracket does not exist (12.0)",associated(left) .and. associated(right),.false.)

  call tree%bracket(0.5d0,left,right)
  call Assert("bracket does not exist ( 0.5)",associated(left) .and. associated(right),.false.)

  deallocate(tree)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Binary_Search_Trees
