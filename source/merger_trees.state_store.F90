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
Contains a module which provides state store/restore functionality for merger trees.
!!}

module Merger_Tree_State_Store
  !!{
  Provides state store/restore functionality for merger trees.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  public

  integer(c_size_t  ) :: treeStateStoreSequence=-1_c_size_t
  !$omp threadprivate(treeStateStoreSequence)

contains

  !![
  <stateStoreTask>
   <unitName>mergerTreeStateStore</unitName>
  </stateStoreTask>
  !!]
  subroutine mergerTreeStateStore(stateFile,gslStateFile,stateOperatorID)
    !!{
    Write the stored snapshot of the random number state to file.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_ptr
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperatorID
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    !$GLC attributes unused :: gslStateFile, stateOperatorID

    write (stateFile) treeStateStoreSequence
    return
  end subroutine mergerTreeStateStore

  !![
  <stateRetrieveTask>
   <unitName>mergerTreeStateRestore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine mergerTreeStateRestore(stateFile,gslStateFile,stateOperatorID)
    !!{
    Write the stored snapshot of the random number state to file.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_ptr
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperatorID
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    !$GLC attributes unused :: gslStateFile, stateOperatorID

    read (stateFile) treeStateStoreSequence
    return
  end subroutine mergerTreeStateRestore

end module Merger_Tree_State_Store
