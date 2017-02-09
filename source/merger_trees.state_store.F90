!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which provides state store/restore functionality for merger trees.

module Merger_Tree_State_Store
  !% Provides state store/restore functionality for merger trees.
  use Galacticus_Nodes
  use Pseudo_Random
  public

  type(mergerTree  ), pointer :: treeStateStore            => null()
  type(pseudoRandom)          :: treeRandomNumberGenerator
  !$omp threadprivate(treeStateStore,treeRandomNumberGenerator)

contains

  !# <galacticusStateSnapshotTask>
  !#  <unitName>mergerTreeStateSnapshot</unitName>
  !# </galacticusStateSnapshotTask> 
  subroutine mergerTreeStateSnapshot()
    !% Store a snapshot of the random number generator internal state.
    use Pseudo_Random
    implicit none

    if (associated(treeStateStore)) treeRandomNumberGenerator=treeStateStore%randomNumberGenerator%clone()
    return
  end subroutine mergerTreeStateSnapshot

  !# <galacticusStateStoreTask>
  !#  <unitName>mergerTreeStateStore</unitName>
  !# </galacticusStateStoreTask>
  subroutine mergerTreeStateStore(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    call treeRandomNumberGenerator%store(stateFile,fgslStateFile)
    return
  end subroutine mergerTreeStateStore

  !# <galacticusStateRetrieveTask>
  !#  <unitName>mergerTreeStateRestore</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine mergerTreeStateRestore(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    call treeStateStore%randomNumberGenerator%restore(stateFile,fgslStateFile)
    return
  end subroutine mergerTreeStateRestore

end module Merger_Tree_State_Store
