!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which stores internal state for the merger tree reading module.

module Merger_Tree_Read_State
  !% Stores internal state for the merger tree reading module.
  public :: Merger_Tree_Read_State_Store, Merger_Tree_Read_State_Retrieve

  ! Current position in the merger tree queue.
  integer, public :: mergerTreeQueuePosition=0  
                                             
contains
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Merger_Tree_Read_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Merger_Tree_Read_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile      
    type   (fgsl_file), intent(in   ) :: fgslStateFile  
                                                     
    write (stateFile) mergerTreeQueuePosition
    return
  end subroutine Merger_Tree_Read_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Merger_Tree_Read_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Merger_Tree_Read_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile      
    type   (fgsl_file), intent(in   ) :: fgslStateFile  
                                                     
    read (stateFile) mergerTreeQueuePosition
    return
  end subroutine Merger_Tree_Read_State_Retrieve

end module Merger_Tree_Read_State
