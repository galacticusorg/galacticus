!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements top-level functions for node components.

module Node_Components
  !% Implements top-level functions for node components.
  private
  public :: Node_Components_Initialize, Node_Components_Thread_Initialize

  logical :: initialized=.false.

contains

  subroutine Node_Components_Initialize(parameters)
    !% Perform initialization tasks for node components.
    use Input_Parameters
    !# <include directive="nodeComponentInitializationTask" type="moduleUse">
    include 'node_components.initialize.moduleUse.inc'
    !# </include>
    implicit none
    type(inputParameters), intent(inout) :: parameters

    if (.not.initialized) then
       !# <include directive="nodeComponentInitializationTask" type="functionCall" functionType="void">
       !#  <functionArgs>parameters</functionArgs>
       include 'node_components.initialize.inc'
       !# </include>
       initialized=.true.
    end if
    return
  end subroutine Node_Components_Initialize
  
  subroutine Node_Components_Thread_Initialize(parameters)
    !% Perform initialization tasks for node components.
    use Input_Parameters
    !# <include directive="nodeComponentThreadInitializationTask" type="moduleUse">
    include 'node_components.threadInitialize.moduleUse.inc'
    !# </include>
    implicit none
    type(inputParameters), intent(inout) :: parameters

    !# <include directive="nodeComponentThreadInitializationTask" type="functionCall" functionType="void">
    !#  <functionArgs>parameters</functionArgs>
    include 'node_components.threadInitialize.inc'
    !# </include>
    return
  end subroutine Node_Components_Thread_Initialize
  
end module Node_Components
