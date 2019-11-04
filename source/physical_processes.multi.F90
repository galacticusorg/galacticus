!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Contains a module which implements a multi physical process class.

  type, public :: multiProcessList
     class(physicalProcessClass), pointer :: process_
     type (multiProcessList    ), pointer :: next     => null()
  end type multiProcessList

  !# <physicalProcess name="physicalProcessMulti">
  !#  <description>A multi physical process property process class.</description>
  !# </physicalProcess>
  type, extends(physicalProcessClass) :: physicalProcessMulti
     !% A multi physical process output process class, which applies multiple physical processes.
     !% processs.
     private
     type(multiProcessList), pointer :: processs => null()
   contains
     final     ::                multiDestructor
     procedure :: nodePromote => multiNodePromote
     procedure :: deepCopy    => multiDeepCopy
  end type physicalProcessMulti

  interface physicalProcessMulti
     !% Constructors for the ``multi'' physical process class.
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface physicalProcessMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !% Constructor for the ``multi'' physical process property process class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (physicalProcessMulti)                :: self
    type   (inputParameters     ), intent(inout) :: parameters
    type   (multiProcessList    ), pointer       :: process_
    integer                                      :: i

    self    %processs => null()
    process_          => null()
    do i=1,parameters%copiesCount('physicalProcessMethod',zeroIfNotPresent=.true.)
       if (associated(process_)) then
          allocate(process_%next)
          process_ => process_%next
       else
          allocate(self%processs)
          process_ => self%processs
       end if
       !# <objectBuilder class="physicalProcess" name="process_%process_" source="parameters" copy="i" />
    end do
    return
  end function multiConstructorParameters

  function multiConstructorInternal(processs) result(self)
    !% Internal constructor for the ``multi'' output process property process class.
    implicit none
    type(physicalProcessMulti)                         :: self
    type(multiProcessList    ), target , intent(in   ) :: processs
    type(multiProcessList    ), pointer                :: process_

    self    %processs => processs
    process_          => processs
    do while (associated(process_))
       !# <referenceCountIncrement owner="process_" object="process_"/>
       process_ => process_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !% Destructor for the {\normalfont \ttfamily multi} output process property process class.
    implicit none
    type(physicalProcessMulti), intent(inout) :: self
    type(multiProcessList    ), pointer       :: process_, processNext

    if (associated(self%processs)) then
       process_ => self%processs
       do while (associated(process_))
          processNext => process_%next
          !# <objectDestructor name="process_%process_"/>
          deallocate(process_)
          process_ => processNext
       end do
    end if
    return
  end subroutine multiDestructor
  
  subroutine multiNodePromote(self,node)
    !% Act on a node promotion event.
    implicit none
    class(physicalProcessMulti), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (multiProcessList    ), pointer       :: process_

    process_ => self%processs
    do while (associated(process_))
       call process_%process_%nodePromote(node)
       process_ => process_%next
    end do
    return
  end subroutine multiNodePromote

  subroutine multiDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily multi} process class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(physicalProcessMulti), intent(inout) :: self
    class(physicalProcessClass), intent(inout) :: destination
    type (multiProcessList    ), pointer       :: process_   , processDestination_, &
         &                                        processNew_

    call self%physicalProcessClass%deepCopy(destination)
    select type (destination)
    type is (physicalProcessMulti)
       ! Copy list of processs.
       destination%processs => null         ()
       processDestination_  => null         ()
       process_             => self%processs
       do while (associated(process_))
          allocate(processNew_)
          if (associated(processDestination_)) then
             processDestination_%next     => processNew_
             processDestination_          => processNew_
          else
             destination        %processs => processNew_
             processDestination_          => processNew_
          end if
          allocate(processNew_%process_,mold=process_%process_)
          !# <deepCopy source="process_%process_" destination="processNew_%process_"/>
          process_ => process_%next
       end do
       class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine multiDeepCopy
