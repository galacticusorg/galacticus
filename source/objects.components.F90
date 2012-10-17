!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module defining the component object type.

module Components
  !% Defines the merger tree object type.
  use Histories
  implicit none
  private
  public :: component, componentArray

  type component
     !% The component object type.
     double precision, allocatable, dimension(:,:) :: properties
     double precision, allocatable, dimension(:)   :: data
     type(history),    allocatable, dimension(:)   :: histories
  end type component

  type componentArray
     !% An array of component objects.
     type(component),  allocatable, dimension(:)   :: instance
     integer                                       :: activeInstanceIndex
   contains
     !@ <objectMethods>
     !@   <object>componentArray</object>
     !@   <objectMethod>
     !@     <method>activeInstance</method>
     !@     <description>Return the index of the active instance in a component array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>activeInstanceSet</method>
     !@     <description>Set the index of the active instance in a component array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>activeInstanceNullify</method>
     !@     <description>Nullfy the active instance in a component array (i.e. make no instance active).</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                     :: activeInstance        => Component_Active_Instance_Get
     procedure                                     :: activeInstanceSet     => Component_Active_Instance_Set
     procedure                                     :: activeInstanceNullify => Component_Active_Instance_Nullify
  end type componentArray
  
  ! Indices for second dimension of properties array used to store value, derivative and scale of each property.
  integer, public, parameter :: labelValue   =1,labelDerivative   =2,labelScale   =3
  integer, public, parameter :: propertyValue=1,propertyDerivative=2,propertyScale=2

  ! Labels for object types.
  integer, public, parameter :: objectTypeProperty= 1
  integer, public, parameter :: objectTypeHistory = 2

  ! Labels for instances.
  integer, public, parameter :: instanceNull      =-1

contains

  subroutine Component_Active_Instance_Set(thisComponentArray,thisInstance)
    !% Set the active instance for a component array.
    use Galacticus_Error
    implicit none
    class(componentArray), intent(inout) :: thisComponentArray
    integer,               intent(in)    :: thisInstance
    
    if (.not.allocated(thisComponentArray%instance)) call Galacticus_Error_Report('Component_Active_Instance_Set','attempt to set active instance in component array with no instances')
    if (thisInstance <= 0 .or. thisInstance > size(thisComponentArray%instance)) call Galacticus_Error_Report('Component_Active_Instance_Set','attempt to set active instance outside of range')
    thisComponentArray%activeInstanceIndex=thisInstance
    return
  end subroutine Component_Active_Instance_Set

  subroutine Component_Active_Instance_Nullify(thisComponentArray)
    !% Nullify the active instance for a component array.
    implicit none
    class(componentArray), intent(inout) :: thisComponentArray

    thisComponentArray%activeInstanceIndex=instanceNull
    return
  end subroutine Component_Active_Instance_Nullify

  integer function Component_Active_Instance_Get(thisComponentArray)
    !% Get the active instance for a component array.
    implicit none
    class(componentArray), intent(in) :: thisComponentArray

    Component_Active_Instance_Get=thisComponentArray%activeInstanceIndex
    return
  end function Component_Active_Instance_Get

end module Components
