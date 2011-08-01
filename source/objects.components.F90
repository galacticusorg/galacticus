!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
