!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  private
  public :: component

  type component
     !% The component object type.
     double precision, allocatable, dimension(:,:) :: properties
     double precision, allocatable, dimension(:)   :: data
     type(history),    allocatable, dimension(:)   :: histories ! memoryManagementIgnore (force memory management system to ignore)
     type(component),  pointer                     :: nextComponentOfType
  end type component

  ! Indices for second dimension of properties array used to store value and derivative of each property.
  integer, public, parameter :: propertyValue=1,propertyDerivative=2

end module Components
