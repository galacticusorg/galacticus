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
