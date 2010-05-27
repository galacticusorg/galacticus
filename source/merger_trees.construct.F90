!% Contains a module which constructs/destructs merger trees.

module Merger_Tree_Construction
  !% Constructs/destructs merger trees.
  use ISO_Varying_String
  use Merger_Trees
  private
  public :: Merger_Tree_Create

  ! Flag to indicate if this module has been initialized.  
  logical              :: mergerTreeConstructInitialized=.false.

  ! Name of tree construction method used.
  type(varying_string) :: mergerTreeConstructMethod
  ! Pointer to the subroutine that tabulates the transfer function and template interface for that subroutine.
  procedure(Merger_Tree_Construct_Template), pointer :: Merger_Tree_Construct => null()
  interface Merger_Tree_Construct_Template
     subroutine Merger_Tree_Construct_Template(thisTree)
       import mergerTree
       type(mergerTree), intent(inout) :: thisTree
     end subroutine Merger_Tree_Construct_Template
  end interface
  
contains

  function Merger_Tree_Create() result(thisTree)
    !% Creates a merger tree.
    use Input_Parameters
    use Galacticus_Error
    !# <include directive="mergerTreeConstructMethod" type="moduleUse">
    include 'merger_trees.construct.modules.inc'
    !# </include>
    implicit none
    type(mergerTree), pointer :: thisTree

    !$omp critical(Merger_Tree_Construct_Initialization) 
    ! Initialize if necessary.
    if (.not.mergerTreeConstructInitialized) then
       !@ <inputParameter>
       !@   <name>mergerTreeConstructMethod</name>
       !@   <defaultValue>build</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Selects the method to be used constructing merger trees.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeConstructMethod',mergerTreeConstructMethod,defaultValue='build')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="mergerTreeConstructMethod" type="code" action="subroutine">
       !#  <subroutineArgs>mergerTreeConstructMethod,Merger_Tree_Construct</subroutineArgs>
       include 'merger_trees.construct.inc'
       !# </include>
       if (.not.associated(Merger_Tree_Construct)) call Galacticus_Error_Report('Merger_Tree_Create','method '&
            &//char(mergerTreeConstructMethod)//' is unrecognized')
       mergerTreeConstructInitialized=.true.
    end if
    !$omp end critical(Merger_Tree_Construct_Initialization)

    ! Create the object.
    allocate(thisTree)

    ! Call the routine to construct the merger tree.
    call Merger_Tree_Construct(thisTree)

    ! Deallocate the tree if no nodes were created.
    if (.not.associated(thisTree%baseNode)) then
       deallocate(thisTree)
       thisTree => null()
    else
       ! Otherwise flag that the tree is uninitialized.
       thisTree%initialized=.false.
    end if
    return
  end function Merger_Tree_Create

end module Merger_Tree_Construction
