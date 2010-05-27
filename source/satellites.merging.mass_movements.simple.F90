!% Contains a module which implements a simple model of mass movements during satellite mergers.

module Satellite_Merging_Mass_Movements_Simple
  !% Implements a simple model of mass movements during satellite mergers.
  private
  public :: Satellite_Merging_Mass_Movements_Simple_Initialize

  ! Mass ratio above which a merger is considered to be "major".
  double precision :: majorMergerMassRatio

contains

  !# <satelliteMergingMassMovementsMethod>
  !#  <unitName>Satellite_Merging_Mass_Movements_Simple_Initialize</unitName>
  !# </satelliteMergingMassMovementsMethod>
  subroutine Satellite_Merging_Mass_Movements_Simple_Initialize(satelliteMergingMassMovementsMethod,Satellite_Merging_Mass_Movement_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingMassMovementsMethod
    procedure(),          pointer, intent(inout) :: Satellite_Merging_Mass_Movement_Get
    
    if (satelliteMergingMassMovementsMethod == 'simple') then
       Satellite_Merging_Mass_Movement_Get => Satellite_Merging_Mass_Movement_Simple
       !@ <inputParameter>
       !@   <name>majorMergerMassRatio</name>
       !@   <defaultValue>0.3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass ratio above which mergers are considered to be ``major'' in the simple merger mass movements method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("majorMergerMassRatio",majorMergerMassRatio,defaultValue=0.3d0)
    end if
    return
  end subroutine Satellite_Merging_Mass_Movements_Simple_Initialize

  subroutine Satellite_Merging_Mass_Movement_Simple(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo)
    !% Return orbital velocities of a satellite selected at random from the fitting function found by \cite{benson_orbital_2005}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Satellite_Merging_Mass_Movements_Descriptors
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(out)             :: gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo
    type(treeNode),                pointer  :: hostNode
    double precision                        :: satelliteMass,hostMass

    ! Get the host node.
    hostNode => thisNode%parentNode

    ! Find the baryonic masses of the two galaxies.
    satelliteMass=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic)
    hostMass     =Galactic_Structure_Enclosed_Mass(hostNode,massType=massTypeGalactic)

    ! Decide if the mass ratio is large enough to trigger a major merger.
    if (satelliteMass >= majorMergerMassRatio*hostMass) then
       gasMovesTo     =movesToSpheroid
       starsMoveTo    =movesToSpheroid
       hostGasMovesTo =movesToSpheroid
       hostStarsMoveTo=movesToSpheroid
    else
       gasMovesTo     =movesToDisk
       starsMoveTo    =movesToSpheroid
       hostGasMovesTo =doesNotMove
       hostStarsMoveTo=doesNotMove
    end if

    return
  end subroutine Satellite_Merging_Mass_Movement_Simple

end module Satellite_Merging_Mass_Movements_Simple
