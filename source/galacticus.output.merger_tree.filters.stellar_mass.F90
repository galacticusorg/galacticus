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

!% Contains a module which filters output on stellar mass.

module Galacticus_Merger_Tree_Output_Filter_Stellar_Masses
  !% Filters output for lightcone geometry.
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output_Filter_Stellar_Mass,Galacticus_Merger_Tree_Output_Filter_Stellar_Mass_Initialize

  ! Flags indicating if the module has been initialized and if this filter is active.
  logical                                      :: stellarMassFilterInitialized=.false.
  logical                                      :: stellarMassFilterActive

  ! The threshold for output.
  double precision                             :: stellarMassFilterThreshold

contains

  !# <mergerTreeOutputFilterInitialize>
  !#   <unitName>Galacticus_Merger_Tree_Output_Filter_Stellar_Mass_Initialize</unitName>
  !# </mergerTreeOutputFilterInitialize>
  subroutine Galacticus_Merger_Tree_Output_Filter_Stellar_Mass_Initialize(filterNames)
    !% Initializes the stellar mass filter module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string), intent(in), dimension(:) :: filterNames

    ! Initialize the filter if necessary.
    if (.not.stellarMassFilterInitialized) then
       ! Determine if this filter has been selected.
       stellarMassFilterActive=any(filterNames == "stellarMass")
       ! If this filter is active, read the minimum stellar mass.
       if (stellarMassFilterActive) then          
          ! Get the minimum stellar mass for output
          !@ <inputParameter>
          !@   <name>stellarMassFilterThreshold</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum stellar mass of a galaxy to pass the {\tt stellarMass} output filter (in units of $M_\odot$).
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('stellarMassFilterThreshold',stellarMassFilterThreshold)
       end if
       ! Flag that this filter is now initialized.
       stellarMassFilterInitialized=.true.
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output_Filter_Stellar_Mass_Initialize
  
  !# <mergerTreeOutputFilter>
  !#   <unitName>Galacticus_Merger_Tree_Output_Filter_Stellar_Mass</unitName>
  !# </mergerTreeOutputFilter>
  subroutine Galacticus_Merger_Tree_Output_Filter_Stellar_Mass(thisNode,doOutput)
    !% Determines whether {\tt thisNode} has sufficient stellar mass to be output.
    use Galacticus_Nodes
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    type(treeNode)  , intent(inout), pointer :: thisNode
    logical         , intent(inout)          :: doOutput
    double precision                         :: stellarMass

    ! Return immediately if this filter is not active.
    if (.not.stellarMassFilterActive) return

    ! Get the total stellar mass of the galaxy.
    stellarMass=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeStellar)

    ! Filter out the galaxy if it is below the stellar mass threshold.
    if (stellarMass < stellarMassFilterThreshold) doOutput=.false.
    return
  end subroutine Galacticus_Merger_Tree_Output_Filter_Stellar_Mass

end module Galacticus_Merger_Tree_Output_Filter_Stellar_Masses
