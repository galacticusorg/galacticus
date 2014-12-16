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

!% Contains a module which filters output on stellar mass spheroid-to-total ratio.

module Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgs
  !% Filters output on stellar mass spheroid-to-total ratio.
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgy,Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgy_Initialize

  ! Flags indicating if the module has been initialized and if this filter is active.
  logical          :: morphologyFilterInitialized=.false.
  logical          :: morphologyFilterActive

  ! The ratio thresholds.
  double precision :: stellarMassMorphologyFilterRatioMinimum, stellarMassMorphologyFilterRatioMaximum

contains

  !# <mergerTreeOutputFilterInitialize>
  !#   <unitName>Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgy_Initialize</unitName>
  !# </mergerTreeOutputFilterInitialize>
  subroutine Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgy_Initialize(filterNames)
    !% Initializes the stellar mass filter module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string), dimension(:), intent(in   ) :: filterNames

    ! Initialize the filter if necessary.
    if (.not.morphologyFilterInitialized) then
       ! Determine if this filter has been selected.
       morphologyFilterActive=any(filterNames == "stellarMassMorphology")
       ! If this filter is active, read the minimum and maximum allowed ratios.
       if (morphologyFilterActive) then
          ! Get the magnitude limits.
          !@ <inputParameter>
          !@   <name>stellarMassMorphologyFilterRatioMinimum</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum spheroid-to-total ratio in stellar mass for which to allow output when the ``{\tt stellarMassMorphology}'' output filter is active.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('stellarMassMorphologyFilterRatioMinimum',stellarMassMorphologyFilterRatioMinimum)
          !@ <inputParameter>
          !@   <name>stellarMassMorphologyFilterRatioMaximum</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum spheroid-to-total ratio in stellar mass for which to allow output when the ``{\tt stellarMassMorphology}'' output filter is active.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('stellarMassMorphologyFilterRatioMaximum',stellarMassMorphologyFilterRatioMaximum)
       end if
       ! Flag that this filter is now initialized.
       morphologyFilterInitialized=.true.
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgy_Initialize

  !# <mergerTreeOutputFilter>
  !#   <unitName>Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgy</unitName>
  !# </mergerTreeOutputFilter>
  subroutine Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgy(thisNode,doOutput)
    !% Determines whether {\tt thisNode} has sufficient stellar mass to be output.
    use Galacticus_Nodes
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    logical                                , intent(inout)          :: doOutput
    class           (nodeComponentDisk    )               , pointer :: thisDisk
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroid
    double precision                                                :: massTotal

    ! Return immediately if this filter is not active.
    if (.not.morphologyFilterActive) return

    ! Get the disk and spheroid component.
    thisDisk     => thisNode%disk    ()
    thisSpheroid => thisNode%spheroid()

    ! Find the total stellar mass.
    massTotal=thisDisk%massStellar()+thisSpheroid%massStellar()
    if (massTotal > 0.0d0) then
       doOutput=                                                                                  &
            &   (                                                                                 &
            &     thisSpheroid%massStellar()/massTotal >= stellarMassMorphologyFilterRatioMinimum &
            &    .and.                                                                            &
            &     thisSpheroid%massStellar()/massTotal <= stellarMassMorphologyFilterRatioMaximum &
            &   )
    else
       doOutput=.false.
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgy

end module Galacticus_Merger_Tree_Output_Stllr_Mss_Mrphlgs
