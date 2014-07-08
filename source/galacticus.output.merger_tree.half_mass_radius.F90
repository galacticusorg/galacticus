!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which handles outputting of galaxy half-mass radii.

module Galacticus_Output_Tree_Half_Mass_Radii
  !% Handles outputting of galaxy half-mass radii.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Half_Mass, Galacticus_Output_Tree_Half_Mass_Property_Count, Galacticus_Output_Tree_Half_Mass_Names

  ! Flag indicating if this module is initialize.
  logical :: outputHalfMassDataInitialized

  ! Flag indicating whether or not half-mass radii are to be output.
  logical :: outputHalfMassData

contains

  subroutine Galacticus_Output_Tree_Half_Mass_Initialize
    !% Initializes the module by determining whether or not half-mass radii should be output.
    use Input_Parameters
    use Stellar_Luminosities_Structure
    implicit none

    if (.not.outputHalfMassDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Half_Mass_Initialize)
       if (.not.outputHalfMassDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputHalfMassData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not half-mass radii should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputHalfMassData',outputHalfMassData,defaultValue=.false.)
          ! Flag that module is now initialized.
          outputHalfMassDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Half_Mass_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Half_Mass_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Half_Mass_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Half_Mass</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Half_Mass_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of half-light properties to be written to the \glc\ output file.
    use ISO_Varying_String
    use Numerical_Constants_Astronomical
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                                          :: iLuminosity

    ! Initialize the module.
    call Galacticus_Output_Tree_Half_Mass_Initialize

    ! Return property names if we are outputting half-light data.
    if (outputHalfMassData) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>halfMassRadius</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Radius enclosing half the galaxy stellar [Mpc]</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='halfMassRadius'
       doublePropertyComments(doubleProperty)='Radius enclosing half the galaxy stellar mass [Mpc]'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
    end if
    return
  end subroutine Galacticus_Output_Tree_Half_Mass_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Half_Mass_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Half_Mass</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Half_Mass_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of half-light properties to be written to the \glc\ output file.
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Half_Mass_Initialize

    ! Increment property count if we are outputting half-mass radii.
    if (outputHalfMassData) doublePropertyCount=doublePropertyCount+1
    return
  end subroutine Galacticus_Output_Tree_Half_Mass_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Half_Mass</unitName>
  !#  <sortName>Galacticus_Output_Tree_Half_Mass</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Half_Mass(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store density contrast properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Stellar_Luminosities_Structure
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)
    double precision                                         :: halfMassRadius

    ! Initialize the module.
    call Galacticus_Output_Tree_Half_Mass_Initialize

    ! Store property data if we are outputting half-light data.
    if (outputHalfMassData) then
       ! Get the half-light radius.
       halfMassRadius=Galactic_Structure_Radius_Enclosing_Mass(thisNode,fractionalMass=0.5d0,massType=massTypeStellar)
       ! Store the resulting radius.
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=halfMassRadius
    end if
    return
  end subroutine Galacticus_Output_Tree_Half_Mass

end module Galacticus_Output_Tree_Half_Mass_Radii
