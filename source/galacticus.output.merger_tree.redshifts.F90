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

!% Contains a module which handles outputting of node redshifts to the \glc\ output file.

module Galacticus_Output_Trees_Redshifts
  !% Handles outputting of node redshift data to the \glc\ output file.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Redshifts, Galacticus_Output_Redshifts_Property_Count, Galacticus_Output_Redshifts_Names

  ! Flag indicating whether redshift output is required.
  logical            :: outputNodeRedshifts

  ! Flag indicating if module is initialized.
  logical            :: redshiftOutputIsInitialized=.false.

  ! Redshift properties.
  logical            :: timeLastIsolatedIsAvailable
  integer            :: redshiftPropertyCount

contains

  subroutine Galacticus_Output_Redshifts_Initalize()
    !% Intialize the ``redshift'' output module.
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.redshiftOutputIsInitialized) then
       !$omp critical(Galacticus_Output_Redshifts_Initalization)
       if (.not.redshiftOutputIsInitialized) then
          ! Read parameter controlling whether or not this module should output.
          !@ <inputParameter>
          !@   <name>outputNodeRedshifts</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Controls whether or not the redshifts corresponding to node times should be output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputNodeRedshifts',outputNodeRedshifts,defaultValue=.false.)
          ! Evaluate which properties we can output.
          redshiftPropertyCount=0
          timeLastIsolatedIsAvailable=.false.
          if (outputNodeRedshifts) then
             timeLastIsolatedIsAvailable=defaultBasicComponent%timeLastIsolatedIsGettable()
             if (timeLastIsolatedIsAvailable) redshiftPropertyCount=redshiftPropertyCount+1
          end if
          ! Flag that the module is now initialized.
          redshiftOutputIsInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Redshifts_Initalization)
    end if

    return
  end subroutine Galacticus_Output_Redshifts_Initalize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Redshifts_Names</unitName>
  !#  <sortName>Galacticus_Output_Redshifts</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Redshifts_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of link properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    call Galacticus_Output_Redshifts_Initalize()
    if (timeLastIsolatedIsAvailable) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>basicRedshiftLastIsolated</name>
       !@   <datatype>double</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The redshift of the epoch at which this node was last isolated.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='basicRedshiftLastIsolated'
       doublePropertyComments(doubleProperty)="The reshift of the epich at which the node was last isolated."
       doublePropertyUnitsSI (doubleProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Redshifts_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Redshifts_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Redshifts</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Redshifts_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    call Galacticus_Output_Redshifts_Initalize()
    doublePropertyCount=doublePropertyCount+redshiftPropertyCount
    return
  end subroutine Galacticus_Output_Redshifts_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Redshifts</unitName>
  !#  <sortName>Galacticus_Output_Redshifts</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Redshifts(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store link properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Cosmology_Functions
    implicit none
    double precision                         , intent(in   )          :: time
    type            (treeNode               ), intent(inout), pointer :: thisNode
    integer                                  , intent(inout)          :: doubleBufferCount            , doubleProperty , &
         &                                                               integerBufferCount           , integerProperty
    integer         (kind=kind_int8         ), intent(inout)          :: integerBuffer           (:,:)
    double precision                         , intent(inout)          :: doubleBuffer            (:,:)
    class           (nodeComponentBasic     )               , pointer :: thisBasic
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctionsDefault

    call Galacticus_Output_Redshifts_Initalize()
    if (timeLastIsolatedIsAvailable) then
       cosmologyFunctionsDefault => cosmologyFunctions      ()
       thisBasic                 => thisNode          %basic()
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=                 &
            & cosmologyFunctionsDefault %redshiftFromExpansionFactor(  &
            &  cosmologyFunctionsDefault%expansionFactor             ( &
            &   thisBasic%timeLastIsolated()                           &
            &  )                                                       &
            & )
    end if
    return
  end subroutine Galacticus_Output_Redshifts

end module Galacticus_Output_Trees_Redshifts
