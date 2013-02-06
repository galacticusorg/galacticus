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

!% Contains a module which handles outputting of galaxy half-light properties (radii and masses).

module Galacticus_Output_Tree_Half_Light_Properties
  !% Handles outputting of galaxy half-light radii and associated masses.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Half_Light, Galacticus_Output_Tree_Half_Light_Property_Count, Galacticus_Output_Tree_Half_Light_Names

  ! Flag indicating if this module is initialize.
  logical :: outputHalfLightDataInitialized

  ! Number of luminosities.
  integer :: luminositiesCount

  ! Number of properties.
  integer :: halfLightPropertyCount

  ! Flag indicating whether or not half-light data is to be output.
  logical :: outputHalfLightData

contains

  subroutine Galacticus_Output_Tree_Half_Light_Initialize
    !% Initializes the module by determining whether or not half-light radius data should be output.
    use Input_Parameters
    use Stellar_Population_Properties_Luminosities
    implicit none

    if (.not.outputHalfLightDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Half_Light_Initialize)
       if (.not.outputHalfLightDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputHalfLightData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not half-light radius data (i.e. radius and mass) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputHalfLightData',outputHalfLightData,defaultValue=.false.)
          
          ! Get the number of luminosities in use.
          if (outputHalfLightData) then
             luminositiesCount=Stellar_Population_Luminosities_Count()
             halfLightPropertyCount=2*luminositiesCount
          end if
          
          ! Flag that module is now initialized.
          outputHalfLightDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Half_Light_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Half_Light_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Half_Light_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Half_Light</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Half_Light_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of half-light properties to be written to the \glc\ output file.
    use ISO_Varying_String 
    use Numerical_Constants_Astronomical
    use Stellar_Population_Properties_Luminosities
    implicit none
    type     (treeNode), intent(inout), pointer      :: thisNode
    double precision   , intent(in   )               :: time
    integer            , intent(inout)               :: integerProperty,doubleProperty
    character(len=*   ), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision   , intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    integer                                          :: iLuminosity

    ! Initialize the module.
    call Galacticus_Output_Tree_Half_Light_Initialize

    ! Return property names if we are outputting half-light data.
    if (outputHalfLightData) then
       do iLuminosity=1,luminositiesCount
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>halfLightRadius</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>Radius enclosing half the galaxy light [Mpc]</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='halfLightRadius'//Stellar_Population_Luminosities_Name(iLuminosity)
          doublePropertyComments(doubleProperty)='Radius enclosing half the galaxy light [Mpc]'
          doublePropertyUnitsSI (doubleProperty)=megaParsec
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>halfLightMass</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>Mass enclosed within the galaxy half-light radius [Solar masses]</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='halfLightMass'//Stellar_Population_Luminosities_Name(iLuminosity)
          doublePropertyComments(doubleProperty)='Mass enclosed within the galaxy half-light radius [Solar masses]'
          doublePropertyUnitsSI (doubleProperty)=massSolar
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Half_Light_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Half_Light_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Half_Light</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Half_Light_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of half-light properties to be written to the \glc\ output file.
    implicit none
    type(treeNode)  , intent(inout), pointer :: thisNode
    double precision, intent(in   )          :: time
    integer         , intent(inout)          :: integerPropertyCount,doublePropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Half_Light_Initialize

    ! Increment property count if we are outputting half-light data.
    if (outputHalfLightData) doublePropertyCount=doublePropertyCount+halfLightPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Half_Light_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Half_Light</unitName>
  !#  <sortName>Galacticus_Output_Tree_Half_Light</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Half_Light(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store density contrast properties in the \glc\ output file buffers.
    use Kind_Numbers 
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)
    integer                                         :: iLuminosity
    double precision                                :: halfLightRadius,massEnclosed

    ! Initialize the module.
    call Galacticus_Output_Tree_Half_Light_Initialize

    ! Store property data if we are outputting half-light data.
    if (outputHalfLightData) then
       
       ! Loop over luminosities.
       do iLuminosity=1,luminositiesCount
          
          ! Get the half-light radius.
          halfLightRadius=Galactic_Structure_Radius_Enclosing_Mass(thisNode,fractionalMass=0.5d0,massType=massTypeStellar&
               &,weightBy=weightByLuminosity,weightIndex=iLuminosity)
          
          ! Find the total mass enclosed.
          massEnclosed   =Galactic_Structure_Enclosed_Mass(thisNode,halfLightRadius,componentType=componentTypeAll,massType&
               &=massTypeAll)

          ! Store the resulting radius.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=halfLightRadius
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=massEnclosed

       end do

    end if
    return
  end subroutine Galacticus_Output_Tree_Half_Light
  
end module Galacticus_Output_Tree_Half_Light_Properties
