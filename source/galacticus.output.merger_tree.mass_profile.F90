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

!% Contains a module which handles outputting of node mass profiles.

module Galacticus_Output_Tree_Mass_Profiles
  !% Handles outputting of node mass profiles.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Mass_Profile, Galacticus_Output_Tree_Mass_Profile_Property_Count, Galacticus_Output_Tree_Mass_Profile_Names

  ! Flag indicating if this module is initialize.
  logical                                     :: outputMassProfileDataInitialized

  ! Number of properties.
  integer                                     :: massProfilePropertyCount

  ! Array of radii.
  double precision, allocatable, dimension(:) :: outputMassProfileRadii

  ! Flag indicating whether or not half-light data is to be output.
  logical                                     :: outputMassProfileData

contains

  subroutine Galacticus_Output_Tree_Mass_Profile_Initialize
    !% Initializes the module by determining whether or not half-light radius data should be output.
    use Input_Parameters
    use Memory_Management
    implicit none

    if (.not.outputMassProfileDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Mass_Profile_Initialize)
       if (.not.outputMassProfileDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputMassProfileData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not half-light radius data (i.e. radius and mass) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputMassProfileData',outputMassProfileData,defaultValue=.false.)

          ! Read radii if necessary.
          if (outputMassProfileData) then
             massProfilePropertyCount=Get_Input_Parameter_Array_Size('outputMassProfileRadii')
             call Alloc_Array(outputMassProfileRadii,[massProfilePropertyCount])
             !@ <inputParameter>
             !@   <name>outputMassProfileRadii</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     A list of radii at which to output the mass profile.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1..*</cardinality>
             !@   <group>output</group>
             !@ </inputParameter>
             call Get_Input_Parameter('outputMassProfileRadii',outputMassProfileRadii)
          end if

          ! Flag that module is now initialized.
          outputMassProfileDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Mass_Profile_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Mass_Profile_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Mass_Profile_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Mass_Profile</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Mass_Profile_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of half-light properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                                          :: iRadius

    ! Initialize the module.
    call Galacticus_Output_Tree_Mass_Profile_Initialize

    ! Return property names if we are outputting half-light data.
    if (outputMassProfileData) then
       do iRadius=1,massProfilePropertyCount
          doubleProperty=doubleProperty+1
          write (doublePropertyNames   (doubleProperty),'(a,e9.3  )') 'massProfile'                      ,outputMassProfileRadii(iRadius)
          write (doublePropertyComments(doubleProperty),'(a,e9.3,a)') 'Mass enclosed within a radius of ',outputMassProfileRadii(iRadius),' Mpc [Solar masses]'
          doublePropertyUnitsSI (doubleProperty)=massSolar
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Mass_Profile_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Mass_Profile_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Mass_Profile</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Mass_Profile_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of half-light properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Mass_Profile_Initialize

    ! Increment property count if we are outputting half-light data.
    if (outputMassProfileData) doublePropertyCount=doublePropertyCount+massProfilePropertyCount
    return
  end subroutine Galacticus_Output_Tree_Mass_Profile_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Mass_Profile</unitName>
  !#  <sortName>Galacticus_Output_Tree_Mass_Profile</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Mass_Profile(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store density contrast properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)
    integer                                                  :: iRadius
    double precision                                         :: massEnclosed

    ! Initialize the module.
    call Galacticus_Output_Tree_Mass_Profile_Initialize

    ! Store property data if we are outputting half-light data.
    if (outputMassProfileData) then

       ! Loop over radii at which to output the mass profile.
       do iRadius=1,massProfilePropertyCount

          ! Find the total mass enclosed.
          massEnclosed=Galactic_Structure_Enclosed_Mass(thisNode,outputMassProfileRadii(iRadius),componentType=componentTypeAll,massType=massTypeAll,haloLoaded=.true.)

          ! Store the resulting mass.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=massEnclosed

       end do

    end if
    return
  end subroutine Galacticus_Output_Tree_Mass_Profile

end module Galacticus_Output_Tree_Mass_Profiles
