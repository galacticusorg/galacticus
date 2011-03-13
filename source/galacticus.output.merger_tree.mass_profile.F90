!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which handles outputting of node mass profiles.

module Galacticus_Output_Tree_Mass_Profiles
  !% Handles outputting of node mass profiles.
  use Tree_Nodes
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

    !$omp critical(Galacticus_Output_Tree_Mass_Profile_Initialize)
    if (.not.outputMassProfileDataInitialized) then
       !@ <inputParameter>
       !@   <name>outputMassProfileData</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not half-light radius data (i.e. radius and mass) should be included in the output.
       !@   </description>
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
          !@ </inputParameter>
          call Get_Input_Parameter('outputMassProfileRadii',outputMassProfileRadii)
       end if

       ! Flag that module is now initialized.
       outputMassProfileDataInitialized=.true.
    end if
    !$omp end critical(Galacticus_Output_Tree_Mass_Profile_Initialize)
    return
  end subroutine Galacticus_Output_Tree_Mass_Profile_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Mass_Profile_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Mass_Profile</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Mass_Profile_Names(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of half-light properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    integer                                       :: iRadius

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
  subroutine Galacticus_Output_Tree_Mass_Profile_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of half-light properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount
    
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
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)
    integer                                         :: iRadius
    double precision                                :: massProfileRadius,massEnclosed

    ! Initialize the module.
    call Galacticus_Output_Tree_Mass_Profile_Initialize

    ! Store property data if we are outputting half-light data.
    if (outputMassProfileData) then
       
       ! Loop over radii at which to output the mass profile.
       do iRadius=1,massProfilePropertyCount
          
          ! Find the total mass enclosed.
          massEnclosed=Galactic_Structure_Enclosed_Mass(thisNode,outputMassProfileRadii(iRadius),componentType=componentTypeAll,massType=massTypeAll)
          
          ! Store the resulting mass.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=massEnclosed
          
       end do

    end if
    return
  end subroutine Galacticus_Output_Tree_Mass_Profile
  
end module Galacticus_Output_Tree_Mass_Profiles
