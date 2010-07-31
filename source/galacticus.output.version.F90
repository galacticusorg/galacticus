!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements writing of the version number and run time to the \glc\ output file.

module Galacticus_Version
  !% Implements writing of the version number and run time to the \glc\ output file.
  private
  public :: Galacticus_Version_Output

contains

  !# <outputFileOpenTask>
  !#  <unitName>Galacticus_Version_Output</unitName>
  !# </outputFileOpenTask>
  subroutine Galacticus_Version_Output
    !% Output version information to the main output file.
    use Galacticus_HDF5_Groups
    use ISO_Varying_String
    use HDF5
    use Galacticus_Error
    use Dates_and_Times
    use File_Utilities
    use FoX_dom
    implicit none
    include 'galacticus.output.version.revision.inc'
    type(varying_string), dimension(1) :: runTime,textBufferVariable
    type(Node),           pointer      :: doc,thisNode,nameNode,emailNode
    type(NodeList),       pointer      :: nodesList
    integer                            :: ioErr
    integer(kind=HID_T)                :: versionGroupID=0,versionMajorID=0,versionMinorID=0,versionRevisionID=0,runTimeID=0&
         &,nameID=0,emailID=0,bazaarRevisionID=0
    type(varying_string)               :: groupName,groupComment
    character(len=128)                 :: textBufferFixed
  
    ! Create a group for version information.
    groupName='Version'
    groupComment='Version and timestamp for this model.'
    versionGroupID=Galacticus_Output_Make_Group(groupName,groupComment)
    call Galacticus_Output_Dataset(versionGroupID,versionMajorID   ,'versionMajor'   ,'Major version number'  ,[0]             )
    call Galacticus_Output_Dataset(versionGroupID,versionMinorID   ,'versionMinor'   ,'Minor version number'  ,[9]             )
    call Galacticus_Output_Dataset(versionGroupID,versionRevisionID,'versionRevision','Revision number'       ,[0]             )
    call Galacticus_Output_Dataset(versionGroupID,bazaarRevisionID ,'bazaarRevision' ,'Bazaar revision number',[bazaarRevision])
    runTime(1)=Formatted_Date_and_Time()
    call Galacticus_Output_Dataset(versionGroupID,runTimeID        ,'runTime'        ,'Time at which model was run',runTime)

    ! Check if a galacticusConfig.xml file exists.
    if (File_Exists("galacticusConfig.xml")) then
       !$omp critical (FoX_DOM_Access)
       doc => parseFile("galacticusConfig.xml",iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Galacticus_Version_Output','Unable to parse config file')
       nodesList => getElementsByTagname(doc,"contact")
       if (getLength(nodesList) >= 0) then
          thisNode => item(nodesList,0)
          nodesList => getElementsByTagname(thisNode,"name")
          if (getLength(nodesList) >= 0) then
             nameNode => item(nodesList,0)
             call extractDataContent(nameNode,textBufferFixed)
             textBufferVariable(1)=trim(textBufferFixed)
             call Galacticus_Output_Dataset(versionGroupID,nameID,'runByName','The name of whosoever ran this model',textBufferVariable)
          end if
          nodesList => getElementsByTagname(thisNode,"email")
          if (getLength(nodesList) >= 0) then
             emailNode => item(nodesList,0)
             call extractDataContent(emailNode,textBufferFixed)
             textBufferVariable(1)=trim(textBufferFixed)
             call Galacticus_Output_Dataset(versionGroupID,emailID,'runByEmail','The e-mail address of whosoever ran this model',textBufferVariable)
          end if
       end if
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine Galacticus_Version_Output
  
end module Galacticus_Version
