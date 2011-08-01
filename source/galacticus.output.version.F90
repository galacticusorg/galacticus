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


!% Contains a module which implements writing of the version number and run time to the \glc\ output file.

module Galacticus_Version
  !% Implements writing of the version number and run time to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Version_Output

contains

  !# <outputFileOpenTask>
  !#  <unitName>Galacticus_Version_Output</unitName>
  !# </outputFileOpenTask>
  subroutine Galacticus_Version_Output
    !% Output version information to the main output file.
    use Galacticus_HDF5
    use IO_HDF5
    use ISO_Varying_String
    use Galacticus_Error
    use Dates_and_Times
    use File_Utilities
    use FoX_dom
    implicit none
    include 'galacticus.output.version.revision.inc'
    type(Node),           pointer      :: doc,thisNode,nameNode,emailNode
    type(NodeList),       pointer      :: nodesList
    integer                            :: ioErr
    character(len=128)                 :: textBufferFixed
    type(hdf5Object)                   :: versionGroup
    type(varying_string)               :: runTime

    ! Create a group for version information.
    versionGroup=IO_HDF5_Open_Group(galacticusOutputFile,'Version','Version and timestamp for this model.')
!    call versionGroup%writeAttribute(0             ,'versionMajor'   )
!    call versionGroup%writeAttribute(9             ,'versionMinor'   )
!    call versionGroup%writeAttribute(0             ,'versionRevision')
!    call versionGroup%writeAttribute(bazaarRevision,'bazaarRevision' )
    runTime=Formatted_Date_and_Time()
    call versionGroup%writeAttribute(runTime       ,'runTime'        )

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
             call versionGroup%writeAttribute(trim(textBufferFixed),'runByName')
          end if
          nodesList => getElementsByTagname(thisNode,"email")
          if (getLength(nodesList) >= 0) then
             emailNode => item(nodesList,0)
             call extractDataContent(emailNode,textBufferFixed)
             call versionGroup%writeAttribute(trim(textBufferFixed),'runByEmail')
          end if
       end if
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
    end if

    ! Close the version group.
    call versionGroup%close()
    return
  end subroutine Galacticus_Version_Output
  
end module Galacticus_Version
