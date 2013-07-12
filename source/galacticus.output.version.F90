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

!% Contains a module which implements writing of the version number and run time to the \glc\ output file.

module Galacticus_Versioning
  !% Implements writing of the version number and run time to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Version_Output, Galacticus_Version

  ! Define the version.
  integer, parameter :: versionMajor   =0
  integer, parameter :: versionMinor   =9
  integer, parameter :: versionRevision=2

  ! Include the automatically generated Bazaar revision number.
  include 'galacticus.output.version.revision.inc'

contains

  function Galacticus_Version()
    !% Returns a string describing the version of \glc.
    use ISO_Varying_String
    use String_Handling
    implicit none
    type(varying_string) :: Galacticus_Version

    Galacticus_Version="v"
    Galacticus_Version=Galacticus_Version//versionMajor//"."//versionMinor//"."//versionRevision//".r"//hgRevision
    return
  end function Galacticus_Version

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
    use FoX_utils
    implicit none
    type     (Node          ), pointer :: doc            , emailNode, nameNode, thisNode
    type     (NodeList      ), pointer :: nodesList
    integer                            :: ioErr
    character(len=128       )          :: textBufferFixed
    type     (hdf5Object    )          :: versionGroup
    type     (varying_string)          :: runTime

	! Write a UUID for this model.
    call galacticusOutputFile%writeAttribute(generate_UUID(4),'UUID')

    ! Create a group for version information.
    versionGroup=galacticusOutputFile%openGroup('Version','Version and timestamp for this model.')
    call versionGroup%writeAttribute(versionMajor   ,'versionMajor'   )
    call versionGroup%writeAttribute(versionMinor   ,'versionMinor'   )
    call versionGroup%writeAttribute(versionRevision,'versionRevision')
    call versionGroup%writeAttribute(hgRevision     ,'hgRevision'     )
    runTime=Formatted_Date_and_Time()
    call versionGroup%writeAttribute(runTime        ,'runTime'        )

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

end module Galacticus_Versioning
