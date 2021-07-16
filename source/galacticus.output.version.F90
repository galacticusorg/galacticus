!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!!{
Contains a module which implements writing of the version number and run time to the \glc\ output file.
!!}

module Galacticus_Versioning
  !!{
  Implements writing of the version number and run time to the \glc\ output file.
  !!}
  implicit none
  private
  public :: Galacticus_Version_Output, Galacticus_Version_String, Galacticus_Version

  ! Include the automatically generated Mercurial revision number.
  include 'galacticus.output.version.revision.inc'

contains

  subroutine Galacticus_Version(gitHash_,gitBranch_,buildTime_)
    !!{
    Return version information
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=)
    implicit none
    character(len=42        ), intent(  out), optional :: gitHash_
    type     (varying_string), intent(  out), optional :: gitBranch_  , buildTime_

    if (present(gitHash_   )) gitHash_   =     gitHash
    if (present(gitBranch_ )) gitBranch_ =trim(gitBranch)
    if (present(buildTime_ )) buildTime_ =trim(buildTime)
    return
  end subroutine Galacticus_Version

  function Galacticus_Version_String()
    !!{
    Returns a string describing the version of \glc.
    !!}
    use :: ISO_Varying_String, only : var_str, varying_string, operator(//)
    implicit none
    type(varying_string) :: Galacticus_Version_String

    Galacticus_Version_String=var_str("revision ")//gitHash//" (branch: "//trim(gitBranch)//"; build time: "//trim(buildTime)//")"
    return
  end function Galacticus_Version_String

  !![
  <outputFileOpenTask>
   <unitName>Galacticus_Version_Output</unitName>
  </outputFileOpenTask>
  !!]
  subroutine Galacticus_Version_Output
    !!{
    Output version information to the main output file.
    !!}
    use :: Dates_and_Times   , only : Formatted_Date_and_Time
    use :: File_Utilities    , only : File_Exists
    use :: FoX_dom           , only : destroy                          , node
    use :: FoX_utils         , only : generate_UUID
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: Galacticus_HDF5   , only : galacticusOutputFile
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Path_Exists, XML_Parse, extractDataContent => extractDataContentTS
    use :: ISO_Varying_String, only : varying_string
    implicit none
    type     (Node          ), pointer :: doc            , emailNode, nameNode
    integer                            :: ioErr
    character(len=128       )          :: textBufferFixed
    type     (hdf5Object    )          :: versionGroup
    type     (varying_string)          :: runTime

    ! Write a UUID for this model.
    !$ call hdf5Access%set()
    call galacticusOutputFile%writeAttribute(generate_UUID(4),'UUID')

    ! Create a group for version information.
    runTime     =Formatted_Date_and_Time()
    versionGroup=galacticusOutputFile%openGroup('Version','Version and timestamp for this model.')
    call versionGroup%writeAttribute(     gitHash   ,'gitHash'  )
    call versionGroup%writeAttribute(trim(gitBranch),'gitBranch')
    call versionGroup%writeAttribute(trim(buildTime),'buildTime')
    call versionGroup%writeAttribute(     runTime   ,'runTime'  )

    ! Check if a galacticusConfig.xml file exists.
    if (File_Exists("galacticusConfig.xml")) then
       !$omp critical (FoX_DOM_Access)
       doc => XML_Parse("galacticusConfig.xml",iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Unable to parse config file'//{introspection:location})
       if (XML_Path_Exists(doc,"contact")) then
          if (XML_Path_Exists(doc,"contact/name")) then
             nameNode => XML_Get_First_Element_By_Tag_Name(doc,"contact/name")
             call extractDataContent(nameNode,textBufferFixed)
             call versionGroup%writeAttribute(trim(textBufferFixed),'runByName')
          end if
          if (XML_Path_Exists(doc,"contact/email")) then
             emailNode => XML_Get_First_Element_By_Tag_Name(doc,"contact/email")
             call extractDataContent(emailNode,textBufferFixed)
             call versionGroup%writeAttribute(trim(textBufferFixed),'runByEmail')
          end if
       end if
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
    end if

    ! Close the version group.
    call versionGroup%close()
    !$ call hdf5Access%unset()
    return
  end subroutine Galacticus_Version_Output

end module Galacticus_Versioning
