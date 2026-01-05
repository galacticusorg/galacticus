!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

! Specify an explicit dependence on the git2.o object file.
!: $(BUILDPATH)/git2.o

module Output_Versioning
  !!{
  Implements writing of the version number and run time to the \glc\ output file.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_char, c_size_t
  implicit none
  private
  public :: Version_String, Version

  ! Include the automatically generated Git revision number.
  include 'output.version.revision.inc'

#ifdef GIT2AVAIL
  interface
     subroutine repoHeadHash(repoPath,hash) bind(c,name='repoHeadHash')
       !!{
       Template for a C function that returns the Git hash of a repo HEAD.
       !!}
       import
       character(kind=c_char) :: repoPath
       character(kind=c_char) :: hash    (41)
     end subroutine repoHeadHash
  end interface
#endif

  ! System clock starting count.
  integer(c_size_t) :: countStartClockSystem
  
contains

  subroutine Version(gitHash_,gitBranch_,buildTime_)
    !!{
    Return version information
    !!}
    use :: ISO_Varying_String, only : assignment(=), varying_string
    implicit none
    character(len=42        ), intent(  out), optional :: gitHash_
    type     (varying_string), intent(  out), optional :: gitBranch_  , buildTime_

    if (present(gitHash_   )) gitHash_   =     gitHash
    if (present(gitBranch_ )) gitBranch_ =trim(gitBranch)
    if (present(buildTime_ )) buildTime_ =trim(buildTime)
    return
  end subroutine Version

  function Version_String()
    !!{
    Returns a string describing the version of \glc.
    !!}
    use :: ISO_Varying_String, only : operator(//), var_str, varying_string
    implicit none
    type(varying_string) :: Version_String

    Version_String=var_str("revision ")//gitHash//" (branch: "//trim(gitBranch)//"; build time: "//trim(buildTime)//")"
    return
  end function Version_String

  !![
  <outputFileOpen function="Version_Output"/>
  !!]
  subroutine Version_Output
    !!{
    Output version information to the main output file.
    !!}
#ifdef GIT2AVAIL
    use, intrinsic :: ISO_C_Binding     , only : c_null_char
#else
    use            :: Input_Paths       , only : pathTypeDataDynamic
    use            :: System_Command    , only : System_Command_Do
    use            :: File_Utilities    , only : File_Name_Temporary              , File_Remove
#endif
    use            :: Dates_and_Times   , only : Formatted_Date_and_Time
    use            :: File_Utilities    , only : File_Exists
    use            :: FoX_dom           , only : destroy                          , node              , extractDataContent
    use            :: FoX_utils         , only : generate_UUID
    use            :: Error             , only : Error_Report
    use            :: Output_HDF5       , only : outputFile
    use            :: HDF5_Access       , only : hdf5Access
    use            :: Input_Paths       , only : inputPath                        , pathTypeDataStatic
    use            :: IO_HDF5           , only : hdf5Object
    use            :: IO_XML            , only : XML_Get_First_Element_By_Tag_Name, XML_Path_Exists   , XML_Parse
    use            :: ISO_Varying_String, only : varying_string                   , char
    use            :: String_Handling   , only : String_C_to_Fortran
    implicit none
    type     (Node          ), pointer :: doc            , emailNode, &
         &                                nameNode
    integer                            :: ioErr
    character(len= 41       )          :: gitHashDatasets
    character(len=128       )          :: textBufferFixed
    type     (hdf5Object    )          :: versionGroup
    type     (varying_string)          :: runTime
#ifndef GIT2AVAIL
    integer                            :: status         , hashUnit
    type     (varying_string)          :: hashFileName
#endif

    ! Record the count of the system clock.
    call System_Clock(count=countStartClockSystem)
    
    ! Write a UUID for this model.
    !$ call hdf5Access%set()
    call outputFile%writeAttribute(generate_UUID(4),'UUID')

    ! Create a group for version information.
    runTime     =Formatted_Date_and_Time()
    versionGroup=outputFile%openGroup('Version','Version and timestamp for this model.')
    call versionGroup%writeAttribute(     gitHash   ,'gitHash'     )
    call versionGroup%writeAttribute(trim(gitBranch),'gitBranch'   )
    call versionGroup%writeAttribute(trim(buildTime),'buildTime'   )
    call versionGroup%writeAttribute(     runTime   ,'runStartTime')
#ifdef GIT2AVAIL
    ! Use the git2 library to get the hash of the datasets repo.
    call repoHeadHash(char(inputPath(pathTypeDataStatic))//c_null_char,gitHashDatasets)
    gitHashDatasets=char(String_C_to_Fortran(gitHashDatasets))
#else
    ! Git2 library is not available. If we have the command line `git` installed, use it insted.
    gitHashDatasets="unknown"
    hashFileName   =File_Name_Temporary("repoHash.txt")
    call System_Command_Do("cd "//char(inputPath(pathTypeDataStatic))//"; if which git > /dev/null && git rev-parse --is-inside-work-tree > /dev/null 2>&1 ; then git rev-parse HEAD > "//char(hashFileName)//"; else echo unknown > "//char(hashFileName)//"; fi",iStatus=status)
    if (status == 0) then
       open(newUnit=hashUnit,file=char(hashFileName),status="old",form="formatted",iostat=ioErr)
       if (ioErr == 0) read (hashUnit,'(a)') gitHashDatasets
       close(hashUnit)
    end if
    call File_Remove(hashFileName)
#endif
    call versionGroup%writeAttribute(trim(gitHashDatasets),'gitHashDatasets')
    
    ! Check if a galacticusConfig.xml file exists.
    if (File_Exists("galacticusConfig.xml")) then
       !$omp critical (FoX_DOM_Access)
       doc => XML_Parse("galacticusConfig.xml",iostat=ioErr)
       if (ioErr /= 0) call Error_Report('Unable to parse config file'//{introspection:location})
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
  end subroutine Version_Output
  
  !![
  <outputFileClose function="Version_Finalize"/>
  !!]
  subroutine Version_Finalize()
    !!{
    Output final version information to the main output file.
    !!}
    use :: Dates_and_Times   , only : Formatted_Date_and_Time
    use :: IO_HDF5           , only : hdf5Object
    use :: HDF5_Access       , only : hdf5Access
    use :: Output_HDF5       , only : outputFile
    use :: ISO_Varying_String, only : varying_string
    implicit none
    type            (hdf5Object    ) :: versionGroup
    type            (varying_string) :: runTime
    integer         (c_size_t      ) :: countEndClockSystem, rateClockSystem
    double precision                 :: timeRunDuration

    !$ call hdf5Access%set()
    call System_Clock(count=countEndClockSystem,count_rate=rateClockSystem)
    runTime        =Formatted_Date_and_Time()
    timeRunDuration=dble(countEndClockSystem-countStartClockSystem)/dble(rateClockSystem)
    versionGroup   =outputFile%openGroup('Version')
    call versionGroup%writeAttribute(runTime        ,'runEndTime' )
    call versionGroup%writeAttribute(timeRunDuration,'runDuration')
    call versionGroup%close         (                             )
    !$ call hdf5Access%unset()
    return
  end subroutine Version_Finalize
    
end module Output_Versioning
