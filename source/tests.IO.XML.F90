!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a program which tests functionality of the XML I/O module.
!!}

program Tests_IO_XML
  !!{
  Tests the XML I/O module.
  !!}
  use            :: Unit_Tests    , only : Assert                        , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  use            :: IO_XML        , only : XML_Count_Elements_By_Tag_Name                        , XML_Array_Read      , XML_Array_Read_Static, XML_Get_First_Element_By_Tag_Name, &
          &                                XML_Parse                                             , XML_Path_Exists     , xmlNodeList          , XML_Get_ELements_By_Tag_Name
  use            :: Display       , only : displayVerbositySet           , verbosityLevelStandard
  use            :: FoX_DOM       , only : destroy                       , node                  , serialize           , extractDataContent
  use            :: Error         , only : Error_Report
  use, intrinsic :: ISO_C_Binding , only : c_size_t
  use            :: System_Command, only : System_Command_Do
  implicit none
  type            (node       )                           , pointer :: doc        , xmlElement
  type            (xmlNodeList), allocatable, dimension(:)          :: xmlElements
  double precision             , allocatable, dimension(:)          :: array1     , array2
  integer                      , allocatable, dimension(:)          :: iarray1
  character       (len=1      ), allocatable, dimension(:)          :: carray1
  integer                                   , dimension(1)          :: iValue
  integer                                                           :: ioErr      , status

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("XML I/O")

  !$omp critical (FoX_DOM_Access)
  ! Parse the XML file.
  doc => XML_Parse("testSuite/data/xmlTest.xml",iostat=ioErr)
  if (ioErr /= 0) call Error_Report("failed to parse 'testSuite/data/xmlTest.xml'"//{introspection:location})

  ! Test XInclude.
  call serialize(doc,"testSuite/outputs/xmlTest.xml")
  call System_Command_Do("diff testSuite/outputs/xmlTest.xml testSuite/data/xmlTestIncluded.xml",status)
  call Assert("Parse XIncludes",status,0)

  ! Test array reading.
  xmlElement => XML_Get_First_Element_By_Tag_Name(doc,"array1")
  call XML_Array_Read       (xmlElement,"datum",array1       )
  call Assert("Read 1D allocatable array from element",array1,[0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0],absTol=1.0d-6)
  call XML_Array_Read_Static(xmlElement,"datum",array1       )
  call Assert("Read 1D static array from element"     ,array1,[0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0],absTol=1.0d-6)
  call XML_Get_Elements_By_Tag_Name(doc,"array",xmlElements)
  call XML_Array_Read       (xmlElements,"value",array1       )
  call Assert("Read 1D allocatable array from list"   ,array1,[0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0],absTol=1.0d-6)
  call XML_Array_Read_Static(xmlElements,"value",array1       )
  call Assert("Read 1D static array from list"        ,array1,[0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0],absTol=1.0d-6)
  allocate(iarray1(10))
  call XML_Get_Elements_By_Tag_Name(doc,"iarray",xmlElements)
  call XML_Array_Read_Static(xmlElements,"value",iarray1       )
  call Assert("Read 1D static integer array from list"        ,iarray1,[0,1,2,3,4,5,6,7,8,9])
  allocate(carray1(10))
  call XML_Get_Elements_By_Tag_Name(doc,"carray",xmlElements)
  call XML_Array_Read_Static(xmlElements,"value",carray1       )
  call Assert("Read 1D static integer array from list"        ,carray1,["0","1","2","3","4","5","6","7","8","9"])
  xmlElement => XML_Get_First_Element_By_Tag_Name(doc,"array2")
  call XML_Array_Read       (xmlElement,"datum",array1,array2)
  call Assert("Read 2D allocatable array from element",[array1,array2],[0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,9.0d0,8.0d0,7.0d0,6.0d0,5.0d0,4.0d0,3.0d0,2.0d0,1.0d0,0.0d0],absTol=1.0d-6)
  xmlElement => XML_Get_First_Element_By_Tag_Name(doc,"test/level1/level2/level3",directChildrenOnly=.true.)
  call extractDataContent(xmlElement,iValue)
  call Assert("Get first element at end of path",iValue,[1234])
  xmlElement => XML_Get_First_Element_By_Tag_Name(doc,"test/some/path/to[@value='scoobydoo']/find",directChildrenOnly=.true.)
  call extractDataContent(xmlElement,iValue)
  call Assert("Process XPath expression",iValue,[9275932])

  ! Test array length functions.
  xmlElement => XML_Get_First_Element_By_Tag_Name(doc,"array1")
  call Assert("Determine array length",XML_Count_Elements_By_Tag_Name(xmlElement,"datum"),10_c_size_t)

  ! Test path detection.
  call Assert("Extant path correctly detected"    ,XML_Path_Exists(doc,"test/level1/level2/level3"),.true. )
  call Assert("Non-extant path correctly detected",XML_Path_Exists(doc,"test/level1/level4/level3"),.false.)

  ! Destroy the XML document.
  call destroy(doc)
  !$omp end critical (FoX_DOM_Access)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_IO_XML
