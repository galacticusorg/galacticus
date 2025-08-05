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

program Tests_IO_HDF5
  !!{
  Tests the HDF5 I/O module.
  !!}
  use :: Display           , only : displayVerbositySet, verbosityLevelStandard
  use :: Error             , only : Error_Report
  use :: HDF5              , only : HSIZE_T
  use :: IO_HDF5           , only : IO_HDF5_Is_HDF5    , hdf5Object            , hdf5VarDouble       , hdf5VarInteger8  , &
       &                            hdf5VarDouble2D
  use :: ISO_Varying_String, only : assignment(=)      , trim                  , varying_string      , var_str
  use :: Kind_Numbers      , only : kind_int8
  use :: System_Command    , only : System_Command_Do
  use :: Unit_Tests        , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                                                                   :: iPass                          , integerValue         , &
       &                                                                       integerValueReread             , i                    , &
       &                                                                       j                              , k
  logical                                                                   :: appendableOK
  integer                                       , dimension(10)             :: integerValueArray
  integer                                       , dimension(10)             :: integerValueArrayRereadStatic
  integer                          , allocatable, dimension( :)             :: integerValueArrayReread
  integer         (kind=kind_int8 )                                         :: integer8Value                  , integer8ValueReread
  integer         (kind=kind_int8 )             , dimension(10)             :: integer8ValueArray
  integer         (kind=kind_int8 )             , dimension(10)             :: integer8ValueArrayRereadStatic
  integer         (kind=kind_int8 ), allocatable, dimension( :)             :: integer8ValueArrayReread       , integerRangeU32
  double precision                                                          :: doubleValue                    , doubleValueReread
  double precision                              , dimension(10)             :: doubleValueArray
  double precision                              , dimension(10)             :: doubleValueArrayRereadStatic
  double precision                 , allocatable, dimension( :)             :: doubleValueArrayReread
  character       (len=32         )                                         :: characterValue                 , characterValueReread
  character       (len=32         )             , dimension(10)             :: characterValueArray
  character       (len=32         )             , dimension(10)             :: characterValueArrayRereadStatic
  character       (len=32         ), allocatable, dimension( :)             :: characterValueArrayReread
  type            (varying_string )                                         :: varStringValue                 , varStringValueReread
  type            (varying_string )             , dimension(10)             :: varStringValueArray
  type            (varying_string )             , dimension(10)             :: varStringValueArrayRereadStatic
  type            (varying_string ), allocatable, dimension( :)             :: varStringValueArrayReread
  double precision                              , dimension(10,10)          :: doubleValueArray2d
  double precision                              , dimension(10,10)          :: doubleValueArray2dRereadStatic
  double precision                 , allocatable, dimension( :, :)          :: doubleValueArray2dReread      , doubleValueArray2dRereadExpect
  double precision                              , dimension(10,10,10)       :: doubleValueArray3d
  double precision                              , dimension(10,10,10)       :: doubleValueArray3dRereadStatic
  double precision                 , allocatable, dimension( :, :, :)       :: doubleValueArray3dReread      , doubleValueArray3dRereadExpect
  double precision                              , dimension(10,10,10,10)    :: doubleValueArray4d
  double precision                              , dimension(10,10,10,10)    :: doubleValueArray4dRereadStatic
  double precision                 , allocatable, dimension( :, :, :, :)    :: doubleValueArray4dReread      , doubleValueArray4dRereadExpect
  double precision                              , dimension(10,10,10,10,10) :: doubleValueArray5d
  double precision                              , dimension(10,10,10,10,10) :: doubleValueArray5dRereadStatic
  double precision                 , allocatable, dimension( :, :, :, :, :) :: doubleValueArray5dReread      , doubleValueArray5dRereadExpect
  type            (varying_string )             , dimension(20)             :: datasetNamesReference
  type            (hdf5VarDouble  ), allocatable, dimension( :)             :: varDoubleArray2D              , varDoubleDataset2DArrayReread
  type            (hdf5VarDouble2D), allocatable, dimension( :)             :: varDoubleArray3D              , varDoubleDataset3DArrayReread
  type            (hdf5VarInteger8), allocatable, dimension( :)             :: varInteger8Array2D            , varInteger8Dataset2dArrayReread

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("HDF5 IO")

  ! Begin passes through the tests.
  do iPass=1,2
     ! Open an HDF5 file.
     block
       type(hdf5Object) :: fileObject, groupObject
       fileObject=hdf5Object("testSuite/outputs/test.IO.HDF5.hdf5",overWrite=.true.,objectsOverwritable=.true.,useLatestFormat=.true.)
       ! Open an HDF5 group.
       select case (iPass)
       case(1)
          call Unit_Tests_Begin_Group("Tests with chunking enabled")
          groupObject=fileObject%openGroup("myGroup",comment="This is my group.",objectsOverwritable=.true.,chunkSize=1024_hsize_t,compressionLevel=+9)
          appendableOK=.true.
       case (2)
          call Unit_Tests_Begin_Group("Tests with chunking disabled")
          groupObject=fileObject%openGroup("myGroup",comment="This is my group.",objectsOverwritable=.true.,chunkSize=  -1_hsize_t,compressionLevel=-1)
          appendableOK=.false.
       case default
          appendableOK=.false.
          call Error_Report('unknown pass'//{introspection:location})
       end select
       
       ! Write a scalar integer attribute to the group.
       integerValue=9
       call groupObject%writeAttribute(integerValue,"integerAttribute")
       ! Overwrite a scalar integer attribute to the group.
       integerValue=10
       call groupObject%writeAttribute(integerValue,"integerAttribute")
       ! Read the scalar integer attribute back.
       call groupObject%readAttribute("integerAttribute",integerValueReread)
       call Assert("re-read scalar integer attribute",integerValue,integerValueReread)

       ! Write a pseudo-scalar integer attribute to the group.
       integerValueArray=7
       call groupObject%writeAttribute(integerValueArray(1:1),"integerAttributePseudoScalar")
       ! Read the pseudo-scalar integer attribute back.
       call groupObject%readAttribute("integerAttributePseudoScalar",integerValueReread,allowPseudoScalar=.true.)
       call Assert("re-read pseudo-scalar integer attribute",integerValueArray(1),integerValueReread)

       ! Write an integer 1-D array attribute to the group.
       integerValueArray=7
       call groupObject%writeAttribute(integerValueArray,"integerAttribute1dArray")
       ! Write an integer 1-D array attribute to the group.
       integerValueArray=6
       call groupObject%writeAttribute(integerValueArray,"integerAttribute1dArray")
       ! Read the integer 1-D array attribute back.
       call groupObject%readAttribute("integerAttribute1dArray",integerValueArrayReread)
       call Assert("re-read 1-D array integer attribute",integerValueArray,integerValueArrayReread)
       ! Read the integer 1-D array attribute back to a static array.
       call groupObject%readAttributeStatic("integerAttribute1dArray",integerValueArrayRereadStatic)
       call Assert("re-read 1-D array integer attribute to static array",integerValueArray,integerValueArrayRereadStatic)

       ! Write a scalar long integer attribute to the group.
       integer8Value=20202020
       call groupObject%writeAttribute(integer8value,"integer8Attribute")
       ! Overwrite a scalar long integer attribute to the group.
       integer8Value=20202021
       call groupObject%writeAttribute(integer8value,"integer8Attribute")
       ! Read the scalar long integer attribute back.
       call groupObject%readAttribute("integer8Attribute",integer8ValueReread)
       call Assert("re-read scalar long integer attribute",integer8Value,integer8ValueReread)

       ! Write a pseudo-scalar long integer attribute to the group.
       integer8ValueArray=20202021
       call groupObject%writeAttribute(integer8ValueArray(1:1),"integer8AttributePseudoScalar")
       ! Read the scalar long integer attribute back.
       call groupObject%readAttribute("integer8AttributePseudoScalar",integer8ValueReread,allowPseudoScalar=.true.)
       call Assert("re-read pseudo-scalar long integer attribute",integer8ValueArray(1),integer8ValueReread)

       ! Write a long integer 1-D array attribute to the group.
       integer8ValueArray=7
       call groupObject%writeAttribute(integer8ValueArray,"integer8Attribute1dArray")
       ! Write a long integer 1-D array attribute to the group.
       integer8ValueArray=6
       call groupObject%writeAttribute(integer8ValueArray,"integer8Attribute1dArray")
       ! Read the long integer 1-D array attribute back.
       call groupObject%readAttribute("integer8Attribute1dArray",integer8ValueArrayReread)
       call Assert("re-read 1-D array long integer attribute",integer8ValueArray,integer8ValueArrayReread)
       ! Read the long integer 1-D array attribute back to a static array.
       call groupObject%readAttributeStatic("integer8Attribute1dArray",integer8ValueArrayRereadStatic)
       call Assert("re-read 1-D array long integer attribute to static array",integer8ValueArray,integer8ValueArrayRereadStatic)

       ! Write a scalar double attribute to the group.
       doubleValue=9.12345d0
       call groupObject%writeAttribute(doubleValue,"doubleAttribute")
       ! Overwrite a scalar double attribute to the group.
       doubleValue=10.9876d0
       call groupObject%writeAttribute(doubleValue,"doubleAttribute")
       ! Read the scalar double attribute back.
       call groupObject%readAttribute("doubleAttribute",doubleValueReread)
       call Assert("re-read scalar double attribute",doubleValue,doubleValueReread)

       ! Write a pseudo-scalar double attribute to the group.
       doubleValueArray=9.12345d0
       call groupObject%writeAttribute(doubleValueArray(1:1),"doubleAttributePseudoScalar")
       ! Read the scalar double attribute back.
       call groupObject%readAttribute("doubleAttributePseudoScalar",doubleValueReread,allowPseudoScalar=.true.)
       call Assert("re-read pseudo-scalar double attribute",doubleValueArray(1),doubleValueReread)

       ! Write an double 1-D array attribute to the group.
       doubleValueArray=7.676767d0
       call groupObject%writeAttribute(doubleValueArray,"doubleAttribute1dArray")
       ! Write an double 1-D array attribute to the group.
       doubleValueArray=6.141414d0
       call groupObject%writeAttribute(doubleValueArray,"doubleAttribute1dArray")
       ! Read the double 1-D array attribute back.
       call groupObject%readAttribute("doubleAttribute1dArray",doubleValueArrayReread)
       call Assert("re-read 1-D array double attribute",doubleValueArray,doubleValueArrayReread)
       ! Read the double 1-D array attribute back to a static array.
       call groupObject%readAttributeStatic("doubleAttribute1dArray",doubleValueArrayRereadStatic)
       call Assert("re-read 1-D array double attribute to static array",doubleValueArray,doubleValueArrayRereadStatic)

       ! Write a scalar character attribute to the group.
       characterValue='abcdefghijklmnopqrstuvwxyz'
       call groupObject%writeAttribute(characterValue,"characterAttribute")
       ! Overwrite a scalar character attribute to the group.
       characterValue='qwertyuiop'
       call groupObject%writeAttribute(characterValue,"characterAttribute")
       ! Read the scalar character attribute back.
       call groupObject%readAttribute("characterAttribute",characterValueReread)
       call Assert("re-read scalar character attribute",characterValue,characterValueReread)

       ! Write a pseudo-scalar character attribute to the group.
       characterValueArray='abcdefghijklmnopqrstuvwxyz'
       call groupObject%writeAttribute(characterValueArray(1:1),"characterAttributePseudoScalar")
       ! Read the scalar character attribute back.
       call groupObject%readAttribute("characterAttributePseudoScalar",characterValueReread,allowPseudoScalar=.true.)
       call Assert("re-read pseudo-scalar character attribute",characterValueArray(1),characterValueReread)

       ! Write a character 1-D array attribute to the group.
       characterValueArray='aAbBcCdDeEfFgGhH'
       call groupObject%writeAttribute(characterValueArray,"characterAttribute1dArray")
       ! Write a character 1-D array attribute to the group.
       characterValueArray='1!2@3#4$5%6^7&8*9(0)'
       call groupObject%writeAttribute(characterValueArray,"characterAttribute1dArray")
       ! Read the character 1-D array attribute back.
       call groupObject%readAttribute("characterAttribute1dArray",characterValueArrayReread)
       call Assert("re-read 1-D array character attribute",characterValueArray,characterValueArrayReread)
       ! Read the character 1-D array attribute back to a static array.
       call groupObject%readAttributeStatic("characterAttribute1dArray",characterValueArrayRereadStatic)
       call Assert("re-read 1-D array character attribute to static array",characterValueArray,characterValueArrayRereadStatic)

       ! Write a scalar varString attribute to the group.
       varStringValue="le sange est dans l'arbre. pad pad pad pad!!"
       call groupObject%writeAttribute(varStringValue,"varStringAttribute")
       ! Overwrite a scalar varString attribute to the group.
       varStringValue="dans l'interior de l'intestine de la mouton!"
       call groupObject%writeAttribute(varStringValue,"varStringAttribute")
       ! Read the scalar varString attribute back.
       call groupObject%readAttribute("varStringAttribute",varStringValueReread)
       call Assert("re-read scalar varString attribute",varStringValue,varStringValueReread)
       ! Write a varying string 1-D array attribute to the group.
       varStringValueArray='qazwsxedcrfvtgbyhnujmikolp'
       call groupObject%writeAttribute(varStringValueArray,"varStringAttribute1dArray")
       ! Read the varying string 1-D array attribute back.
       call groupObject%readAttribute("varStringAttribute1dArray",varStringValueArrayReread)
       call Assert("re-read 1-D array varString attribute",varStringValueArray,varStringValueArrayReread)
       ! Read the varying string 1-D array attribute back to a static array.
       call groupObject%readAttributeStatic("varStringAttribute1dArray",varStringValueArrayRereadStatic)
       call Assert("re-read 1-D array varString attribute to static array",varStringValueArray,varStringValueArrayRereadStatic)

       ! Write a pseudo-scalar varString attribute to the group.
       varStringValueArray="le sange est dans l'arbre. pad pad pad pad!!"
       call groupObject%writeAttribute(varStringValueArray(1:1),"varStringAttributePseudoScalar")
       ! Read the scalar varString attribute back.
       call groupObject%readAttribute("varStringAttributePseudoScalar",varStringValueReread,allowPseudoScalar=.true.)
       call Assert("re-read pseudo-scalar varString attribute",varStringValueArray(1),varStringValueReread)

       ! Write an extensible integer 1-D array dataset to the group.
       if (appendableOK) then
          integerValueArray=[0,11,22,33,44,55,66,77,88,99]
          call groupObject%writeDataset(integerValueArray,"integerDataset1dArrayExtensible","This is an extensible dataset",appendTo=.true.)
          ! Append to the extensible integer 1-D array dataset.
          integerValueArray=[1,-2,3,-4,5,-6,7,-8,9,-10]
          call groupObject%writeDataset(integerValueArray,"integerDataset1dArrayExtensible",appendTo=.true.)
       end if

       ! Write an integer 1-D array dataset to the group.
       integerValueArray=[0,11,22,33,44,55,66,77,88,99]
       call groupObject%writeDataset(integerValueArray,"integerDataset1dArray","This is an example dataset")
       ! Write an integer 1-D array dataset to the group.
       integerValueArray=[1,-2,3,-4,5,-6,7,-8,9,-10]
       call groupObject%writeDataset(integerValueArray,"integerDataset1dArray")
       ! Read an integer 1-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("integerDataset1dArray",integerValueArrayRereadStatic)
       call Assert("re-read 1-D array integer dataset to static array",integerValueArray,integerValueArrayRereadStatic)
       ! Read an integer 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("integerDataset1dArray",integerValueArrayReread)
       call Assert("re-read 1-D array integer dataset to allocatable array",integerValueArray,integerValueArrayReread)
       deallocate(integerValueArrayReread)
       ! Read part of an integer 1-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("integerDataset1dArray",integerValueArrayRereadStatic,int([3],kind=HSIZE_T),int([4],kind=HSIZE_T))
       call Assert("re-read part of a 1-D array integer dataset to static array",integerValueArray(3:6),integerValueArrayRereadStatic)
       ! Read part of an integer 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("integerDataset1dArray",integerValueArrayReread,int([3],kind=HSIZE_T),int([4],kind=HSIZE_T))
       call Assert("re-read part of a 1-D array integer dataset to allocatable array",integerValueArray(3:6),integerValueArrayReread)
       deallocate(integerValueArrayReread)
       block
         type(hdf5Object) :: datasetObject
         ! Open the dataset.
         datasetObject=groupObject%openDataset("integerDataset1dArray")
         ! Create a reference to the dataset object.
         call groupObject%createReference1D(datasetObject,"myReference",int([3],kind=HSIZE_T),int([2],kind=HSIZE_T))
       end block
       ! Read an integer 1-D array dataset from the group into a static array.
       integerValueArrayRereadStatic=0
       call groupObject%readDatasetStatic("myReference",integerValueArrayRereadStatic)
       call Assert("re-read referenced 1-D array integer dataset to static array",integerValueArray(3:4),integerValueArrayRereadStatic(1:2))
       ! Read an integer 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("myReference",integerValueArrayReread)
       call Assert("re-read referenced 1-D array integer dataset to allocatable array",integerValueArray(3:4),integerValueArrayReread(1:2))
       deallocate(integerValueArrayReread)
       ! Read part of an integer 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("myReference",integerValueArrayReread,int([2],kind=HSIZE_T),int([1],kind=HSIZE_T))
       call Assert("re-read part of referenced 1-D array integer dataset to allocatable array",integerValueArray(4:4),integerValueArrayReread(1:1))
       deallocate(integerValueArrayReread)

       ! Write an extensible long integer 1-D array dataset to the group.
       if (appendableOK) then
          integer8ValueArray=[0,11,22,33,44,55,66,77,88,99]
          call groupObject%writeDataset(integer8ValueArray,"integer8Dataset1dArrayExtensible","This is an extensible dataset",appendTo=.true.)
          ! Append to the extensible long integer 1-D array dataset.
          integer8ValueArray=[1,-2,3,-4,5,-6,7,-8,9,-10]
          call groupObject%writeDataset(integer8ValueArray,"integer8Dataset1dArrayExtensible",appendTo=.true.)
       end if

       ! Write a long integer 1-D array dataset to the group.
       integer8ValueArray=[0,11,22,33,44,55,66,77,88,99]
       call groupObject%writeDataset(integer8ValueArray,"integer8Dataset1dArray","This is an example dataset")
       ! Write a long integer 1-D array dataset to the group.
       integer8ValueArray=[1,-2,3,-4,5,-6,7,-8,9,-10]
       call groupObject%writeDataset(integer8ValueArray,"integer8Dataset1dArray")
       ! Read a long integer 1-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("integer8Dataset1dArray",integer8ValueArrayRereadStatic)
       call Assert("re-read 1-D array long integer dataset to static array",integer8ValueArray,integer8ValueArrayRereadStatic)
       ! Read a long integer 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("integer8Dataset1dArray",integer8ValueArrayReread)
       call Assert("re-read 1-D array long integer dataset to allocatable array",integer8ValueArray,integer8ValueArrayReread)
       ! Read part of a long integer 1-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("integer8Dataset1dArray",integer8ValueArrayRereadStatic,int([3],kind=HSIZE_T),int([4],kind=HSIZE_T))
       call Assert("re-read part of a 1-D array long integer dataset to static array",integer8ValueArray(3:6),integer8ValueArrayRereadStatic)
       ! Read part of a long integer 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("integer8Dataset1dArray",integer8ValueArrayReread,int([3],kind=HSIZE_T),int([4],kind=HSIZE_T))
       call Assert("re-read part of a 1-D array long integer dataset to allocatable array",integer8ValueArray(3:6),integer8ValueArrayReread)

       block
         type(hdf5Object) :: datasetObject
         ! Open the dataset.
         datasetObject=groupObject%openDataset("integer8Dataset1dArray")
         ! Create a reference to the dataset object.
         call groupObject%createReference1D(datasetObject,"anotherReference",int([3],kind=HSIZE_T),int([2],kind=HSIZE_T))
       end block
       ! Read a long integer 1-D array dataset from the group into a static array.
       integer8ValueArrayRereadStatic=0
       call groupObject%readDatasetStatic("anotherReference",integer8ValueArrayRereadStatic)
       call Assert("re-read referenced 1-D array long integer dataset to static array",integer8ValueArray(3:4),integer8ValueArrayRereadStatic(1:2))
       ! Read a long integer 1-D array dataset from the group into a static array.
       integer8ValueArrayReread=0
       call groupObject%readDataset("anotherReference",integer8ValueArrayReread)
       call Assert("re-read referenced 1-D array long integer dataset to allocatable array",integer8ValueArray(3:4),integer8ValueArrayReread(1:2))
       ! Read part of a long integer 1-D array dataset from the group into an allocatable array.
       integer8ValueArrayReread=0
       call groupObject%readDataset("anotherReference",integer8ValueArrayReread,int([2],kind=HSIZE_T),int([1],kind=HSIZE_T))
       call Assert("re-read part of referenced 1-D array long integer dataset to allocatable array",integer8ValueArray(4:4),integer8ValueArrayReread(1:1))

       ! Write an extensible double 1-D array dataset to the group.
       if (appendableOK) then
          doubleValueArray=[0.0d0,11.0d0,22.0d0,33.0d0,44.0d0,55.0d0,6.0d06,77.0d0,88.0d0,99.0d0]
          call groupObject%writeDataset(doubleValueArray,"doubleDataset1dArrayExtensible","This is an extensible dataset",appendTo=.true.)
          ! Append to the extensible double 1-D array dataset.
          doubleValueArray=[1.0d0,-.0d02,3.0d0,-4.0d0,5.0d0,-6.0d0,7.0d0,-8.0d0,9.0d0,-10.0d0]
          call groupObject%writeDataset(doubleValueArray,"doubleDataset1dArrayExtensible",appendTo=.true.)
       end if

       ! Write a double 1-D array dataset to the group.
       doubleValueArray=[0.0d0,11.0d0,22.0d0,33.0d0,44.0d0,55.0d0,6.0d06,77.0d0,88.0d0,99.0d0]
       call groupObject%writeDataset(doubleValueArray,"doubleDataset1dArray","This is an example dataset")
       ! Write a double 1-D array dataset to the group.
       doubleValueArray=[1.0d0,-.0d02,3.0d0,-4.0d0,5.0d0,-6.0d0,7.0d0,-8.0d0,9.0d0,-10.0d0]
       call groupObject%writeDataset(doubleValueArray,"doubleDataset1dArray")
       ! Read a double 1-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset1dArray",doubleValueArrayRereadStatic)
       call Assert("re-read 1-D array double dataset to static array",doubleValueArray,doubleValueArrayRereadStatic)
       ! Read a double 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset1dArray",doubleValueArrayReread)
       call Assert("re-read 1-D array double dataset to allocatable array",doubleValueArray,doubleValueArrayReread)
       deallocate(doubleValueArrayReread)
       ! Read part of a double 1-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset1dArray",doubleValueArrayRereadStatic,int([3],kind=HSIZE_T),int([4],kind=HSIZE_T))
       call Assert("re-read part of a 1-D array double dataset to static array",doubleValueArray(3:6),doubleValueArrayRereadStatic)
       ! Read part of a double 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset1dArray",doubleValueArrayReread,int([3],kind=HSIZE_T),int([4],kind=HSIZE_T))
       call Assert("re-read part of a 1-D array double dataset to allocatable array",doubleValueArray(3:6),doubleValueArrayReread)
       deallocate(doubleValueArrayReread)
       block
         type(hdf5Object) :: datasetObject
         ! Open the dataset.
         datasetObject=groupObject%openDataset("doubleDataset1dArray")
         ! Create a reference to the dataset object.
         call groupObject%createReference1D(datasetObject,"doubleReference",int([3],kind=HSIZE_T),int([2],kind=HSIZE_T))
       end block
       ! Read a double 1-D array dataset from the group into a static array.
       doubleValueArrayRereadStatic=0
       call groupObject%readDatasetStatic("doubleReference",doubleValueArrayRereadStatic)
       call Assert("re-read referenced 1-D array double dataset to static array",doubleValueArray(3:4),doubleValueArrayRereadStatic(1:2))
       ! Read a double 1-D array dataset from the group into a static array.
       call groupObject%readDataset("doubleReference",doubleValueArrayReread)
       call Assert("re-read referenced 1-D array double dataset to allocatable array",doubleValueArray(3:4),doubleValueArrayReread(1:2))
       deallocate(doubleValueArrayReread)
       ! Read part of a double 1-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleReference",doubleValueArrayReread,int([2],kind=HSIZE_T),int([1],kind=HSIZE_T))
       call Assert("re-read part of referenced 1-D array double dataset to allocatable array",doubleValueArray(4:4),doubleValueArrayReread(1:1))
       deallocate(doubleValueArrayReread)

       ! Write an extensible double 2-D array dataset to the group.
       if (appendableOK) then
          doubleValueArray2d=3.141d0
          call groupObject%writeDataset(doubleValueArray2d,"doubleDataset2dArrayExtensible","This is an extensible dataset",appendTo=.true.)
          ! Append to the extensible double 2-D array dataset.
          doubleValueArray2d=1.414d0
          call groupObject%writeDataset(doubleValueArray2d,"doubleDataset2dArrayExtensible",appendTo=.true.)
          ! Write a new dataset to which we will then append.
          doubleValueArray2d=3.141d0
          call groupObject%writeDataset(doubleValueArray2d,"doubleDataset2dArrayExtensibleAppendDim","This is an extensible dataset to which we will append",appendTo=.true.)
          ! Append to the extensible double 2-D array dataset.
          doubleValueArray2d=1.414d0
          call groupObject%writeDataset(doubleValueArray2d(1:1,:),"doubleDataset2dArrayExtensibleAppendDim",appendTo=.true.,appendDimension=1)
          ! Re-read the array.
          call groupObject%readDataset("doubleDataset2dArrayExtensibleAppendDim",doubleValueArray2dReread)
          allocate(doubleValueArray2dRereadExpect(11,10))
          doubleValueArray2dRereadExpect( :,:)=3.141d0
          doubleValueArray2dRereadExpect(11,:)=1.414d0
          call Assert("append to arbitrary dimension of a 2-D array double dataset",doubleValueArray2dReread,doubleValueArray2dRereadExpect)
          deallocate(doubleValueArray2dRereadExpect)
          deallocate(doubleValueArray2dReread)
       end if

       ! Write a double 2-D array dataset to the group.
       doubleValueArray2d=3.141d0
       call groupObject%writeDataset(doubleValueArray2d,"doubleDataset2dArray","This is an example dataset")
       ! Write a double 2-D array dataset to the group.
       doubleValueArray2d=reshape([                                                                                               &
            & [0.237788d0 ,0.291066d0,0.845814d0,0.152208d0,0.585537d0,0.193475d0 ,0.810623d0,0.173531d0,0.484983d0 ,0.151863d0], &
            & [0.410484d0 ,0.974563d0,0.785438d0,0.133273d0,0.431549d0,0.348772d0 ,0.906572d0,0.695991d0,0.436928d0 ,0.174665d0], &
            & [0.791948d0 ,0.811843d0,0.190154d0,0.572683d0,0.163506d0,0.160362d0 ,0.735122d0,0.642781d0,0.364369d0 ,0.667146d0], &
            & [0.512184d0 ,0.893128d0,0.541979d0,0.66375d0 ,0.434586d0,0.18809d0  ,0.483802d0,0.307417d0,0.55204d0  ,0.285505d0], &
            & [0.70769d0  ,0.933332d0,0.26389d0 ,0.309795d0,0.619072d0,0.765444d0 ,0.357948d0,0.187444d0,0.249545d0 ,0.443694d0], &
            & [0.0380047d0,0.885859d0,0.863514d0,0.624677d0,0.172947d0,0.700814d0 ,0.187242d0,0.135365d0,0.433696d0 ,0.762279d0], &
            & [0.640738d0 ,0.54239d0 ,0.973174d0,0.04795d0 ,0.85148d0 ,0.0798072d0,0.844763d0,0.622515d0,0.022242d0 ,0.142351d0], &
            & [0.128965d0 ,0.450326d0,0.326675d0,0.662633d0,0.257568d0,0.19853d0  ,0.753418d0,0.41144d0 ,0.490829d0 ,0.702937d0], &
            & [0.40262d0  ,0.508495d0,0.765204d0,0.758824d0,0.55126d0 ,0.373559d0 ,0.223917d0,0.616797d0,0.0957982d0,0.245396d0], &
            & [0.424457d0 ,0.174276d0,0.735554d0,0.312812d0,0.563389d0,0.693084d0 ,0.995854d0,0.549814d0,0.772399d0 ,0.364521d0]  &
            &                     ],                                                                                              &
            & [10,10]                                                                                                             &
            &                    )
       call groupObject%writeDataset(doubleValueArray2d,"doubleDataset2dArray")
       ! Read a double 2-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset2dArray",doubleValueArray2dRereadStatic)
       call Assert("re-read 2-D array double dataset to static array",doubleValueArray2d,doubleValueArray2dRereadStatic)
       ! Read a double 2-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset2dArray",doubleValueArray2dReread)
       call Assert("re-read 2-D array double dataset to allocatable array",doubleValueArray2d,doubleValueArray2dReread)
       deallocate(doubleValueArray2dReread)
       ! Read part of a double 2-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset2dArray",doubleValueArray2dRereadStatic,int([3,6],kind=HSIZE_T),int([4,3],kind=HSIZE_T))
       call Assert("re-read part of a 2-D array double dataset to static array",doubleValueArray2d(3:6,6:8),doubleValueArray2dRereadStatic(1:4,1:3))
       ! Read part of a double 2-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset2dArray",doubleValueArray2dReread,int([3,6],kind=HSIZE_T),int([4,3],kind=HSIZE_T))
       call Assert("re-read part of a 2-D array double dataset to allocatable array",doubleValueArray2d(3:6,6:8),doubleValueArray2dReread(1:4,1:3))
       deallocate(doubleValueArray2dReread)
       block
         type(hdf5Object) :: datasetObject
         ! Open the dataset.
         datasetObject=groupObject%openDataset("doubleDataset2dArray")
         ! Create a reference to the dataset object.
         call groupObject%createReference2D(datasetObject,"double2dReference",int([3,5],kind=HSIZE_T),int([2,3],kind=HSIZE_T))
       end block
       ! Read a double 2-D array dataset from the group into a static array.
       doubleValueArray2dRereadStatic=0
       call groupObject%readDatasetStatic("double2dReference",doubleValueArray2dRereadStatic)
       call Assert("re-read referenced 2-D array double dataset to static array",doubleValueArray2d(3:4,5:7),doubleValueArray2dRereadStatic(1:2,1:3))
       ! Read part of a double 2-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("double2dReference",doubleValueArray2dReread,int([2,2],kind=HSIZE_T),int([1,2],kind=HSIZE_T))
       call Assert("re-read part of referenced 2-D array double dataset to allocatable array",doubleValueArray2d(4:4,6:7),doubleValueArray2dReread(1:1,1:2))
       deallocate(doubleValueArray2dReread)

       ! Read a double 2-D array dataset from the group into a static array.
       call groupObject%readDataset("double2dReference",doubleValueArray2dReread)
       call Assert("re-read referenced 2-D array double dataset to allocatable array",doubleValueArray2d(3:4,5:7),doubleValueArray2dReread(1:2,1:3))
       deallocate(doubleValueArray2dReread)

       ! Write an extensible double 3-D array dataset to the group.
       if (appendableOK) then
          doubleValueArray3d=3.141d0
          call groupObject%writeDataset(doubleValueArray3d,"doubleDataset3dArrayExtensible","This is an extensible dataset",appendTo=.true.)
          ! Append to the extensible double 3-D array dataset.
          doubleValueArray3d=1.414d0
          call groupObject%writeDataset(doubleValueArray3d,"doubleDataset3dArrayExtensible",appendTo=.true.)
          ! Write a new dataset to which we will then append.
          doubleValueArray3d=3.141d0
          call groupObject%writeDataset(doubleValueArray3d,"doubleDataset3dArrayExtensibleAppendDim","This is an extensible dataset to which we will append",appendTo=.true.)
          ! Append to the extensible double 2-D array dataset.
          doubleValueArray3d=1.414d0
          call groupObject%writeDataset(doubleValueArray3d(:,1:1,:),"doubleDataset3dArrayExtensibleAppendDim",appendTo=.true.,appendDimension=2)
          ! Re-read the array.
          call groupObject%readDataset("doubleDataset3dArrayExtensibleAppendDim",doubleValueArray3dReread)
          allocate(doubleValueArray3dRereadExpect(10,11,10))
          doubleValueArray3dRereadExpect(:, :,:)=3.141d0
          doubleValueArray3dRereadExpect(:,11,:)=1.414d0
          call Assert("append to arbitrary dimension of a 3-D array double dataset",doubleValueArray3dReread,doubleValueArray3dRereadExpect)
          deallocate(doubleValueArray3dRereadExpect)
          deallocate(doubleValueArray3dReread)
       end if

       ! Write a double 3-D array dataset to the group.
       doubleValueArray3d=3.141d0
       call groupObject%writeDataset(doubleValueArray3d,"doubleDataset3dArray","This is an example dataset")
       ! Write a double 3-D array dataset to the group.
       doubleValueArray3d=1.414d0
       call groupObject%writeDataset(doubleValueArray3d,"doubleDataset3dArray")
       ! Read a double 3-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset3dArray",doubleValueArray3dRereadStatic)
       call Assert("re-read 3-D array double dataset to static array",doubleValueArray3d,doubleValueArray3dRereadStatic)
       ! Read a double 3-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset3dArray",doubleValueArray3dReread)
       call Assert("re-read 3-D array double dataset to allocatable array",doubleValueArray3d,doubleValueArray3dReread)
       deallocate(doubleValueArray3dReread)
       ! Read part of a double 3-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset3dArray",doubleValueArray3dRereadStatic,int([3,6,2],kind=HSIZE_T),int([4,3,5],kind=HSIZE_T))
       call Assert("re-read part of a 3-D array double dataset to static array",doubleValueArray3d(3:6,6:8,2:6),doubleValueArray3dRereadStatic(1:4,1:3,1:5))
       ! Read part of a double 3-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset3dArray",doubleValueArray3dReread,int([3,6,2],kind=HSIZE_T),int([4,3,5],kind=HSIZE_T))
       call Assert("re-read part of a 3-D array double dataset to allocatable array",doubleValueArray3d(3:6,6:8,2:6),doubleValueArray3dReread(1:4,1:3,1:5))
       deallocate(doubleValueArray3dReread)
       block
         type(hdf5Object) :: datasetObject
         ! Open the dataset.
         datasetObject=groupObject%openDataset("doubleDataset3dArray")
         ! Create a reference to the dataset object.
         call groupObject%createReference3D(datasetObject,"double3dReference",int([3,5,2],kind=HSIZE_T),int([2,3,4],kind=HSIZE_T))
       end block
       ! Read a double 3-D array dataset from the group into a static array.
       doubleValueArray3dRereadStatic=0
       call groupObject%readDatasetStatic("double3dReference",doubleValueArray3dRereadStatic)
       call Assert("re-read referenced 3-D array double dataset to static array",doubleValueArray3d(3:4,5:7,2:5),doubleValueArray3dRereadStatic(1:2,1:3,1:4))
       ! Read a double 3-D array dataset from the group into a static array.
       call groupObject%readDataset("double3dReference",doubleValueArray3dReread)
       call Assert("re-read referenced 3-D array double dataset to allocatable array",doubleValueArray3d(3:4,5:7,2:5),doubleValueArray3dReread(1:2,1:3,1:4))
       deallocate(doubleValueArray3dReread)
       ! Read part of a double 3-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("double3dReference",doubleValueArray3dReread,int([2,2,1],kind=HSIZE_T),int([1,2,3],kind=HSIZE_T))
       call Assert("re-read part of referenced 3-D array double dataset to allocatable array",doubleValueArray3d(4:4,6:7,2:4),doubleValueArray3dReread(1:1,1:2,1:3))
       deallocate(doubleValueArray3dReread)

       ! Write an extensible double 4-D array dataset to the group.
       if (appendableOK) then
          doubleValueArray4d=3.141d0
          call groupObject%writeDataset(doubleValueArray4d,"doubleDataset4dArrayExtensible","This is an extensible dataset",appendTo=.true.)
          ! Append to the extensible double 4-D array dataset.
          doubleValueArray4d=1.414d0
          call groupObject%writeDataset(doubleValueArray4d,"doubleDataset4dArrayExtensible",appendTo=.true.)
          ! Write a new dataset to which we will then append.
          doubleValueArray4d=3.141d0
          call groupObject%writeDataset(doubleValueArray4d,"doubleDataset4dArrayExtensibleAppendDim","This is an extensible dataset to which we will append",appendTo=.true.)
          ! Append to the extensible double 2-D array dataset.
          doubleValueArray4d=1.414d0
          call groupObject%writeDataset(doubleValueArray4d(:,:,1:1,:),"doubleDataset4dArrayExtensibleAppendDim",appendTo=.true.,appendDimension=3)
          ! Re-read the array.
          call groupObject%readDataset("doubleDataset4dArrayExtensibleAppendDim",doubleValueArray4dReread)
          allocate(doubleValueArray4dRereadExpect(10,10,11,10))
          doubleValueArray4dRereadExpect(:,:, :,:)=3.141d0
          doubleValueArray4dRereadExpect(:,:,11,:)=1.414d0
          call Assert("append to arbitrary dimension of a 4-D array double dataset",doubleValueArray4dReread,doubleValueArray4dRereadExpect)
          deallocate(doubleValueArray4dRereadExpect)
          deallocate(doubleValueArray4dReread)
       end if

       ! Write a double 4-D array dataset to the group.
       doubleValueArray4d=3.141d0
       call groupObject%writeDataset(doubleValueArray4d,"doubleDataset4dArray","This is an example dataset")
       ! Write a double 4-D array dataset to the group.
       doubleValueArray4d=1.414d0
       call groupObject%writeDataset(doubleValueArray4d,"doubleDataset4dArray")
       ! Read a double 4-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset4dArray",doubleValueArray4dRereadStatic)
       call Assert("re-read 4-D array double dataset to static array",doubleValueArray4d,doubleValueArray4dRereadStatic)
       ! Read a double 4-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset4dArray",doubleValueArray4dReread)
       call Assert("re-read 4-D array double dataset to allocatable array",doubleValueArray4d,doubleValueArray4dReread)
       deallocate(doubleValueArray4dReread)
       ! Read part of a double 4-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset4dArray",doubleValueArray4dRereadStatic,int([3,6,2,7],kind=HSIZE_T),int([4,3,5,2],kind=HSIZE_T))
       call Assert("re-read part of a 4-D array double dataset to static array",doubleValueArray4d(3:6,6:8,2:6,7:8),doubleValueArray4dRereadStatic(1:4,1:3,1:5,1:2))
       ! Read part of a double 4-D array dataset from the group into a array.
       call groupObject%readDataset("doubleDataset4dArray",doubleValueArray4dReread,int([3,6,2,7],kind=HSIZE_T),int([4,3,5,2],kind=HSIZE_T))
       call Assert("re-read part of a 4-D array double dataset to allocatable array",doubleValueArray4d(3:6,6:8,2:6,7:8),doubleValueArray4dReread(1:4,1:3,1:5,1:2))
       deallocate(doubleValueArray4dReread)
       block
         type(hdf5Object) :: datasetObject
         ! Check the dimensions of the dataset.
         datasetObject=groupObject%openDataset("doubleDataset4dArray")
         call Assert("get dimensions of a dataset",[10,10,10,10],[int(datasetObject%size(1)),int(datasetObject%size(2)),int(datasetObject%size(3)),int(datasetObject%size(4))])
       end block
       block
         type(hdf5Object) :: datasetObject
         ! Open the dataset.
         datasetObject=groupObject%openDataset("doubleDataset4dArray")
         ! Create a reference to the dataset object.
         call groupObject%createReference4D(datasetObject,"double4dReference",int([3,5,2,6],kind=HSIZE_T),int([2,3,4,3],kind=HSIZE_T))
       end block
       ! Read a double 4-D array dataset from the group into a static array.
       doubleValueArray4dRereadStatic=0
       call groupObject%readDatasetStatic("double4dReference",doubleValueArray4dRereadStatic)
       call Assert("re-read referenced 4-D array double dataset to static array",doubleValueArray4d(3:4,5:7,2:5,6:8),doubleValueArray4dRereadStatic(1:2,1:3,1:4,1:3))
       ! Read a double 4-D array dataset from the group into a static array.
       call groupObject%readDataset("double4dReference",doubleValueArray4dReread)
       call Assert("re-read referenced 4-D array double dataset to allocatable array",doubleValueArray4d(3:4,5:7,2:5,6:8),doubleValueArray4dReread(1:2,1:3,1:4,1:3))
       deallocate(doubleValueArray4dReread)
       ! Read part of a double 4-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("double4dReference",doubleValueArray4dReread,int([2,2,1,1],kind=HSIZE_T),int([1,2,3,2],kind=HSIZE_T))
       call Assert("re-read part of referenced 4-D array double dataset to allocatable array",doubleValueArray4d(4:4,6:7,2:4,6:7),doubleValueArray4dReread(1:1,1:2,1:3,1:2))
       deallocate(doubleValueArray4dReread)

       ! Write an extensible double 5-D array dataset to the group.
       if (appendableOK) then
          doubleValueArray5d=3.141d0
          call groupObject%writeDataset(doubleValueArray5d,"doubleDataset5dArrayExtensible","This is an extensible dataset",appendTo=.true.)
          ! Append to the extensible double 5-D array dataset.
          doubleValueArray5d=1.415d0
          call groupObject%writeDataset(doubleValueArray5d,"doubleDataset5dArrayExtensible",appendTo=.true.)
          ! Write a new dataset to which we will then append.
          doubleValueArray5d=3.141d0
          call groupObject%writeDataset(doubleValueArray5d,"doubleDataset5dArrayExtensibleAppendDim","This is an extensible dataset to which we will append",appendTo=.true.)
          ! Append to the extensible double 2-D array dataset.
          doubleValueArray5d=1.415d0
          call groupObject%writeDataset(doubleValueArray5d(:,:,:,1:1,:),"doubleDataset5dArrayExtensibleAppendDim",appendTo=.true.,appendDimension=4)
          ! Re-read the array.
          call groupObject%readDataset("doubleDataset5dArrayExtensibleAppendDim",doubleValueArray5dReread)
          allocate(doubleValueArray5dRereadExpect(10,10,10,11,10))
          doubleValueArray5dRereadExpect(:,:,:, :,:)=3.141d0
          doubleValueArray5dRereadExpect(:,:,:,11,:)=1.415d0
          call Assert("append to arbitrary dimension of a 5-D array double dataset",doubleValueArray5dReread,doubleValueArray5dRereadExpect)
          deallocate(doubleValueArray5dRereadExpect)
          deallocate(doubleValueArray5dReread)
       end if

       ! Write a double 5-D array dataset to the group.
       doubleValueArray5d=3.141d0
       call groupObject%writeDataset(doubleValueArray5d,"doubleDataset5dArray","This is an example dataset")
       ! Write a double 5-D array dataset to the group.
       doubleValueArray5d=1.415d0
       call groupObject%writeDataset(doubleValueArray5d,"doubleDataset5dArray")
       ! Read a double 5-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset5dArray",doubleValueArray5dRereadStatic)
       call Assert("re-read 5-D array double dataset to static array",doubleValueArray5d,doubleValueArray5dRereadStatic)
       ! Read a double 5-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset5dArray",doubleValueArray5dReread)
       call Assert("re-read 5-D array double dataset to allocatable array",doubleValueArray5d,doubleValueArray5dReread)
       deallocate(doubleValueArray5dReread)
       ! Read part of a double 5-D array dataset from the group into a static array.
       call groupObject%readDatasetStatic("doubleDataset5dArray",doubleValueArray5dRereadStatic,int([3,6,2,7,2],kind=HSIZE_T),int([4,3,5,2,6],kind=HSIZE_T))
       call Assert("re-read part of a 5-D array double dataset to static array",doubleValueArray5d(3:6,6:8,2:6,7:8,2:7),doubleValueArray5dRereadStatic(1:4,1:3,1:5,1:2,1:6))
       ! Read part of a double 5-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("doubleDataset5dArray",doubleValueArray5dReread,int([3,6,2,7,2],kind=HSIZE_T),int([4,3,5,2,6],kind=HSIZE_T))
       call Assert("re-read part of a 5-D array double dataset to allocatable array",doubleValueArray5d(3:6,6:8,2:6,7:8,2:7),doubleValueArray5dReread(1:4,1:3,1:5,1:2,1:6))
       deallocate(doubleValueArray5dReread)
       block
         type(hdf5Object) :: datasetObject
         ! Open the dataset.
         datasetObject=groupObject%openDataset("doubleDataset5dArray")
         ! Create a reference to the dataset object.
         call groupObject%createReference5D(datasetObject,"double5dReference",int([3,5,2,6,2],kind=HSIZE_T),int([2,3,4,3,8],kind=HSIZE_T))
       end block
       ! Read a double 5-D array dataset from the group into a static array.
       doubleValueArray5dRereadStatic=0
       call groupObject%readDatasetStatic("double5dReference",doubleValueArray5dRereadStatic)
       call Assert("re-read referenced 5-D array double dataset to static array",doubleValueArray5d(3:4,5:7,2:5,6:8,2:9),doubleValueArray5dRereadStatic(1:2,1:3,1:4,1:3,1:8))
       ! Read a double 5-D array dataset from the group into a static array.
       call groupObject%readDataset("double5dReference",doubleValueArray5dReread)
       call Assert("re-read referenced 5-D array double dataset to allocatable array",doubleValueArray5d(3:4,5:7,2:5,6:8,2:9),doubleValueArray5dReread(1:2,1:3,1:4,1:3,1:8))
       deallocate(doubleValueArray5dReread)
       ! Read part of a double 5-D array dataset from the group into an allocatable array.
       call groupObject%readDataset("double5dReference",doubleValueArray5dReread,int([2,2,1,1,4],kind=HSIZE_T),int([1,2,3,2,4],kind=HSIZE_T))
       call Assert("re-read part of referenced 5-D array double dataset to allocatable array",doubleValueArray5d(4:4,6:7,2:4,6:7,5:8),doubleValueArray5dReread(1:1,1:2,1:3,1:2,1:4))
       deallocate(doubleValueArray5dReread)

       ! Write a character 1-D array dataset to the group.
       characterValueArray='aAbBcCdDeEfFgGhH'
       call groupObject%writeDataset(characterValueArray,"characterDataset1dArray")
       ! Read the character 1-D array dataset back.
       call groupObject%readDataset("characterDataset1dArray",characterValueArrayReread)
       call Assert("re-read 1-D array character dataset",characterValueArray,characterValueArrayReread)
       deallocate(characterValueArrayReread)
       ! Read the character 1-D array dataset back to a static array.
       call groupObject%readDatasetStatic("characterDataset1dArray",characterValueArrayRereadStatic)
       call Assert("re-read 1-D array character dataset to static array",characterValueArray,characterValueArrayRereadStatic)

       ! Write a varying string 1-D array dataset to the group.
       varStringValueArray='qazwsxedcrfvtgbyhnujmikolp'
       call groupObject%writeDataset(varStringValueArray,"varStringDataset1dArray")
       ! Read the varying string 1-D array dataset back.
       call groupObject%readDataset("varStringDataset1dArray",varStringValueArrayReread)
       call Assert("re-read 1-D array varString dataset",varStringValueArray,varStringValueArrayReread)
       deallocate(varStringValueArrayReread)
       ! Read the varying string 1-D array dataset back to a static array.
       call groupObject%readDatasetStatic("varStringDataset1dArray",varStringValueArrayRereadStatic)
       call Assert("re-read 1-D array varString dataset to static array",varStringValueArray,varStringValueArrayRereadStatic)

       ! Write a variable length 2D double dataset to the group.
       allocate(varDoubleArray2D(10))
       do i=1,10
          allocate(varDoubleArray2D(i)%row(11-i))
          varDoubleArray2D(i)%row=float(2*i)
       end do
       call groupObject%writeDataset(varDoubleArray2D,"varDoubleDataset2dArray")
       ! Read the variable-length double array dataset back.
       call groupObject%readDataset("varDoubleDataset2dArray",varDoubleDataset2dArrayReread)
       call Assert(                                                                         &
            &      "re-read 2-D array varDouble dataset"                                  , &
            &      [                                                                        &
            &       all(varDoubleArray2D( 1)%row == varDoubleDataset2dArrayReread( 1)%row), &
            &       all(varDoubleArray2D( 2)%row == varDoubleDataset2dArrayReread( 2)%row), &
            &       all(varDoubleArray2D( 3)%row == varDoubleDataset2dArrayReread( 3)%row), &
            &       all(varDoubleArray2D( 4)%row == varDoubleDataset2dArrayReread( 4)%row), &
            &       all(varDoubleArray2D( 5)%row == varDoubleDataset2dArrayReread( 5)%row), &
            &       all(varDoubleArray2D( 6)%row == varDoubleDataset2dArrayReread( 6)%row), &
            &       all(varDoubleArray2D( 7)%row == varDoubleDataset2dArrayReread( 7)%row), &
            &       all(varDoubleArray2D( 8)%row == varDoubleDataset2dArrayReread( 8)%row), &
            &       all(varDoubleArray2D( 9)%row == varDoubleDataset2dArrayReread( 9)%row), &
            &       all(varDoubleArray2D(10)%row == varDoubleDataset2dArrayReread(10)%row)  &
            &      ]                                                                      , &
            &      spread(.true.,1,10)                                                      &
            &     )
       deallocate(varDoubleDataset2dArrayReread)
       deallocate(varDoubleArray2D             )
     
       ! Write a variable length 3D double dataset to the group.
       allocate(varDoubleArray3D(10))
       do i=1,10
          allocate(varDoubleArray3D(i)%row(11-i,i))
          do j=1,11-i
             do k=1,i
                varDoubleArray3D(i)%row(j,k)=float(2*i-3*j+4*k)
             end do
          end do
       end do
       call groupObject%writeDataset(varDoubleArray3D,"varDoubleDataset3dArray")
       ! Read the variable-length 3D double array dataset back.
       call groupObject%readDataset("varDoubleDataset3dArray",varDoubleDataset3dArrayReread)
       call Assert(                                                                         &
            &      "re-read 3-D array varDouble dataset"                                  , &
            &      [                                                                        &
            &       all(varDoubleArray3D( 1)%row == varDoubleDataset3dArrayReread( 1)%row), &
            &       all(varDoubleArray3D( 2)%row == varDoubleDataset3dArrayReread( 2)%row), &
            &       all(varDoubleArray3D( 3)%row == varDoubleDataset3dArrayReread( 3)%row), &
            &       all(varDoubleArray3D( 4)%row == varDoubleDataset3dArrayReread( 4)%row), &
            &       all(varDoubleArray3D( 5)%row == varDoubleDataset3dArrayReread( 5)%row), &
            &       all(varDoubleArray3D( 6)%row == varDoubleDataset3dArrayReread( 6)%row), &
            &       all(varDoubleArray3D( 7)%row == varDoubleDataset3dArrayReread( 7)%row), &
            &       all(varDoubleArray3D( 8)%row == varDoubleDataset3dArrayReread( 8)%row), &
            &       all(varDoubleArray3D( 9)%row == varDoubleDataset3dArrayReread( 9)%row), &
            &       all(varDoubleArray3D(10)%row == varDoubleDataset3dArrayReread(10)%row)  &
            &      ]                                                                      , &
            &      spread(.true.,1,10)                                                      &
            &     )
       deallocate(varDoubleDataset3DArrayReread)
       deallocate(varDoubleArray3D             )

       ! Write a variable length integer8 dataset to the group.
       allocate(varInteger8Array2D(10))
       do i=1,10
          allocate(varInteger8Array2D(i)%row(11-i))
          varInteger8Array2D(i)%row=2*i
       end do
       call groupObject%writeDataset(varInteger8Array2D,"varInteger8Dataset2dArray")
       ! Read the variable-length integer-8 array dataset back.
       call groupObject%readDataset("varInteger8Dataset2dArray",varInteger8Dataset2dArrayReread)
       call Assert(                                                                             &
            &      "re-read 2-D array varInteger8 dataset"                                    , &
            &      [                                                                            &
            &       all(varInteger8Array2D( 1)%row == varInteger8Dataset2dArrayReread( 1)%row), &
            &       all(varInteger8Array2D( 2)%row == varInteger8Dataset2dArrayReread( 2)%row), &
            &       all(varInteger8Array2D( 3)%row == varInteger8Dataset2dArrayReread( 3)%row), &
            &       all(varInteger8Array2D( 4)%row == varInteger8Dataset2dArrayReread( 4)%row), &
            &       all(varInteger8Array2D( 5)%row == varInteger8Dataset2dArrayReread( 5)%row), &
            &       all(varInteger8Array2D( 6)%row == varInteger8Dataset2dArrayReread( 6)%row), &
            &       all(varInteger8Array2D( 7)%row == varInteger8Dataset2dArrayReread( 7)%row), &
            &       all(varInteger8Array2D( 8)%row == varInteger8Dataset2dArrayReread( 8)%row), &
            &       all(varInteger8Array2D( 9)%row == varInteger8Dataset2dArrayReread( 9)%row), &
            &       all(varInteger8Array2D(10)%row == varInteger8Dataset2dArrayReread(10)%row)  &
            &      ]                                                                          , &
            &      spread(.true.,1,10)                                                          &
            &     )
       deallocate(varInteger8Dataset2dArrayReread)
       deallocate(varInteger8Array2D             )

       ! Write a scalar integer attribute to the group.
       integerValue=2020
       call groupObject%writeAttribute(integerValue,"integerShortAttribute")
       ! Read the scalar integer attribute back into a long integer.
       call groupObject%readAttribute("integerShortAttribute",integer8ValueReread)
       call Assert("read scalar integer attribute to long integer",int(integerValue,kind=kind_int8),integer8ValueReread)

       ! Write an integer 1-D array attribute to the group.
       integerValueArray=7
       call groupObject%writeAttribute(integerValueArray,"integerShortAttribute1dArray")
       ! Read the long integer 1-D array attribute back.
       call groupObject%readAttribute("integerShortAttribute1dArray",integer8ValueArrayReread)
       call Assert("read 1-D array integer attribute to long integer array",int(integerValueArray,kind=kind_int8),integer8ValueArrayReread)
       deallocate(integer8ValueArrayReread)
       ! Read the long integer 1-D array attribute back to a static array.
       call groupObject%readAttributeStatic("integerShortAttribute1dArray",integer8ValueArrayRereadStatic)
       call Assert("read 1-D array integer attribute to static long integer array",int(integerValueArray,kind=kind_int8),integer8ValueArrayRereadStatic)

       ! Write an integer 1-D array dataset.
       integerValueArray=[0,11,22,33,44,55,66,77,88,99]
       call groupObject%writeDataset(integerValueArray,"integerShortDataset1dArray","This is an example dataset")
       ! Read the dataset back into a long integer 1-D array dataset.
       call groupObject%readDataset("integerShortDataset1dArray",integer8ValueArrayReread)
       call Assert("read 1-D array integer dataset to 1-D long integer allocatable array",int(integerValueArray,kind=kind_int8),integer8ValueArrayReread)
       deallocate(integer8ValueArrayReread)
       ! Read the dataset back into a long integer 1-D array dataset.
       call groupObject%readDatasetStatic("integerShortDataset1dArray",integer8ValueArrayRereadStatic)
       call Assert("read 1-D array integer dataset to 1-D long integer allocatable array",int(integerValueArray,kind=kind_int8),integer8ValueArrayRereadStatic)

       ! Retrieve a list of datasets in the group.
       if (iPass == 2) then
          block
            type(varying_string), allocatable, dimension(:) :: datasetNames
            call groupObject%datasets(datasetNames)
            datasetNamesReference=[                              &
                 &                 "anotherReference          ", &
                 &                 "characterDataset1dArray   ", &
                 &                 "double2dReference         ", &
                 &                 "double3dReference         ", &
                 &                 "double4dReference         ", &
                 &                 "double5dReference         ", &
                 &                 "doubleDataset1dArray      ", &
                 &                 "doubleDataset2dArray      ", &
                 &                 "doubleDataset3dArray      ", &
                 &                 "doubleDataset4dArray      ", &
                 &                 "doubleDataset5dArray      ", &
                 &                 "doubleReference           ", &
                 &                 "integer8Dataset1dArray    ", &
                 &                 "integerDataset1dArray     ", &
                 &                 "integerShortDataset1dArray", &
                 &                 "myReference               ", &
                 &                 "varDoubleDataset2dArray   ", &
                 &                 "varDoubleDataset3dArray   ", &
                 &                 "varInteger8Dataset2dArray ", &
                 &                 "varStringDataset1dArray   "  &
                 &                ]
            do i=1,size(datasetNamesReference)
               datasetNamesReference(i)=trim(datasetNamesReference(i))
            end do
            call Assert("recover correct number of datasets in group",size(datasetNames),20)
            call Assert("recover names of datasets in group",datasetNames,datasetNamesReference)
          end block
       end if
       
       ! Write a very large (>4GB) dataset to test that chunking limits
       ! the chunksize to less than the maximum allowed.
       if (iPass == 2) then
          if (allocated(doubleValueArray4dReread)) deallocate(doubleValueArray4dReread)
          allocate(doubleValueArray4dReread(600,100,100,100))
          call fileObject%writeDataset(doubleValueArray4dReread,'bigDataset','A dataset larger than 4GB.',chunkSize=1024_hsize_t)
       end if
       
       ! Write an attribute of length >64KB by forcing dense storage of
       ! attributes in s group.
       block
         type(hdf5Object) :: groupObject2
         groupObject2=fileObject%openGroup("myGroup64k",comment="This is my group for 64k attributes.",objectsOverwritable=.true.,chunkSize=1024_hsize_t&
              &,compressionLevel=9,attributesCompactMaxiumum=0)
         varStringValue=repeat("rain day happy calm dog joy joy over bright falls rain lazy shin ",1025)
         call groupObject2%writeAttribute(varStringValue,"varStringAttribute64k")
       end block

       ! Read a 32-bit unsigned integer 1-D array into a 64-bit signed integer 1-D array.
       if (iPass==2) then
          block
            type(hdf5Object) :: fileObject2, groupObject2
            ! Open the HDF5 file which stores the smallest and the largest 32-bit unsigned integers.
            fileObject2 =hdf5Object("testSuite/data/IntegerRangeU32.hdf5")
            ! Open the root group.
            groupObject2=fileObject2%openGroup("/",comment="Root group.")
            ! Read the dataset.
            call groupObject2%readDataset('IntegerRangeU32',integerRangeU32)
            call Assert("read 32-bit unsigned integers into 64-bit signed integers",[0_kind_int8,4294967295_kind_int8],integerRangeU32)
            deallocate(integerRangeU32)
          end block
       end if
  end do

  ! Test identifying HDF5 file.
  call Assert("test if file is HDF5",IO_HDF5_Is_HDF5('testSuite/outputs/test.IO.HDF5.hdf5'),.true.)

  ! Test of h5py compatibility.
  block
     type(hdf5Object) :: fileObject
     call Unit_Tests_Begin_Group("h5py compatibility")
     call System_Command_Do("./testSuite/scripts/generate_h5py.py")
     fileObject=hdf5Object("testSuite/outputs/h5py.hdf5",overWrite=.false.,objectsOverwritable=.false.)
     call fileObject%readAttribute("stringAttribute"            ,          varStringValueReread                            )
     call fileObject%readAttribute("stringAttribute"            ,          characterValueReread                            )
     call Assert("read h5py string attribute (character)",characterValueReread,"this is a variable length string")
     call Assert("read h5py string attribute (varying_string)",varStringValueReread,var_str("this is a variable length string"))
     call Unit_Tests_End_Group()
   end block

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_IO_HDF5
