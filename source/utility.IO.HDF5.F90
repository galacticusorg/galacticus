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
Contains a module that implements simple and convenient interfaces to a variety of HDF5 functionality.
!!}

! Specify an explicit dependence on the hdf5_cTypes.o object file.
!: $(BUILDPATH)/hdf5_cTypes.o

module IO_HDF5
  !!{
  Implements simple and convenient interfaces to a variety of HDF5 functionality.
  !!}
  use            :: HDF5              , only : hid_t         , hsize_t, size_t
  use, intrinsic :: ISO_C_Binding     , only : c_char        , c_int  , c_ptr , c_size_t
  use            :: ISO_Varying_String, only : varying_string
  use            :: Kind_Numbers      , only : kind_int8
  use            :: Resource_Manager  , only : resourceManager
  implicit none
  private
  public :: hdf5Object          , hdf5VarDouble  , hdf5VarInteger8, ioHDF5AccessInitialize, &
       &    IO_HDF5_Set_Defaults, IO_HDF5_Is_HDF5, hdf5VarDouble2D
#ifdef DEBUGHDF5
  public :: IO_HDF5_Start_Critical, IO_HDF5_End_Critical

  logical                                                :: inCritical=.false.
  !$omp threadprivate(inCritical)
#endif

  ! Record of initialization of this module.
  logical                                                :: hdf5IsInitialized         =.false.
  integer                                                :: initializationsCount      =0

  ! Type of HDF5 objects.
  integer                            , parameter, public :: hdf5ObjectTypeNull        =0
  integer                            , parameter, public :: hdf5ObjectTypeFile        =1
  integer                            , parameter, public :: hdf5ObjectTypeGroup       =2
  integer                            , parameter, public :: hdf5ObjectTypeDataset     =3
  integer                            , parameter, public :: hdf5ObjectTypeAttribute   =4

  ! Data types.
  integer                            , parameter, public :: hdf5DataTypeNull          =0
  integer                            , parameter, public :: hdf5DataTypeInteger       =1
  integer                            , parameter, public :: hdf5DataTypeInteger8      =2
  integer                            , parameter, public :: hdf5DataTypeDouble        =3
  integer                            , parameter, public :: hdf5DataTypeCharacter     =4
  integer                            , parameter, public :: hdf5DataTypeVlenDouble    =5
  integer                            , parameter, public :: hdf5DataTypeVlenVlenDouble=6
  integer                            , parameter, public :: hdf5DataTypeVlenInteger8  =7

  ! Chunking and compression parameters.
  integer                                                :: hdf5CompressionLevel      =-1
  integer(kind=HSIZE_T)                                  :: hdf5ChunkSize             =-1

  ! Arrays of compatible datatypes.
  integer(kind=HID_T  ), dimension(5)           , public :: H5T_NATIVE_DOUBLES
  integer(kind=HID_T  ), dimension(5)           , public :: H5T_NATIVE_INTEGERS
  integer(kind=HID_T  ), dimension(2)           , public :: H5T_NATIVE_UNSIGNED_INTEGERS
  integer(kind=HID_T  ), dimension(3)           , public :: H5T_NATIVE_INTEGER_8S
  integer(kind=HID_T  ), dimension(8)           , public :: H5T_NATIVE_INTEGER_8AS
  integer(kind=HID_T  ), dimension(1)           , public :: H5T_VLEN_DOUBLE             , H5T_VLEN_VLEN_DOUBLE
  integer(kind=HID_T  ), dimension(1)           , public :: H5T_VLEN_INTEGER8
  integer(kind=HID_T  )                         , public :: H5T_INTEGER8

  type hdf5Object
     !!{
     A structure that holds properties of HDF5 objects.
     !!}
     private
     logical                           :: isOpenValue      =  .false.
     logical                           :: isOverwritable
     logical                           :: readOnly
     logical                           :: isTemporary      =  .false.
     integer(hid_t          ), pointer :: objectID         => null() , fileID             => null()
     type   (resourceManager)          :: objectManager              , fileManager
     type   (varying_string )          :: objectFile
     type   (varying_string )          :: objectLocation
     type   (varying_string )          :: objectName
     integer                           :: hdf5ObjectType
     integer(hsize_t        )          :: chunkSize
     integer                           :: compressionLevel
     logical                           :: chunkSizeSet               , compressionLevelSet
     type   (hdf5Object     ), pointer :: parentObject     => null()
   contains
     !![
     <methods>
       <method description="Write an attribute to an HDF5 object." method="writeAttribute" />
       <method description="Write a dataset to an HDF5 group." method="writeDataset" />
       <method description="Read an attribute from an HDF5 object." method="readAttribute" />
       <method description="Read an attribute from an HDF5 object into a static array." method="readAttributeStatic" />
       <method description="Read a dataset from an HDF5 group into an allocatable array." method="readDataset" />
       <method description="Read a dataset from an HDF5 group into a static array." method="readDatasetStatic" />
       <method description="Read a column from an HDF5 table into an allocatable array." method="readTable" />
       <method description="Return the size of a dataset." method="size" />
       <method description="Check if an object has a named attribute." method="hasAttribute" />
       <method description="Check if an object has a named group." method="hasGroup" />
       <method description="Check if an object has a named dataset." method="hasDataset" />
       <method description="Get a list of datasets in a group object." method="datasets" />
       <method description="Check the type and rank of an attribute." method="assertAttributeType" />
       <method description="Check the type and rank of a dataset." method="assertDatasetType" />
       <method description="Return the rank of a dataset." method="rank" />
       <method description="Return true if a dataset is a reference." method="isReference" />
       <method description="Return true if an object is open." method="isOpen" />
       <method description="Remove the named object." method="remove" />
       <method description="Return the object type." method="objectType" />
       <method description="Create a reference to a 1D dataset." method="createReference1D" />
       <method description="Create a reference to a 2D dataset." method="createReference2D" />
       <method description="Create a reference to a 2D dataset." method="createReference3D" />
       <method description="Create a reference to a 2D dataset." method="createReference4D" />
       <method description="Create a reference to a 2D dataset." method="createReference5D" />
       <method description="Flush an HDF5 file to disk." method="flush" />
       <method description="Returns the path to a given object." method="pathTo" />
       <method description="Returns the name of a given object." method="name" />
       <method description="Returns the name of the file containing a given object." method="fileName" />
       <method description="Return a report on the location of an object suitable for inclusion in an error message." method="locationReport" />
       <method description="Open an HDF5 group and return an appropriate HDF5 object." method="openGroup" />
       <method description="Open all HDF5 groups along a path and return the appropriate HDF5 objects." method="openGroupPath" />
       <method description="Open an HDF5 dataset." method="openDataset" />
       <method description="Open an HDF5 attribute." method="openAttribute" />
       <method description="Copy an HDF5 object." method="copy" />
       <method description="Return the parent object." method="parent" />
       <method description="Create a deep copy of the object with a new HDF5 object identifier." method="deepCopy" />
       <method description="Assign HDF5 objects." method="assignment(=)"/>
     </methods>
     !!]
     final     ::                                           IO_HDF5_Finalize
     procedure ::                                           IO_HDF5_Assign
     generic   :: assignment(=)                           =>IO_HDF5_Assign
     procedure :: name                                    =>IO_HDF5_Name
     procedure :: pathTo                                  =>IO_HDF5_Path_To
     procedure :: fileName                                =>IO_HDF5_File_Name
     procedure :: locationReport                          =>IO_HDF5_Location_Report
     procedure :: openGroup                               =>IO_HDF5_Open_Group
     procedure :: openGroupPath                           =>IO_HDF5_Open_Group_Path
     procedure :: openDataset                             =>IO_HDF5_Open_Dataset
     procedure :: openAttribute                           =>IO_HDF5_Open_Attribute
     procedure :: flush                                   =>IO_HDF5_Flush
     procedure :: remove                                  =>IO_HDF5_Remove
     procedure :: IO_HDF5_Write_Attribute_Integer_Scalar
     procedure :: IO_HDF5_Write_Attribute_Integer_1D
     procedure :: IO_HDF5_Write_Attribute_Integer8_Scalar
     procedure :: IO_HDF5_Write_Attribute_Integer8_1D
     procedure :: IO_HDF5_Write_Attribute_Double_Scalar
     procedure :: IO_HDF5_Write_Attribute_Double_1D
     procedure :: IO_HDF5_Write_Attribute_Double_2D
     procedure :: IO_HDF5_Write_Attribute_Character_Scalar
     procedure :: IO_HDF5_Write_Attribute_Character_1D
     procedure :: IO_HDF5_Write_Attribute_VarString_Scalar
     procedure :: IO_HDF5_Write_Attribute_VarString_1D
     procedure :: IO_HDF5_Write_Attribute_Logical_Scalar
     generic   :: writeAttribute      => IO_HDF5_Write_Attribute_Integer_Scalar  , &
          &                              IO_HDF5_Write_Attribute_Integer_1D      , &
          &                              IO_HDF5_Write_Attribute_Integer8_Scalar , &
          &                              IO_HDF5_Write_Attribute_Integer8_1D     , &
          &                              IO_HDF5_Write_Attribute_Double_Scalar   , &
          &                              IO_HDF5_Write_Attribute_Double_1D       , &
          &                              IO_HDF5_Write_Attribute_Double_2D       , &
          &                              IO_HDF5_Write_Attribute_Character_Scalar, &
          &                              IO_HDF5_Write_Attribute_Character_1D    , &
          &                              IO_HDF5_Write_Attribute_VarString_Scalar, &
          &                              IO_HDF5_Write_Attribute_VarString_1D    , &
          &                              IO_HDF5_Write_Attribute_Logical_Scalar
     procedure :: IO_HDF5_Write_Dataset_Integer_1D
     procedure :: IO_HDF5_Write_Dataset_Integer_2D
     procedure :: IO_HDF5_Write_Dataset_Integer_3D
     procedure :: IO_HDF5_Write_Dataset_Integer8_1D
     procedure :: IO_HDF5_Write_Dataset_Integer8_2D
     procedure :: IO_HDF5_Write_Dataset_Integer8_3D
     procedure :: IO_HDF5_Write_Dataset_Double_1D
     procedure :: IO_HDF5_Write_Dataset_Double_2D
     procedure :: IO_HDF5_Write_Dataset_Double_3D
     procedure :: IO_HDF5_Write_Dataset_Double_4D
     procedure :: IO_HDF5_Write_Dataset_Double_5D
     procedure :: IO_HDF5_Write_Dataset_Double_6D
     procedure :: IO_HDF5_Write_Dataset_Character_1D
     procedure :: IO_HDF5_Write_Dataset_VarString_1D
     procedure :: IO_HDF5_Write_Dataset_VarDouble_1D
     procedure :: IO_HDF5_Write_Dataset_VarVarDouble_1D
     procedure :: IO_HDF5_Write_Dataset_VarDouble_2D
     procedure :: IO_HDF5_Write_Dataset_VarInteger8_2D
     generic   :: writeDataset        => IO_HDF5_Write_Dataset_Integer_1D        , &
          &                              IO_HDF5_Write_Dataset_Integer_2D        , &
          &                              IO_HDF5_Write_Dataset_Integer_3D        , &
          &                              IO_HDF5_Write_Dataset_Integer8_1D       , &
          &                              IO_HDF5_Write_Dataset_Integer8_2D       , &
          &                              IO_HDF5_Write_Dataset_Integer8_3D       , &
          &                              IO_HDF5_Write_Dataset_Double_1D         , &
          &                              IO_HDF5_Write_Dataset_Double_2D         , &
          &                              IO_HDF5_Write_Dataset_Double_3D         , &
          &                              IO_HDF5_Write_Dataset_Double_4D         , &
          &                              IO_HDF5_Write_Dataset_Double_5D         , &
          &                              IO_HDF5_Write_Dataset_Double_6D         , &
          &                              IO_HDF5_Write_Dataset_Character_1D      , &
          &                              IO_HDF5_Write_Dataset_VarString_1D      , &
          &                              IO_HDF5_Write_Dataset_VarDouble_1D      , &
          &                              IO_HDF5_Write_Dataset_VarVarDouble_1D   , &
          &                              IO_HDF5_Write_Dataset_VarDouble_2D      , &
          &                              IO_HDF5_Write_Dataset_VarInteger8_2D
     procedure :: IO_HDF5_Read_Attribute_Integer_Scalar
     procedure :: IO_HDF5_Read_Attribute_Integer_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Attribute_Integer_1D_Array_Static
     procedure :: IO_HDF5_Read_Attribute_Integer8_Scalar
     procedure :: IO_HDF5_Read_Attribute_Integer8_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Attribute_Integer8_1D_Array_Static
     procedure :: IO_HDF5_Read_Attribute_Double_Scalar
     procedure :: IO_HDF5_Read_Attribute_Double_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Attribute_Double_1D_Array_Static
     procedure :: IO_HDF5_Read_Attribute_Character_Scalar
     procedure :: IO_HDF5_Read_Attribute_Character_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Attribute_Character_1D_Array_Static
     procedure :: IO_HDF5_Read_Attribute_VarString_Scalar
     procedure :: IO_HDF5_Read_Attribute_VarString_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Attribute_VarString_1D_Array_Static
     generic   :: readAttribute       => IO_HDF5_Read_Attribute_Integer_Scalar                , &
          &                              IO_HDF5_Read_Attribute_Integer_1D_Array_Allocatable  , &
          &                              IO_HDF5_Read_Attribute_Integer8_Scalar               , &
          &                              IO_HDF5_Read_Attribute_Integer8_1D_Array_Allocatable , &
          &                              IO_HDF5_Read_Attribute_Double_Scalar                 , &
          &                              IO_HDF5_Read_Attribute_Double_1D_Array_Allocatable   , &
          &                              IO_HDF5_Read_Attribute_Character_Scalar              , &
          &                              IO_HDF5_Read_Attribute_Character_1D_Array_Allocatable, &
          &                              IO_HDF5_Read_Attribute_VarString_Scalar              , &
          &                              IO_HDF5_Read_Attribute_VarString_1D_Array_Allocatable
     generic   :: readAttributeStatic => IO_HDF5_Read_Attribute_Integer_1D_Array_Static       , &
          &                              IO_HDF5_Read_Attribute_Integer8_1D_Array_Static      , &
          &                              IO_HDF5_Read_Attribute_Double_1D_Array_Static        , &
          &                              IO_HDF5_Read_Attribute_Character_1D_Array_Static     , &
          &                              IO_HDF5_Read_Attribute_VarString_1D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Integer_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Integer_1D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Integer_2D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Integer_2D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Integer8_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Integer8_1D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Integer8_2D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Integer8_2D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Integer8_3D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Double_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Double_1D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Double_2D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Double_2D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Double_3D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Double_3D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Double_4D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Double_4D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Double_5D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Double_5D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Double_6D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Double_6D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_Character_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_Character_1D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_VarString_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_VarString_1D_Array_Static
     procedure :: IO_HDF5_Read_Dataset_VarDouble_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_VarVarDouble_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_VarDouble_2D_Array_Allocatable
     procedure :: IO_HDF5_Read_Dataset_VarInteger8_2D_Array_Allocatable
     generic   :: readDataset         => IO_HDF5_Read_Dataset_Integer_1D_Array_Allocatable     , &
          &                              IO_HDF5_Read_Dataset_Integer_2D_Array_Allocatable     , &
          &                              IO_HDF5_Read_Dataset_Integer8_1D_Array_Allocatable    , &
          &                              IO_HDF5_Read_Dataset_Integer8_2D_Array_Allocatable    , &
          &                              IO_HDF5_Read_Dataset_Integer8_3D_Array_Allocatable    , &
          &                              IO_HDF5_Read_Dataset_Double_1D_Array_Allocatable      , &
          &                              IO_HDF5_Read_Dataset_Double_2D_Array_Allocatable      , &
          &                              IO_HDF5_Read_Dataset_Double_3D_Array_Allocatable      , &
          &                              IO_HDF5_Read_Dataset_Double_4D_Array_Allocatable      , &
          &                              IO_HDF5_Read_Dataset_Double_5D_Array_Allocatable      , &
          &                              IO_HDF5_Read_Dataset_Double_6D_Array_Allocatable      , &
          &                              IO_HDF5_Read_Dataset_Character_1D_Array_Allocatable   , &
          &                              IO_HDF5_Read_Dataset_VarString_1D_Array_Allocatable   , &
          &                              IO_HDF5_Read_Dataset_VarDouble_1D_Array_Allocatable   , &
          &                              IO_HDF5_Read_Dataset_VarVarDouble_1D_Array_Allocatable, &
          &                              IO_HDF5_Read_Dataset_VarDouble_2D_Array_Allocatable   , &
          &                              IO_HDF5_Read_Dataset_VarInteger8_2D_Array_Allocatable
     generic   :: readDatasetStatic   => IO_HDF5_Read_Dataset_Integer_1D_Array_Static          , &
          &                              IO_HDF5_Read_Dataset_Integer_2D_Array_Static          , &
          &                              IO_HDF5_Read_Dataset_Integer8_1D_Array_Static         , &
          &                              IO_HDF5_Read_Dataset_Integer8_2D_Array_Static         , &
          &                              IO_HDF5_Read_Dataset_Double_1D_Array_Static           , &
          &                              IO_HDF5_Read_Dataset_Double_2D_Array_Static           , &
          &                              IO_HDF5_Read_Dataset_Double_3D_Array_Static           , &
          &                              IO_HDF5_Read_Dataset_Double_4D_Array_Static           , &
          &                              IO_HDF5_Read_Dataset_Double_5D_Array_Static           , &
          &                              IO_HDF5_Read_Dataset_Double_6D_Array_Static           , &
          &                              IO_HDF5_Read_Dataset_Character_1D_Array_Static        , &
          &                              IO_HDF5_Read_Dataset_VarString_1D_Array_Static
     procedure :: IO_HDF5_Read_Table_Real_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Table_Integer_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Table_Integer8_1D_Array_Allocatable
     procedure :: IO_HDF5_Read_Table_Character_1D_Array_Allocatable
     generic   :: readTable          =>IO_HDF5_Read_Table_Real_1D_Array_Allocatable    , &
          &                            IO_HDF5_Read_Table_Integer_1D_Array_Allocatable , &
          &                            IO_HDF5_Read_Table_Integer8_1D_Array_Allocatable, &
          &                            IO_HDF5_Read_Table_Character_1D_Array_Allocatable
     procedure :: size               =>IO_HDF5_Dataset_Size
     procedure :: rank               =>IO_HDF5_Dataset_Rank
     procedure :: hasAttribute       =>IO_HDF5_Has_Attribute
     procedure :: hasGroup           =>IO_HDF5_Has_Group
     procedure :: hasDataset         =>IO_HDF5_Has_Dataset
     procedure :: datasets           =>IO_HDF5_Datasets
     procedure :: assertAttributeType=>IO_HDF5_Assert_Attribute_Type
     procedure :: assertDatasetType  =>IO_HDF5_Assert_Dataset_Type
     procedure :: isReference        =>IO_HDF5_Is_Reference
     procedure :: isOpen             =>IO_HDF5_Is_Open
     procedure :: objectType         =>IO_HDF5_Object_Type
     procedure :: createReference1D  =>IO_HDF5_Create_Reference_Scalar_To_1D
     procedure :: createReference2D  =>IO_HDF5_Create_Reference_Scalar_To_2D
     procedure :: createReference3D  =>IO_HDF5_Create_Reference_Scalar_To_3D
     procedure :: createReference4D  =>IO_HDF5_Create_Reference_Scalar_To_4D
     procedure :: createReference5D  =>IO_HDF5_Create_Reference_Scalar_To_5D
     procedure :: copy               =>IO_HDF5_Copy
     procedure :: parent             =>IO_HDF5_Parent
     procedure :: deepCopy           =>IO_HDF5_Deep_Copy
  end type hdf5Object

  interface hdf5Object
     module procedure hdf5FileOpenVarStr
     module procedure hdf5FileOpenChar
  end interface hdf5Object
    
  type :: hdf5VarDouble
     !!{
     Type used for internal storage of variable-length double datasets.
     !!}
     double precision, dimension(:), pointer :: row => null()
   contains
     final :: hdf5VarDoubleDestructor
  end type hdf5VarDouble

  type :: hdf5VarDouble2D
     !!{
     Type used for internal storage of variable-length 2D double datasets.
     !!}
     double precision, dimension(:,:), pointer :: row => null()
   contains
     final :: hdf5VarDouble2DDestructor
  end type hdf5VarDouble2D

  type :: hdf5VarInteger8
     !!{
     Type used for internal storage of variable-length integer-8 datasets.
     !!}
     integer(kind_int8), dimension(:), pointer :: row => null()
   contains
     final :: hdf5VarInteger8Destructor
  end type hdf5VarInteger8
  
  type, bind(C) :: hdf5VlenC
     !!{
     Type used for C-compatible internal storage of variable-length datasets.
     !!}
     integer(c_size_t) :: length
     type   (c_ptr   ) :: p
  end type hdf5VlenC
  
  type :: hdf5VlenVlenC
     !!{
     Type used for C-compatible internal storage of fractal variable-length datasets.
     !!}
     type(hdf5VlenC), allocatable, dimension(:) :: row
  end type hdf5VlenVlenC
  
  ! Interfaces to functions in the HDF5 C API that are required due to the limited datatypes supported by the Fortran API. For
  ! example, h5dread_f() can not read a scalar dataset region reference, so we are forced to go through the h5dread() C function
  ! for this purpose.
  interface
     function H5T_C_S1_Get() bind(c,name='H5T_C_S1_Get')
       !!{
       Template for a C function that returns the {\normalfont \ttfamily H5T\_C\_S1} datatype ID.
       !!}
       import
       integer(kind=hid_t) :: H5T_C_S1_Get
     end function H5T_C_S1_Get
     function H5T_Variable_Get() bind(c,name='H5T_Variable_Get')
       !!{
       Template for a C function that returns the {\normalfont \ttfamily H5T\_C\_S1} datatype ID.
       !!}
       import
       integer(kind=size_t) :: H5T_Variable_Get
     end function H5T_Variable_Get
     function H5Awrite(attr_id,mem_type_id,buf) bind(c, name='H5Awrite')
       !!{
       Template for the HDF5 C API attribute write function.
       !!}
       import
       integer(kind=herr_t)        :: H5Awrite
       integer(kind=hid_t ), value :: attr_id , mem_type_id
       type   (c_ptr      ), value :: buf
     end function H5Awrite
     function H5Aread(attr_id,mem_type_id,buf) bind(c, name='H5Aread')
       !!{
       Template for the HDF5 C API attribute read function.
       !!}
       import
       integer(kind=herr_t)        :: H5Aread
       integer(kind=hid_t ), value :: attr_id, mem_type_id
       type   (c_ptr      ), value :: buf
     end function H5Aread
     function H5Dwrite(dataset_id,mem_type_id,mem_space_id,file_space_id,xfer_plist_id,buf) bind(c, name='H5Dwrite')
       !!{
       Template for the HDF5 C API dataset write function.
       !!}
       import
       integer(kind=herr_t)        :: H5Dwrite
       integer(kind=hid_t ), value :: dataset_id   , file_space_id, mem_space_id, mem_type_id, &
            &                         xfer_plist_id
       type   (c_ptr      ), value :: buf
     end function H5Dwrite
     function H5Dread(dataset_id,mem_type_id,mem_space_id,file_space_id,xfer_plist_id,buf) bind(c, name='H5Dread')
       !!{
       Template for the HDF5 C API dataset read function.
       !!}
       import
       integer(kind=herr_t)        :: H5Dread
       integer(kind=hid_t ), value :: dataset_id   , file_space_id, mem_space_id, mem_type_id, &
            &                         xfer_plist_id
       type   (c_ptr      ), value :: buf
     end function H5Dread
     function H5TBread_fields_name(loc_id,table_name,field_names,start,nrecords,type_size,field_offset,dst_sizes,data) bind(c, name='H5TBread_fields_name')
       !!{
       Template for the HDF5 C API table read fields by name function.
       !!}
       import
       integer  (kind=herr_t )               :: H5TBread_fields_name
       integer  (kind=hid_t  ), value        :: loc_id
       character(kind=c_char )               :: table_name          , field_names(*)
       integer  (kind=hsize_t), value        :: start               , nrecords      , type_size
       integer  (kind=size_t ), dimension(*) :: field_offset        , dst_sizes
       type     (c_ptr       ), value        :: data
     end function H5TBread_fields_name
  end interface

contains

  !! Initialization routines.

  subroutine ioHDF5AccessInitialize()
    !!{
    Initialize the HDF5 access lock.
    !!}
    use :: HDF5_Access, only : hdf5Access, hdf5AccessInitialized
    use :: Locks      , only : ompLock
    implicit none

    if (hdf5AccessInitialized) return
    hdf5Access           =ompLock()
    hdf5AccessInitialized=.true.
    return
  end subroutine ioHDF5AccessInitialize
  
  subroutine IO_HDF5_Initialize
    !!{
    Initialize the HDF5 subsystem.
    !!}
    use :: Error, only : Error_Report
    use :: HDF5 , only : H5T_IEEE_F32BE   , H5T_IEEE_F32LE    , H5T_IEEE_F64BE, H5T_IEEE_F64LE  , &
          &              H5T_NATIVE_DOUBLE, H5T_NATIVE_INTEGER, H5T_STD_I32BE , H5T_STD_U32LE   , &
          &              H5T_STD_I32LE    , H5T_STD_I64BE     , H5T_STD_I64LE , H5T_STD_U32BE   , &
          &              h5tcopy_f        , h5tset_size_f     , h5open_f      , h5tvlen_create_f, &
          &              h5tequal_f
    implicit none
    integer :: errorCode
    logical :: isLittleEndian, isBigEndian

#ifdef DEBUGHDF5
    call IO_HDF5_Assert_In_Critical()
#endif

    if (.not.hdf5IsInitialized) then
       call h5open_f(errorCode)
       if (errorCode < 0) call Error_Report('failed to initialize HDF5 subsystem'//{introspection:location})

       ! Create required datatypes.
       call h5tequal_f(H5T_NATIVE_INTEGER,H5T_STD_I32LE,isLittleEndian,errorCode)
       if (errorCode < 0) call Error_Report('failed to test endianness'//{introspection:location})
       call h5tequal_f(H5T_NATIVE_INTEGER,H5T_STD_I32LE,isBigEndian   ,errorCode)
       if (errorCode < 0) call Error_Report('failed to test endianness'//{introspection:location})
       if (isLittleEndian) then
          call h5tcopy_f(H5T_STD_I64LE,H5T_INTEGER8,errorCode)
       else if (isBigEndian) then
          call h5tcopy_f(H5T_STD_I64BE,H5T_INTEGER8,errorCode)
       else
          call Error_Report('unable to determine native endianness'//{introspection:location})
       end if
       if (errorCode < 0) call Error_Report('failed to copy integer datatype'//{introspection:location})

       ! Ensure native datatype arrays are initialized.
       H5T_NATIVE_DOUBLES          =[H5T_NATIVE_DOUBLE   ,H5T_IEEE_F32BE,H5T_IEEE_F32LE,H5T_IEEE_F64BE,H5T_IEEE_F64LE]
       H5T_NATIVE_INTEGERS         =[H5T_NATIVE_INTEGER  ,H5T_STD_I32BE ,H5T_STD_I32LE ,H5T_STD_I64BE ,H5T_STD_I64LE ]
       H5T_NATIVE_UNSIGNED_INTEGERS=[H5T_STD_U32BE       ,H5T_STD_U32LE                                              ]
       H5T_NATIVE_INTEGER_8S       =[H5T_INTEGER8,H5T_STD_I64BE ,H5T_STD_I64LE                               ]
       H5T_NATIVE_INTEGER_8AS(1:3) =H5T_NATIVE_INTEGERS(1:3)
       H5T_NATIVE_INTEGER_8AS(4:5) =H5T_NATIVE_UNSIGNED_INTEGERS
       H5T_NATIVE_INTEGER_8AS(6:8) =H5T_NATIVE_INTEGER_8S

       ! Create vlen datatypes.
       call h5tvlen_create_f(H5T_NATIVE_DOUBLE   ,H5T_VLEN_DOUBLE     (1),errorCode) 
       if (errorCode < 0) call Error_Report('failed to create vlen double HDF5 datatype'     //{introspection:location})
       call h5tvlen_create_f(H5T_VLEN_DOUBLE  (1),H5T_VLEN_VLEN_DOUBLE(1),errorCode) 
       if (errorCode < 0) call Error_Report('failed to create vlen-veln double HDF5 datatype'//{introspection:location})
       call h5tvlen_create_f(H5T_INTEGER8        ,H5T_VLEN_INTEGER8   (1),errorCode) 
       if (errorCode < 0) call Error_Report('failed to create vlen integer8 HDF5 datatype'   //{introspection:location})

       ! Initialize our OpenMP lock.
       call ioHDF5AccessInitialize()
       
       ! Flag that the hdf5 system is now initialized.
       hdf5IsInitialized=.true.
    end if
    initializationsCount=initializationsCount+1
    return
  end subroutine IO_HDF5_Initialize

  subroutine IO_HDF5_Uninitialize
    !!{
    Uninitialize the HDF5 subsystem.
    !!}
    use :: Error, only : Error_Report
    use :: HDF5 , only : h5close_f   , h5tclose_f
    implicit none
    integer :: errorCode

#ifdef DEBUGHDF5
    call IO_HDF5_Assert_In_Critical()
#endif

    if (hdf5IsInitialized) then
       initializationsCount=initializationsCount-1
       if (initializationsCount == 0) then
          call h5tclose_f(H5T_VLEN_DOUBLE     (1),errorCode)
          if (errorCode < 0) call Error_Report('failed to close vlen double datatype'     //{introspection:location})
          call h5tclose_f(H5T_VLEN_VLEN_DOUBLE(1),errorCode)
          if (errorCode < 0) call Error_Report('failed to close vlen-vlen double datatype'//{introspection:location})
          call h5tclose_f(H5T_VLEN_INTEGER8   (1),errorCode)
          if (errorCode < 0) call Error_Report('failed to close vlen integer8 datatype'   //{introspection:location})
          call h5close_f (                   errorCode)
          if (errorCode < 0) call Error_Report('failed to uninitialize HDF5 subsystem'    //{introspection:location})
          hdf5IsInitialized=.false.
       end if
    end if
    return
  end subroutine IO_HDF5_Uninitialize

  subroutine IO_HDF_Assert_Is_Initialized()
    !!{
    Check if this module has been initialized.
    !!}
    use :: Error, only : Error_Report
    implicit none

#ifdef DEBUGHDF5
    call IO_HDF5_Assert_In_Critical()
#endif

    if (.not.hdf5IsInitialized) call Error_Report('HDF5 IO module has not been initialized'//{introspection:location})
    return
  end subroutine IO_HDF_Assert_Is_Initialized

#ifdef DEBUGHDF5
  subroutine IO_HDF5_Assert_In_Critical()
    !!{
    Assert that we are in an {\normalfont \ttfamily HDF5\_Access} OpenMP critical block.
    !!}
    use :: Error, only : Error_Report
    implicit none

    if (.not.inCritical) call Error_Report('HDF5 functions accessed outside of critical block'//{introspection:location})
    return
  end subroutine IO_HDF5_Assert_In_Critical

  subroutine IO_HDF5_Start_Critical()
    !!{
    Record that we have entered an {\normalfont \ttfamily HDF5\_Access} OpenMP critical block.
    !!}
    implicit none

    inCritical=.true.
    return
  end subroutine IO_HDF5_Start_Critical

  subroutine IO_HDF5_End_Critical()
    !!{
    Record that we have left an {\normalfont \ttfamily HDF5\_Access} OpenMP critical block.
    !!}
    implicit none

    inCritical=.false.
    return
  end subroutine IO_HDF5_End_Critical
#endif

  subroutine IO_HDF5_Set_Defaults(chunkSize,compressionLevel)
    !!{
    Sets the compression level and chunk size for dataset output.
    !!}
    use :: Error      , only : Error_Report
    use :: HDF5       , only : HSIZE_T
    use :: HDF5_Access, only : hdf5Access
    implicit none
    integer(kind=HSIZE_T), intent(in   ), optional :: chunkSize
    integer              , intent(in   ), optional :: compressionLevel

    !$ call hdf5Access%set()
    if (present(chunkSize)) then
       if (chunkSize        ==  0_hsize_t) call Error_Report('zero chunksize is invalid'        //{introspection:location})
       if (chunkSize        <  -1_hsize_t) call Error_Report('chunksize less than -1 is invalid'//{introspection:location})
       hdf5ChunkSize=chunkSize
    end if
    if (present(compressionLevel)) then
       if (compressionLevel <  -1 .or. compressionLevel > 9) call Error_Report('compression level must be in range -1 to 9'//{introspection:location})
       hdf5CompressionLevel=compressionLevel
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine IO_HDF5_Set_Defaults

  !! Utility routines.

  impure elemental subroutine IO_HDF5_Finalize(self)
    !!{
    Finalize an HDF5 object.
    !!}
    use :: File_Utilities    , only : File_Remove
    use :: HDF5_Access       , only : hdf5Access
    use :: Display           , only : displayIndent     , displayMessage  , displayUnindent, verbosityLevelSilent
    use :: Error             , only : Error_Report
    use :: HDF5              , only : h5f_obj_all_f     , h5aclose_f      , h5dclose_f     , h5fclose_f          , &
          &                           h5fget_obj_count_f, h5fget_obj_ids_f, h5gclose_f     , h5iget_name_f       , &
          &                           hid_t             , size_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)    , char
    use :: String_Handling   , only : operator(//)
    implicit none
    type      (hdf5Object               ), intent(inout)               :: self
    integer   (hid_t                    ), allocatable  , dimension(:) :: openObjectIDs
    integer   (size_t                   ), parameter                   :: objectNameSizeMaximum=1024
    integer                                                            :: errorCode
    integer   (size_t                   )                              :: i                         , objectNameSize        , &
         &                                                                openObjectCount           , nonRootOpenObjectCount
    type      (varying_string           )                              :: message
    character (len=objectNameSizeMaximum)                              :: objectName
    !$ logical                                                         :: haveLock
    
    if (self%objectManager%count() == 1) then
       !$ haveLock=hdf5Access%ownedByThread()
       !$ if (.not.haveLock) call hdf5Access%set  ()
       ! Close the object.
       select case (self%hdf5ObjectType)
       case (hdf5ObjectTypeFile     )
          ! Do not close file objects as these are a shared resource.
       case (hdf5ObjectTypeGroup    )
          call h5gclose_f(self%objectID,errorCode)
          if (errorCode /= 0) then
             message="unable to close group object '"//self%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       case (hdf5ObjectTypeDataset  )
          call h5dclose_f(self%objectID,errorCode)
          if (errorCode /= 0) then
             message="unable to close dataset object '"//self%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       case (hdf5ObjectTypeAttribute)
          call h5aclose_f(self%objectID,errorCode)
          if (errorCode /= 0) then
             message="unable to close attribute object '"//self%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end select
       !$ if (.not.haveLock) call hdf5Access%unset()
       nullify(self%parentObject)
    end if
    ! If this is the last reference to the file, close it now.
    if (self%fileManager%count() == 1) then
       !$ haveLock=hdf5Access%ownedByThread()
       !$ if (.not.haveLock) call hdf5Access%set  ()
       ! Check for still-open objects.
       call h5fget_obj_count_f(self%fileID,h5f_obj_all_f,openObjectCount,errorCode)
       if (errorCode /= 0) then
          message="unable to count open objects in file object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       allocate(openObjectIDs(openObjectCount))
       call h5fget_obj_ids_f(self%fileID,h5f_obj_all_f,openObjectCount,openObjectIDs,errorCode)
       if (errorCode /= 0) then
          message="unable to get IDs of open objects in file object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       nonRootOpenObjectCount=0
       if (openObjectCount > 1) then
          do i=1,openObjectCount
             call h5iget_name_f(openObjectIDs(i),objectName,objectNameSizeMaximum,objectNameSize,errorCode)
             if (errorCode /= 0) then
                message="unable to get name of open object in file object '"//self%objectName//"'"
                call Error_Report(message//self%locationReport()//{introspection:location})
             end if
             
             if (trim(objectName) /= "/") nonRootOpenObjectCount=nonRootOpenObjectCount+1
          end do
       end if
       if (nonRootOpenObjectCount > 0 .and. openObjectCount-nonRootOpenObjectCount == 1) then
          message=""
          message=message//nonRootOpenObjectCount//" open object(s) remain in file object '"//self%objectName//"'"
          call displayIndent('Problem closing HDF5 file',verbosityLevelSilent)
          call displayMessage(message,verbosityLevelSilent)
          do i=1,openObjectCount
             call h5iget_name_f(openObjectIDs(i),objectName,objectNameSizeMaximum,objectNameSize,errorCode)
             if (errorCode /= 0) then
                message="unable to get name of open object in file object '"//self%objectName//"'"
                call Error_Report(message//self%locationReport()//{introspection:location})
             end if
             message="Object: "//trim(objectName)//" ["
             message=message//openObjectIDs(i)//"]"
             if (trim(objectName) /= "/") call displayMessage(message,verbosityLevelSilent)
          end do
          call displayUnindent('done',verbosityLevelSilent)
       end if
       call h5fclose_f(self%fileID,errorCode)
       if (errorCode /= 0) then
          message="unable to close file object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       if (self%isTemporary) call File_Remove(char(self%objectName))
       ! Uninitialize the HDF5 library (will only uninitialize if this is the last file to be closed).
       call IO_HDF5_Uninitialize()
       !$ if (.not.haveLock) call hdf5Access%unset()
    end if
    return
  end subroutine IO_HDF5_Finalize
  
  logical function IO_HDF5_Is_Open(self)
    !!{
    Returns true if {\normalfont \ttfamily self} is open.
    !!}
    implicit none
    class(hdf5Object), intent(in   ) :: self

    IO_HDF5_Is_Open=self%isOpenValue
    return
  end function IO_HDF5_Is_Open

  integer function IO_HDF5_Object_Type(self)
    !!{
    Returns the object type for {\normalfont \ttfamily self}.
    !!}
    implicit none
    class(hdf5Object), intent(in   ) :: self

    IO_HDF5_Object_Type=self%hdf5ObjectType
    return
  end function IO_HDF5_Object_Type

  subroutine IO_HDF5_Assign(to,from)
    !!{
    Assignment operator for the {\normalfont \ttfamily hdf5Object} class.
    !!}
    implicit none
    class(hdf5Object), intent(  out) :: to
    class(hdf5Object), intent(in   ) :: from

    to%isOpenValue         =  from%isOpenValue
    to%isOverwritable      =  from%isOverwritable
    to%readOnly            =  from%readOnly
    to%isTemporary         =  from%isTemporary
    to%objectID            => from%objectID
    to%objectManager       =  from%objectManager
    to%fileID              => from%fileID
    to%fileManager         =  from%fileManager
    to%objectLocation      =  from%objectLocation
    to%objectFile          =  from%objectFile
    to%objectName          =  from%objectName
    to%hdf5ObjectType      =  from%hdf5ObjectType
    to%chunkSize           =  from%chunkSize
    to%compressionLevel    =  from%compressionLevel
    to%chunkSizeSet        =  from%chunkSizeSet
    to%compressionLevelSet =  from%compressionLevelSet
    to%parentObject        => from%parentObject
    return
  end subroutine IO_HDF5_Assign
  
  function IO_HDF5_Name(self) result (nameOfObject)
    !!{
    Returns the path to {\normalfont \ttfamily self}.
    !!}
    implicit none
    class(hdf5Object    ), intent(in   ) :: self
    type (varying_string)                :: nameOfObject

    nameOfObject=self%objectName
    return
  end function IO_HDF5_Name

  function IO_HDF5_Path_To(self,includeFileName) result (pathToObject)
    !!{
    Returns the path to {\normalfont \ttfamily self}.
    !!}
    use :: ISO_Varying_String, only : operator(//), assignment(=), operator(/=)
    implicit none
    class  (hdf5Object    ), intent(in   )           :: self
    logical                , intent(in   ), optional :: includeFileName
    type   (varying_string)                          :: pathToObject
    !![
    <optionalArgument name="includeFileName" defaultsTo=".true." />
    !!]

    if (includeFileName_) then
       pathToObject=self%objectFile
    else
       pathToObject=""
    end if
    if (self%objectLocation /= "") then
       if (pathToObject /= "") pathToObject=pathToObject//"/"
       pathToObject=pathToObject//self%objectLocation
    end if
    if (self%objectName     /= "") then
       if (pathToObject /= "") pathToObject=pathToObject//"/"
       pathToObject=pathToObject//self%objectName
    end if
    return
  end function IO_HDF5_Path_To

  function IO_HDF5_File_Name(self) result (fileName)
    !!{
    Returns the name of the file containing {\normalfont \ttfamily self}.
    !!}
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    class(hdf5Object    ), intent(in   ), target  :: self
    type (varying_string)                         :: fileName
    type (hdf5Object    )               , pointer :: parent

    parent  => self
    do while (parent%hdf5ObjectType /= hdf5ObjectTypeFile)
       parent => parent%parentObject
    end do
    fileName=parent%objectName
    return
  end function IO_HDF5_File_Name

  function IO_HDF5_Location_Report(self) result (report)
    !!{
    Returns a report on the location of {\normalfont \ttfamily self} suitable for inclusion in an error message.
    !!}
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    class(hdf5Object    ), intent(in   ), target :: self
    type (varying_string)                        :: report

    report=char(10)//"   object is located in: "//self%pathTo()
    return
  end function IO_HDF5_Location_Report

  subroutine IO_HDF5_Flush(self)
    !!{
    Flush an HDF5 file to disk.
    !!}
    use :: Display           , only : displayMessage
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5F_Scope_Local_F, h5fflush_f
    use :: ISO_Varying_String, only : assignment(=)    , operator(//)
    implicit none
    class  (hdf5Object    ), intent(inout) :: self
    type   (varying_string)                :: message
    integer                                :: errorCode

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the object is open.
    if (.not.self%isOpenValue) then
       message="Attempt to flush unopen HDF5 object '"//self%objectName//"'"
       call displayMessage(message)
       return
    end if

    ! Flush to file.
    call h5fflush_f(self%objectID,H5F_Scope_Local_F,errorCode)
    if (errorCode /= 0) then
       message="unable to flush object '"//self%objectName//"' to file"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Flush

  subroutine IO_HDF5_Remove(self,objectName)
    !!{
    Remove the named object.
    !!}
    use :: Display           , only : displayMessage
    use :: Error             , only : Error_Report
    use :: HDF5              , only : h5ldelete_f
    use :: ISO_Varying_String, only : assignment(=) , operator(//)
    implicit none
    class    (hdf5Object    ), intent(inout) :: self
    character(len=*         ), intent(in   ) :: objectName
    type     (varying_string)                :: message
    integer                                  :: errorCode

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized()
    ! Check that the object is open.
    if (.not.self%isOpenValue) then
       message="Attempt to remove object from unopen HDF5 object '"//self%objectName//"'"
       call displayMessage(message)
       return
    end if
    ! Remove the object.
    call h5ldelete_f(self%objectID,objectName,errorCode)
    if (errorCode /= 0) then
       message="unable to remove '"//objectname//"' from object '"//self%objectName//"' to file"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Remove

  function IO_HDF5_Character_Types(stringLength)
    !!{
    Return datatypes for character data of a given length. Types are for Fortran native and C native types.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_CHARACTER, H5T_STR_NULLPAD_F, H5T_STR_NULLTERM_F, H5T_STR_SPACEPAD_F, &
          &                           H5Tset_strpad_f     , HID_T            , h5tcopy_f         , h5tset_size_f     , &
          &                           size_t
    use :: ISO_Varying_String, only : assignment(=)       , operator(//)
    implicit none
    integer(kind=HID_T    ), dimension(6)  :: IO_HDF5_Character_Types
    integer                , intent(in   ) :: stringLength
    integer                                :: errorCode
    type   (varying_string)                :: message

    call h5tcopy_f(H5T_NATIVE_CHARACTER,IO_HDF5_Character_Types(1),errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_size_f(IO_HDF5_Character_Types(1),int(stringLength,size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_strpad_f(IO_HDF5_Character_Types(1),H5T_STR_SPACEPAD_F,errorCode)
    if (errorCode < 0) then
       message="unable to set padding"
       call Error_Report(message//{introspection:location})
    end if
    call h5tcopy_f(H5T_C_S1_Get(),IO_HDF5_Character_Types(2),errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_size_f(IO_HDF5_Character_Types(2),int(stringLength,size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_strpad_f(IO_HDF5_Character_Types(2),H5T_STR_SPACEPAD_F,errorCode)
    if (errorCode < 0) then
       message="unable to set padding"
       call Error_Report(message//{introspection:location})
    end if
    call h5tcopy_f(H5T_NATIVE_CHARACTER,IO_HDF5_Character_Types(3),errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_size_f(IO_HDF5_Character_Types(3),int(stringLength,size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_strpad_f(IO_HDF5_Character_Types(3),H5T_STR_NULLPAD_F,errorCode)
    if (errorCode < 0) then
       message="unable to set padding"
       call Error_Report(message//{introspection:location})
    end if
    call h5tcopy_f(H5T_C_S1_Get(),IO_HDF5_Character_Types(4),errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_size_f(IO_HDF5_Character_Types(4),int(stringLength,size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_strpad_f(IO_HDF5_Character_Types(4),H5T_STR_NULLPAD_F,errorCode)
    if (errorCode < 0) then
       message="unable to set padding"
       call Error_Report(message//{introspection:location})
    end if
    call h5tcopy_f(H5T_NATIVE_CHARACTER,IO_HDF5_Character_Types(5),errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_size_f(IO_HDF5_Character_Types(5),int(stringLength,size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_strpad_f(IO_HDF5_Character_Types(5),H5T_STR_NULLTERM_F,errorCode)
    if (errorCode < 0) then
       message="unable to set padding"
       call Error_Report(message//{introspection:location})
    end if
    call h5tcopy_f(H5T_C_S1_Get(),IO_HDF5_Character_Types(6),errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_size_f(IO_HDF5_Character_Types(6),int(stringLength,size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size"
       call Error_Report(message//{introspection:location})
    end if
    call h5tset_strpad_f(IO_HDF5_Character_Types(6),H5T_STR_NULLTERM_F,errorCode)
    if (errorCode < 0) then
       message="unable to set padding"
       call Error_Report(message//{introspection:location})
    end if
    return
  end function IO_HDF5_Character_Types

  impure elemental subroutine hdf5VarDoubleDestructor(self)
    !!{
    Destructor for the variable-length double type.
    !!}
    implicit none
    type(hdf5VarDouble), intent(inout) :: self

    if (associated(self%row)) deallocate(self%row)
    return
  end subroutine hdf5VarDoubleDestructor

  impure elemental subroutine hdf5VarDouble2DDestructor(self)
    !!{
    Destructor for the variable-length 2D double type.
    !!}
    implicit none
    type(hdf5VarDouble2D), intent(inout) :: self

    if (associated(self%row)) deallocate(self%row)
    return
  end subroutine hdf5VarDouble2DDestructor

  impure elemental subroutine hdf5VarInteger8Destructor(self)
    !!{
    Destructor for the variable-length integer-8 type.
    !!}
    implicit none
    type(hdf5VarInteger8), intent(inout) :: self

    if (associated(self%row)) deallocate(self%row)
    return
  end subroutine hdf5VarInteger8Destructor
  
  !! File routines.

  function hdf5FileOpenVarStr(fileName,overWrite,readOnly,objectsOverwritable,chunkSize,compressionLevel,sieveBufferSize,useLatestFormat,cacheElementsCount,cacheSizeBytes,isTemporary) result(self)
    !!{
    Constructor for HDF5 object. Will open a file and return an appropriate HDF5 object.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    type   (hdf5Object    )                          :: self
    type   (varying_string), intent(in   )           :: fileName
    logical                , intent(in   ), optional :: objectsOverwritable, overWrite         , readOnly      , isTemporary
    integer(kind=hsize_t  ), intent(in   ), optional :: chunkSize
    integer(kind=size_t   ), intent(in   ), optional :: sieveBufferSize    , cacheElementsCount, cacheSizeBytes
    integer                , intent(in   ), optional :: compressionLevel
    logical                , intent(in   ), optional :: useLatestFormat

    self=hdf5Object(char(fileName),overWrite,readOnly,objectsOverwritable,chunkSize,compressionLevel,sieveBufferSize,useLatestFormat,cacheElementsCount,cacheSizeBytes,isTemporary)
    return
  end function hdf5FileOpenVarStr

  function hdf5FileOpenChar(fileName,overWrite,readOnly,objectsOverwritable,chunkSize,compressionLevel,sieveBufferSize,useLatestFormat,cacheElementsCount,cacheSizeBytes,isTemporary) result(self)
    !!{
    Constructor for HDF5 object. Will open a file and return an appropriate HDF5 object.
    !!}
    use :: File_Utilities    , only : File_Exists        , File_Name_Expand
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5F_ACC_RDONLY_F   , H5F_ACC_RDWR_F        , H5F_ACC_TRUNC_F       , H5F_CLOSE_SEMI_F       , &
         &                            H5F_LIBVER_V18_F   , H5F_LIBVER_LATEST_F   , H5P_FILE_ACCESS_F     , h5fcreate_f            , &
         &                            h5fopen_f          , h5pclose_f            , h5pcreate_f           , h5pset_cache_f         , &
         &                            h5pset_fapl_stdio_f, h5pset_fclose_degree_f, h5pset_libver_bounds_f, h5pset_sieve_buf_size_f, &
         &                            hid_t              , hsize_t               , size_t
    use :: ISO_Varying_String, only : assignment(=)      , len                   , operator(//)          , trim                   , &
         &                            char
    implicit none
    type     (hdf5Object    )                          :: self
    character(len=*         ), intent(in   )           :: fileName
    logical                  , intent(in   ), optional :: objectsOverwritable, overWrite         , readOnly      , isTemporary
    integer  (kind=hsize_t  ), intent(in   ), optional :: chunkSize
    integer  (kind=size_t   ), intent(in   ), optional :: sieveBufferSize    , cacheElementsCount, cacheSizeBytes
    integer                  , intent(in   ), optional :: compressionLevel
    logical                  , intent(in   ), optional :: useLatestFormat
    class    (*             ), pointer                 :: dummyPointer_
    integer                                            :: errorCode          , fileAccess
    logical                                            :: overWriteActual
    type     (varying_string)                          :: message            , fileName_
    integer  (kind=hid_t    )                          :: accessList

    ! Initialize the HDF5 library.
    call IO_HDF5_Initialize()
    ! Expand the file name.
    fileName_=File_Name_Expand(fileName)
    ! Store the location and name of this object.
    self%objectFile    =trim(fileName_)
    self%objectLocation=""
    self%objectName    =""
    ! Mark whether this file is temporary.
    if (present(isTemporary)) then
       self%isTemporary=isTemporary
    else
       self%isTemporary=.false.
    end if
    ! Check if overwriting is allowed.
    if (present(overWrite)) then
       overWriteActual=overWrite
    else
       overWriteActual=.false.
    end if
    ! Create an access list.
    call h5pcreate_f(H5P_FILE_ACCESS_F,accessList,errorCode)
    if (errorCode /= 0) then
       message="failed to create file access list HDF5 file '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5pset_fclose_degree_f(accessList,H5F_CLOSE_SEMI_F,errorCode)
    if (errorCode /= 0) then
       message="failed to set close degree for HDF5 file '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Specify file driver (buffered I/O).
    call h5pset_fapl_stdio_f(accessList,errorCode)
    if (errorCode /= 0) then
       message="failed to set I/O driver for HDF5 file '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Set sieve buffer size.
    if (present(sieveBufferSize)) then
       call h5pset_sieve_buf_size_f(accessList,sieveBufferSize,errorCode)
       if (errorCode /= 0) then
          message="failed to set sieve buffer size for HDF5 file '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Set file format.
    if (present(useLatestFormat)) then
       if (useLatestFormat) then
          call h5pset_libver_bounds_f(accessList,H5F_LIBVER_LATEST_F,H5F_LIBVER_LATEST_F,errorCode)
       else
          call h5pset_libver_bounds_f(accessList,H5F_LIBVER_V18_F   ,H5F_LIBVER_LATEST_F,errorCode)
       end if
       if (errorCode /= 0) then
          message="failed to set file format for HDF5 file '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       call    h5pset_libver_bounds_f(accessList,H5F_LIBVER_V18_F   ,H5F_LIBVER_LATEST_F,errorCode)
    end if
    if (errorCode /= 0) then
       message="failed to set file format for HDF5 file '"//self%objectName//"'"
       call Error_Report(message//{introspection:location})
    end if
    if (present(cacheElementsCount).or.present(cacheSizeBytes)) then
       if (.not.(present(cacheElementsCount).and.present(cacheSizeBytes))) call Error_Report('both or neither of "cacheElementsCount" and "cacheSizeBytes" must be specified'//{introspection:location})
       call h5pset_cache_f(accessList,0,cacheElementsCount,cacheSizeBytes,0.75,errorCode)
    end if
    ! Allocate object and file IDs.
    allocate(self%objectID)
    allocate(self%fileID  )
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%objectID
    self%objectManager =  resourceManager(dummyPointer_)
    dummyPointer_      => self%fileID
    self%fileManager   =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    ! Check if the file exists.
    if (File_Exists(fileName_).and..not.overWriteActual) then
       ! Determine access for file.
       if (present(readOnly)) then
          self%readOnly=readOnly
          if (readOnly) then
             fileAccess=H5F_ACC_RDONLY_F
          else
             fileAccess=H5F_ACC_RDWR_F
          end if
       else
          self%readOnly=.false.
          fileAccess=H5F_ACC_RDWR_F
       end if
       ! Attempt to open the file.
       call h5fopen_f(char(fileName_),fileAccess,self%objectID,errorCode,access_prp=accessList)
       if (errorCode /= 0) then
          message="failed to open HDF5 file '"//self%objectFile//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! If read only was specified, creating the file is not allowed.
       if (present(readOnly)) then
          if (readOnly) then
             message="can not create/overwrite read only file '"//self%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
       ! Attempt to create the file.
       call h5fcreate_f(char(fileName_),H5F_ACC_TRUNC_F,self%objectID,errorCode,access_prp=accessList)
       if (errorCode /= 0) then
          message="failed to create HDF5 file '"//self%objectFile//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Finished with our property list.
    call h5pclose_f(accessList,errorCode)
    if (errorCode /= 0) then
       message="failed to close access property list for HDF5 file '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Store the file ID.
    self%fileID     =self%objectID
    ! Mark this object as open.
    self%isOpenValue=.true.
    ! Object has no parent.
    self%parentObject => null()
    ! Mark this object as a file object.
    self%hdf5ObjectType=hdf5ObjectTypeFile
    ! Set the chunk size if provided.
    if (present(chunkSize)) then
       self%chunkSizeSet=.true.
       self%chunkSize   =chunkSize
    else
       self%chunkSizeSet=.false.
    end if
    ! Set the compression level if provided.
    if (present(compressionLevel)) then
       self%compressionLevelSet=.true.
       self%compressionLevel   =compressionLevel
    else
       self%compressionLevelSet=.false.
    end if
    ! Mark whether objects are overwritable.
    if (present(objectsOverwritable)) then
       self%isOverwritable=objectsOverwritable
    else
       self%isOverwritable=.false.
    end if
    return
  end function hdf5FileOpenChar

  !! Group routines.

  subroutine IO_HDF5_Open_Group_Path(inObject,groupPath,groupObjects)
    !!{
    Open all HDF5 groups in the given path and return objects for all of them.
    !!}
    use :: String_Handling   , only : String_Split_Words, String_Count_Words
    use :: ISO_Varying_String, only : char
    implicit none
    class    (hdf5Object    ), intent(in   ), target                      :: inObject
    character(len=*         ), intent(in   )                              :: groupPath
    type     (hdf5Object    ), intent(inout), allocatable  , dimension(:) :: groupObjects
    type     (varying_string)               , allocatable  , dimension(:) :: groupNames
    integer                                                               :: i           , countGroups

    countGroups=String_Count_Words(groupPath,"/")
    allocate(groupNames  (countGroups))
    allocate(groupObjects(countGroups))
    call String_Split_Words(groupNames,groupPath,"/")
    do i=1,size(groupNames)
       if (i == 1) then
          groupObjects(i)=inObject         %openGroup(char(groupNames(i)))
       else
          groupObjects(i)=groupObjects(i-1)%openGroup(char(groupNames(i)))
       end if
    end do
    return
  end subroutine IO_HDF5_Open_Group_Path

  function IO_HDF5_Open_Group(inObject,groupName,comment,objectsOverwritable,overwriteOverride,chunkSize,compressionLevel,attributesCompactMaxiumum) result (self)
    !!{
    Open an HDF5 group and return an appropriate HDF5 object. The group name can be provided as an input parameter or, if
    not provided, will be taken from the stored object name in {\normalfont \ttfamily self}. The location at which to open the group is
    taken from either {\normalfont \ttfamily inObject} or {\normalfont \ttfamily inPath}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , h5gcreate_f               , h5gopen_f , hsize_t            , &
         &                            h5pcreate_f  , h5pset_attr_phase_change_f, h5pclose_f, H5P_GROUP_CREATE_F
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    implicit none
    type     (hdf5Object    )                          :: self
    character(len=*         ), intent(in   )           :: groupName
    character(len=*         ), intent(in   ), optional :: comment
    logical                  , intent(in   ), optional :: objectsOverwritable, overwriteOverride
    integer  (hsize_t       ), intent(in   ), optional :: chunkSize
    integer                  , intent(in   ), optional :: compressionLevel   , attributesCompactMaxiumum
    class    (hdf5Object    ), intent(in   ), target   :: inObject
    ! <HDF5> Why are "message" and "locationPath" saved? Because if they are not then they get dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (varying_string), save                    :: locationPath       , message
    integer                                            :: errorCode
    integer  (kind=HID_T    )                          :: locationID         , propertyList
    class    (*             ), pointer                 :: dummyPointer_

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Decide where to open the group.
    if (.not.inObject%isOpenValue) then
       message="attempt to open group '"//trim(groupName)//"' in unopen object '"//inObject%objectName//"'"
       call Error_Report(message//inObject%locationReport()//{introspection:location})
    end if
    locationID  =inObject%objectID
    locationPath=inObject%pathTo()

    ! Set the parent for the group.
    select type (inObject)
    type is (hdf5Object)
       self%parentObject => inObject
    end select
    ! Obtain a reference to the file ID.
    self%fileID      => inObject%fileID
    self%fileManager =  inObject%fileManager
    ! Create an ID for this group.
    allocate(self%objectID)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_ => self%objectID
    self%objectManager=resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    ! Check if the group exists.
    if (inObject%hasGroup(groupName)) then
       ! Open the group.
       call h5gopen_f(locationID,trim(groupName),self%objectID,errorCode)
       if (errorCode /= 0) then
          message="failed to open group '"//trim(groupName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
    else
       ! Create a group.
       call h5pcreate_f(H5P_GROUP_CREATE_F,propertyList,errorCode)
       if (errorCode < 0) then
          message="unable to create property list for group '"//trim(groupName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Set the number of attributes allowed in compact
       ! storage. (Setting this to zero will force dense storage for
       ! all attributes, which enable storing of attributes with
       ! length exceeding 64KB.)
       if (present(attributesCompactMaxiumum)) then
          call h5pset_attr_phase_change_f(propertyList,attributesCompactMaxiumum,attributesCompactMaxiumum,errorCode)
          if (errorCode /= 0) then
             message="failed to set attribute phase change for group '"//trim(groupName)
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
       call h5gcreate_f(locationID,trim(groupName),self%objectID,errorCode,gcpl_id=propertyList)
       if (errorCode < 0) then
          message="failed to make group '"//trim(groupName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       call h5pclose_f(propertyList,errorCode)
       if (errorCode /= 0) then
          message="failed to close property list for group '"//trim(groupName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Mark this object as open.
    self%isOpenValue=.true.

    ! Mark this object as a file object.
    self%hdf5ObjectType=hdf5ObjectTypeGroup

    ! Store the name and location of the object.
    self%objectFile    =self%parentObject%objectFile
    self%objectLocation=self%parentObject%pathTo    (includeFileName=.false.)
    self%objectName    =trim(groupName)

    ! Set the chunk size if provided.
    if (present(chunkSize)) then
       self%chunkSizeSet   =.true.
       self%chunkSize      =chunkSize
    else
       ! No chunk size provided. See if we can inherit one from the parent object.
       if (self%parentObject%chunkSizeSet) then
          self%chunkSizeSet=.true.
          self%chunkSize   =self%parentObject%chunkSize
       else
          self%chunkSizeSet=.false.
       end if
    end if

    ! Set the compression level if provided.
    if (present(compressionLevel)) then
       self%compressionLevelSet=.true.
       self%compressionLevel   =compressionLevel
    else
       ! No compression level provided. See if we can inherit one from the parent object.
       if (self%parentObject%compressionLevelSet) then
          self%compressionLevelSet=.true.
          self%compressionLevel   =self%parentObject%compressionLevel
       else
          self%compressionLevelSet=.false.
       end if
    end if

    ! Mark whether objects are overwritable.
    if (present(objectsOverwritable)) then
       if (.not.present(overwriteOverride).or..not.overwriteOverride) then
          if (objectsOverwritable.and..not.self%parentObject%isOverwritable) then
             message="cannot make objects in '"//trim(groupName)//"' overwritable as objects in parent '"&
                  &//self%parentObject%objectName//"' are not overwritable"
             call Error_Report(message//inObject%locationReport()//{introspection:location})
          end if
       end if
       self%isOverwritable=objectsOverwritable
    else
       self%isOverwritable=self%parentObject%isOverwritable
    end if

    ! Set the comment for this group.
    if (present(comment) .and. len_trim(comment) > 0 .and. .not.self%hasAttribute('comment')) call self%writeAttribute(trim(comment),'comment')
    return
  end function IO_HDF5_Open_Group

  logical function IO_HDF5_Has_Group(self,groupName)
    !!{
    Check if {\normalfont \ttfamily self} has a group with the given {\normalfont \ttfamily groupName}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : h5eset_auto_f, h5gget_info_by_name_f
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    implicit none
    class    (hdf5Object    ), intent(in   ) :: self
    character(len=*         ), intent(in   ) :: groupName
    integer                                  :: creationOrderMaximum, errorCode, linkCount, &
         &                                      storageType
    type     (varying_string)                :: message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="object '"//self%objectName//"' in not open"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object exists.
    call h5eset_auto_f(0,errorCode)
    if (errorCode /= 0) then
       message="failed to switch HDF5 error report off"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5gget_info_by_name_f(self%objectID,trim(groupName),storageType,linkCount,creationOrderMaximum,errorCode)
    IO_HDF5_Has_Group=(errorCode == 0)
    call h5eset_auto_f(1,errorCode)
    if (errorCode /= 0) then
       message="failed to switch HDF5 error report on"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end function IO_HDF5_Has_Group

  !! Attribute routines.

  function IO_HDF5_Open_Attribute(inObject,attributeName,attributeDataType,attributeDimensions,isOverwritable,useDataType) result(self)
    !!{
    Open an attribute in {\normalfont \ttfamily inObject}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_CHARACTER, H5T_NATIVE_DOUBLE , H5T_NATIVE_INTEGER, h5screate_simple_f, &
          &                           HID_T               , HSIZE_T           , h5acreate_f       , h5aopen_f         , &
          &                           h5sclose_f
    use :: ISO_Varying_String, only : assignment(=)       , operator(//)      , operator(/=)
    implicit none
    class    (hdf5Object    )              , intent(in   ), target   :: inObject
    type     (hdf5Object    )                                        :: self
    character(len=*         )              , intent(in   )           :: attributeName
    integer                                , intent(in   ), optional :: attributeDataType
    integer  (kind=HSIZE_T  ), dimension(:), intent(in   ), optional :: attributeDimensions
    logical                                , intent(in   ), optional :: isOverwritable
    integer  (kind=HID_T    )              , intent(in   ), optional :: useDataType
    integer                                                          :: attributeRank            , errorCode
    integer  (kind=HID_T    )                                        :: dataSpaceID              , dataTypeID, &
         &                                                              locationID
    integer  (kind=HSIZE_T  ), dimension(7)                          :: attributeDimensionsActual
    type     (varying_string)                                        :: locationPath             , message
    class    (*             ), pointer                               :: dummyPointer_

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Ensure that the object is already open.
    if (inObject%isOpenValue) then
       locationID                                     =inObject    %objectID
       locationPath                                   =inObject    %objectFile
       if (inObject%objectLocation /= "") locationPath=locationPath           //"/"//inObject%objectLocation
       if (inObject%objectName     /= "") locationPath=locationPath           //"/"//inObject%objectName
       select type (inObject)
       type is (hdf5Object)
          self%parentObject => inObject
       end select
    else
       message="attempt to open attribute '"//trim(attributeName)//"' in unopen object '"//inObject%objectName//"'"
       call Error_Report(message//inObject%locationReport()//{introspection:location})
    end if

    ! Determine the rank and dimensions.
    if (present(attributeDimensions)) then
       ! Open data space with the desired dimensions.
       attributeRank=size(attributeDimensions)
       if (attributeRank > 7) then
          message="attributes of rank greater than 7 are not supported - attribute in question is '"//trim(attributeName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       attributeDimensionsActual(1:attributeRank)=attributeDimensions
    else
       ! No dimensions specified, assume a scalar.
       attributeRank=0
    end if
    ! Obtain a reference to the file ID.
    self%fileID      => inObject%fileID
    self%fileManager =  inObject%fileManager
    ! Create an ID for this attribute.
    allocate(self%objectID)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_ => self%objectID
    self%objectManager=resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    ! Check if the attribute exists.
    if (IO_HDF5_Has_Attribute(inObject,attributeName)) then
       ! Open the attribute.
       call h5aopen_f(locationID,trim(attributeName),self%objectID,errorCode)
       if (errorCode /= 0) then
          message="failed to open attribute '"//trim(attributeName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
    else
       ! Ensure that a data type was specified.
       if (.not.present(attributeDataType)) then
          message="no datatype was specified for attribute '"//trim(attributeName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       ! Open a suitable dataspace.
       call h5screate_simple_f(attributeRank,attributeDimensionsActual,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="unable to open dataspace for attribute '"//trim(attributeName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       ! Determine the data type.
       if (present(useDataType)) then
          dataTypeID=useDataType
       else
          select case (attributeDataType)
          case (hdf5DataTypeInteger       )
             dataTypeID=H5T_NATIVE_INTEGER
          case (hdf5DataTypeInteger8      )
             dataTypeID=H5T_INTEGER8
          case (hdf5DataTypeDouble        )
             dataTypeID=H5T_NATIVE_DOUBLE
          case (hdf5DataTypeCharacter     )
             dataTypeID=H5T_NATIVE_CHARACTER
          case (hdf5DataTypeVlenDouble    )
             dataTypeID=H5T_VLEN_DOUBLE     (1)
          case (hdf5DataTypeVlenVlenDouble)
             dataTypeID=H5T_VLEN_VLEN_DOUBLE(1)
          case (hdf5DataTypeVlenInteger8  )
             dataTypeID=H5T_VLEN_INTEGER8   (1)
          end select
       end if
       ! Create the attribute.
       call h5acreate_f(locationID,trim(attributeName),dataTypeID,dataSpaceID,self%objectID,errorCode)
       if (errorCode /= 0) then
          message="failed to create attribute '"//trim(attributeName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       ! Close the dataspace.
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close dataspace for attribute '"//trim(attributeName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
    end if

    ! Mark this object as open.
    self%isOpenValue=.true.

    ! Mark this object as a file object.
    self%hdf5ObjectType=hdf5ObjectTypeAttribute

    ! Store the name and location of the object.
    self%objectFile    =self%parentObject%objectFile
    self%objectLocation=self%parentObject%pathTo    (includeFileName=.false.)
    self%objectName    =trim(attributeName)

    ! Mark whether attribute is overwritable.
    if (present(isOverwritable)) then
       if (isOverwritable.and..not.self%parentObject%isOverwritable) then
          message="cannot make attribute '"//trim(attributeName)//"' overwritable as objects in parent '"&
               &//self%parentObject%objectName//"' are not overwritable"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       self%isOverwritable=isOverwritable
    else
       self%isOverwritable=self%parentObject%isOverwritable
    end if
    return
  end function IO_HDF5_Open_Attribute

  subroutine IO_HDF5_Write_Attribute_Logical_Scalar(self,attributeValue,attributeName)
    !!{
    Open and write a logical scalar attribute in {\normalfont \ttfamily self}.
    !!}
    implicit none
    class    (hdf5Object    ), intent(inout)           :: self
    character(len=*         ), intent(in   ), optional :: attributeName
    logical                  , intent(in   )           :: attributeValue

    if (attributeValue) then
       call self%writeAttribute(1,attributeName)
    else
       call self%writeAttribute(0,attributeName)
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Logical_Scalar

  subroutine IO_HDF5_Write_Attribute_Integer_Scalar(self,attributeValue,attributeName)
    !!{
    Open and write an integer scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_INTEGER, HSIZE_T     , h5awrite_f
    use :: ISO_Varying_String, only : assignment(=)     , operator(//), trim
    implicit none
    class    (hdf5Object    ), intent(inout)           :: self
    character(len=*         ), intent(in   ), optional :: attributeName
    integer                  , intent(in   )           :: attributeValue
    integer  (kind=HSIZE_T  ), dimension(1)            :: attributeDimensions
    integer                                            :: errorCode
    logical                                            :: preExisted
    type     (hdf5Object    )                          :: attributeObject
    type     (varying_string)                          :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a scalar integer.
          call self%assertAttributeType(H5T_NATIVE_INTEGERS,0)
       end if
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeInteger)
       ! Check that pre-existing object is a scalar integer.
       if (preExisted) call attributeObject%assertAttributeType(H5T_NATIVE_INTEGERS,0)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Write the attribute.
    call h5awrite_f(attributeObject%objectID,H5T_NATIVE_INTEGER,attributeValue,attributeDimensions,errorCode)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Integer_Scalar

  subroutine IO_HDF5_Write_Attribute_Integer_1D(self,attributeValue,attributeName)
    !!{
    Open and write an integer 1-D array attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_INTEGER, HSIZE_T     , h5awrite_f
    use :: ISO_Varying_String, only : assignment(=)     , operator(//), trim
    implicit none
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: attributeName
    integer                  , dimension(:), intent(in   )           :: attributeValue
    integer  (kind=HSIZE_T  ), dimension(1)                          :: attributeDimensions
    integer                                                          :: errorCode
    logical                                                          :: preExisted
    type     (hdf5Object    )                                        :: attributeObject
    type     (varying_string)                                        :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D integer.
          call self%assertAttributeType(H5T_NATIVE_INTEGERS,1)
       end if
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (present(attributeName)) then
          attributeNameActual=trim(attributeName)
       else
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeDimensions=shape(attributeValue)
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeInteger,attributeDimensions)
       ! Check that pre-existing object is a 1D integer.
       if (preExisted) call attributeObject%assertAttributeType(H5T_NATIVE_INTEGERS,1)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Write the attribute.
    call h5awrite_f(attributeObject%objectID,H5T_NATIVE_INTEGER,attributeValue,attributeDimensions,errorCode)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Integer_1D

  subroutine IO_HDF5_Write_Attribute_Integer8_Scalar(self,attributeValue,attributeName)
    !!{
    Open and write a long integer scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//), trim
    implicit none
    class    (hdf5Object    ), intent(inout)           :: self
    character(len=*         ), intent(in   ), optional :: attributeName
    integer  (kind=kind_int8), intent(in   ), target   :: attributeValue
    integer                                            :: errorCode
    logical                                            :: preExisted
    type     (hdf5Object    )                          :: attributeObject
    type     (varying_string)                          :: attributeNameActual, message
    type     (c_ptr         )                          :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a scalar integer.
          call self%assertAttributeType(H5T_NATIVE_INTEGER_8S,0)
       end if
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (present(attributeName)) then
          attributeNameActual=trim(attributeName)
       else
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeInteger8)
       ! Check that pre-existing object is a scalar integer.
       if (preExisted) call attributeObject%assertAttributeType(H5T_NATIVE_INTEGER_8S,0)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Write the attribute.
    dataBuffer=c_loc(attributeValue)
    errorCode=H5Awrite(attributeObject%objectID,H5T_INTEGER8,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Integer8_Scalar

  subroutine IO_HDF5_Write_Attribute_Integer8_1D(self,attributeValue,attributeName)
    !!{
    Open and write an integer 1-D array attribute in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : HSIZE_T
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//), trim
    implicit none
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   ), optional :: attributeName
    integer  (kind=kind_int8)             , dimension(:), intent(in   )           :: attributeValue
    integer  (kind=HSIZE_T  )             , dimension(1)                          :: attributeDimensions
    integer  (kind=kind_int8), allocatable, dimension(:)               , target   :: attributeValueContiguous
    integer                                                                       :: errorCode
    logical                                                                       :: preExisted
    type     (hdf5Object    )                                                     :: attributeObject
    type     (varying_string)                                                     :: attributeNameActual     , message
    type     (c_ptr         )                                                     :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D long integer.
          call self%assertAttributeType(H5T_NATIVE_INTEGER_8S,1)
       end if
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (present(attributeName)) then
          attributeNameActual=trim(attributeName)
       else
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeDimensions=shape(attributeValue)
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeInteger8,attributeDimensions)
       ! Check that pre-existing object is a 1D long integer.
       if (preExisted) call attributeObject%assertAttributeType(H5T_NATIVE_INTEGER_8S,1)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Write the attribute.
    ! We're forced to make a copy of attributeValue here because we can't pass attributeValue itself to c_loc()
    ! since it is of assumed shape.
    allocate(attributeValueContiguous,mold=attributeValue)
    attributeValueContiguous=attributeValue
    dataBuffer=c_loc(attributeValueContiguous)
    errorCode=H5Awrite(attributeObject%objectID,H5T_INTEGER8,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    deallocate(attributeValueContiguous)
    return
  end subroutine IO_HDF5_Write_Attribute_Integer8_1D

  subroutine IO_HDF5_Write_Attribute_Double_Scalar(self,attributeValue,attributeName)
    !!{
    Open and write an double scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_DOUBLE, HSIZE_T     , h5awrite_f
    use :: ISO_Varying_String, only : assignment(=)    , operator(//), trim
    implicit none
    class           (hdf5Object    ), intent(inout)           :: self
    character       (len=*         ), intent(in   ), optional :: attributeName
    double precision                , intent(in   )           :: attributeValue
    integer         (kind=HSIZE_T  ), dimension(1)            :: attributeDimensions
    integer                                                   :: errorCode
    logical                                                   :: preExisted
    type            (hdf5Object    )                          :: attributeObject
    type            (varying_string)                          :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a scalar double.
          call self%assertAttributeType(H5T_NATIVE_DOUBLES,0)
       end if
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (present(attributeName)) then
          attributeNameActual=trim(attributeName)
       else
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeDouble)
       ! Check that pre-existing object is a scalar double.
       if (preExisted) call attributeObject%assertAttributeType(H5T_NATIVE_DOUBLES,0)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Write the attribute.
    call h5awrite_f(attributeObject%objectID,H5T_NATIVE_DOUBLE,attributeValue,attributeDimensions,errorCode)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Double_Scalar

  subroutine IO_HDF5_Write_Attribute_Double_1D(self,attributeValue,attributeName)
    !!{
    Open and write an double 1-D array attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_DOUBLE, HSIZE_T     , h5awrite_f
    use :: ISO_Varying_String, only : assignment(=)    , operator(//), trim
    implicit none
    class           (hdf5Object    )              , intent(inout)           :: self
    character       (len=*         )              , intent(in   ), optional :: attributeName
    double precision                , dimension(:), intent(in   )           :: attributeValue
    integer         (kind=HSIZE_T  ), dimension(1)                          :: attributeDimensions
    integer                                                                 :: errorCode
    logical                                                                 :: preExisted
    type            (hdf5Object    )                                        :: attributeObject
    type            (varying_string)                                        :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D double.
          call self%assertAttributeType(H5T_NATIVE_DOUBLES,1)
       end if
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (present(attributeName)) then
          attributeNameActual=trim(attributeName)
       else
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeDimensions=shape(attributeValue)
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeDouble,attributeDimensions)
       ! Check that pre-existing object is a 1D double.
       if (preExisted) call attributeObject%assertAttributeType(H5T_NATIVE_DOUBLES,1)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Write the attribute.
    call h5awrite_f(attributeObject%objectID,H5T_NATIVE_DOUBLE,attributeValue,attributeDimensions,errorCode)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Double_1D

  subroutine IO_HDF5_Write_Attribute_Double_2D(self,attributeValue,attributeName)
    !!{
    Open and write an double 2-D array attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_DOUBLE, HSIZE_T     , h5awrite_f
    use :: ISO_Varying_String, only : assignment(=)    , operator(//), trim
    implicit none
    class           (hdf5Object    )                , intent(inout)           :: self
    character       (len=*         )                , intent(in   ), optional :: attributeName
    double precision                , dimension(:,:), intent(in   )           :: attributeValue
    integer         (kind=HSIZE_T  ), dimension(2)                            :: attributeDimensions
    integer                                                                   :: errorCode
    logical                                                                   :: preExisted
    type            (hdf5Object    )                                          :: attributeObject
    type            (varying_string)                                          :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 2D double.
          call self%assertAttributeType(H5T_NATIVE_DOUBLES,2)
       end if
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (present(attributeName)) then
          attributeNameActual=trim(attributeName)
       else
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeDimensions=shape(attributeValue)
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeDouble,attributeDimensions)
       ! Check that pre-existing object is a 2D double.
       if (preExisted) call attributeObject%assertAttributeType(H5T_NATIVE_DOUBLES,2)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Write the attribute.
    call h5awrite_f(attributeObject%objectID,H5T_NATIVE_DOUBLE,attributeValue,attributeDimensions,errorCode)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Double_2D

  subroutine IO_HDF5_Write_Attribute_Character_Scalar(self,attributeValue,attributeName)
    !!{
    Open and write an character scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_CHARACTER, HID_T       , HSIZE_T      , h5awrite_f, &
          &                           h5tclose_f          , h5tcopy_f   , h5tset_size_f, size_t
    use :: ISO_Varying_String, only : assignment(=)       , operator(//), trim
    implicit none
    class    (hdf5Object    ), intent(inout)           :: self
    character(len=*         ), intent(in   ), optional :: attributeName
    character(len=*         ), intent(in   )           :: attributeValue
    integer  (kind=HSIZE_T  ), dimension(1)            :: attributeDimensions
    integer  (kind=HID_T    )                          :: dataTypeID
    integer                                            :: errorCode
    logical                                            :: preExisted
    type     (hdf5Object    )                          :: attributeObject
    type     (varying_string)                          :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Create a custom datatype.
    call h5tcopy_f(H5T_NATIVE_CHARACTER,dataTypeID,errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype for attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tset_size_f(dataTypeID,int(len(attributeValue),size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size for attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a scalar character.
          call self%assertAttributeType([dataTypeID],0)
       end if
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (present(attributeName)) then
          attributeNameActual=trim(attributeName)
       else
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeCharacter,useDataType=dataTypeID)
       ! Check that pre-existing object is a scalar character.
       if (preExisted) call attributeObject%assertAttributeType([dataTypeID],0)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Write the attribute.
    call h5awrite_f(attributeObject%objectID,dataTypeID,attributeValue,attributeDimensions,errorCode)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID,errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Character_Scalar

  subroutine IO_HDF5_Write_Attribute_Character_1D(self,attributeValue,attributeName)
    !!{
    Open and write an character 1-D array attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_CHARACTER, HID_T        , HSIZE_T, h5awrite_f, &
          &                           h5tcopy_f           , h5tset_size_f, size_t
    use :: ISO_Varying_String, only : assignment(=)       , operator(//) , trim
    implicit none
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: attributeName
    character(len=*         ), dimension(:), intent(in   )           :: attributeValue
    integer  (kind=HSIZE_T  ), dimension(1)                          :: attributeDimensions
    integer  (kind=HID_T    )                                        :: dataTypeID
    integer                                                          :: errorCode
    logical                                                          :: preExisted
    type     (hdf5Object    )                                        :: attributeObject
    type     (varying_string)                                        :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Create a custom datatype.
    call h5tcopy_f(H5T_NATIVE_CHARACTER,dataTypeID,errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype for attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tset_size_f(dataTypeID,int(len(attributeValue),size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size for attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! If this attribute if not overwritable, report an error.
       if (.not.self%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D character.
          call self%assertAttributeType([dataTypeID],1)
       end if
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       attributeNameActual=self%objectName
    else
       ! Check that an attribute name was supplied.
       if (present(attributeName)) then
          attributeNameActual=trim(attributeName)
       else
          message="no name was supplied for attribute in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if attribute already exists.
       preExisted=self%hasAttribute(attributeName)
       ! Open the attribute.
       attributeDimensions=shape(attributeValue)
       attributeObject=self%openAttribute(attributeName,hdf5DataTypeCharacter,attributeDimensions,useDataType=dataTypeID)
       ! Check that pre-existing object is a 1D character.
       if (preExisted) call attributeObject%assertAttributeType([dataTypeID],1)
       ! If this attribute if not overwritable, report an error.
       if (preExisted.and..not.attributeObject%isOverwritable) then
          message="attribute '"//trim(attributeNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Write the attribute.
    call h5awrite_f(attributeObject%objectID,dataTypeID,attributeValue,attributeDimensions,errorCode)
    if (errorCode /= 0) then
       message="unable to write attribute '"//attributeNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Write_Attribute_Character_1D

  subroutine IO_HDF5_Write_Attribute_VarString_Scalar(self,attributeValue,attributeName)
    !!{
    Open and write a varying string scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    class    (hdf5Object    ), intent(inout)           :: self
    character(len=*         ), intent(in   ), optional :: attributeName
    type     (varying_string), intent(in   )           :: attributeValue

    ! Call the character version of this routine to perform the write.
    call IO_HDF5_Write_Attribute_Character_Scalar(self,char(attributeValue),attributeName)
    return
  end subroutine IO_HDF5_Write_Attribute_VarString_Scalar

  subroutine IO_HDF5_Write_Attribute_VarString_1D(self,attributeValue,attributeName)
    !!{
    Open and write a varying string 1-D array attribute in {\normalfont \ttfamily self}.
    !!}
    use :: String_Handling, only : Convert_VarString_To_Char
    implicit none
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: attributeName
    type     (varying_string), dimension(:), intent(in   )           :: attributeValue

    ! Call the character version of this routine to perform the write.
    call IO_HDF5_Write_Attribute_Character_1D(self,Convert_VarString_To_Char(attributeValue),attributeName)

    return
  end subroutine IO_HDF5_Write_Attribute_VarString_1D

  subroutine IO_HDF5_Read_Attribute_Integer_Scalar(self,attributeName,attributeValue,allowPseudoScalar)
    !!{
    Open and read an integer scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_INTEGER, HID_T       , HSIZE_T                    , h5aget_space_f, &
          &                           h5aread_f         , h5sclose_f  , h5sget_simple_extent_dims_f
    use :: ISO_Varying_String, only : assignment(=)     , operator(//), trim
    implicit none
    integer                                , intent(  out)           :: attributeValue
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: attributeName
    logical                                , intent(in   ), optional :: allowPseudoScalar
    integer                  , dimension(1)                          :: pseudoScalarValue
    integer  (kind=HSIZE_T  ), dimension(1)                          :: attributeDimensions    , attributeMaximumDimensions
    integer  (kind=HID_T    )                                        :: attributeDataspaceID
    integer                                                          :: errorCode
    type     (hdf5Object    )                                        :: attributeObject
    type     (varying_string)                                        :: attributeNameActual    , message
    logical                                                          :: allowPseudoScalarActual, matches

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check if pseudo-scalars are allowed.
    if (present(allowPseudoScalar)) then
       allowPseudoScalarActual=allowPseudoScalar
    else
       allowPseudoScalarActual=.false.
    end if

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a scalar integer.
    call attributeObject%assertAttributeType(H5T_NATIVE_INTEGERS,0,matches)
    if (matches) then
       ! Read the scalar attribute.
       call h5aread_f(attributeObject%objectID,H5T_NATIVE_INTEGER,attributeValue,attributeDimensions&
            &,errorCode)
       if (errorCode /= 0) then
          message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else if (allowPseudoScalarActual) then
       ! Attribute is not a scalar. Check if it is a pseudo-scalar.
       call attributeObject%assertAttributeType(H5T_NATIVE_INTEGERS,1,matches)
       if (matches) then
          ! Get the dimensions of the array.
          call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
          if (errorCode /= 0) then
             message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
          if (errorCode < 0) then
             message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          call h5sclose_f(attributeDataspaceID,errorCode)
          if (errorCode /= 0) then
             message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          if (attributeDimensions(1) == 1) then
             call attributeObject%readAttributeStatic(attributeValue=pseudoScalarValue)
             attributeValue=pseudoScalarValue(1)
          else
             call Error_Report("attribute '"//attributeObject%objectName//"' must be an integer scalar or pseudo-scalar"//self%locationReport()//self%locationReport()//{introspection:location})
          end if
       end if
    else
       call       Error_Report("attribute '"//attributeObject%objectName//"' must be an integer scalar"                 //self%locationReport()//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Integer_Scalar

  subroutine IO_HDF5_Read_Attribute_Integer_1D_Array_Allocatable(self,attributeName,attributeValue)
    !!{
    Open and read an integer scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_INTEGER, HID_T       , HSIZE_T                    , h5aget_space_f, &
          &                           h5aread_f         , h5sclose_f  , h5sget_simple_extent_dims_f
    use :: ISO_Varying_String, only : assignment(=)     , operator(//), trim
    implicit none
    integer                  , allocatable, dimension(:), intent(  out)           :: attributeValue
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   ), optional :: attributeName
    integer  (kind=HSIZE_T  )             , dimension(1)                          :: attributeDimensions , attributeMaximumDimensions
    integer                                                                       :: errorCode
    integer  (kind=HID_T    )                                                     :: attributeDataspaceID
    type     (hdf5Object    )                                                     :: attributeObject
    type     (varying_string)                                                     :: attributeNameActual , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a 1D integer array.
    call attributeObject%assertAttributeType(H5T_NATIVE_INTEGERS,1)

    ! Get the dimensions of the array.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(attributeValue)) deallocate(attributeValue)
    !![
    <allocate variable="attributeValue" shape="attributeDimensions"/>
    !!]

    ! Read the attribute.
    call h5aread_f(attributeObject%objectID,H5T_NATIVE_INTEGER,attributeValue,attributeDimensions&
         &,errorCode)
    if (errorCode /= 0) then
       message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Integer_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Attribute_Integer_1D_Array_Static(self,attributeName,attributeValue)
    !!{
    Open and read an integer scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_INTEGER, HID_T       , HSIZE_T                    , h5aget_space_f, &
          &                           h5aread_f         , h5sclose_f  , h5sget_simple_extent_dims_f
    use :: ISO_Varying_String, only : assignment(=)     , operator(//), trim
    implicit none
    integer                  , dimension(:), intent(  out)           :: attributeValue
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: attributeName
    integer  (kind=HSIZE_T  ), dimension(1)                          :: attributeDimensions , attributeMaximumDimensions
    integer                                                          :: errorCode
    integer  (kind=HID_T    )                                        :: attributeDataspaceID
    type     (hdf5Object    )                                        :: attributeObject
    type     (varying_string)                                        :: attributeNameActual , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a 1D integer array.
    call attributeObject%assertAttributeType(H5T_NATIVE_INTEGERS,1)

    ! Get the dimensions of the array.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Ensure that the size of the array is large enough to hold the attributes.
    if (any(shape(attributeValue) < attributeDimensions)) then
       message="array is not large enough to hold attributes from '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the attribute.
    call h5aread_f(attributeObject%objectID,H5T_NATIVE_INTEGER,attributeValue,attributeDimensions&
         &,errorCode)
    if (errorCode /= 0) then
       message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Integer_1D_Array_Static

  subroutine IO_HDF5_Read_Attribute_Integer8_Scalar(self,attributeName,attributeValue,allowPseudoScalar)
    !!{
    Open and read a long integer scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : h5sget_simple_extent_dims_f, HID_T       , HSIZE_T, h5aget_space_f, &
          &                                      h5sclose_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)              , operator(//), trim
    implicit none
    integer  (kind=kind_int8)              , intent(  out)          , target :: attributeValue
    class    (hdf5Object    )              , intent(inout)                   :: self
    character(len=*         )              , intent(in   ), optional         :: attributeName
    logical                                , intent(in   ), optional         :: allowPseudoScalar
    integer  (kind=kind_int8), dimension(1)                                  :: pseudoScalarValue
    integer  (kind=HSIZE_T  ), dimension(1)                                  :: attributeDimensions    , attributeMaximumDimensions
    integer  (kind=HID_T    )                                                :: attributeDataspaceID
    integer                                                                  :: errorCode
    type     (hdf5Object    )                                                :: attributeObject
    type     (varying_string)                                                :: attributeNameActual    , message
    logical                                                                  :: allowPseudoScalarActual, matches
    type     (c_ptr         )                                                :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check if pseudo-scalars are allowed.
    if (present(allowPseudoScalar)) then
       allowPseudoScalarActual=allowPseudoScalar
    else
       allowPseudoScalarActual=.false.
    end if

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a scalar integer.
    call attributeObject%assertAttributeType(H5T_NATIVE_INTEGER_8AS,0,matches)
    if (matches) then
       ! Read the attribute.
       dataBuffer=c_loc(attributeValue)
       errorCode=H5Aread(attributeObject%objectID,H5T_INTEGER8,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else if (allowPseudoScalarActual) then
       ! Attribute is not a scalar. Check if it is a pseudo-scalar.
       call attributeObject%assertAttributeType(H5T_NATIVE_INTEGER_8AS,1,matches)
       if (matches) then
          ! Get the dimensions of the array.
          call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
          if (errorCode /= 0) then
             message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
          if (errorCode < 0) then
             message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          call h5sclose_f(attributeDataspaceID,errorCode)
          if (errorCode /= 0) then
             message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          if (attributeDimensions(1) == 1) then
             call attributeObject%readAttributeStatic(attributeValue=pseudoScalarValue)
             attributeValue=pseudoScalarValue(1)
          else
             call Error_Report("attribute '"//attributeObject%objectName//"' must be a long integer scalar or pseudo-scalar"//self%locationReport()//{introspection:location})
          end if
       end if
    else
       call       Error_Report("attribute '"//attributeObject%objectName//"' must be a long integer scalar"                 //self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Integer8_Scalar

  subroutine IO_HDF5_Read_Attribute_Integer8_1D_Array_Allocatable(self,attributeName,attributeValue)
    !!{
    Open and read an integer scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : h5sget_simple_extent_dims_f, HID_T      , HSIZE_T, h5aget_space_f, &
          &                                      h5sclose_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)             , operator(//), trim
    implicit none
    integer  (kind=kind_int8), allocatable, dimension(:), intent(  out), target   :: attributeValue
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   ), optional :: attributeName
    integer  (kind=HSIZE_T  )             , dimension(1)                          :: attributeDimensions , attributeMaximumDimensions
    integer                                                                       :: errorCode
    integer  (kind=HID_T    )                                                     :: attributeDataspaceID
    type     (hdf5Object    )                                                     :: attributeObject
    type     (varying_string)                                                     :: attributeNameActual , message
    type     (c_ptr         )                                                     :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a 1D long integer array.
    call attributeObject%assertAttributeType(H5T_NATIVE_INTEGER_8AS,1)

    ! Get the dimensions of the array.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(attributeValue)) deallocate(attributeValue)
    !![
    <allocate variable="attributeValue" shape="attributeDimensions"/>
    !!]

    ! Read the attribute.
    dataBuffer=c_loc(attributeValue)
    errorCode=H5Aread(attributeObject%objectID,H5T_INTEGER8,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Integer8_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Attribute_Integer8_1D_Array_Static(self,attributeName,attributeValue)
    !!{
    Open and read an integer scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : h5sget_simple_extent_dims_f, HID_T       , HSIZE_T, h5aget_space_f, &
          &                                      h5sclose_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)              , operator(//), trim
    implicit none
    integer  (kind=kind_int8)             , dimension(:), intent(  out)           :: attributeValue
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   ), optional :: attributeName
    integer  (kind=HSIZE_T  )             , dimension(1)                          :: attributeDimensions     , attributeMaximumDimensions
    integer  (kind=kind_int8), allocatable, dimension(:)               , target   :: attributeValueContiguous
    integer                                                                       :: errorCode
    integer  (kind=HID_T    )                                                     :: attributeDataspaceID
    type     (hdf5Object    )                                                     :: attributeObject
    type     (varying_string)                                                     :: attributeNameActual     , message
    type     (c_ptr         )                                                     :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a 1D long integer array.
    call attributeObject%assertAttributeType(H5T_NATIVE_INTEGER_8AS,1)

    ! Get the dimensions of the array.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Ensure that the size of the array is large enough to hold the attributes.
    if (any(shape(attributeValue) < attributeDimensions)) then
       message="array is not large enough to hold attributes from '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the attribute.
    ! We're forced to make a copy of attributeValue here because we can't pass attributeValue itself to c_loc()
    ! since it is of assumed shape.
    allocate(attributeValueContiguous,mold=attributeValue)
    dataBuffer=c_loc(attributeValueContiguous)
    errorCode=H5Aread(attributeObject%objectID,H5T_INTEGER8,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    attributeValue=attributeValueContiguous    
    return
  end subroutine IO_HDF5_Read_Attribute_Integer8_1D_Array_Static

  subroutine IO_HDF5_Read_Attribute_Double_Scalar(self,attributeName,attributeValue,allowPseudoScalar)
    !!{
    Open and read an double scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_DOUBLE, HID_T       , HSIZE_T                    , h5aget_space_f, &
          &                           h5aread_f        , h5sclose_f  , h5sget_simple_extent_dims_f
    use :: ISO_Varying_String, only : assignment(=)    , operator(//), trim
    implicit none
    double precision                              , intent(  out)           :: attributeValue
    class           (hdf5Object    )              , intent(inout)           :: self
    character       (len=*         )              , intent(in   ), optional :: attributeName
    logical                                       , intent(in   ), optional :: allowPseudoScalar
    integer         (kind=HSIZE_T  ), dimension(1)                          :: attributeDimensions    , attributeMaximumDimensions
    double precision                , dimension(1)                          :: pseudoScalarValue
    integer         (kind=HID_T    )                                        :: attributeDataspaceID
    integer                                                                 :: errorCode
    type            (hdf5Object    )                                        :: attributeObject
    type            (varying_string)                                        :: attributeNameActual    , message
    logical                                                                 :: allowPseudoScalarActual, matches

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check if pseudo-scalars are allowed.
    if (present(allowPseudoScalar)) then
       allowPseudoScalarActual=allowPseudoScalar
    else
       allowPseudoScalarActual=.false.
    end if

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a scalar double.
    call attributeObject%assertAttributeType(H5T_NATIVE_DOUBLES,0,matches)
    if (matches) then
       ! Read the attribute.
       call h5aread_f(attributeObject%objectID,H5T_NATIVE_DOUBLE,attributeValue,attributeDimensions&
            &,errorCode)
       if (errorCode /= 0) then
          message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else if (allowPseudoScalarActual) then
       ! Attribute is not a scalar. Check if it is a pseudo-scalar.
       call attributeObject%assertAttributeType(H5T_NATIVE_DOUBLES,1,matches)
       if (matches) then
          ! Get the dimensions of the array.
          call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
          if (errorCode /= 0) then
             message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
          if (errorCode < 0) then
             message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          call h5sclose_f(attributeDataspaceID,errorCode)
          if (errorCode /= 0) then
             message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          if (attributeDimensions(1) == 1) then
             call attributeObject%readAttributeStatic(attributeValue=pseudoScalarValue)
             attributeValue=pseudoScalarValue(1)
          else
             call Error_Report("attribute '"//attributeObject%objectName//"' must be a double scalar or pseudo-scalar"//self%locationReport()//{introspection:location})
          end if
       else
          call    Error_Report("attribute '"//attributeObject%objectName//"' must be a double scalar or pseudo-scalar"//self%locationReport()//{introspection:location})
       end if
    else
       call       Error_Report("attribute '"//attributeObject%objectName//"' must be a double scalar"                 //self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Double_Scalar

  subroutine IO_HDF5_Read_Attribute_Double_1D_Array_Allocatable(self,attributeName,attributeValue)
    !!{
    Open and read an double scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_DOUBLE, HID_T       , HSIZE_T                    , h5aget_space_f, &
          &                           h5aread_f        , h5sclose_f  , h5sget_simple_extent_dims_f
    use :: ISO_Varying_String, only : assignment(=)    , operator(//), trim
    implicit none
    double precision                , allocatable, dimension(:), intent(  out)           :: attributeValue
    class           (hdf5Object    )                           , intent(inout)           :: self
    character       (len=*         )                           , intent(in   ), optional :: attributeName
    integer         (kind=HSIZE_T  )             , dimension(1)                          :: attributeDimensions , attributeMaximumDimensions
    integer                                                                              :: errorCode
    integer         (kind=HID_T    )                                                     :: attributeDataspaceID
    type            (hdf5Object    )                                                     :: attributeObject
    type            (varying_string)                                                     :: attributeNameActual , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a 1D double array.
    call attributeObject%assertAttributeType(H5T_NATIVE_DOUBLES,1)

    ! Get the dimensions of the array.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(attributeValue)) deallocate(attributeValue)
    !![
    <allocate variable="attributeValue" shape="attributeDimensions"/>
    !!]

    ! Read the attribute.
    call h5aread_f(attributeObject%objectID,H5T_NATIVE_DOUBLE,attributeValue,attributeDimensions&
         &,errorCode)
    if (errorCode /= 0) then
       message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Double_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Attribute_Double_1D_Array_Static(self,attributeName,attributeValue)
    !!{
    Open and read an double scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_NATIVE_DOUBLE, HID_T       , HSIZE_T                    , h5aget_space_f, &
          &                           h5aread_f        , h5sclose_f  , h5sget_simple_extent_dims_f
    use :: ISO_Varying_String, only : assignment(=)    , operator(//), trim
    use :: ISO_Varying_String, only : char
    implicit none
    double precision                , dimension(:), intent(  out)           :: attributeValue
    class           (hdf5Object    )              , intent(inout)           :: self
    character       (len=*         )              , intent(in   ), optional :: attributeName
    integer         (kind=HSIZE_T  ), dimension(1)                          :: attributeDimensions , attributeMaximumDimensions
    integer                                                                 :: errorCode
    integer         (kind=HID_T    )                                        :: attributeDataspaceID
    type            (hdf5Object    )                                        :: attributeObject
    type            (varying_string)                                        :: attributeNameActual , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a 1D double array.
    call attributeObject%assertAttributeType(H5T_NATIVE_DOUBLES,1)

    ! Get the dimensions of the array.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Ensure that the size of the array is large enough to hold the attributes.
    if (any(shape(attributeValue) < attributeDimensions)) then
       message="array is not large enough to hold attributes from '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the attribute.
    call h5aread_f(attributeObject%objectID,H5T_NATIVE_DOUBLE,attributeValue,attributeDimensions&
         &,errorCode)
    if (errorCode /= 0) then
       message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Double_1D_Array_Static

  subroutine IO_HDF5_Read_Attribute_Character_Scalar(self,attributeName,attributeValue,allowPseudoScalar)
    !!{
    Open and read an character scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_loc, c_ptr, c_null_char, c_f_pointer
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : HID_T              , HSIZE_T                    , h5aget_space_f , h5aread_f , &
         &                                       h5sclose_f         , h5sget_simple_extent_dims_f, h5tclose_f     , h5tcopy_f , &
         &                                       h5tset_cset_f      , h5t_cset_utf8_f            , h5tset_size_f  , h5t_string, &
         &                                       H5T_Str_NullTerm_f , H5T_C_S1                   , h5tset_strpad_f
    use            :: ISO_Varying_String, only : assignment(=)      , operator(//)               , trim

    use            :: String_Handling   , only : String_C_to_Fortran
    implicit none
    character(len=*                  )              , intent(  out)           :: attributeValue
    class    (hdf5Object             )              , intent(inout)           :: self
    character(len=*                  )              , intent(in   ), optional :: attributeName
    logical                                         , intent(in   ), optional :: allowPseudoScalar
    integer  (kind=HSIZE_T           ), dimension(1)                          :: attributeDimensions       , attributeMaximumDimensions
    character(len=len(attributeValue)), dimension(1)                          :: pseudoScalarValue
    integer  (kind=HID_T             )                                        :: attributeDataspaceID      , stringType
    integer  (kind=HID_T             )                                        :: dataTypeID             (6)
    integer                                                                   :: errorCode
    type     (hdf5Object             )                                        :: attributeObject
    type     (varying_string         )                                        :: attributeNameActual       , message
    logical                                                                   :: allowPseudoScalarActual   , matches                   , &
         &                                                                       isH5TString
    type     (c_ptr                  )              , target                  :: stringBuffer
    character(c_char                 ), dimension(:), pointer                 :: stringBuffer_   
    
    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check if pseudo-scalars are allowed.
    if (present(allowPseudoScalar)) then
       allowPseudoScalarActual=allowPseudoScalar
    else
       allowPseudoScalarActual=.false.
    end if

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Create a custom datatype. We actually make two types - one is a Fortran native type, the other is a C native type.
    dataTypeID=IO_HDF5_Character_Types(len(attributeValue))

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a scalar character.
    call attributeObject%assertAttributeType (dataTypeID ,0,matches    )
    call attributeObject%assertAttributeType([H5T_String],0,isH5TString)
    if (matches) then
       ! Read the attribute.
       call h5aread_f(attributeObject%objectID,dataTypeID(1),attributeValue,attributeDimensions&
            &,errorCode)
       if (errorCode /= 0) then
          message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else if (isH5TString) then
       ! Attribute is an H5T_String (i.e. a variable length string, as is written by h5py for example). These can be read using
       ! the special datatype size "H5T_Variable" (which is available only in the C API, so we get it via a C function).
       !
       ! First, construct a suitable datatype.
       call H5Tcopy_f(H5T_C_S1, stringType, errorCode)
       if (errorCode /= 0) call Error_Report('unable to copy datatype'       //self%locationReport()//{introspection:location})
       call H5Tset_cset_f  (stringType,H5T_CSET_UTF8_F   ,errorCode)
       if (errorCode /= 0) call Error_Report('unable to set character set'   //self%locationReport()//{introspection:location})
       call H5Tset_size_f  (stringType,h5t_variable_get(),errorCode)            
       if (errorCode /= 0) call Error_Report('unable to set datatype size'   //self%locationReport()//{introspection:location})
       call h5tset_strpad_f(stringType,H5T_STR_NULLTERM_F,errorCode)
       if (errorCode /= 0) call Error_Report('unable to set datatype padding'//self%locationReport()//{introspection:location})
       ! Read the attribute.
       errorCode=H5Aread(attributeObject%objectID,stringType,c_loc(stringBuffer))
       if (errorCode /= 0) then
          message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Extract the attribute from the buffer.
       call c_f_pointer(stringBuffer,stringBuffer_,shape=[len(attributeValue)])
       attributeValue=String_C_to_Fortran(stringBuffer_)
       deallocate(stringBuffer_)
    else if (allowPseudoScalarActual) then
       ! Attribute is not a scalar. Check if it is a pseudo-scalar.
       call attributeObject%assertAttributeType(dataTypeID,1,matches)
       if (matches) then
          ! Get the dimensions of the array.
          call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
          if (errorCode /= 0) then
             message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
          if (errorCode < 0) then
             message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          call h5sclose_f(attributeDataspaceID,errorCode)
          if (errorCode /= 0) then
             message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          if (attributeDimensions(1) == 1) then
             call attributeObject%readAttributeStatic(attributeValue=pseudoScalarValue)
             attributeValue=pseudoScalarValue(1)
          else
             call Error_Report("attribute must be a character scalar, pseudo-scalar, or variable-length string"//self%locationReport()//{introspection:location})
          end if
       end if
    else
       call       Error_Report("attribute must be a character scalar, or variable-length string"               //self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID(1),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tclose_f(dataTypeID(2),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Character_Scalar

  subroutine IO_HDF5_Read_Attribute_Character_1D_Array_Allocatable(self,attributeName,attributeValue)
    !!{
    Open and read an character scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , HSIZE_T                    , h5aget_space_f, h5aread_f, &
          &                           h5sclose_f   , h5sget_simple_extent_dims_f, h5tclose_f
    use :: ISO_Varying_String, only : assignment(=), operator(//)               , trim
    implicit none
    character(len=*         ), allocatable, dimension(:), intent(  out)           :: attributeValue
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   ), optional :: attributeName
    integer  (kind=HSIZE_T  )             , dimension(1)                          :: attributeDimensions , attributeMaximumDimensions
    integer                                                                       :: errorCode
    integer  (kind=HID_T    )                                                     :: attributeDataspaceID, dataTypeID                (6)
    type     (hdf5Object    )                                                     :: attributeObject
    type     (varying_string)                                                     :: attributeNameActual , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Create a custom datatype.
    dataTypeID=IO_HDF5_Character_Types(len(attributeValue))

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a 1D character array.
    call attributeObject%assertAttributeType(dataTypeID,1)

    ! Get the dimensions of the array.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(attributeValue)) deallocate(attributeValue)
    !![
    <allocate variable="attributeValue" shape="attributeDimensions"/>
    !!]

    ! Read the attribute.
    call h5aread_f(attributeObject%objectID,dataTypeID(1),attributeValue,attributeDimensions&
         &,errorCode)
    if (errorCode /= 0) then
       message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID(1),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tclose_f(dataTypeID(2),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Character_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Attribute_Character_1D_Array_Static(self,attributeName,attributeValue)
    !!{
    Open and read an character scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , HSIZE_T                    , h5aget_space_f, h5aread_f, &
          &                           h5sclose_f   , h5sget_simple_extent_dims_f, h5tclose_f
    use :: ISO_Varying_String, only : assignment(=), operator(//)               , trim
    implicit none
    character(len=*         ), dimension(:), intent(  out)           :: attributeValue
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: attributeName
    integer  (kind=HSIZE_T  ), dimension(1)                          :: attributeDimensions , attributeMaximumDimensions
    integer                                                          :: errorCode
    integer  (kind=HID_T    )                                        :: attributeDataspaceID, dataTypeID                (6)
    type     (hdf5Object    )                                        :: attributeObject
    type     (varying_string)                                        :: attributeNameActual , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Create a custom datatype.
    dataTypeID=IO_HDF5_Character_Types(len(attributeValue))

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
          attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Check that the object is a 1D character array.
    call attributeObject%assertAttributeType(dataTypeID,1)

    ! Get the dimensions of the array.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_dims_f(attributeDataspaceID,attributeDimensions,attributeMaximumDimensions,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Ensure that the size of the array is large enough to hold the attributes.
    if (any(shape(attributeValue) < attributeDimensions)) then
       message="array is not large enough to hold attributes from '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the attribute.
    call h5aread_f(attributeObject%objectID,dataTypeID(1),attributeValue,attributeDimensions&
         &,errorCode)
    if (errorCode /= 0) then
       message="unable to read attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID(1),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tclose_f(dataTypeID(2),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(attributeNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Attribute_Character_1D_Array_Static

  subroutine IO_HDF5_Read_Attribute_VarString_Scalar(self,attributeName,attributeValue,allowPseudoScalar)
    !!{
    Open and read an varying string scalar attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , h5aget_type_f, h5tclose_f, h5tget_size_f, &
         &                            H5T_String
    use :: ISO_Varying_String, only : assignment(=), operator(//) , trim      , len_trim
    implicit none
    type     (varying_string), intent(  out)           :: attributeValue
    class    (hdf5Object    ), intent(inout)           :: self
    character(len=*         ), intent(in   ), optional :: attributeName
    logical                  , intent(in   ), optional :: allowPseudoScalar
    integer  (kind=SIZE_T   )                          :: dataTypeSizeMaximum=65536
    integer  (kind=HID_T    )                          :: dataTypeID
    integer  (kind=SIZE_T   )                          :: dataTypeSize             , lengthPrevious
    integer                                            :: errorCode
    logical                                            :: isH5TString
    type     (hdf5Object    )                          :: attributeObject
    type     (varying_string)                          :: attributeNameActual      , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
   if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Get the datatype of this attribute.
    call h5aget_type_f(attributeObject%objectID,dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not get datatype for '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Get the size of the datatype.
    call h5tget_size_f(dataTypeID,dataTypeSize,errorCode)
    if (errorCode /= 0) then
       message="can not get size of datatype for '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not close datatype of '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check for a variable length string.
    call attributeObject%assertAttributeType([H5T_String],0,isH5TString)

    ! Call wrapper routine that will do the remainder of the read.
    if (isH5TString) then
       lengthPrevious=-1
       dataTypeSize  = 1
       do while (dataTypeSize < dataTypeSizeMaximum)
          dataTypeSize=dataTypeSize*2
          call IO_HDF5_Read_Attribute_VarString_Scalar_Do_Read(self,attributeName,attributeValue,dataTypeSize,allowPseudoScalar)
          if (len_trim(attributeValue) == lengthPrevious) then
             exit
          else
             lengthPrevious=len_trim(attributeValue)
          end if
       end do
       if (len_trim(attributeValue) /= lengthPrevious) call Error_Report('variable length HDF5 string is too long'//self%locationReport()//{introspection:location})
attributeValue=trim(attributeValue)
    else
       call IO_HDF5_Read_Attribute_VarString_Scalar_Do_Read(self,attributeName,attributeValue,dataTypeSize,allowPseudoScalar)
    end if

    return
  end subroutine IO_HDF5_Read_Attribute_VarString_Scalar

  subroutine IO_HDF5_Read_Attribute_VarString_Scalar_Do_Read(self,attributeName,attributeValue,dataTypeSize,allowPseudoScalar)
    !!{
    Open and read an varying string scalar attribute in {\normalfont \ttfamily self} by creating a suitably-sized character variable into
    which it can be read.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type     (varying_string  ), intent(  out)           :: attributeValue
    class    (hdf5Object      ), intent(inout)           :: self
    character(len=*           ), intent(in   ), optional :: attributeName
    logical                    , intent(in   ), optional :: allowPseudoScalar
    integer  (kind=SIZE_T     ), intent(in   )           :: dataTypeSize
    character(len=dataTypeSize)                          :: temporaryBuffer

    ! Call the character version of this routine to perform the red.
    call IO_HDF5_Read_Attribute_Character_Scalar(self,attributeName,temporaryBuffer,allowPseudoScalar)

    ! Transfer the results to the varying string variable.
    attributeValue=temporaryBuffer

    return
  end subroutine IO_HDF5_Read_Attribute_VarString_Scalar_Do_Read

  subroutine IO_HDF5_Read_Attribute_VarString_1D_Array_Allocatable(self,attributeName,attributeValue)
    !!{
    Open and read an varying string 1-D array attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , h5aget_type_f, h5tclose_f, h5tget_size_f
    use :: ISO_Varying_String, only : assignment(=), operator(//) , trim
    implicit none
    type     (varying_string), allocatable, dimension(:), intent(  out)           :: attributeValue
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   ), optional :: attributeName
    integer  (kind=HID_T    )                                                     :: dataTypeID
    integer  (kind=SIZE_T   )                                                     :: dataTypeSize
    integer                                                                       :: errorCode
    type     (hdf5Object    )                                                     :: attributeObject
    type     (varying_string)                                                     :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Get the datatype of this attribute.
    call h5aget_type_f(attributeObject%objectID,dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not get datatype for '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Get the size of the datatype.
    call h5tget_size_f(dataTypeID,dataTypeSize,errorCode)
    if (errorCode /= 0) then
       message="can not get size of datatype for '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not close datatype of '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Call wrapper routine that will do the remainder of the read.
    call IO_HDF5_Read_Attribute_VarString_1D_Array_Allocatable_Do_Read(self,attributeName,attributeValue,dataTypeSize)
    return
  end subroutine IO_HDF5_Read_Attribute_VarString_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Attribute_VarString_1D_Array_Allocatable_Do_Read(self,attributeName,attributeValue,dataTypeSize)
    !!{
    Open and read an varying string 1-D array attribute in {\normalfont \ttfamily self} by creating a suitably-sized character variable into
    which it can be read.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type     (varying_string  ), allocatable, dimension(:), intent(  out)           :: attributeValue
    class    (hdf5Object      )                           , intent(inout)           :: self
    character(len=*           )                           , intent(in   ), optional :: attributeName
    integer  (kind=SIZE_T     )                           , intent(in   )           :: dataTypeSize
    character(len=dataTypeSize), allocatable, dimension(:)                          :: temporaryBuffer

    ! Call the character version of this routine to perform the red.
    call IO_HDF5_Read_Attribute_Character_1D_Array_Allocatable(self,attributeName,temporaryBuffer)

    ! Transfer the results to the varying string variable.
    allocate(attributeValue(size(temporaryBuffer)))
    attributeValue=temporaryBuffer
    deallocate(temporaryBuffer)

    return
  end subroutine IO_HDF5_Read_Attribute_VarString_1D_Array_Allocatable_Do_Read

  subroutine IO_HDF5_Read_Attribute_VarString_1D_Array_Static(self,attributeName,attributeValue)
    !!{
    Open and read an varying string 1-D array attribute in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , h5aget_type_f, h5tclose_f, h5tget_size_f
    use :: ISO_Varying_String, only : assignment(=), operator(//) , trim
    implicit none
    type     (varying_string), dimension(:), intent(inout)           :: attributeValue
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: attributeName
    integer  (kind=HID_T    )                                        :: dataTypeID
    integer  (kind=SIZE_T   )                                        :: dataTypeSize
    integer                                                          :: errorCode
    type     (hdf5Object    )                                        :: attributeObject
    type     (varying_string)                                        :: attributeNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the attribute.
    if (present(attributeName)) then
       attributeNameActual=attributeName
    else
       attributeNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read attribute '"//trim(attributeNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an attribute, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeAttribute) then
       ! Object is the attribute.
       select type (self)
       type is (hdf5Object)
       attributeObject=self
       end select
       ! No name should be supplied in this case.
       if (present(attributeName)) then
          message="attribute name was supplied for attribute object '"//trim(attributeName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an attribute name was supplied.
       if (.not.present(attributeName)) then
          message="attribute name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the attribute exists.
       if (.not.self%hasAttribute(attributeName)) then
          message="attribute '"//trim(attributeName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the attribute.
       attributeObject=self%openAttribute(attributeName)
    end if

    ! Get the datatype of this attribute.
    call h5aget_type_f(attributeObject%objectID,dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not get datatype for '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Get the size of the datatype.
    call h5tget_size_f(dataTypeID,dataTypeSize,errorCode)
    if (errorCode /= 0) then
       message="can not get size of datatype for '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not close datatype of '"//trim(attributeNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Call wrapper routine that will do the remainder of the read.
    call IO_HDF5_Read_Attribute_VarString_1D_Array_Static_Do_Read(self,attributeName,attributeValue,dataTypeSize)
    return
  end subroutine IO_HDF5_Read_Attribute_VarString_1D_Array_Static

  subroutine IO_HDF5_Read_Attribute_VarString_1D_Array_Static_Do_Read(self,attributeName,attributeValue,dataTypeSize)
    !!{
    Open and read an varying string 1-D array attribute in {\normalfont \ttfamily self} by creating a suitably-sized character variable into
    which it can be read.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type     (varying_string  ), dimension(:)                   , intent(inout)           :: attributeValue
    class    (hdf5Object      )                                 , intent(inout)           :: self
    character(len=*           )                                 , intent(in   ), optional :: attributeName
    integer  (kind=SIZE_T     )                                 , intent(in   )           :: dataTypeSize
    character(len=dataTypeSize), dimension(size(attributeValue))                          :: temporaryBuffer

    ! Call the character version of this routine to perform the red.
    call IO_HDF5_Read_Attribute_Character_1D_Array_Static(self,attributeName,temporaryBuffer)

    ! Transfer the results to the varying string variable.
    attributeValue=temporaryBuffer

    return
  end subroutine IO_HDF5_Read_Attribute_VarString_1D_Array_Static_Do_Read

  logical function IO_HDF5_Has_Attribute(self,attributeName)
    !!{
    Check if {\normalfont \ttfamily self} has an attribute with the given {\normalfont \ttfamily attributeName}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : h5aexists_f
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    implicit none
    class    (hdf5Object    ), intent(in   ) :: self
    character(len=*         ), intent(in   ) :: attributeName
    integer                                  :: errorCode
    type     (varying_string)                :: message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="object '"//self%objectName//"' in not open"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object exists.
    call h5aexists_f(self%objectID,trim(attributeName),IO_HDF5_Has_Attribute,errorCode)
    if (errorCode /= 0) then
       message="unable to check for attribute '"//trim(attributeName)//"' in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end function IO_HDF5_Has_Attribute

  subroutine IO_HDF5_Assert_Attribute_Type(attributeObject,attributeAssertedType,attributeAssertedRank,matches)
    !!{
    Asserts that an attribute is of a certain type and rank.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T                       , h5aget_space_f, h5aget_type_f, h5sclose_f, &
          &                           h5sget_simple_extent_ndims_f, h5tclose_f    , h5tequal_f
    use :: ISO_Varying_String, only : assignment(=)               , operator(//)
    implicit none
    class  (hdf5Object    )              , intent(in   )           :: attributeObject
    integer                              , intent(in   )           :: attributeAssertedRank
    integer(kind=HID_T    ), dimension(:), intent(in   )           :: attributeAssertedType
    logical                              , intent(  out), optional :: matches
    integer                                                        :: attributeRank        , errorCode
    integer(kind=HID_T    )                                        :: attributeDataspaceID , attributeTypeID
    logical                                                        :: isCorrectType
    integer                                                        :: iType
    type   (varying_string)                                        :: message

    ! Set the return value if present.
    if (present(matches)) matches=.true.

    ! Check the attribute type
    call h5aget_type_f(attributeObject%objectID,attributeTypeID,errorCode)
    if (errorCode /= 0) then
       message="unable to get datatype of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//attributeObject%locationReport()//{introspection:location})
    end if
    isCorrectType=.false. ! Assume that it is of the incorrect type by default.
    do iType=1,size(attributeAssertedType)
       call h5tequal_f(attributeTypeID,attributeAssertedType(iType),isCorrectType,errorCode)
       if (errorCode /= 0) then
          message="unable to test datatype of attribute '"//attributeObject%objectName//"'"
          call Error_Report(message//attributeObject%locationReport()//{introspection:location})
       end if
       ! If a suitable type match has been found, exit the loop.
       if (isCorrectType) exit
    end do
    call h5tclose_f(attributeTypeID,errorCode)
    if (errorCode /= 0) then
       message="unable to close datatype of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//attributeObject%locationReport()//{introspection:location})
    end if
    if (.not.isCorrectType) then
       if (present(matches)) then
          matches=.false.
          return
       else
          message="attribute '"//attributeObject%objectName//"' is of incorrect type"
          call Error_Report(message//attributeObject%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the attribute has the correct rank.
    call h5aget_space_f(attributeObject%objectID,attributeDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//attributeObject%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_ndims_f(attributeDataspaceID,attributeRank,errorCode)
    if (errorCode /= 0) then
       message="unable to get rank of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//attributeObject%locationReport()//{introspection:location})
    end if
    call h5sclose_f(attributeDataspaceID,errorCode)
     if (errorCode /= 0) then
       message="unable to close dataspace of attribute '"//attributeObject%objectName//"'"
       call Error_Report(message//attributeObject%locationReport()//{introspection:location})
    end if
    if (attributeRank /= attributeAssertedRank) then
       if (present(matches)) then
          matches=.false.
          return
       else
          message="attribute '"//attributeObject%objectName//"' has incorrect rank"
          call Error_Report(message//attributeObject%locationReport()//{introspection:location})
       end if
    end if

    return
  end subroutine IO_HDF5_Assert_Attribute_Type

  !! Dataset routines.

  function IO_HDF5_Dataset_Size(datasetObject,dim)
    !!{
    Return the size of the {\normalfont \ttfamily dim}$^\mathrm{th}$ dimension of dataset {\normalfont \ttfamily datasetObject}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T                      , HSIZE_T                     , h5dget_space_f, h5sclose_f, &
          &                           h5sget_simple_extent_dims_f, h5sget_simple_extent_ndims_f
    use :: ISO_Varying_String, only : assignment(=)              , operator(//)
    implicit none
    integer(kind=HSIZE_T  )                              :: IO_HDF5_Dataset_Size
    class  (hdf5Object    ), intent(in   )               :: datasetObject
    integer                , intent(in   )               :: dim
    integer(kind=HSIZE_T  ), allocatable  , dimension(:) :: dimensions          , maximumDimensions
    integer                                              :: datasetRank         , errorCode
    integer(kind=HID_T    )                              :: datasetDataspaceID
    type   (varying_string)                              :: message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Ensure that the object is a dataset.
    if (datasetObject%hdf5ObjectType /= hdf5ObjectTypeDataset) then
       message="object is not a dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if

    ! Ensure that the dataset is open.
    if (.not.datasetObject%isOpenValue) then
       message="attempt to get size of unopen dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if

    ! Get the rank of the dataset
    call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_ndims_f(datasetDataspaceID,datasetRank,errorCode)
    if (errorCode /= 0) then
       message="unable to get rank of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if
    if (datasetRank < dim) then
       message="dataset '"//datasetObject%objectName//"' has rank smaller than the dimension requested"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if

    ! Get the dimensions of the dataspace.
    allocate  (       dimensions(datasetRank))
    allocate  (maximumDimensions(datasetRank))
    call h5sget_simple_extent_dims_f(datasetDataspaceID,dimensions,maximumDimensions,errorCode)
    if (errorCode == -1) then
       message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if
    IO_HDF5_Dataset_Size=dimensions(dim)
    deallocate(       dimensions             )
    deallocate(maximumDimensions             )

    ! Close the dataspace
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if

    return
  end function IO_HDF5_Dataset_Size

  integer function IO_HDF5_Dataset_Rank(datasetObject)
    !!{
    Return the rank of dataset {\normalfont \ttfamily datasetObject}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , h5dget_space_f, h5sget_simple_extent_ndims_f
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    implicit none
    class  (hdf5Object    ), intent(in   )               :: datasetObject
    integer                                              :: datasetRank         , errorCode
    integer(kind=HID_T    )                              :: datasetDataspaceID
    type   (varying_string)                              :: message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Ensure that the object is a dataset.
    if (datasetObject%hdf5ObjectType /= hdf5ObjectTypeDataset) then
       message="object is not a dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if

    ! Ensure that the dataset is open.
    if (.not.datasetObject%isOpenValue) then
       message="attempt to get size of unopen dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if

    ! Get the rank of the dataset
    call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_ndims_f(datasetDataspaceID,datasetRank,errorCode)
    if (errorCode /= 0) then
       message="unable to get rank of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if

    ! Return the rank.
    IO_HDF5_Dataset_Rank=datasetRank
    return
  end function IO_HDF5_Dataset_Rank

  function IO_HDF5_Open_Dataset(inObject,datasetName,comment,datasetDataType,datasetDimensions,isOverwritable,appendTo,appendDimension,useDataType,chunkSize,compressionLevel) result(self)
    !!{
    Open an dataset in {\normalfont \ttfamily inObject}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5P_DATASET_CREATE_F, H5S_UNLIMITED_F      , H5T_NATIVE_CHARACTER, H5T_NATIVE_DOUBLE , &
          &                           H5T_NATIVE_INTEGER  , h5screate_simple_f   , HID_T               , HSIZE_T           , &
          &                           h5dcreate_f         , h5dget_create_plist_f, h5dopen_f           , h5eset_auto_f     , &
          &                           hsize_t             , h5pclose_f           , h5pcreate_f         , h5pget_chunk_f    , &
          &                           h5pset_chunk_f      , h5pset_deflate_f     , h5sclose_f
    use :: ISO_Varying_String, only : assignment(=)       , operator(//)         , operator(/=)
    implicit none
    type     (hdf5Object    )                                        :: self
    character(len=*         )              , intent(in   )           :: datasetName
    character(len=*         )              , intent(in   ), optional :: comment
    integer  (hsize_t       )              , intent(in   ), optional :: chunkSize
    integer                                , intent(in   ), optional :: compressionLevel                           , datasetDataType, &
         &                                                              appendDimension
    integer  (kind=HSIZE_T  ), dimension(:), intent(in   ), optional :: datasetDimensions
    logical                                , intent(in   ), optional :: appendTo                                   , isOverwritable
    integer  (kind=HID_T    )              , intent(in   ), optional :: useDataType
    class    (hdf5Object    )              , intent(in   ), target   :: inObject
    integer  (kind=HSIZE_T  ), parameter                             :: chunkSizeMaximum        =4294967296_hsize_t ! Maximum chunk size of 4GB (see https://support.hdfgroup.org/documentation/hdf5-docs/advanced_topics/chunking_in_hdf5.html).
    integer  (kind=HSIZE_T  ), dimension(7)                          :: chunkDimensions                            , datasetDimensionsActual, &
         &                                                              datasetDimensionsMaximum
    integer  (kind=HSIZE_T  )                                        :: chunkSizeActual
    integer                                                          :: compressionLevelActual                     , datasetRank            , &
         &                                                              errorCode                                  , appendDimensionActual  , &
         &                                                              iDimension
    integer  (kind=HID_T    )                                        :: dataSpaceID                                , dataTypeID             , &
         &                                                              locationID                                 , propertyList
    logical                                                          :: appendToActual
    type     (varying_string)                                        :: locationPath                               , message
    class    (*             ), pointer                               :: dummyPointer_

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Ensure that the object is already open.
    if (inObject%isOpenValue) then
       locationID                                     =inObject    %objectID
       locationPath                                   =inObject    %objectFile
       if (inObject%objectLocation /= "") locationPath=locationPath           //"/"//inObject%objectLocation
       if (inObject%objectName     /= "") locationPath=locationPath           //"/"//inObject%objectName
       select type (inObject)
       type is (hdf5Object)
          self%parentObject => inObject
       end select
    else
       message="attempt to open dataset '"//trim(datasetName)//"' in unopen object '"//inObject%objectName//"'"
       call Error_Report(message//inObject%locationReport()//{introspection:location})
    end if

    ! Determine the rank and dimensions.
    if (present(datasetDimensions)) then
       ! Open data space with the desired dimensions.
       datasetRank=size(datasetDimensions)
       if (datasetRank > 7) then
          message="datasets of rank greater than 7 are not supported - dataset in question is '"//trim(datasetName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       datasetDimensionsActual(1:datasetRank)=datasetDimensions
    else
       ! No dimensions specified, assume a scalar.
       datasetRank=0
    end if

    ! Check whether appending was requested.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    if (present(appendDimension)) then
       appendDimensionActual=appendDimension
    else
       appendDimensionActual=1
    end if
    ! Obtain a reference to the file ID.
    self%fileID      => inObject%fileID
    self%fileManager =  inObject%fileManager
    ! Create an ID for this dataset.
    allocate(self%objectID)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_ => self%objectID
    self%objectManager=resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    ! Store the name and location of the object.
    self%objectFile    =self%parentObject%objectFile
    self%objectLocation=self%parentObject%pathTo    (includeFileName=.false.)
    self%objectName    =trim(datasetName)
    ! Check if the dataset exists.
    if (inObject%hasDataset(datasetName)) then
       ! Open the dataset.
       call h5dopen_f(locationID,trim(datasetName),self%objectID,errorCode)
       if (errorCode /= 0) then
          message="failed to open dataset '"//trim(datasetName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       call h5dget_create_plist_f(self%objectID,propertyList,errorCode)
       if (errorCode /= 0) then
          message="failed to get creation property list for dataset '"//trim(datasetName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       call h5eset_auto_f(0,errorCode)
       if (errorCode /= 0) then
          message="failed to switch HDF5 error report off"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       call h5pget_chunk_f(propertyList,1,chunkDimensions,errorCode)
       if (errorCode < 0) then
          ! Assume that a failed attempt to get chunk size indicates that the dataset is not chunked.
          self%chunkSize=-1
       else
          self%chunkSize=int(chunkDimensions(1))
       end if
       call h5eset_auto_f(1,errorCode)
       if (errorCode /= 0) then
          message="failed to switch HDF5 error report on"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       call h5pclose_f(propertyList,errorCode)
       if (errorCode /= 0) then
          message="failed to close property list for dataset '"//trim(datasetName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
    else
       ! Determine maximum dimensions of this dataset.
       select case (appendToActual)
       case (.true. )
          datasetDimensionsMaximum=H5S_UNLIMITED_F
       case (.false.)
          datasetDimensionsMaximum=datasetDimensionsActual
       end select
       ! Create a suitable dataspace.
       call h5screate_simple_f(datasetRank,datasetDimensionsActual,dataspaceID,errorCode,datasetDimensionsMaximum)
       if (errorCode < 0) then
          message="unable to open dataspace for dataset '"//trim(datasetName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       ! Determine the chunking level.
       if (present(chunkSize)) then
          ! Check that chunk size is valid.
          if (chunkSize == 0 .or. chunkSize < -1) then
             message="invalid chunk size for dataset '"//trim(datasetName)//"' at "//locationPath
             call Error_Report(message//inObject%locationReport()//{introspection:location})
          end if
          chunkSizeActual=chunkSize
       else
          ! No chunk size explicitly provided. Inherit that of parent object if possible, otherwise use default.
          if (inObject%chunkSizeSet) then
             chunkSizeActual=inObject%chunkSize
          else
             chunkSizeActual=hdf5ChunkSize
          end if
       end if
       self%chunkSize=int(chunkSizeActual)
       ! Determine the compression level.
       if (present(compressionLevel)) then
          ! Check that compression level is valid.
          if (compressionLevel == 0 .or. compressionLevel < -1) then
             message="invalid compression level for dataset '"//trim(datasetName)//"' at "//locationPath
             call Error_Report(message//inObject%locationReport()//{introspection:location})
          end if
          compressionLevelActual=compressionLevel
       else
          ! No compression level explicitly provided. Inherit that of parent object if possible, otherwise use default.
          if (inObject%compressionLevelSet) then
             compressionLevelActual=inObject%compressionLevel
          else
             compressionLevelActual=hdf5CompressionLevel
          end if
       end if
       self%compressionLevel=compressionLevelActual
       ! Create a property list for the dataset.
       call h5pcreate_f(H5P_DATASET_CREATE_F,propertyList,errorCode)
       if (errorCode < 0) then
          message="unable to create property list for dataset '"//trim(datasetName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       ! Check if chunk size needs to be set.
       if (chunkSizeActual /= -1) then
          ! It does - determine a suitable chunk size.
          if (appendToActual) then
             ! Extensible dataset, use selected chunk size.
             chunkDimensions(1:datasetRank          )=datasetDimensions
             chunkDimensions(  appendDimensionActual)=chunkSizeActual
          else
             ! Fixed dimension array, use smaller of chunk size and actual size.
             chunkDimensions(1:datasetRank)=min(datasetDimensionsActual(1:datasetRank),chunkSizeActual)
          end if
          ! Reduce chunk size if needed to fit within HDF5's maximum chunk size.
          iDimension=0
          do while (product(chunkDimensions(1:datasetRank))*sizeof(0.0d0) > chunkSizeMaximum)
             iDimension=iDimension+1
             if (iDimension > datasetRank) iDimension=1
             chunkDimensions(iDimension)=max(1_hsize_t,chunkDimensions(iDimension)/2_hsize_t)
          end do
          call h5pset_chunk_f(propertyList,datasetRank,chunkDimensions,errorCode)
          if (errorCode < 0) then
             message="unable to set chunk size for dataset '"//trim(datasetName)//"'"
             call Error_Report(message//inObject%locationReport()//{introspection:location})
          end if
       else
          ! No chunk size was specified. This is problematic if the dataset is appendable.
          if (appendToActual) then
             message="appendable dataset '"//trim(datasetName)//"' requires a chunk size"
             call Error_Report(message//inObject%locationReport()//{introspection:location})
          end if
       end if

       ! Check if compression level should be set.
       if (compressionLevelActual >= 0) then
          call h5pset_deflate_f(propertyList,compressionLevelActual,errorCode)
          if (errorCode < 0) then
             message="could not set compression level for dataset '"//trim(datasetName)//"'"
             call Error_Report(message//inObject%locationReport()//{introspection:location})
          end if
       end if
       ! Ensure that a data type was specified.
       if (.not.present(datasetDataType)) then
          message="no datatype was specified for dataset '"//trim(datasetName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       ! Determine the data type.
       if (present(useDataType)) then
          dataTypeID=useDataType
       else
          select case (datasetDataType)
          case (hdf5DataTypeInteger       )
             dataTypeID=H5T_NATIVE_INTEGER
          case (hdf5DataTypeInteger8      )
             dataTypeID=H5T_INTEGER8
          case (hdf5DataTypeDouble        )
             dataTypeID=H5T_NATIVE_DOUBLE
          case (hdf5DataTypeCharacter     )
             dataTypeID=H5T_NATIVE_CHARACTER
          case (hdf5DataTypeVlenDouble    )
             dataTypeID=H5T_VLEN_DOUBLE     (1)
          case (hdf5DataTypeVlenVlenDouble)
             dataTypeID=H5T_VLEN_VLEN_DOUBLE(1)
          case (hdf5DataTypeVlenInteger8  )
             dataTypeID=H5T_VLEN_INTEGER8   (1)
          end select
       end if
       ! Create the dataset.
       call h5dcreate_f(locationID,trim(datasetName),dataTypeID,dataSpaceID,self%objectID,errorCode,propertyList)
       if (errorCode /= 0) then
          message="failed to create dataset '"//trim(datasetName)//"' at "//locationPath
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       ! Close the dataspace.
       call h5sclose_f(dataSpaceID,errorCode)
       if (errorCode /= 0) then
          message="failed to close dataspace for dataset '"//trim(datasetName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       ! Close the property list.
       call h5pclose_f(propertyList,errorCode)
       if (errorCode /= 0) then
          message="failed to close property list for dataset '"//trim(datasetName)//"'"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
    end if

    ! Mark this object as open.
    self%isOpenValue=.true.

    ! Mark this object as a file object.
    self%hdf5ObjectType=hdf5ObjectTypeDataset

    ! Mark whether dataset is overwritable.
    if (present(isOverwritable)) then
       ! Check overwriting is not requested if parent is not overwritable.
       if (isOverwritable.and..not.self%parentObject%isOverwritable) then
          message="cannot make dataset '"//trim(datasetName)//"' overwritable as objects in parent '"&
               &//self%parentObject%objectName//"' are not overwritable"
          call Error_Report(message//inObject%locationReport()//{introspection:location})
       end if
       self%isOverwritable=isOverwritable
    else
       self%isOverwritable=self%parentObject%isOverwritable
    end if

    ! Set the comment for this dataset.
    if (present(comment) .and. len_trim(comment) > 0 .and. .not.self%hasAttribute('comment')) call self%writeAttribute(trim(comment),'comment')
    return
  end function IO_HDF5_Open_Dataset

  logical function IO_HDF5_Has_Dataset(self,datasetName)
    !!{
    Check if {\normalfont \ttfamily self} has a dataset with the given {\normalfont \ttfamily datasetName}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , h5dclose_f  , h5dopen_f, h5eset_auto_f
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    implicit none
    class    (hdf5Object    ), intent(in   ) :: self
    character(len=*         ), intent(in   ) :: datasetName
    integer                                  :: errorCode
    integer  (kind=HID_T    )                :: datasetID
    type     (varying_string)                :: message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="object '"//self%objectName//"' in not open"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object exists.
    call h5eset_auto_f(0,errorCode)
    if (errorCode /= 0) then
       message="failed to switch HDF5 error report off"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5dopen_f(self%objectID,trim(datasetName),datasetID,errorCode)
    IO_HDF5_Has_Dataset=(errorCode == 0)
    call h5eset_auto_f(1,errorCode)
    if (errorCode /= 0) then
       message="failed to switch HDF5 error report on"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    if (IO_HDF5_Has_Dataset) then
       call h5dclose_f(datasetID,errorCode)
       if (errorCode /= 0) then
          message="failed to close dataset '"//trim(datasetName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    return
  end function IO_HDF5_Has_Dataset

  subroutine IO_HDF5_Datasets(self,datasetNames)
    !!{
    Return a list of all datasets present within {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : h5g_dataset_f, h5gget_obj_info_idx_f, h5gn_members_f, hid_t
    use :: ISO_Varying_String, only : assignment(=), char                 , operator(//)  , trim
    use :: String_Handling   , only : operator(//)
    implicit none
    type     (varying_string), intent(inout), allocatable, dimension(:) :: datasetNames
    class    (hdf5Object    ), intent(in   )                            :: self
    integer                                                             :: errorCode
    integer  (hid_t         )                                           :: locationIdentifier
    type     (varying_string)                                           :: message           , objectName
    character(len=1024      )                                           :: memberName
    integer                                                             :: memberType        , groupMemberCount, &
         &                                                                 i                 , datasetCount

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized()
    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="object '"//self%objectName//"' in not open"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Extract relevant location identifier and object name.
    select case(self%hdf5ObjectType)
    case (hdf5ObjectTypeGroup)
       locationIdentifier=self%parentObject%objectID
       objectName        =self             %objectName
    case (hdf5ObjectTypeFile )
       locationIdentifier=self%objectID
       objectName        ="/"
    case default
       message="object '"//self%objectName//"' is not a group or file"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end select
    ! Get a count of the number of members in the group.
    call h5gn_members_f(locationIdentifier,char(objectName),groupMemberCount,errorCode)
    if (errorCode /= 0) then
       message="failed to get count of members in '"//trim(objectName)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Iterate over members, counting datasets.
    datasetCount=0
    do i=0,groupMemberCount-1
       call h5gget_obj_info_idx_f(locationIdentifier,char(objectName),i,memberName,memberType,errorCode)
       if (errorCode /= 0) then
          message="failed to get info on member "
          message=message//i//" in '"//trim(objectName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Count datasets.
       if (memberType == h5g_dataset_f) datasetCount=datasetCount+1
    end do
    ! Allocate the array of dataset names and retrieve them.
    if (allocated(datasetNames)) deallocate(datasetNames)
    allocate(datasetnames(datasetCount))
    datasetCount=0
    do i=0,groupMemberCount-1
       call h5gget_obj_info_idx_f(locationIdentifier,char(objectName),i,memberName,memberType,errorCode)
       if (errorCode /= 0) then
          message="failed to get info on member "
          message=message//i//" in '"//trim(objectName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Store dataset name.
       if (memberType == h5g_dataset_f) then
          datasetCount=datasetCount+1
          datasetNames(datasetCount)=trim(memberName)
       end if
    end do
    return
  end subroutine IO_HDF5_Datasets

  subroutine IO_HDF5_Assert_Dataset_Type(datasetObject,datasetAssertedType,datasetAssertedRank,status)
    !!{
    Asserts that an dataset is of a certain type and rank.
    !!}
    use :: Error             , only : Error_Report                , errorStatusFail, errorStatusSuccess
    use :: HDF5              , only : HID_T                       , h5dget_space_f , h5dget_type_f     , h5sclose_f, &
          &                           h5sget_simple_extent_ndims_f, h5tclose_f     , h5tequal_f
    use :: ISO_Varying_String, only : assignment(=)               , operator(//)
    implicit none
    class  (hdf5Object    )              , intent(in   ) :: datasetObject
    integer                              , intent(in   ) :: datasetAssertedRank
    integer(kind=HID_T    ), dimension(:), intent(in   ) :: datasetAssertedType
    integer                , optional    , intent(  out) :: status
    integer                                              :: datasetRank        , errorCode
    integer(kind=HID_T    )                              :: datasetDataspaceID , datasetTypeID
    logical                                              :: isCorrectType
    integer                                              :: iType
    type   (varying_string)                              :: message

    ! Set status to success by default.
    if (present(status)) status=errorStatusSuccess
    ! Check the dataset type
    call h5dget_type_f(datasetObject%objectID,datasetTypeID,errorCode)
    if (errorCode /= 0) then
       message="unable to get datatype of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if
    isCorrectType=.false. ! Assume that it is of the incorrect type by default.
    do iType=1,size(datasetAssertedType)
       call h5tequal_f(datasetTypeID,datasetAssertedType(iType),isCorrectType,errorCode)
       if (errorCode /= 0) then
          message="unable to test datatype of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//datasetObject%locationReport()//{introspection:location})
       end if
       ! If a suitable type match has been found, exit the loop.
       if (isCorrectType) exit
    end do
    call h5tclose_f(datasetTypeID,errorCode)
    if (errorCode /= 0) then
       message="unable to close datatype of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if
    if (.not.isCorrectType) then
       if (present(status)) then
          status=errorStatusFail
          return
       else
          message="dataset '"//datasetObject%objectName//"' is of incorrect type"
          call Error_Report(message//datasetObject%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the dataset has the correct rank.
    call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if
    call h5sget_simple_extent_ndims_f(datasetDataspaceID,datasetRank,errorCode)
    if (errorCode /= 0) then
       message="unable to get rank of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if
    if (datasetRank /= datasetAssertedRank) then
       if (present(status)) then
          status=errorStatusFail
          return
       else
          message="dataset '"//datasetObject%objectName//"' has incorrect rank"
          call Error_Report(message//datasetObject%locationReport()//{introspection:location})
       end if
    end if
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//datasetObject%locationReport()//{introspection:location})
    end if

    return
  end subroutine IO_HDF5_Assert_Dataset_Type

  subroutine IO_HDF5_Write_Dataset_Integer_1D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write an integer 1-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_INTEGER         , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: comment                   , datasetName
    integer                  , dimension(:), intent(in   )           :: datasetValue
    logical                                , intent(in   ), optional :: appendTo
    integer  (hsize_t       )              , intent(in   ), optional :: chunkSize
    integer                                , intent(in   ), optional :: compressionLevel
    type     (hdf5Object    )              , intent(  out), optional :: datasetReturned
    integer  (kind=HSIZE_T  ), dimension(1)                          :: datasetDimensions          , hyperslabCount      , &
         &                                                              hyperslabStart             , newDatasetDimensions, &
         &                                                              newDatasetDimensionsMaximum
    integer                                                          :: datasetRank                , errorCode
    integer  (kind=HID_T    )                                        :: dataspaceID                , newDataspaceID
    logical                                                          :: appendToActual             , preExisted
    type     (hdf5Object    )                                        :: datasetObject
    type     (varying_string)                                        :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D integer.
          call self%assertDatasetType(H5T_NATIVE_INTEGERS,1)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeInteger,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 1D integer.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_INTEGERS,1)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=1
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_INTEGER,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Integer_1D

  subroutine IO_HDF5_Write_Dataset_Integer_2D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write an integer 2-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_INTEGER         , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class    (hdf5Object    )                , intent(inout)           :: self
    character(len=*         )                , intent(in   ), optional :: comment                    , datasetName
    integer                  , dimension(:,:), intent(in   )           :: datasetValue
    logical                                  , intent(in   ), optional :: appendTo
    integer  (hsize_t       )                , intent(in   ), optional :: chunkSize
    integer                                  , intent(in   ), optional :: compressionLevel
    type     (hdf5Object    )                , intent(  out), optional :: datasetReturned
    integer  (kind=HSIZE_T  ), dimension(2  )                          :: datasetDimensions          , hyperslabCount      , &
         &                                                                hyperslabStart             , newDatasetDimensions, &
         &                                                                newDatasetDimensionsMaximum
    integer                                                            :: datasetRank                , errorCode
    integer  (kind=HID_T    )                                          :: dataspaceID                , newDataspaceID
    logical                                                            :: appendToActual             , preExisted
    type     (hdf5Object    )                                          :: datasetObject
    type     (varying_string)                                          :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 2D integer.
          call self%assertDatasetType(H5T_NATIVE_INTEGERS,2)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeInteger,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 2D integer.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_INTEGERS,2)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=2
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_INTEGER,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Integer_2D

  subroutine IO_HDF5_Write_Dataset_Integer_3D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write an integer 3-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_INTEGER         , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none  
    class    (hdf5Object    )                  , intent(inout)           :: self
    character(len=*         )                  , intent(in   ), optional :: comment                    , datasetName
    integer                  , dimension(:,:,:), intent(in   )           :: datasetValue
    logical                                    , intent(in   ), optional :: appendTo
    integer  (hsize_t       )                  , intent(in   ), optional :: chunkSize
    integer                                    , intent(in   ), optional :: compressionLevel
    type     (hdf5Object    )                  , intent(  out), optional :: datasetReturned
    integer  (kind=HSIZE_T  ), dimension(3    )                          :: datasetDimensions          , hyperslabCount      , &
         &                                                                  hyperslabStart             , newDatasetDimensions, &
         &                                                                  newDatasetDimensionsMaximum
    integer                                                              :: datasetRank                , errorCode
    integer  (kind=HID_T    )                                            :: dataspaceID                , newDataspaceID
    logical                                                              :: appendToActual             , preExisted
    type     (hdf5Object    )                                            :: datasetObject
    type     (varying_string)                                            :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 3D integer.
          call self%assertDatasetType(H5T_NATIVE_INTEGERS,3)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeInteger,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 3D integer.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_INTEGERS,3)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=3
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_INTEGER,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Integer_3D

  subroutine IO_HDF5_Read_Dataset_Integer_1D_Array_Static(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read an integer scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_INTEGER         , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    integer                     , dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object       )              , intent(inout)           :: self
    character(len=*            )              , intent(in   ), optional :: datasetName
    integer  (kind=HSIZE_T     ), dimension(1), intent(in   ), optional :: readBegin         , readCount
    integer  (kind=HSIZE_T     ), dimension(1)                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                 referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (hdset_reg_ref_t_f), save        , target                  :: referencedRegion
    integer                                                             :: errorCode
    integer  (kind=HID_T       )                                        :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                 memorySpaceID     , storedDatasetID
    logical                                                             :: isReference       , readSubsection
    type     (hdf5Object       )                                        :: datasetObject
    type     (varying_string   )                                        :: datasetNameActual , message
    type     (c_ptr            )                                        :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D integer array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGERS,1)

    ! Get the dimensions of the array to be read.
    storedDatasetID=0
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_INTEGER,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer_1D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Integer_1D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read an integer scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_INTEGER         , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    integer                     , allocatable, dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object       )                           , intent(inout)           :: self
    character(len=*            )                           , intent(in   ), optional :: datasetName
    integer  (kind=HSIZE_T     )             , dimension(1), intent(in   ), optional :: readBegin         , readCount
    integer  (kind=HSIZE_T     )             , dimension(1)                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                              referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (hdset_reg_ref_t_f), save       , target                                :: referencedRegion
    integer                                                                          :: errorCode
    integer  (kind=HID_T       )                                                     :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                              memorySpaceID     , storedDatasetID
    logical                                                                          :: isReference       , readSubsection
    type     (hdf5Object       )                                                     :: datasetObject
    type     (varying_string   )                                                     :: datasetNameActual , message
    type     (c_ptr            )                                                     :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D integer array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGERS,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]
    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_INTEGER,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Dataset_Integer_2D_Array_Static(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read an integer scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_INTEGER         , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    integer                     , dimension(:,:), intent(  out)           :: datasetValue
    class    (hdf5Object       )                , intent(inout)           :: self
    character(len=*            )                , intent(in   ), optional :: datasetName
    integer  (kind=HSIZE_T     ), dimension(2  ), intent(in   ), optional :: readBegin         , readCount
    integer  (kind=HSIZE_T     ), dimension(2  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                   referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (hdset_reg_ref_t_f), save          , target                  :: referencedRegion
    integer                                                               :: errorCode
    integer  (kind=HID_T       )                                          :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                   memorySpaceID     , storedDatasetID
    logical                                                               :: isReference       , readSubsection
    type     (hdf5Object       )                                          :: datasetObject
    type     (varying_string   )                                          :: datasetNameActual , message
    type     (c_ptr            )                                          :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 2D integer array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGERS,2)

    ! Get the dimensions of the array to be read.
    storedDatasetID=0
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_INTEGER,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer_2D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Integer_2D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read an integer scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_INTEGER         , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    integer                     , allocatable, dimension(:,:), intent(  out)           :: datasetValue
    class    (hdf5Object       )                             , intent(inout)           :: self
    character(len=*            )                             , intent(in   ), optional :: datasetName
    integer  (kind=HSIZE_T     )             , dimension(2  ), intent(in   ), optional :: readBegin         , readCount
    integer  (kind=HSIZE_T     )             , dimension(2  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (hdset_reg_ref_t_f), save         , target                                :: referencedRegion
    integer                                                                            :: errorCode
    integer  (kind=HID_T       )                                                       :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                memorySpaceID     , storedDatasetID
    logical                                                                            :: isReference       , readSubsection
    type     (hdf5Object       )                                                       :: datasetObject
    type     (varying_string   )                                                       :: datasetNameActual , message
    type     (c_ptr            )                                                       :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 2D integer array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGERS,2)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]
    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_INTEGER,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer_2D_Array_Allocatable

  subroutine IO_HDF5_Write_Dataset_Integer8_1D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a long integer 1-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F     , H5S_SELECT_SET_F           , h5sselect_hyperslab_f, HID_T     , &
          &                                      HSIZE_T           , h5dget_space_f             , h5dset_extent_f      , h5sclose_f, &
          &                                      h5screate_simple_f, h5sget_simple_extent_dims_f, hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class    (hdf5Object    )                           , intent(inout)                   :: self
    character(len=*         )                           , intent(in   ), optional         :: comment                   , datasetName
    integer  (kind=kind_int8)             , dimension(:), intent(in   )                   :: datasetValue
    logical                                             , intent(in   ), optional         :: appendTo
    integer  (hsize_t       )                           , intent(in   ), optional         :: chunkSize
    integer                                             , intent(in   ), optional         :: compressionLevel
    type     (hdf5Object    )                           , intent(  out), optional         :: datasetReturned
    integer  (kind=kind_int8), allocatable, dimension(:)                         , target :: datasetValueContiguous
    integer  (kind=HSIZE_T  )             , dimension(1)                                  :: datasetDimensions          , hyperslabCount      , &
         &                                                                                   hyperslabStart             , newDatasetDimensions, &
         &                                                                                   newDatasetDimensionsMaximum
    integer                                                                               :: datasetRank                , errorCode
    integer  (kind=HID_T    )                                                             :: dataspaceID                , newDataspaceID
    logical                                                                               :: appendToActual             , preExisted
    type     (hdf5Object    )                                                             :: datasetObject
    type     (varying_string)                                                             :: datasetNameActual          , message
    type     (c_ptr         )                                                             :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D long integer.
          call self%assertDatasetType(H5T_NATIVE_INTEGER_8S,1)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeInteger8,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 1D long integer.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_INTEGER_8S,1)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=1
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    allocate(datasetValueContiguous,mold=datasetValue)
    datasetValueContiguous=datasetValue
    dataBuffer=c_loc(datasetValueContiguous)
    errorCode=h5dwrite(datasetObject%objectID,H5T_INTEGER8,newDataspaceID,dataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    deallocate(datasetValueContiguous)

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Integer8_1D

  subroutine IO_HDF5_Write_Dataset_Integer8_2D(self,datasetValue,datasetName,comment,appendTo,appendDimension,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a long integer 2-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F     , H5S_SELECT_SET_F           , h5sselect_hyperslab_f, HID_T     , &
          &                                      HSIZE_T           , h5dget_space_f             , h5dset_extent_f      , h5sclose_f, &
          &                                      h5screate_simple_f, h5sget_simple_extent_dims_f, hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class    (hdf5Object    )                             , intent(inout)                   :: self
    character(len=*         )                             , intent(in   ), optional         :: comment                     , datasetName
    integer  (kind=kind_int8)             , dimension(:,:), intent(in   )                   :: datasetValue
    logical                                               , intent(in   ), optional         :: appendTo
    integer  (hsize_t       )                             , intent(in   ), optional         :: chunkSize
    integer                                               , intent(in   ), optional         :: appendDimension             , compressionLevel
    type     (hdf5Object    )                             , intent(  out), optional         :: datasetReturned
    integer  (kind=kind_int8), allocatable, dimension(:,:)                         , target :: datasetValueContiguous
    integer  (kind=HSIZE_T  )             , dimension(2  )                                  :: datasetDimensions           , hyperslabCount             , &
         &                                                                                     hyperslabStart              , newDatasetDimensions       , &
         &                                                                                     newDatasetDimensionsFiltered, newDatasetDimensionsMaximum
    integer                                                                                 :: datasetRank                 , errorCode                  , &
         &                                                                                     appendDimensionActual
    integer  (kind=HID_T    )                                                               :: dataspaceID                 , newDataspaceID
    logical                                                                                 :: appendToActual              , preExisted
    type     (hdf5Object    )                                                               :: datasetObject
    type     (varying_string)                                                               :: datasetNameActual           , message
    type     (c_ptr         )                                                               :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 2D long integer.
          call self%assertDatasetType(H5T_NATIVE_INTEGER_8S,2)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeInteger8,datasetDimensions,appendTo&
            &=appendTo,appendDimension=appendDimension,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 2D long integer.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_INTEGER_8S,2)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Determine the dimension for appending.
       appendDimensionActual=1
       if (present(appendDimension)) appendDimensionActual=appendDimension
       ! Ensure that all dimensions other than the one being appended to are of the same size.
       newDatasetDimensionsFiltered                       =newDatasetDimensions
       newDatasetDimensionsFiltered(appendDimensionActual)=dataSetDimensions   (appendDimensionActual)
       if (any(dataSetDimensions /= newDatasetDimensionsFiltered)) then
          message="when appending to dataset '"//trim(datasetNameActual)//"' all dimensions other than that being appended to must be same as original dataset"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Set the hyperslab.
       hyperslabStart                             =0
       hyperslabStart      (appendDimensionActual)=newDatasetDimensions(appendDimensionActual)
       hyperslabCount                             =dataSetDimensions
       newDatasetDimensions(appendDimensionActual)=newDatasetDimensions(appendDimensionActual)+datasetDimensions(appendDimensionActual)
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=2
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset
    allocate(datasetValueContiguous,mold=datasetValue)
    datasetValueContiguous=datasetValue
    dataBuffer=c_loc(datasetValueContiguous)
    errorCode=h5dwrite(datasetObject%objectID,H5T_INTEGER8,newDataspaceID,dataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    deallocate(datasetValueContiguous)

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Integer8_2D

  subroutine IO_HDF5_Write_Dataset_Integer8_3D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a long integer 3-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F     , H5S_SELECT_SET_F           , h5sselect_hyperslab_f, HID_T     , &
          &                                      HSIZE_T           , h5dget_space_f             , h5dset_extent_f      , h5sclose_f, &
          &                                      h5screate_simple_f, h5sget_simple_extent_dims_f, hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class    (hdf5Object    )                               , intent(inout)                   :: self
    character(len=*         )                               , intent(in   ), optional         :: comment                    , datasetName
    integer  (kind=kind_int8)             , dimension(:,:,:), intent(in   )                   :: datasetValue
    logical                                                 , intent(in   ), optional         :: appendTo
    integer  (hsize_t       )                               , intent(in   ), optional         :: chunkSize
    integer                                                 , intent(in   ), optional         :: compressionLevel
    type     (hdf5Object    )                               , intent(  out), optional         :: datasetReturned
    integer  (kind=kind_int8), allocatable, dimension(:,:,:)                         , target :: datasetValueContiguous
    integer  (kind=HSIZE_T  )             , dimension(3    )                                  :: datasetDimensions          , hyperslabCount      , &
         &                                                                                       hyperslabStart             , newDatasetDimensions, &
         &                                                                                       newDatasetDimensionsMaximum
    integer                                                                                   :: datasetRank                , errorCode
    integer  (kind=HID_T    )                                                                 :: dataspaceID                , newDataspaceID
    logical                                                                                   :: appendToActual             , preExisted
    type     (hdf5Object    )                                                                 :: datasetObject
    type     (varying_string)                                                                 :: datasetNameActual          , message
    type     (c_ptr         )                                                                 :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 3D long integer.
          call self%assertDatasetType(H5T_NATIVE_INTEGER_8S,3)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeInteger8,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 3D long integer.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_INTEGER_8S,3)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=3
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    allocate(datasetValueContiguous,mold=datasetValue)
    datasetValueContiguous=datasetValue
    dataBuffer=c_loc(datasetValueContiguous)
    errorCode=h5dwrite(datasetObject%objectID,H5T_INTEGER8,newDataspaceID,dataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    deallocate(datasetValueContiguous)

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Integer8_3D

  subroutine IO_HDF5_Read_Dataset_Integer8_1D_Array_Static(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a long integer scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F             , H5S_SELECT_SET_F           , h5sselect_elements_f, &
          &                                      H5T_STD_REF_DSETREG  , HID_T                 , HSIZE_T                    , h5dclose_f          , &
          &                                      h5dget_space_f       , h5rdereference_f      , h5rget_region_f            , h5sclose_f          , &
          &                                      h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, size_t              , &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)          , trim
    implicit none
    integer  (kind=kind_int8   )             , dimension(:)          , intent(  out)           :: datasetValue
    class    (hdf5Object       )                                     , intent(inout)           :: self
    character(len=*            )                                     , intent(in   ), optional :: datasetName
    integer  (kind=HSIZE_T     )             , dimension(1)          , intent(in   ), optional :: readBegin             , readCount
    integer  (kind=HSIZE_T     )             , dimension(:)          , intent(in   ), optional :: readSelection
    integer  (kind=HSIZE_T     )             , dimension(1)                                    :: datasetDimensions     , datasetMaximumDimensions, &
         &                                                                                        referenceEnd          , referenceStart
    integer  (kind=HSIZE_T     ), allocatable, dimension(:,:)                                  :: readSelectionMap
    integer  (kind=kind_int8   ), allocatable, dimension(:)  , target                          :: datasetValueContiguous
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (hdset_reg_ref_t_f), save                       , target                          :: referencedRegion
    integer                                                                                    :: errorCode             , i
    integer  (kind=HID_T       )                                                               :: datasetDataspaceID    , dereferencedObjectID    , &
         &                                                                                        memorySpaceID         , storedDatasetID
    logical                                                                                    :: isReference           , readSubsection
    type     (hdf5Object       )                                                               :: datasetObject
    type     (varying_string   )                                                               :: datasetNameActual     , message
    type     (c_ptr            )                                                               :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D integer8 array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGER_8AS,1)

    ! Get the dimensions of the array to be read.
    storedDatasetID=0
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    allocate(datasetValueContiguous,mold=datasetValue)
    dataBuffer=c_loc(datasetValueContiguous)
    errorCode=h5dread(datasetObject%objectID,H5T_INTEGER8,memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    datasetValue=datasetValueContiguous
    deallocate(datasetValueContiguous)

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer8_1D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Integer8_1D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a long integer scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F             , H5S_SELECT_SET_F           , hdset_reg_ref_t_f   , &
          &                                      H5T_STD_REF_DSETREG  , HID_T                 , HSIZE_T                    , h5dclose_f          , &
          &                                      h5dget_space_f       , h5rdereference_f      , h5rget_region_f            , h5sclose_f          , &
          &                                      h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, h5sselect_elements_f, &
          &                                      h5sselect_hyperslab_f, size_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)          , trim
    implicit none
    integer  (kind=kind_int8   ), allocatable, dimension(:  ), intent(  out)          , target :: datasetValue
    class    (hdf5Object       )                             , intent(inout)                   :: self
    character(len=*            )                             , intent(in   ), optional         :: datasetName
    integer  (kind=HSIZE_T     )             , dimension(1  ), intent(in   ), optional         :: readBegin         , readCount
    integer  (kind=HSIZE_T     )             , dimension(:  ), intent(in   ), optional         :: readSelection
    integer  (kind=HSIZE_T     )             , dimension(1  )                                  :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                        referenceEnd      , referenceStart
    integer  (kind=HSIZE_T     ), allocatable, dimension(:,:)                                  :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (hdset_reg_ref_t_f), save                                              , target :: referencedRegion
    integer                                                                                  :: errorCode         , i
    integer  (kind=HID_T       )                                                             :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                      memorySpaceID     , storedDatasetID
    logical                                                                                  :: isReference       , readSubsection
    type     (hdf5Object       )                                                             :: datasetObject
    type     (varying_string   )                                                             :: datasetNameActual , message
    type     (c_ptr            )                                                             :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D long integer array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGER_8AS,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]

    ! Read the dataset.
    dataBuffer=c_loc(datasetValue)
    errorCode=h5dread(datasetObject%objectID,H5T_INTEGER8,memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer8_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Dataset_Integer8_2D_Array_Static(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a double scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F             , H5S_SELECT_SET_F           , hdset_reg_ref_t_f   , &
          &                                      H5T_STD_REF_DSETREG  , HID_T                 , HSIZE_T                    , h5dclose_f          , &
          &                                      h5dget_space_f       , h5rdereference_f      , h5rget_region_f            , h5sclose_f          , &
          &                                      h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, h5sselect_elements_f, &
          &                                      h5sselect_hyperslab_f, size_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)          , trim
    implicit none
    integer         (kind_int8        ), dimension(:,:), intent(  out), target   :: datasetValue
    class           (hdf5Object       )                , intent(inout)           :: self
    character       (len=*            )                , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     ), dimension(2  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     ), dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     ), dimension(2  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                          referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), dimension(:,:), allocatable             :: readSelectionMap
    integer         (kind=kind_int8   ), dimension(:,:), allocatable, target     :: datasetValueContiguous
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save          , target                  :: referencedRegion
    integer                                                                      :: errorCode
    integer         (kind=HID_T       )                                          :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                          memorySpaceID     , storedDatasetID
    integer         (kind=HSIZE_T     )                                          :: i                 , j
    logical                                                                      :: isReference       , readSubsection
    type            (hdf5Object       )                                          :: datasetObject
    type            (varying_string   )                                          :: datasetNameActual , message
    type            (c_ptr            )                                          :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 2D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGER_8AS,2)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(2))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(2,size(readSelection)*datasetDimensions(1)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                readSelectionMap(:,(i-1)*datasetDimensions(1)+j)=[j,readSelection(i)]
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,2,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(2)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(2))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
           ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(2,size(readSelection)*datasetDimensions(1)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                readSelectionMap(:,(i-1)*datasetDimensions(1)+j)=[j,readSelection(i)]
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,2,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(2)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    allocate(datasetValueContiguous,mold=datasetValue)
    dataBuffer=c_loc(datasetValueContiguous)
    errorCode=h5dread(datasetObject%objectID,H5T_INTEGER8,memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    datasetValue=datasetValueContiguous
    deallocate(datasetValueContiguous)

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer8_2D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Integer8_2D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a double 2-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F             , H5S_SELECT_SET_F           , hdset_reg_ref_t_f   , &
          &                                      H5T_STD_REF_DSETREG  , HID_T                 , HSIZE_T                    , h5dclose_f          , &
          &                                      h5dget_space_f       , h5rdereference_f      , h5rget_region_f            , h5sclose_f          , &
          &                                      h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, h5sselect_elements_f, &
          &                                      h5sselect_hyperslab_f, size_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)          , trim
    implicit none
    integer         (kind_int8        ), allocatable, dimension(:,:), intent(  out), target   :: datasetValue
    class           (hdf5Object       )                             , intent(inout)           :: self
    character       (len=*            )                             , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(2  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     )             , dimension(2  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                       referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), allocatable, dimension(:,:)                          :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                  :: referencedRegion
    integer                                                                                   :: errorCode
    integer         (kind=HID_T       )                                                       :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                       memorySpaceID     , storedDatasetID
    integer         (kind=HSIZE_T     )                                                       :: i                 , j
    logical                                                                                   :: isReference       , readSubsection
    type            (hdf5Object       )                                                       :: datasetObject
    type            (varying_string   )                                                       :: datasetNameActual , message
    type            (c_ptr            )                                                       :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 2D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGER_8AS,2)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(2))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(2,size(readSelection)*datasetDimensions(1)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                readSelectionMap(:,(i-1)*datasetDimensions(1)+j)=[j,readSelection(i)]
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,2,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(2)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(2))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
           ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(2,size(readSelection)*datasetDimensions(1)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                readSelectionMap(:,(i-1)*datasetDimensions(1)+j)=[j,readSelection(i)]
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,2,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(2)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]

    ! Read the dataset.
    dataBuffer=c_loc(datasetValue)
    errorCode=h5dread(datasetObject%objectID,H5T_INTEGER8,memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer8_2D_Array_Allocatable
  
  subroutine IO_HDF5_Read_Dataset_Integer8_3D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a double 3-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F             , H5S_SELECT_SET_F           , hdset_reg_ref_t_f   , &
          &                                      H5T_STD_REF_DSETREG  , HID_T                 , HSIZE_T                    , h5dclose_f          , &
          &                                      h5dget_space_f       , h5rdereference_f      , h5rget_region_f            , h5sclose_f          , &
          &                                      h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, h5sselect_elements_f, &
          &                                      h5sselect_hyperslab_f, size_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)          , trim
    implicit none
    integer         (kind_int8        ), allocatable, dimension(:,:,:), intent(  out), target   :: datasetValue
    class           (hdf5Object       )                               , intent(inout)           :: self
    character       (len=*            )                               , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(3    ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(:    ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     )             , dimension(3    )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                         referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), allocatable, dimension(:,:  )                          :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                    :: referencedRegion
    integer                                                                                     :: errorCode
    integer         (kind=HID_T       )                                                         :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                         memorySpaceID     , storedDatasetID
    integer         (kind=HSIZE_T     )                                                         :: i                 , j                       , &
         &                                                                                         k
    logical                                                                                     :: isReference       , readSubsection
    type            (hdf5Object       )                                                         :: datasetObject
    type            (varying_string   )                                                         :: datasetNameActual , message
    type            (c_ptr            )                                                         :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 3D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_INTEGER_8AS,3)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(3))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(3,size(readSelection)*datasetDimensions(1)*datasetDimensions(2)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                forall(k=1:datasetDimensions(2))
                   readSelectionMap(:,(i-1)*datasetDimensions(1)*datasetDimensions(2)+(j-1)*datasetDimensions(2)+k)=[j,k,readSelection(i)]
                end forall
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,3,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(3)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(3))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
           ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(3,size(readSelection)*datasetDimensions(1)*datasetDimensions(2)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                forall(k=1:datasetDimensions(2))
                   readSelectionMap(:,(i-1)*datasetDimensions(1)*datasetDimensions(2)+(j-1)*datasetDimensions(2)+k)=[j,k,readSelection(i)]
                end forall
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,3,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(3)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]

    ! Read the dataset.
    dataBuffer=c_loc(datasetValue)
    errorCode=h5dread(datasetObject%objectID,H5T_INTEGER8,memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Integer8_3D_Array_Allocatable
  
  subroutine IO_HDF5_Write_Dataset_Double_1D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a double 1-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_DOUBLE          , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class           (hdf5Object    )              , intent(inout)           :: self
    character       (len=*         )              , intent(in   ), optional :: comment                    , datasetName
    double precision                , dimension(:), intent(in   )           :: datasetValue
    logical                                       , intent(in   ), optional :: appendTo
    integer  (hsize_t       )                     , intent(in   ), optional :: chunkSize
    integer                                       , intent(in   ), optional :: compressionLevel
    type            (hdf5Object    )              , intent(  out), optional :: datasetReturned
    integer         (kind=HSIZE_T  ), dimension(1)                          :: datasetDimensions          , hyperslabCount      , &
         &                                                                     hyperslabStart             , newDatasetDimensions, &
         &                                                                     newDatasetDimensionsMaximum
    integer                                                                 :: datasetRank                , errorCode
    integer         (kind=HID_T    )                                        :: dataspaceID                , newDataspaceID
    logical                                                                 :: appendToActual             , preExisted
    type            (hdf5Object    )                                        :: datasetObject
    type            (varying_string)                                        :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D double.
          call self%assertDatasetType(H5T_NATIVE_DOUBLES,1)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeDouble,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 1D double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,1)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=1
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Double_1D

  subroutine IO_HDF5_Read_Dataset_Double_1D_Array_Static(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a double scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F       , H5S_ALL_F            , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG , HID_T                , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f      , h5dread_f            , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f          , h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_elements_f, h5sselect_hyperslab_f, hdset_reg_ref_t_f     , hsize_t                    , &
          &                                      size_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)         , trim
    implicit none
    double precision                   , dimension(:)  , intent(  out)           :: datasetValue
    class           (hdf5Object       )                , intent(inout)           :: self
    character       (len=*            )                , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     ), dimension(1)  , intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     ), dimension(:)  , intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     ), dimension(1)                            :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                          referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), dimension(:,:), allocatable             :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save          , target                  :: referencedRegion
    integer                                                                      :: errorCode         , i
    integer         (kind=HID_T       )                                          :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                          memorySpaceID     , storedDatasetID
    logical                                                                      :: isReference       , readSubsection
    type            (hdf5Object       )                                          :: datasetObject
    type            (varying_string   )                                          :: datasetNameActual , message
    type            (c_ptr            )                                          :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
                   ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_1D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Double_1D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a double scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F       , H5S_ALL_F            , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG , HID_T                , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f      , h5dread_f            , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f          , h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_elements_f, h5sselect_hyperslab_f, hdset_reg_ref_t_f     , hsize_t                    , &
          &                                      size_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)         , trim
    implicit none
    double precision                   , allocatable, dimension(:  ), intent(  out)           :: datasetValue
    class           (hdf5Object       )                             , intent(inout)           :: self
    character       (len=*            )                             , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(1  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     )             , dimension(1  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                       referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), allocatable, dimension(:,:)                          :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                  :: referencedRegion
    integer                                                                                   :: errorCode         , i
    integer         (kind=HID_T       )                                                       :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                       memorySpaceID     , storedDatasetID
    logical                                                                                   :: isReference       , readSubsection
    type            (hdf5Object       )                                                       :: datasetObject
    type            (varying_string   )                                                       :: datasetNameActual , message
    type            (c_ptr            )                                                       :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_1D_Array_Allocatable

  subroutine IO_HDF5_Write_Dataset_Double_2D(self,datasetValue,datasetName,comment,appendTo,appendDimension,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a double 2-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_DOUBLE          , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class           (hdf5Object    )                , intent(inout)           :: self
    character       (len=*         )                , intent(in   ), optional :: comment                     , datasetName
    double precision                , dimension(:,:), intent(in   )           :: datasetValue
    logical                                         , intent(in   ), optional :: appendTo
    integer  (hsize_t       )                       , intent(in   ), optional :: chunkSize
    integer                                         , intent(in   ), optional :: appendDimension             , compressionLevel
    type            (hdf5Object    )                , intent(  out), optional :: datasetReturned
    integer         (kind=HSIZE_T  ), dimension(2)                            :: datasetDimensions           , hyperslabCount             , &
         &                                                                       hyperslabStart              , newDatasetDimensions       , &
         &                                                                       newDatasetDimensionsFiltered, newDatasetDimensionsMaximum
    integer                                                                   :: appendDimensionActual       , datasetRank                , &
         &                                                                       errorCode
    integer         (kind=HID_T    )                                          :: dataspaceID                 , newDataspaceID
    logical                                                                   :: appendToActual              , preExisted
    type            (hdf5Object    )                                          :: datasetObject
    type            (varying_string)                                          :: datasetNameActual           , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 2D double.
          call self%assertDatasetType(H5T_NATIVE_DOUBLES,2)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeDouble,datasetDimensions,appendTo&
            &=appendTo,appendDimension=appendDimension,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 2D double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,2)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Determine the dimension for appending.
       appendDimensionActual=1
       if (present(appendDimension)) appendDimensionActual=appendDimension
       ! Ensure that all dimensions other than the one being appended to are of the same size.
       newDatasetDimensionsFiltered                       =newDatasetDimensions
       newDatasetDimensionsFiltered(appendDimensionActual)=dataSetDimensions   (appendDimensionActual)
       if (any(dataSetDimensions /= newDatasetDimensionsFiltered)) then
          message="when appending to dataset '"//trim(datasetNameActual)//"' all dimensions other than that being appended to must be same as original dataset"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Set the hyperslab.
       hyperslabStart                             =0
       hyperslabStart      (appendDimensionActual)=newDatasetDimensions(appendDimensionActual)
       hyperslabCount                             =dataSetDimensions
       newDatasetDimensions(appendDimensionActual)=newDatasetDimensions(appendDimensionActual)+datasetDimensions(appendDimensionActual)
    else
       newDatasetDimensions                       =datasetDimensions
       hyperslabStart                             =0
       hyperslabCount                             =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=2
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Double_2D

  subroutine IO_HDF5_Read_Dataset_Double_2D_Array_Static(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a double scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F       , H5S_ALL_F            , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG , HID_T                , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f      , h5dread_f            , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f          , h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_elements_f, h5sselect_hyperslab_f, hdset_reg_ref_t_f     , hsize_t                    , &
          &                                      size_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)         , trim
    implicit none
    double precision                   , dimension(:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                , intent(inout)           :: self
    character       (len=*            )                , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     ), dimension(2  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     ), dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     ), dimension(2  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                          referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), dimension(:,:), allocatable             :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save          , target                  :: referencedRegion
    integer                                                                      :: errorCode
    integer         (kind=HID_T       )                                          :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                          memorySpaceID     , storedDatasetID
    integer         (kind=HSIZE_T     )                                          :: i                 , j
    logical                                                                      :: isReference       , readSubsection
    type            (hdf5Object       )                                          :: datasetObject
    type            (varying_string   )                                          :: datasetNameActual , message
    type            (c_ptr            )                                          :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 2D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,2)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(2))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(2,size(readSelection)*datasetDimensions(1)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                readSelectionMap(:,(i-1)*datasetDimensions(1)+j)=[j,readSelection(i)]
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,2,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(2)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(2))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
           ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(2,size(readSelection)*datasetDimensions(1)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                readSelectionMap(:,(i-1)*datasetDimensions(1)+j)=[j,readSelection(i)]
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,2,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(2)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_2D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Double_2D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a double 2-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F       , H5S_ALL_F            , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG , HID_T                , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f      , h5dread_f            , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f          , h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_elements_f, h5sselect_hyperslab_f, hdset_reg_ref_t_f     , hsize_t                    , &
          &                                      size_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)         , trim
    implicit none
    double precision                   , allocatable, dimension(:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                             , intent(inout)           :: self
    character       (len=*            )                             , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(2  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     )             , dimension(2  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                       referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), allocatable, dimension(:,:)                          :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                  :: referencedRegion
    integer                                                                                   :: errorCode
    integer         (kind=HID_T       )                                                       :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                       memorySpaceID     , storedDatasetID
    integer         (kind=HSIZE_T     )                                                       :: i                 , j
    logical                                                                                   :: isReference       , readSubsection
    type            (hdf5Object       )                                                       :: datasetObject
    type            (varying_string   )                                                       :: datasetNameActual , message
    type            (c_ptr            )                                                       :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 2D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,2)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(2))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(2,size(readSelection)*datasetDimensions(1)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                readSelectionMap(:,(i-1)*datasetDimensions(1)+j)=[j,readSelection(i)]
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,2,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(2)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(2))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
           ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(2,size(readSelection)*datasetDimensions(1)))
          forall(i=1:size(readSelection))
             forall(j=1:datasetDimensions(1))
                readSelectionMap(:,(i-1)*datasetDimensions(1)+j)=[j,readSelection(i)]
             end forall
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,2,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(2)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]
    
    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_2D_Array_Allocatable

  subroutine IO_HDF5_Write_Dataset_Double_3D(self,datasetValue,datasetName,comment,appendTo,appendDimension,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a double 3-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_DOUBLE          , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class           (hdf5Object    )                  , intent(inout)           :: self
    character       (len=*         )                  , intent(in   ), optional :: comment                     , datasetName
    double precision                , dimension(:,:,:), intent(in   )           :: datasetValue
    logical                                           , intent(in   ), optional :: appendTo
    integer         (hsize_t       )                  , intent(in   ), optional :: chunkSize
    integer                                           , intent(in   ), optional :: appendDimension             , compressionLevel
    type            (hdf5Object    )                  , intent(  out), optional :: datasetReturned
    integer         (kind=HSIZE_T  ), dimension(3)                              :: datasetDimensions           , hyperslabCount             , &
         &                                                                         hyperslabStart              , newDatasetDimensions       , &
         &                                                                         newDatasetDimensionsFiltered, newDatasetDimensionsMaximum
    integer                                                                     :: appendDimensionActual       , datasetRank                , &
         &                                                                         errorCode
    integer         (kind=HID_T    )                                            :: dataspaceID                 , newDataspaceID
    logical                                                                     :: appendToActual              , preExisted
    type            (hdf5Object    )                                            :: datasetObject
    type            (varying_string)                                            :: datasetNameActual           , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 3D double.
          call self%assertDatasetType(H5T_NATIVE_DOUBLES,3)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeDouble,datasetDimensions,appendTo&
            &=appendTo,appendDimension=appendDimension,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 3D double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,3)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Determine the dimension for appending.
       appendDimensionActual=1
       if (present(appendDimension)) appendDimensionActual=appendDimension
       ! Ensure that all dimensions other than the one being appended to are of the same size.
       newDatasetDimensionsFiltered                       =newDatasetDimensions
       newDatasetDimensionsFiltered(appendDimensionActual)=dataSetDimensions   (appendDimensionActual)
       if (any(dataSetDimensions /= newDatasetDimensionsFiltered)) then
          message="when appending to dataset '"//trim(datasetNameActual)//"' all dimensions other than that being appended to must be same as original dataset"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Set the hyperslab.
       hyperslabStart                             =0
       hyperslabStart      (appendDimensionActual)=newDatasetDimensions(appendDimensionActual)
       hyperslabCount                             =dataSetDimensions
       newDatasetDimensions(appendDimensionActual)=newDatasetDimensions(appendDimensionActual)+datasetDimensions(appendDimensionActual)
    else
       newDatasetDimensions   =datasetDimensions
       hyperslabStart         =0
       hyperslabCount         =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=3
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Double_3D

  subroutine IO_HDF5_Read_Dataset_Double_3D_Array_Static(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a double scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    double precision                   , dimension(:,:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                  , intent(inout)           :: self
    character       (len=*            )                  , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     ), dimension(3)    , intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     ), dimension(3)                              :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                            referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save            , target                  :: referencedRegion
    integer                                                                        :: errorCode
    integer         (kind=HID_T       )                                            :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                            memorySpaceID     , storedDatasetID
    logical                                                                        :: isReference       , readSubsection
    type            (hdf5Object       )                                            :: datasetObject
    type            (varying_string   )                                            :: datasetNameActual , message
    type            (c_ptr            )                                            :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 3D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,3)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_3D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Double_3D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a double 3-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    double precision                   , allocatable, dimension(:,:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                               , intent(inout)           :: self
    character       (len=*            )                               , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(3)    , intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(3)                              :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                         referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                    :: referencedRegion
    integer                                                                                     :: errorCode
    integer         (kind=HID_T       )                                                         :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                         memorySpaceID     , storedDatasetID
    logical                                                                                     :: isReference       , readSubsection
    type            (hdf5Object       )                                                         :: datasetObject
    type            (varying_string   )                                                         :: datasetNameActual , message
    type            (c_ptr            )                                                         :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 3D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,3)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(3,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_3D_Array_Allocatable

  subroutine IO_HDF5_Write_Dataset_Double_4D(self,datasetValue,datasetName,comment,appendTo,appendDimension,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a double 4-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_DOUBLE          , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class           (hdf5Object    )                    , intent(inout)           :: self
    character       (len=*         )                    , intent(in   ), optional :: comment                     , datasetName
    double precision                , dimension(:,:,:,:), intent(in   )           :: datasetValue
    logical                                             , intent(in   ), optional :: appendTo
    integer         (hsize_t       )                    , intent(in   ), optional :: chunkSize
    integer                                             , intent(in   ), optional :: appendDimension             , compressionLevel
    type            (hdf5Object    )                    , intent(  out), optional :: datasetReturned
    integer         (kind=HSIZE_T  ), dimension(4)                                :: datasetDimensions           , hyperslabCount             , &
         &                                                                           hyperslabStart              , newDatasetDimensions       , &
         &                                                                           newDatasetDimensionsFiltered, newDatasetDimensionsMaximum
    integer                                                                       :: appendDimensionActual       , datasetRank                , &
         &                                                                           errorCode
    integer         (kind=HID_T    )                                              :: dataspaceID                 , newDataspaceID
    logical                                                                       :: appendToActual              , preExisted
    type            (hdf5Object    )                                              :: datasetObject
    type            (varying_string)                                              :: datasetNameActual           , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 4D double.
          call self%assertDatasetType(H5T_NATIVE_DOUBLES,4)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeDouble,datasetDimensions,appendTo&
            &=appendTo,appendDimension=appendDimension,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 4D double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,4)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Determine the dimension for appending.
       appendDimensionActual=1
       if (present(appendDimension)) appendDimensionActual=appendDimension
       ! Ensure that all dimensions other than the one being appended to are of the same size.
       newDatasetDimensionsFiltered                       =newDatasetDimensions
       newDatasetDimensionsFiltered(appendDimensionActual)=dataSetDimensions   (appendDimensionActual)
       if (any(dataSetDimensions /= newDatasetDimensionsFiltered)) then
          message="when appending to dataset '"//trim(datasetNameActual)//"' all dimensions other than that being appended to must be same as original dataset"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Set the hyperslab.
       hyperslabStart                             =0
       hyperslabStart      (appendDimensionActual)=newDatasetDimensions(appendDimensionActual)
       hyperslabCount                             =dataSetDimensions
       newDatasetDimensions(appendDimensionActual)=newDatasetDimensions(appendDimensionActual)+datasetDimensions(appendDimensionActual)
    else
       newDatasetDimensions   =datasetDimensions
       hyperslabStart         =0
       hyperslabCount         =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=4
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Double_4D

  subroutine IO_HDF5_Read_Dataset_Double_4D_Array_Static(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a double scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    double precision                   , dimension(:,:,:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                    , intent(inout)           :: self
    character       (len=*            )                    , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     ), dimension(4)      , intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     ), dimension(4)                                :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                              referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save              , target                  :: referencedRegion
    integer                                                                          :: errorCode
    integer         (kind=HID_T       )                                              :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                              memorySpaceID     , storedDatasetID
    logical                                                                          :: isReference       , readSubsection
    type            (hdf5Object       )                                              :: datasetObject
    type            (varying_string   )                                              :: datasetNameActual , message
    type            (c_ptr            )                                              :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 4D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,4)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(4,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(4,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(4,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_4D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Double_4D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a double 4-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    double precision                   , allocatable, dimension(:,:,:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                                 , intent(inout)           :: self
    character       (len=*            )                                 , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(4)      , intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(4)                                :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                           referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                      :: referencedRegion
    integer                                                                                       :: errorCode
    integer         (kind=HID_T       )                                                           :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                           memorySpaceID     , storedDatasetID
    logical                                                                                       :: isReference       , readSubsection
    type            (hdf5Object       )                                                           :: datasetObject
    type            (varying_string   )                                                           :: datasetNameActual , message
    type            (c_ptr            )                                                           :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 4D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,4)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(4,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(4,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(4,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_4D_Array_Allocatable

  subroutine IO_HDF5_Write_Dataset_Double_5D(self,datasetValue,datasetName,comment,appendTo,appendDimension,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a double 5-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_DOUBLE          , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class           (hdf5Object    )                      , intent(inout)           :: self
    character       (len=*         )                      , intent(in   ), optional :: comment                     , datasetName
    double precision                , dimension(:,:,:,:,:), intent(in   )           :: datasetValue
    logical                                               , intent(in   ), optional :: appendTo
    integer         (hsize_t       )                      , intent(in   ), optional :: chunkSize
    integer                                               , intent(in   ), optional :: appendDimension             , compressionLevel
    type            (hdf5Object    )                      , intent(  out), optional :: datasetReturned
    integer         (kind=HSIZE_T  ), dimension(5)                                  :: datasetDimensions           , hyperslabCount             , &
         &                                                                             hyperslabStart              , newDatasetDimensions       , &
         &                                                                             newDatasetDimensionsFiltered, newDatasetDimensionsMaximum
    integer                                                                         :: appendDimensionActual       , datasetRank                , &
         &                                                                             errorCode
    integer         (kind=HID_T    )                                                :: dataspaceID                 , newDataspaceID
    logical                                                                         :: appendToActual              , preExisted
    type            (hdf5Object    )                                                :: datasetObject
    type            (varying_string)                                                :: datasetNameActual           , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 5D double.
          call self%assertDatasetType(H5T_NATIVE_DOUBLES,5)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeDouble,datasetDimensions,appendTo&
            &=appendTo,appendDimension=appendDimension,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 5D double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,5)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Determine the dimension for appending.
       appendDimensionActual=1
       if (present(appendDimension)) appendDimensionActual=appendDimension
       ! Ensure that all dimensions other than the one being appended to are of the same size.
       newDatasetDimensionsFiltered                       =newDatasetDimensions
       newDatasetDimensionsFiltered(appendDimensionActual)=dataSetDimensions   (appendDimensionActual)
       if (any(dataSetDimensions /= newDatasetDimensionsFiltered)) then
          message="when appending to dataset '"//trim(datasetNameActual)//"' all dimensions other than that being appended to must be same as original dataset"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Set the hyperslab.
       hyperslabStart                             =0
       hyperslabStart      (appendDimensionActual)=newDatasetDimensions(appendDimensionActual)
       hyperslabCount                             =dataSetDimensions
       newDatasetDimensions(appendDimensionActual)=newDatasetDimensions(appendDimensionActual)+datasetDimensions(appendDimensionActual)
    else
       newDatasetDimensions   =datasetDimensions
       hyperslabStart         =0
       hyperslabCount         =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=5
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Double_5D

  subroutine IO_HDF5_Read_Dataset_Double_5D_Array_Static(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a double scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    double precision                   , dimension(:,:,:,:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                      , intent(inout)           :: self
    character       (len=*            )                      , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     ), dimension(5)        , intent(in   ), optional :: readBegin         , readCount                    !   &&
    integer         (kind=HSIZE_T     ), dimension(5)                                  :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save                , target                  :: referencedRegion
    integer                                                                            :: errorCode
    integer         (kind=HID_T       )                                                :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                memorySpaceID     , storedDatasetID
    logical                                                                            :: isReference       , readSubsection
    type            (hdf5Object       )                                                :: datasetObject
    type            (varying_string   )                                                :: datasetNameActual , message
    type            (c_ptr            )                                                :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 5D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,5)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(5,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(5,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(5,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_5D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Double_5D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a double 5-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    double precision                   , allocatable, dimension(:,:,:,:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                                   , intent(inout)           :: self
    character       (len=*            )                                   , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(5)        , intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(5)                                  :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                             referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                        :: referencedRegion
    integer                                                                                         :: errorCode
    integer         (kind=HID_T       )                                                             :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                             memorySpaceID     , storedDatasetID
    logical                                                                                         :: isReference       , readSubsection
    type            (hdf5Object       )                                                             :: datasetObject
    type            (varying_string   )                                                             :: datasetNameActual , message
    type            (c_ptr            )                                                             :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 5D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,5)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(5,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(5,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(5,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_5D_Array_Allocatable

  subroutine IO_HDF5_Write_Dataset_Double_6D(self,datasetValue,datasetName,comment,appendTo,appendDimension,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a double 6-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_DOUBLE          , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, hsize_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class           (hdf5Object    )                        , intent(inout)           :: self
    character       (len=*         )                        , intent(in   ), optional :: comment                     , datasetName
    double precision                , dimension(:,:,:,:,:,:), intent(in   )           :: datasetValue
    logical                                                 , intent(in   ), optional :: appendTo
    integer         (hsize_t       )                        , intent(in   ), optional :: chunkSize
    integer                                                 , intent(in   ), optional :: appendDimension             , compressionLevel
    type            (hdf5Object    )                        , intent(  out), optional :: datasetReturned
    integer         (kind=HSIZE_T  ), dimension(6)                                    :: datasetDimensions           , hyperslabCount             , &
         &                                                                               hyperslabStart              , newDatasetDimensions       , &
         &                                                                               newDatasetDimensionsFiltered, newDatasetDimensionsMaximum
    integer                                                                           :: appendDimensionActual       , datasetRank                , &
         &                                                                               errorCode
    integer         (kind=HID_T    )                                                  :: dataspaceID                 , newDataspaceID
    logical                                                                           :: appendToActual              , preExisted
    type            (hdf5Object    )                                                  :: datasetObject
    type            (varying_string)                                                  :: datasetNameActual           , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 6D double.
          call self%assertDatasetType(H5T_NATIVE_DOUBLES,6)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeDouble,datasetDimensions,appendTo&
            &=appendTo,appendDimension=appendDimension,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 6D double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,6)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Determine the dimension for appending.
       appendDimensionActual=1
       if (present(appendDimension)) appendDimensionActual=appendDimension
       ! Ensure that all dimensions other than the one being appended to are of the same size.
       newDatasetDimensionsFiltered                       =newDatasetDimensions
       newDatasetDimensionsFiltered(appendDimensionActual)=dataSetDimensions   (appendDimensionActual)
       if (any(dataSetDimensions /= newDatasetDimensionsFiltered)) then
          message="when appending to dataset '"//trim(datasetNameActual)//"' all dimensions other than that being appended to must be same as original dataset"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Set the hyperslab.
       hyperslabStart                             =0
       hyperslabStart      (appendDimensionActual)=newDatasetDimensions(appendDimensionActual)
       hyperslabCount                             =dataSetDimensions
       newDatasetDimensions(appendDimensionActual)=newDatasetDimensions(appendDimensionActual)+datasetDimensions(appendDimensionActual)
    else
       newDatasetDimensions   =datasetDimensions
       hyperslabStart         =0
       hyperslabCount         =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=6
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Double_6D

  subroutine IO_HDF5_Read_Dataset_Double_6D_Array_Static(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a double scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    double precision                   , dimension(:,:,:,:,:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                        , intent(inout)           :: self
    character       (len=*            )                        , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     ), dimension(6)          , intent(in   ), optional :: readBegin         , readCount                    !   &&
    integer         (kind=HSIZE_T     ), dimension(6)                                    :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                  referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save                  , target                  :: referencedRegion
    integer                                                                              :: errorCode
    integer         (kind=HID_T       )                                                  :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                  memorySpaceID     , storedDatasetID
    logical                                                                              :: isReference       , readSubsection
    type            (hdf5Object       )                                                  :: datasetObject
    type            (varying_string   )                                                  :: datasetNameActual , message
    type            (c_ptr            )                                                  :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 5D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,6)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(6,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(6,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(6,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_6D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Double_6D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a double 6-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F         , H5S_SELECT_SET_F      , H5T_NATIVE_DOUBLE          , &
          &                                      H5T_STD_REF_DSETREG  , HID_T             , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f       , h5dread_f         , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f           , h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , operator(//)      , trim
    implicit none
    double precision                   , allocatable, dimension(:,:,:,:,:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                                     , intent(inout)           :: self
    character       (len=*            )                                     , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(6)          , intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(6)                                    :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                               referenceEnd      , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                          :: referencedRegion
    integer                                                                                           :: errorCode
    integer         (kind=HID_T       )                                                               :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                               memorySpaceID     , storedDatasetID
    logical                                                                                           :: isReference       , readSubsection
    type            (hdf5Object       )                                                               :: datasetObject
    type            (varying_string   )                                                               :: datasetNameActual , message
    type            (c_ptr            )                                                               :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 6D double array.
    call datasetObject%assertDatasetType(H5T_NATIVE_DOUBLES,6)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(6,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(6,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(6,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,H5T_NATIVE_DOUBLE,datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Double_6D_Array_Allocatable

  subroutine IO_HDF5_Write_Dataset_Character_1D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a character 1-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5S_SELECT_SET_F  , H5T_NATIVE_CHARACTER       , HID_T                , HSIZE_T   , &
          &                           h5dget_space_f    , h5dset_extent_f            , h5dwrite_f           , h5sclose_f, &
          &                           h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, h5tclose_f, &
          &                           h5tcopy_f         , h5tset_size_f              , hsize_t              , size_t
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)               , trim
    implicit none
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: comment                    , datasetName
    character(len=*         ), dimension(:), intent(in   )           :: datasetValue
    logical                                , intent(in   ), optional :: appendTo
    integer  (hsize_t       )              , intent(in   ), optional :: chunkSize
    integer                                , intent(in   ), optional :: compressionLevel
    type     (hdf5Object    )              , intent(  out), optional :: datasetReturned
    integer  (kind=HSIZE_T  ), dimension(1)                          :: datasetDimensions          , hyperslabCount      , &
         &                                                              hyperslabStart             , newDatasetDimensions, &
         &                                                              newDatasetDimensionsMaximum
    integer                                                          :: datasetRank                , errorCode
    integer  (kind=HID_T    )                                        :: dataTypeID                 , dataspaceID         , &
         &                                                              newDataspaceID
    logical                                                          :: appendToActual             , preExisted
    type     (hdf5Object    )                                        :: datasetObject
    type     (varying_string)                                        :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Create a custom datatype.
    call h5tcopy_f(H5T_NATIVE_CHARACTER,dataTypeID,errorCode)
    if (errorCode < 0) then
       message="unable to make custom datatype for attribute '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tset_size_f(dataTypeID,int(len(datasetValue),size_t),errorCode)
    if (errorCode < 0) then
       message="unable to set datatype size for attribute '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is a dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D character.
          call self%assertDatasetType([dataTypeID],1)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject    =self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeCharacter,datasetDimensions,useDataType&
            &=dataTypeID,appendTo =appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 1D integer.
       if (preExisted) call datasetObject%assertDatasetType([dataTypeID],1)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=1
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    call h5dwrite_f(datasetObject%objectID,dataTypeID,datasetValue,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID,errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_Character_1D

  subroutine IO_HDF5_Write_Dataset_VarString_1D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a varying string 1-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: HDF5           , only : hsize_t
    use :: String_Handling, only : Convert_VarString_To_Char
    implicit none
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: comment        , datasetName
    type     (varying_string), dimension(:), intent(in   )           :: datasetValue
    logical                                , intent(in   ), optional :: appendTo
    integer  (hsize_t       )              , intent(in   ), optional :: chunkSize
    integer                                , intent(in   ), optional :: compressionLevel
    type     (hdf5Object    )              , intent(  out), optional :: datasetReturned

    ! Call the character version of this routine to perform the write.
    call IO_HDF5_Write_Dataset_Character_1D(self,Convert_VarString_To_Char(datasetValue),datasetName,comment,appendTo&
         &,chunkSize,compressionLevel,datasetReturned)

    return
  end subroutine IO_HDF5_Write_Dataset_VarString_1D

  subroutine IO_HDF5_Read_Dataset_Character_1D_Array_Static(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read a character scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F     , H5S_ALL_F             , H5S_SELECT_SET_F           , H5T_STD_REF_DSETREG  , &
          &                                      HID_T             , HSIZE_T               , h5dclose_f                 , h5dget_space_f       , &
          &                                      h5dread_f         , h5rdereference_f      , h5rget_region_f            , h5sclose_f           , &
          &                                      h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, &
          &                                      h5tclose_f        , hdset_reg_ref_t_f     , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)     , operator(//)          , trim
    implicit none
    character(len=*            ), dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object       )              , intent(inout)           :: self
    character(len=*            )              , intent(in   ), optional :: datasetName
    integer  (kind=HSIZE_T     ), dimension(1), intent(in   ), optional :: readBegin              , readCount
    integer  (kind=HSIZE_T     ), dimension(1)                          :: datasetDimensions      , datasetMaximumDimensions, &
         &                                                                 referenceEnd           , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (hdset_reg_ref_t_f), save        , target                  :: referencedRegion
    integer                                                             :: errorCode
    integer  (kind=HID_T       )                                        :: dataTypeID          (6), datasetDataspaceID      , &
         &                                                                 dereferencedObjectID   , memorySpaceID           , &
         &                                                                 storedDatasetID
    logical                                                             :: isReference            , readSubsection
    type     (hdf5Object       )                                        :: datasetObject
    type     (varying_string   )                                        :: datasetNameActual      , message
    type     (c_ptr            )                                        :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Create a custom datatype.
    dataTypeID=IO_HDF5_Character_Types(len(datasetValue))

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
   else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D integer array.
    call datasetObject%assertDatasetType(dataTypeID,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Ensure that the size of the array is large enough to hold the datasets.
    if (any(shape(datasetValue) < datasetDimensions)) then
       message="array is not large enough to hold datasets from '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,dataTypeID(1),datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID(1),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tclose_f(dataTypeID(2),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Character_1D_Array_Static

  subroutine IO_HDF5_Read_Dataset_Character_1D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount)
    !!{
    Open and read an integer scalar dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F     , H5S_ALL_F             , H5S_SELECT_SET_F           , H5T_STD_REF_DSETREG  , &
          &                                      HID_T             , HSIZE_T               , h5dclose_f                 , h5dget_space_f       , &
          &                                      h5dread_f         , h5rdereference_f      , h5rget_region_f            , h5sclose_f           , &
          &                                      h5screate_simple_f, h5sget_select_bounds_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f, &
          &                                      h5tclose_f        , hdset_reg_ref_t_f     , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)     , operator(//)          , trim
    implicit none
    character(len=*            ), allocatable, dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object       )                           , intent(inout)           :: self
    character(len=*            )                           , intent(in   ), optional :: datasetName
    integer  (kind=HSIZE_T     )             , dimension(1), intent(in   ), optional :: readBegin              , readCount
    integer  (kind=HSIZE_T     )             , dimension(1)                          :: datasetDimensions      , datasetMaximumDimensions, &
         &                                                                              referenceEnd           , referenceStart
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type     (hdset_reg_ref_t_f), save       , target                                :: referencedRegion
    integer                                                                          :: errorCode
    integer  (kind=HID_T       )                                                     :: dataTypeID          (6), datasetDataspaceID      , &
         &                                                                              dereferencedObjectID   , memorySpaceID           , &
         &                                                                              storedDatasetID
    logical                                                                          :: isReference            , readSubsection
    type     (hdf5Object       )                                                     :: datasetObject
    type     (varying_string   )                                                     :: datasetNameActual      , message
    type     (c_ptr            )                                                     :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Create a custom datatype.
    dataTypeID=IO_HDF5_Character_Types(len(datasetValue))

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
   else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D integer array.
    call datasetObject%assertDatasetType(dataTypeID,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    !![
    <allocate variable="datasetValue" shape="datasetDimensions"/>
    !!]
    ! Read the dataset.
    call h5dread_f(datasetObject%objectID,dataTypeID(1),datasetValue,int(shape(datasetValue),kind=hsize_t),errorCode&
         &,memorySpaceID,datasetDataspaceID)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID(1),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tclose_f(dataTypeID(2),errorCode)
    if (errorCode < 0) then
       message="unable to close custom datatype for attribute '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_Character_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Dataset_VarString_1D_Array_Allocatable(self,datasetName,datasetValue)
    !!{
    Open and read an varying string 1-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , h5dget_type_f, h5tclose_f, h5tget_size_f
    use :: ISO_Varying_String, only : assignment(=), operator(//) , trim
    implicit none
    type     (varying_string), allocatable, dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   ), optional :: datasetName
    integer  (kind=HID_T    )                                                     :: dataTypeID
    integer  (kind=SIZE_T   )                                                     :: dataTypeSize
    integer                                                                       :: errorCode
    type     (hdf5Object    )                                                     :: datasetObject
    type     (varying_string)                                                     :: datasetNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Get the datatype of this dataset.
    call h5dget_type_f(datasetObject%objectID,dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not get datatype for '"//trim(datasetNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Get the size of the datatype.
    call h5tget_size_f(dataTypeID,dataTypeSize,errorCode)
    if (errorCode /= 0) then
       message="can not get size of datatype for '"//trim(datasetNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not close datatype of '"//trim(datasetNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Call wrapper routine that will do the remainder of the read.
    call IO_HDF5_Read_Dataset_VarString_1D_Array_Allocatable_Do_Read(self,datasetName,datasetValue,dataTypeSize)
    return
  end subroutine IO_HDF5_Read_Dataset_VarString_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Dataset_VarString_1D_Array_Allocatable_Do_Read(self,datasetName,datasetValue,dataTypeSize)
    !!{
    Open and read an varying string 1-D array dataset in {\normalfont \ttfamily self} by creating a suitably-sized character variable into
    which it can be read.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type     (varying_string  ), allocatable, dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object      )                           , intent(inout)           :: self
    character(len=*           )                           , intent(in   ), optional :: datasetName
    integer  (kind=SIZE_T     )                           , intent(in   )           :: dataTypeSize
    character(len=dataTypeSize), allocatable, dimension(:)                          :: temporaryBuffer

    ! Call the character version of this routine to perform the red.
    call IO_HDF5_Read_Dataset_Character_1D_Array_Allocatable(self,datasetName,temporaryBuffer)

    ! Transfer the results to the varying string variable.
    allocate(datasetValue(size(temporaryBuffer)))
    datasetValue=temporaryBuffer
    deallocate(temporaryBuffer)

    return
  end subroutine IO_HDF5_Read_Dataset_VarString_1D_Array_Allocatable_Do_Read

  subroutine IO_HDF5_Read_Dataset_VarString_1D_Array_Static(self,datasetName,datasetValue)
    !!{
    Open and read an varying string 1-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : HID_T        , h5dget_type_f, h5tclose_f, h5tget_size_f
    use :: ISO_Varying_String, only : assignment(=), operator(//) , trim
    implicit none
    type     (varying_string), dimension(:), intent(inout)           :: datasetValue
    class    (hdf5Object    )              , intent(inout)           :: self
    character(len=*         )              , intent(in   ), optional :: datasetName
    integer  (kind=HID_T    )                                        :: dataTypeID
    integer  (kind=SIZE_T   )                                        :: dataTypeSize
    integer                                                          :: errorCode
    type     (hdf5Object    )                                        :: datasetObject
    type     (varying_string)                                        :: datasetNameActual, message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetName)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Get the datatype of this dataset.
    call h5dget_type_f(datasetObject%objectID,dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not get datatype for '"//trim(datasetNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Get the size of the datatype.
    call h5tget_size_f(dataTypeID,dataTypeSize,errorCode)
    if (errorCode /= 0) then
       message="can not get size of datatype for '"//trim(datasetNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the datatype.
    call h5tclose_f(dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="can not close datatype of '"//trim(datasetNameActual)//"' located in '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Call wrapper routine that will do the remainder of the read.
    call IO_HDF5_Read_Dataset_VarString_1D_Array_Static_Do_Read(self,datasetName,datasetValue,dataTypeSize)
    return
  end subroutine IO_HDF5_Read_Dataset_VarString_1D_Array_Static

  subroutine IO_HDF5_Read_Dataset_VarString_1D_Array_Static_Do_Read(self,datasetName,datasetValue,dataTypeSize)
    !!{
    Open and read an varying string 1-D array dataset in {\normalfont \ttfamily self} by creating a suitably-sized character variable into
    which it can be read.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type     (varying_string  ), dimension(:)                 , intent(inout)           :: datasetValue
    class    (hdf5Object      )                               , intent(inout)           :: self
    character(len=*           )                               , intent(in   ), optional :: datasetName
    integer  (kind=SIZE_T     )                               , intent(in   )           :: dataTypeSize
    character(len=dataTypeSize), dimension(size(datasetValue))                          :: temporaryBuffer

    ! Call the character version of this routine to perform the red.
    call IO_HDF5_Read_Dataset_Character_1D_Array_Static(self,datasetName,temporaryBuffer)

    ! Transfer the results to the varying string variable.
    datasetValue=temporaryBuffer

    return
  end subroutine IO_HDF5_Read_Dataset_VarString_1D_Array_Static_Do_Read

  subroutine IO_HDF5_Read_Dataset_VarDouble_1D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a varying-length 1D double dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F       , H5S_ALL_F            , H5S_SELECT_SET_F      , size_t                     , &
          &                                      H5T_STD_REF_DSETREG , HID_T                , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f      , h5dread_f            , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f          , h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_elements_f, h5sselect_hyperslab_f, hdset_reg_ref_t_f     , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_f_pointer         , c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)         , trim
    implicit none
    type            (hdf5VarDouble    ), allocatable, dimension(:  ), intent(  out)           :: datasetValue
    class           (hdf5Object       )                             , intent(inout)           :: self
    character       (len=*            )                             , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(1  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     )             , dimension(1  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                       referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), allocatable, dimension(:,:)                          :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                  :: referencedRegion
    integer                                                                                   :: errorCode
    integer         (kind=HSIZE_T     )                                                       :: i
    integer         (kind=HID_T       )                                                       :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                       memorySpaceID     , storedDatasetID
    logical                                                                                   :: isReference       , readSubsection
    type            (hdf5Object       )                                                       :: datasetObject
    type            (varying_string   )                                                       :: datasetNameActual , message
    type            (c_ptr            )                                                       :: dataBuffer
    type            (hdf5VlenC        ), allocatable, dimension(:  ), target                  :: datasetValueC

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=IO_HDF5_Open_Dataset(self,datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D double array.
    call datasetObject%assertDatasetType(H5T_VLEN_DOUBLE,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    allocate(datasetValue (datasetDimensions(1)))
    allocate(datasetValueC(datasetDimensions(1)))
    ! Read the dataset.
    dataBuffer=c_loc(datasetValueC)
    errorCode=h5dread(datasetObject%objectID,H5T_VLEN_DOUBLE(1),memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    do i=1,datasetDimensions(1)
       call c_f_pointer(datasetValueC(i)%p,datasetValue(i)%row,shape=[datasetValueC(i)%length])
    end do
    deallocate(datasetValueC)

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Input was a dataset object. Test if it was a reference.
       if (datasetObject%isReference()) then
          ! It was, so close the referenced dataset.
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Restore the object ID of the original dataset.
          self%objectID=storedDatasetID
       end if
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_VarDouble_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Dataset_VarVarDouble_1D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a varying-length $\times$ varying-length 1D double dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F       , H5S_ALL_F            , H5S_SELECT_SET_F      , size_t                     , &
          &                                      H5T_STD_REF_DSETREG , HID_T                , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f      , h5dread_f            , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f          , h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_elements_f, h5sselect_hyperslab_f, hdset_reg_ref_t_f     , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_f_pointer         , c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)         , trim
    implicit none
    type            (hdf5VarDouble2D  ), allocatable, dimension(:  ), intent(  out)           :: datasetValue
    class           (hdf5Object       )                             , intent(inout)           :: self
    character       (len=*            )                             , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(1  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     )             , dimension(1  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                       referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), allocatable, dimension(:,:)                          :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works. 
    type            (hdset_reg_ref_t_f), save       , target                                  :: referencedRegion
    integer                                                                                   :: errorCode
    integer         (kind=HSIZE_T     )                                                       :: i                 , j
    integer         (kind=HID_T       )                                                       :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                       memorySpaceID     , storedDatasetID
    logical                                                                                   :: isReference       , readSubsection
    type            (hdf5Object       )                                                       :: datasetObject
    type            (varying_string   )                                                       :: datasetNameActual , message
    type            (c_ptr            )                                                       :: dataBuffer
    type            (hdf5VlenC        ), allocatable, dimension(:  ), target                  :: datasetValueC1
    type            (hdf5VlenC        )             , dimension(:  ), pointer                 :: datasetValueC2
    double precision                                , dimension(:  ), pointer                 :: datasetRow

    
    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D double array.
    call datasetObject%assertDatasetType(H5T_VLEN_VLEN_DOUBLE,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    allocate(datasetValue  (datasetDimensions(1)))
    allocate(datasetValueC1(datasetDimensions(1)))
    ! Read the dataset.
    dataBuffer=c_loc(datasetValueC1)
    errorCode=h5dread(datasetObject%objectID,H5T_VLEN_VLEN_DOUBLE(1),memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    do i=1,datasetDimensions(1)
       call c_f_pointer(datasetValueC1(i)%p,datasetValueC2,shape=[datasetValueC1%length])
       allocate(datasetValue(i)%row(datasetValueC2(1)%length,datasetValueC1(i)%length))
       do j=1,datasetValueC1(i)%length
          call c_f_pointer(datasetValueC2(j)%p,datasetRow,shape=[datasetValueC2(j)%length])
          datasetValue(i)%row(:,j)=datasetRow
          deallocate(datasetRow)
          nullify(datasetRow)
       end do
       deallocate(datasetValueC2)
       nullify(datasetValueC2)
    end do
    deallocate(datasetValueC1)
    
    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType /= hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_VarVarDouble_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Dataset_VarDouble_2D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a varying-length 2D double dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F       , H5S_ALL_F            , H5S_SELECT_SET_F      , size_t                     , &
          &                                      H5T_STD_REF_DSETREG , HID_T                , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f      , h5dread_f            , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f          , h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_elements_f, h5sselect_hyperslab_f, hdset_reg_ref_t_f     , hsize_t
    use, intrinsic :: ISO_C_Binding     , only : c_f_pointer         , c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)         , trim
    implicit none
    type            (hdf5VarDouble    ), allocatable, dimension(:,:), intent(  out)           :: datasetValue
    class           (hdf5Object       )                             , intent(inout)           :: self
    character       (len=*            )                             , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(2  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     )             , dimension(2  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                       referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), allocatable, dimension(:,:)                          :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                  :: referencedRegion
    integer                                                                                   :: errorCode
    integer         (kind=HSIZE_T     )                                                       :: i                 , j
    integer         (kind=HID_T       )                                                       :: datasetDataspaceID, dereferencedObjectID, &
         &                                                                                       memorySpaceID     , storedDatasetID
    logical                                                                                   :: isReference       , readSubsection
    type            (hdf5Object       )                                                       :: datasetObject
    type            (varying_string   )                                                       :: datasetNameActual , message
    type            (c_ptr            )                                                       :: dataBuffer
    type            (hdf5VlenC        ), allocatable, dimension(:,:), target                  :: datasetValueC

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a 1D double array.
    call datasetObject%assertDatasetType(H5T_VLEN_DOUBLE,2)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(2,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    allocate(datasetValue (datasetDimensions(1),datasetDimensions(2)))
    allocate(datasetValueC(datasetDimensions(1),datasetDimensions(2)))
    ! Read the dataset.
    dataBuffer=c_loc(datasetValueC)
    errorCode=h5dread(datasetObject%objectID,H5T_VLEN_DOUBLE(1),memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    do i=1,datasetDimensions(1)
       do j=1,datasetDimensions(2)
          call c_f_pointer(datasetValueC(i,j)%p,datasetValue(i,j)%row,shape=[datasetValueC(i,j)%length])
       end do
    end do
    deallocate(datasetValueC)

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_VarDouble_2D_Array_Allocatable

  subroutine IO_HDF5_Write_Dataset_VarDouble_1D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a varying-length double 1-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F, H5S_SELECT_SET_F  , hsize_t                    , HID_T                , &
          &                                      HSIZE_T      , h5dget_space_f    , h5dset_extent_f            , h5dwrite_vl_f        , &
          &                                      h5sclose_f   , h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f
    use            :: ISO_Varying_String, only : assignment(=), operator(//)      , trim
    implicit none
    class           (hdf5Object    ), intent(inout)                       :: self
    character       (len=*         ), intent(in   ), optional             :: comment                    , datasetName
    type            (hdf5VarDouble ), intent(in   ), dimension(:)         :: datasetValue
    logical                         , intent(in   ), optional             :: appendTo
    integer         (hsize_t       ), intent(in   ), optional             :: chunkSize
    integer                         , intent(in   ), optional             :: compressionLevel
    type            (hdf5Object    ), intent(  out), optional             :: datasetReturned
    integer         (kind=HSIZE_T  )               , dimension(1)         :: datasetDimensions          , hyperslabCount      , &
         &                                                                   hyperslabStart             , newDatasetDimensions, &
         &                                                                   newDatasetDimensionsMaximum
    type            (hdf5VlenC     ), allocatable  , dimension(:), target :: datasetValueC
    type            (c_ptr         )                                      :: datasetValueC_
    integer         (kind=HSIZE_T  )                                      :: i
    integer                                                               :: datasetRank                , errorCode
    integer         (kind=HID_T    )                                      :: dataspaceID                , newDataspaceID
    logical                                                               :: appendToActual             , preExisted
    type            (hdf5Object    )                                      :: datasetObject
    type            (varying_string)                                      :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized()

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=     datasetName
    else
       datasetNameActual=self% objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is a dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D vlen double.
          call self%assertDatasetType(H5T_VLEN_DOUBLE,1)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=IO_HDF5_Open_Dataset(self,datasetName,comment,hdf5DataTypeVlenDouble,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 1D double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_VLEN_DOUBLE,1)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=1
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    allocate(datasetValueC(datasetDimensions(1)))
    do i=1,datasetDimensions(1)
       datasetValueC(i)%length=size (datasetValue(i)%row,kind=c_size_t)
       datasetValueC(i)%p     =c_loc(datasetValue(i)%row              )
    end do
    datasetValueC_=c_loc(datasetValueC)
    errorCode     =H5Dwrite(datasetObject%objectID,H5T_VLEN_DOUBLE(1),newDataspaceID,dataspaceID,H5P_DEFAULT_F,datasetValueC_)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    deallocate(datasetValueC)

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_VarDouble_1D
  
  subroutine IO_HDF5_Write_Dataset_VarVarDouble_1D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a varying-length $\times$ varying-length 1D double array dataset in {\normalfont \ttfamily self}.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F, H5S_SELECT_SET_F  , hsize_t                    , HID_T                , &
          &                                      HSIZE_T      , h5dget_space_f    , h5dset_extent_f            , h5dwrite_vl_f        , &
          &                                      h5sclose_f   , h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f
    use            :: ISO_Varying_String, only : assignment(=), operator(//)      , trim
    implicit none
    class           (hdf5Object     ), intent(inout)                       :: self
    character       (len=*          ), intent(in   ), optional             :: comment                    , datasetName
    type            (hdf5VarDouble2D), intent(in   ), dimension(:)         :: datasetValue
    logical                          , intent(in   ), optional             :: appendTo
    integer         (hsize_t        ), intent(in   ), optional             :: chunkSize
    integer                          , intent(in   ), optional             :: compressionLevel
    type            (hdf5Object     ), intent(  out), optional             :: datasetReturned
    integer         (kind=HSIZE_T   )               , dimension(1)         :: datasetDimensions          , hyperslabCount      , &
         &                                                                    hyperslabStart             , newDatasetDimensions, &
         &                                                                    newDatasetDimensionsMaximum
    type            (hdf5VlenC      ), allocatable  , dimension(:), target :: datasetValueC1
    type            (hdf5VlenVlenC  ), allocatable  , dimension(:), target :: datasetValueC2
    type            (c_ptr          )                                      :: datasetValueC_
    integer         (kind=HSIZE_T   )                                      :: i                          , j
    integer                                                                :: datasetRank                , errorCode
    integer         (kind=HID_T     )                                      :: dataspaceID                , newDataspaceID
    logical                                                                :: appendToActual             , preExisted
    type            (hdf5Object     )                                      :: datasetObject
    type            (varying_string )                                      :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized()

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=     datasetName
    else
       datasetNameActual=self% objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is a dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D vlen-vlen double.
          call self%assertDatasetType(H5T_VLEN_VLEN_DOUBLE,1)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=IO_HDF5_Open_Dataset(self,datasetName,comment,hdf5DataTypeVlenVlenDouble,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 1D vlen-vlen double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_VLEN_VLEN_DOUBLE,1)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=1
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    allocate(datasetValueC1(datasetDimensions(1)))
    allocate(datasetValueC2(datasetDimensions(1)))
    do i=1,datasetDimensions(1)
       allocate(datasetValueC2(i)%row(size(datasetValue(i)%row,dim=2,kind=c_size_t)))
       datasetValueC1(i)%length=size (datasetValue  (i)%row,dim=2,kind=c_size_t)
       datasetValueC1(i)%p     =c_loc(datasetValueC2(i)%row                    )
       do j=1,size (datasetValue(i)%row,dim=2,kind=c_size_t)
          datasetValueC2(i)%row(j)%length=size (datasetValue(i)%row     ,dim=1,kind=c_size_t)
          datasetValueC2(i)%row(j)%p     =c_loc(datasetValue(i)%row(1,j)                    )
       end do
    end do
    datasetValueC_=c_loc(datasetValueC1)
    errorCode     =H5Dwrite(datasetObject%objectID,H5T_VLEN_VLEN_DOUBLE(1),newDataspaceID,dataspaceID,H5P_DEFAULT_F,datasetValueC_)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    deallocate(datasetValueC1)
    deallocate(datasetValueC2)

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_VarVarDouble_1D

  subroutine IO_HDF5_Write_Dataset_VarDouble_2D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a varying-length double 2-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F, H5S_SELECT_SET_F  , hsize_t                    , HID_T                , &
          &                                      HSIZE_T      , h5dget_space_f    , h5dset_extent_f            , h5dwrite_vl_f        , &
          &                                      h5sclose_f   , h5screate_simple_f, h5sget_simple_extent_dims_f, h5sselect_hyperslab_f
    use            :: ISO_Varying_String, only : assignment(=), operator(//)      , trim
    implicit none
    class           (hdf5Object    ), intent(inout)                         :: self
    character       (len=*         ), intent(in   ), optional               :: comment                       , datasetName
    type            (hdf5VarDouble ), intent(in   ), dimension(:,:)         :: datasetValue
    logical                         , intent(in   ), optional               :: appendTo
    integer         (hsize_t       ), intent(in   ), optional               :: chunkSize
    integer                         , intent(in   ), optional               :: compressionLevel
    type            (hdf5Object    ), intent(  out), optional               :: datasetReturned
    integer         (kind=HSIZE_T  )               , dimension(2  )         :: datasetDimensions          , hyperslabCount      , &
         &                                                                   hyperslabStart             , newDatasetDimensions, &
         &                                                                   newDatasetDimensionsMaximum
    type            (hdf5VlenC     ), allocatable  , dimension(:,:), target :: datasetValueC
    type            (c_ptr         )                                        :: datasetValueC_
    integer         (kind=HSIZE_T  )                                        :: i                          , j
    integer                                                                 :: datasetRank                , errorCode
    integer         (kind=HID_T    )                                        :: dataspaceID                , newDataspaceID
    logical                                                                 :: appendToActual             , preExisted
    type            (hdf5Object    )                                        :: datasetObject
    type            (varying_string)                                        :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized()

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=     datasetName
    else
       datasetNameActual=self% objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is a dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D vlen double.
          call self%assertDatasetType(H5T_VLEN_DOUBLE,2)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=IO_HDF5_Open_Dataset(self,datasetName,comment,hdf5DataTypeVlenDouble,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a 1D double.
       if (preExisted) call datasetObject%assertDatasetType(H5T_VLEN_DOUBLE,2)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=2
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    allocate(datasetValueC(datasetDimensions(1),datasetDimensions(2)))
    do i=1,datasetDimensions(1)
       do j=1,datasetDimensions(2)
          datasetValueC(i,j)%length=size (datasetValue(i,j)%row,kind=c_size_t)
          datasetValueC(i,j)%p     =c_loc(datasetValue(i,j)%row              )
       end do
    end do
    datasetValueC_=c_loc(datasetValueC)
    errorCode     =H5Dwrite(datasetObject%objectID,H5T_VLEN_DOUBLE(1),newDataspaceID,dataspaceID,H5P_DEFAULT_F,datasetValueC_)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    deallocate(datasetValueC)

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_VarDouble_2D
  
  subroutine IO_HDF5_Read_Dataset_VarInteger8_2D_Array_Allocatable(self,datasetName,datasetValue,readBegin,readCount,readSelection)
    !!{
    Open and read a variable-length integer-8 2D array dataset in {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F       , H5S_ALL_F            , H5S_SELECT_SET_F      , hdset_reg_ref_t_f          , &
          &                                      H5T_STD_REF_DSETREG , HID_T                , HSIZE_T               , h5dclose_f                 , &
          &                                      h5dget_space_f      , h5dread_f            , h5rdereference_f      , h5rget_region_f            , &
          &                                      h5sclose_f          , h5screate_simple_f   , h5sget_select_bounds_f, h5sget_simple_extent_dims_f, &
          &                                      h5sselect_elements_f, h5sselect_hyperslab_f, hsize_t               , size_t
    use, intrinsic :: ISO_C_Binding     , only : c_f_pointer         , c_loc
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)         , trim
    implicit none
    type            (hdf5VarInteger8  ), allocatable, dimension(:  ), intent(  out)           :: datasetValue
    class           (hdf5Object       )                             , intent(inout)           :: self
    character       (len=*            )                             , intent(in   ), optional :: datasetName
    integer         (kind=HSIZE_T     )             , dimension(1  ), intent(in   ), optional :: readBegin         , readCount
    integer         (kind=HSIZE_T     )             , dimension(:  ), intent(in   ), optional :: readSelection
    integer         (kind=HSIZE_T     )             , dimension(1  )                          :: datasetDimensions , datasetMaximumDimensions, &
         &                                                                                       referenceEnd      , referenceStart
    integer         (kind=HSIZE_T     ), allocatable, dimension(:,:)                          :: readSelectionMap
    ! <HDF5> Why is "referencedRegion" saved? Because if it is not then it gets dynamically allocated on the stack, which results
    ! in an invalid pointer error. According to valgrind, this happens because the wrong deallocation function is used (delete
    ! instead of delete[] or vice-verse). Presumably this is an HDF5 library error. Saving the variable prevents it from being
    ! deallocated. This is not an elegant solution, but it works.
    type            (hdset_reg_ref_t_f), save       , target                                  :: referencedRegion
    integer                                                                                   :: errorCode
    integer         (kind=HSIZE_T     )                                                       :: i
    integer         (kind=HID_T       )                                                       :: datasetDataspaceID, dereferencedObjectID    , &
         &                                                                                       memorySpaceID     , storedDatasetID
    logical                                                                                   :: isReference       , readSubsection
    type            (hdf5Object       )                                                       :: datasetObject
    type            (varying_string   )                                                       :: datasetNameActual , message
    type            (c_ptr            )                                                       :: dataBuffer
    type            (hdf5VlenC        ), allocatable, dimension(:  ), target                  :: datasetValueC

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=datasetName
    else
       datasetNameActual=self%objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(datasetNameActual)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Only one of a subsection and a selection can be present.
    if (readSubsection.and.present(readSelection)) then
       message="can not specify both a subsection and selection of dataset '"//trim(datasetNameActual)//"' for reading"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Check if the object is an dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! Object is the dataset.
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select

       ! No name should be supplied in this case.
       if (present(datasetName)) then
          message="dataset name was supplied for dataset object '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Require that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="dataset name was not supplied for object '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Check that the dataset exists.
       if (.not.self%hasDataset(datasetName)) then
          message="dataset '"//trim(datasetName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName)
    end if

    ! Check if the dataset is a reference.
    storedDatasetID=0
    if (datasetObject%isReference()) then
       ! Mark as a reference.
       isReference=.true.
       ! It is, so read the reference.
       dataBuffer=c_loc(referencedRegion)
       errorCode=h5dread(datasetObject%objectID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
       if (errorCode /= 0) then
          message="unable to read reference in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Now dereference the pointer.
       call h5rdereference_f(datasetObject%objectID,referencedRegion,dereferencedObjectID,errorCode)
       if (errorCode < 0) then
          message="unable to dereference pointer in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If the dataset object was opened internally, then close it.
       if (self%hdf5ObjectType /= hdf5ObjectTypeDataset) then
          call h5dclose_f(datasetObject%objectID,errorCode)
          if (errorCode < 0) then
             message="unable to close pointer dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Store the ID of this dataset so that we can replace it later.
          storedDatasetID=datasetObject%objectID
       end if
       ! The dataset object ID is now replaced with the referenced region ID.
       datasetObject%objectID=dereferencedObjectID
       ! Get the dataspace for this referenced region.
       call h5rget_region_f(dereferencedObjectID,referencedRegion,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of referenced region in dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       ! Mark as not reference.
       isReference=.false.
       ! Not a reference, so simply get the dataspace.
       call h5dget_space_f(datasetObject%objectID,datasetDataspaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to get dataspace of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Check that the object is a variable-length 2D integer-8 array.
    call datasetObject%assertDatasetType(H5T_VLEN_INTEGER8,1)

    ! Get the dimensions of the array to be read.
    if (isReference) then
       ! This is a reference, so get the extent of the referenced region.
       call h5sget_select_bounds_f(datasetDataspaceID,referenceStart,referenceEnd,errorCode)
       if (errorCode < 0) then
          message="unable to get bounds of referenced region for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Compute the dimensions of the referenced region.
       datasetDimensions=referenceEnd-referenceStart+1
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,referenceStart-1+readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,datasetDimensions,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       end if
    else
       ! Not a reference, so get the extent of the entire dataset.
       call h5sget_simple_extent_dims_f(datasetDataspaceID,datasetDimensions,datasetMaximumDimensions,errorCode)
       if (errorCode < 0) then
          message="unable to get dimensions of dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! If only a subsection is to be read, then select the appropriate hyperslab.
       if (readSubsection) then
          ! Check that subsection start values are legal.
          if (any(readBegin < 1 .or. readBegin > datasetDimensions)) then
             message="requested subsection begins outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Check that subsection extent is legal.
          if (any(readCount < 1 .or. readBegin+readCount-1 > datasetDimensions)) then
             message="requested subsection count is non-positive or outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab.
          call h5sselect_hyperslab_f(datasetDataspaceID,H5S_SELECT_SET_F,readBegin-1,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select filespace hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Set the size of the data to read in.
          datasetDimensions=readCount
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,readCount,memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,readCount,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else if (present(readSelection)) then
          ! A selection is to be read - create a suitable dataspace selection.
          ! Check that the selection is valid.
          if (any(readSelection < 1 .or. readSelection > datasetDimensions(1))) then
             message="requested selection extends outside of bounds of dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Create a map for selecting elements if necessary.
          allocate(readSelectionMap(1,size(readSelection)))
          forall(i=1:size(readSelection))
             readSelectionMap(:,i)=[readSelection(i)]
          end forall
          ! Create selection.
          call h5sselect_elements_f(datasetDataspaceID,H5S_SELECT_SET_F,1,size(readSelectionMap,dim=2,kind=size_t),readSelectionMap,errorCode)
          if (errorCode < 0) then
             message="could not select filespace selection for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          deallocate(readSelectionMap)
          ! Set the size of the data to read in.
          datasetDimensions(1)=size(readSelection)
          ! Construct a suitable memory space ID to read this data into.
          call h5screate_simple_f(1,int(shape(datasetValue),kind=HSIZE_T),memorySpaceID,errorCode)
          if (errorCode < 0) then
             message="unable to get create memory dataspace for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
          ! Select hyperslab to read to.
          referenceStart=0
          call h5sselect_hyperslab_f(memorySpaceID,H5S_SELECT_SET_F,referenceStart,datasetDimensions,errorCode)
          if (errorCode < 0) then
             message="could not select memory space hyperslab for dataset '"//datasetObject%objectName//"'"
             call Error_Report(message//self%locationReport()//{introspection:location})
          end if
       else
          ! Set the default memory space ID.
          memorySpaceID=H5S_ALL_F
       end if
    end if

    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    allocate(datasetValue (datasetDimensions(1)))
    allocate(datasetValueC(datasetDimensions(1)))
    ! Read the dataset.
    dataBuffer=c_loc(datasetValueC)
    errorCode=h5dread(datasetObject%objectID,H5T_VLEN_INTEGER8(1),memorySpaceID,datasetDataspaceID,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to read dataset '"//trim(datasetNameActual)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    do i=1,datasetDimensions(1)
       call c_f_pointer(datasetValueC(i)%p,datasetValue(i)%row,shape=[datasetValueC(i)%length])
    end do
    deallocate(datasetValueC)

    ! Close the dataspace.
    call h5sclose_f(datasetDataspaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataspace of dataset '"//datasetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Close the memory dataspace if necessary.
    if (memorySpaceID /= H5S_ALL_F) then
       call h5sclose_f(memorySpaceID,errorCode)
       if (errorCode /= 0) then
          message="unable to close memory dataspace for dataset '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! Determine how to close the object.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset .and. datasetObject%isReference()) then
       ! It was, so close the referenced dataset.
       call h5dclose_f(datasetObject%objectID,errorCode)
       if (errorCode < 0) then
          message="unable to close referenced dataset for '"//datasetObject%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Restore the object ID of the original dataset.
       self%objectID=storedDatasetID
    end if
    return
  end subroutine IO_HDF5_Read_Dataset_VarInteger8_2D_Array_Allocatable
  
  subroutine IO_HDF5_Write_Dataset_VarInteger8_2D(self,datasetValue,datasetName,comment,appendTo,chunkSize,compressionLevel,datasetReturned)
    !!{
    Open and write a variable-length integer-8 2-D array dataset in {\normalfont \ttfamily self}.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F, H5S_SELECT_SET_F  , h5sget_simple_extent_dims_f, HID_T        , &
          &                                      HSIZE_T      , h5dget_space_f    , h5dset_extent_f            , h5dwrite_vl_f, &
          &                                      h5sclose_f   , h5screate_simple_f, h5sselect_hyperslab_f      , hsize_t
    use            :: ISO_Varying_String, only : assignment(=), operator(//)      , trim
    implicit none
    class           (hdf5Object     ), intent(inout)                       :: self
    character       (len=*          ), intent(in   ), optional             :: comment                    , datasetName
    type            (hdf5VarInteger8), intent(in   ), dimension(:)         :: datasetValue
    logical                          , intent(in   ), optional             :: appendTo
    integer         (hsize_t        ), intent(in   ), optional             :: chunkSize
    integer                          , intent(in   ), optional             :: compressionLevel
    type            (hdf5Object     ), intent(  out), optional             :: datasetReturned
    integer         (kind=HSIZE_T   )               , dimension(1)         :: datasetDimensions          , hyperslabCount      , &
         &                                                                    hyperslabStart             , newDatasetDimensions, &
         &                                                                    newDatasetDimensionsMaximum
    type            (hdf5VlenC      ), allocatable  , dimension(:), target :: datasetValueC
    type            (c_ptr          )                                      :: datasetValueC_
    integer         (kind=HSIZE_T   )                                      :: i
    integer                                                                :: datasetRank                , errorCode
    integer         (kind=HID_T     )                                      :: dataspaceID                , newDataspaceID
    logical                                                                :: appendToActual             , preExisted
    type            (hdf5Object     )                                      :: datasetObject
    type            (varying_string )                                      :: datasetNameActual          , message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized()

    ! Get the name of the dataset.
    if (present(datasetName)) then
       datasetNameActual=     datasetName
    else
       datasetNameActual=self% objectName
    end if

    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to write dataset '"//trim(datasetNameActual)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Determine append status.
    if (present(appendTo)) then
       appendToActual=appendTo
    else
       appendToActual=.false.
    end if
    ! Determine dataset dimensions
    datasetDimensions=shape(datasetValue)
    ! Check if the object is a dataset, or something else.
    if (self%hdf5ObjectType == hdf5ObjectTypeDataset) then
       ! If this dataset if not overwritable, report an error.
       if (.not.(self%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetNameActual)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       else
          ! Check that the object is a 1D vlen integer8.
          call self%assertDatasetType(H5T_VLEN_INTEGER8,1)
       end if
       select type (self)
       type is (hdf5Object)
          datasetObject=self
       end select
       datasetNameActual=self%objectName
       preExisted       =.true.
    else
       ! Check that an dataset name was supplied.
       if (.not.present(datasetName)) then
          message="no name was supplied for dataset in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       ! Record if dataset already exists.
       preExisted=self%hasDataset(datasetName)
       ! Open the dataset.
       datasetObject=self%openDataset(datasetName,comment,hdf5DataTypeVlenInteger8,datasetDimensions,appendTo&
            &=appendTo,chunkSize=chunkSize,compressionLevel=compressionLevel)
       ! Check that pre-existing object is a variable-length 2D integer-8.
       if (preExisted) call datasetObject%assertDatasetType(H5T_VLEN_INTEGER8,1)
       ! If this dataset if not overwritable, report an error.
       if (preExisted.and..not.(datasetObject%isOverwritable.or.appendToActual)) then
          message="dataset '"//trim(datasetName)//"' is not overwritable"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if

    ! If appending is requested, get the size of the existing dataset.
    if (appendToActual.and.preExisted) then
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sget_simple_extent_dims_f(dataspaceID,newDatasetDimensions,newDatasetDimensionsMaximum,errorCode)
       if (errorCode < 0) then
          message="could not get dataspace extent for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       call h5sclose_f(dataspaceID,errorCode)
       if (errorCode < 0) then
          message="could not close dataspace for dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       hyperslabStart      =newDatasetDimensions
       hyperslabCount      =dataSetDimensions
       newDatasetDimensions=newDatasetDimensions+datasetDimensions
    else
       newDatasetDimensions=datasetDimensions
       hyperslabStart      =0
       hyperslabCount      =datasetDimensions
    end if

    ! Set extent of the dataset.
    if (datasetObject%chunkSize /= -1) then
       call h5dset_extent_f(datasetObject%objectID,newDatasetDimensions,errorCode)
       if (errorCode < 0) then
          message="could not set extent of dataset '"//trim(datasetNameActual)//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    end if
    ! Get the dataspace for the dataset.
    call h5dget_space_f(datasetObject%objectID,dataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not get dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Select hyperslab to write.
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) then
       message="could not select hyperslab for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Create a dataspace for the data to be written.
    datasetRank=1
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="could not create dataspace for data to be written to dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Write the dataset.
    allocate(datasetValueC(datasetDimensions(1)))
    do i=1,datasetDimensions(1)
       datasetValueC(i)%length=size (datasetValue(i)%row,kind=c_size_t)
       datasetValueC(i)%p     =c_loc(datasetValue(i)%row              )
    end do
    datasetValueC_=c_loc(datasetValueC)
    errorCode     =H5Dwrite(datasetObject%objectID,H5T_VLEN_INTEGER8(1),newDataspaceID,dataspaceID,H5P_DEFAULT_F,datasetValueC_)
    if (errorCode /= 0) then
       message="unable to write dataset '"//datasetNameActual//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    deallocate(datasetValueC)

    ! Close the dataspaces.
    call h5sclose_f(dataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) then
       message="unable to close new dataspace for dataset '"//trim(datasetNameActual)//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if

    ! Copy the dataset to return if necessary.
    if (present(datasetReturned)) datasetReturned=datasetObject
    return
  end subroutine IO_HDF5_Write_Dataset_VarInteger8_2D
  
  !! Table routines.

  subroutine IO_HDF5_Read_Table_Real_1D_Array_Allocatable(self,tableName,columnName,datasetValue,readBegin,readCount)
    !!{
    Open and read a real 1D array from a table {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: H5TB              , only : h5tbget_table_info_f, h5tbread_field_name_f
    use :: HDF5              , only : H5T_NATIVE_REAL     , HSIZE_T              , h5tget_size_f
    use :: ISO_Varying_String, only : assignment(=)       , operator(//)
    implicit none
    real                     , allocatable, dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   )           :: tableName      , columnName
    integer  (kind=HSIZE_T  )                           , intent(in   ), optional :: readBegin      , readCount
    integer  (kind=HSIZE_T  )                                                     :: readBeginActual, readCountActual, &
         &                                                                           fieldCount     , recordCount
    integer  (kind=SIZE_T   )                                                     :: recordTypeSize
    integer                                                                       :: errorCode
    logical                                                                       :: readSubsection
    type     (varying_string)                                                     :: message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized
    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read table '"//trim(tableName)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(tableName)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(tableName)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if

    ! Check that the object is a group or a file
    if (self%hdf5ObjectType == hdf5ObjectTypeFile .or. self%hdf5ObjectType == hdf5ObjectTypeGroup) then
       ! Check that the dataset exists.
       if (.not.self%hasDataset(tableName)) then
          message="table '"//trim(tableName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       message="attempt to read table from '"//self%objectName//"' which is neither a file or a group"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Get the table dimensions.
    call h5tbget_table_info_f(self%objectID,tableName,fieldCount,recordCount,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of table '"//tableName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Determine records to read.
    if (readSubsection) then
       readBeginActual=readBegin
       readCountActual=readCount
    else
       readBeginActual=0
       readCountActual=recordCount
    end if
    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    allocate(datasetValue(readCountActual))

    ! Read the column.
    call h5tget_size_f(H5T_NATIVE_REAL,recordTypeSize,errorCode)
    if (errorCode /= 0) then
       message="unable to get real datatype size"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tbread_field_name_f(self%objectID,tableName,columnName,readBeginActual,readCountActual,recordTypeSize,datasetValue,errorCode)
    if (errorCode /= 0) then
       message="unable to read table '"//trim(tableName)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Table_Real_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Table_Integer_1D_Array_Allocatable(self,tableName,columnName,datasetValue,readBegin,readCount)
    !!{
    Open and read an integer 1D array from a table {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: H5TB              , only : h5tbget_table_info_f, h5tbread_field_name_f
    use :: HDF5              , only : H5T_NATIVE_REAL     , HSIZE_T              , h5tget_size_f
    use :: ISO_Varying_String, only : assignment(=)       , operator(//)
    implicit none
    integer                  , allocatable, dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object    )                           , intent(inout)           :: self
    character(len=*         )                           , intent(in   )           :: tableName      , columnName
    integer  (kind=HSIZE_T  )                           , intent(in   ), optional :: readBegin      , readCount
    integer  (kind=HSIZE_T  )                                                     :: readBeginActual, readCountActual, &
         &                                                                           fieldCount     , recordCount
    integer  (kind=SIZE_T   )                                                     :: recordTypeSize
    integer                                                                       :: errorCode
    logical                                                                       :: readSubsection
    type     (varying_string)                                                     :: message

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized
    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read table '"//trim(tableName)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(tableName)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(tableName)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Check that the object is a group or a file
    if (self%hdf5ObjectType == hdf5ObjectTypeFile .or. self%hdf5ObjectType == hdf5ObjectTypeGroup) then
       ! Check that the dataset exists.
       if (.not.self%hasDataset(tableName)) then
          message="table '"//trim(tableName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       message="attempt to read table from '"//self%objectName//"' which is neither a file or a group"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Get the table dimensions.
    call h5tbget_table_info_f(self%objectID,tableName,fieldCount,recordCount,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of table '"//tableName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Determine records to read.
    if (readSubsection) then
       readBeginActual=readBegin
       readCountActual=readCount
    else
       readBeginActual=0
       readCountActual=recordCount
    end if
    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    allocate(datasetValue(readCountActual))
    ! Read the column.
    call h5tget_size_f(H5T_NATIVE_REAL,recordTypeSize,errorCode)
    if (errorCode /= 0) then
       message="unable to get real datatype size"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    call h5tbread_field_name_f(self%objectID,tableName,columnName,readBeginActual,readCountActual,recordTypeSize,datasetValue,errorCode)
    if (errorCode /= 0) then
       message="unable to read table '"//trim(tableName)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Table_Integer_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Table_Integer8_1D_Array_Allocatable(self,tableName,columnName,datasetValue,readBegin,readCount)
    !!{
    Open and read a real scalar from a table {\normalfont \ttfamily self}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: H5TB              , only : h5tbget_table_info_f
    use            :: HDF5              , only : HSIZE_T             , h5tget_size_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc               , c_null_char
    use            :: ISO_Varying_String, only : assignment(=)       , operator(//)
    implicit none
    integer  (kind=kind_int8     ), allocatable, dimension(:), intent(  out), target   :: datasetValue
    class    (hdf5Object         )                           , intent(inout)           :: self
    character(len=*              )                           , intent(in   )           :: tableName      , columnName
    integer  (kind=HSIZE_T  )                                , intent(in   ), optional :: readBegin      , readCount
    integer  (kind=HSIZE_T       )                                                     :: fieldCount     , recordCount    , &
         &                                                                                readBeginActual, readCountActual
    integer  (kind=SIZE_T        )                                                     :: recordTypeSize
    integer                                                                            :: errorCode
    type     (varying_string     )                                                     :: message
    type     (c_ptr              )                                                     :: dataValueC
    logical                                                                            :: readSubsection

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized
    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read table '"//trim(tableName)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(tableName)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(tableName)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Check that the object is a group or a file
    if (self%hdf5ObjectType == hdf5ObjectTypeFile .or. self%hdf5ObjectType == hdf5ObjectTypeGroup) then
       ! Check that the dataset exists.
       if (.not.self%hasDataset(tableName)) then
          message="table '"//trim(tableName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       message="attempt to read table from '"//self%objectName//"' which is neither a file or a group"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Get the table dimensions.
    call h5tbget_table_info_f(self%objectID,tableName,fieldCount,recordCount,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of table '"//tableName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Determine records to read.
    if (readSubsection) then
       readBeginActual=readBegin
       readCountActual=readCount
    else
       readBeginActual=0
       readCountActual=recordCount
    end if
    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    allocate(datasetValue(readCountActual))
    ! Read the column.
    call h5tget_size_f(H5T_INTEGER8,recordTypeSize,errorCode)
    if (errorCode /= 0) then
       message="unable to get long integer datatype size"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    dataValueC=c_loc(datasetValue)
    errorCode=H5TBread_fields_name(self%objectID,trim(tableName)//C_NULL_CHAR,trim(columnName)//C_NULL_CHAR,readBeginActual,readCountActual,recordTypeSize,[0_size_t],[8_size_t],dataValueC)
    if (errorCode /= 0) then
       message="unable to read table '"//trim(tableName)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Read_Table_Integer8_1D_Array_Allocatable

  subroutine IO_HDF5_Read_Table_Character_1D_Array_Allocatable(self,tableName,columnName,datasetValue,readBegin,readCount)
    !!{
    Open and read a real 1D array from a table {\normalfont \ttfamily self}.
    !!}
    use :: Error             , only : Error_Report
    use :: H5TB              , only : h5tbget_table_info_f, h5tbread_field_name_f
    use :: HDF5              , only : HSIZE_T
    use :: ISO_Varying_String, only : assignment(=)       , operator(//)
    implicit none
    character(len=*                ), allocatable, dimension(:), intent(  out)           :: datasetValue
    class    (hdf5Object           )                           , intent(inout)           :: self
    character(len=*                )                           , intent(in   )           :: tableName      , columnName
    integer  (kind=HSIZE_T         )                           , intent(in   ), optional :: readBegin      , readCount
    integer  (kind=HSIZE_T         )                                                     :: readBeginActual, readCountActual, &
         &                                                                                  fieldCount     , recordCount
    integer  (kind=SIZE_T          )                                                     :: recordTypeSize
    integer                                                                              :: errorCode     , i
    logical                                                                              :: readSubsection
    type     (varying_string       )                                                     :: message
    character(len=len(datasetValue))                                                     :: convertValue

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized
    ! Check that the object is already open.
    if (.not.self%isOpenValue) then
       message="attempt to read table '"//trim(tableName)//"' in unopen object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! If a subsection is to be read, we need both start and count values.
    if (present(readBegin)) then
       if (.not.present(readCount)) then
          message="reading a subsection of dataset '"//trim(tableName)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.true.
    else
       if (present(readCount)) then
          message="reading a subsection of dataset '"//trim(tableName)//"' requires both readBegin and readCount to be specified"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
       readSubsection=.false.
    end if
    ! Check that the object is a group or a file
    if (self%hdf5ObjectType == hdf5ObjectTypeFile .or. self%hdf5ObjectType == hdf5ObjectTypeGroup) then
       ! Check that the dataset exists.
       if (.not.self%hasDataset(tableName)) then
          message="table '"//trim(tableName)//"' does not exist in '"//self%objectName//"'"
          call Error_Report(message//self%locationReport()//{introspection:location})
       end if
    else
       message="attempt to read table from '"//self%objectName//"' which is neither a file or a group"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Get the table dimensions.
    call h5tbget_table_info_f(self%objectID,tableName,fieldCount,recordCount,errorCode)
    if (errorCode < 0) then
       message="unable to get dimensions of table '"//tableName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Determine records to read.
    if (readSubsection) then
       readBeginActual=readBegin
       readCountActual=readCount
    else
       readBeginActual=0
       readCountActual=recordCount
    end if
    ! Allocate the array to the appropriate size.
    if (allocated(datasetValue)) deallocate(datasetValue)
    allocate(datasetValue(readCountActual))
    ! Read the column.
    recordTypeSize=len(datasetValue(1))
    call h5tbread_field_name_f(self%objectID,tableName,columnName,readBeginActual,readCountActual,recordTypeSize,datasetValue,errorCode)
    if (errorCode /= 0) then
       message="unable to read table '"//trim(tableName)//"' in object '"//self%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    ! Convert to Fortran form.
    do i=1,size(datasetValue)
       convertValue=datasetValue(i)
       if (index(convertValue,char(0)) /= 0) then
          convertValue(index(convertValue,char(0)):len(convertValue))=repeat(" ",len(convertValue)-index(convertValue,char(0))+1)
          datasetValue(i)=convertValue
       end if
    end do
    return
  end subroutine IO_HDF5_Read_Table_Character_1D_Array_Allocatable

  !! Reference routines.

  subroutine IO_HDF5_Create_Reference_Scalar_To_1D(fromGroup,toDataset,referenceName,referenceStart,referenceCount)
    !!{
    Create a scalar reference to the 1-D {\normalfont \ttfamily toDataset} in the HDF5 group {\normalfont \ttfamily fromGroup}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F        , H5S_SELECT_SET_F, H5T_STD_REF_DSETREG, &
          &                                      HID_T                , HSIZE_T          , h5dclose_f      , h5dcreate_f        , &
          &                                      h5dget_space_f       , h5rcreate_f      , h5sclose_f      , h5screate_simple_f , &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , char             , operator(//)    , trim
    implicit none
    class    (hdf5Object       )              , intent(inout)         :: fromGroup
    type     (hdf5Object       )              , intent(inout)         :: toDataset
    character(len=*            )              , intent(in   )         :: referenceName
    integer  (kind=HSIZE_T     ), dimension(1), intent(in   )         :: referenceCount   , referenceStart
    integer  (kind=HSIZE_T     ), dimension(1)                        :: datasetDimensions, hyperslabCount, hyperslabStart
    type     (hdset_reg_ref_t_f)                             , target :: dataReference
    integer                                                           :: errorCode        , datasetRank
    integer  (kind=HID_T       )                                      :: dataSetID        , dataSpaceID   , dataSubsetSpaceID
    type     (varying_string   )                                      :: message
    type     (c_ptr            )                                      :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the group is already open.
    if (.not.fromGroup%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' in unopen group '"//fromGroup%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check that the dataset is already open.
    if (.not.toDataset%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' to unopen dataset '"//toDataset%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the group really is a group.
    if (fromGroup%hdf5ObjectType /= hdf5ObjectTypeGroup) then
       message="attempt to write reference '"//trim(referenceName)//"' into object '"//fromGroup%objectName//"' which is not a group"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the dataset really is a dataset.
    if (toDataset%hdf5ObjectType /= hdf5ObjectTypeDataset) then
       message="attempt to write reference '"//trim(referenceName)//"' to object '"//toDataset%objectName//"' which is not a dataset"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Get the dataspace of the dataset.
    call h5dget_space_f(toDataset%objectID,dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace for dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Select a hyperslab from this dataspace. Subtract one from the start position since HDF5 uses indexing beginning at 0.
    hyperslabStart=referenceStart-1
    hyperslabCount=referenceCount
    call h5sselect_hyperslab_f(dataSubsetSpaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode /= 0) then
       message="unable to get select hyperslab in dataspace of dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create a dataspace for the reference dataset.
    datasetRank      =0
    datasetDimensions=1
    call h5screate_simple_f(datasetRank,datasetDimensions,dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference dataset.
    call h5dcreate_f(fromGroup%objectID,trim(referenceName),H5T_STD_REF_DSETREG,dataSpaceID,dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference.
    call h5rcreate_f(toDataset%parentObject%objectID,char(toDataset%objectName),dataSubsetSpaceID,dataReference,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Write the reference dataset.
    dataBuffer=c_loc(dataReference)
    errorCode=h5dwrite(dataSetID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset dataspace.
    call h5sclose_f(dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataset dataspace for '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the reference dataset dataspace.
    call h5sclose_f(dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset.
    call h5dclose_f(dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    return
  end subroutine IO_HDF5_Create_Reference_Scalar_To_1D

  subroutine IO_HDF5_Create_Reference_Scalar_To_2D(fromGroup,toDataset,referenceName,referenceStart,referenceCount)
    !!{
    Create a scalar reference to the 2-D {\normalfont \ttfamily toDataset} in the HDF5 group {\normalfont \ttfamily fromGroup}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F        , H5S_SELECT_SET_F, H5T_STD_REF_DSETREG, &
          &                                      HID_T                , HSIZE_T          , h5dclose_f      , h5dcreate_f        , &
          &                                      h5dget_space_f       , h5rcreate_f      , h5sclose_f      , h5screate_simple_f , &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , char             , operator(//)    , trim
    implicit none
    class    (hdf5Object       )              , intent(inout)         :: fromGroup
    type     (hdf5Object       )              , intent(inout)         :: toDataset
    character(len=*            )              , intent(in   )         :: referenceName
    integer  (kind=HSIZE_T     ), dimension(2), intent(in   )         :: referenceCount   , referenceStart
    integer  (kind=HSIZE_T     ), dimension(2)                        :: datasetDimensions, hyperslabCount, hyperslabStart
    type     (hdset_reg_ref_t_f)                             , target :: dataReference
    integer                                                           :: errorCode        , datasetRank
    integer  (kind=HID_T       )                                      :: dataSetID        , dataSpaceID   , dataSubsetSpaceID
    type     (varying_string   )                                      :: message
    type     (c_ptr            )                                      :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the group is already open.
    if (.not.fromGroup%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' in unopen group '"//fromGroup%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check that the dataset is already open.
    if (.not.toDataset%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' to unopen dataset '"//toDataset%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the group really is a group.
    if (fromGroup%hdf5ObjectType /= hdf5ObjectTypeGroup) then
       message="attempt to write reference '"//trim(referenceName)//"' into object '"//fromGroup%objectName//"' which is not a group"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the dataset really is a dataset.
    if (toDataset%hdf5ObjectType /= hdf5ObjectTypeDataset) then
       message="attempt to write reference '"//trim(referenceName)//"' to object '"//toDataset%objectName//"' which is not a dataset"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Get the dataspace of the dataset.
    call h5dget_space_f(toDataset%objectID,dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace for dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Select a hyperslab from this dataspace. Subtract one from the start position since HDF5 uses indexing beginning at 0.
    hyperslabStart=referenceStart-1
    hyperslabCount=referenceCount
    call h5sselect_hyperslab_f(dataSubsetSpaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode /= 0) then
       message="unable to get select hyperslab in dataspace of dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create a dataspace for the reference dataset.
    datasetRank      =0
    datasetDimensions=1
    call h5screate_simple_f(datasetRank,datasetDimensions,dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference dataset.
    call h5dcreate_f(fromGroup%objectID,trim(referenceName),H5T_STD_REF_DSETREG,dataSpaceID,dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference.
    call h5rcreate_f(toDataset%parentObject%objectID,char(toDataset%objectName),dataSubsetSpaceID,dataReference,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Write the reference dataset.
    dataBuffer=c_loc(dataReference)
    errorCode=h5dwrite(dataSetID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset dataspace.
    call h5sclose_f(dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataset dataspace for '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the reference dataset dataspace.
    call h5sclose_f(dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset.
    call h5dclose_f(dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    return
  end subroutine IO_HDF5_Create_Reference_Scalar_To_2D

  subroutine IO_HDF5_Create_Reference_Scalar_To_3D(fromGroup,toDataset,referenceName,referenceStart,referenceCount)
    !!{
    Create a scalar reference to the 3-D {\normalfont \ttfamily toDataset} in the HDF5 group {\normalfont \ttfamily fromGroup}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F        , H5S_SELECT_SET_F, H5T_STD_REF_DSETREG, &
          &                                      HID_T                , HSIZE_T          , h5dclose_f      , h5dcreate_f        , &
          &                                      h5dget_space_f       , h5rcreate_f      , h5sclose_f      , h5screate_simple_f , &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , char             , operator(//)    , trim
    implicit none
    class    (hdf5Object       )              , intent(inout)         :: fromGroup
    type     (hdf5Object       )              , intent(inout)         :: toDataset
    character(len=*            )              , intent(in   )         :: referenceName
    integer  (kind=HSIZE_T     ), dimension(3), intent(in   )         :: referenceCount   , referenceStart
    integer  (kind=HSIZE_T     ), dimension(3)                        :: datasetDimensions, hyperslabCount, hyperslabStart
    type     (hdset_reg_ref_t_f)                             , target :: dataReference
    integer                                                           :: errorCode        , datasetRank
    integer  (kind=HID_T       )                                      :: dataSetID        , dataSpaceID   , dataSubsetSpaceID
    type     (varying_string   )                                      :: message
    type     (c_ptr            )                                      :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the group is already open.
    if (.not.fromGroup%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' in unopen group '"//fromGroup%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check that the dataset is already open.
    if (.not.toDataset%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' to unopen dataset '"//toDataset%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the group really is a group.
    if (fromGroup%hdf5ObjectType /= hdf5ObjectTypeGroup) then
       message="attempt to write reference '"//trim(referenceName)//"' into object '"//fromGroup%objectName//"' which is not a group"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the dataset really is a dataset.
    if (toDataset%hdf5ObjectType /= hdf5ObjectTypeDataset) then
       message="attempt to write reference '"//trim(referenceName)//"' to object '"//toDataset%objectName//"' which is not a dataset"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Get the dataspace of the dataset.
    call h5dget_space_f(toDataset%objectID,dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace for dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Select a hyperslab from this dataspace. Subtract one from the start position since HDF5 uses indexing beginning at 0.
    hyperslabStart=referenceStart-1
    hyperslabCount=referenceCount
    call h5sselect_hyperslab_f(dataSubsetSpaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode /= 0) then
       message="unable to get select hyperslab in dataspace of dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create a dataspace for the reference dataset.
    datasetRank      =0
    datasetDimensions=1
    call h5screate_simple_f(datasetRank,datasetDimensions,dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference dataset.
    call h5dcreate_f(fromGroup%objectID,trim(referenceName),H5T_STD_REF_DSETREG,dataSpaceID,dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference.
    call h5rcreate_f(toDataset%parentObject%objectID,char(toDataset%objectName),dataSubsetSpaceID,dataReference,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Write the reference dataset.
    dataBuffer=c_loc(dataReference)
    errorCode=h5dwrite(dataSetID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset dataspace.
    call h5sclose_f(dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataset dataspace for '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the reference dataset dataspace.
    call h5sclose_f(dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset.
    call h5dclose_f(dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    return
  end subroutine IO_HDF5_Create_Reference_Scalar_To_3D

  subroutine IO_HDF5_Create_Reference_Scalar_To_4D(fromGroup,toDataset,referenceName,referenceStart,referenceCount)
    !!{
    Create a scalar reference to the 4-D {\normalfont \ttfamily toDataset} in the HDF5 group {\normalfont \ttfamily fromGroup}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F        , H5S_SELECT_SET_F, H5T_STD_REF_DSETREG, &
          &                                      HID_T                , HSIZE_T          , h5dclose_f      , h5dcreate_f        , &
          &                                      h5dget_space_f       , h5rcreate_f      , h5sclose_f      , h5screate_simple_f , &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , char             , operator(//)    , trim
    implicit none
    class    (hdf5Object       )              , intent(inout)         :: fromGroup
    type     (hdf5Object       )              , intent(inout)         :: toDataset
    character(len=*            )              , intent(in   )         :: referenceName
    integer  (kind=HSIZE_T     ), dimension(4), intent(in   )         :: referenceCount   , referenceStart
    integer  (kind=HSIZE_T     ), dimension(4)                        :: datasetDimensions, hyperslabCount, hyperslabStart
    type     (hdset_reg_ref_t_f)                             , target :: dataReference
    integer                                                           :: errorCode        , datasetRank
    integer  (kind=HID_T       )                                      :: dataSetID        , dataSpaceID   , dataSubsetSpaceID
    type     (varying_string   )                                      :: message
    type     (c_ptr            )                                      :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the group is already open.
    if (.not.fromGroup%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' in unopen group '"//fromGroup%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check that the dataset is already open.
    if (.not.toDataset%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' to unopen dataset '"//toDataset%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the group really is a group.
    if (fromGroup%hdf5ObjectType /= hdf5ObjectTypeGroup) then
       message="attempt to write reference '"//trim(referenceName)//"' into object '"//fromGroup%objectName//"' which is not a group"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the dataset really is a dataset.
    if (toDataset%hdf5ObjectType /= hdf5ObjectTypeDataset) then
       message="attempt to write reference '"//trim(referenceName)//"' to object '"//toDataset%objectName//"' which is not a dataset"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Get the dataspace of the dataset.
    call h5dget_space_f(toDataset%objectID,dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace for dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Select a hyperslab from this dataspace. Subtract one from the start position since HDF5 uses indexing beginning at 0.
    hyperslabStart=referenceStart-1
    hyperslabCount=referenceCount
    call h5sselect_hyperslab_f(dataSubsetSpaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode /= 0) then
       message="unable to get select hyperslab in dataspace of dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create a dataspace for the reference dataset.
    datasetRank      =0
    datasetDimensions=1
    call h5screate_simple_f(datasetRank,datasetDimensions,dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference dataset.
    call h5dcreate_f(fromGroup%objectID,trim(referenceName),H5T_STD_REF_DSETREG,dataSpaceID,dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference.
    call h5rcreate_f(toDataset%parentObject%objectID,char(toDataset%objectName),dataSubsetSpaceID,dataReference,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Write the reference dataset.
    dataBuffer=c_loc(dataReference)
    errorCode=h5dwrite(dataSetID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset dataspace.
    call h5sclose_f(dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataset dataspace for '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the reference dataset dataspace.
    call h5sclose_f(dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset.
    call h5dclose_f(dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    return
  end subroutine IO_HDF5_Create_Reference_Scalar_To_4D

  subroutine IO_HDF5_Create_Reference_Scalar_To_5D(fromGroup,toDataset,referenceName,referenceStart,referenceCount)
    !!{
    Create a scalar reference to the 5-D {\normalfont \ttfamily toDataset} in the HDF5 group {\normalfont \ttfamily fromGroup}.
    !!}
    use            :: Error             , only : Error_Report
    use            :: HDF5              , only : H5P_DEFAULT_F        , H5S_ALL_F        , H5S_SELECT_SET_F, H5T_STD_REF_DSETREG, &
          &                                      HID_T                , HSIZE_T          , h5dclose_f      , h5dcreate_f        , &
          &                                      h5dget_space_f       , h5rcreate_f      , h5sclose_f      , h5screate_simple_f , &
          &                                      h5sselect_hyperslab_f, hdset_reg_ref_t_f
    use, intrinsic :: ISO_C_Binding     , only : c_loc
    use            :: ISO_Varying_String, only : assignment(=)        , char             , operator(//)    , trim
    implicit none
    class    (hdf5Object       )              , intent(inout)         :: fromGroup
    type     (hdf5Object       )              , intent(inout)         :: toDataset
    character(len=*            )              , intent(in   )         :: referenceName
    integer  (kind=HSIZE_T     ), dimension(5), intent(in   )         :: referenceCount   , referenceStart
    integer  (kind=HSIZE_T     ), dimension(5)                        :: datasetDimensions, hyperslabCount, hyperslabStart
    type     (hdset_reg_ref_t_f)                             , target :: dataReference
    integer                                                           :: errorCode        , datasetRank
    integer  (kind=HID_T       )                                      :: dataSetID        , dataSpaceID   , dataSubsetSpaceID
    type     (varying_string   )                                      :: message
    type     (c_ptr            )                                      :: dataBuffer

    ! Check that this module is initialized.
    call IO_HDF_Assert_Is_Initialized

    ! Check that the group is already open.
    if (.not.fromGroup%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' in unopen group '"//fromGroup%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check that the dataset is already open.
    if (.not.toDataset%isOpenValue) then
       message="attempt to write reference '"//trim(referenceName)//"' to unopen dataset '"//toDataset%objectName//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the group really is a group.
    if (fromGroup%hdf5ObjectType /= hdf5ObjectTypeGroup) then
       message="attempt to write reference '"//trim(referenceName)//"' into object '"//fromGroup%objectName//"' which is not a group"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Check if the dataset really is a dataset.
    if (toDataset%hdf5ObjectType /= hdf5ObjectTypeDataset) then
       message="attempt to write reference '"//trim(referenceName)//"' to object '"//toDataset%objectName//"' which is not a dataset"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Get the dataspace of the dataset.
    call h5dget_space_f(toDataset%objectID,dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to get dataspace for dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Select a hyperslab from this dataspace. Subtract one from the start position since HDF5 uses indexing beginning at 0.
    hyperslabStart=referenceStart-1
    hyperslabCount=referenceCount
    call h5sselect_hyperslab_f(dataSubsetSpaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode /= 0) then
       message="unable to get select hyperslab in dataspace of dataset '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create a dataspace for the reference dataset.
    datasetRank      =0
    datasetDimensions=1
    call h5screate_simple_f(datasetRank,datasetDimensions,dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference dataset.
    call h5dcreate_f(fromGroup%objectID,trim(referenceName),H5T_STD_REF_DSETREG,dataSpaceID,dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Create the reference.
    call h5rcreate_f(toDataset%parentObject%objectID,char(toDataset%objectName),dataSubsetSpaceID,dataReference,errorCode)
    if (errorCode /= 0) then
       message="unable to create reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Write the reference dataset.
    dataBuffer=c_loc(dataReference)
    errorCode=h5dwrite(dataSetID,H5T_STD_REF_DSETREG,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F,dataBuffer)
    if (errorCode /= 0) then
       message="unable to write reference '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset dataspace.
    call h5sclose_f(dataSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to close dataset dataspace for '"//trim(toDataset%objectName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the reference dataset dataspace.
    call h5sclose_f(dataSubsetSpaceID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataspace for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    ! Close the dataset.
    call h5dclose_f(dataSetID,errorCode)
    if (errorCode /= 0) then
       message="unable to reference dataset for '"//trim(referenceName)//"'"
       call Error_Report(message//fromGroup%locationReport()//{introspection:location})
    end if

    return
  end subroutine IO_HDF5_Create_Reference_Scalar_To_5D

  logical function IO_HDF5_Is_Reference(dataset)
    !!{
    Return true if the input dataset is a scalar reference.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : H5T_STD_REF_DSETREG, HID_T       , h5dget_type_f, h5tclose_f, &
          &                           h5tequal_f
    use :: ISO_Varying_String, only : assignment(=)      , operator(//)
    implicit none
    class  (hdf5Object    ), intent(in   ) :: dataset
    integer                                :: errorCode
    integer(kind=HID_T    )                :: dataTypeID
    type   (varying_string)                :: message

    ! Ensure that the dataset is open.
    if (.not.dataset%isOpenValue) then
       message="attempt to check if reference on unopen dataset '"//dataset%objectName//"'"
       call Error_Report(message//dataset%locationReport()//{introspection:location})
    end if

    ! Get the type of the object
    call h5dget_type_f(dataset%objectID,dataTypeID,errorCode)
    if (errorCode < 0) then
       message="unable to get data type for dataset '"//dataset%objectName//"'"
       call Error_Report(message//dataset%locationReport()//{introspection:location})
    end if

    ! Test the type.
    call h5tequal_f(dataTypeID,H5T_STD_REF_DSETREG,IO_HDF5_Is_Reference,errorCode)
    if (errorCode < 0) then
       message="unable to test data type for dataset '"//dataset%objectName//"'"
       call Error_Report(message//dataset%locationReport()//{introspection:location})
    end if

    ! Close the data type.
    call h5tclose_f(dataTypeID,errorCode)
    if (errorCode /= 0) then
       message="unable to close datatype of dataset '"//dataset%objectName//"'"
       call Error_Report(message//dataset%locationReport()//{introspection:location})
    end if

    return
  end function IO_HDF5_Is_Reference

  subroutine IO_HDF5_Copy(self,source,targetObject)
    !!{
    Copy the named object to the target object.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5              , only : h5ocopy_f
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    implicit none
    class    (hdf5Object    ), intent(in   ) :: self
    character(len=*         ), intent(in   ) :: source
    type     (hdf5Object    ), intent(inout) :: targetObject
    integer                                  :: errorCode
    type     (varying_string)                :: message

    call h5ocopy_f(self%objectID,source,targetObject%objectID,source,errorCode)
    if (errorCode < 0) then
       message="unable to copy object '"//source//"' from '"//self%objectName//"' to '"//targetObject%objectName//"'"
       call Error_Report(message//self%locationReport()//{introspection:location})
    end if
    return
  end subroutine IO_HDF5_Copy

  function IO_HDF5_Parent(self) result(parent)
    !!{
    Return the parent object.
    !!}
    implicit none
    class(hdf5Object) , pointer       :: parent
    class(hdf5Object) , intent(in   ) :: self

    parent => self%parentObject
    return
  end function IO_HDF5_Parent

  logical function IO_HDF5_Is_HDF5(fileName)
    !!{
    Return true if the named file is an HDF5 file.
    !!}
    use :: File_Utilities, only : File_Exists
    use :: Error         , only : Error_Report
    use :: HDF5          , only : h5fis_hdf5_f
    implicit none
    character(len=*), intent(in   ) :: fileName
    integer                         :: errorCode

    if (File_Exists(trim(fileName))) then
       call h5fis_hdf5_f(trim(fileName),IO_HDF5_Is_HDF5,errorCode)
       if (errorCode /= 0) call Error_Report('failed to determine nature of file'//{introspection:location})
    else
       IO_HDF5_Is_HDF5=.false.
    end if
    return    
  end function IO_HDF5_Is_HDF5

  subroutine IO_HDF5_Deep_Copy(self,destination)
    !!{
    Make a deep copy of the object, with a new HDF5 object identifier.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char        , var_str, operator(//)
    implicit none
    class(hdf5Object) , intent(in   ) :: self
    type (hdf5Object) , intent(inout) :: destination

    ! Get a new object ID.
    if (self%isOpen()) then
       select case (self%hdf5ObjectType)
       case (hdf5ObjectTypeNull     )
          ! Nothing to do in this case.
       case (hdf5ObjectTypeFile     )
          !![
          <conditionalCall>
           <call>destination=hdf5Object                          (char(self%objectFile),overWrite=.false.,readOnly=self%readOnly,objectsOverwritable=self%isOverwritable{conditions})</call>
           <argument name="compressionLevel" value="self%compressionLevel" condition="self%compressionLevelSet"/>
           <argument name="chunkSize"        value="self%chunkSize"        condition="self%chunkSizeSet"/>
          </conditionalCall>
          !!]
       case (hdf5ObjectTypeGroup    )
          !![
          <conditionalCall>
           <call>     destination=self%parentObject%openGroup    (char(self%objectName)                                         ,objectsOverwritable=self%isOverwritable{conditions})</call>
           <argument name="compressionLevel" value="self%compressionLevel" condition="self%compressionLevelSet"/>
           <argument name="chunkSize"        value="self%chunkSize"        condition="self%chunkSizeSet"/>
          </conditionalCall>
          !!]
       case (hdf5ObjectTypeDataset  )
          !![
          <conditionalCall>
           <call>     destination=self%parentObject%openDataset  (char(self%objectName)                                         ,     isOverwritable=self%isOverwritable{conditions})</call>
           <argument name="compressionLevel" value="self%compressionLevel" condition="self%compressionLevelSet"/>
           <argument name="chunkSize"        value="self%chunkSize"        condition="self%chunkSizeSet"/>
          </conditionalCall>
          !!]
       case (hdf5ObjectTypeAttribute)
          destination               =self%parentObject%openAttribute(char(self%objectName)                                      ,     isOverwritable=self%isOverwritable            )
       case default
          call Error_Report(var_str('unknown HDF5 object type')//self%locationReport()//{introspection:location})
       end select
    end if
    return
  end subroutine IO_HDF5_Deep_Copy

end module IO_HDF5
