!% Contains a module which handles creation of groups and writing of datasets to the \glc\ output file.

module Galacticus_HDF5_Groups
  !% Handles creation of groups and writing of datasets to the \glc\ output file.
  use Galacticus_Error
  use Galacticus_HDF5
  use ISO_Varying_String
  use HDF5
  private
  public :: Galacticus_Output_Make_Group, Galacticus_Output_Dataset

  interface Galacticus_Output_Dataset
     !% Generic interface to the HDF5 dataset writing subroutines.
     module procedure Galacticus_Output_Dataset_Integer
     module procedure Galacticus_Output_Dataset_Double
     module procedure Galacticus_Output_Dataset_Character
     module procedure Galacticus_Output_Dataset_VarString
  end interface
  
contains
  
  subroutine Galacticus_Output_Dataset_Integer(locationID,datasetID,datasetName,commentText,datasetInteger,isExtendable)
    !% Write an integer dataset to the \glc\ output file.
    implicit none
    integer(kind=HID_T),   intent(in)           :: locationID
    integer(kind=HID_T),   intent(inout)        :: datasetID
    integer,               intent(in)           :: datasetInteger(:)
    character(len=*),      intent(in)           :: datasetName,commentText
    logical,               intent(in), optional :: isExtendable
    integer(kind=HID_T),   parameter            :: datasetRank=1
    integer(kind=HSIZE_T), parameter            :: chunkSize=1
    integer                                     :: errorCode
    integer(kind=HSIZE_T)                       :: datasetDimensions(1),datasetDimensionsMaximum(1),hyperslabStart(1),hyperslabCount(1)
    integer(kind=HID_T)                         :: dataspaceID,newDataspaceID,propertyList
    logical                                     :: isExtendableActual

    ! Check if location actually exists.
    if (locationID <= 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','location is null')
    ! Create the dataset if necessary.
    !$omp critical (HDF5_Operation)
    if (datasetID <= 0) then
       datasetDimensions=[size(datasetInteger)]
       if (present(isExtendable)) then
          isExtendableActual=isExtendable
       else
          isExtendableActual=.false.
       end if
       select case (isExtendableActual)
       case (.true.)
          datasetDimensionsMaximum=[H5S_UNLIMITED_F]
       case (.false.)
          datasetDimensionsMaximum=[size(datasetInteger)]
       end select
       call h5screate_simple_f(datasetRank,datasetDimensions,dataspaceID,errorCode,datasetDimensionsMaximum)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to make dataspace&
            & for integer dataset')
       call h5pcreate_f(H5P_DATASET_CREATE_F,propertyList,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to make property list &
            &for integer datset')
       datasetDimensions=[chunkSize]
       call h5pset_chunk_f(propertyList,datasetRank,datasetDimensions,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to set chunk size &
            &for integer dataset')
       call h5dcreate_f(locationID,trim(datasetName),H5T_NATIVE_INTEGER,dataspaceID,datasetID,errorCode,propertyList)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to make dataset &
            &for integer dataset')
       call h5sclose_f(dataspaceID,errorCode)   
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to close dataspace')
       call h5gset_comment_f(datasetID,'.',trim(commentText),errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to set comment &
            &for integer dataset')
       datasetDimensions=[size(datasetInteger)]
       hyperslabStart=[0]
       hyperslabCount=[size(datasetInteger)]
    else
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetID,dataspaceID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to get dataspace')
       call h5sget_simple_extent_dims_f(dataspaceID,datasetDimensions,datasetDimensionsMaximum,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to get dataspace extent')
       call h5sclose_f(dataspaceID,errorCode)   
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to close dataspace')
       hyperslabStart=datasetDimensions
       hyperslabCount=[size(datasetInteger)]
       datasetDimensions=datasetDimensions+[size(datasetInteger)]
    end if
    call h5dset_extent_f(datasetID,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to extend dataset')
    call h5dget_space_f(datasetID,dataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to get new dataspace')
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to select hyperslab')
    datasetDimensions=[size(datasetInteger)]
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to create temporary dataset')
    call h5dwrite_f(datasetID,H5T_NATIVE_INTEGER,datasetInteger,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to write dataset')
    call h5sclose_f(dataspaceID,errorCode)   
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to close dataspace')
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to close dataspace')
    !$omp end critical (HDF5_Operation)
    return
  end subroutine Galacticus_Output_Dataset_Integer

  subroutine Galacticus_Output_Dataset_Double(locationID,datasetID,datasetName,commentText,datasetDouble,isExtendable)
    !% Write an double precision dataset to the \glc\ output file. 
    implicit none
    integer(kind=HID_T),   intent(in)           :: locationID
    integer(kind=HID_T),   intent(inout)        :: datasetID
    double precision,      intent(in)           :: datasetDouble(:)
    character(len=*),      intent(in)           :: datasetName,commentText
    logical,               intent(in), optional :: isExtendable
    integer(kind=HID_T),   parameter            :: datasetRank=1
    integer(kind=HSIZE_T), parameter            :: chunkSize=1
    integer                                     :: errorCode
    integer(kind=HSIZE_T)                       :: datasetDimensions(1),datasetDimensionsMaximum(1),hyperslabStart(1),hyperslabCount(1)
    integer(kind=HID_T)                         :: dataspaceID,newDataspaceID,propertyList
    logical                                     :: isExtendableActual

    ! Check if location actually exists.
    if (locationID <= 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','location is null')

    ! Create the dataset if necessary.
    !$omp critical (HDF5_Operation)
    if (datasetID <= 0) then
       datasetDimensions=[size(datasetDouble)]
       if (present(isExtendable)) then
          isExtendableActual=isExtendable
       else
          isExtendableActual=.false.
       end if
       select case (isExtendableActual)
       case (.true.)
          datasetDimensionsMaximum=[H5S_UNLIMITED_F]
       case (.false.)
          datasetDimensionsMaximum=[size(datasetDouble)]
       end select
       call h5screate_simple_f(datasetRank,datasetDimensions,dataspaceID,errorCode,datasetDimensionsMaximum)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to make dataspace&
            & for double dataset')
       call h5pcreate_f(H5P_DATASET_CREATE_F,propertyList,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to make property list &
            &for double datset')
       datasetDimensions=[chunkSize]
       call h5pset_chunk_f(propertyList,datasetRank,datasetDimensions,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to set chunk size &
            &for double dataset')
       call h5dcreate_f(locationID,trim(datasetName),H5T_NATIVE_DOUBLE,dataspaceID,datasetID,errorCode,propertyList)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to make dataset &
            &for double dataset: '//trim(datasetName))
       call h5sclose_f(dataspaceID,errorCode)   
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to close dataspace')
       call h5gset_comment_f(datasetID,'.',trim(commentText),errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to set comment &
            &for double dataset')
       datasetDimensions=[size(datasetDouble)]
       hyperslabStart=[0]
       hyperslabCount=[size(datasetDouble)]
    else
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetID,dataspaceID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to get dataspace')
       call h5sget_simple_extent_dims_f(dataspaceID,datasetDimensions,datasetDimensionsMaximum,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to get dataspace extent')
       call h5sclose_f(dataspaceID,errorCode)   
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to close dataspace')
       hyperslabStart=datasetDimensions
       hyperslabCount=[size(datasetDouble)]
       datasetDimensions=datasetDimensions+[size(datasetDouble)]
    end if
    call h5dset_extent_f(datasetID,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to extend dataset')
    call h5dget_space_f(datasetID,dataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to get new dataspace')
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to select hyperslab')
    datasetDimensions=[size(datasetDouble)]
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to create temporary dataset')
    call h5dwrite_f(datasetID,H5T_NATIVE_DOUBLE,datasetDouble,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to write dataset')
    call h5sclose_f(dataspaceID,errorCode)   
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to close dataspace')
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to close dataspace')
    !$omp end critical (HDF5_Operation)
    return
  end subroutine Galacticus_Output_Dataset_Double

  subroutine Galacticus_Output_Dataset_Character(locationID,datasetID,datasetName,commentText,datasetCharacter,isExtendable)
    !% Write an character dataset to the \glc\ output file.
    implicit none
    integer(kind=HID_T),   intent(in)           :: locationID
    integer(kind=HID_T),   intent(inout)        :: datasetID
    character(len=*),      intent(in)           :: datasetCharacter(:)
    character(len=*),      intent(in)           :: datasetName,commentText
    logical,               intent(in), optional :: isExtendable
    integer(kind=HID_T),   parameter            :: datasetRank=1
    integer(kind=HSIZE_T), parameter            :: chunkSize=1
    integer                                     :: errorCode
    integer(kind=HSIZE_T)                       :: datasetDimensions(1),datasetDimensionsMaximum(1),hyperslabStart(1),hyperslabCount(1)
    integer(kind=HID_T)                         :: dataspaceID,newDataspaceID,propertyList,dataTypeID
    logical                                     :: isExtendableActual

    ! Check if location actually exists.
    if (locationID <= 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','location is null')

    ! Create a datatype
    !$omp critical (HDF5_Operation)
    call h5tcopy_f(H5T_NATIVE_CHARACTER,dataTypeID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Version_Output','failed to copy datatype')
    call h5tset_size_f(dataTypeID,int(len(datasetCharacter),size_t),errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Version_Output','failed to set datatype size')

    ! Create the dataset if necessary.
    if (datasetID <= 0) then
       datasetDimensions=[size(datasetCharacter)]
       if (present(isExtendable)) then
          isExtendableActual=isExtendable
       else
          isExtendableActual=.false.
       end if
       select case (isExtendableActual)
       case (.true.)
          datasetDimensionsMaximum=[H5S_UNLIMITED_F]
       case (.false.)
          datasetDimensionsMaximum=[size(datasetCharacter)]
       end select
       call h5screate_simple_f(datasetRank,datasetDimensions,dataspaceID,errorCode,datasetDimensionsMaximum)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to make dataspace&
            & for character dataset')
       call h5pcreate_f(H5P_DATASET_CREATE_F,propertyList,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to make property list &
            &for character datset')
       datasetDimensions=[chunkSize]
       call h5pset_chunk_f(propertyList,datasetRank,datasetDimensions,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to set chunk size &
            &for character dataset')
       call h5dcreate_f(locationID,trim(datasetName),datatypeID,dataspaceID,datasetID,errorCode,propertyList)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to make dataset &
            &for character dataset')
       call h5sclose_f(dataspaceID,errorCode)   
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to close dataspace')
       call h5gset_comment_f(datasetID,'.',trim(commentText),errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to set comment &
            &for character dataset')
       datasetDimensions=[size(datasetCharacter)]
       hyperslabStart=[0]
       hyperslabCount=[size(datasetCharacter)]
    else
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetID,dataspaceID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to get dataspace')
       call h5sget_simple_extent_dims_f(dataspaceID,datasetDimensions,datasetDimensionsMaximum,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to get dataspace extent')
       call h5sclose_f(dataspaceID,errorCode)   
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to close dataspace')
       hyperslabStart=datasetDimensions
       hyperslabCount=[size(datasetCharacter)]
       datasetDimensions=datasetDimensions+[size(datasetCharacter)]
    end if
    call h5dset_extent_f(datasetID,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to extend dataset')
    call h5dget_space_f(datasetID,dataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to get new dataspace')
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to select hyperslab')
    datasetDimensions=[size(datasetCharacter)]
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to create temporary dataset')
    call h5dwrite_f(datasetID,datatypeID,datasetCharacter,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to write dataset')
    call h5sclose_f(dataspaceID,errorCode)   
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to close dataspace')
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to close dataspace')
    !$omp end critical (HDF5_Operation)
    return
  end subroutine Galacticus_Output_Dataset_Character

  subroutine Galacticus_Output_Dataset_VarString(locationID,datasetID,datasetName,commentText,datasetVarString,isExtendable)
    !% Write an varying string dataset to the \glc\ output file.
    use ISO_Varying_String
    implicit none
    integer(kind=HID_T),   intent(in)            :: locationID
    integer(kind=HID_T),   intent(inout)         :: datasetID
    type(varying_string),  intent(in)            :: datasetVarString(:)
    character(len=*),      intent(in)            :: datasetName,commentText
    logical,               intent(in), optional  :: isExtendable
    integer(kind=HID_T),   parameter             :: datasetRank=1
    integer(kind=HSIZE_T), parameter             :: chunkSize=1
    integer                                      :: errorCode
    integer(kind=HSIZE_T)                        :: datasetDimensions(1),datasetDimensionsMaximum(1),hyperslabStart(1)&
         &,hyperslabCount(1)
    integer(kind=HID_T)                          :: dataspaceID,newDataspaceID,propertyList,dataTypeID
    logical                                      :: isExtendableActual
    character(len=maxval(len(datasetVarString))) :: datasetCharacter(size(datasetVarString))

    ! Check if location actually exists.
    if (locationID <= 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','location is null')

    ! Copy varying string data into temporary character array.
    datasetCharacter=datasetVarString

    ! Create a datatype
    !$omp critical (HDF5_Operation)
    call h5tcopy_f(H5T_NATIVE_CHARACTER,dataTypeID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Version_Output','failed to copy datatype')
    call h5tset_size_f(dataTypeID,int(len(datasetCharacter),size_t),errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Version_Output','failed to set datatype size')

    ! Create the dataset if necessary.
    if (datasetID <= 0) then
       datasetDimensions=[size(datasetCharacter)]
       if (present(isExtendable)) then
          isExtendableActual=isExtendable
       else
          isExtendableActual=.false.
       end if
       select case (isExtendableActual)
       case (.true.)
          datasetDimensionsMaximum=[H5S_UNLIMITED_F]
       case (.false.)
          datasetDimensionsMaximum=[size(datasetCharacter)]
       end select
       call h5screate_simple_f(datasetRank,datasetDimensions,dataspaceID,errorCode,datasetDimensionsMaximum)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to make dataspace&
            & for varying_string dataset')
       call h5pcreate_f(H5P_DATASET_CREATE_F,propertyList,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to make property list &
            &for varying_string datset')
       datasetDimensions=[chunkSize]
       call h5pset_chunk_f(propertyList,datasetRank,datasetDimensions,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to set chunk size &
            &for varying_string dataset')
       call h5dcreate_f(locationID,trim(datasetName),datatypeID,dataspaceID,datasetID,errorCode,propertyList)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to make dataset &
            &for varying_string dataset')
       call h5sclose_f(dataspaceID,errorCode)   
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to close dataspace')
       call h5gset_comment_f(datasetID,'.',trim(commentText),errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to set comment &
            &for varying_string dataset')
       datasetDimensions=[size(datasetCharacter)]
       hyperslabStart=[0]
       hyperslabCount=[size(datasetCharacter)]
    else
       ! Get size of existing dataset here.
       call h5dget_space_f(datasetID,dataspaceID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to get dataspace')
       call h5sget_simple_extent_dims_f(dataspaceID,datasetDimensions,datasetDimensionsMaximum,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to get dataspace extent')
       call h5sclose_f(dataspaceID,errorCode)   
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to close dataspace')
       hyperslabStart=datasetDimensions
       hyperslabCount=[size(datasetCharacter)]
       datasetDimensions=datasetDimensions+[size(datasetCharacter)]
    end if
    call h5dset_extent_f(datasetID,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to extend dataset')
    call h5dget_space_f(datasetID,dataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to get new dataspace')
    call h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,hyperslabStart,hyperslabCount,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to select hyperslab')
    datasetDimensions=[size(datasetCharacter)]
    call h5screate_simple_f(datasetRank,datasetDimensions,newDataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to create temporary dataset')
    call h5dwrite_f(datasetID,datatypeID,datasetCharacter,datasetDimensions,errorCode,newDataspaceID,dataspaceID)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to write dataset')
    call h5sclose_f(dataspaceID,errorCode)   
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to close dataspace')
    call h5sclose_f(newDataspaceID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to close dataspace')
    !$omp end critical (HDF5_Operation)
    return
  end subroutine Galacticus_Output_Dataset_VarString

  integer function Galacticus_Output_Make_Group(groupName,commentText,locationID)
    !% Create a group within the \glc\ output file.
    use OMP_Lib
    use String_Handling
    implicit none
    type(varying_string), intent(in)           :: groupName
    type(varying_string), intent(in), optional :: commentText
    integer,              intent(in), optional :: locationID
    type(varying_string)                       :: message
    integer                                    :: errorCode,locationIDActual

    ! Decide where to create the group.
    if (present(locationID)) then
       locationIDActual=locationID
    else
       locationIDActual=galacticusOutputID
    end if
    if (locationIDActual <= 0) call Galacticus_Error_Report('Galacticus_Output_Make_Group','location is null')
 
    ! Create a group.
    !$omp critical (HDF5_Operation)
    call h5gcreate_f(locationIDActual,char(groupName),Galacticus_Output_Make_Group,errorCode)
    if (errorCode < 0) then
       message='failed to make group "'//groupName//'" at '//locationID
       call Galacticus_Error_Report('Galacticus_Output_Make_Group',message)
    end if
    if (present(commentText)) then
       call h5gset_comment_f(Galacticus_Output_Make_Group,'.',char(commentText),errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Make_Group','failed to set comment')
    end if
    !$omp end critical (HDF5_Operation)
    return
  end function Galacticus_Output_Make_Group

end module Galacticus_HDF5_Groups
