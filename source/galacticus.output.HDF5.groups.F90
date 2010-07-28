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


!% Contains a module which handles creation of groups and writing of datasets to the \glc\ output file.

module Galacticus_HDF5_Groups
  !% Handles creation of groups and writing of datasets to the \glc\ output file.
  use Galacticus_Error
  use Galacticus_HDF5
  use ISO_Varying_String
  use HDF5
  private
  public :: Galacticus_Output_Make_Group, Galacticus_Output_Close_Group, Galacticus_Output_Dataset

  interface Galacticus_Output_Dataset
     !% Generic interface to the HDF5 dataset writing subroutines.
     module procedure Galacticus_Output_Dataset_Integer
     module procedure Galacticus_Output_Dataset_Double
     module procedure Galacticus_Output_Dataset_Character
     module procedure Galacticus_Output_Dataset_VarString
  end interface

  ! Flag indicating if module is initialized.
  logical                          :: galacticusHDF5GroupsInitialized=.false.
  
contains
  
  subroutine Galacticus_Output_Dataset_Integer(locationID,datasetID,datasetName,commentText,datasetInteger,isExtendable,isNew)
    !% Write an integer dataset to the \glc\ output file.
    implicit none
    integer(kind=HID_T),   intent(in)              :: locationID
    integer(kind=HID_T),   intent(inout)           :: datasetID
    integer,               intent(in)              :: datasetInteger(:)
    character(len=*),      intent(in)              :: datasetName,commentText
    logical,               intent(in),    optional :: isExtendable
    logical,               intent(inout), optional :: isNew
    integer(kind=HID_T),   parameter               :: datasetRank=1
    integer                                        :: errorCode
    integer(kind=HSIZE_T)                          :: datasetDimensions(1),datasetDimensionsMaximum(1),hyperslabStart(1),hyperslabCount(1)
    integer(kind=HID_T)                            :: dataspaceID,newDataspaceID,propertyList
    logical                                        :: isExtendableActual,isNewActual

    ! Check if location actually exists.
    if (locationID <= 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','location is null')

    ! See if the dataset is new.
    if (present(isNew)) then
       isNewActual=isNew
    else
       isNewActual=.true.
    end if

    ! Determine if dataset is extendable.
    if (present(isExtendable)) then
       isExtendableActual=isExtendable
    else
       isExtendableActual=.false.
    end if

    ! Non-new datasets must be extensible.
    if (.not.isNewActual.and..not.isExtendableActual) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','non-new datasets must be extendable')
   
    ! Create the dataset if necessary.
    !$omp critical (HDF5_Operation)
    if (isNewActual) then
       datasetDimensions=[size(datasetInteger)]
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
       datasetDimensions=[min(size(datasetInteger),hdf5ChunkSize)]
       call h5pset_chunk_f(propertyList,datasetRank,datasetDimensions,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to set chunk size &
            &for integer dataset')
       if (hdf5CompressionLevel >= 0) then
          call h5pset_deflate_f(propertyList,hdf5CompressionLevel,errorCode) 
          if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to set compression level &
               &for integer dataset')
       end if
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
       call h5dopen_f(locationID,trim(datasetName),datasetID,errorCode)
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
    call h5dclose_f(datasetID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to close dataset')
    if (isNewActual) then
       call h5pclose_f(propertyList,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Integer','failed to close property list')
    end if
    ! Dataset is no longer new.
    if (present(isNew)) isNew=.false.
    !$omp end critical (HDF5_Operation)
    return
  end subroutine Galacticus_Output_Dataset_Integer

  subroutine Galacticus_Output_Dataset_Double(locationID,datasetID,datasetName,commentText,datasetDouble,isExtendable,isNew)
    !% Write an double precision dataset to the \glc\ output file. 
    implicit none
    integer(kind=HID_T),   intent(in)              :: locationID
    integer(kind=HID_T),   intent(inout)           :: datasetID
    double precision,      intent(in)              :: datasetDouble(:)
    character(len=*),      intent(in)              :: datasetName,commentText
    logical,               intent(in),    optional :: isExtendable
    logical,               intent(inout), optional :: isNew
    integer(kind=HID_T),   parameter               :: datasetRank=1
    integer                                        :: errorCode
    integer(kind=HSIZE_T)                          :: datasetDimensions(1),datasetDimensionsMaximum(1),hyperslabStart(1),hyperslabCount(1)
    integer(kind=HID_T)                            :: dataspaceID,newDataspaceID,propertyList
    logical                                        :: isExtendableActual,isNewActual

    ! Check if location actually exists.
    if (locationID <= 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','location is null')

    ! See if the dataset is new.
    if (present(isNew)) then
       isNewActual=isNew
    else
       isNewActual=.true.
    end if

    ! Determine if dataset is extendable.
    if (present(isExtendable)) then
       isExtendableActual=isExtendable
    else
       isExtendableActual=.false.
    end if

    ! Non-new datasets must be extensible.
    if (.not.isNewActual.and..not.isExtendableActual) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','non-new datasets must be extendable')
   
    ! Create the dataset if necessary.
    !$omp critical (HDF5_Operation)
    if (isNewActual) then
       datasetDimensions=[size(datasetDouble)]
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
       datasetDimensions=[min(size(datasetDouble),hdf5ChunkSize)]
       call h5pset_chunk_f(propertyList,datasetRank,datasetDimensions,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to set chunk size &
            &for double dataset')
       if (hdf5CompressionLevel >= 0) then
          call h5pset_deflate_f(propertyList,hdf5CompressionLevel,errorCode) 
          if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to set compression level &
               &for double dataset')
       end if
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
       call h5dopen_f(locationID,trim(datasetName),datasetID,errorCode)
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
    call h5dclose_f(datasetID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to close dataset')
    if (isNewActual) then
       call h5pclose_f(propertyList,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Double','failed to close property list')
    end if
    ! Dataset is no longer new.
    if (present(isNew)) isNew=.false.
    !$omp end critical (HDF5_Operation)
    return
  end subroutine Galacticus_Output_Dataset_Double

  subroutine Galacticus_Output_Dataset_Character(locationID,datasetID,datasetName,commentText,datasetCharacter,isExtendable,isNew)
    !% Write an character dataset to the \glc\ output file.
    implicit none
    integer(kind=HID_T),   intent(in)              :: locationID
    integer(kind=HID_T),   intent(inout)           :: datasetID
    character(len=*),      intent(in)              :: datasetCharacter(:)
    character(len=*),      intent(in)              :: datasetName,commentText
    logical,               intent(in),    optional :: isExtendable
    logical,               intent(inout), optional :: isNew
    integer(kind=HID_T),   parameter               :: datasetRank=1
    integer                                        :: errorCode
    integer(kind=HSIZE_T)                          :: datasetDimensions(1),datasetDimensionsMaximum(1),hyperslabStart(1),hyperslabCount(1)
    integer(kind=HID_T)                            :: dataspaceID,newDataspaceID,propertyList,dataTypeID
    logical                                        :: isExtendableActual,isNewActual

    ! Check if location actually exists.
    if (locationID <= 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','location is null')

    ! See if the dataset is new.
    if (present(isNew)) then
       isNewActual=isNew
    else
       isNewActual=.true.
    end if

    ! Determine if dataset is extendable.
    if (present(isExtendable)) then
       isExtendableActual=isExtendable
    else
       isExtendableActual=.false.
    end if
    
    ! Non-new datasets must be extensible.
    if (.not.isNewActual.and..not.isExtendableActual) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','non-new datasets must be extendable')
   
    ! Create a datatype
    !$omp critical (HDF5_Operation)
    call h5tcopy_f(H5T_NATIVE_CHARACTER,dataTypeID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Version_Output','failed to copy datatype')
    call h5tset_size_f(dataTypeID,int(len(datasetCharacter),size_t),errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Version_Output','failed to set datatype size')

    ! Create the dataset if necessary.
    if (isNewActual) then
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
       datasetDimensions=[min(size(datasetCharacter),hdf5ChunkSize)]
       call h5pset_chunk_f(propertyList,datasetRank,datasetDimensions,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to set chunk size &
            &for character dataset')
       if (hdf5CompressionLevel >= 0) then
          call h5pset_deflate_f(propertyList,hdf5CompressionLevel,errorCode) 
          if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to set compression level &
               &for character dataset')
       end if
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
       call h5dopen_f(locationID,trim(datasetName),datasetID,errorCode)
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
    call h5dclose_f(datasetID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to close dataset')
    if (isNewActual) then
       call h5pclose_f(propertyList,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_Character','failed to close property list')
    end if
    ! Dataset is no longer new.
    if (present(isNew)) isNew=.false.
    !$omp end critical (HDF5_Operation)
    return
  end subroutine Galacticus_Output_Dataset_Character

  subroutine Galacticus_Output_Dataset_VarString(locationID,datasetID,datasetName,commentText,datasetVarString,isExtendable,isNew)
    !% Write an varying string dataset to the \glc\ output file.
    use ISO_Varying_String
    implicit none
    integer(kind=HID_T),   intent(in)              :: locationID
    integer(kind=HID_T),   intent(inout)           :: datasetID
    type(varying_string),  intent(in)              :: datasetVarString(:)
    character(len=*),      intent(in)              :: datasetName,commentText
    logical,               intent(in),    optional :: isExtendable
    logical,               intent(inout), optional :: isNew
    integer(kind=HID_T),   parameter               :: datasetRank=1
    integer                                        :: errorCode
    integer(kind=HSIZE_T)                          :: datasetDimensions(1),datasetDimensionsMaximum(1),hyperslabStart(1)&
         &,hyperslabCount(1)
    integer(kind=HID_T)                            :: dataspaceID,newDataspaceID,propertyList,dataTypeID
    logical                                        :: isExtendableActual,isNewActual
    character(len=maxval(len(datasetVarString)))   :: datasetCharacter(size(datasetVarString))

    ! Check if location actually exists.
    if (locationID <= 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','location is null')

    ! Copy varying string data into temporary character array.
    datasetCharacter=datasetVarString

    ! See if the dataset is new.
    if (present(isNew)) then
       isNewActual=isNew
    else
       isNewActual=.true.
    end if

    ! Determine if dataset is extendable.
    if (present(isExtendable)) then
       isExtendableActual=isExtendable
    else
       isExtendableActual=.false.
    end if
    
    ! Non-new datasets must be extensible.
    if (.not.isNewActual.and..not.isExtendableActual) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','non-new datasets must be extendable')
   
    ! Create a datatype
    !$omp critical (HDF5_Operation)
    call h5tcopy_f(H5T_NATIVE_CHARACTER,dataTypeID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to copy datatype')
    call h5tset_size_f(dataTypeID,int(len(datasetCharacter),size_t),errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to set datatype size')

    ! Create the dataset if necessary.
    if (isNewActual) then
       datasetDimensions=[size(datasetCharacter)]
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
       datasetDimensions=[min(size(datasetCharacter),hdf5ChunkSize)]
       call h5pset_chunk_f(propertyList,datasetRank,datasetDimensions,errorCode) 
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to set chunk size &
            &for varying_string dataset')
       if (hdf5CompressionLevel >= 0) then
          call h5pset_deflate_f(propertyList,hdf5CompressionLevel,errorCode) 
          if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to set compression level &
               &for varying_string dataset')
       end if
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
       call h5dopen_f(locationID,trim(datasetName),datasetID,errorCode)
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
    call h5dclose_f(datasetID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to close dataset')
    if (isNewActual) then
       call h5pclose_f(propertyList,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Dataset_VarString','failed to close property list')
    end if
    ! Dataset is no longer new.
    if (present(isNew)) isNew=.false.
    !$omp end critical (HDF5_Operation)
    return
  end subroutine Galacticus_Output_Dataset_VarString

  function Galacticus_Output_Make_Group(groupName,commentText,locationID)
    !% Create a group within the \glc\ output file.
    use OMP_Lib
    use String_Handling
    implicit none
    integer(kind=HID_T)                        :: Galacticus_Output_Make_Group
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

  subroutine Galacticus_Output_Close_Group(groupID)
    !% Closes an output group.
    implicit none
    integer(kind=HID_T), intent(in) :: groupID
    integer                         :: errorCode

    call h5gclose_f(groupID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Galacticus_Output_Close_Group','failed to close group')
    return
  end subroutine Galacticus_Output_Close_Group

end module Galacticus_HDF5_Groups
