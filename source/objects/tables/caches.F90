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

!!{RST
Contains a module which implements persistent, mergeable caches of tabulated functions.
!!}

module Table_Caches
  !!{RST
  Implements persistent, mergeable caches of tabulated functions.

  Many tabulations in Galacticus are expensive enough that they are written to a file under ``pathTypeDataDynamic`` and re-read
  by later runs. Because such a table is built on an absolute lattice (see ``Numerical_Ranges``), the abscissae of the cached
  table and of any table built later are identical wherever they overlap, so a cached file whose range does not cover a new
  request need not be discarded: it can be merged with the table in memory, only the genuinely missing points computed, and the
  union written back. Cache files therefore grow monotonically towards the union of every range ever requested, rather than
  thrashing between runs.

  The intended usage is:

  .. code-block:: fortran

     lattice=Range_Pinned(x,pointsPer,scheme,latticeCurrent=table_%lattice)
     if (.not.table_%lattice%covers(lattice)) then
        ! Merge in anything already cached, then re-evaluate what is still needed.
        call Table_Cache_Restore(table_,fileName,status)
        lattice=Range_Pinned(x,pointsPer,scheme,latticeCurrent=table_%lattice)
        call table_%extend(lattice,isComputed)
        ! Compute only the missing points - note that this is done *without* holding any lock.
        do i=1,lattice%count
           if (isComputed(i)) cycle
           call table_%populate(...,i)
        end do
        ! Merge with, and write back, the cache file.
        call Table_Cache_Store(table_,fileName)
     end if

  Locking is handled internally: ``Table_Cache_Restore`` takes a shared lock, and ``Table_Cache_Store`` an exclusive lock which
  it holds only for the duration of the file input/output---never across the tabulation itself. ``Table_Cache_Store`` re-reads
  the file under that exclusive lock before writing, so that a range written by another process since our own read is merged
  rather than overwritten.
  !!}
  use :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: Table_Cache_Restore, Table_Cache_Store, Table_Cache_File_Name, Table_Cache_Lattice_Write, Table_Cache_Lattice_Read

  interface Table_Cache_Restore
     !!{RST
     Generic interface to routines which merge a cached tabulation into a table.
     !!}
     module procedure Table_Cache_Restore_1D
     module procedure Table_Cache_Restore_2D
  end interface Table_Cache_Restore

  interface Table_Cache_Store
     !!{RST
     Generic interface to routines which merge a table with, and write it to, a cache file.
     !!}
     module procedure Table_Cache_Store_1D
     module procedure Table_Cache_Store_2D
  end interface Table_Cache_Store

  ! Version of the cache file format. Files written with a different version are ignored (and so will be recomputed and
  ! overwritten), which allows the format to be changed without a migration.
  integer, parameter :: formatVersionCurrent=1

contains

  function Table_Cache_File_Name(subDirectory,objectType,hashedDescriptor,label) result(fileName)
    !!{RST
    Return the standard name of the cache file for a tabulation. Cache file names take the form

    .. code-block:: none

       {dynamic data path}/{subDirectory}/{objectType}{label}_{hashedDescriptor}.hdf5

    where ``hashedDescriptor`` should be obtained from
    ``hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)`` so that the name changes whenever
    either the parameters of the object or the source code which generates the tabulation changes. Using this function in
    preference to building the name by hand avoids the collisions which arise when a discriminator is accidentally omitted.
    !!}
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , assignment(=)
    implicit none
    type     (varying_string)                          :: fileName
    character(len=*         ), intent(in   )           :: subDirectory    , objectType, &
         &                                                hashedDescriptor
    character(len=*         ), intent(in   ), optional :: label
    type     (varying_string)                          :: label_

    label_=''
    if (present(label)) label_=label
    fileName=inputPath(pathTypeDataDynamic)//trim(subDirectory)//'/'//trim(objectType)//label_//'_'//trim(hashedDescriptor)//'.hdf5'
    return
  end function Table_Cache_File_Name

  subroutine Table_Cache_Lattice_Write(file,axisName,lattice)
    !!{RST
    Record, as attributes of the already-open file ``file``, the ``rangeLattice`` on which an axis of a stored tabulation is
    built. The attributes are named for that axis, so a file may hold the lattices of any number of axes.

    This is for tabulations held in raw arrays rather than in a :galacticus-class:`table1D`; those are written whole by
    ``Table_Cache_Store``, which records their lattices itself.
    !!}
    use :: IO_HDF5         , only : hdf5File
    use :: Numerical_Ranges, only : rangeLattice
    implicit none
    type     (hdf5File    ), intent(inout) :: file
    character(len=*       ), intent(in   ) :: axisName
    type     (rangeLattice), intent(in   ) :: lattice

    call file%writeAttribute(lattice%scheme%ID   ,axisName//'GridScheme'  )
    call file%writeAttribute(lattice%pointsPer   ,axisName//'PointsPer'   )
    call file%writeAttribute(lattice%indexMinimum,axisName//'IndexMinimum')
    call file%writeAttribute(lattice%count       ,axisName//'Count'       )
    return
  end subroutine Table_Cache_Lattice_Write

  subroutine Table_Cache_Lattice_Read(file,axisName,scheme,pointsPer,lattice)
    !!{RST
    Restore, from the already-open file ``file``, the ``rangeLattice`` on which an axis of a stored tabulation was built. The
    lattice is returned undefined---which the caller must treat as the stored tabulation being unusable---unless the file
    records one which is self-consistent and which uses the gridding scheme and density of points that the caller expects. A
    file written before the lattices were recorded, or with a different grid, therefore reports an undefined lattice rather
    than being misread.
    !!}
    use :: IO_HDF5         , only : hdf5File
    use :: Numerical_Ranges, only : rangeLattice, enumerationGridSchemeType
    implicit none
    type     (hdf5File                 ), intent(inout) :: file
    character(len=*                    ), intent(in   ) :: axisName
    type     (enumerationGridSchemeType), intent(in   ) :: scheme
    integer                             , intent(in   ) :: pointsPer
    type     (rangeLattice             ), intent(  out) :: lattice
    integer                                             :: schemeStored, pointsPerStored, &
         &                                                 indexMinimum, countPoints

    lattice=rangeLattice()
    if     (                                                  &
         &   .not.file%hasAttribute(axisName//'GridScheme'  ) &
         &  .or.                                              &
         &   .not.file%hasAttribute(axisName//'PointsPer'   ) &
         &  .or.                                              &
         &   .not.file%hasAttribute(axisName//'IndexMinimum') &
         &  .or.                                              &
         &   .not.file%hasAttribute(axisName//'Count'       ) &
         & ) return
    call file%readAttribute(axisName//'GridScheme'  ,schemeStored   )
    call file%readAttribute(axisName//'PointsPer'   ,pointsPerStored)
    call file%readAttribute(axisName//'IndexMinimum',indexMinimum   )
    call file%readAttribute(axisName//'Count'       ,countPoints    )
    ! Comparing the stored scheme against the one expected is stronger than merely checking that it is a valid member of the
    ! enumeration, so no separate validity test is needed.
    if (enumerationGridSchemeType(schemeStored) /= scheme   ) return
    if (pointsPerStored                         /= pointsPer) return
    lattice=rangeLattice(enumerationGridSchemeType(schemeStored),pointsPerStored,indexMinimum,countPoints)
    if (.not.lattice%isDefined()) lattice=rangeLattice()
    return
  end subroutine Table_Cache_Lattice_Read

  subroutine Table_Cache_Restore_1D(table_,fileName,status)
    !!{RST
    Merge any cached tabulation in ``fileName`` into ``table_``. A shared lock is taken on the file for the duration of the
    read. On return ``status`` is ``errorStatusSuccess`` if the cache contributed anything to ``table_``, and
    ``errorStatusFail`` otherwise (the file does not exist, is unreadable, was written in a different format, uses an
    incommensurate lattice, or adds nothing that ``table_`` does not already contain).
    !!}
    use :: Error             , only : errorStatusFail, errorStatusSuccess
    use :: File_Utilities    , only : File_Exists    , File_Lock         , File_Unlock, lockDescriptor
    use :: ISO_Varying_String, only : char
    use :: Tables            , only : table1D
    implicit none
    class  (table1D       ), allocatable, intent(inout) :: table_
    type   (varying_string)             , intent(in   ) :: fileName
    integer                             , intent(  out) :: status
    class  (table1D       ), allocatable                :: tableCached
    type   (lockDescriptor)                             :: fileLock

    status=errorStatusFail
    if (.not.File_Exists(fileName)) return
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock  (char(fileName),fileLock,lockIsShared=.true.)
    call Table_Cache_Read_1D(tableCached,fileName,status)
    call File_Unlock(fileLock)
    if (status /= errorStatusSuccess) return
    call Table_Cache_Merge_1D(table_,tableCached,status)
    return
  end subroutine Table_Cache_Restore_1D

  subroutine Table_Cache_Store_1D(table_,fileName)
    !!{RST
    Write ``table_`` to the cache file ``fileName``, first merging in anything which the file already contains. An exclusive
    lock is held for the duration, which - since the tabulation itself must already be complete - is only as long as the file
    input/output takes. Re-reading under the exclusive lock ensures that a range written by another process after our own
    ``Table_Cache_Restore`` is merged rather than overwritten. Note that ``table_`` may be extended by this call, since it
    returns holding the union of the cached and in-memory tabulations.
    !!}
    use :: Error             , only : Error_Report  , errorStatusSuccess
    use :: File_Utilities    , only : Directory_Make, File_Exists       , File_Lock  , File_Path, &
         &                            File_Unlock   , lockDescriptor
    use :: ISO_Varying_String, only : char
    use :: Tables            , only : table1D
    implicit none
    class  (table1D       ), allocatable, intent(inout) :: table_
    type   (varying_string)             , intent(in   ) :: fileName
    class  (table1D       ), allocatable                :: tableCached
    type   (lockDescriptor)                             :: fileLock
    integer                                             :: status

    if (.not.allocated(table_)             ) call Error_Report('no table to store'                          //{introspection:location})
    if (.not.table_%lattice%isDefined()    ) call Error_Report('table was not built on an absolute lattice' //{introspection:location})
    call Directory_Make(File_Path(fileName)                              )
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock     (char     (fileName),fileLock,lockIsShared=.false.)
    if (File_Exists(fileName)) then
       call Table_Cache_Read_1D(tableCached,fileName,status)
       if (status == errorStatusSuccess) call Table_Cache_Merge_1D(table_,tableCached,status)
    end if
    call Table_Cache_Write_1D(table_,fileName)
    call File_Unlock   (fileLock)
    return
  end subroutine Table_Cache_Store_1D

  subroutine Table_Cache_Read_1D(tableCached,fileName,status)
    !!{RST
    Read a cached tabulation from ``fileName`` into a newly allocated table. The caller is responsible for holding a suitable
    lock on the file.
    !!}
    use :: Error             , only : errorStatusFail             , errorStatusSuccess
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: Numerical_Ranges  , only : rangeLattice                , enumerationGridSchemeType
    use :: Tables            , only : table1D                     , table1DLinearLinear      , table1DLogarithmicLinear, tableTypeLinearLinear1D, &
         &                            tableTypeLogarithmicLinear1D
    implicit none
    class           (table1D       ), allocatable, intent(  out) :: tableCached
    type            (varying_string)             , intent(in   ) :: fileName
    integer                                      , intent(  out) :: status
    type            (hdf5File      )                             :: file
    type            (rangeLattice  )                             :: lattice
    logical                         , allocatable, dimension(:)  :: isComputed
    double precision                , allocatable, dimension(:,:):: valuesCached
    integer                                                      :: formatVersion, tableType   , &
         &                                                          gridScheme   , pointsPer   , &
         &                                                          indexMinimum , countPoints , &
         &                                                          countTables

    status=errorStatusFail
    !$ call hdf5Access%set()
    file=hdf5File(fileName,readOnly=.true.)
    if (file%hasAttribute('formatVersion')) then
       call file%readAttribute('formatVersion',formatVersion)
       if (formatVersion == formatVersionCurrent) then
          call file%readAttribute('tableType'   ,tableType   )
          call file%readAttribute('gridScheme'  ,gridScheme  )
          call file%readAttribute('pointsPer'   ,pointsPer   )
          call file%readAttribute('indexMinimum',indexMinimum)
          call file%readAttribute('count'       ,countPoints )
          call file%readAttribute('countTables' ,countTables )
          call file%readDataset  ('y'           ,valuesCached)
          lattice=rangeLattice(enumerationGridSchemeType(gridScheme),pointsPer,indexMinimum,countPoints)
          if     (                                                &
               &   lattice%isDefined()                            &
               &  .and.                                           &
               &   size(valuesCached,dim=1) == countPoints        &
               &  .and.                                           &
               &   size(valuesCached,dim=2) == countTables        &
               & ) then
             select case (tableType)
             case (tableTypeLinearLinear1D     %ID)
                allocate(table1DLinearLinear      :: tableCached)
             case (tableTypeLogarithmicLinear1D %ID)
                allocate(table1DLogarithmicLinear :: tableCached)
             end select
             if (allocated(tableCached)) then
                call tableCached%extend(lattice,isComputed,tableCount=countTables)
                tableCached%yv=valuesCached
                status        =errorStatusSuccess
             end if
          end if
       end if
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine Table_Cache_Read_1D

  subroutine Table_Cache_Write_1D(table_,fileName)
    !!{RST
    Write ``table_`` to ``fileName``. The caller is responsible for holding an exclusive lock on the file. The abscissae are
    written purely for the benefit of anyone inspecting the file - they are fully determined by the recorded lattice, and are
    reconstructed from it (not read) when the file is restored.
    !!}
    use :: Error      , only : Error_Report
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5File
    use :: Tables     , only : table1D     , table1DLinearLinear, table1DLogarithmicLinear, tableTypeLinearLinear1D, &
         &                     tableTypeLogarithmicLinear1D
    implicit none
    class  (table1D       ), allocatable, intent(in   ) :: table_
    type   (varying_string)             , intent(in   ) :: fileName
    type   (hdf5File      )                             :: file
    integer                                             :: tableType

    select type (table_)
    type is (table1DLogarithmicLinear)
       tableType=tableTypeLogarithmicLinear1D%ID
    type is (table1DLinearLinear     )
       tableType=tableTypeLinearLinear1D     %ID
    class default
       tableType=-1
       call Error_Report('caching is not supported for this table type'//{introspection:location})
    end select
    !$ call hdf5Access%set()
    file=hdf5File(fileName,overWrite=.true.,readOnly=.false.)
    call file%writeAttribute(formatVersionCurrent       ,'formatVersion')
    call file%writeAttribute(tableType                  ,'tableType'    )
    call file%writeAttribute(table_%lattice%scheme%ID   ,'gridScheme'   )
    call file%writeAttribute(table_%lattice%pointsPer   ,'pointsPer'    )
    call file%writeAttribute(table_%lattice%indexMinimum,'indexMinimum' )
    call file%writeAttribute(table_%lattice%count       ,'count'        )
    call file%writeAttribute(size(table_%yv,dim=2)      ,'countTables'  )
    call file%writeDataset  (table_%xs()                ,'x'            )
    call file%writeDataset  (table_%yv                  ,'y'            )
    !$ call hdf5Access%unset()
    return
  end subroutine Table_Cache_Write_1D

  subroutine Table_Cache_Merge_1D(table_,tableCached,status)
    !!{RST
    Merge the tabulation ``tableCached`` into ``table_``, returning ``errorStatusSuccess`` if anything was contributed.

    Since both tabulations lie on the same absolute lattice, merging is an exact splice of the two value arrays: no
    interpolation or resampling is involved, and any point present in both is identical in both. The cached tabulation is
    ignored (rather than merged) if it is built on an incommensurate lattice, or if the two ranges are disjoint with a gap
    between them---the union would then contain points which neither tabulation has computed, and this code has no means of
    computing them. In practice a gap cannot arise for a tabulation which has a default range, since both tabulations then
    contain that range.
    !!}
    use :: Display         , only : displayMessage, verbosityLevelWorking
    use :: Error           , only : errorStatusFail, errorStatusSuccess
    use :: Numerical_Ranges, only : rangeLattice
    use :: Tables          , only : table1D
    implicit none
    class           (table1D     ), allocatable, intent(inout) :: table_
    class           (table1D     ), allocatable, intent(inout) :: tableCached
    integer                                    , intent(  out) :: status
    type            (rangeLattice)                             :: lattice
    logical                       , allocatable, dimension(:)  :: isComputed
    double precision              , allocatable, dimension(:,:):: valuesCached
    integer                                                    :: offsetCached, i

    status=errorStatusFail
    if (.not.allocated(tableCached)) return
    ! If we have nothing in memory, simply adopt the cached tabulation.
    if (.not.allocated(table_)) then
       call Move_Alloc(tableCached,table_)
       status=errorStatusSuccess
       return
    end if
    if (.not.table_%lattice%isDefined()) then
       deallocate(table_)
       call Move_Alloc(tableCached,table_)
       status=errorStatusSuccess
       return
    end if
    ! Ignore a cached tabulation built on an incommensurate lattice.
    if     (                                                           &
         &   table_%lattice%scheme    /= tableCached%lattice%scheme    &
         &  .or.                                                       &
         &   table_%lattice%pointsPer /= tableCached%lattice%pointsPer &
         & ) then
       call displayMessage('ignoring a cached tabulation built on an incommensurate lattice',verbosityLevelWorking)
       return
    end if
    ! Nothing to contribute if we already contain everything cached.
    if (table_%lattice%covers(tableCached%lattice)) return
    ! Ignore a cached tabulation which neither overlaps nor abuts our own - the union would contain a gap.
    if     (                                                                           &
         &   tableCached%lattice%indexMinimum   > table_     %lattice%indexMaximum()+1 &
         &  .or.                                                                       &
         &   table_     %lattice%indexMinimum   > tableCached%lattice%indexMaximum()+1 &
         & ) then
       call displayMessage('ignoring a cached tabulation whose range is disjoint from that in memory',verbosityLevelWorking)
       return
    end if
    ! Extend to the union of the two ranges, then splice in those cached values which we do not already have.
    valuesCached       =tableCached%yv
    lattice            =table_     %lattice
    lattice%indexMinimum=min(table_%lattice%indexMinimum  ,tableCached%lattice%indexMinimum  )
    lattice%count       =max(table_%lattice%indexMaximum(),tableCached%lattice%indexMaximum())-lattice%indexMinimum+1
    call table_%extend(lattice,isComputed)
    offsetCached=tableCached%lattice%indexMinimum-lattice%indexMinimum
    do i=1,size(valuesCached,dim=1)
       if (.not.isComputed(offsetCached+i)) table_%yv(offsetCached+i,:)=valuesCached(i,:)
    end do
    status=errorStatusSuccess
    return
  end subroutine Table_Cache_Merge_1D

  subroutine Table_Cache_Restore_2D(table_,fileName,status)
    !!{RST
    Merge any cached tabulation in ``fileName`` into the two-dimensional table ``table_``, under a shared lock.

    Merging in two dimensions is more restrictive than in one: the union of two rectangles is a rectangle only if one contains
    the other, so unless the cached tabulation covers that in memory (in which case it is adopted) or is covered by it (in
    which case there is nothing to do), the cached tabulation is reported and ignored. Merging the general case would require
    the union to carry a record of which of its points have been computed, which tables do not presently keep.
    !!}
    use :: Display           , only : displayMessage , verbosityLevelWorking
    use :: Error             , only : errorStatusFail, errorStatusSuccess
    use :: File_Utilities    , only : File_Exists    , File_Lock            , File_Unlock, lockDescriptor
    use :: ISO_Varying_String, only : char
    use :: Tables            , only : table2DLogLogLin
    implicit none
    type   (table2DLogLogLin), intent(inout) :: table_
    type   (varying_string  ), intent(in   ) :: fileName
    integer                  , intent(  out) :: status
    type   (table2DLogLogLin)                :: tableCached
    type   (lockDescriptor  )                :: fileLock

    status=errorStatusFail
    if (.not.File_Exists(fileName)) return
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock  (char(fileName),fileLock,lockIsShared=.true.)
    call Table_Cache_Read_2D(tableCached,fileName,status)
    call File_Unlock(fileLock)
    if (status /= errorStatusSuccess) return
    call Table_Cache_Merge_2D(table_,tableCached,status)
    return
  end subroutine Table_Cache_Restore_2D

  subroutine Table_Cache_Store_2D(table_,fileName)
    !!{RST
    Write the two-dimensional table ``table_`` to the cache file ``fileName``, first merging in anything which the file already
    contains---see ``Table_Cache_Restore_2D`` for the limits of merging in two dimensions. An exclusive lock is held only for
    the duration of the file input/output.
    !!}
    use :: Error             , only : Error_Report  , errorStatusSuccess
    use :: File_Utilities    , only : Directory_Make, File_Exists       , File_Lock, File_Path, &
         &                            File_Unlock   , lockDescriptor
    use :: ISO_Varying_String, only : char
    use :: Tables            , only : table2DLogLogLin
    implicit none
    type   (table2DLogLogLin), intent(inout) :: table_
    type   (varying_string  ), intent(in   ) :: fileName
    type   (table2DLogLogLin)                :: tableCached
    type   (lockDescriptor  )                :: fileLock
    integer                                  :: status

    if (.not.table_%latticeX%isDefined() .or. .not.table_%latticeY%isDefined()) &
         & call Error_Report('table was not built on absolute lattices'//{introspection:location})
    call Directory_Make(File_Path(fileName)                              )
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock     (char     (fileName),fileLock,lockIsShared=.false.)
    if (File_Exists(fileName)) then
       call Table_Cache_Read_2D(tableCached,fileName,status)
       if (status == errorStatusSuccess) call Table_Cache_Merge_2D(table_,tableCached,status)
    end if
    call Table_Cache_Write_2D(table_,fileName)
    call File_Unlock   (fileLock)
    return
  end subroutine Table_Cache_Store_2D

  subroutine Table_Cache_Read_2D(tableCached,fileName,status)
    !!{RST
    Read a cached two-dimensional tabulation from ``fileName``. The caller is responsible for holding a suitable lock on the
    file.
    !!}
    use :: Error           , only : errorStatusFail , errorStatusSuccess
    use :: HDF5_Access     , only : hdf5Access
    use :: IO_HDF5         , only : hdf5File
    use :: Numerical_Ranges, only : rangeLattice    , enumerationGridSchemeType
    use :: Tables          , only : table2DLogLogLin, tableTypeLogLogLin2D
    implicit none
    type            (table2DLogLogLin)                               , intent(  out) :: tableCached
    type            (varying_string  )                               , intent(in   ) :: fileName
    integer                                                          , intent(  out) :: status
    type            (hdf5File        )                                               :: file
    type            (rangeLattice    )                                               :: latticeX     , latticeY
    logical                           , allocatable, dimension(:,:  )                :: isComputed
    double precision                  , allocatable, dimension(:,:,:)                :: valuesCached
    integer                                                                          :: formatVersion, tableType , &
         &                                                                              gridSchemeX  , pointsPerX, &
         &                                                                              indexMinimumX, countX    , &
         &                                                                              gridSchemeY  , pointsPerY, &
         &                                                                              indexMinimumY, countY    , &
         &                                                                              countTables

    status=errorStatusFail
    !$ call hdf5Access%set()
    file=hdf5File(fileName,readOnly=.true.)
    if (file%hasAttribute('formatVersion')) then
       call file%readAttribute('formatVersion',formatVersion)
       if (formatVersion == formatVersionCurrent) then
          call file%readAttribute('tableType'     ,tableType    )
          if (tableType == tableTypeLogLogLin2D%ID) then
             call file%readAttribute('gridSchemeX'   ,gridSchemeX  )
             call file%readAttribute('pointsPerX'    ,pointsPerX   )
             call file%readAttribute('indexMinimumX' ,indexMinimumX)
             call file%readAttribute('countX'        ,countX       )
             call file%readAttribute('gridSchemeY'   ,gridSchemeY  )
             call file%readAttribute('pointsPerY'    ,pointsPerY   )
             call file%readAttribute('indexMinimumY' ,indexMinimumY)
             call file%readAttribute('countY'        ,countY       )
             call file%readAttribute('countTables'   ,countTables  )
             call file%readDataset  ('z'             ,valuesCached )
             latticeX=rangeLattice(enumerationGridSchemeType(gridSchemeX),pointsPerX,indexMinimumX,countX)
             latticeY=rangeLattice(enumerationGridSchemeType(gridSchemeY),pointsPerY,indexMinimumY,countY)
             if     (                                                &
                  &   latticeX%isDefined()                           &
                  &  .and.                                           &
                  &   latticeY%isDefined()                           &
                  &  .and.                                           &
                  &   size(valuesCached,dim=1) == countX             &
                  &  .and.                                           &
                  &   size(valuesCached,dim=2) == countY             &
                  &  .and.                                           &
                  &   size(valuesCached,dim=3) == countTables        &
                  & ) then
                call tableCached%extend(latticeX,latticeY,isComputed,tableCount=countTables)
                tableCached%zv=valuesCached
                status        =errorStatusSuccess
             end if
          end if
       end if
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine Table_Cache_Read_2D

  subroutine Table_Cache_Write_2D(table_,fileName)
    !!{RST
    Write the two-dimensional table ``table_`` to ``fileName``. The caller is responsible for holding an exclusive lock on the
    file.
    !!}
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5File
    use :: Tables     , only : table2DLogLogLin, tableTypeLogLogLin2D
    implicit none
    type(table2DLogLogLin), intent(in   ) :: table_
    type(varying_string  ), intent(in   ) :: fileName
    type(hdf5File        )                :: file

    !$ call hdf5Access%set()
    file=hdf5File(fileName,overWrite=.true.,readOnly=.false.)
    call file%writeAttribute(formatVersionCurrent          ,'formatVersion')
    call file%writeAttribute(tableTypeLogLogLin2D%ID       ,'tableType'    )
    call file%writeAttribute(table_%latticeX%scheme%ID     ,'gridSchemeX'  )
    call file%writeAttribute(table_%latticeX%pointsPer     ,'pointsPerX'   )
    call file%writeAttribute(table_%latticeX%indexMinimum  ,'indexMinimumX')
    call file%writeAttribute(table_%latticeX%count         ,'countX'       )
    call file%writeAttribute(table_%latticeY%scheme%ID     ,'gridSchemeY'  )
    call file%writeAttribute(table_%latticeY%pointsPer     ,'pointsPerY'   )
    call file%writeAttribute(table_%latticeY%indexMinimum  ,'indexMinimumY')
    call file%writeAttribute(table_%latticeY%count         ,'countY'       )
    call file%writeAttribute(size(table_%zv,dim=3)         ,'countTables'  )
    call file%writeDataset  (table_%xs()                   ,'x'            )
    call file%writeDataset  (table_%ys()                   ,'y'            )
    call file%writeDataset  (table_%zv                     ,'z'            )
    !$ call hdf5Access%unset()
    return
  end subroutine Table_Cache_Write_2D

  subroutine Table_Cache_Merge_2D(table_,tableCached,status)
    !!{RST
    Merge the two-dimensional tabulation ``tableCached`` into ``table_``---see ``Table_Cache_Restore_2D`` for why this is
    restricted to the case in which one tabulation contains the other.
    !!}
    use :: Display         , only : displayMessage , verbosityLevelWorking
    use :: Error           , only : errorStatusFail, errorStatusSuccess
    use :: Tables          , only : table2DLogLogLin
    implicit none
    type   (table2DLogLogLin), intent(inout) :: table_
    type   (table2DLogLogLin), intent(inout) :: tableCached
    integer                  , intent(  out) :: status

    status=errorStatusFail
    ! If we have nothing in memory, simply adopt the cached tabulation.
    if (.not.table_%latticeX%isDefined() .or. .not.table_%latticeY%isDefined()) then
       table_=tableCached
       status=errorStatusSuccess
       return
    end if
    ! Ignore a cached tabulation built on incommensurate lattices.
    if     (                                                             &
         &   table_%latticeX%scheme    /= tableCached%latticeX%scheme    &
         &  .or.                                                         &
         &   table_%latticeY%scheme    /= tableCached%latticeY%scheme    &
         &  .or.                                                         &
         &   table_%latticeX%pointsPer /= tableCached%latticeX%pointsPer &
         &  .or.                                                         &
         &   table_%latticeY%pointsPer /= tableCached%latticeY%pointsPer &
         & ) then
       call displayMessage('ignoring a cached tabulation built on incommensurate lattices',verbosityLevelWorking)
       return
    end if
    ! Nothing to contribute if we already contain everything cached.
    if     (                                              &
         &   table_%latticeX%covers(tableCached%latticeX) &
         &  .and.                                         &
         &   table_%latticeY%covers(tableCached%latticeY) &
         & ) return
    ! Adopt the cached tabulation if it contains everything we have.
    if     (                                              &
         &   tableCached%latticeX%covers(table_%latticeX) &
         &  .and.                                         &
         &   tableCached%latticeY%covers(table_%latticeY) &
         & ) then
       table_=tableCached
       status=errorStatusSuccess
       return
    end if
    call displayMessage('ignoring a cached tabulation which neither contains, nor is contained by, that in memory',verbosityLevelWorking)
    return
  end subroutine Table_Cache_Merge_2D

end module Table_Caches
