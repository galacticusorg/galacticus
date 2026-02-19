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
Contains a module which implements a perfect hash algorithm for long integer keys.
!!}

module Hashes_Perfect
  !!{
  Implements a perfect hash algorithm for long integer keys based on methods described by \cite{czech_fundamental_1997}. The
  specific implementation follows the general structure of that given in a Dr. Dobbs
  \href{https://web.archive.org/web/20250613033605/https://www.drdobbs.com/architecture-and-design/generating-perfect-hash-functions/184404506}{article}.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  private
  public :: hashPerfect

  type hashPerfect
     !!{
     A derived type which stores perfect long integer hashes.
     !!}
     private
     logical                                      :: created , hasInverseTable                                                , hasValues
     integer(c_size_t)                            :: hashSize, rowSize
     integer(c_size_t), allocatable, dimension(:) :: r                         !  r(R)=amount row A(R,:) was shifted.
     integer(c_size_t), allocatable, dimension(:) :: C                         !  the shifted rows of A() collapse into C().
     integer(c_size_t), allocatable, dimension(:) :: v                         !  the values corresponding to the keys in C().
   contains
     !![
     <methods>
       <method description="Create a perfect hash." method="create" />
       <method description="Destroy a perfect hash." method="destroy" />
       <method description="Test if a key is present in a perfect hash." method="isPresent" />
       <method description="Return the value corresponding to a key in a perfect hash." method="value" />
       <method description="Return the index corresponding to a key in a perfect hash." method="index" />
       <method description="Return the size of a perfect hash." method="size" />
     </methods>
     !!]
     procedure :: create   =>Hash_Perfect_Create
     procedure :: destroy  =>Hash_Perfect_Destroy
     procedure :: index    =>Hash_Perfect_Index
     procedure :: isPresent=>Hash_Perfect_Is_Present
     procedure :: size     =>Hash_Perfect_Size
     procedure :: value    =>Hash_Perfect_Value
  end type hashPerfect

  type rowStructure
     !!{
     A row structure used in building hashes
     !!}
     integer(c_size_t) :: rowNumber    !  the row number in array A().
     integer(c_size_t) :: rowItemCount !  the number of items in this row of A().
  end type rowStructure

  ! A key value that is impossible
  integer, parameter :: invalidKey=-1

contains

  subroutine Hash_Perfect_Create(hash,keys,values,keepInverseTable)
    !!{
    Create a perfect hash for a given set of keys.
    !!}
    use :: Error            , only : Error_Report
    use :: Kind_Numbers     , only : kind_int8
    implicit none
    class  (hashPerfect   )                             , intent(inout)           :: hash
    integer(kind=kind_int8)             , dimension(:)  , intent(in   )           :: keys
    integer(kind=kind_int8)             , dimension(:)  , intent(in   ), optional :: values
    logical                                             , intent(in   ), optional :: keepInverseTable
    integer(c_size_t      ), allocatable, dimension(:,:)                          :: A                   !  A(i,j)=K (i=K/t, j=K mod t) for each key K.
    integer(c_size_t      ), allocatable, dimension(:,:)                          :: B                   !  B(i,j)=v (i=K/t, j=K mod t) for each key K.
    integer(kind=kind_int8), allocatable, dimension(:)                            :: resizeTemp
    ! row() exists to facilitate sorting the rows of A() by their "fullness".
    type   (rowStructure  ), allocatable, dimension(:)                            :: row                 !  Entry counts for the rows in A().
    integer(c_size_t      )                                                       :: hashTableMax    , i                                               , iColumn, iKey   , &
         &                                                                           iRow            , j                                               , k      , offset
    type   (rowStructure  )                                                       :: tmp

    ! Record options.
    hash%hasInverseTable=.true.
    if (present(keepInverseTable)) hash%hasInverseTable=keepInverseTable
    hash%hasValues=present(values)

    ! Compute values for hash%rowSize and hashTableMax.
    hashTableMax =maxval(keys)
    hash%rowSize=int(sqrt(dble(hashTableMax)))+1

    ! Allocate data structures.
    call hash%destroy()
    allocate   (A     (0:hash%rowSize-1,0:hash%rowSize-1))
    allocate   (hash%r(0:hash%rowSize-1                 ))
    allocate   (hash%C(0:hashTableMax-1                 ))
    allocate   (row   (0:hash%rowSize-1                 ))
    if (hash%hasValues) then
       allocate(B     (0:hash%rowSize-1,0:hash%rowSize-1))
       allocate(hash%v(0:hashTableMax-1                 ))
    else
       allocate(B     (0               ,0               ))
    end if

    ! Record that the hash has been created.
    hash%created=.true.

    ! Initialize data structures.
    ! A row offset may be 0, so the items in r() are set to a negative value to
    ! indicate that the offset for each row is not known yet.
    ! Every item in A() and C() is set to a value that is known to be an invalid
    ! key for the specific application.
    hash%r        =-1         ! Valid offsets are non-negative.
    hash%C        =invalidKey
    if (hash%hasValues) hash%v=0
    A             =invalidKey
    B             =0
    row%RowItemCount=0          ! Indicate that each row is empty.
    forall(iRow=0:hash%rowSize-1)
       row(iRow)%rowNumber=iRow ! Insert the row numbers.
    end forall

    ! The number of items in each row is also computed and stored in Row()%rowItemCount.
    do iKey=1,size(keys)
      iRow                  =    keys(iKey)/hash%rowSize
      iColumn               =mod(keys(iKey),hash%rowSize)
      A  (iRow,iColumn)     =keys(iKey)
      if (hash%hasValues) B(iRow,iColumn)=values(iKey)
      row(iRow)%rowItemCount=row(iRow)%rowItemCount+1
   end do

   ! The algorithm needs to know which row of A() is most full, 2nd most full,
   ! etc. This is most easily done by sorting an array of row-item-counts and
   ! remembering which row the item counts go with. That is what the array
   ! row() does for us.
   ! I saw no point in trying to be clever here, so a simple bubble sort is used.
   do i=0,hash%rowSize-2
      do j=i+1,hash%rowSize-1
         if (row(i)%rowItemCount < row(j)%rowItemCount) then
            tmp   =row(i)
            row(i)=row(j)
            row(j)=tmp
         end if
      end do
   end do

   ! Do the First-Fit Descending Method algorithm.
   ! For each non-empty row:
   !  1. shift the row right until none of its items collide with any of
   !     the items in previous rows.
   !  2. Record the shift amount in array hash%r().
   !  3. Insert this row into the hash table hash%C().
   i=0
   do while (row(i)%rowItemCount > 0)
      iRow=row(i)%rowNumber ! Get the next non-empty row.
      do offset=0,hashTableMax-hash%rowSize-2
         do k=0,hash%rowSize-1 ! Does this offset avoid collisions?
            if ((hash%C(offset+k) /= invalidKey) .and. (A(iRow,k) /= invalidKey)) exit
         end do
         if (k == hash%rowSize) then
            hash%r(iRow)=offset ! Record the shift amount for this row.
            do k=0,hash%rowSize-1 ! Insert this row into the hash table.
               if (A(iRow,k) /= invalidKey) then
                  hash%C(offset+k)=A(iRow,k)
                  if (hash%hasValues) hash%v(offset+k)=B(iRow,k)
               end if
            end do
            exit
         end if
      end do
      if (offset == hashTableMax-hash%rowSize-1) call Error_Report('failed to fit row into hash table - this should not happen'//{introspection:location})
      i=i+1
   end do
   ! Find the size of the resulting hash table.
   hash%hashSize=0_kind_int8
   do k=hashTableMax-1,0,-1
      if (hash%C(k) /= invalidKey) then
         hash%hashSize=k+1
         exit
      end if
   end do

   ! Deallocate data structures.
   deallocate(A  )
   deallocate(B  )
   deallocate(row)

   ! Reduce the size of stored arrays if possible.
   if (hash%hasInverseTable) then
      call move_alloc(hash%C,resizeTemp)
      allocate(hash%C(0:hash%hashSize-1))
      hash%C=resizeTemp(0:hash%hashSize-1)
      deallocate(resizeTemp)
   end if
   if (hash%hasValues) then
      call move_alloc(hash%v,resizeTemp)
      allocate(hash%v(0:hash%hashSize-1))
      hash%v=resizeTemp(0:hash%hashSize-1)
      deallocate(resizeTemp)
   end if

   ! Drop the inverse table if requested.
   if (.not.hash%hasInverseTable) deallocate(hash%C)
   return
 end subroutine Hash_Perfect_Create

 subroutine Hash_Perfect_Destroy(hash)
    !!{
    Destroy a perfect hash.
    !!}
    implicit none
    class(hashPerfect), intent(inout) :: hash

    hash%created=.false.
    if (allocated(hash%r)) deallocate(hash%r)
    if (allocated(hash%C)) deallocate(hash%C)
    if (allocated(hash%v)) deallocate(hash%v)
    return
  end subroutine Hash_Perfect_Destroy

 function Hash_Perfect_Size(hash)
   !!{
   Return the size of the hash table.
   !!}
   use :: Error, only : Error_Report
   implicit none
   integer(c_size_t   )                :: Hash_Perfect_Size
   class  (hashPerfect), intent(in   ) :: hash

   if (.not.hash%created) call Error_Report('hash has not been created'//{introspection:location})
   Hash_Perfect_Size=hash%hashSize
   return
 end function Hash_Perfect_Size

 function Hash_Perfect_Index(hash,key)
   !!{
   Return the index corresponding to a hash key.
   !!}
   use :: Error       , only : Error_Report
   use :: Kind_Numbers, only : kind_int8
   implicit none
   integer(c_size_t      )                :: Hash_Perfect_Index
   class  (hashPerfect   ), intent(in   ) :: hash
   integer(kind=kind_int8), intent(in   ) :: key
   integer(c_size_t      )                :: x                 , y

   if (.not.hash%created) call Error_Report('hash has not been created'//{introspection:location})
   x                 =    key/hash%rowSize
   y                 =mod(key,hash%rowSize)
   Hash_Perfect_Index=hash%r(x)+y
   return
 end function Hash_Perfect_Index

 logical function Hash_Perfect_Is_Present(hash,key)
   !!{
   Returns true if the hash contains the key.
   !!}
   use :: Error       , only : Error_Report
   use :: Kind_Numbers, only : kind_int8
   implicit none
   class  (hashPerfect   ), intent(in   ) :: hash
   integer(kind=kind_int8), intent(in   ) :: key
   integer(c_size_t      )                :: hashIndex

   hashIndex=hash%index(key)
   if (.not.hash%hasInverseTable) call Error_Report('hash does not store inverse table'//{introspection:location})
   Hash_Perfect_Is_Present=(hash%C(hashIndex) == key)
   return
 end function Hash_Perfect_Is_Present

 function Hash_Perfect_Value(hash,key)
   !!{
   Returns the value for a specified key.
   !!}
   use :: Error       , only : Error_Report
   use :: Kind_Numbers, only : kind_int8
   implicit none
   integer(kind=kind_int8)                :: Hash_Perfect_Value
   class  (hashPerfect   ), intent(in   ) :: hash
   integer(kind=kind_int8), intent(in   ) :: key
   integer(c_size_t      )                :: hashIndex

   hashIndex=hash%index(key)
   if (.not.hash%hasValues) call Error_Report('hash does not store values'//{introspection:location})
   Hash_Perfect_Value=hash%v(hashIndex)
   return
 end function Hash_Perfect_Value

end module Hashes_Perfect
