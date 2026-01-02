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
Contains a module defining the history object type.
!!}

module Histories
  !!{
  Defines the history object type.
  !!}
  use :: Kind_Numbers, only : kind_int8
  implicit none
  private
  public :: history, longIntegerHistory, operator(*)

  ! Interface to multiplication operators with history objects as their second argument.
  interface operator(*)
     module procedure History_Multiply_Switched
  end interface operator(*)

  type history
     !!{
     The history object type.
     !!}
     double precision, allocatable, dimension(:  ) :: time
     double precision, allocatable, dimension(:,:) :: data
     integer                                       :: rangeType
   contains
     !![
     <methods>
       <method description="Addition operator." method="operator(+)" />
       <method description="Subtraction operator." method="operator(-)" />
       <method description="Division operator." method="operator(/)" />
       <method description="Multiplication operator." method="operator(*)" />
       <method description="Returns true if the history is entirely zero." method="isZero" />
       <method description="Creates a history object with a specified range of times." method="create" />
       <method description="Build a history object from an XML definition." method="builder" />
       <method description="Dump a history object." method="dump" />
       <method description="Dump a history object in binary." method="dumpRaw" />
       <method description="Read a history object in binary." method="readRaw" />
       <method description="Clone a history object." method="clone" />
       <method description="Destroys a history object." method="destroy" />
       <method description="Removes any times in a history which have become outdated." method="trim" />
       <method description="Removes any times in a history \emph{after} the given time. Optionally returns a history object with the removed history." method="trimForward" />
       <method description="Adds two histories, possibly with different time series." method="increment" />
       <method description="Adds two histories, possibly with different time series, by interpolating the second onto the times of the first and adding the interpolated values." method="interpolatedIncrement" />
       <method description="Extends the time range of a history to encompass the specified limits." method="extend" />
       <method description="Resets all entries in a history to zero." method="reset" />
       <method description="Set all entries in a history to unity." method="setToUnity" />
       <method description="Returns true if the given history has been created." method="exists" />
       <method description="Returns an array with the timesteps (i.e. the intervals between successive times) in the given history." method="timeSteps" />
       <method description="Return a count of the number of properties in a serialized history object." method="serializeCount" />
       <method description="Serialize a history object to an array." method="serialize" />
       <method description="Deserialize a history object from an array." method="deserialize" />
       <method description="Append a history or single instant onto the end of a history." method="append" />
       <method description="Returns the size of any non-static components of the type." method="nonStaticSizeOf" />
     </methods>
     !!]
     procedure ::                          History_Add
     procedure ::                          History_Subtract
     procedure ::                          History_Divide
     procedure ::                          History_Multiply
     generic   :: operator(+)           => History_Add
     generic   :: operator(-)           => History_Subtract
     generic   :: operator(/)           => History_Divide
     generic   :: operator(*)           => History_Multiply
     procedure :: nonStaticSizeOf       => History_Non_Static_Size_Of
     procedure :: isZero                => History_Is_Zero
     procedure :: builder               => History_Builder
     procedure :: dump                  => History_Dump
     procedure :: dumpRaw               => History_Dump_Raw
     procedure :: readRaw               => History_Read_Raw
     procedure :: create                => History_Create
     procedure :: clone                 => History_Clone
     procedure :: destroy               => History_Destroy
     procedure :: trim                  => History_Trim
     procedure :: trimForward           => History_Trim_Forward
     procedure :: extend                => History_Extend
     procedure :: increment             => History_Increment
     procedure :: interpolatedIncrement => History_Interpolated_Increment
     procedure :: reset                 => History_Reset
     procedure :: setToUnity            => History_Set_To_Unity
     procedure :: exists                => History_Exists
     procedure :: timeSteps             => History_Timesteps
     procedure :: serializeCount        => History_Serialize_Count
     procedure :: serialize             => History_Serialize
     procedure :: deserialize           => History_Deserialize
     procedure ::                          History_Append_History
     procedure ::                          History_Append_Epoch
     generic   :: append                => History_Append_History, &
          &                                History_Append_Epoch
  end type history

  type longIntegerHistory
     !!{
     The history object type.
     !!}
     double precision                , allocatable, dimension(:  ) :: time
     integer         (kind=kind_int8), allocatable, dimension(:,:) :: data
     integer                                                       :: rangeType
   contains
     !![
     <methods>
       <method description="Creates a history object with a specified range of times." method="create" />
       <method description="Build a history object from an XML definition." method="builder" />
       <method description="Dump a history object." method="dump" />
       <method description="Dump a history object in binary." method="dumpRaw" />
       <method description="Read a history object in binary." method="readRaw" />
       <method description="Clone a history object." method="clone" />
       <method description="Destroys a history object." method="destroy" />
       <method description="Removes any times in a history which have become outdated." method="trim" />
       <method description="Removes any times in a history \emph{after} the given time. Optionally returns a history object with the removed history." method="trimForward" />
       <method description="Resets all entries in a history to zero." method="reset" />
       <method description="Returns true if the given history has been created." method="exists" />
       <method description="Append a history or single instant onto the end of a history." method="append" />
       <method description="Returns the size of any non-static components of the type." method="nonStaticSizeOf" />
     </methods>
     !!]
     procedure :: builder         => History_Long_Integer_Builder
     procedure :: dump            => History_Long_Integer_Dump
     procedure :: dumpRaw         => History_Long_Integer_Dump_Raw
     procedure :: readRaw         => History_Long_Integer_Read_Raw
     procedure :: create          => History_Long_Integer_Create
     procedure :: clone           => History_Long_Integer_Clone
     procedure :: destroy         => History_Long_Integer_Destroy
     procedure :: trim            => History_Long_Integer_Trim
     procedure :: trimForward     => History_Long_Integer_Trim_Forward
     procedure :: reset           => History_Long_Integer_Reset
     procedure :: exists          => History_Long_Integer_Exists
     procedure ::                    History_Long_Integer_Append_History
     procedure ::                    History_Long_Integer_Append_Epoch
     generic   :: append          => History_Long_Integer_Append_History, &
          &                          History_Long_Integer_Append_Epoch
     procedure :: nonStaticSizeOf => History_Long_Integer_Non_Static_Size_Of
  end type longIntegerHistory

  ! A null history object.
  type            (history)           , public :: nullHistory

  ! Labels for targets when adding to histories.
  integer                  , parameter, public :: historyData               =1
  integer                  , parameter, public :: historyRates              =2
  integer                  , parameter, public :: historyScales             =3

contains

  subroutine History_Create(history_,historyCount,timesCount,timeBegin,timeEnd,rangeType)
    !!{
    Create a history object.
    !!}
    use :: Error            , only : Error_Report
    use :: Numerical_Ranges , only : Make_Range         , rangeTypeLogarithmic, rangeTypeUndefined
    implicit none
    class           (history), intent(inout)           :: history_
    integer                  , intent(in   )           :: historyCount   , timesCount
    double precision         , intent(in   ), optional :: timeBegin      , timeEnd
    integer                  , intent(in   ), optional :: rangeType
    integer                                            :: rangeTypeActual

    if (allocated(history_%time)) then
       call Error_Report('this history appears to have been created already'//{introspection:location})
    else
       allocate(history_%time  (timesCount             ))
       allocate(history_%data  (timesCount,historyCount))
       if (timesCount > 0) then
          if (present(timeBegin)) then
             if (.not.present(timeEnd)) call Error_Report('an end time must be given if a begin time is given'//{introspection:location})
             if (present(rangeType)) then
                rangeTypeActual=rangeType
             else
                rangeTypeActual=rangeTypeLogarithmic
             end if
             history_%time     =Make_Range(timeBegin,timeEnd,timesCount,rangeTypeActual)
             history_%rangeType=rangeTypeActual
          else
             history_%rangeType=rangeTypeUndefined
             history_%time=0.0d0
          end if
          if (historyCount > 0) history_%data=0.0d0
       end if
    end if
    return
  end subroutine History_Create

  subroutine History_Destroy(history_)
    !!{
    Destroy a history.
    !!}
    implicit none
    class  (history), intent(inout) :: history_

    if (allocated(history_%time)) then
       deallocate(history_%time)
       if (allocated(history_%data)) deallocate(history_%data)
    end if
    return
  end subroutine History_Destroy

  subroutine History_Long_Integer_Create(history_,historyCount,timesCount,timeBegin,timeEnd,rangeType)
    !!{
    Create a history object.
    !!}
    use :: Error            , only : Error_Report
    use :: Numerical_Ranges , only : Make_Range         , rangeTypeLogarithmic, rangeTypeUndefined
    implicit none
    class           (longIntegerHistory), intent(inout)           :: history_
    integer                             , intent(in   )           :: historyCount   , timesCount
    double precision                    , intent(in   ), optional :: timeBegin      , timeEnd
    integer                             , intent(in   ), optional :: rangeType
    integer                                                       :: rangeTypeActual

    if (allocated(history_%time)) then
       call Error_Report('this history appears to have been created already'//{introspection:location})
    else
       allocate(history_%time  (timesCount             ))
       allocate(history_%data  (timesCount,historyCount))
       if (timesCount > 0) then
          if (present(timeBegin)) then
             if (.not.present(timeEnd)) call Error_Report('an end time must be given if a begin time is given'//{introspection:location})
             if (present(rangeType)) then
                rangeTypeActual=rangeType
             else
                rangeTypeActual=rangeTypeLogarithmic
             end if
             history_%time     =Make_Range(timeBegin,timeEnd,timesCount,rangeTypeActual)
             history_%rangeType=rangeTypeActual
          else
             history_%rangeType=rangeTypeUndefined
             history_%time=0.0d0
          end if
          if (historyCount > 0) history_%data=0_kind_int8
       end if
    end if
    return
  end subroutine History_Long_Integer_Create

  subroutine History_Long_Integer_Destroy(history_)
    !!{
    Destroy a history.
    !!}
    implicit none
    class  (longIntegerHistory), intent(inout) :: history_

    if (allocated(history_%time)) then
       deallocate(history_%time)
       if (allocated(history_%data)) deallocate(history_%data)
    end if
    return
  end subroutine History_Long_Integer_Destroy

  subroutine History_Builder(self,historyDefinition)
    !!{
    Build a {\normalfont \ttfamily history} object from the given XML {\normalfont \ttfamily historyDefinition}.
    !!}
    use :: FoX_DOM, only : node
    use :: Error  , only : Error_Report
    implicit none
    class(history), intent(inout) :: self
    type (node   ), pointer       :: historyDefinition
    !$GLC attributes unused :: self, historyDefinition

    call Error_Report('building of history objects is not yet supported'//{introspection:location})
    return
  end subroutine History_Builder

  subroutine History_Long_Integer_Builder(self,historyDefinition)
    !!{
    Build a {\normalfont \ttfamily longIntegerHistory} object from the given XML {\normalfont \ttfamily historyDefinition}.
    !!}
    use :: FoX_DOM, only : node
    use :: Error  , only : Error_Report
    implicit none
    class(longIntegerHistory), intent(inout) :: self
    type (node              ), pointer       :: historyDefinition
    !$GLC attributes unused :: self, historyDefinition

    call Error_Report('building of history objects is not yet supported'//{introspection:location})
    return
  end subroutine History_Long_Integer_Builder

  subroutine History_Dump(self,verbosityLevel)
    !!{
    Dumps a history object.
    !!}
    use :: Display           , only : displayMessage, enumerationVerbosityLevelType
    use :: ISO_Varying_String, only : assignment(=) , operator(//)                 , varying_string
    implicit none
    class    (history                      ), intent(in   ) :: self
    type     (enumerationVerbosityLevelType), intent(in   ) :: verbosityLevel
    integer                                                 :: i             , j
    type     (varying_string               )                :: message
    character(len=22                       )                :: label

    if (allocated(self%time)) then
       do i=1,size(self%time)
          write (label,'(i3)') i
          message="("//trim(label)//") "
          write (label,'(e22.16)') self%time(i)
          message=message//label//" :"
          do j=1,size(self%data,dim=2)
             write (label,'(e22.16)') self%data(i,j)
             message=message//" "//label
          end do
          call displayMessage(message,verbosityLevel)
       end do
    end if
    return
  end subroutine History_Dump

  subroutine History_Dump_Raw(self,fileHandle)
    !!{
    Dumps a history object in binary.
    !!}
    implicit none
    class  (history), intent(in   ) :: self
    integer         , intent(in   ) :: fileHandle

    write (fileHandle) self%rangeType
    write (fileHandle) allocated(self%time)
    if (allocated(self%time)) then
       write (fileHandle) shape(self%data)
       write (fileHandle) self%time
       write (fileHandle) self%data
    end if
    return
  end subroutine History_Dump_Raw

  subroutine History_Read_Raw(self,fileHandle)
    !!{
    Read a history object in binary.
    !!}
    implicit none
    class  (history), intent(inout) :: self
    integer         , intent(in   ) :: fileHandle
    logical                         :: isAllocated
    integer         , dimension(2)  :: historyShape

    read (fileHandle) self%rangeType
    read (fileHandle) isAllocated
    if (isAllocated) then
       read (fileHandle) historyShape
       allocate(self%time(historyShape(1)                ))
       allocate(self%data(historyShape(1),historyShape(2)))
       read (fileHandle) self%time
       read (fileHandle) self%data
    end if
    return
  end subroutine History_Read_Raw

  subroutine History_Reset(history_)
    !!{
    Reset a history by zeroing all elements, but leaving the structure (and times) intact.
    !!}
    implicit none
    class(history), intent(inout) :: history_

    if (allocated(history_%time)) history_%data=0.0d0
    return
  end subroutine History_Reset

  subroutine History_Long_Integer_Dump(self,verbosityLevel)
    !!{
    Dumps a history object.
    !!}
    use :: Display           , only : displayMessage, enumerationVerbosityLevelType
    use :: ISO_Varying_String, only : assignment(=) , operator(//)                 , varying_string
    implicit none
    class    (longIntegerHistory           ), intent(in   ) :: self
    type     (enumerationVerbosityLevelType), intent(in   ) :: verbosityLevel
    integer                                                 :: i             , j
    type     (varying_string               )                :: message
    character(len=22                       )                :: label

    if (allocated(self%time)) then
       do i=1,size(self%time)
          write (label,'(i3)') i
          message="("//trim(label)//") "
          write (label,'(e22.16)') self%time(i)
          message=message//label//" :"
          do j=1,size(self%data,dim=2)
             write (label,'(i16)') self%data(i,j)
             message=message//" "//label
          end do
          call displayMessage(message,verbosityLevel)
       end do
    end if
    return
  end subroutine History_Long_Integer_Dump

  subroutine History_Long_Integer_Dump_Raw(self,fileHandle)
    !!{
    Dumps a history object in binary.
    !!}
    implicit none
    class  (longIntegerHistory), intent(in   ) :: self
    integer                    , intent(in   ) :: fileHandle

    write (fileHandle) self%rangeType
    write (fileHandle) allocated(self%time)
    if (allocated(self%time)) then
       write (fileHandle) shape(self%data)
       write (fileHandle) self%time
       write (fileHandle) self%data
    end if
    return
  end subroutine History_Long_Integer_Dump_Raw

  subroutine History_Long_Integer_Read_Raw(self,fileHandle)
    !!{
    Read a history object in binary.
    !!}
    implicit none
    class  (longIntegerHistory), intent(inout) :: self
    integer                    , intent(in   ) :: fileHandle
    logical                                    :: isAllocated
    integer                    , dimension(2)  :: historyShape

    read (fileHandle) self%rangeType
    read (fileHandle) isAllocated
    if (isAllocated) then
       read (fileHandle) historyShape
       allocate(self%time(historyShape(1)                ))
       allocate(self%data(historyShape(1),historyShape(2)))
       read (fileHandle) self%time
       read (fileHandle) self%data
    end if
    return
  end subroutine History_Long_Integer_Read_Raw

  subroutine History_Long_Integer_Reset(history_)
    !!{
    Reset a history by zeroing all elements, but leaving the structure (and times) intact.
    !!}
    implicit none
    class(longIntegerHistory), intent(inout) :: history_

    if (allocated(history_%time)) history_%data=0_kind_int8
    return
  end subroutine History_Long_Integer_Reset

  subroutine History_Set_To_Unity(history_)
    !!{
    Reset a history by zeroing all elements, but leaving the structure (and times) intact.
    !!}
    implicit none
    class(history), intent(inout) :: history_

    if (allocated(history_%time)) history_%data=1.0d0
    return
  end subroutine History_Set_To_Unity

  logical function History_Exists(history_)
    !!{
    Returns true if the history has been created.
    !!}
    implicit none
    class(history), intent(in   ) :: history_

    History_Exists=allocated(history_%time)
    return
  end function History_Exists

  subroutine History_Clone(self,historyToClone)
    !!{
    Clone a history object.
    !!}
    implicit none
    class(history), intent(inout) :: self
    type (history), intent(in   ) :: historyToClone

    if (allocated(self%time)) deallocate(self%time)
    if (allocated(self%data)) deallocate(self%data)
    if (allocated(historyToClone%time)) then
       allocate(self%time(size(historyToClone%time)))
       self%time=historyToClone%time
    end if
    if (allocated(historyToClone%data)) then
       allocate(self%data,mold=historyToClone%data)
       self%data=historyToClone%data
    end if
    self%rangeType=historyToClone%rangeType
    return
  end subroutine History_Clone

  logical function History_Long_Integer_Exists(history_)
    !!{
    Returns true if the history has been created.
    !!}
    implicit none
    class(longIntegerHistory), intent(in   ) :: history_

    History_Long_Integer_Exists=allocated(history_%time)
    return
  end function History_Long_Integer_Exists

  subroutine History_Long_Integer_Clone(self,historyToClone)
    !!{
    Clone a longIntegerHistory object.
    !!}
    implicit none
    class(longIntegerHistory), intent(inout) :: self
    type (longIntegerHistory), intent(in   ) :: historyToClone

    if (allocated(self%time)) deallocate(self%time)
    if (allocated(self%data)) deallocate(self%data)
    if (allocated(historyToClone%time)) then
       allocate(self%time(size(historyToClone%time)))
       self%time=historyToClone%time
    end if
    if (allocated(historyToClone%data)) then
       allocate(self%data,mold=historyToClone%data)
       self%data=historyToClone%data
    end if
    self%rangeType=historyToClone%rangeType
    return
  end subroutine History_Long_Integer_Clone

  logical function History_Is_Zero(self)
    !!{
    Test whether a history object is all zero.
    !!}
    implicit none
    class(history), intent(in   ) :: self

    History_Is_Zero=.true.
    if (allocated(self%data)) then
       if (any(self%data /= 0.0d0)) History_Is_Zero=.false.
    end if
    return
  end function History_Is_Zero

  function History_Add(history1,history2)
    !!{
    Add two history objects.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (history)                          :: History_Add
    class(history), intent(in   )           :: history1
    class(history), intent(in   ), optional :: history2

    ! Clone the first history.
    if (.not.(history1%exists().or.history2%exists())) then
       call History_Add%reset()
    else
       select type (history1)
       type is (history)
          History_Add=history1
          if (present(history2)) then
             if (any(shape(History_Add%data) /= shape(history2%data))) call Error_Report('mismatch in history object shape'//{introspection:location})
             History_Add%data=History_Add%data+history2%data
          end if
       end select
    end if
    return
  end function History_Add

  function History_Subtract(history1,history2)
    !!{
    Subtract two history objects.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (history)                          :: History_Subtract
    class(history), intent(in   )           :: history1
    class(history), intent(in   ), optional :: history2

    ! Clone the first history.
    if (.not.(history1%exists().or.(present(history2).and.history2%exists()))) then
       call History_Subtract%reset()
    else
       select type (history1)
       type is (history)
          History_Subtract=history1
          if (present(history2)) then
             if (any(shape(history1%data) /= shape(history2%data))) call Error_Report('mismatch in history object shape'//{introspection:location})
             History_Subtract%data= history1%data-history2%data
          else
             History_Subtract%data=-history1%data
          end if
       end select
    end if
    return
  end function History_Subtract

  integer function History_Serialize_Count(self)
    !!{
    Return the number of properties required to track a history.
    !!}
    implicit none
    class(history), intent(in   ) :: self

    if (allocated(self%data)) then
       History_Serialize_Count=size(self%data)
    else
       History_Serialize_Count=0
    end if
    return
  end function History_Serialize_Count

  subroutine History_Deserialize(self,historyArray)
    !!{
    Pack history from an array into a history structure.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (history)              , intent(inout) :: self
    double precision         , dimension(:), intent(in   ) :: historyArray

    ! Extract data from array.
    if (allocated(self%data)) then
       self%data=reshape(historyArray,shape(self%data))
    else if (size(historyArray) > 0) then
       call Error_Report('attempt to deserialize into non-existant history'//{introspection:location})
    end if
    return
  end subroutine History_Deserialize

  subroutine History_Serialize(self,historyArray)
    !!{
    Pack history from an array into an history structure.
    !!}
    implicit none
    class           (history)              , intent(in   ) :: self
    double precision         , dimension(:), intent(  out) :: historyArray(:)

    ! Place data into array.
    if (allocated(self%data)) historyArray(:)=reshape(self%data,shape(historyArray))
    return
  end subroutine History_Serialize

  subroutine History_Trim(history_,currentTime,minimumPointsToRemove)
    !!{
    Removes outdated information from ``future histories'' (i.e. histories that store data for future reference). Removes all
    but one entry prior to the given {\normalfont \ttfamily currentTime} (this allows for interpolation of the history to the current
    time). Optionally, the remove is done only if it will remove more than {\normalfont \ttfamily minimumPointsToRemove} entries (since the
    removal can be slow this allows for some optimization).
    !!}
    use            :: Error            , only : Error_Report
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    implicit none
    class           (history), intent(inout)           :: history_
    double precision         , intent(in   )           :: currentTime
    integer                  , intent(in   ), optional :: minimumPointsToRemove
    type            (history)                          :: temporaryHistory
    integer                                            :: currentPointCount    , historyCount               , &
         &                                                iTrim                , minimumPointsToRemoveActual, &
         &                                                newPointCount

    ! Return if no history exists.
    if (.not.allocated(history_%time)) return

    ! Find points to remove.
    currentPointCount=size(history_%time)

    ! Return is nothing to trim.
    if (currentPointCount == 0) return

    ! Decide on the minimum number of points that we will remove.
    if (present(minimumPointsToRemove)) then
       if (minimumPointsToRemove < 1) call Error_Report('minimum number of points to remove must be >= 1'//{introspection:location})
       minimumPointsToRemoveActual=minimumPointsToRemove
    else
       minimumPointsToRemoveActual=1
    end if

    ! Find how much we can trim. Never trim the final two point as they might be needed to extrapolate beyond the end of the
    ! future history. Having found a point which exceeds the current time, pull back two points, so that we leave one point prior
    ! to the current time.
    iTrim=1
    do while (history_%time(iTrim) < currentTime .and. iTrim <= currentPointCount-2)
       iTrim=iTrim+1
    end do
    iTrim=iTrim-2

    ! Check if there are enough removable points to warrant actually doing the removal.
    if (iTrim >= minimumPointsToRemoveActual) then
       ! Move current history to temporary storage.
       call Move_Alloc(history_%time  ,temporaryHistory%time  )
       call Move_Alloc(history_%data  ,temporaryHistory%data  )
       ! Reallocate the history arrays.
       newPointCount=currentPointCount-iTrim
       historyCount =size(temporaryHistory%data,dim=2)
       allocate(history_%time(newPointCount             ))
       allocate(history_%data(newPointCount,historyCount))
       ! Copy the data back into the new arrays.
       history_%time(:  )=temporaryHistory%time(iTrim+1:currentPointCount  )
       history_%data(:,:)=temporaryHistory%data(iTrim+1:currentPointCount,:)
       ! Deallocate the temporary arrays.
       deallocate(temporaryHistory%time  )
       deallocate(temporaryHistory%data  )
    end if
    return
  end subroutine History_Trim

  subroutine History_Trim_Forward(self,time,removedHistory)
    !!{
    Removes all points in a history after the given {\normalfont \ttfamily time}. Optionally, the removed history can be
    returned as {\normalfont \ttfamily removedHistory}.
    !!}
    use            :: Arrays_Search    , only : searchArray
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    implicit none
    class           (history ), intent(inout)           :: self
    double precision          , intent(in   )           :: time
    type            (history ), intent(inout), optional :: removedHistory
    type            (history )                          :: temporaryHistory
    integer         (c_size_t)                          :: trimAt          , trimCount

    ! Ensure the removed history to be returned is initialized.
    if (present(removedHistory).and.allocated(removedHistory%time)) then
       deallocate(removedHistory%time)
       deallocate(removedHistory%data)
    end if
    ! Return if no history exists or if the final time is prior to the trim time.
    if (.not.allocated(self%time).or.self%time(size(self%time)) <= time) return
    ! Find where to trim and number of trimmed points.
    trimAt   =searchArray(self%time,time)+1
    trimCount=size(self%time)-trimAt+1
    ! Transfer data to a temporary history.
    call Move_Alloc(self%time,temporaryHistory%time)
    call Move_Alloc(self%data,temporaryHistory%data)
    ! Reallocate history to trimmed size and populate.
    if (trimAt > 1) then
       allocate(self%time(trimAt-1                                  ))
       allocate(self%data(trimAt-1,size(temporaryHistory%data,dim=2)))
       self%time=temporaryHistory%time(1:trimAt-1  )
       self%data=temporaryHistory%data(1:trimAt-1,:)
    end if
    ! If the trimmed history is to be returned, allocate the arrays and populate.
    if (present(removedHistory)) then
       allocate(removedHistory%time(trimCount                                  ))
       allocate(removedHistory%data(trimCount,size(temporaryHistory%data,dim=2)))
       removedHistory%time=temporaryHistory%time(trimAt:trimAt+trimCount-1  )
       removedHistory%data=temporaryHistory%data(trimAt:trimAt+trimCount-1,:)
    end if
    ! Clean up temporary history.
    deallocate(temporaryHistory%time)
    deallocate(temporaryHistory%data)
    return
  end subroutine History_Trim_Forward

  subroutine History_Long_Integer_Trim(history_,currentTime,minimumPointsToRemove)
    !!{
    Removes outdated information from ``future histories'' (i.e. histories that store data for future reference). Removes all
    but one entry prior to the given {\normalfont \ttfamily currentTime} (this allows for interpolation of the history to the current
    time). Optionally, the remove is done only if it will remove more than {\normalfont \ttfamily minimumPointsToRemove} entries (since the
    removal can be slow this allows for some optimization).
    !!}
    use            :: Error            , only : Error_Report
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    implicit none
    class           (longIntegerHistory), intent(inout)           :: history_
    double precision                    , intent(in   )           :: currentTime
    integer                             , intent(in   ), optional :: minimumPointsToRemove
    type            (longIntegerHistory)                          :: temporaryHistory
    integer                                                       :: currentPointCount    , historyCount               , &
         &                                                           iTrim                , minimumPointsToRemoveActual, &
         &                                                           newPointCount

    ! Return if no history exists.
    if (.not.allocated(history_%time)) return

    ! Find points to remove.
    currentPointCount=size(history_%time)

    ! Return is nothing to trim.
    if (currentPointCount == 0) return

    ! Decide on the minimum number of points that we will remove.
    if (present(minimumPointsToRemove)) then
       if (minimumPointsToRemove < 1) call Error_Report('minimum number of points to remove must be >= 1'//{introspection:location})
       minimumPointsToRemoveActual=minimumPointsToRemove
    else
       minimumPointsToRemoveActual=1
    end if

    ! Find how much we can trim. Never trim the final two point as they might be needed to extrapolate beyond the end of the
    ! future history. Having found a point which exceeds the current time, pull back two points, so that we leave one point prior
    ! to the current time.
    iTrim=1
    do while (history_%time(iTrim) < currentTime .and. iTrim <= currentPointCount-2)
       iTrim=iTrim+1
    end do
    iTrim=iTrim-2

    ! Check if there are enough removable points to warrant actually doing the removal.
    if (iTrim >= minimumPointsToRemoveActual) then
       ! Move current history to temporary storage.
       call Move_Alloc(history_%time  ,temporaryHistory%time  )
       call Move_Alloc(history_%data  ,temporaryHistory%data  )
       ! Reallocate the history arrays.
       newPointCount=currentPointCount-iTrim
       historyCount =size(temporaryHistory%data,dim=2)
       allocate(history_%time(newPointCount             ))
       allocate(history_%data(newPointCount,historyCount))
       ! Copy the data back into the new arrays.
       history_%time(:  )=temporaryHistory%time(iTrim+1:currentPointCount  )
       history_%data(:,:)=temporaryHistory%data(iTrim+1:currentPointCount,:)
       ! Deallocate the temporary arrays.
       deallocate(temporaryHistory%time  )
       deallocate(temporaryHistory%data  )
    end if
    return
  end subroutine History_Long_Integer_Trim

  subroutine History_Long_Integer_Trim_Forward(self,time,removedHistory)
    !!{
    Removes all points in a history after the given {\normalfont \ttfamily time}. Optionally, the removed history can be
    returned as {\normalfont \ttfamily removedHistory}.
    !!}
    use            :: Arrays_Search    , only : searchArray
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    implicit none
    class           (longIntegerHistory), intent(inout)           :: self
    double precision                    , intent(in   )           :: time
    type            (longIntegerHistory), intent(inout), optional :: removedHistory
    type            (longIntegerHistory)                          :: temporaryHistory
    integer         (c_size_t          )                          :: trimAt          , trimCount

    ! Ensure the removed history to be returned is initialized.
    if (present(removedHistory).and.allocated(removedHistory%time)) then
       deallocate(removedHistory%time)
       deallocate(removedHistory%data)
    end if
    ! Return if no history exists or if the final time is prior to the trim time.
    if (.not.allocated(self%time).or.self%time(size(self%time)) <= time) return
    ! Find where to trim and number of trimmed points.
    trimAt   =searchArray(self%time,time)+1
    trimCount=size(self%time)-trimAt+1
    ! Transfer data to a temporary history.
    call Move_Alloc(self%time,temporaryHistory%time)
    call Move_Alloc(self%data,temporaryHistory%data)
    ! Reallocate history to trimmed size and populate.
    if (trimAt > 1) then
       allocate(self%time(trimAt-1                                  ))
       allocate(self%data(trimAt-1,size(temporaryHistory%data,dim=2)))
       self%time=temporaryHistory%time(1:trimAt-1  )
       self%data=temporaryHistory%data(1:trimAt-1,:)
    end if
    ! If the trimmed history is to be returned, allocate the arrays and populate.
    if (present(removedHistory)) then
       allocate(removedHistory%time(trimCount                                  ))
       allocate(removedHistory%data(trimCount,size(temporaryHistory%data,dim=2)))
       removedHistory%time=temporaryHistory%time(trimAt:trimAt+trimCount-1  )
       removedHistory%data=temporaryHistory%data(trimAt:trimAt+trimCount-1,:)
    end if
    ! Clean up temporary history.
    deallocate(temporaryHistory%time)
    deallocate(temporaryHistory%data)
    return
  end subroutine History_Long_Integer_Trim_Forward

  subroutine History_Long_Integer_Append_History(self,append)
    !!{
    Append a history to a long integer history.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class           (longIntegerHistory), intent(inout)                 :: self
    type            (longIntegerHistory), intent(in   )                 :: append
    double precision                    , allocatable  , dimension(:  ) :: timeTmp
    integer         (kind=kind_int8    ), allocatable  , dimension(:,:) :: dataTmp

    if (.not.allocated(self%time)) then
       ! No pre-existing history - simply copy the history to append.
       self%time=append%time
       self%data=append%data
    else
       ! A history already exists. Validate the provided append history.
       if (append%time(1) <= self%time(size(self%time))) call Error_Report('history to append starts before end of history to which it is being appended'//{introspection:location})
       if (size(self%data,dim=2) /= size(append%data,dim=2)) call Error_Report('histories have different cardinalities'//{introspection:location})
       ! Do the append.
       allocate(timeTmp(size(self%time)+size(append%time)                      ))
       allocate(dataTmp(size(self%time)+size(append%time),size(self%data,dim=2)))
       timeTmp(                1:size(self%time)                    )=self  %time
       timeTmp(size(self%time)+1:size(self%time)+size(append%time)  )=append%time
       dataTmp(                1:size(self%time)                  ,:)=self  %data
       dataTmp(size(self%time)+1:size(self%time)+size(append%time),:)=append%data
       deallocate     (        self%time)
       deallocate     (        self%data)
       call move_alloc(timeTmp,self%time)
       call move_alloc(dataTmp,self%data)
    end if
    return
  end subroutine History_Long_Integer_Append_History

  subroutine History_Long_Integer_Append_Epoch(self,time,append)
    !!{
    Append a history to a long integer history.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class           (longIntegerHistory), intent(inout)                 :: self
    double precision                    , intent(in   )                 :: time
    integer         (kind=kind_int8    ), intent(in   ), dimension(:  ) :: append
    double precision                    , allocatable  , dimension(:  ) :: timeTmp
    integer         (kind=kind_int8    ), allocatable  , dimension(:,:) :: dataTmp

    if (.not.allocated(self%time)) then
       ! No pre-existing history - simply copy the history to append.
       allocate(self%time(1              ))
       allocate(self%data(1,size(append)))
       self%time(1  )=time
       self%data(1,:)=append
    else
       ! A history already exists. Validate the provided append history.
       if (time <= self%time(size(self%time))) call Error_Report('history to append starts before end of history to which it is being appended'//{introspection:location})
       if (size(self%data,dim=2) /= size(append)) call Error_Report('histories have different cardinalities'//{introspection:location})
       ! Do the append.
       allocate(timeTmp(size(self%time)+1                      ))
       allocate(dataTmp(size(self%time)+1,size(self%data,dim=2)))
       timeTmp(                1:size(self%time)    )=self   %time
       timeTmp(size(self%time)+1                    )=        time
       dataTmp(                1:size(self%time)  ,:)=self   %data
       dataTmp(size(self%time)+1                  ,:)=        append
       deallocate     (        self%time)
       deallocate     (        self%data)
       call move_alloc(timeTmp,self%time)
       call move_alloc(dataTmp,self%data)
    end if
    return
  end subroutine History_Long_Integer_Append_Epoch

  subroutine History_Append_History(self,append)
    !!{
    Append a history to a long integer history.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class           (history), intent(inout)                 :: self
    type            (history), intent(in   )                 :: append
    double precision         , allocatable  , dimension(:  ) :: timeTmp
    double precision         , allocatable  , dimension(:,:) :: dataTmp

    if (.not.allocated(self%time)) then
       ! No pre-existing history - simply copy the history to append.
       self%time=append%time
       self%data=append%data
    else
       ! A history already exists. Validate the provided append history.
       if (append%time(1) <= self%time(size(self%time))) call Error_Report('history to append starts before end of history to which it is being appended'//{introspection:location})
       if (size(self%data,dim=2) /= size(append%data,dim=2)) call Error_Report('histories have different cardinalities'//{introspection:location})
       ! Do do the append.
       allocate(timeTmp(size(self%time)+size(append%time)                      ))
       allocate(dataTmp(size(self%time)+size(append%time),size(self%data,dim=2)))
       timeTmp(                1:size(self%time)                    )=self  %time
       timeTmp(size(self%time)+1:size(self%time)+size(append%time)  )=append%time
       dataTmp(                1:size(self%time)                  ,:)=self  %data
       dataTmp(size(self%time)+1:size(self%time)+size(append%time),:)=append%data
       deallocate     (        self%time)
       deallocate     (        self%data)
       call move_alloc(timeTmp,self%time)
       call move_alloc(dataTmp,self%data)
    end if
    return
  end subroutine History_Append_History

  subroutine History_Append_Epoch(self,time,append)
    !!{
    Append a history to a long integer history.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class           (history), intent(inout)                 :: self
    double precision         , intent(in   )                 :: time
    double precision         , intent(in   ), dimension(:  ) :: append
    double precision         , allocatable  , dimension(:  ) :: timeTmp
    double precision         , allocatable  , dimension(:,:) :: dataTmp

    if (.not.allocated(self%time)) then
       ! No pre-existing history - simply copy the history to append.
       allocate(self%time(1              ))
       allocate(self%data(1,size(append)))
       self%time(1  )=time
       self%data(1,:)=append
    else
       ! A history already exists. Validate the provided append data.
       if (time <= self%time(size(self%time))) call Error_Report('history to append starts before end of history to which it is being appended'//{introspection:location})
       if (size(self%data,dim=2) /= size(append)) call Error_Report('histories have different cardinalities'//{introspection:location})
       ! Do the append.
       allocate(timeTmp(size(self%time)+1                      ))
       allocate(dataTmp(size(self%time)+1,size(self%data,dim=2)))
       timeTmp(                1:size(self%time)    )=self   %time
       timeTmp(size(self%time)+1                    )=        time
       dataTmp(                1:size(self%time)  ,:)=self   %data
       dataTmp(size(self%time)+1                  ,:)=        append
       deallocate     (        self%time)
       deallocate     (        self%data)
       call move_alloc(timeTmp,self%time)
       call move_alloc(dataTmp,self%data)
    end if
    return
  end subroutine History_Append_Epoch

   subroutine History_Interpolated_Increment(history_,addHistory)
     !!{
     Adds the data in {\normalfont \ttfamily addHistory} to that in {\normalfont \ttfamily history\_}. This function is
     designed for histories that track instantaneous rates. The rates in {\normalfont \ttfamily addHistory} are interpolated to
     the times in {\normalfont \ttfamily history\_} and added to the rates in {\normalfont \ttfamily history\_}.
     !!}
     use            :: Error                  , only : Error_Report
     use, intrinsic :: ISO_C_Binding          , only : c_size_t
     use            :: Numerical_Interpolation, only : interpolator
     implicit none
     class           (history     ), intent(inout) :: history_
     type            (history     ), intent(in   ) :: addHistory
     double precision              , dimension(2)  :: interpolationFactors
     integer                                       :: iPoint              , iHistory
     integer         (c_size_t    )                :: interpolationPoint  , addHistoryPointCount
     type            (interpolator)                :: interpolator_

     select type (history_)
     type is (history)
        ! Return if addHistory does not exist.
        if (.not.allocated(addHistory%time)) return
        ! Get size of addHistory.
        addHistoryPointCount=size(addHistory%time)
        ! Return if addHistory has zero size.
        if (addHistoryPointCount == 0) return
        ! If history_ does not exist, just replace it with addHistory.
        if (.not.allocated(history_%time)) then
           call history_%destroy()
           history_=addHistory
           return
        end if
        ! If history_ has zero size, just replace it with addHistory.
        if (size(history_%time) == 0) then
           call history_%destroy()
           history_=addHistory
           return
        end if
        ! addHistory must have at least two points to permit interpolation.
        if (addHistoryPointCount  < 2) call Error_Report('history to add must have at least two points'//{introspection:location})
        ! The two objects must contain the same number of histories.
        if (size(history_%data,dim=2) /= size(addHistory%data,dim=2)) call Error_Report('two objects contain differing numbers of histories'//{introspection:location})
        ! Loop over each entry in history_.
        interpolationPoint=1
        interpolator_     =interpolator(addHistory%time)
        do iPoint=1,size(history_%time)
           ! If within range of history spanned by addHistory then....
           if (history_%time(iPoint) >= addHistory%time(1) .and. history_%time(iPoint) <= addHistory%time(addHistoryPointCount)) then
              ! Interpolate addHistory to point in history_.
              do while (history_%time(iPoint) > addHistory%time(interpolationPoint) .and. interpolationPoint < addHistoryPointCount-1)
                 interpolationPoint=interpolationPoint+1
              end do
              call interpolator_%linearWeights(history_%time(iPoint),interpolationPoint,interpolationFactors)
              ! Add them.
              forall(iHistory=1:size(history_%data,dim=2))
                 history_%data (iPoint,iHistory)=history_%data (iPoint,iHistory)+addHistory%data(interpolationPoint,iHistory)&
                      &*interpolationFactors(1)+addHistory%data(interpolationPoint+1,iHistory)*interpolationFactors(2)
              end forall
           end if

        end do

     end select

     return
   end subroutine History_Interpolated_Increment

   subroutine History_Increment(history_,addHistory,autoExtend)
     !!{
     Combines the data in {\normalfont \ttfamily addHistory} with that in {\normalfont \ttfamily history\_}. This function is designed for histories that
     track integrated quantities (such as total mass of stars formed in a time interval for example). {\normalfont \ttfamily history\_} will be
     extended if necessary to span the range of {\normalfont \ttfamily addHistory}. Then, the data from {\normalfont \ttfamily addHistory} will be added to
     that in {\normalfont \ttfamily history\_} by finding the fraction of each timestep in {\normalfont \ttfamily addHistory} that overlaps with each timestep
     in {\normalfont \ttfamily history\_} and assuming that the corresponding fraction of the data value should be added to {\normalfont \ttfamily history\_}.
     !!}
     use            :: Arrays_Search   , only : searchArray
     use            :: Error           , only : Error_Report
     use, intrinsic :: ISO_C_Binding   , only : c_size_t
     use            :: Numerical_Ranges, only : rangeTypeUndefined
     implicit none
     class           (history ), intent(inout)           :: history_
     type            (history ), intent(in   )           :: addHistory
     logical                   , intent(in   ), optional :: autoExtend
     integer                                             :: addCount           , addHistoryPointCount
     integer         (c_size_t)                          :: timeBeginIndex     , timeEndIndex        , &
          &                                                 iPoint             , jPoint
     double precision                                    :: fractionContributed, timeBegin           , &
          &                                                 timeEnd            , timeBeginAdd        , &
          &                                                 timeEndAdd
     !![
     <optionalArgument name="autoExtend" defaultsTo=".false." />
     !!]

     select type (history_)
     type is (history)

        ! Return if addHistory does not exist.
        if (.not.allocated(addHistory%time)) return

        ! Get size of addHistory.
        addHistoryPointCount=size(addHistory%time)

        ! Return if addHistory has zero size.
        if (addHistoryPointCount == 0) return

        ! If history_ does not exist, simply replace it with addHistory.
        if (.not.allocated(history_%time)) then
           call history_%destroy()
           history_=addHistory
           return
        end if

        ! If history_ has zero size, simply replace it with addHistory.
        if (size(history_%time) == 0) then
           call history_%destroy()
           history_=addHistory
           return
        end if

        ! The two objects must contain the same number of histories.
        if (size(history_%data,dim=2) /= size(addHistory%data,dim=2)) call Error_Report('two objects contain differing numbers of histories'//{introspection:location})

        ! Determine if we need to extend the time range in history_.
        addCount=size(addHistory%time)
        if (addHistory%time(1) < history_%time(1) .or. addHistory%time(addCount) > history_%time(size(history_%time))) then
           if (.not.autoExtend_) call Error_Report("history needs to be extended, but is not permitted"//{introspection:location})
           if (history_%rangeType == rangeTypeUndefined) then
              ! The history has no defined range type, so pass the time array of the history being added to use as a template for new times.
              call history_%extend(times=addHistory%time)
           else
              ! The history has a defined range type, so simply pass the required extent of the range.
              call history_%extend([addHistory%time(1),addHistory%time(addCount)])
           end if
        end if
        ! Transfer each entry from addHistory to history_.
        do iPoint=1,addCount
           if (iPoint == 1) then
              timeBeginAdd=0.0d0
           else
              timeBeginAdd=addHistory%time(iPoint-1)
           end if
           timeEndAdd     =addHistory%time(iPoint  )
           ! Find indices in history_ spanned by addHistory point.
           if (iPoint > 1) then
              ! Reuse the end index from the previous loop iteration if available.
              timeBeginIndex=timeEndIndex
           else
              timeBeginIndex=1
           end if
           timeEndIndex=min(searchArray(history_%time,addHistory%time(iPoint))+1,size(history_%time))
           ! Loop over all points in history_ to which we need to add this contribution.
           do jPoint=timeBeginIndex,timeEndIndex
              if (jPoint == 1) then
                 timeBegin=                            timeBeginAdd
              else
                 timeBegin=max(history_%time(jPoint-1),timeBeginAdd)
              end if
              timeEnd     =min(history_%time(jPoint  ),timeEndAdd  )
              fractionContributed=max(0.0d0,(timeEnd-timeBegin)/(timeEndAdd-timeBeginAdd))
              history_%data(jPoint,:)=history_%data(jPoint,:)+addHistory%data(iPoint,:)*fractionContributed
           end do
        end do
     end select
     return
   end subroutine History_Increment

   function History_Divide(self,divisor)
     !!{
     Divides history data by a double precision {\normalfont \ttfamily divisor}.
     !!}
     implicit none
     type            (history)                :: History_Divide
     class           (history), intent(in   ) :: self
     double precision         , intent(in   ) :: divisor

     select type(self)
     type is (history)
        History_Divide=self
        if (allocated(History_Divide%data)) History_Divide%data=History_Divide%data/divisor
     end select
     return
   end function History_Divide

   function History_Multiply(self,multiplier)
     !!{
     Multiplies history data by a double precision {\normalfont \ttfamily multiplier}.
     !!}
     implicit none
     type            (history)                :: History_Multiply
     class           (history), intent(in   ) :: self
     double precision         , intent(in   ) :: multiplier

     select type(self)
     type is (history)
        History_Multiply=self
        if (allocated(History_Multiply%data)) History_Multiply%data=History_Multiply%data*multiplier
     end select
     return
   end function History_Multiply

  function History_Multiply_Switched(multiplier,history1)
    !!{
    Multiply a scalar by an history object.
    !!}
    implicit none
    type            (history)                :: History_Multiply_Switched
    type            (history), intent(in   ) :: history1
    double precision         , intent(in   ) :: multiplier

    History_Multiply_Switched=History_Multiply(history1,multiplier)
    return
  end function History_Multiply_Switched

   subroutine History_Extend(history_,timeRange,times)
     !!{
     Extends a history to encompass the given time range.
     !!}
     use :: Error             , only : Error_Report
     use :: ISO_Varying_String, only : assignment(=)  , operator(//)        , varying_string
     use :: Numerical_Ranges  , only : rangeTypeLinear, rangeTypeLogarithmic, rangeTypeUndefined
     use :: String_Handling   , only : operator(//)
     implicit none
     class           (history       )                             , intent(inout)           :: history_
     double precision                             , dimension(2  ), intent(in   ), optional :: timeRange
     double precision                             , dimension(:  ), intent(in   ), optional :: times
     double precision                , allocatable, dimension(:  )                          :: newTimes
     double precision                , allocatable, dimension(:,:)                          :: historyDataTemporary
     double precision                             , dimension(2  )                          :: timeRangeActual
     integer                                                                                :: addCount            , addCountEnd   , addCountStart  , &
          &                                                                                    historyCount        , newTimesAtEnd , newTimesAtStart, &
          &                                                                                    rangeType           , timeBeginIndex, timeCount      , &
          &                                                                                    timeEndIndex
     logical                                                                                :: useRange
     double precision                                                                       :: timeBegin           , timeDelta     , timeEnd
     type            (varying_string)                                                       :: message

     ! Determine the range of times that must be covered.
     if (present(timeRange)) then
        timeRangeActual=timeRange
     else
        if (present(times)) then
           timeRangeActual(1)=times(1          )
           timeRangeActual(2)=times(size(times))
        else
           call Error_Report('either timeRange or times must be specified'//{introspection:location})
        end if
     end if

     ! Determine if we need to extend the time range in history_.
     timeCount     =size(history_%time           )
     timeBegin     =     history_%time(1        )
     timeEnd       =     history_%time(timeCount)
     timeBeginIndex=1
     timeEndIndex  =timeCount
     if (.not.present(times)) then
        select case (history_%rangeType)
        case (rangeTypeLinear     )
           timeDelta=(timeEnd-timeBegin)/dble(timeCount-1)
        case (rangeTypeLogarithmic)
           timeDelta=log(timeEnd/timeBegin)/dble(timeCount-1)
        case default
           timeDelta=0.0d0
           if (history_%rangeType == rangeTypeUndefined) then
              message='undefined range type: '//char(10)
           else
              message='unrecognized range type: '
              message=message//history_%rangeType//char(10)
           end if
           message=message//' -> known types are: '//char(10)//' --> linear      : '//rangeTypeLinear//char(10)//' --> logarithmic : '//rangeTypeLogarithmic
           call Error_Report(message//{introspection:location})
        end select
        if (timeRangeActual(1) < timeBegin) then
           select case (history_%rangeType)
           case (rangeTypeLinear)
              addCountStart =int(    (timeBegin-timeRangeActual(1))/timeDelta)+1
              timeBegin=timeBegin      -dble(addCountStart)*timeDelta
           case (rangeTypeLogarithmic)
              addCountStart =int(log(timeBegin/timeRangeActual(1))/timeDelta)+1
              timeBegin=timeBegin*exp(-dble(addCountStart)*timeDelta)
           case default
              addCountStart=0
              call Error_Report('unknown range type'//{introspection:location})
           end select
           timeBeginIndex=timeBeginIndex+addCountStart
           timeEndIndex  =timeEndIndex  +addCountStart
        else
           addCountStart=0
        end if
        if (timeRangeActual(2) > timeEnd  ) then
           select case (history_%rangeType)
           case (rangeTypeLinear)
              addCountEnd =int(    (timeRangeActual(2)-timeEnd)/timeDelta)+1
              timeEnd  =timeEnd        +dble(addCountEnd)*timeDelta
           case (rangeTypeLogarithmic)
              addCountEnd =int(log(timeRangeActual(2)/timeEnd)/timeDelta)+1
              timeEnd  =timeEnd  *exp(+dble(addCountEnd)*timeDelta)
           case default
              addCountEnd=0
              call Error_Report('unknown range type'//{introspection:location})
           end select
        else
           addCountEnd=0
        end if
        addCount=addCountStart+addCountEnd
        allocate(newTimes(0))
        useRange=.true.
     else
        newTimesAtStart=count(times < history_%time(1                     ))
        newTimesAtEnd  =count(times > history_%time(size(history_%time)))
        timeBeginIndex=timeBeginIndex+newTimesAtStart
        timeEndIndex  =timeEndIndex  +newTimesAtStart
        addCount=newTimesAtStart+newTimesAtEnd
        allocate(newTimes(size(history_%time)+addCount))
        if (newTimesAtStart > 0) newTimes(1             :             newTimesAtStart)=times(1                          :newTimesAtStart)
        if (newTimesAtEnd   > 0) newTimes(timeEndIndex+1:timeEndIndex+newTimesAtEnd  )=times(size(times)-newTimesAtEnd+1:size(times)    )
        newTimes(timeBeginIndex:timeEndIndex)=history_%time
        useRange=.false.
     end if

     ! Create new arrays.
     if (addCount > 0) then
        ! Create copies of current histories.
        call Move_Alloc(history_%data,historyDataTemporary)
        ! Store range type and number of histories.
        rangeType   =history_%rangeType
        historyCount=size(historyDataTemporary,dim=2)
        ! Destroy the history and make a new one.
        call history_%destroy()
        select case (useRange)
        case (.true. )
           call history_%create(historyCount,timeCount+addCount,timeBegin,timeEnd,rangeType)
        case (.false.)
           call history_%create(historyCount,timeCount+addCount)
           history_%time=newTimes
           deallocate(newTimes)
        end select
        ! Copy data back to relevant location.
        history_%data(timeBeginIndex:timeEndIndex,:)=historyDataTemporary
        deallocate(historyDataTemporary)
     end if
     return
   end subroutine History_Extend

  subroutine History_Timesteps(history_,timeSteps)
    !!{
    Return an array of time intervals in {\normalfont \ttfamily history\_}.
    !!}
    use :: Numerical_Ranges , only : rangeTypeLinear, rangeTypeLogarithmic
    implicit none
    class           (history)                           , intent(in   ) :: history_
    double precision         , allocatable, dimension(:), intent(inout) :: timeSteps
    integer                                                             :: iTime
    double precision                                                    :: ratio

    allocate(timeSteps(size(history_%time)))
    select case (history_%rangeType)
    case (rangeTypeLogarithmic)
       ratio=history_%time(2)/history_%time(1)
       forall(iTime=1:size(history_%time))
          timeSteps(iTime)=history_%time(1)*(ratio**(dble(iTime)-0.5d0)-ratio**(dble(iTime)-1.5d0))
       end forall
    case (rangeTypeLinear     )
       timeSteps=history_%time(2)-history_%time(1)
    case default
       timeSteps(1)=history_%time(1)
       forall(iTime=2:size(history_%time))
          timeSteps(iTime)=history_%time(iTime)-history_%time(iTime-1)
       end forall
    end select
    return
  end subroutine History_Timesteps

  function History_Non_Static_Size_Of(self)
    !!{
    Return the size of any non-static components of the object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t)                :: History_Non_Static_Size_Of
    class  (history ), intent(in   ) :: self

    if (allocated(self%time)) then
       History_Non_Static_Size_Of=sizeof(self%time)+sizeof(self%data)
    else
       History_Non_Static_Size_Of=0_c_size_t
    end if
    return
  end function History_Non_Static_Size_Of

  function History_Long_Integer_Non_Static_Size_Of(self)
    !!{
    Return the size of any non-static components of the object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t          )                :: History_Long_Integer_Non_Static_Size_Of
    class  (longIntegerHistory), intent(in   ) :: self

    if (allocated(self%time)) then
       History_Long_Integer_Non_Static_Size_Of=sizeof(self%time)+sizeof(self%data)
    else
       History_Long_Integer_Non_Static_Size_Of=0_c_size_t
    end if
    return
  end function History_Long_Integer_Non_Static_Size_Of

end module Histories
