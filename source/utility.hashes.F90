!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements ``hashes'' (i.e. associative arrays).

module Hashes
  !% Implements ``hashes'' (i.e. associative arrays).
  use ISO_Varying_String
  implicit none

  !# <generic identifier="Type">
  !#  <instance label="integer" intrinsic="integer"          />
  !#  <instance label="double"  intrinsic="double precision" />
  !# </generic>

  private
  public :: {Type¦label}ScalarHash

  type :: {Type¦label}ScalarHash
     !% Derived type for {Type¦label} hashes.
     private
     integer                                       :: allocatedSize   , elementCount
     {Type¦intrinsic}                , allocatable :: hashValues   (:)
     type            (varying_string), allocatable :: hashKeys     (:)
   contains
     !@ <objectMethods>
     !@   <object>{Type¦label}ScalarHash</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <description>Initialize the hash.</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>set</method>
     !@     <description>Set the value of a key in the hash.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless (character(len=*)|varying\_string)\textgreater} key\argin, \intzero\ value\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>delete</method>
     !@     <description>Delete a key from the hash.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless (character(len=*)|varying\_string)\textgreater} key\argin, \intzero\ value\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>value</method>
     !@     <description>Return the value for the given key.</description>
     !@     <type>\intzero</type>
     !@     <arguments>\textcolor{red}{\textless (character(len=*)|varying\_string|\intzero)\textgreater} key\argin, \intzero\ value\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>key</method>
     !@     <description>Return the key of the {\normalfont \ttfamily indexValue}$^{\mathrm th}$ entry in the hash.</description>
     !@     <type>\textcolor{red}{\textless type(varying\_string)\textgreater}</type>
     !@     <arguments>\intzero\ indexValue\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>keys</method>
     !@     <description>Return an array of all keys in the hash.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(varying\_string)[:]\textgreater} keys\arginout</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>values</method>
     !@     <description>Return an array of all values in the hash.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless {Type¦intrinsic}[:]\textgreater} values\arginout</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>exists</method>
     !@     <description>Return true if the specified key exists in the hash.</description>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless (character(len=*)|varying\_string|\intzero)\textgreater} key\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>size</method>
     !@     <description>Return the number of keys in the hash.</description>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <description>Destroy the hash.</description>
     !@     <type>void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                        {Type¦label}ScalarDestructor
     procedure :: initialize           =>Initialize_{Type¦label}_Scalar
     procedure :: Set_{Type¦label}_Scalar_VS
     procedure :: Set_{Type¦label}_Scalar_CH
     generic   :: set        => Set_{Type¦label}_Scalar_VS,Set_{Type¦label}_Scalar_CH
     procedure :: Delete_{Type¦label}_Scalar_VS
     procedure :: Delete_{Type¦label}_Scalar_CH
     generic   :: delete     => Delete_{Type¦label}_Scalar_VS,Delete_{Type¦label}_Scalar_CH
     procedure :: Value_{Type¦label}_Scalar_VS
     procedure :: Value_{Type¦label}_Scalar_CH
     procedure :: Value_{Type¦label}_Scalar_I
     generic   :: value      => Value_{Type¦label}_Scalar_VS,Value_{Type¦label}_Scalar_CH,Value_{Type¦label}_Scalar_I
     procedure :: key                     =>Key_{Type¦label}_Scalar_I
     procedure :: Exists_{Type¦label}_Scalar_VS
     procedure :: Exists_{Type¦label}_Scalar_CH
     procedure :: keys                    =>Keys_{Type¦label}_Scalar
     procedure :: values                  =>Values_{Type¦label}_Scalar
     generic   :: exists     => Exists_{Type¦label}_Scalar_VS,Exists_{Type¦label}_Scalar_CH
     procedure :: size   =>Size_{Type¦label}_Scalar
     procedure :: destroy=>Destroy_{Type¦label}_Scalar
  end type {Type¦label}ScalarHash

  interface {Type¦label}ScalarHash
     module procedure {Type¦label}ScalarHashConstructor
  end interface {Type¦label}ScalarHash
  
  ! The number of new elements by which to extend hashes that need to grow.
  integer, parameter :: hashSizeIncrement=128

contains

  function {Type¦label}ScalarHashConstructor() result(self)
     !% Constructor for scalar hashes.
     implicit none
     type({Type¦label}ScalarHash) :: self

     call Initialize_{Type¦label}_Scalar(self)
     return
   end function {Type¦label}ScalarHashConstructor
  
  subroutine Initialize_{Type¦label}_Scalar(thisHash)
    !% Routine to initialize (or re-initialize) a hash.
    implicit none
    class({Type¦label}ScalarHash), intent(  out) :: thisHash

    thisHash%elementCount =0
    thisHash%allocatedSize=0
    if (allocated(thisHash%hashValues)) deallocate(thisHash%hashValues)
    if (allocated(thisHash%hashKeys  )) deallocate(thisHash%hashKeys  )
    return
  end subroutine Initialize_{Type¦label}_Scalar

  integer function Size_{Type¦label}_Scalar(thisHash)
    !% Returns the number of elements in the specified {\normalfont \ttfamily Hash}.
    implicit none
    class({Type¦label}ScalarHash), intent(in   ) :: thisHash

    Size_{Type¦label}_Scalar=thisHash%elementCount
    return
  end function Size_{Type¦label}_Scalar

  logical function Exists_{Type¦label}_Scalar_CH(thisHash,keyCH)
    !% Returns true if the specified {\normalfont \ttfamily key} exists in the specified {\normalfont \ttfamily thisHash}, false otherwise.
    implicit none
    class    ({Type¦label}ScalarHash), intent(in   ) :: thisHash
    character(len=*            ), intent(in   ) :: keyCH
    type     (varying_string   ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    Exists_{Type¦label}_Scalar_CH=Exists_{Type¦label}_Scalar_VS(thisHash,key)
    return
  end function Exists_{Type¦label}_Scalar_CH

  logical function Exists_{Type¦label}_Scalar_VS(thisHash,key)
    !% Returns true if the specified {\normalfont \ttfamily key} exists in the specified {\normalfont \ttfamily thisHash}, false otherwise.
    implicit none
    class({Type¦label}ScalarHash), intent(in   ) :: thisHash
    type (varying_string   ), intent(in   ) :: key

    if (thisHash%elementCount > 0) then
       Exists_{Type¦label}_Scalar_VS=any(thisHash%hashKeys(1:thisHash%elementCount) == key)
    else
       Exists_{Type¦label}_Scalar_VS=.false.
    end if
    return
  end function Exists_{Type¦label}_Scalar_VS

  subroutine Delete_{Type¦label}_Scalar_CH(thisHash,keyCH)
    !% Deletes entry {\normalfont \ttfamily key} from {\normalfont \ttfamily thisHash}.
    implicit none
    character(len=*            ), intent(in   ) :: keyCH
    class    ({Type¦label}ScalarHash), intent(inout) :: thisHash
    type     (varying_string   ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    call Delete_{Type¦label}_Scalar_VS(thisHash,key)
    return
  end subroutine Delete_{Type¦label}_Scalar_CH

  subroutine Delete_{Type¦label}_Scalar_VS(thisHash,key)
    !% Deletes entry {\normalfont \ttfamily key} from {\normalfont \ttfamily Hash}.
    use Arrays_Search
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
    implicit none
    type   (varying_string   ), intent(in   ) :: key
    class  ({Type¦label}ScalarHash), intent(inout) :: thisHash
    integer(c_size_t         )   , save          :: iKey

    if (Exists_{Type¦label}_Scalar_VS(thisHash,key)) then
       iKey=Search_Array(thisHash%hashKeys(1:thisHash%elementCount),key)
       thisHash%hashKeys  (ikey:thisHash%elementCount-1)=thisHash%hashKeys  (ikey+1:thisHash%elementCount)
       thisHash%hashValues(ikey:thisHash%elementCount-1)=thisHash%hashValues(ikey+1:thisHash%elementCount)
       thisHash%elementCount                        =thisHash%elementCount-1
    else
       call Galacticus_Error_Report('Delete_{Type¦label}_Scalar_VS','key '''//char(key)//''' does not exist in hash')
    end if
    return
  end subroutine Delete_{Type¦label}_Scalar_VS

  function Key_{Type¦label}_Scalar_I(thisHash,indexValue) result (key)
    !% Returns the key of entry number {\normalfont \ttfamily index} in {\normalfont \ttfamily thisHash}.
    implicit none
    type   (varying_string   )                :: key
    integer                   , intent(in   ) :: indexValue
    class  ({Type¦label}ScalarHash), intent(in   ) :: thisHash

    key=thisHash%hashKeys(indexValue)
    return
  end function Key_{Type¦label}_Scalar_I

  subroutine Keys_{Type¦label}_Scalar(thisHash,keys)
    !% Returns an array of all keys in {\normalfont \ttfamily thisHash}.
    implicit none
    type (varying_string   ), allocatable, dimension(:), intent(inout) :: keys
    class({Type¦label}ScalarHash)                           , intent(in   ) :: thisHash

    if (allocated(keys)) deallocate(keys)
    allocate(keys(thisHash%elementCount))
    keys=thisHash%hashKeys(1:thisHash%elementCount)
    return
  end subroutine Keys_{Type¦label}_Scalar

  subroutine Values_{Type¦label}_Scalar(thisHash,values)
    !% Returns an array of all values in {\normalfont \ttfamily thisHash}.
    implicit none
    {Type¦intrinsic}                   , allocatable, dimension(:), intent(inout) :: values
    class  ({Type¦label}ScalarHash)                           , intent(in   ) :: thisHash

    if (allocated(values)) deallocate(values)
    allocate(values(thisHash%elementCount))
    values=thisHash%hashValues(1:thisHash%elementCount)
    return
  end subroutine Values_{Type¦label}_Scalar

   function Value_{Type¦label}_Scalar_I(thisHash,indexValue)
    !% Returns the value of entry number {\normalfont \ttfamily index} in {\normalfont \ttfamily Hash}.
    implicit none
    {Type¦intrinsic}                                        :: Value_{Type¦label}_Scalar_I
    class           ({Type¦label}ScalarHash), intent(in   ) :: thisHash
    integer                                 , intent(in   ) :: indexValue

    Value_{Type¦label}_Scalar_I=thisHash%hashValues(indexValue)
    return
  end function Value_{Type¦label}_Scalar_I

  function Value_{Type¦label}_Scalar_CH(thisHash,keyCH)
    !% Returns the value of {\normalfont \ttfamily Key} in {\normalfont \ttfamily Hash}.
    implicit none
    {Type¦intrinsic}                                        :: Value_{Type¦label}_Scalar_CH
    character       (len=*                 ), intent(in   ) :: keyCH
    class           ({Type¦label}ScalarHash), intent(in   ) :: thisHash
    type            (varying_string        ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    Value_{Type¦label}_Scalar_CH=Value_{Type¦label}_Scalar_VS(thisHash,key)
    return
  end function Value_{Type¦label}_Scalar_CH

  function Value_{Type¦label}_Scalar_VS(thisHash,key)
    !% Returns the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily thisHash}.
    use Arrays_Search
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
    implicit none
    {Type¦intrinsic}                                        :: Value_{Type¦label}_Scalar_VS
    class           ({Type¦label}ScalarHash), intent(in   ) :: thisHash
    type            (varying_string        ), intent(in   ) :: key
    integer         (c_size_t              )                :: iKey

    if (Exists_{Type¦label}_Scalar_VS(thisHash,key)) then
       iKey=Search_Array(thisHash%hashKeys(1:thisHash%elementCount),key)
       Value_{Type¦label}_Scalar_VS=thisHash%hashValues(iKey)
    else
       Value_{Type¦label}_Scalar_VS=0
       call Galacticus_Error_Report('Value_{Type¦label}_Scalar','key '''//char(key)//''' does not exist in hash')
    end if
    return
  end function Value_{Type¦label}_Scalar_VS

  subroutine Set_{Type¦label}_Scalar_CH(thisHash,keyCH,value)
    !% Sets the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily thisHash} to {\normalfont \ttfamily value}.
    implicit none
    {Type¦intrinsic}                        , intent(in   ) :: value
    character       (len=*                 ), intent(in   ) :: keyCH
    class           ({Type¦label}ScalarHash), intent(inout) :: thisHash
    type            (varying_string        ), save          :: key

    key=trim(keyCH)
    call Set_{Type¦label}_Scalar_VS(thisHash,key,value)
    return
  end subroutine Set_{Type¦label}_Scalar_CH

  subroutine Set_{Type¦label}_Scalar_VS(thisHash,key,value)
    !% Sets the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily thisHash} to {\normalfont \ttfamily value}.
    use Arrays_Search
    use, intrinsic :: ISO_C_Binding
    implicit none
    {Type¦intrinsic}                        , intent(in   )               :: Value
    type            (varying_string        ), intent(in   )               :: Key
    class           ({Type¦label}ScalarHash), intent(inout)               :: thisHash
    integer         (c_size_t              )                              :: iKey
    logical                                                               :: keyExists
    {Type¦intrinsic}                        , allocatable  , dimension(:) :: valuesTemporary
    type            (varying_string        ), allocatable  , dimension(:) :: keysTemporary

    ! Check if key already exists.
    if (thisHash%elementCount > 0) then
       keyExists=any(thisHash%hashKeys(1:thisHash%elementCount) == key)
    else
       keyExists=.false.
    end if
    if (keyExists) then
       iKey=Search_Array(thisHash%hashKeys(1:thisHash%elementCount),key)
       thisHash%hashValues(iKey)=value
    else
       ! Increase hash size if necessary.
       if (thisHash%elementCount == thisHash%allocatedSize) then
          if (thisHash%allocatedSize > 0) then
             allocate(valuesTemporary(thisHash%allocatedSize))
             allocate(keysTemporary  (thisHash%allocatedSize))
             valuesTemporary=thisHash%hashValues
             keysTemporary  =thisHash%hashKeys
             deallocate(thisHash%hashValues)
             deallocate(thisHash%hashKeys  )
             thisHash%allocatedSize=thisHash%allocatedSize+hashSizeIncrement
             allocate(thisHash%hashValues(thisHash%allocatedSize))
             allocate(thisHash%hashKeys  (thisHash%allocatedSize))
             thisHash%hashValues(1:size(valuesTemporary))=valuesTemporary
             thisHash%hashKeys  (1:size(valuesTemporary))=keysTemporary
             deallocate(valuesTemporary)
             deallocate(keysTemporary  )
          else
             thisHash%allocatedSize=hashSizeIncrement
             allocate(thisHash%hashValues(thisHash%allocatedSize))
             allocate(thisHash%hashKeys  (thisHash%allocatedSize))
          end if
       end if
       if (thisHash%elementCount > 0) then
          iKey=Search_Array(thisHash%hashKeys(1:thisHash%elementCount),key)
       else
          iKey=1
       end if
       if (iKey > thisHash%elementCount) then
          ! Insert at end.
          thisHash%elementCount                       =thisHash%elementCount+1
          thisHash%hashKeys    (thisHash%elementCount)=key
          thisHash%hashValues  (thisHash%elementCount)=value
       else
          ! Shift array then insert.
          thisHash%hashKeys        (iKey+2:thisHash%elementCount+1)=thisHash%hashKeys  (iKey+1:thisHash%elementCount)
          thisHash%hashValues      (iKey+2:thisHash%elementCount+1)=thisHash%hashValues(iKey+1:thisHash%elementCount)
          thisHash%hashKeys        (iKey+1                        )=key
          thisHash%hashValues      (iKey+1                        )=value
          thisHash%elementCount                                    =thisHash%elementCount+1
       end if
    end if
    return
  end subroutine Set_{Type¦label}_Scalar_VS

  subroutine Destroy_{Type¦label}_Scalar(thisHash)
    !% Destroys {\normalfont \ttfamily thisHash}.
    implicit none
    class  ({Type¦label}ScalarHash), intent(inout) :: thisHash
    integer                                        :: i

    if (allocated(thisHash%hashValues)) deallocate(thisHash%hashValues)
    if (allocated(thisHash%hashKeys  )) then
       do i=1,size(thisHash%hashKeys)
          call thisHash%hashKeys(i)%destroy()
       end do
       deallocate(thisHash%hashKeys)
    end if
    return
  end subroutine Destroy_{Type¦label}_Scalar

  subroutine {Type¦label}ScalarDestructor(self)
    !% Destroys {\normalfont \ttfamily thisHash}.
    implicit none
    type({Type¦label}ScalarHash), intent(inout) :: self

    call Destroy_{Type¦label}_Scalar(self)
    return
  end subroutine {Type¦label}ScalarDestructor

end module Hashes
