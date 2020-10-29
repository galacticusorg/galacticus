!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: ISO_Varying_String, only : varying_string
  implicit none

  !# <generic identifier="Type">
  !#  <instance label="integer"              intrinsic="integer"           attributes=""                              argumentAttributes=""                          assignment="="  null="0"         />
  !#  <instance label="integerSizeT"         intrinsic="integer(c_size_t)" attributes=""                              argumentAttributes=""                          assignment="="  null="0_c_size_t"/>
  !#  <instance label="rank1IntegerSizeT"    intrinsic="integer(c_size_t)" attributes=", allocatable, dimension(:  )" argumentAttributes="         , dimension(:  )" assignment="="  null="0_c_size_t"/>
  !#  <instance label="rank2IntegerSizeT"    intrinsic="integer(c_size_t)" attributes=", allocatable, dimension(:,:)" argumentAttributes="         , dimension(:,:)" assignment="="  null="0_c_size_t"/>
  !#  <instance label="rank1IntegerSizeTPtr" intrinsic="integer(c_size_t)" attributes=", pointer    , dimension(:  )" argumentAttributes=", pointer, dimension(:  )" assignment="=>" null="null()"    />
  !#  <instance label="rank2IntegerSizeTPtr" intrinsic="integer(c_size_t)" attributes=", pointer    , dimension(:,:)" argumentAttributes=", pointer, dimension(:,:)" assignment="=>" null="null()"    />
  !#  <instance label="double"               intrinsic="double precision"  attributes=""                              argumentAttributes=""                          assignment="="  null="0.0d0"     />
  !#  <instance label="rank1Double"          intrinsic="double precision"  attributes=", allocatable, dimension(:  )" argumentAttributes="         , dimension(:  )" assignment="="  null="0.0d0"     />
  !#  <instance label="rank2Double"          intrinsic="double precision"  attributes=", allocatable, dimension(:,:)" argumentAttributes="         , dimension(:,:)" assignment="="  null="0.0d0"     />
  !#  <instance label="rank1DoublePtr"       intrinsic="double precision"  attributes=", pointer    , dimension(:  )" argumentAttributes=", pointer, dimension(:  )" assignment="=>" null="null()"    />
  !#  <instance label="rank2DoublePtr"       intrinsic="double precision"  attributes=", pointer    , dimension(:,:)" argumentAttributes=", pointer, dimension(:,:)" assignment="=>" null="null()"    />
  !#  <instance label="generic"              intrinsic="class(*)"          attributes=", pointer"                     argumentAttributes=", target"                  assignment="=>" null="null()"    />
  !# </generic>

  private
  public :: {Type¦label}Hash

  type :: {Type¦label}Container
  private
     {Type¦intrinsic}{Type¦attributes} :: object
  end type {Type¦label}Container

  type :: {Type¦label}Hash
     !% Derived type for {Type¦label} hashes.
     private
     integer                                                   :: allocatedSize, elementCount
     type   ({Type¦label}Container), allocatable, dimension(:) :: hashValues
     type   (varying_string       ), allocatable, dimension(:) :: hashKeys
   contains
     !@ <objectMethods>
     !@   <object>{Type¦label}Hash</object>
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
     !@     <description>Return the key of the {\normalfont \ttfamily indexValue}$^\mathrm{th}$ entry in the hash.</description>
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
     !@     <arguments>\textcolor{red}{\textless {Type¦intrinsic}{Type¦attributes}[:]\textgreater} values\arginout</arguments>
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
     final     ::                             {Type¦label}Destructor
     procedure :: initialize               => {Type¦label}Initialize
     procedure :: {Type¦label}SetVarStr
     procedure :: {Type¦label}SetChar
     generic   :: set                      => {Type¦label}SetVarStr   , {Type¦label}SetChar
     procedure :: {Type¦label}DeleteVarStr
     procedure :: {Type¦label}DeleteChar
     generic   :: delete                   => {Type¦label}DeleteVarStr, {Type¦label}DeleteChar
     procedure :: {Type¦label}ValueVarStr
     procedure :: {Type¦label}ValueChar
     procedure :: {Type¦label}ValueInt
     generic   :: value                    => {Type¦label}ValueVarStr , {Type¦label}ValueChar , {Type¦label}ValueInt
     procedure :: key                      => {Type¦label}KeyInt
     procedure :: {Type¦label}ExistsVarStr
     procedure :: {Type¦label}ExistsChar
     procedure :: keys                     => {Type¦label}Keys
     procedure :: values                   => {Type¦label}Values
     generic   :: exists                   => {Type¦label}ExistsVarStr, {Type¦label}ExistsChar
     procedure :: size                     => {Type¦label}Size
     procedure :: destroy                  => {Type¦label}Destroy
  end type {Type¦label}Hash

  interface {Type¦label}Hash
     module procedure {Type¦label}HashConstructor
  end interface {Type¦label}Hash

  ! The number of new elements by which to extend hashes that need to grow.
  integer, parameter :: hashSizeIncrement=128

contains

  function {Type¦label}HashConstructor() result(self)
     !% Constructor for scalar hashes.
     implicit none
     type({Type¦label}Hash) :: self

     call {Type¦label}Initialize(self)
     return
   end function {Type¦label}HashConstructor

  subroutine {Type¦label}Initialize(thisHash)
    !% Routine to initialize (or re-initialize) a hash.
    implicit none
    class({Type¦label}Hash), intent(  out) :: thisHash

    thisHash%elementCount =0
    thisHash%allocatedSize=0
    if (allocated(thisHash%hashValues)) deallocate(thisHash%hashValues)
    if (allocated(thisHash%hashKeys  )) deallocate(thisHash%hashKeys  )
    return
  end subroutine {Type¦label}Initialize

  integer function {Type¦label}Size(thisHash)
    !% Returns the number of elements in the specified {\normalfont \ttfamily Hash}.
    implicit none
    class({Type¦label}Hash), intent(in   ) :: thisHash

    {Type¦label}Size=thisHash%elementCount
    return
  end function {Type¦label}Size

  logical function {Type¦label}ExistsChar(thisHash,keyCH)
    !% Returns true if the specified {\normalfont \ttfamily key} exists in the specified {\normalfont \ttfamily thisHash}, false otherwise.
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class    ({Type¦label}Hash), intent(in   ) :: thisHash
    character(len=*           ), intent(in   ) :: keyCH
    type     (varying_string  ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    {Type¦label}ExistsChar={Type¦label}ExistsVarStr(thisHash,key)
    return
  end function {Type¦label}ExistsChar

  logical function {Type¦label}ExistsVarStr(thisHash,key)
    !% Returns true if the specified {\normalfont \ttfamily key} exists in the specified {\normalfont \ttfamily thisHash}, false otherwise.
    use :: ISO_Varying_String, only : operator(==)
    implicit none
    class({Type¦label}Hash), intent(in   ) :: thisHash
    type (varying_string  ), intent(in   ) :: key

    if (thisHash%elementCount > 0) then
       {Type¦label}ExistsVarStr=any(thisHash%hashKeys(1:thisHash%elementCount) == key)
    else
       {Type¦label}ExistsVarStr=.false.
    end if
    return
  end function {Type¦label}ExistsVarStr

  subroutine {Type¦label}DeleteChar(thisHash,keyCH)
    !% Deletes entry {\normalfont \ttfamily key} from {\normalfont \ttfamily thisHash}.
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*           ), intent(in   ) :: keyCH
    class    ({Type¦label}Hash), intent(inout) :: thisHash
    type     (varying_string  ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    call {Type¦label}DeleteVarStr(thisHash,key)
    return
  end subroutine {Type¦label}DeleteChar

  subroutine {Type¦label}DeleteVarStr(thisHash,key)
    !% Deletes entry {\normalfont \ttfamily key} from {\normalfont \ttfamily Hash}.
    use            :: Arrays_Search     , only : searchArray
    use            :: Galacticus_Error  , only : Galacticus_Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : char
    implicit none
    type   (varying_string  ), intent(in   ) :: key
    class  ({Type¦label}Hash), intent(inout) :: thisHash
    integer(c_size_t        ), save          :: iKey
    !$omp threadprivate(iKey)
    integer(c_size_t        )                :: i

    if ({Type¦label}ExistsVarStr(thisHash,key)) then
       iKey=searchArray(thisHash%hashKeys(1:thisHash%elementCount),key)
       do i=iKey,thisHash%elementCount-1
          thisHash%hashKeys  (i)        =                 thisHash%hashKeys  (i+1)
          thisHash%hashValues(i)%object {Type¦assignment} thisHash%hashValues(i+1)%object
       end do
       thisHash%elementCount=thisHash%elementCount-1
    else
       call Galacticus_Error_Report('key '''//char(key)//''' does not exist in hash'//{introspection:location})
    end if
    return
  end subroutine {Type¦label}DeleteVarStr

  function {Type¦label}KeyInt(thisHash,indexValue) result (key)
    !% Returns the key of entry number {\normalfont \ttfamily index} in {\normalfont \ttfamily thisHash}.
    implicit none
    type   (varying_string  )                :: key
    integer                  , intent(in   ) :: indexValue
    class  ({Type¦label}Hash), intent(in   ) :: thisHash

    key=thisHash%hashKeys(indexValue)
    return
  end function {Type¦label}KeyInt

  subroutine {Type¦label}Keys(thisHash,keys)
    !% Returns an array of all keys in {\normalfont \ttfamily thisHash}.
    implicit none
    type (varying_string  ), allocatable, dimension(:), intent(inout) :: keys
    class({Type¦label}Hash)                           , intent(in   ) :: thisHash

    if (allocated(keys)) deallocate(keys)
    allocate(keys(thisHash%elementCount))
    keys=thisHash%hashKeys(1:thisHash%elementCount)
    return
  end subroutine {Type¦label}Keys

  subroutine {Type¦label}Values(thisHash,values)
    !% Returns an array of all values in {\normalfont \ttfamily thisHash}.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
#if {Type¦match¦^(generic|rank\d+[a-zA-Z]+)$¦0¦1}
    {Type¦intrinsic}                  {Type¦attributes}, allocatable, dimension(:), intent(inout) :: values
    class           ({Type¦label}Hash)                                            , intent(in   ) :: thisHash

    if (allocated(values)) deallocate(values)
    allocate(values(thisHash%elementCount))
    values=thisHash%hashValues(1:thisHash%elementCount)%object
#else
    integer                                            , allocatable, dimension(:), intent(inout) :: values
    class           ({Type¦label}Hash)                                            , intent(in   ) :: thisHash
    !$GLC attributes unused :: thisHash, values

    call Galacticus_Error_Report('values method is not supported for generic hashes'//{introspection:location})
#endif
    return
  end subroutine {Type¦label}Values

  function {Type¦label}ValueInt(thisHash,indexValue)
    !% Returns the value of entry number {\normalfont \ttfamily index} in {\normalfont \ttfamily Hash}.
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueInt
    class           ({Type¦label}Hash), intent(in   )   :: thisHash
    integer                           , intent(in   )   :: indexValue

    {Type¦label}ValueInt {Type¦assignment} thisHash%hashValues(indexValue)%object
    return
  end function {Type¦label}ValueInt

  function {Type¦label}ValueChar(thisHash,keyCH)
    !% Returns the value of {\normalfont \ttfamily Key} in {\normalfont \ttfamily Hash}.
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueChar
    character       (len=*           ), intent(in   )   :: keyCH
    class           ({Type¦label}Hash), intent(in   )   :: thisHash
    type            (varying_string  ), save            :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    {Type¦label}ValueChar {Type¦assignment} {Type¦label}ValueVarStr(thisHash,key)
    return
  end function {Type¦label}ValueChar

  function {Type¦label}ValueVarStr(thisHash,key)
    !% Returns the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily thisHash}.
    use            :: Arrays_Search     , only : searchArray
    use            :: Galacticus_Error  , only : Galacticus_Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : char
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueVarStr
    class           ({Type¦label}Hash), intent(in   )   :: thisHash
    type            (varying_string  ), intent(in   )   :: key
    integer         (c_size_t        )                  :: iKey

    if ({Type¦label}ExistsVarStr(thisHash,key)) then
       iKey=searchArray(thisHash%hashKeys(1:thisHash%elementCount),key)
       {Type¦label}ValueVarStr {Type¦assignment} thisHash%hashValues(iKey)%object
    else
       {Type¦label}ValueVarStr {Type¦assignment} {Type¦null}
       call Galacticus_Error_Report('key '''//char(key)//''' does not exist in hash'//{introspection:location})
    end if
    return
  end function {Type¦label}ValueVarStr

  subroutine {Type¦label}SetChar(thisHash,keyCH,value)
    !% Sets the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily thisHash} to {\normalfont \ttfamily value}.
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    {Type¦intrinsic}                  {Type¦argumentAttributes}, intent(in   ) :: value
    character       (len=*           )                         , intent(in   ) :: keyCH
    class           ({Type¦label}Hash)                         , intent(inout) :: thisHash
    type            (varying_string  )                         , save          :: key
    !$omp threadprivate(key)

    key=trim(keyCH)
    call {Type¦label}SetVarStr(thisHash,key,value)
    return
  end subroutine {Type¦label}SetChar

  subroutine {Type¦label}SetVarStr(thisHash,key,value)
    !% Sets the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily thisHash} to {\normalfont \ttfamily value}.
    use            :: Arrays_Search     , only : searchArray
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : operator(==)
    implicit none
    {Type¦intrinsic}                       {Type¦argumentAttributes}, intent(in   )               :: Value
    type            (varying_string       )                         , intent(in   )               :: Key
    class           ({Type¦label}Hash     )                         , intent(inout)               :: thisHash
    integer         (c_size_t             )                                                       :: iKey           , i
    logical                                                                                       :: keyExists
    type            ({Type¦label}Container)                         , allocatable  , dimension(:) :: valuesTemporary
    type            (varying_string       )                         , allocatable  , dimension(:) :: keysTemporary

    ! Check if key already exists.
    if (thisHash%elementCount > 0) then
       keyExists=any(thisHash%hashKeys(1:thisHash%elementCount) == key)
    else
       keyExists=.false.
    end if
    if (keyExists) then
       iKey=searchArray(thisHash%hashKeys(1:thisHash%elementCount),key)
#if {Type¦match¦^rank\d+[a-zA-Z]+$¦1¦0}
       deallocate(thisHash%hashValues(iKey)%object)
#endif
       thisHash%hashValues(iKey)%object {Type¦assignment} value
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
          iKey=searchArray(thisHash%hashKeys(1:thisHash%elementCount),key)
       else
          iKey=1
       end if
       if (iKey > thisHash%elementCount) then
          ! Insert at end.
          thisHash%elementCount                               =                 thisHash%elementCount+1
          thisHash%hashKeys    (thisHash%elementCount)        =                 key
          thisHash%hashValues  (thisHash%elementCount)%object {Type¦assignment} value
       else
          ! Shift array then insert.
          do i=thisHash%elementCount+1,iKey+2,-1
             thisHash%hashKeys    (i     )        =                 thisHash%hashKeys    (i-1)
             thisHash%hashValues  (i     )        =                 thisHash%hashValues  (i-1)
          end do
          thisHash   %hashKeys    (iKey+1)        =                 key
          thisHash   %hashValues  (iKey+1)%object {Type¦assignment} value
          thisHash   %elementCount                =                 thisHash%elementCount     +1
       end if
    end if
    return
  end subroutine {Type¦label}SetVarStr

  subroutine {Type¦label}Destroy(thisHash)
    !% Destroys {\normalfont \ttfamily thisHash}.
    implicit none
    class  ({Type¦label}Hash), intent(inout) :: thisHash
    integer                                  :: i

    if (allocated(thisHash%hashValues)) deallocate(thisHash%hashValues)
    if (allocated(thisHash%hashKeys  )) then
       do i=1,size(thisHash%hashKeys)
          call thisHash%hashKeys(i)%destroy()
       end do
       deallocate(thisHash%hashKeys)
    end if
    return
  end subroutine {Type¦label}Destroy

  subroutine {Type¦label}Destructor(self)
    !% Destroys {\normalfont \ttfamily thisHash}.
    implicit none
    type({Type¦label}Hash), intent(inout) :: self

    call {Type¦label}Destroy(self)
    return
  end subroutine {Type¦label}Destructor

end module Hashes
