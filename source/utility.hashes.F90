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
Contains a module which implements ``hashes'' (i.e. associative arrays).
!!}

module Hashes
  !!{
  Implements ``hashes'' (i.e. associative arrays).
  !!}
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: ISO_Varying_String, only : var_str , varying_string
  implicit none

  !![
  <generic identifier="Type">
   <instance label="integer"              intrinsic="integer"              attributes=""                              argumentAttributes=""                          assignment="="  null="0"           initializor=""         />
   <instance label="integerSizeT"         intrinsic="integer(c_size_t)"    attributes=""                              argumentAttributes=""                          assignment="="  null="0_c_size_t"  initializor=""         />
   <instance label="rank1IntegerSizeT"    intrinsic="integer(c_size_t)"    attributes=", allocatable, dimension(:  )" argumentAttributes="         , dimension(:  )" assignment="="  null="0_c_size_t"  initializor=""         />
   <instance label="rank2IntegerSizeT"    intrinsic="integer(c_size_t)"    attributes=", allocatable, dimension(:,:)" argumentAttributes="         , dimension(:,:)" assignment="="  null="0_c_size_t"  initializor=""         />
   <instance label="rank1IntegerSizeTPtr" intrinsic="integer(c_size_t)"    attributes=", pointer    , dimension(:  )" argumentAttributes=", pointer, dimension(:  )" assignment="=>" null="null()"      initializor="=> null()"/>
   <instance label="rank2IntegerSizeTPtr" intrinsic="integer(c_size_t)"    attributes=", pointer    , dimension(:,:)" argumentAttributes=", pointer, dimension(:,:)" assignment="=>" null="null()"      initializor="=> null()"/>
   <instance label="double"               intrinsic="double precision"     attributes=""                              argumentAttributes=""                          assignment="="  null="0.0d0"       initializor=""         />
   <instance label="rank1Double"          intrinsic="double precision"     attributes=", allocatable, dimension(:  )" argumentAttributes="         , dimension(:  )" assignment="="  null="0.0d0"       initializor=""         />
   <instance label="rank2Double"          intrinsic="double precision"     attributes=", allocatable, dimension(:,:)" argumentAttributes="         , dimension(:,:)" assignment="="  null="0.0d0"       initializor=""         />
   <instance label="rank1DoublePtr"       intrinsic="double precision"     attributes=", pointer    , dimension(:  )" argumentAttributes=", pointer, dimension(:  )" assignment="=>" null="null()"      initializor="=> null()"/>
   <instance label="rank2DoublePtr"       intrinsic="double precision"     attributes=", pointer    , dimension(:,:)" argumentAttributes=", pointer, dimension(:,:)" assignment="=>" null="null()"      initializor="=> null()"/>
   <instance label="varyingString"        intrinsic="type(varying_string)" attributes=""                              argumentAttributes=""                          assignment="="  null="var_str('')" initializor=""         />
   <instance label="generic"              intrinsic="class(*)"             attributes=", pointer"                     argumentAttributes=", target"                  assignment="=>" null="null()"      initializor="=> null()"/>
  </generic>
  !!]

  private
  public :: {Type¦label}Hash

  type :: {Type¦label}Container
  private
     {Type¦intrinsic}{Type¦attributes} :: object {Type¦initializor}
  end type {Type¦label}Container

  type :: {Type¦label}Hash
     !!{
     Derived type for {Type¦label} hashes.
     !!}
     private
     integer                                                   :: allocatedSize=0, elementCount=0
     integer(c_size_t             )                            :: indexPrevious
     type   ({Type¦label}Container), allocatable, dimension(:) :: hashValues
     type   (varying_string       ), allocatable, dimension(:) :: hashKeys
     type   (varying_string       )                            :: keyPrevious
   contains
     !![
     <methods>
       <method description="Initialize the hash."                                                                      method="initialize"/>
       <method description="Set the value of a key in the hash."                                                       method="set"       />
       <method description="Delete a key from the hash."                                                               method="delete"    />
       <method description="Return the value for the given key."                                                       method="value"     />
       <method description="Return the key of the {\normalfont \ttfamily indexValue}$^\mathrm{th}$ entry in the hash." method="key"       />
       <method description="Return an array of all keys in the hash."                                                  method="keys"      />
       <method description="Return an array of all values in the hash."                                                method="values"    />
       <method description="Return true if the specified key exists in the hash."                                      method="exists"    />
       <method description="Return the number of keys in the hash."                                                    method="size"      />
       <method description="Destroy the hash."                                                                         method="destroy"   />
     </methods>
     !!]
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
     !!{
     Constructor for scalar hashes.
     !!}
     implicit none
     type({Type¦label}Hash) :: self

     call {Type¦label}Initialize(self)
     return
   end function {Type¦label}HashConstructor

  subroutine {Type¦label}Initialize(self)
    !!{
    Routine to initialize (or re-initialize) a hash.
    !!}
  use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class({Type¦label}Hash), intent(  out) :: self

    self%elementCount = 0
    self%allocatedSize= 0
    self%indexPrevious=-1
    self%keyPrevious  =''
    if (allocated(self%hashValues)) deallocate(self%hashValues)
    if (allocated(self%hashKeys  )) deallocate(self%hashKeys  )
    return
  end subroutine {Type¦label}Initialize

  integer function {Type¦label}Size(self)
    !!{
    Returns the number of elements in the specified {\normalfont \ttfamily Hash}.
    !!}
    implicit none
    class({Type¦label}Hash), intent(in   ) :: self

    {Type¦label}Size=self%elementCount
    return
  end function {Type¦label}Size

  logical function {Type¦label}ExistsChar(self,keyCH)
    !!{
    Returns true if the specified {\normalfont \ttfamily key} exists in the specified {\normalfont \ttfamily self}, false otherwise.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class    ({Type¦label}Hash), intent(in   ) :: self
    character(len=*           ), intent(in   ) :: keyCH
    type     (varying_string  ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    {Type¦label}ExistsChar={Type¦label}ExistsVarStr(self,key)
    return
  end function {Type¦label}ExistsChar

  logical function {Type¦label}ExistsVarStr(self,key)
    !!{
    Returns true if the specified {\normalfont \ttfamily key} exists in the specified {\normalfont \ttfamily self}, false otherwise.
    !!}
    use :: ISO_Varying_String, only : operator(==)
    implicit none
    class({Type¦label}Hash), intent(in   ) :: self
    type (varying_string  ), intent(in   ) :: key

    if (self%elementCount > 0) then
       {Type¦label}ExistsVarStr=any(self%hashKeys(1:self%elementCount) == key)
    else
       {Type¦label}ExistsVarStr=.false.
    end if
    return
  end function {Type¦label}ExistsVarStr

  subroutine {Type¦label}DeleteChar(self,keyCH)
    !!{
    Deletes entry {\normalfont \ttfamily key} from {\normalfont \ttfamily self}.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*           ), intent(in   ) :: keyCH
    class    ({Type¦label}Hash), intent(inout) :: self
    type     (varying_string  ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    call {Type¦label}DeleteVarStr(self,key)
    return
  end subroutine {Type¦label}DeleteChar

  subroutine {Type¦label}DeleteVarStr(self,key)
    !!{
    Deletes entry {\normalfont \ttfamily key} from {\normalfont \ttfamily Hash}.
    !!}
    use            :: Arrays_Search     , only : searchArray
    use            :: Error             , only : Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : assignment(=), char
    implicit none
    type   (varying_string  ), intent(in   ) :: key
    class  ({Type¦label}Hash), intent(inout) :: self
    integer(c_size_t        ), save          :: iKey
    !$omp threadprivate(iKey)
    integer(c_size_t        )                :: i

    if ({Type¦label}ExistsVarStr(self,key)) then
       iKey=searchArray(self%hashKeys(1:self%elementCount),key)
       do i=iKey,self%elementCount-1
          self%hashKeys  (i)        =                 self%hashKeys  (i+1)
          self%hashValues(i)%object {Type¦assignment} self%hashValues(i+1)%object
       end do
       self%elementCount=self%elementCount-1
       ! Unset memoized key.
       self%  keyPrevious=''
       self%indexPrevious=-1
    else
       call Error_Report('key '''//char(key)//''' does not exist in hash'//{introspection:location})
    end if
    return
  end subroutine {Type¦label}DeleteVarStr

  function {Type¦label}KeyInt(self,indexValue) result (key)
    !!{
    Returns the key of entry number {\normalfont \ttfamily index} in {\normalfont \ttfamily self}.
    !!}
    implicit none
    type   (varying_string  )                :: key
    integer                  , intent(in   ) :: indexValue
    class  ({Type¦label}Hash), intent(in   ) :: self

    key=self%hashKeys(indexValue)
    return
  end function {Type¦label}KeyInt

  subroutine {Type¦label}Keys(self,keys)
    !!{
    Returns an array of all keys in {\normalfont \ttfamily self}.
    !!}
    implicit none
    type (varying_string  ), allocatable, dimension(:), intent(inout) :: keys
    class({Type¦label}Hash)                           , intent(in   ) :: self

    if (allocated(keys)) deallocate(keys)
    allocate(keys(self%elementCount))
    keys=self%hashKeys(1:self%elementCount)
    return
  end subroutine {Type¦label}Keys

  subroutine {Type¦label}Values(self,values)
    !!{
    Returns an array of all values in {\normalfont \ttfamily self}.
    !!}
    use :: Error, only : Error_Report
    implicit none
#if {Type¦match¦^(generic|rank\d+[a-zA-Z]+)$¦0¦1}
    {Type¦intrinsic}                  {Type¦attributes}, allocatable, dimension(:), intent(inout) :: values
    class           ({Type¦label}Hash)                                            , intent(in   ) :: self

    if (allocated(values)) deallocate(values)
    allocate(values(self%elementCount))
    values=self%hashValues(1:self%elementCount)%object
#else
    integer                                            , allocatable, dimension(:), intent(inout) :: values
    class           ({Type¦label}Hash)                                            , intent(in   ) :: self
    !$GLC attributes unused :: self, values

    call Error_Report('values method is not supported for generic hashes'//{introspection:location})
#endif
    return
  end subroutine {Type¦label}Values

  function {Type¦label}ValueInt(self,indexValue)
    !!{
    Returns the value of entry number {\normalfont \ttfamily index} in {\normalfont \ttfamily Hash}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueInt
    class           ({Type¦label}Hash), intent(in   )   :: self
    integer                           , intent(in   )   :: indexValue

    if (indexValue < 1 .or. indexValue > self%size()) call Error_Report('index is out of range'//{introspection:location})
    {Type¦label}ValueInt {Type¦assignment} self%hashValues(indexValue)%object
    return
  end function {Type¦label}ValueInt

  function {Type¦label}ValueChar(self,keyCH)
    !!{
    Returns the value of {\normalfont \ttfamily Key} in {\normalfont \ttfamily Hash}.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueChar
    character       (len=*           ), intent(in   )   :: keyCH
    class           ({Type¦label}Hash), intent(in   )   :: self
    type            (varying_string  ), save            :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    {Type¦label}ValueChar {Type¦assignment} {Type¦label}ValueVarStr(self,key)
    return
  end function {Type¦label}ValueChar

  function {Type¦label}ValueVarStr(self,key)
    !!{
    Returns the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily self}.
    !!}
    use            :: Arrays_Search     , only : searchArray
    use            :: Error             , only : Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : char        , operator(==)
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueVarStr
    class           ({Type¦label}Hash), intent(in   )   :: self
    type            (varying_string  ), intent(in   )   :: key
    integer         (c_size_t        )                  :: iKey

    if (key == self%keyPrevious) then
       {Type¦label}ValueVarStr {Type¦assignment} self%hashValues(self%indexPrevious)%object
    else if ({Type¦label}ExistsVarStr(self,key)) then
       iKey=searchArray(self%hashKeys(1:self%elementCount),key)
       {Type¦label}ValueVarStr {Type¦assignment} self%hashValues(iKey)%object
    else
       {Type¦label}ValueVarStr {Type¦assignment} {Type¦null}
       call Error_Report('key '''//char(key)//''' does not exist in hash'//{introspection:location})
    end if
    return
  end function {Type¦label}ValueVarStr

  subroutine {Type¦label}SetChar(self,keyCH,value)
    !!{
    Sets the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily self} to {\normalfont \ttfamily value}.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    {Type¦intrinsic}                  {Type¦argumentAttributes}, intent(in   ) :: value
    character       (len=*           )                         , intent(in   ) :: keyCH
    class           ({Type¦label}Hash)                         , intent(inout) :: self
    type            (varying_string  )                         , save          :: key
    !$omp threadprivate(key)

    key=trim(keyCH)
    call {Type¦label}SetVarStr(self,key,value)
    return
  end subroutine {Type¦label}SetChar

  subroutine {Type¦label}SetVarStr(self,key,value)
    !!{
    Sets the value of {\normalfont \ttfamily key} in {\normalfont \ttfamily self} to {\normalfont \ttfamily value}.
    !!}
    use            :: Arrays_Search     , only : searchArray
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : assignment(=), char, operator(==)
    implicit none
    {Type¦intrinsic}                       {Type¦argumentAttributes}, intent(in   )               :: Value
    type            (varying_string       )                         , intent(in   )               :: Key
    class           ({Type¦label}Hash     )                         , intent(inout)               :: self
    integer         (c_size_t             )                                                       :: iKey           , i
    logical                                                                                       :: keyExists      , keyChanged
    type            ({Type¦label}Container)                         , allocatable  , dimension(:) :: valuesTemporary
    type            (varying_string       )                         , allocatable  , dimension(:) :: keysTemporary

    ! Check if the key already exists.
    keyChanged=.true.
    if (self%elementCount > 0) then
       if (key == self%keyPrevious) then
          iKey      =self%indexPrevious
          keyExists =.true.
          keyChanged=.false.
       else
          iKey     =searchArray(self%hashKeys(1:self%elementCount),key)
          if (iKey < 1 .or. iKey > self%elementCount) then
             keyExists=.false.
          else
             keyExists=self%hashKeys(iKey) == key
          end if
       end if
    else
       iKey     =-1
       keyExists=.false.
    end if
    if (keyExists) then
#if {Type¦match¦^rank\d+[a-zA-Z]+$¦1¦0}
#if {Type¦match¦Ptr$¦0¦1}
       deallocate(self%hashValues(iKey)%object)
#endif
#endif
       self%hashValues(iKey)%object {Type¦assignment} value
       ! Set memoized key.
       if (keyChanged) then
          self%  keyPrevious=key
          self%indexPrevious=iKey
       end if
    else
       ! Increase hash size if necessary.
       if (self%elementCount == self%allocatedSize) then
          if (self%allocatedSize > 0) then
             allocate(valuesTemporary(self%allocatedSize))
             allocate(keysTemporary  (self%allocatedSize))
             valuesTemporary=self%hashValues
             keysTemporary  =self%hashKeys
             deallocate(self%hashValues)
             deallocate(self%hashKeys  )
             self%allocatedSize=self%allocatedSize+hashSizeIncrement
             allocate(self%hashValues(self%allocatedSize))
             allocate(self%hashKeys  (self%allocatedSize))
             self%hashValues(1:size(valuesTemporary))=valuesTemporary
             self%hashKeys  (1:size(valuesTemporary))=keysTemporary
             deallocate(valuesTemporary)
             deallocate(keysTemporary  )
          else
             self%allocatedSize=hashSizeIncrement
             allocate(self%hashValues(self%allocatedSize))
             allocate(self%hashKeys  (self%allocatedSize))
          end if
       end if
       if (self%elementCount > 0) then
          iKey=searchArray(self%hashKeys(1:self%elementCount),key)
       else
          iKey=1
       end if
       if (iKey > self%elementCount) then
          ! Insert at end.
          self%elementCount                               =                 self%elementCount+1
          self%hashKeys    (self%elementCount)        =                 key
          self%hashValues  (self%elementCount)%object {Type¦assignment} value
          ! Set memoized key.
          self%  keyPrevious=key
          self%indexPrevious=self%elementCount
       else
          ! Shift array then insert.
          do i=self%elementCount+1,iKey+2,-1
             self%hashKeys    (i     )        =                 self%hashKeys    (i-1)
             self%hashValues  (i     )        =                 self%hashValues  (i-1)
          end do
          self   %hashKeys    (iKey+1)        =                 key
          self   %hashValues  (iKey+1)%object {Type¦assignment} value
          self   %elementCount                =                 self%elementCount     +1
          ! Set memoized key.
          self%  keyPrevious=key
          self%indexPrevious=iKey+1
       end if
    end if
    return
  end subroutine {Type¦label}SetVarStr

  subroutine {Type¦label}Destroy(self)
    !!{
    Destroys {\normalfont \ttfamily self}.
    !!}
    implicit none
    class  ({Type¦label}Hash), intent(inout) :: self
    integer                                  :: i

    if (allocated(self%hashValues)) deallocate(self%hashValues)
    if (allocated(self%hashKeys  )) then
       do i=1,size(self%hashKeys)
          call self%hashKeys(i)%destroy()
       end do
       deallocate(self%hashKeys)
    end if
    return
  end subroutine {Type¦label}Destroy

  subroutine {Type¦label}Destructor(self)
    !!{
    Destroys {\normalfont \ttfamily self}.
    !!}
    implicit none
    type({Type¦label}Hash), intent(inout) :: self

    call {Type¦label}Destroy(self)
    return
  end subroutine {Type¦label}Destructor

end module Hashes
