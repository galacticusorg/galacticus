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
Contains a module which implements "dictionaries" (i.e. associative arrays).
!!}

module Dictionaries
  !!{RST
  Implements "dictionaries" (i.e. associative arrays).
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
  public :: {Type¦label}Dictionary

  type :: {Type¦label}Container
     private
     {Type¦intrinsic}{Type¦attributes} :: object {Type¦initializor}
  end type {Type¦label}Container

  type :: {Type¦label}Dictionary
     !!{RST
     Derived type for Type¦label dictionaries.
     !!}
     private
     integer                                                   :: allocatedSize=0, elementCount=0
     integer(c_size_t             )                            :: indexPrevious
     type   ({Type¦label}Container), allocatable, dimension(:) :: dictionaryValues
     type   (varying_string       ), allocatable, dimension(:) :: dictionaryKeys
     type   (varying_string       )                            :: keyPrevious
   contains
     !![
     <methods docformat="rst">
       <method description="Initialize the dictionary."                                                                      method="initialize"   />
       <method description="Set the value of a key in the dictionary."                                                       method="set"          />
       <method description="Delete a key from the dictionary."                                                               method="delete"       />
       <method description="Return the value for the given key."                                                       method="value"        />
       <method description="Return the key of the ``indexValue``\ :math:`^\mathrm{th}` entry in the dictionary."                  method="key"          />
       <method description="Return an array of all keys in the dictionary."                                                  method="keys"         />
       <method description="Return an array of all values in the dictionary."                                                method="values"       />
       <method description="Return true if the specified key exists in the dictionary."                                      method="exists"       />
       <method description="Return the number of keys in the dictionary."                                                    method="size"         />
       <method description="Destroy the dictionary."                                                                         method="destroy"      />
       <method description="Assign dictionary objects."                                                                      method="assignment(=)"/>
     </methods>
     !!]
     final     ::                             {Type¦label}Destructor
     procedure ::                             {Type¦label}Assign
     generic   :: assignment(=)            => {Type¦label}Assign
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
  end type {Type¦label}Dictionary

  interface {Type¦label}Dictionary
     module procedure {Type¦label}DictionaryConstructor
  end interface {Type¦label}Dictionary

  ! The number of elements allocated when a dictionary is first grown; subsequent growth is geometric (doubling).
  integer, parameter :: dictionarySizeIncrement=128

contains

  function {Type¦label}DictionaryConstructor() result(self)
     !!{RST
     Constructor for scalar dictionaries.
     !!}
     implicit none
     type({Type¦label}Dictionary) :: self

     call {Type¦label}Initialize(self)
     return
   end function {Type¦label}DictionaryConstructor

  subroutine {Type¦label}Initialize(self)
    !!{RST
    Routine to initialize (or re-initialize) a dictionary.
    !!}
  use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class({Type¦label}Dictionary), intent(  out) :: self

    self%elementCount = 0
    self%allocatedSize= 0
    self%indexPrevious=-1
    self%keyPrevious  =''
    if (allocated(self%dictionaryValues)) deallocate(self%dictionaryValues)
    if (allocated(self%dictionaryKeys  )) deallocate(self%dictionaryKeys  )
    return
  end subroutine {Type¦label}Initialize

  subroutine {Type¦label}Assign(to,from)
    !!{RST
    Assignment operator for dictionaries.
    !!}
    implicit none
    class  ({Type¦label}Dictionary), intent(  out) :: to
    class  ({Type¦label}Dictionary), intent(in   ) :: from
    integer                                        :: i

    to%allocatedSize=from%allocatedSize
    to%elementCount =from%elementCount
    to%indexPrevious=from%indexPrevious
    to%keyPrevious  =from%keyPrevious
    if (allocated(to  %dictionaryValues)) deallocate(to%dictionaryValues                       )
    if (allocated(to  %dictionaryKeys  )) deallocate(to%dictionaryKeys                         )
    if (allocated(from%dictionaryValues)) then
       allocate(to%dictionaryValues(size(from%dictionaryValues)))
       do i=1,size(from%dictionaryValues)
          to%dictionaryValues(i)=from%dictionaryValues(i)
       end do
    end if
    if (allocated(from%dictionaryKeys  )) then
       allocate(to%dictionaryKeys(size(from%dictionaryKeys)))
       do i=1,size(from%dictionaryKeys)
          to%dictionaryKeys(i)=from%dictionaryKeys(i)
       end do
    end if
    return
  end subroutine {Type¦label}Assign

  integer function {Type¦label}Size(self)
    !!{RST
    Returns the number of elements in the specified ``Dictionary``.
    !!}
    implicit none
    class({Type¦label}Dictionary), intent(in   ) :: self

    {Type¦label}Size=self%elementCount
    return
  end function {Type¦label}Size

  logical function {Type¦label}ExistsChar(self,keyCH)
    !!{RST
    Returns true if the specified ``key`` exists in the specified ``self``, false otherwise.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    class    ({Type¦label}Dictionary), intent(in   ) :: self
    character(len=*                 ), intent(in   ) :: keyCH
    type     (varying_string        ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    {Type¦label}ExistsChar={Type¦label}ExistsVarStr(self,key)
    return
  end function {Type¦label}ExistsChar

  function {Type¦label}Lookup(self,key) result(iKey)
    !!{RST
    Return the index of ``key`` in the sorted keys of ``self``, or zero if ``key`` does not exist. Uses a binary search, so is :math:`\mathcal{O}(\log N)` in the number of entries.
    !!}
    use            :: Arrays_Search     , only : searchArray
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : operator(==)
    implicit none
    integer(c_size_t              )                :: iKey
    class  ({Type¦label}Dictionary), intent(in   ) :: self
    type   (varying_string        ), intent(in   ) :: key

    if (self%elementCount < 1) then
       iKey=0_c_size_t
       return
    end if
    iKey=searchArray(self%dictionaryKeys(1:self%elementCount),key)
    if (iKey < 1_c_size_t .or. iKey > self%elementCount) then
       iKey=0_c_size_t
    else if (.not.(self%dictionaryKeys(iKey) == key)) then
       iKey=0_c_size_t
    end if
    return
  end function {Type¦label}Lookup

  logical function {Type¦label}ExistsVarStr(self,key)
    !!{RST
    Returns true if the specified ``key`` exists in the specified ``self``, false otherwise.
    !!}
    implicit none
    class({Type¦label}Dictionary), intent(in   ) :: self
    type (varying_string        ), intent(in   ) :: key

    {Type¦label}ExistsVarStr={Type¦label}Lookup(self,key) > 0_c_size_t
    return
  end function {Type¦label}ExistsVarStr

  subroutine {Type¦label}DeleteChar(self,keyCH)
    !!{RST
    Deletes entry ``key`` from ``self``.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*                 ), intent(in   ) :: keyCH
    class    ({Type¦label}Dictionary), intent(inout) :: self
    type     (varying_string        ), save          :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    call {Type¦label}DeleteVarStr(self,key)
    return
  end subroutine {Type¦label}DeleteChar

  subroutine {Type¦label}DeleteVarStr(self,key)
    !!{RST
    Deletes entry ``key`` from ``Dictionary``.
    !!}
    use            :: Arrays_Search     , only : searchArray
    use            :: Error             , only : Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : assignment(=), char
    implicit none
    type   (varying_string        ), intent(in   ) :: key
    class  ({Type¦label}Dictionary), intent(inout) :: self
    integer(c_size_t              ), save          :: iKey
    !$omp threadprivate(iKey)
    integer(c_size_t              )                :: i

    if ({Type¦label}ExistsVarStr(self,key)) then
       iKey=searchArray(self%dictionaryKeys(1:self%elementCount),key)
       do i=iKey,self%elementCount-1
          self%dictionaryKeys  (i)        =                 self%dictionaryKeys  (i+1)
          self%dictionaryValues(i)%object {Type¦assignment} self%dictionaryValues(i+1)%object
       end do
       self%elementCount=self%elementCount-1
       ! Unset memoized key.
       self%  keyPrevious=''
       self%indexPrevious=-1
    else
       call Error_Report('key '''//char(key)//''' does not exist in dictionary'//{introspection:location})
    end if
    return
  end subroutine {Type¦label}DeleteVarStr

  function {Type¦label}KeyInt(self,indexValue) result (key)
    !!{RST
    Returns the key of entry number ``index`` in ``self``.
    !!}
    implicit none
    type   (varying_string        )                :: key
    integer                        , intent(in   ) :: indexValue
    class  ({Type¦label}Dictionary), intent(in   ) :: self

    key=self%dictionaryKeys(indexValue)
    return
  end function {Type¦label}KeyInt

  subroutine {Type¦label}Keys(self,keys)
    !!{RST
    Returns an array of all keys in ``self``.
    !!}
    implicit none
    type (varying_string        ), allocatable, dimension(:), intent(inout) :: keys
    class({Type¦label}Dictionary)                           , intent(in   ) :: self

    if (allocated(keys)) deallocate(keys)
    allocate(keys(self%elementCount))
    keys=self%dictionaryKeys(1:self%elementCount)
    return
  end subroutine {Type¦label}Keys

  subroutine {Type¦label}Values(self,values)
    !!{RST
    Returns an array of all values in ``self``.
    !!}
    use :: Error, only : Error_Report
    implicit none
#if {Type¦match¦^(generic|rank\d+[a-zA-Z]+)$¦0¦1}
    {Type¦intrinsic}                        {Type¦attributes}, allocatable, dimension(:), intent(inout) :: values
    class           ({Type¦label}Dictionary)                                            , intent(in   ) :: self

    if (allocated(values)) deallocate(values)
    allocate(values(self%elementCount))
    values=self%dictionaryValues(1:self%elementCount)%object
#else
    integer                                                  , allocatable, dimension(:), intent(inout) :: values
    class           ({Type¦label}Dictionary)                                            , intent(in   ) :: self
    !$GLC attributes unused :: self, values

    call Error_Report('values method is not supported for generic dictionaries'//{introspection:location})
#endif
    return
  end subroutine {Type¦label}Values

  function {Type¦label}ValueInt(self,indexValue)
    !!{RST
    Returns the value of entry number ``index`` in ``Dictionary``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueInt
    class           ({Type¦label}Dictionary), intent(in   )   :: self
    integer                                 , intent(in   )   :: indexValue

    if (indexValue < 1 .or. indexValue > self%size()) call Error_Report('index is out of range'//{introspection:location})
    {Type¦label}ValueInt {Type¦assignment} self%dictionaryValues(indexValue)%object
    return
  end function {Type¦label}ValueInt

  function {Type¦label}ValueChar(self,keyCH)
    !!{RST
    Returns the value of ``Key`` in ``Dictionary``.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueChar
    character       (len=*                 ), intent(in   )   :: keyCH
    class           ({Type¦label}Dictionary), intent(in   )   :: self
    type            (varying_string        ), save            :: key
    !$omp threadprivate(key)
    key=trim(keyCH)
    {Type¦label}ValueChar {Type¦assignment} {Type¦label}ValueVarStr(self,key)
    return
  end function {Type¦label}ValueChar

  function {Type¦label}ValueVarStr(self,key)
    !!{RST
    Returns the value of ``key`` in ``self``. A single binary search (:math:`\mathcal{O}(\log N)`) is used to locate the key, short-circuited by the memoized previous key for repeated access.
    !!}
    use            :: Error             , only : Error_Report
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : char        , operator(==)
    implicit none
    {Type¦intrinsic}                  {Type¦attributes} :: {Type¦label}ValueVarStr
    class           ({Type¦label}Dictionary), intent(in   )   :: self
    type            (varying_string        ), intent(in   )   :: key
    integer         (c_size_t              )                  :: iKey

    if (key == self%keyPrevious) then
       {Type¦label}ValueVarStr {Type¦assignment} self%dictionaryValues(self%indexPrevious)%object
    else
       iKey={Type¦label}Lookup(self,key)
       if (iKey > 0_c_size_t) then
          {Type¦label}ValueVarStr {Type¦assignment} self%dictionaryValues(iKey)%object
       else
          {Type¦label}ValueVarStr {Type¦assignment} {Type¦null}
          call Error_Report('key '''//char(key)//''' does not exist in dictionary'//{introspection:location})
       end if
    end if
    return
  end function {Type¦label}ValueVarStr

  subroutine {Type¦label}SetChar(self,keyCH,value)
    !!{RST
    Sets the value of ``key`` in ``self`` to ``value``.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    {Type¦intrinsic}                        {Type¦argumentAttributes}, intent(in   ) :: value
    character       (len=*                 )                         , intent(in   ) :: keyCH
    class           ({Type¦label}Dictionary)                         , intent(inout) :: self
    type            (varying_string        )                         , save          :: key
    !$omp threadprivate(key)

    key=trim(keyCH)
    call {Type¦label}SetVarStr(self,key,value)
    return
  end subroutine {Type¦label}SetChar

  subroutine {Type¦label}SetVarStr(self,key,value)
    !!{RST
    Sets the value of ``key`` in ``self`` to ``value``.
    !!}
    use            :: Arrays_Search     , only : searchArray
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : assignment(=), char, operator(==)
    implicit none
    {Type¦intrinsic}                        {Type¦argumentAttributes}, intent(in   )               :: Value
    type            (varying_string        )                         , intent(in   )               :: Key
    class           ({Type¦label}Dictionary)                         , intent(inout)               :: self
    integer         (c_size_t              )                                                       :: iKey           , i
    logical                                                                                        :: keyExists      , keyChanged
    type            ({Type¦label}Container )                         , allocatable  , dimension(:) :: valuesTemporary
    type            (varying_string        )                         , allocatable  , dimension(:) :: keysTemporary

    ! Check if the key already exists.
    keyChanged=.true.
    if (self%elementCount > 0) then
       if (key == self%keyPrevious) then
          iKey      =self%indexPrevious
          keyExists =.true.
          keyChanged=.false.
       else
          iKey     =searchArray(self%dictionaryKeys(1:self%elementCount),key)
          if (iKey < 1 .or. iKey > self%elementCount) then
             keyExists=.false.
          else
             keyExists=self%dictionaryKeys(iKey) == key
          end if
       end if
    else
       iKey     =-1
       keyExists=.false.
    end if
    if (keyExists) then
#if {Type¦match¦^rank\d+[a-zA-Z]+$¦1¦0}
#if {Type¦match¦Ptr$¦0¦1}
       deallocate(self%dictionaryValues(iKey)%object)
#endif
#endif
       self%dictionaryValues(iKey)%object {Type¦assignment} value
       ! Set memoized key.
       if (keyChanged) then
          self%  keyPrevious=key
          self%indexPrevious=iKey
       end if
    else
       ! Increase dictionary size if necessary.
       if (self%elementCount == self%allocatedSize) then
          if (self%allocatedSize > 0) then
             allocate(valuesTemporary(self%allocatedSize))
             allocate(keysTemporary  (self%allocatedSize))
             valuesTemporary=self%dictionaryValues
             keysTemporary  =self%dictionaryKeys
             deallocate(self%dictionaryValues)
             deallocate(self%dictionaryKeys  )
             ! Grow geometrically (doubling) to give amortized O(N) build cost.
             self%allocatedSize=self%allocatedSize*2
             allocate(self%dictionaryValues(self%allocatedSize))
             allocate(self%dictionaryKeys  (self%allocatedSize))
             self%dictionaryValues(1:size(valuesTemporary))=valuesTemporary
             self%dictionaryKeys  (1:size(valuesTemporary))=keysTemporary
             deallocate(valuesTemporary)
             deallocate(keysTemporary  )
          else
             self%allocatedSize=dictionarySizeIncrement
             allocate(self%dictionaryValues(self%allocatedSize))
             allocate(self%dictionaryKeys  (self%allocatedSize))
          end if
       end if
       if (self%elementCount > 0) then
          iKey=searchArray(self%dictionaryKeys(1:self%elementCount),key)
       else
          iKey=1
       end if
       if (iKey > self%elementCount) then
          ! Insert at end.
          self%elementCount                                 =                 self%elementCount+1
          self%dictionaryKeys    (self%elementCount)        =                 key
          self%dictionaryValues  (self%elementCount)%object {Type¦assignment} value
          ! Set memoized key.
          self%  keyPrevious=key
          self%indexPrevious=self%elementCount
       else
          ! Shift array then insert.
          do i=self%elementCount+1,iKey+2,-1
             self%dictionaryKeys    (i     )        =                 self%dictionaryKeys    (i-1)
             self%dictionaryValues  (i     )        =                 self%dictionaryValues  (i-1)
          end do
          self   %dictionaryKeys    (iKey+1)        =                 key
          self   %dictionaryValues  (iKey+1)%object {Type¦assignment} value
          self   %elementCount                      =                 self%elementCount     +1
          ! Set memoized key.
          self%  keyPrevious=key
          self%indexPrevious=iKey+1
       end if
    end if
    return
  end subroutine {Type¦label}SetVarStr

  subroutine {Type¦label}Destroy(self)
    !!{RST
    Destroys ``self``.
    !!}
    implicit none
    class  ({Type¦label}Dictionary), intent(inout) :: self
    integer                                        :: i

    if (allocated(self%dictionaryValues)) deallocate(self%dictionaryValues)
    if (allocated(self%dictionaryKeys  )) then
       do i=1,size(self%dictionaryKeys)
          call self%dictionaryKeys(i)%destroy()
       end do
       deallocate(self%dictionaryKeys)
    end if
    return
  end subroutine {Type¦label}Destroy

  subroutine {Type¦label}Destructor(self)
    !!{RST
    Destroys ``self``.
    !!}
    implicit none
    type({Type¦label}Dictionary), intent(inout) :: self

    call {Type¦label}Destroy(self)
    return
  end subroutine {Type¦label}Destructor

end module Dictionaries
