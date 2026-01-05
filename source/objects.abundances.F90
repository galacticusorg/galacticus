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
Contains a module which defines the abundances structure used for describing elemental abundances in \glc.
!!}

module Abundances_Structure
  !!{
  Defines the abundances structure used for describing elemental abundances in \glc.
  !!}
  implicit none
  private
  public :: Abundances_Initialize    , Abundances_Names          , Abundances_Index_From_Name      , Abundances_Atomic_Index, &
       &    Abundances_Property_Count, Abundances_Get_Metallicity, Abundances_Mass_To_Mass_Fraction, abundances             , &
       &    operator(*)              , max                       , abs                             , operator(>)

  ! Interface to multiplication operators with abundances objects as their second argument.
  interface operator(*)
     module procedure Abundances_Multiply_Switched
  end interface operator(*)

  ! Interface to "greater than" operators.
  interface operator(>)
     module procedure Abundances_Greater_Than
  end interface operator(>)

  ! Interface to max() function for abundances objects.
  interface max
     module procedure Abundances_Max
  end interface max

  ! Interface to abs() function for abundances objects.
  interface abs
     module procedure Abundances_Abs
  end interface abs

  type abundances
     !!{
     The abundances structure used for describing elemental abundances in \glc.
     !!}
     private
     double precision                            :: metallicityValue
     double precision, allocatable, dimension(:) :: elementalValue
   contains
     !![
     <methods>
       <method description="Multiply an abundance by a scalar." method="operator(*)" />
       <method description="Divide an abundance by a scalar." method="operator(/)" />
       <method description="Add two abundances." method="operator(+)" />
       <method description="Subtract one abundance from another." method="operator(-)" />
       <method description="Returns the metallicity." method="metallicity" />
       <method description="Sets the metallicity to {\normalfont \ttfamily metallicity}." method="metallicitySet" />
       <method description="Converts abundance masses to mass fractions by dividing by the given {\normalfont \ttfamily mass} while ensuring that fractions are in the range 0--1." method="massToMassFraction" />
       <method description="Increment an abundances object." method="increment" />
       <method description="Return a count of the number of properties in a serialized abundances object." method="serializeCount" />
       <method description="Serialize an abundances object to an array." method="serialize" />
       <method description="Deserialize an abundances object from an array." method="deserialize" />
       <method description="Returns the hydrogen fraction by number." method="hydrogenNumberFraction" />
       <method description="Returns the hydrogen fraction by mass." method="hydrogenMassFraction" />
       <method description="Returns the helium fraction by mass." method="heliumMassFraction" />
       <method description="Returns the helium fraction by number." method="heliumNumberFraction" />
       <method description="Return true if an abundances object is zero." method="isZero" />
       <method description="Destroy an abundances object." method="destroy" />
       <method description="Reset an abundances object." method="reset" />
       <method description="Build an abundances object from a provided XML description." method="builder" />
       <method description="Dump an abundances object." method="dump" />
       <method description="Dump an abundances object to binary." method="dumpRaw" />
       <method description="Read an abundances object from binary." method="readRaw" />
       <method description="Set an abundances object to unity." method="setToUnity" />
       <method description="Store an abundances object in the output buffers." method="output" />
       <method description="Store an abundances object in the output buffers." method="postOutput" />
       <method description="Specify the count of an abundances object for output." method="outputCount" />
       <method description="Specify the names of abundance object properties for output." method="outputNames" />
       <method description="Returns the size of any non-static components of the type." method="nonStaticSizeOf" />
     </methods>
     !!]
     procedure         ::                           Abundances_Add
     procedure         ::                           Abundances_Subtract
     procedure         ::                           Abundances_Multiply
     procedure         ::                           Abundances_Divide
     generic           :: operator(+)            => Abundances_Add
     generic           :: operator(-)            => Abundances_Subtract
     generic           :: operator(*)            => Abundances_Multiply
     generic           :: operator(/)            => Abundances_Divide
     procedure         :: nonStaticSizeOf        => Abundances_Non_Static_Size_Of
     procedure         :: isZero                 => Abundances_Is_Zero
     procedure         :: destroy                => Abundances_Destroy
     procedure         :: reset                  => Abundances_Reset
     procedure         :: builder                => Abundances_Builder
     procedure         :: dump                   => Abundances_Dump
     procedure         :: dumpRaw                => Abundances_Dump_Raw
     procedure         :: readRaw                => Abundances_Read_Raw
     procedure         :: setToUnity             => Abundances_Set_To_Unity
     procedure, nopass :: serializeCount         => Abundances_Property_Count
     procedure         :: serialize              => Abundances_Serialize
     procedure         :: deserialize            => Abundances_Deserialize
     procedure         :: increment              => Abundances_Increment
     procedure         :: metallicity            => Abundances_Get_Metallicity
     procedure         :: metallicitySet         => Abundances_Set_Metallicity
     procedure         :: massToMassFraction     => Abundances_Mass_To_Mass_Fraction_Packed
     procedure         :: hydrogenNumberFraction => Abundances_Hydrogen_Number_Fraction
     procedure         :: hydrogenMassFraction   => Abundances_Hydrogen_Mass_Fraction
     procedure         :: heliumMassFraction     => Abundances_Helium_Mass_Fraction
     procedure         :: heliumNumberFraction   => Abundances_Helium_Number_Fraction
     procedure         :: output                 => Abundances_Output
     procedure         :: postOutput             => Abundances_Post_Output
     procedure         :: outputCount            => Abundances_Output_Count
     procedure         :: outputNames            => Abundances_Output_Names
  end type abundances

  interface abundances
     !!{
     Constructors for the {\normalfont \ttfamily abundances} class.
     !!}
     module procedure abundancesConstructorZero
  end interface abundances

  ! Count of the number of elements being tracked.
  integer                                                                    :: elementsCount     =0
  integer                                                                    :: propertyCount

  ! Names (two-letter labels) of elements to track.
  character       (len=3     )                   , allocatable, dimension(:) :: elementsToTrack

  ! Indices of elements as used in the Atomic_Data module.
  integer                                        , allocatable, dimension(:) :: elementsIndices

  ! Value used to indicate a zero metallicity on logarithmic scales.
  double precision            , parameter, public                            :: logMetallicityZero=-99.0d0

  ! Unit and zero abundances objects.
  type            (abundances)           , public                            :: unitAbundances            , zeroAbundances

  ! Enumeration specifying type of metallicity/abundance measure required.
  !![
  <enumeration>
   <name>metallicityType</name>
   <description>Used to specify the metallicity scale when working with {\normalfont \ttfamily abundances} objects.</description>
   <visibility>public</visibility>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <entry label="linearByMass"                  />
   <entry label="linearByNumber"                />
   <entry label="logarithmicByMassSolar"        />
   <entry label="logarithmicByNumberSolar"      />
   <entry label="linearByMassSolar"             />
   <entry label="linearByNumberSolar"           />
   <entry label="logarithmicByNumberSolarPlus12"/>
  </enumeration>
  !!]

  ! Enumeration used in determining how to update elemental abundances when metallicity is adjusted.
  !![
  <enumeration>
   <name>adjustElements</name>
   <description>Used to specify how elements should be adjusted when the metallicity of an {\normalfont \ttfamily abundances} object is changed.</description>
   <visibility>public</visibility>
   <entry label="none"   />
   <entry label="reset"  />
   <entry label="update" />
  </enumeration>
  !!]

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Abundances_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Abundances_Initialize(parameters_)
    !!{
    Initialize the {\normalfont \ttfamily abundanceStructure} object module. Determines which abundances are to be tracked.
    !!}
    use :: Atomic_Data      , only : Atom_Lookup
    use :: Input_Parameters , only : inputParameters
    implicit none
    type   (inputParameters), intent(inout) :: parameters_
    integer                                 :: iElement
    
    ! Determine how many elements we are required to track.
    elementsCount=parameters_%count('elementsToTrack',zeroIfNotPresent=.true.)
    ! Number of properties to track is one greater, as we always track total metallicity.
    propertyCount=elementsCount+1
    ! If tracking elements, read names of which ones to track.
    if (elementsCount > 0) then
       if (allocated(elementsToTrack)) deallocate(elementsToTrack)
       if (allocated(elementsIndices)) deallocate(elementsIndices)
       allocate(elementsToTrack(elementsCount))
       allocate(elementsIndices(elementsCount))
       !![
       <inputParameter>
         <name>elementsToTrack</name>
         <description>The names of the elements to be tracked.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]
       ! Validate the input names by looking them up in the list of atomic names.
       do iElement=1,elementsCount
          elementsIndices(iElement)=Atom_Lookup(shortLabel=elementsToTrack(iElement))
       end do
    end if
    ! Create zero and unit abundances objects.
    if (allocated(zeroAbundances%elementalValue)) deallocate(zeroAbundances%elementalValue)
    if (allocated(unitAbundances%elementalValue)) deallocate(unitAbundances%elementalValue)
    allocate(zeroAbundances%elementalValue(elementsCount))
    allocate(unitAbundances%elementalValue(elementsCount))
    zeroAbundances%metallicityValue=0.0d0
    zeroAbundances%  elementalValue=0.0d0
    unitAbundances%metallicityValue=1.0d0
    unitAbundances%  elementalValue=1.0d0    
    return
  end subroutine Abundances_Initialize

  function abundancesConstructorZero() result(self)
    !!{
    A constructor for {\normalfont \ttfamily abundances} objects which sets all content to zero.
    !!}
    implicit none
    type(abundances) :: self

    self=zeroAbundances
    return
  end function abundancesConstructorZero

  subroutine Abundances_Destroy(self)
    !!{
    Destroy an abundances object.
    !!}
    implicit none
    class(abundances), intent(inout) :: self

    if (allocated(self%elementalValue)) deallocate(self%elementalValue)
    return
  end subroutine Abundances_Destroy

  subroutine Abundances_Builder(self,abundancesDefinition)
    !!{
    Build a {\normalfont \ttfamily abundances} object from the given XML {\normalfont \ttfamily abundancesDefinition}.
    !!}
    use :: FoX_DOM, only : node                        , extractDataContent
    use :: Error  , only : Error_Report
    use :: IO_XML , only : XML_Get_Elements_By_Tag_Name, xmlNodeList
    implicit none
    class  (abundances ), intent(inout)              :: self
    type   (node       ), intent(in   ), pointer     :: abundancesDefinition
    type   (node       )               , pointer     :: abundance
    type   (xmlNodeList), dimension(:) , allocatable :: abundanceList
    integer                                          :: i

    ! Get the metallicity.
    !$omp critical (FoX_DOM_Access)
    call XML_Get_Elements_By_Tag_Name(abundancesDefinition,'metals',abundanceList)
    !$omp end critical (FoX_DOM_Access)
    if (size(abundanceList) >  1) call Error_Report('multiple metallicity values specified'//{introspection:location})
    if (size(abundanceList) == 1) then
       !$omp critical (FoX_DOM_Access)
       abundance => abundanceList(0)%element
       call extractDataContent(abundance,self%metallicityValue)
       !$omp end critical (FoX_DOM_Access)
    end if
    if (elementsCount > 0) then
       do i=1,elementsCount
          !$omp critical (FoX_DOM_Access)
          call XML_Get_Elements_By_Tag_Name(abundancesDefinition,trim(elementsToTrack(i)),abundanceList)
          !$omp end critical (FoX_DOM_Access)
          if (size(abundanceList) >  1) call Error_Report('multiple '//trim(elementsToTrack(i))//' values specified'//{introspection:location})
          if (size(abundanceList) == 1) then
             !$omp critical (FoX_DOM_Access)
             abundance => abundanceList(0)%element
             call extractDataContent(abundance,self%elementalValue(i))
             !$omp end critical (FoX_DOM_Access)
          end if
       end do
    end if
    return
  end subroutine Abundances_Builder

  subroutine Abundances_Dump(self,verbosityLevel)
    !!{
    Dump properties of an abundances object.
    !!}
    use :: Display           , only : displayMessage, enumerationVerbosityLevelType
    use :: ISO_Varying_String, only : assignment(=) , operator(//)                 , varying_string
    implicit none
    class    (abundances                   ), intent(in   ) :: self
    type     (enumerationVerbosityLevelType), intent(in   ) :: verbosityLevel
    integer                                                 :: i
    character(len=22                       )                :: label
    type     (varying_string               )                :: message

    write (label,'(e22.16)') self%metallicityValue
    message='metallicity: '//label
    call displayMessage(message,verbosityLevel)
    if (elementsCount > 0) then
       do i=1,elementsCount
          write (label,'(e22.16)') self%elementalValue(i)
          message=elementsToTrack(i)//':          '//label
          call displayMessage(message,verbosityLevel)
       end do
    end if
    return
  end subroutine Abundances_Dump

  subroutine Abundances_Dump_Raw(self,fileHandle)
    !!{
    Dump an abundances object to binary.
    !!}
    implicit none
    class  (abundances    ), intent(in   ) :: self
    integer                , intent(in   ) :: fileHandle

    write (fileHandle) self%metallicityValue
    if (elementsCount > 0) write (fileHandle) self%elementalValue
    return
  end subroutine Abundances_Dump_Raw

  subroutine Abundances_Read_Raw(self,fileHandle)
    !!{
    Read an abundances object from binary.
    !!}
    implicit none
    class  (abundances    ), intent(inout) :: self
    integer                , intent(in   ) :: fileHandle

    read (fileHandle) self%metallicityValue
    ! If no individual elements are tracked, our work is done.
    if (elementsCount == 0) return
    ! Ensure elemental values array exists.
    call Abundances_Allocate_Elemental_Values(self)
    ! Read values from file.
    read (fileHandle) self%elementalValue    
    return
  end subroutine Abundances_Read_Raw

  subroutine Abundances_Reset(self)
    !!{
    Reset an abundances object.
    !!}
    implicit none
    class(abundances), intent(inout) :: self

    call self%metallicitySet(0.0d0,adjustElements=adjustElementsReset)
    return
  end subroutine Abundances_Reset

  subroutine Abundances_Set_To_Unity(self)
    !!{
    Set an abundances object to unity.
    !!}
    implicit none
    class(abundances), intent(inout) :: self

    ! Ensure object is initialized.
    call Abundances_Allocate_Elemental_Values(self)
    ! Set values to unity.
    self                       %metallicityValue=1.0d0
    if (elementsCount > 0) self%  elementalValue=1.0d0
    return
  end subroutine Abundances_Set_To_Unity

  logical function Abundances_Is_Zero(self)
    !!{
    Test whether an abundances object is zero.
    !!}
    implicit none
    class(abundances), intent(in   ) :: self

    ! Detect if all abundances are zero.
    Abundances_Is_Zero=.true.
    if (self%metallicityValue /= 0.0d0) Abundances_Is_Zero=.false.
    if (elementsCount > 0 .and. allocated(self%elementalValue)) then
       if (any(self%elementalValue /= 0.0d0)) Abundances_Is_Zero=.false.
    end if
    return
  end function Abundances_Is_Zero

  function Abundances_Add(abundances1,abundances2)
    !!{
    Add two abundances objects.
    !!}
    implicit none
    type (abundances)                          :: Abundances_Add
    class(abundances), intent(in   )           :: abundances1
    class(abundances), intent(in   ), optional :: abundances2

    if (present(abundances2)) then
       Abundances_Add                       %metallicityValue=abundances1%metallicityValue+abundances2%metallicityValue
       if (elementsCount > 0) Abundances_Add%  elementalValue=abundances1%  elementalValue+abundances2%  elementalValue
    else
       Abundances_Add                       %metallicityValue=abundances1%metallicityValue
       if (elementsCount > 0) Abundances_Add%  elementalValue=abundances1%  elementalValue
    end if
    return
  end function Abundances_Add

  subroutine Abundances_Increment(self,increment)
    !!{
    Increment an abundances object.
    !!}
    implicit none
    class(abundances), intent(inout) :: self
    class(abundances), intent(in   ) :: increment

    self                       %metallicityValue=self%metallicityValue+increment%metallicityValue
    if (elementsCount > 0) self%elementalValue  =self%elementalValue  +increment%elementalValue
    return
  end subroutine Abundances_Increment

  function Abundances_Subtract(abundances1,abundances2)
    !!{
    Subtract two abundances objects.
    !!}
    implicit none
    type (abundances)                          :: Abundances_Subtract
    class(abundances), intent(in   )           :: abundances1
    class(abundances), intent(in   ), optional :: abundances2

    if (present(abundances2)) then
       Abundances_Subtract                       %metallicityValue= abundances1%metallicityValue-abundances2%metallicityValue
       if (elementsCount > 0) Abundances_Subtract%  elementalValue= abundances1%  elementalValue-abundances2%  elementalValue
    else
       Abundances_Subtract                       %metallicityValue=-abundances1%metallicityValue
       if (elementsCount > 0) Abundances_Subtract%  elementalValue=-abundances1%  elementalValue
    end if
    return
  end function Abundances_Subtract

  function Abundances_Multiply(abundances1,multiplier)
    !!{
    Multiply an abundances object by a scalar.
    !!}
    implicit none
    type            (abundances)                :: Abundances_Multiply
    class           (abundances), intent(in   ) :: abundances1
    double precision            , intent(in   ) :: multiplier

    Abundances_Multiply                       %metallicityValue=abundances1%metallicityValue*multiplier
    if (elementsCount > 0) Abundances_Multiply%  elementalValue=abundances1%  elementalValue*multiplier
    return
  end function Abundances_Multiply

  function Abundances_Multiply_Switched(multiplier,abundances1)
    !!{
    Multiply a scalar by an abundances object.
    !!}
    implicit none
    type            (abundances)                :: Abundances_Multiply_Switched
    type            (abundances), intent(in   ) :: abundances1
    double precision            , intent(in   ) :: multiplier

    Abundances_Multiply_Switched=Abundances_Multiply(abundances1,multiplier)
    return
  end function Abundances_Multiply_Switched

  logical function Abundances_Greater_Than(abundances1,abundances2)
    !!{
    Return an element-by-element ``$>$'' on two abundances objects.
    !!}
    implicit none
    type(abundances), intent(in   ) :: abundances1, abundances2

    Abundances_Greater_Than=abundances1%metallicityValue > abundances2%metallicityValue
    if (elementsCount > 0) Abundances_Greater_Than=Abundances_Greater_Than .and. all(abundances1%elementalValue > abundances2%elementalValue)
    return
  end function Abundances_Greater_Than

  function Abundances_Max(abundances1,abundances2)
    !!{
    Return an element-by-element {\normalfont \ttfamily max()} on two abundances objects.
    !!}
    implicit none
    type(abundances)                :: Abundances_Max
    type(abundances), intent(in   ) :: abundances1   , abundances2

    Abundances_Max                       %metallicityValue=max(abundances1%metallicityValue,abundances2%metallicityValue)
    if (elementsCount > 0) Abundances_Max%  elementalValue=max(abundances1%  elementalValue,abundances2%  elementalValue)
    return
  end function Abundances_Max

  function Abundances_Abs(abundances1)
    !!{
    Return an element-by-element {\normalfont \ttfamily abs()} on an abundances objects.
    !!}
    implicit none
    type(abundances)                :: Abundances_Abs
    type(abundances), intent(in   ) :: abundances1

    Abundances_Abs                       %metallicityValue=abs(abundances1%metallicityValue)
    if (elementsCount > 0) Abundances_Abs%  elementalValue=abs(abundances1%  elementalValue)
    return
  end function Abundances_Abs

  function Abundances_Divide(abundances1,divisor)
    !!{
    Divide an abundances object by a scalar.
    !!}
    implicit none
    type            (abundances)                :: Abundances_Divide
    class           (abundances), intent(in   ) :: abundances1
    double precision            , intent(in   ) :: divisor

    Abundances_Divide                       %metallicityValue=abundances1%metallicityValue/divisor
    if (elementsCount > 0) Abundances_Divide%  elementalValue=abundances1%  elementalValue/divisor
    return
  end function Abundances_Divide

  integer function Abundances_Property_Count()
    !!{
    Return the number of properties required to track abundances. This is equal to the number of elements tracked, {\normalfont \ttfamily
    elementsCount}, plus one since we always track a total metallicity.
    !!}
    implicit none

    Abundances_Property_Count=propertyCount
    return
  end function Abundances_Property_Count

  function Abundances_Names(index)
    !!{
    Return a name for the specified entry in the abundances structure.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=), varying_string
    implicit none
    type   (varying_string)                :: Abundances_Names
    integer                , intent(in   ) :: index

    select case (index)
    case (1)
       Abundances_Names="Metals"
    case (2:)
       if (index <= elementsCount+1) then
          Abundances_Names=trim(elementsToTrack(index-1))
       else
          call Error_Report('index out of range'//{introspection:location})
       end if
    case default
       Abundances_Names=""
       call Error_Report('index out of range'//{introspection:location})
    end select
    return
  end function Abundances_Names

  integer function Abundances_Index_From_Name(name)
    !!{
    Return the index of an element in the elements array given its name.
    !!}
    implicit none
    character(len=*), intent(in   ) :: name
    integer                         :: i

    Abundances_Index_From_Name=-1
    do i=1,elementsCount
       if (trim(name) == elementsToTrack(i)) then
          Abundances_Index_From_Name=i+1
          exit
       end if
    end do
    return
  end function Abundances_Index_From_Name

  integer function Abundances_Atomic_Index(index)
    !!{
    Return the atomic index for the specified entry in the abundances structure.
    !!}
    use :: Error, only : Error_Report
    implicit none
    integer, intent(in   ) :: index

    select case (index)
    case (1)
       Abundances_Atomic_Index=0 ! Total metallicity.
    case (2:)
       if (index <= elementsCount+1) then
          Abundances_Atomic_Index=elementsIndices(index-1)
       else
          Abundances_Atomic_Index=-1
          call Error_Report('index out of range'//{introspection:location})
       end if
    case default
       Abundances_Atomic_Index=-1
       call Error_Report('index out of range'//{introspection:location})
    end select
    return
  end function Abundances_Atomic_Index

  subroutine Abundances_Allocate_Elemental_Values(self)
    !!{
    Ensure that the {\normalfont \ttfamily elementalValue} array in an {\normalfont \ttfamily abundances} is allocated.
    !!}
    implicit none
    type(abundances), intent(inout) :: self

    if (.not.allocated(self%elementalValue)) allocate(self%elementalValue(elementsCount))
    return
  end subroutine Abundances_Allocate_Elemental_Values

  subroutine Abundances_Deserialize(self,abundancesArray)
    !!{
    Pack abundances from an array into an abundances structure.
    !!}
    implicit none
    class           (abundances)              , intent(inout) :: self
    double precision            , dimension(:), intent(in   ) :: abundancesArray

    ! Extract metallicity from array.
    self%metallicityValue=abundancesArray(1)
    ! If no individual elements are tracked, our work is done.
    if (elementsCount == 0) return
    ! Ensure elemental values array exists.
    call Abundances_Allocate_Elemental_Values(self)
    ! Extract elemental values from array.
    self%elementalValue=abundancesArray(2:elementsCount+1)
    return
  end subroutine Abundances_Deserialize

  subroutine Abundances_Serialize(self,abundancesArray)
    !!{
    Pack abundances from an array into an abundances structure.
    !!}
    implicit none
    double precision            , dimension(:), intent(  out) :: abundancesArray(:)
    class           (abundances)              , intent(in   ) :: self

    ! Place metallicity into array.
    abundancesArray(1)=self%metallicityValue
    ! If no individual elements are tracked, our work is done.
    if (elementsCount == 0) return
    ! Place elemental values into arrays.
    if (allocated(self%elementalValue)) then
       abundancesArray(2:elementsCount+1)=self%elementalValue
    else
       abundancesArray(2:elementsCount+1)=0.0d0
    end if
    return
  end subroutine Abundances_Serialize

  double precision function Abundances_Get_Metallicity(self,metallicityType)
    !!{
    Return the metallicity of the {\normalfont \ttfamily self} structure.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    class  (abundances                    ), intent(in   )           :: self
    type   (enumerationMetallicityTypeType), intent(in   ), optional :: metallicityType
    !![
    <optionalArgument name="metallicityType" defaultsTo="metallicityTypeLinearByMass" />
    !!]
    
    Abundances_Get_Metallicity=self%metallicityValue
    select case (metallicityType_%ID)
    case (metallicityTypeLinearByMass%ID)
       ! Do nothing, this is what we compute by default.
    case (metallicityTypeLogarithmicByMassSolar%ID)
       ! Convert to a logarithmic metallicity by mass relative to Solar.
       if (Abundances_Get_Metallicity > 0.0d0) then
          Abundances_Get_Metallicity=log10(Abundances_Get_Metallicity/metallicitySolar)
       else
          Abundances_Get_Metallicity=logMetallicityZero
       end if
    case (metallicityTypeLinearByMassSolar%ID)
       ! Convert to metallicity by mass relative to Solar.
       Abundances_Get_Metallicity=Abundances_Get_Metallicity/metallicitySolar
    case default
       call Error_Report('metallicity type not supported'//{introspection:location})
    end select
    return
  end function Abundances_Get_Metallicity

  subroutine Abundances_Set_Metallicity(self,metallicity,metallicityType,adjustElements,abundanceIndex)
    !!{
    Set the metallicity of the {\normalfont \ttfamily self} structure to {\normalfont \ttfamily metallicity}.
    !!}
    use :: Atomic_Data                     , only : Atomic_Abundance, normalizationMetals
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    class           (abundances                    ), intent(inout)           :: self
    double precision                                , intent(in   )           :: metallicity
    integer                                         , intent(in   ), optional :: abundanceIndex
    type            (enumerationAdjustElementsType ), intent(in   ), optional :: adjustElements
    type            (enumerationMetallicityTypeType), intent(in   ), optional :: metallicityType
    integer                                                                   :: iElement
    double precision                                                          :: metallicityPrevious
    !![
    <optionalArgument name="adjustElements" defaultsTo="adjustElementsNone" />
    !!]

    ! Store the current metallicity.
    metallicityPrevious        =self%metallicityValue

    ! Store the current metallicity if necessary.
    if (elementsCount > 0 .and. adjustElements_ == adjustElementsUpdate) metallicityPrevious=self%metallicityValue

    ! Update the metallicity.
    self%metallicityValue=metallicity
    if (present(metallicityType)) then
       select case (metallicityType%ID)
       case (metallicityTypeLinearByMass          %ID)
          ! Do nothing, this is how we store metallicity.
       case (metallicityTypeLinearByMassSolar     %ID)
          self%metallicityValue=         self%metallicityValue *metallicitySolar
       case (metallicityTypeLogarithmicByMassSolar%ID)
          self%metallicityValue=(10.0d0**self%metallicityValue)*metallicitySolar
       case default
          call Error_Report('type not supported'//{introspection:location})
       end select
    end if

    ! Determine what we're requested to do with any other elements.
    if (elementsCount > 0) then
       select case (adjustElements_%ID)
       case (adjustElementsNone%ID)
          ! Do nothing to the elemental abundances in this case.
       case (adjustElementsReset%ID)
          ! Ensure that we have an abundanceIndex specified.
          if (self%metallicityValue /= 0.0d0 .and. .not.present(abundanceIndex)) &
               & call Error_Report('an abundance pattern must be specified in order to reset elemental abundances'//{introspection:location})
          ! Ensure elemental values array exists.
          call Abundances_Allocate_Elemental_Values(self)
          if (self%metallicityValue == 0.0d0) then
             self%elementalValue=0.0d0
          else
             do iElement=1,elementsCount
                self%elementalValue(iElement)=self%metallicityValue*Atomic_Abundance(abundanceIndex=abundanceIndex&
                     &,atomIndex=elementsIndices(iElement),normalization=normalizationMetals)
             end do
          end if
       case (adjustElementsUpdate%ID)
          ! Ensure that we have an abundanceIndex specified.
          if (.not.present(abundanceIndex)) call Error_Report('an abundance pattern must be specified in order to reset elemental abundances'//{introspection:location})
          ! Ensure elemental values array exists.
          call Abundances_Allocate_Elemental_Values(self)
          do iElement=1,elementsCount
             self%elementalValue(iElement)=self%elementalValue(iElement)+(self%metallicityValue&
                  &-metallicityPrevious)*Atomic_Abundance(abundanceIndex=abundanceIndex,atomIndex=elementsIndices(iElement)&
                  &,normalization=normalizationMetals)
          end do
        end select
    end if
    return
  end subroutine Abundances_Set_Metallicity

  subroutine Abundances_Mass_To_Mass_Fraction_Packed(self,mass)
    !!{
    Convert abundance masses to mass fractions by dividing by {\normalfont \ttfamily mass} while ensuring that the fractions remain within the range 0--1.
    !!}
    implicit none
    class           (abundances), intent(inout) :: self
    double precision            , intent(in   ) :: mass

    ! Scale metallicity first.
    if      (self%metallicityValue >  mass ) then
       self%metallicityValue=1.0d0
    else if (self%metallicityValue <= 0.0d0) then
       self%metallicityValue=0.0d0
    else
       self%metallicityValue=self%metallicityValue/mass
    end if

    ! Scale elemental abundances.
    if (elementsCount > 0) then
       where     (self%elementalValue >  mass )
          self%elementalValue=1.0d0
       elsewhere (self%elementalValue <= 0.0d0)
          self%elementalValue=0.0d0
       elsewhere
          self%elementalValue=self%elementalValue/mass
       end where
    end if
    return
  end subroutine Abundances_Mass_To_Mass_Fraction_Packed

  subroutine Abundances_Mass_To_Mass_Fraction(self,mass)
    !!{
    Convert abundance masses to mass fractions by dividing by {\normalfont \ttfamily mass} while ensuring that the fractions remain within the range 0--1.
    !!}
    implicit none
    double precision, dimension(:), intent(inout) :: self
    double precision              , intent(in   ) :: mass

    ! Scale abundances.
    where     (self >  mass )
       self=1.0d0
    elsewhere (self <= 0.0d0)
       self=0.0d0
    elsewhere
       self=self/mass
    end where
    return
  end subroutine Abundances_Mass_To_Mass_Fraction

  double precision function Abundances_Hydrogen_Mass_Fraction(self)
    !!{
    Returns the mass fraction of hydrogen.
    !!}
    use :: Numerical_Constants_Astronomical, only : hydrogenByMassPrimordial, hydrogenByMassSolar, metallicitySolar
    implicit none
    class           (abundances), intent(in   ) :: self
    double precision            , parameter     :: massFractionMinimum=0.7d0

    Abundances_Hydrogen_Mass_Fraction=min(                                                      &
         &                                max(                                                  &
         &                                    +self%metallicityValue                            &
         &                                    /     metallicitySolar                            &
         &                                    *(+hydrogenByMassSolar-hydrogenByMassPrimordial)  &
         &                                    +                      hydrogenByMassPrimordial , &
         &                                    +                      massFractionMinimum        &
         &                                   )                                                , &
         &                                    +                      hydrogenByMassPrimordial   &
         &                               )
    return
  end function Abundances_Hydrogen_Mass_Fraction

  double precision function Abundances_Helium_Mass_Fraction(self)
    !!{
    Returns the mass fraction of helium.
    !!}
    use :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial, heliumByMassSolar, metallicitySolar
    implicit none
    class(abundances), intent(in   ) :: self

    Abundances_Helium_Mass_Fraction=min(                                              &
         &                              +self%metallicityValue                        &
         &                              /     metallicitySolar                        &
         &                              *(+heliumByMassSolar-heliumByMassPrimordial)  &
         &                              +                    heliumByMassPrimordial , &
         &                              +                    heliumByMassPrimordial   &
         &                             )
    return
  end function Abundances_Helium_Mass_Fraction

  double precision function Abundances_Hydrogen_Number_Fraction(self)
    !!{
    Returns the number fraction of hydrogen.
    !!}
    use :: Numerical_Constants_Atomic, only : atomicMassHelium, atomicMassHydrogen
    implicit none
    class           (abundances), intent(in   ) :: self
    double precision                            :: numberHelium, numberHydrogen

    numberHydrogen=Abundances_Hydrogen_Mass_Fraction(self)/atomicMassHydrogen
    numberHelium  =Abundances_Helium_Mass_Fraction  (self)/atomicMassHelium
    Abundances_Hydrogen_Number_Fraction=numberHydrogen/(numberHydrogen+numberHelium)
    return
  end function Abundances_Hydrogen_Number_Fraction

  double precision function Abundances_Helium_Number_Fraction(self)
    !!{
    Returns the mass fraction of helium.
    !!}
    use :: Numerical_Constants_Atomic, only : atomicMassHelium, atomicMassHydrogen
    implicit none
    class           (abundances), intent(in   ) :: self
    double precision                            :: numberHelium, numberHydrogen

    numberHydrogen=Abundances_Hydrogen_Mass_Fraction(self)/atomicMassHydrogen
    numberHelium  =Abundances_Helium_Mass_Fraction  (self)/atomicMassHelium
    Abundances_Helium_Number_Fraction=numberHelium/(numberHydrogen+numberHelium)
    return
  end function Abundances_Helium_Number_Fraction

  subroutine Abundances_Output(self,integerProperty,integerBufferCount,integerProperties,doubleProperty,doubleBufferCount,doubleProperties,time,outputInstance)
    !!{
    Store an abundances object in the output buffers.
    !!}
    use :: Kind_Numbers                      , only : kind_int8
    use :: Multi_Counters                    , only : multiCounter
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger , outputPropertyDouble
    implicit none
    class           (abundances           ), intent(in   )               :: self
    double precision                       , intent(in   )               :: time
    integer                                , intent(inout)               :: doubleBufferCount , doubleProperty , &
         &                                                                  integerBufferCount, integerProperty
    type            (outputPropertyInteger), intent(inout), dimension(:) :: integerProperties
    type            (outputPropertyDouble ), intent(inout), dimension(:) :: doubleProperties
    type            (multiCounter         ), intent(in   )               :: outputInstance
    integer                                                              :: i
    !$GLC attributes unused :: time, integerBufferCount, integerProperty, integerProperties, outputInstance

    doubleProperty=doubleProperty+1
    doubleProperties(doubleProperty)%scalar(doubleBufferCount)=self%metallicityValue
    if (elementsCount > 0) then
       do i=1,elementsCount
          doubleProperties(doubleProperty+i)%scalar(doubleBufferCount)=self%elementalValue(i)
       end do
       doubleProperty=doubleProperty+elementsCount
    end if
    return
  end subroutine Abundances_Output

  subroutine Abundances_Post_Output(self,time)
    !!{
    Perform post-output processing of abundances objects.
    !!}
    implicit none
    class           (abundances), intent(in   ) :: self
    double precision            , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    return
  end subroutine Abundances_Post_Output

  subroutine Abundances_Output_Count(self,integerPropertyCount,doublePropertyCount,time)
    !!{
    Increment the output count to account for an abundances object.
    !!}
    implicit none
    class           (abundances), intent(in   ) :: self
    integer                     , intent(inout) :: doublePropertyCount, integerPropertyCount
    double precision            , intent(in   ) :: time
    !$GLC attributes unused :: self, integerPropertyCount, time

    doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Abundances_Output_Count

  subroutine Abundances_Output_Names(self,integerProperty,integerProperties,doubleProperty,doubleProperties,time,prefix,comment,unitsInSI)
    !!{
    Assign names to output buffers for an abundances object.
    !!}
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
    implicit none
    class           (abundances           )              , intent(in   ) :: self
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleProperty    , integerProperty
    type            (outputPropertyInteger), dimension(:), intent(inout) :: integerProperties
    type            (outputPropertyDouble ), dimension(:), intent(inout) :: doubleProperties
    character       (len=*                )              , intent(in   ) :: comment           , prefix
    double precision                                     , intent(in   ) :: unitsInSI
    integer                                                              :: iElement
    !$GLC attributes unused :: self, time, integerProperty, integerProperties

    doubleProperty=doubleProperty+1
    doubleProperties(doubleProperty)%name     =trim(prefix )//  'Metals'
    doubleProperties(doubleProperty)%comment  =trim(comment)//' [Metals]'
    doubleProperties(doubleProperty)%unitsInSI=unitsInSI
    if (elementsCount > 0) then
       do iElement=1,elementsCount
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%name     =trim(prefix )//      trim(elementsToTrack(iElement))
          doubleProperties(doubleProperty)%comment  =trim(comment)//' ['//trim(elementsToTrack(iElement))//']'
          doubleProperties(doubleProperty)%unitsInSI=unitsInSI
       end do
    end if
    return
  end subroutine Abundances_Output_Names

  function Abundances_Non_Static_Size_Of(self)
    !!{
    Return the size of any non-static components of the object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t  )                :: Abundances_Non_Static_Size_Of
    class  (abundances), intent(in   ) :: self

    if (allocated(self%elementalValue)) then
       Abundances_Non_Static_Size_Of=sizeof(self%elementalValue)
    else
       Abundances_Non_Static_Size_Of=0_c_size_t
    end if
    return
  end function Abundances_Non_Static_Size_Of

end module Abundances_Structure
