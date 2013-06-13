!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which defines the abundances structure used for describing elemental abundances in \glc.

module Abundances_Structure
  !% Defines the abundances structure used for describing elemental abundances in \glc.
  use ISO_Varying_String
  use Numerical_Constants_Astronomical
  implicit none
  private
  public :: abundances, Abundances_Names, Abundances_Atomic_Index, Abundances_Property_Count, Abundances_Get_Metallicity&
       &, Abundances_Mass_To_Mass_Fraction, operator(*), max
  
  ! Interface to multiplication operators with abundances objects as their second argument.
  interface operator(*)
     module procedure Abundances_Multiply_Switched
  end interface operator(*)

  ! Interface to max() function for abundances objects.
  interface max
     module procedure Abundances_Max
  end interface max

  type abundances
     !% The abundances structure used for describing elemental abundances in \glc.
     private
     double precision                            :: metallicityValue
     double precision, allocatable, dimension(:) :: elementalValue
   contains
     !@ <objectMethods>
     !@   <object>abundances</object>
     !@   <objectMethod>
     !@     <method>multiply</method>
     !@     <type>\textcolor{red}{\textless type(abundances)\textgreater}</type>
     !@     <arguments>\doublezero\ multiplier\argin</arguments>
     !@     <description>Multiply an abundance by a scalar.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>divide</method>
     !@     <type>\textcolor{red}{\textless type(abundances)\textgreater}</type>
     !@     <arguments>\doublezero\ divisor\argin</arguments>
     !@     <description>Divide an abundance by a scalar.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>add</method>
     !@     <type>\textcolor{red}{\textless type(abundances)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless type(abundances)\textgreater} abundances2\argin</arguments>
     !@     <description>Add two abundances.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>subtract</method>
     !@     <type>\textcolor{red}{\textless type(abundances)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless type(abundances)\textgreater} abundances2\argin</arguments>
     !@     <description>Subtract one abundance from another.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>metallicity</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\enumMetallicityScale\ [metallicityType]\argin</arguments>
     !@     <description>Returns the metallicity.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>metallicitySet</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ metallicity\argin, \enumMetallicityScale\ [metallicityType]\argin, \enumAdjustElements\ [adjustElements]\argin, \intzero\ [abundanceIndex]\argin</arguments>
     !@     <description>Sets the metallicity to {\tt metallicity}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>massToMassFraction</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ mass\argin</arguments>
     !@     <description>Converts abundance masses to mass fractions by dividing by the given {\tt mass} while ensuring that fractions are in the range 0--1.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>increment</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(abundances)\textgreater} addAbundances\argin</arguments>
     !@     <description>Increment an abundances object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serializeCount</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return a count of the number of properties in a serialized abundances object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serialize</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ historyArray\argout</arguments>
     !@     <description>Serialize an abundances object to an array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>deserialize</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ historyArray\argin</arguments>
     !@     <description>Deserialize an abundances object from an array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>hydrogenNumberFraction</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Returns the hydrogen fraction by number.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>hydrogenMassFraction</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Returns the hydrogen fraction by mass.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>heliumMassFraction</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Returns the helium fraction by mass.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>heliumNumberFraction</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Returns the helium fraction by number.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isZero</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if an abundances object is zero.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Destroy an abundances object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Reset an abundances object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>builder</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless *type(node)\textgreater} abundancesDefinition\argin</arguments>
     !@     <description>Build an abundances object from a provided XML description.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dump</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Dump an abundances object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <type>\void</type>
     !@     <method>dumpRaw</method>
     !@     <arguments>\intzero\ fileHandle\argin</arguments>
     !@     <description>Dump an abundances object to binary.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>readRaw</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ fileHandle\argin</arguments>
     !@     <description>Read an abundances object from binary.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setToUnity</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Set an abundances object to unity.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>output</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ integerProperty\arginout, \intzero\ integerBufferCount\arginout, \inttwo\ integerBuffer\arginout, \intzero doubleProperty\arginout, \intzero\ doubleBufferCount\arginout, \doubletwo\ doubleBuffer\arginout, \doublezero\ time\argin, \intzero\ instance\argin</arguments>
     !@     <description>Store an abundances object in the output buffers.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>outputCount</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ integerPropertyCount\arginout, \intzero\ doublePropertyCount\arginout, \doublezero\ time\argin, \intzero\ instance\argin</arguments>
     !@     <description>Specify the count of an abundances object for output.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>outputNames</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ integerProperty\arginout, \textcolor{red}{\textless char[*](:)\textgreater} integerPropertyNames\arginout, \textcolor{red}{\textless char[*](:)\textgreater} integerPropertyComments\arginout, \doubleone\ integerPropertyUnitsSI\arginout, \intzero\ doubleProperty\arginout, \textcolor{red}{\textless char[*](:)\textgreater} doublePropertyNames\arginout, \textcolor{red}{\textless char[*](:)\textgreater} doublePropertyComments\arginout, \doubleone\ doublePropertyUnitsSI\arginout, \doublezero\ time\argin, \intzero\ instance\argin</arguments>
     !@     <description>Specify the names of abundance object properties for output.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: add                    => Abundances_Add
     procedure                 :: subtract               => Abundances_Subtract
     procedure                 :: multiply               => Abundances_Multiply
     procedure                 :: divide                 => Abundances_Divide
     generic                   :: operator(+)            => add
     generic                   :: operator(-)            => subtract
     generic                   :: operator(*)            => multiply
     generic                   :: operator(/)            => divide
     procedure                 :: isZero                 => Abundances_Is_Zero
     procedure                 :: destroy                => Abundances_Destroy
     procedure                 :: reset                  => Abundances_Reset
     procedure                 :: builder                => Abundances_Builder
     procedure                 :: dump                   => Abundances_Dump
     procedure                 :: dumpRaw                => Abundances_Dump_Raw
     procedure                 :: readRaw                => Abundances_Read_Raw
     procedure                 :: setToUnity             => Abundances_Set_To_Unity
     procedure, nopass         :: serializeCount         => Abundances_Property_Count
     procedure                 :: serialize              => Abundances_Serialize
     procedure                 :: deserialize            => Abundances_Deserialize
     procedure                 :: increment              => Abundances_Increment
     procedure                 :: metallicity            => Abundances_Get_Metallicity
     procedure                 :: metallicitySet         => Abundances_Set_Metallicity
     procedure                 :: massToMassFraction     => Abundances_Mass_To_Mass_Fraction_Packed
     procedure                 :: hydrogenNumberFraction => Abundances_Hydrogen_Number_Fraction
     procedure                 :: hydrogenMassFraction   => Abundances_Hydrogen_Mass_Fraction
     procedure                 :: heliumMassFraction     => Abundances_Helium_Mass_Fraction
     procedure                 :: heliumNumberFraction   => Abundances_Helium_Number_Fraction
     procedure                 :: output      => Abundances_Output
     procedure                 :: outputCount => Abundances_Output_Count
     procedure                 :: outputNames => Abundances_Output_Names
  end type abundances

  ! Count of the number of elements being tracked.
  integer                                     :: elementsCount=0
  integer                                     :: propertyCount

  ! Names (two-letter labels) of elements to track.
  character(len=3), allocatable, dimension(:) :: elementsToTrack

  ! Indices of elements as used in the Atomic_Data module.
  integer,          allocatable, dimension(:) :: elementsIndices

  ! Type of metallicity/abundance measure required.
  !@ <enumeration>
  !@  <name>metallicityScale</name>
  !@  <description>Used to specify the metallicity scale when working with {\tt abundances} objects.</description>
  !@  <entry label="linearByMass"             />
  !@  <entry label="linearByNumber"           />
  !@  <entry label="logarithmicByMassSolar"   />
  !@  <entry label="logarithmicByNumberSolar" />
  !@  <entry label="linearByMassSolar"        />
  !@  <entry label="linearByNumberSolar"      />
  !@ </enumeration>
  integer,          parameter, public         :: linearByMass            =0
  integer,          parameter, public         :: linearByNumber          =1
  integer,          parameter, public         :: logarithmicByMassSolar  =2
  integer,          parameter, public         :: logarithmicByNumberSolar=3
  integer,          parameter, public         :: linearByMassSolar       =4
  integer,          parameter, public         :: linearByNumberSolar     =5
  
  ! Value used to indicate a zero metallicity on logarithmic scales.
  double precision, parameter, public         :: logMetallicityZero      =-99.0d0

  ! Flag indicating if this module has been initialized.
  logical                                     :: abundancesInitialized=.false.

  ! Labels used in determining how to update elemental abundances when metallicity is adjusted.
  !@ <enumeration>
  !@  <name>adjustElements</name>
  !@  <description>Used to specify how elements should be adjusted when the metallicity of an {\tt abundances} object is changed.</description>
  !@  <entry label="adjustElementsNone"   />
  !@  <entry label="adjustElementsReset"  />
  !@  <entry label="adjustElementsUpdate" />
  !@ </enumeration>
  integer,          parameter, public         :: adjustElementsNone  =0
  integer,          parameter, public         :: adjustElementsReset =1
  integer,          parameter, public         :: adjustElementsUpdate=2

  ! Unit and zero abundances objects.
  type(abundances),            public         :: zeroAbundances,unitAbundances

contains

  subroutine Abundances_Initialize
    !% Initialize the {\tt abundanceStructure} object module. Determines which abundances are to be tracked.
    use Input_Parameters
    use Memory_Management
    use Atomic_Data
    implicit none
    integer :: iElement

    ! Check if this module has been initialized already.    
    if (.not.abundancesInitialized) then
       !$omp critical (Abundances_Module_Initialize)
       if (.not.abundancesInitialized) then
          ! Determine how many elements we are required to track.
          elementsCount=Get_Input_Parameter_Array_Size('elementsToTrack')
          ! Number of properties to track is one greater, as we always track total metallicity.
          propertyCount=elementsCount+1
          ! If tracking elements, read names of which ones to track.
          if (elementsCount > 0) then
             call Alloc_Array(elementsToTrack,[elementsCount])
             call Alloc_Array(elementsIndices,[elementsCount])
             !@ <inputParameter>
             !@   <name>elementsToTrack</name>
             !@   <defaultValue></defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The names of the elements to be tracked.
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>1..*</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('elementsToTrack',elementsToTrack)
             ! Validate the input names by looking them up in the list of atomic names.
             do iElement=1,elementsCount
                elementsIndices(iElement)=Atom_Lookup(shortLabel=elementsToTrack(iElement))
             end do
          end if
          ! Create zero and unit abundances objects.
          call Alloc_Array(zeroAbundances%elementalValue,[elementsCount])
          call Alloc_Array(unitAbundances%elementalValue,[elementsCount])
          zeroAbundances%metallicityValue=0.0d0
          zeroAbundances%  elementalValue=0.0d0
          unitAbundances%metallicityValue=1.0d0
          unitAbundances%  elementalValue=1.0d0
          ! Flag that this module is now initialized.
          abundancesInitialized=.true.
       end if
       !$omp end critical (Abundances_Module_Initialize)
    end if
    return
  end subroutine Abundances_Initialize

  subroutine Abundances_Destroy(self)
    !% Destroy an abundances object.
    use Memory_Management
    implicit none
    class(abundances), intent(inout) :: self

    if (allocated(self%elementalValue)) call Dealloc_Array(self%elementalValue)
    return
  end subroutine Abundances_Destroy

  subroutine Abundances_Builder(self,abundancesDefinition)
    !% Build a {\tt abundances} object from the given XML {\tt abundancesDefinition}.
    use FoX_DOM
    use Galacticus_Error
    implicit none
    class(abundances), intent(inout)          :: self
    type (node      ), intent(in   ), pointer :: abundancesDefinition
    type (node      ),                pointer :: abundance
    type (nodeList  ),                pointer :: abundanceList
    integer                                   :: i

    ! Get the metallicity.
    abundanceList => getElementsByTagName(abundancesDefinition,'metals')
    if (getLength(abundanceList) >  1) call Galacticus_Error_Report('Abundances_Builder','multiple metallicity values specified')
    if (getLength(abundanceList) == 1) then
       abundance => item(abundanceList,0)
       call extractDataContent(abundance,self%metallicityValue)
    end if
    if (elementsCount > 0) then
       do i=1,elementsCount
          abundanceList => getElementsByTagName(abundancesDefinition,trim(elementsToTrack(i)))
          if (getLength(abundanceList) >  1) call Galacticus_Error_Report('Abundances_Builder','multiple '//trim(elementsToTrack(i))//' values specified')
          if (getLength(abundanceList) == 1) then
             abundance => item(abundanceList,0)
             call extractDataContent(abundance,self%elementalValue(i))
          end if
       end do
    end if
    return
  end subroutine Abundances_Builder

  subroutine Abundances_Dump(self)
    !% Reset an abundances object.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class    (abundances    ), intent(in   ) :: self
    integer                                  :: i
    character(len=12        )                :: label
    type     (varying_string)                :: message

    ! Ensure module is initialized.
    call Abundances_Initialize

    write (label,'(e12.6)') self%metallicityValue
    message='metallicity: '//label
    call Galacticus_Display_Message(message)
    if (elementsCount > 0) then
       do i=1,elementsCount
          write (label,'(e12.6)') self%elementalValue(i)
          message=elementsToTrack(i)//':          '//label
          call Galacticus_Display_Message(message)
       end do
    end if
    return
  end subroutine Abundances_Dump

  subroutine Abundances_Dump_Raw(self,fileHandle)
    !% Dump an abundances object to binary.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class    (abundances    ), intent(in   ) :: self
    integer                  , intent(in   ) :: fileHandle
    integer                                  :: i

    ! Ensure module is initialized.
    call Abundances_Initialize
    ! Dump the content.
    write (fileHandle) self%metallicityValue
    if (elementsCount > 0) write (fileHandle) self%elementalValue
    return
  end subroutine Abundances_Dump_Raw

  subroutine Abundances_Read_Raw(self,fileHandle)
    !% Read an abundances object from binary.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class    (abundances    ), intent(inout) :: self
    integer                  , intent(in   ) :: fileHandle
    integer                                  :: i

    ! Ensure module is initialized.
    call Abundances_Initialize
    ! Read the content.
    read (fileHandle) self%metallicityValue
    if (elementsCount > 0) read (fileHandle) self%elementalValue
    return
  end subroutine Abundances_Read_Raw

  subroutine Abundances_Reset(self)
    !% Reset an abundances object.
    implicit none
    class(abundances), intent(inout) :: self

    ! Ensure module is initialized.
    call Abundances_Initialize

    ! Zero all properties.
    call self%metallicitySet(0.0d0,adjustElements=adjustElementsReset)
    return
  end subroutine Abundances_Reset

  subroutine Abundances_Set_To_Unity(self)
    !% Set an abundances object to unity.
    implicit none
    class(abundances), intent(inout) :: self

    ! Ensure module is initialized.
    call Abundances_Initialize
    ! Ensure object is initialized.
    call Abundances_Allocate_Elemental_Values(self)
    ! Set values to unity.
    self                       %metallicityValue=1.0d0
    if (elementsCount > 0) self%  elementalValue=1.0d0
    return
  end subroutine Abundances_Set_To_Unity

  logical function Abundances_Is_Zero(self)
    !% Test whether an abundances object is zero.
    implicit none
    class(abundances), intent(in) :: self

    ! Ensure module is initialized.
    call Abundances_Initialize

    ! Detect if all abundances are zero.
    Abundances_Is_Zero=.true.
    if (self%metallicityValue /= 0.0d0) Abundances_Is_Zero=.false.
    if (elementsCount > 0 .and. allocated(self%elementalValue)) then
       if (any(self%elementalValue /= 0.0d0)) Abundances_Is_Zero=.false.
    end if
    return
  end function Abundances_Is_Zero

  function Abundances_Add(abundances1,abundances2)
    !% Add two abundances objects.
    implicit none
    type (abundances)                       :: Abundances_Add
    class(abundances), intent(in)           :: abundances1
    class(abundances), intent(in), optional :: abundances2

    ! Ensure module is initialized.
    call Abundances_Initialize
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
    !% Increment an abundances object.
    implicit none
    class(abundances), intent(inout) :: self
    class(abundances), intent(in   ) :: increment

    ! Ensure module is initialized.
    call Abundances_Initialize

    self                       %metallicityValue=self%metallicityValue+increment%metallicityValue
    if (elementsCount > 0) self%elementalValue  =self%elementalValue  +increment%elementalValue
    return
  end subroutine Abundances_Increment

  function Abundances_Subtract(abundances1,abundances2)
    !% Subtract two abundances objects.
    implicit none
    type (abundances)                       :: Abundances_Subtract
    class(abundances), intent(in)           :: abundances1
    class(abundances), intent(in), optional :: abundances2

    ! Ensure module is initialized.
    call Abundances_Initialize
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
    !% Multiply an abundances object by a scalar.
    implicit none
    type (abundances)             :: Abundances_Multiply
    class(abundances), intent(in) :: abundances1
    double precision , intent(in) :: multiplier

    ! Ensure module is initialized.
    call Abundances_Initialize
    Abundances_Multiply                       %metallicityValue=abundances1%metallicityValue*multiplier
    if (elementsCount > 0) Abundances_Multiply%  elementalValue=abundances1%  elementalValue*multiplier
    return
  end function Abundances_Multiply

  function Abundances_Multiply_Switched(multiplier,abundances1)
    !% Multiply a scalar by an abundances object.
    implicit none
    type(abundances)             :: Abundances_Multiply_Switched
    type(abundances), intent(in) :: abundances1
    double precision, intent(in) :: multiplier

    Abundances_Multiply_Switched=Abundances_Multiply(abundances1,multiplier)
    return
  end function Abundances_Multiply_Switched

  function Abundances_Max(abundances1,abundances2)
    !% Return an element-by-element {\tt max()} on two abundances objects.
    implicit none
    type(abundances)             :: Abundances_Max
    type(abundances), intent(in) :: abundances1,abundances2

    Abundances_Max                       %metallicityValue=max(abundances1%metallicityValue,abundances2%metallicityValue)
    if (elementsCount > 0) Abundances_Max%  elementalValue=max(abundances1%  elementalValue,abundances2%  elementalValue)
    return
  end function Abundances_Max

  function Abundances_Divide(abundances1,divisor)
    !% Divide an abundances object by a scalar.
    implicit none
    type (abundances)             :: Abundances_Divide
    class(abundances), intent(in) :: abundances1
    double precision , intent(in) :: divisor

    ! Ensure module is initialized.
    call Abundances_Initialize
    Abundances_Divide                       %metallicityValue=abundances1%metallicityValue/divisor
    if (elementsCount > 0) Abundances_Divide%  elementalValue=abundances1%  elementalValue/divisor
    return
  end function Abundances_Divide

  integer function Abundances_Property_Count()
    !% Return the number of properties required to track abundances. This is equal to the number of elements tracked, {\tt
    !% elementsCount}, plus one since we always track a total metallicity.
    implicit none

    ! Ensure module is initialized.
    call Abundances_Initialize

    Abundances_Property_Count=propertyCount
    return
  end function Abundances_Property_Count

  function Abundances_Names(index)
    !% Return a name for the specified entry in the abundances structure.
    use Galacticus_Error
    implicit none
    type(varying_string)             :: Abundances_Names
    integer,              intent(in) :: index

    ! Ensure module is initialized.
    call Abundances_Initialize

    select case (index)
    case (1)
       Abundances_Names="Metals"
    case (2:)
       if (index <= elementsCount+1) then
          Abundances_Names=trim(elementsToTrack(index-1))
       else
          call Galacticus_Error_Report('Abundances_Names','index out of range')
       end if
    case default
       call Galacticus_Error_Report('Abundances_Names','index out of range')
    end select
    return
  end function Abundances_Names

  integer function Abundances_Atomic_Index(index)
    !% Return the atomic index for the specified entry in the abundances structure.
    use Galacticus_Error
    implicit none
    integer, intent(in) :: index

    ! Ensure module is initialized.
    call Abundances_Initialize

    select case (index)
    case (1)
       Abundances_Atomic_Index=0 ! Total metallicity.
    case (2:)
       if (index <= elementsCount+1) then
          Abundances_Atomic_Index=elementsIndices(index-1)
       else
          call Galacticus_Error_Report('Abundances_Atomic_Index','index out of range')
       end if
    case default
       call Galacticus_Error_Report('Abundances_Atomic_Index','index out of range')
    end select
    return
  end function Abundances_Atomic_Index

  subroutine Abundances_Allocate_Elemental_Values(self)
    !% Ensure that the {\tt elementalValue} array in an {\tt abundances} is allocated.
    use Memory_Management
    implicit none
    type(abundances), intent(inout) :: self

    if (.not.allocated(self%elementalValue)) call Alloc_Array(self%elementalValue,[elementsCount])
    return
  end subroutine Abundances_Allocate_Elemental_Values

  subroutine Abundances_Deserialize(self,abundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    class(abundances), intent(inout)               :: self
    double precision,  intent(in   ), dimension(:) :: abundancesArray

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (self)
    type is (abundances)
    ! Extract metallicity from array.
    self%metallicityValue=abundancesArray(1)
    ! Ensure elemental values array exists.
    call Abundances_Allocate_Elemental_Values(self)
    ! Extract elemental values from array.
    self%elementalValue=abundancesArray(2:elementsCount+1)
    end select
    return
  end subroutine Abundances_Deserialize

  subroutine Abundances_Serialize(self,abundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    double precision,  intent(  out), dimension(:) :: abundancesArray(:)
    class(abundances), intent(in   )               :: self

    ! Ensure module is initialized.
    call Abundances_Initialize

    ! Place metallicity into array.
    abundancesArray(1)=self%metallicityValue
    ! Place elemental values into arrays.
    if (allocated(self%elementalValue)) then
       abundancesArray(2:elementsCount+1)=self%elementalValue
    else
       abundancesArray(2:elementsCount+1)=0.0d0
    end if
    return
  end subroutine Abundances_Serialize

  double precision function Abundances_Get_Metallicity(self,metallicityType)
    !% Return the metallicity of the {\tt self} structure.
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    implicit none
    class(abundances), intent(in)           :: self
    integer,           intent(in), optional :: metallicityType
    integer                                 :: metallicityTypeActual

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (self)
    type is (abundances)
       Abundances_Get_Metallicity=self%metallicityValue
       
       if (present(metallicityType)) then
          metallicityTypeActual=metallicityType
       else
          metallicityTypeActual=linearByMass
       end if
       
       select case (metallicityTypeActual)
       case (linearByMass)
          ! Do nothing, this is what we compute by default.
       case (logarithmicByMassSolar)
          ! Convert to a logarithmic metallicity by mass relative to Solar.
          if (Abundances_Get_Metallicity > 0.0d0) then
             Abundances_Get_Metallicity=log10(Abundances_Get_Metallicity/metallicitySolar)
          else
             Abundances_Get_Metallicity=logMetallicityZero
          end if
       case (linearByMassSolar)
          ! Convert to metallicity by mass relative to Solar.
          Abundances_Get_Metallicity=Abundances_Get_Metallicity/metallicitySolar
       case default
          call Galacticus_Error_Report('Abundances_Get_Metallicity','metallicity type not supported')
       end select
    end select

    return
  end function Abundances_Get_Metallicity

  subroutine Abundances_Set_Metallicity(self,metallicity,metallicityType,adjustElements,abundanceIndex)
    !% Set the metallicity of the {\tt self} structure to {\tt metallicity}.
    use Galacticus_Error
    use Atomic_Data
    implicit none
    class(abundances), intent(inout)           :: self
    double precision,  intent(in   )           :: metallicity
    integer,           intent(in   ), optional :: metallicityType,adjustElements,abundanceIndex
    integer                                    :: adjustElementsActual,iElement
    double precision                           :: metallicityPrevious

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (self)
    type is (abundances)
    ! Store the current metallicity.
    metallicityPrevious        =self%metallicityValue

    ! Determine how elements will be adjusted.
    if (present(adjustElements)) then
       adjustElementsActual=adjustElements
    else
       adjustElementsActual=adjustElementsNone
    end if

    ! Store the current metallicity if necessary.
    if (elementsCount > 0 .and. adjustElementsActual == adjustElementsUpdate) metallicityPrevious=self%metallicityValue

    ! Update the metallicity.
    self%metallicityValue=metallicity
    if (present(metallicityType)) then
       select case (metallicityType)
       case (linearByMass)
          ! Do nothing, this is how we store metallicity.
       case (linearByMassSolar)
          self%metallicityValue=         self%metallicityValue *metallicitySolar
       case (logarithmicByMassSolar)
          self%metallicityValue=(10.0d0**self%metallicityValue)*metallicitySolar
       case default
          call Galacticus_Error_Report('Abundances_Set_Metallicity','type not supported')
       end select
    end if

    ! Determine what we're requested to do with any other elements.
    if (elementsCount > 0) then
       select case (adjustElementsActual)
       case (adjustElementsNone)
          ! Do nothing to the elemental abundances in this case.
       case (adjustElementsReset)
          ! Ensure that we have an abundanceIndex specified.
          if (self%metallicityValue /= 0.0d0 .and. .not.present(abundanceIndex)) call Galacticus_Error_Report('Abundances_Set_Metallicity', &
               & 'an abundance pattern must be specified in order to reset elemental abundances')
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
       case (adjustElementsUpdate)
          ! Ensure that we have an abundanceIndex specified.
          if (.not.present(abundanceIndex)) call Galacticus_Error_Report('Abundances_Set_Metallicity', &
               & 'an abundance pattern must be specified in order to reset elemental abundances')
          ! Ensure elemental values array exists.
          call Abundances_Allocate_Elemental_Values(self)
          do iElement=1,elementsCount
             self%elementalValue(iElement)=self%elementalValue(iElement)+(self%metallicityValue&
                  &-metallicityPrevious)*Atomic_Abundance(abundanceIndex=abundanceIndex,atomIndex=elementsIndices(iElement)&
                  &,normalization=normalizationMetals)
          end do
        end select
    end if
    end select
    return
  end subroutine Abundances_Set_Metallicity

  subroutine Abundances_Mass_To_Mass_Fraction_Packed(self,mass)
    !% Convert abundance masses to mass fractions by dividing by {\tt mass} while ensuring that the fractions remain within the range 0--1.
    implicit none
    class(abundances), intent(inout) :: self
    double precision,  intent(in   ) :: mass

    ! Ensure module is initialized.
    call Abundances_Initialize

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
    !% Convert abundance masses to mass fractions by dividing by {\tt mass} while ensuring that the fractions remain within the range 0--1.
    implicit none
    double precision, intent(inout), dimension(:) :: self
    double precision, intent(in   )               :: mass
    
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
    !% Returns the mass fraction of hydrogen.
    implicit none
    class(abundances), intent(in) :: self
    double precision,  parameter  :: massFractionMinimum=0.7d0

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (self)
    type is (abundances)
    Abundances_Hydrogen_Mass_Fraction=max((self%metallicityValue/metallicitySolar)*(hydrogenByMassSolar&
         &-hydrogenByMassPrimordial)+hydrogenByMassPrimordial,massFractionMinimum)
    end select
    return
  end function Abundances_Hydrogen_Mass_Fraction

  double precision function Abundances_Helium_Mass_Fraction(self)
    !% Returns the mass fraction of helium.
    implicit none
    class(abundances), intent(in) :: self

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (self)
    type is (abundances)
    Abundances_Helium_Mass_Fraction=(self%metallicityValue/metallicitySolar)*(heliumByMassSolar-heliumByMassPrimordial)&
         &+heliumByMassPrimordial
    end select
    return
  end function Abundances_Helium_Mass_Fraction

  double precision function Abundances_Hydrogen_Number_Fraction(self)
    !% Returns the number fraction of hydrogen.
    implicit none
    class(abundances), intent(in) :: self
    double precision              :: numberHydrogen,numberHelium

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (self)
    type is (abundances)
       numberHydrogen=Abundances_Hydrogen_Mass_Fraction(self)/atomicMassHydrogen
       numberHelium  =Abundances_Helium_Mass_Fraction  (self)/atomicMassHelium
       Abundances_Hydrogen_Number_Fraction=numberHydrogen/(numberHydrogen+numberHelium)
    end select
    return
  end function Abundances_Hydrogen_Number_Fraction

  double precision function Abundances_Helium_Number_Fraction(self)
    !% Returns the mass fraction of helium.
    implicit none
    class(abundances), intent(in) :: self
    double precision              :: numberHydrogen,numberHelium

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (self)
    type is (abundances)
       numberHydrogen=Abundances_Hydrogen_Mass_Fraction(self)/atomicMassHydrogen
       numberHelium  =Abundances_Helium_Mass_Fraction  (self)/atomicMassHelium
       Abundances_Helium_Number_Fraction=numberHelium/(numberHydrogen+numberHelium)
    end select
    return
  end function Abundances_Helium_Number_Fraction
  
  subroutine Abundances_Output(self,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount&
       &,doubleBuffer,time)
    !% Store an abundances object in the output buffers.
    implicit none
    class           (abundances    )                 , intent(in   ) :: self
    double precision                                 , intent(in   ) :: time
    integer                                          , intent(inout) :: integerProperty,integerBufferCount&
    &,doubleProperty,doubleBufferCount
    integer         (kind=kind_int8), dimension(:,:), intent(inout) :: integerBuffer
    double precision                , dimension(:,:), intent(inout) :: doubleBuffer

    doubleProperty=doubleProperty+1
    doubleBuffer(doubleBufferCount,doubleProperty)=self%metallicityValue
    if (elementsCount > 0) then
       doubleBuffer(doubleBufferCount,doubleProperty+1:doubleProperty+elementsCount)=self%elementalValue(:)
       doubleProperty=doubleProperty+elementsCount
    end if
    return
  end subroutine Abundances_Output

  subroutine Abundances_Output_Count(self,integerPropertyCount,doublePropertyCount,time)
    !% Increment the output count to account for an abundances object.
    implicit none
    class           (abundances), intent(in   ) :: self
    integer                     , intent(inout) :: integerPropertyCount,doublePropertyCount
    double precision            , intent(in   ) :: time                                    

    doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Abundances_Output_Count
  
  subroutine Abundances_Output_Names(self,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI&
       &,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,prefix,comment,unitsInSI)
    !% Assign names to output buffers for an abundances object.
    implicit none
    class           (abundances       )              , intent(in   ) :: self
    double precision                                 , intent(in   ) :: time
    integer                                          , intent(inout) :: integerProperty,doubleProperty
    character       (len=*            ), dimension(:), intent(inout) :: integerPropertyNames,integerPropertyComments,doublePropertyNames,doublePropertyComments
    double precision                   , dimension(:), intent(inout) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    character       (len=*            )              , intent(in   ) :: prefix,comment
    double precision                                 , intent(in   ) :: unitsInSI
    integer                                                          :: iElement

    doubleProperty=doubleProperty+1
    doublePropertyNames   (doubleProperty)=trim(prefix)//'Metals'
    doublePropertyComments(doubleProperty)=trim(comment)//' [Metals]'
    doublePropertyUnitsSI (doubleProperty)=unitsInSI
    if (elementsCount > 0) then
       do iElement=1,elementsCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)=trim(prefix)//trim(elementsToTrack(iElement))
          doublePropertyComments(doubleProperty)=trim(comment)//' ['//trim(elementsToTrack(iElement))//']'
          doublePropertyUnitsSI (doubleProperty)=unitsInSI
       end do
    end if
    return
  end subroutine Abundances_Output_Names
  
end module Abundances_Structure
