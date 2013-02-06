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

!% Contains a module which defines the structure used for describing chemical abundances in \glc.

module Chemical_Abundances_Structure
  !% Defines the structure used for describing chemical abundances in \glc.
  use ISO_Varying_String
  implicit none
  private
  public :: chemicalAbundances, Chemicals_Names, Chemicals_Index, Chemicals_Property_Count
  
  type chemicalAbundances
     !% The structure used for describing chemical abundances in \glc.
     private
     double precision, allocatable, dimension(:) :: chemicalValue
   contains
     ! Operators.
     procedure                 :: add                    => Chemical_Abundances_Add
     procedure                 :: subtract               => Chemical_Abundances_Subtract
     procedure                 :: multiply               => Chemical_Abundances_Multiply
     procedure                 :: divide                 => Chemical_Abundances_Divide
     generic                   :: operator(+)            => add
     generic                   :: operator(-)            => subtract
     generic                   :: operator(*)            => multiply
     generic                   :: operator(/)            => divide
     ! Serialization methods.
     !@ <objectMethods>
     !@   <object>chemicalAbundances</object>
     !@   <objectMethod>
     !@     <method>serializeCount</method>
     !@     <description>Return a count of the number of properties in a serialized chemical abundances object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serialize</method>
     !@     <description>Serialize a chemical abundances object to an array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>deserialize</method>
     !@     <description>Deserialize a chemical abundances object from an array.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure, nopass         :: serializeCount         => Chemicals_Property_Count
     procedure                 :: serialize              => Chemical_Abundances_Serialize
     procedure                 :: deserialize            => Chemical_Abundances_Deserialize
     !@ <objectMethod>
     !@   <object>chemicalAbundances</object>
     !@   <method>increment</method>
     !@   <description>Increment a chemical abundances object.</description>
     !@ </objectMethod>
     procedure                 :: increment              => Chemical_Abundances_Increment
     ! Abundance methods.
     !@ <objectMethods>
     !@   <object>chemicalAbundances</object>
     !@   <objectMethod>
     !@     <method>abundance</method>
     !@     <description>Returns the abundance of a chemical given its index.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>abundanceSet</method>
     !@     <description>Sets the abundance of a chemical given its index.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <description>Resets abundances to zero.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setToUnity</method>
     !@     <description>Set abundances to unity.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <description>Destroys a chemical abundances object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>numberToMass</method>
     !@     <description>Converts from abundances by number to abundances by mass.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>massToNumber</method>
     !@     <description>Converts from abundances by mass to abundances by number.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>enforcePositive</method>
     !@     <description>Enforces all chemical values to be positive.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>builder</method>
     !@     <description>Build a chemical abundances object from an XML definition.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dump</method>
     !@     <description>Dump a chemical abundances object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dumpRaw</method>
     !@     <description>Dump a chemical abundances object in binary.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>readRaw</method>
     !@     <description>Read a chemical abundances object in binary.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: abundance         => Chemicals_Abundances
     procedure                 :: abundanceSet      => Chemicals_Abundances_Set
     procedure                 :: reset             => Chemicals_Abundances_Reset
     procedure                 :: setToUnity        => Chemicals_Abundances_Set_To_Unity
     procedure                 :: destroy           => Chemicals_Abundances_Destroy
     procedure                 :: numberToMass      => Chemicals_Number_To_Mass
     procedure                 :: massToNumber      => Chemicals_Mass_To_Number
     procedure                 :: enforcePositive   => Chemicals_Enforce_Positive
     procedure                 :: builder           => Chemicals_Builder
     procedure                 :: dump              => Chemicals_Dump
     procedure                 :: dumpRaw           => Chemicals_Dump_Raw
     procedure                 :: readRaw           => Chemicals_Read_Raw
  end type chemicalAbundances

  ! Count of the number of elements being tracked.
  integer                                         :: chemicalsCount=0
  integer                                         :: propertyCount

  ! Names of chemicals to track.
  type(varying_string), allocatable, dimension(:) :: chemicalsToTrack
  integer                                         :: chemicalNameLengthMaximum

  ! Indices of chemicals as used in the Chemical_Structures module.
  integer,              allocatable, dimension(:) :: chemicalsIndices

  ! Net charge and mass (in atomic units) of chemicals.
  double precision,     allocatable, dimension(:) :: chemicalsCharges,chemicalsMasses

  ! Flag indicating if this module has been initialized.
  logical                                         :: chemicalAbundancesInitialized=.false.

  ! Unit and zero chemical abundances objects.
  type(chemicalabundances),           public      :: zeroChemicals,unitChemicals

contains

  subroutine Chemical_Abundances_Initialize
    !% Initialize the {\tt chemicalAbundanceStructure} object module. Determines which chemicals are to be tracked.
    use Input_Parameters
    use Memory_Management
    use Chemical_Structures
    implicit none
    integer                 :: iChemical
    type(chemicalStructure) :: thisChemical

    ! Check if this module has been initialized already.    
    if (.not.chemicalAbundancesInitialized) then
       !$omp critical (Chemical_Abundances_Module_Initialize)
       if (.not.chemicalAbundancesInitialized) then
          
          ! Determine how many elements we are required to track.
          chemicalsCount=Get_Input_Parameter_Array_Size('chemicalsToTrack')
          ! Number of properties to track is the same as the number of chemicals.
          propertyCount=chemicalsCount
          ! If tracking chemicals, read names of which ones to track.
          if (chemicalsCount > 0) then
             allocate(chemicalsToTrack(chemicalsCount))
             call Alloc_Array(chemicalsIndices,[chemicalsCount])
             call Alloc_Array(chemicalsCharges,[chemicalsCount])
             call Alloc_Array(chemicalsMasses ,[chemicalsCount])
             !@ <inputParameter>
             !@   <name>chemicalsToTrack</name>
             !@   <defaultValue></defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The names of the chemicals to be tracked.
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>1..*</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('chemicalsToTrack',chemicalsToTrack)
             ! Validate the input names by looking them up in the list of chemical names.
             chemicalNameLengthMaximum=0
             do iChemical=1,chemicalsCount
                chemicalsIndices(iChemical)=Chemical_Database_Get_Index(char(chemicalsToTrack(iChemical)))
                call thisChemical%retrieve(char(chemicalsToTrack(iChemical)))
                chemicalsCharges(iChemical)=dble(thisChemical%charge())
                chemicalsMasses (iChemical)=     thisChemical%mass  ()
                if (len(chemicalsToTrack(iChemical)) > chemicalNameLengthMaximum) chemicalNameLengthMaximum&
                     &=len(chemicalsToTrack(iChemical))
             end do
          end if
          ! Create zero and unit chemical abundances objects.
          call Alloc_Array(zeroChemicals%chemicalValue,[propertyCount])
          call Alloc_Array(unitChemicals%chemicalValue,[propertyCount])
          zeroChemicals%chemicalValue=0.0d0
          unitChemicals%chemicalValue=1.0d0
          ! Flag that this module is now initialized.
          chemicalAbundancesInitialized=.true.
       end if
       !$omp end critical (Chemical_Abundances_Module_Initialize)
    end if
    return
  end subroutine Chemical_Abundances_Initialize

  integer function Chemicals_Property_Count()
    !% Return the number of properties required to track chemicals. This is equal to the number of chemicals tracked, {\tt
    !% chemicalsCount}.
    implicit none

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize

    Chemicals_Property_Count=propertyCount
    return
  end function Chemicals_Property_Count

  function Chemicals_Names(index)
    !% Return a name for the specified entry in the chemicals structure.
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(varying_string)             :: Chemicals_Names
    integer,              intent(in) :: index

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize

    if (index >= 1 .and. index <= chemicalsCount) then
       Chemicals_Names=trim(chemicalsToTrack(index))
    else
       call Galacticus_Error_Report('Chemicals_Names','index out of range')
    end if
    return
  end function Chemicals_Names

  integer function Chemicals_Index(chemicalName)
    !% Returns the index of a chemical in the chemical abundances structure given the {\tt chemicalName}.
    implicit none
    character(len=*), intent(in) :: chemicalName
    integer                      :: iChemical

    Chemicals_Index=-1 ! Indicates chemical not found.
    do iChemical=1,chemicalsCount
       if (chemicalsToTrack(iChemical) == trim(chemicalName)) then
          Chemicals_Index=iChemical
          return
       end if
    end do
    return
  end function Chemicals_Index
  
  subroutine Chemical_Abundances_Increment(self,increment)
    !% Increment an abundances object.
    implicit none
    class(chemicalAbundances), intent(inout) :: self
    class(chemicalAbundances), intent(in   ) :: increment

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize
    self%chemicalValue=self%chemicalValue+increment%chemicalValue
    return
  end subroutine Chemical_Abundances_Increment

  function Chemical_Abundances_Add(abundances1,abundances2)
    !% Add two abundances objects.
    implicit none
    type (chemicalAbundances)                       :: Chemical_Abundances_Add
    class(chemicalAbundances), intent(in)           :: abundances1
    class(chemicalAbundances), intent(in), optional :: abundances2

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize
    if (chemicalsCount == 0) then
       Chemical_Abundances_Add=zeroChemicals
    else
       if (present(abundances2)) then
          Chemical_Abundances_Add%chemicalValue=abundances1%chemicalValue+abundances2%chemicalValue
       else
          Chemical_Abundances_Add%chemicalValue=abundances1%chemicalValue
       end if
    end if
    return
  end function Chemical_Abundances_Add

  function Chemical_Abundances_Subtract(abundances1,abundances2)
    !% Subtract two abundances objects.
    implicit none
    type (chemicalAbundances)                       :: Chemical_Abundances_Subtract
    class(chemicalAbundances), intent(in)           :: abundances1
    class(chemicalAbundances), intent(in), optional :: abundances2

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize
    if (chemicalsCount == 0) then
       Chemical_Abundances_Subtract=zeroChemicals
    else
       if (present(abundances2)) then
          Chemical_Abundances_Subtract%chemicalValue= abundances1%chemicalValue-abundances2%chemicalValue
       else
          Chemical_Abundances_Subtract%chemicalValue=-abundances1%chemicalValue
       end if
    end if
    return
  end function Chemical_Abundances_Subtract

  function Chemical_Abundances_Multiply(abundances1,multiplier)
    !% Multiply a chemical abundances object by a scalar.
    implicit none
    type (chemicalAbundances)             :: Chemical_Abundances_Multiply
    class(chemicalAbundances), intent(in) :: abundances1
    double precision         , intent(in) :: multiplier

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize
    if (chemicalsCount == 0) then
       Chemical_Abundances_Multiply=zeroChemicals
    else
       Chemical_Abundances_Multiply%chemicalValue=abundances1%chemicalValue*multiplier
    end if
    return
  end function Chemical_Abundances_Multiply

  function Chemical_Abundances_Divide(abundances1,divisor)
    !% Divide a chemical abundances object by a scalar.
    implicit none
    type (chemicalAbundances)             :: Chemical_Abundances_Divide
    class(chemicalAbundances), intent(in) :: abundances1
    double precision         , intent(in) :: divisor

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize
    if (chemicalsCount == 0) then
       Chemical_Abundances_Divide=zeroChemicals
    else
       Chemical_Abundances_Divide%chemicalValue=abundances1%chemicalValue/divisor
    end if
    return
  end function Chemical_Abundances_Divide

  double precision function Chemicals_Abundances(chemicals,moleculeIndex)
    !% Returns the abundance of a molecule in the chemical abundances structure given the {\tt moleculeIndex}.
    implicit none
    class(chemicalAbundances), intent(in) :: chemicals
    integer,                             intent(in) :: moleculeIndex

    Chemicals_Abundances=chemicals%chemicalValue(moleculeIndex)
    return
  end function Chemicals_Abundances

  subroutine Chemicals_Number_To_Mass(chemicals,chemicalsByMass)
    !% Multiply all chemical species by their mass in units of the atomic mass. This converts abundances by number into abundances by mass.
    implicit none
    class(chemicalAbundances), intent(in)    :: chemicals
    type(chemicalAbundances),  intent(inout) :: chemicalsByMass

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicalsByMass)
    
    ! If the input chemical abundances structure is uninitialized, just return zero abundances.
    if (.not.allocated(chemicals%chemicalValue)) then
       call Chemicals_Abundances_Reset(chemicalsByMass)
       return
    end if

    ! Scale by the masses of the chemicals.
    chemicalsByMass%chemicalValue=chemicals%chemicalValue*chemicalsMasses
    return
  end subroutine Chemicals_Number_To_Mass

  subroutine Chemicals_Mass_To_Number(chemicals,chemicalsByNumber)
    !% Divide all chemical species by their mass in units of the atomic mass. This converts abundances by mass into abundances by number.
    implicit none
    class(chemicalAbundances), intent(in)    :: chemicals
    type(chemicalAbundances),  intent(inout) :: chemicalsByNumber

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicalsByNumber)
    
    ! If the input chemical abundances structure is uninitialized, just return zero abundances.
    if (.not.allocated(chemicals%chemicalValue)) then
       call Chemicals_Abundances_Reset(chemicalsByNumber)
       return
    end if

    ! Scale by the masses of the chemicals.
    chemicalsByNumber%chemicalValue=chemicals%chemicalValue/chemicalsMasses
    return
  end subroutine Chemicals_Mass_To_Number

  subroutine Chemicals_Enforce_Positive(chemicals)
    !% Force all chemical values to be positive, by truncating negative values to zero.
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    if (allocated(chemicals%chemicalValue)) then
       where (chemicals%chemicalValue < 0.0d0)
          chemicals%chemicalValue=0.0d0
       end where
    end if
    return
  end subroutine Chemicals_Enforce_Positive

  subroutine Chemicals_Builder(self,chemicalsDefinition)
    !% Build a {\tt chemicalAbundances} object from the given XML {\tt chemicalsDefinition}.
    use FoX_DOM
    use Galacticus_Error
    implicit none
    class(chemicalAbundances), intent(inout) :: self
    type (node              ), pointer       :: chemicalsDefinition

    call Galacticus_Error_Report('Chemicals_Builder','building of chemicalAbundances objects is not yet supported')
    return
  end subroutine Chemicals_Builder

  subroutine Chemicals_Dump(chemicals)
    !% Dump all chemical values.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class(chemicalAbundances), intent(in   ) :: chemicals
    integer                                  :: i
    character(len=12        )                :: label
    type     (varying_string)                :: message

    if (allocated(chemicals%chemicalValue)) then
       do i=1,chemicalsCount
          write (label,'(e12.6)') chemicals%chemicalValue(i)
          message=chemicalsToTrack(i)//': '//repeat(" ",chemicalNameLengthMaximum-len(chemicalsToTrack(i)))//label
          call Galacticus_Display_Message(message)
       end do
    end if
    return
  end subroutine Chemicals_Dump

  subroutine Chemicals_Dump_Raw(chemicals,fileHandle)
    !% Dump all chemical values in binary.
    implicit none
    class(chemicalAbundances), intent(in   ) :: chemicals
    integer                  , intent(in   ) :: fileHandle

    write (fileHandle) allocated(chemicals%chemicalValue)
    if (allocated(chemicals%chemicalValue)) write (fileHandle),chemicals%chemicalValue
    return
  end subroutine Chemicals_Dump_Raw

  subroutine Chemicals_Read_Raw(chemicals,fileHandle)
    !% Read all chemical values in binary.
    use Memory_Management
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals
    integer                  , intent(in   ) :: fileHandle
    logical                                  :: isAllocated

    read (fileHandle) isAllocated
    if (isAllocated) then
       call Alloc_Array(chemicals%chemicalValue,[chemicalsCount])
       read (fileHandle),chemicals%chemicalValue
    end if
    return
  end subroutine Chemicals_Read_Raw

  subroutine Chemicals_Abundances_Set(chemicals,moleculeIndex,abundance)
    !% Sets the abundance of a molecule in the chemical abundances structure given the {\tt moleculeIndex}.
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals
    integer,                   intent(in   ) :: moleculeIndex
    double precision,          intent(in   ) :: abundance
    
    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicals)
    
    select type (chemicals)
    type is (chemicalAbundances)
       chemicals%chemicalValue(moleculeIndex)=abundance
    end select
    return
  end subroutine Chemicals_Abundances_Set

  subroutine Chemicals_Abundances_Reset(chemicals)
    !% Resets all chemical abundances to zero.
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicals)

    ! Reset to zero.
    select type (chemicals)
    type is (chemicalAbundances)
       chemicals%chemicalValue=0.0d0
    end select
    return
  end subroutine Chemicals_Abundances_Reset

  subroutine Chemicals_Abundances_Set_To_Unity(chemicals)
    !% Resets all chemical abundances to unity.
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicals)

    ! Reset to zero.
    select type (chemicals)
    type is (chemicalAbundances)
       chemicals%chemicalValue=1.0d0
    end select
    return
  end subroutine Chemicals_Abundances_Set_To_Unity

  subroutine Chemicals_Abundances_Destroy(chemicals)
    !% Destroy a chemical abundances object.
    use Memory_Management
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    if (allocated(chemicals%chemicalValue)) call Dealloc_Array(chemicals%chemicalValue)
    return
  end subroutine Chemicals_Abundances_Destroy

  subroutine Chemical_Abundances_Allocate_Values(chemicals)
    !% Ensure that the {\tt chemicalValue} array in an {\tt chemicalsStructure} is allocated.
    use Memory_Management
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    select type (chemicals)
    type is (chemicalAbundances)
       if (.not.allocated(chemicals%chemicalValue)) then
          allocate(chemicals%chemicalValue(chemicalsCount))
          call Memory_Usage_Record(sizeof(chemicals%chemicalValue))
       end if
    end select
    return
  end subroutine Chemical_Abundances_Allocate_Values

  subroutine Chemical_Abundances_Deserialize(chemicals,chemicalAbundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    class(chemicalAbundances), intent(inout)            :: chemicals
    double precision,                   intent(in), dimension(:) :: chemicalAbundancesArray

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicals)
    ! Extract chemical values from array.
    select type (chemicals)
    type is (chemicalAbundances)
       chemicals%chemicalValue=chemicalAbundancesArray
    end select
    return
  end subroutine Chemical_Abundances_Deserialize

  subroutine Chemical_Abundances_Serialize(chemicals,chemicalAbundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    double precision,                   intent(out), dimension(:) :: chemicalAbundancesArray(:)
    class(chemicalAbundances), intent(in)                :: chemicals

    ! Ensure module is initialized.
    call Chemical_Abundances_Initialize

    ! Place elemental values into arrays.
    if (allocated(chemicals%chemicalValue)) then
       chemicalAbundancesArray=chemicals%chemicalValue
    else
       chemicalAbundancesArray=0.0d0
    end if
    return
  end subroutine Chemical_Abundances_Serialize

end module Chemical_Abundances_Structure
