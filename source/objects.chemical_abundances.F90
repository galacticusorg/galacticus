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
Contains a module which defines the structure used for describing chemical abundances in \glc.
!!}

module Chemical_Abundances_Structure
  !!{
  Defines the structure used for describing chemical abundances in \glc.
  !!}
  use :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: chemicalAbundances      , Chemical_Abundances_Initialize, Chemicals_Names, Chemicals_Index, &
       &    Chemicals_Property_Count, operator(*)                   , operator(>)

  ! Interface to multiplication operators with chemical abundances objects as their second argument.
  interface operator(*)
     module procedure Chemical_Abundances_Multiply_Switched
  end interface operator(*)

  ! Interface to greater than operators.
  interface operator(>)
     module procedure Chemical_Abundances_Greater_Than
  end interface operator(>)

  type :: chemicalAbundances
     !!{
     The structure used for describing chemical abundances in \glc.
     !!}
     private
     double precision, allocatable, dimension(:) :: chemicalValue
   contains
     !![
     <methods>
       <method description="Multiply a chemical abundance by a scalar." method="operator(*)" />
       <method description="Multiply (in-place) a chemical abundance by a scalar." method="scale" />
       <method description="Divide a chemical abundance by a scalar." method="operator(/)" />
       <method description="Add two chemical abundances." method="operator(+)" />
       <method description="Subtract one chemical abundance from another." method="operator(-)" />
       <method description="Return a count of the number of properties in a serialized chemical abundances object." method="serializeCount" />
       <method description="Serialize a chemical abundances object to an array." method="serialize" />
       <method description="Deserialize a chemical abundances object from an array." method="deserialize" />
       <method description="Increment a chemical abundances object." method="increment" />
       <method description="Returns the abundance of a chemical given its index." method="abundance" />
       <method description="Sets the abundance of a chemical given its index." method="abundanceSet" />
       <method description="Resets abundances to zero." method="reset" />
       <method description="Set abundances to unity." method="setToUnity" />
       <method description="Return true if a chemicals object is zero." method="isZero" />
       <method description="Destroys a chemical abundances object." method="destroy" />
       <method description="Converts from abundances by number to abundances by mass." method="numberToMass" />
       <method description="Converts from abundances by mass to abundances by number." method="massToNumber" />
       <method description="Enforces all chemical values to be positive." method="enforcePositive" />
       <method description="Returns the sum over all chemicals." method="sumOver" />
       <method description="Build a chemical abundances object from an XML definition." method="builder" />
       <method description="Dump a chemical abundances object." method="dump" />
       <method description="Dump a chemical abundances object in binary." method="dumpRaw" />
       <method description="Read a chemical abundances object in binary." method="readRaw" />
       <method description="Returns the size of any non-static components of the type." method="nonStaticSizeOf" />
       <method description="Store a chemical abundances object in the output buffers." method="output" />
       <method description="Store a chemical abundances object in the output buffers." method="postOutput" />
       <method description="Specify the count of a chemical abundances object for output." method="outputCount" />
       <method description="Specify the names of chemical abundance object properties for output." method="outputNames" />
     </methods>
     !!]
     procedure         ::                    Chemical_Abundances_Add
     procedure         ::                    Chemical_Abundances_Subtract
     procedure         ::                    Chemical_Abundances_Multiply
     procedure         ::                    Chemical_Abundances_Divide
     generic           :: operator(+)     => Chemical_Abundances_Add
     generic           :: operator(-)     => Chemical_Abundances_Subtract
     generic           :: operator(*)     => Chemical_Abundances_Multiply
     generic           :: operator(/)     => Chemical_Abundances_Divide
     procedure         :: nonStaticSizeOf => Chemicals_Non_Static_Size_Of
     procedure, nopass :: serializeCount  => Chemicals_Property_Count
     procedure         :: serialize       => Chemical_Abundances_Serialize
     procedure         :: deserialize     => Chemical_Abundances_Deserialize
     procedure         :: increment       => Chemical_Abundances_Increment
     procedure         :: scale           => Chemical_Abundances_Scale
     procedure         :: abundance       => Chemicals_Abundances
     procedure         :: abundanceSet    => Chemicals_Abundances_Set
     procedure         :: reset           => Chemicals_Abundances_Reset
     procedure         :: setToUnity      => Chemicals_Abundances_Set_To_Unity
     procedure         :: isZero          => Chemicals_Abundances_Is_Zero
     procedure         :: destroy         => Chemicals_Abundances_Destroy
     procedure         :: numberToMass    => Chemicals_Number_To_Mass
     procedure         :: massToNumber    => Chemicals_Mass_To_Number
     procedure         :: enforcePositive => Chemicals_Enforce_Positive
     procedure         :: sumOver         => Chemicals_Sum_Over
     procedure         :: builder         => Chemicals_Builder
     procedure         :: dump            => Chemicals_Dump
     procedure         :: dumpRaw         => Chemicals_Dump_Raw
     procedure         :: readRaw         => Chemicals_Read_Raw
     procedure         :: output          => Chemicals_Output
     procedure         :: postOutput      => Chemicals_Post_Output
     procedure         :: outputCount     => Chemicals_Output_Count
     procedure         :: outputNames     => Chemicals_Output_Names
  end type chemicalAbundances

  ! Count of the number of elements being tracked.
  integer                                                         :: chemicalsCount               =0
  integer                                                         :: propertyCount

  ! Names of chemicals to track.
  type            (varying_string    ), allocatable, dimension(:) :: chemicalsToTrack
  integer                                                         :: chemicalNameLengthMaximum

  ! Indices of chemicals as used in the Chemical_Structures module.
  integer                             , allocatable, dimension(:) :: chemicalsIndices

  ! Net charge and mass (in atomic units) of chemicals.
  double precision                    , allocatable, dimension(:) :: chemicalsCharges                     , chemicalsMasses

  ! Unit and zero chemical abundances objects.
  type            (chemicalabundances), public                    :: unitChemicalAbundances               , zeroChemicalAbundances

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Chemical_Abundances_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
   subroutine Chemical_Abundances_Initialize(parameters_)
    !!{
    Initialize the {\normalfont \ttfamily chemicalAbundanceStructure} object module. Determines which chemicals are to be tracked.
    !!}
    use :: Chemical_Structures, only : Chemical_Database_Get_Index, chemicalStructure
    use :: ISO_Varying_String , only : char                       , len
    use :: Input_Parameters   , only : inputParameters
    implicit none
    type   (inputParameters  ), intent(inout) :: parameters_
    integer                                   :: iChemical
    type   (chemicalStructure)                :: chemical

    ! Determine how many elements we are required to track.
    if (parameters_%isPresent('chemicalsToTrack')) then
       chemicalsCount=parameters_%count('chemicalsToTrack')
    else
       chemicalsCount=0
    end if
    ! Number of properties to track is the same as the number of chemicals.
    propertyCount=chemicalsCount
    ! If tracking chemicals, read names of which ones to track.
    if (chemicalsCount > 0) then
       if (allocated(chemicalsToTrack)) deallocate(chemicalsToTrack)
       if (allocated(chemicalsIndices)) deallocate(chemicalsIndices)
       if (allocated(chemicalsCharges)) deallocate(chemicalsCharges)
       if (allocated(chemicalsMasses )) deallocate(chemicalsMasses )
       allocate(chemicalsToTrack(chemicalsCount))
       allocate(chemicalsIndices(chemicalsCount))
       allocate(chemicalsCharges(chemicalsCount))
       allocate(chemicalsMasses (chemicalsCount))
       !![
       <inputParameter>
         <name>chemicalsToTrack</name>
         <description>The names of the chemicals to be tracked.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]
       ! Validate the input names by looking them up in the list of chemical names.
       chemicalNameLengthMaximum=0
       do iChemical=1,chemicalsCount
          chemicalsIndices(iChemical)=Chemical_Database_Get_Index(char(chemicalsToTrack(iChemical)))
          call chemical%retrieve(char(chemicalsToTrack(iChemical)))
          chemicalsCharges(iChemical)=dble(chemical%charge())
          chemicalsMasses (iChemical)=     chemical%mass  ()
          if (len(chemicalsToTrack(iChemical)) > chemicalNameLengthMaximum) chemicalNameLengthMaximum=len(chemicalsToTrack(iChemical))
       end do
    end if
    ! Create zero and unit chemical abundances objects.
    if (allocated(zeroChemicalAbundances%chemicalValue)) deallocate(zeroChemicalAbundances%chemicalValue)
    if (allocated(unitChemicalAbundances%chemicalValue)) deallocate(unitChemicalAbundances%chemicalValue)
    allocate(zeroChemicalAbundances%chemicalValue(propertyCount))
    allocate(unitChemicalAbundances%chemicalValue(propertyCount))
    zeroChemicalAbundances%chemicalValue=0.0d0
    unitChemicalAbundances%chemicalValue=1.0d0
    return
  end subroutine Chemical_Abundances_Initialize

  integer function Chemicals_Property_Count()
    !!{
    Return the number of properties required to track chemicals. This is equal to the number of chemicals tracked, {\normalfont \ttfamily
    chemicalsCount}.
    !!}
    implicit none

    Chemicals_Property_Count=propertyCount
    return
  end function Chemicals_Property_Count

  function Chemicals_Names(index)
    !!{
    Return a name for the specified entry in the chemicals structure.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : trim
    implicit none
    type   (varying_string)                :: Chemicals_Names
    integer                , intent(in   ) :: index

    if (index >= 1 .and. index <= chemicalsCount) then
       Chemicals_Names=trim(chemicalsToTrack(index))
    else
       call Error_Report('index out of range'//{introspection:location})
    end if
    return
  end function Chemicals_Names

  integer function Chemicals_Index(chemicalName,status)
    !!{
    Returns the index of a chemical in the chemical abundances structure given the {\normalfont \ttfamily chemicalName}.
    !!}
    use :: Error             , only : Error_Report, errorStatusFail, errorStatusSuccess
    use :: ISO_Varying_String, only : operator(==)
    implicit none
    character(len=*), intent(in   )           :: chemicalName
    integer         , intent(  out), optional :: status
    integer                                   :: iChemical

    Chemicals_Index=-1 ! Indicates chemical not found.
    do iChemical=1,chemicalsCount
       if (chemicalsToTrack(iChemical) == trim(chemicalName)) then
          Chemicals_Index=iChemical
          if (present(status)) status=errorStatusSuccess
          return
       end if
    end do
    if (present(status)) then
       status=errorStatusFail
    else
       call Error_Report('chemical species "'//trim(chemicalName)//'" is not available - to track this species add it to the <chemicalsToTrack> parameter'//{introspection:location})
    end if
    return
  end function Chemicals_Index

  subroutine Chemical_Abundances_Scale(self,multiplier)
    !!{
    In-place multiplication of a chemical abundances object.
    !!}
    implicit none
    class           (chemicalAbundances), intent(inout) :: self
    double precision                    , intent(in   ) :: multiplier

    if (chemicalsCount == 0) return
    self%chemicalValue= self%chemicalValue &
         &             *     multiplier
    return
  end subroutine Chemical_Abundances_Scale

  subroutine Chemical_Abundances_Increment(self,increment)
    !!{
    Increment an abundances object.
    !!}
    implicit none
    class(chemicalAbundances), intent(inout) :: self
    class(chemicalAbundances), intent(in   ) :: increment

    self%chemicalValue=self%chemicalValue+increment%chemicalValue
    return
  end subroutine Chemical_Abundances_Increment

  logical function Chemicals_Abundances_Is_Zero(self)
    !!{
    Test whether an chemicals object is zero.
    !!}
    implicit none
    class(chemicalAbundances), intent(in   ) :: self

    ! Detect if all chemical abundances are zero.
    Chemicals_Abundances_Is_Zero=all(self%chemicalValue == 0.0d0)
    return
  end function Chemicals_Abundances_Is_Zero

  function Chemical_Abundances_Add(abundances1,abundances2)
    !!{
    Add two abundances objects.
    !!}
    implicit none
    type (chemicalAbundances)                          :: Chemical_Abundances_Add
    class(chemicalAbundances), intent(in   )           :: abundances1
    class(chemicalAbundances), intent(in   ), optional :: abundances2

    if (chemicalsCount == 0) then
       Chemical_Abundances_Add=zeroChemicalAbundances
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
    !!{
    Subtract two abundances objects.
    !!}
    implicit none
    type (chemicalAbundances)                          :: Chemical_Abundances_Subtract
    class(chemicalAbundances), intent(in   )           :: abundances1
    class(chemicalAbundances), intent(in   ), optional :: abundances2

    if (chemicalsCount == 0) then
       Chemical_Abundances_Subtract=zeroChemicalAbundances
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
    !!{
    Multiply a chemical abundances object by a scalar.
    !!}
    implicit none
    type            (chemicalAbundances)                :: Chemical_Abundances_Multiply
    class           (chemicalAbundances), intent(in   ) :: abundances1
    double precision                    , intent(in   ) :: multiplier

    if (chemicalsCount == 0) then
       Chemical_Abundances_Multiply=zeroChemicalAbundances
    else
       Chemical_Abundances_Multiply%chemicalValue=abundances1%chemicalValue*multiplier
    end if
    return
  end function Chemical_Abundances_Multiply

  function Chemical_Abundances_Multiply_Switched(multiplier,abundances1)
    !!{
    Multiply a chemical abundances object by a scalar.
    !!}
    implicit none
    type            (chemicalAbundances)                :: Chemical_Abundances_Multiply_Switched
    double precision                    , intent(in   ) :: multiplier
    class           (chemicalAbundances), intent(in   ) :: abundances1

    if (chemicalsCount == 0) then
       Chemical_Abundances_Multiply_Switched=zeroChemicalAbundances
    else
       Chemical_Abundances_Multiply_Switched%chemicalValue=abundances1%chemicalValue*multiplier
    end if
    return
  end function Chemical_Abundances_Multiply_Switched

  function Chemical_Abundances_Divide(abundances1,divisor)
    !!{
    Divide a chemical abundances object by a scalar.
    !!}
    implicit none
    type            (chemicalAbundances)                :: Chemical_Abundances_Divide
    class           (chemicalAbundances), intent(in   ) :: abundances1
    double precision                    , intent(in   ) :: divisor

    if (chemicalsCount == 0) then
       Chemical_Abundances_Divide=zeroChemicalAbundances
    else
       Chemical_Abundances_Divide%chemicalValue=abundances1%chemicalValue/divisor
    end if
    return
  end function Chemical_Abundances_Divide

  logical function Chemical_Abundances_Greater_Than(abundances1,abundances2)
    !!{
    Multiply a chemical abundances object by a scalar.
    !!}
    implicit none
    class(chemicalAbundances), intent(in   ) :: abundances1, abundances2

    if (chemicalsCount == 0) then
       Chemical_Abundances_Greater_Than=.false.
    else
       Chemical_Abundances_Greater_Than=all(abundances1%chemicalValue > abundances2%chemicalValue)
    end if
    return
  end function Chemical_Abundances_Greater_Than

  double precision function Chemicals_Abundances(chemicals,moleculeIndex)
    !!{
    Returns the abundance of a molecule in the chemical abundances structure given the {\normalfont \ttfamily moleculeIndex}.
    !!}
    implicit none
    class  (chemicalAbundances), intent(in   ) :: chemicals
    integer                    , intent(in   ) :: moleculeIndex

    Chemicals_Abundances=chemicals%chemicalValue(moleculeIndex)
    return
  end function Chemicals_Abundances

  subroutine Chemicals_Number_To_Mass(chemicals,chemicalsByMass)
    !!{
    Multiply all chemical species by their mass in units of the atomic mass. This converts abundances by number into abundances by mass.
    !!}
    implicit none
    class(chemicalAbundances), intent(in   ) :: chemicals
    type (chemicalAbundances), intent(inout) :: chemicalsByMass

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
    !!{
    Divide all chemical species by their mass in units of the atomic mass. This converts abundances by mass into abundances by number.
    !!}
    implicit none
    class(chemicalAbundances), intent(in   ) :: chemicals
    type (chemicalAbundances), intent(inout) :: chemicalsByNumber

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
    !!{
    Force all chemical values to be positive, by truncating negative values to zero.
    !!}
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    if (allocated(chemicals%chemicalValue)) then
       where (chemicals%chemicalValue < 0.0d0)
          chemicals%chemicalValue=0.0d0
       end where
    end if
    return
  end subroutine Chemicals_Enforce_Positive

  double precision function Chemicals_Sum_Over(chemicals)
    !!{
    Return the sum over all chemicals.
    !!}
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    if (allocated(chemicals%chemicalValue)) then
       Chemicals_Sum_Over=sum(chemicals%chemicalValue)
    else
       Chemicals_Sum_Over=0.0d0
    end if
    return
  end function Chemicals_Sum_Over

  subroutine Chemicals_Builder(self,chemicalsDefinition)
    !!{
    Build a {\normalfont \ttfamily chemicalAbundances} object from the given XML {\normalfont \ttfamily chemicalsDefinition}.
    !!}
    use :: FoX_DOM           , only : node                        , extractDataContent
    use :: Error             , only : Error_Report
    use :: IO_XML            , only : XML_Get_Elements_By_Tag_Name, xmlNodeList
    use :: ISO_Varying_String, only : char
    implicit none
    class  (chemicalAbundances), intent(inout)              :: self
    type   (node              ), pointer                    :: chemicalsDefinition
    type   (node              )               , pointer     :: chemical
    type   (xmlNodeList       ), dimension(:) , allocatable :: chemicalAbundanceList
    integer                                                 :: i

    if (chemicalsCount > 0) then
       do i=1,chemicalsCount
          !$omp critical (FoX_DOM_Access)
          call XML_Get_Elements_By_Tag_Name(chemicalsDefinition,char(chemicalsToTrack(i)),chemicalAbundanceList)
          !$omp end critical (FoX_DOM_Access)
          if (size(chemicalAbundanceList) >  1) call Error_Report('multiple '//char(chemicalsToTrack(i))//' values specified'//{introspection:location})
          if (size(chemicalAbundanceList) == 1) then
             !$omp critical (FoX_DOM_Access)
             chemical => chemicalAbundanceList(0)%element
             call extractDataContent(chemical,self%chemicalValue(i))
             !$omp end critical (FoX_DOM_Access)
          end if
       end do
    end if
    return
  end subroutine Chemicals_Builder

  subroutine Chemicals_Dump(chemicals,verbosityLevel)
    !!{
    Dump all chemical values.
    !!}
    use :: Display           , only : displayMessage, enumerationVerbosityLevelType
    use :: ISO_Varying_String, only : len           , operator(//)
    implicit none
    class    (chemicalAbundances           ), intent(in   ) :: chemicals
    type     (enumerationVerbosityLevelType), intent(in   ) :: verbosityLevel
    integer                                                 :: i
    character(len=22                       )                :: label
    type     (varying_string               )                :: message

    if (allocated(chemicals%chemicalValue)) then
       do i=1,chemicalsCount
          write (label,'(e22.16)') chemicals%chemicalValue(i)
          message=chemicalsToTrack(i)//': '//repeat(" ",chemicalNameLengthMaximum-len(chemicalsToTrack(i)))//label
          call displayMessage(message,verbosityLevel)
       end do
    end if
    return
  end subroutine Chemicals_Dump

  subroutine Chemicals_Dump_Raw(chemicals,fileHandle)
    !!{
    Dump all chemical values in binary.
    !!}
    implicit none
    class  (chemicalAbundances), intent(in   ) :: chemicals
    integer                    , intent(in   ) :: fileHandle

    write (fileHandle) allocated(chemicals%chemicalValue)
    if (allocated(chemicals%chemicalValue)) write (fileHandle) chemicals%chemicalValue
    return
  end subroutine Chemicals_Dump_Raw

  subroutine Chemicals_Read_Raw(chemicals,fileHandle)
    !!{
    Read all chemical values in binary.
    !!}
    implicit none
    class  (chemicalAbundances), intent(inout) :: chemicals
    integer                    , intent(in   ) :: fileHandle
    logical                                    :: isAllocated

    read (fileHandle) isAllocated
    if (isAllocated) then
       allocate(chemicals%chemicalValue(chemicalsCount))
       read (fileHandle) chemicals%chemicalValue
    end if
    return
  end subroutine Chemicals_Read_Raw

  subroutine Chemicals_Abundances_Set(chemicals,moleculeIndex,abundance)
    !!{
    Sets the abundance of a molecule in the chemical abundances structure given the {\normalfont \ttfamily moleculeIndex}.
    !!}
    implicit none
    class           (chemicalAbundances), intent(inout) :: chemicals
    integer                             , intent(in   ) :: moleculeIndex
    double precision                    , intent(in   ) :: abundance

    ! Ensure values array exists.
    call Chemical_Abundances_Allocate_Values(chemicals)

    select type (chemicals)
    type is (chemicalAbundances)
       chemicals%chemicalValue(moleculeIndex)=abundance
    end select
    return
  end subroutine Chemicals_Abundances_Set

  subroutine Chemicals_Abundances_Reset(chemicals)
    !!{
    Resets all chemical abundances to zero.
    !!}
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
    !!{
    Resets all chemical abundances to unity.
    !!}
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
    !!{
    Destroy a chemical abundances object.
    !!}
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    if (allocated(chemicals%chemicalValue)) deallocate(chemicals%chemicalValue)
    return
  end subroutine Chemicals_Abundances_Destroy

  subroutine Chemical_Abundances_Allocate_Values(chemicals)
    !!{
    Ensure that the {\normalfont \ttfamily chemicalValue} array in an {\normalfont \ttfamily chemicalsStructure} is allocated.
    !!}
    implicit none
    class(chemicalAbundances), intent(inout) :: chemicals

    select type (chemicals)
    type is (chemicalAbundances)
       if (.not.allocated(chemicals%chemicalValue)) then
          allocate(chemicals%chemicalValue(chemicalsCount))
       end if
    end select
    return
  end subroutine Chemical_Abundances_Allocate_Values

  subroutine Chemical_Abundances_Deserialize(chemicals,chemicalAbundancesArray)
    !!{
    Pack abundances from an array into an abundances structure.
    !!}
    implicit none
    class           (chemicalAbundances)              , intent(inout) :: chemicals
    double precision                    , dimension(:), intent(in   ) :: chemicalAbundancesArray

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
    !!{
    Pack abundances from an array into an abundances structure.
    !!}
    implicit none
    double precision                    , dimension(:), intent(  out) :: chemicalAbundancesArray(:)
    class           (chemicalAbundances)              , intent(in   ) :: chemicals

    ! Place elemental values into arrays.
    if (allocated(chemicals%chemicalValue)) then
       chemicalAbundancesArray=chemicals%chemicalValue
    else
       chemicalAbundancesArray=0.0d0
    end if
    return
  end subroutine Chemical_Abundances_Serialize

  function Chemicals_Non_Static_Size_Of(self)
    !!{
    Return the size of any non-static components of the object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t          )                :: Chemicals_Non_Static_Size_Of
    class  (chemicalAbundances), intent(in   ) :: self

    if (allocated(self%chemicalValue)) then
       Chemicals_Non_Static_Size_Of=sizeof(self%chemicalValue)
    else
       Chemicals_Non_Static_Size_Of=0_c_size_t
    end if
    return
  end function Chemicals_Non_Static_Size_Of

  subroutine Chemicals_Output(self,integerProperty,integerBufferCount,integerProperties,doubleProperty,doubleBufferCount,doubleProperties,time,outputInstance)
    !!{
    Store an abundances object in the output buffers.
    !!}
    use :: Kind_Numbers                      , only : kind_int8
    use :: Multi_Counters                    , only : multiCounter
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger , outputPropertyDouble
    implicit none
    class           (chemicalAbundances   ), intent(in   )               :: self
    double precision                       , intent(in   )               :: time
    integer                                , intent(inout)               :: doubleBufferCount , doubleProperty , &
         &                                                                  integerBufferCount, integerProperty
    type            (outputPropertyInteger), intent(inout), dimension(:) :: integerProperties
    type            (outputPropertyDouble ), intent(inout), dimension(:) :: doubleProperties
    type            (multiCounter         ), intent(in   )               :: outputInstance
    integer                                                              :: i
    !$GLC attributes unused :: time, integerBufferCount, integerProperty, integerProperties, outputInstance

    if (chemicalsCount > 0) then
       do i=1,chemicalsCount
          doubleProperties(doubleProperty+i)%scalar(doubleBufferCount)=self%chemicalValue(i)
       end do
       doubleProperty=doubleProperty+chemicalsCount
    end if
    return
  end subroutine Chemicals_Output

  subroutine Chemicals_Post_Output(self,time)
    !!{
    Perform post-output processing of abundances objects.
    !!}
    implicit none
    class           (chemicalAbundances), intent(in   ) :: self
    double precision                    , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    return
  end subroutine Chemicals_Post_Output

  subroutine Chemicals_Output_Count(self,integerPropertyCount,doublePropertyCount,time)
    !!{
    Increment the output count to account for an abundances object.
    !!}
    implicit none
    class           (chemicalAbundances), intent(in   ) :: self
    integer                             , intent(inout) :: doublePropertyCount, integerPropertyCount
    double precision                    , intent(in   ) :: time
    !$GLC attributes unused :: self, integerPropertyCount, time

    doublePropertyCount=doublePropertyCount+chemicalsCount
    return
  end subroutine Chemicals_Output_Count

  subroutine Chemicals_Output_Names(self,integerProperty,integerProperties,doubleProperty,doubleProperties,time,prefix,comment,unitsInSI)
    !!{
    Assign names to output buffers for an abundances object.
    !!}
    use :: ISO_Varying_String                , only : char
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
    implicit none
    class           (chemicalAbundances   )              , intent(in   ) :: self
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleProperty    , integerProperty
    type            (outputPropertyInteger), dimension(:), intent(inout) :: integerProperties
    type            (outputPropertyDouble ), dimension(:), intent(inout) :: doubleProperties
    character       (len=*                )              , intent(in   ) :: comment           , prefix
    double precision                                     , intent(in   ) :: unitsInSI
    integer                                                              :: i
    !$GLC attributes unused :: self, time, integerProperty, integerProperties

    if (chemicalsCount > 0) then
       do i=1,chemicalsCount
          doubleProperty=doubleProperty+1
          doubleProperties(doubleProperty)%name     =trim(prefix )//      char(chemicalsToTrack(i))
          doubleProperties(doubleProperty)%comment  =trim(comment)//' ['//char(chemicalsToTrack(i))//']'
          doubleProperties(doubleProperty)%unitsInSI=unitsInSI
       end do
    end if
    return
  end subroutine Chemicals_Output_Names

end module Chemical_Abundances_Structure
