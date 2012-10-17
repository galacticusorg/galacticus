!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  public :: abundancesStructure, Abundances_Names, Abundances_Atomic_Index, Abundances_Property_Count, Abundances_Get_Metallicity&
       &, Abundances_Mass_To_Mass_Fraction
  
  type abundancesStructure
     !% The abundances structure used for describing elemental abundances in \glc.
     private
     double precision                            :: metallicityValue
     double precision, allocatable, dimension(:) :: elementalValue
   contains
     ! Pack/unpack methods.
     !@ <objectMethods>
     !@   <object>abundancesStructure</object>
     !@   <objectMethod>
     !@     <method>pack</method>
     !@     <description>Packs abundance data into a 1-D array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unpack(array)</method>
     !@     <description>Unpacks abundance data from a 1-D array.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: pack                   => Abundances_Pack
     procedure                 :: unpack                 => Abundances_Unpack
     ! Metallicity methods.
     !@ <objectMethods>
     !@   <object>abundancesStructure</object>
     !@   <objectMethod>
     !@     <method>metallicity</method>
     !@     <description>Returns the metallicity.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>metallicitySet(metallicity)</method>
     !@     <description>Sets the metallicity to {\tt metallicity}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>massToMassFraction(mass)</method>
     !@     <description>Converts abundance masses to mass fractions by dividing by the given {\tt mass} while ensuring that fractions are in the range 0--1.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: metallicity            => Abundances_Get_Metallicity
     procedure                 :: metallicitySet         => Abundances_Set_Metallicity
     procedure                 :: massToMassFraction     => Abundances_Mass_To_Mass_Fraction_Packed
     ! Hydrogen/helium methods.
     !@ <objectMethods>
     !@   <object>abundancesStructure</object>
     !@   <objectMethod>
     !@     <method>hydrogenNumberFraction</method>
     !@     <description>Returns the hydrogen fraction by number.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>hydrogenMassFraction</method>
     !@     <description>Returns the hydrogen fraction by mass.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>heliumMassFraction</method>
     !@     <description>Returns the helium fraction by mass.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>heliumNumberFraction</method>
     !@     <description>Returns the helium fraction by number.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                 :: hydrogenNumberFraction => Abundances_Hydrogen_Number_Fraction
     procedure                 :: hydrogenMassFraction   => Abundances_Hydrogen_Mass_Fraction
     procedure                 :: heliumMassFraction     => Abundances_Helium_Mass_Fraction
     procedure                 :: heliumNumberFraction   => Abundances_Helium_Number_Fraction
  end type abundancesStructure

  ! Count of the number of elements being tracked.
  integer                                     :: elementsCount=0
  integer                                     :: propertyCount

  ! Names (two-letter labels) of elements to track.
  character(len=3), allocatable, dimension(:) :: elementsToTrack

  ! Indices of elements as used in the Atomic_Data module.
  integer,          allocatable, dimension(:) :: elementsIndices

  ! Type of metallicity/abundance measure required.
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
  integer,          parameter, public         :: adjustElementsNone  =0
  integer,          parameter, public         :: adjustElementsReset =1
  integer,          parameter, public         :: adjustElementsUpdate=2

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
          ! Flag that this module is now initialized.
          abundancesInitialized=.true.
       end if
       !$omp end critical (Abundances_Module_Initialize)
    end if
    return
  end subroutine Abundances_Initialize

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

  subroutine Abundances_Allocate_Elemental_Values(abundances)
    !% Ensure that the {\tt elementalValue} array in an {\tt abundancesStructure} is allocated.
    use Memory_Management
    implicit none
    type(abundancesStructure), intent(inout) :: abundances

    if (.not.allocated(abundances%elementalValue)) call Alloc_Array(abundances%elementalValue,[elementsCount])
    return
  end subroutine Abundances_Allocate_Elemental_Values

  subroutine Abundances_Pack(abundances,abundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    class(abundancesStructure), intent(inout)              :: abundances
    double precision,           intent(in),   dimension(:) :: abundancesArray

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (abundances)
    type is (abundancesStructure)
    ! Extract metallicity from array.
    abundances%metallicityValue=abundancesArray(1)
    ! Ensure elemental values array exists.
    call Abundances_Allocate_Elemental_Values(abundances)
    ! Extract elemental values from array.
    abundances%elementalValue=abundancesArray(2:elementsCount+1)
    end select
    return
  end subroutine Abundances_Pack

  subroutine Abundances_Unpack(abundances,abundancesArray)
    !% Pack abundances from an array into an abundances structure.
    implicit none
    double precision,           intent(out), dimension(:) :: abundancesArray(:)
    class(abundancesStructure), intent(in)                :: abundances

    ! Ensure module is initialized.
    call Abundances_Initialize

    ! Place metallicity into array.
    abundancesArray(1)=abundances%metallicityValue
    ! Place elemental values into arrays.
    if (allocated(abundances%elementalValue)) then
       abundancesArray(2:elementsCount+1)=abundances%elementalValue
    else
       abundancesArray(2:elementsCount+1)=0.0d0
    end if
    return
  end subroutine Abundances_Unpack

  double precision function Abundances_Get_Metallicity(abundances,metallicityType)
    !% Return the metallicity of the {\tt abundances} structure.
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    implicit none
    class(abundancesStructure), intent(in)           :: abundances
    integer,                    intent(in), optional :: metallicityType
    integer                                          :: metallicityTypeActual

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (abundances)
    type is (abundancesStructure)
       Abundances_Get_Metallicity=abundances%metallicityValue
       
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
             Abundances_Get_Metallicity=dlog10(Abundances_Get_Metallicity/metallicitySolar)
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

  subroutine Abundances_Set_Metallicity(abundances,metallicity,metallicityType,adjustElements,abundanceIndex)
    !% Set the metallicity of the {\tt abundances} structure to {\tt metallicity}.
    use Galacticus_Error
    use Atomic_Data
    implicit none
    class(abundancesStructure), intent(inout)        :: abundances
    double precision,           intent(in)           :: metallicity
    integer,                    intent(in), optional :: metallicityType,adjustElements,abundanceIndex
    integer                                          :: adjustElementsActual,iElement
    double precision                                 :: metallicityPrevious

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (abundances)
    type is (abundancesStructure)
    ! Store the current metallicity.
    metallicityPrevious        =abundances%metallicityValue

    ! Determine how elements will be adjusted.
    if (present(adjustElements)) then
       adjustElementsActual=adjustElements
    else
       adjustElementsActual=adjustElementsNone
    end if

    ! Store the current metallicity if necessary.
    if (elementsCount > 0 .and. adjustElementsActual == adjustElementsUpdate) metallicityPrevious=abundances%metallicityValue

    ! Update the metallicity.
    abundances%metallicityValue=metallicity
    if (present(metallicityType)) then
       select case (metallicityType)
       case (linearByMass)
          ! Do nothing, this is how we store metallicity.
       case (linearByMassSolar)
          abundances%metallicityValue=         abundances%metallicityValue *metallicitySolar
       case (logarithmicByMassSolar)
          abundances%metallicityValue=(10.0d0**abundances%metallicityValue)*metallicitySolar
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
          if (.not.present(abundanceIndex)) call Galacticus_Error_Report('Abundances_Set_Metallicity', &
               & 'an abundance pattern must be specified in order to reset elemental abundances')
          ! Ensure elemental values array exists.
          call Abundances_Allocate_Elemental_Values(abundances)
          do iElement=1,elementsCount
             abundances%elementalValue(iElement)=abundances%metallicityValue*Atomic_Abundance(abundanceIndex=abundanceIndex&
                  &,atomIndex=elementsIndices(iElement),normalization=normalizationMetals)
          end do
       case (adjustElementsUpdate)
          ! Ensure that we have an abundanceIndex specified.
          if (.not.present(abundanceIndex)) call Galacticus_Error_Report('Abundances_Set_Metallicity', &
               & 'an abundance pattern must be specified in order to reset elemental abundances')
          ! Ensure elemental values array exists.
          call Abundances_Allocate_Elemental_Values(abundances)
          do iElement=1,elementsCount
             abundances%elementalValue(iElement)=abundances%elementalValue(iElement)+(abundances%metallicityValue&
                  &-metallicityPrevious)*Atomic_Abundance(abundanceIndex=abundanceIndex,atomIndex=elementsIndices(iElement)&
                  &,normalization=normalizationMetals)
          end do
        end select
    end if
    end select
    return
  end subroutine Abundances_Set_Metallicity

  subroutine Abundances_Mass_To_Mass_Fraction_Packed(abundances,mass)
    !% Convert abundance masses to mass fractions by dividing by {\tt mass} while ensuring that the fractions remain within the range 0--1.
    implicit none
    class(abundancesStructure), intent(inout) :: abundances
    double precision,           intent(in)    :: mass

    ! Ensure module is initialized.
    call Abundances_Initialize

    ! Scale metallicity first.
    if      (abundances%metallicityValue >  mass ) then
       abundances%metallicityValue=1.0d0
    else if (abundances%metallicityValue <= 0.0d0) then
       abundances%metallicityValue=0.0d0
    else
       abundances%metallicityValue=abundances%metallicityValue/mass
    end if
    
    ! Scale elemental abundances.
    if (elementsCount > 0) then
       where     (abundances%elementalValue >  mass )
          abundances%elementalValue=1.0d0
       elsewhere (abundances%elementalValue <= 0.0d0)
          abundances%elementalValue=0.0d0
       elsewhere
          abundances%elementalValue=abundances%elementalValue/mass
       end where
    end if
    return
  end subroutine Abundances_Mass_To_Mass_Fraction_Packed

  subroutine Abundances_Mass_To_Mass_Fraction(abundances,mass)
    !% Convert abundance masses to mass fractions by dividing by {\tt mass} while ensuring that the fractions remain within the range 0--1.
    implicit none
    double precision, intent(inout), dimension(:) :: abundances
    double precision, intent(in)                  :: mass
    
    ! Scale abundances.
    where     (abundances >  mass )
       abundances=1.0d0
    elsewhere (abundances <= 0.0d0)
       abundances=0.0d0
    elsewhere
       abundances=abundances/mass
    end where
    return
  end subroutine Abundances_Mass_To_Mass_Fraction

  double precision function Abundances_Hydrogen_Mass_Fraction(abundances)
    !% Returns the mass fraction of hydrogen.
    implicit none
    class(abundancesStructure), intent(in) :: abundances
    double precision,           parameter  :: massFractionMinimum=0.7d0

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (abundances)
    type is (abundancesStructure)
    Abundances_Hydrogen_Mass_Fraction=max((abundances%metallicityValue/metallicitySolar)*(hydrogenByMassSolar&
         &-hydrogenByMassPrimordial)+hydrogenByMassPrimordial,massFractionMinimum)
    end select
    return
  end function Abundances_Hydrogen_Mass_Fraction

  double precision function Abundances_Helium_Mass_Fraction(abundances)
    !% Returns the mass fraction of helium.
    implicit none
    class(abundancesStructure), intent(in) :: abundances

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (abundances)
    type is (abundancesStructure)
    Abundances_Helium_Mass_Fraction=(abundances%metallicityValue/metallicitySolar)*(heliumByMassSolar-heliumByMassPrimordial)&
         &+heliumByMassPrimordial
    end select
    return
  end function Abundances_Helium_Mass_Fraction

  double precision function Abundances_Hydrogen_Number_Fraction(abundances)
    !% Returns the number fraction of hydrogen.
    implicit none
    class(abundancesStructure), intent(in) :: abundances
    double precision                       :: numberHydrogen,numberHelium

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (abundances)
    type is (abundancesStructure)

    numberHydrogen=Abundances_Hydrogen_Mass_Fraction(abundances)/atomicMassHydrogen
    numberHelium  =Abundances_Helium_Mass_Fraction  (abundances)/atomicMassHelium
    Abundances_Hydrogen_Number_Fraction=numberHydrogen/(numberHydrogen+numberHelium)

    end select

    return
  end function Abundances_Hydrogen_Number_Fraction

  double precision function Abundances_Helium_Number_Fraction(abundances)
    !% Returns the mass fraction of helium.
    implicit none
    class(abundancesStructure), intent(in) :: abundances
    double precision                       :: numberHydrogen,numberHelium

    ! Ensure module is initialized.
    call Abundances_Initialize

    select type (abundances)
    type is (abundancesStructure)

    numberHydrogen=Abundances_Hydrogen_Mass_Fraction(abundances)/atomicMassHydrogen
    numberHelium  =Abundances_Helium_Mass_Fraction  (abundances)/atomicMassHelium
    Abundances_Helium_Number_Fraction=numberHelium/(numberHydrogen+numberHelium)

    end select

    return
  end function Abundances_Helium_Number_Fraction

end module Abundances_Structure
