!% Contains a module which provides various atomic data.

module Atomic_Data
  !% Provides various atomic data.
  use Atomic_Data_Type
  private
  public :: Atom_Lookup, Abundance_Pattern_Lookup, Atomic_Mass, Atomic_Abundance

  ! Flag indicating if module has been initialized.
  logical                                                       :: atomicDataInitialized=.false.

  ! Array used to store atomic data.
  type(atomicData), allocatable, dimension(:)                   :: atoms

  ! Array used to index atomic data by atomic number.
  integer,          allocatable, dimension(:)                   :: atomicNumberIndex
  integer                                                       :: atomicNumberMaximum

  ! Abundance pattern information.
  integer,          parameter                                   :: abundancePatternCount=1
  character(len=*), parameter, dimension(abundancePatternCount) :: abundancePatternFiles=           &
       &                                                            ["data/Solar_Composition.xml"], & 
       &                                                           abundancePatternNames=           &
       &                                                            ["solar"]

contains

  double precision function Atomic_Mass(atomicNumber,shortLabel,name)
    !% Returns the atomic mass of an element specified by atomic number, name or short label.
    implicit none
    integer,          intent(in), optional :: atomicNumber
    character(len=*), intent(in), optional :: shortLabel,name
    integer                                :: iAtom

    ! Ensure the module is initialized.
    call Atomic_Data_Initialize
    
    ! Look up the index of this atom in the array.
    iAtom=Atom_Lookup(atomicNumber,shortLabel,name)
    
    ! Return the atomic mass.
    Atomic_Mass=atoms(iAtom)%atomicMass

    return
  end function Atomic_Mass

  double precision function Atomic_Abundance(abundanceIndex,abundanceName,atomicNumber,shortLabel,name)
    !% Returns the abundance by mass of a given atom in a given abundance pattern.
    implicit none
    integer,          intent(in), optional :: abundanceIndex,atomicNumber
    character(len=*), intent(in), optional :: abundanceName,shortLabel,name
    integer                                :: iAbundancePattern,iAtom

    ! Ensure the module is initialized.
    call Atomic_Data_Initialize
    
    ! Look up the index of this atom in the array.
    iAbundancePattern=Abundance_Pattern_Lookup(abundanceIndex,abundanceName)
    
    ! Look up the index of this atom in the array.
    iAtom=Atom_Lookup(atomicNumber,shortLabel,name)
    
    ! Return the atomic mass.
    Atomic_Abundance=atoms(iAtom)%abundanceByMass(iAbundancePattern)

    return
  end function Atomic_Abundance

  subroutine Atomic_Data_Initialize
    !% Ensure that the module is initialized by reading in data.
    use FoX_dom
    use Memory_Management
    use Galacticus_Error
    use String_Handling
    implicit none
    type(Node),       pointer      :: doc,thisElement,abundanceTypeElement
    type(NodeList),   pointer      :: elementList
    integer,          dimension(1) :: elementValueInteger
    double precision, dimension(1) :: elementValueDouble
    integer                        :: ioErr,iAtom,iAbundancePattern,atomsIndex,atomicNumber
    double precision               :: totalMass,abundance
    character(len=100)             :: abundanceType
  
    ! Check if module is initialized.
    if (.not.atomicDataInitialized) then

       ! Read in the atomic data.
       !$omp critical (FoX_DOM_Access)
       doc => parseFile("./data/Atomic_Data.xml",iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Atomic_Data_Initialize','Unable to parse data file')

       ! Get list of all element elements.
       elementList => getElementsByTagname(doc,"element")
       
       ! Allocate storage space.
       call Alloc_Array(atoms,getLength(elementList),'atoms')

       ! Extract the data into our array.
       atomicNumberMaximum=0
       do iAtom=1,getLength(elementList)
          ! Allocate abundance pattern array for this element.
          call Alloc_Array(atoms(iAtom)%abundanceByMass,abundancePatternCount,"abundanceByMass")
          atoms(iAtom)%abundanceByMass=0.0d0
          ! Get atomic number.
          thisElement => item(getElementsByTagname(item(elementList,iAtom-1),"atomicNumber"),0)
          call extractDataContent(thisElement,elementValueInteger    )
          atoms(iAtom)%atomicNumber=elementValueInteger(1)
          ! Get atomic mass.
          thisElement => item(getElementsByTagname(item(elementList,iAtom-1),"atomicMass"  ),0)
          call extractDataContent(thisElement,elementValueDouble     )
          atoms(iAtom)%atomicMass  =elementValueDouble(1)
          ! Get short label.
          thisElement => item(getElementsByTagname(item(elementList,iAtom-1),"shortLabel"  ),0)
          call extractDataContent(thisElement,atoms(iAtom)%shortLabel)
          ! Get name.
          thisElement => item(getElementsByTagname(item(elementList,iAtom-1),"name"        ),0)
          call extractDataContent(thisElement,atoms(iAtom)%name      )
          ! Convert name to lower case.
          atoms(iAtom)%name=String_Lower_Case(atoms(iAtom)%name)
          ! Record the maximum atomic number found.
          atomicNumberMaximum=max(atomicNumberMaximum,atoms(iAtom)%atomicNumber)
       end do

       ! Allocate space for atomic number lookup array.
       call Alloc_Array(atomicNumberIndex,atomicNumberMaximum,'atomicNumberIndex')

       ! Create lookup array by atomic number.
       forall(iAtom=1:size(atoms))
          atomicNumberIndex(atoms(iAtom)%atomicNumber)=iAtom
       end forall

       ! Destroy the document.
       call destroy(doc)

       ! Load tables of abundance patterns.
       do iAbundancePattern=1,abundancePatternCount

          ! Parse the abundance pattern file.
          doc => parseFile(abundancePatternFiles(iAbundancePattern),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('Atomic_Data_Initialize','Unable to parse data file')

          ! Get list of all element elements.
          elementList => getElementsByTagname(doc,"element")

          ! Loop over elements.
          do iAtom=1,getLength(elementList)
             ! Get atomic number.
             thisElement => item(getElementsByTagname(item(elementList,iAtom-1),"atomicNumber"),0)
             call extractDataContent(thisElement,elementValueInteger)
             atomicNumber=elementValueInteger(1)
             ! Get the abundance.
             thisElement => item(getElementsByTagname(item(elementList,iAtom-1),"abundance"),0)
             call extractDataContent(thisElement,elementValueDouble )
             abundance=elementValueDouble(1)
             ! Store in the atoms array.
             atoms(atomicNumberIndex(atomicNumber))%abundanceByMass(iAbundancePattern)=abundance
          end do

          ! Determine the type of abundance just loaded.
          abundanceTypeElement => item(getElementsByTagname(doc,"abundanceType"),0)
          call extractDataContent(abundanceTypeElement,abundanceType)
          if (trim(abundanceType) == "number relative to hydrogen") then
             ! Convert to abundances by mass.
             totalMass=0.0d0
             do iAtom=1,size(atoms)
                totalMass=totalMass+atoms(iAtom)%abundanceByMass(iAbundancePattern)*atoms(iAtom)%atomicMass
             end do
             do iAtom=1,size(atoms)
                atoms(iAtom)%abundanceByMass(iAbundancePattern)=atoms(iAtom)%abundanceByMass(iAbundancePattern)&
                     &*atoms(iAtom)%atomicMass/totalMass
             end do
          else
             call Galacticus_Error_Report("Atomic_Data_Initialize","unrecognized abundance type")
          end if

          ! Destroy the document.
          call destroy(doc)
       end do
       !$omp end critical (FoX_DOM_Access)
       

       ! Mark module as initialized.
       atomicDataInitialized=.true.
    end if
    
    return
  end subroutine Atomic_Data_Initialize

  integer function Atom_Lookup(atomicNumber,shortLabel,name)
    !% Returns the position in the {\tt atoms()} array of an element specified by atomic number, name or short label.
    use String_Handling
    use Galacticus_Error
    implicit none
    integer,          intent(in), optional :: atomicNumber
    character(len=*), intent(in), optional :: shortLabel,name
    integer                                :: iAtom
    character(len=20)                      :: nameLowerCase
 
    ! Look up by atomic number if present.
    if (present(atomicNumber)) then
       Atom_Lookup=atomicNumberIndex(atomicNumber)
       return
    end if

    ! Look up by short label if present.
    if (present(shortLabel)) then
       do iAtom=1,size(atoms)
          if (trim(atoms(iAtom)%shortLabel) == trim(shortLabel)) then
             Atom_Lookup=iAtom
             return
          end if
       end do
    end if

    ! Look up by name if present.
    if (present(name)) then
       nameLowerCase=String_Lower_Case(name)
       do iAtom=1,size(atoms)
          if (trim(atoms(iAtom)%name) == trim(nameLowerCase)) then
             Atom_Lookup=iAtom
             return
          end if
       end do
    end if

    ! Element was not found, report an error.
    call Galacticus_Error_Report('Atom_Lookup','could not find this element')
    return
  end function Atom_Lookup

  integer function Abundance_Pattern_Lookup(abundanceIndex,abundanceName)
    !% Returns the position in the {\tt atoms()} array of an element specified by atomic number, name or short label.
    use String_Handling
    use Galacticus_Error
    implicit none
    integer,          intent(in), optional :: abundanceIndex
    character(len=*), intent(in), optional :: abundanceName
    integer                                :: iAbundancePattern
    character(len=100)                     :: nameLowerCase
 
    ! Look up by abundance index if present.
    if (present(abundanceIndex)) then
       Abundance_Pattern_Lookup=abundanceIndex
       return
    end if

    ! Look up by name if present.
    if (present(abundanceName)) then
       nameLowerCase=String_Lower_Case(abundanceName)
       do iAbundancePattern=1,abundancePatternCount
          if (trim(abundancePatternNames(iAbundancePattern)) == trim(nameLowerCase)) then
             Abundance_Pattern_Lookup=iAbundancePattern
             return
          end if
       end do
    end if

    ! Abundance pattern was not found, report an error.
    call Galacticus_Error_Report('Abundance_Pattern_Lookup','could not find this abundance pattern')
    return
  end function Abundance_Pattern_Lookup

end module Atomic_Data
