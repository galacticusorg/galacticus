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

!% Contains a module which implements structures that describe chemicals.

module Chemical_Structures
  !% Implements structures that describe chemicals.
  use ISO_Varying_String
  use Numerical_Constants_Physical
  use Numerical_Constants_Atomic
  implicit none
  private
  public :: chemicalStructure, atomicStructure,atomicBond, Chemical_Database_Get_Index

  type atomicStructure
     !% A type that defines an atom within a chemical.
     character(len=10) :: name
     character(len=2 ) :: shortLabel
     double precision  :: mass
  end type atomicStructure

  type atomicBond
     !% A type that defines an atomic bond within a chemical.
     integer :: atom(2)
  end type atomicBond

  type chemicalStructure
     !% A type that defines a chemical.
     integer                                          :: index
     type(varying_string)                             :: name
     integer                                          :: chargeValue
     double precision                                 :: massValue
     type(atomicStructure), allocatable, dimension(:) :: atom
     type(atomicBond),      allocatable, dimension(:) :: bond
   contains
     ! Data methods.
     !@ <objectMethods>
     !@   <object>chemicalStructure</object>
     !@   <objectMethod>
     !@     <method>retrieve</method>
     !@     <description>Get a chemical from the database.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>export</method>
     !@     <description>Write a chemical structure to a CML file.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                        :: retrieve => Chemical_Database_Get
     procedure                                        :: export   => Chemical_Structure_Export
     !@ <objectMethod>
     !@   <object>chemicalStructure</object>
     !@   <method>charge</method>
     !@   <description>Return the charge of a chemical.</description>
     !@ </objectMethod>
     procedure                                        :: charge   => Chemical_Structure_Charge
      !@ <objectMethod>
     !@   <object>chemicalStructure</object>
     !@   <method>mass</method>
     !@   <description>Return the mass of a chemical in atomic mass units.</description>
     !@ </objectMethod>
     procedure                                        :: mass     => Chemical_Structure_Mass
  end type chemicalStructure

  ! Atoms (we include an electron here for convenience).
  type(atomicStructure),   parameter                 :: atoms(2)=[ &
       &                                                          atomicStructure("electron","e",electronMass/atomicMassUnit), &
       &                                                          atomicStructure("hydrogen","H",atomicMassHydrogen         )  &
       &                                                         ]

  ! Chemicals.
  type(chemicalStructure), allocatable, dimension(:) :: chemicals

  ! Flag indicating if the database has been initialized.
  logical                                            :: chemicalDatabaseInitialized=.false.

contains

  subroutine Chemical_Structure_Initialize
    !% Initialize the chemical structure database by reading the atomic structure database. Note: this implementation is not
    !% fully compatible with chemical markup language (CML), but only a limited subset of it.
    use FoX_dom
    use Galacticus_Error
    use Galacticus_Input_Paths
    implicit none
    type(Node),        pointer :: doc,thisChemical,thisElement,thisAtom,thisBond
    type(NodeList),    pointer :: chemicalList,atomList,bondList,thisList
    integer                    :: iChemical,ioErr,iAtom,jAtom,iBond
    character(len=128)         :: name

    ! Check if the chemical database is initialized.
    if (.not.chemicalDatabaseInitialized) then
       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(Galacticus_Input_Path())//'data/abundances/Chemical_Database.cml',iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Chemical_Structure_Initialize','Unable to find chemical database file')
       ! Get a list of all chemicals.
       chemicalList => getElementsByTagname(doc,"chemical")
       ! Allocate the array of chemicals.
       allocate(chemicals(getLength(chemicalList)))
       ! Loop over all chemicals.
       do iChemical=0,getLength(chemicalList)-1
       ! Set the index for this chemical.
          chemicals(iChemical+1)%index=iChemical+1
          ! Get the chemical.
          thisChemical => item(chemicalList,iChemical)
          ! Get the name of the chemical.
          thisList    => getElementsByTagname(thisChemical,"id")
          thisElement => item(thisList,0)
          call extractDataContent(thisElement,name)
          chemicals(iChemical+1)%name=trim(name)
          ! Get the charge of the chemical.
          thisList    => getElementsByTagname(thisChemical,"formalCharge")
          thisElement => item(thisList,0)
          call extractDataContent(thisElement,chemicals(iChemical+1)%chargeValue)
          ! Get a list of atoms in the chemical.
          atomList => getElementsByTagname(thisChemical,"atom")
          ! Allocate array for atoms
          allocate(chemicals(iChemical+1)%atom(getLength(atomList)))
          ! Loop over atoms.
          do iAtom=0,getLength(atomList)-1
             ! Get the atom.
             thisAtom => item(atomList,iAtom)
             ! Get the element type.
             thisList    => getElementsByTagname(thisAtom,"elementType")
             thisElement => item(thisList,0)
             call extractDataContent(thisElement,name)
             ! Find the element in the list of atoms.
             do jAtom=1,size(atoms)
                if (trim(name) == trim(atoms(jAtom)%shortLabel)) chemicals(iChemical+1)%atom(iAtom+1)=atoms(jAtom)
             end do
          end do
          ! Compute the mass of the chemical.
          chemicals(iChemical+1)%massValue=sum(chemicals(iChemical+1)%atom(:)%mass)
          ! Get a list of bonds in the chemical.
          bondList => getElementsByTagname(thisChemical,"bond")
          ! Retrieve bonds if any were found.
          if (getLength(bondList) > 0) then
             ! Allocate array for bonds.
             allocate(chemicals(iChemical+1)%bond(getLength(bondList)))
             ! Loop over bond.
             do iBond=0,getLength(bondList)-1
                ! Get the bond.
                thisBond => item(bondList,iBond)
                ! Get the atom references.
                thisList    => getElementsByTagname(thisBond,"atomRefs2")
                thisElement => item(thisList,0)
                call extractDataContent(thisElement,chemicals(iChemical+1)%bond(iBond+1)%atom)
             end do
          end if
       end do
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
       ! Flag that the database is now initialized.
       chemicalDatabaseInitialized=.true.
    end if
    return
  end subroutine Chemical_Structure_Initialize

  subroutine Chemical_Structure_Export(thisChemical,outputFile)
    !% Export a chemical structure to a chemical markup language (CML) file.
    use FoX_wxml
    use FoX_dom
    use String_Handling
    implicit none
    class(chemicalStructure), intent(in) :: thisChemical
    character(len=*),         intent(in) :: outputFile
    type(xmlf_t)                         :: exportedChemicalDoc
    character(len=10)                    :: label
    integer                              :: iAtom,iBond

    ! Create the output file.
    call xml_OpenFile(outputFile,exportedChemicalDoc)
    ! Begin the chemical.
    call xml_NewElement   (exportedChemicalDoc,"chemical")
    ! Write the chemical name.
    call xml_NewElement   (exportedChemicalDoc,"id")
    call xml_AddCharacters(exportedChemicalDoc,char(thisChemical%name))
    call xml_EndElement   (exportedChemicalDoc,"id")
    ! Write the charge.
    call xml_NewElement   (exportedChemicalDoc,"formalCharge")
    call xml_AddCharacters(exportedChemicalDoc,thisChemical%chargeValue)
    call xml_EndElement   (exportedChemicalDoc,"formalCharge")
    ! Begin the array of atoms.
    call xml_NewElement   (exportedChemicalDoc,"atomArray")
    ! Loop over all atoms.
    do iAtom=1,size(thisChemical%atom)
       ! Begin the atom.
       call xml_NewElement   (exportedChemicalDoc,"atom")
       ! Write atom ID.
       write (label,'(a,i1)') "a",iAtom
       call xml_NewElement   (exportedChemicalDoc,"id")
       call xml_AddCharacters(exportedChemicalDoc,trim(label))
       call xml_EndElement   (exportedChemicalDoc,"id")
       ! Write atom element type.
       call xml_NewElement   (exportedChemicalDoc,"elementType")
       call xml_AddCharacters(exportedChemicalDoc,trim(thisChemical%atom(iAtom)%shortLabel))
       call xml_EndElement   (exportedChemicalDoc,"elementType")
       ! Close the atom.
       call xml_EndElement   (exportedChemicalDoc,"atom")
    end do
    ! End the atom array.
    call xml_EndElement   (exportedChemicalDoc,"atomArray")
    ! Check if we have bonds to output
    if (allocated(thisChemical%bond)) then
       ! Begin the array of bonds.
       call xml_NewElement   (exportedChemicalDoc,"bondArray")
       ! Loop over all bonds.
       do iBond=1,size(thisChemical%bond)
          ! Begin the bond.
          call xml_NewElement   (exportedChemicalDoc,"bond")
          ! Write the bond atomic references.
          write (label,'(a,i1,x,a,i1)') "a",thisChemical%bond(iBond)%atom(1),"a",thisChemical%bond(iBond)%atom(2)
          call xml_NewElement   (exportedChemicalDoc,"atomRefs2")
          call xml_AddCharacters(exportedChemicalDoc,trim(label))
          call xml_EndElement   (exportedChemicalDoc,"atomRefs2")
          ! Write the bond order.
          call xml_NewElement   (exportedChemicalDoc,"order")
          call xml_AddCharacters(exportedChemicalDoc,"1")
          call xml_EndElement   (exportedChemicalDoc,"order")
          ! Close the bond.
          call xml_EndElement   (exportedChemicalDoc,"bond")
       end do
       ! Close the array of bonds.
       call xml_EndElement   (exportedChemicalDoc,"bondArray")
    end if
    ! Close the chemical.
    call xml_EndElement   (exportedChemicalDoc,"chemical")
    ! Close the file.
    call xml_Close        (exportedChemicalDoc)
    return
  end subroutine Chemical_Structure_Export

  integer function Chemical_Database_Get_Index(chemicalName)
    !% Find a chemical in the database and return it.
    use Galacticus_Error
    implicit none
    character(len=*), intent(in)  :: chemicalName
    integer                       :: iChemical

    ! Initialize the database.
    call Chemical_Structure_Initialize

    ! Scan through chemicals searching for that requested.
    do iChemical=1,size(chemicals)
       if (chemicals(iChemical)%name == chemicalName) then
          Chemical_Database_Get_Index=iChemical
          return
       end if
    end do
    call Galacticus_Error_Report('Chemical_Database_Get_Index','chemical was not found in database')
    return
  end function Chemical_Database_Get_Index

  subroutine Chemical_Database_Get(thisChemical,chemicalName)
    !% Find a chemical in the database and return it.
    use Galacticus_Error
    implicit none
    class(chemicalStructure), intent(inout) :: thisChemical

    character(len=*),         intent(in)    :: chemicalName

    select type (thisChemical)
    type is (chemicalStructure)
       thisChemical=chemicals(Chemical_Database_Get_Index(chemicalName))
    end select
    return
  end subroutine Chemical_Database_Get

  integer function Chemical_Structure_Charge(thisChemical)
    !% Return the charge on a chemical.
    use Galacticus_Error
    implicit none
    class(chemicalStructure), intent(in) :: thisChemical

    Chemical_Structure_Charge=thisChemical%chargeValue
    return
  end function Chemical_Structure_Charge

  double precision function Chemical_Structure_Mass(thisChemical)
    !% Return the mass of a chemical.
    use Galacticus_Error
    implicit none
    class(chemicalStructure), intent(in) :: thisChemical

    Chemical_Structure_Mass=thisChemical%massValue
    return
  end function Chemical_Structure_Mass

end module Chemical_Structures
