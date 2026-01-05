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
Contains a module which implements structures that describe chemicals.
!!}

module Chemical_Structures
  !!{
  Implements structures that describe chemicals.
  !!}
  use :: ISO_Varying_String          , only : varying_string
  use :: Numerical_Constants_Atomic  , only : atomicMassHydrogen, atomicMassUnit
  use :: Numerical_Constants_Physical, only : electronMass
  implicit none
  private
  public :: chemicalStructure, atomicStructure,atomicBond, Chemical_Database_Get_Index

  type atomicStructure
     !!{
     A type that defines an atom within a chemical.
     !!}
     character       (len=10) :: name
     character       (len=2 ) :: shortLabel
     double precision         :: mass
  end type atomicStructure

  type atomicBond
     !!{
     A type that defines an atomic bond within a chemical.
     !!}
     integer :: atom(2)
  end type atomicBond

  type chemicalStructure
     !!{
     A type that defines a chemical.
     !!}
     integer                                                      :: index
     type            (varying_string )                            :: name
     integer                                                      :: chargeValue
     double precision                                             :: massValue
     type            (atomicStructure), allocatable, dimension(:) :: atom
     type            (atomicBond     ), allocatable, dimension(:) :: bond
   contains
     ! Data methods.
     !![
     <methods>
       <method description="Get a chemical from the database."                   method="retrieve"/>
       <method description="Write a chemical structure to a CML file."           method="export"  />
       <method description="Return the charge of a chemical."                    method="charge"  />
       <method description="Return the mass of a chemical in atomic mass units." method="mass"    />
     </methods>
     !!]
     procedure :: retrieve => Chemical_Database_Get
     procedure :: export   => Chemical_Structure_Export
     procedure :: charge   => Chemical_Structure_Charge
     procedure :: mass     => Chemical_Structure_Mass
  end type chemicalStructure

  ! Atoms (we include an electron here for convenience).
  type(atomicStructure),   parameter                 :: atoms(2)=[                                                             &
       &                                                          atomicStructure("electron","e",electronMass/atomicMassUnit), &
       &                                                          atomicStructure("hydrogen","H",atomicMassHydrogen         )  &
       &                                                         ]

  ! Chemicals.
  type   (chemicalStructure), allocatable, dimension(:) :: chemicals

  ! Flag indicating if the database has been initialized.
  logical                                               :: chemicalDatabaseInitialized=.false.

contains

  subroutine Chemical_Structure_Initialize
    !!{
    Initialize the chemical structure database by reading the atomic structure database. Note: this implementation is not
    fully compatible with chemical markup language (CML), but only a limited subset of it.
    !!}
    use :: FoX_dom           , only : Node                        , destroy           , extractDataContent
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath                   , pathTypeDataStatic
    use :: IO_XML            , only : XML_Get_Elements_By_Tag_Name, xmlNodeList       , XML_Parse
    use :: ISO_Varying_String, only : char                        , assignment(=)     , varying_string
    implicit none
    type     (Node          ), pointer                   :: doc     , atom    , bond        , chemical, element
    type     (xmlNodeList   ), allocatable, dimension(:) :: atomList, bondList, chemicalList, list
    integer                                              :: iAtom   , iBond   , iChemical   , ioErr   , jAtom
    character(len=128       )                            :: name
    type     (varying_string)                           :: fileName

    ! Check if the chemical database is initialized.
    if (.not.chemicalDatabaseInitialized) then
       fileName=char(inputPath(pathTypeDataStatic))//'abundances/Chemical_Database.cml'
       !$omp critical (FoX_DOM_Access)
       doc => XML_Parse(char(fileName),iostat=ioErr)
       if (ioErr /= 0) call Error_Report('Unable to find chemical database file'//{introspection:location})
       ! Get a list of all chemicals.
       call XML_Get_Elements_By_Tag_Name(doc,"chemical",chemicalList)
       ! Allocate the array of chemicals.
       allocate(chemicals(size(chemicalList)))
       ! Loop over all chemicals.
       do iChemical=0,size(chemicalList)-1
       ! Set the index for this chemical.
          chemicals(iChemical+1)%index=iChemical+1
          ! Get the chemical.
          chemical => chemicalList(iChemical)%element
          ! Get the name of the chemical.
          call XML_Get_Elements_By_Tag_Name(chemical,"id",list)
          element => list(0)%element
          call extractDataContent(element,name)
          chemicals(iChemical+1)%name=trim(name)
          ! Get the charge of the chemical.
          call XML_Get_Elements_By_Tag_Name(chemical,"formalCharge",list)
          element => list(0)%element
          call extractDataContent(element,chemicals(iChemical+1)%chargeValue)
          ! Get a list of atoms in the chemical.
          call XML_Get_Elements_By_Tag_Name(chemical,"atom",atomList)
          ! Allocate array for atoms
          allocate(chemicals(iChemical+1)%atom(size(atomList)))
          ! Loop over atoms.
          do iAtom=0,size(atomList)-1
             ! Get the atom.
             atom => atomList(iAtom)%element
             ! Get the element type.
             call XML_Get_Elements_By_Tag_Name(atom,"elementType",list)
             element => list(0)%element
             call extractDataContent(element,name)
             ! Find the element in the list of atoms.
             do jAtom=1,size(atoms)
                if (trim(name) == trim(atoms(jAtom)%shortLabel)) chemicals(iChemical+1)%atom(iAtom+1)=atoms(jAtom)
             end do
          end do
          ! Compute the mass of the chemical.
          chemicals(iChemical+1)%massValue=sum(chemicals(iChemical+1)%atom(:)%mass)
          ! Get a list of bonds in the chemical.
          call XML_Get_Elements_By_Tag_Name(chemical,"bond",bondList)
          ! Retrieve bonds if any were found.
          if (size(bondList) > 0) then
             ! Allocate array for bonds.
             allocate(chemicals(iChemical+1)%bond(size(bondList)))
             ! Loop over bond.
             do iBond=0,size(bondList)-1
                ! Get the bond.
                bond => bondList(iBond)%element
                ! Get the atom references.
                call XML_Get_Elements_By_Tag_Name(bond,"atomRefs2",list)
                element => list(0)%element
                call extractDataContent(element,chemicals(iChemical+1)%bond(iBond+1)%atom)
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

  subroutine Chemical_Structure_Export(chemical,outputFile)
    !!{
    Export a chemical structure to a chemical markup language (CML) file.
    !!}
    use :: FoX_wxml          , only : xml_AddCharacters, xml_Close, xml_EndElement, xml_NewElement, &
          &                           xml_OpenFile     , xmlf_t
    use :: ISO_Varying_String, only : char
    implicit none
    class    (chemicalStructure), intent(in   ) :: chemical
    character(len=*            ), intent(in   ) :: outputFile
    type     (xmlf_t           )                :: exportedChemicalDoc
    character(len=10           )                :: label
    integer                                     :: iAtom              , iBond

    ! Create the output file.
    call xml_OpenFile(outputFile,exportedChemicalDoc)
    ! Begin the chemical.
    call xml_NewElement   (exportedChemicalDoc,"chemical")
    ! Write the chemical name.
    call xml_NewElement   (exportedChemicalDoc,"id")
    call xml_AddCharacters(exportedChemicalDoc,char(chemical%name))
    call xml_EndElement   (exportedChemicalDoc,"id")
    ! Write the charge.
    call xml_NewElement   (exportedChemicalDoc,"formalCharge")
    call xml_AddCharacters(exportedChemicalDoc,chemical%chargeValue)
    call xml_EndElement   (exportedChemicalDoc,"formalCharge")
    ! Begin the array of atoms.
    call xml_NewElement   (exportedChemicalDoc,"atomArray")
    ! Loop over all atoms.
    do iAtom=1,size(chemical%atom)
       ! Begin the atom.
       call xml_NewElement   (exportedChemicalDoc,"atom")
       ! Write atom ID.
       write (label,'(a,i1)') "a",iAtom
       call xml_NewElement   (exportedChemicalDoc,"id")
       call xml_AddCharacters(exportedChemicalDoc,trim(label))
       call xml_EndElement   (exportedChemicalDoc,"id")
       ! Write atom element type.
       call xml_NewElement   (exportedChemicalDoc,"elementType")
       call xml_AddCharacters(exportedChemicalDoc,trim(chemical%atom(iAtom)%shortLabel))
       call xml_EndElement   (exportedChemicalDoc,"elementType")
       ! Close the atom.
       call xml_EndElement   (exportedChemicalDoc,"atom")
    end do
    ! End the atom array.
    call xml_EndElement   (exportedChemicalDoc,"atomArray")
    ! Check if we have bonds to output
    if (allocated(chemical%bond)) then
       ! Begin the array of bonds.
       call xml_NewElement   (exportedChemicalDoc,"bondArray")
       ! Loop over all bonds.
       do iBond=1,size(chemical%bond)
          ! Begin the bond.
          call xml_NewElement   (exportedChemicalDoc,"bond")
          ! Write the bond atomic references.
          write (label,'(a,i1,x,a,i1)') "a",chemical%bond(iBond)%atom(1),"a",chemical%bond(iBond)%atom(2)
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
    !!{
    Find a chemical in the database and return it.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : operator(==)
    implicit none
    character(len=*), intent(in   ) :: chemicalName
    integer                         :: iChemical

    ! Initialize the database.
    call Chemical_Structure_Initialize()

    ! Scan through chemicals searching for that requested.
    do iChemical=1,size(chemicals)
       if (chemicals(iChemical)%name == chemicalName) then
          Chemical_Database_Get_Index=iChemical
          return
       end if
    end do
    Chemical_Database_Get_Index=-1
    call Error_Report('chemical was not found in database'//{introspection:location})
    return
  end function Chemical_Database_Get_Index

  subroutine Chemical_Database_Get(chemical,chemicalName)
    !!{
    Find a chemical in the database and return it.
    !!}
    implicit none
    class    (chemicalStructure), intent(inout) :: chemical
    character(len=*            ), intent(in   ) :: chemicalName

    select type (chemical)
    type is (chemicalStructure)
       chemical=chemicals(Chemical_Database_Get_Index(chemicalName))
    end select
    return
  end subroutine Chemical_Database_Get

  integer function Chemical_Structure_Charge(chemical)
    !!{
    Return the charge on a chemical.
    !!}
    implicit none
    class(chemicalStructure), intent(in   ) :: chemical

    Chemical_Structure_Charge=chemical%chargeValue
    return
  end function Chemical_Structure_Charge

  double precision function Chemical_Structure_Mass(chemical)
    !!{
    Return the mass of a chemical.
    !!}
    implicit none
    class(chemicalStructure), intent(in   ) :: chemical

    Chemical_Structure_Mass=chemical%massValue
    return
  end function Chemical_Structure_Mass

end module Chemical_Structures
