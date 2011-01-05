!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements structures that describe molecules.

module Molecular_Structures
  !% Implements structures that describe molecules.
  use ISO_Varying_String
  use Numerical_Constants_Physical
  use Numerical_Constants_Atomic
  private
  public :: molecularStructure, atomicStructure,atomicBond, Molecular_Database_Get_Index

  type atomicStructure
     !% A type that defines an atom within a molecule.
     character(len=10) :: name
     character(len=2 ) :: shortLabel
     double precision  :: mass
  end type atomicStructure

  type atomicBond
     !% A type that defines an atomic bond within a molecule.
     integer :: atom(2)
  end type atomicBond

  type molecularStructure
     !% A type that defines a molecule.
     integer                                          :: index
     type(varying_string)                             :: name
     integer                                          :: chargeValue
     double precision                                 :: massValue
     type(atomicStructure), allocatable, dimension(:) :: atom
     type(atomicBond),      allocatable, dimension(:) :: bond
   contains
     ! Data methods.
     !@ <objectMethods>
     !@   <object>molecularStructure</object>
     !@   <objectMethod>
     !@     <method>retrieve</method>
     !@     <description>Get a molecule from the database.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>export</method>
     !@     <description>Write a molecular structure to a CML file.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                        :: retrieve => Molecular_Database_Get
     procedure                                        :: export   => Molecular_Structure_Export
     !@ <objectMethod>
     !@   <object>molecularStructure</object>
     !@   <method>charge</method>
     !@   <description>Return the charge of a molecule.</description>
     !@ </objectMethod>
     procedure                                        :: charge   => Molecular_Structure_Charge
      !@ <objectMethod>
     !@   <object>molecularStructure</object>
     !@   <method>mass</method>
     !@   <description>Return the mass of a molecule in atomic mass units.</description>
     !@ </objectMethod>
     procedure                                        :: mass     => Molecular_Structure_Mass
  end type molecularStructure

  ! Atoms (we include an electron here for convenience).
  type(atomicStructure),    parameter                 :: atoms(2)=[ &
       &                                                           atomicStructure("electron","e",electronMass/atomicMassUnit), &
       &                                                           atomicStructure("hydrogen","H",atomicMassHydrogen         )  &
       &                                                          ]

  ! Molecules.
  type(molecularStructure), allocatable, dimension(:) :: molecules

  ! Flag indicating if the database has been initialized.
  logical                                             :: molecularDatabaseInitialized=.false.

contains

  subroutine Molecular_Structure_Initialize
    !% Initialize the molecular structure database by reading the atomic structure database. Note: this implementation is not
    !% fully compatible with chemical markup language (CML), but only a limited subset of it.
    use FoX_dom
    use Galacticus_Error
    implicit none
    type(Node),        pointer :: doc,thisMolecule,thisElement,thisAtom,thisBond
    type(NodeList),    pointer :: moleculeList,atomList,bondList,thisList
    integer                    :: iMolecule,ioErr,iAtom,jAtom,iBond
    character(len=128)         :: name

    ! Check if the molecular database is initialized.
    if (.not.molecularDatabaseInitialized) then
       !$omp critical (FoX_DOM_Access)
       doc => parseFile('data/Molecular_Database.cml',iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Molecular_Structure_Initialize','Unable to find molecular database file')
       ! Get a list of all molecules.
       moleculeList => getElementsByTagname(doc,"molecule")
       ! Allocate the array of molecules.
       allocate(molecules(getLength(moleculeList)))
       ! Loop over all molecules.
       do iMolecule=0,getLength(moleculeList)-1
       ! Set the index for this molecule.
          molecules(iMolecule+1)%index=iMolecule+1
          ! Get the molecule.
          thisMolecule => item(moleculeList,iMolecule)
          ! Get the name of the molecule.
          thisList    => getElementsByTagname(thisMolecule,"id")
          thisElement => item(thisList,0)
          call extractDataContent(thisElement,name)
          molecules(iMolecule+1)%name=trim(name)
          ! Get the charge of the molecule.
          thisList    => getElementsByTagname(thisMolecule,"formalCharge")
          thisElement => item(thisList,0)
          call extractDataContent(thisElement,molecules(iMolecule+1)%chargeValue)
          ! Get a list of atoms in the molecule.
          atomList => getElementsByTagname(thisMolecule,"atom")
          ! Allocate array for atoms
          allocate(molecules(iMolecule+1)%atom(getLength(atomList)))
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
                if (trim(name) == trim(atoms(jAtom)%shortLabel)) molecules(iMolecule+1)%atom(iAtom+1)=atoms(jAtom)
             end do
          end do
          ! Compute the mass of the molecule.
          molecules(iMolecule+1)%massValue=sum(molecules(iMolecule+1)%atom(:)%mass)
          ! Get a list of bonds in the molecule.
          bondList => getElementsByTagname(thisMolecule,"bond")
          ! Retrieve bonds if any were found.
          if (getLength(bondList) > 0) then
             ! Allocate array for bonds.
             allocate(molecules(iMolecule+1)%bond(getLength(bondList)))
             ! Loop over bond.
             do iBond=0,getLength(bondList)-1
                ! Get the bond.
                thisBond => item(bondList,iBond)
                ! Get the atom references.
                thisList    => getElementsByTagname(thisBond,"atomRefs2")
                thisElement => item(thisList,0)
                call extractDataContent(thisElement,molecules(iMolecule+1)%bond(iBond+1)%atom)
             end do
          end if
       end do
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
       ! Flag that the database is now initialized.
       molecularDatabaseInitialized=.true.
    end if
    return
  end subroutine Molecular_Structure_Initialize

  subroutine Molecular_Structure_Export(thisMolecule,outputFile)
    !% Export a molecular structure to a chemical markup language (CML) file.
    use FoX_wxml
    use FoX_dom
    use String_Handling
    implicit none
#ifdef GCC45
    class(molecularStructure), intent(in) :: thisMolecule
#else
    type(molecularStructure),  intent(in) :: thisMolecule
#endif
    character(len=*),          intent(in) :: outputFile
    type(xmlf_t)                          :: exportedMoleculeDoc
    character(len=10)                     :: label
    integer                               :: iAtom,iBond

    ! Create the output file.
    call xml_OpenFile(outputFile,exportedMoleculeDoc)
    ! Begin the molecule.
    call xml_NewElement   (exportedMoleculeDoc,"molecule")
    ! Write the molecule name.
    call xml_NewElement   (exportedMoleculeDoc,"id")
    call xml_AddCharacters(exportedMoleculeDoc,char(thisMolecule%name))
    call xml_EndElement   (exportedMoleculeDoc,"id")
    ! Write the charge.
    call xml_NewElement   (exportedMoleculeDoc,"formalCharge")
    call xml_AddCharacters(exportedMoleculeDoc,thisMolecule%chargeValue)
    call xml_EndElement   (exportedMoleculeDoc,"formalCharge")
    ! Begin the array of atoms.
    call xml_NewElement   (exportedMoleculeDoc,"atomArray")
    ! Loop over all atoms.
    do iAtom=1,size(thisMolecule%atom)
       ! Begin the atom.
       call xml_NewElement   (exportedMoleculeDoc,"atom")
       ! Write atom ID.
       write (label,'(a,i1)') "a",iAtom
       call xml_NewElement   (exportedMoleculeDoc,"id")
       call xml_AddCharacters(exportedMoleculeDoc,trim(label))
       call xml_EndElement   (exportedMoleculeDoc,"id")
       ! Write atom element type.
       call xml_NewElement   (exportedMoleculeDoc,"elementType")
       call xml_AddCharacters(exportedMoleculeDoc,trim(thisMolecule%atom(iAtom)%shortLabel))
       call xml_EndElement   (exportedMoleculeDoc,"elementType")
       ! Close the atom.
       call xml_EndElement   (exportedMoleculeDoc,"atom")
    end do
    ! End the atom array.
    call xml_EndElement   (exportedMoleculeDoc,"atomArray")
    ! Check if we have bonds to output
    if (allocated(thisMolecule%bond)) then
       ! Begin the array of bonds.
       call xml_NewElement   (exportedMoleculeDoc,"bondArray")
       ! Loop over all bonds.
       do iBond=1,size(thisMolecule%bond)
          ! Begin the bond.
          call xml_NewElement   (exportedMoleculeDoc,"bond")
          ! Write the bond atomic references.
          write (label,'(a,i1,x,a,i1)') "a",thisMolecule%bond(iBond)%atom(1),"a",thisMolecule%bond(iBond)%atom(2)
          call xml_NewElement   (exportedMoleculeDoc,"atomRefs2")
          call xml_AddCharacters(exportedMoleculeDoc,trim(label))
          call xml_EndElement   (exportedMoleculeDoc,"atomRefs2")
          ! Write the bond order.
          call xml_NewElement   (exportedMoleculeDoc,"order")
          call xml_AddCharacters(exportedMoleculeDoc,"1")
          call xml_EndElement   (exportedMoleculeDoc,"order")
          ! Close the bond.
          call xml_EndElement   (exportedMoleculeDoc,"bond")
       end do
       ! Close the array of bonds.
       call xml_EndElement   (exportedMoleculeDoc,"bondArray")
    end if
    ! Close the molecule.
    call xml_EndElement   (exportedMoleculeDoc,"molecule")
    ! Close the file.
    call xml_Close        (exportedMoleculeDoc)
    return
  end subroutine Molecular_Structure_Export

  integer function Molecular_Database_Get_Index(moleculeName)
    !% Find a molecule in the database and return it.
    use Galacticus_Error
    implicit none
    character(len=*), intent(in)  :: moleculeName
    integer                       :: iMolecule

    ! Initialize the database.
    call Molecular_Structure_Initialize

    ! Scan through molecules searching for that requested.
    do iMolecule=1,size(molecules)
       if (molecules(iMolecule)%name == moleculeName) then
          Molecular_Database_Get_Index=iMolecule
          return
       end if
    end do
    call Galacticus_Error_Report('Molecular_Database_Get_Index','molecule was not found in database')
    return
  end function Molecular_Database_Get_Index

  subroutine Molecular_Database_Get(thisMolecule,moleculeName)
    !% Find a molecule in the database and return it.
    use Galacticus_Error
    implicit none
#ifdef GCC45
    class(molecularStructure), intent(inout) :: thisMolecule
#else
    type(molecularStructure),  intent(inout) :: thisMolecule
#endif

    character(len=*),         intent(in)    :: moleculeName

#ifdef GCC45
    select type (thisMolecule)
    type is (molecularStructure)
#endif
       thisMolecule=molecules(Molecular_Database_Get_Index(moleculeName))
#ifdef GCC45
    end select
#endif
    return
  end subroutine Molecular_Database_Get

  integer function Molecular_Structure_Charge(thisMolecule)
    !% Return the charge on a molecule.
    use Galacticus_Error
    implicit none
#ifdef GCC45
    class(molecularStructure), intent(in) :: thisMolecule
#else
    type(molecularStructure),  intent(in) :: thisMolecule
#endif

    Molecular_Structure_Charge=thisMolecule%chargeValue
    return
  end function Molecular_Structure_Charge

  double precision function Molecular_Structure_Mass(thisMolecule)
    !% Return the mass of a molecule.
    use Galacticus_Error
    implicit none
#ifdef GCC45
    class(molecularStructure), intent(in) :: thisMolecule
#else
    type(molecularStructure),  intent(in) :: thisMolecule
#endif

    Molecular_Structure_Mass=thisMolecule%massValue
    return
  end function Molecular_Structure_Mass

end module Molecular_Structures
