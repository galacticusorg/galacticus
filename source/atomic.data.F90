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


!% Contains a module which provides various atomic data.

module Atomic_Data
  !% Provides various atomic data.
  implicit none
  private
  public :: Atom_Lookup, Abundance_Pattern_Lookup, Atomic_Mass, Atomic_Abundance, Atomic_Data_Atoms_Count, Atomic_Short_Label

  type atomicData
     !% Data type for storing atomic data.
     integer                                     :: atomicNumber
     double precision                            :: atomicMass
     double precision, allocatable, dimension(:) :: abundanceByMass
     character(len=3)                            :: shortLabel
     character(len=20)                           :: name
  end type atomicData

  ! Flag indicating if module has been initialized.
  logical                                                       :: atomicDataInitialized=.false.

  ! Array used to store atomic data.
  type(atomicData), allocatable, dimension(:)                   :: atoms

  ! Metal mass normalizations.
  double precision, allocatable, dimension(:)                   :: metalMassNormalization

  ! Array used to index atomic data by atomic number.
  integer,          allocatable, dimension(:)                   :: atomicNumberIndex
  integer                                                       :: atomicNumberMaximum

  ! Abundance pattern information.
  integer,          parameter                                   :: abundancePatternCount=1
  character(len=*), parameter, dimension(abundancePatternCount) :: abundancePatternFiles=           &
       &                                                            ["data/Solar_Composition.xml"], & 
       &                                                           abundancePatternNames=           &
       &                                                            ["solar"]

  ! Mass normalization options.
  integer,          parameter, public                           :: normalizationTotal =0
  integer,          parameter, public                           :: normalizationMetals=1

contains

  character(len=3) function Atomic_Short_Label(atomIndex,atomicNumber,name)
    !% Return the short label for an atom.
    implicit none
    integer,          intent(in), optional :: atomIndex,atomicNumber
    character(len=*), intent(in), optional :: name
    integer                                :: iAtom
    
    ! Ensure the module is initialized.
    call Atomic_Data_Initialize
    
    ! Look up the index of this atom in the array if necessary.
    if (present(atomIndex)) then
       iAtom=atomIndex
    else
       iAtom=Atom_Lookup(atomicNumber=atomicNumber,name=name)
    end if

    Atomic_Short_Label=atoms(iAtom)%shortLabel
    return
  end function Atomic_Short_Label

  double precision function Atomic_Mass(atomIndex,atomicNumber,shortLabel,name)
    !% Returns the atomic mass of an element specified by atomic number, name or short label.
    implicit none
    integer,          intent(in), optional :: atomIndex,atomicNumber
    character(len=*), intent(in), optional :: shortLabel,name
    integer                                :: iAtom

    ! Ensure the module is initialized.
    call Atomic_Data_Initialize
    
    ! Look up the index of this atom in the array if necessary.
    if (present(atomIndex)) then
       iAtom=atomIndex
    else
       iAtom=Atom_Lookup(atomicNumber,shortLabel,name)
    end if

    ! Return the atomic mass.
    Atomic_Mass=atoms(iAtom)%atomicMass

    return
  end function Atomic_Mass

  double precision function Atomic_Abundance(abundanceIndex,abundanceName,atomIndex,atomicNumber,shortLabel,name,normalization)
    !% Returns the abundance by mass of a given atom in a given abundance pattern.
    implicit none
    integer,          intent(in), optional :: abundanceIndex,atomIndex,atomicNumber,normalization
    character(len=*), intent(in), optional :: abundanceName,shortLabel,name
    integer                                :: iAbundancePattern,iAtom

    ! Ensure the module is initialized.
    call Atomic_Data_Initialize
    
    ! Look up the index of this atom in the array.
    iAbundancePattern=Abundance_Pattern_Lookup(abundanceIndex,abundanceName)
    
    ! Look up the index of this atom in the array if necessary.
    if (present(atomIndex)) then
       iAtom=atomIndex
    else
       iAtom=Atom_Lookup(atomicNumber,shortLabel,name)
    end if
    
    ! Return the atomic mass.
    Atomic_Abundance=atoms(iAtom)%abundanceByMass(iAbundancePattern)

    ! Normalize the abundance as requested.
    if (present(normalization)) then
       select case (normalization)
       case (normalizationMetals)
          Atomic_Abundance=Atomic_Abundance*metalMassNormalization(iAbundancePattern)
       end select
    end if

    return
  end function Atomic_Abundance

  integer function Atomic_Data_Atoms_Count()
    !% Return the number of atomic species known in this module.
    implicit none
    
    ! Ensure the module is initialized.
    call Atomic_Data_Initialize
    
    ! Return the number of atomic species stored.
    Atomic_Data_Atoms_Count=size(atoms)
    return
  end function Atomic_Data_Atoms_Count

  subroutine Atomic_Data_Initialize
    !% Ensure that the module is initialized by reading in data.
    use FoX_dom
    use Memory_Management
    use Galacticus_Error
    use String_Handling
    implicit none
    type(Node),       pointer      :: doc,thisElement,abundanceTypeElement,thisAtom
    type(NodeList),   pointer      :: elementList
    integer,          dimension(1) :: elementValueInteger
    double precision, dimension(1) :: elementValueDouble
    integer                        :: ioErr,iAtom,iAbundancePattern,atomicNumber
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
       allocate(atoms(getLength(elementList)))
       call Memory_Usage_Record(sizeof(atoms))

       ! Extract the data into our array.
       atomicNumberMaximum=0
       do iAtom=1,getLength(elementList)
          ! Allocate abundance pattern array for this element.
          call Alloc_Array(atoms(iAtom)%abundanceByMass,[abundancePatternCount])
          atoms(iAtom)%abundanceByMass=0.0d0
          ! Get atom.
          thisAtom => item(elementList,iAtom-1)
          ! Get atomic number.
          thisElement => item(getElementsByTagname(thisAtom,"atomicNumber"),0)
          call extractDataContent(thisElement,elementValueInteger    )
          atoms(iAtom)%atomicNumber=elementValueInteger(1)
          ! Get atomic mass.
          thisElement => item(getElementsByTagname(thisAtom,"atomicMass"  ),0)
          call extractDataContent(thisElement,elementValueDouble     )
          atoms(iAtom)%atomicMass  =elementValueDouble(1)
          ! Get short label.
          thisElement => item(getElementsByTagname(thisAtom,"shortLabel"  ),0)
          call extractDataContent(thisElement,atoms(iAtom)%shortLabel)
          ! Get name.
          thisElement => item(getElementsByTagname(thisAtom,"name"        ),0)
          call extractDataContent(thisElement,atoms(iAtom)%name      )
          ! Convert name to lower case.
          atoms(iAtom)%name=String_Lower_Case(atoms(iAtom)%name)
          ! Record the maximum atomic number found.
          atomicNumberMaximum=max(atomicNumberMaximum,atoms(iAtom)%atomicNumber)
       end do

       ! Allocate space for atomic number lookup array.
       call Alloc_Array(atomicNumberIndex,[atomicNumberMaximum])

       ! Create lookup array by atomic number.
       forall(iAtom=1:size(atoms))
          atomicNumberIndex(atoms(iAtom)%atomicNumber)=iAtom
       end forall

       ! Destroy the document.
       call destroy(doc)

       ! Allocate metal mass normalizations array.
       call Alloc_Array(metalMassNormalization,[abundancePatternCount])

       ! Load tables of abundance patterns.
       do iAbundancePattern=1,abundancePatternCount

          ! Parse the abundance pattern file.
          doc => parseFile(abundancePatternFiles(iAbundancePattern),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('Atomic_Data_Initialize','Unable to parse data file')

          ! Get list of all element elements.
          elementList => getElementsByTagname(doc,"element")

          ! Loop over elements.
          do iAtom=1,getLength(elementList)
             ! Get atom.
             thisAtom => item(elementList,iAtom-1)
             ! Get atomic number.
             thisElement => item(getElementsByTagname(thisAtom,"atomicNumber"),0)
             call extractDataContent(thisElement,elementValueInteger)
             atomicNumber=elementValueInteger(1)
             ! Get the abundance.
             thisElement => item(getElementsByTagname(thisAtom,"abundance"),0)
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

          ! Compute the normalization for unit metal mass in this abundance pattern.
          metalMassNormalization(iAbundancePattern)=0.0d0
          do iAtom=1,getLength(elementList)
             if (atoms(iAtom)%atomicNumber > 2) metalMassNormalization(iAbundancePattern)&
                  &=metalMassNormalization(iAbundancePattern)+atoms(iAtom)%abundanceByMass(iAbundancePattern)
          end do
          metalMassNormalization(iAbundancePattern)=1.0d0/metalMassNormalization(iAbundancePattern)

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

    ! Ensure the module is initialized.
    call Atomic_Data_Initialize
 
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

    ! Ensure the module is initialized.
    call Atomic_Data_Initialize
 
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
