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
Contains a module which provides various atomic data.
!!}

module Atomic_Data
  !!{
  Provides various atomic data.
  !!}
  implicit none
  private
  public :: Atom_Lookup, Abundance_Pattern_Lookup, Atomic_Number, Atomic_Mass, Atomic_Abundance, Atomic_Data_Atoms_Count, Atomic_Short_Label

  type atomicData
     !!{
     Data type for storing atomic data.
     !!}
     integer                                             :: atomicNumber
     double precision                                    :: atomicMass
     double precision        , allocatable, dimension(:) :: abundanceByMass
     character       (len=3 )                            :: shortLabel
     character       (len=20)                            :: name
  end type atomicData

  ! Flag indicating if module has been initialized.
  logical                                                                   :: atomicDataInitialized =.false.

  ! Array used to store atomic data.
  type            (atomicData), allocatable, dimension(:                  ) :: atoms

  ! Metal mass normalizations.
  double precision            , allocatable, dimension(:                  ) :: metalMassNormalization

  ! Array used to index atomic data by atomic number.
  integer                     , allocatable, dimension(:                  ) :: atomicNumberIndex
  integer                                                                   :: atomicNumberMaximum

  ! Abundance pattern information.
  integer                     , parameter                                   :: abundancePatternCount =1
  character       (len=*     ), parameter, dimension(abundancePatternCount) :: abundancePatternFiles=["abundances/Solar_Composition_Cloudy_08.00.xml"], &
       &                                                                       abundancePatternNames=["solar"]

  ! Mass normalization options.
  integer                     , parameter, public                           :: normalizationTotal   =0
  integer                     , parameter, public                           :: normalizationMetals  =1

contains

  character(len=3) function Atomic_Short_Label(atomIndex,atomicNumber,name)
    !!{
    Return the short label for an atom.
    !!}
    implicit none
    integer         , intent(in   ), optional :: atomIndex, atomicNumber
    character(len=*), intent(in   ), optional :: name
    integer                                   :: iAtom

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

  integer function Atomic_Number(atomIndex,shortLabel,name)
    !!{
    Returns the atomic number of an element specified by name or short label.
    !!}
    implicit none
    integer         , intent(in   ), optional :: atomIndex
    character(len=*), intent(in   ), optional :: name     , shortLabel
    integer                                   :: iAtom

    ! Ensure the module is initialized.
    call Atomic_Data_Initialize

    ! Look up the index of this atom in the array if necessary.
    if (present(atomIndex)) then
       iAtom=atomIndex
    else
       iAtom=Atom_Lookup(shortLabel=shortLabel,name=name)
    end if

    ! Return the atomic number.
    Atomic_Number=atoms(iAtom)%atomicNumber

    return
  end function Atomic_Number

  double precision function Atomic_Mass(atomIndex,atomicNumber,shortLabel,name)
    !!{
    Returns the atomic mass of an element specified by atomic number, name or short label.
    !!}
    implicit none
    integer         , intent(in   ), optional :: atomIndex, atomicNumber
    character(len=*), intent(in   ), optional :: name     , shortLabel
    integer                                   :: iAtom

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
    !!{
    Returns the abundance by mass of a given atom in a given abundance pattern.
    !!}
    implicit none
    integer         , intent(in   ), optional :: abundanceIndex   , atomIndex, atomicNumber, &
         &                                       normalization
    character(len=*), intent(in   ), optional :: abundanceName    , name     , shortLabel
    integer                                   :: iAbundancePattern, iAtom

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
    !!{
    Return the number of atomic species known in this module.
    !!}
    implicit none

    ! Ensure the module is initialized.
    call Atomic_Data_Initialize

    ! Return the number of atomic species stored.
    Atomic_Data_Atoms_Count=size(atoms)
    return
  end function Atomic_Data_Atoms_Count

  subroutine Atomic_Data_Initialize
    !!{
    Ensure that the module is initialized by reading in data.
    !!}
    use :: FoX_dom           , only : destroy              , getElementsByTagname             , node                        , extractDataContent
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath            , pathTypeDataStatic
    use :: IO_XML            , only : XML_Array_Read_Static, XML_Get_First_Element_By_Tag_Name, XML_Get_Elements_By_Tag_Name, XML_Parse         , &
         &                            xmlNodeList
    use :: ISO_Varying_String, only : char                 , varying_string                   , assignment(=)               , operator(//)
    use :: String_Handling   , only : String_Lower_Case    , char
    implicit none
    type            (Node          )              , pointer     :: abundanceTypeElement, doc                , &
         &                                                         atom                , element
    type            (xmlNodeList   ), dimension(:), allocatable :: elementList
    integer                         , dimension(1)              :: elementValueInteger
    double precision                , dimension(1)              :: elementValueDouble
    integer                                                     :: atomicNumber        , iAbundancePattern  , &
         &                                                         iAtom               , ioErr
    double precision                                            :: abundance           , totalMass
    character       (len=100       )                            :: abundanceType
    type            (varying_string)                            :: fileName            , abundancePatternFile

    ! Check if module is initialized.
    !$omp critical(atomicDataInitialize)
    if (.not.atomicDataInitialized) then

       ! Read in the atomic data.
       fileName=inputPath(pathTypeDataStatic)//"abundances/Atomic_Data.xml"
       !$omp critical (FoX_DOM_Access)
       doc => XML_Parse(fileName,iostat=ioErr)
       if (ioErr /= 0) call Error_Report('Unable to parse data file'//{introspection:location})

       ! Get list of all element elements.
       call XML_Get_Elements_By_Tag_Name(doc,"element",elementList)

       ! Allocate storage space.
       allocate(atoms(size(elementList)))
       ! Allocate abundance pattern array for elements.
       do iAtom=1,size(elementList)
          allocate(atoms(iAtom)%abundanceByMass(abundancePatternCount))
          atoms(iAtom)%abundanceByMass=0.0d0
       end do
       ! Get atom properties.
       call XML_Array_Read_Static(elementList,"atomicNumber",atoms(:)%atomicNumber)
       call XML_Array_Read_Static(elementList,"atomicMass"  ,atoms(:)%atomicMass  )
       call XML_Array_Read_Static(elementList,"shortLabel"  ,atoms(:)%shortLabel  )
       call XML_Array_Read_Static(elementList,"name"        ,atoms(:)%name        )
       ! Destroy the document.
       call destroy(doc)
       ! Convert name to lower case.
       atoms%name=String_Lower_Case(atoms%name)
       ! Record the maximum atomic number found.
       atomicNumberMaximum=maxval(atoms%atomicNumber)

       ! Allocate space for atomic number lookup array.
       allocate(atomicNumberIndex(atomicNumberMaximum))
       ! Create lookup array by atomic number.
       forall(iAtom=1:size(atoms))
          atomicNumberIndex(atoms(iAtom)%atomicNumber)=iAtom
       end forall
       ! Allocate metal mass normalizations array.
       allocate(metalMassNormalization(abundancePatternCount))
       ! Load tables of abundance patterns.
       do iAbundancePattern=1,abundancePatternCount

          ! Parse the abundance pattern file.
          abundancePatternFile=inputPath(pathTypeDataStatic)//abundancePatternFiles(iAbundancePattern)
          doc => XML_Parse(abundancePatternFile,iostat=ioErr)
          if (ioErr /= 0) call Error_Report('Unable to parse data file'//{introspection:location})

          ! Get list of all element elements.
          call XML_Get_Elements_By_Tag_Name(doc,"element",elementList)

          ! Loop over elements.
          do iAtom=1,size(elementList)
             ! Get atom.
             atom => elementList(iAtom-1)%element
             ! Get atomic number.
             element => XML_Get_First_Element_By_Tag_Name(atom,"atomicNumber")
             call extractDataContent(element,elementValueInteger)
             atomicNumber=elementValueInteger(1)
             ! Get the abundance.
             element => XML_Get_First_Element_By_Tag_Name(atom,"abundance"   )
             call extractDataContent(element,elementValueDouble )
             abundance=elementValueDouble    (1)
             ! Store in the atoms array.
             atoms(atomicNumberIndex(atomicNumber))%abundanceByMass(iAbundancePattern)=abundance
          end do

          ! Determine the type of abundance just loaded.
          abundanceTypeElement => XML_Get_First_Element_By_Tag_Name(doc,"abundanceType")
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
             call Error_Report("unrecognized abundance type"//{introspection:location})
          end if

          ! Compute the normalization for unit metal mass in this abundance pattern.
          metalMassNormalization(iAbundancePattern)=0.0d0
          do iAtom=1,size(elementList)
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
    !$omp end critical(atomicDataInitialize)
    return
  end subroutine Atomic_Data_Initialize

  integer function Atom_Lookup(atomicNumber,shortLabel,name)
    !!{
    Returns the position in the {\normalfont \ttfamily atoms()} array of an element specified by atomic number, name or short label.
    !!}
    use :: Error          , only : Error_Report
    use :: String_Handling, only : String_Lower_Case
    implicit none
    integer          , intent(in   ), optional :: atomicNumber
    character(len=* ), intent(in   ), optional :: name         , shortLabel
    integer                                    :: iAtom
    character(len=20)                          :: nameLowerCase

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
    Atom_Lookup=-1
    call Error_Report('could not find this element'//{introspection:location})
    return
  end function Atom_Lookup

  integer function Abundance_Pattern_Lookup(abundanceIndex,abundanceName)
    !!{
    Returns the position in the {\normalfont \ttfamily atoms()} array of an element specified by atomic number, name or short label.
    !!}
    use :: Error          , only : Error_Report
    use :: String_Handling, only : String_Lower_Case
    implicit none
    integer           , intent(in   ), optional :: abundanceIndex
    character(len=*  ), intent(in   ), optional :: abundanceName
    integer                                     :: iAbundancePattern
    character(len=100)                          :: nameLowerCase

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
    Abundance_Pattern_Lookup=-1
    call Error_Report('could not find this abundance pattern'//{introspection:location})
    return
  end function Abundance_Pattern_Lookup

end module Atomic_Data
