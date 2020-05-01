!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements an N-body data importer which merges data from other importers.
  
  type, public :: nbodyImporterList
     class(nbodyImporterClass), pointer :: importer_
     type (nbodyImporterList ), pointer :: next      => null()
  end type nbodyImporterList

  !# <nbodyImporter name="nbodyImporterMerge">
  !#  <description>An importer which merges data from other importers.</description>
  !#  <deepCopy>
  !#   <linkedList type="nbodyImporterList" variable="importers" next="next" object="importer_" objectType="nbodyImporterClass"/>
  !#  </deepCopy>
  !# </nbodyImporter>
  type, extends(nbodyImporterClass) :: nbodyImporterMerge
     !% An importer which merges data from other importers.
     private
     type(nbodyImporterList), pointer :: importers => null()
   contains
     final     ::           mergeDestructor
     procedure :: import => mergeImport
  end type nbodyImporterMerge

  interface nbodyImporterMerge
     !% Constructors for the {\normalfont \ttfamily merge} N-body importer class.
     module procedure mergeConstructorParameters
     module procedure mergeConstructorInternal
  end interface nbodyImporterMerge

contains

  function mergeConstructorParameters(parameters) result (self)
    !% Constructor for the {\normalfont \ttfamily merge} N-body importer class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyImporterMerge)                :: self
    type   (inputParameters   ), intent(inout) :: parameters
    type   (nbodyImporterList ), pointer       :: importer_
    integer                                    :: i

    self     %importers => null()
    importer_           => null()
    do i=1,parameters%copiesCount('nbodyImporterMethod',zeroIfNotPresent=.true.)
       if (associated(importer_)) then
          allocate(importer_%next)
          importer_ => importer_%next
       else
          allocate(self%importers)
          importer_ => self%importers
       end if
       !# <objectBuilder class="nbodyImporter" name="importer_%importer_" source="parameters" copy="i" />
    end do
    return
  end function mergeConstructorParameters

  function mergeConstructorInternal(importers) result (self)
    !% Internal constructor for the {\normalfont \ttfamily merge} N-body importer class.
    implicit none
    type(nbodyImporterMerge)                        :: self
    type(nbodyImporterList ), target, intent(in   ) :: importers
    type(nbodyImporterList ), pointer               :: importer_

    self     %importers => importers
    importer_           => importers
    do while (associated(importer_))
       !# <referenceCountIncrement owner="importer_" object="importer_"/>
       importer_ => importer_%next
    end do
    return
  end function mergeConstructorInternal

  subroutine mergeDestructor(self)
    !% Destructor for {\normalfont \ttfamily merge} importer class.
    implicit none
    type(nbodyImporterMerge), intent(inout) :: self
    type(nbodyImporterList ), pointer       :: importer_, importerNext

    if (associated(self%importers)) then
       importer_ => self%importers
       do while (associated(importer_))
          importerNext => importer_%next
          !# <objectDestructor name="importer_%importer_"/>
          deallocate(importer_)
          importer_ => importerNext
       end do
    end if
    return
  end subroutine mergeDestructor

  function mergeImport(self) result(simulation)
    !% Merge data from multiple importers.
    use :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, verbosityStandard
    use :: Hashes            , only : rank1IntegerSizeTHash    , rank1DoubleHash
    implicit none
    type            (nBodyData         )                            :: simulation
    class           (nbodyImporterMerge), intent(inout)             :: self
    type            (nbodyImporterList ), pointer                   :: importer_
    type            (nBodyData         ), allocatable, dimension(:) :: simulations
    integer         (c_size_t          ), allocatable, dimension(:) :: propertyInteger
    double precision                    , allocatable, dimension(:) :: propertyReal
    integer                                                         :: countImporters , i, &
         &                                                             j
    integer         (c_size_t          )                            :: countObjects

    call Galacticus_Display_Indent('merging imported data',verbosityStandard)
    countImporters =  0
    importer_      => self%importers
    do while (associated(importer_))
       countImporters =  countImporters+1
       importer_      => importer_%next
    end do
    allocate(simulations(countImporters))
    i            = 0
    countObjects = 0_c_size_t
    importer_    => self%importers
    do while (associated(importer_))
       i              =  i+1
       simulations(i) =  importer_%importer_%import()
       importer_      => importer_%next
       countObjects   =  countObjects+size(simulations(i)%particleIDs)
    end do
    allocate(simulation%particleIDs    (  countObjects))
    allocate(simulation%position       (3,countObjects))
    allocate(simulation%velocity       (3,countObjects))
    allocate(           propertyInteger(  countObjects))
    allocate(           propertyReal   (  countObjects))
    countObjects=0_c_size_t
    ! Default properties.
    do i=1,countImporters
       simulation%particleIDs(countObjects+1:countObjects+size(simulations(i)%particleIDs))=simulations(i)%particleIDs
       do j=1,3
          simulation%position(j,countObjects+1:countObjects+size(simulations(i)%particleIDs))=simulations(i)%position(j,:)
          simulation%velocity(j,countObjects+1:countObjects+size(simulations(i)%particleIDs))=simulations(i)%velocity(j,:)
       end do
       countObjects=countObjects+size(simulations(i)%particleIDs)
    end do
    ! Extra properties.
    !! Integer properties.
    simulation%propertiesInteger=rank1IntegerSizeTHash()
    if (simulations(1)%propertiesInteger%size() > 0) then
       do j=1,simulations(1)%propertiesInteger%size()
          countObjects=0_c_size_t
          do i=1,countImporters
             propertyInteger(countObjects+1:countObjects+size(simulations(i)%particleIDs))=simulations(i)%propertiesInteger%value(j)
             countObjects=countObjects+size(simulations(i)%particleIDs)
          end do
          call simulation%propertiesInteger%set(simulations(i)%propertiesInteger%key(j),propertyInteger)
       end do
    end if
    !! Real properties.
    simulation%propertiesReal   =rank1DoubleHash      ()
    if (simulations(1)%propertiesReal   %size() > 0) then
       do j=1,simulations(1)%propertiesReal   %size()
          countObjects=0_c_size_t
          do i=1,countImporters
             propertyReal   (countObjects+1:countObjects+size(simulations(i)%particleIDs))=simulations(i)%propertiesReal   %value(j)
             countObjects=countObjects+size(simulations(i)%particleIDs)
          end do
          call simulation%propertiesReal   %set(simulations(i)%propertiesReal   %key(j),propertyReal   )
       end do
    end if
    call Galacticus_Display_Unindent('done',verbosityStandard)
    return
  end function mergeImport
