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
     type(varying_string   )          :: label
   contains
     final     ::           mergeDestructor
     procedure :: import => mergeImport
     procedure :: isHDF5 => mergeIsHDF5
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

    !# <inputParameter>
    !#   <name>label</name>
    !#   <source>parameters</source>
    !#   <variable>self%label</variable>
    !#   <description>A label for the simulation</description>
    !#   <defaultValue>var_str('*')</defaultValue>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
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

  function mergeConstructorInternal(label,importers) result (self)
    !% Internal constructor for the {\normalfont \ttfamily merge} N-body importer class.
    implicit none
    type(nbodyImporterMerge)                         :: self
    type(nbodyImporterList ), target , intent(in   ) :: importers
    type(varying_string    )         , intent(in   ) :: label
    type(nbodyImporterList ), pointer                :: importer_
    !# <constructorAssign variables="label"/>
    
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

  subroutine mergeImport(self,simulations)
    !% Merge data from multiple importers.
    use :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, verbosityStandard
    use :: Hashes            , only : rank1IntegerSizeTHash    , rank1DoubleHash
    implicit none
    class           (nbodyImporterMerge), intent(inout)                            :: self
    type            (nBodyData         ), intent(  out), allocatable, dimension(:) :: simulations
    type            (nbodyImporterList )               , pointer                   :: importer_
    integer         (c_size_t          )               , allocatable, dimension(:) :: propertyInteger
    double precision                                   , allocatable, dimension(:) :: propertyReal
    integer                                                                        :: j              , k
    integer         (c_size_t          )                                           :: countObjects

    call Galacticus_Display_Indent('merging imported data',verbosityStandard)
    allocate(simulations(1))
    countObjects          =  0_c_size_t
    importer_             => self%importers
    do while (associated(importer_))
       call importer_%importer_%import(importer_%simulations)
       do j=1,size(importer_%simulations)
          countObjects =  countObjects+size(importer_%simulations(j)%particleIDs)
       end do
       importer_ => importer_%next
    end do
    allocate(simulations(1)%particleIDs    (  countObjects))
    allocate(simulations(1)%position       (3,countObjects))
    allocate(simulations(1)%velocity       (3,countObjects))
    allocate(               propertyInteger(  countObjects))
    allocate(               propertyReal   (  countObjects))
    ! Set a label.
    if (self%label == "*") then
       simulations(1)%label=self%importers%simulations(1)%label
    else
       simulations(1)%label=self%label
    end if
    ! Default properties.
    countObjects =  0_c_size_t
    importer_    => self%importers
    do while (associated(importer_))
       do j=1,size(importer_%simulations)
          simulations   (1)%particleIDs(  countObjects+1:countObjects+size(importer_%simulations(j)%particleIDs))=importer_%simulations(j)%particleIDs
          do k=1,3
             simulations(1)%position   (k,countObjects+1:countObjects+size(importer_%simulations(j)%particleIDs))=importer_%simulations(j)%position   (k,:)
             simulations(1)%velocity   (k,countObjects+1:countObjects+size(importer_%simulations(j)%particleIDs))=importer_%simulations(j)%velocity   (k,:)
          end do
          countObjects=countObjects+size(importer_%simulations(j)%particleIDs)
       end do
       importer_ => importer_%next
    end do
    ! Extra properties.
    !! Integer properties.
    simulations(1)%propertiesInteger=rank1IntegerSizeTHash()
    if (self%importers%simulations(1)%propertiesInteger%size() > 0) then
       do k=1,self%importers%simulations(1)%propertiesInteger%size()
          countObjects = 0_c_size_t
          importer_    => self%importers
          do while (associated(importer_))
             do j=1,size(importer_%simulations)
                propertyInteger(countObjects+1:countObjects+size(importer_%simulations(k)%particleIDs))=importer_%simulations(j)%propertiesInteger%value(k)
                countObjects=countObjects+size(importer_%simulations(j)%particleIDs)
             end do
             importer_ => importer_%next
          end do
          call simulations(1)%propertiesInteger%set(self%importers%simulations(1)%propertiesInteger%key(k),propertyInteger)
       end do
    end if
    !! Real properties.
    simulations(1)%propertiesReal   =rank1DoubleHash      ()
    if (self%importers%simulations(1)%propertiesReal   %size() > 0) then
       do k=1,self%importers%simulations(1)%propertiesReal   %size()
          countObjects = 0_c_size_t
          importer_    => self%importers
          do while (associated(importer_))
             do j=1,size(importer_%simulations)
                propertyReal   (countObjects+1:countObjects+size(importer_%simulations(k)%particleIDs))=importer_%simulations(j)%propertiesReal   %value(k)
                countObjects=countObjects+size(importer_%simulations(j)%particleIDs)
             end do
             importer_ => importer_%next
          end do
          call simulations(1)%propertiesReal   %set(self%importers%simulations(1)%propertiesReal   %key(k),propertyReal   )
       end do
    end if
    ! Close analysis groups in merged importers.
    importer_ => self%importers
    do while (associated(importer_))
       do j=1,size(importer_%simulations)
          if (importer_%simulations(j)%analysis%isOpen()) call importer_%simulations(j)%analysis%close()
      end do
       importer_ => importer_%next
    end do
    call Galacticus_Display_Unindent('done',verbosityStandard)
    return
  end subroutine mergeImport

  logical function mergeIsHDF5(self)
    !% Return whether or not the imported data is from an HDF5 file.
    implicit none
    class(nbodyImporterMerge), intent(inout) :: self

    mergeIsHDF5=.false.
    return
  end function mergeIsHDF5
