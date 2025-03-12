!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements an N-body data importer which imports using multiple other importers.
!!}
  
  !![
  <nbodyImporter name="nbodyImporterMultiple">
    <description>An importer which imports using multiple other importers.</description>
    <linkedList type="nbodyImporterList" variable="importers" next="next" object="importer_" objectType="nbodyImporterClass"/>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterMultiple
     !!{
     An importer which imports using multiple other importers.
     !!}
     private
     type   (nbodyImporterList), pointer :: importers => null()
     logical                             :: allHDF5   =  .true.
   contains
     final     ::           multipleDestructor
     procedure :: import => multipleImport
     procedure :: isHDF5 => multipleIsHDF5
  end type nbodyImporterMultiple

  interface nbodyImporterMultiple
     !!{
     Constructors for the {\normalfont \ttfamily multiple} N-body importer class.
     !!}
     module procedure multipleConstructorParameters
     module procedure multipleConstructorInternal
  end interface nbodyImporterMultiple

contains

  function multipleConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily multiple} N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyImporterMultiple)                :: self
    type   (inputParameters      ), intent(inout) :: parameters
    type   (nbodyImporterList    ), pointer       :: importer_
    integer                                       :: i

    self     %allHDF5   =  .true.
    self     %importers => null()
    importer_           => null()
    do i=1,parameters%copiesCount('nbodyImporter',zeroIfNotPresent=.true.)
       if (associated(importer_)) then
          allocate(importer_%next)
          importer_ => importer_%next
       else
          allocate(self%importers)
          importer_ => self%importers
       end if
       !![
       <objectBuilder class="nbodyImporter" name="importer_%importer_" source="parameters" copy="i" />
       !!]
       self%allHDF5= self               %allHDF5   &
            &       .and.                          &
            &        importer_%importer_% isHDF5()
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="nbodyImporter"/>
    !!]
    return
  end function multipleConstructorParameters

  function multipleConstructorInternal(importers) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily multiple} N-body importer class.
    !!}
    implicit none
    type(nbodyImporterMultiple)                        :: self
    type(nbodyImporterList    ), target, intent(in   ) :: importers
    type(nbodyImporterList    ), pointer               :: importer_

    self     %allHDF5   =  .true.
    self     %importers => importers
    importer_           => importers
    do while (associated(importer_))
       !![
       <referenceCountIncrement owner="importer_" object="importer_"/>
       !!]
       self     %allHDF5 =   self               %allHDF5   &
            &               .and.                          &
            &                importer_%importer_% isHDF5()
       importer_         =>  importer_%next
    end do
    return
  end function multipleConstructorInternal

  subroutine multipleDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily multiple} importer class.
    !!}
    implicit none
    type(nbodyImporterMultiple), intent(inout) :: self
    type(nbodyImporterList    ), pointer       :: importer_, importerNext

    if (associated(self%importers)) then
       importer_ => self%importers
       do while (associated(importer_))
          importerNext => importer_%next
          !![
          <objectDestructor name="importer_%importer_"/>
          !!]
          deallocate(importer_)
          importer_ => importerNext
       end do
    end if
    return
  end subroutine multipleDestructor

  subroutine multipleImport(self,simulations)
    !!{
    Import data using multiple importers.
    !!}
    use :: Display, only : displayIndent     , displayUnindent         , verbosityLevelStandard
    use :: Hashes , only : rank1DoublePtrHash, rank1IntegerSizeTPtrHash, rank2DoublePtrHash    , rank2IntegerSizeTPtrHash, &
         &                 doubleHash        , integerSizeTHash        , varyingStringHash     , genericHash
    implicit none
    class  (nbodyImporterMultiple), intent(inout)                            :: self
    type   (nBodyData            ), intent(  out), allocatable, dimension(:) :: simulations
    type   (nbodyImporterList    )               , pointer                   :: importer_
    integer                                                                  :: i
    integer(c_size_t             )                                           :: countSimulations

    call displayIndent('merging imported data',verbosityLevelStandard)
    countSimulations =  0_c_size_t
    importer_        => self%importers
    do while (associated(importer_))
       call importer_%importer_%import(importer_%simulations)
       countSimulations=countSimulations+size(importer_%simulations)
       importer_ => importer_%next
    end do
    allocate(simulations(countSimulations))
    ! Default properties.
    countSimulations =  0_c_size_t
    importer_        => self%importers
    do while (associated(importer_))
       do i=1,size(importer_%simulations)
          countSimulations                  =countSimulations+1_c_size_t
          simulations     (countSimulations)=importer_%simulations(i)
       end do
       importer_ => importer_%next
    end do
    ! Remove pointers to simulation data in combined importers.
    importer_ => self%importers
    do while (associated(importer_))
       do i=1,size(importer_%simulations)
          importer_%simulations(i)%propertiesInteger     =rank1IntegerSizeTPtrHash()
          importer_%simulations(i)%propertiesIntegerRank1=rank2IntegerSizeTPtrHash()
          importer_%simulations(i)%propertiesReal        =rank1DoublePtrHash      ()
          importer_%simulations(i)%propertiesRealRank1   =rank2DoublePtrHash      ()
          importer_%simulations(1)%attributesInteger     =integerSizeTHash        ()
          importer_%simulations(1)%attributesReal        =doubleHash              ()
          importer_%simulations(1)%attributesText        =varyingStringHash       ()
          importer_%simulations(1)%attributesGeneric     =genericHash             ()
       end do
       importer_ => importer_%next
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine multipleImport

  logical function multipleIsHDF5(self)
    !!{
    Return whether or not the imported data is from an HDF5 file.
    !!}
    implicit none
    class(nbodyImporterMultiple), intent(inout) :: self

    multipleIsHDF5=self%allHDF5
    return
  end function multipleIsHDF5
