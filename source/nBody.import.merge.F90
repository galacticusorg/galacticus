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
Implements an N-body data importer which merges data from other importers.
!!}
  
  !![
  <nbodyImporter name="nbodyImporterMerge">
   <description>An importer which merges data from other importers.</description>
   <linkedList type="nbodyImporterList" variable="importers" next="next" object="importer_" objectType="nbodyImporterClass"/>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterMerge
     !!{
     An importer which merges data from other importers.
     !!}
     private
     type(nbodyImporterList), pointer :: importers => null()
     type(varying_string   )          :: label
   contains
     final     ::           mergeDestructor
     procedure :: import => mergeImport
     procedure :: isHDF5 => mergeIsHDF5
  end type nbodyImporterMerge

  interface nbodyImporterMerge
     !!{
     Constructors for the \refClass{nbodyImporterMerge} N-body importer class.
     !!}
     module procedure mergeConstructorParameters
     module procedure mergeConstructorInternal
  end interface nbodyImporterMerge

contains

  function mergeConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyImporterMerge} N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyImporterMerge)                :: self
    type   (inputParameters   ), intent(inout) :: parameters
    type   (nbodyImporterList ), pointer       :: importer_
    integer                                    :: i

    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <variable>self%label</variable>
      <description>A label for the simulation</description>
      <defaultValue>var_str('*')</defaultValue>
    </inputParameter>
    !!]
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
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="nbodyImporter"/>
    !!]
    return
  end function mergeConstructorParameters

  function mergeConstructorInternal(label,importers) result (self)
    !!{
    Internal constructor for the \refClass{nbodyImporterMerge} N-body importer class.
    !!}
    implicit none
    type(nbodyImporterMerge)                         :: self
    type(nbodyImporterList ), target , intent(in   ) :: importers
    type(varying_string    )         , intent(in   ) :: label
    type(nbodyImporterList ), pointer                :: importer_
    !![
    <constructorAssign variables="label"/>
    !!]
    
    self     %importers => importers
    importer_           => importers
    do while (associated(importer_))
       !![
       <referenceCountIncrement owner="importer_" object="importer_"/>
       !!]
       importer_ => importer_%next
    end do
    return
  end function mergeConstructorInternal

  subroutine mergeDestructor(self)
    !!{
    Destructor for the \refClass{nbodyImporterMerge} importer class.
    !!}
    implicit none
    type(nbodyImporterMerge), intent(inout) :: self
    type(nbodyImporterList ), pointer       :: importer_, importerNext

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
  end subroutine mergeDestructor

  subroutine mergeImport(self,simulations)
    !!{
    Merge data from multiple importers.
    !!}
    use :: Display, only : displayIndent     , displayUnindent         , verbosityLevelStandard
    use :: Error  , only : Error_Report
    use :: Hashes , only : doubleHash        , integerSizeTHash        , rank1DoublePtrHash    , rank1IntegerSizeTPtrHash, &
          &                rank2DoublePtrHash, rank2IntegerSizeTPtrHash, varyingStringHash     , genericHash
    implicit none
    class           (nbodyImporterMerge), intent(inout)                              :: self
    type            (nBodyData         ), intent(  out), allocatable, dimension(  :) :: simulations
    type            (nbodyImporterList )               , pointer                     :: importer_
    integer         (c_size_t          )               , pointer    , dimension(  :) :: propertyInteger     , propertyIntegerMerged
    double precision                                   , pointer    , dimension(  :) :: propertyReal        , propertyRealMerged
    integer         (c_size_t          )               , pointer    , dimension(:,:) :: propertyIntegerRank1, propertyIntegerRank1Merged
    double precision                                   , pointer    , dimension(:,:) :: propertyRealRank1   , propertyRealRank1Merged
    integer                                                                          :: j                   , k
    integer         (c_size_t          )                                             :: countObjects        , countObjectsMerged
    
    call displayIndent('merging imported data',verbosityLevelStandard)
    allocate(simulations(1))
    countObjectsMerged    =  0_c_size_t
    importer_             => self%importers
    do while (associated(importer_))
       call importer_%importer_%import(importer_%simulations)
       do j=1,size(importer_%simulations)
          if (importer_%simulations(j)%propertiesInteger          %size() > 0) then
             propertyInteger      =>  importer_%simulations(j)%propertiesInteger     %value(1)
             countObjectsMerged   =  +countObjectsMerged               &
                  &                  +size(propertyInteger     ,dim=1)
          else if (importer_%simulations(j)%propertiesReal        %size() > 0) then
             propertyReal         =>  importer_%simulations(j)%propertiesReal        %value(1)
             countObjectsMerged   =  +countObjectsMerged               &
                  &                  +size(propertyReal        ,dim=1)
          else if (importer_%simulations(j)%propertiesIntegerRank1%size() > 0) then
             propertyIntegerRank1 =>  importer_%simulations(j)%propertiesIntegerRank1%value(1)
             countObjectsMerged   =  +countObjectsMerged               &
                  &                  +size(propertyIntegerRank1,dim=2)
          else if (importer_%simulations(j)%propertiesRealRank1   %size() > 0) then
             propertyRealRank1    =>  importer_%simulations(j)%propertiesRealRank1   %value(1)
             countObjectsMerged   =  +countObjectsMerged               &
                  &                  +size(propertyRealRank1   ,dim=2)
          else
             call Error_Report('no properties are available in the simulations'//{introspection:location})
          end if
       end do
       importer_ => importer_%next
    end do    
    ! Set a label.
    if (self%label == "*") then
       simulations(1)%label=self%importers%simulations(1)%label
    else
       simulations(1)%label=self%label
    end if
    ! Extra attributes/properties.
    !! Integer attributes.
    simulations(1)%attributesInteger=integerSizeTHash     ()
    if (self%importers%simulations(1)%attributesInteger%size() > 0) then
       do k=1,self%importers%simulations(1)%attributesInteger%size()
          call simulations(1)%attributesInteger%set(self%importers%simulations(1)%attributesInteger%key(k),self%importers%simulations(1)%attributesInteger%value(k))
       end do
    end if
    !! Real attributes.
    simulations(1)%attributesReal   =doubleHash           ()
    if (self%importers%simulations(1)%attributesReal   %size() > 0) then
       do k=1,self%importers%simulations(1)%attributesReal   %size()
          call simulations(1)%attributesReal   %set(self%importers%simulations(1)%attributesReal   %key(k),self%importers%simulations(1)%attributesReal   %value(k))
       end do
    end if
    !! Text attributes.
    simulations(1)%attributesText   =varyingStringHash    ()
    if (self%importers%simulations(1)%attributesText   %size() > 0) then
       do k=1,self%importers%simulations(1)%attributesText   %size()
          call simulations(1)%attributesText   %set(self%importers%simulations(1)%attributesText   %key(k),self%importers%simulations(1)%attributesText   %value(k))
       end do
    end if
    !! Generic attributes.
    simulations(1)%attributesGeneric=genericHash          ()
    if (self%importers%simulations(1)%attributesGeneric%size() > 0) then
       do k=1,self%importers%simulations(1)%attributesGeneric%size()
          call simulations(1)%attributesGeneric%set(self%importers%simulations(1)%attributesGeneric%key(k),self%importers%simulations(1)%attributesGeneric%value(k))
       end do
    end if
    !! Scalar integer properties.
    simulations(1)%propertiesInteger=rank1IntegerSizeTPtrHash()
    if (self%importers%simulations(1)%propertiesInteger%size() > 0) then
       do k=1,self%importers%simulations(1)%propertiesInteger%size()
          countObjects          = 0_c_size_t
          importer_             => self%importers
          propertyIntegerMerged => null()
          do while (associated(importer_))
             do j=1,size(importer_%simulations)
                propertyInteger                                                                =>  importer_%simulations(j)%propertiesInteger%value(k)
                if (.not.associated(propertyIntegerMerged)) allocate(propertyIntegerMerged(countObjectsMerged))
                propertyIntegerMerged(countObjects+1:countObjects+size(propertyInteger,dim=1)) =   propertyInteger
                countObjects                                                                   =  +countObjects                &
                     &                                                                            +size(propertyInteger,dim=1)
             end do
             importer_ => importer_%next
          end do
          call simulations(1)%propertiesInteger%set(self%importers%simulations(1)%propertiesInteger%key(k),propertyIntegerMerged)
          nullify(propertyIntegerMerged)
       end do
    end if
    !! Scalar real properties.
    simulations(1)%propertiesReal=rank1DoublePtrHash()
    if (self%importers%simulations(1)%propertiesReal%size() > 0) then
       do k=1,self%importers%simulations(1)%propertiesReal%size()
          countObjects       = 0_c_size_t
          importer_          => self%importers
          propertyRealMerged => null()
          do while (associated(importer_))
             do j=1,size(importer_%simulations)
                propertyReal                                                             =>  importer_%simulations(j)%propertiesReal%value(k)
                if (.not.associated(propertyRealMerged)) allocate(propertyRealMerged(countObjectsMerged))
                propertyRealMerged(countObjects+1:countObjects+size(propertyReal,dim=1)) =   propertyReal
                countObjects                                                             =  +countObjects                &
                     &                                                                      +size(propertyReal,dim=1)
             end do
             importer_ => importer_%next
          end do
          call simulations(1)%propertiesReal%set(self%importers%simulations(1)%propertiesReal%key(k),propertyRealMerged)
          nullify(propertyRealMerged)
       end do
    end if
    !! Rank-1 integer properties.
    simulations(1)%propertiesIntegerRank1=rank2IntegerSizeTPtrHash()
    if (self%importers%simulations(1)%propertiesIntegerRank1%size() > 0) then
       do k=1,self%importers%simulations(1)%propertiesIntegerRank1%size()
          countObjects               = 0_c_size_t
          importer_                  => self%importers
          propertyIntegerRank1Merged => null()
          do while (associated(importer_))
             do j=1,size(importer_%simulations)
                propertyIntegerRank1                                                                       =>  importer_%simulations(j)%propertiesIntegerRank1%value(k)
                if (.not.associated(propertyIntegerRank1Merged)) allocate(propertyIntegerRank1Merged(size(propertyIntegerRank1,dim=1),countObjectsMerged))
                propertyIntegerRank1Merged(:,countObjects+1:countObjects+size(propertyIntegerRank1,dim=2)) =   propertyIntegerRank1
                countObjects                                                                               =  +countObjects                     &
                     &                                                                                        +size(propertyIntegerRank1,dim=2)
             end do
             importer_ => importer_%next
          end do
          call simulations(1)%propertiesIntegerRank1%set(self%importers%simulations(1)%propertiesIntegerRank1%key(k),propertyIntegerRank1Merged)
          nullify(propertyIntegerRank1Merged)
       end do
    end if
    !! Rank-1 real properties.
    simulations(1)%propertiesRealRank1=rank2DoublePtrHash()
    if (self%importers%simulations(1)%propertiesRealRank1%size() > 0) then
       do k=1,self%importers%simulations(1)%propertiesRealRank1%size()
          countObjects            = 0_c_size_t
          importer_               => self%importers
          propertyRealRank1Merged => null()
          do while (associated(importer_))
             do j=1,size(importer_%simulations)
                propertyRealRank1                                                                    =>  importer_%simulations(j)%propertiesRealRank1%value(k)
                if (.not.associated(propertyRealRank1Merged)) allocate(propertyRealRank1Merged(size(propertyRealRank1,dim=1),countObjectsMerged))
                propertyRealRank1Merged(:,countObjects+1:countObjects+size(propertyRealRank1,dim=2)) =   propertyRealRank1
                countObjects                                                                         =  +countObjects                     &
                     &                                                                                  +size(propertyRealRank1,dim=2)
             end do
             importer_ => importer_%next
          end do
          call simulations(1)%propertiesRealRank1%set(self%importers%simulations(1)%propertiesRealRank1%key(k),propertyRealRank1Merged)
          nullify(propertyRealRank1Merged)
       end do
    end if
    ! Close analysis groups in merged importers and deallocate their simulation data.
    importer_ => self%importers
    do while (associated(importer_))
       deallocate(importer_%simulations)
       importer_ => importer_%next
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine mergeImport

  logical function mergeIsHDF5(self)
    !!{
    Return whether or not the imported data is from an HDF5 file.
    !!}
    implicit none
    class(nbodyImporterMerge), intent(inout) :: self

    mergeIsHDF5=.false.
    return
  end function mergeIsHDF5
