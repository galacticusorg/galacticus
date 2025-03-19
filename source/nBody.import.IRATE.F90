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
Implements an N-body data importer for IRATE files.
!!}

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: IO_HDF5             , only : hdf5Object

  !![
  <nbodyImporter name="nbodyImporterIRATE">
   <description>An importer for IRATE files.</description>
   <runTimeFileDependencies paths="fileName"/>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterIRATE
     !!{
     An importer for IRATE files.
     !!}
     private
     class  (cosmologyParametersClass), pointer                   :: cosmologyParameters_ => null()
     class  (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_  => null()
     type   (varying_string          )                            :: fileName                      , label
     type   (hdf5Object              )                            :: file
     integer                                                      :: snapshot
     logical                                                      :: haveProperties
     type   (varying_string          ), allocatable, dimension(:) :: properties
   contains
     final     ::           irateDestructor
     procedure :: import => irateImport
     procedure :: isHDF5 => irateIsHDF5
  end type nbodyImporterIRATE

  interface nbodyImporterIRATE
     !!{
     Constructors for the {\normalfont \ttfamily irate} N-body importer class.
     !!}
     module procedure irateConstructorParameters
     module procedure irateConstructorInternal
  end interface nbodyImporterIRATE

contains

  function irateConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily irate} N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyImporterIRATE      )                            :: self
    type   (inputParameters         ), intent(inout)             :: parameters
    class  (cosmologyParametersClass), pointer                   :: cosmologyParameters_
    class  (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_
    type   (varying_string          ), allocatable, dimension(:) :: properties
    type   (varying_string          )                            :: fileName            , label
    integer                                                      :: snapshot

    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file to read.</description>
    </inputParameter>
    <inputParameter>
      <name>snapshot</name>
      <source>parameters</source>
      <description>The snapshot number to read.</description>
    </inputParameter>
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <defaultValue>var_str('primary')</defaultValue>
      <description>A label for the simulation.</description>
    </inputParameter>
    !!]
    allocate(properties(parameters%count('properties',zeroIfNotPresent=.true.)))
    if (size(properties) > 0) then
       !![
       <inputParameter>
         <name>properties</name>
         <source>parameters</source>
         <description>Properties to read from the simulation.</description>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=nbodyImporterIRATE(fileName,label,snapshot,properties,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function irateConstructorParameters

  function irateConstructorInternal(fileName,label,snapshot,properties,cosmologyParameters_,cosmologyFunctions_) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily irate} N-body importer class.
    !!}
    implicit none
    type   (nbodyImporterIRATE      )                              :: self
    type   (varying_string          ), intent(in   )               :: fileName            , label
    integer                          , intent(in   )               :: snapshot
    type   (varying_string          ), intent(in   ), dimension(:) :: properties
    class  (cosmologyParametersClass), intent(in   ), target       :: cosmologyParameters_
    class  (cosmologyFunctionsClass ), intent(in   ), target       :: cosmologyFunctions_
    !![
    <constructorAssign variables="fileName, label, snapshot, properties, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    self%haveProperties=size(properties) > 0
    return
  end function irateConstructorInternal

  subroutine irateDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily irate} importer class.
    !!}
    implicit none
    type(nbodyImporterIRATE), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    if (self%file%isOpen()) call self%file%close()
    return
  end subroutine irateDestructor

  subroutine irateImport(self,simulations)
    !!{
    Import data from a IRATE file.
    !!}
    use :: Display         , only : displayIndent     , displayUnindent         , verbosityLevelStandard
    use :: Error, only : errorStatusSuccess
    use :: Hashes          , only : doubleHash        , integerSizeTHash        , rank1DoublePtrHash    , rank1IntegerSizeTPtrHash, &
          &                         rank2DoublePtrHash, rank2IntegerSizeTPtrHash, varyingStringHash     , genericHash
    use :: HDF5_Access     , only : hdf5Access
    use :: IO_HDF5         , only : H5T_NATIVE_DOUBLES, H5T_NATIVE_INTEGERS     , hdf5Object
    use :: IO_IRATE        , only : irate
    implicit none
    class           (nbodyImporterIRATE), intent(inout)                              :: self
    type            (nBodyData         ), intent(  out), dimension(:  ), allocatable :: simulations
    type            (varying_string    )               , dimension(:  ), allocatable :: datasetNames
    integer         (c_size_t          )               , dimension(:  ), pointer     :: propertyInteger, particleIDs
    double precision                                   , dimension(:  ), pointer     :: propertyReal
    double precision                                   , dimension(:,:), pointer     :: position       , velocity
    type            (irate             )                                             :: irate_
    type            (hdf5Object        )                                             :: snapshotGroup  , dataset
    character       (len=13            )                                             :: snapshotLabel
    integer                                                                          :: i              , status
    double precision                                                                 :: boxSize

    call displayIndent('import simulation from IRATE file',verbosityLevelStandard)
    allocate(simulations(1))
    simulations(1)%label                 =self%label
    simulations(1)%propertiesInteger     =rank1IntegerSizeTPtrHash()
    simulations(1)%propertiesReal        =rank1DoublePtrHash      ()
    simulations(1)%propertiesIntegerRank1=rank2IntegerSizeTPtrHash()
    simulations(1)%propertiesRealRank1   =rank2DoublePtrHash      ()
    simulations(1)%attributesInteger     =integerSizeTHash        ()
    simulations(1)%attributesReal        =doubleHash              ()
    simulations(1)%attributesText        =varyingStringHash       ()
    simulations(1)%attributesGeneric     =genericHash             ()
    irate_=irate(char(self%fileName),self%cosmologyParameters_,self%cosmologyFunctions_)
    call irate_        %readSimulation    (          boxSize)
    call simulations(1)%attributesReal%set('boxSize',boxSize)
    !![
    <conditionalCall>
     <call>call irate_%readHalos(self%snapshot{conditions})</call>
     <argument name="center"   value="position"    condition=".not.self%haveProperties .or. any(self%properties == 'position'  )"/>
     <argument name="velocity" value="velocity"    condition=".not.self%haveProperties .or. any(self%properties == 'velocity'  )"/>
     <argument name="IDs"      value="particleIDs" condition=".not.self%haveProperties .or. any(self%properties == 'particleID')"/>
    </conditionalCall>
    !!]
    if (.not.self%haveProperties .or. any(self%properties == 'particleID')) call simulations(1)%propertiesInteger  %set('particleID',particleIDs)
    if (.not.self%haveProperties .or. any(self%properties == 'position'  )) call simulations(1)%propertiesRealRank1%set('position'  ,position   )
    if (.not.self%haveProperties .or. any(self%properties == 'velocity'  )) call simulations(1)%propertiesRealRank1%set('velocity'  ,velocity   )
    write (snapshotLabel,'(a,i5.5)') 'Snapshot',self%snapshot
    !$ call hdf5Access%set()
    call self%file%openFile(char(self%fileName),readOnly=.false.,objectsOverwritable=.true.)
    snapshotGroup            =self%file         %openGroup(snapshotLabel)
    simulations  (1)%analysis=     snapshotGroup%openGroup('HaloCatalog')
    call simulations(1)%analysis%datasets(datasetNames)
    do i=1,size(datasetNames)
       if     (                               &
            &   datasetNames(i) == "Center"   &
            &  .or.                           &
            &   datasetNames(i) == "Velocity" &
            &  .or.                           &
            &   datasetNames(i) == "HaloID"   &
            & ) cycle
       if (self%haveProperties .and. .not.any(self%properties == datasetNames(i))) cycle
       dataset=simulations(1)%analysis%openDataset(char(datasetNames(i)))
       call dataset%assertDatasetType(H5T_NATIVE_INTEGERS,1,status)
       if (status == errorStatusSuccess) then
          allocate(propertyInteger(dataset%size(1)))
          call dataset       %readDatasetStatic    (                datasetValue=propertyInteger)
          call simulations(1)%propertiesInteger%set(datasetNames(i),             propertyInteger)
          nullify(propertyInteger)
       end if
       call dataset%assertDatasetType(H5T_NATIVE_DOUBLES ,1,status)
       if (status == errorStatusSuccess) then
          allocate(propertyReal   (dataset%size(1)))
          call dataset       %readDatasetStatic    (                datasetValue=propertyReal   )
          call simulations(1)%propertiesReal   %set(datasetNames(i),             propertyReal   )
          nullify(propertyReal   )
       end if
       call dataset%close()
    end do
    call snapshotGroup%close()
    !$ call hdf5Access%unset()
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine irateImport

  logical function irateIsHDF5(self)
    !!{
    Return whether or not the imported data is from an HDF5 file.
    !!}
    implicit none
    class(nbodyImporterIRATE), intent(inout) :: self

    irateIsHDF5=.true.
    return
  end function irateIsHDF5
