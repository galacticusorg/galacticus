!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements an N-body data importer for Gadget HDF5 files.
!!}

  use :: IO_HDF5, only : hdf5Object

  !![
  <nbodyImporter name="nbodyImporterGadgetHDF5">
   <description>An importer for Gadget HDF5 files.</description>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterGadgetHDF5
     !!{
     An importer for Gadget HDF5 files.
     !!}
     private
     type            (varying_string) :: fileName       , label
     type            (hdf5Object    ) :: file
     integer                          :: particleType
     double precision                 :: lengthSoftening, unitMassInSI    , &
          &                              unitLengthInSI , unitVelocityInSI
   contains
     final     ::           gadgetHDF5Destructor
     procedure :: import => gadgetHDF5Import
     procedure :: isHDF5 => gadgetHDF5IsHDF5
  end type nbodyImporterGadgetHDF5

  interface nbodyImporterGadgetHDF5
     !!{
     Constructors for the ``gadgetHDF5'' N-body importer class.
     !!}
     module procedure gadgetHDF5ConstructorParameters
     module procedure gadgetHDF5ConstructorInternal
  end interface nbodyImporterGadgetHDF5

contains

  function gadgetHDF5ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``gadgetHDF5'' N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyImporterGadgetHDF5)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    integer                                                  :: particleType
    double precision                                         :: lengthSoftening, unitMassInSI    , &
         &                                                      unitLengthInSI , unitVelocityInSI
    type            (varying_string         )                :: fileName       , label

    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file to read.</description>
    </inputParameter>
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <description>A label for the simulation</description>
      <defaultValue>var_str('primary')</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>particleType</name>
      <source>parameters</source>
      <description>The particle type to read from the Gadget HDF5 file.</description>
    </inputParameter>
    <inputParameter>
      <name>lengthSoftening</name>
      <source>parameters</source>
      <description>The softening length.</description>
    </inputParameter>
    <inputParameter>
      <name>unitMassInSI</name>
      <source>parameters</source>
      <description>The mass unit expressed in the SI system.</description>
    </inputParameter>
    <inputParameter>
      <name>unitLengthInSI</name>
      <source>parameters</source>
      <description>The length unit expressed in the SI system.</description>
    </inputParameter>
    <inputParameter>
      <name>unitVelocityInSI</name>
      <source>parameters</source>
      <description>The velocity unit expressed in the SI system.</description>
    </inputParameter>
    !!]
    self=nbodyImporterGadgetHDF5(fileName,label,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gadgetHDF5ConstructorParameters

  function gadgetHDF5ConstructorInternal(fileName,label,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI) result (self)
    !!{
    Internal constructor for the ``gadgetHDF5'' N-body importer class.
    !!}
    implicit none
    type            (nbodyImporterGadgetHDF5)                :: self
    type            (varying_string         ), intent(in   ) :: fileName       , label
    integer                                  , intent(in   ) :: particleType
    double precision                         , intent(in   ) :: lengthSoftening, unitMassInSI   , &
         &                                                      unitLengthInSI ,unitVelocityInSI
    !![
    <constructorAssign variables="fileName, label, particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI"/>
    !!]

    return
  end function gadgetHDF5ConstructorInternal

  subroutine gadgetHDF5Destructor(self)
    !!{
    Destructor for Gadget HF5 importer class.
    !!}
    implicit none
    type(nbodyImporterGadgetHDF5), intent(inout) :: self

    if (self%file%isOpen()) call self%file%close()
    return
  end subroutine gadgetHDF5Destructor

  subroutine gadgetHDF5Import(self,simulations)
    !!{
    Import data from a Gadget HDF5 file.
    !!}
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Hashes                          , only : rank1IntegerSizeTPtrHash, rank2IntegerSizeTPtrHash, rank1DoublePtrHash, rank2DoublePtrHash, &
         &                                          doubleHash
    use :: Numerical_Constants_Astronomical, only : massSolar               , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (nbodyImporterGadgetHDF5), intent(inout)                              :: self
    type            (nBodyData              ), intent(  out), allocatable, dimension(:  ) :: simulations
    double precision                                                     , dimension(6  ) :: massParticleType
    double precision                                        , pointer    , dimension(:,:) :: position         , velocity            , &
         &                                                                                   sampleWeight
    integer         (c_size_t               )               , pointer    , dimension(:,:) :: boundStatus
    integer         (c_size_t               )               , pointer    , dimension(:  ) :: particleID
    integer         (c_size_t               )                                             :: countParticles   , countBootstrapSample
    character       (len=9                  )                                             :: particleGroupName
    type            (hdf5Object             )                                             :: header           , dataset
    double precision                                                                      :: lengthSoftening  , massParticle

    allocate(simulations(1))
    simulations(1)%label=self%label
    ! Open the data file of the current snapshot.
    call self%file%openFile(char(self%fileName),objectsOverwritable=.true.)
    ! Construct the particle type group to read and verify that it exists.
    write (particleGroupName,'(a8,i1)') "PartType",self%particleType
    if (.not.self%file%hasGroup(particleGroupName)) call Galacticus_Error_Report('particle group does not exist'//{introspection:location})
    ! Read values from header.
    header=self%file%openGroup('Header')
    call header%readAttributeStatic('MassTable',massParticleType)
    call header%close()
    ! Compute softening length and particle mass.
    lengthSoftening=+self%lengthSoftening                  &
         &          *self%unitLengthInSI                   &
         &          /megaParsec
    massParticle   =+massParticleType(self%particleType+1) &
         &          *self%unitMassInSI                     &
         &          /massSolar
    ! Open the particle group - this group will be used for analysis output.
    simulations(1)%analysis=self%file%openGroup(particleGroupName)
    ! Import the particle postions, velocities and IDs. Optionally import also
    ! the bound status of particles.
    dataset=simulations(1)%analysis%openDataset('ParticleIDs')
    countParticles=dataset%size(1)
    call dataset%close()
    allocate(particleID(  countParticles))
    allocate(position  (3,countParticles))
    allocate(velocity  (3,countParticles))
    call simulations(1)%analysis%readDatasetStatic('Coordinates',position  )
    call simulations(1)%analysis%readDatasetStatic('Velocities' ,velocity  )
    call simulations(1)%analysis%readDatasetStatic('ParticleIDs',particleID)
    ! Convert position and velocities to internal units.
    position=+position                    &
         &   *self      %unitLengthInSI   &
         &   /megaParsec
    velocity=+velocity                    &
         &   *self      %unitVelocityInSI &
         &   /kilo
    if (simulations(1)%analysis%hasDataset('selfBoundStatus')) then
       dataset=simulations(1)%analysis%openDataset('selfBoundStatus')
       countBootstrapSample=dataset%size(2)
       call dataset%close()
       allocate(boundStatus (countParticles,countBootstrapSample))
       allocate(sampleWeight(countParticles,countBootstrapSample))
       call simulations(1)%analysis%readDatasetStatic('selfBoundStatus',boundStatus )
       call simulations(1)%analysis%readDatasetStatic('weight'         ,sampleWeight)
    end if
    ! Store the data.
    simulations(1)%propertiesInteger     =rank1IntegerSizeTPtrHash()
    simulations(1)%propertiesIntegerRank1=rank2IntegerSizeTPtrHash()
    simulations(1)%propertiesReal        =rank1DoublePtrHash      ()
    simulations(1)%propertiesRealRank1   =rank2DoublePtrHash      ()
    simulations(1)%attributesReal        =doubleHash              ()
    call simulations(1)%propertiesRealRank1%set('position'       ,position       )
    call simulations(1)%propertiesRealRank1%set('velocity'       ,velocity       )
    call simulations(1)%propertiesInteger  %set('particleID'     ,particleID     )
    call simulations(1)%attributesReal     %set('massParticle'   ,massParticle   )
    call simulations(1)%attributesReal     %set('lengthSoftening',lengthSoftening)
    if (simulations(1)%analysis%hasDataset('selfBoundStatus')) then
       call simulations(1)%propertiesIntegerRank1%set('isBound'     ,boundStatus )
       call simulations(1)%propertiesRealRank1   %set('sampleWeight',sampleWeight)
    end if
    return
  end subroutine gadgetHDF5Import

  logical function gadgetHDF5IsHDF5(self)
    !!{
    Return whether or not the imported data is from an HDF5 file.
    !!}
    implicit none
    class(nbodyImporterGadgetHDF5), intent(inout) :: self

    gadgetHDF5IsHDF5=.true.
    return
  end function gadgetHDF5IsHDF5
