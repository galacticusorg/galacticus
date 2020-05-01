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

!% Contains a module which implements an N-body data importer for Gadget HDF5 files.

  use :: IO_HDF5, only : hdf5Object

  !# <nbodyImporter name="nbodyImporterGadgetHDF5">
  !#  <description>An importer for Gadget HDF5 files.</description>
  !# </nbodyImporter>
  type, extends(nbodyImporterClass) :: nbodyImporterGadgetHDF5
     !% An importer for Gadget HDF5 files.
     private
     type            (varying_string) :: fileName
     type            (hdf5Object    ) :: file
     integer                          :: particleType
     double precision                 :: lengthSoftening, unitMassInSI    , &
          &                              unitLengthInSI , unitVelocityInSI
   contains
     final     ::           gadgetHDF5Destructor
     procedure :: import => gadgetHDF5Import
  end type nbodyImporterGadgetHDF5

  interface nbodyImporterGadgetHDF5
     !% Constructors for the ``gadgetHDF5'' N-body importer class.
     module procedure gadgetHDF5ConstructorParameters
     module procedure gadgetHDF5ConstructorInternal
  end interface nbodyImporterGadgetHDF5

contains

  function gadgetHDF5ConstructorParameters(parameters) result (self)
    !% Constructor for the ``gadgetHDF5'' N-body importer class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyImporterGadgetHDF5)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    integer                                                  :: particleType
    double precision                                         :: lengthSoftening, unitMassInSI    , &
         &                                                      unitLengthInSI , unitVelocityInSI
    type            (varying_string         )                :: fileName

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of the file to read.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>particleType</name>
    !#   <source>parameters</source>
    !#   <description>The particle type to read from the Gadget HDF5 file.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>lengthSoftening</name>
    !#   <source>parameters</source>
    !#   <description>The softening length.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>unitMassInSI</name>
    !#   <source>parameters</source>
    !#   <description>The mass unit expressed in the SI system.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>unitLengthInSI</name>
    !#   <source>parameters</source>
    !#   <description>The length unit expressed in the SI system.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>unitVelocityInSI</name>
    !#   <source>parameters</source>
    !#   <description>The velocity unit expressed in the SI system.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=nbodyImporterGadgetHDF5(fileName,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI)
    !# <inputParametersValidate source="parameters"/>
    return
  end function gadgetHDF5ConstructorParameters

  function gadgetHDF5ConstructorInternal(fileName,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI) result (self)
    !% Internal constructor for the ``gadgetHDF5'' N-body importer class.
    implicit none
    type            (nbodyImporterGadgetHDF5)                :: self
    type            (varying_string         ), intent(in   ) :: fileName
    integer                                  , intent(in   ) :: particleType
    double precision                         , intent(in   ) :: lengthSoftening, unitMassInSI   , &
         &                                                      unitLengthInSI ,unitVelocityInSI
    !# <constructorAssign variables="fileName, particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI"/>

    return
  end function gadgetHDF5ConstructorInternal

  subroutine gadgetHDF5Destructor(self)
    !% Destructor for Gadget HF5 importer class.
    implicit none
    type(nbodyImporterGadgetHDF5), intent(inout) :: self

    if (self%file%isOpen()) call self%file%close()
    return
  end subroutine gadgetHDF5Destructor

  function gadgetHDF5Import(self)
    !% Import data from a Gadget HDF5 file.
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Numerical_Constants_Astronomical, only : massSolar              , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    type            (nBodyData              )                :: gadgetHDF5Import
    class           (nbodyImporterGadgetHDF5), intent(inout) :: self
    double precision                         , dimension(6)  :: massParticleType
    character       (len=9                  )                :: particleGroupName
    type            (hdf5Object             )                :: header

    ! Open the data file of the current snapshot.
    call self%file%openFile(char(self%fileName),objectsOverwritable=.true.)
    ! Construct the particle type group to read and verify that it exists.
    write (particleGroupName,'(a8,i1)') "PartType",self%particleType
    if (.not.self%file%hasGroup(particleGroupName)) call Galacticus_Error_Report('particle group does not exist'//{introspection:location})
    ! Read values from header.
    header=self%file%openGroup('Header')
    call header%readAttributeStatic('MassTable',massParticleType)
    call header%close()
    ! Set softening length and particle mass.
    gadgetHDF5Import%lengthSoftening=+self%lengthSoftening                  &
         &                           *self%unitLengthInSI                   &
         &                           /megaParsec
    gadgetHDF5Import%massParticle   =+massParticleType(self%particleType+1) &
         &                           *self%unitMassInSI                     &
         &                           /massSolar
    ! Open the particle group - this group will be used for analysis output.
    gadgetHDF5Import%analysis=self%file%openGroup(particleGroupName)
    ! Import the particle postions, velocities and IDs.
    call gadgetHDF5Import%analysis%readDataset('Coordinates',gadgetHDF5Import%position   )
    call gadgetHDF5Import%analysis%readDataset('Velocities' ,gadgetHDF5Import%velocity   )
    call gadgetHDF5Import%analysis%readDataset('ParticleIDs',gadgetHDF5Import%particleIDs)
    ! Convert position and velocities to internal units.
    gadgetHDF5Import%position=+gadgetHDF5Import%position &
         &                    *self%unitLengthInSI       &
         &                    /megaParsec
    gadgetHDF5Import%velocity=+gadgetHDF5Import%velocity &
         &                    *self%unitVelocityInSI     &
         &                    /kilo
    return
  end function gadgetHDF5Import
