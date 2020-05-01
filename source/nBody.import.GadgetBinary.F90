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

!% Contains a module which implements an N-body data importer for Gadget binary files.

  use :: ISO_Varying_String, only : varying_string

  !# <nbodyImporter name="nbodyImporterGadgetBinary">
  !#  <description>An importer for Gadget binary files.</description>
  !# </nbodyImporter>
  type, extends(nbodyImporterClass) :: nbodyImporterGadgetBinary
     !% An importer for Gadget HDF5 files.
     private
     type            (varying_string) :: outputFileName , fileName
     integer                          :: particleType
     double precision                 :: lengthSoftening, unitMassInSI    , &
          &                              unitLengthInSI , unitVelocityInSI
   contains
     procedure :: import => gadgetBinaryImport
  end type nbodyImporterGadgetBinary

  interface nbodyImporterGadgetBinary
     !% Constructors for the ``gadgetBinary'' N-body importer class.
     module procedure gadgetBinaryConstructorParameters
     module procedure gadgetBinaryConstructorInternal
  end interface nbodyImporterGadgetBinary

contains

  function gadgetBinaryConstructorParameters(parameters) result (self)
    !% Constructor for the ``gadgetBinary'' N-body importer class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyImporterGadgetBinary)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    type            (varying_string           )                :: outputFileName , fileName
    integer                                                    :: particleType
    double precision                                           :: lengthSoftening, unitMassInSI    , &
         &                                                        unitLengthInSI , unitVelocityInSI

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
    !#   <description>The particle type to read from the Gadget bianry file.</description>
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
    !# <inputParameter>
    !#   <name>unitVelocityInSI</name>
    !#   <source>parameters</source>
    !#   <description>The velocity unit expressed in the SI system.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputFileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of the file to which any output analysis should be written.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=nbodyImporterGadgetBinary(fileName,outputFileName,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI)
    !# <inputParametersValidate source="parameters"/>
    return
  end function gadgetBinaryConstructorParameters

  function gadgetBinaryConstructorInternal(fileName,outputFileName,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI) result (self)
    !% Internal constructor for the ``gadgetBinary'' N-body importer class.
    implicit none
    type            (nbodyImporterGadgetBinary)                :: self
    type            (varying_string           ), intent(in   ) :: outputFileName , fileName
    integer                                    , intent(in   ) :: particleType
    double precision                           , intent(in   ) :: lengthSoftening, unitMassInSI   , &
         &                                                        unitLengthInSI ,unitVelocityInSI
    !# <constructorAssign variables="fileName, outputFileName,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI"/>

    return
  end function gadgetBinaryConstructorInternal

  function gadgetBinaryImport(self)
    !% Import data from a Gadget HDF5 file.
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Memory_Management               , only : allocateArray
    use :: Numerical_Constants_Astronomical, only : massSolar              , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    type            (nBodyData                )                                :: gadgetBinaryImport
    class           (nbodyImporterGadgetBinary), intent(inout)                 :: self
    integer                                                   , dimension(6  ) :: numberParticleType    , numberParticleTypeFile, &
         &                                                                        numberParticleTypeRead
    real                                       , allocatable  , dimension(:,:) :: position              , velocity
    double precision                                          , dimension(6  ) :: massParticleType
    double precision                                                           :: time                  , redshift
    integer                                                                    :: file                  , flagSFR               , &
         &                                                                        flagFeedback          , flagCooling           , &
         &                                                                        numberFiles           , fileNumber
    character       (len=6                    )                                :: fileNumberText

    ! Open the file.
    numberFiles           =huge(0)
    fileNumber            =     0
    numberParticleTypeRead=     0
    do while (fileNumber < numberFiles)
       write (fileNumberText,'(i4)') fileNumber
       open(newUnit=file,file=char(self%fileName)//"."//trim(adjustl(fileNumberText)),status='old',form='unformatted')
       read (file) numberParticleTypeFile, massParticleType, time, redshift, flagSFR, flagFeedback, numberParticleType, flagCooling, numberFiles
       if (fileNumber == 0) then
          call allocateArray(gadgetBinaryImport%position  ,[3,numberParticleType(self%particleType+1)])
          call allocateArray(gadgetBinaryImport%velocity  ,[3,numberParticleType(self%particleType+1)])
          call allocateArray(gadgetBinaryImport%identifier,[  numberParticleType(self%particleType+1)])
       end if
       allocate(position(3,numberParticleTypeFile(self%particleType+1)))
       allocate(velocity(3,numberParticleTypeFile(self%particleType+1)))
       read (file) position
       read (file) velocity
       gadgetBinaryImport%position(:,numberParticleTypeRead(self%particleType+1)+1:numberParticleTypeRead(self%particleType+1)+numberParticleTypeFile(self%particleType+1))=position
       gadgetBinaryImport%velocity(:,numberParticleTypeRead(self%particleType+1)+1:numberParticleTypeRead(self%particleType+1)+numberParticleTypeFile(self%particleType+1))=velocity
       deallocate(position)
       deallocate(velocity)
       read (file) gadgetBinaryImport%identifier(numberParticleTypeRead(self%particleType+1)+1:numberParticleTypeRead(self%particleType+1)+numberParticleTypeFile(self%particleType+1))
       close(file)
       numberParticleTypeRead=numberParticleTypeRead+numberParticleTypeFile
       fileNumber            =fileNumber            +1
    end do
    ! Set particle mass.
    gadgetBinaryImport%massParticle=+massParticleType(self%particleType+1) &
         &                          *self%unitMassInSI                     &
         &                          /massSolar
    ! Convert position and velocities to internal units.
    gadgetBinaryImport%position   =+gadgetBinaryImport%position &
         &                         *self%unitLengthInSI         &
         &                         /megaParsec
    gadgetBinaryImport%velocity   =+gadgetBinaryImport%velocity &
         &                         *self%unitVelocityInSI       &
         &                         /kilo
    ! Open the particle group - this group will be used for analysis output.
    call gadgetBinaryImport%analysis%openFile    (char(self%outputFileName)                 )
    call gadgetBinaryImport%analysis%writeDataset(gadgetBinaryImport%identifier,'identifier')
    return
  end function gadgetBinaryImport
