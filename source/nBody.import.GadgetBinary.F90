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
     type            (varying_string) :: fileName       , label
     integer                          :: particleType
     double precision                 :: lengthSoftening, unitMassInSI    , &
          &                              unitLengthInSI , unitVelocityInSI
   contains
     procedure :: import => gadgetBinaryImport
     procedure :: isHDF5 => gadgetBinaryIsHDF5
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
    type            (varying_string           )                :: fileName       , label
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
    !#   <name>label</name>
    !#   <source>parameters</source>
    !#   <description>A label for the simulation</description>
    !#   <defaultValue>var_str('primary')</defaultValue>
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
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>unitMassInSI</name>
    !#   <source>parameters</source>
    !#   <description>The mass unit expressed in the SI system.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>unitLengthInSI</name>
    !#   <source>parameters</source>
    !#   <description>The length unit expressed in the SI system.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>unitVelocityInSI</name>
    !#   <source>parameters</source>
    !#   <description>The velocity unit expressed in the SI system.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>unitVelocityInSI</name>
    !#   <source>parameters</source>
    !#   <description>The velocity unit expressed in the SI system.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=nbodyImporterGadgetBinary(fileName,label,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI)
    !# <inputParametersValidate source="parameters"/>
    return
  end function gadgetBinaryConstructorParameters

  function gadgetBinaryConstructorInternal(fileName,label,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI) result (self)
    !% Internal constructor for the ``gadgetBinary'' N-body importer class.
    implicit none
    type            (nbodyImporterGadgetBinary)                :: self
    type            (varying_string           ), intent(in   ) :: fileName       , label
    integer                                    , intent(in   ) :: particleType
    double precision                           , intent(in   ) :: lengthSoftening, unitMassInSI   , &
         &                                                        unitLengthInSI ,unitVelocityInSI
    !# <constructorAssign variables="fileName, label, particleType, lengthSoftening, unitMassInSI, unitLengthInSI, unitVelocityInSI"/>

    return
  end function gadgetBinaryConstructorInternal

  subroutine gadgetBinaryImport(self,simulations)
    !% Import data from a Gadget HDF5 file.
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Memory_Management               , only : allocateArray
    use :: Numerical_Constants_Astronomical, only : massSolar              , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (nbodyImporterGadgetBinary), intent(inout)                                :: self
    type            (nBodyData                ), intent(  out), allocatable  , dimension(:  ) :: simulations
    integer                                                                  , dimension(6  ) :: numberParticleType    , numberParticleTypeFile, &
         &                                                                                       numberParticleTypeRead
    real                                                      , allocatable  , dimension(:,:) :: position              , velocity
    double precision                                                         , dimension(6  ) :: massParticleType
    double precision                                                                          :: time                  , redshift
    integer                                                                                   :: file                  , flagSFR               , &
         &                                                                                       flagFeedback          , flagCooling           , &
         &                                                                                       numberFiles           , fileNumber
    character       (len=6                    )                                               :: fileNumberText

    ! Open the file.
    allocate(simulations(1))
    simulations(1)%label  =self%label
    numberFiles           =huge(0)
    fileNumber            =     0
    numberParticleTypeRead=     0
    do while (fileNumber < numberFiles)
       write (fileNumberText,'(i4)') fileNumber
       open(newUnit=file,file=char(self%fileName)//"."//trim(adjustl(fileNumberText)),status='old',form='unformatted')
       read (file) numberParticleTypeFile, massParticleType, time, redshift, flagSFR, flagFeedback, numberParticleType, flagCooling, numberFiles
       if (fileNumber == 0) then
          call allocateArray(simulations(1)%position,[3,numberParticleType(self%particleType+1)])
          call allocateArray(simulations(1)%velocity,[3,numberParticleType(self%particleType+1)])
       end if
       allocate(position(3,numberParticleTypeFile(self%particleType+1)))
       allocate(velocity(3,numberParticleTypeFile(self%particleType+1)))
       read (file) position
       read (file) velocity
       close(file)
       simulations(1)%position(:,numberParticleTypeRead(self%particleType+1)+1:numberParticleTypeRead(self%particleType+1)+numberParticleTypeFile(self%particleType+1))=position
       simulations(1)%velocity(:,numberParticleTypeRead(self%particleType+1)+1:numberParticleTypeRead(self%particleType+1)+numberParticleTypeFile(self%particleType+1))=velocity
       deallocate(position)
       deallocate(velocity)
       numberParticleTypeRead=numberParticleTypeRead+numberParticleTypeFile
       fileNumber            =fileNumber            +1
    end do
    ! Set particle mass.
    simulations(1)%massParticle=+massParticleType(self%particleType+1) &
         &                      *self%unitMassInSI                     &
         &                      /massSolar
    ! Convert position and velocities to internal units.
    simulations(1)%position   =+simulations(1)%position         &
         &                     *self         %unitLengthInSI    &
         &                     /megaParsec
    simulations(1)%velocity   =+simulations(1)%velocity         &
         &                     *self          %unitVelocityInSI &
         &                     /kilo
    return
  end subroutine gadgetBinaryImport

  logical function gadgetBinaryIsHDF5(self)
    !% Return whether or not the imported data is from an HDF5 file.
    implicit none
    class(nbodyImporterGadgetBinary), intent(inout) :: self

    gadgetBinaryIsHDF5=.false.
    return
  end function gadgetBinaryIsHDF5
