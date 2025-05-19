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
Implements an N-body data importer for Gadget binary files.
!!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: ISO_Varying_String  , only : varying_string

  !![
  <nbodyImporter name="nbodyImporterGadgetBinary">
   <description>An importer for Gadget binary files.</description>
   <runTimeFileDependencies paths="fileName"/>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterGadgetBinary
     !!{
     An importer for Gadget HDF5 files.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     type            (varying_string          )          :: fileName                      , label
     integer                                             :: particleType
     double precision                                    :: lengthSoftening               , unitMassInSI    , &
          &                                                 unitLengthInSI                , unitVelocityInSI
     logical                                             :: isCosmological                , setParticleType
   contains
     final     ::           gadgetBinaryDestructor
     procedure :: import => gadgetBinaryImport
     procedure :: isHDF5 => gadgetBinaryIsHDF5
  end type nbodyImporterGadgetBinary

  interface nbodyImporterGadgetBinary
     !!{
     Constructors for the {\normalfont \ttfamily gadgetBinary} N-body importer class.
     !!}
     module procedure gadgetBinaryConstructorParameters
     module procedure gadgetBinaryConstructorInternal
  end interface nbodyImporterGadgetBinary
  
contains

  function gadgetBinaryConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily gadgetBinary} N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyImporterGadgetBinary)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    class           (cosmologyParametersClass ), pointer       :: cosmologyParameters_
    type            (varying_string           )                :: fileName            , label
    integer                                                    :: particleType
    double precision                                           :: lengthSoftening     , unitMassInSI    , &
         &                                                        unitLengthInSI      , unitVelocityInSI
    logical                                                    :: isCosmological      , setParticleType

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
      <description>The particle type to read from the Gadget binary file.</description>
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
    <inputParameter>
      <name>isCosmological</name>
      <source>parameters</source>
      <description>Set to true if this is a cosmological simulation, false otherwise.</description>
    </inputParameter>
    <inputParameter>
      <name>setParticleType</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, particle type values will be set.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=nbodyImporterGadgetBinary(fileName,label,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI,isCosmological,setParticleType,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function gadgetBinaryConstructorParameters

  function gadgetBinaryConstructorInternal(fileName,label,particleType,lengthSoftening,unitMassInSI,unitLengthInSI,unitVelocityInSI,isCosmological,setParticleType,cosmologyParameters_) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily gadgetBinary} N-body importer class.
    !!}
    implicit none
    type            (nbodyImporterGadgetBinary)                        :: self
    type            (varying_string           ), intent(in   )         :: fileName            , label
    integer                                    , intent(in   )         :: particleType
    double precision                           , intent(in   )         :: lengthSoftening     , unitMassInSI    , &
         &                                                                unitLengthInSI      , unitVelocityInSI
    logical                                    , intent(in   )         :: isCosmological      , setParticleType
    class           (cosmologyParametersClass ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="fileName, label, particleType, lengthSoftening, unitMassInSI, unitLengthInSI, unitVelocityInSI, isCosmological, setParticleType, *cosmologyParameters_"/>
    !!]

    return
  end function gadgetBinaryConstructorInternal

  subroutine gadgetBinaryDestructor(self)
    !!{
    Destructor for the Gadget binary N-body importer class.
    !!}
    implicit none
    type(nbodyImporterGadgetBinary), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine gadgetBinaryDestructor

  subroutine gadgetBinaryImport(self,simulations)
    !!{
    Import data from a Gadget HDF5 file.
    !!}
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Display                         , only : displayIndent           , displayUnindent         , verbosityLevelStandard
    use            :: Cosmology_Parameters            , only : hubbleUnitsLittleH
    use            :: File_Utilities                  , only : File_Exists
    use            :: Hashes                          , only : rank1IntegerSizeTPtrHash, rank2IntegerSizeTPtrHash, rank1DoublePtrHash    , rank2DoublePtrHash, &
         &                                                     doubleHash              , varyingStringHash       , integerSizeTHash      , genericHash
    use            :: Numerical_Constants_Astronomical, only : massSolar               , megaParsec
    use            :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (nbodyImporterGadgetBinary), intent(inout)                                :: self
    type            (nBodyData                ), intent(  out), allocatable  , dimension(:  ) :: simulations
    integer                                                                  , dimension(6  ) :: numberParticleType     , numberParticleTypeFile
    integer                                                                                   :: numberParticleTotal    , numberParticleTotalFile, &
         &                                                                                       numberParticleTotalRead, numberParticleStartRead, &
         &                                                                                       numberMassStartRead    , numberTypeAssigned
    real                                                      , allocatable  , dimension(:,:) :: position               , velocity
    double precision                                          , pointer      , dimension(:,:) :: position_              , velocity_
    double precision                                                         , dimension(6  ) :: massParticleType
    real                                                      , allocatable  , dimension(:  ) :: mass
    double precision                                          , pointer      , dimension(:  ) :: mass_
    integer                                                   , allocatable  , dimension(:  ) :: id
    integer         (c_size_t                 )               , pointer      , dimension(:  ) :: id_                    , type_
    double precision                                                                          :: time                   , redshift               , &
         &                                                                                       hubbleConstantLittleH  , boxSize
    integer                                                                                   :: file                   , flagSFR                , &
         &                                                                                       flagFeedback           , flagCooling            , &
         &                                                                                       numberFiles            , fileNumber             , &
         &                                                                                       i                      , j
    character       (len=6                    )                                               :: fileNumberText
    logical                                                                                   :: particleTypeAll        , massesDiffer           , &
         &                                                                                       readMasses

    call displayIndent('import simulation from Gadget binary file',verbosityLevelStandard)
    allocate(simulations(1))
    simulations(1)%label   =self%label
    particleTypeAll        =self%particleType < 0
    numberFiles            =huge(0)
    fileNumber             =     0
    numberParticleTotalRead=     0
    numberParticleTotal    =     0
    massesDiffer           =.false.
    readMasses             =.false.
    if (self%isCosmological) then
       hubbleConstantLittleH=self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    else
       hubbleConstantLittleH=1.0d0
    end if
    do while (fileNumber < numberFiles)
       ! Check for the existence of the named file with no subfile suffix. This can occur, for example, if reading initial conditions files.
       if (fileNumber == 0 .and. File_Exists(self%fileName)) then
          ! Open the file with no subfile suffix.
          open(newUnit=file,file=char(self%fileName)                                    ,status='old',form='unformatted')
       else
          ! Append the appropriate subfile suffix and open the file.
          write (fileNumberText,'(i4)') fileNumber
          open(newUnit=file,file=char(self%fileName)//"."//trim(adjustl(fileNumberText)),status='old',form='unformatted')
       end if
       read (file) numberParticleTypeFile, massParticleType, time, redshift, flagSFR, flagFeedback, numberParticleType, flagCooling, numberFiles, boxSize
       if (fileNumber == 0) then
          if (particleTypeAll) then
             numberParticleTotal=sum(numberParticleType                     )
             massesDiffer       =any(massParticleType                      == 0.0d0 .and. numberParticleType > 0) .or. any(massParticleType /= massParticleType(1))
          else
             numberParticleTotal=    numberParticleType(self%particleType+1)
             massesDiffer       =    massParticleType(self%particleType+1) == 0.0d0 
          end if
          readMasses=massesDiffer .and. any(massParticleType == 0.0d0)
          allocate       (position_(3,numberParticleTotal))
          allocate       (velocity_(3,numberParticleTotal))
          allocate       (id_      (  numberParticleTotal))
          if (     massesDiffer   )                         &
               & allocate(mass_    (  numberParticleTotal))
          if (self%setParticleType) then
             allocate    (type_    (  numberParticleTotal))
          else
             allocate    (type_    (                    0))
          end if
       end if
       if (particleTypeAll) then
          numberParticleTotalFile   =sum(numberParticleTypeFile                       )
          numberParticleStartRead   =0
       else
          numberParticleTotalFile   =    numberParticleTypeFile(  self%particleType+1)
          if (self%particleType > 0) then
             numberParticleStartRead=sum(numberParticleTypeFile(1:self%particleType  ))
          else
             numberParticleStartRead=0
          end if
       end if
       allocate       (position(3,sum(numberParticleTypeFile                               )))
       allocate       (velocity(3,sum(numberParticleTypeFile                               )))
       allocate       (id      (  sum(numberParticleTypeFile                               )))
       if (readMasses)                                                                         &
            & allocate(mass    (  sum(numberParticleTypeFile,mask=massParticleType == 0.0d0)))
       read        (file) position
       read        (file) velocity
       read        (file) id
       if (readMasses .and. size(mass) > 0) &
            & read (file) mass
       close(file)
       position_(:,numberParticleTotalRead+1:numberParticleTotalRead+numberParticleTotalFile)=position(:,numberParticleStartRead+1:numberParticleStartRead+numberParticleTotalFile)
       velocity_(:,numberParticleTotalRead+1:numberParticleTotalRead+numberParticleTotalFile)=velocity(:,numberParticleStartRead+1:numberParticleStartRead+numberParticleTotalFile)
       id_      (  numberParticleTotalRead+1:numberParticleTotalRead+numberParticleTotalFile)=id      (  numberParticleStartRead+1:numberParticleStartRead+numberParticleTotalFile)
       if (massesDiffer) then
          if (particleTypeAll) then
             j                  =0
             numberMassStartRead=0
             do i=1,6
                if (numberParticleTypeFile(i) == 0) cycle
                if (readMasses .and. massParticleType(i) == 0.0d0) then
                   mass_              (numberParticleTotalRead+j+1:numberParticleTotalRead+j+numberParticleTypeFile(i))=                    mass                  (numberMassStartRead+1:numberMassStartRead+numberParticleTypeFile  (i))
                   numberMassStartRead                                                                                 =numberMassStartRead+numberParticleTypeFile(i)
                else
                   mass_              (numberParticleTotalRead+j+1:numberParticleTotalRead+j+numberParticleTypeFile(i))=                    massParticleType      (                                                                   i )
                end if
                j                                                                                                      =j                  +numberParticleTypeFile(                                                                   i )
             end do
          else
             if (readMasses) then
                if (self%particleType > 0) then
                   numberMassStartRead=sum(numberParticleTypeFile(1:self%particleType),mask=massParticleType(1:self%particleType) == 0.0d0)
                else
                   numberMassStartRead=0
                end if
                mass_                 (numberParticleTotalRead  +1:numberParticleTotalRead  +numberParticleTotalFile   )=                    mass                  (numberMassStartRead+1:numberMassStartRead+numberParticleTotalFile   )
             else
                mass_                 (numberParticleTotalRead  +1:numberParticleTotalRead  +numberParticleTotalFile   )=                    massParticleType      (self%particleType  +1                                               )
             end if
          end if
       end if
       if (self%setParticleType) then
          if (particleTypeAll) then
             numberTypeAssigned=0
             do i=1,6
                type_(numberParticleTotalRead+numberTypeAssigned+1:numberParticleTotalRead+numberTypeAssigned+numberParticleTypeFile (i))=i-1
                numberTypeAssigned=numberTypeAssigned+numberParticleTypeFile(i)
             end do
          else
             type_   (numberParticleTotalRead                   +1:numberParticleTotalRead                   +numberParticleTotalFile   )=self%particleType
          end if
       end if
       deallocate     (position)
       deallocate     (velocity)
       deallocate     (id      )
       if (readMasses) &
            deallocate(mass    )
       if (particleTypeAll) then
          numberParticleTotalRead=numberParticleTotalRead+sum(numberParticleTypeFile                     )
       else
          numberParticleTotalRead=numberParticleTotalRead+    numberParticleTypeFile(self%particleType+1)
       end if
       fileNumber=fileNumber+1
    end do
    ! Convert position and velocities to internal units.
    position_=+position_                   &
         &    *self      %unitLengthInSI   &
         &    /hubbleConstantLittleH       &
         &    /megaParsec
    velocity_=+velocity_                   &
         &    *self      %unitVelocityInSI &
         &    /kilo
    boxSize  =+boxSize                     &
         &    *self      %unitLengthInSI   &
         &    /hubbleConstantLittleH       &
         &    /megaParsec
    !! For cosmological simulations velocities are in internal Gadget units and must be multiplied by √a to get peculiar velocity.
    if (self%isCosmological)         &
         velocity_=+velocity_        &
         &         *sqrt(            &
         &               +  1.0d0    &
         &               /(          &
         &                 +1.0d0    &
         &                 +redshift &
         &                )          &
         &              )
    ! Set the properties.
    simulations(1)%propertiesInteger     =rank1IntegerSizeTPtrHash()
    simulations(1)%propertiesIntegerRank1=rank2IntegerSizeTPtrHash()
    simulations(1)%propertiesReal        =rank1DoublePtrHash      ()
    simulations(1)%propertiesRealRank1   =rank2DoublePtrHash      ()
    simulations(1)%attributesReal        =doubleHash              ()
    simulations(1)%attributesText        =varyingStringHash       ()
    simulations(1)%propertiesRealRank1   =rank2DoublePtrHash      ()
    simulations(1)%attributesGeneric     =genericHash             ()
    call simulations(1)%attributesReal     %set('boxSize'     ,boxSize  )
    call simulations(1)%attributesReal     %set('redshift'    ,redshift )
    call simulations(1)%propertiesRealRank1%set('position'    ,position_)
    call simulations(1)%propertiesRealRank1%set('velocity'    ,velocity_)
    call simulations(1)%propertiesInteger  %set('particleID'  ,id_      )
    call simulations(1)%propertiesInteger  %set('particleType',type_    )
    ! Set particle mass.
    if (massesDiffer) then
       mass_  =+mass_                  &
            &  *self     %unitMassInSI &
            &  /hubbleConstantLittleH  &
            &  /massSolar
       call simulations(1)%propertiesReal%set(                                        &
            &                                 'massParticle'                        , &
            &                                  mass_                                  &
            &                                )
       call simulations(1)%attributesReal%set(                                        &
            &                                 'massParticle'                        , &
            &                                 -1.0d0                                  &
            &                                )
    else
       call simulations(1)%attributesReal%set(                                        &
            &                                 'massParticle'                        , &
            &                                 +massParticleType(self%particleType+1)  &
            &                                 *self%unitMassInSI                      &
            &                                 /hubbleConstantLittleH                  &
            &                                 /massSolar                              &
            &                                )
    end if
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine gadgetBinaryImport

  logical function gadgetBinaryIsHDF5(self)
    !!{
    Return whether or not the imported data is from an HDF5 file.
    !!}
    implicit none
    class(nbodyImporterGadgetBinary), intent(inout) :: self

    gadgetBinaryIsHDF5=.false.
    return
  end function gadgetBinaryIsHDF5
