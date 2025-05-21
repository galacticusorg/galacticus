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
Implements an N-body data importer for Millennium database CSV files.
!!}

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass
  
  !![
  <nbodyImporter name="nbodyImporterMillenniumCSV">
   <description>An importer for Millennium database CSV files.</description>
   <runTimeFileDependencies paths="fileName"/>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterMillenniumCSV
     !!{
     An importer for Millennium database CSV files.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     type            (varying_string          )          :: fileName                      , label
     double precision                                    :: time                          , redshift
   contains
     final     ::           millenniumCSVDestructor
     procedure :: import => millenniumCSVImport
     procedure :: isHDF5 => millenniumCSVIsHDF5
  end type nbodyImporterMillenniumCSV

  interface nbodyImporterMillenniumCSV
     !!{
     Constructors for the \refClass{nbodyImporterMillenniumCSV} N-body importer class.
     !!}
     module procedure millenniumCSVConstructorParameters
     module procedure millenniumCSVConstructorInternal
  end interface nbodyImporterMillenniumCSV

  type columnType
     !!{
     Type used to describe column types.
     !!}
     logical                 :: isPosition=.false., isVelocity=.false., &
          &                     isReal    =.false., isInteger =.false., &
          &                     isID
     integer                 :: index
     type   (varying_string) :: label
  end type columnType
  
contains

  function millenniumCSVConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyImporterMillenniumCSV} N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyImporterMillenniumCSV)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    type            (varying_string            )                :: fileName            , label
    double precision                                            :: redshift            , time
    
    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file to read.</description>
    </inputParameter>
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <description>A label for the simulation.</description>
      <defaultValue>var_str('primary')</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <description>The redshift of the data.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    self=nbodyImporterMillenniumCSV(fileName,label,time,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function millenniumCSVConstructorParameters

  function millenniumCSVConstructorInternal(fileName,label,time,cosmologyParameters_,cosmologyFunctions_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyImporterMillenniumCSV} N-body importer class.
    !!}
    implicit none
    type            (nbodyImporterMillenniumCSV)                        :: self
    class           (cosmologyParametersClass  ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), intent(in   ), target :: cosmologyFunctions_
    type            (varying_string            ), intent(in   )         :: fileName            , label
    double precision                            , intent(in   )         :: time
    !![
    <constructorAssign variables="fileName, label, time, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    self%redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    return
  end function millenniumCSVConstructorInternal

  subroutine millenniumCSVDestructor(self)
    !!{
    Destructor for the \refClass{nbodyImporterMillenniumCSV} importer class.
    !!}
    implicit none
    type(nbodyImporterMillenniumCSV), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine millenniumCSVDestructor

  subroutine millenniumCSVImport(self,simulations)
    !!{
    Import data from a MillenniumCSV file.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Display             , only : displayCounter        , displayCounterClear     , displayIndent     , displayUnindent         , &
          &                             verbosityLevelStandard
    use :: Error               , only : Error_Report
    use :: File_Utilities      , only : Count_Lines_in_File
    use :: Hashes              , only : rank1DoublePtrHash    , rank1IntegerSizeTPtrHash, rank2DoublePtrHash, rank2IntegerSizeTPtrHash, &
         &                              doubleHash            , integerSizeTHash        , varyingStringHash , genericHash
    use :: String_Handling     , only : String_Count_Words    , String_Split_Words
    implicit none
    class           (nbodyImporterMillenniumCSV), intent(inout)                              :: self
    type            (nBodyData                 ), intent(  out), dimension(  :), allocatable :: simulations
    type            (columnType                )               , dimension(  :), allocatable :: columns
    double precision                                           , dimension(:,:), pointer     :: position                , velocity
    type            (nbodyPropertiesRealList   )               , dimension(  :), allocatable :: propertiesReal
    type            (nbodyPropertiesIntegerList)               , dimension(  :), allocatable :: propertiesInteger
    integer         (c_size_t                  )               , dimension(  :), pointer     :: particleID
    character       (len=  32                  )               , dimension(  :), allocatable :: columnsText
    double precision                            , parameter                                  :: massUnit         =1.0d10
    integer         (c_size_t                  )                                             :: countPoints             , i
    integer                                                                                  :: status                  , file        , &
         &                                                                                      j                       , countColumns, &
         &                                                                                      countReal               , countInteger
    character       (len=1024                  )                                             :: line
    logical                                                                                  :: isComment               , readHeader  , &
         &                                                                                      gotPositionX            , gotPositionY, &
         &                                                                                      gotPositionZ            , gotVelocityX, &
         &                                                                                      gotVelocityY            , gotVelocityZ, &
         &                                                                                      gotID
    !$GLC attributes initialized :: columns

    call displayIndent('import simulation from MillenniumCSV file',verbosityLevelStandard)
    allocate(simulations(1))
    simulations(1)%label=self%label
    ! Count lines in file. (Subtract 1 since one lines gives the column headers.)
    countPoints=Count_Lines_In_File(self%fileName,comment_char="#")-1_c_size_t
    ! Allocate storage
    allocate(position  (3_c_size_t,countPoints))
    allocate(velocity  (3_c_size_t,countPoints))
    allocate(particleID(           countPoints))
    allocate(propertiesReal   (0))
    allocate(propertiesInteger(0))
    ! Initialize status.
    gotPositionX=.false.
    gotPositionY=.false.
    gotPositionZ=.false.
    gotVelocityX=.false.
    gotVelocityY=.false.
    gotVelocityZ=.false.
    gotID       =.false.
    ! Open the file and read.
    readHeader=.false.
    i=0_c_size_t
    open(newUnit=file,file=char(self%fileName),status='old',form='formatted',iostat=status)
    do while (status == 0)
       read (file,'(a)',iostat=status) line
       if (status /= 0) exit
       isComment=line(1:1) == "#"
       if (.not.isComment) then
          if (.not.readHeader) then
             countColumns=String_Count_Words(line,",")
             allocate(columns    (countColumns))
             allocate(columnsText(countColumns))
             call String_Split_Words(columnsText,line,",")
             countReal   =0
             countInteger=0
             do j=1,countColumns
                select case (trim(columnsText(j)))
                case ('id','haloID')
                   columns     (j)%isID      =.true.
                   gotID                     =.true.
                case ('x')
                   columns     (j)%isPosition=.true.
                   columns     (j)%index     =1
                   gotPositionX=.true.
                case ('y')
                   columns     (j)%isPosition=.true.
                   columns     (j)%index     =2
                   gotPositionY=.true.
                case ('z')
                   columns     (j)%isPosition=.true.
                   columns     (j)%index     =3
                   gotPositionZ=.true.
                case ('vx','velX')
                   columns     (j)%isVelocity=.true.
                   columns     (j)%index     =1
                   gotVelocityX=.true.
                case ('vy','velY')
                   columns     (j)%isVelocity=.true.
                   columns     (j)%index     =2
                   gotVelocityY=.true.
                case ('vz','velZ')
                   columns     (j)%isVelocity=.true.
                   columns     (j)%index     =3
                   gotVelocityZ=.true.
                case ('m_TopHat')
                   countReal             =countReal+1
                   columns     (j)%isReal=.true.
                   columns     (j)%label ='massVirial'
                   columns     (j)%index =countReal
                end select
             end do
             if (.not.gotID       ) call Error_Report('particle ID is missing'//{introspection:location})
             if (.not.gotPositionX) call Error_Report('x position is missing' //{introspection:location})
             if (.not.gotPositionY) call Error_Report('y position is missing' //{introspection:location})
             if (.not.gotPositionZ) call Error_Report('z position is missing' //{introspection:location})
             if (.not.gotVelocityX) call Error_Report('x velocity is missing' //{introspection:location})
             if (.not.gotVelocityY) call Error_Report('y velocity is missing' //{introspection:location})
             if (.not.gotVelocityZ) call Error_Report('z velocity is missing' //{introspection:location})
             deallocate(propertiesReal   )
             deallocate(propertiesInteger)
             allocate(propertiesReal   (countReal   ))
             allocate(propertiesInteger(countInteger))
             do j=1,countReal
                allocate(propertiesReal   (j)%property(countPoints))
             end do
             do j=1,countInteger
                allocate(propertiesInteger(j)%property(countPoints))
             end do
             readHeader=.true.
          else
             ! Parse the line.
             i=i+1_c_size_t
             call String_Split_Words(columnsText,line,",")
             do j=1,countColumns
                if (columns(j)%isID) then                
                   read (columnsText(j),*) particleID       (                           i)
                else if (columns(j)%isPosition) then                
                   read (columnsText(j),*) position         (columns(j)%index ,         i)
                else if (columns(j)%isVelocity) then                
                   read (columnsText(j),*) velocity         (columns(j)%index ,         i)
                else  if (columns(j)%isReal) then
                   read (columnsText(j),*) propertiesReal   (columns(j)%index)%property(i)
                else  if (columns(j)%isInteger) then
                   read (columnsText(j),*) propertiesInteger(columns(j)%index)%property(i)
                end if
             end do
          end if
       end if
       call displayCounter(int(100.0d0*dble(i)/dble(countPoints)),verbosity=verbosityLevelStandard,isNew=i == 1_c_size_t)
    end do
    call displayCounterClear(verbosityLevelStandard)
    close(file)
    ! Convert position to internal units (physical Mpc).
    do i=1_c_size_t,countPoints
       position(:,i)=+position                                     (:,i               ) &
            &        *self    %cosmologyFunctions_ %expansionFactor(self%time         ) &
            &        /self    %cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH)
    end do
    ! Set positions, velocities, and particleIDs.
    simulations(1)%propertiesInteger     =rank1IntegerSizeTPtrHash()
    simulations(1)%propertiesIntegerRank1=rank2IntegerSizeTPtrHash()
    simulations(1)%propertiesReal        =rank1DoublePtrHash      ()
    simulations(1)%propertiesRealRank1   =rank2DoublePtrHash      ()
    simulations(1)%attributesInteger     =integerSizeTHash        ()
    simulations(1)%attributesReal        =doubleHash              ()
    simulations(1)%attributesText        =varyingStringHash       ()
    simulations(1)%attributesGeneric     =genericHash             ()
    call simulations(1)%propertiesRealRank1%set('position'  ,position  )
    call simulations(1)%propertiesRealRank1%set('velocity'  ,velocity  )
    call simulations(1)%propertiesInteger  %set('particleID',particleID)
    ! Convert and set other properties.
    do j=1,countColumns
       if (columns(j)%isReal) then
          select case (char(columns(j)%label))
          case ('massVirial')
             propertiesReal(columns(j)%index)%property=+propertiesReal(columns(j)%index)%property                    &
                  &                                    *massUnit                                                     &
                  &                                    /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
          end select
          call simulations(1)%propertiesReal    %set(columns(j)%label,propertiesReal   (columns(j)%index)%property)
       else if (columns(j)%isInteger) then
           call simulations(1)%propertiesInteger%set(columns(j)%label,propertiesInteger(columns(j)%index)%property)
      end if
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine millenniumCSVImport

  logical function millenniumCSVIsHDF5(self)
    !!{
    Return whether or not the imported data is from an HDF5 file.
    !!}
    implicit none
    class(nbodyImporterMillenniumCSV), intent(inout) :: self

    millenniumCSVIsHDF5=.false.
    return
  end function millenniumCSVIsHDF5
