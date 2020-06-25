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

!% Contains a module which implements an N-body data importer for Rockstar files.

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  
  !# <nbodyImporter name="nbodyImporterRockstar">
  !#  <description>An importer for Rockstar files.</description>
  !# </nbodyImporter>
  type, extends(nbodyImporterClass) :: nbodyImporterRockstar
     !% An importer for Rockstar files.
     private
     class  (cosmologyParametersClass)              , pointer     :: cosmologyParameters_    => null()
     integer                          , dimension(:), allocatable :: readColumns                      , readColumnsType
     integer                                                      :: readColumnsIntegerCount          , readColumnsRealCount
     type   (varying_string          )                            :: fileName                         , label
   contains
     final     ::           rockstarDestructor
     procedure :: import => rockstarImport
     procedure :: isHDF5 => rockstarIsHDF5
  end type nbodyImporterRockstar

  interface nbodyImporterRockstar
     !% Constructors for the {\normalfont \ttfamily rockstar} N-body importer class.
     module procedure rockstarConstructorParameters
     module procedure rockstarConstructorInternal
  end interface nbodyImporterRockstar

  ! Enumeration of Rockstar columns.
  !# <enumeration>
  !#  <name>rockstarColumn</name>
  !#  <description>Enumeration of columns in Rockstar output files.</description>
  !#  <indexing>0</indexing>
  !#  <visibility>public</visibility>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <validator>yes</validator>
  !#  <entry label="scale"     />
  !#  <entry label="id"        />
  !#  <entry label="desc_scale"/>
  !#  <entry label="desc_id"   />
  !#  <entry label="num_prog"  />
  !#  <entry label="pid"       />
  !#  <entry label="upid"      />
  !#  <entry label="desc_pid"  />
  !#  <entry label="phantom"   />
  !#  <entry label="Mvir"      />
  !# </enumeration>

  !# <enumeration>
  !#  <name>columnType</name>
  !#  <description>Enumeration of columns types in Rockstar output files.</description>
  !#  <visibility>private</visibility>
  !#  <entry label="integer"/>
  !#  <entry label="real"   />
  !# </enumeration>
  
contains

  function rockstarConstructorParameters(parameters) result (self)
    !% Constructor for the {\normalfont \ttfamily rockstar} N-body importer class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyImporterRockstar   )                              :: self
    type   (inputParameters         ), intent(inout)               :: parameters
    class  (cosmologyParametersClass), pointer                     :: cosmologyParameters_
    type   (varying_string          ), allocatable  , dimension(:) :: readColumnsText
    integer                          , allocatable  , dimension(:) :: readColumns
    type   (varying_string          )                              :: fileName            , label
    integer                                                        :: i
    
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
    if (parameters%isPresent('readColumns')) then
       allocate(readColumnsText(parameters%count('readColumns')))
       allocate(readColumns    (parameters%count('readColumns')))
       !# <inputParameter>
       !#   <name>readColumns</name>
       !#   <source>parameters</source>
       !#   <variable>readColumnsText</variable>
       !#   <description>The names of additional columns to read.</description>
       !#   <type>string</type>
       !#   <cardinality>0..1</cardinality>
       !# </inputParameter>
       do i=1,size(readColumns)
          readColumns(i)=enumerationRockstarColumnEncode(char(readColumnsText(i)),includesPrefix=.false.)
       end do
    end if
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <conditionalCall>
    !# <call>
    !#  self=nbodyImporterRockstar(fileName,label,cosmologyParameters_{conditions})
    !# </call>
    !# <argument name="readColumns" value="readColumns" parameterPresent="parameters" />
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    return
  end function rockstarConstructorParameters

  function rockstarConstructorInternal(fileName,label,cosmologyParameters_,readColumns) result (self)
    !% Internal constructor for the {\normalfont \ttfamily rockstar} N-body importer class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type   (nbodyImporterRockstar   )                                        :: self
    class  (cosmologyParametersClass), intent(in   ), target                 :: cosmologyParameters_
    integer                          , intent(in   ), dimension(:), optional :: readColumns
    type   (varying_string          ), intent(in   )                         :: fileName            , label
    integer                                                                  :: i
    !# <constructorAssign variables="fileName, label, *cosmologyParameters_, readColumns"/>

    ! If extra read columns are requested, determine the number of integer and real columns.
    if (allocated(self%readColumns)) then
       self%readColumnsIntegerCount=0
       self%readColumnsRealCount   =0
       allocate(self%readColumnsType(size(self%readColumns)))
       do i=1,size(self%readColumns)
          select case (self%readColumns(i))
          case   (                          &
               &  rockstarColumnScale     , &
               &  rockstarColumnDesc_scale, &
               &  rockstarColumnMvir        &
               & )
             self%readColumnsRealCount      =self%readColumnsRealCount   +1
             self%readColumnsType        (i)=columnTypeReal
          case   (                          &
               &  rockstarColumnId        , &
               &  rockstarColumnDesc_id   , &
               &  rockstarColumnNum_prog  , &
               &  rockstarColumnPid       , &
               &  rockstarColumnUpid      , &
               &  rockstarColumnDesc_pid  , &
               &  rockstarColumnPhantom     &
               & )
             self%readColumnsIntegerCount   =self%readColumnsIntegerCount+1
             self%readColumnsType        (i)=columnTypeInteger
          case default
             call Galacticus_Error_Report('unknown column'//{introspection:location})
          end select
       end do
    end if
    return
  end function rockstarConstructorInternal

  subroutine rockstarDestructor(self)
    !% Destructor for {\normalfont \ttfamily rockstar} importer class.
    implicit none
    type(nbodyImporterRockstar), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    return
  end subroutine rockstarDestructor

  subroutine rockstarImport(self,simulations)
    !% Import data from a Rockstar file.
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Galacticus_Display  , only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Counter, Galacticus_Display_Counter_Clear, &
         &                              verbosityStandard
    use :: Memory_Management   , only : allocateArray            , deallocateArray
    use :: File_Utilities      , only : Count_Lines_in_File
    use :: Hashes              , only : rank1IntegerSizeTHash    , rank1DoubleHash
    implicit none
    class           (nbodyImporterRockstar), intent(inout)                                 :: self
    type            (nBodyData            ), intent(  out), dimension( :    ), allocatable :: simulations
    double precision                                      , dimension( :    ), allocatable :: expansionFactor
    double precision                                      , dimension( :  ,:), allocatable :: propertiesReal
    integer         (c_size_t             )               , dimension( :  ,:), allocatable :: propertiesInteger
    double precision                                      , dimension(0:22  )              :: columnsReal
    integer         (c_size_t             ), dimension(0:22  )                             :: columnsInteger
    integer         (c_size_t             )                                                :: countHalos       , countTrees, &
         &                                                                                    i
    integer                                                                                :: status           , j         , &
         &                                                                                    file             , jInteger  , &
         &                                                                                    jReal
    character       (len=1024             )                                                :: line
    character       (len=64               )                                                :: columnName
    logical                                                                                :: isComment

    call Galacticus_Display_Indent('import simulation from Rockstar file',verbosityStandard)
    allocate(simulations(1))
    simulations(1)%label=self%label
    ! Count lines in file. (Subtract 1 since one lines gives the number of distinct trees.)
    countHalos=Count_Lines_In_File(self%fileName,comment_char="#")-1_c_size_t
    countTrees=                                                    0_c_size_t
    ! Allocate storage
    call allocateArray(simulations(1)%position       ,[3_c_size_t,countHalos])
    call allocateArray(simulations(1)%velocity       ,[3_c_size_t,countHalos])
    call allocateArray(simulations(1)%particleIDs    ,[           countHalos])
    call allocateArray(               expansionFactor,[           countHalos])
    if (allocated(self%readColumns)) then
       allocate(propertiesInteger(countHalos,self%readColumnsIntegerCount))
       allocate(propertiesReal   (countHalos,self%readColumnsRealCount   ))
    else
       allocate(propertiesInteger(         0,                           0))
       allocate(propertiesReal   (         0,                           0))
    end if
    ! Open the file and read.
    i=0_c_size_t
    open(newUnit=file,file=char(self%fileName),status='old',form='formatted',iostat=status)
    do while (status == 0)
       read (file,'(a)',iostat=status) line
       isComment=line(1:1) == "#"
       if (.not.isComment) then
          if (countTrees > 0_c_size_t) then
             ! We have already determined the number of trees.
             i=i+1_c_size_t
             read (line,*) columnsReal   (0   ), &
                  &        columnsInteger(1   ), &
                  &        columnsReal   (2   ), &
                  &        columnsInteger(3: 8), &
                  &        columnsReal   (9:22)
             expansionFactor            (  i)=columnsReal   ( 0   )
             simulations(1) %particleIDs(  i)=columnsInteger( 1   )
             simulations(1) %position   (:,i)=columnsReal   (17:19)
             simulations(1) %velocity   (:,i)=columnsReal   (20:22)
             ! Read any extra columns.
             if (allocated(self%readColumns)) then
                jInteger=0
                jReal   =0
                do j=1,size(self%readColumns)
                   select case (self%readColumnsType(j))
                   case (columnTypeInteger)
                      jInteger                     =jInteger                           +1
                      propertiesInteger(i,jInteger)=columnsInteger(self%readColumns(j))
                   case (columnTypeReal   )
                      jReal                        =jReal                              +1
                      propertiesReal   (i,jReal   )=columnsReal   (self%readColumns(j))
                   end select
                end do
             end if
          else
             ! We have not yet determined the number of trees - read that now.
             read (line,*) countTrees
          end if
          ! Exit once all halos have been read.
          if (i == countHalos) exit
       end if
       call Galacticus_Display_Counter(int(100.0d0*dble(i)/dble(countHalos)),verbosity=verbosityStandard,isNew=i == 1_c_size_t)
    end do
    call Galacticus_Display_Counter_Clear(verbosityStandard)
    close(file)
    ! Convert position to internal units (physical Mpc).
    do i=1_c_size_t,countHalos
       simulations(1)%position(:,i)=+simulations(1)           %position       (:,i                 ) &
            &                       *                          expansionFactor(  i                 ) &
            &                       /self%cosmologyParameters_%HubbleConstant (  hubbleUnitsLittleH)
    end do
    call deallocateArray(expansionFactor)
    ! Add any additional properties.
    simulations(1)%propertiesInteger=rank1IntegerSizeTHash()
    simulations(1)%propertiesReal   =rank1DoubleHash      ()
    if (allocated(self%readColumns)) then
       jInteger                        =0
       jReal                           =0
       do j=1,size(self%readColumns)
          select case (self%readColumnsType(j))
          case (columnTypeInteger)
             jInteger=jInteger+1
             select case (self%readColumns(j))
             case (rockstarColumnId        )
                columnName='particleID'
             case (rockstarColumnDesc_id   )
                columnName='descendentID'
             case (rockstarColumnNum_prog  )
                columnName='progenitorCount'
             case (rockstarColumnPid       )
                columnName='hostID'
             case (rockstarColumnUpid      )
                columnName='hostRootID'
             case (rockstarColumnDesc_pid  )
                columnName='descendentHostID'
             case (rockstarColumnPhantom   )
                columnName='isPhantom'
             end select
             call simulations(1)%propertiesInteger%set(columnName,propertiesInteger(:,jInteger))
          case (columnTypeReal   )
             jReal   =jReal   +1
             select case (self%readColumns(j))
             case (rockstarColumnScale     )
                columnName='expansionFactor'
             case (rockstarColumnDesc_scale)
                columnName='descendentExpansionFactor'
             case (rockstarColumnMvir      )
                columnName='massVirial'
                propertiesReal(:,jReal)=+propertiesReal                                    (:,jReal             ) &
                     &                  /self          %cosmologyParameters_%HubbleConstant(  hubbleUnitsLittleH)
             end select
             call simulations(1)%propertiesReal   %set(columnName,propertiesReal   (:,jReal   ))
          end select
       end do
    end if
    call Galacticus_Display_Unindent('done',verbosityStandard)
    return
  end subroutine rockstarImport

  logical function rockstarIsHDF5(self)
    !% Return whether or not the imported data is from an HDF5 file.
    implicit none
    class(nbodyImporterRockstar), intent(inout) :: self

    rockstarIsHDF5=.false.
    return
  end function rockstarIsHDF5
