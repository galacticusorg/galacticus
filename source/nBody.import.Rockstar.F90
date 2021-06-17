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
Contains a module which implements an N-body data importer for Rockstar files.
!!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  
  !![
  <nbodyImporter name="nbodyImporterRockstar">
   <description>An importer for Rockstar files.</description>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterRockstar
     !!{
     An importer for Rockstar files.
     !!}
     private
     class  (cosmologyParametersClass)              , pointer     :: cosmologyParameters_    => null()
     integer                          , dimension(:), allocatable :: readColumns                      , readColumnsType
     integer                                                      :: readColumnsIntegerCount          , readColumnsRealCount
     type   (varying_string          )                            :: fileName                         , label
     logical                                                      :: havePosition                     , haveVelocity        , &
          &                                                          expansionFactorNeeded
   contains
     final     ::           rockstarDestructor
     procedure :: import => rockstarImport
     procedure :: isHDF5 => rockstarIsHDF5
  end type nbodyImporterRockstar

  interface nbodyImporterRockstar
     !!{
     Constructors for the {\normalfont \ttfamily rockstar} N-body importer class.
     !!}
     module procedure rockstarConstructorParameters
     module procedure rockstarConstructorInternal
  end interface nbodyImporterRockstar

  ! Enumeration of Rockstar columns.
  !![
  <enumeration>
   <name>rockstarColumn</name>
   <description>Enumeration of columns in Rockstar output files.</description>
   <indexing>0</indexing>
   <visibility>public</visibility>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <entry label="scale"                          />
   <entry label="id"                             />
   <entry label="desc_scale"                     />
   <entry label="desc_id"                        />
   <entry label="num_prog"                       />
   <entry label="pid"                            />
   <entry label="upid"                           />
   <entry label="desc_pid"                       />
   <entry label="phantom"                        />
   <entry label="sam_Mvir"                       />
   <entry label="Mvir"                           />
   <entry label="Rvir"                           />
   <entry label="rs"                             />
   <entry label="vrms"                           />
   <entry label="mmp"                            />
   <entry label="scale_of_last_MM"               />
   <entry label="Vmax"                           />
   <entry label="X"                              />
   <entry label="Y"                              />
   <entry label="Z"                              />
   <entry label="VX"                             />
   <entry label="VY"                             />
   <entry label="VZ"                             />
   <entry label="JX"                             />
   <entry label="JY"                             />
   <entry label="JZ"                             />
   <entry label="Spin"                           />
   <entry label="Breadth_first_ID"               />
   <entry label="Depth_first_ID"                 />
   <entry label="Tree_root_ID"                   />
   <entry label="Orig_halo_ID"                   />
   <entry label="Snap_num"                       />
   <entry label="Next_coprogenitor_depthfirst_ID"/>
   <entry label="Last_progenitor_depthfirst_ID"  />
   <entry label="Last_mainleaf_depthfirst_ID"    />
   <entry label="Tidal_Force"                    />
   <entry label="Tidal_ID"                       />
   <entry label="Rs_Klypin"                      />
   <entry label="Mmvir_all"                      />
   <entry label="M200b"                          />
   <entry label="M200c"                          />
   <entry label="M500c"                          />
   <entry label="M2500c"                         />
   <entry label="Xoff"                           />
   <entry label="Voff"                           />
   <entry label="Spin_Bullock"                   />
   <entry label="b_to_a"                         />
   <entry label="c_to_a"                         />
   <entry label="Ax"                             />
   <entry label="Ay"                             />
   <entry label="Az"                             />
   <entry label="b_to_a500c"                     />
   <entry label="c_to_a500c"                     />
   <entry label="Ax500c"                         />
   <entry label="Ay500c"                         />
   <entry label="Az500c"                         />
   <entry label="TU"                             />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>columnType</name>
   <description>Enumeration of columns types in Rockstar output files.</description>
   <visibility>private</visibility>
   <entry label="integer"/>
   <entry label="real"   />
  </enumeration>
  !!]
  
contains

  function rockstarConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily rockstar} N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyImporterRockstar   )                              :: self
    type   (inputParameters         ), intent(inout)               :: parameters
    class  (cosmologyParametersClass), pointer                     :: cosmologyParameters_
    type   (varying_string          ), allocatable  , dimension(:) :: readColumnsText
    integer                          , allocatable  , dimension(:) :: readColumns
    type   (varying_string          )                              :: fileName            , label
    integer                                                        :: i
    
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
    !!]
    if (parameters%isPresent('readColumns')) then
       allocate(readColumnsText(parameters%count('readColumns')))
       allocate(readColumns    (parameters%count('readColumns')))
       !![
       <inputParameter>
         <name>readColumns</name>
         <source>parameters</source>
         <variable>readColumnsText</variable>
         <description>The names of additional columns to read.</description>
       </inputParameter>
       !!]
       do i=1,size(readColumns)
          readColumns(i)=enumerationRockstarColumnEncode(char(readColumnsText(i)),includesPrefix=.false.)
       end do
    end if
    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <conditionalCall>
    <call>
     self=nbodyImporterRockstar(fileName,label,cosmologyParameters_{conditions})
    </call>
    <argument name="readColumns" value="readColumns" parameterPresent="parameters" />
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function rockstarConstructorParameters

  function rockstarConstructorInternal(fileName,label,cosmologyParameters_,readColumns) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily rockstar} N-body importer class.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type   (nbodyImporterRockstar   )                                        :: self
    class  (cosmologyParametersClass), intent(in   ), target                 :: cosmologyParameters_
    integer                          , intent(in   ), dimension(:), optional :: readColumns
    type   (varying_string          ), intent(in   )                         :: fileName            , label
    integer                                                                  :: i
    !![
    <constructorAssign variables="fileName, label, *cosmologyParameters_, readColumns"/>
    !!]

    ! If extra read columns are requested, determine the number of integer and real columns.
    if (allocated(self%readColumns)) then
       self%readColumnsIntegerCount=0
       self%readColumnsRealCount   =0
       allocate(self%readColumnsType(size(self%readColumns)))
       do i=1,size(self%readColumns)
          select case (self%readColumns(i))
          case   (                                               &
               &  rockstarColumnScale                          , &
               &  rockstarColumnDesc_scale                     , &
               &  rockstarColumnSam_Mvir                       , &
               &  rockstarColumnMvir                           , &
               &  rockstarColumnRvir                           , &
               &  rockstarColumnRs                             , &
               &  rockstarColumnVrms                           , &
               &  rockstarColumnscale_of_last_MM               , &
               &  rockstarColumnVmax                           , &
               &  rockstarColumnX                              , &
               &  rockstarColumnY                              , &
               &  rockstarColumnZ                              , &
               &  rockstarColumnVX                             , &
               &  rockstarColumnVY                             , &
               &  rockstarColumnVZ                             , &
               &  rockstarColumnJX                             , &
               &  rockstarColumnJY                             , &
               &  rockstarColumnJZ                             , &
               &  rockstarColumnSpin                           , &
               &  rockstarColumnTidal_Force                    , &
               &  rockstarColumnRs_Klypin                      , &
               &  rockstarColumnMmvir_all                      , &
               &  rockstarColumnM200b                          , &
               &  rockstarColumnM200c                          , &
               &  rockstarColumnM500c                          , &
               &  rockstarColumnM2500c                         , &
               &  rockstarColumnXoff                           , &
               &  rockstarColumnVoff                           , &
               &  rockstarColumnSpin_Bullock                   , &
               &  rockstarColumnb_to_a                         , &
               &  rockstarColumnc_to_a                         , &
               &  rockstarColumnAx                             , &
               &  rockstarColumnAy                             , &
               &  rockstarColumnAz                             , &
               &  rockstarColumnb_to_a500c                     , &
               &  rockstarColumnc_to_a500c                     , &
               &  rockstarColumnAx500c                         , &
               &  rockstarColumnAy500c                         , &
               &  rockstarColumnAz500c                         , &
               &  rockstarColumnTU                               &
               & )
             self%readColumnsRealCount      =self%readColumnsRealCount   +1
             self%readColumnsType        (i)=columnTypeReal
          case   (                                               &
               &  rockstarColumnId                             , &
               &  rockstarColumnDesc_id                        , &
               &  rockstarColumnNum_prog                       , &
               &  rockstarColumnPid                            , &
               &  rockstarColumnUpid                           , &
               &  rockstarColumnDesc_pid                       , &
               &  rockstarColumnPhantom                        , &
               &  rockstarColumnMmp                            , &
               &  rockstarColumnBreadth_first_ID               , &
               &  rockstarColumnDepth_first_ID                 , &
               &  rockstarColumnTree_root_ID                   , &
               &  rockstarColumnOrig_halo_ID                   , &
               &  rockstarColumnSnap_num                       , &
               &  rockstarColumnNext_coprogenitor_depthfirst_ID, &
               &  rockstarColumnLast_progenitor_depthfirst_ID  , &
               &  rockstarColumnLast_mainleaf_depthfirst_ID    , &
               &  rockstarColumnTidal_ID                         &
               & )
             self%readColumnsIntegerCount   =self%readColumnsIntegerCount+1
             self%readColumnsType        (i)=columnTypeInteger
          case default
             call Galacticus_Error_Report('unknown column'//{introspection:location})
          end select
       end do
       self%havePosition         =any(self%readColumns == rockstarColumnX ) .and. any(self%readColumns == rockstarColumnY ) .and. any(self%readColumns == rockstarColumnZ )
       self%haveVelocity         =any(self%readColumns == rockstarColumnVX) .and. any(self%readColumns == rockstarColumnVY) .and. any(self%readColumns == rockstarColumnVZ)
       self%expansionFactorNeeded=any(self%readColumns == rockstarColumnX ) .or.  any(self%readColumns == rockstarColumnY ) .or.  any(self%readColumns == rockstarColumnZ )
    else
       self%havePosition         =.false.
       self%haveVelocity         =.false.
       self%expansionFactorNeeded=.false.
    end if
    return
  end function rockstarConstructorInternal

  subroutine rockstarDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily rockstar} importer class.
    !!}
    implicit none
    type(nbodyImporterRockstar), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine rockstarDestructor

  subroutine rockstarImport(self,simulations)
    !!{
    Import data from a Rockstar file.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Display             , only : displayCounter        , displayCounterClear     , displayIndent     , displayUnindent         , &
          &                             verbosityLevelStandard
    use :: File_Utilities      , only : Count_Lines_in_File
    use :: Hashes              , only : doubleHash            , integerSizeTHash        , rank1DoublePtrHash, rank1IntegerSizeTPtrHash, &
          &                             rank2DoublePtrHash    , rank2IntegerSizeTPtrHash, varyingStringHash
    use :: Memory_Management   , only : allocateArray         , deallocateArray
    implicit none
    class           (nbodyImporterRockstar     ), intent(inout)                                 :: self
    type            (nBodyData                 ), intent(  out), dimension( :    ), allocatable :: simulations
    double precision                                           , dimension( :    ), allocatable :: expansionFactor
    type            (nbodyPropertiesRealList   )               , dimension( :    ), allocatable :: propertiesReal
    type            (nbodyPropertiesIntegerList)               , dimension( :    ), allocatable :: propertiesInteger
    double precision                                           , dimension( :  ,:), pointer     :: position         , velocity
    double precision                                           , dimension(0:56  )              :: columnsReal
    integer         (c_size_t                  )               , dimension(0:56  )              :: columnsInteger
    integer         (c_size_t                  )                                                :: countHalos       , countTrees, &
         &                                                                                         i
    integer                                                                                     :: status           , j         , &
         &                                                                                         file             , jInteger  , &
         &                                                                                         jReal
    character       (len=1024                  )                                                :: line
    character       (len=64                    )                                                :: columnName
    logical                                                                                     :: isComment
    double precision                                                                            :: boxSize

    call displayIndent('import simulation from Rockstar file',verbosityLevelStandard)
    allocate(simulations(1))
    simulations(1)%label=self%label
    ! Count lines in file. (Subtract 1 since one lines gives the number of distinct trees.)
    countHalos=Count_Lines_In_File(self%fileName,comment_char="#")-1_c_size_t
    countTrees=                                                    0_c_size_t    
    ! Allocate storage
    if (self%expansionFactorNeeded) allocate(expansionFactor(countHalos))
    if (allocated(self%readColumns)) then
       allocate(propertiesInteger(self%readColumnsIntegerCount))
       allocate(propertiesReal   (self%readColumnsRealCount   ))
       do i=1,self%readColumnsIntegerCount
          allocate(propertiesInteger(i)%property(countHalos))
       end do
       do i=1,self%readColumnsRealCount
          allocate(propertiesReal   (i)%property(countHalos))
       end do
    else
       allocate(propertiesInteger(0                           ))
       allocate(propertiesReal   (0                           ))
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
             read (line,*) columnsReal   ( 0   ), &
                  &        columnsInteger( 1   ), &
                  &        columnsReal   ( 2   ), &
                  &        columnsInteger( 3: 8), &
                  &        columnsReal   ( 9:13), &
                  &        columnsInteger(14:14), &
                  &        columnsReal   (15:26), &
                  &        columnsInteger(27:34), &
                  &        columnsReal   (35   ), &
                  &        columnsInteger(36   ), &
                  &        columnsReal   (37:56)
             if (self%expansionFactorNeeded) expansionFactor(i)=columnsReal(0)
             ! Read any extra columns.
             if (allocated(self%readColumns)) then
                jInteger=0
                jReal   =0
                do j=1,size(self%readColumns)
                   select case (self%readColumnsType(j))
                   case (columnTypeInteger)
                      jInteger                               =jInteger                           +1
                      propertiesInteger(jInteger)%property(i)=columnsInteger(self%readColumns(j))
                   case (columnTypeReal   )
                      jReal                                  =jReal                              +1
                      propertiesReal   (jReal   )%property(i)=columnsReal   (self%readColumns(j))
                   end select
                end do
             end if
          else
             ! We have not yet determined the number of trees - read that now.
             read (line,*) countTrees
          end if
          ! Exit once all halos have been read.
          if (i == countHalos) exit
       else if (line(1:16) == "#Full box size =") then
          ! Extract box size.
          read (line(17:len_trim(line)),*) boxSize
       end if
       call displayCounter(int(100.0d0*dble(i)/dble(countHalos)),verbosity=verbosityLevelStandard,isNew=i == 1_c_size_t)
    end do
    call displayCounterClear(verbosityLevelStandard)
    close(file)
    ! Convert positions to internal units (physical Mpc).
    jReal=0
    do j=1,size(self%readColumns)
       if (self%readColumnsType(j) /= columnTypeReal) cycle
       jReal=jReal+1
       select case (self%readColumns(j))
       case (rockstarColumnX,rockstarColumnY,rockstarColumnZ)
          propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
            &                            *                                           expansionFactor                    &
            &                            /self                 %cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH)
       end select
    end do
    if (allocated(expansionFactor)) deallocate(expansionFactor)
    ! Convert box size to internal units (comoving Mpc).
    boxSize=+boxSize                                                      &
         &  /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    ! Store attribues.
    simulations(1)%attributesInteger=integerSizeTHash ()
    simulations(1)%attributesReal   =doubleHash       ()
    simulations(1)%attributesText   =varyingStringHash()
    call simulations(1)%attributesReal%set('boxSize',boxSize)
    ! Add any additional properties.
    simulations(1)%propertiesInteger     =rank1IntegerSizeTPtrHash()
    simulations(1)%propertiesReal        =rank1DoublePtrHash      ()
    simulations(1)%propertiesIntegerRank1=rank2IntegerSizeTPtrHash()
    simulations(1)%propertiesRealRank1   =rank2DoublePtrHash      ()
    if (allocated(self%readColumns)) then
       if (self%havePosition) allocate(position(3,countHalos))
       if (self%haveVelocity) allocate(velocity(3,countHalos))
       jInteger                        =0
       jReal                           =0
       do j=1,size(self%readColumns)
          select case (self%readColumnsType(j))
          case (columnTypeInteger)
             jInteger=jInteger+1
             select case (self%readColumns(j))
             case (rockstarColumnId          )
                columnName='particleID'
             case (rockstarColumnDesc_id     )
                columnName='descendentID'
             case (rockstarColumnNum_prog    )
                columnName='progenitorCount'
             case (rockstarColumnPid         )
                columnName='hostID'
             case (rockstarColumnUpid        )
                columnName='isolatedHostID'
             case (rockstarColumnDesc_pid    )
                columnName='descendentHostID'
             case (rockstarColumnMmp         )
                columnName='isMostMassiveProgenitor'
             case (rockstarColumnPhantom     )
                columnName='isPhantom'
             case (rockstarColumnSnap_num    )
                columnName='snapshotID'
             case (rockstarColumnTree_root_ID)
                columnName='treeID'
             end select
             call simulations(1)%propertiesInteger%set(columnName,propertiesInteger(jInteger)%property)
          case (columnTypeReal   )
             jReal   =jReal   +1
             select case (self%readColumns(j))
             case (rockstarColumnScale     )
                columnName='expansionFactor'
             case (rockstarColumnDesc_scale)
                columnName='descendentExpansionFactor'
             case (rockstarColumnMvir      )
                columnName='massVirial'
                propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
                     &                         /self                 %cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
             case (rockstarColumnRvir      )
                columnName='radiusVirial'
                propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
                     &                         /self                 %cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
             case (rockstarColumnSpin      )
                columnName='spin'
             case (rockstarColumnrs        )
                columnName='radiusScale'
                propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
                     &                         /self                 %cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
             case (rockstarColumnTU        )
                columnName='virialRatio'
                propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
                     &                         *2.0d0
             end select
             if      (self%havePosition.and.self%readColumns(j) == rockstarColumnX ) then
                position(1,:)=propertiesReal(jReal)%property
             else if (self%havePosition.and.self%readColumns(j) == rockstarColumnY ) then
                position(2,:)=propertiesReal(jReal)%property
             else if (self%havePosition.and.self%readColumns(j) == rockstarColumnZ ) then
                position(3,:)=propertiesReal(jReal)%property
             else if (self%haveVelocity.and.self%readColumns(j) == rockstarColumnVX) then
                velocity(1,:)=propertiesReal(jReal)%property
             else if (self%haveVelocity.and.self%readColumns(j) == rockstarColumnVY) then
                velocity(2,:)=propertiesReal(jReal)%property
             else if (self%haveVelocity.and.self%readColumns(j) == rockstarColumnVZ) then
                velocity(3,:)=propertiesReal(jReal)%property
             else
                call simulations(1)%propertiesReal%set(columnName,propertiesReal(jReal)%property)
             end if
          end select
       end do
       if (self%havePosition) call simulations(1)%propertiesRealRank1%set('position',position)
       if (self%haveVelocity) call simulations(1)%propertiesRealRank1%set('velocity',velocity)
    end if
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine rockstarImport

  logical function rockstarIsHDF5(self)
    !!{
    Return whether or not the imported data is from an HDF5 file.
    !!}
    implicit none
    class(nbodyImporterRockstar), intent(inout) :: self

    rockstarIsHDF5=.false.
    return
  end function rockstarIsHDF5
