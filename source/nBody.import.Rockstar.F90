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
Implements an N-body data importer for Rockstar files.
!!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  ! Enumeration of Rockstar columns.
  !![
  <enumeration>
   <name>rockstarColumn</name>
   <description>Enumeration of columns in Rockstar output files.</description>
   <indexing>0</indexing>
   <visibility>public</visibility>
   <decodeFunction>yes</decodeFunction>
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
   <entry label="M_pe_Behroozi"                  />
   <entry label="M_pe_Diemer"                    />
   <entry label="Halfmass_Radius"                />
   <entry label="RVmax"                          />
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
  
  !![
  <nbodyImporter name="nbodyImporterRockstar">
   <description>An importer for Rockstar files.</description>
   <runTimeFileDependencies paths="fileName"/>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterRockstar
     !!{
     An importer for Rockstar files.
     !!}
     private
     class  (cosmologyParametersClass     )              , pointer     :: cosmologyParameters_    => null()
     type   (enumerationRockstarColumnType), dimension(:), allocatable :: readColumns                      , readColumnsMapped
     type   (enumerationColumnTypeType    ), dimension(:), allocatable :: readColumnsType
     integer                                                           :: readColumnsIntegerCount          , readColumnsRealCount
     type   (varying_string               )                            :: fileName                         , label
     logical                                                           :: havePosition                     , haveVelocity        , &
          &                                                               expansionFactorNeeded
   contains
     final     ::           rockstarDestructor
     procedure :: import => rockstarImport
     procedure :: isHDF5 => rockstarIsHDF5
  end type nbodyImporterRockstar

  interface nbodyImporterRockstar
     !!{
     Constructors for the \refClass{nbodyImporterRockstar} N-body importer class.
     !!}
     module procedure rockstarConstructorParameters
     module procedure rockstarConstructorInternal
  end interface nbodyImporterRockstar

contains

  function rockstarConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyImporterRockstar} N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyImporterRockstar        )                              :: self
    type   (inputParameters              ), intent(inout)               :: parameters
    class  (cosmologyParametersClass     ), pointer                     :: cosmologyParameters_
    type   (varying_string               ), allocatable  , dimension(:) :: readColumns
    type   (enumerationRockstarColumnType), allocatable  , dimension(:) :: readColumns_
    type   (varying_string               )                              :: fileName            , label
    integer                                                             :: i
    
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
       allocate(readColumns (parameters%count('readColumns')))
       allocate(readColumns_(parameters%count('readColumns')))
       !![
       <inputParameter>
         <name>readColumns</name>
         <source>parameters</source>
         <description>The names of additional columns to read.</description>
       </inputParameter>
       !!]
       do i=1,size(readColumns)
          readColumns_(i)=enumerationRockstarColumnEncode(char(readColumns(i)),includesPrefix=.false.)
       end do
    end if
    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <conditionalCall>
    <call>
     self=nbodyImporterRockstar(fileName,label,cosmologyParameters_{conditions})
    </call>
    <argument name="readColumns" value="readColumns_" parameterPresent="parameters" />
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function rockstarConstructorParameters

  function rockstarConstructorInternal(fileName,label,cosmologyParameters_,readColumns) result (self)
    !!{
    Internal constructor for the \refClass{nbodyImporterRockstar} N-body importer class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (nbodyImporterRockstar        )                                        :: self
    class  (cosmologyParametersClass     ), intent(in   ), target                 :: cosmologyParameters_
    type   (enumerationRockstarColumnType), intent(in   ), dimension(:), optional :: readColumns
    type   (varying_string               ), intent(in   )                         :: fileName            , label
    integer                                                                  :: i
    !![
    <constructorAssign variables="fileName, label, *cosmologyParameters_, readColumns"/>
    !!]

    ! If extra read columns are requested, determine the number of integer and real columns.
    if (allocated(self%readColumns)) then
       self%readColumnsIntegerCount=0
       self%readColumnsRealCount   =0
       allocate(self%readColumnsType  (size(self%readColumns)))
       allocate(self%readColumnsMapped(size(self%readColumns)))
       self%readColumnsMapped=self%readColumns
       do i=1,size(self%readColumns)
          select case (self%readColumns(i)%ID)
          case   (                                                  &
               &  rockstarColumnScale                          %ID, &
               &  rockstarColumnDesc_scale                     %ID, &
               &  rockstarColumnSam_Mvir                       %ID, &
               &  rockstarColumnMvir                           %ID, &
               &  rockstarColumnRvir                           %ID, &
               &  rockstarColumnRs                             %ID, &
               &  rockstarColumnVrms                           %ID, &
               &  rockstarColumnscale_of_last_MM               %ID, &
               &  rockstarColumnVmax                           %ID, &
               &  rockstarColumnX                              %ID, &
               &  rockstarColumnY                              %ID, &
               &  rockstarColumnZ                              %ID, &
               &  rockstarColumnVX                             %ID, &
               &  rockstarColumnVY                             %ID, &
               &  rockstarColumnVZ                             %ID, &
               &  rockstarColumnJX                             %ID, &
               &  rockstarColumnJY                             %ID, &
               &  rockstarColumnJZ                             %ID, &
               &  rockstarColumnSpin                           %ID, &
               &  rockstarColumnTidal_Force                    %ID, &
               &  rockstarColumnRs_Klypin                      %ID, &
               &  rockstarColumnMmvir_all                      %ID, &
               &  rockstarColumnM200b                          %ID, &
               &  rockstarColumnM200c                          %ID, &
               &  rockstarColumnM500c                          %ID, &
               &  rockstarColumnM2500c                         %ID, &
               &  rockstarColumnXoff                           %ID, &
               &  rockstarColumnVoff                           %ID, &
               &  rockstarColumnSpin_Bullock                   %ID, &
               &  rockstarColumnb_to_a                         %ID, &
               &  rockstarColumnc_to_a                         %ID, &
               &  rockstarColumnAx                             %ID, &
               &  rockstarColumnAy                             %ID, &
               &  rockstarColumnAz                             %ID, &
               &  rockstarColumnb_to_a500c                     %ID, &
               &  rockstarColumnc_to_a500c                     %ID, &
               &  rockstarColumnAx500c                         %ID, &
               &  rockstarColumnAy500c                         %ID, &
               &  rockstarColumnAz500c                         %ID, &
               &  rockstarColumnTU                             %ID, &
               &  rockstarColumnM_pe_Behroozi                  %ID, &
               &  rockstarColumnM_pe_Diemer                    %ID, &
               &  rockstarColumnHalfmass_Radius                %ID, &
               &  rockstarColumnRVmax                          %ID  &
               & )
             self%readColumnsRealCount      =self%readColumnsRealCount   +1
             self%readColumnsType        (i)=columnTypeReal
          case   (                                                  &
               &  rockstarColumnId                             %ID, &
               &  rockstarColumnDesc_id                        %ID, &
               &  rockstarColumnNum_prog                       %ID, &
               &  rockstarColumnPid                            %ID, &
               &  rockstarColumnUpid                           %ID, &
               &  rockstarColumnDesc_pid                       %ID, &
               &  rockstarColumnPhantom                        %ID, &
               &  rockstarColumnMmp                            %ID, &
               &  rockstarColumnBreadth_first_ID               %ID, &
               &  rockstarColumnDepth_first_ID                 %ID, &
               &  rockstarColumnTree_root_ID                   %ID, &
               &  rockstarColumnOrig_halo_ID                   %ID, &
               &  rockstarColumnSnap_num                       %ID, &
               &  rockstarColumnNext_coprogenitor_depthfirst_ID%ID, &
               &  rockstarColumnLast_progenitor_depthfirst_ID  %ID, &
               &  rockstarColumnLast_mainleaf_depthfirst_ID    %ID, &
               &  rockstarColumnTidal_ID                       %ID  &
               & )
             self%readColumnsIntegerCount   =self%readColumnsIntegerCount+1
             self%readColumnsType        (i)=columnTypeInteger
          case default
             call Error_Report('unknown column'//{introspection:location})
          end select
       end do
       self%havePosition         =     any(self%readColumns == rockstarColumnX   ) .and. any(self%readColumns == rockstarColumnY ) .and. any(self%readColumns == rockstarColumnZ   )
       self%haveVelocity         =     any(self%readColumns == rockstarColumnVX  ) .and. any(self%readColumns == rockstarColumnVY) .and. any(self%readColumns == rockstarColumnVZ  )
       self%expansionFactorNeeded=     any(self%readColumns == rockstarColumnX   ) .or.  any(self%readColumns == rockstarColumnY ) .or.  any(self%readColumns == rockstarColumnZ   ) &
            &                     .or. any(self%readColumns == rockstarColumnRvir) .or.  any(self%readColumns == rockstarColumnRs) .or.  any(self%readColumns == rockstarColumnXoff) &
            &                     .or. any(self%readColumns == rockstarColumnAx  ) .or.  any(self%readColumns == rockstarColumnAy) .or.  any(self%readColumns == rockstarColumnAz  )
    else
       self%havePosition         =.false.
       self%haveVelocity         =.false.
       self%expansionFactorNeeded=.false.
    end if
    return
  end function rockstarConstructorInternal

  subroutine rockstarDestructor(self)
    !!{
    Destructor for the \refClass{nbodyImporterRockstar} importer class.
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
    use :: Cosmology_Parameters        , only : hubbleUnitsLittleH
    use :: Display                     , only : displayCounter        , displayCounterClear     , displayIndent     , displayUnindent         , &
          &                                     verbosityLevelStandard
    use :: Error                       , only : Error_Report
    use :: File_Utilities              , only : Count_Lines_in_File
    use :: Hashes                      , only : doubleHash            , integerSizeTHash        , rank1DoublePtrHash, rank1IntegerSizeTPtrHash, &
          &                                     rank2DoublePtrHash    , rank2IntegerSizeTPtrHash, varyingStringHash , genericHash
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: String_Handling             , only : String_Split_Words    , String_Count_Words
    implicit none
    class           (nbodyImporterRockstar     ), intent(inout)                                 :: self
    type            (nBodyData                 ), intent(  out), dimension( :    ), allocatable :: simulations
    double precision                                           , dimension( :    ), allocatable :: expansionFactor
    type            (nbodyPropertiesRealList   )               , dimension( :    ), allocatable :: propertiesReal
    type            (nbodyPropertiesIntegerList)               , dimension( :    ), allocatable :: propertiesInteger
    double precision                                           , dimension( :  ,:), pointer     :: position         , velocity
    double precision                                           , dimension(0:60  )              :: columnsReal
    integer         (c_size_t                  )               , dimension(0:60  )              :: columnsInteger
    type            (varying_string            )               , dimension( :    ), allocatable :: columnNames
    integer         (c_size_t                  )                                                :: countHalos       , countTrees, &
         &                                                                                         i
    integer                                                                                     :: status           , j         , &
         &                                                                                         file             , jInteger  , &
         &                                                                                         jReal            , columnMap , &
         &                                                                                         lineStatus
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
    if (self%expansionFactorNeeded) then
       allocate(expansionFactor(countHalos))
    else
       allocate(expansionFactor(         0))
    end if
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
    if (countHalos > 0_c_size_t) then
       i        =0_c_size_t
       columnMap=-1
       open(newUnit=file,file=char(self%fileName),status='old',form='formatted',iostat=status)
       do while (status == 0)
          read (file,'(a)',iostat=status) line
          isComment=line(1:1) == "#"
          if (isComment .and. columnMap < 0) then
             ! Different versions of Rockstar have slightly different column layouts. Figure out which is used in this file.
             allocate(columnNames(String_Count_Words(line)))
             call String_Split_Words(columnNames,line)
             if (columnNames(36) == "Mvir_all"       ) then
                ! This version has `Mvir_all` in column 35 (note Rockstar columns are zero-indexed), and a total of 57 columns.
                !
                ! Example:
                !  #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Mmvir_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius
               columnMap=1
             else if      (columnNames(36) == "Rs_Klypin"      ) then
                ! This version has `Rs_Klypin` in column 35 (note Rockstar columns are zero-indexed), and a total of 58
                ! columns. Compared to columnMap=1, it added `Rs_Klypin` in column 36, shifting later columns.
                !
                ! Example:
                !  #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Rs_Klypin Mmvir_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius
                columnMap=2
             else if (columnNames(36) == "Tidal_Force(35)") then
                ! These versions have `Tidal_Force(35)` in column 35 (note Rockstar columns are zero-indexed).
                if      (size(columnNames) == 59) then
                   ! This version has 59 columns. Compared to columnMap=2, it added `Tidal_Force(35) Tidal_ID(36)` prior to
                   ! `Rs_Klypin`, shifting later columns.
                   !
                   ! Example:
                   !  #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) Rs_Klypin Mmvir_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer
                   columnMap=3
                else if (size(columnNames) == 60) then
                   ! This version has 60 columns. Compared to columnMap=3, it added `Halfmass_Radius` as a final column.
                   !
                   ! Note that there are a few different versions of this, where the names of mass columns differ (e.g. "Mmvir_all" vs. "Mvir_all").
                   !
                   ! Example:
                   ! #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_Mmvir(9) Mmvir(10) Rmvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_idx(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) Rs_Klypin Mmvir_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius
                   columnMap=4
                else if (size(columnNames) == 61) then
                   ! This version has 61 columns. Compared to columnMap=4, it added `RVmax` as a final column.
                   ! #scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_Mvir(9) Mvir(10) Rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_idx(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Tidal_Force(35) Tidal_ID(36) Rs_Klypin Mvir_all M200b M200c M500c M2500c Xoff Voff Spin_Bullock b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius RVmax
                   columnMap=5
                else
                   call Error_Report('unrecognized column layout'//{introspection:location})
                end if
             else
                call Error_Report   ('unrecognized column layout'//{introspection:location})
             end if
             ! Adjust column numbers to be read.
             select case (columnMap)
             case (1)
                do j=1,size(self%readColumns)
                   if     (                                                  &
                        &   self%readColumns(j) == rockstarColumnTidal_Force &
                        &  .or.                                              &
                        &   self%readColumns(j) == rockstarColumnTidal_ID    &
                        & ) call Error_Report('tidal properties not available'//{introspection:location})
                   if (self%readColumns(j)%ID > 33)                  &
                        & call self%readColumnsMapped(j)%subtract(3)
                   if (self%readColumnsMapped(j)%ID > 54)                                                                      &
                        & call Error_Report(                                                                                   &
                        &                   "property '"                                                                    // &
                        &                   enumerationRockstarColumnDecode(self%readColumnsMapped(j),includePrefix=.false.)// &
                        &                   "' not available"                                                               // &
                        &                   {introspection:location}                                                           &
                        &                  )
                end do
             case (2)
                do j=1,size(self%readColumns)
                   if     (                                                  &
                        &   self%readColumns(j) == rockstarColumnTidal_Force &
                        &  .or.                                              &
                        &   self%readColumns(j) == rockstarColumnTidal_ID    &
                        & ) call Error_Report('tidal properties not available'//{introspection:location})
                   if (self%readColumns(j)%ID > 34)                  &
                        & call self%readColumnsMapped(j)%subtract(2)
                   if (self%readColumnsMapped(j)%ID > 54)                                                                      &
                        & call Error_Report(                                                                                   &
                        &                   "property '"                                                                    // &
                        &                   enumerationRockstarColumnDecode(self%readColumnsMapped(j),includePrefix=.false.)// &
                        &                   "' not available"                                                               // &
                        &                   {introspection:location}                                                           &
                        &                  )
                end do
             case (3)
                do j=1,size(self%readColumns)
                   if (self%readColumnsMapped(j)%ID > 59)                                                                      &
                        & call Error_Report(                                                                                   &
                        &                   "property '"                                                                    // &
                        &                   enumerationRockstarColumnDecode(self%readColumnsMapped(j),includePrefix=.false.)// &
                        &                   "' not available"                                                               // &
                        &                   {introspection:location}                                                           &
                        &                  )
                end do
             case (4)
                do j=1,size(self%readColumns)
                   if (self%readColumnsMapped(j)%ID > 60)                                                                      &
                        & call Error_Report(                                                                                   &
                        &                   "property '"                                                                    // &
                        &                   enumerationRockstarColumnDecode(self%readColumnsMapped(j),includePrefix=.false.)// &
                        &                   "' not available"                                                               // &
                        &                   {introspection:location}                                                           &
                        &                  )
                end do
             case (5)
                do j=1,size(self%readColumns)
                   if (self%readColumnsMapped(j)%ID > 61)                                                                      &
                        & call Error_Report(                                                                                   &
                        &                   "property '"                                                                    // &
                        &                   enumerationRockstarColumnDecode(self%readColumnsMapped(j),includePrefix=.false.)// &
                        &                   "' not available"                                                               // &
                        &                   {introspection:location}                                                           &
                        &                  )
                end do
             end select
          end if
          if (.not.isComment) then
             if (countTrees > 0_c_size_t) then
                ! We have already determined the number of trees.
                i=i+1_c_size_t
                select case (columnMap)
                case (1)
                   ! Older layout with "Mvir_all" in column 35.
                   read (line,*,ioStat=lineStatus) columnsReal   ( 0   ), &
                        &                          columnsInteger( 1   ), &
                        &                          columnsReal   ( 2   ), &
                        &                          columnsInteger( 3: 8), &
                        &                          columnsReal   ( 9:13), &
                        &                          columnsInteger(14:14), &
                        &                          columnsReal   (15:26), &
                        &                          columnsInteger(27:33), &
                        &                          columnsReal   (35:54)
                case (2)
                   ! Older layout with "Rs_Klypin" in column 35.
                   read (line,*,ioStat=lineStatus) columnsReal   ( 0   ), &
                        &                          columnsInteger( 1   ), &
                        &                          columnsReal   ( 2   ), &
                        &                          columnsInteger( 3: 8), &
                        &                          columnsReal   ( 9:13), &
                        &                          columnsInteger(14:14), &
                        &                          columnsReal   (15:26), &
                        &                          columnsInteger(27:34), &
                        &                          columnsReal   (35:54)
                case (3)
                   ! Newer layout with "Tidal_Force" in column 35 and 58 columns total.
                   read (line,*,ioStat=lineStatus) columnsReal   ( 0   ), &
                        &                          columnsInteger( 1   ), &
                        &                          columnsReal   ( 2   ), &
                        &                          columnsInteger( 3: 8), &
                        &                          columnsReal   ( 9:13), &
                        &                          columnsInteger(14:14), &
                        &                          columnsReal   (15:26), &
                        &                          columnsInteger(27:34), &
                        &                          columnsReal   (35   ), &
                        &                          columnsInteger(36   ), &
                        &                          columnsReal   (37:58)
                case (4)
                   ! Newer layout with "Tidal_Force" in column 35 and Halfmass_Radius in column 59.
                   read (line,*,ioStat=lineStatus) columnsReal   ( 0   ), &
                        &                          columnsInteger( 1   ), &
                        &                          columnsReal   ( 2   ), &
                        &                          columnsInteger( 3: 8), &
                        &                          columnsReal   ( 9:13), &
                        &                          columnsInteger(14:14), &
                        &                          columnsReal   (15:26), &
                        &                          columnsInteger(27:34), &
                        &                          columnsReal   (35   ), &
                        &                          columnsInteger(36   ), &
                        &                          columnsReal   (37:59)
                case (5)
                   ! Newer layout with "Tidal_Force" in column 35 and RVmax in column 60.
                   read (line,*,ioStat=lineStatus) columnsReal   ( 0   ), &
                        &                          columnsInteger( 1   ), &
                        &                          columnsReal   ( 2   ), &
                        &                          columnsInteger( 3: 8), &
                        &                          columnsReal   ( 9:13), &
                        &                          columnsInteger(14:14), &
                        &                          columnsReal   (15:26), &
                        &                          columnsInteger(27:34), &
                        &                          columnsReal   (35   ), &
                        &                          columnsInteger(36   ), &
                        &                          columnsReal   (37:60)
                case default
                   call Error_Report('unknown column layout'//{introspection:location})
                end select
                if (lineStatus /= 0) then
                   block
                     character(len=1) :: columnMapLabel
                     write (columnMapLabel,'(i1)') columnMap
                     call Error_Report('failed to parse line (column map: '//columnMapLabel//'):'//char(10)//'"'//trim(line)//'"'//{introspection:location})
                   end block
                end if
                if (self%expansionFactorNeeded) expansionFactor(i)=columnsReal(0)
                ! Read any extra columns.
                if (allocated(self%readColumns)) then
                   jInteger=0
                   jReal   =0
                   do j=1,size(self%readColumns)
                      select case (self%readColumnsType(j)%ID)
                      case (columnTypeInteger%ID)
                         jInteger                               =jInteger                                    +1
                         propertiesInteger(jInteger)%property(i)=columnsInteger(self%readColumnsMapped(j)%ID)
                      case (columnTypeReal   %ID)
                         jReal                                  =jReal                                       +1
                         propertiesReal   (jReal   )%property(i)=columnsReal   (self%readColumnsMapped(j)%ID)
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
    end if
    ! Convert positions to internal units (physical Mpc).
    jReal=0
    do j=1,size(self%readColumns)
       if (self%readColumnsType(j) /= columnTypeReal) cycle
       jReal=jReal+1
       select case (self%readColumns(j)%ID)
       case (rockstarColumnX%ID,rockstarColumnY%ID,rockstarColumnZ%ID)
          propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
            &                            *                                           expansionFactor                    &
            &                            /self                 %cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH)
       end select
    end do
    ! Convert box size to internal units (comoving Mpc).
    boxSize=+boxSize                                                      &
         &  /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    ! Store attributes.
    simulations(1)%attributesInteger=integerSizeTHash ()
    simulations(1)%attributesReal   =doubleHash       ()
    simulations(1)%attributesText   =varyingStringHash()
    simulations(1)%attributesGeneric=genericHash      ()
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
          select case (self%readColumnsType(j)%ID)
          case (columnTypeInteger%ID)
             jInteger=jInteger+1
             select case (self%readColumns(j)%ID)
             case (rockstarColumnId          %ID)
                columnName='particleID'
             case (rockstarColumnDesc_id     %ID)
                columnName='descendantID'
             case (rockstarColumnNum_prog    %ID)
                columnName='progenitorCount'
             case (rockstarColumnPid         %ID)
                columnName='hostID'
             case (rockstarColumnUpid        %ID)
                columnName='isolatedHostID'
             case (rockstarColumnDesc_pid    %ID)
                columnName='descendantHostID'
             case (rockstarColumnMmp         %ID)
                columnName='isMostMassiveProgenitor'
             case (rockstarColumnPhantom     %ID)
                columnName='isPhantom'
             case (rockstarColumnSnap_num    %ID)
                columnName='snapshotID'
             case (rockstarColumnTree_root_ID%ID)
                columnName='treeID'
             end select
             call simulations(1)%propertiesInteger%set(columnName,propertiesInteger(jInteger)%property)
          case (columnTypeReal   %ID)
             jReal   =jReal   +1
             select case (self%readColumns(j)%ID)
             case (rockstarColumnScale     %ID)
                columnName='expansionFactor'
             case (rockstarColumnDesc_scale%ID)
                columnName='descendantExpansionFactor'
             case (rockstarColumnMvir      %ID)
                columnName='massVirial'
                propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
                     &                         /self                 %cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
             case (rockstarColumnRvir      %ID)
                columnName='radiusVirial'
                propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
                     &                         *                                           expansionFactor                    &
                     &                         /self                 %cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH) &
                     &                         /kilo
             case (rockstarColumnSpin      %ID)
                columnName='spin'
             case (rockstarColumnrs        %ID)
                columnName='radiusScale'
                propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
                     &                         *                                           expansionFactor                    &
                     &                         /self                 %cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH) &
                     &                         /kilo
             case (rockstarColumnVmax      %ID)
                columnName='velocityMaximum'
             case (rockstarColumnTU        %ID)
                columnName='virialRatio'
                propertiesReal(jReal)%property=+propertiesReal(jReal)                     %property                           &
                     &                         *2.0d0
             case (rockstarColumnRVmax     %ID)
                columnName='radiusVelocityMaximum'
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
    if (allocated(expansionFactor)) deallocate(expansionFactor)
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
