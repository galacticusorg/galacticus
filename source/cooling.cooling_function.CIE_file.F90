!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% Implements a cooling function class which interpolates in a tabulated cooling function read from file.
  
  use FGSL

  !# <coolingFunction name="coolingFunctionCIEFile" defaultThreadPrivate="yes">
  !#  <description>
  !#   Class providing a cooling function interpolated from a table read from file.  The XML file containing the table should have the following form:
  !#   \begin{verbatim}
  !#    &lt;coolingFunctions&gt;
  !#    &lt;coolingFunction&gt;
  !#      &lt;temperature&gt;
  !#        &lt;datum&gt;10000.0&lt;/datum&gt;
  !#        &lt;datum&gt;15000.0&lt;/datum&gt;
  !#        .
  !#        .
  !#        .
  !#      &lt;/temperature&gt;
  !#      &lt;coolingRate&gt;
  !#        &lt;datum&gt;1.0e-23&lt;/datum&gt;
  !#        &lt;datum&gt;1.7e-23&lt;/datum&gt;
  !#        .
  !#        .
  !#        .
  !#      &lt;/coolingRate&gt;
  !#      &lt;metallicity&gt;-4.0&lt;/metallicity&gt;
  !#    &lt;/coolingFunction&gt;
  !#    &lt;coolingFunction&gt;
  !#    .
  !#    .
  !#    .
  !#    &lt;/coolingFunction&gt;
  !#    &lt;description&gt;Some description of what this cooling function is.&lt;/description&gt;
  !#    &lt;extrapolation&gt;
  !#     &lt;metallicity&gt;
  !#       &lt;limit&gt;low&lt;/limit&gt;
  !#       &lt;method&gt;powerLaw&lt;/method&gt;
  !#     &lt;/metallicity&gt;
  !#     &lt;metallicity&gt;
  !#       &lt;limit&gt;high&lt;/limit&gt;
  !#       &lt;method&gt;powerLaw&lt;/method&gt;
  !#     &lt;/metallicity&gt;
  !#     &lt;temperature&gt;
  !#       &lt;limit&gt;low&lt;/limit&gt;
  !#       &lt;method&gt;powerLaw&lt;/method&gt;
  !#    &lt;/temperature&gt;
  !#    &lt;temperature&gt;
  !#       &lt;limit&gt;high&lt;/limit&gt;
  !#       &lt;method&gt;powerLaw&lt;/method&gt;
  !#     &lt;/temperature&gt;
  !#   &lt;/extrapolation&gt;
  !#   &lt;/coolingFunctions&gt;
  !#   \end{verbatim}
  !#   Each {\normalfont \ttfamily coolingFunction} element should contain two lists (inside
  !#   {\normalfont \ttfamily temperature} and {\normalfont \ttfamily coolingRate} tags) of
  !#   {\normalfont \ttfamily datum} elements which specify temperature (in Kelvin) and cooling
  !#   function (in ergs cm$^3$ s$^{-1}$ computed for a hydrogen density of 1 cm$^{-3}$)
  !#   respectively, and a {\normalfont \ttfamily metallicity} element which gives the
  !#   logarithmic metallcity relative to Solar (a value of -999 or less is taken to imply zero
  !#   metallicity). Any number of {\normalfont \ttfamily coolingFunction} elements may appear,
  !#   but they must be in order of increasing metallicity and must all contain the same set of
  !#   temperatures. The {\normalfont \ttfamily extrapolation} element defines how the table is
  !#   to be extrapolated in the {\normalfont \ttfamily low} and {\normalfont \ttfamily high}
  !#   limits of {\normalfont \ttfamily temperature} and {\normalfont \ttfamily
  !#   metallicity}. The {\normalfont \ttfamily method} elements can take the following values:
  !#   \begin{description}
  !#    \item[{\normalfont \ttfamily zero}] The cooling function is set to zero beyond the relevant limit.
  !#    \item[{\normalfont \ttfamily fixed}] The cooling function is held fixed at the value at the relevant limit.
  !#    \item[{\normalfont \ttfamily powerLaw}] The cooling function is extrapolated assuming a
  !#    power-law dependence beyond the relevant limit. This option is only allowed if the
  !#    cooling function is everywhere positive.  
  !#   \end{description}
  !#   If the cooling function is everywhere positive the interpolation will be done in the
  !#   logarithmic of temperature, metallicity\footnote{The exception is if the first cooling
  !#   function is tabulated for zero metallicity. In that case, a linear interpolation in
  !#   metallicity is always used between zero and the first non-zero tabulated metallicity.}
  !#   and cooling function. Otherwise, interpolation is linear in these quantities. The cooling
  !#   function is scaled assuming a quadratic dependence on hydrogen density.
  !#  </description>
  !# </coolingFunction>
  type, extends(coolingFunctionClass) :: coolingFunctionCIEFile
     !% A cooling function class which interpolates in a tabulated cooling function read from file.
     private
     type            (varying_string   )                              :: fileName
     double precision                                                 :: metallicityMaximum        , metallicityMinimum          , &
          &                                                              temperatureMaximum        , temperatureMinimum
     integer                                                          :: extrapolateMetallicityHigh, extrapolateMetallicityLow   , &
          &                                                              extrapolateTemperatureHigh, extrapolateTemperatureLow
     logical                                                          :: firstMetallicityIsZero    , logarithmicTable
     integer                                                          :: metallicityCount          , temperatureCount
     double precision                                                 :: firstNonZeroMetallicity
     double precision                   , allocatable, dimension(:)   :: metallicities             , temperatures
     double precision                   , allocatable, dimension(:,:) :: coolingFunctionTable
     logical                                                          :: resetMetallicity          , resetTemperature                      
     type            (fgsl_interp_accel)                              :: acceleratorMetallicity    , acceleratorTemperature
     double precision                                                 :: temperaturePrevious       , metallicityPrevious         , &
          &                                                              temperatureSlopePrevious  , metallicitySlopePrevious    , &
          &                                                              coolingFunctionPrevious   , coolingFunctionSlopePrevious
   contains
     !@ <objectMethods>
     !@   <object>coolingFunctionCIEFile</object>
     !@   <objectMethod>
     !@     <method>readFile</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless char(len=*)\textgreater} fileName\argin</arguments>
     !@     <description>Read the named cooling function file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>interpolatingFactors</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ temperature\argin, \doublezero\ metallicity\argin, \textcolor{red}{\textless integer(c\_size\_t)\textgreater} iTemperature\argout, \doublezero\ hTemperature\argout, \textcolor{red}{\textless integer(c\_size\_t)\textgreater} iMetallicity\argout, \doublezero\ hMetallicity\argout</arguments>
     !@     <description>Compute interpolating factors in a CIE cooling function file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>interpolate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless integer(c\_size\_t)\textgreater} iTemperature\argin, \doublezero\ hTemperature\argin, \textcolor{red}{\textless integer(c\_size\_t)\textgreater} iMetallicity\argin, \doublezero\ hMetallicity\argin</arguments>
     !@     <description>Interpolate in the cooling function.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                       cieFileDestructor
     procedure :: readFile                           => cieFileReadFile
     procedure :: interpolatingFactors               => cieFileInterpolatingFactors
     procedure :: interpolate                        => cieFileInterpolate
     procedure :: coolingFunction                    => cieFileCoolingFunction
     procedure :: coolingFunctionTemperatureLogSlope => cieFileCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => cieFileCoolingFunctionDensityLogSlope
     procedure :: descriptor                         => cieFileDescriptor
  end type coolingFunctionCIEFile

  interface coolingFunctionCIEFile
     !% Constructors for the ``CIE file'' cooling function class.
     module procedure cieFileConstructorParameters
     module procedure cieFileConstructorInternal
  end interface coolingFunctionCIEFile

  ! Current file format version for CIE cooling function files.
  integer, parameter :: cieFileFormatVersionCurrent=1

  ! Initialization status and default parameters.
  logical                 :: cieFileInitialized         =.false.
  type   (varying_string) :: cieFileFileName

contains

  function cieFileConstructorParameters(parameters)
    !% Constructor for the ``CIE file'' cooling function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(coolingFunctionCIEFile)                :: cieFileConstructorParameters
    type(inputParameters       ), intent(in   ) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    if (.not.cieFileInitialized) then
       !$omp critical(cieFileInitialize)
       if (.not.cieFileInitialized) then
          call parameters%checkParameters(allowedParameterNames)    
          !# <inputParameter>
          !#   <name>fileName</name>
          !#   <source>parameters</source>
          !#   <variable>cieFileFileName</variable>
          !#   <description>The name of the file containing a tabulation of the collisional ionization equilibrium cooling function.</description>
          !#   <type>string</type>
          !#   <cardinality>1</cardinality>
          !# </inputParameter>
          cieFileInitialized=.true.
       end if
       !$omp end critical(cieFileInitialize)
    end if
    ! Construct the instance.    
    cieFileConstructorParameters=cieFileConstructorInternal(char(cieFileFileName))
    return
  end function cieFileConstructorParameters
  
  function cieFileConstructorInternal(fileName)
    !% Internal constructor for the ``CIE file'' cooling function class.
    implicit none
    type     (coolingFunctionCIEFile)                :: cieFileConstructorInternal
    character(len=*                 ), intent(in   ) :: fileName
    
    ! Read the file.
    cieFileConstructorInternal%fileName=fileName
    call cieFileConstructorInternal%readFile(fileName)
    !initialize
    cieFileConstructorInternal%temperaturePrevious     =-1.0d0
    cieFileConstructorInternal%metallicityPrevious     =-1.0d0
    cieFileConstructorInternal%temperatureSlopePrevious=-1.0d0
    cieFileConstructorInternal%metallicitySlopePrevious=-1.0d0
    cieFileConstructorInternal%resetMetallicity        =.true.
    cieFileConstructorInternal%resetTemperature        =.true.
    return
  end function cieFileConstructorInternal
  
  subroutine cieFileDestructor(self)
    !% Destructor for the ``CIE file'' cooling function class.
    use Numerical_Interpolation
    implicit none
    type(coolingFunctionCIEFile), intent(inout) :: self

    ! Free all FGSL objects.
    call Interpolate_Done(                                                      &
         &                interpolationAccelerator=self%acceleratorMetallicity, &
         &                reset                   =self%      resetMetallicity  &
         &               )
    call Interpolate_Done(                                                      &
         &                interpolationAccelerator=self%acceleratorTemperature, &
         &                reset                   =self%      resetTemperature  &
         &               )
    return
  end subroutine cieFileDestructor

  double precision function cieFileCoolingFunction(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the cooling function by interpolating in tabulated CIE data read from a file.
    use, intrinsic :: ISO_C_Binding
    use               Abundances_Structure
    use               Radiation_Structure
    use               Chemical_Abundances_Structure
    use               Table_Labels
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    double precision                        , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances            ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances    ), intent(in   ) :: chemicalDensities
    type            (radiationStructure    ), intent(in   ) :: radiation
    integer         (c_size_t              )                :: iMetallicity         , iTemperature
    double precision                                        :: hMetallicity         , hTemperature  , &
         &                                                     metallicityUse       , temperatureUse

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow)
       case (extrapolationTypeZero)
          cieFileCoolingFunction=0.0d0
          return
       case (extrapolationTypeFix,extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh)
       case (extrapolationTypeZero)
          cieFileCoolingFunction=0.0d0
          return
       case (extrapolationTypeFix,extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMaximum
       end select
    end if
    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(gasAbundances,metallicityType=metallicityTypeLinearByMassSolar)
    if (metallicityUse < self%metallicityMinimum) then
       select case (self%extrapolateMetallicityLow)
       case (extrapolationTypeZero)
          cieFileCoolingFunction=0.0d0
          return
       case (extrapolationTypeFix)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh)
       case (extrapolationTypeZero)
          cieFileCoolingFunction=0.0d0
          return
       case (extrapolationTypeFix)
          metallicityUse=self%metallicityMaximum
       end select
    end if
    ! Check if we need to recompute the cooling function.
    if     (                                            &
         &   temperatureUse /= self%temperaturePrevious &
         &  .or.                                        &
         &   metallicityUse /= self%metallicityPrevious &
         & ) then
       ! Get the interpolation.
       call self%interpolatingFactors(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Do the interpolation.
       self%coolingFunctionPrevious=self%interpolate(iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Store the temperature and metallicity for which calculation was performed.
       self%temperaturePrevious=temperatureUse
       self%metallicityPrevious=metallicityUse
    end if
    ! Scale to the specified density assuming all processes are proportional to hydrogen density squared.
    cieFileCoolingFunction=+self%coolingFunctionPrevious    &
         &                 *numberDensityHydrogen       **2
    return
  end function cieFileCoolingFunction

  double precision function cieFileCoolingFunctionTemperatureLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the slope of the cooling function with respect to temperature by interpolating in tabulated CIE data
    !% read from a file.
    use, intrinsic :: ISO_C_Binding
    use               Abundances_Structure
    use               Chemical_Abundances_Structure
    use               Radiation_Structure
    use               Table_Labels
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    double precision                        , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances            ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances    ), intent(in   ) :: chemicalDensities
    type            (radiationStructure    ), intent(in   ) :: radiation
    double precision                                        :: coolingFunction
    integer         (c_size_t              )                :: iMetallicity         , iTemperature
    double precision                                        :: hMetallicity         , hTemperature  , &
         &                                                     metallicityUse       , temperatureUse

    ! Get the cooling function.
    coolingFunction=self%coolingFunction(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < self%temperatureMinimum) then
       select case (self%extrapolateTemperatureLow)
       case (extrapolationTypeZero,extrapolationTypeFix)
          cieFileCoolingFunctionTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMinimum
       end select
    end if
    if (temperatureUse > self%temperatureMaximum) then
       select case (self%extrapolateTemperatureHigh)
       case (extrapolationTypeZero,extrapolationTypeFix)
          cieFileCoolingFunctionTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypePowerLaw)
          temperatureUse=self%temperatureMaximum
       end select
    end if
    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(gasAbundances,metallicityType=metallicityTypeLinearByMassSolar)
    if (metallicityUse < self%metallicityMinimum) then
       select case (self%extrapolateMetallicityLow)
       case (extrapolationTypeZero)
          cieFileCoolingFunctionTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeFix)
          metallicityUse=self%metallicityMinimum
       end select
    end if
    if (metallicityUse > self%metallicityMaximum) then
       select case (self%extrapolateMetallicityHigh)
       case (extrapolationTypeZero)
          cieFileCoolingFunctionTemperatureLogSlope=0.0d0
          return
       case (extrapolationTypeFix)
          metallicityUse=self%metallicityMaximum
       end select
    end if
    ! Check if we need to recompute the cooling function.
    if     (                                                 &
         &   temperatureUse /= self%temperatureSlopePrevious &
         &  .or.                                             &
         &   metallicityUse /= self%metallicitySlopePrevious &
         & ) then

       ! Get the interpolation.
       call self%interpolatingFactors(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Do the interpolation.
       self%coolingFunctionSlopePrevious=+(                                                            &
            &                              +(+self%coolingFunctionTable(iTemperature+1,iMetallicity  ) &
            &                                -self%coolingFunctionTable(iTemperature  ,iMetallicity  ) &
            &                               )                                                          &
            &                              *(1.0d0-hMetallicity)                                       &
            &                              +(+self%coolingFunctionTable(iTemperature+1,iMetallicity+1) &
            &                                -self%coolingFunctionTable(iTemperature  ,iMetallicity+1) &
            &                               )                                                          &
            &                              *(      hMetallicity)                                       &
            &                              )                                                           &
            &                              /(                                                          &
            &                                +self%temperatures   (iTemperature+1                    ) &
            &                                -self%temperatures   (iTemperature                      ) &
            &                               )

       ! Convert to logarithmic gradient if table was not stored logarithmically.
       if (.not.self%logarithmicTable)                                                &
            & self%coolingFunctionSlopePrevious=                                      &
            &  +self%coolingFunctionSlopePrevious                                     &
            &  *temperature                                                           &
            &  /self%interpolate(iTemperature,hTemperature,iMetallicity,hMetallicity)
       ! Store the temperature and metallicity for which calculation was performed.
       self%temperatureSlopePrevious=temperatureUse
       self%metallicitySlopePrevious=metallicityUse
    end if
    ! Return the stored value.
    cieFileCoolingFunctionTemperatureLogSlope=self%coolingFunctionSlopePrevious
    return
  end function cieFileCoolingFunctionTemperatureLogSlope

  double precision function cieFileCoolingFunctionDensityLogSlope(self,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic slope of the cooling function with respect to density.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    double precision                        , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances            ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances    ), intent(in   ) :: chemicalDensities
    type            (radiationStructure    ), intent(in   ) :: radiation

    ! Logarithmic slope is always 2 for a CIE cooling function.
    cieFileCoolingFunctionDensityLogSlope=2.0d0
    return
  end function cieFileCoolingFunctionDensityLogSlope

  subroutine cieFileReadFile(self,fileName)
    !% Read in data from a cooling function file.
    use Galacticus_Error
    use FoX_dom
    use Memory_Management
    use Numerical_Comparison
    use Galacticus_Display
    use IO_XML
    use Table_Labels
    implicit none
    class           (coolingFunctionCIEFile), intent(inout)             :: self
    character       (len=*                 ), intent(in   )             :: fileName
    double precision                        , allocatable, dimension(:) :: temperaturesReference
    type            (Node                  ), pointer                   :: doc                                  , extrapolation               , &
         &                                                                 extrapolationElement                 , metallicityElement          , &
         &                                                                 thisCoolingFunction                  , thisCoolingRate             , &
         &                                                                 thisTemperature                      , version
    type            (NodeList              ), pointer                   :: coolingFunctionList                  , metallicityExtrapolationList, &
         &                                                                 temperatureExtrapolationList
    double precision                        , parameter                 :: metallicityLogarithmicZero  =-999.0d0
    integer                                                             :: extrapolationMethod                  , fileFormatVersion           , &
         &                                                                 iCoolingFunction                     , iExtrapolation              , &
         &                                                                 ioErr
    character       (len=32                )                            :: limitType
    
    !$omp critical (FoX_DOM_Access)
    ! Parse the XML file.
    call Galacticus_Display_Indent('Parsing file: '//fileName,verbosityWorking)
    call Galacticus_Display_Counter(0,.true.,verbosityWorking)
    doc => parseFile(fileName,iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('cieFileReadFile','unable to find cooling function file')
    ! Check the file format version of the file.
    version => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
    call extractDataContent(version,fileFormatVersion)
    if (fileFormatVersion /= cieFileFormatVersionCurrent) call Galacticus_Error_Report('cieFileReadFile','file format version is out of date')
    call Galacticus_Display_Counter(50,.false.,verbosityWorking)
    ! Get a list of all <coolingFunction> elements.
    coolingFunctionList   => getElementsByTagname(doc,"coolingFunction")
    self%metallicityCount =getLength(coolingFunctionList)
    ! Extract data from first cooling function and count number of temperatures present.
    thisCoolingFunction   => item(coolingFunctionList,0)
    thisTemperature       => XML_Get_First_Element_By_Tag_Name(thisCoolingFunction,"temperature")
    self%temperatureCount =  XML_Array_Length                 (thisTemperature    ,"datum"      )
    ! Allocate space for the table.
    if (allocated(self%metallicities       )) call Dealloc_Array(self%metallicities       )
    if (allocated(self%temperatures        )) call Dealloc_Array(self%temperatures        )
    if (allocated(self%coolingFunctionTable)) call Dealloc_Array(self%coolingFunctionTable)
    call Alloc_Array(self%metallicities       ,[                      self%metallicityCount])
    call Alloc_Array(self%temperatures        ,[self%temperatureCount                      ])
    call Alloc_Array(self%coolingFunctionTable,[self%temperatureCount,self%metallicityCount])
    ! Extract data from the cooling functions and populate metallicity and temperature arrays.
    do iCoolingFunction=0,self%metallicityCount-1
       ! Get required cooling function.
       thisCoolingFunction => item(coolingFunctionList,iCoolingFunction)
       ! Extract the metallicity from the <metallicity> element.
       metallicityElement  => XML_Get_First_Element_By_Tag_Name(thisCoolingFunction,"metallicity")
       call extractDataContent(metallicityElement,self%metallicities(iCoolingFunction+1))
       ! Extract the data.
       thisTemperature => XML_Get_First_Element_By_Tag_Name(thisCoolingFunction,"temperature")
       thisCoolingRate => XML_Get_First_Element_By_Tag_Name(thisCoolingFunction,"coolingRate")
       ! Check that number of temperatures is consistent.
       if (XML_Array_Length(thisTemperature,"datum") /= self%temperatureCount) &
            & call Galacticus_Error_Report('cieFileReadFile','sizes of temperature grids must be the same for all metallicities')
       ! Check that number of cooling rates matches number of temperatures.
       if (XML_Array_Length(thisTemperature,"datum") /= XML_Array_Length(thisCoolingRate,"datum")) &
            & call Galacticus_Error_Report('cieFileReadFile','sizes of temperature and cooling rate arrays must match')
       ! Extract data.
       call XML_Array_Read_Static(thisTemperature,"datum",self%temperatures                              )
       call XML_Array_Read_Static(thisCoolingRate,"datum",self%coolingFunctionTable(:,iCoolingFunction+1))
       call Galacticus_Display_Counter(                                                   &
            &                           int(                                              &
            &                                50.0d0                                       &
            &                               +50.0d0                                       &
            &                               *dble(iCoolingFunction     )                  &
            &                               /dble(self%metallicityCount)                  &
            &                              )                                              &
            &                          ,.false.                                           &
            &                          ,verbosityWorking                                  &
            &                         )
       if (iCoolingFunction == 0) then
          ! Make a copy of the temperatures to use as a reference for future temperature reads.
          temperaturesReference=self%temperatures
       else
          ! Check that temperature grids are aligned.
          if (any(Values_Differ(self%temperatures,temperaturesReference,relTol=1.0d-6))) &
               & call Galacticus_Error_Report('cieFileReadFile','temperature grids mismatch')
       end if
    end do
    deallocate(temperaturesReference)
    where (self%metallicities > metallicityLogarithmicZero)
       self%metallicities=10.0d0**self%metallicities
    elsewhere
       self%metallicities=0.0d0
    end where
    ! Extract extrapolation methods from the file.
    extrapolationElement         => XML_Get_First_Element_By_Tag_Name(doc,"extrapolation")
    metallicityExtrapolationList => getElementsByTagname(extrapolationElement,"metallicity")
    do iExtrapolation=0,getLength(metallicityExtrapolationList)-1
       extrapolation => item(metallicityExtrapolationList,iExtrapolation)
       call XML_Extrapolation_Element_Decode(extrapolation,limitType,extrapolationMethod,allowedMethods=[extrapolationTypeZero,extrapolationTypeFix,extrapolationTypePowerLaw])
       select case (trim(limitType))
       case ('low')
          self%extrapolateMetallicityLow=extrapolationMethod
       case ('high')
          self%extrapolateMetallicityHigh=extrapolationMethod
       case default
         call Galacticus_Error_Report('cieFileReadFile','unrecognized extrapolation limit')
       end select
    end do
    temperatureExtrapolationList => getElementsByTagname(extrapolationElement,"temperature")
    do iExtrapolation=0,getLength(temperatureExtrapolationList)-1
       extrapolation => item(temperatureExtrapolationList,iExtrapolation)
       call XML_Extrapolation_Element_Decode(extrapolation,limitType,extrapolationMethod,allowedMethods=[extrapolationTypeZero,extrapolationTypeFix,extrapolationTypePowerLaw])
       select case (trim(limitType))
       case ('low')
          self%extrapolateTemperatureLow=extrapolationMethod
       case ('high')
          self%extrapolateTemperatureHigh=extrapolationMethod
       case default
          call Galacticus_Error_Report('cieFileReadFile','unrecognized extrapolation limit')
       end select
    end do
    ! Destroy the document.
    call destroy(doc)
    call Galacticus_Display_Counter_Clear(       verbosityWorking)
    call Galacticus_Display_Unindent     ('done',verbosityWorking)
    !$omp end critical (FoX_DOM_Access)
    ! Store table ranges for convenience.
    self%metallicityMinimum=self%metallicities(                    1)
    self%metallicityMaximum=self%metallicities(self%metallicityCount)
    self%temperatureMinimum=self%temperatures (                    1)
    self%temperatureMaximum=self%temperatures (self%temperatureCount)
    ! Decide whether or not to make the tables logarithmic.
    self%logarithmicTable=all(self%coolingFunctionTable > 0.0d0)
    if (self%logarithmicTable) then
       self%firstMetallicityIsZero=(self%metallicities(1) == 0.0d0)
       if (self%firstMetallicityIsZero) self%firstNonZeroMetallicity=self%metallicities(2)
       where (self%metallicities > 0.0d0)
          self%metallicities=log(self%metallicities)
       elsewhere
          self%metallicities=metallicityLogarithmicZero
       end where
       self%temperatures        =log(self%temperatures        )
       self%coolingFunctionTable=log(self%coolingFunctionTable)
    else
       if     (                                                             &
            &  self%extrapolateTemperatureLow  == extrapolationTypePowerLaw &
            &  .or.                                                         &
            &  self%extrapolateTemperatureHigh == extrapolationTypePowerLaw &
            & )                                                             &
            & call Galacticus_Error_Report('cieFileReadFile','power law extrapolation allowed only in loggable tables')
    end if
    if     (                                                             &
         &  self%extrapolateMetallicityLow  == extrapolationTypePowerLaw &
         &   .or.                                                        &
         &  self%extrapolateMetallicityHigh == extrapolationTypePowerLaw &
         & )                                                             &
         & call Galacticus_Error_Report('cieFileReadFile','power law extrapolation not allowed in metallicity')
    ! Force interpolation accelerators to be reset.
    self%resetTemperature=.true.
    self%resetMetallicity=.true.
    return
  end subroutine cieFileReadFile

  subroutine cieFileInterpolatingFactors(self,temperature,metallicity,iTemperature,hTemperature,iMetallicity,hMetallicity)
    !% Determine the interpolating paramters.
    use, intrinsic :: ISO_C_Binding
    use               Numerical_Interpolation
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    double precision                        , intent(in   ) :: metallicity , temperature
    integer         (c_size_t              ), intent(  out) :: iMetallicity  , iTemperature
    double precision                        , intent(  out) :: hMetallicity  , hTemperature
    double precision                                        :: metallicityUse, temperatureUse

    ! Copy the input parameters.
    temperatureUse=    temperature
    metallicityUse=max(metallicity,0.0d0)
    ! Get interpolation in temperature.
    if (self%logarithmicTable) temperatureUse=log(temperatureUse)
    iTemperature=max(                                                    &
         &           min(                                                &
         &               Interpolate_Locate(                             &
         &                                  self%temperatures          , &
         &                                  self%acceleratorTemperature, &
         &                                  temperatureUse             , &
         &                                  self%resetTemperature        &
         &                                 )                           , &
         &               self%temperatureCount-1                         &
         &              )                                              , &
         &               1                                               &
         &          )
    hTemperature=+(     temperatureUse                -self%temperatures(iTemperature)) &
         &       /(self%temperatures  (iTemperature+1)-self%temperatures(iTemperature))
    ! Get interpolation in metallicity.
    if (self%firstMetallicityIsZero .and. metallicityUse < self%firstNonZeroMetallicity) then
       iMetallicity=1
       hMetallicity=metallicityUse/self%firstNonZeroMetallicity
    else
       if (self%logarithmicTable) metallicityUse=log(metallicityUse)
       iMetallicity=max(                                                    &
            &           min(                                                &
            &               Interpolate_Locate(                             &
            &                                  self%metallicities         , &
            &                                  self%acceleratorMetallicity, &
            &                                  metallicityUse             , &
            &                                  self%resetMetallicity        &
            &                                 )                           , &
            &               self%metallicityCount-1                         &
            &              )                                              , &
            &               1                                               &
            &          )
       hMetallicity=+(     metallicityUse                -self%metallicities(iMetallicity)) &
            &       /(self%metallicities (iMetallicity+1)-self%metallicities(iMetallicity))
    end if
    return
  end subroutine cieFileInterpolatingFactors

  double precision function cieFileInterpolate(self,iTemperature,hTemperature,iMetallicity,hMetallicity)
    !% Perform the interpolation.
    use, intrinsic :: ISO_C_Binding
    implicit none
    class           (coolingFunctionCIEFile), intent(inout) :: self
    integer         (c_size_t              ), intent(in   ) :: iMetallicity, iTemperature
    double precision                        , intent(in   ) :: hMetallicity, hTemperature

    ! Do the interpolation.
    cieFileInterpolate=+self%coolingFunctionTable(iTemperature  ,iMetallicity  )*(1.0d0-hTemperature)*(1.0d0-hMetallicity) &
         &             +self%coolingFunctionTable(iTemperature  ,iMetallicity+1)*(1.0d0-hTemperature)*(      hMetallicity) &
         &             +self%coolingFunctionTable(iTemperature+1,iMetallicity  )*(      hTemperature)*(1.0d0-hMetallicity) &
         &             +self%coolingFunctionTable(iTemperature+1,iMetallicity+1)*(      hTemperature)*(      hMetallicity)
    ! Exponentiate the result if the table was stored as the log.
    if (self%logarithmicTable) cieFileInterpolate=exp(cieFileInterpolate)
    return
  end function cieFileInterpolate

  subroutine cieFileDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class(coolingFunctionCIEFile), intent(inout) :: self
    type (inputParameters       ), intent(inout) :: descriptor
    type (node                  ), pointer       :: parameterNode
    type (inputParameters       )                :: subParameters

    call descriptor%addParameter("coolingFunctionMethod","cieFile")
    parameterNode => descriptor%node("coolingFunctionMethod")
    subParameters=inputParameters(parameterNode)
    call subParameters%addParameter("fileName",char(self%fileName))
    return
  end subroutine cieFileDescriptor
