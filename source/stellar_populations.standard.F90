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
  Implements a standard stellar population class.
  !!}

  use :: Numerical_Interpolation                   , only : interpolator
  use :: Stellar_Astrophysics                      , only : stellarAstrophysics     , stellarAstrophysicsClass
  use :: Stellar_Feedback                          , only : stellarFeedback         , stellarFeedbackClass
  use :: Stellar_Population_Spectra                , only : stellarPopulationSpectra, stellarPopulationSpectraClass
  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunction     , initialMassFunctionClass
  use :: Supernovae_Type_Ia                        , only : supernovaeTypeIa        , supernovaeTypeIaClass

  abstract interface
     !!{
     Interface for stellar population property integrands.
     !!}
     double precision function integrandTemplate(massInitial)
       double precision, intent(in   ) :: massInitial
     end function integrandTemplate
  end interface

  type :: populationTable
     !!{
     Table used to store properties of the stellar population as a function of age and metallicity.
     !!}
     type            (varying_string   )                              :: label
     procedure       (integrandTemplate), pointer    , nopass         :: integrand
     double precision                   , allocatable, dimension(:  ) :: age              , metallicity
     double precision                   , allocatable, dimension(:,:) :: property
     double precision                                                 :: toleranceAbsolute, toleranceRelative
     type            (interpolator     )                              :: interpolatorAge  , interpolatorMetallicity
     logical                                                          :: computed         , instantaneousApproximation
  end type populationTable

  interface populationTable
     !!{
     Constructors for the {\normalfont \ttfamily populationTable} class.
     !!}
     module procedure populationTableConstructor
  end interface populationTable

  !![
  <stellarPopulation name="stellarPopulationStandard">
   <description>A standard stellar population class.</description>
  </stellarPopulation>
  !!]
  type, extends(stellarPopulationClass) :: stellarPopulationStandard
     !!{
     A standard stellar population class.
     !!}
     private
     logical                                                                    :: instantaneousRecyclingApproximation          , instantaneousYieldApproximation, &
         &                                                                         instantaneousEnergyInputApproximation        , metalYieldInitialized          , &
          &                                                                        recycledFractionInitialized
     double precision                                                           :: massLongLived                                , ageEffective                   , &
          &                                                                        recycledFraction_                            , metalYield_                    , &
          &                                                                        recycledFraction                             , metalYield
     class           (stellarAstrophysicsClass     ), pointer                   :: stellarAstrophysics_                => null()
     class           (initialMassFunctionClass     ), pointer                   :: initialMassFunction_                => null()
     class           (stellarFeedbackClass         ), pointer                   :: stellarFeedback_                    => null()
     class           (supernovaeTypeIaClass        ), pointer                   :: supernovaeTypeIa_                   => null()
     class           (stellarPopulationSpectraClass), pointer                   :: stellarPopulationSpectra_           => null()
     type            (populationTable              )                            :: recycleFraction
     type            (populationTable              )                            :: energyOutput
     type            (populationTable              ), allocatable, dimension(:) :: yield
   contains
     !![
     <methods>
       <method description="Return true if the star of given initial mass and metallicity has evolved off of the main sequence by the given age." method="starIsEvolved" />
       <method description="Interpolate in the given property to return the mean rate of production of that property from the stellar population between the given minimum and maximum ages." method="interpolate" />
     </methods>
     !!]
     final     ::                                  standardDestructor
     procedure :: rateRecycling                 => standardRateRecycling
     procedure :: rateYield                     => standardRateYield
     procedure :: rateEnergy                    => standardRateEnergy
     procedure :: recycledFractionInstantaneous => standardRecycledFractionInstantaneous
     procedure :: yieldInstantaneous            => standardYieldInstantaneous
     procedure :: starIsEvolved                 => standardStarIsEvolved
     procedure :: interpolate                   => standardInterpolate
     procedure :: spectra                       => standardSpectra
  end type stellarPopulationStandard

  interface stellarPopulationStandard
     !!{
     Constructors for the \refClass{stellarPopulationStandard} stellar population class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface stellarPopulationStandard

  ! Module-scope variables used in integrations.
  class           (stellarPopulationStandard), pointer   :: self_
  double precision                                       :: lifetime_                     , metallicity_
  integer                                                :: indexElement_
  class           (stellarAstrophysicsClass ), pointer   :: stellarAstrophysics_
  class           (initialMassFunctionClass ), pointer   :: initialMassFunction_
  class           (stellarFeedbackClass     ), pointer   :: stellarFeedback_
  class           (supernovaeTypeIaClass    ), pointer   :: supernovaeTypeIa_
  logical                                                :: instantaneousApproximation
  !$omp threadprivate(indexElement_,self_,lifetime_,metallicity_,stellarAstrophysics_,initialMassFunction_,stellarFeedback_,supernovaeTypeIa_,instantaneousApproximation)

  ! Tabulation resolution.
  integer                                    , parameter :: tableMetallicityCount  =10
  integer                                    , parameter :: tableAgeCount          =50
  double precision                           , parameter :: tableMetallicityMinimum=1.0d-4
  double precision                           , parameter :: tableMetallicityMaximum=0.6d-1
  double precision                           , parameter :: tableAgeMinimum        =1.0d-3
  double precision                           , parameter :: tableAgeMaximum        =1.0d+2

  ! File format version.
  integer                                    , parameter :: fileFormatCurrent      =1

contains

  function populationTableConstructor(label,integrand,toleranceAbsolute,toleranceRelative,instantaneousApproximation) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationStandard} stellar population class which takes a parameter list as input.
    !!}
    implicit none
    type            (populationTable)                :: self
    character       (len=*          ), intent(in   ) :: label
    double precision                 , external      :: integrand
    double precision                 , intent(in   ) :: toleranceAbsolute         , toleranceRelative
    logical                          , intent(in   ) :: instantaneousApproximation
    !![
    <constructorAssign variables="label, *integrand, toleranceAbsolute, toleranceRelative, instantaneousApproximation"/>
    !!]

    self%computed=.false.
    return
  end function populationTableConstructor

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationStandard} stellar population class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarPopulationStandard    )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (initialMassFunctionClass     ), pointer       :: initialMassFunction_
    class           (stellarAstrophysicsClass     ), pointer       :: stellarAstrophysics_
    class           (stellarFeedbackClass         ), pointer       :: stellarFeedback_
    class           (supernovaeTypeIaClass        ), pointer       :: supernovaeTypeIa_
    class           (stellarPopulationSpectraClass), pointer       :: stellarPopulationSpectra_
    logical                                                        :: instantaneousRecyclingApproximation  , instantaneousYieldApproximation, &
         &                                                            instantaneousEnergyInputApproximation
    double precision                                               :: massLongLived                        , ageEffective                   , &
         &                                                            recycledFraction                     , metalYield
    !$GLC attributes initialized :: self

    !![
    <inputParameter>
      <name>instantaneousRecyclingApproximation</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, then use an instantaneous recycling approximation when computing recycling rates.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>instantaneousYieldApproximation</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, then use an instantaneous recycling approximation when computing yield rates.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>instantaneousEnergyInputApproximation</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, then use an instantaneous recycling approximation when computing energy input rates.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massLongLived</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass below which stars are assumed to be infinitely long-lived in the instantaneous approximation for stellar evolution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ageEffective</name>
      <defaultValue>13.8d0</defaultValue>
      <description>The effective age to use for computing SNeIa yield when using the instantaneous stellar evolution approximation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>recycledFraction</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The recycled fraction to use in the instantaneous stellar evolution approximation. (If not specified it will be computed internally.)</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>metalYield</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The metal yield to use in the instantaneous stellar evolution approximation. (If not specified it will be computed internally.)</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="initialMassFunction"      name="initialMassFunction_"      source="parameters"/>
    <objectBuilder class="stellarAstrophysics"      name="stellarAstrophysics_"      source="parameters"/>
    <objectBuilder class="stellarFeedback"          name="stellarFeedback_"          source="parameters"/>
    <objectBuilder class="supernovaeTypeIa"         name="supernovaeTypeIa_"         source="parameters"/>
    <objectBuilder class="stellarPopulationSpectra" name="stellarPopulationSpectra_" source="parameters"/>
    <conditionalCall>
    <call>
     self=stellarPopulationStandard(                                                                  &amp;
      &amp;                                                    instantaneousRecyclingApproximation  , &amp;
      &amp;                                                    instantaneousYieldApproximation      , &amp;
      &amp;                                                    instantaneousEnergyInputApproximation, &amp;
      &amp;                                                    massLongLived                        , &amp;
      &amp;                                                    ageEffective                         , &amp;
      &amp;                          initialMassFunction_     =initialMassFunction_                 , &amp;
      &amp;                          stellarAstrophysics_     =stellarAstrophysics_                 , &amp;
      &amp;                          stellarFeedback_         =stellarFeedback_                     , &amp;
      &amp;                          supernovaeTypeIa_        =supernovaeTypeIa_                    , &amp;
      &amp;                          stellarPopulationSpectra_=stellarPopulationSpectra_              &amp;
      &amp;                          {conditions}                                                     &amp;
      &amp;                         )
     </call>
     <argument name="recycledFraction" value="recycledFraction" parameterPresent="parameters"/>
     <argument name="metalYield"       value="metalYield"       parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="initialMassFunction_"     />
    <objectDestructor name="stellarAstrophysics_"     />
    <objectDestructor name="stellarFeedback_"         />
    <objectDestructor name="supernovaeTypeIa_"        />
    <objectDestructor name="stellarPopulationSpectra_"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(instantaneousRecyclingApproximation,instantaneousYieldApproximation,instantaneousEnergyInputApproximation,massLongLived,ageEffective,recycledFraction,metalYield,initialMassFunction_,stellarAstrophysics_,stellarFeedback_,supernovaeTypeIa_,stellarPopulationSpectra_) result(self)
    !!{
    Internal constructor for the \refClass{stellarPopulationStandard} stellar population.
    !!}
    use :: Abundances_Structure, only : Abundances_Names, Abundances_Property_Count
    implicit none
    type            (stellarPopulationStandard    )                          :: self
    logical                                        , intent(in   )           :: instantaneousRecyclingApproximation  , instantaneousYieldApproximation, &
         &                                                                      instantaneousEnergyInputApproximation
    double precision                               , intent(in   )           :: massLongLived                        , ageEffective
    double precision                               , intent(in   ), optional :: recycledFraction                     , metalYield
    class           (initialMassFunctionClass     ), intent(in   ), target   :: initialMassFunction_
    class           (stellarAstrophysicsClass     ), intent(in   ), target   :: stellarAstrophysics_
    class           (stellarFeedbackClass         ), intent(in   ), target   :: stellarFeedback_
    class           (supernovaeTypeIaClass        ), intent(in   ), target   :: supernovaeTypeIa_
    class           (stellarPopulationSpectraClass), intent(in   ), target   :: stellarPopulationSpectra_
    integer                                                                  :: i                                    , countElements
    !![
    <constructorAssign variables="instantaneousRecyclingApproximation, instantaneousYieldApproximation, instantaneousEnergyInputApproximation, massLongLived, ageEffective, *initialMassFunction_, *stellarAstrophysics_, *stellarFeedback_, *supernovaeTypeIa_, *stellarPopulationSpectra_"/>
    !!]

    self%recycledFractionInitialized=present(recycledFraction)
    self%metalYieldInitialized      =present(metalYield      )
    if (self%recycledFractionInitialized) then
       self%recycledFraction_=recycledFraction
       self%recycledFraction =recycledFraction
    else
       self%recycledFraction_=-huge(0.0d0)
       self%recycledFraction =-huge(0.0d0)
    end if
    if (self%metalYieldInitialized      ) then
       self%metalYield_      =metalYield
       self%metalYield       =metalYield
    else
       self%metalYield_      =-huge(0.0d0)
       self%metalYield       = huge(0.0d0)
    end if
    countElements=Abundances_Property_Count()
    self   %recycleFraction   =populationTable('recycledFraction'                ,standardIntegrandRecycledFraction,toleranceAbsolute=1.0d-3,toleranceRelative=1.0d-4,instantaneousApproximation=instantaneousRecyclingApproximation  )
    self   %energyOutput      =populationTable('energyOutput'                    ,standardIntegrandEnergyOutput    ,toleranceAbsolute=0.0d+0,toleranceRelative=1.0d-3,instantaneousApproximation=instantaneousEnergyInputApproximation)
    allocate(self%yield(countElements))
    do i=1,countElements
       self%yield          (i)=populationTable('yield'//char(Abundances_Names(i)),standardIntegrandYield           ,toleranceAbsolute=1.0d-4,toleranceRelative=1.0d-5,instantaneousApproximation=instantaneousYieldApproximation      )
    end do
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{stellarPopulationStandard} stellar population class.
    !!}
    implicit none
    type(stellarPopulationStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%initialMassFunction_"     />
    <objectDestructor name="self%stellarAstrophysics_"     />
    <objectDestructor name="self%stellarFeedback_"         />
    <objectDestructor name="self%supernovaeTypeIa_"        />
    <objectDestructor name="self%stellarPopulationSpectra_"/>
    !!]
    return
  end subroutine standardDestructor

  double precision function standardRateRecycling(self,abundances_,ageMinimum,ageMaximum)
    !!{
    Return the rate at which mass is being recycled from this stellar population. The mean recycling rate (i.e. the fraction of
    the population's mass returned to the \gls{ism} per Gyr) is computed between the given {\normalfont \ttfamily ageMinimum}
    and {\normalfont \ttfamily ageMaximum} (in Gyr).
    !!}
    implicit none
    class           (stellarPopulationStandard), intent(inout) :: self
    double precision                           , intent(in   ) :: ageMinimum , ageMaximum
    type            (abundances               ), intent(in   ) :: abundances_

    standardRateRecycling=self%interpolate(abundances_,ageMinimum,ageMaximum,self%recycleFraction)
    return
  end function standardRateRecycling

  double precision function standardRateYield(self,abundances_,ageMinimum,ageMaximum,elementIndex)
    !!{
    Return the rate at which mass is being recycled from this stellar population. The mean recycling rate (i.e. the fraction of
    the population's mass returned to the \gls{ism} per Gyr) is computed between the given {\normalfont \ttfamily ageMinimum}
    and {\normalfont \ttfamily ageMaximum} (in Gyr).
    !!}
    use :: Abundances_Structure, only : Abundances_Atomic_Index
    implicit none
    class           (stellarPopulationStandard), intent(inout)           :: self
    double precision                           , intent(in   )           :: ageMinimum  , ageMaximum
    type            (abundances               ), intent(in   )           :: abundances_
    integer                                    , intent(in   ), optional :: elementIndex
    !![
    <optionalArgument name="elementIndex" defaultsTo="1"/>
    !!]

    indexElement_    =Abundances_Atomic_Index(elementIndex_)
    standardRateYield=self%interpolate(abundances_,ageMinimum,ageMaximum,self%yield(elementIndex_))
    return
  end function standardRateYield

  double precision function standardRateEnergy(self,abundances_,ageMinimum,ageMaximum)
    !!{
    Return the rate at which energy is being output by this stellar population in (km/s)$^2$ Gyr$^{-1}$. The mean energy output
    rate per Gyr is computed between the given {\normalfont \ttfamily ageMinimum} and {\normalfont \ttfamily ageMaximum} (in
    Gyr).
    !!}
    implicit none
    class           (stellarPopulationStandard), intent(inout) :: self
    double precision                           , intent(in   ) :: ageMinimum , ageMaximum
    type            (abundances               ), intent(in   ) :: abundances_

    standardRateEnergy=self%interpolate(abundances_,ageMinimum,ageMaximum,self%energyOutput)
    return
  end function standardRateEnergy

  double precision function standardInterpolate(self,abundances_,ageMinimum,ageMaximum,property)
    !!{
    Return the rate at which at cumulative property is being produced from this stellar population. The cumulative property is
    computed on a grid of age and metallicity. This is stored to file and will be read back in on subsequent runs. This is
    useful as computation of the table is relatively slow.
    !!}
    use            :: Abundances_Structure            , only : Abundances_Get_Metallicity
    use            :: Dates_and_Times                 , only : Formatted_Date_and_Time
    use            :: Display                         , only : displayCounter                   , displayCounterClear , displayIndent, displayUnindent, &
          &                                                    verbosityLevelWorking
    use            :: File_Utilities                  , only : Directory_Make                   , File_Exists         , File_Lock    , File_Path      , &
          &                                                    File_Unlock                      , lockDescriptor
    use            :: Error                           , only : Error_Report
    use            :: Input_Paths                     , only : inputPath                        , pathTypeDataDynamic
    use            :: HDF5_Access                     , only : hdf5Access
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Input_Parameters                , only : inputParameters
    use            :: Numerical_Constants_Astronomical, only : metallicitySolar
    use            :: Numerical_Integration2          , only : integratorCompositeGaussKronrod1D
    use            :: Numerical_Ranges                , only : Make_Range                       , rangeTypeLogarithmic
    use            :: Table_Labels                    , only : extrapolationTypeExtrapolate
    implicit none
    class           (stellarPopulationStandard        ), intent(inout), target :: self
    double precision                                   , intent(in   )         :: ageMinimum        , ageMaximum
    type            (abundances                       ), intent(in   )         :: abundances_
    type            (populationTable                  ), intent(inout)         :: property
    double precision                                   , dimension(2)          :: metallicityFactors, rate            , &
         &                                                                        propertyCumulative, age
    type            (inputParameters                  ), save                  :: descriptor
    !$omp threadprivate(descriptor)
    integer         (c_size_t                         )                        :: metallicityIndex
    type            (integratorCompositeGaussKronrod1D)                        :: integrator_
    integer                                                                    :: fileFormat        , iAge            , &
         &                                                                        iMetallicity      , loopCount       , &
         &                                                                        loopCountTotal    , i
    double precision                                                           :: maximumMass       , minimumMass     , &
         &                                                                        metallicity
    type            (hdf5Object                       )                        :: file
    character       (len=20                           )                        :: progressMessage
    type            (varying_string                   )                        :: fileName          , descriptorString
    logical                                                                    :: makeFile
    type            (lockDescriptor                   )                        :: lock

    ! Compute the property if not already done.
    if (.not.property%computed) then
       ! Check for previously computed data.
       makeFile=.false.
       fileName=char(inputPath(pathTypeDataDynamic))//'stellarPopulations/'//property%label//'_'//self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)//'.hdf5'
       call Directory_Make(File_Path(fileName))
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),lock,lockIsShared=.false.)
       if (File_Exists(fileName)) then
          ! Open the file containing cumulative property data.
          call displayIndent('Reading file: '//fileName,verbosityLevelWorking)
          !$ call hdf5Access%set          (                         )
          file=hdf5Object(char(fileName))
          call file%readAttribute('fileFormat',fileFormat)
          if (fileFormat /= fileFormatCurrent) then
             makeFile=.true.
             call displayUnindent('done',verbosityLevelWorking)
          end if
          !$ call hdf5Access%unset()
       else
          makeFile=.true.
       end if
       if (.not.makeFile) then
          ! Read the cumulative property data from file.
          !$ call hdf5Access%set        (                                         )
          call    file      %readDataset("age"               ,property%age        )
          call    file      %readDataset("metallicity"       ,property%metallicity)
          call    file      %readDataset(char(property%label),property%property   )
          !$ call hdf5Access%unset      (                                         )
          call displayUnindent('done',verbosityLevelWorking)
       else
          allocate(property%age        (tableAgeCount                      ))
          allocate(property%metallicity(              tableMetallicityCount))
          allocate(property%property   (tableAgeCount,tableMetallicityCount))
          property%metallicity(1                              )=0.0d0
          property%metallicity(2:tableMetallicityCount)=Make_Range(tableMetallicityMinimum,tableMetallicityMaximum,tableMetallicityCount-1,rangeType=rangeTypeLogarithmic)
          property%age                                 =Make_Range(tableAgeMinimum        ,tableAgeMaximum        ,tableAgeCount          ,rangeType=rangeTypeLogarithmic)
          ! Open file to output the data to.
          descriptor=inputParameters()
          call self%descriptor(descriptor,includeClass=.true.)
          descriptorString=descriptor%serializeToString()
          call descriptor%destroy()
          !$ call hdf5Access%set()
          block
            type(hdf5Object) :: dataset
            file=hdf5Object(char(fileName))
            call file   %writeAttribute(fileFormatCurrent                                                        ,'fileFormat'                         )
            call file   %writeAttribute(char(property%label           )                                          ,'description'                        )
            call file   %writeAttribute('Computed by Galacticus'                                                 ,'source'                             )
            call file   %writeAttribute(char(Formatted_Date_and_Time())                                          ,'date'                               )
            call file   %writeAttribute(char(descriptorString         )                                          ,'parameters'                         )
            call file   %writeDataset  (property%age                                                             ,'age'        ,datasetReturned=dataset)
            call dataset%writeAttribute('Age of the stellar population in Gyr'                                   ,'description'                        )
            call file   %writeDataset  (property%metallicity                                                     ,'metallicity',datasetReturned=dataset)
            call dataset%writeAttribute('Metallicity (fractional mass of total metals) of the stellar population','description'                        )
          end block
          !$ call hdf5Access%unset()
          ! Loop over ages and metallicities and compute the property.
          call displayIndent('Tabulating property: '//char(property%label),verbosityLevelWorking)
          call displayCounter(0,.true.,verbosityLevelWorking)
          loopCountTotal             =  +tableMetallicityCount &
               &                        *tableAgeCount
          loopCount                  =   0
          self_                      =>  self
          instantaneousApproximation =   property%instantaneousApproximation
          !$omp parallel private (iAge,iMetallicity,progressMessage,minimumMass,maximumMass,integrator_) copyin(self_,indexElement_)
          allocate(stellarAstrophysics_,mold=self%stellarAstrophysics_)
          allocate(initialMassFunction_,mold=self%initialMassFunction_)
          allocate(stellarFeedback_    ,mold=self%stellarFeedback_    )
          allocate(supernovaeTypeIa_   ,mold=self%supernovaeTypeIa_   )
          !$omp critical(stellarPopulationsStandardDeepCopy)
          !![
          <deepCopyReset variables="self%stellarAstrophysics_ self%initialMassFunction_ self%stellarFeedback_ self%supernovaeTypeIa_"/>
          <deepCopy source="self%stellarAstrophysics_" destination="stellarAstrophysics_"/>
          <deepCopy source="self%initialMassFunction_" destination="initialMassFunction_"/>
          <deepCopy source="self%stellarFeedback_"     destination="stellarFeedback_"    />
          <deepCopy source="self%supernovaeTypeIa_"    destination="supernovaeTypeIa_"   />
          <deepCopyFinalize variables="stellarAstrophysics_ initialMassFunction_ stellarFeedback_ supernovaeTypeIa_"/>
          !!]
          !$omp end critical(stellarPopulationsStandardDeepCopy)
          call integrator_%initialize  (24                        ,61                        )
          call integrator_%toleranceSet(property%toleranceAbsolute,property%toleranceRelative)
          call integrator_%integrandSet(property%integrand                                   )
          !$omp do schedule(dynamic)
          do i=0,loopCountTotal-1
             iMetallicity=mod( i                  ,tableMetallicityCount)+1
             iAge        =    (i-(iMetallicity-1))/tableMetallicityCount +1
             lifetime_   =property%age(iAge)
             ! Set the metallicity. If using the instantaneous recycling approximation, assume Solar metallicity always.
             if (instantaneousApproximation) then
                metallicity_=metallicitySolar
             else
                metallicity_=property%metallicity(iMetallicity)
             end if
             ! Find the minimum and maximum masses to integrate over for this IMF.
             minimumMass=self%initialMassFunction_%massMinimum()
             maximumMass=self%initialMassFunction_%massMaximum()
             ! Integrate ejected mass over the IMF between these limits.
             property%property(iAge,iMetallicity)=integrator_%evaluate(minimumMass,maximumMass)
             ! Update the counter.
             !$omp atomic
             loopCount=loopCount+1
             call displayCounter(                                                   &
                  &              int(100.0d0*dble(loopCount)/dble(loopCountTotal)), &
                  &              .false.                                          , &
                  &              verbosityLevelWorking                              &
                  &             )
          end do
          !$omp end do
          !![
          <objectDestructor name="stellarAstrophysics_"/>
          <objectDestructor name="initialMassFunction_"/>
          <objectDestructor name="stellarFeedback_"    />
	  <objectDestructor name="supernovaeTypeIa_"   />
	  !!]
          !$omp end parallel
          do iAge=1,tableAgeCount
             do iMetallicity=1,tableMetallicityCount
                ! Enforce monotonicity in the cumulative property. Non-monotonicity can arise due to the vagaries of interpolating
                ! stellar lifetimes in an irregular grid of stellar models.
                if (iAge > 1 )                                        &
                     & property%property       (iAge  ,iMetallicity)  &
                     &   =max(                                        &
                     &        property%property(iAge  ,iMetallicity), &
                     &        property%property(iAge-1,iMetallicity)  &
                     &       )
             end do
          end do
          call displayCounterClear(           verbosityLevelWorking)
          call displayUnindent     ('finished',verbosityLevelWorking)
          !$ call hdf5Access%set         (                                      )
          call    file      %writeDataset(property%property,char(property%label))
          !$ call hdf5Access%unset       (                                      )
          call displayIndent('Storing to file: '//fileName,verbosityLevelWorking)
       end if
       call File_Unlock(lock)
       ! Build interpolators.
       property%interpolatorMetallicity=interpolator(property%metallicity                                               )
       property%interpolatorAge        =interpolator(property%age        ,extrapolationType=extrapolationTypeExtrapolate)
       ! Record that this IMF has now been tabulated.
       property%computed=.true.
    end if
    ! Interpolate to get the derivative in the property at two adjacent metallicities.
    metallicity=max(Abundances_Get_Metallicity(abundances_),0.0d0)
    if (metallicity > property%metallicity(size(property%metallicity))) then
       metallicityIndex=size(property%metallicity)
       metallicityFactors=[1.0d0,0.0d0]
    else
       call property%interpolatorMetallicity%linearFactors(metallicity,metallicityIndex,metallicityFactors)
    end if
    ! Interpolate in age at both metallicities.
    age=[ageMinimum,ageMaximum]
    do iMetallicity=0,1
       if (metallicityFactors(iMetallicity+1) /= 0.0d0) then
          ! Get average property between ageMinimum and ageMaximum.
          do iAge=1,2
             if (age(iAge) > 0.0d0) then
                propertyCumulative(iAge)=property%interpolatorAge%interpolate(age(iAge),property%property(:,metallicityIndex+iMetallicity))
             else
                propertyCumulative(iAge)=0.0d0
             end if
          end do
          rate(iMetallicity+1)=+(+propertyCumulative(2)-propertyCumulative(1)) &
               &               /(+               age(2)-               age(1))
       else
          rate(iMetallicity+1)=0.0d0
       end if
    end do
    ! Interpolate in metallicity to get the actual rate.
    standardInterpolate=sum(metallicityFactors*rate)
    return
  end function standardInterpolate

  double precision function standardIntegrandRecycledFraction(massInitial)
    !!{
    Integrand used in evaluating recycled fractions.
    !!}
    implicit none
    double precision, intent(in   ) :: massInitial

    if (self_%starIsEvolved(massInitial,metallicity_,lifetime_)) then
       standardIntegrandRecycledFraction=+initialMassFunction_%phi        (massInitial             ) &
            &                            *stellarAstrophysics_%massEjected(massInitial,metallicity_)
    else
       standardIntegrandRecycledFraction=+0.0d0
    end if
    return
  end function standardIntegrandRecycledFraction

  double precision function standardIntegrandYield(massInitial)
    !!{
    Integrand used in evaluating metal yields.
    !!}
    implicit none
    double precision, intent(in   ) :: massInitial
    double precision                :: sneIaLifetime, yieldMass

    ! Include yields from isolated stars.
    if (self_%starIsEvolved(massInitial,metallicity_,lifetime_)) then
       select case (indexElement_)
       case (0)
          ! Total metallicity required.
          yieldMass=stellarAstrophysics_%massYield(massInitial,metallicity_              )
       case default
          ! Individual element required.
          yieldMass=stellarAstrophysics_%massYield(massInitial,metallicity_,indexElement_)
       end select
       standardIntegrandYield=+initialMassFunction_%phi(massInitial) &
            &                 *yieldMass
    else
       standardIntegrandYield=0.0d0
    end if
    ! Include yield from Type Ia supernovae.
    if (self_%instantaneousYieldApproximation) then
       ! In the instantaneous stellar evolution approximation use the effective age to compute the SneIa yield.
       sneIaLifetime=self_%ageEffective
    else
       ! In the standard calculation simply use the current age.
       sneIaLifetime=lifetime_
    end if
    select case (indexElement_)
    case (0)
       ! Total metallicity required.
       yieldMass=supernovaeTypeIa_%yield(initialMassFunction_,massInitial,sneIaLifetime,metallicity_              )
    case default
       yieldMass=supernovaeTypeIa_%yield(initialMassFunction_,massInitial,sneIaLifetime,metallicity_,indexElement_)
    end select
    standardIntegrandYield=+standardIntegrandYield                &
         &                 +initialMassFunction_%phi(massInitial) &
         &                 *yieldMass
    return
  end function standardIntegrandYield

  double precision function standardIntegrandEnergyOutput(massInitial)
    !!{
    Integrand used in evaluating energy output.
    !!}
    implicit none
    double precision, intent(in   ) :: massInitial
    double precision                :: lifetime

    if (self_%instantaneousEnergyInputApproximation) then
       ! In the instantaneous stellar evolution approximation, assume stars more massive than the long-lived star cut off
       ! contribute to the energy input (with an age equal to the specified effective age), while less massive stars contribute
       ! nothing.
       if (massInitial > self_%massLongLived) then
          lifetime=self_%ageEffective
       else
          lifetime=0.0d0
       end if
    else
       ! In the standard calculation, simply use the current lifetime.
       lifetime=lifetime_
    end if
    standardIntegrandEnergyOutput=+initialMassFunction_%phi                  (                     massInitial                      ) &
         &                        *stellarFeedback_    %energyInputCumulative(initialMassFunction_,massInitial,lifetime,metallicity_)
    return
  end function standardIntegrandEnergyOutput

  logical function standardStarIsEvolved(self,massInitial,metallicity,age)
    !!{
    Returns true if the specified star is evolved by the given {\normalfont \ttfamily age}.
    !!}
    implicit none
    class           (stellarPopulationStandard), intent(inout) :: self
    double precision                           , intent(in   ) :: age        , massInitial, &
         &                                                        metallicity

    if (instantaneousApproximation) then
       ! Instantaneous calculation - star is evolved if it is more massive that the specified mass of long-lived stars.
       standardStarIsEvolved=                                      massInitial              > self%massLongLived
    else
       ! Standard calculation - star is evolved if its lifetime is less than the supplied age.
       standardStarIsEvolved=stellarAstrophysics_%lifetime(massInitial,metallicity) < age
    end if
    return
  end function standardStarIsEvolved

  double precision function standardRecycledFractionInstantaneous(self)
    !!{
    Return the recycled fraction from the stellar population in the instantaneous approximation.
    !!}
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    class(stellarPopulationStandard), intent(inout) :: self
    type (abundances               )                :: abundances_

    if (.not.self%recycledFractionInitialized) then
       call abundances_%metallicitySet(metallicitySolar)
       self%recycledFraction_          =self%rateRecycling(abundances_,ageMinimum=0.0d0,ageMaximum=self%ageEffective)*self%ageEffective
       self%recycledFractionInitialized=.true.
    end if
    standardRecycledFractionInstantaneous=self%recycledFraction_
    return
  end function standardRecycledFractionInstantaneous

  double precision function standardYieldInstantaneous(self)
    !!{
    Return the metal yield from the stellar population in the instantaneous approximation.
    !!}
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    class(stellarPopulationStandard), intent(inout) :: self
    type (abundances               )                :: abundances_

    if (.not.self%metalYieldInitialized) then
       call abundances_%metallicitySet(metallicitySolar)
       self%metalYield_          =self%rateYield(abundances_,ageMinimum=0.0d0,ageMaximum=self%ageEffective)*self%ageEffective
       self%metalYieldInitialized=.true.
    end if
    standardYieldInstantaneous=self%metalYield_
    return
  end function standardYieldInstantaneous

  !![
  <workaround type="gfortran" PR="93422" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=93422">
   <description>
    If the function name is used as the result variable, instead of using "result(spectra)", this PR is triggered.
   </description>
  </workaround>
  !!]
  function standardSpectra(self) result(spectra)
    !!{
    Return the stellar spectra associated with this population.
    !!}
    implicit none
    class(stellarPopulationSpectraClass), pointer       :: spectra
    class(stellarPopulationStandard    ), intent(inout) :: self

    spectra => self%stellarPopulationSpectra_
    return
  end function standardSpectra
