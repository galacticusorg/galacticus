!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  !+    Contributions to this file made by: Sachi Weerasooriya, Andrew Benson

  !!{
  Implements an emission line luminosity for AGN node property extractor class.
  !!}

  use    :: Numerical_Interpolation             , only : interpolator
  use    :: ISO_Varying_String                  , only : varying_string
  !$ use :: OMP_Lib                             , only : omp_lock_kind
  use    :: Output_Times                        , only : outputTimesClass
  use    :: Black_Hole_Accretion_Rates          , only : blackHoleAccretionRateClass
  use    :: Accretion_Disks                     , only : accretionDisksClass
  use    :: Atomic_Rates_Recombination_Radiative, only : atomicRecombinationRateRadiativeClass
  !![
  <nodePropertyExtractor name="nodePropertyExtractorLmnstyEmssnLineAGN">
    <description>
      An emission line luminosity property extractor class for AGN narrow line regions. The luminosity of the named emission lines
      (given by the {\normalfont \ttfamily lineNames} parameter are computed, largely following the model of
      \cite{feltre_nuclear_2016}.
    </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLmnstyEmssnLineAGN
     !!{
     A stellar luminosity output analysis property extractor class.
     !!}
     private
     class           (accretionDisksClass                  ), pointer                           :: accretionDisks_                   => null()
     class           (blackHoleAccretionRateClass          ), pointer                           :: blackHoleAccretionRate_           => null()
     class           (outputTimesClass                     ), pointer                           :: outputTimes_                      => null()
     class           (atomicRecombinationRateRadiativeClass), pointer                           :: atomicRecombinationRateRadiative_ => null()
     type            (varying_string                       ), allocatable, dimension(:        ) :: lineNames                                  , names_                  , &
          &                                                                                        descriptions_
     integer                                                                                    :: countLines
     integer         (c_size_t                             )                                    :: indexDensityHydrogen                       , indexIonizationParameter, &
          &                                                                                        indexMetallicity                           , indexSpectralIndex
     double precision                                       , allocatable, dimension(:        ) :: metallicity                                , densityHydrogen         , &
          &                                                                                        spectralIndex                              , ionizationParameter     , &
          &                                                                                        wavelengths
     double precision                                       , allocatable, dimension(:,:,:,:,:) :: luminosity
     type            (interpolator                         ), allocatable, dimension(:        ) :: interpolator_
     double precision                                                                           :: indexSpectralShortWavelength               , factorFillingVolume, densityHydrogenFixed, temperature
     !$ integer      (omp_lock_kind                        )                                    :: interpolateLock
   contains
     final     ::                 lmnstyEmssnLineAGNDestructor
     procedure :: elementCount => lmnstyEmssnLineAGNElementCount
     procedure :: extract      => lmnstyEmssnLineAGNExtract
     procedure :: names        => lmnstyEmssnLineAGNNames
     procedure :: descriptions => lmnstyEmssnLineAGNDescriptions
     procedure :: unitsInSI    => lmnstyEmssnLineAGNUnitsInSI
  end type nodePropertyExtractorLmnstyEmssnLineAGN

  interface nodePropertyExtractorLmnstyEmssnLineAGN
     !!{
     Constructors for the ``lmnstyEmssnLineAGN'' output analysis class.
     !!}
     module procedure lmnstyEmssnLineAGNConstructorParameters
     module procedure lmnstyEmssnLineAGNConstructorInternal
  end interface nodePropertyExtractorLmnstyEmssnLineAGN

  ! Enumeration for interpolants in the AGN emission line table.
  !![
  <enumeration>
   <name>interpolants</name>
   <description>Specifies the different interpolants for AGN emission line calculations.</description>
   <indexing>1</indexing>
   <entry label="density"            />
   <entry label="ionizationParameter"/>
   <entry label="metallicity"        />
   <entry label="spectralIndex"      />
  </enumeration>
  !!]
contains
  function lmnstyEmssnLineAGNConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``lmnstyEmssnLineAGN'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorLmnstyEmssnLineAGN)                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    type            (varying_string                         ), allocatable  , dimension(:) :: lineNames
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    class           (accretionDisksClass                    ), pointer                     :: accretionDisks_
    class           (blackHoleAccretionRateClass            ), pointer                     :: blackHoleAccretionRate_
    class           (atomicRecombinationRateRadiativeClass  ), pointer                     :: atomicRecombinationRateRadiative_
    double precision                                                                       :: indexSpectralShortWavelength     , factorFillingVolume, densityHydrogenFixed, temperature

    allocate(lineNames(parameters%count('lineNames')))
    !![
    <inputParameter>
      <name>lineNames</name>
      <source>parameters</source>
      <description>The names of the emission lines to extract.</description>
    </inputParameter>
    <inputParameter>
      <name>indexSpectralShortWavelength</name>
      <defaultValue>-1.7d0</defaultValue>
      <source>parameters</source>
      <description>The index, $\alpha$, of the power-law spectrum at wavelengths shortward of 0.25$\mu$m: $S_\nu \propto \nu^\alpha$ \citep{feltre_nuclear_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>factorFillingVolume</name>
      <defaultValue>0.01d0</defaultValue>
      <source>parameters</source>
      <description>The volume-filling factor, i.e. the ratio of the volume-averaged hydrogen density to the hydrogen density.</description>
    </inputParameter>
    <inputParameter>
      <name>densityHydrogenFixed</name>
      <defaultValue>1.0d+3</defaultValue>
      <source>parameters</source>
      <description>Density of hydrogen in cm⁻³ </description>
    </inputParameter>
    <inputParameter>
      <name>temperature</name>
      <defaultValue>+1.00d+4</defaultValue>
      <source>parameters</source>
      <description> emperature in K </description>
    </inputParameter>
    <objectBuilder class="accretionDisks"                   name="accretionDisks_"                   source="parameters"/>
    <objectBuilder class="blackHoleAccretionRate"           name="blackHoleAccretionRate_"           source="parameters"/>
    <objectBuilder class="outputTimes"                      name="outputTimes_"                      source="parameters"/>
    <objectBuilder class="atomicRecombinationRateRadiative" name="atomicRecombinationRateRadiative_" source="parameters"/>
    !!]
    self=nodePropertyExtractorLmnstyEmssnLineAGN(accretionDisks_,blackHoleAccretionRate_,outputTimes_,atomicRecombinationRateRadiative_,lineNames,indexSpectralShortWavelength,factorFillingVolume,densityHydrogenFixed, temperature)
    !![
    <inputParametersValidate source="parameters"              />
    <objectDestructor name="accretionDisks_"                  />
    <objectDestructor name="blackHoleAccretionRate_"          />
    <objectDestructor name="outputTimes_"                     />
    <objectDestructor name="atomicRecombinationRateRadiative_"/>
    !!]
    return
  end function lmnstyEmssnLineAGNConstructorParameters

  function lmnstyEmssnLineAGNConstructorInternal(accretionDisks_,blackHoleAccretionRate_,outputTimes_,atomicRecombinationRateRadiative_,lineNames,indexSpectralShortWavelength,factorFillingVolume,densityHydrogenFixed,temperature,outputMask) result(self)
    !!{
    Internal constructor for the ``lmnstyEmssnLineAGN'' output analysis property extractor class.
    !!}
    use            :: Error                         , only : Error_Report
    use            :: Input_Paths                   , only : inputPath             , pathTypeDataStatic
    use            :: HDF5_Access                   , only : hdf5Access
    use            :: IO_HDF5                       , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Instruments_Filters           , only : Filter_Extent         , Filter_Get_Index
    use            :: Output_Times                  , only : outputTimesClass
    use            :: String_Handling               , only : String_Join           , char
    use            :: Galacticus_Nodes              , only : nodeComponentBlackHole
     use           :: Table_Labels                  , only : extrapolationTypeFix
    implicit none
    type            (nodePropertyExtractorLmnstyEmssnLineAGN)                                                :: self
    double precision                                         , intent(in   )                                 :: indexSpectralShortWavelength,     factorFillingVolume, densityHydrogenFixed, temperature
    type            (varying_string                         ), intent(in   ), dimension(:        )           :: lineNames
    logical                                                  , intent(in   ), dimension(:        ), optional :: outputMask
    class           (accretionDisksClass                    ), intent(in   ), target                         :: accretionDisks_
    class           (blackHoleAccretionRateClass            ), intent(in   ), target                         :: blackHoleAccretionRate_
    class           (outputTimesClass                       ), intent(in   ), target                         :: outputTimes_
    class           (atomicRecombinationRateRadiativeClass  ), intent(in   ), target                         :: atomicRecombinationRateRadiative_
    type            (hdf5Object                             )                                                :: emissionLinesFile                , lines             , &
         &                                                                                                      lineDataset                      , dataset
    integer                                                                                                  :: i
    integer         (c_size_t)                                              , dimension(5        )           :: shapeLines                       , permutation
    double precision                                         , allocatable  , dimension(:,:,:,:,:)           :: luminosity
    !![
    <constructorAssign variables="lineNames, indexSpectralShortWavelength, factorFillingVolume, densityHydrogenFixed, temperature, *accretionDisks_, *blackHoleAccretionRate_, *outputTimes_, *atomicRecombinationRateRadiative_"/>
    !!]
    
    ! Read the table of emission line luminosities.
    !$ call hdf5Access%set()
    call emissionLinesFile%openFile(char(inputPath(pathTypeDataStatic))//"hiiRegions/emissionLineLuminosities_AGN.hdf5",readOnly=.true.)
    lines=emissionLinesFile%openGroup('lines')
    do i=1,size(lineNames)
       if (.not.lines%hasDataset(char(self%lineNames(i)))) call Error_Report('line "'//char(self%lineNames(i))//'" not found'//{introspection:location})
    end do
    call emissionLinesFile%readDataset('densityHydrogen'    ,self%densityHydrogen    )
    call emissionLinesFile%readDataset('ionizationParameter',self%ionizationParameter)
    call emissionLinesFile%readDataset('metallicity'        ,self%metallicity        )
    call emissionLinesFile%readDataset('spectralIndex'      ,self%spectralIndex      )
    ! Extract indexing into the lines arrays.
    dataset=emissionLinesFile%openDataset('densityHydrogen'    )
    call dataset%readAttribute('index',self%indexDensityHydrogen    )
    call dataset%close        (                                     )
    dataset=emissionLinesFile%openDataset('ionizationParameter')
    call dataset%readAttribute('index',self%indexIonizationParameter)
    call dataset%close        (                                     )
    dataset=emissionLinesFile%openDataset('metallicity'        )
    call dataset%readAttribute('index',self%indexMetallicity        )
    call dataset%close        (                                     )
    dataset=emissionLinesFile%openDataset('spectralIndex'      )
    call dataset%readAttribute('index',self%indexSpectralIndex      )
    call dataset%close        (                                     )
    ! Offset indexing to Fortran standard (i.e. starting from 1 instead of 0).
    self%indexDensityHydrogen    =self%indexDensityHydrogen    +1
    self%indexIonizationParameter=self%indexIonizationParameter+1
    self%indexMetallicity        =self%indexMetallicity        +1
    self%indexSpectralIndex      =self%indexSpectralIndex      +1
    ! Establish arrays.
    shapeLines(self%indexDensityHydrogen    )=size(self%densityHydrogen    )
    shapeLines(self%indexIonizationParameter)=size(self%ionizationParameter)
    shapeLines(self%indexMetallicity        )=size(self%metallicity        )
    shapeLines(self%indexSpectralIndex      )=size(self%spectralIndex      )
    shapeLines(     5                       )=size(     lineNames          )
    ! Allocate a temporary luminosities array into which we will read data from the table. The dimension ordering here is whatever
    ! ordering was used in the table file. This will be reordered into our preferred, internal order later.
    !![
    <allocate variable="luminosity" shape="shapeLines"/>
    !!]
    allocate(self%wavelengths(size(self%lineNames)))
    do i=1,size(self%lineNames)
       call lines%readDatasetStatic(char(self%lineNames(i)),luminosity(:,:,:,:,i))
       lineDataset=lines%openDataset(char(self%lineNames(i)))
       call lineDataset%readAttribute('wavelength',self%wavelengths(i))
       call lineDataset%close        (                                )
    end do
    call lines            %close()
    call emissionLinesFile%close()
    !$ call hdf5Access%unset()
    ! Re-order the luminosities table into our preferred order.
    !! First, allocate our final table array with our preferred ordering of dimensions.
    allocate(                                 &
         &   self%luminosity                  &
         &   (                                &
         &    size(self%spectralIndex      ), &
         &    size(self%metallicity        ), &
         &    size(self%ionizationParameter), &
         &    size(self%densityHydrogen    ), &
         &    size(self%lineNames          )  &
         &   )                                &
         &  )
    !! Construct a permutation - mapping the indices of dimensions in the file into the order we want internally.
    permutation=[                               &
         &       self%indexSpectralIndex      , &
         &       self%indexMetallicity        , &
         &       self%indexIonizationParameter, &
         &       self%indexDensityHydrogen    , &
         &       5_c_size_t                     &
         &      ]
    !! Reorder the table read from file into our internal table.
    self%luminosity=reshape(luminosity,shape(self%luminosity),order=permutation)
    ! Convert parameters and luminosities to logarithmic form.
    self%densityHydrogen    =log10(self%densityHydrogen    )
    self%ionizationParameter=log10(self%ionizationParameter)
    self%metallicity        =log10(self%metallicity        )
    ! Initialize interpolators.
    allocate(self%interpolator_(4))
    self%interpolator_(interpolantsDensity            %ID)=interpolator(self%densityHydrogen,extrapolationType=extrapolationTypeFix)
    self%interpolator_(interpolantsIonizationParameter%ID)=interpolator(self%ionizationParameter,extrapolationType=extrapolationTypeFix)
    self%interpolator_(interpolantsMetallicity        %ID)=interpolator(self%metallicity,extrapolationType=extrapolationTypeFix)
    self%interpolator_(interpolantsSpectralIndex      %ID)=interpolator(self%spectralIndex,extrapolationType=extrapolationTypeFix)
    !$ call OMP_Init_Lock(self%interpolateLock)
    ! Construct names and descriptions.
    allocate(self%names_       (size(lineNames)))
    allocate(self%descriptions_(size(lineNames)))
    do i=1,size(lineNames)
       self%names_       (i)="luminosityEmissionLineAGN:"//lineNames(i)
       self%descriptions_(i)="Luminosity of the "        //lineNames(i)//" AGN emission line [ergs/s]"
    end do
    self%countLines=size(lineNames)
    return
  end function lmnstyEmssnLineAGNConstructorInternal

  subroutine lmnstyEmssnLineAGNDestructor(self)
    !!{
    Destructor for the ``lmnstyEmssnLineAGN'' output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLmnstyEmssnLineAGN), intent(inout) :: self

    !$ call OMP_Destroy_Lock(self%interpolateLock)
    !![
    <objectDestructor name="self%accretionDisks_"                  />
    <objectDestructor name="self%blackHoleAccretionRate_"          />
    <objectDestructor name="self%outputTimes_"                     />
    <objectDestructor name="self%atomicRecombinationRateRadiative_"/>
    !!]
    return
  end subroutine lmnstyEmssnLineAGNDestructor

  function lmnstyEmssnLineAGNExtract(self,node,time,instance)
    !!{
    Implement an emission line output analysis property extractor.
    !!}
    use            :: Atomic_Rates_Recombination_Radiative, only : recombinationCaseB
    use            :: Abundances_Structure                , only : abundances           , max                  , metallicityTypeLogarithmicByMassSolar
    use            :: Galacticus_Nodes                    , only : nodeComponentBasic   , nodeComponentDisk    , nodeComponentSpheroid                    , nodeComponentBlackHole, &
         &                                                         treeNode
    use, intrinsic :: ISO_C_Binding                       , only : c_size_t
    use            :: Numerical_Constants_Astronomical    , only : hydrogenByMassSolar  , luminosityZeroPointAB, massSolar                                , megaParsec            , &
         &                                                         metallicitySolar     , parsec               , massSolar                                , gigaYear
    use            :: Numerical_Constants_Atomic          , only : atomicMassHydrogen   , atomicMassUnit       , lymanSeriesLimitWavelengthHydrogen_atomic
    use            :: Numerical_Constants_Math            , only : Pi
    use            :: Numerical_Constants_Physical        , only : plancksConstant      , speedLight
    use            :: Numerical_Constants_Prefixes        , only : centi                , hecto                , mega                                     , micro
    use            :: Numerical_Constants_Units           , only : metersToAngstroms
    use            :: Galactic_Structure_Options          , only : massTypeStellar
    implicit none
    double precision                                                 , dimension(:) , allocatable :: lmnstyEmssnLineAGNExtract
    class           (nodePropertyExtractorLmnstyEmssnLineAGN        ), intent(inout), target      :: self
    type            (treeNode                                       ), intent(inout), target      :: node
    double precision                                                 , intent(in   )              :: time
    type            (multiCounter                                   ), intent(inout), optional    :: instance
    class           (nodeComponentBasic                             ), pointer                    :: basic
    class           (nodeComponentDisk                              ), pointer                    :: disk
    class           (nodeComponentSpheroid                          ), pointer                    :: spheroid
    class           (nodeComponentBlackHole                         ), pointer                    :: blackHole
    double precision                                                 , parameter                  :: massMinimum              =+1.00d-06
    double precision                                                 , parameter                  :: radiusMinimum            =+1.00d-06
    double precision                                                 , parameter                  :: metallicityISMLocal      =+2.00d-02 ! Metallicity in the local ISM.
    double precision                                                 , parameter                  :: wavelengthZeroPoint      =+5.50d+03 ! Angstroms.
    double precision                                                 , parameter                  :: frequency0p001Microns    =speedLight/( 0.001d0*micro) ! Frequency at  0.001μm in Hz.
    double precision                                                 , parameter                  :: frequency0p250Microns    =speedLight/( 0.250d0*micro) ! Frequency at  0.250μm in Hz.
    double precision                                                 , parameter                  :: frequency10p00Microns    =speedLight/(10.000d0*micro) ! Frequency at 10.000μm in Hz.
    double precision                                                 , parameter                  :: frequencyLymanLimit      =speedLight/(lymanSeriesLimitWavelengthHydrogen_atomic/metersToAngstroms) ! Frequency at the Lyman limit in Hz.
    type            (abundances                                     ),                            :: abundancesGas
    logical                                                          ,                            :: isPhysical
    integer         (c_size_t                                       ), dimension(0:1,4)           :: interpolateIndex
    double precision                                                 , dimension(0:1,4)           :: interpolateFactor
    double precision                                                                              :: weight                                                 
    double precision                                                                              :: luminosityBolometricUnnormalized, recombinationCoefficient, &
         &                                                                                           rateMassAccretionSpheroid       , rateMassAccretionHotHalo, &
         &                                                                                           rateAccretionNuclearStarCluster , luminosityBolometricAGN , &
         &                                                                                           normalization                                             , &
         &                                                                                           radiativeEfficiency             , rateAccretionBlackHole  , &
         &                                                                                           densityHydrogenLogarithmic      , rateIonizingPhotons     , &
         &                                                                                           ionizationParameterLogarithmic  , massGas                 , &
         &                                                                                           radiusGalaxy                    , metallicityGas          , &
         &                                                                                           radiusStromgren                 , massHydrogen
    integer         (c_size_t                                       )                             :: output
    integer                                                                                       :: i                               , j                       , &
         &                                                                                           k                               , l                       , &
         &                                                                                           line
    !$GLC attributes unused :: instance

    ! Retrieve components.
    basic     => node%basic    ()
    disk      => node%disk     ()
    spheroid  => node%spheroid ()
    blackHole => node%blackHole()
    ! Determine output index.
    output   =  self%outputTimes_%index(basic%time(),findClosest=.true.)
    ! Extract all required properties.
    abundancesGas  =(disk%abundancesGas()+spheroid%abundancesGas())
    massGas        =(disk%massGas      ()+spheroid%massGas      ())
    radiusGalaxy   =(disk%radius       ()+spheroid%radius       ())*megaParsec ! SI units.
    ! Determine if component is physically reasonable.
    isPhysical= massGas      > massMinimum   &
         &     .and.                         &
         &      radiusGalaxy > radiusMinimum             
    ! Compute the logarithmic metallicity of the gas in each component in Solar units.
    call abundancesGas%massToMassFraction(massGas)
    if (isPhysical) then
       metallicityGas=abundancesGas%metallicity(metallicityTypeLogarithmicByMassSolar)
    else
       metallicityGas=0.0d0
    end if
    ! Find the hydrogen mass in SI units.
    massHydrogen=+              massGas                &
         &       *              massSolar              &
         &       *abundancesGas%hydrogenMassFraction()
    ! Compute the integral over all frequencies of the unnormalized AGN spectrum. The spectrum is taken from Feltre, Gutkin &
    ! Charlot (2016; MNRAS; 456; 3354; https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.3354F), their equation 5.
    luminosityBolometricUnnormalized=+  frequency0p250Microns**(self%indexSpectralShortWavelength+0.50d0) & ! ⎫ Match normalization of piecewise-power-law at 10.0μm.
         &                           *  frequency10p00Microns**                                  -2.50d0  & ! ⎭ 
         &                           *  frequency10p00Microns**                                  +3.00d0  & ! ⎫ Integral of ν²...
         &                           /                                                            3.00d0  & ! ⎭ ...from 0μm to 10.0μm.
         &                           +  frequency0p250Microns**(self%indexSpectralShortWavelength+0.25d0) & ! } Match normalization of piecewise-power-law at 0.25μm.
         &                           *(                                                                   & ! ⎫ Integral of ν^{-0.5}...
         &                             +frequency0p250Microns**                                   +0.5d0  & ! ⎪ ...from 0.250μm...
         &                             -frequency10p00Microns**                                   +0.5d0  & ! ⎬ ...to 10.0μm.
         &                            )                                                                   & ! ⎪
         &                           /                                                             0.5d0  & ! ⎭
         &                           +(                                                                   & ! ⎫ Integral of ν^α...
         &                             +frequency0p001Microns**(self%indexSpectralShortWavelength+1.00d0) & ! ⎪ ...from 0.001μm...
         &                             -frequency0p250Microns**(self%indexSpectralShortWavelength+1.00d0) & ! ⎬ ...to 0.250μm.
         &                            )                                                                   & ! ⎪
         &                           /                         (self%indexSpectralShortWavelength+1.00d0)   ! ⎭
    ! Find the corresponding normalization constant to unit bolometric luminosity
    normalization=+1.0d0                            &
         &        /luminosityBolometricUnnormalized
    ! Get the hydrogen recombination coefficient (in m³ s⁻¹).
    recombinationCoefficient=+self%atomicRecombinationRateRadiative_%rate(1,1,self%temperature,recombinationCaseB) &
         &                   *micro

    ! Get black hole accretion rates.
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo,rateAccretionNuclearStarCluster)
    ! Find the total black hole accretion rate in kg s⁻¹.
    rateAccretionBlackHole=+massSolar                         &
         &                 /gigaYear                          &
         &                 *(                                 &
         &                   +rateMassAccretionHotHalo        &
         &                   +rateMassAccretionSpheroid       &
         &                   +rateAccretionNuclearStarCluster &
         &                 )
    ! Get the radiative efficiency of black hole.
    radiativeEfficiency=self%accretionDisks_%efficiencyRadiative(blackHole,rateMassAccretionSpheroid+rateMassAccretionHotHalo+rateAccretionNuclearStarCluster)
    ! Compute the bolometric luminosity of AGN in W.
    luminosityBolometricAGN=+rateAccretionBlackHole    &
         &                  *speedLight            **2 &
         &                  *radiativeEfficiency
    ! Find the emission rate of ionizing photons, Q_H, in s⁻¹. This is done by integrating over the AGN spectrum: 
    ! ∫_vLy^ν(0.001μm) S_ν/hν dν, where νLy is the frequency at the Lyman limit.
    rateIonizingPhotons=+normalization                                              &
         &              *luminosityBolometricAGN                                    &
         &              /plancksConstant                                            &
         &              *(                                                          &
         &                +frequency0p001Microns**self%indexSpectralShortWavelength &
         &                -frequencyLymanLimit  **self%indexSpectralShortWavelength &
         &               )                                                          &
         &              /                         self%indexSpectralShortWavelength
    ! Calculate the Strömgren radius (equation 3 of Feltre, Gutkin & Charlot (2016; MNRAS; 456; 3354;
    ! https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.3354F).
    radiusStromgren=+(                             &
         &            +3.0d0                       &
         &            *rateIonizingPhotons         &
         &            /(                           &
         &              +4.0d0                     &
         &              *Pi                        &
         &              *(self%densityHydrogenFixed*mega)**2 &
         &              *self%factorFillingVolume  &
         &              *recombinationCoefficient  &
         &             )                           &
         &           )**(1.0d0/3.0d0)
    ! Compute the ionization parameter.
    if (luminosityBolometricAGN > 0.0d0) then
       ionizationParameterLogarithmic=log10(                      &
            &                               +rateIonizingPhotons  &
            &                               /4.0d0                &
            &                               /Pi                   &
            &                               /radiusStromgren**2   &
            &                               /self%densityHydrogenFixed &
            &                               /mega                 &
            &                               /speedLight           &
            &                              )
    else
      ionizationParameterLogarithmic=0.0d0
    end if
    ! Find the logarithmic density.
    densityHydrogenLogarithmic=log10(self%densityHydrogenFixed)
    ! Truncate properties to table bounds where necessary to avoid unphysical extrapolations.
    if      (metallicityGas                  < self%metallicity       (1                             )) then
      metallicityGas                 =self%metallicity        (1                             )
    else if (metallicityGas                  > self%metallicity       (size(self%metallicity        ))) then
       metallicityGas                =self%metallicity        (size(self%metallicity        ))
    end if
    if      (ionizationParameterLogarithmic < self%ionizationParameter(1                             )) then
       ionizationParameterLogarithmic=self%ionizationParameter(1                             )
    else if (ionizationParameterLogarithmic > self%ionizationParameter(size(self%ionizationParameter))) then
       ionizationParameterLogarithmic=self%ionizationParameter(size(self%ionizationParameter))
    end if
    ! Find interpolating factors in all four interpolants, preventing extrapolation beyond the tabulated ranges.
    !$ call OMP_Set_Lock  (self%interpolateLock)
    call self%interpolator_(interpolantsDensity            %ID)%linearFactors(densityHydrogenLogarithmic       ,interpolateIndex(0,interpolantsDensity            %ID),interpolateFactor(:,interpolantsDensity            %ID))
    call self%interpolator_(interpolantsIonizationParameter%ID)%linearFactors(ionizationParameterLogarithmic   ,interpolateIndex(0,interpolantsIonizationParameter%ID),interpolateFactor(:,interpolantsIonizationParameter%ID))
    call self%interpolator_(interpolantsMetallicity        %ID)%linearFactors(metallicityGas                   ,interpolateIndex(0,interpolantsMetallicity        %ID),interpolateFactor(:,interpolantsMetallicity        %ID))
    call self%interpolator_(interpolantsSpectralIndex      %ID)%linearFactors(self%indexSpectralShortWavelength,interpolateIndex(0,interpolantsSpectralIndex      %ID),interpolateFactor(:,interpolantsSpectralIndex      %ID))
    !$ call OMP_Unset_Lock(self%interpolateLock)
    interpolateIndex(1,:)=interpolateIndex(0,:)+1
    interpolateFactor=max(min(interpolateFactor,1.0d0),0.0d0)
    ! Iterate over lines.
    allocate(lmnstyEmssnLineAGNExtract(self%countLines))
    lmnstyEmssnLineAGNExtract=0.0d0
    do line=1,size(self%luminosity,dim=5)
       ! Interpolate in all four interpolants.
       do i=0,1
          do j=0,1
             do k=0,1
                do l=0,1
                   weight                         =+                interpolateFactor(i   ,1)  &
                        &                          *                interpolateFactor(j   ,2)  &
                        &                          *                interpolateFactor(k   ,3)  &
                        &                          *                interpolateFactor(l   ,4)  
                   lmnstyEmssnLineAGNExtract(line)=+lmnstyEmssnLineAGNExtract        (line  )  &
                        &                          +weight                                     &
                        &                          *(radiusStromgren*hecto)**2                 & ! The quantity given by Cloudy is (4π*intensity in cm⁻²). Multiply by the Strömgren radius squared to convert to luminosity.
                        &                          *self%luminosity(                           &
                        &                                           interpolateIndex (l   ,4), &
                        &                                           interpolateIndex (k   ,3), &
                        &                                           interpolateIndex (j   ,2), &
                        &                                           interpolateIndex (i   ,1), &
                        &                                           line                       &
                        &                                          )
                end do
             end do
          end do
       end do
    end do
    return
  end function lmnstyEmssnLineAGNExtract

  integer function lmnstyEmssnLineAGNElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily lmnstyEmssnLineAGN} property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorLmnstyEmssnLineAGN), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: time

    lmnstyEmssnLineAGNElementCount=self%countLines
    return
  end function lmnstyEmssnLineAGNElementCount

  subroutine lmnstyEmssnLineAGNNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily emissionLines}.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class           (nodePropertyExtractorLmnstyEmssnLineAGN), intent(inout)                            :: self
    double precision                                         , intent(in   )                            :: time
    type            (varying_string                         ), intent(inout), dimension(:), allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%countLines))
    names=self%names_
    return
  end subroutine lmnstyEmssnLineAGNNames

  subroutine lmnstyEmssnLineAGNDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily emission line luminosity} property.
    !!}
    implicit none
    class           (nodePropertyExtractorLmnstyEmssnLineAGN), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(self%countLines))
    descriptions=self%descriptions_
    return
  end subroutine lmnstyEmssnLineAGNDescriptions

  function lmnstyEmssnLineAGNUnitsInSI(self,time) result(unitsInSI)
  !!{
    Return the units of the {\normalfont \ttfamily lmnstyEmssnLineAGN} properties in the SI system.
    !!}
    use :: Numerical_Constants_Units, only : ergs
    implicit none
    double precision                                         , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorLmnstyEmssnLineAGN), intent(inout)               :: self
    double precision                                         , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(unitsInSI(self%countLines))
    unitsInSI=ergs
    return
  end function lmnstyEmssnLineAGNUnitsInSI
 
