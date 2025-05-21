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

! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
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

!!{
Implements an emission line luminosity node property extractor class.
!!}

  use    :: Numerical_Interpolation          , only : interpolator
  use    :: ISO_Varying_String               , only : varying_string
  !$ use :: OMP_Lib                          , only : omp_lock_kind
  use    :: Output_Times                     , only : outputTimesClass
  use    :: Star_Formation_Rates_Disks       , only : starFormationRateDisksClass
  use    :: Star_Formation_Rates_Spheroids   , only : starFormationRateSpheroidsClass
  use    :: Stellar_Spectra_Dust_Attenuations, only : stellarSpectraDustAttenuationClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorLmnstyEmssnLinePanuzzo2003">
    <description>
      An emission line luminosity property extractor class. The luminosity of the named emission line (given by the {\normalfont
      \ttfamily lineNames} parameter: if multiple lines are named, the sum of their luminosities) is computed. Additional dust
      attenuation for emission line luminosities can be specified via the {\normalfont \ttfamily depthOpticalISMCoefficient}
      parameter.
    </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorLmnstyEmssnLinePanuzzo2003
     !!{
     A stellar luminosity output analysis property extractor class.
     !!}
     private
     class           (starFormationRateDisksClass       ), pointer                             :: starFormationRateDisks_        => null()
     class           (starFormationRateSpheroidsClass   ), pointer                             :: starFormationRateSpheroids_    => null()
     class           (stellarSpectraDustAttenuationClass), pointer                             :: stellarSpectraDustAttenuation_ => null()
     class           (outputTimesClass                  ), pointer                             :: outputTimes_                   => null()
     type            (varying_string                    )                                      :: name_                                   , description_
     type            (varying_string                    ), allocatable, dimension(:          ) :: lineNames
     double precision                                    , allocatable, dimension(:          ) :: metallicity                             , densityHydrogen             , &
          &                                                                                       ionizingFluxHydrogen                    , ionizingFluxHeliumToHydrogen, &
          &                                                                                       ionizingFluxOxygenToHelium              , wavelength
     double precision                                    , allocatable, dimension(:,:,:,:,:,:) :: luminosity
     integer                                             , allocatable, dimension(:,:        ) :: ionizingContinuumIndex
     double precision                                                 , dimension(2,3        ) :: filterExtent
     type            (interpolator                      ), allocatable, dimension(:          ) :: interpolator_
     double precision                                                                          :: depthOpticalISMCoefficient
     !$ integer      (omp_lock_kind                     )                                      :: interpolateLock
   contains
     final     ::                lmnstyEmssnLinePanuzzo2003Destructor
     procedure :: extract     => lmnstyEmssnLinePanuzzo2003Extract
     procedure :: quantity    => lmnstyEmssnLinePanuzzo2003Quantity
     procedure :: name        => lmnstyEmssnLinePanuzzo2003Name
     procedure :: description => lmnstyEmssnLinePanuzzo2003Description
     procedure :: unitsInSI   => lmnstyEmssnLinePanuzzo2003UnitsInSI
  end type nodePropertyExtractorLmnstyEmssnLinePanuzzo2003

  interface nodePropertyExtractorLmnstyEmssnLinePanuzzo2003
     !!{
     Constructors for the \refClass{nodePropertyExtractorLmnstyEmssnLinePanuzzo2003} output analysis class.
     !!}
     module procedure lmnstyEmssnLinePanuzzo2003ConstructorParameters
     module procedure lmnstyEmssnLinePanuzzo2003ConstructorInternal
  end interface nodePropertyExtractorLmnstyEmssnLinePanuzzo2003

  ! Enumerations for galactic components and ionizing continuua.
  !![
  <enumeration>
   <name>component</name>
   <description>Specifies the galactic component for emission line calculations.</description>
   <indexing>1</indexing>
   <entry label="disk"    />
   <entry label="spheroid"/>
  </enumeration>
  <enumeration>
   <name>ionizingContinuum</name>
   <description>Specifies the ionizing continuum for emission line calculations.</description>
   <indexing>1</indexing>
   <entry label="Hydrogen"/>
   <entry label="Helium"  />
   <entry label="Oxygen"  />
  </enumeration>
  <enumeration>
   <name>interpolant</name>
   <description>Specifies the different interpolants for emission line calculations.</description>
   <indexing>1</indexing>
   <entry label="metallicity"/>
   <entry label="density"    />
   <entry label="hydrogen"   />
   <entry label="helium"     />
   <entry label="oxygen"     />
  </enumeration>
  !!]

contains

  function lmnstyEmssnLinePanuzzo2003ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorLmnstyEmssnLinePanuzzo2003} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorLmnstyEmssnLinePanuzzo2003)                              :: self
    type            (inputParameters                                ), intent(inout)               :: parameters
    type            (varying_string                                 ), allocatable  , dimension(:) :: lineNames
    class           (starFormationRateDisksClass                    ), pointer                     :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                ), pointer                     :: starFormationRateSpheroids_
    class           (stellarSpectraDustAttenuationClass             ), pointer                     :: stellarSpectraDustAttenuation_
    class           (outputTimesClass                               ), pointer                     :: outputTimes_
    double precision                                                                               :: depthOpticalISMCoefficient

    allocate(lineNames(parameters%count('lineNames')))
    !![
    <inputParameter>
      <name>lineNames</name>
      <source>parameters</source>
      <description>The emission lines to extract.</description>
    </inputParameter>
    <inputParameter>
      <name>depthOpticalISMCoefficient</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>Multiplicative coefficient for optical depth in the ISM.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateDisks"        name="starFormationRateDisks_"        source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"    name="starFormationRateSpheroids_"    source="parameters"/>
    <objectBuilder class="stellarSpectraDustAttenuation" name="stellarSpectraDustAttenuation_" source="parameters"/>
    <objectBuilder class="outputTimes"                   name="outputTimes_"                   source="parameters"/>
    !!]
    self=nodePropertyExtractorLmnstyEmssnLinePanuzzo2003(starFormationRateDisks_,starFormationRateSpheroids_,stellarSpectraDustAttenuation_,outputTimes_,lineNames,depthOpticalISMCoefficient)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"       />
    <objectDestructor name="starFormationRateSpheroids_"   />
    <objectDestructor name="stellarSpectraDustAttenuation_"/>
    <objectDestructor name="outputTimes_"                  />
    !!]
    return
  end function lmnstyEmssnLinePanuzzo2003ConstructorParameters

  function lmnstyEmssnLinePanuzzo2003ConstructorInternal(starFormationRateDisks_,starFormationRateSpheroids_,stellarSpectraDustAttenuation_,outputTimes_,lineNames,depthOpticalISMCoefficient,outputMask) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorLmnstyEmssnLinePanuzzo2003} output analysis property extractor class.
    !!}
    use            :: Error                         , only : Error_Report
    use            :: Input_Paths                   , only : inputPath              , pathTypeDataStatic
    use            :: HDF5_Access                   , only : hdf5Access
    use            :: IO_HDF5                       , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Instruments_Filters           , only : Filter_Extent          , Filter_Get_Index
    use            :: Output_Times                  , only : outputTimesClass
    use            :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    use            :: String_Handling               , only : String_Join            , char
    implicit none
    type            (nodePropertyExtractorLmnstyEmssnLinePanuzzo2003)                                        :: self
    double precision                                                 , intent(in   )                         :: depthOpticalISMCoefficient
    type            (varying_string                                 ), intent(in   ), dimension(:)           :: lineNames
    logical                                                          , intent(in   ), dimension(:), optional :: outputMask
    class           (starFormationRateDisksClass                    ), intent(in   ), target                 :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass                ), intent(in   ), target                 :: starFormationRateSpheroids_
    class           (stellarSpectraDustAttenuationClass             ), intent(in   ), target                 :: stellarSpectraDustAttenuation_
    class           (outputTimesClass                               ), intent(in   ), target                 :: outputTimes_
    type            (hdf5Object                                     )                                        :: emissionLinesFile             , lines, &
         &                                                                                                      lineDataset
    integer         (c_size_t                                       )                                        :: i
    !![
    <constructorAssign variables="lineNames, depthOpticalISMCoefficient, *starFormationRateDisks_, *starFormationRateSpheroids_, *stellarSpectraDustAttenuation_, *outputTimes_"/>
    !!]

    ! Read the table of emission line luminosities.
    !$ call hdf5Access%set()
    call emissionLinesFile%openFile(char(inputPath(pathTypeDataStatic))//"hiiRegions/emissionLinesPanuzzo2003.hdf5",readOnly=.true.)
    lines=emissionLinesFile%openGroup('lines')
    do i=1,size(lineNames)
       if (.not.lines%hasDataset(char(self%lineNames(i)))) call Error_Report('line "'//char(self%lineNames(i))//'" not found'//{introspection:location})
    end do
    call emissionLinesFile%readDataset('metallicity'                  ,self%metallicity                 )
    call emissionLinesFile%readDataset('densityHydrogen'              ,self%densityHydrogen             )
    call emissionLinesFile%readDataset('ionizingFluxHydrogen'         ,self%ionizingFluxHydrogen        )
    call emissionLinesFile%readDataset('ionizingFluxHeliumToHydrogen' ,self%ionizingFluxHeliumToHydrogen)
    call emissionLinesFile%readDataset('ionizingFluxOxygenToHelium'   ,self%ionizingFluxOxygenToHelium  )
    allocate(                                          &
         &   self%luminosity                           &
         &   (                                         &
         &    size(self%ionizingFluxOxygenToHelium  ), &
         &    size(self%ionizingFluxHeliumToHydrogen), &
         &    size(self%ionizingFluxHydrogen        ), &
         &    size(self%densityHydrogen             ), &
         &    size(self%metallicity                 ), &
         &    size(self%lineNames                   )  &
         &   )                                         &
         &  )
    allocate(                                          &
         &   self%wavelength                           &
         &   (                                         &
         &    size(self%lineNames                   )  &
         &   )                                         &
         &  )
    do i=1,size(lineNames)
       call lines      %readDatasetStatic(char(self%lineNames(i)),self%luminosity(:,:,:,:,:,i))
       lineDataset=lines%openDataset(char(self%lineNames(i)))
       call lineDataset%readAttribute('wavelength',self%wavelength(i))
       call lineDataset%close        (                               )
    end do
    call lines            %close      (                                                                 )
    call emissionLinesFile%close      (                                                                 )
    !$ call hdf5Access%unset()
    ! Convert parameters and luminosities to log form.
    self%metallicity                 =log10(self%metallicity                 )
    self%densityHydrogen             =log10(self%densityHydrogen             )
    self%ionizingFluxHydrogen        =log10(self%ionizingFluxHydrogen        )
    self%ionizingFluxHeliumToHydrogen=log10(self%ionizingFluxHeliumToHydrogen)
    self%ionizingFluxOxygenToHelium  =log10(self%ionizingFluxOxygenToHelium  )
    self%luminosity                  =log10(self%luminosity                  )
    ! Find indices of ionizing continuua filters.
    allocate(self%ionizingContinuumIndex(self%outputTimes_%count(),3_c_size_t))
    do i=1,self%outputTimes_%count()
       if (present(outputMask).and..not.outputMask(i)) then
          self%ionizingContinuumIndex(i,:                        )=-1
       else
          self%ionizingContinuumIndex(i,ionizingContinuumHydrogen%ID)=unitStellarLuminosities%index('Lyc'            ,'rest',self%outputTimes_%redshift(i))
          self%ionizingContinuumIndex(i,ionizingContinuumHelium  %ID)=unitStellarLuminosities%index('HeliumContinuum','rest',self%outputTimes_%redshift(i))
          self%ionizingContinuumIndex(i,ionizingContinuumOxygen  %ID)=unitStellarLuminosities%index('OxygenContinuum','rest',self%outputTimes_%redshift(i))
       end if
    end do
    ! Read wavelength intervals of ionizing continuum filters.
    self%filterExtent(:,ionizingContinuumHydrogen%ID)=Filter_Extent(Filter_Get_Index(var_str('Lyc'            )))
    self%filterExtent(:,ionizingContinuumHelium  %ID)=Filter_Extent(Filter_Get_Index(var_str('HeliumContinuum')))
    self%filterExtent(:,ionizingContinuumOxygen  %ID)=Filter_Extent(Filter_Get_Index(var_str('OxygenContinuum')))
    ! Initialize interpolators.
    allocate(self%interpolator_(5))
    self%interpolator_(interpolantMetallicity%ID)=interpolator(self%metallicity                 )
    self%interpolator_(interpolantDensity    %ID)=interpolator(self%densityHydrogen             )
    self%interpolator_(interpolantHydrogen   %ID)=interpolator(self%ionizingFluxHydrogen        )
    self%interpolator_(interpolantHelium     %ID)=interpolator(self%ionizingFluxHeliumToHydrogen)
    self%interpolator_(interpolantOxygen     %ID)=interpolator(self%ionizingFluxOxygenToHelium  )
    !$ call OMP_Init_Lock(self%interpolateLock)
    ! Construct name and description.
    self%name_       ="luminosityEmissionLine:"//String_Join(lineNames,"+")
    self%description_="Luminosity of the "     //String_Join(lineNames,"+")//" emission line"
    if (size(lineNames) > 1) self%description_=self%description_//"s"
    self%description_=self%description_//" [ergs/s]"
    return
  end function lmnstyEmssnLinePanuzzo2003ConstructorInternal

  subroutine lmnstyEmssnLinePanuzzo2003Destructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorLmnstyEmssnLinePanuzzo2003} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLmnstyEmssnLinePanuzzo2003), intent(inout) :: self

    !$ call OMP_Destroy_Lock(self%interpolateLock)
    !![
    <objectDestructor name="self%starFormationRateDisks_"       />
    <objectDestructor name="self%starFormationRateSpheroids_"   />
    <objectDestructor name="self%stellarSpectraDustAttenuation_"/>
    <objectDestructor name="self%outputTimes_"                  />
    !!]
    return
  end subroutine lmnstyEmssnLinePanuzzo2003Destructor

  double precision function lmnstyEmssnLinePanuzzo2003Extract(self,node,instance)
    !!{
    Implement an emission line output analysis property extractor.
    !!}
    use            :: Abundances_Structure            , only : abundances         , max                  , metallicityTypeLogarithmicByMassSolar
    use            :: Galacticus_Nodes                , only : nodeComponentBasic , nodeComponentDisk    , nodeComponentSpheroid                , treeNode
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : hydrogenByMassSolar, luminosityZeroPointAB, massSolar                            , megaParsec, &
          &                                                    metallicitySolar   , parsec
    use            :: Numerical_Constants_Atomic      , only : atomicMassHydrogen , atomicMassUnit
    use            :: Numerical_Constants_Math        , only : Pi
    use            :: Numerical_Constants_Physical    , only : plancksConstant
    use            :: Numerical_Constants_Prefixes    , only : centi              , hecto                , mega
    use            :: Stellar_Luminosities_Structure  , only : max                , stellarLuminosities
    implicit none
    class           (nodePropertyExtractorLmnstyEmssnLinePanuzzo2003), intent(inout), target   :: self
    type            (treeNode                                       ), intent(inout), target   :: node
    type            (multiCounter                                   ), intent(inout), optional :: instance
    class           (nodeComponentBasic                             ), pointer                 :: basic
    class           (nodeComponentDisk                              ), pointer                 :: disk
    class           (nodeComponentSpheroid                          ), pointer                 :: spheroid
    double precision                                                 , parameter               :: massMinimum                   =1.0d-06
    double precision                                                 , parameter               :: radiusMinimum                 =1.0d-06
    double precision                                                 , parameter               :: rateStarFormationMinimum      =1.0d-06
    double precision                                                 , parameter               :: luminosityIonizingMinimum     =1.0d-20
    double precision                                                 , parameter               :: massHIIRegion                 =7.5d+03                     ! Mass of gas in HII region; M☉.
    double precision                                                 , parameter               :: massGMC                       =3.7d+07                     ! Mass of a giant molecular cloud at critical surface density; M☉.
    double precision                                                 , parameter               :: lifetimeHIIRegion             =1.0d-03                     ! Lifetime of HII region; Gyr.
    double precision                                                 , parameter               :: efficiencyHIIRegion           =1.0d-02                     ! Efficiency of HII region (fraction of mass turned into stars).
    double precision                                                 , parameter               :: densitySurfaceCritical        =8.5d+13                     ! Critical surface density for molecular clouds; M☉ Mpc⁻².
    double precision                                                 , parameter               :: metallicityISMLocal           =+2.00d-02                   ! Metallicity in the local ISM.
    double precision                                                 , parameter               :: AVToEBV                       =+3.10d+00                   ! (A_V/E(B-V); Savage & Mathis 1979)
    double precision                                                 , parameter               :: NHToEBV                       =+5.80d+21                   ! (N_H/E(B-V); atoms/cm²/mag; Savage & Mathis 1979)
    double precision                                                 , parameter               :: wavelengthZeroPoint           =+5.50d+03                   ! Angstroms
    double precision                                                 , parameter               :: depthOpticalToMagnitudes      =+2.50d+00                 & ! Conversion factor from optical depth to magnitudes of extinction.
         &                                                                                                                       *log10(                   &
         &                                                                                                                              +exp(              &
         &                                                                                                                                   +1.0d0        &
         &                                                                                                                                  )              &
         &                                                                                                                             )
    double precision                                                 , parameter               :: depthOpticalNormalization     =+AVToEBV                  &
         &                                                                                                                       /NHToEBV                  &
         &                                                                                                                       *hydrogenByMassSolar      &
         &                                                                                                                       /atomicMassUnit*massSolar &
         &                                                                                                                       /(                        &
         &                                                                                                                         +parsec                 &
         &                                                                                                                         *hecto                  &
         &                                                                                                                       )**2                      &
         &                                                                                                                       /metallicityISMLocal      &
         &                                                                                                                       /depthOpticalToMagnitudes
    type            (stellarLuminosities                            ), dimension(  2  )        :: luminositiesStellar
    type            (abundances                                     ), dimension(  2  )        :: abundancesGas
    double precision                                                 , dimension(3,2  )        :: luminosityIonizing
    double precision                                                 , dimension(  2  )        :: massGas                                                 , radius                       , &
         &                                                                                        rateStarFormation                                       , metallicityGas               , &
         &                                                                                        densityHydrogen                                         , luminosityLymanContinuum     , &
         &                                                                                        ratioLuminosityHeliumToHydrogen                         , ratioLuminosityOxygenToHelium, &
         &                                                                                        countHIIRegion                                          , densitySurfaceGas            , &
         &                                                                                        massClouds                                              , densitySurfaceClouds         , &
         &                                                                                        depthOpticalDiffuse                                     , densitySurfaceMetals         , &
         &                                                                                        ionizingFluxMultiplier
    logical                                                          , dimension(  2  )        :: isPhysical
    integer         (c_size_t                                       ), dimension(0:1,5)        :: interpolateIndex
    double precision                                                 , dimension(0:1,5)        :: interpolateFactor
    double precision                                                                           :: weight                                                  , luminosityLinePerHIIRegion
    integer         (c_size_t                                       )                          :: output
    integer                                                                                    :: component                                               , continuum                    , &
         &                                                                                        i                                                       , j                            , &
         &                                                                                        k                                                       , l                            , &
         &                                                                                        m                                                       , line
    !$GLC attributes unused :: instance

    ! Retrieve components.
    basic    => node%basic   ()
    disk     => node%disk    ()
    spheroid => node%spheroid()
    ! Determine output index.
    output   =  self%outputTimes_%index(basic%time(),findClosest=.true.)
    ! Extract all required properties.
    luminositiesStellar(componentDisk    %ID)=disk    %luminositiesStellar             (    )
    luminositiesStellar(componentSpheroid%ID)=spheroid%luminositiesStellar             (    )
    abundancesGas      (componentDisk    %ID)=disk    %abundancesGas                   (    )
    abundancesGas      (componentSpheroid%ID)=spheroid%abundancesGas                   (    )
    massGas            (componentDisk    %ID)=disk    %massGas                         (    )
    massGas            (componentSpheroid%ID)=spheroid%massGas                         (    )
    radius             (componentDisk    %ID)=disk    %radius                          (    )
    radius             (componentSpheroid%ID)=spheroid%radius                          (    )
    rateStarFormation  (componentDisk    %ID)=self    %starFormationRateDisks_    %rate(node)
    rateStarFormation  (componentSpheroid%ID)=self    %starFormationRateSpheroids_%rate(node)
    ! Extract ionizing continuum luminosities.
    do component=1,2
       do continuum=1,3
          luminosityIonizing(continuum,component)=luminositiesStellar(component)%luminosity(self%ionizingContinuumIndex(output,continuum))
       end do
    end do
    ! Determine if component is physically reasonable.
    isPhysical= massGas                                            > massMinimum               &
         &     .and.                                                                           &
         &      radius                                             > radiusMinimum             &
         &     .and.                                                                           &
         &      rateStarFormation                                  > rateStarFormationMinimum  &
         &     .and.                                                                           &
         &      luminosityIonizing(ionizingContinuumHydrogen%ID,:) > luminosityIonizingMinimum
    ! Convert ionizing continuum luminosities from AB units to units of photons s⁻¹.
    forall(continuum=1:3)
       luminosityIonizing(continuum,:)=+luminosityIonizing(continuum,:)     &
            &                          *luminosityZeroPointAB               &
            &                          /plancksConstant                     &
            &                          *log(                                &
            &                               +self%filterExtent(2,continuum) &
            &                               /self%filterExtent(1,continuum) &
            &                              )
    end forall
    !  Compute the logarithmic metallicity of the gas in each component in Solar units.
    do component=1,2
       if (isPhysical(component)) then
          call abundancesGas(component)%massToMassFraction(massGas(component))
          metallicityGas(component)=abundancesGas(component)%metallicity(metallicityTypeLogarithmicByMassSolar)
       else
          metallicityGas(component)=0.0d0
       end if
    end do
    ! Compute the (logarithm of) hydrogen density, based on the model of Krumholz, McKee, & Tumlinson (2009) for molecular cloud
    ! properties. The assumption is that clouds have a mean surface density, Σ, which is the larger of the gas surface density and
    ! the critical density. Assuming spherical clouds, with a mass, M, that scales with the gas density, then we can determine
    ! their radius, r, from that fact that Σ=M/πr², and from this determine their volume density ρ=3M/4πr³.
    densitySurfaceGas=+0.0d0
    where (isPhysical)
       densitySurfaceGas              =+massGas                            &
            &                          /2.0d0                              &
            &                          /Pi                                 &
            &                          /radius                     **2
       massClouds                     =+massGMC                            &
            &                          *densitySurfaceGas                  &
            &                          /densitySurfaceCritical
       densitySurfaceClouds           =+max(                               &
            &                               +densitySurfaceCritical,       &
            &                               +densitySurfaceGas             &
            &                              )
       densityHydrogen                =+log10(                             &
            &                                 +3.0d0                       &
            &                                 /4.0d0                       &
            &                                 *sqrt(Pi       )             &
            &                                 /sqrt(massClouds)            &
            &                                 *densitySurfaceClouds**1.5d0 &
            &                                 /megaParsec          **3     &
            &                                 *centi               **3     &
            &                                 *hydrogenByMassSolar         &
            &                                 *massSolar                   &
            &                                 /atomicMassUnit              &
            &                                 /atomicMassHydrogen          &
            &                                )
       ! Compute logarithm of Lyman continuum luminosity.
       luminosityLymanContinuum       =+log10(                                                                                    luminosityIonizing(ionizingContinuumHydrogen%ID,:)                           )
       ! Compute helium to Lyman continuum luminosity logarithmic ratio.
       ratioLuminosityHeliumToHydrogen=+log10(max(luminosityIonizing(ionizingContinuumHelium%ID,:),luminosityIonizingMinimum)/    luminosityIonizing(ionizingContinuumHydrogen%ID,:)                           )
       ! Compute oxygen to helium continuum luminosity logarithmic ratio.
       ratioLuminosityOxygenToHelium  =+log10(max(luminosityIonizing(ionizingContinuumOxygen%ID,:),luminosityIonizingMinimum)/max(luminosityIonizing(ionizingContinuumHelium  %ID,:),luminosityIonizingMinimum))
       ! Compute number of HII regions.
       countHIIRegion                 =+rateStarFormation   &
            &                          *lifetimeHIIRegion   &
            &                          /massHIIRegion       &
            &                          /efficiencyHIIRegion
       !  Convert the hydrogen ionizing luminosity to be per HII region.
       luminosityLymanContinuum       =+luminosityLymanContinuum &
            &                          -log10(countHIIRegion)
    elsewhere
       ! For non-physical components set the properties to zero. This avoids attempts to use uninitialized values in what follows.
       luminosityLymanContinuum       =0.0d0
       densityHydrogen                =0.0d0
       metallicityGas                 =0.0d0
       ratioLuminosityHeliumToHydrogen=0.0d0
       ratioLuminosityOxygenToHelium  =0.0d0
    end where
    ! Truncate properties to table bounds where necessary to avoid unphysical extrapolations.
    !
    !! Hydrogen density and H-ionizing flux are truncated at both the lower and upper extent
    !! of the tabulated range. The assumption is that the table covers the plausible
    !! physical range for these values. In the case of H-ionizing flux, if we truncate the
    !! value we must then apply a multiplicative correction to the line luminosity to ensure
    !! that we correctly account for all ionizing photons produced.
    !!
    !! For metallicity and the He/H and O/He ionizing flux ratios we truncate only at the
    !! upper extent of the tabulated range. The table is assumed to be tabulated up to the
    !! maximum physically plausible extent for these quantities. Extrapolation to lower
    !! values should be reasonably robust (and the table is assumed to extend to
    !! sufficiently low values that the consequences of extrapolation are unlikely to be
    !! observationally relevant anyway).
    where     (luminosityLymanContinuum        < self%ionizingFluxHydrogen        (                                     1 ))
       ionizingFluxMultiplier         =10.0d0**(                                                            &
            &                                   -self%ionizingFluxHydrogen(                             1 ) &
            &                                   +luminosityLymanContinuum                                   &
            &                                  )
       luminosityLymanContinuum       =self%ionizingFluxHydrogen(                                             1 )
    elsewhere (luminosityLymanContinuum        < self%ionizingFluxHydrogen        (size(self%ionizingFluxHydrogen        )))
       ionizingFluxMultiplier         =10.0d0**(                                                            &
            &                                   -self%ionizingFluxHydrogen(size(self%ionizingFluxHydrogen)) &
            &                                   +luminosityLymanContinuum                                   &
            &                                  )
       luminosityLymanContinuum       =self%ionizingFluxHydrogen        (size(self%ionizingFluxHydrogen        ))
    elsewhere
       ionizingFluxMultiplier=1.0d0
    end where
    where     (densityHydrogen                 < self%densityHydrogen             (                                     1 ))
       densityHydrogen                =self%densityHydrogen             (                                     1 )
    elsewhere (densityHydrogen                < self%densityHydrogen              (size(self%densityHydrogen             )))
       densityHydrogen                =self%densityHydrogen             (size(self%densityHydrogen             ))
    end where
    where     (metallicityGas                  > self%metallicity                 (size(self%metallicity                 )))
       metallicityGas                 =self%metallicity                 (size(self%metallicity                 ))
    end where
    where     (ratioLuminosityHeliumToHydrogen > self%ionizingFluxHeliumToHydrogen(size(self%ionizingFluxHeliumToHydrogen)))
       ratioLuminosityHeliumToHydrogen=self%ionizingFluxHeliumToHydrogen(size(self%ionizingFluxHeliumToHydrogen))
    end where
    where     (ratioLuminosityOxygenToHelium   > self%ionizingFluxOxygenToHelium  (size(self%ionizingFluxOxygenToHelium  )))
       ratioLuminosityOxygenToHelium  =self%ionizingFluxOxygenToHelium  (size(self%ionizingFluxOxygenToHelium  ))
    end where
    ! Perform dust calculation if necessary.
    if (self%depthOpticalISMCoefficient > 0.0d0) then
       where (isPhysical)
          ! Compute surface densities of metals in units of M☉ pc⁻².
          densitySurfaceMetals           =+10.0d0**metallicityGas       &
               &                          *        metallicitySolar     &
               &                          *        densitySurfaceGas    &
               &                          /        mega             **2
          ! Compute optical depth of diffuse dust.
          depthOpticalDiffuse            =+self%depthOpticalISMCoefficient &
               &                          *     depthOpticalNormalization  &
               &                          *     densitySurfaceMetals
       end where
    else
       depthOpticalDiffuse            =+0.0d0
    end if
    ! Iterate over components.
    lmnstyEmssnLinePanuzzo2003Extract=0.0d0
    do component=1,2
       if (.not.isPhysical(component)) cycle
       ! Find interpolating factors in all five interpolants, preventing extrapolation beyond the tabulated ranges.
       !$ call OMP_Set_Lock  (self%interpolateLock)
       call self%interpolator_(interpolantMetallicity%ID)%linearFactors(metallicityGas                 (component),interpolateIndex(0,interpolantMetallicity%ID),interpolateFactor(:,interpolantMetallicity%ID))
       call self%interpolator_(interpolantDensity    %ID)%linearFactors(densityHydrogen                (component),interpolateIndex(0,interpolantDensity    %ID),interpolateFactor(:,interpolantDensity    %ID))
       call self%interpolator_(interpolantHydrogen   %ID)%linearFactors(luminosityLymanContinuum       (component),interpolateIndex(0,interpolantHydrogen   %ID),interpolateFactor(:,interpolantHydrogen   %ID))
       call self%interpolator_(interpolantHelium     %ID)%linearFactors(ratioLuminosityHeliumToHydrogen(component),interpolateIndex(0,interpolantHelium     %ID),interpolateFactor(:,interpolantHelium     %ID))
       call self%interpolator_(interpolantOxygen     %ID)%linearFactors(ratioLuminosityOxygenToHelium  (component),interpolateIndex(0,interpolantOxygen     %ID),interpolateFactor(:,interpolantOxygen     %ID))
       !$ call OMP_Unset_Lock(self%interpolateLock)
       interpolateIndex (1,:                     )=interpolateIndex(0,:)+1
       interpolateFactor=max(min(interpolateFactor,1.0d0),0.0d0)
       ! Iterate over lines.
       do line=1,size(self%luminosity,dim=6)
          ! Interpolate in all five interpolants.
          luminosityLinePerHIIRegion=0.0d0
          do i=0,1
             do j=0,1
                do k=0,1
                   do l=0,1
                      do m=0,1
                         weight                    =+                interpolateFactor(i,1)  &
                              &                     *                interpolateFactor(j,2)  &
                              &                     *                interpolateFactor(k,3)  &
                              &                     *                interpolateFactor(l,4)  &
                              &                     *                interpolateFactor(m,5)
                         luminosityLinePerHIIRegion=+luminosityLinePerHIIRegion              &
                              &                     +weight                                  &
                              &                     *self%luminosity(                        &
                              &                                      interpolateIndex (m,5), &
                              &                                      interpolateIndex (l,4), &
                              &                                      interpolateIndex (k,3), &
                              &                                      interpolateIndex (j,2), &
                              &                                      interpolateIndex (i,1), &
                              &                                      line                    &
                              &                                     )
                      end do
                   end do
                end do
             end do
          end do
          ! Compute the final luminosity in ergs s⁻¹.
          lmnstyEmssnLinePanuzzo2003Extract=+        lmnstyEmssnLinePanuzzo2003Extract                                                                                   &
               &                 +10.0d0**luminosityLinePerHIIRegion                                                                               &
               &                 *        ionizingFluxMultiplier                                                                      (component)  &
               &                 *        countHIIRegion                                                                              (component)  &
               &                 *exp(                                                                                                             &
               &                      -self%stellarSpectraDustAttenuation_%attenuation(                                                            &
               &                                                                       wavelength      =self               %wavelength(     line), &
               &                                                                       age             =0.0d0                                    , &
               &                                                                       vBandAttenuation=depthOpticalDiffuse           (component)  &
               &                                                                      )                                                            &
               &                     )
       end do
    end do
    return
  end function lmnstyEmssnLinePanuzzo2003Extract


  function lmnstyEmssnLinePanuzzo2003Quantity(self)
    !!{
    Return the class of the emission line luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityLuminosity
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType  )                :: lmnstyEmssnLinePanuzzo2003Quantity
    class(nodePropertyExtractorLmnstyEmssnLinePanuzzo2003), intent(inout) :: self
    !$GLC attributes unused :: self

    lmnstyEmssnLinePanuzzo2003Quantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function lmnstyEmssnLinePanuzzo2003Quantity

  function lmnstyEmssnLinePanuzzo2003Name(self)
    !!{
    Return the name of the lmnstyEmssnLinePanuzzo2003 property.
    !!}
    implicit none
    type (varying_string                                 )                :: lmnstyEmssnLinePanuzzo2003Name
    class(nodePropertyExtractorLmnstyEmssnLinePanuzzo2003), intent(inout) :: self

    lmnstyEmssnLinePanuzzo2003Name=self%name_
    return
  end function lmnstyEmssnLinePanuzzo2003Name

  function lmnstyEmssnLinePanuzzo2003Description(self)
    !!{
    Return a description of the lmnstyEmssnLinePanuzzo2003 property.
    !!}
    implicit none
    type (varying_string                                 )                :: lmnstyEmssnLinePanuzzo2003Description
    class(nodePropertyExtractorLmnstyEmssnLinePanuzzo2003), intent(inout) :: self

    lmnstyEmssnLinePanuzzo2003Description=self%description_
    return
  end function lmnstyEmssnLinePanuzzo2003Description

  double precision function lmnstyEmssnLinePanuzzo2003UnitsInSI(self)
    !!{
    Return the units of the lmnstyEmssnLinePanuzzo2003 property in the SI system.
    !!}
    use :: Numerical_Constants_Units, only : ergs
    implicit none
    class(nodePropertyExtractorLmnstyEmssnLinePanuzzo2003), intent(inout) :: self
    !$GLC attributes unused :: self

    lmnstyEmssnLinePanuzzo2003UnitsInSI=ergs
    return
  end function lmnstyEmssnLinePanuzzo2003UnitsInSI
