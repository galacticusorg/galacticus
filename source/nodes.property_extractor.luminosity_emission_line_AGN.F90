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
  use    :: Numerical_Interpolation          , only : interpolator
  use    :: ISO_Varying_String               , only : varying_string
  !$ use :: OMP_Lib                          , only : omp_lock_kind
  use    :: Output_Times                     , only : outputTimesClass
  use    :: Stellar_Spectra_Dust_Attenuations, only : stellarSpectraDustAttenuationClass
  use    :: Black_Hole_Accretion_Rates       , only : blackHoleAccretionRateClass
  use    :: Accretion_Disks                  , only : accretionDisksClass
  !![
  <nodePropertyExtractor name="nodePropertyExtractorLmnstyEmssnLineAGN">
    <description>
      An emission line luminosity property extractor class. The luminosity of the named emission line (given by the {\normalfont
      \ttfamily lineNames} parameter: if multiple lines are named, the sum of their luminosities) is computed. }
      parameter.
    </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLmnstyEmssnLineAGN
     !!{
     A stellar luminosity output analysis property extractor class.
     !!}
     private
     class(accretionDisksClass        ), pointer :: accretionDisks_         => null()
     class(blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
     class(outputTimesClass           ), pointer :: outputTimes_            => null()
     type            (varying_string                    ), allocatable, dimension(:          ) :: lineNames, names_, descriptions_
     integer                                                                                   :: countLines
     integer         (c_size_t                         )                                       :: indexDensityHydrogen, indexIonizationParameter, indexMetallicity,indexSpectralIndex
     double precision                                    , allocatable, dimension(:          ) :: metallicity                             , densityHydrogen             , &
          &                                                                                       spectralIndex                    , ionizationParameter, &
          &                                                                                       wavelengths
     double precision                                    , allocatable, dimension(:,:,:,:,:) :: luminosity
     double precision                                                 , dimension(2,3        ) :: filterExtent
     type            (interpolator                      ), allocatable, dimension(:          ) :: interpolator_
     double precision                                                                          :: alpha, volume_filling_factor
     !$ integer      (omp_lock_kind                     )                                      :: interpolateLock
   contains
     final     ::                lmnstyEmssnLineAGNDestructor
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

  ! Enumerations for galactic components and spectral index.
  !![
  <enumeration>
   <name>interpolants</name>
   <description>Specifies the different interpolants for AGN emission line calculations.</description>
   <indexing>1</indexing>
   <entry label="density"/>
   <entry label="ionizationParameter"/>
   <entry label="metallicity"/>
   <entry label="spectralIndex"/>
  </enumeration>
  !!]
contains
  function lmnstyEmssnLineAGNConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``lmnstyEmssnLineAGN'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorLmnstyEmssnLineAGN        )                              :: self
    type            (inputParameters                                ), intent(inout)               :: parameters
    type            (varying_string                                 ), allocatable  , dimension(:) :: lineNames
    class           (outputTimesClass                               ), pointer                     :: outputTimes_
    class           (accretionDisksClass                            ), pointer                     :: accretionDisks_
    class(blackHoleAccretionRateClass                               ), pointer       :: blackHoleAccretionRate_
    double precision                                                                               :: alpha, volume_filling_factor
    allocate(lineNames(parameters%count('lineNames')))
    !![
    <inputParameter>
      <name>lineNames</name>
      <source>parameters</source>
      <description>The emission lines to extract.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha</name>
      <defaultValue>-1.7d0</defaultValue>
      <source>parameters</source>
      <description>Multiplicative coefficient for optical depth in the ISM.</description>
    </inputParameter>
    <inputParameter>
      <name>volume_filling_factor</name>
      <defaultValue>0.01d0</defaultValue>
      <source>parameters</source>
      <description> Ratio of the volume-averaged hydrogen density to hydrogen density </description>
    </inputParameter>
    <objectBuilder class="accretionDisks"        name="accretionDisks_"         source="parameters"/>
   <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"            source="parameters"/>
    !!]
    self=nodePropertyExtractorLmnstyEmssnLineAGN(accretionDisks_,blackHoleAccretionRate_,outputTimes_,lineNames,alpha,volume_filling_factor)
    !![
    <inputParametersValidate source="parameters"    />
    <objectDestructor name="accretionDisks_"        />
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="outputTimes_"           />
    !!]
    return
  end function lmnstyEmssnLineAGNConstructorParameters

  function lmnstyEmssnLineAGNConstructorInternal(accretionDisks_,blackHoleAccretionRate_,outputTimes_,lineNames,alpha,volume_filling_factor,outputMask) result(self)
    !!{
    Internal constructor for the ``lmnstyEmssnLineAGN'' output analysis property extractor class.
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
    use            :: Galacticus_Nodes              , only : nodeComponentBlackHole
     use           :: Table_Labels                  , only : extrapolationTypeFix
    implicit none
    type            (nodePropertyExtractorLmnstyEmssnLineAGN)                                                :: self
    double precision                                                 , intent(in   )                         :: alpha, volume_filling_factor
    type            (varying_string                                 ), intent(in   ), dimension(:)           :: lineNames
    logical                                                          , intent(in   ), dimension(:), optional :: outputMask
    class           (accretionDisksClass                            ), intent(in   ), target                 :: accretionDisks_
    class           (blackHoleAccretionRateClass                    ), intent(in   ), target                 :: blackHoleAccretionRate_
    class           (outputTimesClass                               ), intent(in   ), target                 :: outputTimes_
    type            (hdf5Object                                 )                                            :: emissionLinesFile, lines, lineDataset, dataset
    integer                                                                                                  :: i,k, countBlackHoles
    integer         (c_size_t)                                               , dimension(5    )                        :: shapeLines, permutation
    double precision                                               , allocatable  , dimension(:,:,:,:,:)     :: luminosity
    !![
    <constructorAssign variables="lineNames, alpha, volume_filling_factor, *accretionDisks_, *blackHoleAccretionRate_, *outputTimes_"/>
    !!]
    ! Read the table of emission line luminosities.
    !$ call hdf5Access%set()
    call emissionLinesFile%openFile(char(inputPath(pathTypeDataStatic))//"hiiRegions/emissionLineLuminosities_AGN.hdf5",readOnly=.true.)
    lines=emissionLinesFile%openGroup('lines')
    do i=1,size(lineNames)
       if (.not.lines%hasDataset(char(self%lineNames(i)))) call Error_Report('line "'//char(self%lineNames(i))//'" not found'//{introspection:location})
    end do
    call emissionLinesFile%readDataset('densityHydrogen'              ,self%densityHydrogen             )
    call emissionLinesFile%readDataset('ionizationParameter'          ,self%ionizationParameter         )
    call emissionLinesFile%readDataset('metallicity'                  ,self%metallicity                 )
    call emissionLinesFile%readDataset('spectralIndex'                ,self%spectralIndex               )
    ! Extract indexing into the lines arrays.
    dataset=emissionLinesFile%openDataset('densityHydrogen'               )
    call dataset%readAttribute('index',self%indexDensityHydrogen                                        )
    call dataset%close        (                     )
    dataset=emissionLinesFile%openDataset('ionizationParameter'                                         )
    call dataset%readAttribute('index',self%indexIonizationParameter                                    )
    call dataset%close        (                     )
    dataset=emissionLinesFile%openDataset('metallicity'                                                 )
    call dataset%readAttribute('index',self%indexMetallicity                                            )
    call dataset%close        (                     )
    dataset=emissionLinesFile%openDataset('spectralIndex'                                               )
    call dataset%readAttribute('index',self%indexSpectralIndex                                          )
    call dataset%close        (                     )
    ! Offset indexing to Fortran standard (i.e. starting from 1 instead of 0).
    self%indexDensityHydrogen           =self%indexDensityHydrogen               +1
    self%indexIonizationParameter       =self%indexIonizationParameter           +1
    self%indexMetallicity               =self%indexMetallicity                   +1
    self%indexSpectralIndex             =self%indexSpectralIndex                 +1
     ! Establish arrays.
    shapeLines(self%indexDensityHydrogen           )=size(self%densityHydrogen           )
    shapeLines(self%indexIonizationParameter       )=size(self%ionizationParameter       )
    shapeLines(self%indexMetallicity               )=size(self%metallicity               )
    shapeLines(self%indexSpectralIndex             )=size(self%spectralIndex             )
    shapeLines(     5                              )=size(     lineNames                 )
    ! Allocate a temporary luminosities array into which we will read data from the table. The dimension ordering here is whatever
    ! ordering was used in the table file. This will be reordered into our preferred, internal order later.
    !![
    <allocate variable="luminosity" shape="shapeLines"/>
    !!]
    allocate(                                          &
         &   self%wavelengths                          &
         &   (                                         &
         &    size(self%lineNames                   )  &
         &   )                                         &
         &  )
    do i=1,size(self%lineNames)
       call lines%readDatasetStatic(char(self%lineNames(i)),luminosity(:,:,:,:,i))
       lineDataset=lines%openDataset(char(self%lineNames(i)))
       call lineDataset%readAttribute('wavelength',self%wavelengths(i))
       call lineDataset%close        (                               )
    end do
    call lines            %close      (                                                                 )
    call emissionLinesFile%close      (                                                                 )
    !$ call hdf5Access%unset()

    ! Re-order the luminosities table into our preferred order.
    !! First, allocate our final table array with our preferred ordering of dimensions.
    allocate(                                          &
         &   self%luminosity                           &
         &   (                                         &
         &    size(self%spectralIndex               ), &
         &    size(self%metallicity                 ), &
         &    size(self%ionizationParameter         ), &
         &    size(self%densityHydrogen             ), &
         &    size(self%lineNames                   )  &
         &   )                                         &
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
    
    ! Convert parameters and luminosities to log form.
    self%densityHydrogen             =log10(self%densityHydrogen             )
    self%ionizationParameter         =log10(self%ionizationParameter         )
    ! Cloudy table gives metallicity Z (not in log)
    self%metallicity                 =log10(self%metallicity                 )
    
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
       self%  names_(i)="luminosityEmissionLineAGN:"//lineNames(i)
       self   %descriptions_(i)="Luminosity of the "             //lineNames(i)//" AGN emission line [ergs/s]"
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
    <objectDestructor name="self%accretionDisks_"       />
    <objectDestructor name="self%blackHoleAccretionRate_"   />
    <objectDestructor name="self%outputTimes_"                  />
    !!]
    return
  end subroutine lmnstyEmssnLineAGNDestructor

  function lmnstyEmssnLineAGNExtract(self,node,time,instance)
    !!{
    Implement an emission line output analysis property extractor.
    !!}
    use            :: Abundances_Structure            , only : abundances         , max                  , metallicityTypeLogarithmicByMassSolar
    use            :: Galacticus_Nodes                , only : nodeComponentBasic , nodeComponentDisk    , nodeComponentSpheroid, treeNode, nodeComponentBlackHole
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : hydrogenByMassSolar, luminosityZeroPointAB, massSolar                            , megaParsec, &
          &                                                    metallicitySolar   , parsec
    use            :: Numerical_Constants_Atomic      , only : atomicMassHydrogen , atomicMassUnit
    use            :: Numerical_Constants_Math        , only : Pi
    use            :: Numerical_Constants_Physical    , only : plancksConstant, speedLight
    use            :: Numerical_Constants_Prefixes    , only : centi              , hecto                , mega, micro
    use            :: Numerical_Constants_Astronomical, only : massSolar, gigaYear, parsec
    use            :: Stellar_Luminosities_Structure  , only : max                , stellarLuminosities
    use :: Galactic_Structure_Options, only : massTypeStellar
    use :: Mass_Distributions            , only : massDistributionClass
    implicit none
    double precision                                              , dimension(:) , allocatable :: lmnstyEmssnLineAGNExtract
    class           (nodePropertyExtractorLmnstyEmssnLineAGN        ), intent(inout), target   :: self
    type            (treeNode                                       ), intent(inout), target   :: node
    double precision                                                 , intent(in   )           :: time
    type            (multiCounter                                   ), intent(inout), optional :: instance
    class           (nodeComponentBasic                             ), pointer                 :: basic
    class           (nodeComponentDisk                              ), pointer                 :: disk
    class           (nodeComponentSpheroid                          ), pointer                 :: spheroid
    class           (nodeComponentBlackHole                         ), pointer                 :: blackHole
    class           (massDistributionClass                          ), pointer                 :: massDistribution_              
    double precision                                                 , parameter               :: massMinimum                   =1.0d-06
    double precision                                                 , parameter               :: radiusMinimum                 =1.0d-06
    double precision                                                 , parameter               :: rateStarFormationMinimum      =1.0d-06
    double precision                                                 , parameter               :: luminosityIonizingMinimum     =1.0d-20
    double precision                                                 , parameter               :: metallicityISMLocal           =+2.00d-02  ! Metallicity in the local ISM.
    double precision                                                 , parameter               :: wavelengthZeroPoint           =+5.50d+03  ! Angstroms
    type            (stellarLuminosities                            ), dimension(  2  )        :: luminositiesStellar
    type            (abundances                                     ),                         :: abundancesGas
    double precision                                                 , dimension(3,2  )        :: luminosityIonizing                                                                                                                                      
    logical                                                          ,                         :: isPhysical
    integer         (c_size_t                                       ), dimension(0:1,4)        :: interpolateIndex
    double precision                                                 , dimension(0:1,4)        :: interpolateFactor
    double precision                                                                           :: weight                                                 
    double precision                                                                           :: integral, recombinationCoefficient, rateMassAccretionSpheroid, rateMassAccretionHotHalo, wavelength_,nu_001,nu_091,nu_25,nu_10,L_agn,&
                    &                                                                             hydrogenDensity,spectralIndex ,normalizationConstant,radiativeEfficiency,blackHoleAccretionRate,rateIonizingPhotonsNormalized, rateIonizingPhotons, ionizationParam,&
                    &                                                                             massGas, radius, rateStarFormation, metallicityGas, stromgren_radius, denom , halfMassRadius, hydrogen_mass                  
    !$GLC attributes unused :: instance
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
    abundancesGas      =disk    %abundancesGas                   (    ) + spheroid%abundancesGas                   (    )
    massGas            =disk    %massGas                         (    ) + spheroid%massGas                         (    )
    radius             =disk    %radius                          (    ) + spheroid%radius                          (    )
    massDistribution_            => node             %massDistribution   (massType      =massTypeStellar)
    halfMassRadius                              =  mega*parsec*massDistribution_%radiusEnclosingMass(massFractional=0.5d0          )
    radius=mega*parsec*radius
    hydrogen_mass = massGas*massSolar
    hydrogenDensity=1000.d0
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Determine if component is physically reasonable.
    isPhysical= massGas                                            > massMinimum               &
         &     .and.                                                                           &
         &      radius                                             > radiusMinimum             
    !  Compute the logarithmic metallicity of the gas in each component in Solar units.
    
    
    call abundancesGas%massToMassFraction(massGas)
    ! Galacticus galaxies gives the metallicity in logrithmic scale
    metallicityGas=abundancesGas%metallicity(metallicityTypeLogarithmicByMassSolar)
    if (isnan(metallicityGas)) then
        metallicityGas=0.0d0
    end if            

    blackHole                => node     %blackHole()
    !limit frequencies in SI units (per seconds)
    nu_001=speedLight/(0.001d0*micro)
    nu_091=speedLight/(0.0912d0*micro)
    nu_25=speedLight/(0.25d0*micro)
    nu_10=speedLight/(10.0d0*micro)

    ! Calculate integral for AGN luminosity from 0 to infinity
    integral=(                                                        &
             &     (nu_25**(self%alpha+0.5d0)) * (nu_10**-2.5d0)      &
             &    *((nu_10**3.0d0)/3.0d0)                             &
             &     )                                                  &
             &    +(                                                  &
             &        (nu_25**(self%alpha+0.25d0)                     &
             &     )                                                  &
             &    *(2*(SQRT(nu_25)   -  SQRT(nu_10)) ) )             &
             &    + ((1.0d0/(self%alpha+1.0d0))* ( (nu_001**(self%alpha+1.0d0)) - (nu_25**(self%alpha+1.0d0)) ) ) 
    
    ! recombination coefficient
    recombinationCoefficient=micro*2.6d-13


    !calculate black hole accretion rate
    ! change to call  self%blackHoleAccretionRate_%rateAccretion(blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo)
    call  self%blackHoleAccretionRate_%rateAccretion(blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo)

    !Black hole accretion rate in kg/s
    blackHoleAccretionRate =(massSolar/gigaYear)* (rateMassAccretionHotHalo + rateMassAccretionSpheroid)

    !Radiative effciency of black hole
    radiativeEfficiency =  self%accretionDisks_%efficiencyRadiative(blackHole,rateMassAccretionSpheroid+rateMassAccretionHotHalo)


    !Luminosity of AGN in J/s
    L_agn=blackHoleAccretionRate*speedLight*speedLight*radiativeEfficiency
    
    !normalization constant for the integral over ionizing spectrum
    normalizationConstant=1/integral
    
    ! Rate of ionizing photons Q_H
    rateIonizingPhotons=L_agn*((nu_001**self%alpha) - (nu_091**self%alpha)) &
            &                     /(plancksConstant*self%alpha)
    
    ! Normalized rate of ionizing photons
    rateIonizingPhotonsNormalized=rateIonizingPhotons*normalizationConstant
    
    ! Calculate stromgren radius
    stromgren_radius = (3*rateIonizingPhotonsNormalized/(4*Pi*((hydrogenDensity*mega)**2.0d0)*self%volume_filling_factor*recombinationCoefficient))**(1.0d0/3.0d0)

    if(L_agn>0.0d0) then
      ionizationParam = rateIonizingPhotonsNormalized/(4*Pi*stromgren_radius*stromgren_radius*hydrogenDensity*mega*speedLight)
      ionizationParam=log10(ionizationParam)
    else
      rateIonizingPhotonsNormalized=0.0d0
      ionizationParam=0.0d0
    end if

    if (isPhysical) then
      hydrogenDensity=log10(1000.0000000000000d0)
    else 
      hydrogenDensity=0.0d0   
      ionizationParam=0.0d0        
    end if                                                        
    ! Truncate properties to table bounds where necessary to avoid unphysical extrapolations.
    if (metallicityGas                  < self%metallicity (1)) then
      metallicityGas = self%metallicity(1)
    else if     (metallicityGas                  > self%metallicity                 (size(self%metallicity                 ))) then
       metallicityGas                 =self%metallicity                 (size(self%metallicity                 ))
    end if

    if    (ionizationParam < self%ionizationParameter(1)) then
       ionizationParam=self%ionizationParameter(1)
    else if     (ionizationParam > self%ionizationParameter(size(self%ionizationParameter))) then
       ionizationParam=self%ionizationParameter(size(self%ionizationParameter))
    end if
    

    ! Find interpolating factors in all four interpolants, preventing extrapolation beyond the tabulated ranges.
    !$ call OMP_Set_Lock  (self%interpolateLock)
    call self%interpolator_(interpolantsDensity    %ID)%linearFactors(hydrogenDensity  ,interpolateIndex(0,interpolantsDensity    %ID),interpolateFactor(:,interpolantsDensity    %ID))
    call self%interpolator_(interpolantsIonizationParameter%ID)%linearFactors(ionizationParam ,interpolateIndex(0,interpolantsIonizationParameter%ID),interpolateFactor(:,interpolantsIonizationParameter     %ID))
    call self%interpolator_(interpolantsMetallicity%ID)%linearFactors(metallicityGas ,interpolateIndex(0,interpolantsMetallicity%ID),interpolateFactor(:,interpolantsMetallicity%ID))
    call self%interpolator_(interpolantsSpectralIndex%ID)%linearFactors(self%alpha,interpolateIndex(0,interpolantsSpectralIndex%ID),interpolateFactor(:,interpolantsSpectralIndex   %ID))
  
    !$ call OMP_Unset_Lock(self%interpolateLock)
    interpolateIndex (1,:                     )=interpolateIndex(0,:)+1
    interpolateFactor=max(min(interpolateFactor,1.0d0),0.0d0)
    ! Iterate over lines.
    allocate(lmnstyEmssnLineAGNExtract(self%countLines))
    lmnstyEmssnLineAGNExtract=0.0
    do line=1,size(self%luminosity,dim=5)
        ! Interpolate in all four interpolants.
        do i=0,1
          do j=0,1
              do k=0,1
                do l=0,1
                    weight                    =+                interpolateFactor(i,1)   &
                          &                     *                interpolateFactor(j,2)  &
                          &                     *                interpolateFactor(k,3)  &
                          &                     *                interpolateFactor(l,4)  
                    lmnstyEmssnLineAGNExtract(line)=+lmnstyEmssnLineAGNExtract(line)     &
                          &                     +weight                                  &
                          ! intensity given by emission line models of Cloudy (intensity*4pi) is multiplied by stromgren radius**2 (cm^2)
                          &                     *stromgren_radius*stromgren_radius*1.0d4 &
                          &                     *self%luminosity(                        &
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
   return
  end function lmnstyEmssnLineAGNExtract


  integer function lmnstyEmssnLineAGNElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily lmnstyEmssnLineAGN} property extractors.
    !!}
    implicit none
    class     (nodePropertyExtractorLmnstyEmssnLineAGN), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    lmnstyEmssnLineAGNElementCount=self%countLines
    return
  end function lmnstyEmssnLineAGNElementCount

  subroutine lmnstyEmssnLineAGNNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily emissionLines}.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class           (nodePropertyExtractorLmnstyEmssnLineAGN    ), intent(inout)                            :: self
    double precision                                             , intent(in   )                            :: time
    type            (varying_string                             ), intent(inout), dimension(:), allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(self%countLines))
    names=self%names_
    return
  end subroutine lmnstyEmssnLineAGNNames

  subroutine lmnstyEmssnLineAGNDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily emission line luminosity} property.
    !!}
    implicit none
    class           (nodePropertyExtractorLmnstyEmssnLineAGN), intent(inout)                            :: self
    double precision                            , intent(in   )                             :: time
    type            (varying_string            ), intent(inout), dimension(:) , allocatable :: descriptions
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
    double precision                                             , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorLmnstyEmssnLineAGN    ), intent(inout)               :: self
    double precision                                             , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(unitsInSI(self%countLines))
    unitsInSI=ergs
    return
  end function lmnstyEmssnLineAGNUnitsInSI
 
