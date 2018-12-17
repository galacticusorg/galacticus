!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a stellar mass output analysis property extractor class.

  !$ use OMP_Lib
  use    ISO_Varying_String
  use    FGSL                             , only : fgsl_interp_accel
  use    Stellar_Spectra_Dust_Attenuations
  use    Output_Times

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorLmnstyEmssnLine" defaultThreadPrivate="yes">
  !#  <description>A stellar luminosity output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorLmnstyEmssnLine
     !% A stellar luminosity output analysis property extractor class.
     private
     class           (stellarSpectraDustAttenuationClass), pointer                             :: stellarSpectraDustAttenuation_
     class           (outputTimesClass                  ), pointer                             :: outputTimes_
     type            (varying_string                    ), allocatable, dimension(:          ) :: lineNames
     double precision                                    , allocatable, dimension(:          ) :: metallicity                   , densityHydrogen             , &
          &                                                                                       ionizingFluxHydrogen          , ionizingFluxHeliumToHydrogen, &
          &                                                                                       ionizingFluxOxygenToHydrogen  , wavelength
     double precision                                    , allocatable, dimension(:,:,:,:,:,:) :: luminosity
     integer                                             , allocatable, dimension(:,:        ) :: ionizingContinuumIndex
     double precision                                                 , dimension(2,3        ) :: filterExtent
     type            (fgsl_interp_accel                  )            , dimension(5          ) :: accelerator
     logical                                                          , dimension(5          ) :: interpolateReset
     double precision                                                                          :: depthOpticalISMCoefficient
     !$ integer      (omp_lock_kind                      )                                     :: interpolateLock
   contains
     final     ::             lmnstyEmssnLineDestructor
     procedure :: extract  => lmnstyEmssnLineExtract
     procedure :: type     => lmnstyEmssnLineType
     procedure :: quantity => lmnstyEmssnLineQuantity
  end type outputAnalysisPropertyExtractorLmnstyEmssnLine

  interface outputAnalysisPropertyExtractorLmnstyEmssnLine
     !% Constructors for the ``lmnstyEmssnLine'' output analysis class.
     module procedure lmnstyEmssnLineConstructorParameters
     module procedure lmnstyEmssnLineConstructorInternal
  end interface outputAnalysisPropertyExtractorLmnstyEmssnLine

  ! Enumerations for galactic components and ionizing continuua.
  !# <enumeration>
  !#  <name>galacticComponent</name>
  !#  <description>Specifies the galactic component for emission line calculations.</description>
  !#  <indexing>1</indexing>
  !#  <entry label="disk"    />
  !#  <entry label="spheroid"/>
  !# </enumeration>
  !# <enumeration>
  !#  <name>ionizingContinuum</name>
  !#  <description>Specifies the ionizing continuum for emission line calculations.</description>
  !#  <indexing>1</indexing>
  !#  <entry label="Hydrogen"/>
  !#  <entry label="Helium"  />
  !#  <entry label="Oxygen"  />
  !# </enumeration>
  !# <enumeration>
  !#  <name>interpolant</name>
  !#  <description>Specifies the different interpolants for emission line calculations.</description>
  !#  <indexing>1</indexing>
  !#  <entry label="metallicity"/>
  !#  <entry label="density"  />
  !#  <entry label="hydrogen"  />
  !#  <entry label="helium"  />
  !#  <entry label="oxygen"  />
  !# </enumeration>

contains
  
  function lmnstyEmssnLineConstructorParameters(parameters) result(self)
    !% Constructor for the ``lmnstyEmssnLine'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisPropertyExtractorLmnstyEmssnLine)                              :: self
    type            (inputParameters                               ), intent(inout)               :: parameters
    type            (varying_string                                ), allocatable  , dimension(:) :: lineNames
    class           (stellarSpectraDustAttenuationClass            ), pointer                     :: stellarSpectraDustAttenuation_
    class           (outputTimesClass                              ), pointer                     :: outputTimes_
    double precision                                                                              :: depthOpticalISMCoefficient
    
    allocate(lineNames(parameters%count('lineNames')))
    !# <inputParameter>
    !#   <name>lineNames</name>
    !#   <source>parameters</source>
    !#   <description>The emission lines to extract.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>depthOpticalISMCoefficient</name>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Multiplicative coefficient for optical depth in the ISM.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="stellarSpectraDustAttenuation" name="stellarSpectraDustAttenuation_" source="parameters"/>
    !# <objectBuilder class="outputTimes"                   name="outputTimes_"                   source="parameters"/>
    self=outputAnalysisPropertyExtractorLmnstyEmssnLine(stellarSpectraDustAttenuation_,outputTimes_,lineNames,depthOpticalISMCoefficient)
    !# <inputParametersValidate source="parameters"/>
    return
  end function lmnstyEmssnLineConstructorParameters

  function lmnstyEmssnLineConstructorInternal(stellarSpectraDustAttenuation_,outputTimes_,lineNames,depthOpticalISMCoefficient,outputMask) result(self)
    !% Internal constructor for the ``lmnstyEmssnLine'' output analysis property extractor class.
    use, intrinsic :: ISO_C_Binding
    use               Instruments_Filters
    use               Output_Times
    use               Galacticus_Paths
    use               Memory_Management
    use               Stellar_Luminosities_Structure
    use               IO_HDF5
    use               Galacticus_Error
    implicit none
    type            (outputAnalysisPropertyExtractorLmnstyEmssnLine)                                        :: self
    double precision                                                , intent(in   )                         :: depthOpticalISMCoefficient
    type            (varying_string                                ), intent(in   ), dimension(:)           :: lineNames
    logical                                                         , intent(in   ), dimension(:), optional :: outputMask
    class           (stellarSpectraDustAttenuationClass            ), intent(in   ), target                 :: stellarSpectraDustAttenuation_
    class           (outputTimesClass                              ), intent(in   ), target                 :: outputTimes_
    type            (hdf5Object                                    )                                        :: emissionLinesFile             , lines, &
         &                                                                                                     lineDataset
    integer         (c_size_t                                      )                                        :: i
    !# <constructorAssign variables="lineNames, depthOpticalISMCoefficient, *stellarSpectraDustAttenuation_, *outputTimes_"/>

    ! Read the table of emission line luminosities.
    !$ call hdf5Access%set()
    call emissionLinesFile%openFile(char(galacticusPath(pathTypeDataStatic))//"hiiRegions/emissionLines.hdf5",readOnly=.true.)
    lines=emissionLinesFile%openGroup('lines')
    do i=1,size(lineNames)
       if (.not.lines%hasDataset(char(self%lineNames(i)))) call Galacticus_Error_Report('line "'//char(self%lineNames(i))//'" not found'//{introspection:location})
    end do
    call emissionLinesFile%readDataset('metallicity'                  ,self%metallicity                 )
    call emissionLinesFile%readDataset('densityHydrogen'              ,self%densityHydrogen             )
    call emissionLinesFile%readDataset('ionizingFluxHydrogen'         ,self%ionizingFluxHydrogen        )
    call emissionLinesFile%readDataset('ionizingFluxHeliumToHydrogen' ,self%ionizingFluxHeliumToHydrogen)
    call emissionLinesFile%readDataset('ionizingFluxOxygenToHydrogen' ,self%ionizingFluxOxygenToHydrogen)
    call allocateArray(                                          &
         &             self%luminosity,                          &
         &             [                                         &
         &              size(self%ionizingFluxOxygenToHydrogen), &
         &              size(self%ionizingFluxHeliumToHydrogen), &
         &              size(self%ionizingFluxHydrogen        ), &
         &              size(self%densityHydrogen             ), &
         &              size(self%metallicity                 ), &
         &              size(self%lineNames                   )  &
         &             ]                                         &
         &            )
    call allocateArray(                                          &
         &             self%wavelength,                          &
         &             [                                         &
         &              size(self%lineNames                   )  &
         &             ]                                         &
         &            )
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
    self%ionizingFluxOxygenToHydrogen=log10(self%ionizingFluxOxygenToHydrogen)
    self%luminosity                  =log10(self%luminosity                  )
    ! Find indices of ionizing continuua filters.
    call allocateArray(self%ionizingContinuumIndex,[self%outputTimes_%count(),3_c_size_t])
    do i=1,self%outputTimes_%count()
       if (present(outputMask).and..not.outputMask(i)) then
          self%ionizingContinuumIndex(i,:                        )=-1
       else
          self%ionizingContinuumIndex(i,ionizingContinuumHydrogen)=unitStellarLuminosities%index('Lyc'            ,'rest',self%outputTimes_%redshift(i))
          self%ionizingContinuumIndex(i,ionizingContinuumHelium  )=unitStellarLuminosities%index('HeliumContinuum','rest',self%outputTimes_%redshift(i))
          self%ionizingContinuumIndex(i,ionizingContinuumOxygen  )=unitStellarLuminosities%index('OxygenContinuum','rest',self%outputTimes_%redshift(i))
       end if
    end do
    ! Read wavelength intervals of ionizing continuum filters.
    self%filterExtent(:,ionizingContinuumHydrogen)=Filter_Extent(Filter_Get_Index(var_str('Lyc'            )))
    self%filterExtent(:,ionizingContinuumHelium  )=Filter_Extent(Filter_Get_Index(var_str('HeliumContinuum')))
    self%filterExtent(:,ionizingContinuumOxygen  )=Filter_Extent(Filter_Get_Index(var_str('OxygenContinuum')))
    ! Initialize interpolations.
    self%interpolateReset=.true.
    !$ call OMP_Init_Lock(self%interpolateLock)
    return
  end function lmnstyEmssnLineConstructorInternal

  subroutine lmnstyEmssnLineDestructor(self)
    !% Destructor for the ``lmnstyEmssnLine'' output analysis property extractor class.
    use Numerical_Interpolation
    implicit none
    type   (outputAnalysisPropertyExtractorLmnstyEmssnLine), intent(inout) :: self
    integer                                                                :: i

    do i=1,5
       call Interpolate_Done(                                                   &
            &                interpolationAccelerator=self%accelerator     (i), &
            &                reset                   =self%interpolateReset(i)  &
            &               )
    end do
    !$ call OMP_Destroy_Lock(self%interpolateLock)
    !# <objectDestructor name="self%stellarSpectraDustAttenuation_"/>
    !# <objectDestructor name="self%outputTimes_"                  />
    return
  end subroutine lmnstyEmssnLineDestructor

  double precision function lmnstyEmssnLineExtract(self,node)
    !% Implement an emission line output analysis property extractor.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Nodes                , only : nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroid
    use               Stellar_Luminosities_Structure
    use               Numerical_Constants_Physical
    use               Numerical_Constants_Astronomical
    use               Numerical_Constants_Atomic
    use               Numerical_Constants_Prefixes
    use               Numerical_Interpolation
    use               Abundances_Structure
    implicit none
    class           (outputAnalysisPropertyExtractorLmnstyEmssnLine), intent(inout)    :: self
    type            (treeNode                                      ), intent(inout)    :: node
    class           (nodeComponentBasic                            ), pointer          :: basic
    class           (nodeComponentDisk                             ), pointer          :: disk
    class           (nodeComponentSpheroid                         ), pointer          :: spheroid
    double precision                                                , parameter        :: massMinimum                   =1.0d-06
    double precision                                                , parameter        :: radiusMinimum                 =1.0d-06
    double precision                                                , parameter        :: rateStarFormationMinimum      =1.0d-06
    double precision                                                , parameter        :: luminosityIonizingMinimum     =1.0d-20
    double precision                                                , parameter        :: massHIIRegion                 =7.5d+03                     ! Mass of gas in HII region; M☉.
    double precision                                                , parameter        :: massGMC                       =3.7d+07                     ! Mass of a giant molecular cloud at critical surface density; M☉.
    double precision                                                , parameter        :: lifetimeHIIRegion             =1.0d-03                     ! Lifetime of HII region; Gyr.
    double precision                                                , parameter        :: densitySurfaceCritical        =8.5d+13                     ! Critical surface density for molecular clouds; M☉ Mpc⁻³
    double precision                                                , parameter        :: metallicityISMLocal           =+2.00d-02                   ! Metallicity in the local ISM.
    double precision                                                , parameter        :: AVToEBV                       =+3.10d+00                   ! (A_V/E(B-V); Savage & Mathis 1979)
    double precision                                                , parameter        :: NHToEBV                       =+5.80d+21                   ! (N_H/E(B-V); atoms/cm^2/mag; Savage & Mathis 1979)
    double precision                                                , parameter        :: wavelengthZeroPoint           =+5.50d+03                   ! Angstroms
    double precision                                                , parameter        :: depthOpticalToMagnitudes      =+2.50d+00                 & ! Conversion factor from optical depth to magnitudes of extinction.
         &                                                                                                               *log10(                   &
         &                                                                                                                      +exp(              &
         &                                                                                                                           +1.0d0        &
         &                                                                                                                          )              &
         &                                                                                                                     )
    double precision                                                , parameter        :: depthOpticalNormalization     =+AVToEBV                  &
         &                                                                                                               /NHToEBV                  &
         &                                                                                                               *hydrogenByMassSolar      &
         &                                                                                                               /atomicMassUnit*massSolar &
         &                                                                                                               /(                        &
         &                                                                                                                 +parsec                 &
         &                                                                                                                 *hecto                  &
         &                                                                                                               )**2                      &
         &                                                                                                               /metallicityISMLocal      &
         &                                                                                                               /depthOpticalToMagnitudes
    type            (stellarLuminosities                           ), dimension(  2  ) :: luminositiesStellar
    type            (abundances                                    ), dimension(  2  ) :: abundancesGas
    double precision                                                , dimension(3,2  ) :: luminosityIonizing
    double precision                                                , dimension(  2  ) :: massGas                                                 , radius                         , &
         &                                                                                rateStarFormation                                       , metallicityGas                 , &
         &                                                                                densityHydrogen                                         , luminosityLymanContinuum       , &
         &                                                                                ratioLuminosityHeliumToHydrogen                         , ratioLuminosityOxygenToHydrogen, &
         &                                                                                countHIIRegion                                          , densitySurfaceGas              , &
         &                                                                                massClouds                                              , densitySurfaceClouds           , &
         &                                                                                depthOpticalDiffuse                                     , densitySurfaceMetals
    logical                                                         , dimension(  2  ) :: isPhysical
    integer         (c_size_t                                      ), dimension(0:1,5) :: interpolateIndex
    double precision                                                , dimension(0:1,5) :: interpolateFactor
    double precision                                                                   :: weight                                                  , luminosityLinePerHIIRegion
    integer         (c_size_t                                      )                   :: output
    integer                                                                            :: component                                               , continuum                      , &
         &                                                                                i                                                       , j                              , &
         &                                                                                k                                                       , l                              , &
         &                                                                                m                                                       , line

    ! Retrieve components.
    basic    => node%basic   ()
    disk     => node%disk    ()
    spheroid => node%spheroid()
    ! Determine output index.
    output   =  self%outputTimes_%index(basic%time())
    ! Extract all required properties.
    luminositiesStellar(galacticComponentDisk    )=disk    %luminositiesStellar()
    luminositiesStellar(galacticComponentSpheroid)=spheroid%luminositiesStellar()
    abundancesGas      (galacticComponentDisk    )=disk    %abundancesGas      ()
    abundancesGas      (galacticComponentSpheroid)=spheroid%abundancesGas      ()
    massGas            (galacticComponentDisk    )=disk    %massGas            ()
    massGas            (galacticComponentSpheroid)=spheroid%massGas            ()
    radius             (galacticComponentDisk    )=disk    %radius             ()
    radius             (galacticComponentSpheroid)=spheroid%radius             ()
    rateStarFormation  (galacticComponentDisk    )=disk    %starFormationRate  ()
    rateStarFormation  (galacticComponentSpheroid)=spheroid%starFormationRate  ()
    ! Extract ionizing continuum luminosities.
    do component=1,2
       do continuum=1,3
          luminosityIonizing(continuum,component)=luminositiesStellar(component)%luminosity(self%ionizingContinuumIndex(output,continuum))
       end do
    end do
    ! Determine if component is physically reasonable.
    isPhysical= massGas                                         > massMinimum               &
         &     .and.                                                                        &
         &      radius                                          > radiusMinimum             &
         &     .and.                                                                        &
         &      rateStarFormation                               > rateStarFormationMinimum  &
         &     .and.                                                                        &
         &      luminosityIonizing(ionizingContinuumHydrogen,:) > luminosityIonizingMinimum
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
    ! Compute the (logarithm of) hydrogen density, based on the model of Krumholz, McKee, & Tumlinson (2009) for molecular cloud properties.
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
       luminosityLymanContinuum       =+log10(                                                                             luminosityIonizing(ionizingContinuumHydrogen,:))
       ! Compute helium to Lyman continuum luminosity logarithmic ratio.
       ratioLuminosityHeliumToHydrogen=+log10(max(luminosityIonizing(ionizingContinuumHelium,:),luminosityIonizingMinimum)/luminosityIonizing(ionizingContinuumHydrogen,:))
       ! Compute oxygen to Lyman continuum luminosity logarithmic ratio.
       ratioLuminosityOxygenToHydrogen=+log10(max(luminosityIonizing(ionizingContinuumOxygen,:),luminosityIonizingMinimum)/luminosityIonizing(ionizingContinuumHydrogen,:))
       ! Compute number of HII regions.
       countHIIRegion                 =+rateStarFormation &
            &                          *lifetimeHIIRegion &
            &                          /massHIIRegion
       !  Convert the hydrogen ionizing luminosity to be per HII region.
       luminosityLymanContinuum       =+luminosityLymanContinuum &
            &                          -log10(countHIIRegion)
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
    lmnstyEmssnLineExtract=0.0d0
    do component=1,2
       if (.not.isPhysical(component)) cycle
       ! Find interpolating factors in all five interpolants, preventing extrapolation beyond the tabulated ranges.
       !$ call OMP_Set_Lock  (self%interpolateLock)
       interpolateIndex (0,interpolantMetallicity)=Interpolate_Locate                 (self%metallicity                 ,self            %accelerator(  interpolantMetallicity),metallicityGas                 (component),reset=self%interpolateReset(interpolantMetallicity))
       interpolateIndex (0,interpolantDensity    )=Interpolate_Locate                 (self%densityHydrogen             ,self            %accelerator(  interpolantDensity    ),densityHydrogen                (component),reset=self%interpolateReset(interpolantDensity    ))
       interpolateIndex (0,interpolantHydrogen   )=Interpolate_Locate                 (self%ionizingFluxHydrogen        ,self            %accelerator(  interpolantHydrogen   ),luminosityLymanContinuum       (component),reset=self%interpolateReset(interpolantHydrogen   ))
       interpolateIndex (0,interpolantHelium     )=Interpolate_Locate                 (self%ionizingFluxHeliumToHydrogen,self            %accelerator(  interpolantHelium     ),ratioLuminosityHeliumToHydrogen(component),reset=self%interpolateReset(interpolantHelium     ))
       interpolateIndex (0,interpolantOxygen     )=Interpolate_Locate                 (self%ionizingFluxOxygenToHydrogen,self            %accelerator(  interpolantOxygen     ),ratioLuminosityOxygenToHydrogen(component),reset=self%interpolateReset(interpolantOxygen     ))
       !$ call OMP_Unset_Lock(self%interpolateLock)
       interpolateFactor(:,interpolantMetallicity)=Interpolate_Linear_Generate_Factors(self%metallicity                 ,interpolateIndex            (0,interpolantMetallicity),metallicityGas                 (component)                                                    )
       interpolateFactor(:,interpolantDensity    )=Interpolate_Linear_Generate_Factors(self%densityHydrogen             ,interpolateIndex            (0,interpolantDensity    ),densityHydrogen                (component)                                                    )
       interpolateFactor(:,interpolantHydrogen   )=Interpolate_Linear_Generate_Factors(self%ionizingFluxHydrogen        ,interpolateIndex            (0,interpolantHydrogen   ),luminosityLymanContinuum       (component)                                                    )
       interpolateFactor(:,interpolantHelium     )=Interpolate_Linear_Generate_Factors(self%ionizingFluxHeliumToHydrogen,interpolateIndex            (0,interpolantHelium     ),ratioLuminosityHeliumToHydrogen(component)                                                    )
       interpolateFactor(:,interpolantOxygen     )=Interpolate_Linear_Generate_Factors(self%ionizingFluxOxygenToHydrogen,interpolateIndex            (0,interpolantOxygen     ),ratioLuminosityOxygenToHydrogen(component)                                                    )
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
          lmnstyEmssnLineExtract=+        lmnstyEmssnLineExtract                                                                                   &
               &                 +10.0d0**luminosityLinePerHIIRegion                                                                               &
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
  end function lmnstyEmssnLineExtract
  
  integer function lmnstyEmssnLineType(self)
    !% Return the type of the emission line luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorLmnstyEmssnLine), intent(inout) :: self
    !GCC$ attributes unused :: self

    lmnstyEmssnLineType=outputAnalysisPropertyTypeLinear
    return
  end function lmnstyEmssnLineType

  integer function lmnstyEmssnLineQuantity(self)
    !% Return the class of the emission line luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorLmnstyEmssnLine), intent(inout) :: self
    !GCC$ attributes unused :: self

    lmnstyEmssnLineQuantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function lmnstyEmssnLineQuantity
