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
Implements a stellar luminosity output analysis property extractor class which applies the dust model of \cite{charlot_simple_2000}.
!!}

  use :: ISO_Varying_String, only : varying_string
  use :: Output_Times      , only : outputTimesClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorLmnstyStllrCF2000">
   <description>A stellar luminosity output analysis property extractor class which applies the dust model of \cite{charlot_simple_2000}.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorLmnstyStllrCF2000
     !!{
     A stellar luminosity output analysis property extractor class which applies the dust model of \cite{charlot_simple_2000}.
     !!}
     private
     type            (varying_string  )                            :: filterName                          , filterType                   , &
          &                                                           name_                               , description_
     double precision                                              :: redshiftBand                        , wavelengthFilterEffective    , &
          &                                                           depthOpticalISMCoefficient          , depthOpticalCloudsCoefficient, &
          &                                                           wavelengthExponent
     integer                           , allocatable, dimension(:) :: luminosityIndex                     , luminosityRecentIndex
     class           (outputTimesClass), pointer                   :: outputTimes_               => null()
   contains
     final     ::                lmnstyStllrChrltFll2000Destructor
     procedure :: extract     => lmnstyStllrChrltFll2000Extract
     procedure :: quantity    => lmnstyStllrChrltFll2000Quantity
     procedure :: name        => lmnstyStllrCF2000Name
     procedure :: description => lmnstyStllrCF2000Description
     procedure :: unitsInSI   => lmnstyStllrCF2000UnitsInSI
  end type nodePropertyExtractorLmnstyStllrCF2000

  interface nodePropertyExtractorLmnstyStllrCF2000
     !!{
     Constructors for the \refClass{nodePropertyExtractorLmnstyStllrCF2000} output analysis class.
     !!}
     module procedure lmnstyStllrChrltFll2000ConstructorParameters
     module procedure lmnstyStllrChrltFll2000ConstructorInternal
  end interface nodePropertyExtractorLmnstyStllrCF2000

contains

  function lmnstyStllrChrltFll2000ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorLmnstyStllrCF2000} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorLmnstyStllrCF2000)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (outputTimesClass                      ), pointer       :: outputTimes_
    type            (varying_string                        )                :: filterName                   , filterType
    double precision                                                        :: redshiftBand                 , depthOpticalISMCoefficient, &
         &                                                                     depthOpticalCloudsCoefficient, wavelengthExponent
    logical                                                                 :: redshiftBandIsPresent

    redshiftBandIsPresent=parameters%isPresent('redshiftBand'    )
    !![
    <inputParameter>
      <name>filterName</name>
      <source>parameters</source>
      <description>The filter to select.</description>
    </inputParameter>
    <inputParameter>
      <name>filterType</name>
      <source>parameters</source>
      <description>The filter type (rest or observed) to select.</description>
    </inputParameter>
    !!]
    if (redshiftBandIsPresent) then
       !![
       <inputParameter>
         <name>redshiftBand</name>
         <source>parameters</source>
         <description>The redshift of the band (if not the output redshift).</description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>depthOpticalISMCoefficient</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>Multiplicative coefficient for optical depth in the ISM.</description>
    </inputParameter>
    <inputParameter>
      <name>depthOpticalCloudsCoefficient</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>Multiplicative coefficient for optical depth in birth clouds.</description>
    </inputParameter>
    <inputParameter>
      <name>wavelengthExponent</name>
      <defaultValue>0.7d0</defaultValue>
      <source>parameters</source>
      <description>Exponent of wavelength in the optical depth.</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    if (redshiftBandIsPresent) then
       self=nodePropertyExtractorLmnstyStllrCF2000(char(filterName),char(filterType),depthOpticalISMCoefficient,depthOpticalCloudsCoefficient,wavelengthExponent,outputTimes_,redshiftBand=redshiftBand)
    else
       self=nodePropertyExtractorLmnstyStllrCF2000(char(filterName),char(filterType),depthOpticalISMCoefficient,depthOpticalCloudsCoefficient,wavelengthExponent,outputTimes_                          )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function lmnstyStllrChrltFll2000ConstructorParameters

  function lmnstyStllrChrltFll2000ConstructorInternal(filterName,filterType,depthOpticalISMCoefficient,depthOpticalCloudsCoefficient,wavelengthExponent,outputTimes_,redshiftBand,outputMask) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorLmnstyStllrCF2000} output analysis property extractor class.
    !!}
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Instruments_Filters           , only : Filter_Get_Index       , Filter_Wavelength_Effective
    use            :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    type            (nodePropertyExtractorLmnstyStllrCF2000)                                        :: self
    character       (len=*                                 ), intent(in   )                         :: filterName                , filterType
    double precision                                        , intent(in   )                         :: depthOpticalISMCoefficient, depthOpticalCloudsCoefficient, &
         &                                                                                             wavelengthExponent
    class           (outputTimesClass                      ), intent(in   ), target                 :: outputTimes_
    double precision                                        , intent(in   ),               optional :: redshiftBand
    logical                                                 , intent(in   ), dimension(:), optional :: outputMask
    integer         (c_size_t                              )                                        :: i
    character       (len=7                                 )                                        :: label
    !![
    <constructorAssign variables="filterName, filterType, redshiftBand, depthOpticalISMCoefficient, depthOpticalCloudsCoefficient, wavelengthExponent, *outputTimes_"/>
    !!]

    self%wavelengthFilterEffective=Filter_Wavelength_Effective(Filter_Get_Index(var_str(filterName)))
    allocate(self%luminosityIndex      (self%outputTimes_%count()))
    allocate(self%luminosityRecentIndex(self%outputTimes_%count()))
    do i=1,self%outputTimes_%count()
       if (present(outputMask).and..not.outputMask(i)) then
          self%luminosityIndex      (i)=-1
          self%luminosityRecentIndex(i)=-1
       else
          self%luminosityIndex      (i)=unitStellarLuminosities%index(filterName,filterType,self%outputTimes_%redshift(i),redshiftBand,postprocessChain='default')
          self%luminosityRecentIndex(i)=unitStellarLuminosities%index(filterName,filterType,self%outputTimes_%redshift(i),redshiftBand,postprocessChain='recent' )
       end if
    end do
    self%name_       ="luminosityStellar:"//filterName//":"//filterType
    self%description_="Total stellar luminosity luminosity in the "//filterType//"-frame "//filterName//" filter"
    if (present(redshiftBand)) then
       write (label,'(f7.3)') redshiftBand
       self%name_      =self%name_        //":z"            //trim(adjustl(label))
       self%description_=self%description_//" shifted to z="//trim(adjustl(label))
    end if
    self%name_       =self%name_       //":dustCharlotFall2000"
    self%description_=self%description_//" with Charlot & Fall (2000) dust extinction."
    return
  end function lmnstyStllrChrltFll2000ConstructorInternal

  subroutine lmnstyStllrChrltFll2000Destructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorLmnstyStllrCF2000} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLmnstyStllrCF2000), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine lmnstyStllrChrltFll2000Destructor

  double precision function lmnstyStllrChrltFll2000Extract(self,node,instance)
    !!{
    Implement a stellar luminosity output analysis property extractor.
    !!}
    use            :: Abundances_Structure            , only : abundances         , metallicityTypeLinearByMass
    use            :: Galacticus_Nodes                , only : nodeComponentBasic , nodeComponentDisk          , nodeComponentSpheroid, treeNode
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : hydrogenByMassSolar, massSolar                  , parsec
    use            :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use            :: Numerical_Constants_Math        , only : Pi
    use            :: Numerical_Constants_Prefixes    , only : hecto              , mega
    use            :: Stellar_Luminosities_Structure  , only : stellarLuminosities
    implicit none
    class           (nodePropertyExtractorLmnstyStllrCF2000), intent(inout), target   :: self
    type            (treeNode                              ), intent(inout), target   :: node
    type            (multiCounter                          ), intent(inout), optional :: instance
    class           (nodeComponentBasic                    ), pointer                 :: basic
    class           (nodeComponentDisk                     ), pointer                 :: disk
    class           (nodeComponentSpheroid                 ), pointer                 :: spheroid
    double precision                                        , parameter               :: metallicityISMLocal      =+2.00d-02            ! Metallicity in the local ISM.
    double precision                                        , parameter               :: AVToEBV                  =+3.10d+00            ! (A_V/E(B-V); Savage & Mathis 1979)
    double precision                                        , parameter               :: NHToEBV                  =+5.80d+21            ! (N_H/E(B-V); atoms/cm²/mag; Savage & Mathis 1979)
    double precision                                        , parameter               :: wavelengthZeroPoint      =+5.50d+03            ! Angstroms
    double precision                                        , parameter               :: depthOpticalToMagnitudes =+2.50d+00          & ! Conversion factor from optical depth to magnitudes of extinction.
         &                                                                                                         *log10(            &
         &                                                                                                                +exp(       &
         &                                                                                                                     +1.0d0 &
         &                                                                                                                    )       &
         &                                                                                                               )
    double precision                                        , parameter               :: depthOpticalNormalization=+AVToEBV                  &
         &                                                                                                         /NHToEBV                  &
         &                                                                                                         *hydrogenByMassSolar      &
         &                                                                                                         /atomicMassUnit*massSolar &
         &                                                                                                         /(                        &
         &                                                                                                           +parsec                 &
         &                                                                                                           *hecto                  &
         &                                                                                                         )**2                      &
         &                                                                                                         /metallicityISMLocal      &
         &                                                                                                         /depthOpticalToMagnitudes
    integer         (c_size_t                              )                          :: i
    type            (stellarLuminosities                   )                          :: luminositiesStellar
    type            (abundances                            )                          :: abundancesGas
    double precision                                                                  :: luminosityDisk                                     , luminosityDiskRecent        , &
         &                                                                               luminositySpheroid                                 , luminositySpheroidRecent    , &
         &                                                                               metallicityDisk                                    , metallicitySpheroid         , &
         &                                                                               densitySurfaceMetalsDisk                           , densitySurfaceMetalsSpheroid, &
         &                                                                               depthOpticalDiffuseDisk                            , depthOpticalDiffuseSpheroid , &
         &                                                                               depthOpticalCloudsDisk                             , depthOpticalCloudsSpheroid
    !$GLC attributes unused :: instance

    ! Extract luminosities and metallicities of disk and spheroid.
    basic                          =>                         node               %basic              (                             )
    disk                           =>                         node               %disk               (                             )
    spheroid                       =>                         node               %spheroid           (                             )
    i                              =  self%outputTimes_%index(basic              %time               (                             ),findClosest=.true.)
    luminositiesStellar            =                          disk               %luminositiesStellar(                             )
    abundancesGas                  =                          disk               %abundancesGas      (                             )
    call abundancesGas%massToMassFraction(disk    %massGas())
    metallicityDisk                =                          abundancesGas      %metallicity        (metallicityTypeLinearByMass  )
    abundancesGas                  =                          spheroid           %abundancesGas      (                             )
    call abundancesGas%massToMassFraction(spheroid%massGas())
    metallicitySpheroid            =                          abundancesGas      %metallicity        (metallicityTypeLinearByMass  )
    luminosityDisk                 =                          luminositiesStellar%luminosity         (self%luminosityIndex      (i))
    luminosityDiskRecent           =                          luminositiesStellar%luminosity         (self%luminosityRecentIndex(i))
    luminositiesStellar            =                          spheroid           %luminositiesStellar(                             )
    luminositySpheroid             =                          luminositiesStellar%luminosity         (self%luminosityIndex      (i))
    luminositySpheroidRecent       =                          luminositiesStellar%luminosity         (self%luminosityRecentIndex(i))
    ! Compute surface densities of metals in units of M☉/pc².
    if     (                            &
         &   disk    %radius () > 0.0d0 &
         &  .and.                       &
         &   disk    %massGas() > 0.0d0 &
         & ) then
       densitySurfaceMetalsDisk    =+metallicityDisk      &
            &                       *  disk    %massGas() &
            &                       /2.0d0                &
            &                       /Pi                   &
            &                       /(                    &
            &                         +mega               &
            &                         *disk    %radius () &
            &                        )**2
    else
       densitySurfaceMetalsDisk    =+0.0d0
    end if
    if     (                            &
         &   spheroid%radius () > 0.0d0 &
         &  .and.                       &
         &   spheroid%massGas() > 0.0d0 &
         & ) then
       densitySurfaceMetalsSpheroid=+metallicitySpheroid  &
            &                       *  spheroid%massGas() &
            &                       /2.0d0                &
            &                       /Pi                   &
            &                       /(                    &
            &                         +mega               &
            &                         *spheroid%radius () &
            &                        )**2
    else
       densitySurfaceMetalsSpheroid=+0.0d0
    end if
    ! Compute optical depth of diffuse dust.
    depthOpticalDiffuseDisk    =+  self%depthOpticalISMCoefficient   &
         &                      *depthOpticalNormalization           &
         &                      *densitySurfaceMetalsDisk            &
         &                      /(                                   &
         &                        +self%wavelengthFilterEffective    &
         &                        /wavelengthZeroPoint               &
         &                      )**self%wavelengthExponent
    depthOpticalDiffuseSpheroid=+  self%depthOpticalISMCoefficient   &
         &                      *depthOpticalNormalization           &
         &                      *densitySurfaceMetalsSpheroid        &
         &                      /(                                   &
         &                        +self%wavelengthFilterEffective    &
         &                        /wavelengthZeroPoint               &
         &                      )**self%wavelengthExponent
    ! Compute optical depth in clouds.
    depthOpticalCloudsDisk    =+  self%depthOpticalCloudsCoefficient &
         &                     *metallicityDisk                      &
         &                     /metallicityISMLocal                  &
         &                     /(                                    &
         &                       +self%wavelengthFilterEffective     &
         &                       /wavelengthZeroPoint                &
         &                     )**self%wavelengthExponent
    depthOpticalCloudsSpheroid=+  self%depthOpticalCloudsCoefficient &
         &                     *metallicitySpheroid                  &
         &                     /metallicityISMLocal                  &
         &                     /(                                    &
         &                       +self%wavelengthFilterEffective     &
         &                       /wavelengthZeroPoint                &
         &                     )**self%wavelengthExponent
    ! Compute the attenuated luminosity.
    lmnstyStllrChrltFll2000Extract=+(                                   &
         &                           +(                                 &
         &                             +luminosityDisk                  &
         &                             -luminosityDiskRecent            &
         &                            )                                 &
         &                           +  luminosityDiskRecent            &
         &                           *exp(-depthOpticalCloudsDisk     ) &
         &                          )                                   &
         &                         *  exp(-depthOpticalDiffuseDisk    ) &
         &                         +(                                   &
         &                           +(                                 &
         &                             +luminositySpheroid              &
         &                             -luminositySpheroidRecent        &
         &                            )                                 &
         &                           +  luminositySpheroidRecent        &
         &                           *exp(-depthOpticalCloudsSpheroid ) &
         &                          )                                   &
         &                         *  exp(-depthOpticalDiffuseSpheroid)
    return
  end function lmnstyStllrChrltFll2000Extract


  function lmnstyStllrChrltFll2000Quantity(self)
    !!{
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityLuminosity
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: lmnstyStllrChrltFll2000Quantity
    class(nodePropertyExtractorLmnstyStllrCF2000       ), intent(inout) :: self
    !$GLC attributes unused :: self

    lmnstyStllrChrltFll2000Quantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function lmnstyStllrChrltFll2000Quantity

  function lmnstyStllrCF2000Name(self)
    !!{
    Return the name of the lmnstyStllrCF2000 property.
    !!}
    implicit none
    type (varying_string                        )                :: lmnstyStllrCF2000Name
    class(nodePropertyExtractorLmnstyStllrCF2000), intent(inout) :: self

    lmnstyStllrCF2000Name=self%name_
    return
  end function lmnstyStllrCF2000Name

  function lmnstyStllrCF2000Description(self)
    !!{
    Return a description of the lmnstyStllrCF2000 property.
    !!}
    implicit none
    type (varying_string                        )                :: lmnstyStllrCF2000Description
    class(nodePropertyExtractorLmnstyStllrCF2000), intent(inout) :: self

    lmnstyStllrCF2000Description=self%description_
    return
  end function lmnstyStllrCF2000Description

  double precision function lmnstyStllrCF2000UnitsInSI(self)
    !!{
    Return the units of the lmnstyStllrCF2000 property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminosityZeroPointAB
    implicit none
    class(nodePropertyExtractorLmnstyStllrCF2000), intent(inout) :: self
    !$GLC attributes unused :: self

    lmnstyStllrCF2000UnitsInSI=luminosityZeroPointAB
    return
  end function lmnstyStllrCF2000UnitsInSI
