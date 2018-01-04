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

!% Contains a module which implements a stellar luminosity output analysis property extractor class which applies the dust model of \cite{charlot_simple_2000}.

  use ISO_Varying_String

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorLmnstyStllrCF2000">
  !#  <description>A stellar luminosity output analysis property extractor class which applies the dust model of \cite{charlot_simple_2000}.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorLmnstyStllrCF2000
     !% A stellar luminosity output analysis property extractor class which applies the dust model of \cite{charlot_simple_2000}.
     private
     type            (varying_string)                            :: filterName                , filterType
     double precision                                            :: redshiftBand              , wavelengthFilterEffective    , &
          &                                                         depthOpticalISMCoefficient, depthOpticalCloudsCoefficient, &
          &                                                         wavelengthExponent
     integer                         , allocatable, dimension(:) :: luminosityIndex           , luminosityRecentIndex
   contains
     procedure :: extract  => lmnstyStllrChrltFll2000Extract
     procedure :: type     => lmnstyStllrChrltFll2000Type
     procedure :: quantity => lmnstyStllrChrltFll2000Quantity
  end type outputAnalysisPropertyExtractorLmnstyStllrCF2000

  interface outputAnalysisPropertyExtractorLmnstyStllrCF2000
     !% Constructors for the ``lmnstyStllrChrltFll2000'' output analysis class.
     module procedure lmnstyStllrChrltFll2000ConstructorParameters
     module procedure lmnstyStllrChrltFll2000ConstructorInternal
  end interface outputAnalysisPropertyExtractorLmnstyStllrCF2000

contains

  function lmnstyStllrChrltFll2000ConstructorParameters(parameters) result(self)
    !% Constructor for the ``lmnstyStllrChrltFll2000'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisPropertyExtractorLmnstyStllrCF2000)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    type            (varying_string                                  )                :: filterName                   , filterType
    double precision                                                                  :: redshiftBand                 , depthOpticalISMCoefficient, &
         &                                                                               depthOpticalCloudsCoefficient, wavelengthExponent
    logical                                                                           :: redshiftBandIsPresent

    redshiftBandIsPresent=parameters%isPresent('redshiftBand'    )
    !# <inputParameter>
    !#   <name>filterName</name>
    !#   <source>parameters</source>
    !#   <description>The filter to select.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>filterType</name>
    !#   <source>parameters</source>
    !#   <description>The filter type (rest or observed) to select.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    if (redshiftBandIsPresent) then
       !# <inputParameter>
       !#   <name>redshiftBand</name>
       !#   <source>parameters</source>
       !#   <description>The redshift of the band (if not the output redshift).</description>
       !#   <type>float</type>
       !#   <cardinality>0..1</cardinality>
       !# </inputParameter>
    end if
    !# <inputParameter>
    !#   <name>depthOpticalISMCoefficient</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Multiplicative coefficient for optical depth in the ISM.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>depthOpticalCloudsCoefficient</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Multiplicative coefficient for optical depth in birth clouds.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>wavelengthExponent</name>
    !#   <defaultValue>0.7d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>Exponent of wavelength in the optical depth.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    if (redshiftBandIsPresent) then
       self=outputAnalysisPropertyExtractorLmnstyStllrCF2000(char(filterName),char(filterType),depthOpticalISMCoefficient,depthOpticalCloudsCoefficient,wavelengthExponent,redshiftBand=redshiftBand)
    else
       self=outputAnalysisPropertyExtractorLmnstyStllrCF2000(char(filterName),char(filterType),depthOpticalISMCoefficient,depthOpticalCloudsCoefficient,wavelengthExponent                          )
    end if
    !# <inputParametersValidate source="parameters"/>
    return
  end function lmnstyStllrChrltFll2000ConstructorParameters

  function lmnstyStllrChrltFll2000ConstructorInternal(filterName,filterType,depthOpticalISMCoefficient,depthOpticalCloudsCoefficient,wavelengthExponent,redshiftBand,outputMask) result(self)
    !% Internal constructor for the ``lmnstyStllrChrltFll2000'' output analysis property extractor class.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Output_Times
    use               Stellar_Luminosities_Structure
    use               Memory_Management
    use               Instruments_Filters
    implicit none
    type            (outputAnalysisPropertyExtractorLmnstyStllrCF2000)                                        :: self
    character       (len=*                                           ), intent(in   )                         :: filterName                , filterType
    double precision                                                  , intent(in   )                         :: depthOpticalISMCoefficient, depthOpticalCloudsCoefficient, &
         &                                                                                                       wavelengthExponent
    double precision                                                  , intent(in   ),               optional :: redshiftBand
    logical                                                           , intent(in   ), dimension(:), optional :: outputMask
    integer         (c_size_t                                        )                                        :: i
    !# <constructorAssign variables="filterName, filterType, redshiftBand, depthOpticalISMCoefficient, depthOpticalCloudsCoefficient, wavelengthExponent"/>

    self%wavelengthFilterEffective=Filter_Wavelength_Effective(Filter_Get_Index(var_str(filterName)))
    call allocateArray(self%luminosityIndex,[Galacticus_Output_Time_Count()])
    call allocateArray(self%luminosityRecentIndex,[Galacticus_Output_Time_Count()])
    do i=1,Galacticus_Output_Time_Count()
       if (present(outputMask).and..not.outputMask(i)) then
          self%luminosityIndex      (i)=-1
          self%luminosityRecentIndex(i)=-1
       else
          self%luminosityIndex      (i)=unitStellarLuminosities%index(filterName,filterType,Galacticus_Output_Redshift(i),redshiftBand,postprocessChain='default')
          self%luminosityRecentIndex(i)=unitStellarLuminosities%index(filterName,filterType,Galacticus_Output_Redshift(i),redshiftBand,postprocessChain='recent' )
       end if
    end do
    return
  end function lmnstyStllrChrltFll2000ConstructorInternal
  
  double precision function lmnstyStllrChrltFll2000Extract(self,node)
    !% Implement a stellar luminosity output analysis property extractor.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Output_Times
    use               Galacticus_Nodes
    use               Abundances_Structure
    use               Stellar_Luminosities_Structure
    use               Numerical_Constants_Atomic
    use               Numerical_Constants_Astronomical
    use               Numerical_Constants_Prefixes
    use               Numerical_Constants_Math
    implicit none
    class           (outputAnalysisPropertyExtractorLmnstyStllrCF2000), intent(inout) :: self
    type            (treeNode                                        ), intent(inout) :: node
    class           (nodeComponentBasic                              ), pointer       :: basic
    class           (nodeComponentDisk                               ), pointer       :: disk
    class           (nodeComponentSpheroid                           ), pointer       :: spheroid
    double precision                                                  , parameter     :: metallicityISMLocal      =+2.00d-02            ! Metallicity in the local ISM.
    double precision                                                  , parameter     :: AVToEBV                  =+3.10d+00            ! (A_V/E(B-V); Savage & Mathis 1979)
    double precision                                                  , parameter     :: NHToEBV                  =+5.80d+21            ! (N_H/E(B-V); atoms/cm^2/mag; Savage & Mathis 1979)
    double precision                                                  , parameter     :: wavelengthZeroPoint      =+5.50d+03            ! Angstroms
    double precision                                                  , parameter     :: depthOpticalToMagnitudes =+2.50d+00          & ! Conversion factor from optical depth to magnitudes of extinction.
         &                                                                                                         *log10(            &
         &                                                                                                                +exp(       &
         &                                                                                                                     +1.0d0 &
         &                                                                                                                    )       &
         &                                                                                                               )
    double precision                                                  , parameter     :: depthOpticalNormalization=+AVToEBV                  &
         &                                                                                                         /NHToEBV                  &
         &                                                                                                         *hydrogenByMassSolar      &
         &                                                                                                         /atomicMassUnit*massSolar &
         &                                                                                                         /(                        &
         &                                                                                                           +parsec                 &
         &                                                                                                           *hecto                  &
         &                                                                                                         )**2                      &
         &                                                                                                         /metallicityISMLocal      &
         &                                                                                                         /depthOpticalToMagnitudes
    integer         (c_size_t                                        )                :: i
    type            (stellarLuminosities                             )                :: luminositiesStellar
    type            (abundances                                      )                :: abundancesGas
    double precision                                                                  :: luminosityDisk                                     , luminosityDiskRecent        , &
         &                                                                               luminositySpheroid                                 , luminositySpheroidRecent    , &
         &                                                                               metallicityDisk                                    , metallicitySpheroid         , &
         &                                                                               densitySurfaceMetalsDisk                           , densitySurfaceMetalsSpheroid, &
         &                                                                               depthOpticalDiffuseDisk                            , depthOpticalDiffuseSpheroid , &
         &                                                                               depthOpticalCloudsDisk                             , depthOpticalCloudsSpheroid

    ! Extract luminosities and metallicities of disk and spheroid.
    basic                          =>                                  node               %basic              (                             )
    disk                           =>                                  node               %disk               (                             )
    spheroid                       =>                                  node               %spheroid           (                             )
    i                              =  Galacticus_Output_Time_Index    (basic              %time               (                             ))
    luminositiesStellar            =                                   disk               %luminositiesStellar(                             )
    abundancesGas                  =                                   disk               %abundancesGas      (                             )
    call abundancesGas%massToMassFraction(disk    %massGas())
    metallicityDisk                =                                   abundancesGas      %metallicity        (metallicityTypeLinearByMass  )
    abundancesGas                  =                                   spheroid           %abundancesGas      (                             )
    call abundancesGas%massToMassFraction(spheroid%massGas())
    metallicitySpheroid            =                                   abundancesGas      %metallicity        (metallicityTypeLinearByMass  )
    luminosityDisk                 =                                   luminositiesStellar%luminosity         (self%luminosityIndex      (i))
    luminosityDiskRecent           =                                   luminositiesStellar%luminosity         (self%luminosityRecentIndex(i))
    luminositiesStellar            =                                   spheroid           %luminositiesStellar(                             )
    luminositySpheroid             =                                   luminositiesStellar%luminosity         (self%luminosityIndex      (i))
    luminositySpheroidRecent       =                                   luminositiesStellar%luminosity         (self%luminosityRecentIndex(i))
    ! Compute surface densities of metals in units of Msun/pc^2.
    if (disk%radius() > 0.0d0) then
       densitySurfaceMetalsDisk    =+metallicityDisk     &
            &                            * disk    %massGas() &
            &                            /2.0d0               &
            &                            /Pi                  &
            &                            /(                   &
            &                             +mega               &
            &                             *disk    %radius () &
            &                            )**2
    else
       densitySurfaceMetalsDisk    =+0.0d0
    end if
    if (spheroid%radius() > 0.0d0) then
       densitySurfaceMetalsSpheroid=+metallicitySpheroid &
            &                            * spheroid%massGas() &
            &                            /2.0d0               &
            &                            /Pi                  &
            &                            /(                   &
            &                             +mega               &
            &                             *spheroid%radius()  &
            &                            )**2
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

  integer function lmnstyStllrChrltFll2000Type(self)
    !% Return the type of the stellar luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorLmnstyStllrCF2000), intent(inout) :: self
    !GCC$ attributes unused :: self

    lmnstyStllrChrltFll2000Type=outputAnalysisPropertyTypeLinear
    return
  end function lmnstyStllrChrltFll2000Type

  integer function lmnstyStllrChrltFll2000Quantity(self)
    !% Return the class of the stellar luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorLmnstyStllrCF2000), intent(inout) :: self
    !GCC$ attributes unused :: self

    lmnstyStllrChrltFll2000Quantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function lmnstyStllrChrltFll2000Quantity
