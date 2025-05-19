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
Implements an intracluster medium cooling power in band property extractor class.
!!}

  use :: Cooling_Functions      , only : coolingFunction    , coolingFunctionClass
  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorICMCoolingPowerInBand">
   <description>An intracluster medium cooling power in band property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorICMCoolingPowerInBand
     !!{
     A property extractor class which extracts the fraction of the ICM cooling power due to emission in a given energy band.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     class           (coolingFunctionClass    ), pointer :: coolingFunction_     => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     double precision                                    :: energyLow                     , energyHigh
     type            (varying_string          )          :: label
   contains
     final     ::                icmCoolingPowerInBandDestructor
     procedure :: extract     => icmCoolingPowerInBandExtract
     procedure :: name        => icmCoolingPowerInBandName
     procedure :: description => icmCoolingPowerInBandDescription
     procedure :: unitsInSI   => icmCoolingPowerInBandUnitsInSI
  end type nodePropertyExtractorICMCoolingPowerInBand

  interface nodePropertyExtractorICMCoolingPowerInBand
     !!{
     Constructors for the {\normalfont \ttfamily icmCoolingPowerInBand} output analysis class.
     !!}
     module procedure icmCoolingPowerInBandConstructorParameters
     module procedure icmCoolingPowerInBandConstructorInternal
  end interface nodePropertyExtractorICMCoolingPowerInBand

contains

  function icmCoolingPowerInBandConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily icmCoolingPowerInBand} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorICMCoolingPowerInBand)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                  ), pointer       :: darkMatterHaloScale_
    class           (coolingFunctionClass                      ), pointer       :: coolingFunction_
    class           (cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_
    double precision                                                            :: energyLow           , energyHigh
    type            (varying_string                            )                :: label


    !![
    <inputParameter>
      <name>energyLow</name>
      <source>parameters</source>
      <description>The minimum energy (in units of keV) for the band.</description>
    </inputParameter>
    <inputParameter>
      <name>energyHigh</name>
      <source>parameters</source>
      <description>The maximum energy (in units of keV) for the band.</description>
    </inputParameter>
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <description>A label to use as a suffix for this property.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="coolingFunction"     name="coolingFunction_"     source="parameters"/>
    !!]
    self=nodePropertyExtractorICMCoolingPowerInBand(energyLow,energyHigh,label,cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="coolingFunction_"    />
    !!]
    return
  end function icmCoolingPowerInBandConstructorParameters

  function icmCoolingPowerInBandConstructorInternal(energyLow,energyHigh,label,cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily icmCoolingPowerInBand} property extractor class.
    !!}
    implicit none
    type            (nodePropertyExtractorICMCoolingPowerInBand)                        :: self
    double precision                                            , intent(in   )         :: energyLow           , energyHigh
    type            (varying_string                            ), intent(in   )         :: label
    class           (cosmologyFunctionsClass                   ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_
    class           (coolingFunctionClass                      ), intent(in   ), target :: coolingFunction_
    !![
    <constructorAssign variables="energyLow, energyHigh, label, *cosmologyFunctions_, *darkMatterHaloScale_, *coolingFunction_"/>
    !!]

    return
  end function icmCoolingPowerInBandConstructorInternal

  subroutine icmCoolingPowerInBandDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily icmCoolingPowerInBand} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorICMCoolingPowerInBand), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%coolingFunction_"    />
    !!]
    return
  end subroutine icmCoolingPowerInBandDestructor

  double precision function icmCoolingPowerInBandExtract(self,node,instance)
    !!{
    Implement an ICM X-ray properties extractor.
    !!}
    use :: Numerical_Constants_Physical, only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Numerical_Constants_Units   , only : electronVolt
    use :: Numerical_Integration       , only : integrator
    use :: Radiation_Fields            , only : radiationFieldCosmicMicrowaveBackground
    use :: Mass_Distributions          , only : massDistributionClass                  , kinematicsDistributionClass
    use :: Galactic_Structure_Options  , only : componentTypeHotHalo                   , massTypeGaseous
    implicit none
    class           (nodePropertyExtractorICMCoolingPowerInBand), intent(inout), target   :: self
    type            (treeNode                                  ), intent(inout), target   :: node
    type            (multiCounter                              ), intent(inout), optional :: instance
    type            (radiationFieldCosmicMicrowaveBackground   ), pointer                 :: radiation_
    class           (massDistributionClass                     ), pointer                 :: massDistribution_
    class           (kinematicsDistributionClass               ), pointer                 :: kinematicsDistribution_
    type            (integrator                                )                          :: integratorTotal        , integratorInBand
    double precision                                                                      :: luminosityTotal        , luminosityInBand
    !$GLC attributes unused :: instance

    ! Initialize radiation field.
    allocate(radiation_)
    !![
    <referenceConstruct object="radiation_" constructor="radiationFieldCosmicMicrowaveBackground(self%cosmologyFunctions_)"/>
    !!]
    ! Get the mass distribution.
    massDistribution_       => node             %massDistribution      (componentTypeHotHalo,massTypeGaseous)
    kinematicsDistribution_ => massDistribution_%kinematicsDistribution(                                    )
    ! Compute luminosity and temperature.
    integratorTotal =integrator                (integrandLuminosityTotal ,toleranceRelative                           =1.0d-3)
    integratorInBand=integrator                (integrandLuminosityInBand,toleranceRelative                           =1.0d-3)
    luminosityTotal =integratorTotal %integrate(0.0d0                    ,self%darkMatterHaloScale_%radiusVirial(node)       )
    luminosityInBand=integratorInBand%integrate(0.0d0                    ,self%darkMatterHaloScale_%radiusVirial(node)       )
    if (luminosityTotal > 0.0d0) then
       icmCoolingPowerInBandExtract=+luminosityInBand &
            &                       /luminosityTotal
    else
       icmCoolingPowerInBandExtract=+0.0d0
    end if
    !![
    <objectDestructor name="radiation_"             />
    <objectDestructor name="massDistribution_"      />
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return

  contains
    
    double precision function integrandLuminosityTotal(radius)
      !!{
      Integrand function used for computing ICM X-ray luminosities.
      !!}
      use :: Abundances_Structure         , only : abundances
      use :: Chemical_Abundances_Structure, only : chemicalAbundances
      implicit none
      double precision                      , intent(in   ) :: radius
      type            (abundances          )                :: abundancesICM
      type            (chemicalAbundances  )                :: densityChemicalICM
      double precision                                      :: numberDensityHydrogen, temperature

      call icmProperties(radius,numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM)
      ! Do not include factors of 4π and unit conversions, as these would cancel out of the final result.
      integrandLuminosityTotal=+radius**2                                                                                                                 &
           &                   *self%coolingFunction_%coolingFunction(node,numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM,radiation_)
      return
    end function integrandLuminosityTotal
    
    double precision function integrandLuminosityInBand(radius)
      !!{
      Integrand function used for computing ICM X-ray luminosities.
      !!}
      use :: Abundances_Structure         , only : abundances
      use :: Chemical_Abundances_Structure, only : chemicalAbundances
      implicit none
      double precision                      , intent(in   ) :: radius
      type            (abundances          )                :: abundancesICM
      type            (chemicalAbundances  )                :: densityChemicalICM
      double precision                                      :: numberDensityHydrogen, temperature
      
      call icmProperties(radius,numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM)
      ! Do not include factors of 4π and unit conversions, as these would cancel out of the final result.
      integrandLuminosityInBand=+radius**2                                                                                                                                                              &
           &                    *self%coolingFunction_%coolingFunction              (node,numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM,radiation_                               ) &
           &                    *self%coolingFunction_%coolingFunctionFractionInBand(node,numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM,radiation_,self%energyLow,self%energyHigh)
      return
    end function integrandLuminosityInBand

    subroutine icmProperties(radius,numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM)
      use :: Abundances_Structure             , only : abundances
      use :: Chemical_Abundances_Structure    , only : chemicalAbundances
      use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Fraction_Conversion
      use :: Galacticus_Nodes                 , only : nodeComponentHotHalo
      use :: Numerical_Constants_Atomic       , only : massHydrogenAtom
      use :: Numerical_Constants_Prefixes     , only : hecto
      use :: Numerical_Constants_Astronomical , only : massSolar                            , megaParsec
      use :: Coordinates                      , only : coordinateSpherical                  , assignment(=)
      implicit none
      double precision                      , intent(in   ) :: radius
      double precision                      , intent(  out) :: numberDensityHydrogen  , temperature
      type            (abundances          ), intent(  out) :: abundancesICM
      type            (chemicalAbundances  ), intent(  out) :: densityChemicalICM
      class           (nodeComponentHotHalo), pointer       :: hotHalo
      type            (chemicalAbundances  )                :: massChemicalICM
      type            (coordinateSpherical )                :: coordinates
      double precision                                      :: density                , massICM     , &
           &                                                   massToDensityConversion

      ! Set the coordinates.
      coordinates     =  [radius,0.0d0,0.0d0]
      ! Get the density of the ICM.
      density         =  massDistribution_      %density    (coordinates)
      ! Get the temperature of the ICM.
      temperature     =  kinematicsDistribution_%temperature(coordinates)
      ! Get abundances and chemistry of the ICM.
      hotHalo         => node   %hotHalo   ()
      massICM         =  hotHalo%mass      ()
      abundancesICM   =  hotHalo%abundances()
      massChemicalICM =  hotHalo%chemicals ()
      call abundancesICM  %massToMassFraction(           massICM)
      call massChemicalICM%massToNumber      (densityChemicalICM)
      ! Compute factor converting mass of chemicals in (M☉) to number density in cm⁻³ per total mass density.
      if (hotHalo%mass() > 0.0d0) then
         massToDensityConversion=Chemicals_Mass_To_Fraction_Conversion(hotHalo%mass())
      else
         massToDensityConversion=0.0d0
      end if
      ! Convert to number density.
      densityChemicalICM= densityChemicalICM      &
           &             *massToDensityConversion &
           &             *density
      ! Compute number density of hydrogen (in cm⁻³).
      numberDensityHydrogen=+density                                    &
           &                *abundancesICM   %hydrogenMassFraction()    &
           &                *massSolar                                  &
           &                /massHydrogenAtom                           &
           &                /hecto                                  **3 &
           &                /megaParsec                             **3
      return
    end subroutine icmProperties
    
  end function icmCoolingPowerInBandExtract

  function icmCoolingPowerInBandName(self)
    !!{
    Return the name of the cooling power in band property.
    !!}
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    type (varying_string                            )                :: icmCoolingPowerInBandName
    class(nodePropertyExtractorICMCoolingPowerInBand), intent(inout) :: self

    icmCoolingPowerInBandName=var_str('icmCoolingPowerFractionInBand')//self%label
    return
  end function icmCoolingPowerInBandName

  function icmCoolingPowerInBandDescription(self)
    !!{
    Return a description of the cooling power in band property.
    !!}
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    type     (varying_string                            )                :: icmCoolingPowerInBandDescription
    class    (nodePropertyExtractorICMCoolingPowerInBand), intent(inout) :: self
    character(len=16                                    )                :: labelLow                        , labelHigh

    write (labelLow ,'(f16.12)') self%energyLow
    write (labelHigh,'(f16.12)') self%energyHigh
    icmCoolingPowerInBandDescription=var_str('The fraction of the cooling power of the hot halo emitted in the ')//trim(adjustl(labelLow))//"-"//trim(adjustl(labelHigh))//' keV band.'
    return
  end function icmCoolingPowerInBandDescription

  double precision function icmCoolingPowerInBandUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily icmCoolingPowerInBand} properties in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorICMCoolingPowerInBand), intent(inout) :: self
    !$GLC attributes unused :: self

    icmCoolingPowerInBandUnitsInSI=1.0d0
    return
  end function icmCoolingPowerInBandUnitsInSI


