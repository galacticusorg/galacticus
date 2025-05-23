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
Implements an intracluster medium X-ray luminosity-weighted temperature property extractor class.
!!}

  use :: Cooling_Functions      , only : coolingFunction    , coolingFunctionClass
  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorICMXRayTemperature">
   <description>An intracluster medium X-ray luminosity-weighted temperature property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorICMXRayTemperature
     !!{
     An ICM luminosity-weighted temperature property extractor class.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     class(coolingFunctionClass    ), pointer :: coolingFunction_     => null()
     class(cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
   contains
     final     ::                icmXRayTemperatureDestructor
     procedure :: extract     => icmXRayTemperatureExtract
     procedure :: name        => icmXRayTemperatureName
     procedure :: description => icmXRayTemperatureDescription
     procedure :: unitsInSI   => icmXRayTemperatureUnitsInSI
  end type nodePropertyExtractorICMXRayTemperature

  interface nodePropertyExtractorICMXRayTemperature
     !!{
     Constructors for the \refClass{nodePropertyExtractorICMXRayTemperature} output analysis class.
     !!}
     module procedure icmXRayTemperatureConstructorParameters
     module procedure icmXRayTemperatureConstructorInternal
  end interface nodePropertyExtractorICMXRayTemperature

contains

  function icmXRayTemperatureConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorICMXRayTemperature} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorICMXRayTemperature)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    class(coolingFunctionClass                   ), pointer       :: coolingFunction_
    class(cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="coolingFunction"     name="coolingFunction_"     source="parameters"/>
    !!]
    self=nodePropertyExtractorICMXRayTemperature(cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="coolingFunction_"    />
    !!]
    return
  end function icmXRayTemperatureConstructorParameters

  function icmXRayTemperatureConstructorInternal(cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorICMXRayTemperature} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorICMXRayTemperature)                        :: self
    class(cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class(darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    class(coolingFunctionClass                   ), intent(in   ), target :: coolingFunction_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *darkMatterHaloScale_, *coolingFunction_"/>
    !!]

    return
  end function icmXRayTemperatureConstructorInternal

  subroutine icmXRayTemperatureDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorICMXRayTemperature} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorICMXRayTemperature), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%coolingFunction_"    />
    !!]
    return
  end subroutine icmXRayTemperatureDestructor

  double precision function icmXRayTemperatureExtract(self,node,instance)
    !!{
    Implement an ICM X-ray properties extractor.
    !!}
    use :: Galacticus_Nodes            , only : nodeComponentHotHalo                   , treeNode
    use :: Numerical_Constants_Physical, only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Numerical_Constants_Units   , only : electronVolt
    use :: Numerical_Integration       , only : integrator
    use :: Radiation_Fields            , only : radiationFieldCosmicMicrowaveBackground
    use :: Mass_Distributions          , only : massDistributionClass                  , kinematicsDistributionClass
    use :: Galactic_Structure_Options  , only : componentTypeHotHalo                   , massTypeGaseous
    implicit none
    class           (nodePropertyExtractorICMXRayTemperature), intent(inout), target   :: self
    type            (treeNode                               ), intent(inout), target   :: node
    type            (multiCounter                           ), intent(inout), optional :: instance
    type            (radiationFieldCosmicMicrowaveBackground), pointer                 :: radiation_
    type            (integrator                             )                          :: integratorLuminosity   , integratorTemperature
    class           (massDistributionClass                  ), pointer                 :: massDistribution_
    class           (kinematicsDistributionClass            ), pointer                 :: kinematicsDistribution_
    double precision                                                                   :: luminosity             , temperature
    !$GLC attributes unused :: self, instance

    ! Initialize radiation field.
    allocate(radiation_)
    !![
    <referenceConstruct object="radiation_" constructor="radiationFieldCosmicMicrowaveBackground(self%cosmologyFunctions_)"/>
    !!]
    ! Get the mass distribution.
    massDistribution_       => node             %massDistribution      (componentTypeHotHalo,massTypeGaseous)
    kinematicsDistribution_ => massDistribution_%kinematicsDistribution(                                    )
    ! Compute luminosity and temperature.
    integratorLuminosity =integrator                     (integrandLuminosityXray ,toleranceRelative                           =1.0d-3)
    integratorTemperature=integrator                     (integrandTemperatureXray,toleranceRelative                           =1.0d-3)
    luminosity           =integratorLuminosity %integrate(0.0d0                   ,self%darkMatterHaloScale_%radiusVirial(node)       )
    temperature          =integratorTemperature%integrate(0.0d0                   ,self%darkMatterHaloScale_%radiusVirial(node)       )
    if (luminosity > 0.0d0) then
       temperature=+temperature        &
            &      /luminosity         &
            &      *boltzmannsConstant &
            &      /kilo               &
            &      /electronVolt
    else
       luminosity =+0.0d0
       temperature=+0.0d0
    end if
    !![
    <objectDestructor name="radiation_"             />
    <objectDestructor name="massDistribution_"      />
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    icmXRayTemperatureExtract=temperature
    return

  contains

    double precision function integrandLuminosityXray(radius)
      !!{
      Integrand function used for computing ICM X-ray luminosities.
      !!}
      use :: Abundances_Structure             , only : abundances
      use :: Chemical_Abundances_Structure    , only : chemicalAbundances
      use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Fraction_Conversion
      use :: Numerical_Constants_Astronomical , only : massSolar                            , megaParsec
      use :: Numerical_Constants_Atomic       , only : massHydrogenAtom
      use :: Numerical_Constants_Math         , only : Pi
      use :: Numerical_Constants_Prefixes     , only : centi                                , hecto
      use :: Coordinates                      , only : coordinateSpherical                  , assignment(=)
      implicit none
      double precision                      , intent(in   ) :: radius
      class           (nodeComponentHotHalo), pointer       :: hotHalo
      type            (coordinateSpherical )                :: coordinates
      double precision                                      :: density                , temperature        , &
           &                                                   numberDensityHydrogen  , massICM            , &
           &                                                   massToDensityConversion
      type            (abundances          )                :: abundancesICM
      type            (chemicalAbundances  )                :: massChemicalICM        , fractionChemicalICM

      ! Set the coordinates.
      coordinates             =  [radius,0.0d0,0.0d0]
      ! Get the density of the ICM.
      density                 =  massDistribution_      %density    (coordinates)
      ! Get the temperature of the ICM.
      temperature             =  kinematicsDistribution_%temperature(coordinates)
      ! Get abundances and chemistry of the ICM.
      hotHalo         => node   %hotHalo   ()
      massICM         =  hotHalo%mass      ()
      abundancesICM   =  hotHalo%abundances()
      massChemicalICM =  hotHalo%chemicals ()
      call abundancesICM  %massToMassFraction(            massICM)
      call massChemicalICM%massToNumber      (fractionChemicalICM)
      ! Compute factor converting mass of chemicals in (M☉) to number density in cm⁻³ per total mass density.
      if (hotHalo%mass() > 0.0d0) then
         massToDensityConversion=Chemicals_Mass_To_Fraction_Conversion(hotHalo%mass())
      else
         massToDensityConversion=0.0d0
      end if
      ! Convert to number density per unit total mass density.
      fractionChemicalICM=fractionChemicalICM*massToDensityConversion
      ! Compute number density of hydrogen (in cm⁻³).
      numberDensityHydrogen  =+density                                    &
           &                  *abundancesICM   %hydrogenMassFraction()    &
           &                  *massSolar                                  &
           &                  /massHydrogenAtom                           &
           &                  /hecto                                  **3 &
           &                  /megaParsec                             **3
      ! Evaluate the integrand.
      integrandLuminosityXray=+4.0d0                                                                                                                              &
           &                  *Pi                                                                                                                                 &
           &                  *radius**2                                                                                                                          &
           &                  *self%coolingFunction_%coolingFunction(node,numberDensityHydrogen,temperature,abundancesICM,fractionChemicalICM*density,radiation_) &
           &                  *(                                                                                                                                  &
           &                    +megaParsec                                                                                                                       &
           &                    /centi                                                                                                                            &
           &                   )**3
      return
    end function integrandLuminosityXray

    double precision function integrandTemperatureXray(radius)
      !!{
      Integrand function used for computing ICM X-ray luminosity-weighted temperatures.
      !!}
      use :: Coordinates, only : coordinateSpherical, assignment(=)
      implicit none
      double precision                     , intent(in   ) :: radius
      type            (coordinateSpherical)                :: coordinates

      coordinates             =[radius,0.0d0,0.0d0]
      integrandTemperatureXray=+integrandLuminosityXray            (     radius) &
           &                   *kinematicsDistribution_%temperature(coordinates)
      return
    end function integrandTemperatureXray

  end function icmXRayTemperatureExtract

  function icmXRayTemperatureName(self)
    !!{
    Return the names of the {\normalfont \ttfamily icmXRayTemperature} properties.
    !!}
    implicit none
    type (varying_string                         )                :: icmXRayTemperatureName
    class(nodePropertyExtractorICMXRayTemperature), intent(inout) :: self
    !$GLC attributes unused :: self

    icmXRayTemperatureName=var_str('icmXrayTemperature')
    return
  end function icmXRayTemperatureName

  function icmXRayTemperatureDescription(self)
    !!{
    Return descriptions of the {\normalfont \ttfamily icmXRayTemperature} properties.
    !!}
    implicit none
    type (varying_string                         )                :: icmXRayTemperatureDescription
    class(nodePropertyExtractorICMXRayTemperature), intent(inout) :: self
    !$GLC attributes unused :: self

    icmXRayTemperatureDescription=var_str('X-ray luminosity-weighted temperature of the ICM [keV]')
    return
  end function icmXRayTemperatureDescription

  double precision function icmXRayTemperatureUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily icmXRayTemperature} properties in the SI system.
    !!}
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Numerical_Constants_Units   , only : electronVolt
    implicit none
    class(nodePropertyExtractorICMXRayTemperature), intent(inout) :: self
    !$GLC attributes unused :: self

    icmXRayTemperatureUnitsInSI=kilo*electronVolt
    return
  end function icmXRayTemperatureUnitsInSI

