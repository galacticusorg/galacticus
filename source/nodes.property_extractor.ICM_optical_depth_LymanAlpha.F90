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
  
  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorICMOpticalDepthLymanAlpha">
    <description>An intracluster medium cooling power in band property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorICMOpticalDepthLymanAlpha
     !!{
     A property extractor class which extracts the fraction of the ICM cooling power due to emission in a given energy band.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_   => null()
     class  (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_    => null()
     integer                                    :: speciesHydrogenNeutral          , speciesHydrogenIonized
   contains
     final     ::                icmOpticalDepthLymanAlphaDestructor
     procedure :: extract     => icmOpticalDepthLymanAlphaExtract
     procedure :: name        => icmOpticalDepthLymanAlphaName
     procedure :: description => icmOpticalDepthLymanAlphaDescription
     procedure :: unitsInSI   => icmOpticalDepthLymanAlphaUnitsInSI
  end type nodePropertyExtractorICMOpticalDepthLymanAlpha

  interface nodePropertyExtractorICMOpticalDepthLymanAlpha
     !!{
     Constructors for the {\normalfont \ttfamily icmOpticalDepthLymanAlpha} output analysis class.
     !!}
     module procedure icmOpticalDepthLymanAlphaConstructorParameters
     module procedure icmOpticalDepthLymanAlphaConstructorInternal
  end interface nodePropertyExtractorICMOpticalDepthLymanAlpha

contains

  function icmOpticalDepthLymanAlphaConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily icmOpticalDepthLymanAlpha} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorICMOpticalDepthLymanAlpha)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass                      ), pointer       :: darkMatterHaloScale_
    class(cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorICMOpticalDepthLymanAlpha(cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function icmOpticalDepthLymanAlphaConstructorParameters

  function icmOpticalDepthLymanAlphaConstructorInternal(cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily icmOpticalDepthLymanAlpha} property extractor class.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    implicit none
    type (nodePropertyExtractorICMOpticalDepthLymanAlpha)                        :: self
    class(cosmologyFunctionsClass                       ), intent(in   ), target :: cosmologyFunctions_
    class(darkMatterHaloScaleClass                      ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]

    self%speciesHydrogenNeutral=Chemicals_Index("AtomicHydrogen"      )
    self%speciesHydrogenIonized=Chemicals_Index("AtomicHydrogenCation")
    return
  end function icmOpticalDepthLymanAlphaConstructorInternal

  subroutine icmOpticalDepthLymanAlphaDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily icmOpticalDepthLymanAlpha} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorICMOpticalDepthLymanAlpha), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine icmOpticalDepthLymanAlphaDestructor

  double precision function icmOpticalDepthLymanAlphaExtract(self,node,instance)
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
    class(nodePropertyExtractorICMOpticalDepthLymanAlpha), intent(inout), target   :: self
    type (treeNode                                      ), intent(inout), target   :: node
    type (multiCounter                                  ), intent(inout), optional :: instance
    type (radiationFieldCosmicMicrowaveBackground       ), pointer                 :: radiation_
    class(massDistributionClass                         ), pointer                 :: massDistribution_
    class(kinematicsDistributionClass                   ), pointer                 :: kinematicsDistribution_
    type (integrator                                    )                          :: integrator_
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
    integrator_                     =integrator           (integrandOpticalDepth,toleranceRelative                           =1.0d-3)
    icmOpticalDepthLymanAlphaExtract=integrator_%integrate(0.0d0                ,self%darkMatterHaloScale_%radiusVirial(node)       )    
    !![
    <objectDestructor name="radiation_"             />
    <objectDestructor name="massDistribution_"      />
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return

  contains
    
    double precision function integrandOpticalDepth(radius)
      !!{
      Integrand function used for computing ICM optical depth to Lyman-$\alpha$.
      !!}
      use :: Abundances_Structure            , only : abundances
      use :: Chemical_Abundances_Structure   , only : chemicalAbundances
      use :: Numerical_Constants_Prefixes    , only : centi
      use :: Numerical_Constants_Astronomical, only : megaParsec
      implicit none
      double precision                      , intent(in   ) :: radius
      type            (abundances          )                :: abundancesICM
      type            (chemicalAbundances  )                :: densityChemicalICM
      double precision                                      :: numberDensityHydrogen  , temperature, &
           &                                                   fractionHydrogenNeutral

      call icmProperties(radius,numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM)
      if (numberDensityHydrogen > 0.0d0) then
         fractionHydrogenNeutral=+  densityChemicalICM%abundance(self%speciesHydrogenNeutral) &
              &                  /(                                                           &
              &                    +densityChemicalICM%abundance(self%speciesHydrogenNeutral) &
              &                    +densityChemicalICM%abundance(self%speciesHydrogenIonized) &
              &                   )
         integrandOpticalDepth  =+crossSectionLymanAlphaLineCenter(temperature) &
              &                  *numberDensityHydrogen                         &
              &                  *fractionHydrogenNeutral                       &
              &                  *megaParsec                                    &
              &                  /centi
      else
         integrandOpticalDepth  =+0.0d0
      end if
      return
    end function integrandOpticalDepth

    double precision function crossSectionLymanAlphaLineCenter(temperature)
      !!{
      The cross-section at the center of the Lyman-$\alpha$ line at finite temperature \cite[][eqn.~54]{dijkstra_saas-fee_2017}.
      !!}
      implicit none
      double precision, intent(in   ) :: temperature

      crossSectionLymanAlphaLineCenter=+5.9d-14           & ! cm²
           &                           /sqrt(             &
           &                                 +temperature &
           &                                 /1.0d4       &
           &                                )
      return
    end function crossSectionLymanAlphaLineCenter
    
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
    
  end function icmOpticalDepthLymanAlphaExtract

  function icmOpticalDepthLymanAlphaName(self)
    !!{
    Return the name of the cooling power in band property.
    !!}
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    type (varying_string                                )                :: icmOpticalDepthLymanAlphaName
    class(nodePropertyExtractorICMOpticalDepthLymanAlpha), intent(inout) :: self

    icmOpticalDepthLymanAlphaName=var_str('icmOpticalDepthLymanAlpha')
    return
  end function icmOpticalDepthLymanAlphaName

  function icmOpticalDepthLymanAlphaDescription(self)
    !!{
    Return a description of the cooling power in band property.
    !!}
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    type (varying_string                                )                :: icmOpticalDepthLymanAlphaDescription
    class(nodePropertyExtractorICMOpticalDepthLymanAlpha), intent(inout) :: self

    icmOpticalDepthLymanAlphaDescription=var_str('The optical depth through the ICM at the center of the Lyman-α line.')
    return
  end function icmOpticalDepthLymanAlphaDescription

  double precision function icmOpticalDepthLymanAlphaUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily icmOpticalDepthLymanAlpha} properties in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorICMOpticalDepthLymanAlpha), intent(inout) :: self
    !$GLC attributes unused :: self

    icmOpticalDepthLymanAlphaUnitsInSI=1.0d0
    return
  end function icmOpticalDepthLymanAlphaUnitsInSI


