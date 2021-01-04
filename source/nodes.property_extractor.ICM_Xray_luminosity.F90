!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements an intracluster medium X-ray luminosity property extractor class.

  use :: Cooling_Functions            , only : coolingFunction          , coolingFunctionClass
  use :: Cosmology_Functions          , only : cosmologyFunctions       , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales      , only : darkMatterHaloScale      , darkMatterHaloScaleClass
  use :: Hot_Halo_Mass_Distributions  , only : hotHaloMassDistribution  , hotHaloMassDistributionClass
  use :: Hot_Halo_Temperature_Profiles, only : hotHaloTemperatureProfile, hotHaloTemperatureProfileClass

  !# <nodePropertyExtractor name="nodePropertyExtractorICMXRayLuminosity">
  !#  <description>An intracluster medium X-ray luminosity property extractor class.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorICMXRayLuminosity
     !% A icmXRayLuminosity property extractor class.
     private
     class(darkMatterHaloScaleClass      ), pointer :: darkMatterHaloScale_       => null()
     class(hotHaloMassDistributionClass  ), pointer :: hotHaloMassDistribution_   => null()
     class(hotHaloTemperatureProfileClass), pointer :: hotHaloTemperatureProfile_ => null()
     class(coolingFunctionClass          ), pointer :: coolingFunction_           => null()
     class(cosmologyFunctionsClass       ), pointer :: cosmologyFunctions_        => null()
   contains
     final     ::                 icmXRayLuminosityDestructor
     procedure :: elementCount => icmXRayLuminosityElementCount
     procedure :: extract      => icmXRayLuminosityExtract
     procedure :: names        => icmXRayLuminosityNames
     procedure :: descriptions => icmXRayLuminosityDescriptions
     procedure :: unitsInSI    => icmXRayLuminosityUnitsInSI
     procedure :: type         => icmXRayLuminosityType
  end type nodePropertyExtractorICMXRayLuminosity

  interface nodePropertyExtractorICMXRayLuminosity
     !% Constructors for the ``icmXRayLuminosity'' output analysis class.
     module procedure icmXRayLuminosityConstructorParameters
     module procedure icmXRayLuminosityConstructorInternal
  end interface nodePropertyExtractorICMXRayLuminosity

contains

  function icmXRayLuminosityConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily icmXRayLuminosity} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorICMXRayLuminosity)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    class(hotHaloMassDistributionClass          ), pointer       :: hotHaloMassDistribution_
    class(hotHaloTemperatureProfileClass        ), pointer       :: hotHaloTemperatureProfile_
    class(coolingFunctionClass                  ), pointer       :: coolingFunction_
    class(cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_

    !# <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"       name="darkMatterHaloScale_"       source="parameters"/>
    !# <objectBuilder class="hotHaloMassDistribution"   name="hotHaloMassDistribution_"   source="parameters"/>
    !# <objectBuilder class="hotHaloTemperatureProfile" name="hotHaloTemperatureProfile_" source="parameters"/>
    !# <objectBuilder class="coolingFunction"           name="coolingFunction_"           source="parameters"/>
    self=nodePropertyExtractorICMXRayLuminosity(cosmologyFunctions_,darkMatterHaloScale_,hotHaloMassDistribution_,hotHaloTemperatureProfile_,coolingFunction_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"       />
    !# <objectDestructor name="darkMatterHaloScale_"      />
    !# <objectDestructor name="hotHaloMassDistribution_"  />
    !# <objectDestructor name="hotHaloTemperatureProfile_"/>
    !# <objectDestructor name="coolingFunction_"          />
    return
  end function icmXRayLuminosityConstructorParameters

  function icmXRayLuminosityConstructorInternal(cosmologyFunctions_,darkMatterHaloScale_,hotHaloMassDistribution_,hotHaloTemperatureProfile_,coolingFunction_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily icmXRayLuminosity} property extractor class.
    implicit none
    type (nodePropertyExtractorICMXRayLuminosity)                        :: self
    class(cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class(darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    class(hotHaloMassDistributionClass          ), intent(in   ), target :: hotHaloMassDistribution_
    class(hotHaloTemperatureProfileClass        ), intent(in   ), target :: hotHaloTemperatureProfile_
    class(coolingFunctionClass                  ), intent(in   ), target :: coolingFunction_
    !# <constructorAssign variables="*cosmologyFunctions_, *darkMatterHaloScale_, *hotHaloMassDistribution_, *hotHaloTemperatureProfile_, *coolingFunction_"/>

    return
  end function icmXRayLuminosityConstructorInternal

  subroutine icmXRayLuminosityDestructor(self)
    !% Destructor for the {\normalfont \ttfamily icmXRayLuminosity} property extractor class.
    implicit none
    type(nodePropertyExtractorICMXRayLuminosity), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"       />
    !# <objectDestructor name="self%darkMatterHaloScale_"      />
    !# <objectDestructor name="self%hotHaloMassDistribution_"  />
    !# <objectDestructor name="self%hotHaloTemperatureProfile_"/>
    !# <objectDestructor name="self%coolingFunction_"          />
    return
  end subroutine icmXRayLuminosityDestructor

  integer function icmXRayLuminosityElementCount(self,time)
    !% Return the number of elements in the lightconeple property extractors.
    implicit none
    class           (nodePropertyExtractorICMXRayLuminosity), intent(inout) :: self
    double precision                                        , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    icmXRayLuminosityElementCount=2
    return
  end function icmXRayLuminosityElementCount

  function icmXRayLuminosityExtract(self,node,time,instance)
    !% Implement an ICM X-ray properties extractor.
    use :: Galacticus_Nodes            , only : nodeComponentHotHalo                   , treeNode
    use :: Numerical_Constants_Physical, only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Numerical_Constants_Units   , only : electronVolt
    use :: Numerical_Integration       , only : integrator
    use :: Radiation_Fields            , only : radiationFieldCosmicMicrowaveBackground
    implicit none
    double precision                                         , dimension(:) , allocatable :: icmXRayLuminosityExtract
    class           (nodePropertyExtractorICMXRayLuminosity ), intent(inout), target      :: self
    type            (treeNode                               ), intent(inout), target      :: node
    double precision                                         , intent(in   )              :: time
    type            (multiCounter                           ), intent(inout), optional    :: instance
    type            (radiationFieldCosmicMicrowaveBackground), pointer                    :: radiation_
    type            (integrator                             )                             :: integratorLuminosity    , integratorTemperature
    double precision                                                                      :: luminosity              , temperature
    !$GLC attributes unused :: self, time, instance

    allocate(icmXRayLuminosityExtract(2))
    ! Initialize radiation field.
    allocate(radiation_)
    !# <referenceConstruct object="radiation_" constructor="radiationFieldCosmicMicrowaveBackground(self%cosmologyFunctions_)"/>
    ! Compute luminosity and temperature.
    integratorLuminosity =integrator                     (integrandLuminosityXray ,toleranceRelative                           =1.0d-3)
    integratorTemperature=integrator                     (integrandTemperatureXray,toleranceRelative                           =1.0d-3)
    luminosity           =integratorLuminosity %integrate(0.0d0                   ,self%darkMatterHaloScale_%virialRadius(node)       )
    temperature          =integratorTemperature%integrate(0.0d0                   ,self%darkMatterHaloScale_%virialRadius(node)       )
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
    !# <objectDestructor name="radiation_"/>
    icmXRayLuminosityExtract=[luminosity,temperature]
    return

  contains

    double precision function integrandLuminosityXray(radius)
      !% Integrand function used for computing ICM X-ray luminosities.
      use :: Abundances_Structure             , only : abundances
      use :: Chemical_Abundances_Structure    , only : chemicalAbundances
      use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
      use :: Numerical_Constants_Astronomical , only : massSolar                           , megaParsec
      use :: Numerical_Constants_Atomic       , only : massHydrogenAtom
      use :: Numerical_Constants_Math         , only : Pi
      use :: Numerical_Constants_Prefixes     , only : centi                               , hecto
      implicit none
      double precision                      , intent(in   ) :: radius
      class           (nodeComponentHotHalo), pointer       :: hotHalo
      double precision                                      :: density                , temperature       , &
           &                                                   numberDensityHydrogen  , massICM           , &
           &                                                   massToDensityConversion
      type            (abundances          )                :: abundancesICM
      type            (chemicalAbundances  )                :: massChemicalICM        , densityChemicalICM

      ! Get the density of the ICM.
      density    =self%hotHaloMassDistribution_  %density    (node,radius)
      ! Get the temperature of the ICM.
      temperature=self%hotHaloTemperatureProfile_%temperature(node,radius)
      ! Get abundances and chemistry of the ICM.
      hotHalo         => node   %hotHalo   ()
      massICM         =  hotHalo%mass      ()
      abundancesICM   =  hotHalo%abundances()
      massChemicalICM =  hotHalo%chemicals ()
      call abundancesICM  %massToMassFraction(           massICM)
      call massChemicalICM%massToNumber      (densityChemicalICM)
      ! Compute factor converting mass of chemicals in (M☉/Mₐₘᵤ) to number density in cm⁻³.
      if (hotHalo%outerRadius() > 0.0d0) then
         massToDensityConversion=Chemicals_Mass_To_Density_Conversion(hotHalo%outerRadius())
      else
         massToDensityConversion=0.0d0
      end if
      ! Convert to number density.
      densityChemicalICM=densityChemicalICM*massToDensityConversion
      ! Compute number density of hydrogen (in cm⁻³).
      numberDensityHydrogen  =+density                                    &
           &                  *abundancesICM   %hydrogenMassFraction()    &
           &                  *massSolar                                  &
           &                  /massHydrogenAtom                           &
           &                  /hecto                                  **3 &
           &                  /megaParsec                             **3
      ! Evaluate the integrand.
      integrandLuminosityXray=+4.0d0                                                                                                                &
           &                  *Pi                                                                                                                   &
           &                  *radius**2                                                                                                            &
           &                  *self%coolingFunction_%coolingFunction(numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM,radiation_) &
           &                  *(                                                                                                                    &
           &                    +megaParsec                                                                                                         &
           &                    /centi                                                                                                              &
           &                   )**3
      return
    end function integrandLuminosityXray

    double precision function integrandTemperatureXray(radius)
      !% Integrand function used for computing ICM X-ray luminosity-weighted temperatures.
      implicit none
      double precision, intent(in   ) :: radius

      integrandTemperatureXray=+integrandLuminosityXray                     (     radius) &
           &                    *self%hotHaloTemperatureProfile_%temperature(node,radius)
      return
    end function integrandTemperatureXray

  end function icmXRayLuminosityExtract

  function icmXRayLuminosityNames(self,time)
    !% Return the names of the {\normalfont \ttfamily icmXRayLuminosity} properties.
    implicit none
    type            (varying_string                        ), dimension(:) , allocatable :: icmXRayLuminosityNames
    class           (nodePropertyExtractorICMXRayLuminosity), intent(inout)              :: self
    double precision                                        , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(icmXRayLuminosityNames(2))
    icmXRayLuminosityNames=[var_str('icmXrayLuminosity'),var_str('icmXrayTemperature')]
    return
  end function icmXRayLuminosityNames

  function icmXRayLuminosityDescriptions(self,time)
    !% Return descriptions of the {\normalfont \ttfamily icmXRayLuminosity} properties.
    implicit none
    type            (varying_string                        ), dimension(:) , allocatable :: icmXRayLuminosityDescriptions
    class           (nodePropertyExtractorICMXRayLuminosity), intent(inout)              :: self
    double precision                                        , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(icmXRayLuminosityDescriptions(2))
    icmXRayLuminosityDescriptions=[var_str('X-ray luminosity of the ICM [ergs/s]'),var_str('X-ray luminosity-weighted temperature of the ICM [keV]')]
    return
  end function icmXRayLuminosityDescriptions

  function icmXRayLuminosityUnitsInSI(self,time)
    !% Return the units of the {\normalfont \ttfamily icmXRayLuminosity} properties in the SI system.
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Numerical_Constants_Units   , only : electronVolt, ergs
    implicit none
    double precision                                        , allocatable  , dimension(:) :: icmXRayLuminosityUnitsInSI
    class           (nodePropertyExtractorICMXRayLuminosity), intent(inout)               :: self
    double precision                                        , intent(in   )               :: time
    !$GLC attributes unused :: self, time

    allocate(icmXRayLuminosityUnitsInSI(2))
    icmXRayLuminosityUnitsInSI=[ergs,kilo*electronVolt]
    return
  end function icmXRayLuminosityUnitsInSI

  integer function icmXRayLuminosityType(self)
    !% Return the type of the {\normalfont \ttfamily icmXRayLuminosity} properties.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorICMXRayLuminosity), intent(inout) :: self
    !$GLC attributes unused :: self

    icmXRayLuminosityType=outputAnalysisPropertyTypeLinear
    return
  end function icmXRayLuminosityType

