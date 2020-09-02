!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
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

!% Contains a module which implements an intracluster medium Sunyaev-Zeldovich property extractor class.
  use :: Chemical_States              , only : chemicalState            , chemicalStateClass
  use :: Cosmology_Functions          , only : cosmologyFunctions       , cosmologyFunctionsClass
  use :: Cosmology_Parameters         , only : cosmologyParameters      , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales      , only : darkMatterHaloScale      , darkMatterHaloScaleClass
  use :: Hot_Halo_Mass_Distributions  , only : hotHaloMassDistribution  , hotHaloMassDistributionClass
  use :: Hot_Halo_Temperature_Profiles, only : hotHaloTemperatureProfile, hotHaloTemperatureProfileClass

  !# <nodePropertyExtractor name="nodePropertyExtractorICMSZ">
  !#  <description>An intracluster medium Sunyaev-Zeldovich property extractor class.</description>
  !#  <deepCopy>
  !#   <functionClass variables="densityContrastExtractor_"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="densityContrastExtractor_"/>
  !#  </stateStorable>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorICMSZ
     !% An intracluster medium Sunyaev-Zeldovich property extractor class.
     private
     class           (cosmologyParametersClass             ), pointer :: cosmologyParameters_       => null()
     class           (cosmologyFunctionsClass              ), pointer :: cosmologyFunctions_        => null()
     class           (darkMatterHaloScaleClass             ), pointer :: darkMatterHaloScale_       => null()
     class           (hotHaloMassDistributionClass         ), pointer :: hotHaloMassDistribution_   => null()
     class           (hotHaloTemperatureProfileClass       ), pointer :: hotHaloTemperatureProfile_ => null()
     class           (chemicalStateClass                   ), pointer :: chemicalState_             => null()
     type            (nodePropertyExtractorDensityContrasts), pointer :: densityContrastExtractor_  => null()
     double precision                                                 :: densityContrast
     logical                                                          :: useDensityContrast
     type            (varying_string                       )          :: name_
   contains
     final     ::                icmSZDestructor
     procedure :: extract     => icmSZExtract
     procedure :: name        => icmSZName
     procedure :: description => icmSZDescription
     procedure :: unitsInSI   => icmSZUnitsInSI
     procedure :: type        => icmSZType
  end type nodePropertyExtractorICMSZ

  interface nodePropertyExtractorICMSZ
     !% Constructors for the ``icmSZ'' output analysis class.
     module procedure icmSZConstructorParameters
     module procedure icmSZConstructorInternal
  end interface nodePropertyExtractorICMSZ

contains

  function icmSZConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily icmSZ} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorICMSZ    )                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class           (hotHaloMassDistributionClass  ), pointer       :: hotHaloMassDistribution_
    class           (hotHaloTemperatureProfileClass), pointer       :: hotHaloTemperatureProfile_
    class           (chemicalStateClass            ), pointer       :: chemicalState_
    double precision                                                :: densityContrast

    !# <objectBuilder class="cosmologyParameters"       name="cosmologyParameters_"       source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"       name="darkMatterHaloScale_"       source="parameters"/>
    !# <objectBuilder class="hotHaloMassDistribution"   name="hotHaloMassDistribution_"   source="parameters"/>
    !# <objectBuilder class="hotHaloTemperatureProfile" name="hotHaloTemperatureProfile_" source="parameters"/>
    !# <objectBuilder class="chemicalState"             name="chemicalState_"             source="parameters"/>
    if (parameters%isPresent('densityContrast')) then
       !# <inputParameter>
       !#   <name>densityContrast</name>
       !#   <description>The density contrast within which to compute the Sunyaev-Zeldovich parameter.</description>
       !#   <source>parameters</source>
       !# </inputParameter>
    end if
    !# <conditionalCall>
    !#  <call>self=nodePropertyExtractorICMSZ(cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,hotHaloMassDistribution_,hotHaloTemperatureProfile_,chemicalState_{conditions})</call>
    !#  <argument name="densityContrast" value="densityContrast" parameterPresent="parameters"/>
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"      />
    !# <objectDestructor name="cosmologyFunctions_"       />
    !# <objectDestructor name="darkMatterHaloScale_"      />
    !# <objectDestructor name="hotHaloMassDistribution_"  />
    !# <objectDestructor name="hotHaloTemperatureProfile_"/>
    !# <objectDestructor name="chemicalState_"            />
    return
  end function icmSZConstructorParameters

  function icmSZConstructorInternal(cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,hotHaloMassDistribution_,hotHaloTemperatureProfile_,chemicalState_,densityContrast) result(self)
    !% Internal constructor for the {\normalfont \ttfamily icmSZ} property extractor class.
    implicit none
    type            (nodePropertyExtractorICMSZ    )                          :: self
    class           (cosmologyParametersClass      ), intent(in   ), target   :: cosmologyParameters_
    class           (cosmologyFunctionsClass       ), intent(in   ), target   :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass      ), intent(in   ), target   :: darkMatterHaloScale_
    class           (hotHaloMassDistributionClass  ), intent(in   ), target   :: hotHaloMassDistribution_
    class           (hotHaloTemperatureProfileClass), intent(in   ), target   :: hotHaloTemperatureProfile_
    class           (chemicalStateClass            ), intent(in   ), target   :: chemicalState_
    double precision                                , intent(in   ), optional :: densityContrast
    character       (len=8                         )                          :: label
    !# <constructorAssign variables="densityContrast, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *hotHaloMassDistribution_, *hotHaloTemperatureProfile_, *chemicalState_"/>

    self%useDensityContrast=present(densityContrast)
    if (self%useDensityContrast) then
       allocate(self%densityContrastExtractor_)
       !# <referenceConstruct owner="self" isResult="yes" object="densityContrastExtractor_">
       !#  <constructor>
       !#   nodePropertyExtractorDensityContrasts(                                           &amp;
       !#    &amp;                                densityContrasts    =[densityContrast]   , &amp;
       !#    &amp;                                darkMatterOnly      =.true.              , &amp;
       !#    &amp;                                cosmologyparameters_=cosmologyParameters_, &amp;
       !#    &amp;                                cosmologyFunctions_ =cosmologyFunctions_ , &amp;
       !#    &amp;                                darkMatterHaloScale_=darkMatterHaloScale_  &amp;
       !#    &amp;                               )
       !#  </constructor>
       !# </referenceConstruct>
       write (label,'(f7.2)') densityContrast
       self%name_="icmSZComptonYMeanR"//trim(adjustl(label))
    else
       self%name_="icmSZComptonYMeanVirial"
    end if
    return
  end function icmSZConstructorInternal

  subroutine icmSZDestructor(self)
    !% Destructor for the {\normalfont \ttfamily icmSZ} property extractor class.
    implicit none
    type(nodePropertyExtractorICMSZ), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"      />
    !# <objectDestructor name="self%cosmologyFunctions_"       />
    !# <objectDestructor name="self%darkMatterHaloScale_"      />
    !# <objectDestructor name="self%hotHaloMassDistribution_"  />
    !# <objectDestructor name="self%hotHaloTemperatureProfile_"/>
    !# <objectDestructor name="self%chemicalState_"            />
    if (self%useDensityContrast) then
       !# <objectDestructor name="self%densityContrastExtractor_"/>
    end if
    return
  end subroutine icmSZDestructor

  double precision function icmSZExtract(self,node,instance)
    !% Implement a Sunyaev-Zeldovich effect property extractor.
    use Numerical_Integration   , only : integrator
    use Galacticus_Nodes        , only : nodeComponentHotHalo                   , nodeComponentBasic
    use Numerical_Constants_Math, only : Pi
    use Radiation_Fields        , only : radiationFieldCosmicMicrowaveBackground
    implicit none
    class           (nodePropertyExtractorICMSZ             ), intent(inout)           :: self
    type            (treeNode                               ), intent(inout), target   :: node
    type            (multiCounter                           ), intent(inout), optional :: instance
    type            (radiationFieldCosmicMicrowaveBackground), pointer                 :: radiation_
    double precision                                         , dimension(2)            :: densityContrastProperties
    class           (nodeComponentBasic                     ), pointer                 :: basic
    type            (integrator                             )                          :: integrator_
    double precision                                                                   :: radiusOuter              , time

    ! Initialize radiation field.
    allocate(radiation_)
    !# <referenceConstruct object="radiation_" constructor="radiationFieldCosmicMicrowaveBackground(self%cosmologyFunctions_)"/>
    ! Get upper limit for integral.
    if (self%useDensityContrast) then
       basic                     => node %basic                                 (                  )
       time                      =  basic%time                                  (                  )
       densityContrastProperties =  self %densityContrastExtractor_%extract     (node,time,instance)
       radiusOuter               =        densityContrastProperties             (1                 )
    else
       radiusOuter              =   self %darkMatterHaloScale_     %virialRadius(node              )
    end if
    ! Compute mean Compton-y parameter within this radius.
    integrator_ = integrator           (integrandComptionY,toleranceRelative=1.0d-3)
    icmSZExtract=+integrator_%integrate(0.0d0             ,radiusOuter             ) &
         &       /Pi                                                                 &
         &       /radiusOuter**2
    !# <objectDestructor name="radiation_"/>
    return

  contains

    double precision function integrandComptionY(radius)
      !% Integrand function used for computing ICM SZ properties.
      use :: Abundances_Structure            , only : abundances
      use :: Numerical_Constants_Astronomical, only : massSolar         , megaParsec
      use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
      use :: Numerical_Constants_Physical    , only : boltzmannsConstant, electronMass, speedLight, thomsonCrossSection
      use :: Numerical_Constants_Prefixes    , only : centi             , hecto
      implicit none
      double precision                      , intent(in   ) :: radius
      class           (nodeComponentHotHalo), pointer       :: hotHalo
      double precision                                      :: density              , temperature, &
           &                                                   numberDensityHydrogen, massICM
      type            (abundances          )                :: abundancesICM

      ! Get the density of the ICM.
      density    =self%hotHaloMassDistribution_  %density    (node,radius)
      ! Get the temperature of the ICM.
      temperature=self%hotHaloTemperatureProfile_%temperature(node,radius)
      ! Get abundances and chemistry of the ICM.
      hotHalo         => node   %hotHalo   ()
      massICM         =  hotHalo%mass      ()
      abundancesICM   =  hotHalo%abundances()
      call abundancesICM%massToMassFraction(massICM)
      ! Compute number density of hydrogen (in cm⁻³).
      numberDensityHydrogen  =+density                                    &
           &                  *abundancesICM   %hydrogenMassFraction()    &
           &                  *massSolar                                  &
           &                  /massHydrogenAtom                           &
           &                  /hecto                                  **3 &
           &                  /megaParsec                             **3
      ! Evaluate the integrand.
      integrandComptionY=+4.0d0                                                                                              &
           &             *Pi                                                                                                 &
           &             *radius                                                                                         **2 &
           &             *boltzmannsConstant                                                                                 &
           &             *thomsonCrossSection                                                                                &
           &             /electronMass                                                                                       &
           &             /speedLight                                                                                     **2 &
           &             *self%chemicalState_%electronDensity(numberDensityHydrogen,temperature,abundancesICM,radiation_)    &
           &             *temperature                                                                                        &
           &             *megaParsec                                                                                         &
           &             /centi                                                                                          **3
      return
    end function integrandComptionY

  end function icmSZExtract

  function icmSZName(self)
    !% Return the name of the last isolated redshift property.
    implicit none
    type (varying_string            )                :: icmSZName
    class(nodePropertyExtractorICMSZ), intent(inout) :: self

    icmSZName=self%name_
    return
  end function icmSZName

  function icmSZDescription(self)
    !% Return a description of the intracluster medium Sunyaev-Zeldovich property.
    implicit none
    type (varying_string            )                :: icmSZDescription
    class(nodePropertyExtractorICMSZ), intent(inout) :: self
    !$GLC attributes unused :: self

    icmSZDescription=var_str('Mean thermal Sunyaev-Zeldovich Compton y-parameter of the ICM within the virial radius.')
    return
  end function icmSZDescription

  double precision function icmSZUnitsInSI(self)
    !% Return the units of the last isolated redshift property in the SI system.
    implicit none
    class(nodePropertyExtractorICMSZ), intent(inout) :: self
    !$GLC attributes unused :: self

    icmSZUnitsInSI=0.0d0
    return
  end function icmSZUnitsInSI

  integer function icmSZType(self)
    !% Return the type of the last isolated redshift property.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorICMSZ), intent(inout) :: self
    !$GLC attributes unused :: self

    icmSZType=outputAnalysisPropertyTypeLinear
    return
  end function icmSZType

