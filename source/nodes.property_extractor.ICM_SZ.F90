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
Implements an intracluster medium Sunyaev-Zeldovich Compton-y parameter property extractor class.
!!}
  use :: Chemical_States        , only : chemicalState      , chemicalStateClass
  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass , enumerationDensityCosmologicalType
  use :: Cosmology_Parameters   , only : cosmologyParameters, cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorICMSZ">
   <description>
    An intracluster medium Sunyaev-Zeldovich Compton-$y$ parameter property extractor class. Specifically, the quantity
    extracted is
    \begin{equation}
     Y = {\sigma_\mathrm{T} \over \mathrm{m}_\mathrm{e} \mathrm{c}^2} \int_0^{R_\mathrm{outer}} n_\mathrm{e}(R) \mathrm{k}_\mathrm{B} T(r) {4 \pi R^2 \mathrm{d} R \over D_\mathrm{A}^2},
    \end{equation}
    where $D_\mathrm{A}$ is the angular diameter distance to the halo, and the result is expressed in units of square
    arcminutes. The angular diameter distance is, by default, computed from the epoch of the halo. Alternatively, a fixed
    angular diameter distance can be specified via the {\normalfont \ttfamily [distanceAngular]} parameter. The outer radius,
    $R_\mathrm{out}$, is either the halo virial radius (by default), or the radius enclosing the density contrast specified by
    the optional {\normalfont \ttfamily [densityContrast]} parameter. This density contrast is relative to either {\normalfont
    \ttfamily mean} or {\normalfont \ttfamily critical} density as specified by the {\normalfont \ttfamily
    densityContrastRelativeTo} parameter.
   </description>
   <deepCopy>
    <functionClass variables="densityContrastExtractor_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="densityContrastExtractor_"/>
   </stateStorable>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorICMSZ
     !!{
     An intracluster medium Sunyaev-Zeldovich Compton-y parameter property extractor class.
     !!}
     private
     class           (cosmologyParametersClass             ), pointer :: cosmologyParameters_       => null()
     class           (cosmologyFunctionsClass              ), pointer :: cosmologyFunctions_        => null()
     class           (darkMatterHaloScaleClass             ), pointer :: darkMatterHaloScale_       => null()
     class           (chemicalStateClass                   ), pointer :: chemicalState_             => null()
     type            (nodePropertyExtractorDensityContrasts), pointer :: densityContrastExtractor_  => null()
     double precision                                                 :: densityContrast                     , distanceAngular
     logical                                                          :: useDensityContrast                  , useFixedDistance
     type            (enumerationDensityCosmologicalType   )          :: densityContrastRelativeTo
     type            (varying_string                       )          :: name_
   contains
     final     ::                icmSZDestructor
     procedure :: extract     => icmSZExtract
     procedure :: name        => icmSZName
     procedure :: description => icmSZDescription
     procedure :: unitsInSI   => icmSZUnitsInSI
  end type nodePropertyExtractorICMSZ

  interface nodePropertyExtractorICMSZ
     !!{
     Constructors for the \refClass{nodePropertyExtractorICMSZ} output analysis class.
     !!}
     module procedure icmSZConstructorParameters
     module procedure icmSZConstructorInternal
  end interface nodePropertyExtractorICMSZ

contains

  function icmSZConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorICMSZ} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters   , only : inputParameter                      , inputParameters
    use :: Cosmology_Functions, only : enumerationDensityCosmologicalEncode
    implicit none
    type            (nodePropertyExtractorICMSZ)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass  ), pointer       :: darkMatterHaloScale_
    class           (chemicalStateClass        ), pointer       :: chemicalState_
    double precision                                            :: densityContrast           , distanceAngular
    type            (varying_string            )                :: densityContrastRelativeTo

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="chemicalState"       name="chemicalState_"       source="parameters"/>
    !!]
    if (parameters%isPresent('densityContrast')) then
       !![
       <inputParameter>
         <name>densityContrast</name>
         <description>The density contrast within which to compute the Sunyaev-Zeldovich parameter.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>densityContrastRelativeTo</name>
         <description>The density ({\normalfont \ttfamily mean} or {\normalfont \ttfamily critical}) used in defining the density contrast.</description>
         <source>parameters</source>
         <defaultValue>var_str('mean')</defaultValue>
       </inputParameter>
       !!]
    end if
    if (parameters%isPresent('distanceAngular')) then
       !![
       <inputParameter>
         <name>distanceAngular</name>
         <description>The fixed angular diameter distance at which to compute the Sunyaev-Zeldovich parameter.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <conditionalCall>
     <call>self=nodePropertyExtractorICMSZ(cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,chemicalState_{conditions})</call>
     <argument name="densityContrast"           value="densityContrast"                                                                              parameterPresent="parameters"                                />
     <argument name="densityContrastRelativeTo" value="enumerationDensityCosmologicalEncode(char(densityContrastRelativeTo),includesPrefix=.false.)" parameterPresent="parameters" parameterName="densityContrast"/>
     <argument name="distanceAngular"           value="distanceAngular"                                                                              parameterPresent="parameters"                                />
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="chemicalState_"      />
    !!]
    return
  end function icmSZConstructorParameters

  function icmSZConstructorInternal(cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,chemicalState_,densityContrast,densityContrastRelativeTo,distanceAngular) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorICMSZ} property extractor class.
    !!}
    use :: Cosmology_Functions, only : densityCosmologicalMean, enumerationDensityCosmologicalDecode, enumerationDensityCosmologicalType
    use :: ISO_Varying_String , only : char
    use :: String_Handling    , only : String_Upper_Case_First
    implicit none
    type            (nodePropertyExtractorICMSZ        )                          :: self
    class           (cosmologyParametersClass          ), intent(in   ), target   :: cosmologyParameters_
    class           (cosmologyFunctionsClass           ), intent(in   ), target   :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass          ), intent(in   ), target   :: darkMatterHaloScale_
    class           (chemicalStateClass                ), intent(in   ), target   :: chemicalState_
    double precision                                    , intent(in   ), optional :: densityContrast           , distanceAngular
    type            (enumerationDensityCosmologicalType), intent(in   ), optional :: densityContrastRelativeTo
    character       (len=8                             )                          :: label
    !![
    <constructorAssign variables="densityContrast, densityContrastRelativeTo, distanceAngular, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *chemicalState_"/>
    !!]

    self%useDensityContrast=present(densityContrast)
    self%useFixedDistance  =present(distanceAngular)
    if (.not.present(densityContrastRelativeTo)) self%densityContrastRelativeTo=densityCosmologicalMean
    if (self%useDensityContrast) then
       allocate(self%densityContrastExtractor_)
       !![
       <referenceConstruct owner="self" isResult="yes" object="densityContrastExtractor_">
        <constructor>
         nodePropertyExtractorDensityContrasts(                                                     &amp;
          &amp;                                densityContrasts         =[densityContrast]        , &amp;
          &amp;                                darkMatterOnly           =.true.                   , &amp;
          &amp;                                densityContrastRelativeTo=densityContrastRelativeTo, &amp;
          &amp;                                cosmologyparameters_     =cosmologyParameters_     , &amp;
          &amp;                                cosmologyFunctions_      =cosmologyFunctions_      , &amp;
          &amp;                                darkMatterHaloScale_     =darkMatterHaloScale_       &amp;
          &amp;                               )
        </constructor>
       </referenceConstruct>
       !!]
       write (label,'(f7.2)') densityContrast
       self%name_="icmSZComptonYR"//trim(adjustl(label))//String_Upper_Case_First(char(enumerationDensityCosmologicalDecode(densityContrastRelativeTo,includePrefix=.false.)))
    else
       self%name_="icmSZComptonYVirial"
    end if
    return
  end function icmSZConstructorInternal

  subroutine icmSZDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorICMSZ} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorICMSZ), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%chemicalState_"      />
    !!]
    if (self%useDensityContrast) then
       !![
       <objectDestructor name="self%densityContrastExtractor_"/>
       !!]
    end if
    return
  end subroutine icmSZDestructor

  double precision function icmSZExtract(self,node,instance)
    !!{
    Implement a Sunyaev-Zeldovich effect property extractor.
    !!}
    use :: Numerical_Integration           , only : integrator
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentHotHalo                   , nodeComponentBasic
    use :: Numerical_Constants_Astronomical, only : degreesToRadians                       , arcminutesToDegrees
    use :: Numerical_Constants_Math        , only : Pi
    use :: Radiation_Fields                , only : radiationFieldCosmicMicrowaveBackground
    implicit none
    class           (nodePropertyExtractorICMSZ             ), intent(inout), target   :: self
    type            (treeNode                               ), intent(inout), target   :: node
    type            (multiCounter                           ), intent(inout), optional :: instance
    type            (radiationFieldCosmicMicrowaveBackground), pointer                 :: radiation_
    double precision                                         , dimension(2)            :: densityContrastProperties
    class           (nodeComponentBasic                     ), pointer                 :: basic
    type            (integrator                             )                          :: integrator_
    double precision                                                                   :: radiusOuter              , time, &
         &                                                                                distanceAngular
    
    ! Initialize radiation field.
    allocate(radiation_)
    !![
    <referenceConstruct object="radiation_" constructor="radiationFieldCosmicMicrowaveBackground(self%cosmologyFunctions_)"/>
    !!]
    ! Extract the time at which this node exists.
    basic => node %basic()
    time  =  basic%time ()
    ! Get upper limit for integral.
    if (self%useDensityContrast) then
       densityContrastProperties=reshape(self%densityContrastExtractor_%extract     (node,time,instance),[2])
       radiusOuter              =             densityContrastProperties             (1                 )
    else
       radiusOuter              =        self%darkMatterHaloScale_     %radiusVirial(node              )
    end if
    ! Find the angular diameter distance to use.
    if (self%useFixedDistance) then
       distanceAngular=self%                    distanceAngular
    else
       distanceAngular=self%cosmologyFunctions_%distanceAngular(time)
    end if
    if (distanceAngular <= 0.0d0) call Error_Report('non-positive angular diameter distance'//{introspection:location})
    ! Compute the integrated Compton-y parameter within this radius, divided by the angular diameter distance squared, and
    ! converted to units of arcminutes².
    integrator_ = integrator           (integrandComptonY,toleranceRelative=1.0d-3)
    icmSZExtract=+integrator_%integrate(0.0d0            ,radiusOuter             )    &
         &       /distanceAngular                                                  **2 &
         &       /degreesToRadians                                                 **2 &
         &       /arcminutesToDegrees                                              **2
    !![
    <objectDestructor name="radiation_"/>
    !!]
    return

  contains

    double precision function integrandComptonY(radius)
      !!{
      Integrand function used for computing ICM SZ properties.
      !!}
      use :: Abundances_Structure            , only : abundances
      use :: Numerical_Constants_Astronomical, only : massSolar            , megaParsec
      use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
      use :: Numerical_Constants_Physical    , only : boltzmannsConstant   , electronMass               , speedLight, thomsonCrossSection
      use :: Numerical_Constants_Prefixes    , only : centi                , hecto
      use :: Mass_Distributions              , only : massDistributionClass, kinematicsDistributionClass
      use :: Coordinates                     , only : coordinateSpherical  , assignment(=)
      use :: Galactic_Structure_Options      , only : componentTypeHotHalo , massTypeGaseous
      implicit none
      double precision                             , intent(in   ) :: radius
      class           (nodeComponentHotHalo       ), pointer       :: hotHalo
      class           (massDistributionClass      ), pointer       :: massDistribution_
      class           (kinematicsDistributionClass), pointer       :: kinematicsDistribution_
      type            (coordinateSpherical        )                :: coordinates
      double precision                                             :: density                , temperature, &
           &                                                          numberDensityHydrogen  , massICM
      type            (abundances                 )                :: abundancesICM

      ! Get the mass distribution.
      coordinates             =  [radius,0.0d0,0.0d0]
      massDistribution_       => node                   %massDistribution      (componentTypeHotHalo,massTypeGaseous)
      kinematicsDistribution_ => massDistribution_      %kinematicsDistribution(                                    )      
      ! Get the density of the ICM.
      density                 =  massDistribution_      %density               (coordinates                         )
      ! Get the temperature of the ICM.
      temperature             =  kinematicsDistribution_%temperature           (coordinates                         )
      !![
      <objectDestructor name="massDistribution_"      />
      <objectDestructor name="kinematicsDistribution_"/>
      !!]          
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
      ! Evaluate the integrand. This gives a result in units of Mpc² - we will divide by the angular diameter distance (in Mpc)
      ! squared later.
      integrandComptonY=+4.0d0                                                                                              &
           &            *Pi                                                                                                 &
           &            *radius                                                                                         **2 &
           &            *boltzmannsConstant                                                                                 &
           &            *thomsonCrossSection                                                                                &
           &            /electronMass                                                                                       &
           &            /speedLight                                                                                     **2 &
           &            *self%chemicalState_%electronDensity(numberDensityHydrogen,temperature,abundancesICM,radiation_)    &
           &            *temperature                                                                                        &
           &            *megaParsec                                                                                         &
           &            /centi                                                                                          **3
      return
    end function integrandComptonY

  end function icmSZExtract

  function icmSZName(self)
    !!{
    Return the name of the last isolated redshift property.
    !!}
    implicit none
    type (varying_string            )                :: icmSZName
    class(nodePropertyExtractorICMSZ), intent(inout) :: self

    icmSZName=self%name_
    return
  end function icmSZName

  function icmSZDescription(self)
    !!{
    Return a description of the intracluster medium Sunyaev-Zeldovich property.
    !!}
    implicit none
    type (varying_string            )                :: icmSZDescription
    class(nodePropertyExtractorICMSZ), intent(inout) :: self
    !$GLC attributes unused :: self

    icmSZDescription=var_str('Mean thermal Sunyaev-Zeldovich Compton y-parameter of the ICM within the virial radius.')
    return
  end function icmSZDescription

  double precision function icmSZUnitsInSI(self)
    !!{
    Return the units of the last isolated redshift property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorICMSZ), intent(inout) :: self
    !$GLC attributes unused :: self

    icmSZUnitsInSI=0.0d0
    return
  end function icmSZUnitsInSI


