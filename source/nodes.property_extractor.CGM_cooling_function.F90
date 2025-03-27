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
  Implements a property extractor class for the CGM cooling function at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleClass
  use :: Cosmology_Functions                 , only : cosmologyFunctionsClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier
  use :: Cooling_Functions                   , only : coolingFunctionClass
  use :: Radiation_Fields                    , only : radiationFieldCosmicMicrowaveBackground
  !![
  <nodePropertyExtractor name="nodePropertyExtractorCGMCoolingFunction">
   <description>A property extractor class for the CGM cooling function at a set of radii.</description>
   <deepCopy>
    <functionClass variables="radiation"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="radiation"/>
   </stateStorable>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorCGMCoolingFunction
     !!{
     A property extractor class for the CGM cooling function at a set of radii.
     !!}
     private
     class  (cosmologyFunctionsClass                ), pointer                   :: cosmologyFunctions_           => null()
     class  (darkMatterHaloScaleClass               ), pointer                   :: darkMatterHaloScale_          => null()
     class  (coolingFunctionClass                   ), pointer                   :: coolingFunction_              => null()
     type   (radiationFieldCosmicMicrowaveBackground), pointer                   :: radiation                     => null()
     integer                                                                     :: radiiCount                             , elementCount_       , &
          &                                                                         abundancesCount                        , chemicalsCount      , &
          &                                                                         indexRadii                             , indexDensity
     logical                                                                     :: includeRadii                           , includeDensity
     type   (varying_string                         ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusSpecifier                        ), allocatable, dimension(:) :: radii
     logical                                                                     :: darkMatterScaleRadiusIsNeeded          , diskIsNeeded        , &
          &                                                                         spheroidIsNeeded                       , virialRadiusIsNeeded, &
          &                                                                         nuclearStarClusterIsNeeded             , satelliteIsNeeded
     type   (varying_string                         )                            :: label
   contains
     final     ::                       cgmCoolingFunctionDestructor
     procedure :: columnDescriptions => cgmCoolingFunctionColumnDescriptions
     procedure :: size               => cgmCoolingFunctionSize
     procedure :: elementCount       => cgmCoolingFunctionElementCount
     procedure :: extract            => cgmCoolingFunctionExtract
     procedure :: names              => cgmCoolingFunctionNames
     procedure :: descriptions       => cgmCoolingFunctionDescriptions
     procedure :: unitsInSI          => cgmCoolingFunctionUnitsInSI
  end type nodePropertyExtractorCGMCoolingFunction

  interface nodePropertyExtractorCGMCoolingFunction
     !!{
     Constructors for the {\normalfont \ttfamily cgmCoolingFunction} output analysis class.
     !!}
     module procedure cgmCoolingFunctionConstructorParameters
     module procedure cgmCoolingFunctionConstructorInternal
  end interface nodePropertyExtractorCGMCoolingFunction

contains

  function cgmCoolingFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily cgmCoolingFunction} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorCGMCoolingFunction)                              :: self
    type   (inputParameters                        ), intent(inout)               :: parameters
    type   (varying_string                         ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass               ), pointer                     :: darkMatterHaloScale_
    class  (coolingFunctionClass                   ), pointer                     :: coolingFunction_
    class  (cosmologyFunctionsClass                ), pointer                     :: cosmologyFunctions_
    logical                                                                       :: includeRadii        , includeDensity
    type   (varying_string                         )                              :: label
    
    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter>
      <name>radiusSpecifiers</name>
      <description>A list of radius specifiers at which to output the CGM cooling function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the radii at which density data are output should also be included in the output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeDensity</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the total hydrogen densities ($n_\mathrm{H}$) at each radius should be included in the output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>label</name>
      <defaultValue>var_str('')</defaultValue>
      <description>A label to distinguish this cooling function from others.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="coolingFunction"     name="coolingFunction_"     source="parameters"/>
    !!]
    self=nodePropertyExtractorCGMCoolingFunction(radiusSpecifiers,includeRadii,includeDensity,label,cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="coolingFunction_"    />
    !!]
    return
  end function cgmCoolingFunctionConstructorParameters

  function cgmCoolingFunctionConstructorInternal(radiusSpecifiers,includeRadii,includeDensity,label,cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily cgmCoolingFunction} property extractor class.
    !!}
    use :: Abundances_Structure                , only : Abundances_Property_Count
    use :: Chemical_Abundances_Structure       , only : Chemicals_Property_Count
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    use :: String_Handling                     , only : String_Upper_Case_First
    implicit none
    type   (nodePropertyExtractorCGMCoolingFunction)                              :: self
    type   (varying_string                         ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (cosmologyFunctionsClass                ), intent(in   ), target       :: cosmologyFunctions_
    class  (darkMatterHaloScaleClass               ), intent(in   ), target       :: darkMatterHaloScale_
    class  (coolingFunctionClass                   ), intent(in   ), target       :: coolingFunction_
    logical                                         , intent(in   )               :: includeRadii        , includeDensity
    type   (varying_string                         ), intent(in   )               :: label
    !![
    <constructorAssign variables="radiusSpecifiers, includeRadii, includeDensity, *cosmologyFunctions_, *darkMatterHaloScale_, *coolingFunction_"/>
    !!]

    ! Decode radii specifiers.
    self   %elementCount_=                   1
    if (includeRadii  ) then
       self%elementCount_=self%elementCount_+1
       self%indexRadii   =self%elementCount_
    end if
    if (includeDensity) then
       self%elementCount_=self%elementCount_+1
       self%indexDensity =self%elementCount_
    end if
    self%radiiCount      =size(radiusSpecifiers)
    call Galactic_Structure_Radii_Definition_Decode(                                    &
         &                                          radiusSpecifiers                  , &
         &                                          self%radii                        , &
         &                                          self%diskIsNeeded                 , &
         &                                          self%spheroidIsNeeded             , &
         &                                          self%nuclearStarClusterIsNeeded   , &
         &                                          self%satelliteIsNeeded            , &
         &                                          self%virialRadiusIsNeeded         , &
         &                                          self%darkMatterScaleRadiusIsNeeded  &
         &                                         )
    ! Get a count of the number of abundances and chemicals properties.
    self%abundancesCount=Abundances_Property_Count()
    self%chemicalsCount =Chemicals_Property_Count ()
    ! Initialize radiation field.
    allocate(self%radiation)
    !![
    <referenceConstruct isResult="yes" owner="self" object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
    !!]
    ! Upper case the label.
    self%label=String_Upper_Case_First(char(label))
    return
  end function cgmCoolingFunctionConstructorInternal

  subroutine cgmCoolingFunctionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily cgmCoolingFunction} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorCGMCoolingFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%coolingFunction_"    />
    <objectDestructor name="self%radiation"           />
    !!]
    return
  end subroutine cgmCoolingFunctionDestructor

  integer function cgmCoolingFunctionElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily cgmCoolingFunction} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorCGMCoolingFunction), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: time

    cgmCoolingFunctionElementCount=self%elementCount_
    return
  end function cgmCoolingFunctionElementCount

  function cgmCoolingFunctionSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily cgmCoolingFunction} property extractors.
    !!}
    implicit none
    integer         (c_size_t                               )                :: cgmCoolingFunctionSize
    class           (nodePropertyExtractorCGMCoolingFunction), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: time

    cgmCoolingFunctionSize=self%radiiCount
    return
  end function cgmCoolingFunctionSize

  function cgmCoolingFunctionExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily cgmCoolingFunction} property extractor.
    !!}
    use :: Abundances_Structure                , only : abundances
    use :: Chemical_Abundances_Structure       , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities   , only : Chemicals_Mass_To_Fraction_Conversion
    use :: Galactic_Structure_Options          , only : componentTypeAll                     , componentTypeHotHalo        , massTypeGaseous                   , massTypeGalactic                          , &
         &                                              massTypeStellar
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius      , radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius              , radiusTypeGalacticLightFraction           , &
         &                                              radiusTypeGalacticMassFraction       , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius  , radiusTypeSpheroidRadius                  , &
         &                                              radiusTypeStellarMassFraction        , radiusTypeVirialRadius      , radiusTypeNuclearStarClusterRadius, radiusTypeNuclearStarClusterHalfMassRadius    
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile       , nodeComponentDisk           , nodeComponentSpheroid             , treeNode                                  , &
         &                                              nodeComponentBasic                   , nodeComponentHotHalo        , nodeComponentNSC
    use :: Mass_Distributions                  , only : massDistributionClass                , kinematicsDistributionClass
    use :: Coordinates                         , only : coordinateSpherical                  , assignment(=)
    use :: Numerical_Constants_Astronomical    , only : massSolar                            , megaParsec
    use :: Numerical_Constants_Atomic          , only : massHydrogenAtom
    use :: Numerical_Constants_Prefixes        , only : hecto
    use :: Error                               , only : Error_Report
    implicit none
    double precision                                         , dimension(:,:), allocatable :: cgmCoolingFunctionExtract
    class           (nodePropertyExtractorCGMCoolingFunction), intent(inout) , target      :: self
    type            (treeNode                               ), intent(inout) , target      :: node
    double precision                                         , intent(in   )               :: time
    type            (multiCounter                           ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic                     )                , pointer     :: basic
    class           (nodeComponentHotHalo                   )                , pointer     :: hotHalo
    class           (nodeComponentDisk                      )                , pointer     :: disk
    class           (nodeComponentSpheroid                  )                , pointer     :: spheroid
    class           (nodeComponentNSC                       )                , pointer     :: nuclearStarCluster
    class           (nodeComponentDarkMatterProfile         )                , pointer     :: darkMatterProfile
    class           (massDistributionClass                  )                , pointer     :: massDistribution_
    class           (kinematicsDistributionClass            )                , pointer     :: kinematicsDistribution_
    type            (coordinateSpherical                    )                              :: coordinates
    integer                                                                                :: i
    double precision                                                                       :: radius                   , radiusVirial         , &
         &                                                                                    density                  , temperature          , &
         &                                                                                    massToDensityConversion  , numberDensityHydrogen
    type            (abundances                             )                              :: abundancesGas
    type            (chemicalAbundances                     )                              :: chemicalMasses           , fractionsChemical
    !$GLC attributes unused :: time, instance

    allocate(cgmCoolingFunctionExtract(self%radiiCount,self%elementCount_))
    radiusVirial                                               =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial       =  self%darkMatterHaloScale_%radiusVirial(node                    )
    if (self%                 diskIsNeeded) disk               =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid           =>                                        node%spheroid         ()
    if (self%   nuclearStarClusterIsNeeded) nuclearStarCluster =>                                        node%NSC              ()
    if (self%darkMatterScaleRadiusIsNeeded) darkMatterProfile  =>                                        node%darkMatterProfile()
    do i=1,self%radiiCount
       radius=self%radii(i)%value
       select case (self%radii(i)%type%ID)
       case   (radiusTypeRadius                          %ID)
          ! Nothing to do.
       case   (radiusTypeVirialRadius                    %ID)
          radius=+radius*radiusVirial
       case   (radiusTypeDarkMatterScaleRadius           %ID)
          radius=+radius*darkMatterProfile %         scale()
       case   (radiusTypeDiskRadius                      %ID)
          radius=+radius*disk              %        radius()
       case   (radiusTypeSpheroidRadius                  %ID)
          radius=+radius*spheroid          %        radius()
       case   (radiusTypeNuclearStarClusterRadius        %ID)
          radius=+radius*nuclearStarCluster%        radius()
       case   (radiusTypeDiskHalfMassRadius              %ID)
          radius=+radius*disk              %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius          %ID)
          radius=+radius*spheroid          %halfMassRadius()
       case   (radiusTypeNuclearStarClusterHalfMassRadius%ID)
          radius=+radius*nuclearStarCluster%halfMassRadius()
       case   (radiusTypeGalacticMassFraction            %ID,  &
            &  radiusTypeGalacticLightFraction           %ID)
          massDistribution_ =>  node             %massDistribution   (                                                &
               &                                                      massType      =              massTypeStellar ,  &
               &                                                      componentType =              componentTypeAll,  &
               &                                                      weightBy      =self%radii(i)%weightBy        ,  &
               &                                                      weightIndex   =self%radii(i)%weightByIndex      &
               &                                                     )
          radius            =  +radius                                                                                &
               &               *massDistribution_%radiusEnclosingMass(                                                &
               &                                                      massFractional=self%radii(i)%fraction           &
               &                                                     )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
        case   (radiusTypeStellarMassFraction  %ID)
           massDistribution_ =>  node             %massDistribution  (                                                &
               &                                                      massType      =              massTypeStellar ,  &
               &                                                      componentType =              componentTypeAll,  &
               &                                                      weightBy      =self%radii(i)%weightBy        ,  &
               &                                                      weightIndex   =self%radii(i)%weightByIndex      &
               &                                                     )
          radius            =  +radius                                                                                &
               &               *massDistribution_%radiusEnclosingMass(                                                &
               &                                                      massFractional=self%radii(i)%fraction           &
               &                                                     )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       case default
          call Error_Report('unrecognized radius type'//{introspection:location})
       end select
       ! Extract properties needed for the cooling function.       
       basic   => node%basic  ()
       hotHalo => node%hotHalo()
       ! Set epoch for radiation field.
       call self%radiation%timeSet(basic%time())
       ! Compute metal abundances.
       abundancesGas=hotHalo%abundances()
       call abundancesGas%massToMassFraction(hotHalo%mass())
       ! Compute the chemicals for this node.
       if (self%chemicalsCount > 0) then
          chemicalMasses=hotHalo%chemicals()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(fractionsChemical)
          ! Compute factor converting mass of chemicals in (M☉) to number density per unit total mass density (in cm⁻³ / M☉
          ! Mpc⁻³).
          if (hotHalo%mass() > 0.0d0) then
             massToDensityConversion=Chemicals_Mass_To_Fraction_Conversion(hotHalo%mass())
          else
             massToDensityConversion=0.0d0
          end if          
          ! Convert to number density per unit total mass density.
          call fractionsChemical%scale(massToDensityConversion)
       end if
       ! Get density and temperature.
       coordinates             =  [radius,0.0d0,0.0d0]
       massDistribution_       => node                   %massDistribution      (componentTypeHotHalo,massTypeGaseous)
       kinematicsDistribution_ => massDistribution_      %kinematicsDistribution(                                    )
       density                 =  massDistribution_      %density               (coordinates                         )
       temperature             =  kinematicsDistribution_%temperature           (coordinates                         )
       !![
       <objectDestructor name="massDistribution_"      />
       <objectDestructor name="kinematicsDistribution_"/>
       !!]
       ! Compute number density of hydrogen (in cm⁻³).
       numberDensityHydrogen=+density                                    &
            &                *abundancesGas   %hydrogenMassFraction()    &
            &                *massSolar                                  &
            &                /massHydrogenAtom                           &
            &                /hecto                                  **3 &
            &                /megaParsec                             **3
       if (numberDensityHydrogen > 0.0d0) then
          ! Extract the cooling function. Note that we must convert to the standard definition of "per hydrogen atom" here.
          cgmCoolingFunctionExtract    (i,                1)=+self%coolingFunction_%coolingFunction(                              &
               &                                                                                    node                     ,    &
               &                                                                                    numberDensityHydrogen    ,    &
               &                                                                                    temperature              ,    &
               &                                                                                    abundancesGas            ,    &
               &                                                                                    fractionsChemical*density,    &
               &                                                                                    self%radiation                &
               &                                                                                   )                              &
               &                                             /                                      numberDensityHydrogen     **2
       else
          cgmCoolingFunctionExtract    (i,                1)=+0.0d0
       end if
       if (self%includeRadii  )                                                                                              &
            & cgmCoolingFunctionExtract(i,self%indexRadii  )=                                       radius
       if (self%includeDensity)                                                                                              &
            & cgmCoolingFunctionExtract(i,self%indexDensity)=                                       numberDensityHydrogen       
    end do
    return
  end function cgmCoolingFunctionExtract

  subroutine cgmCoolingFunctionNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily cgmCoolingFunction} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorCGMCoolingFunction), intent(inout)                             :: self
    double precision                                         , intent(in   ), optional                   :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names(1)="cgm"//self%label//"CoolingFunction"
    if (self%includeRadii  ) names(self%indexRadii  )="cgm"//self%label//"CoolingFunctionRadius"
    if (self%includeDensity) names(self%indexDensity)="cgm"//self%label//"CoolingFunctionDensityHydrogen"
    return
  end subroutine cgmCoolingFunctionNames

  subroutine cgmCoolingFunctionDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily cgmCoolingFunction} property.
    !!}
    implicit none
    class           (nodePropertyExtractorCGMCoolingFunction), intent(inout)                             :: self
    double precision                                         , intent(in   ), optional                   :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions       (1)="CGM cooling function at a given radius [ergs cm³ s¯¹]."
    if (self%includeRadii  )                                                                                                          &
         & descriptions(self%indexRadii  )="Radius at which the CGM cooling function is output [Mpc]."
    if (self%includeDensity)                                                                                                          &
         & descriptions(self%indexDensity)="Total hydrogen density at the radius at which the CGM cooling function is output [cm¯³]."
    return
  end subroutine cgmCoolingFunctionDescriptions

  subroutine cgmCoolingFunctionColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily cgmCoolingFunction} property.
    !!}
    implicit none
    class           (nodePropertyExtractorCGMCoolingFunction), intent(inout)                            :: self
    double precision                                         , intent(in   ), optional                  :: time
    type            (varying_string                         ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                         , intent(inout), dimension(:), allocatable :: values 
    type            (varying_string                         ), intent(  out)                            :: valuesDescription
    double precision                                         , intent(  out)                            :: valuesUnitsInSI
    !$GLC attributes unused :: time

    allocate(descriptions(self%radiiCount))
    allocate(values      (              0))
    valuesDescription=var_str('')
    valuesUnitsInSI  =0.0d0
    descriptions     =self%radii%name
    return
  end subroutine cgmCoolingFunctionColumnDescriptions

  function cgmCoolingFunctionUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily cgmCoolingFunction} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : centi
    use :: Numerical_Constants_Units       , only : ergs
    implicit none
    double precision                                         , allocatable  , dimension(:) :: cgmCoolingFunctionUnitsInSI
    class           (nodePropertyExtractorCGMCoolingFunction), intent(inout)               :: self
    double precision                                         , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(cgmCoolingFunctionUnitsInSI(self%elementCount_))
    cgmCoolingFunctionUnitsInSI       (                1)=ergs *centi     **3
    if (self%includeRadii  )                                                  &
         & cgmCoolingFunctionUnitsInSI(self%indexRadii  )=      megaParsec
    if (self%includeDensity)                                                  &
         & cgmCoolingFunctionUnitsInSI(self%indexDensity)=1.0d0/centi     **3
    return
  end function cgmCoolingFunctionUnitsInSI

