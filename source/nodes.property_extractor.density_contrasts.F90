!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Implements a property extractor class for the mass and radii of spheres are specified density contrast.
!!}

  use :: Cosmology_Functions       , only : cosmologyFunctions     , cosmologyFunctionsClass , enumerationDensityCosmologicalType
  use :: Cosmology_Parameters      , only : cosmologyParameters    , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale    , darkMatterHaloScaleClass
  use :: Galactic_Structure_Options, only : enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Root_Finder               , only : rootFinder

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDensityContrasts">
   <description>
    A property extractor class for the mass and radii of spheres of specified density contrast. A list of density contrasts,
    $\Delta$ (defined in units of the mean density of the Universe), is specified via the {\normalfont \ttfamily
    [densityContrasts]} parameter. For each specified density contrast, two properties are output for each node: {\normalfont
    \ttfamily nodeRadius}$\Delta$ and {\normalfont \ttfamily nodeMass}$\Delta$ which give the radius enclosing a mean density
    contrast of $\Delta$ and the mass enclosed within that radius. The parameter {\normalfont \ttfamily [darkMatterOnly]}
    controls whether density contrasts are measured for total mass ({\normalfont \ttfamily false}) or dark matter mass only
    ({\normalfont \ttfamily true}). In the latter case, density contrasts are defined relative to the mean dark matter density
    of the Universe.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorDensityContrasts
     !!{
     A property extractor class for the mass and radii of spheres are specified density contrast.
     !!}
     private
     class           (cosmologyParametersClass          ), pointer                   :: cosmologyParameters_      => null()
     class           (cosmologyFunctionsClass           ), pointer                   :: cosmologyFunctions_       => null()
     class           (darkMatterHaloScaleClass          ), pointer                   :: darkMatterHaloScale_      => null()
     type            (rootFinder                        )                            :: finder
     integer                                                                         :: elementCount_                      , countDensityContrasts
     type            (enumerationMassTypeType           )                            :: massTypeSelected
     type            (enumerationDensityCosmologicalType)                            :: densityContrastRelativeTo
     logical                                                                         :: darkMatterOnly
     double precision                                    , allocatable, dimension(:) :: densityContrasts
   contains
     final     ::                       densityContrastsDestructor
     procedure :: size               => densityContrastsSize
     procedure :: elementCount       => densityContrastsElementCount
     procedure :: extract            => densityContrastsExtract
     procedure :: names              => densityContrastsNames
     procedure :: descriptions       => densityContrastsDescriptions
     procedure :: columnDescriptions => densityContrastsColumnDescriptions
     procedure :: unitsInSI          => densityContrastsUnitsInSI
  end type nodePropertyExtractorDensityContrasts

  interface nodePropertyExtractorDensityContrasts
     !!{
     Constructors for the \refClass{nodePropertyExtractorDensityContrasts} output analysis class.
     !!}
     module procedure densityContrastsConstructorParameters
     module procedure densityContrastsConstructorInternal
  end interface nodePropertyExtractorDensityContrasts

  ! Module-scope variables used in root finding.
  class           (massDistributionClass), pointer :: massDistribution_
  double precision                                 :: densityTarget
  !$omp threadprivate(massDistribution_,densityTarget)

contains

  function densityContrastsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorDensityContrasts} property extractor class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : enumerationDensityCosmologicalEncode
    use :: Input_Parameters   , only : inputParameter                      , inputParameters
    implicit none
    type            (nodePropertyExtractorDensityContrasts)                              :: self
    type            (inputParameters                      ), intent(inout)               :: parameters
    class           (cosmologyParametersClass             ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass              ), pointer                     :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass             ), pointer                     :: darkMatterHaloScale_
    double precision                                       , allocatable  , dimension(:) :: densityContrasts
    logical                                                                              :: darkMatterOnly
    type            (varying_string                       )                              :: densityContrastRelativeTo
    
    allocate(densityContrasts(parameters%count('densityContrasts')))
    !![
    <inputParameter>
      <name>densityContrasts</name>
      <description>A list of density contrasts at which to output data.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityContrastRelativeTo</name>
      <description>The density ({\normalfont \ttfamily mean} or {\normalfont \ttfamily critical}) used in defining the density contrast.</description>
      <source>parameters</source>
      <defaultValue>var_str('mean')</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>darkMatterOnly</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not density contrast data should be computed using the dark matter component alone.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorDensityContrasts(densityContrasts,darkMatterOnly,enumerationDensityCosmologicalEncode(char(densityContrastRelativeTo),includesPrefix=.false.),cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function densityContrastsConstructorParameters

  function densityContrastsConstructorInternal(densityContrasts,darkMatterOnly,densityContrastRelativeTo,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorDensityContrasts} property extractor class.
    !!}
    use :: Galactic_Structure_Options, only : massTypeAll              , massTypeDark
    use :: Root_Finder               , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type            (nodePropertyExtractorDensityContrasts)                              :: self
    class           (cosmologyFunctionsClass              ), intent(in   ), target       :: cosmologyFunctions_
    class           (cosmologyParametersClass             ), intent(in   ), target       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass             ), intent(in   ), target       :: darkMatterHaloScale_
    double precision                                       , intent(in   ), dimension(:) :: densityContrasts
    logical                                                , intent(in   )               :: darkMatterOnly
    type            (enumerationDensityCosmologicalType   ), intent(in   )               :: densityContrastRelativeTo
    double precision                                       , parameter                   :: toleranceAbsolute        =0.0d0, toleranceRelative=1.0d-3
    !![
    <constructorAssign variables="densityContrasts, darkMatterOnly, densityContrastRelativeTo, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]

    self%countDensityContrasts=size(densityContrasts)
    self%elementCount_        =2
    select case (darkMatterOnly)
    case (.true.)
       self%massTypeSelected=massTypeDark
    case (.false.)
       self%massTypeSelected=massTypeAll
    end select
    self%finder=rootFinder(                                                             &
         &                 rootFunction                 =densityContrastsRoot         , &
         &                 rangeExpandDownward          =0.5d0                        , &
         &                 rangeExpandUpward            =2.0d0                        , &
         &                 rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
         &                 rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                 rangeExpandType              =rangeExpandMultiplicative    , &
         &                 toleranceAbsolute            =toleranceAbsolute            , &
         &                 toleranceRelative            =toleranceRelative              &
         &                )
    return
  end function densityContrastsConstructorInternal

  subroutine densityContrastsDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorDensityContrasts} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDensityContrasts), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine densityContrastsDestructor

  integer function densityContrastsElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily densityContrasts} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityContrasts), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    densityContrastsElementCount=self%elementCount_
    return
  end function densityContrastsElementCount

  function densityContrastsSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily densityContrasts} property extractors.
    !!}
    implicit none
    integer         (c_size_t                             )                :: densityContrastsSize
    class           (nodePropertyExtractorDensityContrasts), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    densityContrastsSize=self%countDensityContrasts
    return
  end function densityContrastsSize

  function densityContrastsExtract(self,node,time,instance)
    !!{
    Implement a last isolated redshift output analysis.
    !!}
    use :: Cosmology_Functions       , only : densityCosmologicalMean, densityCosmologicalCritical
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure_Options, only : componentTypeAll
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    implicit none
    double precision                                        , dimension(:,:), allocatable :: densityContrastsExtract
    class           (nodePropertyExtractorDensityContrasts ), intent(inout) , target      :: self
    type            (treeNode                              ), intent(inout) , target      :: node
    double precision                                        , intent(in   )               :: time
    type            (multiCounter                          ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic                    )                , pointer     :: basic
    double precision                                        , parameter                   :: radiusTiny             =1.0d-12
    integer                                                                               :: i
    double precision                                                                      :: enclosedMass                   , radius, &
         &                                                                                   densityReference
    !$GLC attributes unused :: time, instance

    allocate(densityContrastsExtract(self%countDensityContrasts,self%elementCount_))
    ! Find the reference density at this epoch.
    basic            => node                    %basic               (            )
    densityReference =  self%cosmologyFunctions_%matterDensityEpochal(basic%time())
    select case (self%densityContrastRelativeTo%ID)
    case (densityCosmologicalMean    %ID)
       ! No modification required.
    case (densityCosmologicalCritical%ID)
       ! Modify reference density to be the critical density.
       densityReference=+densityReference                                                          &
            &           /self%cosmologyFunctions_%OmegaMatterEpochal(basic%time())
    case default
       call Error_Report('unknown cosmological density'//{introspection:location})
    end select
    ! If dark matter only is used, multiply the reference density by the dark matter fraction.
    if (self%darkMatterOnly) densityReference=+ densityReference                                                                 &
         &                                    *(self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon()) &
         &                                    / self%cosmologyParameters_%OmegaMatter()
    ! Iterate over density contrasts.
    massDistribution_ => node%massDistribution(massType=self%massTypeSelected)
    do i=1,self%countDensityContrasts
       densityTarget=self%densityContrasts(i)*densityReference
       if (densityContrastsRoot(radiusTiny*self%darkMatterHaloScale_%radiusVirial(node)) < 0.0d0) then
          ! The target density contrast is not reached even at this tiny radius. This happens in cored density profiles. Return
          ! zero mass and radius.
          radius      =0.0d0
          enclosedMass=0.0d0
       else
          ! The target density is reached, so find the exact radius at which it occurs.
          radius      =self             %finder%find                (rootGuess=self%darkMatterHaloScale_%radiusVirial(node))
          enclosedMass=massDistribution_       %massEnclosedBySphere(                                    radius            )
       end if
       densityContrastsExtract(i,:)=[radius,enclosedMass]
    end do
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function densityContrastsExtract

  subroutine densityContrastsNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily densityContrasts} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityContrasts), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names(1)='nodeRadius'
    names(2)='nodeMass'
    return
  end subroutine densityContrastsNames

  subroutine densityContrastsDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily densityContrasts} property.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityContrasts), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions(1)='Radius enclosing a given density contrast [Mpc].'
    descriptions(2)='Mass within a given density contrast [M☉].'
    return
  end subroutine densityContrastsDescriptions

  subroutine densityContrastsColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily densityContrasts} property.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityContrasts), intent(inout)                            :: self
    double precision                                       , intent(in   ), optional                  :: time
    type            (varying_string                       ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                       , intent(inout), dimension(:), allocatable :: values
    type            (varying_string                       ), intent(  out)                            :: valuesDescription
    double precision                                       , intent(  out)                            :: valuesUnitsInSI
    character       (len=32                               )                                           :: label
    integer         (c_size_t                             )                                           :: i
    !$GLC attributes unused :: time

    allocate(descriptions(self%countDensityContrasts))
    allocate(values      (                         0))
    do i=1,self%countDensityContrasts
       write (label,'(a2,f9.4)') 'Δ=',self%densityContrasts(i)
       descriptions(i)=trim(label)
    end do
    valuesDescription=var_str('')
    valuesUnitsInSI  =0.0d0
    return
  end subroutine densityContrastsColumnDescriptions

  function densityContrastsUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily densityContrasts} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                       , allocatable  , dimension(:) :: densityContrastsUnitsInSI
    class           (nodePropertyExtractorDensityContrasts), intent(inout)               :: self
    double precision                                       , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(densityContrastsUnitsInSI(self%elementCount_))
    densityContrastsUnitsInSI(1:2)=[megaParsec,massSolar]
    return
  end function densityContrastsUnitsInSI


  double precision function densityContrastsRoot(radius)
    !!{
    Root function used in finding the radius that encloses a given density contrast.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    double precision, intent(in   ) :: radius
    
    densityContrastsRoot=+3.0d0                                          &
         &               *massDistribution_%massEnclosedBySphere(radius) &
         &               /4.0d0                                          &
         &               /Pi                                             &
         &               /radius**3                                      &
         &               -densityTarget
    return
  end function densityContrastsRoot
