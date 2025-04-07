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
  Implements a property extractor class for the velocity dispersion at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale    , darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier
  use :: Galactic_Structure_Options          , only : enumerationMassTypeType, enumerationComponentTypeType, enumerationWeightByType
  use :: Mass_Distributions                  , only : massDistributionClass  , kinematicsDistributionClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorVelocityDispersion">
   <description>
    A property extractor class for the velocity dispersion at a set of radii. The radii and types of projected density to output
    is specified by the {\normalfont \ttfamily radiusSpecifiers} parameter. This parameter's value can contain multiple
    entries, each of which should be a valid
    \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#sec.radiusSpecifiers}{radius
    specifier}, but with an additional, colon-separated, value at the end indicating the direction in which the velocity
    dispersion should be computed. This direction should be one of {\normalfont \ttfamily radial} (computes the radial
    component of velocity dispersion), {\normalfont \ttfamily lineOfSight\{\textless luminosity\textgreater\}} (computes the
    line-of-sight velocity dispersion), {\normalfont \ttfamily lineOfSightInteriorAverage\{\textless luminosity\textgreater\}}
    (computes the line-of-sight velocity dispersion averaged interior to the given radius), or {\normalfont \ttfamily
    lambdaR\{\textless luminosity\textgreater\}} (computes the $\lambda_\mathrm{R}$ statistic of
    \citealt{cappellari_sauron_2007})---in the latter three cases {\normalfont \ttfamily \{\textless luminosity\textgreater\}}
    specifies which band should be used to weight the velocity dispersion, alternatively setting {\normalfont \ttfamily
    \{\textless luminosity\textgreater\}}$=${\normalfont \ttfamily mass} (or just leaving off this specifier entirely) will use
    mass weighting instead.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorVelocityDispersion
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_          => null()
     integer                                                      :: radiiCount                             , elementCount_
     logical                                                      :: includeRadii                           , integrationFailureIsFatal
     double precision                                             :: toleranceRelative
     type   (varying_string          ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusSpecifier         ), allocatable, dimension(:) :: radii
     logical                                                      :: darkMatterScaleRadiusIsNeeded          , diskIsNeeded        , &
          &                                                          spheroidIsNeeded                       , virialRadiusIsNeeded, &
          &                                                          nuclearStarClusterIsNeeded             , satelliteIsNeeded   , &
          &                                                          hotHaloIsNeeded
   contains
     final     ::                       velocityDispersionDestructor
     procedure :: columnDescriptions => velocityDispersionColumnDescriptions
     procedure :: size               => velocityDispersionSize
     procedure :: elementCount       => velocityDispersionElementCount
     procedure :: extract            => velocityDispersionExtract
     procedure :: names              => velocityDispersionNames
     procedure :: descriptions       => velocityDispersionDescriptions
     procedure :: unitsInSI          => velocityDispersionUnitsInSI
  end type nodePropertyExtractorVelocityDispersion

  interface nodePropertyExtractorVelocityDispersion
     !!{
     Constructors for the {\normalfont \ttfamily velocityDispersion} output analysis class.
     !!}
     module procedure velocityDispersionConstructorParameters
     module procedure velocityDispersionConstructorInternal
  end interface nodePropertyExtractorVelocityDispersion

  ! Module-scope variables used in integrands.
  class           (massDistributionClass                  ), pointer :: massDistribution_                     , massDistributionStellarDisk_      , &
       &                                                                massDistributionWeighted_             , massDistributionStellarSpheroid_  , &
       &                                                                massDistributionTotal_
  class           (kinematicsDistributionClass            ), pointer :: kinematicsDistribution_               , kinematicsDistributionStellarDisk_, &
       &                                                                kinematicsDistributionStellarSpheroid_
  class           (nodePropertyExtractorVelocityDispersion), pointer :: self_
  type            (enumerationMassTypeType                )          :: massType_
  type            (enumerationComponentTypeType           )          :: componentType_
  type            (enumerationWeightByType                )          :: weightBy_
  integer                                                            :: weightIndex_
  double precision                                                   :: radiusImpact_ , radiusOuter_
  !$omp threadprivate(massDistribution_,massDistributionWeighted_,massDistributionStellarDisk_,massDistributionStellarSpheroid_,massDistributionTotal_,kinematicsDistribution_,kinematicsDistributionStellarDisk_,kinematicsDistributionStellarSpheroid_,self_,weightBy_,componentType_,massType_,weightIndex_,radiusImpact_,radiusOuter_)

contains

  function velocityDispersionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily velocityDispersion} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorVelocityDispersion)                              :: self
    type            (inputParameters                        ), intent(inout)               :: parameters
    type            (varying_string                         ), allocatable  , dimension(:) :: radiusSpecifiers
    class           (darkMatterHaloScaleClass               ), pointer                     :: darkMatterHaloScale_
    double precision                                                                       :: toleranceRelative
    logical                                                                                :: includeRadii        , integrationFailureIsFatal

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter>
      <name>radiusSpecifiers</name>
      <description>A list of radius specifiers at which to output the velocity dispersion.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the radii at which velocity dispersion data are output should also be included in the output file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>integrationFailureIsFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, failure of line-of-sight integrals is fatal. Otherwise, such errors are tolerated.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelative</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The relative tolerance to use in integrals.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorVelocityDispersion(radiusSpecifiers,includeRadii,integrationFailureIsFatal,toleranceRelative,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function velocityDispersionConstructorParameters

  function velocityDispersionConstructorInternal(radiusSpecifiers,includeRadii,integrationFailureIsFatal,toleranceRelative,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily velocityDispersion} property extractor class.
    !!}
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type            (nodePropertyExtractorVelocityDispersion)                              :: self
    type            (varying_string                         ), intent(in   ), dimension(:) :: radiusSpecifiers
    class           (darkMatterHaloScaleClass               ), intent(in   ), target       :: darkMatterHaloScale_
    logical                                                  , intent(in   )               :: includeRadii        , integrationFailureIsFatal
    double precision                                         , intent(in   )               :: toleranceRelative
    !![
    <constructorAssign variables="radiusSpecifiers, includeRadii, integrationFailureIsFatal, toleranceRelative, *darkMatterHaloScale_"/>
    !!]

    if (includeRadii) then
       self%elementCount_=2
    else
       self%elementCount_=1
    end if
    self%radiiCount      =size(radiusSpecifiers)
    call Galactic_Structure_Radii_Definition_Decode(                                    &
         &                                          radiusSpecifiers                  , &
         &                                          self%radii                        , &
         &                                          self%hotHaloIsNeeded              , &
         &                                          self%diskIsNeeded                 , &
         &                                          self%spheroidIsNeeded             , &
         &                                          self%nuclearStarClusterIsNeeded   , &
         &                                          self%satelliteIsNeeded            , &
         &                                          self%virialRadiusIsNeeded         , &
         &                                          self%darkMatterScaleRadiusIsNeeded  &
         &                                         )
    return
  end function velocityDispersionConstructorInternal

  subroutine velocityDispersionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily velocityDispersion} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorVelocityDispersion), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine velocityDispersionDestructor

  integer function velocityDispersionElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily velocityDispersion} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorVelocityDispersion), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: time

    velocityDispersionElementCount=self%elementCount_
    return
  end function velocityDispersionElementCount

  function velocityDispersionSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily velocityDispersion} property extractors.
    !!}
    implicit none
    integer         (c_size_t                               )                :: velocityDispersionSize
    class           (nodePropertyExtractorVelocityDispersion), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: time

    velocityDispersionSize=self%radiiCount
    return
  end function velocityDispersionSize

  function velocityDispersionExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily velocityDispersion} property extractor.
    !!}
    use :: Galactic_Structure_Options          , only : componentTypeAll               , componentTypeDisk           , componentTypeSpheroid                     , massTypeGalactic                  , &
          &                                             massTypeStellar                , massTypeAll
    use :: Galactic_Structure_Radii_Definitions, only : directionLambdaR               , directionLineOfSight        , directionLineOfSightInteriorAverage       , directionRadial                   , &
          &                                             radiusTypeDarkMatterScaleRadius, radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius                      , radiusTypeGalacticLightFraction   , &
          &                                             radiusTypeGalacticMassFraction , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius          , radiusTypeSpheroidRadius          , &
          &                                             radiusTypeStellarMassFraction  , radiusTypeVirialRadius      , radiusTypeNuclearStarClusterHalfMassRadius, radiusTypeNuclearStarClusterRadius, &
          &                                             radiusTypeHotHaloOuterRadius
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile , nodeComponentDisk           , nodeComponentSpheroid                     , nodeComponentNSC                  , &
          &                                             nodeComponentHotHalo           , treeNode
    use :: Coordinates                         , only : coordinateSpherical            , assignment(=)
    use :: Numerical_Integration               , only : integrator
    use :: Error                               , only : Error_Report
    implicit none
    double precision                                         , dimension(:,:), allocatable :: velocityDispersionExtract
    class           (nodePropertyExtractorVelocityDispersion), intent(inout) , target      :: self
    type            (treeNode                               ), intent(inout) , target      :: node
    double precision                                         , intent(in   )               :: time
    type            (multiCounter                           ), intent(inout) , optional    :: instance
    class           (nodeComponentHotHalo                   ), pointer                     :: hotHalo
    class           (nodeComponentDisk                      ), pointer                     :: disk
    class           (nodeComponentSpheroid                  ), pointer                     :: spheroid
    class           (nodeComponentNSC                       ), pointer                     :: nuclearStarCluster
    class           (nodeComponentDarkMatterProfile         ), pointer                     :: darkMatterProfile
    double precision                                         , parameter                   :: outerRadiusMultiplier           =10.0d0
    integer                                                                                :: i
    double precision                                                                       :: radius                                 , radiusVirial            , &
         &                                                                                    radiusFromFraction                     , densityIntegrand        , &
         &                                                                                    radiusZero                             , velocityDensityIntegrand, &
         &                                                                                    numerator                              , denominator             , &
         &                                                                                    massDisk                               , massSpheroid
    logical                                                                                :: scaleIsZero
    type            (integrator                             )                              :: integratorVelocitySurfaceDensity       , integratorSurfaceDensity, &
         &                                                                                    integratorLambdaR2                     , integratorLambdaR1
    type            (coordinateSpherical                    )                              :: coordinates
    !$GLC attributes unused :: time, instance

    integratorVelocitySurfaceDensity=integrator(velocityDispersionVelocitySurfaceDensityIntegrand,toleranceRelative=1.0d-3)
    integratorSurfaceDensity        =integrator(velocityDispersionSurfaceDensityIntegrand        ,toleranceRelative=1.0d-3)
    integratorLambdaR1              =integrator(velocityDispersionLambdaRIntegrand1              ,toleranceRelative=1.0d-2)
    integratorLambdaR2              =integrator(velocityDispersionLambdaRIntegrand2              ,toleranceRelative=1.0d-2)
    allocate(velocityDispersionExtract(self%radiiCount,self%elementCount_))
    radiusVirial                                         =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial       =  self%darkMatterHaloScale_%radiusVirial(node                    )
    if (self%              hotHaloIsNeeded) hotHalo            =>                                        node%hotHalo          ()
    if (self%                 diskIsNeeded) disk               =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid           =>                                        node%spheroid         ()
    if (self%   nuclearStarClusterIsNeeded) nuclearStarCluster =>                                        node%NSC              ()
    if (self%darkMatterScaleRadiusIsNeeded) darkMatterProfile  =>                                        node%darkMatterProfile()
    do i=1,self%radiiCount
       scaleIsZero=.false.
       radius     =self%radii(i)%value
       select case (self%radii(i)%type%ID)
       case   (radiusTypeRadius                          %ID)
          radiusOuter_=    radius                                    *outerRadiusMultiplier
       case   (radiusTypeVirialRadius                    %ID)
          radius                       =    radius*radiusVirial
          radiusOuter_=max(radius,radiusVirial                      )*outerRadiusMultiplier
       case   (radiusTypeDarkMatterScaleRadius           %ID)
          radius                       =    radius*darkMatterProfile %         scale()
          radiusOuter_=max(radius,darkMatterProfile%         scale())*outerRadiusMultiplier
       case   (radiusTypeHotHaloOuterRadius              %ID)
          radius                       =    radius*hotHalo           %   outerRadius()
          radiusOuter_=max(radius,hotHalo          %   outerRadius())*outerRadiusMultiplier
          scaleIsZero                  =(hotHalo                    %   outerRadius() <= 0.0d0)
       case   (radiusTypeDiskRadius                      %ID)
          radius                       =    radius*disk              %        radius()
          radiusOuter_=max(radius,disk             %        radius())*outerRadiusMultiplier
          scaleIsZero                  =(disk                       %        radius() <= 0.0d0)
       case   (radiusTypeSpheroidRadius                  %ID)
          radius                       =    radius*spheroid          %        radius()
          radiusOuter_=max(radius,spheroid         %        radius())*outerRadiusMultiplier
          scaleIsZero                  =(spheroid                   %        radius() <= 0.0d0)
       case   (radiusTypeNuclearStarClusterRadius        %ID)
          radius                       =    radius*nuclearStarCluster%        radius()
          radiusOuter_=max(radius,nuclearStarCluster              %        radius())*outerRadiusMultiplier
          scaleIsZero                  =(nuclearStarCluster                        %        radius() <= 0.0d0)
       case   (radiusTypeDiskHalfMassRadius              %ID)
          radius                       =    radius*disk              %halfMassRadius()
          radiusOuter_=max(radius,disk             %halfMassRadius())*outerRadiusMultiplier
          scaleIsZero                  =(disk                       %halfMassRadius() <= 0.0d0)
       case   (radiusTypeSpheroidHalfMassRadius          %ID)
          radius                       =    radius*spheroid          %halfMassRadius()
          radiusOuter_=max(radius,spheroid         %halfMassRadius())*outerRadiusMultiplier
          scaleIsZero                  =(spheroid                   %halfMassRadius() <= 0.0d0)
       case   (radiusTypeNuclearStarClusterHalfMassRadius%ID)
          radius                       =    radius*nuclearStarCluster              %halfMassRadius()
          radiusOuter_=max(radius,nuclearStarCluster                %halfMassRadius())*outerRadiusMultiplier
          scaleIsZero                  =(nuclearStarCluster                        %halfMassRadius() <= 0.0d0)
       case   (radiusTypeGalacticMassFraction            %ID,  &
            &  radiusTypeGalacticLightFraction           %ID )
          massDistribution_  =>  node             %massDistribution   (                                                &
               &                                                       massType      =              massTypeGalactic,  &
               &                                                       componentType =              componentTypeAll,  &
               &                                                       weightBy      =self%radii(i)%weightBy        ,  &
               &                                                       weightIndex   =self%radii(i)%weightByIndex      &
               &                                                      )
          radiusFromFraction =  +massDistribution_%radiusEnclosingMass(                                                &
               &                                                       massFractional=self%radii(i)%fraction           &
               &                                                      )
          radius             =  +radius*radiusFromFraction
          radiusOuter_       =  max(radius,radiusFromFraction)*outerRadiusMultiplier
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       case   (radiusTypeStellarMassFraction            %ID)
          massDistribution_  =>  node             %massDistribution   (                                                &
               &                                                       massType      =              massTypeStellar ,  &
               &                                                       componentType =              componentTypeAll,  &
               &                                                       weightBy      =self%radii(i)%weightBy        ,  &
               &                                                       weightIndex   =self%radii(i)%weightByIndex      &
               &                                                      )
          radiusFromFraction =  +massDistribution_%radiusEnclosingMass(                                                &
               &                                                       massFractional=self%radii(i)%fraction           &
               &                                                      )
          radius             =  +radius*radiusFromFraction
          radiusOuter_       =  max(radius,radiusFromFraction)*outerRadiusMultiplier
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       case default
          call Error_Report('unrecognized radius type'//{introspection:location})
       end select
       if (scaleIsZero) then
          ! Do not compute dispersions if the component scale is zero.
          velocityDispersionExtract(i,1)=0.0d0
       else
          massDistribution_         => node             %      massDistribution(componentType=self%radii(i)%component       ,massType=self%radii(i)%mass                                                                                )
          massDistributionWeighted_ => node             %      massDistribution(componentType=self%radii(i)%component       ,massType=self%radii(i)%mass        ,weightBy=self%radii(i)%weightBy,weightIndex=self%radii(i)%weightByIndex)
          massDistributionTotal_    => node             %      massDistribution(componentType=              componentTypeAll,massType=              massTypeAll                                                                         )
          kinematicsDistribution_   => massDistribution_%kinematicsDistribution(                                                                                                                                                        ) 
          select case (self%radii(i)%direction%ID)
          case (directionRadial                    %ID)
             ! Radial velocity dispersion.
             coordinates                   =[radius,0.0d0,0.0d0]
             velocityDispersionExtract(i,1)=kinematicsDistribution_%velocityDispersion1D(coordinates,massDistribution_,massDistributionTotal_)
          case (directionLineOfSight               %ID)
             ! Line-of-sight velocity dispersion.
             self_               => self
             massType_           =  self%radii(i)%mass
             componentType_      =  self%radii(i)%component
             weightBy_           =  self%radii(i)%integralWeightBy
             weightIndex_        =  self%radii(i)%integralWeightByIndex
             radiusImpact_       =  radius
             velocityDispersionExtract      (i,1) =  velocityDispersionLineOfSightVelocityDispersionIntegrand(radius)
          case (directionLineOfSightInteriorAverage%ID)
             ! Average over the line-of-sight velocity dispersion within the radius.
             self_          => self
             massType_      =  self%radii(i)%mass
             componentType_ =  self%radii(i)%component
             weightBy_      =  self%radii(i)%integralWeightBy
             weightIndex_   =  self%radii(i)%integralWeightByIndex
             radiusZero                      =  0.0d0
             radiusImpact_  =  radius
             velocityDensityIntegrand        =integratorVelocitySurfaceDensity%integrate(radiusZero,radiusOuter_)
             densityIntegrand                =integratorSurfaceDensity        %integrate(radiusZero,radiusOuter_)
             if (velocityDensityIntegrand <= 0.0d0) then
                velocityDispersionExtract(i,1)=0.0d0
             else
                velocityDispersionExtract(i,1)=sqrt(velocityDensityIntegrand/densityIntegrand)
             end if
          case (directionLambdaR                   %ID)
             ! The "lambdaR" parameter of Cappellari et al. (2007; MNRAS; 379; 418)
             self_          => self
             massType_      =  self%radii(i)%mass
             componentType_ =  self%radii(i)%component
             weightBy_      =  self%radii(i)%integralWeightBy
             weightIndex_   =  self%radii(i)%integralWeightByIndex
             ! Check the total masses of the disk and spheroid components. If either is zero we can use the solutions for the
             ! appropriate limiting case.
             massDistributionStellarDisk_     => node%massDistribution(componentType=componentTypeDisk    ,massType=massTypeStellar,weightBy=weightBy_,weightIndex=weightIndex_)
             massDistributionStellarSpheroid_ => node%massDistribution(componentType=componentTypeSpheroid,massType=massTypeStellar,weightBy=weightBy_,weightIndex=weightIndex_)
             massSpheroid=massDistributionStellarSpheroid_%massTotal()
             massDisk    =massDistributionStellarDisk_    %massTotal()
             if      (massDisk     <= 0.0d0) then
                velocityDispersionExtract(i,1)=0.0d0
             else if (massSpheroid <= 0.0d0) then
                velocityDispersionExtract(i,1)=1.0d0
             else
                ! Full calculation is required.
                radiusZero             =  0.0d0
                numerator              =  integratorLambdaR2%integrate(radiusZero,radius)
                denominator            =  integratorLambdaR1%integrate(radiusZero,radius)
                if (denominator <= 0.0d0) then
                   velocityDispersionExtract(i,1)=0.0d0
                else
                   velocityDispersionExtract(i,1)=numerator/denominator
                end if
             end if
             !![
	     <objectDestructor name="massDistributionStellarDisk_"    />
	     <objectDestructor name="massDistributionStellarSpheroid_"/>
	     !!]
          end select
          !![
	  <objectDestructor name="massDistribution_"        />
	  <objectDestructor name="massDistributionWeighted_"/>
	  <objectDestructor name="kinematicsDistribution_"  />
	  <objectDestructor name="massDistributionTotal_"   />
	  !!]
       end if
       if (self%includeRadii)                       &
            & velocityDispersionExtract(i,2)=radius
    end do
    return
  end function velocityDispersionExtract

  subroutine velocityDispersionNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily velocityDispersion} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorVelocityDispersion), intent(inout)                             :: self
    double precision                                         , intent(in   ), optional                   :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time
    
    allocate(names(self%elementCount_))
    names       (1)="velocityDispersion"
    if (self%includeRadii)                                       &
         & names(2)="velocityDispersionRadius"
    return
  end subroutine velocityDispersionNames

  subroutine velocityDispersionDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily velocityDispersion} property.
    !!}
    implicit none
    class           (nodePropertyExtractorVelocityDispersion), intent(inout)                             :: self
    double precision                                         , intent(in   ), optional                   :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions       (1)="Velocity dispersion at a given radius [km s⁻¹]."
    if (self%includeRadii)                                                                          &
         & descriptions(2)="Radius at which velocity dispersion is output [Mpc]."
    return
  end subroutine velocityDispersionDescriptions

  subroutine velocityDispersionColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily velocityDispersion} property.
    !!}
    implicit none
    class           (nodePropertyExtractorVelocityDispersion), intent(inout)                            :: self
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
  end subroutine velocityDispersionColumnDescriptions

  function velocityDispersionUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily velocityDispersion} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                         , allocatable  , dimension(:) :: velocityDispersionUnitsInSI
    class           (nodePropertyExtractorVelocityDispersion), intent(inout)               :: self
    double precision                                         , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(velocityDispersionUnitsInSI(self%elementCount_))
    velocityDispersionUnitsInSI       (1)=kilo
    if (self%includeRadii)                           &
         & velocityDispersionUnitsInSI(2)=megaParsec
    return
  end function velocityDispersionUnitsInSI

  double precision function velocityDispersionLambdaRIntegrand1(radius)
    !!{
    Integrand function used for integrating the $\lambda_\mathrm{R}$ statistic of \cite{cappellari_sauron_2007}. In this case we
    want to evaluate
    \begin{equation}
    \int_0^r 2 \pi r^\prime \Sigma(r^\prime) \sqrt{\sigma^2(r^\prime)+V^2(r^\prime)} \mathrm{d}r^\prime,
    \end{equation}
    where $\Sigma(r)$ is the projected surface density (in mass or light) of the galaxy at radius $r$, $\sigma^2(r)$ is the
    measured velocity dispersion and $V(r)$ the measured rotation speed. Assuming that the selected component is purely
    dispersion dominated with velocity dispersion $\sigma_\mathrm{s}(r)$, and that rotation is present in only the disk component
    with rotation curve $V_\mathrm{d}(r)$ then we can model the velocity distribution, $P(V)$, at $r$ as the sum of a Gaussian of
    width $\sigma_\mathrm{s}(r)$ and normalized area $\Sigma_\mathrm{s}(r)$, and a delta function at $V_\mathrm{d}(r)$ with normalized
    area $\Sigma_\mathrm{d}(r)$. The measured rotation speed is then:
    \begin{equation}
    V(r) = \left. \int_{-\infty}^{+\infty} P(V) V \mathrm{d}V \right/ \int_{-\infty}^{+\infty} P(V) \mathrm{d}V = {\Sigma_\mathrm{d}(r) V_\mathrm{d}(r) \over [\Sigma_\mathrm{d}(r)+\Sigma_\mathrm{s}(r)]},
    \end{equation}
    and the measured velocity dispersion is:
    \begin{equation}
    \sigma^2(r) = \left. \int_{-\infty}^{+\infty} P(V) [V-V(r)]^2 \mathrm{d}V \right/ \int_{-\infty}^{+\infty} P(V) \mathrm{d}V = {  \Sigma_\mathrm{s}(r) [\sigma_\mathrm{s}^2(r)] + \Sigma_\mathrm{d}(r) [V_\mathrm{d}(r)-V(r)]^2  \over [\Sigma_\mathrm{d}(r)+\Sigma_\mathrm{s}(r)]}.
    \end{equation}
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Coordinates             , only : coordinateCylindrical, assignment(=)
    implicit none
    double precision                       , intent(in   ) :: radius
    double precision                       , parameter     :: fractionSmall                         =1.0d-3
    type            (integrator           )                :: integratorDensity                            , integratorVelocityDensity
    double precision                                       :: sigmaLineOfSightSquaredSpheroidDensity       , densitySpheroid          , &
         &                                                    densityDisk                                  , velocityDisk             , &
         &                                                    velocityMean                                 , sigmaLineOfSightSquared
    type            (coordinateCylindrical)                :: coordinates

    if (radius <= 0.0d0) then
       velocityDispersionLambdaRIntegrand1=0.0d0
    else
       radiusImpact_                 = radius
       coordinates                   =[radius,0.0d0,0.0d0]
       integratorDensity             =integrator                                 (velocityDispersionDensityIntegrand,toleranceRelative=1.0d-2)
       densitySpheroid               =integratorDensity           %integrate     (radius                            ,radiusOuter_            )
       densityDisk                   =massDistributionStellarDisk_%surfaceDensity(coordinates                                                )
       velocityDisk                  =massDistributionTotal_      %rotationCurve (radius                                                     )
       ! Test if the spheroid density is significant....
       if (densitySpheroid < fractionSmall*densityDisk) then
          ! ...it is not, so we can avoid computing the spheroid velocity dispersion.
          velocityDispersionLambdaRIntegrand1=+2.0d0         &
               &                               *Pi           &
               &                               *radius       &
               &                               *densityDisk  &
               &                               *velocityDisk
       else
          ! ...it is, so we must do the full calculation.
          integratorVelocityDensity             =integrator                         (velocityDispersionVelocityDensityIntegrand,toleranceRelative=1.0d-2)
          sigmaLineOfSightSquaredSpheroidDensity=integratorVelocityDensity%integrate(radius                                    ,radiusOuter_            )
          velocityMean                       =+densityDisk                              &
               &                              *velocityDisk                             &
               &                              /(densityDisk+densitySpheroid)
          sigmaLineOfSightSquared            =+(                                        &
               &                                +sigmaLineOfSightSquaredSpheroidDensity &
               &                                +densitySpheroid                        &
               &                                *              velocityMean **2         &
               &                                +densityDisk                            &
               &                                *(velocityDisk-velocityMean)**2         &
               &                               )                                        &
               & /(densityDisk+densitySpheroid)
          velocityDispersionLambdaRIntegrand1=+2.0d0                                    &
               &                              *Pi                                       &
               &                              *radius                                   &
               &                              *(densityDisk+densitySpheroid)            &
               &                              *sqrt(                                    &
               &                                    +sigmaLineOfSightSquared            &
               &                                    +velocityMean**2                    &
               &                                   )
       end if
    end if
    return
  end function velocityDispersionLambdaRIntegrand1

  double precision function velocityDispersionLambdaRIntegrand2(radius)
    !!{
    Integrand function used for integrating the $\lambda_\mathrm{R}$ statistic of \cite{cappellari_sauron_2007}. In this case we
    want to evaluate
    \begin{equation}
    \int_0^r 2 \pi r^\prime \Sigma(r^\prime) V(r^\prime) \mathrm{d}r^\prime,
    \end{equation}
    where $\Sigma(r)$ is the projected surface density (in mass or light) of the galaxy at radius $r$, and $V(r)$ the measured
    rotation speed. Assuming that the selected component is purely dispersion dominated with velocity dispersion $\sigma_\mathrm{
    s}(r)$, and that rotation is present in only the disk component with rotation curve $V_\mathrm{d}(r)$ then we can model the
    velocity distribution, $P(V)$, at $r$ as the sum of a Gaussian of width $\sigma_\mathrm{s}(r)$ and normalized area
    $\Sigma_\mathrm{s}(r)$, and a delta function at $V_\mathrm{d}(r)$ with normalized area $\Sigma_\mathrm{d}(r)$. The measured rotation
    speed is then:
    \begin{equation}
    V(r) = \left. \int_{-\infty}^{+\infty} P(V) V \mathrm{d}V \right/ \int_{-\infty}^{+\infty} P(V) \mathrm{d}V = {\Sigma_\mathrm{d}(r) V_\mathrm{d}(r) \over [\Sigma_\mathrm{d}(r)+\Sigma_\mathrm{s}(r)]}.
    \end{equation}
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Coordinates             , only : coordinateCylindrical, assignment(=)
    implicit none
    double precision                       , intent(in   ) :: radius
    double precision                                       :: densityDisk, velocityDisk
    type            (coordinateCylindrical)                :: coordinates

    if (radius <= 0.0d0) then
       velocityDispersionLambdaRIntegrand2=0.0d0
    else
       coordinates                        =[radius,0.0d0,0.0d0]
       densityDisk                        =massDistributionStellarDisk_%surfaceDensity(coordinates)
       velocityDisk                       =massDistributionTotal_      %rotationCurve (radius     )
       velocityDispersionLambdaRIntegrand2=+2.0d0        &
            &                              *Pi           &
            &                              *radius       &
            &                              *densityDisk  &
            &                              *velocityDisk
    end if
    return
  end function velocityDispersionLambdaRIntegrand2

  double precision function velocityDispersionVelocitySurfaceDensityIntegrand(radius)
    !!{
    Integrand function used for integrating line-of-sight velocity dispersion over surface density.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    double precision                     , intent(in   ) :: radius
    type            (coordinateSpherical)                :: coordinates
    
    if (radius <= 0.0d0) then
       velocityDispersionVelocitySurfaceDensityIntegrand=0.0d0
    else
       coordinates                                      =[radius,0.0d0,0.0d0]
       velocityDispersionVelocitySurfaceDensityIntegrand=+                          velocityDispersionSolidAngleInCylinder(radius                                              )    &
            &                                            *                                                                 radius                                               **2 &
            &                                            *massDistributionWeighted_%density                               (coordinates                                         )    &
            &                                            *kinematicsDistribution_  %velocityDispersion1D                  (coordinates,massDistribution_,massDistributionTotal_)**2
    end if
   return
  end function velocityDispersionVelocitySurfaceDensityIntegrand

  double precision function velocityDispersionSurfaceDensityIntegrand(radius)
    !!{
    Integrand function used for integrating line-of-sight surface density dispersion over area.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    double precision                     , intent(in   ) :: radius
    type            (coordinateSpherical)                :: coordinates

    if (radius <= 0.0d0) then
       velocityDispersionSurfaceDensityIntegrand=+0.0d0
    else
       coordinates                              =[radius,0.0d0,0.0d0]
       velocityDispersionSurfaceDensityIntegrand=+                          velocityDispersionSolidAngleInCylinder(radius     )    &
            &                                    *                                                                 radius      **2 &
            &                                    *massDistributionWeighted_%density                               (coordinates)
    end if
    return
  end function velocityDispersionSurfaceDensityIntegrand

  double precision function velocityDispersionSolidAngleInCylinder(radius)
    !!{
    Computes the solid angle of a spherical shell of given {\normalfont \ttfamily radius} that lies within a cylinder of radius {\normalfont \ttfamily
    radiusImpact}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    ! Test size of sphere relative to cylinder.
    if (radius <= radiusImpact_) then
       ! Sphere is entirely within the cylinder.
       velocityDispersionSolidAngleInCylinder=+4.0d0*Pi
    else
       ! Some of sphere lies outside of cylinder.
       velocityDispersionSolidAngleInCylinder=+4.0d0*Pi*(1.0d0-cos(asin(radiusImpact_/radius)))
    end if
    return
  end function velocityDispersionSolidAngleInCylinder

  double precision function velocityDispersionLineOfSightVelocityDispersionIntegrand(radius)
    !!{
    Compute the line-of-sight velocity dispersion at the given {\normalfont \ttfamily radius}.
    !!}
    use :: Error                , only : Error_Report, errorStatusSuccess
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: radius
    type            (integrator)                :: integratorVelocityDensity, integratorDensity
    double precision                            :: densityIntegral          , velocityDensityIntegral
    integer                                     :: statusVelocityDensity    , statusDensity

    integratorVelocityDensity=integrator                         (velocityDispersionVelocityDensityIntegrand,toleranceRelative=self_%toleranceRelative                             )
    integratorDensity        =integrator                         (velocityDispersionDensityIntegrand        ,toleranceRelative=self_%toleranceRelative                             )
    velocityDensityIntegral  =integratorVelocityDensity%integrate(radius                                    ,radiusOuter_                             ,status=statusVelocityDensity)
    densityIntegral          =integratorDensity        %integrate(radius                                    ,radiusOuter_                             ,status=statusDensity        )
    if     (                                                                                                   &
         &   (                                                                                                 &
         &     statusVelocityDensity /= errorStatusSuccess                                                     &
         &    .or.                                                                                             &
         &     statusDensity         /= errorStatusSuccess                                                     &
         &   )                                                                                                 &
         &  .and.                                                                                              &
         &   self_%integrationFailureIsFatal                                                                   &
         & ) call Error_Report('line-of-sight velocity dispersion integral failure'//{introspection:location})
    if (velocityDensityIntegral <= 0.0d0) then
       velocityDispersionLineOfSightVelocityDispersionIntegrand=0.0d0
    else
       velocityDispersionLineOfSightVelocityDispersionIntegrand=sqrt(velocityDensityIntegral/densityIntegral)
    end if
    return
  end function velocityDispersionLineOfSightVelocityDispersionIntegrand

  double precision function velocityDispersionDensityIntegrand(radius)
    !!{
    Integrand function used for computing line-of-sight velocity dispersions.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    if (radius <= radiusImpact_) then
       velocityDispersionDensityIntegrand=+0.0d0
    else
        velocityDispersionDensityIntegrand=+massDistributionWeighted_%densitySphericalAverage(radius) &
            &                             *     radius                                               &
            &                             /sqrt(radius**2-radiusImpact_**2)
    end if
    return
  end function velocityDispersionDensityIntegrand

  double precision function velocityDispersionVelocityDensityIntegrand(radius)
    !!{
    Integrand function used for computing line-of-sight velocity dispersions. Specifically, we wish to evaluate the integral:
    \begin{equation}
    \int_{r_\mathrm{i}}^{r_\mathrm{o}} \sigma^2(r) \rho(r) {r \over \sqrt{r^2-r_\mathrm{i}^2}} \mathrm{d}r,
    \label{eq:velocityDispersionDensityIntegral}
    \end{equation}
    where $r_\mathrm{i}$ is the impact parameter, $r_\mathrm{o}$ is an outer radius at which we assume $\rho(r_\mathrm{
    o})\sigma^2(r_\mathrm{o}) = 0$ (i.e. it is the radius at which we begin integrating the Jeans equation), $\rho(r)$ is density,
    and $\sigma(r)$ is the velocity dispersion at radius $r$. Assuming spherical symmetry and isotropic velocity dispersion,
    the Jeans equation tells us
    \begin{equation}
    \rho(r) \sigma^2(r) = \int^{r_\mathrm{o}}_r {\mathrm{G} M(<r^\prime) \over r^{\prime 2}} \rho(r^\prime) \mathrm{d}r^\prime,
    \label{eq:sphericalIsotropicJeans}
    \end{equation}
    where $\mathrm{G}$ is the gravitational constant, and $M(<r)$ is the total mass contained within radius
    $r$. Equation~(\ref{eq:velocityDispersionDensityIntegral}) can then be simplified using integration by parts to give:
    \begin{equation}
    \left[ \sigma^2(r)\rho(r)\sqrt{r^2-r_\mathrm{i}^2}\right]_{r_\mathrm{i}}^{r_\mathrm{o}} + \int_{r_\mathrm{i}}^{r_\mathrm{o}} {\mathrm{d}\over \mathrm{d}r} \left[ \sigma^2(r) \rho(r) \right] \sqrt{r^2-r_\mathrm{i}^2} \mathrm{d}r.
    \end{equation}
    The first term is zero at both limits (due to the constraint $\rho(r_\mathrm{o})\sigma^2(r_\mathrm{o}) = 0$ at $r_\mathrm{o}$ and
    due to $\sqrt{r^2-r_\mathrm{i}^2}=0$ at $r_\mathrm{i}$), and the second term can be simplified using
    eqn.~(\ref{eq:sphericalIsotropicJeans}) to give
    \begin{equation}
    \int_{r_\mathrm{i}}^{r_\mathrm{o}} {\mathrm{G} M(<r) \over r^2} \rho(r) \sqrt{r^2-r_\mathrm{i}^2} \mathrm{d}r.
    \end{equation}
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    double precision, intent(in   ) :: radius

    if (radius <= radiusImpact_) then
       velocityDispersionVelocityDensityIntegrand=+0.0d0
    else
       velocityDispersionVelocityDensityIntegrand=+gravitationalConstant_internal                            &
            &                                     *massDistributionWeighted_%densitySphericalAverage(radius) &
            &                                     *massDistributionTotal_   %massEnclosedBySphere   (radius) &
            &                                     /     radius**2                                            &
            &                                     *sqrt(radius**2-radiusImpact_**2)
    end if
    return
  end function velocityDispersionVelocityDensityIntegrand
