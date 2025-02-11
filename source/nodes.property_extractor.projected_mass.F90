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
  Contains a module which implements a property extractor class for the projected density at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale, darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier

  !![
  <nodePropertyExtractor name="nodePropertyExtractorProjectedMass">
   <description>
    A property extractor class for the projected mass at a set of radii. The radii and types of projected mass to output
    is specified by the {\normalfont \ttfamily radiusSpecifiers} parameter. This parameter's value can contain multiple
    entries, each of which should be a valid
    \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#sec.radiusSpecifiers}{radius
    specifier}.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorProjectedMass
     !!{
     A property extractor class for the projected density at a set of radii.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_          => null()
     integer                                                      :: radiiCount                             , elementCount_
     logical                                                      :: includeRadii
     type   (varying_string          ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusSpecifier         ), allocatable, dimension(:) :: radii
     logical                                                      :: darkMatterScaleRadiusIsNeeded          , diskIsNeeded        , &
          &                                                          spheroidIsNeeded                       , virialRadiusIsNeeded, &
          &                                                          nuclearStarClusterIsNeeded             , satelliteIsNeeded
   contains
     final     ::                       projectedMassDestructor
     procedure :: columnDescriptions => projectedMassColumnDescriptions
     procedure :: size               => projectedMassSize
     procedure :: elementCount       => projectedMassElementCount
     procedure :: extract            => projectedMassExtract
     procedure :: names              => projectedMassNames
     procedure :: descriptions       => projectedMassDescriptions
     procedure :: unitsInSI          => projectedMassUnitsInSI
  end type nodePropertyExtractorProjectedMass

  interface nodePropertyExtractorProjectedMass
     !!{
     Constructors for the ``projectedMass'' output analysis class.
     !!}
     module procedure projectedMassConstructorParameters
     module procedure projectedMassConstructorInternal
  end interface nodePropertyExtractorProjectedMass

  ! Module-scope variables used in integrands.
  double precision :: radius_
  !$omp threadprivate(radius_)

contains

  function projectedMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily projectedMass} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorProjectedMass)                              :: self
    type   (inputParameters                   ), intent(inout)               :: parameters
    type   (varying_string                    ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass          ), pointer                     :: darkMatterHaloScale_
    logical                                                                  :: includeRadii

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter>
      <name>radiusSpecifiers</name>
      <description>A list of radius specifiers at which to output the projected mass profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the radii at which projected mass data are output should also be included in the output file.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorProjectedMass(radiusSpecifiers,includeRadii,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function projectedMassConstructorParameters

  function projectedMassConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily projectedMass} property extractor class.
    !!}
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorProjectedMass)                              :: self
    type   (varying_string                    ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass          ), intent(in   ), target       :: darkMatterHaloScale_
    logical                                    , intent(in   )               :: includeRadii
    !![
    <constructorAssign variables="radiusSpecifiers, includeRadii, *darkMatterHaloScale_"/>
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
         &                                          self%diskIsNeeded                 , &
         &                                          self%spheroidIsNeeded             , &
         &                                          self%nuclearStarClusterIsNeeded   , &
         &                                          self%satelliteIsNeeded            , &
         &                                          self%virialRadiusIsNeeded         , &
         &                                          self%darkMatterScaleRadiusIsNeeded  &
         &                                         )
    return
  end function projectedMassConstructorInternal

  subroutine projectedMassDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily projectedMass} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorProjectedMass), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine projectedMassDestructor

  integer function projectedMassElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily projectedMass} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedMass), intent(inout) :: self
    double precision                                    , intent(in   ) :: time
    !$GLC attributes unused :: time

    projectedMassElementCount=self%elementCount_
    return
  end function projectedMassElementCount

  function projectedMassSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily projectedMass} property extractors.
    !!}
    implicit none
    integer         (c_size_t                          )                :: projectedMassSize
    class           (nodePropertyExtractorProjectedMass), intent(inout) :: self
    double precision                                    , intent(in   ) :: time
    !$GLC attributes unused :: time

    projectedMassSize=self%radiiCount
    return
  end function projectedMassSize

  function projectedMassExtract(self,node,time,instance) result(massProjected)
    !!{
    Implement a {\normalfont \ttfamily projectedMass} property extractor.
    !!}
    use :: Galactic_Structure_Options          , only : componentTypeAll                          , massTypeGalactic                  , massTypeStellar
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius           , radiusTypeDiskHalfMassRadius      , radiusTypeDiskRadius            , radiusTypeGalacticLightFraction, &
          &                                             radiusTypeGalacticMassFraction            , radiusTypeRadius                  , radiusTypeSpheroidHalfMassRadius, radiusTypeSpheroidRadius       , &
          &                                             radiusTypeNuclearStarClusterHalfMassRadius, radiusTypeNuclearStarClusterRadius, radiusTypeStellarMassFraction   , radiusTypeVirialRadius         
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile            , nodeComponentDisk                 , nodeComponentSpheroid           , nodeComponentNSC               , &
          &                                             treeNode
    use :: Numerical_Integration               , only : integrator, GSL_Integ_Gauss15
    use :: Numerical_Comparison                , only : Values_Agree
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Error                               , only : Error_Report
    implicit none
    double precision                                    , dimension(:,:), allocatable :: massProjected
    class           (nodePropertyExtractorProjectedMass), intent(inout) , target      :: self
    type            (treeNode                          ), intent(inout) , target      :: node
    double precision                                    , intent(in   )               :: time
    type            (multiCounter                      ), intent(inout) , optional    :: instance
    class           (nodeComponentDisk                 ), pointer                     :: disk
    class           (nodeComponentSpheroid             ), pointer                     :: spheroid
    class           (nodeComponentNSC                  ), pointer                     :: nuclearStarCluster
    class           (nodeComponentDarkMatterProfile    ), pointer                     :: darkMatterProfile
    class           (massDistributionClass             ), pointer                     :: massDistribution_
    double precision                                    , parameter                   :: toleranceRelative   =1.0d-2
    type            (integrator                        )                              :: integrator_
    integer                                                                           :: i
    double precision                                                                  :: radiusVirial               , radiusOuter          , &
         &                                                                               massProjectedCurrent       , massProjectedPrevious
    logical                                                                           :: converged
    !$GLC attributes unused :: time, instance

    allocate(massProjected(self%radiiCount,self%elementCount_))
    radiusVirial                                               =  self%darkMatterHaloScale_%radiusVirial(node                    )
    if (self%                 diskIsNeeded) disk               =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid           =>                                        node%spheroid         ()
    if (self%   nuclearStarClusterIsNeeded) nuclearStarCluster =>                                        node%NSC              ()
    if (self%darkMatterScaleRadiusIsNeeded) darkMatterProfile  =>                                        node%darkMatterProfile()
    integrator_=integrator(projectedMassIntegrand,toleranceRelative=1.0d-3,hasSingularities=.true.,integrationRule=GSL_Integ_Gauss15)
    do i=1,self%radiiCount
       radius_=self%radii(i)%value
       select case (self%radii(i)%type%ID)
       case   (radiusTypeRadius                          %ID)
          ! Nothing to do.
       case   (radiusTypeVirialRadius                    %ID)
          radius_=+radius_*radiusVirial
       case   (radiusTypeDarkMatterScaleRadius           %ID)
          radius_=+radius_*darkMatterProfile %         scale()
       case   (radiusTypeDiskRadius                      %ID)
          radius_=+radius_*disk              %        radius()
       case   (radiusTypeSpheroidRadius                  %ID)
          radius_=+radius_*spheroid          %        radius()
       case   (radiusTypeNuclearStarClusterRadius        %ID)
          radius_=+radius_*nuclearStarCluster%        radius() 
       case   (radiusTypeDiskHalfMassRadius              %ID)
          radius_=+radius_*disk              %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius          %ID)
          radius_=+radius_*spheroid          %halfMassRadius()
       case   (radiusTypeNuclearStarClusterHalfMassRadius%ID)
          radius_=+radius_*nuclearStarCluster%halfMassRadius()
       case   (radiusTypeGalacticMassFraction            %ID,  &
            &  radiusTypeGalacticLightFraction           %ID)
          massDistribution_ =>  node             %massDistribution   (                                                &
               &                                                      massType      =              massTypeStellar ,  &
               &                                                      componentType =              componentTypeAll,  &
               &                                                      weightBy      =self%radii(i)%weightBy        ,  &
               &                                                      weightIndex   =self%radii(i)%weightByIndex      &
               &                                                     )
          radius_           =  +radius_                                                                               &
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
          radius_           =  +radius_                                                                               &
               &               *massDistribution_%radiusEnclosingMass(                                                &
               &                                                      massFractional=self%radii(i)%fraction           &
               &                                                     )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       case default
          call Error_Report('unrecognized radius type'//{introspection:location})
       end select
       massProjectedPrevious=0.0d0
       radiusOuter          =max(radius_*2.0d0,radiusVirial)
       ! Evaluate the integral, then add on the mass of the sphere entirely enclosed inside the cylinder.
       massDistribution_ => node%massDistribution(componentType=self%radii(i)%component,massType=self%radii(i)%mass)
       converged         =  .false.
       do while (.not.converged)
          massProjectedCurrent=integrator_%integrate(log(radius_),log(radiusOuter))
          converged           =Values_Agree(massProjectedCurrent,massProjectedPrevious,relTol=toleranceRelative)
          if (.not.converged) then
             radiusOuter          =2.0d0*radiusOuter
             massProjectedPrevious=      massProjectedCurrent
          end if
       end do
       massProjected(i,1)=+                  massProjectedCurrent          &
            &             +massDistribution_%massEnclosedBySphere(radius_)
       if (self%includeRadii) massProjected(i,2)=radius_
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end do
    return

  contains

    double precision function projectedMassIntegrand(radiusLogarithmic)
      !!{
      Integrand function used for computing projected masses.
      !!}
      use :: Coordinates             , only : coordinateSpherical, assignment(=)
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision                     , intent(in   ) :: radiusLogarithmic
      double precision                                     :: radius
      type            (coordinateSpherical)                :: coordinates

      radius=exp(radiusLogarithmic)
      if (radius <= radius_) then
         projectedMassIntegrand=+0.0d0
      else
         coordinates           =[radius,0.0d0,0.0d0]
         projectedMassIntegrand=+4.0d0                                  & ! ⎫
              &                 *Pi                                     & ! ⎬ Surface area of the spherical shell.
              &                 *radius**2                              & ! ⎭
              &                 *(                                      & ! ⎫
              &                   +1.0d0                                & ! ⎪
              &                   -sqrt(                                & ! ⎪
              &                         +1.0d0                          & ! ⎪
              &                         -(                              & ! ⎪
              &                           +radius_                      & ! ⎬ Fraction of shell solid angle lying
              &                           /radius                       & ! ⎪ inside the cylinder
              &                          )**2                           & ! ⎪
              &                        )                                & ! ⎪
              &                  )                                      & ! ⎭
              &                 *massDistribution_%density(coordinates) & ! } Density of the spherical shell.
              &                 *radius                                   ! } Account for logarithmic integration variable. 
      end if
      return
    end function projectedMassIntegrand

  end function projectedMassExtract

  subroutine projectedMassNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily projectedMass} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedMass), intent(inout)                             :: self
    double precision                                    , intent(in   ), optional                   :: time
    type            (varying_string                    ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names       (1)="projectedMass"
    if (self%includeRadii)                &
         & names(2)="projectedMassRadius"
    return
  end subroutine projectedMassNames

  subroutine projectedMassDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily projectedMass} property.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedMass), intent(inout)                             :: self
    double precision                                    , intent(in   ), optional                   :: time
    type            (varying_string                    ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions       (1)="Projected mass at a given radius [M☉]."
    if (self%includeRadii)                                                      &
         & descriptions(2)="Radius at which projected density is output [Mpc]."
    return
  end subroutine projectedMassDescriptions

  subroutine projectedMassColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily projectedMass} property.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedMass), intent(inout)                            :: self
    double precision                                    , intent(in   ), optional                  :: time
    type            (varying_string                    ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                    , intent(inout), dimension(:), allocatable :: values
    type            (varying_string                    ), intent(  out)                            :: valuesDescription
    double precision                                    , intent(  out)                            :: valuesUnitsInSI

    !$GLC attributes unused :: time

    allocate(descriptions(self%radiiCount))
    allocate(values      (              0))
    valuesDescription=var_str('')
    valuesUnitsInSI  =0.0d0
    descriptions     =self%radii%name
    return
  end subroutine projectedMassColumnDescriptions

  function projectedMassUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily projectedMass} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                    , allocatable  , dimension(:) :: projectedMassUnitsInSI
    class           (nodePropertyExtractorProjectedMass), intent(inout)               :: self
    double precision                                    , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(projectedMassUnitsInSI(self%elementCount_))
    projectedMassUnitsInSI       (1)=massSolar
    if (self%includeRadii)                      &
         & projectedMassUnitsInSI(2)=megaParsec
    return
  end function projectedMassUnitsInSI

