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
Implements an output analysis property extractor class that extracts the Einstein radius of the halo.
!!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Root_Finder            , only : rootFinder
  use :: Mass_Distributions     , only : massDistributionClass
  use :: Numerical_Integration  , only : integrator

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusEinstein">
   <description>An output analysis property extractor class that extracts the Einstein radius of the halo.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusEinstein
     !!{
     A property extractor output analysis class that extracts the Einstein radius of the halo.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_       => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_      => null()
     double precision                                    :: redshiftSource                     , timeSource
     type            (rootFinder              )          :: finder
     type            (integrator              )          :: integratorImpactParameter          , integratorLineOfSight
   contains
     final     ::                radiusEinsteinDestructor
     procedure :: extract     => radiusEinsteinExtract
     procedure :: name        => radiusEinsteinName
     procedure :: description => radiusEinsteinDescription
     procedure :: unitsInSI   => radiusEinsteinUnitsInSI
  end type nodePropertyExtractorRadiusEinstein

  interface nodePropertyExtractorRadiusEinstein
     !!{
     Constructors for the {\normalfont \ttfamily radiusEinstein} output analysis class.
     !!}
     module procedure radiusEinsteinConstructorParameters
     module procedure radiusEinsteinConstructorInternal
  end interface nodePropertyExtractorRadiusEinstein

  ! Submodule-scope variables used in root-finding.
  class           (nodePropertyExtractorRadiusEinstein), pointer :: self_
  class           (massDistributionClass              ), pointer :: massDistribution_
  double precision                                               :: densitySurfaceCritical     , radiusImpact_, &
       &                                                            distanceLineOfSightMaximum_
  !$omp threadprivate(self_,massDistribution_,densitySurfaceCritical,radiusImpact_,distanceLineOfSightMaximum_)
  
contains

  function radiusEinsteinConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiusEinstein} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodePropertyExtractorRadiusEinstein)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    double precision                                                     :: redshiftSource
        
     !![
    <inputParameter>
     <name>redshiftSource</name>
     <source>parameters</source>
     <description>The source redshift to using in Einstein radius calculations.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusEinstein(                                                                                                 &
         &                                   cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftSource)), &
         &                                   cosmologyFunctions_                                                                            , &
         &                                   darkMatterHaloScale_                                                                             &
         &                                  )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function radiusEinsteinConstructorParameters

  function radiusEinsteinConstructorInternal(timeSource,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily radiusEinstein} node property extractor.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type            (nodePropertyExtractorRadiusEinstein)                        :: self
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                     , intent(in   )         :: timeSource
    !![
    <constructorAssign variables="timeSource, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]

    ! Compute corresponding redshift - these are needed for the descriptor.
    self%redshiftSource=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeSource))
    ! Build required objects.
    self%finder=rootFinder                   (                                                                                  &
         &                                    rootFunction                 =radiusEinsteinProjectedDensityRoot                , &
         &                                    rangeExpandDownward          =0.5d0                                             , &
         &                                    rangeExpandUpward            =2.0d0                                             , &
         &                                    rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive                     , &
         &                                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative                     , &
         &                                    rangeExpandType              =rangeExpandMultiplicative                         , &
         &                                    toleranceRelative            =1.0d-3                                              &
         &                                   )
    self%integratorImpactParameter=integrator(                                                                                  &
         &                                    integrand                    =radiusEinsteinProjectedDensityIntegrandImpact     , &
         &                                    toleranceRelative            =1.0d-3                                              &
         &                                   )
    self%integratorLineOfSight    =integrator(                                                                                  &
         &                                    integrand                    =radiusEinsteinProjectedDensityIntegrandLineOfSight, &
         &                                    toleranceRelative            =1.0d-3                                              &
         &                                   )
    return
  end function radiusEinsteinConstructorInternal

  subroutine radiusEinsteinDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily radiusEinstein} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusEinstein), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine radiusEinsteinDestructor

  double precision function radiusEinsteinExtract(self,node,instance)
    !!{
    Implement an Einstein radius calculation.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo 
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, arcsecondsToDegrees   , degreesToRadians
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , nodeComponentSatellite
    implicit none
    class           (nodePropertyExtractorRadiusEinstein), intent(inout), target   :: self
    type            (treeNode                           ), intent(inout), target   :: node
    type            (multiCounter                       ), intent(inout), optional :: instance
    class           (nodeComponentBasic                 )               , pointer  :: basic
    class           (nodeComponentSatellite             )               , pointer  :: satellite
    double precision                                     , parameter               :: radiusEinsteinTiny       =+1.0d-3                                    &
         &                                                                                                      *arcsecondstoDegrees                       &
         &                                                                                                      *degreesToRadians
    double precision                                                               :: distanceAngularLensSource                     , distanceAngularLens, &
         &                                                                            distanceAngularSource
    !$GLC attributes unused :: instance

    radiusEinsteinExtract=-1.0d0
    ! Find the angular diameter distances.
    basic => node%basic()
    if (basic%time() <= self%timeSource) return
    distanceAngularLens      =self%cosmologyFunctions_%distanceAngular(                basic%time())
    distanceAngularSource    =self%cosmologyFunctions_%distanceAngular(self%timeSource             )
    distanceAngularLensSource=self%cosmologyFunctions_%distanceAngular(self%timeSource,basic%time())
    if (distanceAngularLens <= 0.0d0 .or. distanceAngularLensSource <= 0.0d0) return
    ! Compute the critical surface density.
    densitySurfaceCritical=+(                              &
         &                   +speedLight                   &
         &                   /kilo                         &
         &                  )**2                           &
         &                 /4.0d0                          &
         &                 /Pi                             &
         &                 /gravitationalConstant_internal &
         &                 *distanceAngularSource          &
         &                 /distanceAngularLens            &
         &                 /distanceAngularLensSource
    ! Find the outer radius of the halo.
    massDistribution_ => node%massDistribution()
    satellite                   => node             %satellite          (                                       )
    distanceLineOfSightMaximum_ =  massDistribution_%radiusEnclosingMass(min(satellite%boundMass(),basic%mass()))
    ! Find the radius within which the mean projected surface density equals the critical density.
    self_ => self
    if (radiusEinsteinProjectedDensityRoot(radiusEinsteinTiny*distanceAngularLens) < 0.0d0) then
       ! For extremely tiny Einstein radii we simply return zero, to avoid unnecessary computation.
       radiusEinsteinExtract=+0.0d0
    else
       radiusEinsteinExtract=+self%finder%find(rootGuess=self%darkMatterHaloScale_%radiusVirial(node)) &
            &                /distanceAngularLens                                                      &
            &                /degreesToRadians                                                         &
            &                /arcsecondstoDegrees
    end if
    return
  end function radiusEinsteinExtract

  double precision function radiusEinsteinProjectedDensityRoot(radius)
    !!{
    Root function used to find the Einstein radius.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    radiusEinsteinProjectedDensityRoot=+self_%integratorImpactParameter%integrate(0.0d0,radius) &
         &                             /Pi                                                      &
         &                             /radius**2                                               &
         &                             -densitySurfaceCritical    
    return
  end function radiusEinsteinProjectedDensityRoot
  
  double precision function radiusEinsteinProjectedDensityIntegrandImpact(radiusImpact)
    !!{
    Integrand function used in finding the mean enclosed projected density.
      !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radiusImpact
    double precision                :: distanceLineOfSightMaximum

    if (radiusImpact < distanceLineOfSightMaximum_) then
       radiusImpact_                                =+radiusImpact    
       distanceLineOfSightMaximum                   =+sqrt(                                &
            &                                              +distanceLineOfSightMaximum_**2 &
            &                                              -radiusImpact               **2 &
            &                                             )
       radiusEinsteinProjectedDensityIntegrandImpact=+4.0d0                                                                   & ! Factor 2π for circumference, factor 2 since
            &                                        *Pi                                                                      & ! we integrate from 0 to the maximum distance,
            &                                        *radiusImpact                                                            & ! and assume symmetry at negative distances,
            &                                        *self_%integratorLineOfSight%integrate(0.0d0,distanceLineOfSightMaximum)
    else
       radiusEinsteinProjectedDensityIntegrandImpact=+0.0d0
    end if
    return
  end function radiusEinsteinProjectedDensityIntegrandImpact
  
  double precision function radiusEinsteinProjectedDensityIntegrandLineOfSight(distance)
    !!{
    Integrand function used in finding the mean enclosed projected density.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll   , massTypeAll
    use :: Coordinates               , only : coordinateSpherical, assignment(=)
    implicit none
    double precision                     , intent(in   ) :: distance
    double precision                                     :: radius
    type            (coordinateSpherical)                :: coordinates

    radius                                            =+sqrt(                  &
         &                                                   +distance     **2 &
         &                                                   +radiusImpact_**2 &
         &                                                  )
    coordinates                                       = [radius,0.0d0,0.0d0]
    radiusEinsteinProjectedDensityIntegrandLineOfSight=+massDistribution_%density(coordinates)
    return
  end function radiusEinsteinProjectedDensityIntegrandLineOfSight
  
  function radiusEinsteinName(self)
    !!{
    Return the name of the Einstein radius property.
    !!}
    implicit none
    type (varying_string                     )                :: radiusEinsteinName
    class(nodePropertyExtractorRadiusEinstein), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusEinsteinName=var_str('radiusEinstein')
    return
  end function radiusEinsteinName

  function radiusEinsteinDescription(self) result(description)
    !!{
    Return a description of the Einstein radius property.
    !!}
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    implicit none
    type     (varying_string                     )                :: description
    class    (nodePropertyExtractorRadiusEinstein), intent(inout) :: self
    character(len=8                              )                :: label
    
    description='The Einstein radius of the node ('
    write (label,'(f7.4)') self%redshiftSource
    description=description//"zₛ = "//trim(adjustl(label))//") [arcseconds]"
    return
  end function radiusEinsteinDescription

  double precision function radiusEinsteinUnitsInSI(self)
    !!{
    Return the units of the Einstein radius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : degreesToRadians, arcsecondsToDegrees
    implicit none
    class(nodePropertyExtractorRadiusEinstein), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusEinsteinUnitsInSI=+degreesToRadians    &
         &                  *arcsecondsToDegrees
    return
  end function radiusEinsteinUnitsInSI
