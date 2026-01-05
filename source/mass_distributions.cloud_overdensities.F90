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
  Implementation of a mass distribution class which overlays clouds on another mass distribution.
  !!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Nearest_Neighbors       , only : nearestNeighbors
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <massDistribution name="massDistributionCloudOverdensities">
   <description>A mass distribution class which overlays clouds on another mass distribution.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionClass) :: massDistributionCloudOverdensities
     !!{
     A mass distribution class which overlays clouds on another mass distribution.
     !!}
     class           (randomNumberGeneratorClass), pointer                     :: randomNumberGenerator_    => null()
     class           (massDistributionClass     ), pointer                     :: massDistribution_         => null()
     double precision                            , allocatable, dimension(:,:) :: positions
     type            (nearestNeighbors          )                              :: neighbors
     double precision                                                          :: radius                             , densityContrast, &
          &                                                                       volumeFillingFactor                , radiusBoundary , &
          &                                                                       densityContrastIntercloud          , radiusSquared
     integer         (c_size_t                  )                              :: countClouds
   contains
     final     ::            cloudOverdensitiesDestructor
     procedure :: density => cloudOverdensitiesDensity
  end type massDistributionCloudOverdensities

  interface massDistributionCloudOverdensities
     !!{
     Constructors for the \refClass{massDistributionCloudOverdensities} mass distribution class.
     !!}
     module procedure cloudOverdensitiesConstructorParameters
     module procedure cloudOverdensitiesConstructorInternal
  end interface massDistributionCloudOverdensities

contains

  function cloudOverdensitiesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionCloudOverdensities} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionCloudOverdensities)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (massDistributionClass             ), pointer       :: massDistribution_
    class           (randomNumberGeneratorClass        ), pointer       :: randomNumberGenerator_
    double precision                                                    :: radius                , densityContrast, &
          &                                                                volumeFillingFactor   , radiusBoundary
    logical                                                             :: dimensionless
    type            (varying_string                    )                :: componentType
    type            (varying_string                    )                :: massType

    !![
    <inputParameter>
      <name>radius</name>
      <description>The cloud radius.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityContrast</name>
      <description>The cloud density contrast.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>volumeFillingFactor</name>
      <description>The cloud volume filling factor.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusBoundary</name>
      <description>The boundary radius within which to populate clouds.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the cloud overdensities profile is considered to be dimensionless.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="massDistribution"      name="massDistribution_"      source="parameters"/>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=massDistributionCloudOverdensities(radius,densityContrast,volumeFillingFactor,radiusBoundary,massDistribution_,randomNumberGenerator_,dimensionless,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <objectDestructor name="massDistribution_"     />
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cloudOverdensitiesConstructorParameters
  
  function cloudOverdensitiesConstructorInternal(radius,densityContrast,volumeFillingFactor,radiusBoundary,massDistribution_,randomNumberGenerator_,dimensionless,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionCloudOverdensities} mass distribution class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Sorting                 , only : sortIndex, sortByIndex
    implicit none
    type            (massDistributionCloudOverdensities)                              :: self
    double precision                                    , intent(in   )               :: radius                , densityContrast, &
          &                                                                              volumeFillingFactor   , radiusBoundary
    class           (massDistributionClass             ), intent(in   ), target       :: massDistribution_
    class           (randomNumberGeneratorClass        ), intent(in   ), target       :: randomNumberGenerator_
    logical                                             , intent(in   ), optional     :: dimensionless
    type            (enumerationComponentTypeType      ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional     :: massType
    integer         (c_size_t                          ), allocatable  , dimension(:) :: order
    double precision                                                                  :: positionRadius        , positionTheta  , &
         &                                                                               positionPhi
    integer         (c_size_t                          )                              :: i
    !![
    <constructorAssign variables="radius, densityContrast, volumeFillingFactor, radiusBoundary, componentType, massType, *massDistribution_, *randomNumberGenerator_"/>
    !!]

    ! Determine if profile is dimensionless.
    if (present(dimensionless)) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! Compute square of cloud radius - used for determining if points are within a cloud.
    self%radiusSquared=self%radius**2
    ! Generate number of clouds.
    self%countClouds=int(                      &
         &               +volumeFillingFactor  &
         &               *(                    &
         &                 +radiusBoundary     &
         &                 /radius             &
         &                )**3               , &
         &               c_size_t              &
         &              )
    ! Determine intercloud density reduction factor.
    self%densityContrastIntercloud=+1.0d0                          &
         &                         /(                              &
         &                           +        volumeFillingFactor  &
         &                           *        densityContrast      &
         &                           +(+1.0d0-volumeFillingFactor) &
         &                          )
    ! Sample cloud positions.
    allocate(self%positions(self%countClouds,3))
    do i=1,self%countClouds
       positionRadius=              +self%radiusBoundary                                               &
            &                       *self%randomNumberGenerator_%uniformSample()**(1.0d0/3.0d0)
       positionTheta =acos(+2.0d0   *self%randomNumberGenerator_%uniformSample()               -1.0d0)
       positionPhi   =     +2.0d0*Pi*self%randomNumberGenerator_%uniformSample()
       self%positions(i,:)= [                                     &
            &                sin(positionTheta)*cos(positionPhi), &
            &                sin(positionTheta)*sin(positionPhi), &
            &                cos(positionTheta)                   &
            &               ]                                     &
            &              *positionRadius
    end do
    ! Sort the clouds by x-coordinate.
    order=sortIndex(self%positions(:,1))
    do i=1,3
       call sortByIndex(self%positions(:,i),order)
    end do
    self%neighbors=nearestNeighbors(self%positions)
    return
  end function cloudOverdensitiesConstructorInternal
  
  subroutine cloudOverdensitiesDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionCloudOverdensities} mass distribution class.
    !!}
    implicit none
    type(massDistributionCloudOverdensities), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"     />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine cloudOverdensitiesDestructor

  double precision function cloudOverdensitiesDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a cloud overdensities mass distribution.
    !!}
    use :: Coordinates  , only : assignment(=), coordinateCartesian
    implicit none
    class           (massDistributionCloudOverdensities), intent(inout) :: self
    class           (coordinate                        ), intent(in   ) :: coordinates
    type            (coordinateCartesian               )                :: position
    double precision                                    , dimension(3)  :: positionComponents
    double precision                                                    :: densityContrast
    integer                                                             :: neighborCount
    
    ! Extract the position.
    position          =coordinates
    positionComponents=position
    ! Determine if this point is within a cloud.
    call self%neighbors%searchFixedRadius(positionComponents,self%radius,0.0d0,neighborCount)    
    ! Determine density contrast.
    if (neighborCount > 0) then
       densityContrast=self%densityContrast*self%densityContrastIntercloud
    else
       densityContrast=                     self%densityContrastIntercloud
    end if
    ! Compute the final density.
    cloudOverdensitiesDensity=+self%massDistribution_%density        (coordinates) &
         &                    *                       densityContrast
    return
  end function cloudOverdensitiesDensity
