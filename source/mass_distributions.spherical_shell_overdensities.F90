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
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <massDistribution name="massDistributionSphericalShellOverdensities">
   <description>
    A mass distribution class which overlays overdense spherical shells on another mass distribution.
    
    The intent is to mimic the effects of a 3-D distribution of spherical clouds, but along a single sight-line from the center of a
    spherically symmetric mass distribution. This is useful in computing radiative transfer through cloudy media for spherically
    symmetric systems.
    
    In the 3-D case clouds are defined by a radius, $r_\mathrm{c}$, a volume filling factor, $f_\mathrm{v}$, and a density contrast,
    $\Delta_\mathrm{c}$. For this case of spherical shells the same quantities are used, except that the radius is referred to as
    the ``half-width'' of the shell, but is still labeled $r_\mathrm{c}$.
    
    In the 3-D case the number density of clouds is
    \begin{equation}
     n_\mathrm{c} = {f_\mathrm{v} \over (4 \pi / 3 ) r_\mathrm{c}^3}.
    \end{equation}
    Along a sightline of length $l$ (specified by the {\normalfont \ttfamily [radiusBoundary]} parameter) the number of clouds
    intersected is
    \begin{equation}
     N_\mathrm{c} = n_\mathrm{c} l 4 \pi r_\mathrm{c}^2 = 3 f_\mathrm{v} {l \over r_\mathrm{c}}.
    \end{equation}
    This last relation is used to determine the number of spherical shells to generate. These shells are then placed randomly
    in radius between $0$ and $l$. Each shell is also assigned an impact parameter, $b$, meant to represent the distance of the
    center of the notional spherical cloud from the line of sight. The effective half-width of the shell is then
    $\sqrt{r_\mathrm{c}^2-b^2}$.
   </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionClass) :: massDistributionSphericalShellOverdensities
     !!{
     A mass distribution class which overlays clouds on another mass distribution.
     !!}
     class           (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
     class           (massDistributionClass     ), pointer                   :: massDistribution_      => null()
     double precision                            , allocatable, dimension(:) :: radii                           , impactParameter
     double precision                                                        :: halfWidth                       , densityContrast          , &
          &                                                                     volumeFillingFactor             , densityContrastIntershell, &
          &                                                                     radiusBoundary
     integer         (c_size_t                  )                            :: countShells
   contains
     final     ::            sphericalShellOverdensitiesDestructor
     procedure :: density => sphericalShellOverdensitiesDensity
  end type massDistributionSphericalShellOverdensities

  interface massDistributionSphericalShellOverdensities
     !!{
     Constructors for the \refClass{massDistributionSphericalShellOverdensities} mass distribution class.
     !!}
     module procedure sphericalShellOverdensitiesConstructorParameters
     module procedure sphericalShellOverdensitiesConstructorInternal
  end interface massDistributionSphericalShellOverdensities

contains

  function sphericalShellOverdensitiesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalShellOverdensities} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalShellOverdensities)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (massDistributionClass                      ), pointer       :: massDistribution_
    class           (randomNumberGeneratorClass                 ), pointer       :: randomNumberGenerator_
    double precision                                                             :: halfWidth             , densityContrast, &
          &                                                                         volumeFillingFactor   , radiusBoundary
    logical                                                                      :: dimensionless
    type            (varying_string                             )                :: componentType
    type            (varying_string                             )                :: massType

    !![
    <inputParameter>
      <name>halfWidth</name>
      <description>The shell half-width.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityContrast</name>
      <description>The shell density contrast.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>volumeFillingFactor</name>
      <description>The shell volume filling factor.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusBoundary</name>
      <description>The boundary radius within which to populate shells.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the shell overdensities profile is considered to be dimensionless.</description>
      <source>parameters</source>
    </inputParameter>
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
    <objectBuilder class="massDistribution"      name="massDistribution_"      source="parameters"/>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=massDistributionSphericalShellOverdensities(halfWidth,densityContrast,volumeFillingFactor,radiusBoundary,massDistribution_,randomNumberGenerator_,dimensionless,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <objectDestructor name="massDistribution_"     />
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sphericalShellOverdensitiesConstructorParameters
  
  function sphericalShellOverdensitiesConstructorInternal(halfWidth,densityContrast,volumeFillingFactor,radiusBoundary,massDistribution_,randomNumberGenerator_,dimensionless,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalShellOverdensities} mass distribution class.
    !!}
    use :: Sorting      , only : sort
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf, mpiBarrier
#endif    
    implicit none
    type            (massDistributionSphericalShellOverdensities)                          :: self
    double precision                                             , intent(in   )           :: halfWidth             , densityContrast, &
          &                                                                                   volumeFillingFactor   , radiusBoundary
    class           (massDistributionClass                      ), intent(in   ), target   :: massDistribution_
    class           (randomNumberGeneratorClass                 ), intent(in   ), target   :: randomNumberGenerator_
    logical                                                      , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType               ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                    ), intent(in   ), optional :: massType
    integer         (c_size_t                                   )                          :: i
    double precision                                                                       :: countShellsMean
    !![
    <constructorAssign variables="halfWidth, densityContrast, volumeFillingFactor, radiusBoundary, componentType, massType, *massDistribution_, *randomNumberGenerator_"/>
    !!]

    ! Determine if profile is dimensionless.
    if (present(dimensionless)) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! Generate number of shells. If running MPI parallelized, generate the number of shells on the master process, and broadcast
    ! this to other processes, to ensure that they all have identical numbers.
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       countShellsMean=+3.0d0                &
            &          *volumeFillingFactor  &
            &          *radiusBoundary       &
            &          /halfWidth
       self%countShells=self%randomNumberGenerator_%poissonSample(countShellsMean)
#ifdef USEMPI
    end if
    call mpiSelf%broadcastData(0,self%countShells)
    call mpiBarrier()
#endif
    ! Determine intershell density reduction factor.
    self%densityContrastIntershell=+1.0d0                          &
         &                         /(                              &
         &                           +        volumeFillingFactor  &
         &                           *        densityContrast      &
         &                           +(+1.0d0-volumeFillingFactor) &
         &                          )
    ! Sample shell radii. If running MPI parallelized, generate the shell radii on the master process, and broadcast this to other
    ! processes, to ensure that they all have identical radii.
    allocate(self%radii          (self%countShells))
    allocate(self%impactParameter(self%countShells))
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       do i=1,self%countShells
          self%radii          (i)=+     self%radiusBoundary                          &
               &                  *     self%randomNumberGenerator_%uniformSample()
          self%impactParameter(i)=+     self%halfWidth                               &
               &                  *sqrt(self%randomNumberGenerator_%uniformSample())
       end do
       ! Sort the clouds by radial coordinate. (We don't need to sort impact parameters as they're independent of the radii.)
       call sort(self%radii)
#ifdef USEMPI
    end if
    call mpiSelf%broadcastData(0,self%radii          )
    call mpiBarrier()
    call mpiSelf%broadcastData(0,self%impactParameter)
    call mpiBarrier()
#endif
    return
  end function sphericalShellOverdensitiesConstructorInternal
  
  subroutine sphericalShellOverdensitiesDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionSphericalShellOverdensities} mass distribution class.
    !!}
    implicit none
    type(massDistributionSphericalShellOverdensities), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"     />
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine sphericalShellOverdensitiesDestructor

  double precision function sphericalShellOverdensitiesDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a cloud overdensities mass distribution.
    !!}
    use :: Arrays_Search, only : searchArrayClosest
    use :: Coordinates  , only : assignment(=)     , coordinateSpherical
    implicit none
    class           (massDistributionSphericalShellOverdensities), intent(inout) :: self
    class           (coordinate                                 ), intent(in   ) :: coordinates
    integer         (c_size_t                                   )                :: i
    type            (coordinateSpherical                        )                :: position
    double precision                                                             :: densityContrast, radius

    ! Extract the position.
    position=coordinates
    radius  =position   %r()
    ! Determine if this point is within a cloud.
    i=searchArrayClosest(self%radii,radius)
    ! Determine density contrast.
    if (abs(radius-self%radii(i)) <= sqrt(self%halfWidth**2-self%impactParameter(i)**2)) then
       densityContrast=self%densityContrast*self%densityContrastIntershell
    else
       densityContrast=                     self%densityContrastIntershell
    end if
    ! Compute the final density.
    sphericalShellOverdensitiesDensity=+self%massDistribution_%density        (coordinates) &
         &                             *                       densityContrast
    return
  end function sphericalShellOverdensitiesDensity
