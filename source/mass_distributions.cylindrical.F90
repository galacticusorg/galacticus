!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of an abstract mass distribution class for cylindrically symmetric distributions.
  !!}

  !![
  <massDistribution name="massDistributionCylindrical" abstract="yes">
   <description>An abstract mass distribution class for cylindrically symmetric distributions.</description>
  </massDistribution>
  !!]
  type, extends(massDistributionClass), abstract :: massDistributionCylindrical
     !!{
     Implementation of an abstract mass distribution class for cylindrically symmetric distributions.
     !!}
     private
   contains
     !![
     <methods>
       <method description="Returns the cylindrical radius enclosing half of the mass of the mass distribution." method="radiusHalfMass"/>
     </methods>
     !!]
     procedure                                      :: symmetry              => cylindricalSymmetry
     procedure                                      :: chandrasekharIntegral => cylindricalChandrasekharIntegral
     procedure(cylindricalRadiusHalfMass), deferred :: radiusHalfMass
  end type massDistributionCylindrical

  abstract interface

     double precision function cylindricalRadiusHalfMass(self,componentType,massType)
       !!{
       Interface for cylindrically symmetric mass distribution half mass radii functions.
       !!}
       import massDistributionCylindrical, enumerationComponentTypeType, enumerationMassTypeType
       class(massDistributionCylindrical ), intent(inout)           :: self
       type (enumerationComponentTypeType), intent(in   ), optional :: componentType
       type (enumerationMassTypeType     ), intent(in   ), optional :: massType
     end function cylindricalRadiusHalfMass

  end interface

contains

  function cylindricalSymmetry(self)
    !!{
    Returns symmetry label for mass distributions with cylindrical symmetry.
    !!}
    implicit none
    type (enumerationMassDistributionSymmetryType)                :: cylindricalSymmetry
    class(massDistributionCylindrical            ), intent(inout) :: self
    !$GLC attributes unused :: self

    cylindricalSymmetry=massDistributionSymmetryCylindrical
    return
  end function cylindricalSymmetry

  function cylindricalChandrasekharIntegral(self,massDistributionEmbedding,massDistributionPerturber,coordinates,velocity,componentType,massType)
    !!{
    Compute the Chandrasekhar integral at the specified {\normalfont \ttfamily coordinates} in a spherical mass distribution.
    !!}
    use :: Coordinates                     , only : coordinateCartesian            , coordinateSpherical, coordinateCylindrical       , assignment(=)
    use :: Galactic_Structure_Options      , only : componentTypeAll               , massTypeAll        , enumerationComponentTypeType, enumerationMassTypeType
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Linear_Algebra                  , only : vector                         , matrix             , assignment(=)


    !! WORKAROUND
    use    :: Display      , only : displayBold       , displayRed     , displayReset



    implicit none


    !! WORKAROUND
    logical :: warned=.false.
    

    double precision                                   , dimension(3)            :: cylindricalChandrasekharIntegral
    class           (massDistributionCylindrical      ), intent(inout)           :: self
    class           (massDistributionClass            ), intent(inout)           :: massDistributionEmbedding                           , massDistributionPerturber
    class           (coordinate                       ), intent(in   )           :: coordinates                                         , velocity
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType
    double precision                                   , dimension(3)            :: velocityCartesian_
    double precision                                   , parameter               :: toomreQRadiusHalfMass                        =1.50d0  ! The Toomre Q-parameter at the disk half-mass radius (Benson et al.,
    ! 2004 , https://ui.adsabs.harvard.edu/abs/2004MNRAS.351.1215B, Appendix A).
    double precision                                   , parameter               :: toomreQFactor                                =3.36d0  ! The factor appearing in the definition of the Toomre Q-parameter for
    ! a stellar disk (Binney & Tremaine, eqn. 6.71).
    double precision                                   , dimension(3)            :: velocityDisk                                         , velocityRelative                , &
         &                                                                          positionCartesianMidplane                                 , &
         &                                                                          positionCylindricalHalfMass                          , positionCartesian
    type            (massDistributionGaussianEllipsoid), save                    :: velocityDistribution
    logical                                            , save                    :: velocityDistributionInitialized              =.false.
    !$omp threadprivate(velocityDistribution,velocityDistributionInitialized)
    type            (coordinateCartesian              )                          :: coordinates_                                         , coordinatesMidplane             , &
         &                                                                          coordinatesMidplaneHalfMass                          , velocityCartesian
    double precision                                                             :: velocityDispersionRadial                             , velocityDispersionAzimuthal     , &
         &                                                                          velocityDispersionVertical                           , velocityCircular                , &
         &                                                                          velocityCircularHalfMassRadius                       , velocityCircularSquaredGradient , &
         &                                                                          velocityCircularSquaredGradientHalfMassRadius        , density                         , &
         &                                                                          densityMidPlane                                      , densitySurface                  , &
         &                                                                          heightScale                                          , radiusMidplane                  , &
         &                                                                          frequencyCircular                                    , frequencyEpicyclic              , &
         &                                                                          frequencyCircularHalfMassRadius                      , frequencyEpicyclicHalfMassRadius, &
         &                                                                          densitySurfaceRadiusHalfMass                         , velocityDispersionRadialHalfMass, &
         &                                                                          velocityDispersionMaximum                            , velocityRelativeMagnitude       , &
         &                                                                          factorSuppressionExtendedMass                        , extentPerturber
    type            (matrix                           )                          :: rotation

    cylindricalChandrasekharIntegral=0.0d0
    if (.not.self%matches(componentType,massType)) return
    coordinates_                                  =   coordinates
    positionCartesian                             =   coordinates_
    positionCartesianMidplane                     =  [     positionCartesian(1),positionCartesian(2),0.0d0]
    positionCylindricalHalfMass                   =  [self%radiusHalfMass   ( ),0.0d0               ,0.0d0]
    coordinatesMidplane                           =   positionCartesianMidplane
    coordinatesMidplaneHalfMass                   =   positionCylindricalHalfMass
    radiusMidplane                                =   coordinatesMidplane%rSpherical()
    velocityCircular                              =   massDistributionEmbedding%rotationCurve        (     radiusMidplane  ,componentType,massType)
    velocityCircularSquaredGradient               =   massDistributionEmbedding%rotationCurveGradient(     radiusMidplane  ,componentType,massType)
    velocityCircularHalfMassRadius                =   massDistributionEmbedding%rotationCurve        (self%radiusHalfMass(),componentType,massType)
    velocityCircularSquaredGradientHalfMassRadius =   massDistributionEmbedding%rotationCurveGradient(self%radiusHalfMass(),componentType,massType)
    velocityDisk                                  = +[positionCartesianMidplane(2),-positionCartesianMidplane(1),0.0d0] &
         &                                          /radiusMidplane                                                     &
         &                                          *velocityCircular
    ! Compute epicyclic frequency.
    frequencyCircular               =velocityCircular              /     radiusMidplane
    frequencyCircularHalfMassRadius =velocityCircularHalfMassRadius/self%radiusHalfMass()
    frequencyEpicyclic              =sqrt(velocityCircularSquaredGradient              /     radiusMidplane  +2.0d0*frequencyCircular              **2)
    frequencyEpicyclicHalfMassRadius=sqrt(velocityCircularSquaredGradientHalfMassRadius/self%radiusHalfMass()+2.0d0*frequencyCircularHalfMassRadius**2)
    ! Get disk structural properties.
    density                     =+self%density       (coordinates                ,componentType,massType)
    densityMidPlane             =+self%density       (coordinatesMidplane        ,componentType,massType)
    densitySurface              =+self%surfaceDensity(coordinatesMidplane        ,componentType,massType)
    densitySurfaceRadiusHalfMass=+self%surfaceDensity(coordinatesMidplaneHalfMass,componentType,massType)
    if (density <= 0.0d0) return
    heightScale                 =+0.5d0           &
         &                       *densitySurface  &
         &                       /densityMidPlane
    ! Compute normalization of the radial velocity dispersion.
    velocityDispersionRadialHalfMass=+toomreQFactor                    &
         &                           *gravitationalConstantGalacticus  &
         &                           *densitySurfaceRadiusHalfMass     &
         &                           *toomreQRadiusHalfMass            &
         &                           /frequencyEpicyclicHalfMassRadius
    ! Find the velocity dispersion components of the disk.
    velocityDispersionRadial   =+velocityDispersionRadialHalfMass      &
         &                      *sqrt(                                 &
         &                            +densitySurface                  &
         &                            /densitySurfaceRadiusHalfMass    &
         &                           )
    velocityDispersionAzimuthal=+velocityDispersionRadial*frequencyEpicyclic/2.0d0/frequencyCircular
    velocityDispersionVertical =+sqrt(Pi*gravitationalConstantGalacticus*densitySurface*heightScale)
    velocityDispersionMaximum  =+maxval([velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical])
    velocityDispersionRadial   =+velocityDispersionRadial   /velocityDispersionMaximum
    velocityDispersionAzimuthal=+velocityDispersionAzimuthal/velocityDispersionMaximum
    velocityDispersionVertical =+velocityDispersionVertical /velocityDispersionMaximum
    if (any([velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical] <= 0.0d0)) return
    ! Find the relative velocity of the perturber and the disk.
    velocityCartesian          = velocity
    velocityCartesian_         = velocityCartesian
    velocityRelative           =(velocityCartesian_-velocityDisk)/velocityDispersionMaximum
    ! Handle limiting case of large relative velocity.
    velocityRelativeMagnitude  =sqrt(sum(velocityRelative**2))
    ! Initialize the velocity distribution.
    rotation=reshape(                                                                               &
         &            [                                                                             &
         &             +positionCartesianMidplane(1),-positionCartesianMidplane(2),+0.0d0         , &
         &             +positionCartesianMidplane(2),+positionCartesianMidplane(1),+0.0d0         , &
         &             +0.0d0                       ,+0.0d0                       ,+radiusMidplane  &
         &            ]                                                                             &
         &           /radiusMidplane                                                              , &
         &           [3,3]                                                                          &
         &          )
    coordinates_=velocityRelative
    if (.not.velocityDistributionInitialized) then
       velocityDistribution           =massDistributionGaussianEllipsoid(scaleLength=[1.0d0,1.0d0,1.0d0],rotation=rotation,mass=1.0d0,dimensionless=.true.)
       velocityDistributionInitialized=.true.
    end if
    call velocityDistribution%initialize(scaleLength=[velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical],rotation=rotation)
    ! Compute suppression factor due to satellite being an extended mass distribution. This is largely untested - it is meant to
    ! simply avoid extremely large accelerations for subhalo close to the disk plane when that subhalo is much more extended than
    ! the disk.

    !! WORKAROUND - this is disabled until we have all components implementing massDistributions (as we need to include the dark matter here)
    ! extentPerturber=massDistributionPerturber%radiusEnclosingMass(                                 &
    !      &                                                        massFractional=0.5d0           , &
    !      &                                                        componentType =componentTypeAll, &
    !      &                                                        massType      =massTypeAll       &
    !      &                                                       )    
    extentPerturber=1.0d-2
    if (.not.warned) then
       write (0,*) displayRed()//displayBold()//'WORKAROUND - can not compute extent of perturber for Chandrasekhar integral'//displayReset()
       call sleep(10)
       warned=.true.
    end if
    
    if (extentPerturber > heightScale) then
       factorSuppressionExtendedMass=+heightScale     &
            &                        /extentPerturber
    else
       factorSuppressionExtendedMass=+1.0d0
    end if
    ! Evaluate the integral.
    cylindricalChandrasekharIntegral=+density                                                     &
         &                           *velocityDistribution         %acceleration(coordinates_)    &
         &                           /velocityDispersionMaximum                               **2 &
         &                           *factorSuppressionExtendedMass
    return
  end function cylindricalChandrasekharIntegral
