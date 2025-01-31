!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implementation of a SIDM Parametric profile mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionSIDMParametricProfile">
   <description>An mass distribution class for SIDM Parametric profile distributions.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionSIDMParametricProfile
     !!{
     The SIDM PArametric profile: $\rho(r)=$
     !!}
     double precision :: beta                  , scaleRadius          , coreRadius           , densityNormalization  , &
          &              momentRadial2Previous , momentRadial3Previous, momentRadial2XPrevious, &
          &              momentRadial3XPrevious, outerRadius
   contains
     procedure :: density               => SIDMParametricProfileDensity
     procedure :: densityGradientRadial => SIDMParametricProfileDensityGradientRadial
!     procedure :: densityRadialMoment   => betaProfileDensityRadialMoment
!     procedure :: densitySquareIntegral => betaProfileDensitySquareIntegral
!     procedure :: massEnclosedBySphere  => betaProfileMassEnclosedBySphere
!     procedure :: potentialIsAnalytic   => betaProfilePotentialIsAnalytic
!     procedure :: potential             => betaProfilePotential
!     procedure :: descriptor            => betaProfileDescriptor
  end type massDistributionSIDMParametricProfile

  interface massDistributionSIDMParametricProfile
     !!{
     Constructors for the {\normalfont \ttfamily SIDMParametricProfile} mass distribution class.
     !!}
     module procedure SIDMParametricProfileConstructorParameters
     module procedure SIDMParametricProfileConstructorInternal
  end interface massDistributionSIDMParametricProfile

contains

  function SIDMParametricProfileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily SIDMParametricProfile} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSIDMParametricProfile)      :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: beta         , densityNormalization , &
         &                                                          mass         , scaleRadius          , &
         &                                                          coreRadius
    logical                                                      :: dimensionless, truncateAtOuterRadius
    type            (varying_string             )                :: componentType
    type            (varying_string             )                :: massType

    !![
    <inputParameter>
      <name>beta</name>
      <defaultValue>4.0d0</defaultValue>
      <description>The value $\beta$ in a SIDMParametric-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The density normalization of a SIDMParametric-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleRadius</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The scale of a SIDMParametric-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>coreRadius</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The core radius of a SIDMParametric-model mass distribution.</description>
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
    !!]
    self=massDistributionSIDMParametricProfile(beta,densityNormalization,scaleRadius,coreRadius,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.))

    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function SIDMParametricProfileConstructorParameters

  function SIDMParametricProfileConstructorInternal(beta,densityNormalization,scaleRadius,coreRadius,componentType,massType) result(self)
    !!{
    Constructor for ``SIDMParametricProfile'' convergence class.
    !!}
    use :: Display                 , only : displayIndent      , displayMessage, displayUnindent, displayVerbosity, &
          &                                 verbosityLevelDebug
    use :: Error                   , only : Error_Report
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Comparison    , only : Values_Agree       , Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionSIDMParametricProfile )                :: self
    double precision                              , intent(in   )           :: beta
    double precision                              , intent(in   ), optional :: densityNormalization              , &
         &                                                                     scaleRadius                       , coreRadius
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="beta, densityNormalization, scaleRadius, coreRadius, componentType, massType"/>
    !!]

    return
  end function SIDMParametricProfileConstructorInternal

  double precision function SIDMParametricProfileDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a SIDMParametric-profile mass distribution.
    !!}
    implicit none
    class           (massDistributionSIDMParametricProfile ), intent(inout) :: self
    class           (coordinate                  ), intent(in   ) :: coordinates
    double precision                                              :: r, term1, term2

!    if (.not.self%matches(componentType,massType)) then
!       SIDMParametricProfileDensity=0.0d0
!       return
!    end if
    
    r=coordinates%rSpherical()
    ! Compute density.
!    print *,'r, self%scaleRadius: ', r, self%scaleRadius
    term1 = ((r ** self%beta + self%coreRadius ** self%beta) ** (1.0d0 / self%beta)) / self%scaleRadius
    term2 = (1.0d0 + r / self%scaleRadius) ** 2

    SIDMParametricProfileDensity=self%densityNormalization/(term1 * term2)

    return
  end function SIDMParametricProfileDensity

  double precision function SIDMParametricProfileDensityGradientRadial(self,coordinates,logarithmic)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a SIDMParametric-profile mass distribution.
    !!}
    implicit none
    class           (massDistributionSIDMParametricProfile ), intent(inout), target   :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    logical                                       , intent(in   ), optional :: logarithmic
    double precision                                                        :: r, term1, term2, term3

    ! Set default options.
!    logarithmicActual=.false.
!    if (present(logarithmic)) logarithmicActual=logarithmic
    ! Get position in spherical coordinate system.
    r=coordinates%rSpherical()

    ! Compute density gradient.
    term1 = (r ** self%beta + self%coreRadius ** self%beta) ** ((-1.0d0-self%beta)/self%beta)
    term2 = 3.0d0 * r**(1.0d0+self%beta) + 2.0d0 * r * self%coreRadius**self%beta + self%scaleRadius * r**self%beta
    term3 = r * (r+self%scaleRadius)**3.0d0

!    print *,'r, self%scaleRadius, term3, self%coreRadius: ', r, self%scaleRadius, term3, self%coreRadius
    if (r>0.0d0) then
       SIDMParametricProfileDensityGradientRadial= (-1.0d0) * term1 * term2 * self%scaleRadius**3.0d0 *self%densityNormalization / term3
    else
!       print *,'correct gradient is used!'
       SIDMParametricProfileDensityGradientRadial= 0.0d0
!       SIDMParametricProfileDensityGradientRadial= (-2.0d0)/(self%coreRadius * self%scaleRadius**3)
    end if

    return
  end function SIDMParametricProfileDensityGradientRadial

