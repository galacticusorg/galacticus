!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% Implementation of a \cite{hernquist_analytical_1990} mass distribution class.

  !# <massDistribution name="massDistributionHernquist">
  !#  <description>A \cite{hernquist_analytical_1990} mass distribution class.</description>
  !# </massDistribution>
  type, public, extends(massDistributionSpherical) :: massDistributionHernquist
     !% The Hernquist \citep{hernquist_analytical_1990} mass distribution.
     private
     double precision :: densityNormalization, mass, &
          &              scaleLength
   contains
     procedure :: density              => hernquistDensity
     procedure :: densityRadialMoment  => hernquistDensityRadialMoment
     procedure :: massEnclosedBySphere => hernquistMassEnclosedBySphere
     procedure :: potential            => hernquistPotential
     procedure :: radiusHalfMass       => hernquistRadiusHalfMass
  end type massDistributionHernquist

  interface massDistributionHernquist
     !% Constructors for the {\normalfont \ttfamily hernquist} mass distribution class.
     module procedure hernquistConstructorParameters
     module procedure hernquistConstructorInternal
  end interface massDistributionHernquist

contains

  function hernquistConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily hernquist} mass distribution class which builds the object from a parameter
    !% set.
    use :: Input_Parameters        , only : inputParameter, inputParameters
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionHernquist)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: mass                , scaleLength, &
         &                                                        densityNormalization
    logical                                                    :: dimensionless

    !# <inputParameter>
    !#   <name>densityNormalization</name>
    !#   <defaultValue>0.5d0/Pi</defaultValue>
    !#   <description>The density normalization of the Hernquist profile.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>scaleLength</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The scale radius of the Hernquist profile.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>mass</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The mass of the Hernquist profile.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>dimensionless</name>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If true the Hernquist profile is considered to be dimensionless.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <conditionalCall>
    !#  <call>self=massDistributionHernquist({conditions})</call>
    !#  <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
    !#  <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
    !#  <argument name="scaleLength"          value="scaleLength"          parameterPresent="parameters"/>
    !#  <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    return
  end function hernquistConstructorParameters

  function hernquistConstructorInternal(densityNormalization,mass,scaleLength,dimensionless) result(self)
    !% Internal constructor for ``hernquist'' mass distribution class.
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: Numerical_Comparison    , only : Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionHernquist)                          :: self
    double precision                           , intent(in   ), optional :: densityNormalization, mass, scaleLength
    logical                                    , intent(in   ), optional :: dimensionless

    ! Determine if profile is dimensionless.
    self%dimensionless=.false.
    if (present(dimensionless)) self%dimensionless=dimensionless
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (present(scaleLength         )) then
          if (Values_Differ(scaleLength         ,1.0d0   ,absTol=1.0d-6)) call Galacticus_Error_Report('scaleLength should be unity for a dimensionless profile (or simply do not specify a scale length)'                //{introspection:location})
       end if
       if (present(mass                )) then
          if (Values_Differ(mass                ,1.0d0   ,absTol=1.0d-6)) call Galacticus_Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'                               //{introspection:location})
       end if
       if (present(densityNormalization)) then
          if (Values_Differ(densityNormalization,0.5d0/Pi,absTol=1.0d-6)) call Galacticus_Error_Report('densityNormalization should be 1/2π for a dimensionless profile (or simply do not specify a densityNormalization)'//{introspection:location})
       end if
       self%scaleLength         =1.0d0
       self%mass                =1.0d0
       self%densityNormalization=0.5d0/Pi
    else
       if      (present(scaleLength         )) then
          self%scaleLength=scaleLength
       else
          call Galacticus_Error_Report('"scaleLength" must be specified'//{introspection:location})
       end if
       if      (present(densityNormalization)) then
          self%densityNormalization=densityNormalization
          self%mass                =densityNormalization*2.0d0*Pi*scaleLength**3
       else if (present(mass                )) then
          self%densityNormalization=mass                /2.0d0/Pi/scaleLength**3
          self%mass                =mass
       else
          call Galacticus_Error_Report('one of "densityNormalization" or "mass" must be specified'//{introspection:location})
       end if
    end if
    return
  end function hernquistConstructorInternal

  double precision function hernquistDensity(self,coordinates)
    !% Return the density at the specified {\normalfont \ttfamily coordinates} in a Hernquist mass distribution.
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class           (massDistributionHernquist), intent(inout) :: self
    class           (coordinate               ), intent(in   ) :: coordinates
    type            (coordinateSpherical      )                :: position
    double precision                                           :: r

    ! Get position in spherical coordinate system.
    position        = coordinates
    ! Compute the density at this position.
    r               =+position%r                   () &
         &           /self    %scaleLength
    hernquistDensity=+self    %densityNormalization   &
         &           /        r                       &
         &           /(+1.0d0+r)**3
    return
  end function hernquistDensity

  double precision function hernquistDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !% Returns a radial density moment for the Hernquist mass distribution.
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: Numerical_Comparison    , only : Values_Agree
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionHernquist), intent(inout)           :: self
    double precision                           , intent(in   )           :: moment
    double precision                           , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                    , intent(  out), optional :: isInfinite

    ! Abort on limited ranges.
    if (present(radiusMinimum).or.present(radiusMaximum)) call Galacticus_Error_Report('ranges are not supported'//{introspection:location})
    if (moment <= 0.0d0 .or. moment >= 3.0d0) then
       ! Handle cases where the moment is infinite.
       if (present(isInfinite)) then
          hernquistDensityRadialMoment=0.0d0
          isInfinite=.true.
          return
       else
          hernquistDensityRadialMoment=0.0d0
          call Galacticus_Error_Report('requested radial density moment is infinite'//{introspection:location})
       end if
    else if (                                              &
         &   Values_Agree(moment,1.0d0,absTol=1.0d-3) .or. &
         &   Values_Agree(moment,2.0d0,absTol=1.0d-3)      &
         &  ) then
       ! Handle special cases.
       if (present(isInfinite)) isInfinite=.false.
       hernquistDensityRadialMoment=0.5d0
    else
       ! Compute the general case.
       if (present(isInfinite)) isInfinite=.false.
       hernquistDensityRadialMoment=0.5d0*Pi*(moment-1.0d0)*(moment-2.0d0)/sin(Pi*moment)
    end if
    hernquistDensityRadialMoment=+hernquistDensityRadialMoment         &
         &                       *self%densityNormalization            &
         &                       *self%scaleLength            **moment
    return
  end function hernquistDensityRadialMoment

  double precision function hernquistMassEnclosedBySphere(self,radius)
    !% Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for Hernquist mass distributions.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionHernquist), intent(inout), target :: self
    double precision                           , intent(in   )         :: radius
    double precision                           , parameter             :: fractionalRadiusLarge=1.0d6
    double precision                                                   :: fractionalRadius

    fractionalRadius=radius/self%scaleLength
    if (fractionalRadius > fractionalRadiusLarge) then
       ! For very large radius approximate the mass enclosed as the total mass.
       hernquistMassEnclosedBySphere=2.0d0*Pi*self%densityNormalization
    else
       hernquistMassEnclosedBySphere=2.0d0*Pi*self%densityNormalization*fractionalRadius**2/(1.0d0+fractionalRadius)**2
    end if
    return
  end function hernquistMassEnclosedBySphere

  double precision function hernquistPotential(self,coordinates)
    !% Return the potential at the specified {\normalfont \ttfamily coordinates} in a Hernquist mass distribution.
    use :: Coordinates                 , only : assignment(=)                  , coordinateSpherical
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class(massDistributionHernquist), intent(inout) :: self
    class(coordinate               ), intent(in   ) :: coordinates
    type (coordinateSpherical      )                :: position

    ! Get position in spherical coordinate system.
    position=coordinates
    ! Compute the potential at this position.
    hernquistPotential=-self%mass/(self%scaleLength+position%r())
    if (.not.self%isDimensionless()) hernquistPotential=+gravitationalConstantGalacticus &
         &                                              *hernquistPotential
    return
  end function hernquistPotential

  double precision function hernquistRadiusHalfMass(self)
    !% Return the half-mass radius of a Hernquist mass distribution.
    implicit none
    class           (massDistributionHernquist), intent(inout) :: self
    double precision                           , parameter     :: radiusHalfMassToScaleRadius=1.0d0/(sqrt(2.0d0)-1.0d0)

    hernquistRadiusHalfMass=+radiusHalfMassToScaleRadius &
         &                  *self%scaleLength
    return
  end function hernquistRadiusHalfMass
