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
  Implementation of a mass distribution class that mimics the ``hydrostatic'' profile used by the Enzo simulation code.
  !!}

  !![
  <massDistribution name="massDistributionEnzoHydrostatic">
    <description>
      A hot halo mass distribution class which adopts a spherically symmetric density profile for the hot halo motivated by the
      ``hydrostatic'' profile available in the \gls{enzo} code. Specifically,
      \begin{equation}
      \rho_\mathrm{hot halo}(r) \propto \left\{ \begin{array}{ll} T^{-1} r^{-1} &amp; \hbox{ if } r &gt; r_\mathrm{core} \\ T^{-1}
      r_\mathrm{core}^{-1} &amp; \hbox{ if } r \le r_\mathrm{core}, \end{array} \right.
      \end{equation}
      where the core radius, $r_\mathrm{core}$, is set using the selected cored profile core radius method (see
      \refPhysics{hotHaloMassDistributionCoreRadius}). The profile is normalized such that the current mass in the
      hot gas profile is contained within the outer radius of the hot halo, $r_\mathrm{hot, outer}$. Note that the \gls{enzo}
      hydrostatic profile does not include this core, but without introducing this the profile mass can be divergent at small
      radii.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionEnzoHydrostatic
     !!{
     The Enzo hydrostatic mass distribution.
     !!}
     private
     double precision :: radiusScale          , radiusOuter                 , &
          &              mass                 , normalizationDensity_
     logical          :: truncateAtOuterRadius, normalizationDensityComputed
   contains
     !![
     <methods>
       <method description="Return the normalization of the density profile." method="normalizationDensity" />
     </methods>
     !!]
     procedure :: normalizationDensity  => enzoHydrostaticNormalizationDensity
     procedure :: density               => enzoHydrostaticDensity
     procedure :: densityGradientRadial => enzoHydrostaticDensityGradientRadial
  end type massDistributionEnzoHydrostatic

  interface massDistributionEnzoHydrostatic
     !!{
     Constructors for the \refClass{massDistributionEnzoHydrostatic} mass distribution class.
     !!}
     module procedure enzoHydrostaticConstructorParameters
     module procedure enzoHydrostaticConstructorInternal
  end interface massDistributionEnzoHydrostatic

contains

  function enzoHydrostaticConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionEnzoHydrostatic} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionEnzoHydrostatic)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    type            (varying_string                 )                :: componentType
    type            (varying_string                 )                :: massType
    double precision                                                 :: radiusScale          , radiusOuter, &
         &                                                              mass
    logical                                                          :: truncateAtOuterRadius

    !![
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
    <inputParameter>
      <name>mass</name>
      <description>The mass within the outer radius.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusOuter</name>
      <description>The outer radius of the mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusScale</name>
      <description>The core radius of the mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>truncateAtOuterRadius</name>
      <defaultValue>.false.</defaultValue>
      <description>If true then the mass distribution is truncated beyond the outer radius.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=massDistributionEnzoHydrostatic(mass,radiusOuter,radiusScale,truncateAtOuterRadius,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function enzoHydrostaticConstructorParameters

  function enzoHydrostaticConstructorInternal(mass,radiusOuter,radiusScale,truncateAtOuterRadius,componentType,massType) result(self)
    !!{
    Internal constructor for ``enzoHydrostatic'' mass distribution class.
    !!}
    implicit none
    type            (massDistributionEnzoHydrostatic)                          :: self
    double precision                                 , intent(in   )           :: radiusScale          , radiusOuter, &
         &                                                                        mass
    logical                                          , intent(in   ), optional :: truncateAtOuterRadius
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="radiusScale, radiusOuter, mass, truncateAtOuterRadius, componentType, massType"/>
    !!]

    ! This distribution profile is never dimensionless.
    self%dimensionless               =.false.
    self%normalizationDensityComputed=.false.
    return
  end function enzoHydrostaticConstructorInternal

  double precision function enzoHydrostaticNormalizationDensity(self) result(normalizationDensity)
    !!{
    Return the density normalization in a {\normalfont \ttfamily enzoHydrostatic} mass distribution.
    !!}
    implicit none
    class           (massDistributionEnzoHydrostatic), intent(inout) :: self
    double precision                                                 :: massEnclosed

    if (.not.self%normalizationDensityComputed) then
       self%normalizationDensityComputed=.true.
       if     (                           &
            &   self%mass        <= 0.0d0 &
            &  .or.                       &
            &   self%radiusOuter <= 0.0d0 &
            & ) then
          normalizationDensity=+0.0d0
       else
          ! Compute the normalization. Initially set this to 1, then compute the enclosed mass, then compute the multiplicative
          ! factor needed to get the required mass.
          self%normalizationDensity_=+1.0d0
          massEnclosed              =+self%massEnclosedBySphere(self%radiusOuter)
          self%normalizationDensity_=+self%mass                                   &
               &                     /     massEnclosed
       end if
    end if
    normalizationDensity=self%normalizationDensity_
    return
  end function enzoHydrostaticNormalizationDensity

  double precision function enzoHydrostaticDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an Enzo hydrostatic mass distribution.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class           (massDistributionEnzoHydrostatic), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    type            (coordinateSpherical            )                :: coordinatesEffective
    double precision                                                 :: radiusEffective

    radiusEffective     = max(coordinates%rSpherical(),self%radiusScale)
    coordinatesEffective= coordinates                   &
         &               *            radiusEffective   &
         &               /coordinates%rSpherical     ()
    density             =+self                        %normalizationDensity(                    )    &
         &               /self%kinematicsDistribution_%temperature         (coordinatesEffective)    &
         &               /radiusEffective                                                        **3

    return
  end function enzoHydrostaticDensity

  double precision function enzoHydrostaticDensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an EnzoHydrostatic \citep{navarro_structure_1996} mass distribution.
    !!}
    implicit none
    class  (massDistributionEnzoHydrostatic), intent(inout), target   :: self
    class  (coordinate                     ), intent(in   )           :: coordinates
    logical                                 , intent(in   ), optional :: logarithmic
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    if (coordinates%rSpherical() > self%radiusScale) then
       densityGradientRadial                       =-self       %kinematicsDistribution_%temperatureGradientLogarithmic(coordinates) &
            &                                       -3.0d0
       if (.not.logarithmic_) densityGradientRadial=+                                    densityGradientRadial                       &
            &                                       *self                               %density                       (coordinates) &
            &                                       /coordinates                        %rSpherical                    (           )
    else
       densityGradientRadial                       =+0.0d0
    end if    
    return
  end function enzoHydrostaticDensityGradientRadial
