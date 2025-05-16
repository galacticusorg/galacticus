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
  Implementation of a mass distribution for accretion flow using the fitting function of \cite{diemer_dependence_2014}.
  !!}

  !![
  <massDistribution name="massDistributionDiemerKravtsov2014">
    <description>
      A mass distribution class for accretion flows which models the accretion flow using the fitting function of
      \cite{diemer_dependence_2014}. Specifically, the density profile of the accretion flow is modeled using their equation~(4).
   </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionDiemerKravtsov2014
     !!{
     A mass distribution for accretion flow using the fitting function of \cite{diemer_dependence_2014}.
     !!}
     private
     double precision :: radius200Mean, densityMean, &
          &              b            , s
     logical          :: includeMean
   contains
     procedure :: density               => diemerKravtsov2014Density
     procedure :: densityGradientRadial => diemerKravtsov2014DensityGradientRadial
  end type massDistributionDiemerKravtsov2014

  interface massDistributionDiemerKravtsov2014
     !!{
     Constructors for the \refClass{massDistributionDiemerKravtsov2014} mass distribution class.
     !!}
     module procedure massDistributionDiemerKravtsov2014ConstructorParameters
     module procedure massDistributionDiemerKravtsov2014ConstructorInternal
  end interface massDistributionDiemerKravtsov2014

contains

  function massDistributionDiemerKravtsov2014ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionDiemerKravtsov2014} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionDiemerKravtsov2014)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    double precision                                                    :: radius200Mean, densityMean, &
         &                                                                 b            , s
    logical                                                             :: includeMean
    type            (varying_string                    )                :: componentType
    type            (varying_string                    )                :: massType

    !![
    <inputParameter>
      <name>densityMean</name>
      <description>The mean density of the universe in the \cite{diemer_dependence_2014} accretion flow mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radius200Mean</name>
      <description>The radius enclosing a density of 200 times the mean density of the universe in the \cite{diemer_dependence_2014} accretion flow mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeMean</name>
      <description>If true, include the mean density of the universe in the profile, otherwise, subtract off that mean density.</description>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>b</name>
      <description>The coefficient $b$ in the \cite{diemer_dependence_2014} accretion flow mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>s</name>
      <description>The exponent $s$ in the \cite{diemer_dependence_2014} accretion flow mass distribution.</description>
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
    self=massDistributionDiemerKravtsov2014(densityMean,radius200Mean,includeMean,b,s,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massDistributionDiemerKravtsov2014ConstructorParameters

  function massDistributionDiemerKravtsov2014ConstructorInternal(densityMean,radius200Mean,includeMean,b,s,componentType,massType) result(self)
    !!{
    Internal constructor for ``diemerKravtsov2014'' mass distribution class.
    !!}
    implicit none
    type            (massDistributionDiemerKravtsov2014)                          :: self
    double precision                                    , intent(in   )           :: densityMean  , radius200Mean, &
         &                                                                           b            , s
    logical                                             , intent(in   )           :: includeMean
    type            (enumerationComponentTypeType      ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType           ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="densityMean, radius200Mean, includeMean, b, s, componentType, massType"/>
    !!]

    self%dimensionless=.false.
    return
  end function massDistributionDiemerKravtsov2014ConstructorInternal

  double precision function diemerKravtsov2014Density(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a \cite{diemer_dependence_2014} mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class(massDistributionDiemerKravtsov2014), intent(inout) :: self
    class(coordinate                        ), intent(in   ) :: coordinates

    density=+self%densityMean              &
         &  *self%b                        &
         &  /(                             &
         &    +coordinates%rSpherical   () &
         &    /5.0d0                       &
         &    /self       %radius200Mean   &
         &   )**self%s
    if (self%includeMean)            &
         & density=+     density     &
         &         +self%densityMean
    return
  end function diemerKravtsov2014Density

  double precision function diemerKravtsov2014DensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density gradient at the specified {\normalfont \ttfamily coordinates} in a \cite{diemer_dependence_2014} mass distribution.
    !!}
    implicit none
    class  (massDistributionDiemerKravtsov2014), intent(inout), target   :: self
    class  (coordinate                        ), intent(in   )           :: coordinates
    logical                                    , intent(in   ), optional :: logarithmic
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    densityGradient=-self%s                          &
         &          /(                               &
         &            +1.0d0                         &
         &            +(                             &
         &              +coordinates%rSpherical   () &
         &              /5.0d0                       &
         &              /self       %radius200Mean   &
         &             )**self%s                     &
         &            /self%b                        &
         &          )
    if (.not.logarithmic_)                                            &
         & densityGradient=+             densityGradient              &
         &                 *self        %density        (coordinates) &
         &                 /coordinates%rSpherical      (           )
    return
  end function diemerKravtsov2014DensityGradientRadial
