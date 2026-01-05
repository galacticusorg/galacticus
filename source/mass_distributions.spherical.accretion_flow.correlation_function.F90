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
  Implementation of a mass distribution for accretion flow using the 2-halo correlation function.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Numerical_Interpolation, only : interpolator
  
  !![
  <massDistribution name="massDistributionCorrelationFunction">
    <description>
      An accretion flow class which models the accretion flow using the 2-halo correlation function.
   </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionCorrelationFunction
     !!{
     A mass distribution for accretion flow using the 2-halo correlation function.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_  => null()
     type            (interpolator           )                            :: correlationFunction_
     double precision                                                     :: mass                          , time               , &
          &                                                                  redshift
     double precision                         , allocatable, dimension(:) :: radius                        , correlationFunction
   contains
     procedure :: density               => correlationFunctionDensity
     procedure :: densityGradientRadial => correlationFunctionDensityGradientRadial 
 end type massDistributionCorrelationFunction

  interface massDistributionCorrelationFunction
     !!{
     Constructors for the \refClass{massDistributionCorrelationFunction} mass distribution class.
     !!}
     module procedure massDistributionCorrelationFunctionConstructorParameters
     module procedure massDistributionCorrelationFunctionConstructorInternal
  end interface massDistributionCorrelationFunction

contains

  function massDistributionCorrelationFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionCorrelationFunction} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    use :: Numerical_Interpolation   , only : interpolator
    implicit none
    type            (massDistributionCorrelationFunction)                            :: self
    type            (inputParameters                    ), intent(inout)             :: parameters
    class           (cosmologyFunctionsClass            ), pointer                   :: cosmologyFunctions_
    double precision                                     , allocatable, dimension(:) :: radius             , correlationFunction
    double precision                                                                 :: mass               , redshift
    type            (varying_string                     )                            :: componentType
    type            (varying_string                     )                            :: massType

    !![
    <inputParameter>
      <name>mass</name>
      <description>The mass of the halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <description>The redshift of the halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radius</name>
      <description>The radius in the tabulated correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>correlationFunction</name>
      <description>The correlation in the tabulated correlation function.</description>
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
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=massDistributionCorrelationFunction(mass,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),radius,correlationFunction,cosmologyFunctions_,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function massDistributionCorrelationFunctionConstructorParameters

  function massDistributionCorrelationFunctionConstructorInternal(mass,time,radius,correlationFunction,cosmologyFunctions_,componentType,massType) result(self)
    !!{
    Internal constructor for ``correlationFunction'' mass distribution class.
    !!}
    implicit none
    type            (massDistributionCorrelationFunction)                              :: self
    class           (cosmologyFunctionsClass            ), intent(in   ), target       :: cosmologyFunctions_
    double precision                                     , intent(in   )               :: mass               , time
    double precision                                     , intent(in   ), dimension(:) :: radius             , correlationFunction
    type            (enumerationComponentTypeType       ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType            ), intent(in   ), optional     :: massType
    !![
    <constructorAssign variables="mass, time, radius, correlationFunction, *cosmologyFunctions_, componentType, massType"/>
    !!]

    self%dimensionless       =.false.
    self%redshift            =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    self%correlationFunction_=interpolator(radius,correlationFunction)
    return
  end function massDistributionCorrelationFunctionConstructorInternal

  subroutine correlationFunctionDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionCorrelationFunction} accretion flow mass distribution class.
    !!}
    implicit none
    type(massDistributionCorrelationFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine correlationFunctionDestructor

  double precision function correlationFunctionDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a accretion flow modeled on the 2-halo correlation function.
    !!}
    implicit none
    class(massDistributionCorrelationFunction), intent(inout) :: self
    class(coordinate                         ), intent(in   ) :: coordinates

    density=+(                                                                          &
         &    +1.0d0                                                                    &
         &    +self%correlationFunction_%interpolate         (coordinates%rSpherical()) &
         &   )                                                                          &
         &  *  self%cosmologyFunctions_ %matterDensityEpochal(self       %time        )
    return
  end function correlationFunctionDensity

  double precision function correlationFunctionDensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a accretion flow modeled on the 2-halo correlation function.
    !!}
    implicit none
    class  (massDistributionCorrelationFunction), intent(inout), target   :: self
    class  (coordinate                         ), intent(in   )           :: coordinates
    logical                                     , intent(in   ), optional :: logarithmic
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    densityGradientRadial=+self%correlationFunction_%derivative          (coordinates%rSpherical()) &
         &                *self%cosmologyFunctions_ %matterDensityEpochal(self       %time        )
    return
  end function correlationFunctionDensityGradientRadial
