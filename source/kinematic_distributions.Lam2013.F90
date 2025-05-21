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
  Implementation of a kinematic distribution class for the \cite{lam_modeling_2013} model of halo accretion flows.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Numerical_Interpolation, only : interpolator

  !![
  <kinematicsDistribution name="kinematicsDistributionLam2013">
   <description>A kinematic distribution class for the \cite{lam_modeling_2013} model of halo accretion flows.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionCollisionless) :: kinematicsDistributionLam2013
     !!{
     A kinematics distribution for the \cite{lam_modeling_2013} model of halo accretion flows.
     !!}
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_                => null()
     type            (interpolator           )                            :: correlationFunctionVolumeAveraged_
     double precision                         , allocatable, dimension(:) :: radius                                      , correlationFunctionVolumeAveraged
     double precision                                                     :: overdensityCritical                         , radiusVirial                     , &
          &                                                                  massVirial                                  , scaleFactorVelocity              , &
          &                                                                  redshift                                    , time                             , &
          &                                                                  rateLinearGrowth
   contains
     final     ::                   lam2013Destructor
     procedure :: velocityRadial => lam2013VelocityRadial
  end type kinematicsDistributionLam2013

  interface kinematicsDistributionLam2013
     !!{
     Constructors for the \refClass{kinematicsDistributionLam2013} kinematic distribution class.
     !!}
     module procedure lam2013ConstructorParameters
     module procedure lam2013ConstructorInternal
  end interface kinematicsDistributionLam2013

contains

  function lam2013ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionLam2013} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionLam2013)                              :: self
    type            (inputParameters              ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass      ), pointer                     :: cosmologyFunctions_
    double precision                               , allocatable  , dimension(:) :: radius             , correlationFunctionVolumeAveraged
    double precision                                                             :: radiusVirial       , massVirial                       , &
         &                                                                          overdensityCritical, redshift                         , &
         &                                                                          scaleFactorVelocity, rateLinearGrowth

    !![
    <inputParameter>
      <name>scaleFactorVelocity</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>A scale factor to be applied to inflow velocities.</description>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <description>The redshift of the halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massVirial</name>
      <description>The virial mass of the halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusVirial</name>
      <description>The virial radius of the halo.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>overdensityCritical</name>
      <description>The critical overdensity.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>rateLinearGrowth</name>
      <description>The logarithmic derivative of the linear growth factor with respect to expansion factor.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radius</name>
      <description>The radius in the tabulated volume-averaged correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>correlationFunctionVolumeAveraged</name>
      <description>The correlation in the tabulated volume-averaged correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=kinematicsDistributionLam2013(massVirial,radiusVirial,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),overdensityCritical,rateLinearGrowth,scaleFactorVelocity,radius,correlationFunctionVolumeAveraged,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function lam2013ConstructorParameters
  
  function lam2013ConstructorInternal(massVirial,radiusVirial,time,overdensityCritical,rateLinearGrowth,scaleFactorVelocity,radius,correlationFunctionVolumeAveraged,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionLam2013} kinematic distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionLam2013)                              :: self
    class           (cosmologyFunctionsClass      ), intent(in   ), target       :: cosmologyFunctions_
    double precision                               , intent(in   )               :: radiusVirial       , massVirial                       , &
         &                                                                          overdensityCritical, time                             , &
         &                                                                          scaleFactorVelocity, rateLinearGrowth
    double precision                               , intent(in   ), dimension(:) :: radius             , correlationFunctionVolumeAveraged
    !![
    <constructorAssign variables="massVirial, radiusVirial, overdensityCritical, rateLinearGrowth, time, scaleFactorVelocity, radius, correlationFunctionVolumeAveraged, *cosmologyFunctions_"/>
    !!]
    
    self%redshift                          =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    self%correlationFunctionVolumeAveraged_=interpolator(radius,correlationFunctionVolumeAveraged)
    return
  end function lam2013ConstructorInternal
  
  subroutine lam2013Destructor(self)
    !!{
    Destructor for the \refClass{kinematicsDistributionLam2013} accretion flow mass distribution class.
    !!}
    implicit none
    type(kinematicsDistributionLam2013), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine lam2013Destructor
  
  double precision function lam2013VelocityRadial(self,coordinates,massDistributionEmbedding) result(velocityRadial)
    !!{
     Return the radial velocity at the specified {\normalfont \ttfamily coordinates} in the \cite{lam_modeling_2013} model for the accretion flow around a halo.
     !!}
     implicit none
    class           (kinematicsDistributionLam2013), intent(inout) :: self
    class           (coordinate                   ), intent(in   ) :: coordinates
    class           (massDistributionClass        ), intent(inout) :: massDistributionEmbedding
    double precision                                               :: massShell                , radius, &
         &                                                            densityContrastNonLinear

    radius=coordinates%rSpherical()
    ! Evaluate the mass in the shell outside the halo virial radius using equation (B4) of Lam et al. (2013).
    if (radius > self%radiusVirial) then
       massShell            =+self %cosmologyFunctions_ %matterDensityEpochal              (self%time)                              &
            &                *4.0d0                                                                                                 &
            &                *Pi                                                                                                    &
            &                /3.0d0                                                                                                 &
            &                *(                                                                                                     &
            &                  +     radius      **3*(1.0d0+self%correlationFunctionVolumeAveraged_%interpolate(     radius      )) &
            &                  -self%radiusVirial**3*(1.0d0+self%correlationFunctionVolumeAveraged_%interpolate(self%radiusVirial)) &
            &                 )
    else
       massShell            =+0.0d0
    end if
    ! Compute the nonlinear density contrast using equation (B1) of Lam et al. (2013).
    densityContrastNonlinear=+(                                                                                                     &
         &                     +self%massVirial                                                                                     &
         &                     +      massShell                                                                                     &
         &                    )                                                                                                     &
         &                   /self%cosmologyFunctions_%matterDensityEpochal                (self%time)                              &
         &                   /(                                                                                                     &
         &                     +4.0d0                                                                                               &
         &                     *Pi                                                                                                  &
         &                     /3.0d0                                                                                               &
         &                     *radius**3                                                                                           &
         &                    )
    ! Evaluate the inflow velocity in the spherical collapse model using equation (B2) of Lam et al. (2013).
    velocityRadial          =-self%scaleFactorVelocity                                                                              &
         &                   *self%cosmologyFunctions_%hubbleParameterEpochal              (self%time)                              &
         &                   *radius                                                                                                &
         &                   *self%cosmologyFunctions_%expansionFactor                     (self%time)                              &
         &                   *self%rateLinearGrowth                                                                                 &
         &                   /3.0d0                                                                                                 &
         &                   *                                   self%overdensityCritical                                           &
         &                   *(                                                                                                     &
         &                     +densityContrastNonLinear**(1.0d0/self%overdensityCritical)                                          &
         &                     -1.0d0                                                                                               &
         &                   )
    return
  end function lam2013VelocityRadial
