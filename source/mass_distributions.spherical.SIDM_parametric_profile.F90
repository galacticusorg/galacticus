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

  !+    Contributions to this file made by: Niusha Ahvazi
  
  !!{RST
  Implementation of a mass distribution class for the SIDM parametric profile of :cite:t:`yang_parametric_2024`.
  !!}

  !![
  <massDistribution name="massDistributionSIDMParametricProfile" docformat="rst">
    <description>
    A mass distribution class for the SIDM parametric profile of :cite:t:`yang_parametric_2024`. The density profile is given by:

    .. math::

        \rho(r) = \rho_\mathrm{s} \left[ \left( \left[\frac{r}{r_\mathrm{s}}\right]^\beta +
       \left[\frac{r_\mathrm{c}}{r_\mathrm{s}}\right]^\beta \right)^{1/\beta} \left( 1 + \frac{r}{r_\mathrm{s}} \right)^2
       \right]^{-1}.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSphericalTabulated) :: massDistributionSIDMParametricProfile
     !!{RST
     A mass distribution class for the SIDM parametric profile of :cite:t:`yang_parametric_2024`.
     !!}
     double precision :: beta          , radiusScale         , &
          &              radiusCore    , densityNormalization, &
          &              densityCentral
   contains
     procedure :: density                         => SIDMParametricProfileDensity
     procedure :: densityGradientRadial           => SIDMParametricProfileDensityGradientRadial
     procedure :: radiusEnclosingDensity          => SIDMParametricProfileRadiusEnclosingDensity
     procedure :: radiusEnclosingDensityNumerical => SIDMParametricProfileRadiusEnclosingDensityNumerical
     procedure :: parameters                      => SIDMParametricProfileParameters
     procedure :: factoryTabulation               => SIDMParametricProfileFactoryTabulation
     procedure :: descriptor                      => SIDMParametricProfileDescriptor
     procedure :: suffix                          => SIDMParametricProfileSuffix
  end type massDistributionSIDMParametricProfile
  
  interface massDistributionSIDMParametricProfile
     !!{RST
     Constructors for the :galacticus-class:`massDistributionSIDMParametricProfile` mass distribution class.
     !!}
     module procedure SIDMParametricProfileConstructorParameters
     module procedure SIDMParametricProfileConstructorInternal
  end interface massDistributionSIDMParametricProfile

  ! Tabulated solutions.
  logical                                     :: containerSIDMParametricProfileInitialized=.false.
  type   (massDistributionContainer), pointer :: containerSIDMParametricProfile
  !$omp threadprivate(containerSIDMParametricProfile,containerSIDMParametricProfileInitialized)
  
  ! Generate a source digest.
  !![
  <sourceDigest name="massDistributionSIDMParametricProfileSourceDigest"/>
  !!]

contains

  function SIDMParametricProfileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionSIDMParametricProfile` mass distribution class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSIDMParametricProfile)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: beta         , densityNormalization, &
         &                                                                    radiusScale  , radiusCore
    type            (varying_string                       )                :: componentType
    type            (varying_string                       )                :: massType

    !![
    <inputParameter docformat="rst">
      <name>beta</name>
      <defaultValue>4.0d0</defaultValue>
      <description>
      The value :math:`\beta` in a SIDM parametric mass distribution.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>densityNormalization</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The density normalization of a SIDM parametric mass distribution.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusScale</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The scale of a SIDM parametric mass distribution.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusCore</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The core radius of a SIDM parametric mass distribution.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>
      The component type that this mass distribution represents.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>
      The mass type that this mass distribution represents.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=massDistributionSIDMParametricProfile(beta,densityNormalization,radiusScale,radiusCore,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.))

    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function SIDMParametricProfileConstructorParameters

  function SIDMParametricProfileConstructorInternal(beta,densityNormalization,radiusScale,radiusCore,componentType,massType) result(self)
    !!{RST
    Internal constructor for :galacticus-class:`massDistributionSIDMParametricProfile` mass distribution class.
    !!}
    implicit none
    type            (massDistributionSIDMParametricProfile)                          :: self
    double precision                                       , intent(in   )           :: beta
    double precision                                       , intent(in   ), optional :: densityNormalization, radiusScale, &
         &                                                                              radiusCore
    type            (enumerationComponentTypeType         ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType              ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="beta, densityNormalization, radiusScale, radiusCore, componentType, massType"/>
    !!]

    self%dimensionless=.false.
    self%densityCentral=+densityNormalization &
         &              *radiusScale          &
         &              /radiusCore
    return
  end function SIDMParametricProfileConstructorInternal

  function SIDMParametricProfileFactoryTabulation(self,parameters) result(instance)
    !!{RST
    Construct an instance of this class using tabulation parameters.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated   )               , pointer      :: instance
    class           (massDistributionSIDMParametricProfile), intent(inout)               :: self
    double precision                                       , intent(in   ), dimension(:) :: parameters
    !$GLC attributes unused :: self
    
    allocate(massDistributionSIDMParametricProfile :: instance)
    select type(instance)
    type is (massDistributionSIDMParametricProfile)
       instance=massDistributionSIDMParametricProfile(radiusScale=1.0d0,densityNormalization=1.0d0,radiusCore=parameters(1),beta=parameters(2))
    end select
    return
  end function SIDMParametricProfileFactoryTabulation
  
  double precision function SIDMParametricProfileDensity(self,coordinates) result(density)
    !!{RST
    Return the density at the specified ``coordinates`` in a SIDM parametric mass distribution.
    !!}
    implicit none
    class           (massDistributionSIDMParametricProfile), intent(inout) :: self
    class           (coordinate                           ), intent(in   ) :: coordinates
    double precision                                                       :: radius
    
    radius =+coordinates%rSpherical()
    density=+self%densityNormalization                       &
         &  /(                                               &
         &    +1.0d0                                         &
         &    +      radius    /self%radiusScale             &
         &   )** 2                                           &
         &  /(                                               &
         &    +(     radius    /self%radiusScale)**self%beta &
         &    +(self%radiusCore/self%radiusScale)**self%beta &
         &   )**(1.0d0/self%beta)
    return
  end function SIDMParametricProfileDensity

  double precision function SIDMParametricProfileDensityGradientRadial(self,coordinates,logarithmic) result(densityGradientRadial)
    !!{RST
    Return the density at the specified ``coordinates`` in a SIDM parametric mass distribution.
    !!}
    implicit none
    class           (massDistributionSIDMParametricProfile), intent(inout), target   :: self
    class           (coordinate                           ), intent(in   )           :: coordinates
    logical                                                , intent(in   ), optional :: logarithmic
    double precision                                                                 :: radius
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    ! Get position in spherical coordinate system.
    radius=coordinates%rSpherical()
    ! Compute the density gradient.
    if (radius > 0.0d0) then
       ! Evaluate the logarithmic density profile gradient.
       densityGradientRadial=-1.0d0                                              &
            &                +1.0d0/(1.0d0+(radius/self%radiusCore )**self%beta) &
            &                -2.0d0*radius/(radius+self%radiusScale)
       ! Convert to non-logarithmic form if necessary.
       if (.not.logarithmic_)                                                &
            & densityGradientRadial=+     densityGradientRadial              &
            &                       *self%density              (coordinates) &
            &                       /     radius
    else
       densityGradientRadial=+0.0d0
    end if
    return
  end function SIDMParametricProfileDensityGradientRadial
  
  double precision function SIDMParametricProfileRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{RST
    Computes the radius enclosing a given mean density for SIDM parametric mass distributions.
    !!}
    implicit none
    class           (massDistributionSIDMParametricProfile), intent(inout), target   :: self
    double precision                                       , intent(in   )           :: density
    double precision                                       , intent(in   ), optional :: radiusGuess
    double precision                                       , parameter               :: epsilonDeltaDensityFractional=1.0d-3

    if (density >= self%densityCentral) then
       ! Above the central density, return a radius of zero.
       radius=0.0d0
    else if (density > self%densityCentral*(1.0d0-epsilonDeltaDensityFractional)) then
       ! For densities close to, but below the central density, use a series solution.
       radius=+0.5d0                       &
            & *(                           &
            &   +self%radiusScale          &
            &   -self%radiusCore           &
            &   *     density              &
            &   /self%densityNormalization &
            &  )
    else
       radius=sphericalTabulatedRadiusEnclosingDensity(self,density,radiusGuess)
    end if
    return
  end function SIDMParametricProfileRadiusEnclosingDensity

  double precision function SIDMParametricProfileRadiusEnclosingDensityNumerical(self,density,radiusGuess) result(radius)
    !!{RST
    Computes the radius enclosing a given mean density for SIDM parametric mass distributions.
    !!}
    implicit none
    class           (massDistributionSIDMParametricProfile), intent(inout), target   :: self
    double precision                                       , intent(in   )           :: density
    double precision                                       , intent(in   ), optional :: radiusGuess
    double precision                                       , parameter               :: epsilonDeltaDensityFractional=1.0d-3

    if (density >= self%densityCentral) then
       ! Above the central density, return a radius of zero.
       radius=+0.0d0
    else if (density >  self%densityCentral*(1.0d0-epsilonDeltaDensityFractional)) then
       ! For densities close to, but below the central density, use a series solution.
       radius=+0.5d0                       &
            & *(                           &
            &   +self%radiusScale          &
            &   -self%radiusCore           &
            &   *     density              &
            &   /self%densityNormalization &
            &  )
    else
       radius=massDistributionRadiusEnclosingDensityNumerical(self,density,radiusGuess)
    end if
    return
  end function SIDMParametricProfileRadiusEnclosingDensityNumerical

  subroutine SIDMParametricProfileParameters(self,densityNormalization,radiusNormalization,parameters,container)
    !!{RST
    Establish parameters for tabulation.
    !!}
    implicit none
    class           (massDistributionSIDMParametricProfile ), intent(inout)                              :: self
    double precision                                        , intent(  out)                              :: densityNormalization, radiusNormalization
    double precision                                        , intent(inout), allocatable, dimension(:  ) :: parameters
    type            (massDistributionContainer             ), intent(  out), pointer                     :: container

    if (.not.containerSIDMParametricProfileInitialized) then
       ! Allocate the table and initialize.
       allocate(containerSIDMParametricProfile)
       call containerSIDMParametricProfile%initialize(2)
       ! Specify the number of tabulation points per interval in radius and each parameter.
       containerSIDMParametricProfile%mass                      %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%mass                      %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%radiusEnclosingDensity    %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%radiusEnclosingDensity    %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%potential                 %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%potential                 %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%velocityDispersion1D      %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%velocityDispersion1D      %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%energy                    %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%energy                    %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%radiusFreefall            %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%radiusFreefall            %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%radiusFreefallIncreaseRate%radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%radiusFreefallIncreaseRate%parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%densityRadialMoment0      %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%densityRadialMoment0      %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%densityRadialMoment1      %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%densityRadialMoment1      %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%densityRadialMoment2      %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%densityRadialMoment2      %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%densityRadialMoment3      %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%densityRadialMoment3      %parametersCountPer   =+20_c_size_t
       containerSIDMParametricProfile%fourierTransform          %radiusCountPer       =+20_c_size_t
       containerSIDMParametricProfile%fourierTransform          %parametersCountPer   =+20_c_size_t
       ! Specify names and descriptions of the parameters.
       containerSIDMParametricProfile%nameParameters                               (1)='radiusCoreOverRadiusScale'
       containerSIDMParametricProfile%descriptionParameters                        (1)='The ratio of core to scale radii.'
       containerSIDMParametricProfile%nameParameters                               (2)='beta'
       containerSIDMParametricProfile%descriptionParameters                        (2)='The β parameter.'
       ! Record that the table is now initialized.
       containerSIDMParametricProfileInitialized                                      =.true.
    end if
    allocate(parameters(2))
    densityNormalization    =  +self%densityNormalization
    radiusNormalization     =  +self%radiusScale
    parameters          (1) =  +self%radiusCore                     &
         &                     /self%radiusScale
    parameters          (2) =  +self%beta
    container               =>       containerSIDMParametricProfile
    return
  end subroutine SIDMParametricProfileParameters

  subroutine SIDMParametricProfileDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionSIDMParametricProfile), intent(inout)           :: self
    type     (inputParameters                      ), intent(inout)           :: descriptor
    logical                                         , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                               )                          :: parameterLabel
    type     (inputParameters                      )                          :: parameters
    !$GLC attributes unused :: includeFileModificationTimes
    
    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','SIDMParametricProfile')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%radiusScale
    call parameters%addParameter('radiusScale'         ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%radiusCore
    call parameters%addParameter('radiusCore'          ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%beta
    call parameters%addParameter('beta'                ,trim(adjustl(parameterLabel)))
    return
  end subroutine SIDMParametricProfileDescriptor

  function SIDMParametricProfileSuffix(self) result(suffix)
    !!{RST
    Return a suffix for tabulated file names.
    !!}
    use :: String_Handling, only : String_C_To_Fortran
    implicit none
    type (varying_string                       )                :: suffix
    class(massDistributionSIDMParametricProfile), intent(inout) :: self
    !$GLC attributes unused :: self

    suffix=String_C_To_Fortran(massDistributionSIDMParametricProfileSourceDigest)
    return
  end function SIDMParametricProfileSuffix
