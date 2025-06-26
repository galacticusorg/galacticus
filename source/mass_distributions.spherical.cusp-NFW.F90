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
  Implementation of the cusp-NFW \citep{delos_cusp-halo_2025} mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionCuspNFW">
    <description>
      The cusp-NFW mass distribution \citep{delos_cusp-halo_2025}. The density profile is given by:
      \begin{equation}
       \rho_\mathrm{dark matter}(r) = \rho_\mathrm{s} \left(y^2+{r\over r_\mathrm{s}}\right)^{1/2} \left({r\over r_\mathrm{s}}\right)^{-3/2} \left[1 + \left({r\over r_\mathrm{s}}\right) \right]^{-2},
      \end{equation}      
      where $\rho_\mathrm{s}$ and $r_\mathrm{s}$ are the usual NFW density normalization and scale length, and $y = A/\rho_s
      r_\mathrm{s}^{3/2}$ characterizes the amplitude of the cusp, with $A$ being the ``cusp coefficient'' as defined by
      \cite{delos_cusp-halo_2025}.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSphericalTabulated) :: massDistributionCuspNFW
     !!{
     The cusp-NFW \citep{delos_cusp-halo_2025} mass distribution.
     !!}
     private
     double precision :: densityNormalization, radiusScale, &
          &              y
   contains
     procedure :: density               => cuspNFWDensity
     procedure :: densityGradientRadial => cuspNFWDensityGradientRadial
     procedure :: massEnclosedBySphere  => cuspNFWMassEnclosedBySphere
     procedure :: parameters            => cuspNFWParameters
     procedure :: factoryTabulation     => cuspNFWFactoryTabulation
     procedure :: descriptor            => cuspNFWDescriptor
     procedure :: suffix                => cuspNFWSuffix
  end type massDistributionCuspNFW
  
  interface massDistributionCuspNFW
     !!{
     Constructors for the \refClass{massDistributionCuspNFW} mass distribution class.
     !!}
     module procedure cuspNFWConstructorParameters
     module procedure cuspNFWConstructorInternal
  end interface massDistributionCuspNFW

  ! Tabulated solutions.
  logical                                     :: containerCuspNFWInitialized=.false.
  type   (massDistributionContainer), pointer :: containerCuspNFW
  !$omp threadprivate(containerCuspNFW,containerCuspNFWInitialized)
  
  ! Generate a source digest.
  !![
  <sourceDigest name="massDistributionCuspNFWSourceDigest"/>
  !!]

contains

  function cuspNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionCuspNFW} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    type            (massDistributionCuspNFW)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    double precision                                          :: mass                      , radiusScale  , &
         &                                                       densityNormalization      , concentration, &
         &                                                       radiusVirial              , y            , &
         &                                                       toleranceRelativePotential
    logical                                                   :: dimensionless
    type            (varying_string          )                :: componentType
    type            (varying_string          )                :: massType

    !![
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>1.0d0/2.0d0/Pi/(log(4.0d0)-1.0d0)</defaultValue>
      <description>The density normalization of the cusp-NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusScale</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the cusp-NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>y</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The cusp amplitude parameter the cusp-NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the cored profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>concentration</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The concentration of the cusp-NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusVirial</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The virial radius of the cusp-NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the cusp-NFW profile is considered to be dimensionless.</description>
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
    <inputParameter>
      <name>toleranceRelativePotential</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the gravitational potential.</description>
    </inputParameter>
    <conditionalCall>
     <call>self=massDistributionCuspNFW(y=y,toleranceRelativePotential=toleranceRelativePotential,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
     <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
     <argument name="radiusScale"          value="radiusScale"          parameterPresent="parameters"/>
     <argument name="radiusVirial"         value="radiusVirial"         parameterPresent="parameters"/>
     <argument name="concentration"        value="concentration"        parameterPresent="parameters"/>
     <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cuspNFWConstructorParameters

  function cuspNFWConstructorInternal(radiusScale,y,concentration,densityNormalization,mass,radiusVirial,dimensionless,componentType,massType,toleranceRelativePotential) result(self)
    !!{
    Internal constructor for \refClass{massDistributionCuspNFW} mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionCuspNFW     )                          :: self
    double precision                              , intent(in   ), optional :: radiusScale               , concentration      , &
         &                                                                     densityNormalization      , mass               , &
         &                                                                     radiusVirial              , y                  , &
         &                                                                     toleranceRelativePotential
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    !![
    <constructorAssign variables="componentType, massType, toleranceRelativePotential"/>
    !!]

    ! Determine y parameter.
    if      (                            &
         &   present(y            )      &
         &  ) then
       self%y                   =y
    else
       call Error_Report('no means to determinecusp amplitude' //{introspection:location})
    end if
    ! Determine scale radius.
    if      (                            &
         &   present(radiusScale  )      &
         &  ) then
       self%radiusScale         =radiusScale
    else if (                            &
         &   present(concentration).and. &
         &   present(radiusVirial )      &
         &  ) then
       self%radiusScale         =radiusVirial/concentration
    else
       call Error_Report('no means to determine scale radius'//{introspection:location})
    end if
    ! Determine density normalization.
    if      (                                   &
         &   present(densityNormalization)      &
         &  ) then
       self%densityNormalization=densityNormalization
    else if (                                   &
         &   present(mass                ).and. &
         &   present(radiusVirial        )      &
         &  ) then
       radiusScaleFree          =+     radiusVirial/self%radiusScale
       self%densityNormalization=+mass                                                                                                                 &
            &                    /self%radiusScale**3                                                                                                  &
            &                    /4.0d0                                                                                                                &
            &                    /Pi                                                                                                                   &
            &                    /(                                                                                                                    &
            &                      + 2.0d0                                 *asinh(sqrt(radiusScaleFree                 )/                 self%y     ) &
            &                      -(2.0d0-self%y**2)/sqrt(1.0d0-self%y**2)*atanh(sqrt(radiusScaleFree*(1.0d0-self%y**2)/(radiusScaleFree+self%y**2))) &
            &                      -sqrt(radiusScaleFree*(radiusScaleFree+self%y**2))/(1.0d0+radiusScaleFree)                                          &
            &                     )
    else
       call Error_Report('either "densityNormalization", or "mass" and "radiusVirial" must be specified'//{introspection:location})
    end if
    ! Determine if profile is dimensionless.
    if      (present(dimensionless     )) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    return
  end function cuspNFWConstructorInternal

  function cuspNFWFactoryTabulation(self,parameters) result(instance)
    !!{
    Construct an instance of this class using tabulation parameters.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated)               , pointer      :: instance
    class           (massDistributionCuspNFW           ), intent(inout)               :: self
    double precision                                    , intent(in   ), dimension(:) :: parameters
    !$GLC attributes unused :: self
    
    allocate(massDistributionCuspNFW :: instance)
    select type(instance)
    type is (massDistributionCuspNFW)
       instance=massDistributionCuspNFW(radiusScale=1.0d0,densityNormalization=1.0d0,y=parameters(1),toleranceRelativePotential=self%toleranceRelativePotential)
    end select
    return
  end function cuspNFWFactoryTabulation
  
  double precision function cuspNFWDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a cusp-NFW mass distribution.
    !!}
    implicit none
    class           (massDistributionCuspNFW), intent(inout) :: self
    class           (coordinate             ), intent(in   ) :: coordinates
    double precision                                         :: radiusScaleFree

    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %radiusScale
    density        =+self       %densityNormalization       &
         &          *sqrt(self%y**2+radiusScaleFree)        &
         &          /               radiusScaleFree **1.5d0 &
         &          /    (1.0d0    +radiusScaleFree)**2
    return
  end function cuspNFWDensity
  
  double precision function cuspNFWDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a cusp-NFW mass distribution.
    !!}
    implicit none
    class           (massDistributionCuspNFW), intent(inout), target   :: self
    class           (coordinate             ), intent(in   )           :: coordinates
    logical                                  , intent(in   ), optional :: logarithmic
    double precision                                                   :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %radiusScale
    densityGradient=-3.5d0                                                 &
         &          +2.0d0                /(     1.0d0   +radiusScaleFree) &
         &          +0.5d0*radiusScaleFree/(self%y    **2+radiusScaleFree)
    if (.not.logarithmic_)                                    &
         & densityGradient=+     densityGradient              &
         &                 *self%density        (coordinates) &
         &                 /self%radiusScale                  &
         &                 /     radiusScaleFree
    return
  end function cuspNFWDensityGradientRadial

  double precision function cuspNFWMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Return the mass enclosed by a sphere of the specified {\normalfont \ttfamily radius} in a cusp-NFW mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionCuspNFW), intent(inout), target :: self
    double precision                         , intent(in   )         :: radius
    double precision                         , parameter             :: fractionSmall  =1.0d-3
    double precision                                                 :: radiusScaleFree

    radiusScaleFree=+     radius      &
         &          /self%radiusScale
    if (radiusScaleFree < fractionSmall*self%y**2) then
       ! Use a series expansion for small radii.
       mass  =+8.0d0                            &
            & *Pi                               &
            & /3.0d0                            &
            & *self%densityNormalization        &
            & *self%radiusScale         **3     &
            & *     radiusScaleFree     **1.5d0 &
            & *self%y                           &
            & *(                                &
            &   + 1.0d0                         &
            &   + 3.0d0                         &
            &   /10.0d0                         &
            &   *radiusScaleFree                &
            &   *(+1.0d0-4.0d0*self%y**2)       &
            &   /              self%y**2        &
            & )
    else
       ! Use the full solution for larger radii.
       mass  =+4.0d0                                                                                                                &
            & *Pi                                                                                                                   &
            & *self%radiusScale**3                                                                                                  &
            & *self%densityNormalization                                                                                            &
            & *(                                                                                                                    &
            &   + 2.0d0                                 *asinh(sqrt(radiusScaleFree                 )/                 self%y     ) &
            &   -(2.0d0-self%y**2)/sqrt(1.0d0-self%y**2)*atanh(sqrt(radiusScaleFree*(1.0d0-self%y**2)/(radiusScaleFree+self%y**2))) &
            &   -sqrt(radiusScaleFree*(radiusScaleFree+self%y**2))/(1.0d0+radiusScaleFree)                                          &
            &  )
    end if
    return
  end function cuspNFWMassEnclosedBySphere

  subroutine cuspNFWParameters(self,densityNormalization,radiusNormalization,parameters,container)
    !!{
    Establish parameters for tabulation.
    !!}
    implicit none
    class           (massDistributionCuspNFW  ), intent(inout)                              :: self
    double precision                           , intent(  out)                              :: densityNormalization, radiusNormalization
    double precision                           , intent(inout), allocatable, dimension(:  ) :: parameters
    type            (massDistributionContainer), intent(  out), pointer                     :: container

    if (.not.containerCuspNFWInitialized) then
       ! Allocate the table and initialize.
       allocate(containerCuspNFW)
       call containerCuspNFW%initialize(1)
       ! Specify the number of tabulation points per interval in radius and each parameter.
       containerCuspNFW%mass                      %radiusCountPer       =+20_c_size_t
       containerCuspNFW%mass                      %parametersCountPer   =+20_c_size_t
       containerCuspNFW%potential                 %radiusCountPer       =+20_c_size_t
       containerCuspNFW%potential                 %parametersCountPer   =+20_c_size_t
       containerCuspNFW%velocityDispersion1D      %radiusCountPer       =+20_c_size_t
       containerCuspNFW%velocityDispersion1D      %parametersCountPer   =+20_c_size_t
       containerCuspNFW%energy                    %radiusCountPer       =+20_c_size_t
       containerCuspNFW%energy                    %parametersCountPer   =+20_c_size_t
       containerCuspNFW%radiusFreefall            %radiusCountPer       =+20_c_size_t
       containerCuspNFW%radiusFreefall            %parametersCountPer   =+20_c_size_t
       containerCuspNFW%radiusFreefallIncreaseRate%radiusCountPer       =+20_c_size_t
       containerCuspNFW%radiusFreefallIncreaseRate%parametersCountPer   =+20_c_size_t
       containerCuspNFW%densityRadialMoment0      %radiusCountPer       =+20_c_size_t
       containerCuspNFW%densityRadialMoment0      %parametersCountPer   =+20_c_size_t
       containerCuspNFW%densityRadialMoment1      %radiusCountPer       =+20_c_size_t
       containerCuspNFW%densityRadialMoment1      %parametersCountPer   =+20_c_size_t
       containerCuspNFW%densityRadialMoment2      %radiusCountPer       =+20_c_size_t
       containerCuspNFW%densityRadialMoment2      %parametersCountPer   =+20_c_size_t
       containerCuspNFW%densityRadialMoment3      %radiusCountPer       =+20_c_size_t
       containerCuspNFW%densityRadialMoment3      %parametersCountPer   =+20_c_size_t
       containerCuspNFW%fourierTransform          %radiusCountPer       =+20_c_size_t
       containerCuspNFW%fourierTransform          %parametersCountPer   =+20_c_size_t
       ! Specify names and descriptions of the parameters.
       containerCuspNFW%nameParameters                               (1)='y'
       containerCuspNFW%descriptionParameters                        (1)='The cusp amplitude.'
       ! Record that the table is now initialized.
       containerCuspNFWInitialized                                      =.true.
    end if
    allocate(parameters(1))
    densityNormalization    =  +self%densityNormalization
    radiusNormalization     =  +self%radiusScale
    parameters          (1) =  +self%y
    container               =>       containerCuspNFW
    return
  end subroutine cuspNFWParameters

  subroutine cuspNFWDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionCuspNFW), intent(inout)           :: self
    type     (inputParameters        ), intent(inout)           :: descriptor
    logical                           , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                 )                          :: parameterLabel
    type     (inputParameters        )                          :: parameters
    !$GLC attributes unused :: includeFileModificationTimes
    
    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','cuspNFW')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%radiusScale
    call parameters%addParameter('radiusScale'         ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%y
    call parameters%addParameter('y'                   ,trim(adjustl(parameterLabel)))
    return
  end subroutine cuspNFWDescriptor

  function cuspNFWSuffix(self) result(suffix)
    !!{
    Return a suffix for tabulated file names.
    !!}
    use :: String_Handling, only : String_C_To_Fortran
    implicit none
    type (varying_string         )                :: suffix
    class(massDistributionCuspNFW), intent(inout) :: self
    !$GLC attributes unused :: self

    suffix=String_C_To_Fortran(massDistributionCuspNFWSourceDigest)
    return
  end function cuspNFWSuffix
