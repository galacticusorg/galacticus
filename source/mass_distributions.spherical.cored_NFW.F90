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
  Implementation of a cored NFW \citep{navarro_structure_1996} mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionCoredNFW">
    <description>
      A cored NFW \citep{navarro_structure_1996} mass distribution class. The density profile is given by:
      \begin{equation}
       \rho_\mathrm{dark matter}(r) \propto \left({r_\mathrm{c}\over r_\mathrm{s}}+{r\over r_\mathrm{s}}\right)^{-1} \left[1 + \left({r\over r_\mathrm{s}}\right)
      \right]^{-2}.
      \end{equation}
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSphericalTabulated) :: massDistributionCoredNFW
     !!{
     The cored NFW \citep{navarro_structure_1996} mass distribution.
     !!}
     private
     double precision :: densityNormalization, radiusScale        , &
          &              radiusCore          , radiusCoreScaleFree
   contains
     procedure :: density               => coredNFWDensity
     procedure :: densityGradientRadial => coredNFWDensityGradientRadial
     procedure :: parameters            => coredNFWParameters
     procedure :: factoryTabulation     => coredNFWFactoryTabulation
     procedure :: descriptor            => coredNFWDescriptor
     procedure :: suffix                => coredNFWSuffix
  end type massDistributionCoredNFW
  
  interface massDistributionCoredNFW
     !!{
     Constructors for the {\normalfont \ttfamily coredNFW} mass distribution class.
     !!}
     module procedure coredNFWConstructorParameters
     module procedure coredNFWConstructorInternal
  end interface massDistributionCoredNFW

  ! Tabulated solutions.
  logical                                     :: containerCoredNFWInitialized=.false.
  type   (massDistributionContainer), pointer :: containerCoredNFW
  !$omp threadprivate(containerCoredNFW,containerCoredNFWInitialized)
  
  ! Generate a source digest.
  !![
  <sourceDigest name="massDistributionCoredNFWSourceDigest"/>
  !!]

contains

  function coredNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily coredNFW} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    type            (massDistributionCoredNFW)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    double precision                                          :: mass                      , radiusScale  , &
         &                                                       densityNormalization      , concentration, &
         &                                                       radiusVirial              , radiusCore   , &
         &                                                       toleranceRelativePotential
    logical                                                   :: dimensionless
    type            (varying_string          )                :: componentType
    type            (varying_string          )                :: massType

    !![
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>1.0d0/2.0d0/Pi/(log(4.0d0)-1.0d0)</defaultValue>
      <description>The density normalization of the cored NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusScale</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the cored NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusCore</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The core radius of the cored NFW profile.</description>
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
      <description>The concentration of the cored NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusVirial</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The virial radius of the cored NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the cored NFW profile is considered to be dimensionless.</description>
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
     <call>self=massDistributionCoredNFW(radiusCore=radiusCore,toleranceRelativePotential=toleranceRelativePotential,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
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
  end function coredNFWConstructorParameters

  function coredNFWConstructorInternal(radiusScale,radiusCore,concentration,densityNormalization,mass,radiusVirial,dimensionless,componentType,massType,toleranceRelativePotential) result(self)
    !!{
    Internal constructor for ``coreNFW'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionCoredNFW    )                          :: self
    double precision                              , intent(in   ), optional :: radiusScale               , concentration      , &
         &                                                                     densityNormalization      , mass               , &
         &                                                                     radiusVirial              , radiusCore         , &
         &                                                                     toleranceRelativePotential
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree           , radiusCoreScaleFree
    !![
    <constructorAssign variables="componentType, massType, toleranceRelativePotential"/>
    !!]

    ! Determine core radius.
    if      (                            &
         &   present(radiusCore   )      &
         &  ) then
       self%radiusCore          =radiusCore
    else
       call Error_Report('no means to determine core radius' //{introspection:location})
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
       self%radiusScale=radiusVirial/concentration
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
       radiusCoreScaleFree      =+self%radiusCore  /self%radiusScale
       self%densityNormalization=+mass                                                                                                    &
            &                    /self%radiusScale**3                                                                                     &
            &                    /4.0d0                                                                                                   &
            &                    /Pi                                                                                                      &
            &                    *(1.0d0-radiusCoreScaleFree)**2                                                                          &
            &                    /(                                                                                                       &
            &                      +(1.0d0-2.0d0*radiusCoreScaleFree   )*log(  +1.0d0              +radiusScaleFree                     ) &
            &                      +             radiusCoreScaleFree**2 *log(+(+radiusCoreScaleFree+radiusScaleFree)/radiusCoreScaleFree) &
            &                      -(1.0d0-      radiusCoreScaleFree   )*radiusScaleFree/(1.0d0+radiusScaleFree)                          &
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
    ! Compute scale-free core radius.
    self%radiusCoreScaleFree=+self%radiusCore  &
         &                   /self%radiusScale
    return
  end function coredNFWConstructorInternal

  function coredNFWFactoryTabulation(self,parameters) result(instance)
    !!{
    Construct an instance of this class using tabulation parameters.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated)               , pointer      :: instance
    class           (massDistributionCoredNFW          ), intent(inout)               :: self
    double precision                                    , intent(in   ), dimension(:) :: parameters
    !$GLC attributes unused :: self
    
    allocate(massDistributionCoredNFW :: instance)
    select type(instance)
    type is (massDistributionCoredNFW)
       instance=massDistributionCoredNFW(radiusScale=1.0d0,densityNormalization=1.0d0,radiusCore=parameters(1),toleranceRelativePotential=self%toleranceRelativePotential)
    end select
    return
  end function coredNFWFactoryTabulation
  
  double precision function coredNFWDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a cored NFW mass distribution.
    !!}
    implicit none
    class           (massDistributionCoredNFW), intent(inout) :: self
    class           (coordinate              ), intent(in   ) :: coordinates
    double precision                                          :: radiusScaleFree

    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %radiusScale
    density        =+self       %densityNormalization              &
         &          /(self%radiusCoreScaleFree+radiusScaleFree)    &
         &          /(1.0d0                   +radiusScaleFree)**2
    return
  end function coredNFWDensity
  
  double precision function coredNFWDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a cored NFW mass distribution.
    !!}
    implicit none
    class           (massDistributionCoredNFW), intent(inout), target   :: self
    class           (coordinate              ), intent(in   )           :: coordinates
    logical                                   , intent(in   ), optional :: logarithmic
    double precision                                                    :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %radiusScale
    if (logarithmic_) then
       densityGradient=-     3.0d0                                                          &
            &          +     2.0d0              /(     1.0d0              +radiusScaleFree) &
            &          +self%radiusCoreScaleFree/(self%radiusCoreScaleFree+radiusScaleFree)
    else
       densityGradient=-self%densityNormalization                                        &
            &          /self%radiusScale                                                 &
            &          *(+1.0d0+2.0d0*self%radiusCoreScaleFree+3.0d0*radiusScaleFree)    &
            &          /(      +      self%radiusCoreScaleFree+      radiusScaleFree)**2 &
            &          /(+1.0d0                               +      radiusScaleFree)**3
    end if
    return
  end function coredNFWDensityGradientRadial

  subroutine coredNFWParameters(self,densityNormalization,radiusNormalization,parameters,container)
    !!{
    Establish parameters for tabulation.
    !!}
    implicit none
    class           (massDistributionCoredNFW ), intent(inout)                              :: self
    double precision                           , intent(  out)                              :: densityNormalization, radiusNormalization
    double precision                           , intent(inout), allocatable, dimension(:  ) :: parameters
    type            (massDistributionContainer), intent(  out), pointer                     :: container

    if (.not.containerCoredNFWInitialized) then
       ! Allocate the table and initialize.
       allocate(containerCoredNFW)
       call containerCoredNFW%initialize(1)
       ! Specify the number of tabulation points per interval in radius and each parameter.
       containerCoredNFW%mass                      %radiusCountPer       =+20_c_size_t
       containerCoredNFW%mass                      %parametersCountPer   =+20_c_size_t
       containerCoredNFW%potential                 %radiusCountPer       =+20_c_size_t
       containerCoredNFW%potential                 %parametersCountPer   =+20_c_size_t
       containerCoredNFW%velocityDispersion1D      %radiusCountPer       =+20_c_size_t
       containerCoredNFW%velocityDispersion1D      %parametersCountPer   =+20_c_size_t
       containerCoredNFW%energy                    %radiusCountPer       =+20_c_size_t
       containerCoredNFW%energy                    %parametersCountPer   =+20_c_size_t
       containerCoredNFW%radiusFreefall            %radiusCountPer       =+20_c_size_t
       containerCoredNFW%radiusFreefall            %parametersCountPer   =+20_c_size_t
       containerCoredNFW%radiusFreefallIncreaseRate%radiusCountPer       =+20_c_size_t
       containerCoredNFW%radiusFreefallIncreaseRate%parametersCountPer   =+20_c_size_t
       containerCoredNFW%densityRadialMoment0      %radiusCountPer       =+20_c_size_t
       containerCoredNFW%densityRadialMoment0      %parametersCountPer   =+20_c_size_t
       containerCoredNFW%densityRadialMoment1      %radiusCountPer       =+20_c_size_t
       containerCoredNFW%densityRadialMoment1      %parametersCountPer   =+20_c_size_t
       containerCoredNFW%densityRadialMoment2      %radiusCountPer       =+20_c_size_t
       containerCoredNFW%densityRadialMoment2      %parametersCountPer   =+20_c_size_t
       containerCoredNFW%densityRadialMoment3      %radiusCountPer       =+20_c_size_t
       containerCoredNFW%densityRadialMoment3      %parametersCountPer   =+20_c_size_t
       containerCoredNFW%fourierTransform          %radiusCountPer       =+20_c_size_t
       containerCoredNFW%fourierTransform          %parametersCountPer   =+20_c_size_t
       ! Specify names and descriptions of the parameters.
       containerCoredNFW%nameParameters                               (1)='radiusCoreOverRadiusScale'
       containerCoredNFW%descriptionParameters                        (1)='The ratio of core to scale radii.'
       ! Record that the table is now initialized.
       containerCoredNFWInitialized                                =.true.
    end if
    allocate(parameters(1))
    densityNormalization    =  +self%densityNormalization
    radiusNormalization     =  +self%radiusScale
    parameters          (1) =  +self%radiusCoreScaleFree
    container               =>       containerCoredNFW
    return
  end subroutine coredNFWParameters

  subroutine coredNFWDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionCoredNFW), intent(inout)           :: self
    type     (inputParameters         ), intent(inout)           :: descriptor
    logical                            , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                  )                          :: parameterLabel
    type     (inputParameters         )                          :: parameters
    !$GLC attributes unused :: includeFileModificationTimes
    
    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','coredNFW')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%radiusScale
    call parameters%addParameter('radiusScale'         ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%radiusCore
    call parameters%addParameter('radiusCore'          ,trim(adjustl(parameterLabel)))
    return
  end subroutine coredNFWDescriptor

  function coredNFWSuffix(self) result(suffix)
    !!{
    Return a suffix for tabulated file names.
    !!}
    use :: String_Handling, only : String_C_To_Fortran
    implicit none
    type (varying_string          )                :: suffix
    class(massDistributionCoredNFW), intent(inout) :: self
    !$GLC attributes unused :: self

    suffix=String_C_To_Fortran(massDistributionCoredNFWSourceDigest)
    return
  end function coredNFWSuffix
