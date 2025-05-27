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
   Implementation of a Soliton and NFW \citep{Schive2014} mass distribution class.
   !!}

   !![
   <massDistribution name="massDistributionSolitonNFW">
      <description>
         A Soliton and NFW mass distribution.
      </description>
   </massDistribution>
   !!]

   type, public, extends(massDistributionSphericalTabulated) :: massDistributionSolitonNFW
         !!{
         The Soliton + NFW mass distribution.
         This model consists of a soliton core for small radii and transitions to an NFW profile at larger radii.
         !!}
         private
         double precision :: densityNormalization, radiusScale, radiusCore, radiusSol,     &
                          &  solitonCentralDensity,  concentration, radiusVirial,          &
                          &  radiusCoreScaleFree,  radiusSolScaleFree, densitySolScaleFree
      contains
         procedure :: density               => solitonNFWDensity
         procedure :: densityGradientRadial => solitonNFWDensityGradientRadial
         procedure :: parameters            => solitonNFWParameters
         procedure :: factoryTabulation     => solitonNFWFactoryTabulation
         procedure :: descriptor            => solitonNFWDescriptor
         procedure :: suffix                => solitonNFWSuffix
   end type massDistributionSolitonNFW

   interface massDistributionSolitonNFW
     !!{
     Constructors for the {\normalfont \ttfamily solitonNFW} mass distribution class.
     !!}
     module procedure solitonNFWConstructorParameters
     module procedure solitonNFWConstructorInternal
   end interface massDistributionSolitonNFW

   logical                                     :: containerSolitonNFWInitialized=.false.
   type   (massDistributionContainer), pointer :: containerSolitonNFW
   !$omp threadprivate(containerSolitonNFW,containerSolitonNFWInitialized)

   !![
   <sourceDigest name="massDistributionSolitonNFWSourceDigest"/>
   !!]

contains

   function solitonNFWConstructorParameters(parameters) result(self)
      !!{
      Constructor for the soliton and NFW mass distribution class which builds the object from a parameter set.
      !!}
      use :: Input_Parameters          , only : inputParameter                , inputParameters
      use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
      implicit none
      type          (massDistributionSolitonNFW)               :: self
      type          (inputParameters), intent(inout)           :: parameters
      double precision                                         :: radiusScale,          solitonCentralDensity, &
         &                                                        densityNormalization, radiusCore, radiusSol, &
         &                                                        concentration,        radiusVirial         , &
         &                                                        toleranceRelativePotential
      logical                                                  :: dimensionless
      type(varying_string)                                     :: componentType
      type(varying_string)                                     :: massType

      !![
      <inputParameter>
         <name>radiusScale</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The scale radius of the Soliton and NFW profile (r_s).</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>radiusCore</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The soliton core radius (r_c).</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>radiusSol</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The soliton transition radius (r_sol).</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>solitonCentralDensity</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The central density of the soliton (rho_0).</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>concentration</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The concentration parameter for the NFW profile.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>radiusVirial</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The virial radius of the NFW profile.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>toleranceRelativePotential</name>
         <defaultValue>1.0d-3</defaultValue>
         <description>The relative tolerance for numerical solutions.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>dimensionless</name>
         <defaultValue>.true.</defaultValue>
         <description>If true the profile is dimensionless.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>componentType</name>
         <defaultValue>var_str('unknown')</defaultValue>
         <description>The component type for the profile.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>massType</name>
         <defaultValue>var_str('unknown')</defaultValue>
         <description>The mass type for the profile.</description>
         <source>parameters</source>
      </inputParameter>
      <conditionalCall>
       <call>self = massDistributionSolitonNFW(toleranceRelativePotential=toleranceRelativePotential,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType, includesPrefix=.false.){conditions})</call>
       <argument name="radiusCore"             value="radiusCore"             parameterPresent="parameters"/>
       <argument name="radiusSol"              value="radiusSol"              parameterPresent="parameters"/>
       <argument name="solitonCentralDensity"  value="solitonCentralDensity"  parameterPresent="parameters"/>
       <argument name="densityNormalization"   value="densityNormalization"   parameterPresent="parameters"/>
       <argument name="radiusScale"            value="radiusScale"            parameterPresent="parameters"/>
       <argument name="radiusVirial"           value="radiusVirial"           parameterPresent="parameters"/>
       <argument name="concentration"          value="concentration"          parameterPresent="parameters"/>
       <argument name="dimensionless"          value="dimensionless"          parameterPresent="parameters"/>
       </conditionalCall>
      <inputParametersValidate source="parameters"/>
      !!]
      return
   end function solitonNFWConstructorParameters

   function solitonNFWConstructorInternal(radiusScale, radiusCore, radiusSol, solitonCentralDensity, densityNormalization, radiusVirial, concentration, dimensionless, componentType, massType, toleranceRelativePotential) result(self)
      !!{
      Internal constructor for ``soliton and NFW'' mass distribution class.
      !!}
      use :: Error                   , only : Error_Report
      implicit none
      type      (massDistributionSolitonNFW) :: self
      double precision                  , intent(in), optional :: radiusScale          , radiusCore                , &
            &                                                     solitonCentralDensity, densityNormalization      , &
            &                                                     radiusSol            , radiusVirial              , &
            &                                                     concentration        , toleranceRelativePotential
      logical                           , intent(in), optional :: dimensionless
      type(enumerationComponentTypeType), intent(in), optional :: componentType
      type(enumerationMassTypeType)     , intent(in), optional :: massType
      !![
      <constructorAssign variables="componentType, massType, toleranceRelativePotential"/>
      !!]

      if (present(radiusCore)) then
         self%radiusCore = radiusCore
      else
         call Error_Report('no means to determine core radius')
      end if

      if (present(radiusSol)) then
         self%radiusSol = radiusSol
      else
         call Error_Report('no means to determine Soliton radius')
      end if

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
          call Error_Report('no means to determine scale radius')
      end if

      if (present(densityNormalization)) then
         self%densityNormalization = densityNormalization
      else
         call Error_Report('densityNormalization must be specified')
      end if

      if (present(solitonCentralDensity)) then
         self%solitonCentralDensity = solitonCentralDensity
      else
         call Error_Report('solitonCentralDensity must be specified')
      end if

      if (present(dimensionless)) then
         self%dimensionless = dimensionless
      else
         self%dimensionless = .false.
      end if

      self%radiusCoreScaleFree = self%radiusCore / self%radiusScale
      self%radiusSolScaleFree  = self%radiusSol  / self%radiusScale
      self%densitySolScaleFree = self%solitonCentralDensity / self%densityNormalization
      return
   end function solitonNFWConstructorInternal

   function solitonNFWFactoryTabulation(self,parameters) result(instance)
      !!{
      Construct an instance of this class using tabulation parameters.
      !!}
      implicit none
      class           (massDistributionSphericalTabulated)    , pointer                  :: instance
      class           (massDistributionSolitonNFW)            , intent(inout)            :: self
      double precision                                        , intent(in), dimension(:) :: parameters
      !$GLC attributes unused :: self
    
      allocate(massDistributionSolitonNFW :: instance)
      select type(instance)
      type is (massDistributionSolitonNFW)
         instance = massDistributionSolitonNFW(             &
         radiusScale               = 1.0d0                , &
         densityNormalization      = 1.0d0                , &
         radiusCore                = parameters(1)        , &
         radiusSol                 = parameters(2)        , &
         solitonCentralDensity     = parameters(3)        , &
         toleranceRelativePotential=self%toleranceRelativePotential)
      end select
      return
   end function solitonNFWFactoryTabulation

   double precision function solitonNFWDensity(self,coordinates) result(density)
      !!{
      Return the density at the specified {\normalfont \ttfamily coordinates} in a Soliton and NFW mass distribution.
      !!}
      implicit none
      class           (massDistributionSolitonNFW) ,intent(inout) :: self
      class           (coordinate              )   ,intent(in   ) :: coordinates
      double precision                                            :: radiusScaleFree
      double precision                                            :: radiusCoreFree

      radiusScaleFree = +coordinates%rSpherical() / self%radiusScale
      radiusCoreFree  = +coordinates%rSpherical() / self%radiusCore

      if (coordinates%rSpherical() < self%radiusSol) then
         density = self%solitonCentralDensity / (1.0d0 + 0.091d0 * radiusCoreFree**2)**8
      else
         density = self%densityNormalization / radiusScaleFree / (1.0d0 + radiusScaleFree)**2
      end if
      return
   end function solitonNFWDensity

   double precision function solitonNFWDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
      !!{
      Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a Soliton and NFW mass distribution.
      !!}
      implicit none
      class           (massDistributionSolitonNFW) , intent(inout), target   :: self
      class           (coordinate              )   , intent(in   )           :: coordinates
      logical                                      , intent(in   ), optional :: logarithmic
      double precision                                                       :: radiusScaleFree
      double precision                                                       :: radiusCoreFree
      !![
      <optionalArgument name="logarithmic" defaultsTo=".false."/>
      !!]

      radiusScaleFree = +coordinates%rSpherical() / self%radiusScale
      radiusCoreFree  = +coordinates%rSpherical() / self%radiusCore

      if (coordinates%rSpherical() < self%radiusSol) then
         if (logarithmic) then
            densityGradient = -16.0d0 * 0.091d0 * radiusScaleFree**2 / (1.0d0 + 0.091d0 * radiusCoreFree**2)
         else
            densityGradient = -16.0d0 * 0.091d0 * self%solitonCentralDensity * radiusScaleFree/ self%radiusScale / (1.0d0 + 0.091d0 * radiusCoreFree**2)**9
         end if
      else
         if (logarithmic) then
            densityGradient = -(2.0d0 + radiusScaleFree) / (radiusScaleFree * (1.0d0 + radiusScaleFree))
         else
            densityGradient = -self%densityNormalization / (radiusScaleFree**2 * (1.0d0 + radiusScaleFree)**3)
         end if
      end if
      return
   end function solitonNFWDensityGradientRadial

   subroutine solitonNFWParameters(self,densityNormalization,radiusNormalization,parameters,container)
      !!{
      Establish parameters for tabulation.
      !!}
      implicit none
      class       (massDistributionSolitonNFW), intent(inout)                            :: self
      double precision                        , intent(out)                              :: densityNormalization, radiusNormalization
      double precision                        , intent(inout), allocatable, dimension(:) :: parameters
      type         (massDistributionContainer), intent(out)  , pointer                   :: container

      if (.not.containerSolitonNFWInitialized) then
         allocate(containerSolitonNFW)
         call containerSolitonNFW%initialize(3)
         containerSolitonNFW%mass                      %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%mass                      %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%potential                 %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%potential                 %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%velocityDispersion1D      %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%velocityDispersion1D      %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%energy                    %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%energy                    %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%radiusFreefall            %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%radiusFreefall            %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%radiusFreefallIncreaseRate%radiusCountPer       =+20_c_size_t
         containerSolitonNFW%radiusFreefallIncreaseRate%parametersCountPer   =+20_c_size_t
         containerSolitonNFW%densityRadialMoment0      %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%densityRadialMoment0      %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%densityRadialMoment1      %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%densityRadialMoment1      %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%densityRadialMoment2      %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%densityRadialMoment2      %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%densityRadialMoment3      %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%densityRadialMoment3      %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%fourierTransform          %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%fourierTransform          %parametersCountPer   =+20_c_size_t
         
         containerSolitonNFW%nameParameters        (1) = 'radiusCoreOverradiusScale'
         containerSolitonNFW%descriptionParameters (1) = 'The ratio of core to scale radii(r_c/r_s).'

         containerSolitonNFW%nameParameters        (2) = 'radiusSolitonOverradiusScale'
         containerSolitonNFW%descriptionParameters (2) = 'The ratio of soliton to scale radii(r_sol/r_s).'

         containerSolitonNFW%nameParameters        (3) = 'densitySolitonCentralOverdensityScale'
         containerSolitonNFW%descriptionParameters (3) = 'The ratio of soliton central to scale densities(rho_0/rho_s).'
         
         containerSolitonNFWInitialized                                      = .true.
      end if
      allocate(parameters(3))
      densityNormalization = self%densityNormalization
      radiusNormalization  = self%radiusScale
      parameters(1)        = self%radiusCoreScaleFree
      parameters(2)        = self%radiusSolScaleFree
      parameters(3)        = self%densitySolScaleFree
      container            => containerSolitonNFW
      return
   end subroutine solitonNFWParameters

   subroutine solitonNFWDescriptor(self, descriptor, includeClass, includeFileModificationTimes)
      !!{
      Return an input parameter list descriptor which could be used to recreate this object.
      !!}
      use :: Input_Parameters, only : inputParameters
      implicit none
      class    (massDistributionSolitonNFW), intent(inout)           :: self
      type     (inputParameters           ), intent(inout)           :: descriptor
      logical                              , intent(in   ), optional :: includeClass, includeFileModificationTimes
      character(len=18)                                              :: parameterLabel
      type     (inputParameters           )                          :: parameters
      !$GLC attributes unused :: includeFileModificationTimes

      if (.not.present(includeClass) .or. includeClass) call descriptor%addParameter('massDistribution','solitonNFW')
      parameters = descriptor%subparameters('massDistribution')
      write(parameterLabel,'(e17.10)') self%densityNormalization
      call parameters%addParameter('densityNormalization', trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%radiusScale
      call parameters%addParameter('radiusScale', trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%radiusCore
      call parameters%addParameter('radiusCore', trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%radiusSol
      call parameters%addParameter('radiusSol', trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%solitonCentralDensity
      call parameters%addParameter('solitonCentralDensity', trim(adjustl(parameterLabel)))
      return
   end subroutine solitonNFWDescriptor

   function solitonNFWSuffix(self) result(suffix)
      !!{
      Return a suffix for tabulated file names.
      !!}
      use :: String_Handling, only : String_C_To_Fortran
      implicit none
      type (varying_string          )                  :: suffix
      class(massDistributionSolitonNFW), intent(inout) :: self
      !$GLC attributes unused :: self

      suffix = String_C_To_Fortran(massDistributionSolitonNFWSourceDigest)
      return
   end function solitonNFWSuffix
