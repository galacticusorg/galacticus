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

  !+    Contributions to this file made by: Yu Zhao

  !!{
  Implementation of a mass distribution class for fuzzy dark matter halos consisting of soliton profiles \citep{schive_understanding_2014}.   
  !!}
  
  !![
  <massDistribution name="massDistributionSoliton">
    <description>
      A mass distribution class for fuzzy dark matter halos \citep{schive_understanding_2014} consisting of soliton core.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSphericalTabulated) :: massDistributionSoliton
     !!{
     A mass distribution class for fuzzy dark matter halos \citep{schive_understanding_2014} consisting of soliton core.     
     !!}
     private
     double precision :: radiusCore, densitySolitonCentral
   contains
     procedure :: massEnclosedBySphere   => solitonMassEnclosedBySphere
     procedure :: density                => solitonDensity
     procedure :: densityGradientRadial  => solitonDensityGradientRadial
     procedure :: radiusEnclosingDensity => solitonRadiusEnclosingDensity
     procedure :: parameters             => solitonParameters
     procedure :: factoryTabulation      => solitonFactoryTabulation
     procedure :: descriptor             => solitonDescriptor
     procedure :: suffix                 => solitonSuffix
  end type massDistributionSoliton
  
   interface massDistributionSoliton
     !!{
     Constructors for the {\normalfont \ttfamily soliton} mass distribution class.
     !!}
     module procedure solitonConstructorParameters
     module procedure solitonConstructorInternal
   end interface massDistributionSoliton

   logical                                       :: containerSolitonInitialized=.false.
   type   (massDistributionContainer), pointer   :: containerSoliton
   !$omp threadprivate(containerSoliton,containerSolitonInitialized)

   ! Coefficient of the dimensionless radius in the soliton profile.
   double precision                  , parameter :: coefficientCore               =0.091d0 ! Schive et al. (2014; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S; equation 3).

   !![
   <sourceDigest name="massDistributionSolitonSourceDigest"/>
   !!]

contains

   function solitonConstructorParameters(parameters) result(self)
      !!{
      Constructor for the soliton mass distribution class which builds the object from a parameter set.
      !!}
      use :: Input_Parameters          , only : inputParameter                , inputParameters
      use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
      implicit none
      type            (massDistributionSoliton)                :: self
      type            (inputParameters        ), intent(inout) :: parameters
      double precision                                         :: radiusCore                , &
         &                                                        densitySolitonCentral     , &
         &                                                        toleranceRelativePotential
      logical                                                  :: dimensionless
      type            (varying_string         )                :: componentType
      type            (varying_string         )                :: massType

      !![
      <inputParameter>
         <name>radiusCore</name>
         <description>The soliton core radius.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>densitySolitonCentral</name>
         <description>The central density of the soliton.</description>
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
	<call>
	  self=massDistributionSoliton(toleranceRelativePotential=toleranceRelativePotential,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})
	</call>
       <argument name="radiusCore"              value="radiusCore"              parameterPresent="parameters"/>
       <argument name="densitySolitonCentral"   value="densitySolitonCentral"   parameterPresent="parameters"/>
       <argument name="dimensionless"           value="dimensionless"           parameterPresent="parameters"/>
      </conditionalCall>
      <inputParametersValidate source="parameters"/>
      !!]
      return
   end function solitonConstructorParameters

   function solitonConstructorInternal(radiusCore,densitySolitonCentral,dimensionless,componentType,massType,toleranceRelativePotential) result(self)
     !!{
     Internal constructor for ``soliton'' mass distribution class.
     !!}
     use :: Error, only : Error_Report
     implicit none
     type            (massDistributionSoliton     )                       :: self
     double precision                              , intent(in), optional :: radiusCore                , &
          &                                                                  densitySolitonCentral     , &
          &                                                                  toleranceRelativePotential
     logical                                       , intent(in), optional :: dimensionless
     type            (enumerationComponentTypeType), intent(in), optional :: componentType
     type            (enumerationMassTypeType     ), intent(in), optional :: massType
     !![
     <constructorAssign variables="componentType, massType, toleranceRelativePotential"/>
     !!]
     
     if (present(radiusCore   )) then
        self%radiusCore             =+radiusCore
     else
        call Error_Report('no means to determine core radius')
     end if
     if (present(densitySolitonCentral)) then
        self%densitySolitonCentral  =+densitySolitonCentral
     else
        call Error_Report('densitySolitonCentral must be specified')
     end if
     if (present(dimensionless)) then
        self%dimensionless          =dimensionless
     else
        self%dimensionless          =.false.
     end if
     return
   end function solitonConstructorInternal

   function solitonFactoryTabulation(self,parameters) result(instance)
      !!{
      Construct an instance of this class using tabulation parameters.
      !!}
      implicit none
      class           (massDistributionSphericalTabulated), pointer                     :: instance
      class           (massDistributionSoliton           ), intent(inout)               :: self
      double precision                                    , intent(in   ), dimension(:) :: parameters
    
      allocate(massDistributionSoliton :: instance)
      select type(instance)
      type is (massDistributionSoliton)
         instance= massDistributionSoliton(                                                             &
              &                             radiusCore                =     1.0d0                     , &
              &                             densitySolitonCentral     =     1.0d0                     , &
              &                             toleranceRelativePotential=self%toleranceRelativePotential  &
              &                              )
      end select
      return
   end function solitonFactoryTabulation

   double precision function solitonMassEnclosedBySphere(self,radius) result(mass)
      !!{
      Return the mass at the specified {\normalfont \ttfamily coordinates} in a soliton mass distribution.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      class           (massDistributionSoliton), intent(inout), target :: self
      double precision                         , intent(in   )         :: radius
      double precision                         , parameter             :: fractionRadiusSmall=0.1d0
      double precision                                                 :: termArcTangent           , termPrefactor , &
             &                                                            termCore                 , termPolynomial

      if (radius <= 0.0d0) then
         ! Zero radius.
         mass           =+0.0d0
      else if (radius < fractionRadiusSmall*self%radiusCore) then
           ! At small radii use a Taylor series for numerical accuracy.
         mass          =+  4.0d0/3.0d0*Pi                   *self%densitySolitonCentral*radius**3                    &
              &         - 32.0d0/5.0d0*Pi*coefficientCore   *self%densitySolitonCentral*radius**5/self%radiusCore**2 &
              &         +144.0d0/7.0d0*Pi*coefficientCore**2*self%densitySolitonCentral*radius**7/self%radiusCore**4
      else
         ! Use the exact solution.
         termCore      =+sqrt(coefficientCore)                                        &
              &         *                 radius                                      &
              &         *                            self%radiusCore                  &
              &         /(coefficientCore*radius**2+self%radiusCore**2)**7
         termPolynomial=+  3465.0d0*coefficientCore**6*radius**12                     &
              &         + 23100.0d0*coefficientCore**5*radius**10*self%radiusCore** 2 &
              &         + 65373.0d0*coefficientCore**4*radius** 8*self%radiusCore** 4 &
              &         +101376.0d0*coefficientCore**3*radius** 6*self%radiusCore** 6 &
              &         + 92323.0d0*coefficientCore**2*radius** 4*self%radiusCore** 8 &
              &         + 48580.0d0*coefficientCore*   radius** 2*self%radiusCore**10 &
              &         -  3465.0d0                               *self%radiusCore**12
         termArcTangent=+3465.0d0                                                     &
              &         *atan(sqrt(coefficientCore)*radius/self%radiusCore)
         termPrefactor =+Pi                                                           &
              &         *self%densitySolitonCentral                                   &
              &         *self%radiusCore           **3                                &
              &         /(                                                            &
              &           +53760.0d0                                                  &
              &           *coefficientCore**1.5d0                                     &
              &         )
         mass          = + termPrefactor                                              &
              &         *(                                                            &
              &           +termArcTangent                                             &
              &           +termCore                                                   &
              &           *termPolynomial                                             &
              &          )
      end if
      return
    end function solitonMassEnclosedBySphere

    double precision function solitonDensity(self,coordinates) result(density)
      !!{
      Return the density at the specified {\normalfont \ttfamily coordinates} in a soliton mass distribution.
      !!}
      implicit none
      class           (massDistributionSoliton), intent(inout) :: self
      class           (coordinate             ), intent(in   ) :: coordinates
      double precision                                         :: radiusCoreFree
      
      radiusCoreFree =+coordinates%rSpherical()/self%radiusCore
      
      if (coordinates%rSpherical() <= 0.0d0) then
         ! Zero radius.
         density=+self%densitySolitonCentral
      else
         density=+self%densitySolitonCentral/(+1.0d0+coefficientCore*radiusCoreFree**2)**8
      end if
      return
   end function solitonDensity

   double precision function solitonDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
      !!{
      Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a soliton mass distribution.
      !!}
      implicit none
      class           (massDistributionSoliton) , intent(inout), target   :: self
      class           (coordinate             ) , intent(in   )           :: coordinates
      logical                                   , intent(in   ), optional :: logarithmic
      double precision                                                    :: radiusCoreFree
      !![
      <optionalArgument name="logarithmic" defaultsTo=".false."/>
      !!]

      radiusCoreFree =+coordinates%rSpherical()/self%radiusCore

      if (coordinates%rSpherical() <= 0.0d0) then
         ! Zero radius.
         densityGradient=+0.0d0
      else if (logarithmic) then
         densityGradient=-16.0d0                                       &
                 &       *       coefficientCore*radiusCoreFree**2     &
                 &       /(1.0d0+coefficientCore*radiusCoreFree**2)
      else
         densityGradient=-16.0d0                                       &
                 &       *self%densitySolitonCentral                   &
                 &       /self%radiusCore                              &
                 &       *       coefficientCore*radiusCoreFree        &
                 &       /(1.0d0+coefficientCore*radiusCoreFree**2)**9
      end if
      return
   end function solitonDensityGradientRadial

   double precision function solitonRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for soliton mass distributions.
    !!}
    implicit none
    class           (massDistributionSoliton), intent(inout), target   :: self
    double precision                         , intent(in   )           :: density
    double precision                         , intent(in   ), optional :: radiusGuess

    if (density >= self%densitySolitonCentral) then
       radius=0.0d0
    else
       radius=sphericalTabulatedRadiusEnclosingDensity(self,density,radiusGuess)
    end if
    return
  end function solitonRadiusEnclosingDensity

   subroutine solitonParameters(self,densityNormalization,radiusNormalization,parameters,container)
      !!{
      Establish parameters for tabulation.
      !!}
      implicit none
      class           (massDistributionSoliton   ), intent(inout)                            :: self
      double precision                            , intent(  out)                            :: densityNormalization, radiusNormalization
      double precision                            , intent(inout), allocatable, dimension(:) :: parameters
      type             (massDistributionContainer), intent(  out), pointer                   :: container

      if (.not.containerSolitonInitialized) then
         allocate(containerSoliton)
         call containerSoliton%initialize(0)
         containerSoliton%mass                      %radiusCountPer       =+20_c_size_t
         containerSoliton%mass                      %parametersCountPer   =+20_c_size_t
         containerSoliton%radiusEnclosingDensity    %radiusCountPer       =+20_c_size_t
         containerSoliton%radiusEnclosingDensity    %parametersCountPer   =+20_c_size_t
         containerSoliton%potential                 %radiusCountPer       =+20_c_size_t
         containerSoliton%potential                 %parametersCountPer   =+20_c_size_t
         containerSoliton%velocityDispersion1D      %radiusCountPer       =+20_c_size_t
         containerSoliton%velocityDispersion1D      %parametersCountPer   =+20_c_size_t
         containerSoliton%energy                    %radiusCountPer       =+20_c_size_t
         containerSoliton%energy                    %parametersCountPer   =+20_c_size_t
         containerSoliton%radiusFreefall            %radiusCountPer       =+20_c_size_t
         containerSoliton%radiusFreefall            %parametersCountPer   =+20_c_size_t
         containerSoliton%radiusFreefallIncreaseRate%radiusCountPer       =+20_c_size_t
         containerSoliton%radiusFreefallIncreaseRate%parametersCountPer   =+20_c_size_t
         containerSoliton%densityRadialMoment0      %radiusCountPer       =+20_c_size_t
         containerSoliton%densityRadialMoment0      %parametersCountPer   =+20_c_size_t
         containerSoliton%densityRadialMoment1      %radiusCountPer       =+20_c_size_t
         containerSoliton%densityRadialMoment1      %parametersCountPer   =+20_c_size_t
         containerSoliton%densityRadialMoment2      %radiusCountPer       =+20_c_size_t
         containerSoliton%densityRadialMoment2      %parametersCountPer   =+20_c_size_t
         containerSoliton%densityRadialMoment3      %radiusCountPer       =+20_c_size_t
         containerSoliton%densityRadialMoment3      %parametersCountPer   =+20_c_size_t
         containerSoliton%fourierTransform          %radiusCountPer       =+20_c_size_t
         containerSoliton%fourierTransform          %parametersCountPer   =+20_c_size_t
         containerSolitonInitialized                                      =.true.
      end if
      densityNormalization =  self%densitySolitonCentral
      radiusNormalization  =  self%radiusCore
      container            => containerSoliton
      return
   end subroutine solitonParameters

   subroutine solitonDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
      !!{
      Return an input parameter list descriptor which could be used to recreate this object.
      !!}
      use :: Input_Parameters, only : inputParameters
      implicit none
      class    (massDistributionSoliton), intent(inout)           :: self
      type     (inputParameters        ), intent(inout)           :: descriptor
      logical                           , intent(in   ), optional :: includeClass  , includeFileModificationTimes
      character(len=18)                                           :: parameterLabel
      type     (inputParameters        )                          :: parameters
      !$GLC attributes unused :: includeFileModificationTimes

      if (.not.present(includeClass) .or. includeClass) call descriptor%addParameter('massDistribution','soliton')
      parameters = descriptor%subparameters('massDistribution')
      write(parameterLabel,'(e17.10)') self%radiusCore
      call parameters%addParameter('radiusCore'             ,trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%densitySolitonCentral
      call parameters%addParameter('densitySolitonCentral'  ,trim(adjustl(parameterLabel)))
      return
   end subroutine solitonDescriptor

   function solitonSuffix(self) result(suffix)
      !!{
      Return a suffix for tabulated file names.
      !!}
      use :: String_Handling, only : String_C_To_Fortran
      implicit none
      type (varying_string         )                :: suffix
      class(massDistributionSoliton), intent(inout) :: self
      !$GLC attributes unused :: self

      suffix=String_C_To_Fortran(massDistributionSolitonSourceDigest)
      return
   end function solitonSuffix
