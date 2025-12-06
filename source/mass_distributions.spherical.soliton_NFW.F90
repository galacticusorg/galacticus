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
  Implementation of a mass distribution class for fuzzy dark matter halos consisting of soliton and NFW profiles
  \citep{schive_understanding_2014}.   
  !!}
  
  !![
  <massDistribution name="massDistributionSolitonNFW">
    <description>
      A mass distribution class for fuzzy dark matter halos \citep{schive_understanding_2014} consisting of soliton core for small
      radii, transitioning to an \gls{nfw} profile at larger radii.
    </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSphericalTabulated) :: massDistributionSolitonNFW
     !!{
     A mass distribution class for fuzzy dark matter halos \citep{schive_understanding_2014} consisting of soliton core for small
     radii, transitioning to an \gls{nfw} profile at larger radii.     
     !!}
     private
     double precision :: densityNormalizationNFW, radiusScale            , &
          &              radiusCore             , radiusSoliton          , &
          &              densitySolitonCentral  , concentration          , &
          &              radiusVirial           , radiusCoreScaleFree    , &
          &              radiusSolitonScaleFree , densitySolitonScaleFree, &
          &              massSolitionTransition , massNFWTransition      , &
          &              densityTransition
   contains
     procedure :: massEnclosedBySphere            => solitonNFWMassEnclosedBySphere
     procedure :: density                         => solitonNFWDensity
     procedure :: densityGradientRadial           => solitonNFWDensityGradientRadial
     procedure :: radiusEnclosingDensity          => solitonNFWRadiusEnclosingDensity
     procedure :: radiusEnclosingDensityNumerical => solitonNFWRadiusEnclosingDensityNumerical
     procedure :: parameters                      => solitonNFWParameters
     procedure :: factoryTabulation               => solitonNFWFactoryTabulation
     procedure :: descriptor                      => solitonNFWDescriptor
     procedure :: suffix                          => solitonNFWSuffix
  end type massDistributionSolitonNFW
  
   interface massDistributionSolitonNFW
     !!{
     Constructors for the {\normalfont \ttfamily solitonNFW} mass distribution class.
     !!}
     module procedure solitonNFWConstructorParameters
     module procedure solitonNFWConstructorInternal
   end interface massDistributionSolitonNFW

   logical                                       :: containerSolitonNFWInitialized=.false.
   type   (massDistributionContainer), pointer   :: containerSolitonNFW
   !$omp threadprivate(containerSolitonNFW,containerSolitonNFWInitialized)

   ! Coefficient of the dimensionless radius in the soliton profile.
   double precision                  , parameter :: coefficientCore               =0.091d0 ! Schive et al. (2014; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S; equation 3).

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
      type            (massDistributionSolitonNFW)                :: self
      type            (inputParameters           ), intent(inout) :: parameters
      double precision                                            :: radiusScale          , radiusCore                , &
         &                                                           densitySolitonCentral, densityNormalizationNFW   , &
         &                                                           radiusSoliton        , radiusVirial              , &
         &                                                           concentration        , toleranceRelativePotential
      logical                                                     :: dimensionless
      type            (varying_string            )                :: componentType
      type            (varying_string            )                :: massType

      !![
      <inputParameter>
         <name>radiusScale</name>
         <description>The scale radius of the NFW component of the mass distribution.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>radiusCore</name>
         <description>The soliton core radius.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>radiusSoliton</name>
         <description>The soliton transition radius.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>densitySolitonCentral</name>
         <description>The central density of the soliton.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>densityNormalizationNFW</name>
         <description>The density normalization of the NFW profile.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>concentration</name>
         <description>The concentration of the NFW profile.</description>
         <source>parameters</source>
      </inputParameter>
      <inputParameter>
         <name>radiusVirial</name>
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
	<call>
	  self=massDistributionSolitonNFW(toleranceRelativePotential=toleranceRelativePotential,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})
	</call>
       <argument name="radiusCore"              value="radiusCore"              parameterPresent="parameters"/>
       <argument name="radiusSoliton"           value="radiusSoliton"           parameterPresent="parameters"/>
       <argument name="densitySolitonCentral"   value="densitySolitonCentral"   parameterPresent="parameters"/>
       <argument name="densityNormalizationNFW" value="densityNormalizationNFW" parameterPresent="parameters"/>
       <argument name="radiusScale"             value="radiusScale"             parameterPresent="parameters"/>
       <argument name="radiusVirial"            value="radiusVirial"            parameterPresent="parameters"/>
       <argument name="concentration"           value="concentration"           parameterPresent="parameters"/>
       <argument name="dimensionless"           value="dimensionless"           parameterPresent="parameters"/>
      </conditionalCall>
      <inputParametersValidate source="parameters"/>
      !!]
      return
   end function solitonNFWConstructorParameters

   function solitonNFWConstructorInternal(radiusScale,radiusCore,radiusSoliton,densitySolitonCentral,densityNormalizationNFW,radiusVirial,concentration,dimensionless,componentType,massType,toleranceRelativePotential) result(self)
     !!{
     Internal constructor for ``soliton and NFW'' mass distribution class.
     !!}
     use :: Error      , only : Error_Report
     use :: Coordinates, only : coordinateSpherical, assignment(=)
     implicit none
     type            (massDistributionSolitonNFW)                         :: self
     double precision                              , intent(in), optional :: radiusScale          , radiusCore                , &
          &                                                                  densitySolitonCentral, densityNormalizationNFW   , &
          &                                                                  radiusSoliton        , radiusVirial              , &
          &                                                                  concentration        , toleranceRelativePotential
     logical                                       , intent(in), optional :: dimensionless
     type            (enumerationComponentTypeType), intent(in), optional :: componentType
     type            (enumerationMassTypeType     ), intent(in), optional :: massType
     type            (coordinateSpherical         )                       :: coordinates
     !![
     <constructorAssign variables="componentType, massType, toleranceRelativePotential"/>
     !!]
     
     if (present(radiusCore   )) then
        self%radiusCore             =+radiusCore
     else
        self%radiusCore             =+0.0d0
        call Error_Report('no means to determine core radius')
     end if
     if (present(radiusSoliton)) then
        self%radiusSoliton          =+radiusSoliton
     else
        self%radiusSoliton          =+0.0d0
        call Error_Report('no means to determine Soliton radius')
     end if
     if      (                            &
          &   present(radiusScale  )      &
          &  ) then
        self%radiusScale            =+radiusScale
     else if (                            &
          &   present(concentration).and. &
          &   present(radiusVirial )      &
          &  ) then
        self%radiusScale            =+radiusVirial  &
             &                       /concentration
     else
        self%radiusScale            =+0.0d0
        call Error_Report('no means to determine scale radius')
     end if
     if (present(densityNormalizationNFW)) then
        self%densityNormalizationNFW=+densityNormalizationNFW
     else
        self%densityNormalizationNFW=+0.0d0
        call Error_Report('densityNormalizationNFW must be specified')
     end if
     if (present(densitySolitonCentral)) then
        self%densitySolitonCentral  =+densitySolitonCentral
     else
        self%densitySolitonCentral  =+0.0d0
        call Error_Report('densitySolitonCentral must be specified')
     end if
     if (present(dimensionless)) then
        self%dimensionless          =dimensionless
     else
        self%dimensionless          =.false.
     end if
     self%radiusCoreScaleFree    =+self%radiusCore           /self%radiusScale
     self%radiusSolitonScaleFree =+self%radiusSoliton        /self%radiusScale
     self%densitySolitonScaleFree=+self%densitySolitonCentral/self%densityNormalizationNFW
     ! Determine the masses of NFW and soliton profiles at the soliton radius.
     coordinates                 =[self%radiusSoliton,0.0d0,0.0d0]
     self%massSolitionTransition =+self%massEnclosedBySphere(self%radiusSoliton)
     self%massNFWTransition      =+4.0d0                                     &
          &                       *Pi                                        &
          &                       *self%densityNormalizationNFW              &
          &                       *self%radiusScale**3                       &
          &                       *(                                         &
          &                         +log(+1.0d0+self%radiusSolitonScaleFree) &
          &                         -           self%radiusSolitonScaleFree  &
          &                         /   (+1.0d0+self%radiusSolitonScaleFree) &
          &                       )
     self%densityTransition      =+self%density(coordinates)
     return
   end function solitonNFWConstructorInternal

   function solitonNFWFactoryTabulation(self,parameters) result(instance)
      !!{
      Construct an instance of this class using tabulation parameters.
      !!}
      implicit none
      class           (massDistributionSphericalTabulated), pointer                  :: instance
      class           (massDistributionSolitonNFW)        , intent(inout)            :: self
      double precision                                    , intent(in), dimension(:) :: parameters
    
      allocate(massDistributionSolitonNFW :: instance)
      select type(instance)
      type is (massDistributionSolitonNFW)
         instance= massDistributionSolitonNFW(                                                            &
              &                               radiusScale               =     1.0d0                     , &
              &                               densityNormalizationNFW   =     1.0d0                     , &
              &                               radiusCore                =     parameters(1)             , &
              &                               radiusSoliton             =     parameters(2)             , &
              &                               densitySolitonCentral     =     parameters(3)             , &
              &                               toleranceRelativePotential=self%toleranceRelativePotential  &
              &                              )
      end select
      return
   end function solitonNFWFactoryTabulation

   double precision function solitonNFWMassEnclosedBySphere(self,radius) result(mass)
      !!{
      Return the mass at the specified {\normalfont \ttfamily coordinates} in a soliton and NFW mass distribution.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      class           (massDistributionSolitonNFW), intent(inout), target :: self
      double precision                            , intent(in   )         :: radius
      double precision                                                    :: termLogarithmic, termInverse

      if      (radius <= 0.0d0             ) then
         ! Zero radius.
         mass           =+0.0d0
      else if (radius <= self%radiusSoliton) then
         ! Soliton regime.
         mass           =+massSoliton(radius)
      else
         ! NFW regime.
         termInverse    =+self%radiusScale                             &
              &          *(                                            &
              &           +1.0d0/(     radius     +self%radiusScale  ) &
              &           -1.0d0/(self%radiusScale+self%radiusSoliton) &
              &          )
         termLogarithmic=+log(     radius     +self%radiusScale  ) &
              &          -log(self%radiusScale+self%radiusSoliton)
         mass           =+massSoliton(self%radiusSoliton) &
              &          +4.0d0                           &
              &          *Pi                              &
              &          *self%densityNormalizationNFW    &
              &          *self%radiusScale**3             &
              &          *(                               &
              &            +termInverse                   &
              &            +termLogarithmic               &
              &          )
      end if
      return

    contains

      double precision function massSoliton(radius_) result(mass)
        !!{
        Return the mass at the specified {\normalfont \ttfamily coordinates} in the soliton region.
        !!}
        implicit none
        double precision, intent(in   ) :: radius_
        double precision, parameter     :: fractionRadiusSmall=0.1d0
        double precision                :: termArcTangent           , termPrefactor , &
             &                             termCore                 , termPolynomial

        if (radius_ < fractionRadiusSmall*self%radiusCore) then
           ! At small radii use a Taylor series for numerical accuracy.
           mass          =+  4.0d0/3.0d0*Pi                   *self%densitySolitonCentral*radius_**3                    &
                &         - 32.0d0/5.0d0*Pi*coefficientCore   *self%densitySolitonCentral*radius_**5/self%radiusCore**2 &
                &         +144.0d0/7.0d0*Pi*coefficientCore**2*self%densitySolitonCentral*radius_**7/self%radiusCore**4
        else
           ! Use the exact solution.
           termCore      =+sqrt(coefficientCore)                              &
                &         *                 radius_                           &
                &         *                            self%radiusCore        &
                &         /(coefficientCore*radius_**2+self%radiusCore**2)**7
           termPolynomial=+  3465.0d0*coefficientCore**6*radius_**12                     &
                &         + 23100.0d0*coefficientCore**5*radius_**10*self%radiusCore** 2 &
                &         + 65373.0d0*coefficientCore**4*radius_** 8*self%radiusCore** 4 &
                &         +101376.0d0*coefficientCore**3*radius_** 6*self%radiusCore** 6 &
                &         + 92323.0d0*coefficientCore**2*radius_** 4*self%radiusCore** 8 &
                &         + 48580.0d0*coefficientCore*   radius_** 2*self%radiusCore**10 &
                &         -  3465.0d0                               *self%radiusCore**12
           termArcTangent=+3465.0d0                                            &
                &         *atan(sqrt(coefficientCore)*radius_/self%radiusCore)
           termPrefactor =+Pi                            &
                &         *self%densitySolitonCentral    &
                &         *self%radiusCore           **3 &
                &         /(                             &
                &           +53760.0d0                   &
                &           *coefficientCore**1.5d0      &
                &         )
           mass          = + termPrefactor               &
                &         *(                             &
                &           +termArcTangent              &
                &           +termCore                    &
                &           *termPolynomial              &
                &          )
        end if
        return
      end function massSoliton
      
    end function solitonNFWMassEnclosedBySphere

    double precision function solitonNFWDensity(self,coordinates) result(density)
      !!{
      Return the density at the specified {\normalfont \ttfamily coordinates} in a soliton and NFW mass distribution.
      !!}
      implicit none
      class           (massDistributionSolitonNFW), intent(inout) :: self
      class           (coordinate                ), intent(in   ) :: coordinates
      double precision                                            :: radiusScaleFree, radiusCoreFree

      radiusScaleFree=+coordinates%rSpherical()/self%radiusScale
      radiusCoreFree =+coordinates%rSpherical()/self%radiusCore
      
      if (coordinates%rSpherical() <= 0.0d0) then
         ! Zero radius.
         density=+self%densitySolitonCentral
      else if (coordinates%rSpherical() < self%radiusSoliton) then
         ! Soliton regime.
         density=+self%densitySolitonCentral                  /(+1.0d0+coefficientCore*radiusCoreFree**2)**8
      else
         ! NFW regime.
         density=+self%densityNormalizationNFW/radiusScaleFree/(+1.0d0+                radiusScaleFree  )**2
      end if
      return
   end function solitonNFWDensity

   double precision function solitonNFWDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
      !!{
      Return the radial density gradient at the specified {\normalfont \ttfamily coordinates} in a soliton and NFW mass distribution.
      !!}
      implicit none
      class           (massDistributionSolitonNFW) , intent(inout), target   :: self
      class           (coordinate              )   , intent(in   )           :: coordinates
      logical                                      , intent(in   ), optional :: logarithmic
      double precision                                                       :: radiusScaleFree, radiusCoreFree
      !![
      <optionalArgument name="logarithmic" defaultsTo=".false."/>
      !!]

      radiusScaleFree=+coordinates%rSpherical()/self%radiusScale
      radiusCoreFree =+coordinates%rSpherical()/self%radiusCore

      if (coordinates%rSpherical() <= 0.0d0) then
         ! Zero radius.
         densityGradient   =0.0d0
      else if (coordinates%rSpherical() < self%radiusSoliton) then
         ! Soliton regime.
         if (logarithmic) then
            densityGradient=-16.0d0                                       &
                 &          *       coefficientCore*radiusCoreFree**2     &
                 &          /(1.0d0+coefficientCore*radiusCoreFree**2)
         else
            densityGradient=-16.0d0                                       &
                 &          *self%densitySolitonCentral                   &
                 &          /self%radiusCore                              &
                 &          *       coefficientCore*radiusCoreFree        &
                 &          /(1.0d0+coefficientCore*radiusCoreFree**2)**9
         end if
      else
         ! NFW regime.
         if (logarithmic) then
            densityGradient=  -1.0d0                  &
                 &            -2.0d0*radiusScaleFree  &
                 &          /(+1.0d0+radiusScaleFree)
         else
            densityGradient=-self%densityNormalizationNFW      &
                 &          /self%radiusScale                  &
                 &          *(+1.0d0+3.0d0*radiusScaleFree)    &
                 &          /(+1.0d0+      radiusScaleFree)**3 &
                 &          /              radiusScaleFree **2
         end if
      end if
      return
   end function solitonNFWDensityGradientRadial
   
   double precision function solitonNFWRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
     !!{
     Computes the radius enclosing a given mean density for soliton NFW mass distributions.
     !!}
     implicit none
     class           (massDistributionSolitonNFW), intent(inout), target   :: self
     double precision                            , intent(in   )           :: density
     double precision                            , intent(in   ), optional :: radiusGuess
     double precision                            , parameter               :: epsilonDeltaDensityFractional=4.0d-2
     
     if      (density >= self%densitySolitonCentral                                      ) then
        ! Above the central density, return a radius of zero.
        radius=0.0d0
     else if (density >  self%densitySolitonCentral*(1.0d0-epsilonDeltaDensityFractional)) then
        ! For densities close to, but below the central density, use a series solution.
        radius=+0.5d0                              &
             & *self%radiusCore                    &
             & *sqrt(                              &
             &       +5.0d0                        &
             &       /6.0d0                        &
             &       /coefficientCore              &
             &       *(                            &
             &         +1.0d0                      &
             &         -     density               &
             &         /self%densitySolitonCentral &
             &        )                            &
             &      )
     else
        radius=sphericalTabulatedRadiusEnclosingDensity(self,density,radiusGuess)
     end if
     return
   end function solitonNFWRadiusEnclosingDensity
   
   double precision function solitonNFWRadiusEnclosingDensityNumerical(self,density,radiusGuess) result(radius)
     !!{
     Computes the radius enclosing a given mean density for soliton NFW mass distributions.
     !!}
     implicit none
     class           (massDistributionSolitonNFW), intent(inout), target   :: self
     double precision                            , intent(in   )           :: density
     double precision                            , intent(in   ), optional :: radiusGuess
     double precision                            , parameter               :: epsilonDeltaDensityFractional=1.0d-2

     if      (density >= self%densitySolitonCentral                                      ) then
        radius=0.0d0
     else if (density >  self%densitySolitonCentral*(1.0d0-epsilonDeltaDensityFractional)) then
        ! For densities close to, but below the central density, use a series solution.
        radius=+0.5d0                              &
             & *self%radiusCore                    &
             & *sqrt(                              &
             &       +5.0d0                        &
             &       /6.0d0                        &
             &       /coefficientCore              &
             &       *(                            &
             &         +1.0d0                      &
             &         -     density               &
             &         /self%densitySolitonCentral &
             &        )                            &
             &      )
     else
        radius=massDistributionRadiusEnclosingDensityNumerical(self,density,radiusGuess)
     end if
     return
   end function solitonNFWRadiusEnclosingDensityNumerical
   
   subroutine solitonNFWParameters(self,densityNormalization,radiusNormalization,parameters,container)
      !!{
      Establish parameters for tabulation.
      !!}
      implicit none
      class           (massDistributionSolitonNFW), intent(inout)                            :: self
      double precision                            , intent(  out)                            :: densityNormalization, radiusNormalization
      double precision                            , intent(inout), allocatable, dimension(:) :: parameters
      type             (massDistributionContainer), intent(  out), pointer                   :: container

      if (.not.containerSolitonNFWInitialized) then
         allocate(containerSolitonNFW)
         call containerSolitonNFW%initialize(3)
         containerSolitonNFW%mass                      %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%mass                      %parametersCountPer   =+20_c_size_t
         containerSolitonNFW%radiusEnclosingDensity    %radiusCountPer       =+20_c_size_t
         containerSolitonNFW%radiusEnclosingDensity    %parametersCountPer   =+20_c_size_t
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
         containerSolitonNFW%nameParameters                               (1)='radiusCoreOverRadiusScale'
         containerSolitonNFW%descriptionParameters                        (1)='The ratio of core to scale radii(r_c/rₛ).'
         containerSolitonNFW%nameParameters                               (2)='radiusSolitonOverRadiusScale'
         containerSolitonNFW%descriptionParameters                        (2)='The ratio of soliton to scale radii(rₛₒₗ/rₛ).'
         containerSolitonNFW%nameParameters                               (3)='densitySolitonCentralOverDensityScale'
         containerSolitonNFW%descriptionParameters                        (3)='The ratio of soliton central to scale densities(ρ₀/ρₛ).'
         containerSolitonNFWInitialized                                      =.true.
      end if
      allocate(parameters(3))
      densityNormalization =  self%densityNormalizationNFW
      radiusNormalization  =  self%radiusScale
      parameters(1)        =  self%radiusCoreScaleFree
      parameters(2)        =  self%radiusSolitonScaleFree
      parameters(3)        =  self%densitySolitonScaleFree
      container            => containerSolitonNFW
      return
   end subroutine solitonNFWParameters

   subroutine solitonNFWDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
      !!{
      Return an input parameter list descriptor which could be used to recreate this object.
      !!}
      use :: Input_Parameters, only : inputParameters
      implicit none
      class    (massDistributionSolitonNFW), intent(inout)           :: self
      type     (inputParameters           ), intent(inout)           :: descriptor
      logical                              , intent(in   ), optional :: includeClass  , includeFileModificationTimes
      character(len=18)                                              :: parameterLabel
      type     (inputParameters           )                          :: parameters
      !$GLC attributes unused :: includeFileModificationTimes

      if (.not.present(includeClass) .or. includeClass) call descriptor%addParameter('massDistribution','solitonNFW')
      parameters = descriptor%subparameters('massDistribution')
      write(parameterLabel,'(e17.10)') self%densityNormalizationNFW
      call parameters%addParameter('densityNormalizationNFW',trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%radiusScale
      call parameters%addParameter('radiusScale'            ,trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%radiusCore
      call parameters%addParameter('radiusCore'             ,trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%radiusSoliton
      call parameters%addParameter('radiusSoliton'          ,trim(adjustl(parameterLabel)))
      write(parameterLabel,'(e17.10)') self%densitySolitonCentral
      call parameters%addParameter('densitySolitonCentral'  ,trim(adjustl(parameterLabel)))
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

      suffix=String_C_To_Fortran(massDistributionSolitonNFWSourceDigest)
      return
   end function solitonNFWSuffix
