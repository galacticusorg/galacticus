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
  Implements a finite resolution NFW spherical mass distribution.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <massDistribution name="massDistributionSphericalFiniteResolutionNFW">
   <description>
    A mass distribution class which applies a finite resolution to an NFW density profile, typically to mimic the effects
    of finite resolution in an N-body simulation. Specifically, the density profile is given by
    \begin{equation}
    \rho(r) = \rho_\mathrm{NFW}(r) \left( 1 + \left[ \frac{\Delta x}{r} \right]^2 \right)^{-1/2},
    \end{equation}
    where $\Delta x$ is the larger of the resolution length, {\normalfont \ttfamily [lengthResolution]}, and the radius in the
    original profile enclosing the mass resolution, {\normalfont \ttfamily [massResolution]}.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSpherical) :: massDistributionSphericalFiniteResolutionNFW
     !!{
     Implementation of a finite resolution spherical mass distribution.
     !!}
     private
     double precision                                            :: lengthResolution                                       , radiusScale                                   , &
          &                                                         radiusVirial                                           , mass                                          , &
          &                                                         densityNormalization                                   , lengthResolutionScaleFree
     double precision                                            :: potentialRadiusPrevious                                , potentialPrevious                             , &
          &                                                         massEnclosedMassPrevious                               , massEnclosedRadiusPrevious                    , &
          &                                                         densityRadiusPrevious                                  , densityPrevious                               , &
          &                                                         densityNormalizationPrevious                           , radiusEnclosingDensityDensityPrevious         , &
          &                                                         radiusEnclosingDensityPrevious                         , radiusEnclosingMassMassPrevious               , &
          &                                                         radiusEnclosingMassPrevious                            , energyPrevious                                , &
          &                                                         lengthResolutionScaleFreePotentialPrevious             , lengthResolutionPotentialTerm1                , &
          &                                                         lengthResolutionPotentialTerm2                         , lengthResolutionPotentialSquare               , &
          &                                                         lengthResolutionPotentialRootSquare                    , lengthResolutionPotentialOnePlusSquare        , &
          &                                                         lengthResolutionPotentialSqrtOnePlusSquare             , lengthResolutionPotentialOnePlusTwoSquare     , &
          &                                                         lengthResolutionPotentialOnePlusSquareP1p5             , lengthResolutionPotentialAtanh                , &
          &                                                         concentrationPotentialTerm
     ! Radius-enclosing-density tabulation.
     logical                                                     :: radiusEnclosingDensityTableInitialized
     integer                                                     :: radiusEnclosingDensityTableLengthResolutionCount       , radiusEnclosingDensityTableDensityCount
     double precision              , allocatable, dimension(:  ) :: radiusEnclosingDensityTableLengthResolution            , radiusEnclosingDensityTableDensity
     double precision              , allocatable, dimension(:,:) :: radiusEnclosingDensityTable
     type            (interpolator), allocatable                 :: radiusEnclosingDensityTableLengthResolutionInterpolator, radiusEnclosingDensityTableDensityInterpolator
     double precision                                            :: radiusEnclosingDensityDensityMinimum                   , radiusEnclosingDensityDensityMaximum          , &
          &                                                         radiusEnclosingDensityLengthResolutionMinimum          , radiusEnclosingDensityLengthResolutionMaximum
     ! Radius-enclosing-mass tabulation.
     logical                                                     :: radiusEnclosingMassTableInitialized
     integer                                                     :: radiusEnclosingMassTableLengthResolutionCount          , radiusEnclosingMassTableMassCount
     double precision              , allocatable, dimension(:  ) :: radiusEnclosingMassTableLengthResolution               , radiusEnclosingMassTableMass
     double precision              , allocatable, dimension(:,:) :: radiusEnclosingMassTable
     type            (interpolator), allocatable                 :: radiusEnclosingMassTableLengthResolutionInterpolator   , radiusEnclosingMassTableMassInterpolator
     double precision                                            :: radiusEnclosingMassMassMinimum                         , radiusEnclosingMassMassMaximum                , &
          &                                                         radiusEnclosingMassLengthResolutionMinimum             , radiusEnclosingMassLengthResolutionMaximum
     ! Energy tabulation.
     logical                                                     :: energyTableInitialized
     integer                                                     :: energyTableLengthResolutionCount                       , energyTableRadiusOuterCount
     double precision              , allocatable, dimension(:  ) :: energyTableLengthResolution                            , energyTableRadiusOuter
     double precision              , allocatable, dimension(:,:) :: energyTable
     type            (interpolator), allocatable                 :: energyTableLengthResolutionInterpolator                , energyTableRadiusOuterInterpolator
     double precision                                            :: energyRadiusOuterMinimum                               , energyRadiusOuterMaximum                      , &
          &                                                         energyLengthResolutionMinimum                          , energyLengthResolutionMaximum
     ! Enclosed mass quantities.
     double precision                                            :: lengthResolutionScaleFreeLowerTerm                     , lengthResolutionScaleFreeSquared              , &
          &                                                         lengthResolutionScaleFreeCubed                         , lengthResolutionScaleFreeOnePlusTerm          , &
          &                                                         lengthResolutionScaleFreeOnePlus2Term                  , lengthResolutionScaleFreeSqrtTerm             , &
          &                                                         lengthResolutionScaleFreeSqrt2Term                     , lengthResolutionScaleFreeSqrtCubedTerm        , &
          &                                                         lengthResolutionScaleFreePrevious
   contains
     !![
     <methods>
       <method method="radiusEnclosingDensityTabulate" description="Tabulate the radius enclosing a given density as a function of density and core radius."           />
       <method method="radiusEnclosingMassTabulate"    description="Tabulate the radius enclosing a given mass as a function of density and core radius."              />
       <method method="energyTabulate"                 description="Tabulate the energy as a function of concentration and core radius."                               />
       <method method="densityScaleFree"               description="The density of the profile in units where the mass and scale length are both 1."                   />
       <method method="massEnclosedScaleFree"          description="The mass enclosed of the profile in units where the mass and scale length are both 1."             />
       <method method="storeDensityTable"              description="Store the tabulated radius-enclosing-density to file."                                             />
       <method method="restoreDensityTable"            description="Attempt to restore the tabulated radius-enclosing-density from file, returning true if successful."/>
       <method method="storeMassTable"                 description="Store the tabulated radius-enclosing-mass to file."                                                />
       <method method="restoreMassTable"               description="Attempt to restore the tabulated radius-enclosing-mass from file, returning true if successful."   />
       <method method="storeEnergyTable"               description="Store the tabulated energy to file."                                                               />
       <method method="restoreEnergyTable"             description="Attempt to restore the tabulated energy from file, returning true if successful."                  />
     </methods>
     !!]
     procedure :: density                        => sphericalFiniteResolutionNFWDensity
     procedure :: densityGradientRadial          => sphericalFiniteResolutionNFWDensityGradientRadial
     procedure :: massEnclosedBySphere           => sphericalFiniteResolutionNFWMassEnclosedBySphere
     procedure :: potentialIsAnalytic            => sphericalFiniteResolutionNFWPotentialIsAnalytic
     procedure :: potential                      => sphericalFiniteResolutionNFWPotential
     procedure :: radiusEnclosingMass            => sphericalFiniteResolutionNFWRadiusEnclosingMass
     procedure :: radiusEnclosingDensity         => sphericalFiniteResolutionNFWRadiusEnclosingDensity
     procedure :: energy                         => sphericalFiniteResolutionNFWEnergy
     procedure :: radiusEnclosingDensityTabulate => sphericalFiniteResolutionNFWRadiusEnclosingDensityTabulate
     procedure :: radiusEnclosingMassTabulate    => sphericalFiniteResolutionNFWRadiusEnclosingMassTabulate
     procedure :: energyTabulate                 => sphericalFiniteResolutionNFWEnergyTabulate
     procedure :: densityScaleFree               => sphericalFiniteResolutionNFWDensityScaleFree
     procedure :: massEnclosedScaleFree          => sphericalFiniteResolutionNFWMassEnclosedScaleFree
     procedure :: storeDensityTable              => sphericalFiniteResolutionNFWStoreDensityTable
     procedure :: restoreDensityTable            => sphericalFiniteResolutionNFWRestoreDensityTable
     procedure :: storeMassTable                 => sphericalFiniteResolutionNFWStoreMassTable
     procedure :: restoreMassTable               => sphericalFiniteResolutionNFWRestoreMassTable
     procedure :: storeEnergyTable               => sphericalFiniteResolutionNFWStoreEnergyTable
     procedure :: restoreEnergyTable             => sphericalFiniteResolutionNFWRestoreEnergyTable
  end type massDistributionSphericalFiniteResolutionNFW

  interface massDistributionSphericalFiniteResolutionNFW
     !!{
     Constructors for the {\normalfont \ttfamily sphericalFiniteResolutionNFW} mass distribution class.
     !!}
     module procedure sphericalFiniteResolutionNFWConstructorParameters
     module procedure sphericalFiniteResolutionNFWConstructorInternal
  end interface massDistributionSphericalFiniteResolutionNFW

  ! Tabulation resolution parameters.
  integer                                                       , parameter :: radiusEnclosingDensityTableDensityPointsPerDecade         =100
  integer                                                       , parameter :: radiusEnclosingDensityTableLengthResolutionPointsPerDecade=100
  integer                                                       , parameter :: radiusEnclosingMassTableMassPointsPerDecade               =100
  integer                                                       , parameter :: radiusEnclosingMassTableLengthResolutionPointsPerDecade   =100
  integer                                                       , parameter :: energyTableRadiusOuterPointsPerDecade                     =100
  integer                                                       , parameter :: energyTableLengthResolutionPointsPerDecade                =100

  ! Sub-module-scope variables used in integrations.
  class           (massDistributionSphericalFiniteResolutionNFW), pointer   :: self_
  integer                                                                   :: iLengthResolution_                                               , iDensity_, &
       &                                                                       iMass_
  !$omp threadprivate(self_,iLengthResolution_,iDensity_,iMass_)

  ! Largest radius for precise arctanh() evaluation.
  double precision                                              , parameter :: radiusScaleFreeLargeATanh                                 =1.0d+6

contains

  function sphericalFiniteResolutionNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sphericalFiniteResolutionNFW} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalFiniteResolutionNFW)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: lengthResolution, radiusScale, &
         &                                                                           radiusVirial    , mass
    type            (varying_string                              )                :: componentType   , massType

    !![
    <inputParameter>
      <name>lengthResolution</name>
      <source>parameters</source>
      <description>The resolution length scale.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusScale</name>
      <source>parameters</source>
      <description>The NFW scale radius.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusVirial</name>
      <source>parameters</source>
      <description>The virial radius.</description>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <source>parameters</source>
      <description>The mass within the virial radius.</description>
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
    self=massDistributionSphericalFiniteResolutionNFW(lengthResolution,radiusScale,radiusVirial,mass,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sphericalFiniteResolutionNFWConstructorParameters
  
  function sphericalFiniteResolutionNFWConstructorInternal(lengthResolution,radiusScale,radiusVirial,mass,componentType,massType) result(self)
    !!{
    Constructor for {\normalfont \ttfamily sphericalFiniteResolutionNFW} mass distribution class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionSphericalFiniteResolutionNFW)                          :: self
    double precision                                              , intent(in   )           :: lengthResolution , radiusScale, &
         &                                                                                     radiusVirial     , mass
    type            (enumerationComponentTypeType                ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType                     ), intent(in   ), optional :: massType
    double precision                                                                        :: radiusScaleFree
    !![
    <constructorAssign variables="lengthResolution, radiusScale, radiusVirial, mass, componentType, massType"/>
    !!]

    self%dimensionless                                 =.false.
    self%lengthResolutionScalefreePrevious             =-huge(0.0d0)
    self%massEnclosedMassPrevious                      =-huge(0.0d0)
    self%massEnclosedRadiusPrevious                    =-huge(0.0d0)
    self%potentialPrevious                             =-huge(0.0d0)
    self%potentialRadiusPrevious                       =-huge(0.0d0)
    self%lengthResolutionScaleFreePotentialPrevious    =-huge(0.0d0)
    self%concentrationPotentialTerm                    =-huge(0.0d0)
    self%densityRadiusPrevious                         =-huge(0.0d0)
    self%densityPrevious                               =-huge(0.0d0)
    self%densityNormalizationPrevious                  =-huge(0.0d0)
    self%radiusEnclosingDensityDensityPrevious         =-huge(0.0d0)
    self%radiusEnclosingDensityPrevious                =-huge(0.0d0)
    self%radiusEnclosingMassMassPrevious               =-huge(0.0d0)
    self%radiusEnclosingMassPrevious                   =-huge(0.0d0)
    self%energyPrevious                                =+huge(0.0d0)
    ! Radius enclosing density table initialization.
    self%radiusEnclosingDensityDensityMinimum          =+huge(0.0d0)
    self%radiusEnclosingDensityDensityMaximum          =-huge(0.0d0)
    self%radiusEnclosingDensityLengthResolutionMinimum =+huge(0.0d0)
    self%radiusEnclosingDensityLengthResolutionMaximum =-huge(0.0d0)
    self%radiusEnclosingDensityTableInitialized        =.false.
    ! Radius enclosing mass table initialization.
    self%radiusEnclosingMassMassMinimum                =+huge(0.0d0)
    self%radiusEnclosingMassMassMaximum                =-huge(0.0d0)
    self%radiusEnclosingMassLengthResolutionMinimum    =+huge(0.0d0)
    self%radiusEnclosingMassLengthResolutionMaximum    =-huge(0.0d0)
    self%radiusEnclosingMassTableInitialized           =.false.
    ! Energy table initialization.
    self%energyRadiusOuterMinimum                      =+huge(0.0d0)
    self%energyRadiusOuterMaximum                      =-huge(0.0d0)
    self%energyLengthResolutionMinimum                 =+huge(0.0d0)
    self%energyLengthResolutionMaximum                 =-huge(0.0d0)
    self%energyTableInitialized                        =.false.
    ! Construct profile quantities.
    radiusScaleFree               =+    radiusVirial/radiusScale
    self%lengthResolutionScaleFree=+lengthResolution/radiusScale
    self%densityNormalization     =+mass/4.0d0/Pi/radiusScale**3/(log(1.0d0+radiusScaleFree)-radiusScaleFree/(1.0d0+radiusScaleFree))
    return
  end function sphericalFiniteResolutionNFWConstructorInternal

  double precision function sphericalFiniteResolutionNFWDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    class           (coordinate                                  ), intent(in   ) :: coordinates
    double precision                                                              :: radiusScaleFree
    
    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %radiusScale
    density        =+self       %densityNormalization                             &
         &          /sqrt(+self%lengthResolutionScaleFree**2+radiusScaleFree **2) &
         &          /    (+1.0d0                            +radiusScaleFree)**2
    return
  end function sphericalFiniteResolutionNFWDensity

  double precision function sphericalFiniteResolutionNFWDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a finiteResolution spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target   :: self
    class           (coordinate                                  ), intent(in   )           :: coordinates
    logical                                                       , intent(in   ), optional :: logarithmic
    double precision                                                                        :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    radiusScaleFree=+coordinates%rSpherical () &
         &          /self       %radiusScale
    densityGradient=-3.0d0                                                             &
         &          +2.0d0/(1.0d0+ radiusScaleFree                                   ) &
         &          +1.0d0/(1.0d0+(radiusScaleFree/self%lengthResolutionScaleFree)**2)
    if (.not.logarithmic_) &
         densityGradient=+            densityGradient              &
         &               *self       %density        (coordinates) &
         &               /coordinates%rSpherical     (           )
   return
  end function sphericalFiniteResolutionNFWDensityGradientRadial

  double precision function sphericalFiniteResolutionNFWMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Returns the enclosed mass (in $M_\odot$) at the given {\normalfont \ttfamily radius} (given in units of Mpc). The analytic
    solution (computed using Mathematica) is
    \begin{equation}
    M(x) = 4 \pi \rho_0 r_\mathrm{s}^3 \left[ -\frac{\sqrt{x^2+X^2}}{(1+x) \left(1+X^2\right)}+\tanh ^{-1}\left(\frac{x}{\sqrt{x^2+X^2}}\right)+\frac{\left(1+2X^2\right) \tanh ^{-1}\left(\frac{X^2-x}{\sqrt{1+X^2} \sqrt{x^2+X^2}}\right)}{\left(1+X^2\right)^{3/2}} -\frac{\left(1 + 2 X^2\right) \tanh ^{-1}\left(\sqrt{\frac{X^2}{1 + X^2}}\right)}{\left(1+ X^2\right)^{3/2}}+\frac{\sqrt{X^2}}{1 + X^2} \right],
    \end{equation}
    where $x=r/r_\mathrm{s}$, $X = \Delta x/r_\mathrm{s}$, and $r_\mathrm{s}$ is the NFW scale length.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target :: self
    double precision                                              , intent(in   )         :: radius
    double precision                                                                      :: radiusScaleFree
    
    radiusScaleFree=+     radius                                                                   &
         &          /self%radiusScale
    mass           =+self%densityNormalization                                                     &
         &          *self%radiusScale                                                          **3 &
         &          *self%massEnclosedScaleFree(radiusScaleFree,self%lengthResolutionScaleFree)
    return
  end function sphericalFiniteResolutionNFWMassEnclosedBySphere

  double precision function sphericalFiniteResolutionNFWMassEnclosedScaleFree(self,radiusScaleFree,lengthResolutionScaleFree) result(mass)
    !!{
    Returns the scale-free enclosed mass at the given scale-free radius. The analytic solution (computed using Mathematica) is
    \begin{equation}
    M(x) = 4 \pi \left[ -\frac{\sqrt{x^2+X^2}}{(1+x) \left(1+X^2\right)}+\tanh ^{-1}\left(\frac{x}{\sqrt{x^2+X^2}}\right)+\frac{\left(1+2X^2\right) \tanh ^{-1}\left(\frac{X^2-x}{\sqrt{1+X^2} \sqrt{x^2+X^2}}\right)}{\left(1+X^2\right)^{3/2}} -\frac{\left(1 + 2 X^2\right) \tanh ^{-1}\left(\sqrt{\frac{X^2}{1 + X^2}}\right)}{\left(1+ X^2\right)^{3/2}}+\frac{\sqrt{X^2}}{1 + X^2} \right],
    \end{equation}
    where $x=r/r_\mathrm{s}$, $X = \Delta x/r_\mathrm{s}$, and $r_\mathrm{s}$ is the NFW scale length.
    !!}
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    double precision                                              , intent(in   ) :: radiusScaleFree                       , lengthResolutionScaleFree
    double precision                                              , parameter     :: radiusScaleFreeSmall           =1.0d-3, radiusScaleFreeLarge      =1.0d4
    double precision                                                              :: radiusScaleFreeEffective              , arctanhTerm1                    , &
         &                                                                           arctanhTerm                           , radiusNormTerm                  , &
         &                                                                           radiusScaleFreeEffectiveSquared

    
    if (radiusScaleFree /= self%massEnclosedRadiusPrevious) then
       self%massEnclosedRadiusPrevious=+radiusScaleFree
       if (lengthResolutionScaleFree /= self%lengthResolutionScaleFreePrevious) then
          ! Construct quantities used in the mass enclosed within a sphere.
          self%lengthResolutionScaleFreePrevious     =     lengthResolutionScaleFree
          self%lengthResolutionScaleFreeSquared      =self%lengthResolutionScaleFreePrevious**2
          self%lengthResolutionScaleFreeCubed        =self%lengthResolutionScaleFreePrevious**3
          self%lengthResolutionScaleFreeOnePlusTerm  =+1.0d0+      self%lengthResolutionScaleFreeSquared
          self%lengthResolutionScaleFreeOnePlus2Term =+1.0d0+2.0d0*self%lengthResolutionScaleFreeSquared
          self%lengthResolutionScaleFreeSqrtTerm     =sqrt(self%lengthResolutionScaleFreeOnePlusTerm )
          self%lengthResolutionScaleFreeSqrt2Term    =sqrt(self%lengthResolutionScaleFreeOnePlus2Term)
          self%lengthResolutionScaleFreeSqrtCubedTerm=self%lengthResolutionScaleFreeSqrtTerm**3
          ! For large values of the argument to arctanh(), use a series solution to avoiding floating point errors.
          if (self%lengthResolutionScaleFreePrevious > radiusScaleFreeLargeATanh) then
             arctanhTerm=-log(                                        &
                  &           +2.0d0                                  &
                  &           *self%lengthResolutionScaleFreePrevious &
                  &          )                                        &
                  &      /2.0d0                                       &
                  &      +1.0d0                                       &
                  &      /2.0d0                                       &
                  &      /self%lengthResolutionScaleFreePrevious      &
                  &      +1.0d0                                       &
                  &      /8.0d0                                       &
                  &      /self%lengthResolutionScaleFreePrevious**2
          else
             arctanhTerm=+atanh(                                                 &
                  &             +(+1.0d0-self%lengthResolutionScaleFreePrevious) &
                  &             /self%lengthResolutionScaleFreeSqrtTerm          &
                  &            ) 
          end if
          self%lengthResolutionScaleFreeLowerTerm=+self%lengthResolutionScaleFreePrevious      &
               &                                  /self%lengthResolutionScaleFreeOnePlusTerm   &
               &                                  +2.0d0                                       &
               &                                  *self%lengthResolutionScaleFreeOnePlus2Term  &
               &                                  *arctanhTerm                                 &
               &                                  /self%lengthResolutionScaleFreeSqrtCubedTerm
       end if       
       if (radiusScaleFree < radiusScaleFreeSmall) then
          ! Series expansion for small radii.
          self%massEnclosedMassPrevious=+  radiusScaleFree**3                                                                                                             &
            &                           *(                                                                                                                                &
            &                                                                                                       +1.0d0 /self%lengthResolutionScaleFreePrevious/ 3.0d0 &
            &                             +radiusScaleFree   *(                                                     +1.0d0 /self%lengthResolutionScaleFreePrevious/ 2.0d0 &
            &                             +radiusScaleFree   * ( 1.0d0+(+6.0d0*self%lengthResolutionScaleFreeSquared-1.0d0)/self%lengthResolutionScaleFreeCubed   /10.0d0 &
            &                             +radiusScaleFree   *  (1.0d0-(+4.0d0*self%lengthResolutionScaleFreeSquared-1.0d0)/self%lengthResolutionScaleFreeCubed   / 6.0d0 &
            &                                                   )                                                                                                         &
            &                                                  )                                                                                                          &
            &                                                 )                                                                                                           &
            &                            )
       else
          ! Full analytic solution.
          !! Limit the evaluation to some large radius.
          radiusScaleFreeEffective=min(radiusScaleFree,radiusScaleFreeLarge)
          radiusScaleFreeEffective       =min(radiusScaleFree,radiusScaleFreeLarge)
          radiusScaleFreeEffectiveSquared=radiusScaleFreeEffective**2
          radiusNormTerm                 =sqrt(+radiusScaleFreeEffectiveSquared+self%lengthResolutionScaleFreeSquared)
          if (radiusScaleFreeEffective > radiusScaleFreeLargeATanh*self%lengthResolutionScaleFreePrevious) then
             arctanhTerm1=+log  (                                       &
                  &              +4.0d0                                 &
                  &              *     radiusScaleFreeEffectiveSquared  &
                  &              /self%lengthResolutionScaleFreeSquared &
                  &             )                                       &
                  &       /2.0d0                                        &
                  &       -self%lengthResolutionScaleFreeSquared        &
                  &       /8.0d0                                        &
                  &       /     radiusScaleFreeEffectiveSquared
          else
             arctanhTerm1=+atanh(                          &
                  &              +radiusScaleFreeEffective &
                  &              /radiusNormTerm           &
                  &             )
          end if
          self%massEnclosedMassPrevious=-        radiusNormTerm                        &
               &                        /(+1.0d0+radiusScaleFreeEffective)             &
               &                        /self%lengthResolutionScaleFreeOnePlusTerm     &
               &                        -2.0d0                                         &
               &                        *self%lengthResolutionScaleFreeOnePlus2Term    &
               &                        *atanh(                                        &
               &                               +(                                      &
               &                                 +1.0d0                                &
               &                                 +radiusScaleFreeEffective             &
               &                                 -radiusNormTerm                       &
               &                                )                                      &
               &                               /self%lengthResolutionScaleFreeSqrtTerm &
               &                              )                                        &
               &                        /self%lengthResolutionScaleFreeSqrtCubedTerm   &
               &                        +arctanhTerm1                                  &
               &                        +self%lengthResolutionScaleFreeLowerTerm
          !! Beyond the limiting radius assume logarithmic growth in mass as appropriate for an r⁻³ profile.
          if (radiusScaleFree > radiusScaleFreeEffective)                     &
               & self%massEnclosedMassPrevious=+self%massEnclosedMassPrevious &
               &                               *log(                          &
               &                                    +radiusScaleFree          &
               &                                    /radiusScaleFreeEffective &
               &                                   )
       end if
       self%massEnclosedMassPrevious=+4.0d0                         &
            &                        *Pi                            &
            &                        *self%massEnclosedMassPrevious
    end if
    mass=self%massEnclosedMassPrevious
    return
  end function sphericalFiniteResolutionNFWMassEnclosedScaleFree
  
  logical function sphericalFiniteResolutionNFWPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self

    isAnalytic=.true.
    return
  end function sphericalFiniteResolutionNFWPotentialIsAnalytic

  double precision function sphericalFiniteResolutionNFWPotential(self,coordinates,status) result(potential)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc). The analytic solution (computed using Mathematica) is
    \begin{eqnarray}
    \Phi(x) &=& -\frac{\mathrm{G} M}{r_\mathrm{s}}  \nonumber \\
            & & \left\{ +\frac{\sqrt{x^2+X^2}}{x \left(X^2+1\right)} \right. \nonumber \\
            & & -\frac{X^2 \log \left(\sqrt{X^2+1} \sqrt{x^2+X^2}-x+X^2\right)}{\left(X^2+1\right)^{3/2}} \nonumber \\
            & & -\frac{\tanh ^{-1}\left(\frac{x}{\sqrt{x^2+X^2}}\right)}{x} \nonumber \\
            & & -\frac{\left(2 X^2+1\right) \tanh ^{-1}\left(\frac{X^2-x}{\sqrt{X^2+1} \sqrt{x^2+X^2}}\right)}{x \left(X^2+1\right)^{3/2}} \nonumber \\
            & & -\frac{\sqrt{X^2}}{x \left(X^2+1\right)}+\frac{X^2 \log (x+1)}{\left(X^2+1\right)^{3/2}} \nonumber \\
            & & +\frac{\left(2 X^2+1\right) \tanh ^{-1}\left(\sqrt{\frac{X^2}{X^2+1}}\right)}{x \left(X^2+1\right)^{3/2}} \nonumber \\
            & & \left. +\frac{ \left(\sqrt{X^2+1}-X^2 \log \left(\sqrt{X^2+1}-1\right)\right)}{\left(X^2+1\right)^{3/2}} \right\} \nonumber \\
            & & /\left[\log (1+c)-\frac{c}{1+c}\right]
    \end{eqnarray}
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess     , structureErrorCodeInfinite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Error                           , only : Error_Report
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target   :: self
    class           (coordinate                                  ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType           ), intent(  out), optional :: status
    double precision                                              , parameter               :: radiusScaleFreeSmall     =1.0d-3
    double precision                                              , parameter               :: radiusScaleFreeLarge     =1.0d+5
    double precision                                                                        :: lengthResolutionScaleFree       , radiusScaleFree, &
         &                                                                                     concentration
    
    if (present(status)) status=structureErrorCodeSuccess
    if (coordinates%rSpherical() /= self%potentialRadiusPrevious) then
       self%potentialRadiusPrevious=+coordinates%rSpherical ()
       radiusScaleFree             =+coordinates%rSpherical () &
            &                       /self       %radiusScale
       lengthResolutionScaleFree   =+self       %lengthResolution &
            &                       /self       %radiusScale
       if (lengthResolutionScaleFree /= self%lengthResolutionScaleFreePotentialPrevious) then
          ! Recompute terms that depend only on the scale free resolution length.
          self%lengthResolutionScaleFreePotentialPrevious=lengthResolutionScaleFree
          self%lengthResolutionPotentialSquare           =lengthResolutionScaleFree**2
          self%lengthResolutionPotentialRootSquare       =sqrt(self%lengthResolutionPotentialSquare)
          self%lengthResolutionPotentialOnePlusSquare    =1.0d0+self%lengthResolutionPotentialSquare
          self%lengthResolutionPotentialSqrtOnePlusSquare=sqrt(self%lengthResolutionPotentialOnePlusSquare)
          self%lengthResolutionPotentialOnePlusTwoSquare =1.0d0+2.0d0*self%lengthResolutionPotentialSquare
          self%lengthResolutionPotentialOnePlusSquareP1p5=self%lengthResolutionPotentialOnePlusSquare**1.5d0
          self%lengthResolutionPotentialAtanh            =atanh(                                                                                         &
               &                                                +sqrt(+self%lengthResolutionPotentialSquare/self%lengthResolutionPotentialOnePlusSquare) &
               &                                               )
          self%lengthResolutionPotentialTerm1            =+(+1.0d0-lengthResolutionScaleFree   )                            &
               &                                          /self%lengthResolutionPotentialOnePlusSquare                      &
               &                                          +self%lengthResolutionPotentialSquare                             &
               &                                          *(                                                                &
               &                                            +asinh(lengthResolutionScaleFree   )                            &
               &                                            +log  (                                                         &
               &                                                   +(1.0d0+self%lengthResolutionPotentialSqrtOnePlusSquare) &
               &                                                   /            lengthResolutionScaleFree                   &
               &                                                  )                                                         &
               &                                           )                                                                &
               &                                          /self%lengthResolutionPotentialOnePlusSquareP1p5
          self%lengthResolutionPotentialTerm2            =+(                                                      &
               &                                            +self%lengthResolutionPotentialSqrtOnePlusSquare      &
               &                                            -self%lengthResolutionPotentialSquare                 &
               &                                            *log(                                                 &
               &                                                 -1.0d0                                           &
               &                                                 +self%lengthResolutionPotentialSqrtOnePlusSquare &
               &                                                )                                                 &
               &                                           )                                                      &
               &                                          /self%lengthResolutionPotentialOnePlusSquareP1p5
          concentration                                  =+self%radiusVirial &
               &                                          /self%radiusScale
          self%concentrationPotentialTerm                =1.0d0/(                          &
               &                                                 -          concentration  &
               &                                                 /   (1.0d0+concentration) &
               &                                                 +log(1.0d0+concentration) &
               &                                                )
       end if
       ! Evaluate the potential.
       if      (radiusScaleFree > radiusScaleFreeLarge) then
          ! Truncate to zero at very large radii.
          self%potentialPrevious=0.0d0
       else if (radiusScaleFree < radiusScaleFreeSmall) then
          ! Series expansion for small radii.
          self%potentialPrevious       =  -gravitationalConstant_internal         &
               &                          *self%mass                              &
               &                          /self%radiusScale                       &
               &                          *(                                      &
               &                            +self%lengthResolutionPotentialTerm1  &
               &                            -        radiusScaleFree**2           &
               &                            *(+1.0d0-radiusScaleFree            ) &
               &                            /        lengthResolutionScaleFree    &
               &                            /6.0d0                                &
               &                           )                                      &
               &                          *self%concentrationPotentialTerm
       else
          self%potentialPrevious       =  -gravitationalConstant_internal                                                                                                &
               &                          *self%mass                                                                                                                     &
               &                          /self%radiusScale                                                                                                              &
               &                          *(                                                                                                                             &
               &                            +     self%lengthResolutionPotentialRootSquare                                                                               &
               &                            /                   radiusScaleFree                                        /self%lengthResolutionPotentialOnePlusSquare      &
               &                            -       sqrt(      +radiusScaleFree**2+self%lengthResolutionPotentialSquare                                           )      &
               &                            /                   radiusScaleFree                                        /self%lengthResolutionPotentialOnePlusSquare      &
               &                            -                                                                           self%lengthResolutionPotentialOnePlusTwoSquare   &
               &                            *self%lengthResolutionPotentialAtanh                                                                                         &
               &                            /                   radiusScaleFree                                        /self%lengthResolutionPotentialOnePlusSquareP1p5  &
               &                            +atanh(                                                                                                                      &
               &                                   +            radiusScaleFree                                                                                          &
               &                                   /sqrt(      +radiusScaleFree**2+self%lengthResolutionPotentialSquare                                                ) &
               &                                  )                                                                                                                      &
               &                            /                   radiusScaleFree                                                                                          &
               &                            +                                                                           self%lengthResolutionPotentialOnePlusTwoSquare   &
               &                            *atanh(                                                                                                                      &
               &                                        (      -radiusScaleFree   +self%lengthResolutionPotentialSquare                                                ) &
               &                                   /                                                                    self%lengthResolutionPotentialSqrtOnePlusSquare  &
               &                                   /sqrt(      +radiusScaleFree**2+self%lengthResolutionPotentialSquare                                                ) &
               &                                  )                                                                                                                      &
               &                            /                   radiusScaleFree                                        /self%lengthResolutionPotentialOnePlusSquareP1p5  &
               &                            -                                      self%lengthResolutionPotentialSquare                                                  &
               &                            *log  (                                                                                                                      &
               &                                         +1.0d0+radiusScaleFree                                                                                          &
               &                                  )                                                                                                                      &
               &                            /                                                                           self%lengthResolutionPotentialOnePlusSquareP1p5  &
               &                            +                                      self%lengthResolutionPotentialSquare/self%lengthResolutionPotentialOnePlusSquareP1p5  &
               &                            *log  (                                                                                                                      &
               &                                               -radiusScaleFree   +self%lengthResolutionPotentialSquare                                                  &
               &                                   +self%lengthResolutionPotentialSqrtOnePlusSquare                                                                      &
               &                                   *sqrt(      +radiusScaleFree**2+self%lengthResolutionPotentialSquare                                                ) &
               &                                  )                                                                                                                      &
               &                            +self%lengthResolutionPotentialTerm2                                                                                         &
               &                           )                                                                                                                             &
               &                          *self%concentrationPotentialTerm
       end if
    end if
    potential=self%potentialPrevious
    return
  end function sphericalFiniteResolutionNFWPotential
  
  double precision function sphericalFiniteResolutionNFWRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for finite-resolution NFW distributions.
    !!}    
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target   :: self
    double precision                                              , intent(in   ), optional :: mass             , massFractional
    integer         (c_size_t                                    ), dimension(0:1)          :: jLengthResolution
    double precision                                              , dimension(0:1)          :: hLengthResolution
    integer                                                                                 :: iLengthResolution
    double precision                                                                        :: mass_            , massScaleFree

    mass_=0.0d0
    if (present(mass)) then
       mass_=mass
    else if (present(massFractional)) then
       call Error_Report('mass is unbounded, so mass fraction is undefined'//{introspection:location})
    else
       call Error_Report('either mass or massFractional must be supplied'  //{introspection:location})
    end if
    if (mass /= self%radiusEnclosingMassMassPrevious) then
       self%radiusEnclosingMassMassPrevious=mass
       ! Find scale free mass, and the maximum such mass reached in the profile.
       massScaleFree=+     mass                    &
            &        /self%densityNormalization    &
            &        /self%radiusScale         **3
       ! Ensure table is sufficiently extensive.
       call self%radiusEnclosingMassTabulate(massScaleFree,self%lengthResolutionScaleFree)
       ! Interpolate to get the scale free radius enclosing the scale free mass.
       call self%radiusEnclosingMassTableLengthResolutionInterpolator%linearFactors(self%lengthResolutionScaleFree,jLengthResolution(0),hLengthResolution)
       jLengthResolution(1)=jLengthResolution(0)+1
       self%radiusEnclosingMassPrevious=0.0d0
       do iLengthResolution=0,1
          self%radiusEnclosingMassPrevious=+self%radiusEnclosingMassPrevious                                                                                                               &
               &                           +self%radiusEnclosingMassTableMassInterpolator%interpolate(massScaleFree,self%radiusEnclosingMassTable(:,jLengthResolution(iLengthResolution))) &
               &                           *                                                                                                        hLengthResolution(iLengthResolution)
       end do
       self%radiusEnclosingMassPrevious=+self%radiusEnclosingMassPrevious &
            &                           *self%radiusScale
    end if
    radius=self%radiusEnclosingMassPrevious
    return
  end function sphericalFiniteResolutionNFWRadiusEnclosingMass
  
  subroutine sphericalFiniteResolutionNFWRadiusEnclosingMassTabulate(self,mass,lengthResolution)
    !!{
    Tabulates the radius enclosing a given mass for finite resolution NFW mass profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Make_Range               , rangeTypeLogarithmic
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target :: self
    double precision                                              , intent(in   )         :: mass                   , lengthResolution
    double precision                                              , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-9
    logical                                                                               :: retabulate
    integer                                                                               :: iLengthResolution      , iMass                   , &
         &                                                                                   i
    type            (rootFinder                                  )                        :: finder
    
    do i=1,2
       retabulate=.false.
       if (.not.self%radiusEnclosingMassTableInitialized) then
          retabulate=.true.
       else if (                                                                    &
            &    mass             < self%radiusEnclosingMassMassMinimum             &
            &   .or.                                                                &
            &    mass             > self%radiusEnclosingMassMassMaximum             &
            &   .or.                                                                &
            &    lengthResolution < self%radiusEnclosingMassLengthResolutionMinimum &
            &   .or.                                                                &
            &    lengthResolution > self%radiusEnclosingMassLengthResolutionMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (retabulate     .and.i==1) call self%restoreMassTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%radiusEnclosingMassMassMinimum               =min(0.5d0*mass            ,self%radiusEnclosingMassMassMinimum            )
       self%radiusEnclosingMassMassMaximum               =max(2.0d0*mass            ,self%radiusEnclosingMassMassMaximum            )
       self%radiusEnclosingMassLengthResolutionMinimum   =min(0.5d0*lengthResolution,self%radiusEnclosingMassLengthResolutionMinimum)
       self%radiusEnclosingMassLengthResolutionMaximum   =max(2.0d0*lengthResolution,self%radiusEnclosingMassLengthResolutionMaximum)
       self%radiusEnclosingMassTableMassCount            =int(log10(self%radiusEnclosingMassMassMaximum            /self%radiusEnclosingMassMassMinimum            )*dble(radiusEnclosingMassTableMassPointsPerDecade            ))+1
       self%radiusEnclosingMassTableLengthResolutionCount=int(log10(self%radiusEnclosingMassLengthResolutionMaximum/self%radiusEnclosingMassLengthResolutionMinimum)*dble(radiusEnclosingMassTableLengthResolutionPointsPerDecade))+1
       if (allocated(self%radiusEnclosingMassTableMass)) then
          deallocate(self%radiusEnclosingMassTableLengthResolution)
          deallocate(self%radiusEnclosingMassTableMass            )
          deallocate(self%radiusEnclosingMassTable                )
       end if
       allocate(self%radiusEnclosingMassTableLengthResolution(                                       self%radiusEnclosingMassTableLengthResolutionCount))
       allocate(self%radiusEnclosingMassTableMass            (self%radiusEnclosingMassTableMassCount                                                   ))
       allocate(self%radiusEnclosingMassTable                (self%radiusEnclosingMassTableMassCount,self%radiusEnclosingMassTableLengthResolutionCount))
       ! Create a range of radii and core radii.
       self%radiusEnclosingMassTableMass            =Make_Range(self%radiusEnclosingMassMassMinimum            ,self%radiusEnclosingMassMassMaximum            ,self%radiusEnclosingMassTableMassCount            ,rangeType=rangeTypeLogarithmic)
       self%radiusEnclosingMassTableLengthResolution=Make_Range(self%radiusEnclosingMassLengthResolutionMinimum,self%radiusEnclosingMassLengthResolutionMaximum,self%radiusEnclosingMassTableLengthResolutionCount,rangeType=rangeTypeLogarithmic)
       ! Initialize our root finder.
       finder=rootFinder(                                                             &
            &            rootFunction                 =rootMass                     , &
            &            toleranceAbsolute            =toleranceAbsolute            , &
            &            toleranceRelative            =toleranceRelative            , &
            &            rangeExpandDownward          =0.5d0                        , &
            &            rangeExpandUpward            =2.0d0                        , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
            &           )
       ! Loop over mass and core radius and populate tables.
       self_ => self
       do iLengthResolution=1,self%radiusEnclosingMassTableLengthResolutionCount
          iLengthResolution_=iLengthResolution
          do iMass=1,self%radiusEnclosingMassTableMassCount
             iMass_=iMass
             ! Check that the root condition is satisfied at infinitely large radius. If it is not, then no radius encloses the
             ! required mass. Simply set the radius to an infinitely large value in such case.
             if (rootMass(radius=huge(0.0d0)) < 0.0d0) then
                self%radiusEnclosingMassTable(iMass,iLengthResolution)=huge(0.0d0)
             else
                self%radiusEnclosingMassTable(iMass,iLengthResolution)=finder%find(rootGuess=1.0d0)
             end if
          end do
       end do
       ! Build interpolators.
       if (allocated(self%radiusEnclosingMassTableLengthResolutionInterpolator)) deallocate(self%radiusEnclosingMassTableLengthResolutionInterpolator)
       if (allocated(self%radiusEnclosingMassTableMassInterpolator            )) deallocate(self%radiusEnclosingMassTableMassInterpolator            )
       allocate(self%radiusEnclosingMassTableLengthResolutionInterpolator)
       allocate(self%radiusEnclosingMassTableMassInterpolator            )
       self%radiusEnclosingMassTableLengthResolutionInterpolator=interpolator(self%radiusEnclosingMassTableLengthResolution)
       self%radiusEnclosingMassTableMassInterpolator            =interpolator(self%radiusEnclosingMassTableMass            )
       ! Specify that tabulation has been made.
       self%radiusEnclosingMassTableInitialized=.true.
       call self%storeMassTable()
    end if
    return
  end subroutine sphericalFiniteResolutionNFWRadiusEnclosingMassTabulate

  double precision function rootMass(radius)
    !!{
    Root function used in finding the radius enclosing a given mean mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    rootMass=+self_%massEnclosedScaleFree       (radius,self_%radiusEnclosingMassTableLengthResolution(iLengthResolution_)) &
         &   -self_%radiusEnclosingMassTableMass(                                                      iMass_             )
    return
  end function rootMass

  subroutine sphericalFiniteResolutionNFWStoreMassTable(self)
    !!{
    Store the tabulated radius-enclosing-mass data to file.
    !!}
    use :: File_Utilities    , only : File_Lock     , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , char
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5Object                                  )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Mass_'                                          // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    call file%openFile(char(fileName),overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    call file%writeDataset(self%radiusEnclosingMassTableLengthResolution,'lengthResolution')
    call file%writeDataset(self%radiusEnclosingMassTableMass            ,'mass'            )
    call file%writeDataset(self%radiusEnclosingMassTable                ,'radius'          )
    call file%close()
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine sphericalFiniteResolutionNFWStoreMassTable

  subroutine sphericalFiniteResolutionNFWRestoreMassTable(self)
    !!{
    Restore the tabulated radius-enclosing-mass data from file, returning true if successful.
    !!}
    use :: File_Utilities    , only : File_Exists    , File_Lock         , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5Object                                  )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Mass_'                                          // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    if (File_Exists(fileName)) then
       if (allocated(self%radiusEnclosingMassTableMass)) then
          deallocate(self%radiusEnclosingMassTableLengthResolution)
          deallocate(self%radiusEnclosingMassTableMass            )
          deallocate(self%radiusEnclosingMassTable                )
       end if
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName))
       call file%readDataset('lengthResolution',self%radiusEnclosingMassTableLengthResolution)
       call file%readDataset('mass'            ,self%radiusEnclosingMassTableMass            )
       call file%readDataset('radius'          ,self%radiusEnclosingMassTable                )
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%radiusEnclosingMassTableMassCount            =size(self%radiusEnclosingMassTableMass            )
       self%radiusEnclosingMassTableLengthResolutionCount=size(self%radiusEnclosingMassTableLengthResolution)
       self%radiusEnclosingMassMassMinimum               =self%radiusEnclosingMassTableMass            (                                                 1)
       self%radiusEnclosingMassMassMaximum               =self%radiusEnclosingMassTableMass            (self%radiusEnclosingMassTableMassCount            )
       self%radiusEnclosingMassLengthResolutionMinimum   =self%radiusEnclosingMassTableLengthResolution(                                                 1)
       self%radiusEnclosingMassLengthResolutionMaximum   =self%radiusEnclosingMassTableLengthResolution(self%radiusEnclosingMassTableLengthResolutionCount)
       if (allocated(self%radiusEnclosingMassTableLengthResolutionInterpolator)) deallocate(self%radiusEnclosingMassTableLengthResolutionInterpolator)
       if (allocated(self%radiusEnclosingMassTableMassInterpolator            )) deallocate(self%radiusEnclosingMassTableMassInterpolator            )
       allocate(self%radiusEnclosingMassTableLengthResolutionInterpolator)
       allocate(self%radiusEnclosingMassTableMassInterpolator            )
       self%radiusEnclosingMassTableLengthResolutionInterpolator=interpolator(self%radiusEnclosingMassTableLengthResolution)
       self%radiusEnclosingMassTableMassInterpolator            =interpolator(self%radiusEnclosingMassTableMass            )
       self%radiusEnclosingMassTableInitialized                 =.true.
    end if    
    return
  end subroutine sphericalFiniteResolutionNFWRestoreMassTable

  double precision function sphericalFiniteResolutionNFWRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given mean density for finite-resolution NFW mass distributions.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target   :: self
    double precision                                              , intent(in   )           :: density
    double precision                                              , intent(in   ), optional :: radiusGuess
    double precision                                              , parameter               :: epsilonDensity         =1.0d-3
    double precision                                                                        :: densityScaleFreeMaximum       , densityScaleFree
    integer         (c_size_t                                    ), dimension(0:1)          :: jLengthResolution
    double precision                                              , dimension(0:1)          :: hLengthResolution
    integer                                                                                 :: iLengthResolution
    
    if (density /= self%radiusEnclosingDensityDensityPrevious) then
       self%radiusEnclosingDensityDensityPrevious=density
       ! Find scale free density, and the maximum such density reached in the profile.
       densityScaleFree       =+     density                   &
            &                  /self%densityNormalization
       densityScaleFreeMaximum=+1.0d0                          &
            &                  /self%lengthResolutionScaleFree
       if      (densityScaleFree >= densityScaleFreeMaximum) then
          ! Maximum density is exceeded - return zero radius.
          self%radiusEnclosingDensityPrevious=0.0d0
       else if (densityScaleFree >= densityScaleFreeMaximum*(1.0d0-epsilonDensity)) then
          ! For densities close to the maximum density, use a series solution.
          self%radiusEnclosingDensityPrevious=+0.5d0                     &
               &                              *(                         &
               &                                +1.0d0                   &
               &                                -densityScaleFree        &
               &                                /densityScaleFreeMaximum &
               &                               )                         &
               &                              *self%radiusScale
       else
          ! Use a tabulated solution in other regimes.   
          ! Ensure table is sufficiently extensive.
          call self%radiusEnclosingDensityTabulate(densityScaleFree,self%lengthResolutionScaleFree)
          ! Interpolate to get the scale free radius enclosing the scale free density.
          call self%radiusEnclosingDensityTableLengthResolutionInterpolator%linearFactors(self%lengthResolutionScaleFree,jLengthResolution(0),hLengthResolution)
          jLengthResolution(1)=jLengthResolution(0)+1
          self%radiusEnclosingDensityPrevious=0.0d0
          do iLengthResolution=0,1
             self%radiusEnclosingDensityPrevious=+self%radiusEnclosingDensityPrevious                                                                                                                        &
                  &                              +self%radiusEnclosingDensityTableDensityInterpolator%interpolate(densityScaleFree,self%radiusEnclosingDensityTable(:,jLengthResolution(iLengthResolution))) &
                  &                              *                                                                                                                    hLengthResolution(iLengthResolution)
         end do
          self%radiusEnclosingDensityPrevious=+self%radiusEnclosingDensityPrevious &
               &                              *self%radiusScale
       end if
    end if
    radius=self%radiusEnclosingDensityPrevious
    return
  end function sphericalFiniteResolutionNFWRadiusEnclosingDensity
  
  subroutine sphericalFiniteResolutionNFWRadiusEnclosingDensityTabulate(self,density,lengthResolution)
    !!{
    Tabulates the radius enclosing a given density for finite resolution NFW density profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Make_Range               , rangeTypeLogarithmic
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target :: self
    double precision                                              , intent(in   )         :: density                , lengthResolution
    double precision                                              , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-9
    logical                                                                               :: retabulate
    integer                                                                               :: iLengthResolution      , iDensity                , &
         &                                                                                   i
    type            (rootFinder                                  )                        :: finder

    do i=1,2
       retabulate=.false.
       if (.not.self%radiusEnclosingDensityTableInitialized) then
          retabulate=.true.
       else if (                                                                       &
            &    density          < self%radiusEnclosingDensityDensityMinimum          &
            &   .or.                                                                   &
            &    density          > self%radiusEnclosingDensityDensityMaximum          &
            &   .or.                                                                   &
            &    lengthResolution < self%radiusEnclosingDensityLengthResolutionMinimum &
            &   .or.                                                                   &
            &    lengthResolution > self%radiusEnclosingDensityLengthResolutionMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (retabulate     .and.i==1) call self%restoreDensityTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%radiusEnclosingDensityDensityMinimum            =min(0.5d0*density         ,self%radiusEnclosingDensityDensityMinimum         )
       self%radiusEnclosingDensityDensityMaximum            =max(2.0d0*density         ,self%radiusEnclosingDensityDensityMaximum         )
       self%radiusEnclosingDensityLengthResolutionMinimum   =min(0.5d0*lengthResolution,self%radiusEnclosingDensityLengthResolutionMinimum)
       self%radiusEnclosingDensityLengthResolutionMaximum   =max(2.0d0*lengthResolution,self%radiusEnclosingDensityLengthResolutionMaximum)
       self%radiusEnclosingDensityTableDensityCount         =int(log10(self%radiusEnclosingDensityDensityMaximum         /self%radiusEnclosingDensityDensityMinimum         )*dble(radiusEnclosingDensityTableDensityPointsPerDecade         ))+1
       self%radiusEnclosingDensityTableLengthResolutionCount=int(log10(self%radiusEnclosingDensityLengthResolutionMaximum/self%radiusEnclosingDensityLengthResolutionMinimum)*dble(radiusEnclosingDensityTableLengthResolutionPointsPerDecade))+1
       if (allocated(self%radiusEnclosingDensityTableDensity)) then
          deallocate(self%radiusEnclosingDensityTableLengthResolution)
          deallocate(self%radiusEnclosingDensityTableDensity         )
          deallocate(self%radiusEnclosingDensityTable                )
       end if
       allocate(self%radiusEnclosingDensityTableLengthResolution(                                             self%radiusEnclosingDensityTableLengthResolutionCount))
       allocate(self%radiusEnclosingDensityTableDensity         (self%radiusEnclosingDensityTableDensityCount                                                      ))
       allocate(self%radiusEnclosingDensityTable                (self%radiusEnclosingDensityTabledensityCount,self%radiusEnclosingDensityTableLengthResolutionCount))
       ! Create a range of radii and core radii.
       self%radiusEnclosingDensityTableDensity         =Make_Range(self%radiusEnclosingDensityDensityMinimum         ,self%radiusEnclosingDensityDensityMaximum         ,self%radiusEnclosingDensityTableDensityCount         ,rangeType=rangeTypeLogarithmic)
       self%radiusEnclosingDensityTableLengthResolution=Make_Range(self%radiusEnclosingDensityLengthResolutionMinimum,self%radiusEnclosingDensityLengthResolutionMaximum,self%radiusEnclosingDensityTableLengthResolutionCount,rangeType=rangeTypeLogarithmic)
       ! Initialize our root finder.
       finder=rootFinder(                                                             &
            &            rootFunction                 =rootDensity                  , &
            &            toleranceAbsolute            =toleranceAbsolute            , &
            &            toleranceRelative            =toleranceRelative            , &
            &            rangeExpandDownward          =0.5d0                        , &
            &            rangeExpandUpward            =2.0d0                        , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
            &           )
       ! Loop over density and core radius and populate tables.
       self_ => self
       do iLengthResolution=1,self%radiusEnclosingDensityTableLengthResolutionCount
          iLengthResolution_=iLengthResolution
          do iDensity=1,self%radiusEnclosingDensityTableDensityCount
             iDensity_=iDensity
             if (self%radiusEnclosingDensityTableDensity(iDensity) > 1.0d0/self%radiusEnclosingDensityTableLengthResolution(iLengthResolution)) then
                ! Density exceeds the maximum density in the profile - so set zero radius.
                self%radiusEnclosingDensityTable(iDensity,iLengthResolution)=0.0d0
             else
                self%radiusEnclosingDensityTable(iDensity,iLengthResolution)=finder%find(rootGuess=1.0d0)
             end if
          end do
       end do
       ! Build interpolators.
       if (allocated(self%radiusEnclosingDensityTableLengthResolutionInterpolator)) deallocate(self%radiusEnclosingDensityTableLengthResolutionInterpolator)
       if (allocated(self%radiusEnclosingDensityTableDensityInterpolator         )) deallocate(self%radiusEnclosingDensityTableDensityInterpolator         )
       allocate(self%radiusEnclosingDensityTableLengthResolutionInterpolator)
       allocate(self%radiusEnclosingDensityTableDensityInterpolator         )
       self%radiusEnclosingDensityTableLengthResolutionInterpolator=interpolator(self%radiusEnclosingDensityTableLengthResolution)
       self%radiusEnclosingDensityTableDensityInterpolator         =interpolator(self%radiusEnclosingDensityTableDensity         )
       ! Specify that tabulation has been made.
       self%radiusEnclosingDensityTableInitialized=.true.
       call self%storeDensityTable()
    end if
    return
  end subroutine sphericalFiniteResolutionNFWRadiusEnclosingDensityTabulate

  double precision function rootDensity(radius)
    !!{
    Root function used in finding the radius enclosing a given mean density.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    rootDensity=+3.0d0                                                                                                                     &
         &      *self_%massEnclosedScaleFree             (radius,self_%radiusEnclosingDensityTableLengthResolution(iLengthResolution_))    &
         &      /4.0d0                                                                                                                     &
         &      /Pi                                                                                                                        &
         &      /                                         radius                                                                       **3 &
         &      -self_%radiusEnclosingDensityTableDensity(                                                         iDensity_          )
    return
  end function rootDensity

  subroutine sphericalFiniteResolutionNFWStoreDensityTable(self)
    !!{
    Store the tabulated radius-enclosing-density data to file.
    !!}
    use :: File_Utilities    , only : File_Lock     , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , char
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5Object                                  )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Density_'                                       // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    call file%openFile(char(fileName),overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    call file%writeDataset(self%radiusEnclosingDensityTableLengthResolution,'lengthResolution')
    call file%writeDataset(self%radiusEnclosingDensityTableDensity         ,'density'         )
    call file%writeDataset(self%radiusEnclosingDensityTable                ,'radius'          )
    call file%close()
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine sphericalFiniteResolutionNFWStoreDensityTable

  subroutine sphericalFiniteResolutionNFWRestoreDensityTable(self)
    !!{
    Restore the tabulated radius-enclosing-density data from file, returning true if successful.
    !!}
    use :: File_Utilities    , only : File_Exists    , File_Lock         , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5Object                                  )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Density_'                                       // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    if (File_Exists(fileName)) then
       if (allocated(self%radiusEnclosingDensityTableDensity)) then
          deallocate(self%radiusEnclosingDensityTableLengthResolution)
          deallocate(self%radiusEnclosingDensityTableDensity   )
          deallocate(self%radiusEnclosingDensityTable          )
       end if
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName))
       call file%readDataset('lengthResolution',self%radiusEnclosingDensityTableLengthResolution)
       call file%readDataset('density'   ,self%radiusEnclosingDensityTableDensity   )
       call file%readDataset('radius'    ,self%radiusEnclosingDensityTable          )
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%radiusEnclosingDensityTableDensityCount         =size(self%radiusEnclosingDensityTableDensity         )
       self%radiusEnclosingDensityTableLengthResolutionCount=size(self%radiusEnclosingDensityTableLengthResolution)
       self%radiusEnclosingDensityDensityMinimum            =self%radiusEnclosingDensityTableDensity         (                                                    1)
       self%radiusEnclosingDensityDensityMaximum            =self%radiusEnclosingDensityTableDensity         (self%radiusEnclosingDensityTableDensityCount         )
       self%radiusEnclosingDensityLengthResolutionMinimum   =self%radiusEnclosingDensityTableLengthResolution(                                                    1)
       self%radiusEnclosingDensityLengthResolutionMaximum   =self%radiusEnclosingDensityTableLengthResolution(self%radiusEnclosingDensityTableLengthResolutionCount)
       if (allocated(self%radiusEnclosingDensityTableLengthResolutionInterpolator)) deallocate(self%radiusEnclosingDensityTableLengthResolutionInterpolator)
       if (allocated(self%radiusEnclosingDensityTableDensityInterpolator         )) deallocate(self%radiusEnclosingDensityTableDensityInterpolator         )
       allocate(self%radiusEnclosingDensityTableLengthResolutionInterpolator)
       allocate(self%radiusEnclosingDensityTableDensityInterpolator         )
       self%radiusEnclosingDensityTableLengthResolutionInterpolator=interpolator(self%radiusEnclosingDensityTableLengthResolution)
       self%radiusEnclosingDensityTableDensityInterpolator         =interpolator(self%radiusEnclosingDensityTableDensity         )
       self%radiusEnclosingDensityTableInitialized                 =.true.
    end if    
    return
  end subroutine sphericalFiniteResolutionNFWRestoreDensityTable

  double precision function sphericalFiniteResolutionNFWEnergy(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{
    Compute the energy within a given {\normalfont \ttfamily radius} in a finite-resolution NFW mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout) , target :: self
    double precision                                              , intent(in   )          :: radiusOuter
    class           (massDistributionClass                       ), intent(inout) , target :: massDistributionEmbedding
    integer         (c_size_t                                    ), dimension(0:1)         :: jLengthResolution
    double precision                                              , dimension(0:1)         :: hLengthResolution
    integer                                                                                :: iLengthResolution

    if (self%energyPrevious > 0.0d0) then
       ! Ensure table is sufficiently extensive.
       call self%energyTabulate(self%lengthResolutionScaleFree,radiusOuter/self%radiusScale)
       ! Interpolate to get the scale free energy.
       call self%energyTableLengthResolutionInterpolator%linearFactors(self%lengthResolutionScaleFree,jLengthResolution(0),hLengthResolution)
       jLengthResolution(1)=jLengthResolution(0)+1
       self%energyPrevious=0.0d0
       do iLengthResolution=0,1
          self%energyPrevious=+self%energyPrevious                                                                                                                        &
               &              +self%energyTableRadiusOuterInterpolator%interpolate(radiusOuter/self%radiusScale,self%energyTable(:,jLengthResolution(iLengthResolution))) &
               &              *                                                                                                    hLengthResolution(iLengthResolution)
       end do
       self%energyPrevious=+self             %energyPrevious          &
            &              *gravitationalConstant_internal            &
            &              *self             %densityNormalization**2 &
            &              *self             %radiusScale         **5
    end if
    energy=self%energyPrevious
    return
  end function sphericalFiniteResolutionNFWEnergy
  
  subroutine sphericalFiniteResolutionNFWEnergyTabulate(self,lengthResolution,radiusOuter)
    !!{
    Tabulates the energy for finite resolution NFW mass profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range, rangeTypeLogarithmic
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target :: self
    double precision                                              , intent(in   )         :: radiusOuter                , lengthResolution
    double precision                                              , parameter             :: multiplierRadius   =100.0d0
    type            (integrator                                  )                        :: integratorPotential        , integratorKinetic, &
         &                                                                                   integratorPressure
    double precision                                                                      :: pseudoPressure             , energyKinetic    , &
         &                                                                                   energyPotential            , radiusOuter_
    logical                                                                               :: retabulate
    integer                                                                               :: iLengthResolution          , iRadiusOuter     , &
         &                                                                                   i

    do i=1,2
       retabulate=.false.
       if (.not.self%energyTableInitialized) then
          retabulate=.true.
       else if (                                                       &
            &    radiusOuter      < self%energyRadiusOuterMinimum      &
            &   .or.                                                   &
            &    radiusOuter      > self%energyRadiusOuterMaximum      &
            &   .or.                                                   &
            &    lengthResolution < self%energyLengthResolutionMinimum &
            &   .or.                                                   &
            &    lengthResolution > self%energyLengthResolutionMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (     retabulate.and.i==1) call self%restoreEnergyTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       self%energyRadiusOuterMinimum        =min(0.5d0*radiusOuter     ,self%energyRadiusOuterMinimum     )
       self%energyRadiusOuterMaximum        =max(2.0d0*radiusOuter     ,self%energyRadiusOuterMaximum     )
       self%energyLengthResolutionMinimum   =min(0.5d0*lengthResolution,self%energyLengthResolutionMinimum)
       self%energyLengthResolutionMaximum   =max(2.0d0*lengthResolution,self%energyLengthResolutionMaximum)
       self%energyTableRadiusOuterCount     =int(log10(self%energyRadiusOuterMaximum     /self%energyRadiusOuterMinimum     )*dble(energyTableRadiusOuterPointsPerDecade     ))+1
       self%energyTableLengthResolutionCount=int(log10(self%energyLengthResolutionMaximum/self%energyLengthResolutionMinimum)*dble(energyTableLengthResolutionPointsPerDecade))+1
       if (allocated(self%energyTableRadiusOuter)) then
          deallocate(self%energyTableLengthResolution)
          deallocate(self%energyTableRadiusOuter     )
          deallocate(self%energyTable                )
       end if
       allocate(self%energyTableLengthResolution(                                 self%energyTableLengthResolutionCount))
       allocate(self%energyTableRadiusOuter     (self%energyTableRadiusOuterCount                                      ))
       allocate(self%energyTable                (self%energyTableradiusOuterCount,self%energyTableLengthResolutionCount))
       ! Create a range of radii and core radii.
       self%energyTableRadiusOuter     =Make_Range(self%energyRadiusOuterMinimum     ,self%energyRadiusOuterMaximum     ,self%energyTableRadiusOuterCount     ,rangeType=rangeTypeLogarithmic)
       self%energyTableLengthResolution=Make_Range(self%energyLengthResolutionMinimum,self%energyLengthResolutionMaximum,self%energyTableLengthResolutionCount,rangeType=rangeTypeLogarithmic)
       ! Initialize integrators.
       integratorPotential=integrator(integrandEnergyPotential,toleranceRelative=1.0d-3)
       integratorKinetic  =integrator(integrandEnergyKinetic  ,toleranceRelative=1.0d-3)
       integratorPressure =integrator(integrandPseudoPressure ,toleranceRelative=1.0d-3)
       ! Loop over radiusOuter and core radius and populate tables.
       self_ => self
       do iLengthResolution=1,self%energyTableLengthResolutionCount
          iLengthResolution_=iLengthResolution
          do iRadiusOuter=1,self%energyTableRadiusOuterCount
             radiusOuter_                                    =self%energyTableRadiusOuter(iRadiusOuter)
             energyPotential                                 =+integratorPotential%integrate(       0.0d0,                 radiusOuter_)
             energyKinetic                                   =+integratorKinetic  %integrate(       0.0d0,                 radiusOuter_)
             pseudoPressure                                  =+integratorPressure %integrate(radiusOuter_,multiplierRadius*radiusOuter_)
             self%energyTable(iRadiusOuter,iLengthResolution)=-0.5d0                                                                                             &
                  &                                           *(                                                                                                 &
                  &                                             +energyPotential                                                                                 &
                  &                                             +self%massEnclosedScaleFree(radiusOuter_,self%energyTableLengthResolution(iLengthResolution))**2 &
                  &                                             /radiusOuter_                                                                                    &
                  &                                            )                                                                                                 &
                  &                                           +2.0d0                                                                                             &
                  &                                           *Pi                                                                                                &
                  &                                           *(                                                                                                 &
                  &                                             +radiusOuter_                                                                                **3 &
                  &                                             *pseudoPressure                                                                                  &
                  &                                             +energyKinetic                                                                                   &
                  &                                            )
            end do
       end do
       ! Build interpolators.
       if (allocated(self%energyTableLengthResolutionInterpolator)) deallocate(self%energyTableLengthResolutionInterpolator)
       if (allocated(self%energyTableRadiusOuterInterpolator     )) deallocate(self%energyTableRadiusOuterInterpolator     )
       allocate(self%energyTableLengthResolutionInterpolator)
       allocate(self%energyTableRadiusOuterInterpolator     )
       self%energyTableLengthResolutionInterpolator=interpolator(self%energyTableLengthResolution)
       self%energyTableRadiusOuterInterpolator     =interpolator(self%energyTableRadiusOuter     )
       ! Specify that tabulation has been made.
       self%energyTableInitialized=.true.
       call self%storeEnergyTable()
    end if
    return
  end subroutine sphericalFiniteResolutionNFWEnergyTabulate

  double precision function integrandEnergyPotential(radius)
    !!{
    Integrand for potential energy of the halo.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       integrandEnergyPotential=(                                                                                           &
            &                    +self_%massEnclosedScaleFree(radius,self_%energyTableLengthResolution(iLengthResolution_)) &
            &                    /                            radius                                                        &
            &                   )**2
    else
       integrandEnergyPotential=+0.0d0
    end if
    return
  end function integrandEnergyPotential
  
  double precision function integrandEnergyKinetic(radius)
    !!{
    Integrand for kinetic energy of the halo.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       integrandEnergyKinetic=+self_%massEnclosedScaleFree(radius,self_%energyTableLengthResolution(iLengthResolution_)) &
            &                 *self_%densityScaleFree     (radius,self_%energyTableLengthResolution(iLengthResolution_)) &
            &                 *                            radius
    else
       integrandEnergyKinetic=+0.0d0
    end if
    return
  end function integrandEnergyKinetic
  
  double precision function integrandPseudoPressure(radius)
    !!{
    Integrand for pseudo-pressure ($\rho(r) \sigma^2(r)$) of the halo.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       integrandPseudoPressure=+self_%massEnclosedScaleFree(radius,self_%energyTableLengthResolution(iLengthResolution_))    &
            &                  *self_%densityScaleFree     (radius,self_%energyTableLengthResolution(iLengthResolution_))    &
            &                  /                            radius                                                       **2
    else
       integrandPseudoPressure=+0.0d0
    end if
    return
  end function integrandPseudoPressure

  double precision function sphericalFiniteResolutionNFWDensityScaleFree(self,radius,radiusCore) result(densityScaleFree)
    !!{
    Returns the scale-free density in the dark matter profile at the given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    double precision                                              , intent(in   ) :: radius, radiusCore
    !$GLC attributes unused :: self
    
    densityScaleFree=1.0d0/(1.0d0+radius)**2/sqrt(radius**2+radiusCore**2)
    return
  end function sphericalFiniteResolutionNFWDensityScaleFree
    
  subroutine sphericalFiniteResolutionNFWStoreEnergyTable(self)
    !!{
    Store the tabulated energy data to file.
    !!}
    use :: File_Utilities    , only : File_Lock     , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)       , char
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5Object                                  )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Energy_'                                        // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    call Directory_Make(char(File_Path(char(fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    call file%openFile(char(fileName),overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    call file%writeDataset(self%energyTableLengthResolution,'lengthResolution')
    call file%writeDataset(self%energyTableRadiusOuter     ,'radiusOuter'     )
    call file%writeDataset(self%energyTable                ,'energy'          )
    call file%close()
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine sphericalFiniteResolutionNFWStoreEnergyTable

  subroutine sphericalFiniteResolutionNFWRestoreEnergyTable(self)
    !!{
    Restore the tabulated radius-enclosing-mass data from file, returning true if successful.
    !!}
    use :: File_Utilities    , only : File_Exists    , File_Lock         , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5Object                                  )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)                   // &
         &   'darkMatter/'                                    // &
         &   self%objectType      (                          )// &
         &   'Energy_'                                        // &
         &   self%hashedDescriptor(includeSourceDigest=.true.)// &
         &   '.hdf5'
    if (File_Exists(fileName)) then
       if (allocated(self%energyTableRadiusOuter)) then
          deallocate(self%energyTableLengthResolution   )
          deallocate(self%energyTableRadiusOuter)
          deallocate(self%energyTable             )
       end if
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName))
       call file%readDataset('lengthResolution',self%energyTableLengthResolution)
       call file%readDataset('radiusOuter'     ,self%energyTableRadiusOuter     )
       call file%readDataset('energy'          ,self%energyTable                )
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%energyTableRadiusOuterCount     =size(self%energyTableRadiusOuter      )
       self%energyTableLengthResolutionCount=size(self%energyTableLengthResolution)
       self%energyRadiusOuterMinimum        =self%energyTableRadiusOuter     (                                    1)
       self%energyRadiusOuterMaximum        =self%energyTableRadiusOuter     (self%energyTableRadiusOuterCount     )
       self%energyLengthResolutionMinimum   =self%energyTableLengthResolution(                                    1)
       self%energyLengthResolutionMaximum   =self%energyTableLengthResolution(self%energyTableLengthResolutionCount)
       if (allocated(self%energyTableLengthResolutionInterpolator)) deallocate(self%energyTableLengthResolutionInterpolator)
       if (allocated(self%energyTableRadiusOuterInterpolator     )) deallocate(self%energyTableRadiusOuterInterpolator     )
       allocate(self%energyTableLengthResolutionInterpolator)
       allocate(self%energyTableRadiusOuterInterpolator     )
       self%energyTableLengthResolutionInterpolator=interpolator(self%energyTableLengthResolution)
       self%energyTableRadiusOuterInterpolator     =interpolator(self%energyTableRadiusOuter     )
       self%energyTableInitialized                 =.true.
    end if    
    return
  end subroutine sphericalFiniteResolutionNFWRestoreEnergyTable
