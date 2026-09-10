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

  !+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a finite resolution NFW spherical mass distribution.
  !!}

  use :: Numerical_Interpolation, only : interpolator
  use :: Numerical_Ranges       , only : rangeLattice

  !![
  <massDistribution name="massDistributionSphericalFiniteResolutionNFW" docformat="rst">
   <description>
   A mass distribution class which applies a finite resolution to an NFW density profile, typically to mimic the effects of finite resolution in an N-body simulation. Specifically, the density profile is given by

   .. math::

      \rho(r) = \rho_\mathrm{NFW}(r) \left( 1 + \left[ \frac{\Delta x}{r} \right]^2 \right)^{-1/2},

   where :math:`\Delta x` is the larger of the resolution length, ``[lengthResolution]``, and the radius in the original profile enclosing the mass resolution, ``[massResolution]``.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSpherical) :: massDistributionSphericalFiniteResolutionNFW
     !!{RST
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
     ! Enclosed mass quantities.
     double precision                                            :: lengthResolutionScaleFreeLowerTerm                     , lengthResolutionScaleFreeSquared              , &
          &                                                         lengthResolutionScaleFreeCubed                         , lengthResolutionScaleFreeOnePlusTerm          , &
          &                                                         lengthResolutionScaleFreeOnePlus2Term                  , lengthResolutionScaleFreeSqrtTerm             , &
          &                                                         lengthResolutionScaleFreeSqrt2Term                     , lengthResolutionScaleFreeSqrtCubedTerm        , &
          &                                                         lengthResolutionScaleFreePrevious
   contains
     !![
     <methods docformat="rst">
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
       <method method="suffix"                         description="Return a file name suffix (containing a source code digest."                                       />
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
     procedure :: suffix                         => finiteResolutionNFWSuffix
  end type massDistributionSphericalFiniteResolutionNFW

  interface massDistributionSphericalFiniteResolutionNFW
     !!{RST
     Constructors for the :galacticus-class:`massDistributionSphericalFiniteResolutionNFW` mass distribution class.
     !!}
     module procedure sphericalFiniteResolutionNFWConstructorParameters
     module procedure sphericalFiniteResolutionNFWConstructorInternal
  end interface massDistributionSphericalFiniteResolutionNFW

  ! Tabulation resolution parameters.
  integer                                                       , parameter                   :: radiusEnclosingDensityTableDensityPointsPerDecade                =100
  integer                                                       , parameter                   :: radiusEnclosingDensityTableLengthResolutionPointsPerDecade       =100
  integer                                                       , parameter                   :: radiusEnclosingMassTableMassPointsPerDecade                      =100
  integer                                                       , parameter                   :: radiusEnclosingMassTableLengthResolutionPointsPerDecade          =100
  integer                                                       , parameter                   :: energyTableRadiusOuterPointsPerDecade                            =100
  integer                                                       , parameter                   :: energyTableLengthResolutionPointsPerDecade                       =100

  ! Interval, in lattice steps, to which the bounds of every axis of these tabulations are pinned - that is, every bound is
  ! rounded outward to a tenth of a decade. Whole decades are not used because each of these tabulations is two-dimensional and
  ! its cost is the product of its two extents, and because the axes here reach a long way: measured against the extents an
  ! accumulated velocity dispersion tabulation actually reached (nine decades of radius by 1.3 decades of length resolution),
  ! whole decades would inflate the number of points by 73%, where tenths of a decade cost 3%. The margin applied to a request
  ! is a factor of two either side, so a tenth of a decade is also small compared with the width of the window being pinned.
  integer                                                       , parameter                   :: tabulationAnchorsPerDecade                                       = 10

  ! Largest radius for precise arctanh() evaluation.
  double precision                                              , parameter                   :: radiusScaleFreeLargeATanh                                        =1.0d+6

  ! Radius-enclosing-density tabulation.
  logical                                                                                     :: radiusEnclosingDensityTableInitialized                           =.false.
  integer                                                                                     :: radiusEnclosingDensityTableLengthResolutionCount                              , radiusEnclosingDensityTableDensityCount
  double precision                                              , allocatable, dimension(:  ) :: radiusEnclosingDensityTableLengthResolution                                   , radiusEnclosingDensityTableDensity
  double precision                                              , allocatable, dimension(:,:) :: radiusEnclosingDensityTable
  type            (interpolator                                ), allocatable                 :: radiusEnclosingDensityTableLengthResolutionInterpolator                       , radiusEnclosingDensityTableDensityInterpolator
  double precision                                                                            :: radiusEnclosingDensityDensityMinimum                             =+huge(0.0d0), radiusEnclosingDensityDensityMaximum          =-huge(0.0d0), &
       &                                                                                         radiusEnclosingDensityLengthResolutionMinimum                    =+huge(0.0d0), radiusEnclosingDensityLengthResolutionMaximum =-huge(0.0d0)
  type            (rangeLattice                                )                             :: radiusEnclosingDensityLatticeDensity
  type            (rangeLattice                                )                             :: radiusEnclosingDensityLatticeLengthResolution
  !$omp threadprivate(radiusEnclosingDensityLatticeDensity,radiusEnclosingDensityLatticeLengthResolution,radiusEnclosingDensityTableInitialized,radiusEnclosingDensityTableLengthResolutionCount,radiusEnclosingDensityTableDensityCount,radiusEnclosingDensityTableLengthResolution, radiusEnclosingDensityTableDensity,radiusEnclosingDensityTable,radiusEnclosingDensityTableLengthResolutionInterpolator,radiusEnclosingDensityTableDensityInterpolator,radiusEnclosingDensityDensityMinimum,radiusEnclosingDensityDensityMaximum,radiusEnclosingDensityLengthResolutionMinimum,radiusEnclosingDensityLengthResolutionMaximum)
  
  ! Radius-enclosing-mass tabulation.
  logical                                                                                     :: radiusEnclosingMassTableInitialized                              =.false.
  integer                                                                                     :: radiusEnclosingMassTableLengthResolutionCount                                 , radiusEnclosingMassTableMassCount
  double precision                                              , allocatable, dimension(:  ) :: radiusEnclosingMassTableLengthResolution                                      , radiusEnclosingMassTableMass
  double precision                                              , allocatable, dimension(:,:) :: radiusEnclosingMassTable
  type            (interpolator                                ), allocatable                 :: radiusEnclosingMassTableLengthResolutionInterpolator                          , radiusEnclosingMassTableMassInterpolator
  double precision                                                                            :: radiusEnclosingMassMassMinimum                                   =+huge(0.0d0), radiusEnclosingMassMassMaximum                =-huge(0.0d0), &
       &                                                                                         radiusEnclosingMassLengthResolutionMinimum                       =+huge(0.0d0), radiusEnclosingMassLengthResolutionMaximum    =-huge(0.0d0)
  type            (rangeLattice                                )                             :: radiusEnclosingMassLatticeMass
  type            (rangeLattice                                )                             :: radiusEnclosingMassLatticeLengthResolution
  !$omp threadprivate(radiusEnclosingMassLatticeMass,radiusEnclosingMassLatticeLengthResolution,radiusEnclosingMassTableInitialized,radiusEnclosingMassTableLengthResolutionCount,radiusEnclosingMassTableMassCount,radiusEnclosingMassTableLengthResolution,radiusEnclosingMassTableMass,radiusEnclosingMassTable,radiusEnclosingMassTableLengthResolutionInterpolator,radiusEnclosingMassTableMassInterpolator,radiusEnclosingMassMassMinimum,radiusEnclosingMassMassMaximum,radiusEnclosingMassLengthResolutionMinimum,radiusEnclosingMassLengthResolutionMaximum)
  
  ! Energy tabulation.
  logical                                                                                     :: energyTableInitialized                                           =.false.
  integer                                                                                     :: energyTableLengthResolutionCount                                              , energyTableRadiusOuterCount
  double precision                                              , allocatable, dimension(:  ) :: energyTableLengthResolution                                                   , energyTableRadiusOuter
  double precision                                              , allocatable, dimension(:,:) :: energyTable
  type            (interpolator                                ), allocatable                 :: energyTableLengthResolutionInterpolator                                       , energyTableRadiusOuterInterpolator
  double precision                                                                            :: energyRadiusOuterMinimum                                         =+huge(0.0d0), energyRadiusOuterMaximum                      =-huge(0.0d0), &
       &                                                                                         energyLengthResolutionMinimum                                    =+huge(0.0d0), energyLengthResolutionMaximum                 =-huge(0.0d0)
  type            (rangeLattice                                )                             :: energyLatticeRadiusOuter
  type            (rangeLattice                                )                             :: energyLatticeLengthResolution
  !$omp threadprivate(energyLatticeRadiusOuter,energyLatticeLengthResolution,energyTableInitialized,energyTableLengthResolutionCount,energyTableRadiusOuterCount,energyTableLengthResolution,energyTableRadiusOuter,energyTable,energyTableLengthResolutionInterpolator,energyTableRadiusOuterInterpolator,energyRadiusOuterMinimum,energyRadiusOuterMaximum,energyLengthResolutionMinimum,energyLengthResolutionMaximum) 
  
  ! Sub-module-scope variables used in integrations.
  class           (massDistributionSphericalFiniteResolutionNFW), pointer                     :: self_
  integer                                                                                     :: iLengthResolution_                                               , iDensity_, &
       &                                                                                         iMass_
  !$omp threadprivate(self_,iLengthResolution_,iDensity_,iMass_)

  ! Generate a source digest.
  !![
  <sourceDigest name="massDistributionFiniteResolutionNFWSourceDigest"/>
  !!]

contains

  function sphericalFiniteResolutionNFWConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionSphericalFiniteResolutionNFW` mass distribution class which builds the object from a parameter set.
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
    <inputParameter docformat="rst">
      <name>lengthResolution</name>
      <source>parameters</source>
      <description>
      The spatial resolution length scale (in Mpc) of the N-body simulation being modeled; sets the minimum effective radius below which the NFW density profile is softened.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusScale</name>
      <source>parameters</source>
      <description>
      The NFW scale radius (in Mpc) at which the density profile transitions from the inner :math:`\rho \propto r^{-1}` slope to the outer :math:`\rho \propto r^{-3}` slope.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusVirial</name>
      <source>parameters</source>
      <description>
      The virial radius (in Mpc) of the halo, defining the outer boundary of the NFW profile at which the mean enclosed density equals the virial overdensity threshold.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>mass</name>
      <source>parameters</source>
      <description>
      The total mass (in :math:`\mathrm{M}_\odot`) enclosed within the virial radius, used together with ``radiusScale`` and ``radiusVirial`` to normalize the NFW density profile.
      </description>
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
    self=massDistributionSphericalFiniteResolutionNFW(lengthResolution,radiusScale,radiusVirial,mass,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sphericalFiniteResolutionNFWConstructorParameters
  
  function sphericalFiniteResolutionNFWConstructorInternal(lengthResolution,radiusScale,radiusVirial,mass,componentType,massType) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionSphericalFiniteResolutionNFW` mass distribution class.
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

    self%dimensionless                             =.false.
    self%lengthResolutionScalefreePrevious         =-huge(0.0d0)
    self%massEnclosedMassPrevious                  =-huge(0.0d0)
    self%massEnclosedRadiusPrevious                =-huge(0.0d0)
    self%potentialPrevious                         =-huge(0.0d0)
    self%potentialRadiusPrevious                   =-huge(0.0d0)
    self%lengthResolutionScaleFreePotentialPrevious=-huge(0.0d0)
    self%concentrationPotentialTerm                =-huge(0.0d0)
    self%densityRadiusPrevious                     =-huge(0.0d0)
    self%densityPrevious                           =-huge(0.0d0)
    self%densityNormalizationPrevious              =-huge(0.0d0)
    self%radiusEnclosingDensityDensityPrevious     =-huge(0.0d0)
    self%radiusEnclosingDensityPrevious            =-huge(0.0d0)
    self%radiusEnclosingMassMassPrevious           =-huge(0.0d0)
    self%radiusEnclosingMassPrevious               =-huge(0.0d0)
    self%energyPrevious                            =+huge(0.0d0)
    ! Construct profile quantities.
    radiusScaleFree                                =+    radiusVirial/radiusScale
    self%lengthResolutionScaleFree                 =+lengthResolution/radiusScale
    self%densityNormalization                      =+mass/4.0d0/Pi/radiusScale**3/(log(1.0d0+radiusScaleFree)-radiusScaleFree/(1.0d0+radiusScaleFree))
    return
  end function sphericalFiniteResolutionNFWConstructorInternal

  double precision function sphericalFiniteResolutionNFWDensity(self,coordinates) result(density)
    !!{RST
    Return the density at the specified ``coordinates`` in a scaled spherical mass distribution.
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
    !!{RST
    Return the density at the specified ``coordinates`` in a finiteResolution spherical mass distribution.
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
    !!{RST
    Returns the enclosed mass (in :math:`\mathrm{M}_\odot`) at the given ``radius`` (given in units of Mpc). The analytic solution (computed using Mathematica) is

    .. math::

       M(x) = 4 \pi \rho_0 r_\mathrm{s}^3 \left[ -\frac{\sqrt{x^2+X^2}}{(1+x) \left(1+X^2\right)}+\tanh ^{-1}\left(\frac{x}{\sqrt{x^2+X^2}}\right)+\frac{\left(1+2X^2\right) \tanh ^{-1}\left(\frac{X^2-x}{\sqrt{1+X^2} \sqrt{x^2+X^2}}\right)}{\left(1+X^2\right)^{3/2}} -\frac{\left(1 + 2 X^2\right) \tanh ^{-1}\left(\sqrt{\frac{X^2}{1 + X^2}}\right)}{\left(1+ X^2\right)^{3/2}}+\frac{\sqrt{X^2}}{1 + X^2} \right],

    where :math:`x=r/r_\mathrm{s}`, :math:`X = \Delta x/r_\mathrm{s}`, and :math:`r_\mathrm{s}` is the NFW scale length.
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
    !!{RST
    Returns the scale-free enclosed mass at the given scale-free radius. The analytic solution (computed using Mathematica) is

    .. math::

       M(x) = 4 \pi \left[ -\frac{\sqrt{x^2+X^2}}{(1+x) \left(1+X^2\right)}+\tanh ^{-1}\left(\frac{x}{\sqrt{x^2+X^2}}\right)+\frac{\left(1+2X^2\right) \tanh ^{-1}\left(\frac{X^2-x}{\sqrt{1+X^2} \sqrt{x^2+X^2}}\right)}{\left(1+X^2\right)^{3/2}} -\frac{\left(1 + 2 X^2\right) \tanh ^{-1}\left(\sqrt{\frac{X^2}{1 + X^2}}\right)}{\left(1+ X^2\right)^{3/2}}+\frac{\sqrt{X^2}}{1 + X^2} \right],

    where :math:`x=r/r_\mathrm{s}`, :math:`X = \Delta x/r_\mathrm{s}`, and :math:`r_\mathrm{s}` is the NFW scale length.
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
    !!{RST
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self

    isAnalytic=.true.
    return
  end function sphericalFiniteResolutionNFWPotentialIsAnalytic

  double precision function sphericalFiniteResolutionNFWPotential(self,coordinates,status) result(potential)
    !!{RST
    Returns the potential (in (km/s)\ :math:`^2`) in the dark matter profile of ``node`` at the given ``radius`` (given in units of Mpc). The analytic solution (computed using Mathematica) is

    .. math::

       \Phi(x) &=& -\frac{\mathrm{G} M}{r_\mathrm{s}}  \nonumber \\
       & & \left\{ +\frac{\sqrt{x^2+X^2}}{x \left(X^2+1\right)} \right. \nonumber \\
       & & -\frac{X^2 \log \left(\sqrt{X^2+1} \sqrt{x^2+X^2}-x+X^2\right)}{\left(X^2+1\right)^{3/2}} \nonumber \\
       & & -\frac{\tanh ^{-1}\left(\frac{x}{\sqrt{x^2+X^2}}\right)}{x} \nonumber \\
       & & -\frac{\left(2 X^2+1\right) \tanh ^{-1}\left(\frac{X^2-x}{\sqrt{X^2+1} \sqrt{x^2+X^2}}\right)}{x \left(X^2+1\right)^{3/2}} \nonumber \\
       & & -\frac{\sqrt{X^2}}{x \left(X^2+1\right)}+\frac{X^2 \log (x+1)}{\left(X^2+1\right)^{3/2}} \nonumber \\
       & & +\frac{\left(2 X^2+1\right) \tanh ^{-1}\left(\sqrt{\frac{X^2}{X^2+1}}\right)}{x \left(X^2+1\right)^{3/2}} \nonumber \\
       & & \left. +\frac{ \left(\sqrt{X^2+1}-X^2 \log \left(\sqrt{X^2+1}-1\right)\right)}{\left(X^2+1\right)^{3/2}} \right\} \nonumber \\
       & & /\left[\log (1+c)-\frac{c}{1+c}\right]
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
    !!{RST
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
       call radiusEnclosingMassTableLengthResolutionInterpolator%linearFactors(self%lengthResolutionScaleFree,jLengthResolution(0),hLengthResolution)
       jLengthResolution(1)=jLengthResolution(0)+1
       self%radiusEnclosingMassPrevious=0.0d0
       do iLengthResolution=0,1
          self%radiusEnclosingMassPrevious=+self%radiusEnclosingMassPrevious                                                                                                     &
               &                           +radiusEnclosingMassTableMassInterpolator%interpolate(massScaleFree,radiusEnclosingMassTable(:,jLengthResolution(iLengthResolution))) &
               &                           *                                                                                              hLengthResolution(iLengthResolution)
       end do
       self%radiusEnclosingMassPrevious=+self%radiusEnclosingMassPrevious &
            &                           *self%radiusScale
    end if
    radius=self%radiusEnclosingMassPrevious
    return
  end function sphericalFiniteResolutionNFWRadiusEnclosingMass
  
  subroutine sphericalFiniteResolutionNFWRadiusEnclosingMassTabulate(self,mass,lengthResolution)
    !!{RST
    Tabulates the radius enclosing a given mass for finite resolution NFW mass profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Range_Pinned             , Range_Lattice_Extend_2D      , gridSchemePerDecade
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target       :: self
    double precision                                              , intent(in   )               :: mass                   , lengthResolution
    double precision                                              , parameter                   :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-9
    logical                                                                                     :: retabulate
    integer                                                                                     :: iLengthResolution      , iMass                   , &
         &                                                                                         i
    type            (rootFinder                                  )                              :: finder
    type            (rangeLattice                                )                              :: latticeMass            , latticeLengthResolution
    logical                                                       , allocatable, dimension(:,:) :: isComputed

    do i=1,2
       retabulate=.false.
       if (.not.radiusEnclosingMassTableInitialized) then
          retabulate=.true.
       else if (                                                               &
            &    mass             < radiusEnclosingMassMassMinimum             &
            &   .or.                                                           &
            &    mass             > radiusEnclosingMassMassMaximum             &
            &   .or.                                                           &
            &    lengthResolution < radiusEnclosingMassLengthResolutionMinimum &
            &   .or.                                                           &
            &    lengthResolution > radiusEnclosingMassLengthResolutionMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (retabulate     .and.i==1) call self%restoreMassTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Find the range to tabulate, pinning each axis to an absolute lattice so that the points evaluated - and therefore every
       ! value interpolated between them - depend only on which lattice points are spanned, and not on the sequence of values
       ! which happened to be requested. Each request is passed as the target and the range already tabulated is unioned in
       ! through `latticeCurrent`; folding the latter into the target instead - as the `min`/`max` against the current bounds
       ! formerly did - would apply the factor-of-two margin to an already-margined bound and so ratchet the range outward on
       ! every retabulation.
       latticeMass            =Range_Pinned(                                                                          &
            &                                              [mass                                                   ], &
            &                                               radiusEnclosingMassTableMassPointsPerDecade             , &
            &                                               gridSchemePerDecade                                     , &
            &                               marginFactor  = 2.0d0                                                   , &
            &                               anchorEvery   =+radiusEnclosingMassTableMassPointsPerDecade               &
            &                                              /tabulationAnchorsPerDecade                              , &
            &                               latticeCurrent= radiusEnclosingMassLatticeMass                            &
            &                              )
       latticeLengthResolution=Range_Pinned(                                                                          &
            &                                              [lengthResolution                                       ], &
            &                                               radiusEnclosingMassTableLengthResolutionPointsPerDecade , &
            &                                               gridSchemePerDecade                                     , &
            &                               marginFactor  = 2.0d0                                                   , &
            &                               anchorEvery   =+radiusEnclosingMassTableLengthResolutionPointsPerDecade   &
            &                                              /tabulationAnchorsPerDecade                              , &
            &                               latticeCurrent= radiusEnclosingMassLatticeLengthResolution                &
            &                              )
       ! Extend the table onto the new lattices, carrying over the solutions already found. Every offset is computed in exact
       ! integer arithmetic from the lattice indices, so no abscissa is compared.
       call Range_Lattice_Extend_2D(radiusEnclosingMassLatticeMass,latticeMass,radiusEnclosingMassLatticeLengthResolution,latticeLengthResolution,radiusEnclosingMassTable,isComputed)
       radiusEnclosingMassLatticeMass               =latticeMass
       radiusEnclosingMassLatticeLengthResolution   =latticeLengthResolution
       radiusEnclosingMassTableMassCount            =latticeMass            %count
       radiusEnclosingMassTableLengthResolutionCount=latticeLengthResolution%count
       radiusEnclosingMassMassMinimum               =latticeMass            %minimum()
       radiusEnclosingMassMassMaximum               =latticeMass            %maximum()
       radiusEnclosingMassLengthResolutionMinimum   =latticeLengthResolution%minimum()
       radiusEnclosingMassLengthResolutionMaximum   =latticeLengthResolution%maximum()
       ! Take the abscissae from the lattices. They must come from there, and never from a range laid out across the current
       ! extent: the lattice evaluates them through a single, deliberately un-inlined path, so that a given lattice point is
       ! bit-identical between one tabulation and another regardless of how many points each spans.
       if (allocated(radiusEnclosingMassTableMass            )) deallocate(radiusEnclosingMassTableMass            )
       if (allocated(radiusEnclosingMassTableLengthResolution)) deallocate(radiusEnclosingMassTableLengthResolution)
       radiusEnclosingMassTableMass            =latticeMass            %values()
       radiusEnclosingMassTableLengthResolution=latticeLengthResolution%values()
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
       do iLengthResolution=1,radiusEnclosingMassTableLengthResolutionCount
          iLengthResolution_=iLengthResolution
          do iMass=1,radiusEnclosingMassTableMassCount
             ! Skip the solutions carried over from the tabulation already in hand - finding them again would merely reproduce
             ! them, at the cost of a root find apiece.
             if (isComputed(iMass,iLengthResolution)) cycle
             iMass_=iMass
             ! Check that the root condition is satisfied at infinitely large radius. If it is not, then no radius encloses the
             ! required mass. Simply set the radius to an infinitely large value in such case.
             if (rootMass(radius=huge(0.0d0)) < 0.0d0) then
                radiusEnclosingMassTable(iMass,iLengthResolution)=huge(0.0d0)
             else
                radiusEnclosingMassTable(iMass,iLengthResolution)=finder%find(rootGuess=1.0d0)
             end if
          end do
       end do
       ! Build interpolators.
       if (allocated(radiusEnclosingMassTableLengthResolutionInterpolator)) deallocate(radiusEnclosingMassTableLengthResolutionInterpolator)
       if (allocated(radiusEnclosingMassTableMassInterpolator            )) deallocate(radiusEnclosingMassTableMassInterpolator            )
       allocate(radiusEnclosingMassTableLengthResolutionInterpolator)
       allocate(radiusEnclosingMassTableMassInterpolator            )
       radiusEnclosingMassTableLengthResolutionInterpolator=interpolator(radiusEnclosingMassTableLengthResolution)
       radiusEnclosingMassTableMassInterpolator            =interpolator(radiusEnclosingMassTableMass            )
       ! Specify that tabulation has been made.
       radiusEnclosingMassTableInitialized=.true.
       call self%storeMassTable()
    end if
    return
  end subroutine sphericalFiniteResolutionNFWRadiusEnclosingMassTabulate

  double precision function rootMass(radius)
    !!{RST
    Root function used in finding the radius enclosing a given mean mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    rootMass=+self_%massEnclosedScaleFree       (radius,radiusEnclosingMassTableLengthResolution(iLengthResolution_)) &
         &   -      radiusEnclosingMassTableMass(                                                iMass_             )
    return
  end function rootMass

  subroutine sphericalFiniteResolutionNFWStoreMassTable(self)
    !!{RST
    Store the tabulated radius-enclosing-mass data to file.
    !!}
    use :: File_Utilities    , only : File_Lock                , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: Table_Caches      , only : Table_Cache_Lattice_Write
    use :: Input_Paths       , only : inputPath                , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string           , operator(//)       , char
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5File                                    )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)// &
         &   'darkMatter/'                 // &
         &   self%objectType(             )// &
         &   '_mass_'                      // &
         &   self%suffix    (             )// &
         &   '.hdf5'
    call Directory_Make(File_Path(fileName))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(fileName,fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    file=hdf5File(fileName,overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    ! Record the lattices on which the two axes are built, so that a restored tabulation is recognized as lying on the
    ! same absolute lattice as one built here.
    call Table_Cache_Lattice_Write(file,'mass',radiusEnclosingMassLatticeMass)
    call Table_Cache_Lattice_Write(file,'lengthResolution',radiusEnclosingMassLatticeLengthResolution)
    call file%writeDataset(radiusEnclosingMassTableLengthResolution,'lengthResolution')
    call file%writeDataset(radiusEnclosingMassTableMass            ,'mass'            )
    call file%writeDataset(radiusEnclosingMassTable                ,'radius'          )
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine sphericalFiniteResolutionNFWStoreMassTable

  subroutine sphericalFiniteResolutionNFWRestoreMassTable(self)
    !!{RST
    Restore the tabulated radius-enclosing-mass data from file.

    The stored tabulation is adopted only if the file records, for both axes, a lattice which is self-consistent and which uses
    the density of points that this object would use, and if the array stored alongside them has the extent they imply. It is
    also declined if it does not *contain* the tabulation already in hand: this is called only when the latter has been found
    insufficient, so adopting a narrower tabulation in its place would discard solutions which must then be found again.
    !!}
    use :: File_Utilities    , only : File_Exists             , File_Lock          , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: Input_Paths       , only : inputPath               , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string          , operator(//)       , char
    use :: Numerical_Ranges  , only : gridSchemePerDecade
    use :: Table_Caches      , only : Table_Cache_Lattice_Read
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout)               :: self
    type            (lockDescriptor                              )                              :: fileLock
    type            (hdf5File                                    )                              :: file
    type            (varying_string                              )                              :: fileName
    type            (rangeLattice                                )                              :: latticeMass, latticeLengthResolution
    double precision                                              , allocatable, dimension(:,:) :: tableStored
    logical                                                                                     :: isUsable

    fileName=inputPath(pathTypeDataDynamic)// &
         &   'darkMatter/'                 // &
         &   self%objectType(             )// &
         &   '_mass_'                      // &
         &   self%suffix    (             )// &
         &   '.hdf5'
    if (.not.File_Exists(fileName)) return
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads. Note that the file is
    ! opened read-only: opening it for writing would take an exclusive lock and abort any concurrent reader.
    call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
    !$ call hdf5Access%set()
    file=hdf5File(fileName,readOnly=.true.)
    call Table_Cache_Lattice_Read(file,'mass',gridSchemePerDecade,radiusEnclosingMassTableMassPointsPerDecade,latticeMass)
    call Table_Cache_Lattice_Read(file,'lengthResolution',gridSchemePerDecade,radiusEnclosingMassTableLengthResolutionPointsPerDecade,latticeLengthResolution)
    isUsable=latticeMass%isDefined() .and. latticeLengthResolution%isDefined()
    if (isUsable) call file%readDataset('radius',tableStored)
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    if (isUsable)                                                                 &
         & isUsable=      size(tableStored,dim=1) == latticeMass%count            &
         &          .and. size(tableStored,dim=2) == latticeLengthResolution%count
    if (isUsable .and. radiusEnclosingMassTableInitialized)                                                                                                             &
         & isUsable=      (.not.radiusEnclosingMassLatticeMass            %isDefined() .or. latticeMass            %covers(radiusEnclosingMassLatticeMass            )) &
         &          .and. (.not.radiusEnclosingMassLatticeLengthResolution%isDefined() .or. latticeLengthResolution%covers(radiusEnclosingMassLatticeLengthResolution))
    if (.not.isUsable) return
    ! Adopt the stored tabulation, recovering everything which describes its extent from the lattices rather than from the
    ! stored abscissae, so that a restored tabulation cannot come to be described differently from a freshly built one.
    if (allocated(radiusEnclosingMassTableMass            )) deallocate(radiusEnclosingMassTableMass            )
    if (allocated(radiusEnclosingMassTableLengthResolution)) deallocate(radiusEnclosingMassTableLengthResolution)
    if (allocated(radiusEnclosingMassTable                )) deallocate(radiusEnclosingMassTable                )
    call Move_Alloc(tableStored,radiusEnclosingMassTable)
    radiusEnclosingMassLatticeMass               =latticeMass
    radiusEnclosingMassLatticeLengthResolution   =latticeLengthResolution
    radiusEnclosingMassTableMass                 =latticeMass            %values ()
    radiusEnclosingMassTableLengthResolution     =latticeLengthResolution%values ()
    radiusEnclosingMassTableMassCount            =latticeMass            %count
    radiusEnclosingMassTableLengthResolutionCount=latticeLengthResolution%count
    radiusEnclosingMassMassMinimum               =latticeMass            %minimum()
    radiusEnclosingMassMassMaximum               =latticeMass            %maximum()
    radiusEnclosingMassLengthResolutionMinimum   =latticeLengthResolution%minimum()
    radiusEnclosingMassLengthResolutionMaximum   =latticeLengthResolution%maximum()
    if (allocated(radiusEnclosingMassTableLengthResolutionInterpolator)) deallocate(radiusEnclosingMassTableLengthResolutionInterpolator)
    if (allocated(radiusEnclosingMassTableMassInterpolator            )) deallocate(radiusEnclosingMassTableMassInterpolator            )
    allocate(radiusEnclosingMassTableLengthResolutionInterpolator)
    allocate(radiusEnclosingMassTableMassInterpolator            )
    radiusEnclosingMassTableLengthResolutionInterpolator=interpolator(radiusEnclosingMassTableLengthResolution)
    radiusEnclosingMassTableMassInterpolator            =interpolator(radiusEnclosingMassTableMass            )
    radiusEnclosingMassTableInitialized                 =.true.
    return
  end subroutine sphericalFiniteResolutionNFWRestoreMassTable

  double precision function sphericalFiniteResolutionNFWRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{RST
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
          call radiusEnclosingDensityTableLengthResolutionInterpolator%linearFactors(self%lengthResolutionScaleFree,jLengthResolution(0),hLengthResolution)
          jLengthResolution(1)=jLengthResolution(0)+1
          self%radiusEnclosingDensityPrevious=0.0d0
          do iLengthResolution=0,1
             self%radiusEnclosingDensityPrevious=+self%radiusEnclosingDensityPrevious                                                                                                              &
                  &                              +radiusEnclosingDensityTableDensityInterpolator%interpolate(densityScaleFree,radiusEnclosingDensityTable(:,jLengthResolution(iLengthResolution))) &
                  &                              *                                                                                                          hLengthResolution(iLengthResolution)
         end do
          self%radiusEnclosingDensityPrevious=+self%radiusEnclosingDensityPrevious &
               &                              *self%radiusScale
       end if
    end if
    radius=self%radiusEnclosingDensityPrevious
    return
  end function sphericalFiniteResolutionNFWRadiusEnclosingDensity
  
  subroutine sphericalFiniteResolutionNFWRadiusEnclosingDensityTabulate(self,density,lengthResolution)
    !!{RST
    Tabulates the radius enclosing a given density for finite resolution NFW density profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Range_Pinned             , Range_Lattice_Extend_2D      , gridSchemePerDecade
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target       :: self
    double precision                                              , intent(in   )               :: density                , lengthResolution
    double precision                                              , parameter                   :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-9
    logical                                                                                     :: retabulate
    integer                                                                                     :: iLengthResolution      , iDensity                , &
         &                                                                                         i
    type            (rootFinder                                  )                              :: finder
    type            (rangeLattice                                )                              :: latticeDensity         , latticeLengthResolution
    logical                                                       , allocatable, dimension(:,:) :: isComputed

    do i=1,2
       retabulate=.false.
       if (.not.radiusEnclosingDensityTableInitialized) then
          retabulate=.true.
       else if (                                                                  &
            &    density          < radiusEnclosingDensityDensityMinimum          &
            &   .or.                                                              &
            &    density          > radiusEnclosingDensityDensityMaximum          &
            &   .or.                                                              &
            &    lengthResolution < radiusEnclosingDensityLengthResolutionMinimum &
            &   .or.                                                              &
            &    lengthResolution > radiusEnclosingDensityLengthResolutionMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (retabulate     .and.i==1) call self%restoreDensityTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Find the range to tabulate, pinning each axis to an absolute lattice - see the corresponding comment in
       ! `sphericalFiniteResolutionNFWRadiusEnclosingMassTabulate` for why the request, and not the range already tabulated, is
       ! passed as the target.
       latticeDensity         =Range_Pinned(                                                                             &
            &                                              [density                                                   ], &
            &                                               radiusEnclosingDensityTableDensityPointsPerDecade          , &
            &                                               gridSchemePerDecade                                        , &
            &                               marginFactor  = 2.0d0                                                      , &
            &                               anchorEvery   =+radiusEnclosingDensityTableDensityPointsPerDecade            &
            &                                              /tabulationAnchorsPerDecade                                 , &
            &                               latticeCurrent= radiusEnclosingDensityLatticeDensity                         &
            &                              )
       latticeLengthResolution=Range_Pinned(                                                                             &
            &                                              [lengthResolution                                          ], &
            &                                               radiusEnclosingDensityTableLengthResolutionPointsPerDecade , &
            &                                               gridSchemePerDecade                                        , &
            &                               marginFactor  = 2.0d0                                                      , &
            &                               anchorEvery   =+radiusEnclosingDensityTableLengthResolutionPointsPerDecade   &
            &                                              /tabulationAnchorsPerDecade                                 , &
            &                               latticeCurrent= radiusEnclosingDensityLatticeLengthResolution                &
            &                              )
       ! Extend the table onto the new lattices, carrying over the solutions already found.
       call Range_Lattice_Extend_2D(radiusEnclosingDensityLatticeDensity,latticeDensity,radiusEnclosingDensityLatticeLengthResolution,latticeLengthResolution,radiusEnclosingDensityTable,isComputed)
       radiusEnclosingDensityLatticeDensity            =latticeDensity
       radiusEnclosingDensityLatticeLengthResolution   =latticeLengthResolution
       radiusEnclosingDensityTableDensityCount         =latticeDensity         %count
       radiusEnclosingDensityTableLengthResolutionCount=latticeLengthResolution%count
       radiusEnclosingDensityDensityMinimum            =latticeDensity         %minimum()
       radiusEnclosingDensityDensityMaximum            =latticeDensity         %maximum()
       radiusEnclosingDensityLengthResolutionMinimum   =latticeLengthResolution%minimum()
       radiusEnclosingDensityLengthResolutionMaximum   =latticeLengthResolution%maximum()
       ! Take the abscissae from the lattices.
       if (allocated(radiusEnclosingDensityTableDensity         )) deallocate(radiusEnclosingDensityTableDensity         )
       if (allocated(radiusEnclosingDensityTableLengthResolution)) deallocate(radiusEnclosingDensityTableLengthResolution)
       radiusEnclosingDensityTableDensity         =latticeDensity         %values()
       radiusEnclosingDensityTableLengthResolution=latticeLengthResolution%values()
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
       do iLengthResolution=1,radiusEnclosingDensityTableLengthResolutionCount
          iLengthResolution_=iLengthResolution
          do iDensity=1,radiusEnclosingDensityTableDensityCount
             ! Skip the solutions carried over from the tabulation already in hand.
             if (isComputed(iDensity,iLengthResolution)) cycle
             iDensity_=iDensity
             if (radiusEnclosingDensityTableDensity(iDensity) > 1.0d0/radiusEnclosingDensityTableLengthResolution(iLengthResolution)) then
                ! Density exceeds the maximum density in the profile - so set zero radius.
                radiusEnclosingDensityTable(iDensity,iLengthResolution)=0.0d0
             else
                radiusEnclosingDensityTable(iDensity,iLengthResolution)=finder%find(rootGuess=1.0d0)
             end if
          end do
       end do
       ! Build interpolators.
       if (allocated(radiusEnclosingDensityTableLengthResolutionInterpolator)) deallocate(radiusEnclosingDensityTableLengthResolutionInterpolator)
       if (allocated(radiusEnclosingDensityTableDensityInterpolator         )) deallocate(radiusEnclosingDensityTableDensityInterpolator         )
       allocate(radiusEnclosingDensityTableLengthResolutionInterpolator)
       allocate(radiusEnclosingDensityTableDensityInterpolator         )
       radiusEnclosingDensityTableLengthResolutionInterpolator=interpolator(radiusEnclosingDensityTableLengthResolution)
       radiusEnclosingDensityTableDensityInterpolator         =interpolator(radiusEnclosingDensityTableDensity         )
       ! Specify that tabulation has been made.
       radiusEnclosingDensityTableInitialized=.true.
       call self%storeDensityTable()
    end if
    return
  end subroutine sphericalFiniteResolutionNFWRadiusEnclosingDensityTabulate

  double precision function rootDensity(radius)
    !!{RST
    Root function used in finding the radius enclosing a given mean density.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    rootDensity=+3.0d0                                                                                                               &
         &      *self_%massEnclosedScaleFree             (radius,radiusEnclosingDensityTableLengthResolution(iLengthResolution_))    &
         &      /4.0d0                                                                                                               &
         &      /Pi                                                                                                                  &
         &      /                                         radius                                                                 **3 &
         &      -      radiusEnclosingDensityTableDensity(                                                   iDensity_          )
    return
  end function rootDensity

  subroutine sphericalFiniteResolutionNFWStoreDensityTable(self)
    !!{RST
    Store the tabulated radius-enclosing-density data to file.
    !!}
    use :: File_Utilities    , only : File_Lock                , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: Table_Caches      , only : Table_Cache_Lattice_Write
    use :: Input_Paths       , only : inputPath                , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string           , operator(//)       , char
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5File                                    )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)// &
         &   'darkMatter/'                 // &
         &   self%objectType(             )// &
         &   '_density_'                   // &
         &   self%suffix    (             )// &
         &   '.hdf5'
    call Directory_Make(File_Path(fileName))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(fileName,fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    file=hdf5File(fileName,overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    ! Record the lattices on which the two axes are built, so that a restored tabulation is recognized as lying on the
    ! same absolute lattice as one built here.
    call Table_Cache_Lattice_Write(file,'density'         ,radiusEnclosingDensityLatticeDensity         )
    call Table_Cache_Lattice_Write(file,'lengthResolution',radiusEnclosingDensityLatticeLengthResolution)
    call file%writeDataset(radiusEnclosingDensityTableLengthResolution,'lengthResolution')
    call file%writeDataset(radiusEnclosingDensityTableDensity         ,'density'         )
    call file%writeDataset(radiusEnclosingDensityTable                ,'radius'          )
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine sphericalFiniteResolutionNFWStoreDensityTable

  subroutine sphericalFiniteResolutionNFWRestoreDensityTable(self)
    !!{RST
    Restore the tabulated radius-enclosing-density data from file.

    The stored tabulation is adopted only if the file records, for both axes, a lattice which is self-consistent and which uses
    the density of points that this object would use, and if the array stored alongside them has the extent they imply. It is
    also declined if it does not *contain* the tabulation already in hand: this is called only when the latter has been found
    insufficient, so adopting a narrower tabulation in its place would discard solutions which must then be found again.
    !!}
    use :: File_Utilities    , only : File_Exists             , File_Lock          , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: Input_Paths       , only : inputPath               , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string          , operator(//)       , char
    use :: Numerical_Ranges  , only : gridSchemePerDecade
    use :: Table_Caches      , only : Table_Cache_Lattice_Read
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout)               :: self
    type            (lockDescriptor                              )                              :: fileLock
    type            (hdf5File                                    )                              :: file
    type            (varying_string                              )                              :: fileName
    type            (rangeLattice                                )                              :: latticeDensity, latticeLengthResolution
    double precision                                              , allocatable, dimension(:,:) :: tableStored
    logical                                                                                     :: isUsable

    fileName=inputPath(pathTypeDataDynamic)// &
         &   'darkMatter/'                 // &
         &   self%objectType(             )// &
         &   '_density_'                   // &
         &   self%suffix    (             )// &
         &   '.hdf5'
    if (.not.File_Exists(fileName)) return
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads. Note that the file is
    ! opened read-only: opening it for writing would take an exclusive lock and abort any concurrent reader.
    call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
    !$ call hdf5Access%set()
    file=hdf5File(fileName,readOnly=.true.)
    call Table_Cache_Lattice_Read(file,'density'         ,gridSchemePerDecade,radiusEnclosingDensityTableDensityPointsPerDecade         ,latticeDensity         )
    call Table_Cache_Lattice_Read(file,'lengthResolution',gridSchemePerDecade,radiusEnclosingDensityTableLengthResolutionPointsPerDecade,latticeLengthResolution)
    isUsable=latticeDensity%isDefined() .and. latticeLengthResolution%isDefined()
    if (isUsable) call file%readDataset('radius',tableStored)
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    if (isUsable)                                                                 &
         & isUsable=      size(tableStored,dim=1) == latticeDensity%count         &
         &          .and. size(tableStored,dim=2) == latticeLengthResolution%count
    if (isUsable .and. radiusEnclosingDensityTableInitialized)                                                                                                                &
         & isUsable=      (.not.radiusEnclosingDensityLatticeDensity         %isDefined() .or. latticeDensity         %covers(radiusEnclosingDensityLatticeDensity         )) &
         &          .and. (.not.radiusEnclosingDensityLatticeLengthResolution%isDefined() .or. latticeLengthResolution%covers(radiusEnclosingDensityLatticeLengthResolution))
    if (.not.isUsable) return
    ! Adopt the stored tabulation, recovering everything which describes its extent from the lattices rather than from the
    ! stored abscissae, so that a restored tabulation cannot come to be described differently from a freshly built one.
    if (allocated(radiusEnclosingDensityTableDensity         )) deallocate(radiusEnclosingDensityTableDensity         )
    if (allocated(radiusEnclosingDensityTableLengthResolution)) deallocate(radiusEnclosingDensityTableLengthResolution)
    if (allocated(radiusEnclosingDensityTable                )) deallocate(radiusEnclosingDensityTable                )
    call Move_Alloc(tableStored,radiusEnclosingDensityTable)
    radiusEnclosingDensityLatticeDensity            =latticeDensity
    radiusEnclosingDensityLatticeLengthResolution   =latticeLengthResolution
    radiusEnclosingDensityTableDensity              =latticeDensity         %values ()
    radiusEnclosingDensityTableLengthResolution     =latticeLengthResolution%values ()
    radiusEnclosingDensityTableDensityCount         =latticeDensity         %count
    radiusEnclosingDensityTableLengthResolutionCount=latticeLengthResolution%count
    radiusEnclosingDensityDensityMinimum            =latticeDensity         %minimum()
    radiusEnclosingDensityDensityMaximum            =latticeDensity         %maximum()
    radiusEnclosingDensityLengthResolutionMinimum   =latticeLengthResolution%minimum()
    radiusEnclosingDensityLengthResolutionMaximum   =latticeLengthResolution%maximum()
    if (allocated(radiusEnclosingDensityTableLengthResolutionInterpolator)) deallocate(radiusEnclosingDensityTableLengthResolutionInterpolator)
    if (allocated(radiusEnclosingDensityTableDensityInterpolator         )) deallocate(radiusEnclosingDensityTableDensityInterpolator         )
    allocate(radiusEnclosingDensityTableLengthResolutionInterpolator)
    allocate(radiusEnclosingDensityTableDensityInterpolator         )
    radiusEnclosingDensityTableLengthResolutionInterpolator=interpolator(radiusEnclosingDensityTableLengthResolution)
    radiusEnclosingDensityTableDensityInterpolator         =interpolator(radiusEnclosingDensityTableDensity         )
    radiusEnclosingDensityTableInitialized                 =.true.
    return
  end subroutine sphericalFiniteResolutionNFWRestoreDensityTable

  double precision function sphericalFiniteResolutionNFWEnergy(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{RST
    Compute the energy within a given ``radius`` in a finite-resolution NFW mass distribution.
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
       call energyTableLengthResolutionInterpolator%linearFactors(self%lengthResolutionScaleFree,jLengthResolution(0),hLengthResolution)
       jLengthResolution(1)=jLengthResolution(0)+1
       self%energyPrevious=0.0d0
       do iLengthResolution=0,1
          self%energyPrevious=+self%energyPrevious                                                                                                              &
               &              +energyTableRadiusOuterInterpolator%interpolate(radiusOuter/self%radiusScale,energyTable(:,jLengthResolution(iLengthResolution))) &
               &              *                                                                                          hLengthResolution(iLengthResolution)
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
    !!{RST
    Tabulates the energy for finite resolution NFW mass profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Range_Pinned, Range_Lattice_Extend_2D, gridSchemePerDecade
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout), target       :: self
    double precision                                              , intent(in   )               :: radiusOuter                , lengthResolution
    double precision                                              , parameter                   :: multiplierRadius   =100.0d0
    type            (integrator                                  )                              :: integratorPotential        , integratorKinetic      , &
         &                                                                                         integratorPressure
    double precision                                                                            :: pseudoPressure             , energyKinetic          , &
         &                                                                                         energyPotential            , radiusOuter_
    logical                                                                                     :: retabulate
    integer                                                                                     :: iLengthResolution          , iRadiusOuter           , &
         &                                                                                         i
    type            (rangeLattice                                )                              :: latticeRadiusOuter         , latticeLengthResolution
    logical                                                       , allocatable, dimension(:,:) :: isComputed

    do i=1,2
       retabulate=.false.
       if (.not.energyTableInitialized) then
          retabulate=.true.
       else if (                                                  &
            &    radiusOuter      < energyRadiusOuterMinimum      &
            &   .or.                                              &
            &    radiusOuter      > energyRadiusOuterMaximum      &
            &   .or.                                              &
            &    lengthResolution < energyLengthResolutionMinimum &
            &   .or.                                              &
            &    lengthResolution > energyLengthResolutionMaximum &
            &  ) then
          retabulate=.true.
       end if
       if (     retabulate.and.i==1) call self%restoreEnergyTable()
       if (.not.retabulate         ) exit
    end do
    if (retabulate) then
       ! Find the range to tabulate, pinning each axis to an absolute lattice - see the corresponding comment in
       ! `sphericalFiniteResolutionNFWRadiusEnclosingMassTabulate`.
       latticeRadiusOuter              =Range_Pinned(                                                             &
            &                                                       [radiusOuter                               ], &
            &                                                        energyTableRadiusOuterPointsPerDecade      , &
            &                                                        gridSchemePerDecade                        , &
            &                                        marginFactor  = 2.0d0                                      , &
            &                                        anchorEvery   = energyTableRadiusOuterPointsPerDecade        &
            &                                                       /tabulationAnchorsPerDecade                 , &
            &                                        latticeCurrent= energyLatticeRadiusOuter                     &
            &                                       )
       latticeLengthResolution         =Range_Pinned(                                                             &
            &                                                       [lengthResolution                          ], &
            &                                                        energyTableLengthResolutionPointsPerDecade , &
            &                                                        gridSchemePerDecade                        , &
            &                                        marginFactor  = 2.0d0                                      , &
            &                                        anchorEvery   = energyTableLengthResolutionPointsPerDecade   &
            &                                                       /tabulationAnchorsPerDecade                 , &
            &                                        latticeCurrent= energyLatticeLengthResolution                &
            &                                       )
       ! Extend the table onto the new lattices, carrying over the solutions already found.
       call Range_Lattice_Extend_2D(energyLatticeRadiusOuter,latticeRadiusOuter,energyLatticeLengthResolution,latticeLengthResolution,energyTable,isComputed)
       energyLatticeRadiusOuter        =latticeRadiusOuter
       energyLatticeLengthResolution   =latticeLengthResolution
       energyTableRadiusOuterCount     =latticeRadiusOuter     %count
       energyTableLengthResolutionCount=latticeLengthResolution%count
       energyRadiusOuterMinimum        =latticeRadiusOuter     %minimum()
       energyRadiusOuterMaximum        =latticeRadiusOuter     %maximum()
       energyLengthResolutionMinimum   =latticeLengthResolution%minimum()
       energyLengthResolutionMaximum   =latticeLengthResolution%maximum()
       ! Take the abscissae from the lattices.
       if (allocated(energyTableRadiusOuter     )) deallocate(energyTableRadiusOuter     )
       if (allocated(energyTableLengthResolution)) deallocate(energyTableLengthResolution)
       energyTableRadiusOuter     =latticeRadiusOuter     %values()
       energyTableLengthResolution=latticeLengthResolution%values()
       ! Initialize integrators.
       integratorPotential=integrator(integrandEnergyPotential,toleranceRelative=1.0d-3)
       integratorKinetic  =integrator(integrandEnergyKinetic  ,toleranceRelative=1.0d-3)
       integratorPressure =integrator(integrandPseudoPressure ,toleranceRelative=1.0d-3)
       ! Loop over radiusOuter and core radius and populate tables.
       self_ => self
       do iLengthResolution=1,energyTableLengthResolutionCount
          iLengthResolution_=iLengthResolution
          do iRadiusOuter=1,energyTableRadiusOuterCount
             ! Skip the solutions carried over from the tabulation already in hand - evaluating them again would merely
             ! reproduce them, at the cost of three numerical integrals apiece.
             if (isComputed(iRadiusOuter,iLengthResolution)) cycle
             radiusOuter_                                    =energyTableRadiusOuter(iRadiusOuter)
             energyPotential                                 =+integratorPotential%integrate(       0.0d0,                 radiusOuter_)
             energyKinetic                                   =+integratorKinetic  %integrate(       0.0d0,                 radiusOuter_)
             pseudoPressure                                  =+integratorPressure %integrate(radiusOuter_,multiplierRadius*radiusOuter_)
             energyTable(iRadiusOuter,iLengthResolution)=-0.5d0                                                                                             &
                  &                                           *(                                                                                            &
                  &                                             +energyPotential                                                                            &
                  &                                             +self%massEnclosedScaleFree(radiusOuter_,energyTableLengthResolution(iLengthResolution))**2 &
                  &                                             /radiusOuter_                                                                               &
                  &                                            )                                                                                            &
                  &                                           +2.0d0                                                                                        &
                  &                                           *Pi                                                                                           &
                  &                                           *(                                                                                            &
                  &                                             +radiusOuter_                                                                           **3 &
                  &                                             *pseudoPressure                                                                             &
                  &                                             +energyKinetic                                                                              &
                  &                                            )
            end do
       end do
       ! Build interpolators.
       if (allocated(energyTableLengthResolutionInterpolator)) deallocate(energyTableLengthResolutionInterpolator)
       if (allocated(energyTableRadiusOuterInterpolator     )) deallocate(energyTableRadiusOuterInterpolator     )
       allocate(energyTableLengthResolutionInterpolator)
       allocate(energyTableRadiusOuterInterpolator     )
       energyTableLengthResolutionInterpolator=interpolator(energyTableLengthResolution)
       energyTableRadiusOuterInterpolator     =interpolator(energyTableRadiusOuter     )
       ! Specify that tabulation has been made.
       energyTableInitialized=.true.
       call self%storeEnergyTable()
    end if
    return
  end subroutine sphericalFiniteResolutionNFWEnergyTabulate

  double precision function integrandEnergyPotential(radius)
    !!{RST
    Integrand for potential energy of the halo.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       integrandEnergyPotential=(                                                                                     &
            &                    +self_%massEnclosedScaleFree(radius,energyTableLengthResolution(iLengthResolution_)) &
            &                    /                            radius                                                  &
            &                   )**2
    else
       integrandEnergyPotential=+0.0d0
    end if
    return
  end function integrandEnergyPotential
  
  double precision function integrandEnergyKinetic(radius)
    !!{RST
    Integrand for kinetic energy of the halo.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       integrandEnergyKinetic=+self_%massEnclosedScaleFree(radius,energyTableLengthResolution(iLengthResolution_)) &
            &                 *self_%densityScaleFree     (radius,energyTableLengthResolution(iLengthResolution_)) &
            &                 *                            radius
    else
       integrandEnergyKinetic=+0.0d0
    end if
    return
  end function integrandEnergyKinetic
  
  double precision function integrandPseudoPressure(radius)
    !!{RST
    Integrand for pseudo-pressure (:math:`\rho(r) \sigma^2(r)`) of the halo.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    
    if (radius > 0.0d0) then
       integrandPseudoPressure=+self_%massEnclosedScaleFree(radius,energyTableLengthResolution(iLengthResolution_))    &
            &                  *self_%densityScaleFree     (radius,energyTableLengthResolution(iLengthResolution_))    &
            &                  /                            radius                                                 **2
    else
       integrandPseudoPressure=+0.0d0
    end if
    return
  end function integrandPseudoPressure

  double precision function sphericalFiniteResolutionNFWDensityScaleFree(self,radius,radiusCore) result(densityScaleFree)
    !!{RST
    Returns the scale-free density in the dark matter profile at the given ``radius``.
    !!}
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    double precision                                              , intent(in   ) :: radius, radiusCore
    !$GLC attributes unused :: self
    
    densityScaleFree=1.0d0/(1.0d0+radius)**2/sqrt(radius**2+radiusCore**2)
    return
  end function sphericalFiniteResolutionNFWDensityScaleFree
    
  subroutine sphericalFiniteResolutionNFWStoreEnergyTable(self)
    !!{RST
    Store the tabulated energy data to file.
    !!}
    use :: File_Utilities    , only : File_Lock                , File_Unlock        , lockDescriptor, Directory_Make, &
         &                            File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: Table_Caches      , only : Table_Cache_Lattice_Write
    use :: Input_Paths       , only : inputPath                , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string           , operator(//)       , char
    implicit none
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    type (lockDescriptor                              )                :: fileLock
    type (hdf5File                                    )                :: file
    type (varying_string                              )                :: fileName

    fileName=inputPath(pathTypeDataDynamic)// &
         &   'darkMatter/'                 // &
         &   self%objectType(             )// &
         &   '_energy_'                    // &
         &   self%suffix    (             )// &
         &   '.hdf5'
    call Directory_Make(File_Path(fileName))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(fileName,fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    file=hdf5File(fileName,overWrite=.true.,objectsOverwritable=.true.,readOnly=.false.)
    ! Record the lattices on which the two axes are built, so that a restored tabulation is recognized as lying on the
    ! same absolute lattice as one built here.
    call Table_Cache_Lattice_Write(file,'radiusOuter'     ,energyLatticeRadiusOuter     )
    call Table_Cache_Lattice_Write(file,'lengthResolution',energyLatticeLengthResolution)
    call file%writeDataset(energyTableLengthResolution,'lengthResolution')
    call file%writeDataset(energyTableRadiusOuter     ,'radiusOuter'     )
    call file%writeDataset(energyTable                ,'energy'          )
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine sphericalFiniteResolutionNFWStoreEnergyTable

  subroutine sphericalFiniteResolutionNFWRestoreEnergyTable(self)
    !!{RST
    Restore the tabulated energy data from file.

    The stored tabulation is adopted only if the file records, for both axes, a lattice which is self-consistent and which uses
    the density of points that this object would use, and if the array stored alongside them has the extent they imply. It is
    also declined if it does not *contain* the tabulation already in hand: this is called only when the latter has been found
    insufficient, so adopting a narrower tabulation in its place would discard solutions which must then be found again.
    !!}
    use :: File_Utilities    , only : File_Exists             , File_Lock          , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: Input_Paths       , only : inputPath               , pathTypeDataDynamic
    use :: ISO_Varying_String, only : varying_string          , operator(//)       , char
    use :: Numerical_Ranges  , only : gridSchemePerDecade
    use :: Table_Caches      , only : Table_Cache_Lattice_Read
    implicit none
    class           (massDistributionSphericalFiniteResolutionNFW), intent(inout)               :: self
    type            (lockDescriptor                              )                              :: fileLock
    type            (hdf5File                                    )                              :: file
    type            (varying_string                              )                              :: fileName
    type            (rangeLattice                                )                              :: latticeRadiusOuter, latticeLengthResolution
    double precision                                              , allocatable, dimension(:,:) :: tableStored
    logical                                                                                     :: isUsable

    fileName=inputPath(pathTypeDataDynamic)// &
         &   'darkMatter/'                 // &
         &   self%objectType(             )// &
         &   '_energy_'                    // &
         &   self%suffix    (             )// &
         &   '.hdf5'
    if (.not.File_Exists(fileName)) return
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads. Note that the file is
    ! opened read-only: opening it for writing would take an exclusive lock and abort any concurrent reader.
    call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
    !$ call hdf5Access%set()
    file=hdf5File(fileName,readOnly=.true.)
    call Table_Cache_Lattice_Read(file,'radiusOuter',gridSchemePerDecade,energyTableRadiusOuterPointsPerDecade,latticeRadiusOuter)
    call Table_Cache_Lattice_Read(file,'lengthResolution',gridSchemePerDecade,energyTableLengthResolutionPointsPerDecade,latticeLengthResolution)
    isUsable=latticeRadiusOuter%isDefined() .and. latticeLengthResolution%isDefined()
    if (isUsable) call file%readDataset('radius',tableStored)
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    if (isUsable)                                                                  &
         & isUsable=      size(tableStored,dim=1) == latticeRadiusOuter     %count &
         &          .and. size(tableStored,dim=2) == latticeLengthResolution%count
    if (isUsable .and. energyTableInitialized)                                                                                                &
         & isUsable=      (.not.energyLatticeRadiusOuter     %isDefined() .or. latticeRadiusOuter     %covers(energyLatticeRadiusOuter     )) &
         &          .and. (.not.energyLatticeLengthResolution%isDefined() .or. latticeLengthResolution%covers(energyLatticeLengthResolution))
    if (.not.isUsable) return
    ! Adopt the stored tabulation, recovering everything which describes its extent from the lattices rather than from the
    ! stored abscissae, so that a restored tabulation cannot come to be described differently from a freshly built one.
    if (allocated(energyTableRadiusOuter     )) deallocate(energyTableRadiusOuter     )
    if (allocated(energyTableLengthResolution)) deallocate(energyTableLengthResolution)
    if (allocated(energyTable                )) deallocate(energyTable                )
    call Move_Alloc(tableStored,energyTable)
    energyLatticeRadiusOuter        =latticeRadiusOuter
    energyLatticeLengthResolution   =latticeLengthResolution
    energyTableRadiusOuter          =latticeRadiusOuter     %values ()
    energyTableLengthResolution     =latticeLengthResolution%values ()
    energyTableRadiusOuterCount     =latticeRadiusOuter     %count
    energyTableLengthResolutionCount=latticeLengthResolution%count
    energyRadiusOuterMinimum        =latticeRadiusOuter     %minimum()
    energyRadiusOuterMaximum        =latticeRadiusOuter     %maximum()
    energyLengthResolutionMinimum   =latticeLengthResolution%minimum()
    energyLengthResolutionMaximum   =latticeLengthResolution%maximum()
    if (allocated(energyTableLengthResolutionInterpolator)) deallocate(energyTableLengthResolutionInterpolator)
    if (allocated(energyTableRadiusOuterInterpolator     )) deallocate(energyTableRadiusOuterInterpolator     )
    allocate(energyTableLengthResolutionInterpolator)
    allocate(energyTableRadiusOuterInterpolator     )
    energyTableLengthResolutionInterpolator=interpolator(energyTableLengthResolution)
    energyTableRadiusOuterInterpolator     =interpolator(energyTableRadiusOuter     )
    energyTableInitialized                 =.true.
    return
  end subroutine sphericalFiniteResolutionNFWRestoreEnergyTable

  function finiteResolutionNFWSuffix(self) result(suffix)
    !!{RST
    Return a suffix for tabulated file names.
    !!}
    use :: String_Handling, only : String_C_To_Fortran
    implicit none
    type (varying_string                              )                :: suffix
    class(massDistributionSphericalFiniteResolutionNFW), intent(inout) :: self
    !$GLC attributes unused :: self

    suffix=String_C_To_Fortran(massDistributionFiniteResolutionNFWSourceDigest)
    return
  end function finiteResolutionNFWSuffix
