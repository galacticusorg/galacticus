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
  Implementation of a mass distribution class for exponential disks which have been tidally truncated and heated.
  !!}

  !![
  <massDistribution name="massDistributionExponentialDiskTidal" docformat="rst">
   <description>
   A mass distribution class for a dimensionless exponential disk (see :galacticus-class:`massDistributionExponentialDisk`)
   which has been truncated by tidal stripping and heated by tidal shocks. Lengths are in units of the scale length of the
   disk, :math:`R_\mathrm{d}`.

   **Truncation.** The degree of truncation is specified by :math:`x=`\ ``[truncation]``, such that the disk retains a fraction
   :math:`f=\exp(-x)` of the mass of the untruncated disk. Beyond a truncation radius, :math:`R_\mathrm{t}`, the surface
   density of the exponential disk, :math:`\Sigma_\mathrm{e}(R)`, is multiplied by

   .. math::

      g(u) = u^{1/y_\mathrm{d}} \exp\left(-\frac{u-1}{y_\mathrm{d}}\right),

   where :math:`u=R/R_\mathrm{t}` and :math:`y_\mathrm{d}=`\ ``[radiusDecayOverRadiusTruncate]``. Since :math:`g(1)=1` and
   :math:`g^\prime(1)=0` the surface density and its gradient are continuous at :math:`R_\mathrm{t}`, while :math:`g(u)&lt;1`
   for all :math:`u&gt;1`, such that the truncated disk always has less mass beyond :math:`R_\mathrm{t}` than the untruncated
   disk, and :math:`x` increases monotonically as :math:`R_\mathrm{t}` decreases. (The spherical truncation of
   :galacticus-class:`massDistributionSphericalTruncatedExponential` instead joins a power law to the density at the truncation
   radius. For an exponential disk that power law becomes shallower than the disk itself for :math:`R_\mathrm{t} \gtrsim
   R_\mathrm{d}/y_\mathrm{d}`, so is not used here.) :math:`R_\mathrm{t}` is chosen such that the truncated disk retains a
   fraction :math:`f` of the mass of the untruncated disk.

   **Heating.** Tidal shocks impart a change in specific angular momentum to each star which, to first order in the impulse
   approximation (e.g. :cite:t:`carr_migration_2022`), averages to zero around every ring of the disk. The disk therefore
   retains its angular momentum, while the guiding radii of its stars are spread out. Assuming the self-gravity of the disk to
   be negligible, and the rotation curve to be flat, a ring at radius :math:`R` is spread in radius with a dispersion
   :math:`\sigma_R = s R^2/R_\mathrm{d}`, where :math:`s=`\ ``[heating]``. The ring is spread with a log-normal distribution
   in :math:`R`, of width :math:`\sigma_{\ln R} = s R/R_\mathrm{d}`, shifted such that its mean radius (and so, for a flat
   rotation curve, its mean specific angular momentum) is conserved. Heating is applied after truncation.

   **Gravitational field.** The mid-plane rotation curve is computed exactly, for the :math:`\mathrm{sech}^2` vertical structure
   of the exponential disk, as a Hankel transform of the modified surface density, evaluated using :term:`FFTLog`. It is tabulated
   as its ratio to the rotation curve of the untruncated, unheated disk, which is then multiplied by the rotation curve of
   :galacticus-class:`massDistributionExponentialDisk` (such that any errors in the transform cancel to first order). The
   mid-plane potential is computed exactly, as :math:`\Phi(R,0)=-\int_R^\infty V_\mathrm{c}^2(R^\prime) \mathrm{d}\ln R^\prime`,
   and is tabulated in the same way. Away from the mid-plane, the potential and acceleration are *approximated*. Close to the disk
   they are those of the untruncated, unheated disk, multiplied by the ratio of the mid-plane potentials, or of the squares of the
   mid-plane rotation curves, respectively, evaluated at the spherical radius of the point considered. This retains the flattening
   of the untruncated disk, which a truncated disk does not have. So, further out, they are blended smoothly into the multipole
   expansion of the potential to quadrupole order, with quadrupole moment :math:`Q=M(\langle z^2 \rangle - \langle R^2
   \rangle/2)`, normalized at each spherical radius such that it reproduces the mid-plane potential and rotation curve. The weight
   given to the multipole expansion rises from zero at 2, to unity at 5, times the root-mean-square radius, :math:`\sqrt{\langle
   R^2 \rangle + \langle z^2 \rangle}`, of the disk. Since the multipole expansion cannot describe the field close to material at
   large radii (such as the exponential tail of a weakly truncated disk), this weight is further multiplied by
   :math:`(1-\sigma)^2`, where :math:`\sigma` is the ratio of the surface density to that of the untruncated, unheated disk at the
   spherical radius of the point (or zero if :math:`\sigma>1`). The tidal tensor is the (symmetrized) gradient of the approximate
   acceleration field, with its vertical-vertical component chosen such that its trace satisfies Poisson's equation for the
   density of the disk. The potential, acceleration, and tidal tensor are therefore exact in the mid-plane, and approximate away
   from it. Compared with direct summation over the disk (for :math:`h=0.137 R_\mathrm{d}` and :math:`y_\mathrm{d}=0.1`), beyond
   five root-mean-square radii the potential away from the mid-plane is accurate to better than 0.5%, and the vertical
   acceleration to about 2% (and to better than 1% beyond six). Within a few root-mean-square radii of a truncated disk, but
   beyond its truncation radius, the errors can reach tens of percent in the potential, and factors of several in the vertical
   acceleration close to the mid-plane. Beyond the tabulated range of radii (see below) the multipole expansion of the potential
   to quadrupole order is used without normalization.  Quantities which depend on the mass distribution in three dimensions (such
   as the mass enclosed by a sphere) use the thin-disk approximation, as does :galacticus-class:`massDistributionExponentialDisk`.

   **Tabulation.** The rotation curve, the mid-plane potential, the mass enclosed within a cylinder, and the surface density are
   tabulated as their ratios to those of the untruncated, unheated disk at the same radius (ratios which, unlike the quantities
   themselves, vary smoothly with :math:`R_\mathrm{t}` and :math:`s` close to the center of the disk), as functions of
   :math:`R/R_\mathrm{t}` (for :math:`10^{-3} \le R/R_\mathrm{t} \lesssim 300`) on a grid of :math:`R_\mathrm{t}` and
   :math:`s`, and interpolated (using four-point Lagrange interpolation in :math:`\ln R_\mathrm{t}`,
   :math:`\sinh^{-1}(s/s_0)`, and :math:`\ln(R/R_\mathrm{t})`) thereafter. The tables are stored in the dynamic data path, and
   depend only on the dimensionless scale height of the disk and :math:`y_\mathrm{d}`. Truncations with :math:`x &lt; 10^{-6}`
   are treated as :math:`x=10^{-6}`. Truncations stronger than that of the smallest tabulated truncation radius,
   :math:`10^{-2}R_\mathrm{d}` (for which :math:`x \approx 9`), retain the shape of the disk truncated at that radius, with
   its mass reduced to :math:`\exp(-x)` of that of the untruncated disk. Heating parameters :math:`s&gt;0.3` are limited to
   that value. The width of the log-normal distribution used for heating is limited to :math:`\sigma_{\ln R} \le 1`: rings for
   which :math:`s R &gt; R_\mathrm{d}` are far beyond the regime in which the first-order treatment of heating applies, and are
   treated as fully mixed in radius.
   </description>
  </massDistribution>
  !!]

  type :: exponentialDiskTidalTable
     !!{RST
     Tables of the properties of tidally truncated and heated exponential disks of a given dimensionless scale height and
     truncation decay length.
     !!}
     double precision                                                                       :: scaleHeight                               , radiusDecayOverRadiusTruncate
     integer                                                                                :: countTruncation                           , countHeating                     , &
          &                                                                                    countRadius
     double precision                                                                       :: radiusTruncateLogarithmicMinimum          , radiusTruncateLogarithmicInterval, &
          &                                                                                    heatingInterval                           , radiusLogarithmicMinimum         , &
          &                                                                                    radiusLogarithmicInterval
     logical                                                                                :: monotonic
     double precision                           , allocatable, dimension(:    )             :: truncationLogarithmic
     double precision                           , allocatable, dimension(:,:  )             :: moment1                                   , moment2                          , &
          &                                                                                    moment3
     double precision                           , allocatable, dimension(:,:,:)             :: rotationCurveRatioLogarithmic             , massRatioLogarithmic             , &
          &                                                                                    surfaceDensityRatioLogarithmic            , potentialRatioLogarithmic
     type            (exponentialDiskTidalTable), pointer                                   :: next                             => null()
  end type exponentialDiskTidalTable

  type, public, extends(massDistributionCylindrical) :: massDistributionExponentialDiskTidal
     !!{RST
     A mass distribution class for tidally truncated and heated exponential disks.
     !!}
     private
     class           (massDistributionExponentialDisk), pointer         :: massDistribution_             => null()
     type            (exponentialDiskTidalTable      ), pointer         :: table                         => null()
     double precision                                                   :: truncation                             , heating                  , &
          &                                                                radiusDecayOverRadiusTruncate          , radiusTruncateLogarithmic, &
          &                                                                massTotal_                             , radiusHalfMass_          , &
          &                                                                massScaleLogarithmic                   , quadrupole               , &
          &                                                                radiusBlendMinimum                     , radiusBlendMaximum
     integer                                          , dimension(4  ) :: indexTruncation                         , indexHeating
     double precision                                 , dimension(4,4) :: weight
   contains
     !![
     <methods docformat="rst">
       <method method="interpolate"        description="Interpolate a tabulated quantity at the given dimensionless radius, optionally returning its logarithmic gradient."                         />
       <method method="ratioRotationCurve" description="Return the ratio of the square of the mid-plane rotation curve to that of the untruncated, unheated disk at the given dimensionless radius."/>
     </methods>
     !!]
     final     ::                                            exponentialDiskTidalDestructor
     procedure :: interpolate                             => exponentialDiskTidalInterpolate
     procedure :: ratioRotationCurve                      => exponentialDiskTidalRatioRotationCurve
     procedure :: isDimensionless                         => exponentialDiskTidalIsDimensionless
     procedure :: assumeMonotonicDecreasingSurfaceDensity => exponentialDiskTidalAssumeMonotonicDecreasingSurfaceDensity
     procedure :: massTotal                               => exponentialDiskTidalMassTotal
     procedure :: density                                 => exponentialDiskTidalDensity
     procedure :: densityGradientRadial                   => exponentialDiskTidalDensityGradientRadial
     procedure :: densitySlopeLogarithmicCentral          => exponentialDiskTidalDensitySlopeLogarithmicCentral
     procedure :: densitySphericalAverage                 => exponentialDiskTidalDensitySphericalAverage
     procedure :: surfaceDensity                          => exponentialDiskTidalSurfaceDensity
     procedure :: radiusHalfMass                          => exponentialDiskTidalRadiusHalfMass
     procedure :: massEnclosedBySphere                    => exponentialDiskTidalMassEnclosedBySphere
     procedure :: massEnclosedByCylinder                  => exponentialDiskTidalMassEnclosedByCylinder
     procedure :: potentialIsAnalytic                     => exponentialDiskTidalPotentialIsAnalytic
     procedure :: potential                               => exponentialDiskTidalPotential
     procedure :: rotationCurve                           => exponentialDiskTidalRotationCurve
     procedure :: rotationCurveGradient                   => exponentialDiskTidalRotationCurveGradient
     procedure :: surfaceDensityRadialMoment              => exponentialDiskTidalSurfaceDensityRadialMoment
     procedure :: acceleration                            => exponentialDiskTidalAcceleration
     procedure :: tidalTensor                             => exponentialDiskTidalTidalTensor
     procedure :: positionSample                          => exponentialDiskTidalPositionSample
  end type massDistributionExponentialDiskTidal

  interface massDistributionExponentialDiskTidal
     !!{RST
     Constructors for the :galacticus-class:`massDistributionExponentialDiskTidal` mass distribution class.
     !!}
     module procedure exponentialDiskTidalConstructorParameters
     module procedure exponentialDiskTidalConstructorInternal
  end interface massDistributionExponentialDiskTidal

  ! The list of tables of tidally-modified exponential disks. Tables are shared between all threads, and are created (or read
  ! from file) on first use.
  type            (exponentialDiskTidalTable), pointer   :: exponentialDiskTidalTables                  => null()

  ! Grid of truncation radii (in units of the disk scale length) and degree of truncation.
  double precision                           , parameter :: exponentialDiskTidalRadiusTruncateMinimum   =  1.0d-2
  integer                                    , parameter :: exponentialDiskTidalRadiusTruncatePerDecade =  20
  double precision                           , parameter :: exponentialDiskTidalTruncationNegligible    =  1.0d-6
  ! Grid of heating parameter, s, which is uniform in asinh(s/s₀).
  double precision                           , parameter :: exponentialDiskTidalHeatingScale            =  1.0d-3
  integer                                    , parameter :: exponentialDiskTidalHeatingPerDecade        =  20
  double precision                           , parameter :: exponentialDiskTidalHeatingMaximum          =  3.0d-1
  ! Grid of radius in units of the truncation radius, u=R/Rₜ.
  double precision                           , parameter :: exponentialDiskTidalRadiusMinimum           =  1.0d-3, exponentialDiskTidalRadiusMaximum         =3.0d2
  integer                                    , parameter :: exponentialDiskTidalRadiusPerDecade         =  60
  ! Grids (in units of the disk scale length) used in the FFTLog transforms and in the convolution which applies heating.
  integer                                    , parameter :: exponentialDiskTidalCountTransform          =  16384 , exponentialDiskTidalCountConvolution      =2000
  double precision                           , parameter :: exponentialDiskTidalTransformRadiusMinimum  =  1.0d-8, exponentialDiskTidalTransformRadiusMaximum=1.0d+6
  ! Radius (in units of the disk scale length) below which the ratio of rotation curves is held constant. The transforms ring near
  ! the lower end of their grid of radii, where the rotation curves (which scale as R² at small radii) are tiny. Well inside the
  ! disk both rotation curves scale as R², so their ratio is constant.
  double precision                           , parameter :: exponentialDiskTidalRadiusRatioFrozen       =  1.0d-4
  ! Radius (in units of the disk scale length) beyond which the mid-plane potential is taken to be Keplerian. This lies well
  ! beyond the largest radius tabulated, and well inside the upper end of the grid used in the transforms.
  double precision                           , parameter :: exponentialDiskTidalRadiusPotentialKeplerian=  1.0d+4
  ! Indices of tabulated quantities.
  integer                                    , parameter :: exponentialDiskTidalQuantityRotationCurve   =  1     , exponentialDiskTidalQuantityMass          =2     , &
       &                                                    exponentialDiskTidalQuantitySurfaceDensity  =  3     , exponentialDiskTidalQuantityPotential     =4
  ! Version of the format of the tabulation file. Increment if the grids above, or the content of the file, are changed.
  integer                                    , parameter :: exponentialDiskTidalFileFormatVersion       =  2
  ! Radii, in units of the root-mean-square radius of the disk, between which the field away from the mid-plane is blended from that
  ! of the scaled untruncated, unheated disk to that of the normalized multipole expansion.
  double precision                           , parameter :: exponentialDiskTidalRadiusBlendMinimum      =  2.0d0 , exponentialDiskTidalRadiusBlendMaximum    =5.0d0

contains

  function exponentialDiskTidalConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionExponentialDiskTidal` mass distribution class which builds the object
    from a parameter set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionExponentialDiskTidal)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (massDistributionClass               ), pointer       :: massDistribution_
    double precision                                                      :: truncation                   , heating, &
         &                                                                   radiusDecayOverRadiusTruncate

    !![
    <inputParameter docformat="rst">
      <name>truncation</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The degree of truncation, :math:`x=-\log f`, where :math:`f` is the fraction of the mass of the untruncated disk which is
      retained.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>heating</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The dimensionless heating parameter, :math:`s`, such that a ring of the disk at radius :math:`R` is spread in radius with a
      dispersion :math:`s R^2/R_\mathrm{d}`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusDecayOverRadiusTruncate</name>
      <defaultValue>0.1d0</defaultValue>
      <description>
      The ratio of the decay length of the truncation to the truncation radius.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="massDistribution" name="massDistribution_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionExponentialDisk)
       self=massDistributionExponentialDiskTidal(truncation,heating,radiusDecayOverRadiusTruncate,massDistribution_)
    class default
       call Error_Report('an exponential disk mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function exponentialDiskTidalConstructorParameters

  function exponentialDiskTidalConstructorInternal(truncation,heating,radiusDecayOverRadiusTruncate,massDistribution_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`massDistributionExponentialDiskTidal` mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionExponentialDiskTidal)                        :: self
    class           (massDistributionExponentialDisk     ), intent(in   ), target :: massDistribution_
    double precision                                      , intent(in   )         :: truncation                   , heating        , &
         &                                                                           radiusDecayOverRadiusTruncate
    double precision                                      , dimension(4)          :: weightTruncation             , weightHeating
    double precision                                                              :: positionTruncation           , positionHeating     , &
         &                                                                           truncationLogarithmic        , radiusSquaredMean   , &
         &                                                                           heightSquaredMean            , radiusTableMaximum
    integer                                                                       :: i                            , j
    !![
    <constructorAssign variables="truncation, heating, radiusDecayOverRadiusTruncate, *massDistribution_"/>
    !!]

    if (.not.self%massDistribution_%isDimensionless()) call Error_Report('the exponential disk must be dimensionless'            //{introspection:location})
    if (truncation < 0.0d0 .or. heating < 0.0d0      ) call Error_Report('truncation and heating must be non-negative'           //{introspection:location})
    if (self%massDistribution_%scaleHeight <= 0.0d0  ) call Error_Report('the exponential disk must have a positive scale height'//{introspection:location})
    self%componentType   =  massDistribution_%componentType
    self%     massType   =  massDistribution_%     massType
    self%massTotal_      =  exp(-truncation)*self%massDistribution_%massTotal()
    self%radiusHalfMass_ =  -1.0d0
    self%table           => exponentialDiskTidalTableGet(massDistribution_%scaleHeight/massDistribution_%radiusScale,radiusDecayOverRadiusTruncate)
    ! Find the position of this disk in the grid of truncation radius. The degree of truncation decreases monotonically with
    ! increasing truncation radius.
    truncationLogarithmic=log(max(truncation,exponentialDiskTidalTruncationNegligible))
    self%massScaleLogarithmic=0.0d0
    if      (truncationLogarithmic >= self%table%truncationLogarithmic(                         1)) then
       ! The truncation is stronger than that of the smallest tabulated truncation radius. Retain the shape of the disk truncated at
       ! that radius, but scale its mass to that required.
       positionTruncation       =0.0d0
       self%massScaleLogarithmic=exp(self%table%truncationLogarithmic(1))-truncation
    else if (truncationLogarithmic <= self%table%truncationLogarithmic(self%table%countTruncation)) then
       positionTruncation=dble(self%table%countTruncation-1)
    else
       i=1
       j=self%table%countTruncation
       do while (j-i > 1)
          if (self%table%truncationLogarithmic((i+j)/2) > truncationLogarithmic) then
             i=(i+j)/2
          else
             j=(i+j)/2
          end if
       end do
       positionTruncation=+dble(i-1)                                                                 &
            &             +(self%table%truncationLogarithmic(i)-           truncationLogarithmic   ) &
            &             /(self%table%truncationLogarithmic(i)-self%table%truncationLogarithmic(j))
    end if
    self%radiusTruncateLogarithmic=+self%table%radiusTruncateLogarithmicMinimum  &
         &                         +self%table%radiusTruncateLogarithmicInterval &
         &                         *positionTruncation
    ! Find the position of this disk in the grid of heating parameter.
    positionHeating=min(                                                                                         &
         &              asinh(min(heating,exponentialDiskTidalHeatingMaximum)/exponentialDiskTidalHeatingScale)  &
         &              /self%table%heatingInterval                                                            , &
         &              dble(self%table%countHeating-1)                                                          &
         &             )
    ! Evaluate the interpolating weights.
    call exponentialDiskTidalWeights(positionTruncation,self%table%countTruncation,self%indexTruncation,weightTruncation)
    call exponentialDiskTidalWeights(positionHeating   ,self%table%countHeating   ,self%indexHeating   ,weightHeating   )
    do i=1,4
       do j=1,4
          self%weight(i,j)=weightTruncation(i)*weightHeating(j)
       end do
    end do
    ! Evaluate the quadrupole moment of the disk, Q = M(⟨z²⟩-⟨R²⟩/2), used beyond the tabulated range of radii. For the sech²
    ! vertical profile ⟨z²⟩=π²h²/12, while ⟨R²⟩=2π∫Σ R³ dR/M is interpolated from the tabulated moments.
    radiusSquaredMean=0.0d0
    do i=1,4
       do j=1,4
          radiusSquaredMean=radiusSquaredMean+2.0d0*Pi*self%weight(i,j)*self%table%moment3(self%indexTruncation(i),self%indexHeating(j))
       end do
    end do
    heightSquaredMean=Pi**2*self%massDistribution_%scaleHeight**2/12.0d0
    self%quadrupole  =+self%massTotal_              &
         &            *(                            &
         &              +      heightSquaredMean    &
         &              -0.5d0*radiusSquaredMean    &
         &             )
    ! Find the radii between which the field away from the mid-plane is blended from that of the scaled untruncated, unheated disk
    ! to that of the normalized multipole expansion. These are set by the root-mean-square radius of the disk, but limited to lie
    ! within the tabulated range.
    radiusTableMaximum     =+exp(                                                                     &
         &                       +self%radiusTruncateLogarithmic                                      &
         &                       +self%table%radiusLogarithmicMinimum                                 &
         &                       +self%table%radiusLogarithmicInterval*dble(self%table%countRadius-1) &
         &                      )
    self%radiusBlendMaximum=+min(                                                                                   &
         &                       +exponentialDiskTidalRadiusBlendMaximum*sqrt(radiusSquaredMean+heightSquaredMean), &
         &                       +radiusTableMaximum                                                                &
         &                      )                                                                                   &
         &                  *self%massDistribution_%radiusScale
    self%radiusBlendMinimum=+self%radiusBlendMaximum                &
         &                  *exponentialDiskTidalRadiusBlendMinimum &
         &                  /exponentialDiskTidalRadiusBlendMaximum
    return
  end function exponentialDiskTidalConstructorInternal

  subroutine exponentialDiskTidalDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`massDistributionExponentialDiskTidal` mass distribution class.
    !!}
    implicit none
    type(massDistributionExponentialDiskTidal), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine exponentialDiskTidalDestructor

  subroutine exponentialDiskTidalWeights(position,count,indices,weights,derivative)
    !!{RST
    Compute the indices and weights for four-point Lagrange interpolation at the given ``position`` (zero-based, in units of
    the grid spacing) in a uniform grid of ``count`` points. The position is limited to lie within the grid. If ``derivative``
    is true, the weights for the derivative with respect to position are returned instead.
    !!}
    implicit none
    double precision              , intent(in   )           :: position
    integer                       , intent(in   )           :: count
    integer         , dimension(4), intent(  out)           :: indices
    double precision, dimension(4), intent(  out)           :: weights
    logical                       , intent(in   ), optional :: derivative
    double precision                                        :: p
    integer                                                 :: i
    !![
    <optionalArgument name="derivative" defaultsTo=".false."/>
    !!]

    p      =max(0.0d0,min(position,dble(count-1)))
    i      =min(max(int(p)+1,2),count-2)
    p      =p-dble(i-1)
    indices=[i-1,i,i+1,i+2]
    if (derivative_) then
       weights=[                                   &
            &   -(3.0d0*p**2-6.0d0*p+2.0d0)/6.0d0, &
            &   +(3.0d0*p**2-4.0d0*p-1.0d0)/2.0d0, &
            &   -(3.0d0*p**2-2.0d0*p-2.0d0)/2.0d0, &
            &   +(3.0d0*p**2        -1.0d0)/6.0d0  &
            &  ]
    else
       weights=[                                      &
            &   -(p      )*(p-1.0d0)*(p-2.0d0)/6.0d0, &
            &   +(p+1.0d0)*(p-1.0d0)*(p-2.0d0)/2.0d0, &
            &   -(p+1.0d0)*(p      )*(p-2.0d0)/2.0d0, &
            &   +(p+1.0d0)*(p      )*(p-1.0d0)/6.0d0  &
            &  ]
    end if
    return
  end subroutine exponentialDiskTidalWeights

  function exponentialDiskTidalTableGet(scaleHeight,radiusDecayOverRadiusTruncate) result(table)
    !!{RST
    Return a pointer to the table for the given dimensionless scale height and truncation decay length, creating it if
    necessary.
    !!}
    implicit none
    type            (exponentialDiskTidalTable), pointer       :: table
    double precision                           , intent(in   ) :: scaleHeight, radiusDecayOverRadiusTruncate

    !$omp critical (massDistributionExponentialDiskTidalTables)
    table => exponentialDiskTidalTables
    do while (associated(table))
       if (table%scaleHeight == scaleHeight .and. table%radiusDecayOverRadiusTruncate == radiusDecayOverRadiusTruncate) exit
       table => table%next
    end do
    if (.not.associated(table)) then
       allocate(table)
       table%scaleHeight                  =  scaleHeight
       table%radiusDecayOverRadiusTruncate=  radiusDecayOverRadiusTruncate
       call exponentialDiskTidalTabulate(table)
       table%next                         => exponentialDiskTidalTables
       exponentialDiskTidalTables         => table
    end if
    !$omp end critical (massDistributionExponentialDiskTidalTables)
    return
  end function exponentialDiskTidalTableGet

  subroutine exponentialDiskTidalTabulate(table)
    !!{RST
    Build (or read from file) the tables for a tidally truncated and heated exponential disk.
    !!}
    use :: Display                                  , only : displayIndent                       , displayUnindent                       , verbosityLevelWorking
    use :: Error                                    , only : Error_Report
    use :: FFTLogs                                  , only : FFTLog                              , fftLogForward
    use :: File_Utilities                           , only : Directory_Make                      , File_Exists                           , File_Lock               , File_Path, &
         &                                                   File_Unlock                         , lockDescriptor
    use :: HDF5_Access                              , only : hdf5Access
    use :: IO_HDF5                                  , only : hdf5File
    use :: ISO_Varying_String                       , only : varying_string                      , operator(//)                          , char
    use :: Input_Paths                              , only : inputPath                           , pathTypeDataDynamic
    use :: Mass_Distributions_Cylindrical_Utilities , only : Kernel_Vertical_Integrators_Allocate, Kernel_Vertical_Integrators_Deallocate, kernelVerticalSechSquared
    use :: Numerical_Constants_Math                 , only : Pi
    use :: String_Handling                          , only : operator(//)
    implicit none
    type            (exponentialDiskTidalTable), intent(inout)               :: table
    double precision                           , allocatable  , dimension(:) :: radiusTransform           , wavenumberTransform         , &
         &                                                                      kernel                    , rotationCurveSquaredParent  , &
         &                                                                      truncationLogarithmic     , massCumulativeParent        , &
         &                                                                      potentialParent
    type            (kernelVerticalSechSquared)                              :: kernel_
    type            (varying_string           )                              :: fileName
    type            (lockDescriptor           )                              :: fileLock
    character       (len=8                    )                              :: labelHeight               , labelDecay
    double precision                                                         :: radiusTruncate            , intervalLogarithmicTransform, &
         &                                                                      massParent
    integer                                                                  :: i                         , iTruncation                 , &
         &                                                                      iHeating                  , iPlane                      , &
         &                                                                      countTruncationMaximum
    logical                                                                  :: tableRead                 , monotonic

    ! Set the grids which do not depend on the degree of truncation.
    table%heatingInterval                  =log(10.0d0)/dble(exponentialDiskTidalHeatingPerDecade)
    table%countHeating                     =int(asinh(exponentialDiskTidalHeatingMaximum/exponentialDiskTidalHeatingScale)/table%heatingInterval)+2
    table%radiusLogarithmicMinimum         =log(exponentialDiskTidalRadiusMinimum)
    table%radiusLogarithmicInterval        =log(10.0d0)/dble(exponentialDiskTidalRadiusPerDecade)
    table%countRadius                      =int(log(exponentialDiskTidalRadiusMaximum/exponentialDiskTidalRadiusMinimum)/table%radiusLogarithmicInterval)+1
    table%radiusTruncateLogarithmicMinimum =log(exponentialDiskTidalRadiusTruncateMinimum)
    table%radiusTruncateLogarithmicInterval=log(10.0d0)/dble(exponentialDiskTidalRadiusTruncatePerDecade)
    ! Construct a file name for the table.
    write (labelHeight,'(f8.6)') table%scaleHeight
    write (labelDecay ,'(f8.6)') table%radiusDecayOverRadiusTruncate
    fileName=inputPath(pathTypeDataDynamic)                                                                    // &
         &   'galacticStructure/massDistributionExponentialDiskTidal_h'//trim(adjustl(labelHeight))            // &
         &   '_d'                                                      //trim(adjustl(labelDecay ))            // &
         &   '_v'                                                      //exponentialDiskTidalFileFormatVersion // &
         &   '.hdf5'
    call Directory_Make(File_Path(fileName))
    ! Read the table from file if it exists, holding a shared lock while doing so. Otherwise, obtain an exclusive lock before
    ! building and writing the table, and check for the file again once that lock is obtained (as another process may have built
    ! it in the meantime). Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    tableRead=.false.
    if (File_Exists(fileName)) then
       call File_Lock  (fileName,fileLock,lockIsShared=.true.)
       call tableReadFile()
       call File_Unlock(         fileLock                    )
       tableRead=.true.
    end if
    if (.not.tableRead) call File_Lock(fileName,fileLock,lockIsShared=.false.)
    if (tableRead) then
       ! Table was read from file - nothing more to do.
    else if (File_Exists(fileName)) then
       call tableReadFile()
    else
       call displayIndent("tabulating tidally-modified exponential disk",verbosityLevelWorking)
       ! Construct the grid of radii used in the transforms. The transform of the surface density yields a grid of wavenumbers,
       ! and the transform back yields the original radii.
       allocate(radiusTransform           (exponentialDiskTidalCountTransform))
       allocate(wavenumberTransform       (exponentialDiskTidalCountTransform))
       allocate(rotationCurveSquaredParent(exponentialDiskTidalCountTransform))
       allocate(kernel                    (exponentialDiskTidalCountTransform))
       intervalLogarithmicTransform=log(exponentialDiskTidalTransformRadiusMaximum/exponentialDiskTidalTransformRadiusMinimum)/dble(exponentialDiskTidalCountTransform-1)
       do i=1,exponentialDiskTidalCountTransform
          radiusTransform(i)=exp(log(exponentialDiskTidalTransformRadiusMinimum)+dble(i-1)*intervalLogarithmicTransform)
       end do
       ! Find the wavenumbers of the transform, and evaluate the vertical kernel at those wavenumbers.
       call FFTLog(radiusTransform,wavenumberTransform,exp(-radiusTransform),kernel,0.0d0,fftLogForward)
       kernel_=kernelVerticalSechSquared(table%scaleHeight,[0.0d0])
       call Kernel_Vertical_Integrators_Allocate  ()
       do i=1,exponentialDiskTidalCountTransform
          kernel(i)=kernel_%value(wavenumberTransform(i),1,1)
       end do
       call Kernel_Vertical_Integrators_Deallocate()
       ! Compute the mass, cumulative mass, and mid-plane rotation curve, of the untruncated, unheated disk.
       allocate(massCumulativeParent(exponentialDiskTidalCountTransform))
       call massCumulative_(surfaceDensityTruncated(radiusTransform,huge(0.0d0)),massCumulativeParent)
       massParent=massCumulativeParent(exponentialDiskTidalCountTransform)
       call rotationCurveSquared(surfaceDensityTruncated(radiusTransform,huge(0.0d0)),rotationCurveSquaredParent)
       allocate(potentialParent(exponentialDiskTidalCountTransform))
       call potentialMidplane(rotationCurveSquaredParent,massCumulativeParent,potentialParent)
       ! Find the grid of truncation radii, extending it until the degree of truncation becomes negligible.
       countTruncationMaximum=int(log(1.0d4)/table%radiusTruncateLogarithmicInterval)
       allocate(truncationLogarithmic(countTruncationMaximum))
       table%countTruncation=0
       do while (table%countTruncation < countTruncationMaximum)
          table%countTruncation=table%countTruncation+1
          radiusTruncate       =exp(table%radiusTruncateLogarithmicMinimum+dble(table%countTruncation-1)*table%radiusTruncateLogarithmicInterval)
          truncationLogarithmic(table%countTruncation)=log(-log(1.0d0-massRemoved(radiusTruncate)/massParent))
          if (table%countTruncation > 1) then
             if (truncationLogarithmic(table%countTruncation) >= truncationLogarithmic(table%countTruncation-1)) &
                  & call Error_Report('degree of truncation is not monotonic in truncation radius'//{introspection:location})
          end if
          ! Continue for three further points beyond the point at which truncation becomes negligible, so that four-point
          ! interpolation is possible throughout the range of interest.
          if (table%countTruncation > 3) then
             if (truncationLogarithmic(table%countTruncation-3) < log(exponentialDiskTidalTruncationNegligible)) exit
          end if
       end do
       if (table%countTruncation >= countTruncationMaximum) call Error_Report('failed to reach negligible truncation'//{introspection:location})
       table%truncationLogarithmic=truncationLogarithmic(1:table%countTruncation)
       allocate(table%moment1                       (                  table%countTruncation,table%countHeating))
       allocate(table%moment2                       (                  table%countTruncation,table%countHeating))
       allocate(table%moment3                       (                  table%countTruncation,table%countHeating))
       allocate(table%rotationCurveRatioLogarithmic (table%countRadius,table%countTruncation,table%countHeating))
       allocate(table%massRatioLogarithmic          (table%countRadius,table%countTruncation,table%countHeating))
       allocate(table%surfaceDensityRatioLogarithmic(table%countRadius,table%countTruncation,table%countHeating))
       allocate(table%potentialRatioLogarithmic     (table%countRadius,table%countTruncation,table%countHeating))
       ! Tabulate each plane of truncation radius and heating parameter.
       monotonic=.true.
       !$omp parallel do schedule(dynamic) private(iTruncation,iHeating) reduction(.and.:monotonic)
       do iPlane=1,table%countTruncation*table%countHeating
          iTruncation=(iPlane-1)/table%countHeating+1
          iHeating   =mod(iPlane-1,table%countHeating)+1
          call tabulatePlane(iTruncation,iHeating,monotonic)
       end do
       !$omp end parallel do
       table%monotonic=monotonic
       ! Write the table to file.
       !$ call hdf5Access%set()
       hdf5WriteScope: block
         type(hdf5File) :: file
         file=hdf5File(char(fileName),overWrite=.true.,readOnly=.false.)
         call file%writeAttribute(exponentialDiskTidalFileFormatVersion  ,'fileFormat'                       )
         call file%writeAttribute(table%countTruncation                  ,'countTruncation'                  )
         call file%writeAttribute(table%countHeating                     ,'countHeating'                     )
         call file%writeAttribute(table%countRadius                      ,'countRadius'                      )
         call file%writeAttribute(table%radiusTruncateLogarithmicMinimum ,'radiusTruncateLogarithmicMinimum' )
         call file%writeAttribute(table%radiusTruncateLogarithmicInterval,'radiusTruncateLogarithmicInterval')
         call file%writeAttribute(table%heatingInterval                  ,'heatingInterval'                  )
         call file%writeAttribute(table%radiusLogarithmicMinimum         ,'radiusLogarithmicMinimum'         )
         call file%writeAttribute(table%radiusLogarithmicInterval        ,'radiusLogarithmicInterval'        )
         if (table%monotonic) then
            call file%writeAttribute(1,'monotonic')
         else
            call file%writeAttribute(0,'monotonic')
         end if
         call file%writeDataset(table%truncationLogarithmic         ,'truncationLogarithmic'         )
         call file%writeDataset(table%moment1                       ,'moment1'                       )
         call file%writeDataset(table%moment2                       ,'moment2'                       )
         call file%writeDataset(table%moment3                       ,'moment3'                       )
         call file%writeDataset(table%rotationCurveRatioLogarithmic ,'rotationCurveRatioLogarithmic' )
         call file%writeDataset(table%massRatioLogarithmic          ,'massRatioLogarithmic'          )
         call file%writeDataset(table%surfaceDensityRatioLogarithmic,'surfaceDensityRatioLogarithmic')
         call file%writeDataset(table%potentialRatioLogarithmic     ,'potentialRatioLogarithmic'     )
       end block hdf5WriteScope
       !$ call hdf5Access%unset()
       call displayUnindent("done",verbosityLevelWorking)
    end if
    if (.not.tableRead) call File_Unlock(fileLock)
    return

  contains

    subroutine tableReadFile()
      !!{RST
      Read the table from file. The caller must hold a lock on the file.
      !!}
      implicit none
      integer :: fileFormat, monotonic

      !$ call hdf5Access%set()
      hdf5ReadScope: block
        type(hdf5File) :: file
        file=hdf5File(char(fileName),readOnly=.true.)
        call file%readAttribute('fileFormat'                       ,fileFormat                             )
        if (fileFormat /= exponentialDiskTidalFileFormatVersion) call Error_Report('incompatible file format'//{introspection:location})
        call file%readAttribute('countTruncation'                  ,table%countTruncation                  )
        call file%readAttribute('countHeating'                     ,table%countHeating                     )
        call file%readAttribute('countRadius'                      ,table%countRadius                      )
        call file%readAttribute('radiusTruncateLogarithmicMinimum' ,table%radiusTruncateLogarithmicMinimum )
        call file%readAttribute('radiusTruncateLogarithmicInterval',table%radiusTruncateLogarithmicInterval)
        call file%readAttribute('heatingInterval'                  ,table%heatingInterval                  )
        call file%readAttribute('radiusLogarithmicMinimum'         ,table%radiusLogarithmicMinimum         )
        call file%readAttribute('radiusLogarithmicInterval'        ,table%radiusLogarithmicInterval        )
        call file%readAttribute('monotonic'                        ,monotonic                              )
        table%monotonic=monotonic == 1
        call file%readDataset  ('truncationLogarithmic'         ,table%truncationLogarithmic         )
        call file%readDataset  ('moment1'                       ,table%moment1                       )
        call file%readDataset  ('moment2'                       ,table%moment2                       )
        call file%readDataset  ('moment3'                       ,table%moment3                       )
        call file%readDataset  ('rotationCurveRatioLogarithmic' ,table%rotationCurveRatioLogarithmic )
        call file%readDataset  ('massRatioLogarithmic'          ,table%massRatioLogarithmic          )
        call file%readDataset  ('surfaceDensityRatioLogarithmic',table%surfaceDensityRatioLogarithmic)
        call file%readDataset  ('potentialRatioLogarithmic'     ,table%potentialRatioLogarithmic     )
      end block hdf5ReadScope
      !$ call hdf5Access%unset()
      return
    end subroutine tableReadFile

    function surfaceDensityTruncatedLogarithmic(radius,radiusTruncate) result(surfaceDensityLogarithmic)
      !!{RST
      Return the logarithm of the surface density of the truncated (but unheated) disk of unit mass and scale length at the
      given radii. Evaluating the logarithm avoids underflow far beyond the truncation radius.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:           ) :: radius
      double precision, intent(in   )                          :: radiusTruncate
      double precision               , dimension(size(radius)) :: surfaceDensityLogarithmic
      integer                                                  :: i

      do i=1,size(radius)
         surfaceDensityLogarithmic(i)=-log(2.0d0*Pi)-radius(i)
         if (radius(i) > radiusTruncate)                                          &
              & surfaceDensityLogarithmic(i)=+surfaceDensityLogarithmic(i)        &
              &                              +(                                   &
              &                                +log(radius(i)/radiusTruncate)     &
              &                                -(radius(i)/radiusTruncate-1.0d0)  &
              &                               )                                   &
              &                              /table%radiusDecayOverRadiusTruncate
      end do
      return
    end function surfaceDensityTruncatedLogarithmic

    function surfaceDensityTruncated(radius,radiusTruncate) result(surfaceDensity)
      !!{RST
      Return the surface density of the truncated (but unheated) disk of unit mass and scale length at the given radii.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:           ) :: radius
      double precision, intent(in   )                          :: radiusTruncate
      double precision               , dimension(size(radius)) :: surfaceDensity

      surfaceDensity=exp(surfaceDensityTruncatedLogarithmic(radius,radiusTruncate))
      return
    end function surfaceDensityTruncated

    double precision function massRemoved(radiusTruncate) result(mass)
      !!{RST
      Return the mass removed by truncation from the disk of unit (untruncated) mass and scale length, computed on the grid of
      radii used in the transforms. Computing the removed mass directly (rather than as the difference between the masses of the
      untruncated and truncated disks) retains precision when the truncation is small.
      !!}
      implicit none
      double precision, intent(in   ) :: radiusTruncate

      mass=+sum(                                                                           &
           &    +2.0d0*Pi*radiusTransform**2                                               &
           &    *(                                                                         &
           &      +exp(surfaceDensityTruncatedLogarithmic(radiusTransform,huge(0.0d0)   )) &
           &      -exp(surfaceDensityTruncatedLogarithmic(radiusTransform,radiusTruncate)) &
           &     )                                                                         &
           &   )                                                                           &
           & *intervalLogarithmicTransform
      return
    end function massRemoved

    subroutine massCumulative_(surfaceDensity,massCumulative)
      !!{RST
      Compute the cumulative mass of a disk with the given surface density, tabulated on the grid of radii used in the transforms,
      by trapezoidal integration in the logarithm of radius.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:) :: surfaceDensity
      double precision, intent(  out), dimension(:) :: massCumulative
      integer                                       :: j

      massCumulative(1)=0.5d0*Pi*radiusTransform(1)**2*surfaceDensity(1)
      do j=2,size(surfaceDensity)
         massCumulative(j)=+massCumulative(j-1)                           &
              &            +Pi                                            &
              &            *(                                             &
              &              +radiusTransform(j-1)**2*surfaceDensity(j-1) &
              &              +radiusTransform(j  )**2*surfaceDensity(j  ) &
              &             )                                             &
              &            *intervalLogarithmicTransform
      end do
      return
    end subroutine massCumulative_

    subroutine rotationCurveSquared(surfaceDensity,rotationCurveSquared_)
      !!{RST
      Return the square of the mid-plane rotation curve (for unit gravitational constant) of a disk with the given surface density
      (tabulated on the grid of radii used in the transforms) and the vertical structure of the exponential disk. This is
      :math:`V_\mathrm{c}^2(R) = 2\pi R \int_0^\infty k J_1(kR) \tilde{\Sigma}(k) K(k,0) \mathrm{d}k`, where
      :math:`\tilde{\Sigma}(k)=\int_0^\infty \Sigma(R) J_0(kR) R \mathrm{d}R` and :math:`K(k,0)` is the mid-plane vertical kernel.
      Since ``FFTLog()`` computes :math:`\int_0^\infty f(r) (kr)^{1/2} J_\mu(kr) \mathrm{d}r`, the transforms are applied to
      :math:`\Sigma(R) R^{1/2}` (giving :math:`\tilde{\Sigma}(k) k^{1/2}`), and to :math:`\tilde{\Sigma}(k) k^{1/2} K(k,0)`
      (giving :math:`V_\mathrm{c}^2(R) R^{-1/2}/2\pi`).
      !!}
      implicit none
      double precision, intent(in   ), dimension(:) :: surfaceDensity
      double precision, intent(  out), dimension(:) :: rotationCurveSquared_
      double precision, allocatable  , dimension(:) :: transform            , wavenumber, &
           &                                           radius

      allocate(transform (size(surfaceDensity)))
      allocate(wavenumber(size(surfaceDensity)))
      allocate(radius    (size(surfaceDensity)))
      call FFTLog(radiusTransform,wavenumber,surfaceDensity*sqrt(radiusTransform),transform            ,0.0d0,fftLogForward)
      call FFTLog(wavenumber     ,radius    ,transform     *kernel               ,rotationCurveSquared_,1.0d0,fftLogForward)
      rotationCurveSquared_=2.0d0*Pi*sqrt(radius)*rotationCurveSquared_
      return
    end subroutine rotationCurveSquared

    subroutine potentialMidplane(rotationCurveSquared_,massCumulative,potential)
      !!{RST
      Return the mid-plane potential (for unit gravitational constant) of a disk with the given square of the mid-plane rotation
      curve and cumulative mass (both tabulated on the grid of radii used in the transforms). Since :math:`\partial \Phi/\partial
      R = V_\mathrm{c}^2/R` in the mid-plane, the potential is :math:`\Phi(R,0) = -\int_R^\infty V_\mathrm{c}^2(R^\prime)
      \mathrm{d}\ln R^\prime`. This is evaluated by trapezoidal integration in the logarithm of radius, inward from a radius,
      :math:`R_\mathrm{K}`, beyond which the disk is treated as a point mass, contributing :math:`-M(R_\mathrm{K})/R_\mathrm{K}`.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:) :: rotationCurveSquared_, massCumulative
      double precision, intent(  out), dimension(:) :: potential
      integer                                       :: j                    , jKeplerian

      jKeplerian=int(log(exponentialDiskTidalRadiusPotentialKeplerian/exponentialDiskTidalTransformRadiusMinimum)/intervalLogarithmicTransform)+1
      do j=size(potential),jKeplerian,-1
         potential(j)=-massCumulative(j)/radiusTransform(j)
      end do
      do j=jKeplerian-1,1,-1
         potential(j)=+potential(j+1)                &
              &       -0.5d0                         &
              &       *(                             &
              &         +rotationCurveSquared_(j  )  &
              &         +rotationCurveSquared_(j+1)  &
              &        )                             &
              &       *intervalLogarithmicTransform
      end do
      return
    end subroutine potentialMidplane

    subroutine tabulatePlane(iTruncation,iHeating,monotonic)
      !!{RST
      Tabulate the properties of the disk with the given indices of truncation radius and heating parameter.
      !!}
      implicit none
      integer                  , intent(in   )               :: iTruncation                             , iHeating
      logical                  , intent(inout)               :: monotonic
      double precision         , parameter                   :: monotonicTolerance               =1.0d-6
      double precision         , allocatable  , dimension(:) :: surfaceDensity                          , surfaceDensityLogarithmic        , &
           &                                                    massCumulative                          , rotationCurveSquaredRatio        , &
           &                                                    potentialRatio
      double precision                                       :: radiusTruncate                          , heating                          , &
           &                                                    massTotal                               , radiusLogarithmic                , &
           &                                                    position                                , massRatio                        , &
           &                                                    surfaceDensityLogarithmicModified       , surfaceDensityLogarithmicPrevious
      integer                                                :: iRadius                                 , j
      character       (len=256)                              :: message

      allocate(surfaceDensity           (exponentialDiskTidalCountTransform))
      allocate(surfaceDensityLogarithmic(exponentialDiskTidalCountTransform))
      allocate(massCumulative           (exponentialDiskTidalCountTransform))
      allocate(rotationCurveSquaredRatio(exponentialDiskTidalCountTransform))
      allocate(potentialRatio           (exponentialDiskTidalCountTransform))
      radiusTruncate=exp(table%radiusTruncateLogarithmicMinimum+dble(iTruncation-1)*table%radiusTruncateLogarithmicInterval)
      heating       =exponentialDiskTidalHeatingScale*sinh(dble(iHeating-1)*table%heatingInterval)
      ! Find the surface density of the truncated disk, and apply heating.
      surfaceDensityLogarithmic=max(surfaceDensityTruncatedLogarithmic(radiusTransform,radiusTruncate),log(tiny(0.0d0)))
      if (iHeating > 1) call heat(heating,radiusTruncate,surfaceDensityLogarithmic)
      surfaceDensity=exp(surfaceDensityLogarithmic)
      ! Find the cumulative mass, and the radial moments.
      call massCumulative_(surfaceDensity,massCumulative)
      massTotal=massCumulative(exponentialDiskTidalCountTransform)
      table%moment1(iTruncation,iHeating)=sum(surfaceDensity*radiusTransform**2)*intervalLogarithmicTransform/massTotal
      table%moment2(iTruncation,iHeating)=sum(surfaceDensity*radiusTransform**3)*intervalLogarithmicTransform/massTotal
      table%moment3(iTruncation,iHeating)=sum(surfaceDensity*radiusTransform**4)*intervalLogarithmicTransform/massTotal
      ! Find the rotation curve relative to that of the untruncated, unheated disk.
      call rotationCurveSquared(surfaceDensity,rotationCurveSquaredRatio)
      rotationCurveSquaredRatio=rotationCurveSquaredRatio/rotationCurveSquaredParent
      j=int(log(exponentialDiskTidalRadiusRatioFrozen/exponentialDiskTidalTransformRadiusMinimum)/intervalLogarithmicTransform)+1
      rotationCurveSquaredRatio(1:j-1)=rotationCurveSquaredRatio(j)
      ! Find the mid-plane potential relative to that of the untruncated, unheated disk. The rotation curve used here is that
      ! reconstructed from the ratio, in which ringing at small radii has been removed.
      call potentialMidplane(rotationCurveSquaredRatio*rotationCurveSquaredParent,massCumulative,potentialRatio)
      potentialRatio=potentialRatio/potentialParent
      ! Interpolate onto the grid of radii in units of the truncation radius.
      do iRadius=1,table%countRadius
         radiusLogarithmic=+table%radiusLogarithmicMinimum                  &
              &            +table%radiusLogarithmicInterval*dble(iRadius-1) &
              &            +log(radiusTruncate)
         position         =(radiusLogarithmic-log(exponentialDiskTidalTransformRadiusMinimum))/intervalLogarithmicTransform
         j                =min(max(int(position)+1,1),exponentialDiskTidalCountTransform-1)
         position         =position-dble(j-1)
         if (rotationCurveSquaredRatio(j) <= 0.0d0 .or. rotationCurveSquaredRatio(j+1) <= 0.0d0) then
            write (message,'(a,e12.6,a,e12.6,a,e12.6,a,e12.6)') 'non-positive rotation curve in tidally-modified exponential disk: Rₜ=',radiusTruncate,' s=',heating,' R=',exp(radiusLogarithmic),' V²/V²₀=',min(rotationCurveSquaredRatio(j),rotationCurveSquaredRatio(j+1))
            call Error_Report(trim(message)//{introspection:location})
         end if
         table%rotationCurveRatioLogarithmic(iRadius,iTruncation,iHeating)=  +(1.0d0-position)*log(rotationCurveSquaredRatio(j  )) &
              &                                                              +       position *log(rotationCurveSquaredRatio(j+1))
         table%potentialRatioLogarithmic    (iRadius,iTruncation,iHeating)=  +(1.0d0-position)*log(potentialRatio           (j  )) &
              &                                                              +       position *log(potentialRatio           (j+1))
         massRatio                                                        =+(                                                      &
              &                                                              +(1.0d0-position)*    massCumulative           (j  )  &
              &                                                              +       position *    massCumulative           (j+1)  &
              &                                                             )                                                      &
              &                                                            /(                                                      &
              &                                                              +(1.0d0-position)*    massCumulativeParent     (j  )  &
              &                                                              +       position *    massCumulativeParent     (j+1)  &
              &                                                             )
         surfaceDensityLogarithmicModified                                =  +(1.0d0-position)*    surfaceDensityLogarithmic(j  ) &
              &                                                              +       position *    surfaceDensityLogarithmic(j+1)
         table%massRatioLogarithmic         (iRadius,iTruncation,iHeating)=log(massRatio)
         if (surfaceDensityLogarithmicModified <= log(tiny(0.0d0))) then
            ! The surface density has underflowed - record a ratio which will also underflow.
            table%surfaceDensityRatioLogarithmic(iRadius,iTruncation,iHeating)=log(tiny(0.0d0))
         else
            table%surfaceDensityRatioLogarithmic(iRadius,iTruncation,iHeating)=+surfaceDensityLogarithmicModified &
                 &                                                             +log(2.0d0*Pi)                     &
                 &                                                             +exp(radiusLogarithmic)
         end if
         if (iRadius > 1) then
            if (surfaceDensityLogarithmicModified > surfaceDensityLogarithmicPrevious+monotonicTolerance) monotonic=.false.
         end if
         surfaceDensityLogarithmicPrevious=surfaceDensityLogarithmicModified
      end do
      return
    end subroutine tabulatePlane

    subroutine heat(heating,radiusTruncate,surfaceDensityLogarithmic)
      !!{RST
      Apply heating to the surface density (given as its logarithm on the grid of radii used in the transforms). The mass of each
      cell of a coarse logarithmic grid in radius is spread with a log-normal distribution in radius, of width
      :math:`\sigma_{\ln R}=sR`, and mean :math:`\ln R-\sigma_{\ln R}^2/2` (such that the mean radius is conserved). The change in
      the logarithm of the surface density on the coarse grid is then interpolated onto the grid used in the transforms (and added
      to the logarithm of the unheated surface density there), such that features of the unheated surface density on scales
      smaller than the coarse grid are preserved.
      !!}
      implicit none
      double precision, intent(in   )                                                  :: heating                              , radiusTruncate
      double precision, intent(inout), dimension(exponentialDiskTidalCountTransform  ) :: surfaceDensityLogarithmic
      double precision, parameter                                                      :: widthKernel                    =8.0d0, massNegligible                 =1.0d-20, &
           &                                                                              changeMaximum                  =1.0d1, widthMaximum                   =1.0d+00
      double precision               , dimension(exponentialDiskTidalCountConvolution) :: radius                               , massCell                               , &
           &                                                                              massCellHeated                       , changeLogarithmic                      , &
           &                                                                              surfaceDensityLogarithmicCoarse      , surfaceDensityLogarithmicHeated
      logical                        , dimension(exponentialDiskTidalCountConvolution) :: isEmpty
      double precision                                                                 :: interval                            , width                                  , &
           &                                                                              mean                                , lower                                  , &
           &                                                                              upper                               , position
      integer                                                                          :: i                                   , k                                      , &
           &                                                                              kMinimum                            , kMaximum                               , &
           &                                                                              j

      ! Construct the coarse grid, and find the mass in each cell.
      interval=log(exponentialDiskTidalTransformRadiusMaximum/exponentialDiskTidalTransformRadiusMinimum)/dble(exponentialDiskTidalCountConvolution-1)
      do i=1,exponentialDiskTidalCountConvolution
         radius(i)=exp(log(exponentialDiskTidalTransformRadiusMinimum)+dble(i-1)*interval)
      end do
      surfaceDensityLogarithmicCoarse=surfaceDensityTruncatedLogarithmic(radius,radiusTruncate)
      massCell                       =2.0d0*Pi*radius**2*exp(surfaceDensityLogarithmicCoarse)*interval
      ! Spread the mass of each cell.
      massCellHeated=0.0d0
      do i=1,exponentialDiskTidalCountConvolution
         width=min(heating*radius(i),widthMaximum)
         if (massCell(i) < massNegligible .or. width < 0.25d0*interval) then
            ! The cell has negligible mass, or the kernel is narrow compared to the grid - leave the mass in its cell.
            massCellHeated(i)=massCellHeated(i)+massCell(i)
            cycle
         end if
         mean    =log(radius(i))-0.5d0*width**2
         kMinimum=max(1                                   ,int((mean-widthKernel*width-log(radius(1)))/interval)+1)
         kMaximum=min(exponentialDiskTidalCountConvolution,int((mean+widthKernel*width-log(radius(1)))/interval)+2)
         do k=kMinimum,kMaximum
            lower=(log(radius(k))-0.5d0*interval-mean)/width/sqrt(2.0d0)
            upper=(log(radius(k))+0.5d0*interval-mean)/width/sqrt(2.0d0)
            ! Evaluate the fraction of the mass in this cell using complementary error functions on the appropriate side of the
            ! mean, to avoid loss of precision far from the mean.
            if (lower > 0.0d0) then
               massCellHeated(k)=massCellHeated(k)+0.5d0*massCell(i)*(      erfc( lower)-erfc( upper))
            else if (upper < 0.0d0) then
               massCellHeated(k)=massCellHeated(k)+0.5d0*massCell(i)*(      erfc(-upper)-erfc(-lower))
            else
               massCellHeated(k)=massCellHeated(k)+0.5d0*massCell(i)*(2.0d0-erfc( upper)-erfc(-lower))
            end if
         end do
      end do
      ! Find the logarithm of the heated surface density, and its change due to heating, on the coarse grid.
      isEmpty=massCellHeated <= 0.0d0
      where (.not.isEmpty)
         surfaceDensityLogarithmicHeated=log(massCellHeated/(2.0d0*Pi*radius**2*interval))
         changeLogarithmic              =surfaceDensityLogarithmicHeated-surfaceDensityLogarithmicCoarse
      elsewhere
         surfaceDensityLogarithmicHeated=log(tiny(0.0d0))
         changeLogarithmic              =0.0d0
      end where
      ! Interpolate onto the grid used in the transforms. Where heating changes the surface density only modestly, the change is
      ! interpolated and added to the unheated surface density, preserving features (such as the truncation) on scales smaller than
      ! the coarse grid. Elsewhere (e.g. far beyond the truncation radius, where the unheated surface density falls by many orders of
      ! magnitude between points of the coarse grid, while the heated surface density varies smoothly) the logarithm of the heated
      ! surface density is interpolated directly.
      do i=1,exponentialDiskTidalCountTransform
         position=(log(radiusTransform(i))-log(radius(1)))/interval
         j       =min(max(int(position)+1,1),exponentialDiskTidalCountConvolution-1)
         position=position-dble(j-1)
         if (isEmpty(j) .or. isEmpty(j+1)) then
            surfaceDensityLogarithmic(i)=log(tiny(0.0d0))
         else if (abs(changeLogarithmic(j)) < changeMaximum .and. abs(changeLogarithmic(j+1)) < changeMaximum) then
            surfaceDensityLogarithmic(i)=+surfaceDensityLogarithmic(i)                          &
                 &                       +(1.0d0-position)*changeLogarithmic              (j  ) &
                 &                       +       position *changeLogarithmic              (j+1)
         else
            surfaceDensityLogarithmic(i)=+(1.0d0-position)*surfaceDensityLogarithmicHeated(j  ) &
                 &                       +       position *surfaceDensityLogarithmicHeated(j+1)
         end if
         surfaceDensityLogarithmic(i)=max(surfaceDensityLogarithmic(i),log(tiny(0.0d0)))
      end do
      return
    end subroutine heat

  end subroutine exponentialDiskTidalTabulate

  double precision function exponentialDiskTidalInterpolate(self,radius,quantity,gradientLogarithmic) result(ratio)
    !!{RST
    Interpolate the ratio of the tabulated ``quantity`` (the square of the mid-plane rotation curve, the mid-plane potential, the
    mass enclosed within a cylinder, or the surface density) to that of the untruncated, unheated exponential disk at the given dimensionless ``radius``
    (in units of the disk scale length). The logarithm of the ratio is interpolated. If present, ``gradientLogarithmic`` returns
    the derivative of the logarithm of the ratio with respect to :math:`\ln R`. Below the tabulated range the ratio is taken to
    be constant (as all of these quantities are unaffected by truncation and heating close to the center of the disk). Above the
    tabulated range the ratio at the upper limit of the table is returned, and callers must instead apply the appropriate
    asymptotic behavior.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout)           :: self
    double precision                                      , intent(in   )           :: radius
    integer                                               , intent(in   )           :: quantity
    double precision                                      , intent(  out), optional :: gradientLogarithmic
    integer                                               , dimension(4)            :: indexRadius
    double precision                                      , dimension(4)            :: weightRadius       , weightRadiusDerivative
    double precision                                                                :: position
    logical                                                                         :: inRange
    integer                                                                         :: i                  , j

    if (radius > 0.0d0) then
       position=(log(radius)-self%radiusTruncateLogarithmic-self%table%radiusLogarithmicMinimum)/self%table%radiusLogarithmicInterval
    else
       position=-1.0d0
    end if
    inRange =position >= 0.0d0 .and. position <= dble(self%table%countRadius-1)
    position=max(0.0d0,min(position,dble(self%table%countRadius-1)))
    call exponentialDiskTidalWeights(position,self%table%countRadius,indexRadius,weightRadius)
    ratio=0.0d0
    do i=1,4
       do j=1,4
          select case (quantity)
          case (exponentialDiskTidalQuantityRotationCurve )
             ratio=ratio+self%weight(i,j)*sum(weightRadius*self%table%rotationCurveRatioLogarithmic (indexRadius,self%indexTruncation(i),self%indexHeating(j)))
          case (exponentialDiskTidalQuantityMass          )
             ratio=ratio+self%weight(i,j)*sum(weightRadius*self%table%massRatioLogarithmic          (indexRadius,self%indexTruncation(i),self%indexHeating(j)))
          case (exponentialDiskTidalQuantitySurfaceDensity)
             ratio=ratio+self%weight(i,j)*sum(weightRadius*self%table%surfaceDensityRatioLogarithmic(indexRadius,self%indexTruncation(i),self%indexHeating(j)))
          case (exponentialDiskTidalQuantityPotential     )
             ratio=ratio+self%weight(i,j)*sum(weightRadius*self%table%potentialRatioLogarithmic     (indexRadius,self%indexTruncation(i),self%indexHeating(j)))
          end select
       end do
    end do
    ! All tabulated quantities scale linearly with the mass of the disk, so apply any scaling of the mass.
    ratio=exp(max(ratio+self%massScaleLogarithmic,log(tiny(0.0d0))))
    if (present(gradientLogarithmic)) then
       gradientLogarithmic=0.0d0
       if (inRange) then
          call exponentialDiskTidalWeights(position,self%table%countRadius,indexRadius,weightRadiusDerivative,derivative=.true.)
          do i=1,4
             do j=1,4
                select case (quantity)
                case (exponentialDiskTidalQuantityRotationCurve )
                   gradientLogarithmic=gradientLogarithmic+self%weight(i,j)*sum(weightRadiusDerivative*self%table%rotationCurveRatioLogarithmic (indexRadius,self%indexTruncation(i),self%indexHeating(j)))
                case (exponentialDiskTidalQuantityMass          )
                   gradientLogarithmic=gradientLogarithmic+self%weight(i,j)*sum(weightRadiusDerivative*self%table%massRatioLogarithmic          (indexRadius,self%indexTruncation(i),self%indexHeating(j)))
                case (exponentialDiskTidalQuantitySurfaceDensity)
                   gradientLogarithmic=gradientLogarithmic+self%weight(i,j)*sum(weightRadiusDerivative*self%table%surfaceDensityRatioLogarithmic(indexRadius,self%indexTruncation(i),self%indexHeating(j)))
                case (exponentialDiskTidalQuantityPotential     )
                   gradientLogarithmic=gradientLogarithmic+self%weight(i,j)*sum(weightRadiusDerivative*self%table%potentialRatioLogarithmic     (indexRadius,self%indexTruncation(i),self%indexHeating(j)))
                end select
             end do
          end do
          gradientLogarithmic=gradientLogarithmic/self%table%radiusLogarithmicInterval
       end if
    end if
    return
  end function exponentialDiskTidalInterpolate

  logical function exponentialDiskTidalBeyondTable(self,radius) result(beyond)
    !!{RST
    Return true if the given dimensionless ``radius`` (in units of the disk scale length) lies beyond the tabulated range.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius

    beyond=radius > 0.0d0
    if (beyond)                                                                                                            &
         & beyond=+log(radius)                                                                                             &
         &        -self%radiusTruncateLogarithmic                                                                          &
         &        >                                                                                                        &
         &        +self%table%radiusLogarithmicMinimum+self%table%radiusLogarithmicInterval*dble(self%table%countRadius-1)
    return
  end function exponentialDiskTidalBeyondTable

  double precision function exponentialDiskTidalRatioRotationCurve(self,radius) result(ratio)
    !!{RST
    Return the ratio of the square of the mid-plane rotation curve to that of the untruncated, unheated disk at the given
    dimensionless ``radius`` (in units of the disk scale length). Beyond the tabulated range the rotation curve of the
    tidally-modified disk is that of its multipole expansion to quadrupole order, :math:`V_\mathrm{c}^2 = GM/R - 3GQ/2R^3`.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: radius_

    if (exponentialDiskTidalBeyondTable(self,radius)) then
       radius_=radius*self%massDistribution_%radiusScale
       ratio  =+(                                          &
            &    +      self%massTotal_/radius_            &
            &    -1.5d0*self%quadrupole/radius_**3         &
            &   )                                          &
            &  /self%massDistribution_%rotationCurve(radius_)**2
    else
       ratio=self%interpolate(radius,exponentialDiskTidalQuantityRotationCurve)
    end if
    return
  end function exponentialDiskTidalRatioRotationCurve

  logical function exponentialDiskTidalIsDimensionless(self)
    !!{RST
    Return the dimensional status.
    !!}
    implicit none
    class(massDistributionExponentialDiskTidal), intent(inout) :: self
    !$GLC attributes unused :: self

    exponentialDiskTidalIsDimensionless=.true.
    return
  end function exponentialDiskTidalIsDimensionless

  logical function exponentialDiskTidalAssumeMonotonicDecreasingSurfaceDensity(self) result(assumeMonotonicDecreasingSurfaceDensity)
    !!{RST
    Return true if the tabulated surface densities are all monotonically decreasing.
    !!}
    implicit none
    class(massDistributionExponentialDiskTidal), intent(inout) :: self

    assumeMonotonicDecreasingSurfaceDensity=self%table%monotonic
    return
  end function exponentialDiskTidalAssumeMonotonicDecreasingSurfaceDensity

  double precision function exponentialDiskTidalMassTotal(self)
    !!{RST
    Return the total mass of a tidally-modified exponential disk.
    !!}
    implicit none
    class(massDistributionExponentialDiskTidal), intent(inout) :: self

    exponentialDiskTidalMassTotal=self%massTotal_
    return
  end function exponentialDiskTidalMassTotal

  double precision function exponentialDiskTidalSurfaceDensity(self,coordinates)
    !!{RST
    Return the surface density at the specified ``coordinates`` in a tidally-modified exponential disk.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCylindrical
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout) :: self
    class           (coordinate                          ), intent(in   ) :: coordinates
    type            (coordinateCylindrical               )                :: position
    double precision                                                      :: radius

    position=coordinates
    radius  =position%r()/self%massDistribution_%radiusScale
    if (exponentialDiskTidalBeyondTable(self,radius)) then
       exponentialDiskTidalSurfaceDensity=0.0d0
    else
       exponentialDiskTidalSurfaceDensity=+self%interpolate(radius,exponentialDiskTidalQuantitySurfaceDensity) &
            &                             *self%massDistribution_%surfaceDensity(coordinates)
    end if
    return
  end function exponentialDiskTidalSurfaceDensity

  double precision function exponentialDiskTidalVerticalFactor(self,height,gradientLogarithmic) result(factor)
    !!{RST
    Return the factor, :math:`\mathrm{sech}^2(z/h)/2h`, which converts surface density to density at the given ``height``
    for the :math:`\mathrm{sech}^2` vertical profile of the exponential disk. If present, ``gradientLogarithmic`` returns
    :math:`\mathrm{d}\ln \rho/\mathrm{d}\ln z`.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout)           :: self
    double precision                                      , intent(in   )           :: height
    double precision                                      , intent(  out), optional :: gradientLogarithmic
    double precision                                      , parameter               :: coshArgumentMaximum=50.0d0
    double precision                                                                :: z

    z=abs(height)/self%massDistribution_%scaleHeight
    if (z > coshArgumentMaximum) then
       factor=(2.0d0*exp(-z)/(1.0d0+exp(-2.0d0*z)))**2
    else
       factor=1.0d0/cosh(z)**2
    end if
    factor=factor/2.0d0/self%massDistribution_%scaleHeight
    if (present(gradientLogarithmic)) gradientLogarithmic=-2.0d0*z*tanh(z)
    return
  end function exponentialDiskTidalVerticalFactor

  double precision function exponentialDiskTidalDensity(self,coordinates)
    !!{RST
    Return the density at the specified ``coordinates`` in a tidally-modified exponential disk.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCylindrical
    implicit none
    class(massDistributionExponentialDiskTidal), intent(inout) :: self
    class(coordinate                          ), intent(in   ) :: coordinates
    type (coordinateCylindrical               )                :: position

    position                   =coordinates
    exponentialDiskTidalDensity=+self%surfaceDensity               (     coordinates    ) &
         &                      *exponentialDiskTidalVerticalFactor(self,position   %z())
    return
  end function exponentialDiskTidalDensity

  double precision function exponentialDiskTidalDensityGradientRadial(self,coordinates,logarithmic)
    !!{RST
    Return the gradient of the density with respect to spherical radius at the specified ``coordinates`` in a tidally-modified
    exponential disk.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCylindrical
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout), target   :: self
    class           (coordinate                          ), intent(in   )           :: coordinates
    logical                                               , intent(in   ), optional :: logarithmic
    type            (coordinateCylindrical               )                          :: position
    double precision                                                                :: gradientRadial, gradientVertical, &
         &                                                                             density_      , ratio
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    position=coordinates
    ratio   =self%interpolate(position%r()/self%massDistribution_%radiusScale,exponentialDiskTidalQuantitySurfaceDensity,gradientRadial)
    density_=exponentialDiskTidalVerticalFactor(self,position%z(),gradientVertical)
    ! The logarithmic gradient with respect to spherical radius is the sum of those with respect to cylindrical radius and height.
    ! The gradient with respect to cylindrical radius is that of the ratio, plus that of the exponential disk.
    exponentialDiskTidalDensityGradientRadial=+gradientRadial                                  &
         &                                    -position%r()/self%massDistribution_%radiusScale &
         &                                    +gradientVertical
    if (.not.logarithmic_)                                                                                   &
         & exponentialDiskTidalDensityGradientRadial=+exponentialDiskTidalDensityGradientRadial              &
         &                                           *self    %density                         (coordinates) &
         &                                           /position%rSpherical                      (           )
    return
  end function exponentialDiskTidalDensityGradientRadial

  double precision function exponentialDiskTidalDensitySlopeLogarithmicCentral(self) result(slope)
    !!{RST
    Return the central logarithmic slope of the density profile in a tidally-modified exponential disk. Neither truncation nor
    heating changes the disk at its center, so this is the slope of the exponential disk.
    !!}
    implicit none
    class(massDistributionExponentialDiskTidal), intent(inout) :: self

    slope=self%massDistribution_%densitySlopeLogarithmicCentral()
    return
  end function exponentialDiskTidalDensitySlopeLogarithmicCentral

  double precision function exponentialDiskTidalDensitySphericalAverage(self,radius)
    !!{RST
    Return the spherically-averaged density at the specified ``radius`` in a tidally-modified exponential disk, using the
    thin-disk approximation.
    !!}
    use :: Coordinates, only : coordinateCylindrical
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius

    exponentialDiskTidalDensitySphericalAverage=+0.5d0                                                          &
         &                                      *self%surfaceDensity(coordinateCylindrical(radius,0.0d0,0.0d0)) &
         &                                      /     radius
    return
  end function exponentialDiskTidalDensitySphericalAverage

  double precision function exponentialDiskTidalRadiusHalfMass(self)
    !!{RST
    Return the half-mass radius of a tidally-modified exponential disk.
    !!}
    implicit none
    class(massDistributionExponentialDiskTidal), intent(inout) :: self

    if (self%radiusHalfMass_ < 0.0d0) self%radiusHalfMass_=self%radiusCylindricalEnclosingMass(massFractional=0.5d0)
    exponentialDiskTidalRadiusHalfMass=self%radiusHalfMass_
    return
  end function exponentialDiskTidalRadiusHalfMass

  double precision function exponentialDiskTidalMassEnclosedByCylinder(self,radius) result(mass)
    !!{RST
    Return the mass enclosed within a cylinder of the given ``radius`` in a tidally-modified exponential disk.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout), target :: self
    double precision                                      , intent(in   )         :: radius

    if (radius <= 0.0d0) then
       mass=0.0d0
    else if (exponentialDiskTidalBeyondTable(self,radius/self%massDistribution_%radiusScale)) then
       mass   = self                  %massTotal_
    else
       mass   =+self                  %interpolate           (radius/self%massDistribution_%radiusScale,exponentialDiskTidalQuantityMass) &
            &  *self%massDistribution_%massEnclosedByCylinder(radius                                                                    )
    end if
    return
  end function exponentialDiskTidalMassEnclosedByCylinder

  double precision function exponentialDiskTidalMassEnclosedBySphere(self,radius)
    !!{RST
    Return the mass enclosed within a sphere of the given ``radius`` in a tidally-modified exponential disk, using the thin-disk
    approximation.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout), target :: self
    double precision                                      , intent(in   )         :: radius

    exponentialDiskTidalMassEnclosedBySphere=self%massEnclosedByCylinder(radius)
    return
  end function exponentialDiskTidalMassEnclosedBySphere

  logical function exponentialDiskTidalPotentialIsAnalytic(self) result(isAnalytic)
    !!{RST
    Return false, as the potential of a tidally-modified exponential disk is not analytic.
    !!}
    implicit none
    class(massDistributionExponentialDiskTidal), intent(inout) :: self
    !$GLC attributes unused :: self

    isAnalytic=.false.
    return
  end function exponentialDiskTidalPotentialIsAnalytic

  double precision function exponentialDiskTidalRotationCurve(self,radius)
    !!{RST
    Return the mid-plane rotation curve of a tidally-modified exponential disk.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius

    if (radius <= 0.0d0) then
       exponentialDiskTidalRotationCurve=0.0d0
    else
       exponentialDiskTidalRotationCurve=+sqrt(self                  %ratioRotationCurve(radius/self%massDistribution_%radiusScale)) &
            &                            *     self%massDistribution_%rotationCurve     (radius                                    )
    end if
    return
  end function exponentialDiskTidalRotationCurve

  double precision function exponentialDiskTidalRotationCurveGradient(self,radius)
    !!{RST
    Return the gradient, :math:`\mathrm{d}V_\mathrm{c}^2/\mathrm{d}R`, of the mid-plane rotation curve of a tidally-modified
    exponential disk.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: ratio  , gradientLogarithmic, &
         &                                                                   radius_

    if (radius <= 0.0d0) then
       exponentialDiskTidalRotationCurveGradient=0.0d0
       return
    end if
    radius_=radius/self%massDistribution_%radiusScale
    if (exponentialDiskTidalBeyondTable(self,radius_)) then
       ! Beyond the tabulated range use the multipole expansion to quadrupole order.
       exponentialDiskTidalRotationCurveGradient=-      self%massTotal_/radius**2 &
            &                                    +4.5d0*self%quadrupole/radius**4
    else
       ratio                                    =+self%interpolate(radius_,exponentialDiskTidalQuantityRotationCurve,gradientLogarithmic)
       exponentialDiskTidalRotationCurveGradient=+ratio                                                     &
            &                                    *(                                                         &
            &                                      +self%massDistribution_%rotationCurveGradient(radius)    &
            &                                      +self%massDistribution_%rotationCurve        (radius)**2 &
            &                                      *gradientLogarithmic                                     &
            &                                      /radius                                                  &
            &                                     )
    end if
    return
  end function exponentialDiskTidalRotationCurveGradient

  subroutine exponentialDiskTidalMultipole(self,positionCartesian,potential,acceleration,tidalTensor)
    !!{RST
    Evaluate the potential, acceleration, and tidal tensor of the disk at ``positionCartesian`` from its multipole expansion to
    quadrupole order, :math:`\Phi = -GM/r - G Q f`, where :math:`f = (3z^2-r^2)/2r^5` and :math:`Q=M(\langle z^2 \rangle -
    \langle R^2 \rangle/2)` is the quadrupole moment of the disk. This is used beyond the tabulated range of radii.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout)                           :: self
    double precision                                      , intent(in   ), dimension(3  )           :: positionCartesian
    double precision                                      , intent(  out)                , optional :: potential
    double precision                                      , intent(  out), dimension(3  ), optional :: acceleration
    double precision                                      , intent(  out), dimension(3,3), optional :: tidalTensor
    double precision                                                                                :: radius           , height
    integer                                                                                         :: i                , j

    radius=sqrt(sum(positionCartesian**2))
    height=positionCartesian(3)
    if (present(potential   )) potential   =-self%massTotal_                                   &
         &                                  /radius                                            &
         &                                  -self%quadrupole                                   &
         &                                  *(3.0d0*height**2-radius**2)/2.0d0/radius**5
    if (present(acceleration)) then
       acceleration   =-self%massTotal_*positionCartesian/radius**3                            &
            &          +self%quadrupole                                                        &
            &          *(                                                                      &
            &            +1.5d0*positionCartesian/radius**5                                    &
            &            -7.5d0*height**2*positionCartesian/radius**7                          &
            &           )
       acceleration(3)=+acceleration(3)                                                        &
            &          +self%quadrupole*3.0d0*height/radius**5
    end if
    if (present(tidalTensor )) then
       do i=1,3
          do j=1,3
             tidalTensor(i,j)=+self%massTotal_*3.0d0*positionCartesian(i)*positionCartesian(j)/radius**5 &
                  &           +self%quadrupole                                                           &
                  &           *(                                                                         &
                  &             +52.5d0*height**2*positionCartesian(i)*positionCartesian(j)/radius**9    &
                  &             - 7.5d0          *positionCartesian(i)*positionCartesian(j)/radius**7    &
                  &            )
          end do
       end do
       do i=1,3
          tidalTensor(i,i)=+tidalTensor(i,i)            &
               &           -self%massTotal_/radius**3   &
               &           +self%quadrupole             &
               &           *(                           &
               &             +1.5d0          /radius**5 &
               &             -7.5d0*height**2/radius**7 &
               &            )
          tidalTensor(i,3)=tidalTensor(i,3)-self%quadrupole*15.0d0*height*positionCartesian(i)/radius**7
          tidalTensor(3,i)=tidalTensor(3,i)-self%quadrupole*15.0d0*height*positionCartesian(i)/radius**7
       end do
       tidalTensor(3,3)=tidalTensor(3,3)+self%quadrupole*3.0d0/radius**5
       ! Products such as xᵢxⱼ and xⱼxᵢ may differ by rounding, so copy the upper triangle to the lower to make the tensor exactly
       ! symmetric.
       do i=2,3
          do j=1,i-1
             tidalTensor(i,j)=tidalTensor(j,i)
          end do
       end do
    end if
    return
  end subroutine exponentialDiskTidalMultipole

  subroutine exponentialDiskTidalMultipoleNormalized(self,positionCartesian,potential,acceleration,tidalTensor)
    !!{RST
    Evaluate the potential, acceleration, and tidal tensor of the disk at ``positionCartesian`` (which must lie within the
    tabulated range of radii) from its multipole expansion to quadrupole order, normalized such that, at the same spherical
    radius, the mid-plane potential and rotation curve match those tabulated. The potential is multiplied by :math:`p_\Phi(r) =
    \Phi(r,0)/\Phi_\mathrm{m}(r,0)`, and the acceleration by :math:`p(r) = V_\mathrm{c}^2(r)/V_\mathrm{c,m}^2(r)`, where the
    subscript "m" indicates the multipole expansion. The tidal tensor is the (symmetrized) gradient of :math:`p(r)
    \mathbf{a}_\mathrm{m}(\mathbf{x})`. In the mid-plane these are therefore exact, while far from the disk they approach the
    multipole expansion.
    !!}
    use :: Coordinates, only : coordinateCylindrical
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout)                           :: self
    double precision                                      , intent(in   ), dimension(3  )           :: positionCartesian
    double precision                                      , intent(  out)                , optional :: potential
    double precision                                      , intent(  out), dimension(3  ), optional :: acceleration
    double precision                                      , intent(  out), dimension(3,3), optional :: tidalTensor
    double precision                                                     , dimension(3  )           :: acceleration_
    double precision                                                                                :: radius                       , ratio               , &
         &                                                                                             rotationCurveSquaredMultipole, rotationCurveSquared, &
         &                                                                                             gradientLogarithmic
    integer                                                                                         :: i                            , j

    radius=sqrt(sum(positionCartesian**2))
    if (present(potential)) then
       call exponentialDiskTidalMultipole(self,positionCartesian,potential=potential)
       potential=+potential                                                                                                                 &
            &    *exponentialDiskTidalRatioAtPosition(self,coordinateCylindrical(radius,0.0d0,0.0d0),exponentialDiskTidalQuantityPotential) &
            &    *self%massDistribution_%potential   (     coordinateCylindrical(radius,0.0d0,0.0d0)                                      ) &
            &    /(                                                                                                                         &
            &      -      self%massTotal_/radius                                                                                            &
            &      +0.5d0*self%quadrupole/radius**3                                                                                         &
            &     )
    end if
    if (present(acceleration) .or. present(tidalTensor)) then
       rotationCurveSquaredMultipole=+      self%massTotal_/radius    &
            &                        -1.5d0*self%quadrupole/radius**3
       rotationCurveSquared         =+self%rotationCurve(radius)**2
       ratio                        =+rotationCurveSquared            &
            &                        /rotationCurveSquaredMultipole
       call exponentialDiskTidalMultipole(self,positionCartesian,acceleration=acceleration_,tidalTensor=tidalTensor)
       if (present(acceleration)) acceleration=ratio*acceleration_
       if (present(tidalTensor )) then
          ! The logarithmic gradient of the ratio is that of the tabulated V², less that of the multipole V².
          gradientLogarithmic=+radius*self%rotationCurveGradient(radius)/rotationCurveSquared                                    &
               &              -radius*(-self%massTotal_/radius**2+4.5d0*self%quadrupole/radius**4)/rotationCurveSquaredMultipole
          do i=1,3
             do j=1,3
                tidalTensor(i,j)=+ratio                                                                           &
                     &           *(                                                                               &
                     &             +tidalTensor(i,j)                                                              &
                     &             +gradientLogarithmic                                                           &
                     &             *(acceleration_(i)*positionCartesian(j)+acceleration_(j)*positionCartesian(i)) &
                     &             /2.0d0                                                                         &
                     &             /radius**2                                                                     &
                     &            )
             end do
          end do
       end if
    end if
    return
  end subroutine exponentialDiskTidalMultipoleNormalized

  subroutine exponentialDiskTidalBlend(self,radius,weight,weightDerivative)
    !!{RST
    Return the ``weight`` given to the normalized multipole expansion (as opposed to the scaled field of the untruncated, unheated
    disk) away from the mid-plane at spherical ``radius``, and, optionally, its derivative with respect to radius. The weight is
    :math:`w = w_r(r) g(\sigma)`. Here :math:`w_r` rises smoothly (as :math:`3s^2-2s^3`, with :math:`s` linear in :math:`\ln r`)
    from zero at ``radiusBlendMinimum`` to unity at ``radiusBlendMaximum``. The multipole expansion cannot describe the field
    close to material at large radii (such as the exponential tail of a weakly truncated disk), while the scaled field of the
    untruncated disk does, so :math:`g = (1-\min[\sigma,1])^2`, where :math:`\sigma(r)` is the ratio of the surface density to
    that of the untruncated, unheated disk at radius :math:`r`, limits the weight to regions from which material has been
    removed.
    !!}
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout)           :: self
    double precision                                      , intent(in   )           :: radius
    double precision                                      , intent(  out)           :: weight
    double precision                                      , intent(  out), optional :: weightDerivative
    double precision                                                                :: position              , weightRadial, &
         &                                                                             weightRadialDerivative, ratio       , &
         &                                                                             gradientLogarithmic   , factor      , &
         &                                                                             factorDerivative

    if (radius <= self%radiusBlendMinimum) then
       weight=0.0d0
       if (present(weightDerivative)) weightDerivative=0.0d0
       return
    else if (radius >= self%radiusBlendMaximum) then
       weightRadial          =1.0d0
       weightRadialDerivative=0.0d0
    else
       position              =+log(radius                 /self%radiusBlendMinimum) &
            &                 /log(self%radiusBlendMaximum/self%radiusBlendMinimum)
       weightRadial          =+position**2                                          &
            &                 *(3.0d0-2.0d0*position)
       weightRadialDerivative=+6.0d0                                                &
            &                 *position                                             &
            &                 *(1.0d0-position)                                     &
            &                 /radius                                               &
            &                 /log(self%radiusBlendMaximum/self%radiusBlendMinimum)
    end if
    ratio=self%interpolate(radius/self%massDistribution_%radiusScale,exponentialDiskTidalQuantitySurfaceDensity,gradientLogarithmic)
    if (ratio < 1.0d0) then
       factor          =+(1.0d0-ratio)**2
       ! d g/d r = -2 (1-σ) σ (d ln σ/d ln r)/r.
       factorDerivative=-2.0d0               &
            &           *(1.0d0-ratio)       &
            &           *ratio               &
            &           *gradientLogarithmic &
            &           /radius
    else
       factor          =0.0d0
       factorDerivative=0.0d0
    end if
    weight=weightRadial*factor
    if (present(weightDerivative))                  &
         & weightDerivative=+weightRadialDerivative &
         &                  *factor                 &
         &                  +weightRadial           &
         &                  *factorDerivative
    return
  end subroutine exponentialDiskTidalBlend

  double precision function exponentialDiskTidalRatioAtPosition(self,coordinates,quantity,gradientLogarithmic) result(ratio)
    !!{RST
    Return the ratio of the tabulated ``quantity`` to that of the untruncated, unheated disk at the spherical radius of the given
    ``coordinates`` (which must lie within the tabulated range). This ratio is used to approximate quantities away from the
    mid-plane. If present, ``gradientLogarithmic`` returns the derivative of the logarithm of the ratio with respect to the
    logarithm of the spherical radius.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout)           :: self
    class           (coordinate                          ), intent(in   )           :: coordinates
    integer                                               , intent(in   )           :: quantity
    double precision                                      , intent(  out), optional :: gradientLogarithmic
    type            (coordinateSpherical                 )                          :: position
    double precision                                                                :: radius

    position=coordinates
    radius  =position%r()/self%massDistribution_%radiusScale
    if (radius <= 0.0d0) radius=exp(self%radiusTruncateLogarithmic+self%table%radiusLogarithmicMinimum)
    ratio   =self%interpolate(radius,quantity,gradientLogarithmic)
    return
  end function exponentialDiskTidalRatioAtPosition

  logical function exponentialDiskTidalBeyondTableAtPosition(self,coordinates) result(beyond)
    !!{RST
    Return true if the spherical radius of the given ``coordinates`` lies beyond the tabulated range, where the multipole
    expansion of the potential is used.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class(massDistributionExponentialDiskTidal), intent(inout) :: self
    class(coordinate                          ), intent(in   ) :: coordinates
    type (coordinateSpherical                 )                :: position

    position=coordinates
    beyond  =position%r() > 0.0d0
    if (beyond) beyond=exponentialDiskTidalBeyondTable(self,position%r()/self%massDistribution_%radiusScale)
    return
  end function exponentialDiskTidalBeyondTableAtPosition

  double precision function exponentialDiskTidalPotential(self,coordinates,status)
    !!{RST
    Return the gravitational potential of a tidally-modified exponential disk. Close to the disk this is the potential of the
    untruncated, unheated disk multiplied by the ratio of the mid-plane potentials at the spherical radius of the point. Further
    out it is blended smoothly into the normalized multipole expansion (see
    ``exponentialDiskTidalMultipoleNormalized``). Both are exact in the mid-plane. Beyond the tabulated range of radii the
    multipole expansion of the potential to quadrupole order is used.
    !!}
    use :: Coordinates               , only : assignment(=)            , coordinateCartesian
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout), target   :: self
    class           (coordinate                          ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType   ), intent(  out), optional :: status
    type            (coordinateCartesian                 )                          :: position
    double precision                                      , dimension(3)            :: positionCartesian
    double precision                                                                :: weight           , potentialMultipole

    if (present(status)) status=structureErrorCodeSuccess
    position         =coordinates
    positionCartesian=position
    if (exponentialDiskTidalBeyondTableAtPosition(self,coordinates)) then
       call exponentialDiskTidalMultipole(self,positionCartesian,potential=exponentialDiskTidalPotential)
    else
       call exponentialDiskTidalBlend(self,sqrt(sum(positionCartesian**2)),weight)
       exponentialDiskTidalPotential=0.0d0
       if (weight < 1.0d0)                                                                                                               &
            & exponentialDiskTidalPotential=+(1.0d0-weight)                                                                              &
            &                               *exponentialDiskTidalRatioAtPosition(self,coordinates,exponentialDiskTidalQuantityPotential) &
            &                               *self%massDistribution_%potential   (     coordinates,status                               )
       if (weight > 0.0d0) then
          call exponentialDiskTidalMultipoleNormalized(self,positionCartesian,potential=potentialMultipole)
          exponentialDiskTidalPotential=+exponentialDiskTidalPotential &
               &                        +weight                        &
               &                        *potentialMultipole
       end if
    end if
    return
  end function exponentialDiskTidalPotential

  function exponentialDiskTidalAcceleration(self,coordinates)
    !!{RST
    Return the gravitational acceleration due to a tidally-modified exponential disk. Close to the disk this is the acceleration
    due to the untruncated, unheated disk multiplied by the ratio of the squares of the mid-plane rotation curves at the
    spherical radius of the point. Further out it is blended smoothly into the normalized multipole expansion (see
    ``exponentialDiskTidalMultipoleNormalized``). Both are exact in the mid-plane. Beyond the tabulated range of radii the
    multipole expansion of the potential to quadrupole order is used.
    !!}
    use :: Coordinates, only : coordinateCartesian, assignment(=)
    implicit none
    type            (coordinateCartesian                 )                :: exponentialDiskTidalAcceleration
    class           (massDistributionExponentialDiskTidal), intent(inout) :: self
    class           (coordinate                          ), intent(in   ) :: coordinates
    type            (coordinateCartesian                 )                :: position
    double precision                                      , dimension(3)  :: acceleration                    , positionCartesian, &
         &                                                                   accelerationMultipole
    double precision                                                      :: weight

    position         = coordinates
    positionCartesian= position
    if (exponentialDiskTidalBeyondTableAtPosition(self,coordinates)) then
       call exponentialDiskTidalMultipole(self,positionCartesian,acceleration=acceleration)
    else
       call exponentialDiskTidalBlend(self,sqrt(sum(positionCartesian**2)),weight)
       acceleration=0.0d0
       if (weight < 1.0d0) then
          acceleration=self%massDistribution_%acceleration(coordinates)
          acceleration=+(1.0d0-weight)                                                                                  &
               &       *exponentialDiskTidalRatioAtPosition(self,coordinates,exponentialDiskTidalQuantityRotationCurve) &
               &       *acceleration
       end if
       if (weight > 0.0d0) then
          call exponentialDiskTidalMultipoleNormalized(self,positionCartesian,acceleration=accelerationMultipole)
          acceleration=+acceleration          &
               &       +weight                &
               &       *accelerationMultipole
       end if
    end if
    exponentialDiskTidalAcceleration=acceleration
    return
  end function exponentialDiskTidalAcceleration

  function exponentialDiskTidalTidalTensor(self,coordinates)
    !!{RST
    Return the gravitational tidal tensor due to a tidally-modified exponential disk. Close to the disk, the acceleration is
    approximated as :math:`\mathbf{a}_\mathrm{s}(\mathbf{x}) = q(r) \mathbf{a}_0(\mathbf{x})`, where :math:`\mathbf{a}_0` is
    the acceleration due to the untruncated, unheated disk, and :math:`q(r)` is the ratio of the squares of the mid-plane rotation
    curves at spherical radius :math:`r`. Its tidal tensor is the symmetrized gradient of this acceleration,

    .. math::

       T_{\mathrm{s},ij} = q T_{0,ij} + \frac{q^\prime}{2r} (a_{0,i} x_j + a_{0,j} x_i).

    Further out this is blended into the tidal tensor, :math:`T_{\mathrm{m},ij}`, of the normalized multipole expansion,
    :math:`\mathbf{a}_\mathrm{m}` (see ``exponentialDiskTidalMultipoleNormalized``), with weight :math:`w(r)`, as the
    symmetrized gradient of :math:`(1-w)\mathbf{a}_\mathrm{s}+w\mathbf{a}_\mathrm{m}`:

    .. math::

       T_{ij} = (1-w) T_{\mathrm{s},ij} + w T_{\mathrm{m},ij} + \frac{w^\prime}{2r} (\Delta_i x_j + \Delta_j x_i),

    where :math:`\boldsymbol{\Delta}=\mathbf{a}_\mathrm{m}-\mathbf{a}_\mathrm{s}` (which vanishes in the mid-plane).
    :math:`T_{zz}` is then adjusted such that :math:`\mathrm{Tr}(T) = -4\pi G\rho` (Poisson's equation). In the mid-plane this
    is exact. Beyond the tabulated range of radii the multipole expansion of the potential to quadrupole order is used.
    !!}
    use :: Coordinates             , only : assignment(=), coordinateCartesian
    use :: Numerical_Constants_Math, only : Pi
    use :: Tensors                 , only : assignment(=)
    implicit none
    type            (tensorRank2Dimension3Symmetric      )                  :: exponentialDiskTidalTidalTensor
    class           (massDistributionExponentialDiskTidal), intent(inout)   :: self
    class           (coordinate                          ), intent(in   )   :: coordinates
    type            (coordinateCartesian                 )                  :: position
    double precision                                      , dimension(3  )  :: positionCartesian              , acceleration          , &
         &                                                                     accelerationMultipole          , accelerationDifference
    double precision                                      , dimension(3,3)  :: tidalTensor                    , tidalTensorMultipole
    double precision                                                        :: radius                         , ratio                 , &
         &                                                                     gradientLogarithmic            , weight                , &
         &                                                                     weightDerivative
    integer                                                                 :: i                              , j

    position         =coordinates
    positionCartesian=position
    radius           =sqrt(sum(positionCartesian**2))
    if (exponentialDiskTidalBeyondTableAtPosition(self,coordinates)) then
       call exponentialDiskTidalMultipole(self,positionCartesian,tidalTensor=tidalTensor)
    else
       call exponentialDiskTidalBlend(self,radius,weight,weightDerivative)
       tidalTensor =0.0d0
       acceleration=0.0d0
       if (weight < 1.0d0) then
          ! The scaled field of the untruncated, unheated disk.
          ratio      =exponentialDiskTidalRatioAtPosition(self,coordinates,exponentialDiskTidalQuantityRotationCurve,gradientLogarithmic)
          tidalTensor=self%massDistribution_%tidalTensor (coordinates)
          tidalTensor=ratio*tidalTensor
          if (radius > 0.0d0) then
             acceleration=self%massDistribution_%acceleration(coordinates)
             ! The gradient of the ratio is ∇q = q (d ln q/d ln r) x/r².
             do i=1,3
                do j=1,3
                   tidalTensor(i,j)=+tidalTensor(i,j)                                                            &
                        &           +ratio                                                                       &
                        &           *gradientLogarithmic                                                         &
                        &           *(acceleration(i)*positionCartesian(j)+acceleration(j)*positionCartesian(i)) &
                        &           /2.0d0                                                                       &
                        &           /radius**2
                end do
             end do
             acceleration=ratio*acceleration
          end if
       end if
       if (weight > 0.0d0) then
          ! Blend in the normalized multipole expansion.
          call exponentialDiskTidalMultipoleNormalized(self,positionCartesian,acceleration=accelerationMultipole,tidalTensor=tidalTensorMultipole)
          accelerationDifference=accelerationMultipole-acceleration
          do i=1,3
             do j=1,3
                tidalTensor(i,j)=+(1.0d0-weight)*tidalTensor         (i,j)                                                        &
                     &           +       weight *tidalTensorMultipole(i,j)                                                        &
                     &           +weightDerivative                                                                                &
                     &           *(accelerationDifference(i)*positionCartesian(j)+accelerationDifference(j)*positionCartesian(i)) &
                     &           /2.0d0                                                                                           &
                     &           /radius
             end do
          end do
       end if
       ! Adjust the vertical-vertical component to satisfy Poisson's equation.
       tidalTensor(3,3)=+tidalTensor(3,3)                   &
            &           -4.0d0*Pi*self%density(coordinates) &
            &           -(                                  &
            &             +tidalTensor(1,1)                 &
            &             +tidalTensor(2,2)                 &
            &             +tidalTensor(3,3)                 &
            &            )
    end if
    exponentialDiskTidalTidalTensor=tidalTensor
    return
  end function exponentialDiskTidalTidalTensor

  double precision function exponentialDiskTidalSurfaceDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{RST
    Return the radial moment, :math:`\int \Sigma(R) R^m \mathrm{d}R`, of the surface density of a tidally-modified exponential
    disk. The first and second moments over all radii are interpolated from tabulated values. Other moments are computed by
    numerical integration.
    !!}
    use :: Coordinates          , only : coordinateCylindrical
    use :: Error                , only : Error_Report         , errorStatusSuccess, errorStatusRound
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (massDistributionExponentialDiskTidal), intent(inout)           :: self
    double precision                                      , intent(in   )           :: moment
    double precision                                      , intent(in   ), optional :: radiusMinimum           , radiusMaximum
    logical                                               , intent(  out), optional :: isInfinite
    type            (integrator                          )                          :: integrator_
    double precision                                                                :: radiusLogarithmicMinimum, radiusLogarithmicMaximum, &
         &                                                                             radiusLogarithmicInner  , radiusLogarithmicUpper  , &
         &                                                                             radiusInnerMinimum
    logical                                                                         :: includeInner
    integer                                                                         :: i                       , j                       , &
         &                                                                             status

    if (present(isInfinite)) isInfinite=.false.
    if (.not.present(radiusMinimum) .and. .not.present(radiusMaximum) .and. (moment == 1.0d0 .or. moment == 2.0d0)) then
       exponentialDiskTidalSurfaceDensityRadialMoment=0.0d0
       do i=1,4
          do j=1,4
             if (moment == 1.0d0) then
                exponentialDiskTidalSurfaceDensityRadialMoment=exponentialDiskTidalSurfaceDensityRadialMoment+self%weight(i,j)*self%table%moment1(self%indexTruncation(i),self%indexHeating(j))
             else
                exponentialDiskTidalSurfaceDensityRadialMoment=exponentialDiskTidalSurfaceDensityRadialMoment+self%weight(i,j)*self%table%moment2(self%indexTruncation(i),self%indexHeating(j))
             end if
          end do
       end do
    else
       ! Integrate in the logarithm of radius over the requested range, limited to the tabulated range.
       radiusLogarithmicMinimum=self%radiusTruncateLogarithmic+self%table%radiusLogarithmicMinimum
       radiusLogarithmicMaximum=self%radiusTruncateLogarithmic+self%table%radiusLogarithmicMinimum+self%table%radiusLogarithmicInterval*dble(self%table%countRadius-1)
       if (present(radiusMinimum)) then
          if (radiusMinimum > 0.0d0) radiusLogarithmicMinimum=max(radiusLogarithmicMinimum,log(radiusMinimum/self%massDistribution_%radiusScale))
       end if
       if (present(radiusMaximum)) radiusLogarithmicMaximum=min(radiusLogarithmicMaximum,log(radiusMaximum/self%massDistribution_%radiusScale))
       if (radiusLogarithmicMaximum <= radiusLogarithmicMinimum) then
          exponentialDiskTidalSurfaceDensityRadialMoment=0.0d0
       else
          integrator_                                   =integrator(momentIntegrand,toleranceRelative=1.0d-6)
          exponentialDiskTidalSurfaceDensityRadialMoment=integrator_%integrate(radiusLogarithmicMinimum,radiusLogarithmicMaximum,status=status)
          ! The interpolated surface density has discontinuous derivatives at the tabulated radii, so the integration can be limited
          ! by round-off error before the required tolerance is reached. The result is then as accurate as round-off allows, so is
          ! accepted.
          if (status /= errorStatusSuccess .and. status /= errorStatusRound) call Error_Report('integration of the radial moment failed'//{introspection:location})
       end if
       ! Below the tabulated range the ratio of surface densities is constant. Add the contribution from that region if it is
       ! included in the requested range.
       radiusLogarithmicInner=self%radiusTruncateLogarithmic+self%table%radiusLogarithmicMinimum
       includeInner          =.true.
       if (present(radiusMinimum)) includeInner=radiusMinimum <= 0.0d0 .or. log(radiusMinimum/self%massDistribution_%radiusScale) < radiusLogarithmicInner
       if (includeInner) then
          radiusLogarithmicUpper=radiusLogarithmicInner
          if (present(radiusMaximum)) radiusLogarithmicUpper=min(radiusLogarithmicUpper,log(radiusMaximum/self%massDistribution_%radiusScale))
          radiusInnerMinimum=0.0d0
          if (present(radiusMinimum)) radiusInnerMinimum=max(radiusMinimum,0.0d0)
          if (exp(radiusLogarithmicUpper)*self%massDistribution_%radiusScale > radiusInnerMinimum)                                                                                 &
               & exponentialDiskTidalSurfaceDensityRadialMoment=+exponentialDiskTidalSurfaceDensityRadialMoment                                                                    &
               &                                                +self%interpolate(exp(radiusLogarithmicInner),exponentialDiskTidalQuantitySurfaceDensity)                          &
               &                                                *self%massDistribution_%surfaceDensityRadialMoment(                                                                &
               &                                                                                                   moment                                                        , &
               &                                                                                                       radiusInnerMinimum                                        , &
               &                                                                                                   exp(radiusLogarithmicUpper)*self%massDistribution_%radiusScale  &
               &                                                                                                  )
       end if
       return
    end if
    ! Scale to the mass and scale length of the disk.
    exponentialDiskTidalSurfaceDensityRadialMoment=+exponentialDiskTidalSurfaceDensityRadialMoment     &
         &                                         *self%massTotal_                                    &
         &                                         *self%massDistribution_%radiusScale**(moment-1.0d0)
    return

  contains

    double precision function momentIntegrand(radiusLogarithmic)
      !!{RST
      Integrand for radial moments of the surface density, with respect to the logarithm of radius.
      !!}
      implicit none
      double precision, intent(in   ) :: radiusLogarithmic

      momentIntegrand=+self%surfaceDensity(coordinateCylindrical(exp(radiusLogarithmic)*self%massDistribution_%radiusScale,0.0d0,0.0d0)) &
           &          *(exp(radiusLogarithmic)*self%massDistribution_%radiusScale)**(moment+1.0d0)
      return
    end function momentIntegrand

  end function exponentialDiskTidalSurfaceDensityRadialMoment

  function exponentialDiskTidalPositionSample(self,randomNumberGenerator_)
    !!{RST
    Sample a position from a tidally-modified exponential disk. The height and azimuth are sampled from the untruncated,
    unheated disk, and the radius by inverting the tabulated enclosed mass.
    !!}
    use :: Coordinates, only : coordinateCartesian, assignment(=)
    implicit none
    type            (coordinateCartesian                 )                :: exponentialDiskTidalPositionSample
    class           (massDistributionExponentialDiskTidal), intent(inout) :: self
    class           (randomNumberGeneratorClass          ), intent(inout) :: randomNumberGenerator_
    double precision                                      , dimension(3)  :: position
    integer                                               , parameter     :: countBisection                    =60
    double precision                                                      :: radiusLogarithmicLower               , radiusLogarithmicUpper, &
         &                                                                   radiusLogarithmic                    , fraction              , &
         &                                                                   radiusParent
    integer                                                               :: i

    position=self%massDistribution_%positionSample(randomNumberGenerator_)
    fraction=randomNumberGenerator_%uniformSample()
    ! Bisect for the radius enclosing the sampled fraction of the mass.
    radiusLogarithmicLower=self%radiusTruncateLogarithmic+self%table%radiusLogarithmicMinimum-log(1.0d6)
    radiusLogarithmicUpper=self%radiusTruncateLogarithmic+self%table%radiusLogarithmicMinimum+self%table%radiusLogarithmicInterval*dble(self%table%countRadius-1)
    do i=1,countBisection
       radiusLogarithmic=0.5d0*(radiusLogarithmicLower+radiusLogarithmicUpper)
       if (self%massEnclosedByCylinder(exp(radiusLogarithmic)*self%massDistribution_%radiusScale) < fraction*self%massTotal_) then
          radiusLogarithmicLower=radiusLogarithmic
       else
          radiusLogarithmicUpper=radiusLogarithmic
       end if
    end do
    ! Replace the cylindrical radius of the sampled position, retaining its azimuth and height.
    radiusParent=sqrt(position(1)**2+position(2)**2)
    if (radiusParent > 0.0d0) then
       position(1:2)=position(1:2)*exp(radiusLogarithmic)*self%massDistribution_%radiusScale/radiusParent
    else
       position(1:2)=[exp(radiusLogarithmic)*self%massDistribution_%radiusScale,0.0d0]
    end if
    exponentialDiskTidalPositionSample=position
    return
  end function exponentialDiskTidalPositionSample
