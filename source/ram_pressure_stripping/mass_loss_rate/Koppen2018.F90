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
  Implementation of a ram pressure stripping class for disks following :cite:t:`koppen_ram_2018`.
  !!}

  use :: Galactic_Component_Truncations , only : galacticComponentTruncationClass
  use :: Kind_Numbers                   , only : kind_int8
  use :: Hot_Halo_Ram_Pressure_Histories, only : hotHaloRamPressureHistoryClass
  use :: Mass_Distributions             , only : massDistributionClass
  use :: Root_Finder                    , only : rootFinder

  !![
  <ramPressureStripping name="ramPressureStrippingKoppen2018" docformat="rst">
   <description>
   A ram pressure stripping class for disks which follows :cite:t:`koppen_ram_2018`, in which the outcome of stripping depends on
   both the maximum ram pressure, :math:`P`, experienced by the satellite and its time integral, :math:`I` (both supplied by the
   selected :galacticus-class:`hotHaloRamPressureHistoryClass`). Gas at radius :math:`r` in the disk is stripped by a long pulse
   of ram pressure if :math:`P` exceeds the maximum vertical restoring force per unit area,

   .. math::

      p_0(r) = f(r) \Sigma_\mathrm{gas}(r) \left| \frac{\partial \Phi}{\partial z} \right|_\mathrm{max},

   and by a short pulse if :math:`I` exceeds the momentum per unit area needed for the gas to escape,

   .. math::

      I_0(r) = \Sigma_\mathrm{gas}(r) \sqrt{-2\Phi(r,0)-V_\mathrm{c}^2(r)},

   where :math:`\Phi` is the gravitational potential of the satellite, :math:`V_\mathrm{c}` its rotation curve, and the maximum of
   the vertical acceleration is found numerically. These limits are joined using the interpolation of :cite:t:`koppen_ram_2018`
   (their eqn. 11), :math:`P \ge p_0 [1 + 0.5/(I/I_0-1)]`, i.e. gas is stripped if :math:`I &gt; I_0` and :math:`P (I-I_0) \ge p_0
   (I-I_0/2)`. The factor :math:`f(r) = 1 + a \exp(-r/R_0)` (their eqn. 21), with :math:`a=`\ ``[enhancementMolecular]`` and
   :math:`R_0=`\ ``[radiusEnhancementMolecular]`` (in units of the disk scale length), accounts for molecular gas being harder
   to strip than atomic gas, when :math:`\Sigma_\mathrm{gas}` refers to atomic gas only. Since the gas surface density of the
   disk includes all gas, the default is :math:`a=0`.

   The stripping radius, :math:`r_\mathrm{s}`, beyond which all gas is stripped, is found as the root of this criterion. The gas of
   the disk (which must be truncated, see :galacticus-class:`galacticComponentTruncationClass`) then relaxes toward a truncation at
   this radius: if :math:`x_\mathrm{s}=-\log f_\mathrm{s}`, where :math:`f_\mathrm{s}` is the fraction of its mass retained by the
   untruncated gas distribution truncated at :math:`r_\mathrm{s}`, exceeds the current degree of truncation of the gas,
   :math:`x_\mathrm{gas}`, the mass loss rate is

   .. math::

      \dot{M}_\mathrm{gas} = M_\mathrm{gas} \frac{x_\mathrm{s} - x_\mathrm{gas}}{\tau},

   such that (see :galacticus-class:`nodeOperatorRamPressureMassLossDisks`) :math:`\dot{x}_\mathrm{gas} = (x_\mathrm{s} -
   x_\mathrm{gas})/\tau`. The timescale is :math:`\tau = \alpha T_\mathrm{vert}(r_\mathrm{s})`, where :math:`\alpha=`\
   ``[factorPeriodVertical]`` and :math:`T_\mathrm{vert} = 2 \pi/\sqrt{\partial^2\Phi/\partial z^2|_{z=0}}` is the period of
   vertical oscillations in the disk (their eqn. 1). If gas is stripped at all radii (i.e. within the smallest truncation
   radius modeled), the gas relaxes toward complete truncation, with :math:`\tau` evaluated at that radius.

   The stripping radius is found for the actual (truncated, and, if :galacticus-class:`galacticDiskHeatingClass` is active,
   heated) gas distribution, but :math:`x_\mathrm{s}` is found from the geometry of the truncated, *unheated* distribution. For a
   heated disk, which has gas beyond its truncation radius, this is an approximation: gas of the heated distribution remaining
   beyond :math:`r_\mathrm{s}` once :math:`x_\mathrm{gas}=x_\mathrm{s}` is not stripped.

   The vertical restoring force, :math:`|\partial \Phi/\partial z|`, of a truncated disk is taken from
   :galacticus-class:`massDistributionExponentialDiskTidal`, whose field away from the mid-plane is accurate to about 1.5% for
   truncated disks, and to about 2% for truncated and heated disks (see the documentation of that class).

   Only disk components are supported, and the disk must be truncated (i.e. ``[galacticComponentTruncation]`` must not be
   ``none``), with the same ``[galacticComponentTruncation]`` as the disk.
   </description>
  </ramPressureStripping>
  !!]
  type, extends(ramPressureStrippingClass) :: ramPressureStrippingKoppen2018
     !!{RST
     Implementation of a model of ram pressure stripping of disks following :cite:t:`koppen_ram_2018`.
     !!}
     private
     class           (hotHaloRamPressureHistoryClass  ), pointer :: hotHaloRamPressureHistory_   => null()
     class           (galacticComponentTruncationClass), pointer :: galacticComponentTruncation_ => null()
     double precision                                            :: factorPeriodVertical                        , enhancementMolecular, &
          &                                                         radiusEnhancementMolecular
     integer                                                     :: truncationGasID
     type            (rootFinder                      )          :: finder
     integer         (kind_int8                       )          :: uniqueIDLast                 =  -1_kind_int8
     double precision                                            :: radiusStrippingLast          =  -1.0d0
   contains
     final     ::                 koppen2018Destructor
     procedure :: rateMassLoss => koppen2018RateMassLoss
  end type ramPressureStrippingKoppen2018

  interface ramPressureStrippingKoppen2018
     !!{RST
     Constructors for the :galacticus-class:`ramPressureStrippingKoppen2018` model of ram pressure stripping of disks class.
     !!}
     module procedure koppen2018ConstructorParameters
     module procedure koppen2018ConstructorInternal
  end interface ramPressureStrippingKoppen2018

  ! The smallest radius (in units of the disk scale length) considered in finding the stripping radius. This is the smallest
  ! truncation radius modeled by galacticComponentTruncationExponential, within which the disk is regarded as completely truncated.
  double precision, parameter :: koppen2018RadiusMinimum   =1.0d-2
  ! The largest radius (in units of the disk scale length) considered in finding the stripping radius.
  double precision, parameter :: koppen2018RadiusMaximum   =3.0d+1
  ! The degree of truncation toward which gas stripped at all radii relaxes.
  double precision, parameter :: koppen2018TruncationLarge =2.5d+1

  ! Module-scope variables used in root finding.
  class           (massDistributionClass), pointer :: massDistribution_          , massDistributionGas_
  double precision                                 :: radiusDisk_                , forceMaximum_        , &
       &                                              momentumIntegrated_        , enhancementMolecular_, &
       &                                              radiusEnhancementMolecular_
  !$omp threadprivate(massDistribution_,massDistributionGas_,radiusDisk_,forceMaximum_,momentumIntegrated_,enhancementMolecular_,radiusEnhancementMolecular_)

contains

  function koppen2018ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`ramPressureStrippingKoppen2018` model of ram pressure stripping of disks class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (ramPressureStrippingKoppen2018  )                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (hotHaloRamPressureHistoryClass  ), pointer       :: hotHaloRamPressureHistory_
    class           (galacticComponentTruncationClass), pointer       :: galacticComponentTruncation_
    double precision                                                  :: factorPeriodVertical        , enhancementMolecular, &
         &                                                               radiusEnhancementMolecular

    !![
    <inputParameter docformat="rst">
      <name>factorPeriodVertical</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The timescale on which the gas of the disk relaxes toward truncation at the stripping radius, in units of the period of
      vertical oscillations at that radius.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>enhancementMolecular</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The factor, :math:`a`, by which the presence of molecular gas enhances the restoring force at the center of the disk
      :cite:p:`koppen_ram_2018`. Following :cite:t:`vollmer_ram_2001`, :cite:t:`koppen_ram_2018` use :math:`a=15`, with the gas
      surface density being that of atomic gas only. Since the gas surface density of the disk includes all gas, the default is zero.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusEnhancementMolecular</name>
      <defaultValue>0.67d0</defaultValue>
      <description>
      The radial scale, :math:`R_0`, of the enhancement of the restoring force by molecular gas, in units of the disk scale length.
      The default corresponds to the value of 2 kpc used by :cite:t:`koppen_ram_2018` for a disk with a scale length of 3 kpc.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="hotHaloRamPressureHistory"   name="hotHaloRamPressureHistory_"   source="parameters"/>
    <objectBuilder class="galacticComponentTruncation" name="galacticComponentTruncation_" source="parameters"/>
    !!]
    self=ramPressureStrippingKoppen2018(factorPeriodVertical,enhancementMolecular,radiusEnhancementMolecular,hotHaloRamPressureHistory_,galacticComponentTruncation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloRamPressureHistory_"  />
    <objectDestructor name="galacticComponentTruncation_"/>
    !!]
    return
  end function koppen2018ConstructorParameters

  function koppen2018ConstructorInternal(factorPeriodVertical,enhancementMolecular,radiusEnhancementMolecular,hotHaloRamPressureHistory_,galacticComponentTruncation_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`ramPressureStrippingKoppen2018` model of ram pressure stripping of disks
    class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (ramPressureStrippingKoppen2018  )                        :: self
    double precision                                  , intent(in   )         :: factorPeriodVertical        , enhancementMolecular, &
         &                                                                       radiusEnhancementMolecular
    class           (hotHaloRamPressureHistoryClass  ), intent(in   ), target :: hotHaloRamPressureHistory_
    class           (galacticComponentTruncationClass), intent(in   ), target :: galacticComponentTruncation_
    !![
    <constructorAssign variables="factorPeriodVertical, enhancementMolecular, radiusEnhancementMolecular, *hotHaloRamPressureHistory_, *galacticComponentTruncation_"/>
    !!]

    if (.not.self%galacticComponentTruncation_%isActive())                                                                        &
         & call Error_Report(                                                                                                     &
         &                   'the `koppen2018` ram pressure stripping class requires that the gas of the disk be truncated - ' // &
         &                   'set [galacticComponentTruncation] (e.g. to `exponential`) at the top level of the parameter file'// &
         &                   {introspection:location}                                                                             &
         &                  )
    !![
    <addMetaProperty component="disk" name="truncationGas" id="self%truncationGasID" isEvolvable="yes" isCreator="no"/>
    !!]
    self%finder=rootFinder(                                    &
         &                 rootFunction     =koppen2018Margin, &
         &                 toleranceAbsolute=0.0d+0          , &
         &                 toleranceRelative=1.0d-6            &
         &                )
    return
  end function koppen2018ConstructorInternal

  subroutine koppen2018Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`ramPressureStrippingKoppen2018` model of ram pressure stripping of disks class.
    !!}
    implicit none
    type(ramPressureStrippingKoppen2018), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloRamPressureHistory_"  />
    <objectDestructor name="self%galacticComponentTruncation_"/>
    !!]
    return
  end subroutine koppen2018Destructor

  double precision function koppen2018RateMassLoss(self,component) result(rateMassLoss)
    !!{RST
    Computes the mass loss rate from disks due to ram pressure stripping following :cite:t:`koppen_ram_2018`.
    !!}
    use :: Coordinates                     , only : coordinateCartesian      , assignment(=)
    use :: Display                         , only : displayGreen             , displayBlue                 , displayMagenta                , displayReset
    use :: Error                           , only : Error_Report             , errorStatusSuccess
    use :: Galactic_Structure_Options      , only : componentTypeDisk        , massTypeGaseous
    use :: Galacticus_Nodes                , only : nodeComponentDisk        , treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear                 , megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Root_Finder                     , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: Tensors                         , only : assignment(=)
    implicit none
    class           (ramPressureStrippingKoppen2018), intent(inout)  :: self
    class           (nodeComponent                 ), intent(inout)  :: component
    type            (treeNode                      ), pointer        :: node
    double precision                                , dimension(3,3) :: tidalTensor
    type            (coordinateCartesian           )                 :: coordinates
    double precision                                                 :: massGas        , truncationGas      , &
         &                                                              radiusStripping, truncationStripping, &
         &                                                              timescale      , curvatureVertical  , &
         &                                                              marginMinimum  , marginMaximum
    integer                                                          :: status
    logical                                                          :: foundRoot

    ! Assume no mass loss rate due to ram pressure by default.
    rateMassLoss=0.0d0
    ! Return immediately if this is not a satellite.
    node => component%hostNode
    if (.not.node%isSatellite()) return
    ! Get required properties of the disk.
    select type (component)
    class is (nodeComponentDisk)
       massGas      =component%massGas()
       radiusDisk_  =component%radius ()
       truncationGas=0.0d0
       if (component%floatRank0MetaPropertyIsCreated(self%truncationGasID)) truncationGas=max(0.0d0,component%floatRank0MetaPropertyGet(self%truncationGasID))
    class default
       massGas      =0.0d0
       truncationGas=0.0d0
       call Error_Report(                                                                                                                                                                                                                                                                     &
            &            'only "'//displayBlue()//'disk'//displayReset()//'" components are supported by the "'//displayGreen()//'koppen2018'//displayReset()//'" '//displayBlue()//'ramPressureStripping'//displayReset()//' class'//char(10)//                                              &
            &            displayGreen()//'HELP:'//displayReset()//' see '//displayMagenta()//'https://galacticus.readthedocs.io/en/latest/manuals/user-guide/troubleshooting/run-time-errors.html#inconsistent-assumptions-for-ram-pressure-models'//displayReset()//{introspection:location} &
            &           )
    end select
    if (massGas <= 0.0d0 .or. radiusDisk_ <= 0.0d0) return
    ! Get the maximum and time-integrated ram pressure. The time integral is converted to units of (km/s) M☉/Mpc², those of
    ! momentum per unit area.
    forceMaximum_      =+self%hotHaloRamPressureHistory_%forceMaximum       (node) 
    momentumIntegrated_=+self%hotHaloRamPressureHistory_%forceTimeIntegrated(node) &
         &              *kilo                                                      &
         &              *gigaYear                                                  &
         &              /megaParsec
    if (forceMaximum_ <= 0.0d0 .or. momentumIntegrated_ <= 0.0d0) return
    ! Get the mass distributions of the satellite, and of the gas of its disk.
    enhancementMolecular_      =  self%enhancementMolecular
    radiusEnhancementMolecular_=  self%radiusEnhancementMolecular
    massDistribution_    => node%massDistribution(                                                 )
    massDistributionGas_ => node%massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)
    ! Find the stripping radius.
    foundRoot=.false.
    if (node%uniqueID() == self%uniqueIDLast .and. self%radiusStrippingLast > 0.0d0) then
       ! The stripping radius was previously found for this node. Use it as an initial guess, expanding the range as necessary.
       ! (Note that the previous radius is not reset between ODE steps, as we specifically want to retain knowledge from the
       ! previous step.) If no root is found within the range considered the search is repeated below, after checking whether gas
       ! is stripped at all, or at all radii.
       call self%finder%rangeExpand(                                                                   &
            &                       rangeExpandDownward          =0.9d0                              , &
            &                       rangeExpandUpward            =1.1d0                              , &
            &                       rangeExpandType              =rangeExpandMultiplicative          , &
            &                       rangeDownwardLimit           =koppen2018RadiusMinimum*radiusDisk_, &
            &                       rangeUpwardLimit             =koppen2018RadiusMaximum*radiusDisk_, &
            &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative      , &
            &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive        &
            &                      )
       radiusStripping=self%finder%find(rootGuess=min(max(self%radiusStrippingLast,koppen2018RadiusMinimum*radiusDisk_),koppen2018RadiusMaximum*radiusDisk_),status=status)
       foundRoot      =status == errorStatusSuccess
    end if
    if (.not.foundRoot) then
       marginMinimum=koppen2018Margin(koppen2018RadiusMinimum*radiusDisk_)
       marginMaximum=koppen2018Margin(koppen2018RadiusMaximum*radiusDisk_)
       if      (marginMaximum <  0.0d0) then
          ! Gas is retained at all radii - no stripping occurs.
          radiusStripping=-1.0d0
       else if (marginMinimum >= 0.0d0) then
          ! Gas is stripped at all radii - relax toward complete truncation.
          radiusStripping=koppen2018RadiusMinimum*radiusDisk_
       else
          call self%finder%rangeExpand()
          radiusStripping=self%finder%find(rootRange=[koppen2018RadiusMinimum,koppen2018RadiusMaximum]*radiusDisk_,rootRangeValues=[marginMinimum,marginMaximum])
          foundRoot      =.true.
       end if
    end if
    if (foundRoot) then
       self%uniqueIDLast       =node%uniqueID()
       self%radiusStrippingLast=radiusStripping
       truncationStripping     =min(                                                                                                         &
            &                       -log(max(self%galacticComponentTruncation_%fractionMassRetainedCylindrical(radiusStripping/radiusDisk_), &
            &                                exp(-koppen2018TruncationLarge)                                                                 &
            &                               )                                                                                                &
            &                           )                                                                                                  , &
            &                       koppen2018TruncationLarge                                                                                &
            &                      )
    else if (radiusStripping > 0.0d0) then
       truncationStripping     =koppen2018TruncationLarge
    else
       truncationStripping     =0.0d0
    end if
    if (truncationStripping > truncationGas) then
       ! Find the period of vertical oscillations at the stripping radius. The tidal tensor is -∂²Φ/∂xᵢ∂xⱼ, in units of
       ! (km/s/Mpc)².
       coordinates      = [radiusStripping,0.0d0,0.0d0]
       tidalTensor      = massDistribution_%tidalTensor(coordinates)
       curvatureVertical=-tidalTensor(3,3)
       if (curvatureVertical > 0.0d0) then
          timescale   =+self%factorPeriodVertical &
               &       *2.0d0                     &
               &       *Pi                        &
               &       /sqrt(curvatureVertical)   &
               &       *megaParsec                &
               &       /kilo                      &
               &       /gigaYear
          rateMassLoss=+massGas                             &
               &       *(truncationStripping-truncationGas) &
               &       /timescale
       end if
    end if
    !![
    <objectDestructor name="massDistribution_"   />
    <objectDestructor name="massDistributionGas_"/>
    !!]
    return
  end function koppen2018RateMassLoss

  double precision function koppen2018Margin(radius) result(margin)
    !!{RST
    Return the margin by which the stripping criterion of :cite:t:`koppen_ram_2018` is satisfied at the given radius in the disk.
    This is positive if gas at that radius is stripped, and negative otherwise. If the time-integrated ram pressure, :math:`I`,
    exceeds :math:`I_0`, the margin is :math:`(1-I_0/I)-(p_0/P)(1-I_0/2I)`, which is zero where :math:`P (I-I_0) = p_0
    (I-I_0/2)`, and tends to :math:`1-p_0/P` (the criterion of Gunn & Gott 1972) for long pulses, and to
    :math:`1-I_0/I` for short pulses. Otherwise, it is :math:`(1-I_0/I)-p_0/2P`, which is negative, and continuous with the
    former at :math:`I=I_0`.
    !!}
    use :: Coordinates                     , only : coordinateCylindrical, assignment(=)
    use :: Numerical_Constants_Astronomical, only : gigaYear             , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                       , intent(in   ) :: radius
    type            (coordinateCylindrical)                :: coordinates
    double precision                                       :: surfaceDensityGas, forceRestoring, &
         &                                                    velocitySquared  , momentumEscape

    coordinates      =[radius,0.0d0,0.0d0]
    surfaceDensityGas=massDistributionGas_%surfaceDensity(coordinates)
    ! Maximum restoring force per unit area, converting the acceleration from km/s/Gyr to (km/s)²/Mpc.
    forceRestoring   =+(                                                      &
         &              +1.0d0                                                &
         &              +enhancementMolecular_                                &
         &              *exp(-radius/radiusEnhancementMolecular_/radiusDisk_) &
         &             )                                                      &
         &            *surfaceDensityGas                                      &
         &            *koppen2018AccelerationVerticalMaximum(radius)          &
         &            *megaParsec                                             &
         &            /kilo                                                   &
         &            /gigaYear
    ! Momentum per unit area needed for gas to escape.
    velocitySquared  =max(                                                        &
         &                -2.0d0*massDistribution_%potential    (coordinates)     &
         &                -      massDistribution_%rotationCurve(radius     )**2, &
         &                +0.0d0                                                  &
         &               )
    momentumEscape   =+surfaceDensityGas     &
         &            *sqrt(velocitySquared)
    if (momentumIntegrated_ > momentumEscape) then
       margin=+(1.0d0-      momentumEscape/momentumIntegrated_) &
            & -forceRestoring/forceMaximum_                     &
            & *(1.0d0-0.5d0*momentumEscape/momentumIntegrated_)
    else
       margin=+(1.0d0-      momentumEscape/momentumIntegrated_) &
            & -forceRestoring/forceMaximum_                     &
            & *0.5d0
    end if
    return
  end function koppen2018Margin

  double precision function koppen2018AccelerationVerticalMaximum(radius) result(accelerationMaximum)
    !!{RST
    Return the maximum, over height above the mid-plane, of the magnitude of the vertical gravitational acceleration (in units
    of km/s/Gyr) at the given radius in the disk. The acceleration is evaluated on a logarithmic grid of heights, from
    :math:`10^{-3}` disk scale lengths to the larger of 10 disk scale lengths and three times the radius (since, where the
    potential is dominated by a spherical halo, the maximum occurs at a height comparable to the radius), and the maximum refined
    by golden section search in the logarithm of height. Since the acceleration is stationary at its maximum, the tolerance of
    :math:`10^{-2}` in the logarithm of height gives the maximum to a precision of order :math:`10^{-4}`.
    !!}
    use :: Coordinates, only : coordinateCartesian, assignment(=)
    implicit none
    double precision                                , intent(in   )           :: radius
    integer                                         , parameter               :: countHeights          =9
    double precision                                , parameter               :: heightMinimum         =1.0d-3              , heightMaximum          =1.0d+1, &
         &                                                                       toleranceLogarithmic  =1.0d-2              , factorRadius           =3.0d+0
    double precision                                , parameter               :: ratioGolden           =0.3819660112501051d0
    double precision                                , dimension(countHeights) :: accelerations
    double precision                                                          :: heightLogarithmicLower                      , heightLogarithmicUpper       , &
         &                                                                       heightLogarithmic1                          , heightLogarithmic2           , &
         &                                                                       acceleration1                               , acceleration2                , &
         &                                                                       intervalLogarithmic                         , heightMaximumActual
    integer                                                                   :: i                                           , iMaximum

    ! Evaluate the acceleration on a grid of heights (in units of the disk scale length).
    heightMaximumActual=max(heightMaximum,factorRadius*radius/radiusDisk_)
    intervalLogarithmic=log(heightMaximumActual/heightMinimum)/dble(countHeights-1)
    do i=1,countHeights
       accelerations(i)=accelerationVertical(log(heightMinimum)+dble(i-1)*intervalLogarithmic)
    end do
    iMaximum=maxloc(accelerations,dim=1)
    ! Refine the maximum by golden section search between the neighboring grid points.
    heightLogarithmicLower=log(heightMinimum)+dble(max(iMaximum-1,           1)-1)*intervalLogarithmic
    heightLogarithmicUpper=log(heightMinimum)+dble(min(iMaximum+1,countHeights)-1)*intervalLogarithmic
    heightLogarithmic1    =heightLogarithmicUpper-(1.0d0-ratioGolden)*(heightLogarithmicUpper-heightLogarithmicLower)
    heightLogarithmic2    =heightLogarithmicLower+(1.0d0-ratioGolden)*(heightLogarithmicUpper-heightLogarithmicLower)
    acceleration1         =accelerationVertical(heightLogarithmic1)
    acceleration2         =accelerationVertical(heightLogarithmic2)
    do while (heightLogarithmicUpper-heightLogarithmicLower > toleranceLogarithmic)
       if (acceleration1 > acceleration2) then
          heightLogarithmicUpper=heightLogarithmic2
          heightLogarithmic2    =heightLogarithmic1
          acceleration2         =acceleration1
          heightLogarithmic1    =heightLogarithmicUpper-(1.0d0-ratioGolden)*(heightLogarithmicUpper-heightLogarithmicLower)
          acceleration1         =accelerationVertical(heightLogarithmic1)
       else
          heightLogarithmicLower=heightLogarithmic1
          heightLogarithmic1    =heightLogarithmic2
          acceleration1         =acceleration2
          heightLogarithmic2    =heightLogarithmicLower+(1.0d0-ratioGolden)*(heightLogarithmicUpper-heightLogarithmicLower)
          acceleration2         =accelerationVertical(heightLogarithmic2)
       end if
    end do
    accelerationMaximum=max(maxval(accelerations),acceleration1,acceleration2)
    return

  contains

    double precision function accelerationVertical(heightLogarithmic)
      !!{RST
      Return the magnitude of the vertical acceleration at the given logarithm of height (in units of the disk scale length).
      !!}
      implicit none
      double precision                      , intent(in   ) :: heightLogarithmic
      type            (coordinateCartesian)                :: coordinates
      double precision, dimension(3)                       :: acceleration

      coordinates         =[radius,0.0d0,radiusDisk_*exp(heightLogarithmic)]
      acceleration        =massDistribution_%acceleration(coordinates)
      accelerationVertical=abs(acceleration(3))
      return
    end function accelerationVertical

  end function koppen2018AccelerationVerticalMaximum
