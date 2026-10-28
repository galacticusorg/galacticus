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
  An implementation of virial orbits in which satellites fall in preferentially along the cosmic web filaments feeding their
  host, applied to another virial orbit class.
  !!}

  use :: Cosmological_Density_Field, only : criticalOverdensityClass
  use :: Dark_Matter_Halo_Filaments, only : darkMatterHaloFilamentsClass

  !![
  <virialOrbit name="virialOrbitFilamentary" docformat="rst">
   <description>
   A virial orbit class which assigns infall directions to the orbits returned by another virial orbit class such that
   satellites fall in preferentially along the cosmic web filaments feeding their host. The directions of those filaments,
   which persist along each branch of the merger tree, are given by the :galacticus-class:`darkMatterHaloFilamentsClass`
   object.

   The direction of infall, :math:`\hat{\mathbf{n}}`, is drawn from a mixture of an isotropic distribution and a von
   Mises--Fisher distribution about each filament axis, :math:`\hat{\mathbf{f}}_i`:

   .. math::

      p(\hat{\mathbf{n}}) = \frac{1-f}{4\pi} + \frac{f}{N} \sum_{i=1}^N \frac{\kappa}{4\pi \sinh\kappa}
      \exp\left(\kappa\, \hat{\mathbf{n}}\cdot\hat{\mathbf{f}}_i\right),

   where :math:`N` is the number of filaments feeding the host, :math:`f` is the fraction of satellites accreted along filaments,
   and :math:`\kappa` sets the angular width of each filament. The tangential velocity is oriented isotropically, so that the mean
   vector tangential velocity and the mean vector angular momentum both remain zero.

   The fraction accreted along filaments depends on the mass of the host, :math:`M`, relative to the mass, :math:`M_*`, of
   halos just collapsing at the time of accretion, :math:`\sigma(M_*,t) = \delta_\mathrm{c}(t)`:

   .. math::

      f = f_* + \alpha \log_{10}(M/M_*),

   limited to the range :math:`0 \le f \le 1`, where :math:`f_*` is given by ``fractionFilamentaryCharacteristic`` and
   :math:`\alpha` by ``fractionFilamentaryLogarithmicSlope``. The width of a filament is specified, by
   ``radiusFilamentFractional``, as its radius, :math:`x`, in units of the virial radius of the host. Its angular half-width
   on the virial sphere is then :math:`\theta_\mathrm{half} = \sin^{-1}[\min(1,x)]`, and :math:`\kappa = 1/\theta_\mathrm{half}^2`.

   Two consequences of this convention should be noted. First, :math:`\theta_\mathrm{half}` is the dispersion of the angle
   from the axis in the limit of large :math:`\kappa` (for which that distribution is Gaussian, with
   :math:`\sigma=\kappa^{-1/2}`), not the half-width at half maximum, which is larger (:math:`31^\circ` for the default
   :math:`\theta_\mathrm{half}=26^\circ`). Second, :math:`\kappa` does not vanish for wide filaments: any :math:`x \ge 1`
   gives :math:`\theta_\mathrm{half}=\pi/2` and :math:`\kappa=4/\pi^2\approx 0.41`, for which infall directions of the
   filamentary component are still :math:`\cosh(4/\pi^2)\approx 1.08` times more probable along the filament than
   perpendicular to it. Increasing the width of the filaments therefore never makes infall fully isotropic. Isotropic infall
   is obtained only with :math:`f=0`, i.e. with ``fractionFilamentaryCharacteristic`` and
   ``fractionFilamentaryLogarithmicSlope`` both set to zero.

   The default parameters are from a joint fit of this model to the distributions of the angles between the points at which
   subhalos enter their hosts and both the axes of slowest (:math:`\mathbf{e}_3`) and fastest (:math:`\mathbf{e}_1`) collapse of
   the velocity shear tensor (smoothed on a scale of four virial radii) measured by :cite:t:`libeskind_universal_2014`, in four
   bins of :math:`M/M_*`. The model was evaluated with the default :galacticus-class:`darkMatterHaloFilamentsRandomSequential`
   filaments (three filaments per halo: one passing through the halo, plus one further filament at least :math:`60^\circ` from
   it), identifying :math:`\mathbf{e}_3` and :math:`\mathbf{e}_1` with the major and minor axes of the tensor
   :math:`\sum_i \hat{\mathbf{f}}_i \hat{\mathbf{f}}_i^\mathrm{T}` of the filament directions. With a single :math:`\kappa` the fit gives
   :math:`f_*=0.291`, :math:`\alpha=0.102`, and :math:`\kappa=5.0`. It matches both distributions in each mass bin with root mean
   square residuals of 0.009--0.033. The number of filaments is constrained by :math:`\mathbf{e}_1`: infall avoids the axis of
   fastest collapse in all bins above :math:`M/M_*=0.1`, which a single filament (being axisymmetric) cannot reproduce, but which
   follows from filaments lying preferentially in a plane. Fewer than about 2.8 filaments per halo on average are strongly
   disfavored. The concentration shows no trend with mass, consistent with the profiles of filaments being universal when scaled
   to the virial radius of the halo they feed :cite:p:`xu_universal_2026`, which is why the width is specified in units of the
   virial radius.  Since the measurement averages over all accretion since :math:`z \approx 5`, the fit describes accretion
   averaged over redshift.

   The calibration depends on the filaments: with more filaments per halo, a larger :math:`f` is needed to give the same
   anisotropy relative to any one axis. For example, for a single filament passing through each halo (a mean number of two
   filaments in :galacticus-class:`darkMatterHaloFilamentsRandomSequential`) a fit to the distribution relative to
   :math:`\mathbf{e}_3` alone gives :math:`f_*=0.184`, :math:`\alpha=0.055`, and :math:`\kappa=4.49`
   (``radiusFilamentFractional`` :math:`=0.4546`), but reproduces only about half of the avoidance of :math:`\mathbf{e}_1` in
   the most massive hosts.

   This class sets the angles :math:`\theta`, :math:`\phi`, and :math:`\epsilon` of the orbit, so wrapping another class which
   sets them (such as :galacticus-class:`virialOrbitIsotropic` or :galacticus-class:`virialOrbitSpinCorrelated`) is not
   permitted. Conversely, wrapping this class in one of those classes would discard the angles that it sets.
   </description>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitFilamentary
     !!{RST
     A virial orbit class in which satellites fall in preferentially along the cosmic web filaments feeding their host.
     !!}
     private
     class           (virialOrbitClass            ), pointer :: virialOrbit_                      => null()
     class           (criticalOverdensityClass    ), pointer :: criticalOverdensity_              => null()
     class           (darkMatterHaloFilamentsClass), pointer :: darkMatterHaloFilaments_          => null()
     double precision                                        :: fractionFilamentaryCharacteristic          , fractionFilamentaryLogarithmicSlope, &
          &                                                     radiusFilamentFractional                   , concentration                      , &
          &                                                     timePrevious                      =-1.0d0  , massCollapsingPrevious
   contains
     final     ::                                    filamentaryDestructor
     procedure :: orbit                           => filamentaryOrbit
     procedure :: densityContrastDefinition       => filamentaryDensityContrastDefinition
     procedure :: velocityTangentialMagnitudeMean => filamentaryVelocityTangentialMagnitudeMean
     procedure :: velocityTangentialVectorMean    => filamentaryVelocityTangentialVectorMean
     procedure :: angularMomentumMagnitudeMean    => filamentaryAngularMomentumMagnitudeMean
     procedure :: angularMomentumVectorMean       => filamentaryAngularMomentumVectorMean
     procedure :: velocityTotalRootMeanSquared    => filamentaryVelocityTotalRootMeanSquared
     procedure :: energyMean                      => filamentaryEnergyMean
     procedure :: isAngularlyResolved             => filamentaryIsAngularlyResolved
     procedure :: fractionFilamentary             => filamentaryFractionFilamentary
  end type virialOrbitFilamentary

  interface virialOrbitFilamentary
     !!{RST
     Constructors for the :galacticus-class:`virialOrbitFilamentary` virial orbits class.
     !!}
     module procedure filamentaryConstructorParameters
     module procedure filamentaryConstructorInternal
  end interface virialOrbitFilamentary

contains

  function filamentaryConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`virialOrbitFilamentary` virial orbits class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (virialOrbitFilamentary      )                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (virialOrbitClass            ), pointer       :: virialOrbit_
    class           (criticalOverdensityClass    ), pointer       :: criticalOverdensity_
    class           (darkMatterHaloFilamentsClass), pointer       :: darkMatterHaloFilaments_
    double precision                                              :: fractionFilamentaryCharacteristic, fractionFilamentaryLogarithmicSlope, &
         &                                                           radiusFilamentFractional

    !![
    <inputParameter docformat="rst">
      <name>fractionFilamentaryCharacteristic</name>
      <defaultValue>0.291d0</defaultValue>
      <defaultSource>Fit to :cite:t:`libeskind_universal_2014`, for the default filaments.</defaultSource>
      <source>parameters</source>
      <description>
      The fraction, :math:`f_*`, of satellites accreted along filaments by hosts of mass equal to the mass of halos just
      collapsing, :math:`M_*`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>fractionFilamentaryLogarithmicSlope</name>
      <defaultValue>0.102d0</defaultValue>
      <defaultSource>Fit to :cite:t:`libeskind_universal_2014`, for the default filaments.</defaultSource>
      <source>parameters</source>
      <description>
      The rate of change, :math:`\alpha=\mathrm{d}f/\mathrm{d}\log_{10}(M/M_*)`, of the fraction of satellites accreted along
      filaments with the mass of the host.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusFilamentFractional</name>
      <defaultValue>0.4325d0</defaultValue>
      <defaultSource>Fit to :cite:t:`libeskind_universal_2014`, corresponding to a concentration :math:`\kappa=5.0`.</defaultSource>
      <source>parameters</source>
      <description>
      The radius of a filament, :math:`x`, in units of the virial radius of the host, which sets the angular half-width,
      :math:`\theta_\mathrm{half} = \sin^{-1}[\min(1,x)]`, of the filament on the virial sphere. All values :math:`x \ge 1`
      give the same, weak but non-zero, concentration of infall about the filament, :math:`\kappa=4/\pi^2`. Infall is
      isotropic only if the filamentary fraction is zero.
      </description>
    </inputParameter>
    <objectBuilder class="virialOrbit"             name="virialOrbit_"             source="parameters"/>
    <objectBuilder class="criticalOverdensity"     name="criticalOverdensity_"     source="parameters"/>
    <objectBuilder class="darkMatterHaloFilaments" name="darkMatterHaloFilaments_" source="parameters"/>
    !!]
    self=virialOrbitFilamentary(fractionFilamentaryCharacteristic,fractionFilamentaryLogarithmicSlope,radiusFilamentFractional,virialOrbit_,criticalOverdensity_,darkMatterHaloFilaments_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialOrbit_"            />
    <objectDestructor name="criticalOverdensity_"    />
    <objectDestructor name="darkMatterHaloFilaments_"/>
    !!]
    return
  end function filamentaryConstructorParameters

  function filamentaryConstructorInternal(fractionFilamentaryCharacteristic,fractionFilamentaryLogarithmicSlope,radiusFilamentFractional,virialOrbit_,criticalOverdensity_,darkMatterHaloFilaments_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`virialOrbitFilamentary` virial orbits class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (virialOrbitFilamentary      )                        :: self
    double precision                              , intent(in   )         :: fractionFilamentaryCharacteristic, fractionFilamentaryLogarithmicSlope, &
         &                                                                   radiusFilamentFractional
    class           (virialOrbitClass            ), intent(in   ), target :: virialOrbit_
    class           (criticalOverdensityClass    ), intent(in   ), target :: criticalOverdensity_
    class           (darkMatterHaloFilamentsClass), intent(in   ), target :: darkMatterHaloFilaments_
    !![
    <constructorAssign variables="fractionFilamentaryCharacteristic, fractionFilamentaryLogarithmicSlope, radiusFilamentFractional, *virialOrbit_, *criticalOverdensity_, *darkMatterHaloFilaments_"/>
    !!]

    if (radiusFilamentFractional <= 0.0d0) call Error_Report('[radiusFilamentFractional] must be positive'//{introspection:location})
    ! The wrapped class must not itself set the angles of the orbit, as they would be silently overwritten.
    if (self%virialOrbit_%isAngularlyResolved()) call Error_Report('the wrapped virial orbit class sets the angles of the orbit, which this class would overwrite'//{introspection:location})
    ! Find the concentration of the von Mises-Fisher distribution from the angular half-width of the filament.
    self%concentration=1.0d0/asin(min(1.0d0,radiusFilamentFractional))**2
    return
  end function filamentaryConstructorInternal

  subroutine filamentaryDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`virialOrbitFilamentary` virial orbits class.
    !!}
    implicit none
    type(virialOrbitFilamentary), intent(inout) :: self

    !![
    <objectDestructor name="self%virialOrbit_"            />
    <objectDestructor name="self%criticalOverdensity_"    />
    <objectDestructor name="self%darkMatterHaloFilaments_"/>
    !!]
    return
  end subroutine filamentaryDestructor

  function filamentaryOrbit(self,node,host,acceptUnboundOrbits)
    !!{RST
    Return orbital parameters for a satellite, with the direction of infall drawn preferentially along the filaments feeding
    the host.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Vectors                 , only : Vector_Product
    implicit none
    type            (keplerOrbit           )                                :: filamentaryOrbit
    class           (virialOrbitFilamentary), intent(inout), target         :: self
    type            (treeNode              ), intent(inout)                 :: host               , node
    logical                                 , intent(in   )                 :: acceptUnboundOrbits
    double precision                        , allocatable  , dimension(:,:) :: axes
    double precision                                       , dimension(3  ) :: axis               , direction, &
         &                                                                     tangent            , binormal
    double precision                                                        :: cosineTheta        , phi      , &
         &                                                                     uniformDeviate
    integer                                                                 :: i

    ! Get the underlying orbit.
    filamentaryOrbit=self%virialOrbit_%orbit(node,host,acceptUnboundOrbits)
    ! Choose between accretion along a filament and isotropic accretion.
    if (node%hostTree%randomNumberGenerator_%uniformSample() < self%fractionFilamentary(host)) then
       ! Choose one of the filaments of the host at random.
       axes=self%darkMatterHaloFilaments_%axes(host)
       i   =min(int(dble(size(axes,dim=2))*node%hostTree%randomNumberGenerator_%uniformSample())+1,size(axes,dim=2))
       axis=axes(:,i)
       ! Draw the angle from the axis from the von Mises-Fisher distribution. Its cumulative distribution inverts analytically.
       uniformDeviate=node%hostTree%randomNumberGenerator_%uniformSample()
       cosineTheta   =+1.0d0                                                                     &
            &         +log(uniformDeviate+(1.0d0-uniformDeviate)*exp(-2.0d0*self%concentration)) &
            &         /                                                     self%concentration
       cosineTheta   =max(-1.0d0,min(+1.0d0,cosineTheta))
       phi           =2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()
       ! Construct an orthonormal basis about the axis, using whichever Cartesian axis is far from it.
       if (abs(axis(1)) < 0.9d0) then
          tangent=Vector_Product(axis,[1.0d0,0.0d0,0.0d0])
       else
          tangent=Vector_Product(axis,[0.0d0,1.0d0,0.0d0])
       end if
       tangent  =tangent/sqrt(sum(tangent**2))
       binormal =Vector_Product(axis,tangent)
       direction=+           cosineTheta                               *axis      &
            &    +sqrt(1.0d0-cosineTheta**2)*(cos(phi)*tangent+sin(phi)*binormal)
    else
       ! Draw from an isotropic distribution.
       cosineTheta=+2.0d0   *node%hostTree%randomNumberGenerator_%uniformSample()-1.0d0
       phi        =+2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()
       direction  =[                                     &
            &       sqrt(1.0d0-cosineTheta**2)*cos(phi), &
            &       sqrt(1.0d0-cosineTheta**2)*sin(phi), &
            &                  cosineTheta               &
            &      ]
    end if
    call filamentaryOrbit%thetaSet(acos  (max(-1.0d0,min(+1.0d0,direction(3)))))
    call filamentaryOrbit%phiSet  (modulo(atan2(direction(2),direction(1)),2.0d0*Pi))
    ! The tangential velocity is oriented isotropically.
    call filamentaryOrbit%epsilonSet(2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample())
    return
  end function filamentaryOrbit

  double precision function filamentaryFractionFilamentary(self,host)
    !!{RST
    Return the fraction of satellites accreted along filaments by the given host.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(virialOrbitFilamentary), intent(inout) :: self
    type (treeNode              ), intent(inout) :: host
    class(nodeComponentBasic    ), pointer       :: basic

    basic => host%basic()
    ! The collapsing mass depends only on time, so retain the last value computed.
    if (basic%time() /= self%timePrevious) then
       self%timePrevious          =                                              basic%time()
       self%massCollapsingPrevious=self%criticalOverdensity_%collapsingMass(time=basic%time())
    end if
    if (self%fractionFilamentaryLogarithmicSlope == 0.0d0) then
       filamentaryFractionFilamentary=+self%fractionFilamentaryCharacteristic
    else if (self%massCollapsingPrevious <= 0.0d0) then
       ! The collapsing mass is too small to be found, so the host is effectively of infinite mass in units of it.
       filamentaryFractionFilamentary=merge(1.0d0,0.0d0,self%fractionFilamentaryLogarithmicSlope > 0.0d0)
    else
       filamentaryFractionFilamentary=+self%fractionFilamentaryCharacteristic   &
            &                         +self%fractionFilamentaryLogarithmicSlope &
            &                         *log10(                                   &
            &                                +basic%mass                  ()    &
            &                                /self %massCollapsingPrevious      &
            &                               )
    end if
    filamentaryFractionFilamentary=max(0.0d0,min(1.0d0,filamentaryFractionFilamentary))
    return
  end function filamentaryFractionFilamentary

  function filamentaryDensityContrastDefinition(self)
    !!{RST
    Return a virial density contrast object defining that used in the definition of virial orbits.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: filamentaryDensityContrastDefinition
    class(virialOrbitFilamentary    ), intent(inout) :: self

    filamentaryDensityContrastDefinition => self%virialOrbit_%densityContrastDefinition()
    return
  end function filamentaryDensityContrastDefinition

  double precision function filamentaryVelocityTangentialMagnitudeMean(self,node,host)
    !!{RST
    Return the mean magnitude of the tangential velocity.
    !!}
    implicit none
    class(virialOrbitFilamentary), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node, host

    filamentaryVelocityTangentialMagnitudeMean=self%virialOrbit_%velocityTangentialMagnitudeMean(node,host)
    return
  end function filamentaryVelocityTangentialMagnitudeMean

  function filamentaryVelocityTangentialVectorMean(self,node,host)
    !!{RST
    Return the mean of the vector tangential velocity.
    !!}
    implicit none
    double precision                        , dimension(3)  :: filamentaryVelocityTangentialVectorMean
    class           (virialOrbitFilamentary), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node                                   , host
    !$GLC attributes unused :: self, node, host

    ! The tangential velocity is oriented isotropically in the plane tangent to the virial sphere, so its mean is zero at any
    ! point of infall, and therefore zero however anisotropic the distribution of infall directions.
    filamentaryVelocityTangentialVectorMean=0.0d0
    return
  end function filamentaryVelocityTangentialVectorMean

  double precision function filamentaryAngularMomentumMagnitudeMean(self,node,host)
    !!{RST
    Return the mean magnitude of the angular momentum.
    !!}
    implicit none
    class(virialOrbitFilamentary), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node, host

    filamentaryAngularMomentumMagnitudeMean=self%virialOrbit_%angularMomentumMagnitudeMean(node,host)
    return
  end function filamentaryAngularMomentumMagnitudeMean

  function filamentaryAngularMomentumVectorMean(self,node,host)
    !!{RST
    Return the mean of the vector angular momentum.
    !!}
    implicit none
    double precision                        , dimension(3)  :: filamentaryAngularMomentumVectorMean
    class           (virialOrbitFilamentary), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node                                , host
    !$GLC attributes unused :: self, node, host

    ! The tangential velocity is oriented isotropically in the plane tangent to the virial sphere, so the mean angular momentum
    ! is zero at any point of infall, and therefore zero however anisotropic the distribution of infall directions.
    filamentaryAngularMomentumVectorMean=0.0d0
    return
  end function filamentaryAngularMomentumVectorMean

  double precision function filamentaryVelocityTotalRootMeanSquared(self,node,host)
    !!{RST
    Return the root mean squared of the total velocity.
    !!}
    implicit none
    class(virialOrbitFilamentary), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node, host

    filamentaryVelocityTotalRootMeanSquared=self%virialOrbit_%velocityTotalRootMeanSquared(node,host)
    return
  end function filamentaryVelocityTotalRootMeanSquared

  double precision function filamentaryEnergyMean(self,node,host)
    !!{RST
    Return the mean of the total energy.
    !!}
    implicit none
    class(virialOrbitFilamentary), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node, host

    filamentaryEnergyMean=self%virialOrbit_%energyMean(node,host)
    return
  end function filamentaryEnergyMean

  logical function filamentaryIsAngularlyResolved(self)
    !!{RST
    Return true indicating that orbits are angularly-resolved.
    !!}
    implicit none
    class(virialOrbitFilamentary), intent(inout) :: self
    !$GLC attributes unused :: self

    filamentaryIsAngularlyResolved=.true.
    return
  end function filamentaryIsAngularlyResolved
