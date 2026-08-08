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

  !!{RST
  An implementation of dark matter halo splashback radii using the fitting functions of :cite:t:`shi_outer_2016`.
  !!}

  !+ Implemented by Andrew Benson with assistance from Claude.

  use :: Cosmology_Functions                        , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                       , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Mass_Accretion_Histories  , only : darkMatterHaloMassAccretionHistoryClass
  use :: Virial_Density_Contrast                    , only : virialDensityContrastClass
  use :: Dark_Matter_Halo_Splashback_Radii_Reference, only : splashbackReferenceScales              , splashbackDensityContrastReference

  ! The smallest mass accretion rate for which the Shi (2016) fitting functions will be evaluated. The fitting function for
  ! the splashback radius involves ln(s), and so diverges as s→0.
  double precision, parameter, public :: shi2016AccretionRateMinimum=1.0d-3
  ! The fitting functions are made public so that they can be validated directly against reference values by the
  ! "tests.splashback_radii.exe" unit test.
  public :: shi2016RadiusRatioEvaluate, shi2016DensityContrastEvaluate

  !![
  <darkMatterHaloSplashbackRadius name="darkMatterHaloSplashbackRadiusShi2016" docformat="rst">
   <description>
   A splashback radius class implementing the fitting functions of :cite:t:`shi_outer_2016`, which are calibrated against
   the self-similar spherical collapse solution for the accretion flow around a halo. Specifically,

   .. math::

      {R_\mathrm{sp} \over R_\mathrm{200m}} = \exp\left[ \left( 0.24 + 0.074 \ln s \right) \ln \Omega_\mathrm{M} + 0.55 - 0.15 s \right],

   and the mean overdensity enclosed by the splashback radius is

   .. math::

      \Delta_\mathrm{sp} = 33 \Omega_\mathrm{M}^{-0.45} \exp\left[ \left( 0.88 - 0.14 \ln \Omega_\mathrm{M} \right) s^{0.6} \right],

   from which the splashback mass follows as :math:`M_\mathrm{sp}/M_\mathrm{200m} = (\Delta_\mathrm{sp}/200)
   (R_\mathrm{sp}/R_\mathrm{200m})^3`. Here :math:`s = \mathrm{d}\log M / \mathrm{d}\log a` is the *instantaneous* mass
   accretion rate---note that this differs from the :math:`\Gamma_\mathrm{dyn}` averaged over a dynamical time which is used
   by the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2017` and
   :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2020` classes.

   These fitting functions were calibrated over the range :math:`0.5 \le s \le 5`; results outside of that range should be
   treated as extrapolations. This class makes no prediction for the scatter in :math:`R_\mathrm{sp}` or
   :math:`M_\mathrm{sp}`, so attempting to evaluate those will result in an error.
   </description>
  </darkMatterHaloSplashbackRadius>
  !!]
  type, extends(darkMatterHaloSplashbackRadiusClass) :: darkMatterHaloSplashbackRadiusShi2016
     !!{RST
     A dark matter halo splashback radius class implementing the fitting functions of :cite:t:`shi_outer_2016`.
     !!}
     private
     class(cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_                 => null()
     class(cosmologyParametersClass               ), pointer :: cosmologyParameters_                => null()
     class(virialDensityContrastClass             ), pointer :: virialDensityContrast_              => null()
     class(darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
   contains
     !![
     <methods docformat="rst">
       <method method="propertiesGet" description="Return the properties of the halo needed to evaluate the splashback fitting functions."/>
     </methods>
     !!]
     final     ::                  shi2016Destructor
     procedure :: radius        => shi2016Radius
     procedure :: radiusRatio   => shi2016RadiusRatio
     procedure :: mass          => shi2016Mass
     procedure :: massRatio     => shi2016MassRatio
     procedure :: propertiesGet => shi2016PropertiesGet
  end type darkMatterHaloSplashbackRadiusShi2016

  interface darkMatterHaloSplashbackRadiusShi2016
     !!{RST
     Constructors for the :galacticus-class:`darkMatterHaloSplashbackRadiusShi2016` splashback radius class.
     !!}
     module procedure shi2016ConstructorParameters
     module procedure shi2016ConstructorInternal
  end interface darkMatterHaloSplashbackRadiusShi2016

contains

  function shi2016ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusShi2016` splashback radius class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (darkMatterHaloSplashbackRadiusShi2016  )                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class(virialDensityContrastClass             ), pointer       :: virialDensityContrast_
    class(darkMatterHaloMassAccretionHistoryClass), pointer       :: darkMatterHaloMassAccretionHistory_

    !![
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"/>
    <objectBuilder class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"/>
    <objectBuilder class="virialDensityContrast"              name="virialDensityContrast_"              source="parameters"/>
    <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    !!]
    self=darkMatterHaloSplashbackRadiusShi2016(cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,darkMatterHaloMassAccretionHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="cosmologyParameters_"               />
    <objectDestructor name="virialDensityContrast_"             />
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
  end function shi2016ConstructorParameters

  function shi2016ConstructorInternal(cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,darkMatterHaloMassAccretionHistory_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusShi2016` splashback radius class.
    !!}
    implicit none
    type (darkMatterHaloSplashbackRadiusShi2016  )                        :: self
    class(cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologyParametersClass               ), intent(in   ), target :: cosmologyParameters_
    class(virialDensityContrastClass             ), intent(in   ), target :: virialDensityContrast_
    class(darkMatterHaloMassAccretionHistoryClass), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_, *virialDensityContrast_, *darkMatterHaloMassAccretionHistory_"/>
    !!]

    return
  end function shi2016ConstructorInternal

  subroutine shi2016Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusShi2016` splashback radius class.
    !!}
    implicit none
    type(darkMatterHaloSplashbackRadiusShi2016), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%cosmologyParameters_"               />
    <objectDestructor name="self%virialDensityContrast_"             />
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
  end subroutine shi2016Destructor

  subroutine shi2016PropertiesGet(self,node,accretionRate,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    !!{RST
    Return the properties of the halo in ``node`` needed to evaluate the :cite:t:`shi_outer_2016` splashback fitting
    functions---namely the instantaneous mass accretion rate, :math:`s = \mathrm{d}\log M / \mathrm{d}\log a`, the matter
    density parameter, :math:`\Omega_\mathrm{M}`, and :math:`R_\mathrm{200m}` and :math:`M_\mathrm{200m}`.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloSplashbackRadiusShi2016), intent(inout)           :: self
    type            (treeNode                             ), intent(inout)           :: node
    double precision                                       , intent(  out)           :: accretionRate    , omegaMatter, &
         &                                                                              radius200Mean    , mass200Mean
    class           (massDistributionClass                ), intent(inout), optional :: massDistribution_
    class           (nodeComponentBasic                   ), pointer                 :: basic
    double precision                                                                 :: time

    basic        =>  node                     %basic             (    )
    time         =   basic                    %time              (    )
    omegaMatter  =   self%cosmologyFunctions_ %omegaMatterEpochal(time)
    call splashbackReferenceScales(                                                     &
         &                                                node                        , &
         &                         cosmologyParameters_  =self%cosmologyParameters_   , &
         &                         cosmologyFunctions_   =self%cosmologyFunctions_    , &
         &                         virialDensityContrast_=self%virialDensityContrast_ , &
         &                         mass                  =     mass200Mean            , &
         &                         radius                =     radius200Mean          , &
         &                         massDistribution_     =     massDistribution_        &
         &                        )
    ! Compute the instantaneous logarithmic growth rate of the halo mass, s = dlogM/dloga = (dM/dt)/(M H).
    if (basic%mass() > 0.0d0) then
       accretionRate=+self%darkMatterHaloMassAccretionHistory_%massAccretionRate(node,time)                                      &
            &        /basic                                   %mass             (         )                                      &
            &        /self%cosmologyFunctions_                %expansionRate    (self%cosmologyFunctions_%expansionFactor(time))
    else
       accretionRate=+0.0d0
    end if
    ! The fitting functions involve ln(s), and so must not be evaluated at, or below, zero accretion rate.
    accretionRate   =max(accretionRate,shi2016AccretionRateMinimum)
    return
  end subroutine shi2016PropertiesGet

  double precision function shi2016RadiusRatio(self,node,massDistribution_) result(radiusRatio)
    !!{RST
    Return the ratio :math:`R_\mathrm{sp}/R_\mathrm{200m}` for the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusShi2016), intent(inout)           :: self
    type            (treeNode                             ), intent(inout)           :: node
    class           (massDistributionClass                ), intent(inout), optional :: massDistribution_
    double precision                                                                 :: accretionRate    , omegaMatter, &
         &                                                                              radius200Mean    , mass200Mean

    call self%propertiesGet(node,accretionRate,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    radiusRatio=shi2016RadiusRatioEvaluate(accretionRate,omegaMatter)
    return
  end function shi2016RadiusRatio

  double precision function shi2016Radius(self,node,massDistribution_) result(radius)
    !!{RST
    Return the splashback radius (in Mpc) of the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusShi2016), intent(inout)           :: self
    type            (treeNode                             ), intent(inout)           :: node
    class           (massDistributionClass                ), intent(inout), optional :: massDistribution_
    double precision                                                                 :: accretionRate    , omegaMatter, &
         &                                                                              radius200Mean    , mass200Mean

    call self%propertiesGet(node,accretionRate,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    radius=+shi2016RadiusRatioEvaluate(accretionRate,omegaMatter) &
         & *radius200Mean
    return
  end function shi2016Radius

  double precision function shi2016MassRatio(self,node,massDistribution_) result(massRatio)
    !!{RST
    Return the ratio :math:`M_\mathrm{sp}/M_\mathrm{200m}` for the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusShi2016), intent(inout)           :: self
    type            (treeNode                             ), intent(inout)           :: node
    class           (massDistributionClass                ), intent(inout), optional :: massDistribution_
    double precision                                                                 :: accretionRate    , omegaMatter, &
         &                                                                              radius200Mean    , mass200Mean

    call self%propertiesGet(node,accretionRate,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    massRatio=+shi2016DensityContrastEvaluate(accretionRate,omegaMatter)         &
         &    /                               splashbackDensityContrastReference &
         &    *shi2016RadiusRatioEvaluate    (accretionRate,omegaMatter)**3
    return
  end function shi2016MassRatio

  double precision function shi2016Mass(self,node,massDistribution_) result(mass)
    !!{RST
    Return the splashback mass (in :math:`\mathrm{M}_\odot`) of the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusShi2016), intent(inout)           :: self
    type            (treeNode                             ), intent(inout)           :: node
    class           (massDistributionClass                ), intent(inout), optional :: massDistribution_
    double precision                                                                 :: accretionRate    , omegaMatter, &
         &                                                                              radius200Mean    , mass200Mean

    call self%propertiesGet(node,accretionRate,omegaMatter,radius200Mean,mass200Mean,massDistribution_)
    mass=+shi2016DensityContrastEvaluate(accretionRate,omegaMatter)         &
         & /                             splashbackDensityContrastReference &
         & *shi2016RadiusRatioEvaluate  (accretionRate,omegaMatter)**3      &
         & *mass200Mean
    return
  end function shi2016Mass

  double precision function shi2016RadiusRatioEvaluate(accretionRate,omegaMatter) result(ratio)
    !!{RST
    Evaluate the :cite:t:`shi_outer_2016` fitting function for :math:`R_\mathrm{sp}/R_\mathrm{200m}`,

    .. math::

       {R_\mathrm{sp} \over R_\mathrm{200m}} = \exp\left[ \left( 0.24 + 0.074 \ln s \right) \ln \Omega_\mathrm{M} + 0.55 - 0.15 s \right],

    where :math:`s` is the instantaneous mass accretion rate.
    !!}
    implicit none
    double precision, intent(in   ) :: accretionRate, omegaMatter

    ratio=exp(                              &
         &    +(                            &
         &      +0.240d0                    &
         &      +0.074d0*log(accretionRate) &
         &     )        *log(omegaMatter  ) &
         &    +0.550d0                      &
         &    -0.150d0  *    accretionRate  &
         &   )
    return
  end function shi2016RadiusRatioEvaluate

  double precision function shi2016DensityContrastEvaluate(accretionRate,omegaMatter) result(densityContrast)
    !!{RST
    Evaluate the :cite:t:`shi_outer_2016` fitting function for the mean overdensity enclosed by the splashback radius, in
    units of the mean density of the universe,

    .. math::

       \Delta_\mathrm{sp} = 33 \Omega_\mathrm{M}^{-0.45} \exp\left[ \left( 0.88 - 0.14 \ln \Omega_\mathrm{M} \right) s^{0.6} \right],

    where :math:`s` is the instantaneous mass accretion rate.
    !!}
    implicit none
    double precision, intent(in   ) :: accretionRate, omegaMatter

    densityContrast=+33.0d0                         &
         &          *omegaMatter**(-0.45d0)         &
         &          *exp(                           &
         &               +(                         &
         &                 +0.88d0                  &
         &                 -0.14d0*log(omegaMatter) &
         &                )                         &
         &               *accretionRate**0.6d0      &
         &              )
    return
  end function shi2016DensityContrastEvaluate
