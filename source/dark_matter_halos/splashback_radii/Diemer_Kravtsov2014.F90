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
  An implementation of dark matter halo splashback radii using the truncation radius of :cite:t:`diemer_dependence_2014`.
  !!}

  !+ Contributions to this file made by: Andrew Benson, Claude.

  use :: Cosmology_Functions                        , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                       , only : cosmologyParametersClass
  use :: Cosmological_Density_Field                 , only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Virial_Density_Contrast                    , only : virialDensityContrastClass
  use :: Dark_Matter_Halo_Splashback_Radii_Reference, only : splashbackReferenceScales

  !![
  <darkMatterHaloSplashbackRadius name="darkMatterHaloSplashbackRadiusDiemerKravtsov2014" docformat="rst">
   <description>
   A splashback radius class which uses the truncation radius of the outer density profile of
   :cite:t:`diemer_dependence_2014`, given by their equation (6),

   .. math::

      {R_\mathrm{t} \over R_\mathrm{200m}} = \alpha - \beta \nu,

   where :math:`\alpha=`\ ``[alpha]``, :math:`\beta=`\ ``[beta]``, and :math:`\nu` is the peak height computed from
   :math:`M_\mathrm{200m}`.

   **Note that this is not a splashback radius.** :math:`R_\mathrm{t}` is a fitting parameter of the
   :cite:t:`diemer_dependence_2014` outer density profile which controls where that profile steepens; it correlates with, but
   is not equal to, the splashback radius. It is provided by this class both because it is the quantity historically used by
   the accretion flow mass distribution classes in Galacticus, and because it is a useful point of comparison for the models
   which do predict a true splashback radius. For physical applications prefer
   :galacticus-class:`darkMatterHaloSplashbackRadiusDiemer2020`.

   Since this model makes no prediction for the splashback mass, or for the scatter in either quantity, attempting to
   evaluate those will result in an error.
   </description>
  </darkMatterHaloSplashbackRadius>
  !!]
  type, extends(darkMatterHaloSplashbackRadiusClass) :: darkMatterHaloSplashbackRadiusDiemerKravtsov2014
     !!{RST
     A dark matter halo splashback radius class using the truncation radius of :cite:t:`diemer_dependence_2014`.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_      => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (virialDensityContrastClass   ), pointer :: virialDensityContrast_    => null()
     double precision                                         :: alpha                              , beta
   contains
     !![
     <methods docformat="rst">
       <method method="propertiesGet" description="Return the properties of the halo needed to evaluate the truncation radius."/>
     </methods>
     !!]
     final     ::                  diemerKravtsov2014Destructor
     procedure :: radius        => diemerKravtsov2014Radius
     procedure :: radiusRatio   => diemerKravtsov2014RadiusRatio
     procedure :: propertiesGet => diemerKravtsov2014PropertiesGet
  end type darkMatterHaloSplashbackRadiusDiemerKravtsov2014

  interface darkMatterHaloSplashbackRadiusDiemerKravtsov2014
     !!{RST
     Constructors for the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemerKravtsov2014` splashback radius class.
     !!}
     module procedure diemerKravtsov2014ConstructorParameters
     module procedure diemerKravtsov2014ConstructorInternal
  end interface darkMatterHaloSplashbackRadiusDiemerKravtsov2014

contains

  function diemerKravtsov2014ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemerKravtsov2014` splashback radius class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterHaloSplashbackRadiusDiemerKravtsov2014)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                         ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass                        ), pointer       :: cosmologyParameters_
    class           (criticalOverdensityClass                        ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                   ), pointer       :: cosmologicalMassVariance_
    class           (virialDensityContrastClass                      ), pointer       :: virialDensityContrast_
    double precision                                                                  :: alpha                    , beta

    !![
    <inputParameter docformat="rst">
      <name>alpha</name>
      <defaultValue>1.90d0</defaultValue>
      <defaultSource>:cite:t:`diemer_dependence_2014`, equation (6).</defaultSource>
      <description>
      The parameter :math:`\alpha` in the truncation radius :math:`R_\mathrm{t}/R_\mathrm{200m} = \alpha - \beta \nu` of
      :cite:t:`diemer_dependence_2014`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>beta</name>
      <defaultValue>0.18d0</defaultValue>
      <defaultSource>:cite:t:`diemer_dependence_2014`, equation (6).</defaultSource>
      <description>
      The parameter :math:`\beta` in the truncation radius :math:`R_\mathrm{t}/R_\mathrm{200m} = \alpha - \beta \nu` of
      :cite:t:`diemer_dependence_2014`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="virialDensityContrast"    name="virialDensityContrast_"    source="parameters"/>
    !!]
    self=darkMatterHaloSplashbackRadiusDiemerKravtsov2014(alpha,beta,cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="virialDensityContrast_"   />
    !!]
    return
  end function diemerKravtsov2014ConstructorParameters

  function diemerKravtsov2014ConstructorInternal(alpha,beta,cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemerKravtsov2014` splashback radius class.
    !!}
    implicit none
    type            (darkMatterHaloSplashbackRadiusDiemerKravtsov2014)                        :: self
    double precision                                                  , intent(in   )         :: alpha                    , beta
    class           (cosmologyFunctionsClass                         ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass                        ), intent(in   ), target :: cosmologyParameters_
    class           (criticalOverdensityClass                        ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                   ), intent(in   ), target :: cosmologicalMassVariance_
    class           (virialDensityContrastClass                      ), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="alpha, beta, *cosmologyFunctions_, *cosmologyParameters_, *criticalOverdensity_, *cosmologicalMassVariance_, *virialDensityContrast_"/>
    !!]

    return
  end function diemerKravtsov2014ConstructorInternal

  subroutine diemerKravtsov2014Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`darkMatterHaloSplashbackRadiusDiemerKravtsov2014` splashback radius class.
    !!}
    implicit none
    type(darkMatterHaloSplashbackRadiusDiemerKravtsov2014), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%virialDensityContrast_"   />
    !!]
    return
  end subroutine diemerKravtsov2014Destructor

  subroutine diemerKravtsov2014PropertiesGet(self,node,peakHeight,radius200Mean,massDistribution_)
    !!{RST
    Return the peak height, :math:`\nu`, computed from :math:`M_\mathrm{200m}`, and :math:`R_\mathrm{200m}` itself, for the
    halo in ``node``.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemerKravtsov2014), intent(inout)           :: self
    type            (treeNode                                        ), intent(inout)           :: node
    double precision                                                  , intent(  out)           :: peakHeight       , radius200Mean
    class           (massDistributionClass                           ), intent(inout), optional :: massDistribution_
    class           (nodeComponentBasic                              ), pointer                 :: basic
    double precision                                                                            :: time             , rootVariance , &
         &                                                                                         mass200Mean

    basic       =>  node %basic()
    time        =   basic%time ()
    call splashbackReferenceScales(                                                    &
         &                                                node                        , &
         &                         cosmologyParameters_  =self%cosmologyParameters_   , &
         &                         cosmologyFunctions_   =self%cosmologyFunctions_    , &
         &                         virialDensityContrast_=self%virialDensityContrast_ , &
         &                         mass                  =     mass200Mean            , &
         &                         radius                =     radius200Mean          , &
         &                         massDistribution_     =     massDistribution_        &
         &                        )
    rootVariance=self%cosmologicalMassVariance_%rootVariance(time=time,mass=mass200Mean          )
    if (rootVariance > 0.0d0) then
       peakHeight=+self%criticalOverdensity_   %value       (time=time,mass=mass200Mean,node=node) &
            &     /                             rootVariance
    else
       peakHeight=+0.0d0
    end if
    return
  end subroutine diemerKravtsov2014PropertiesGet

  double precision function diemerKravtsov2014RadiusRatio(self,node,massDistribution_) result(radiusRatio)
    !!{RST
    Return the ratio :math:`R_\mathrm{t}/R_\mathrm{200m}` for the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemerKravtsov2014), intent(inout)           :: self
    type            (treeNode                                        ), intent(inout)           :: node
    class           (massDistributionClass                           ), intent(inout), optional :: massDistribution_
    double precision                                                                            :: peakHeight       , radius200Mean

    call self%propertiesGet(node,peakHeight,radius200Mean,massDistribution_)
    radiusRatio=+self%alpha            &
         &      -self%beta *peakHeight
    return
  end function diemerKravtsov2014RadiusRatio

  double precision function diemerKravtsov2014Radius(self,node,massDistribution_) result(radius)
    !!{RST
    Return the truncation radius (in Mpc) of the halo in ``node``.
    !!}
    implicit none
    class           (darkMatterHaloSplashbackRadiusDiemerKravtsov2014), intent(inout)           :: self
    type            (treeNode                                        ), intent(inout)           :: node
    class           (massDistributionClass                           ), intent(inout), optional :: massDistribution_
    double precision                                                                            :: peakHeight       , radius200Mean

    call self%propertiesGet(node,peakHeight,radius200Mean,massDistribution_)
    radius=+(                       &
         &   +self%alpha            &
         &   -self%beta *peakHeight &
         &  )                       &
         & *radius200Mean
    return
  end function diemerKravtsov2014Radius
