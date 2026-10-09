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
  Implements a dust screen whose attenuation depends on galaxy stellar mass and redshift.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <dustAttenuation name="dustAttenuationStellarMassRedshift" docformat="rst">
   <description>
   A uniform dust screen whose attenuation at a reference wavelength (by default that of H\ :math:`\alpha`) depends on the
   stellar mass and redshift of the galaxy, following the relation of :cite:t:`garn_predicting_2010` between H\
   :math:`\alpha` attenuation and stellar mass, generalized with free offsets and slopes:

   .. math::

      \bar{A} = A_\mathrm{GB10}(X) + \delta_0 + \delta_\mathrm{M} X + \delta_z u + \delta_{\mathrm{M}z} X u,

   where :math:`X = \log_{10}(M_\star/10^{10}\mathrm{M}_\odot)`, with :math:`M_\star` limited to the range
   [``[massStellarMinimum]``, ``[massStellarMaximum]``], :math:`u = \ln[(1+z)/(1+z_\mathrm{p})]` with :math:`z_\mathrm{p}`
   given by ``[redshiftPivot]``, and

   .. math::

      A_\mathrm{GB10}(X) = 0.91 + 0.77 X + 0.11 X^2 - 0.09 X^3.

   The parameters :math:`\delta_0`, :math:`\delta_\mathrm{M}`, :math:`\delta_z`, and :math:`\delta_{\mathrm{M}z}` are given by
   ``[delta0]``, ``[deltaMass]``, ``[deltaRedshift]``, and ``[deltaMassRedshift]``. The attenuation applied is
   :math:`\max(\bar{A},0)` (in magnitudes) at the wavelength ``[wavelengthReference]``, and scales with wavelength following
   the ``[dustExtinctionCurve]``, so that the same screen may attenuate other emission lines consistently. Any scatter in
   the attenuation about this mean relation is applied separately, by the
   :galacticus-class:`outputAnalysisDistributionOperatorAttenuationScatter` distribution operator.

   The attenuation depends only on the properties of the galaxy as a whole, so is applied equally to all components
   (including a central black hole), and to their combined emission.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationScreen) :: dustAttenuationStellarMassRedshift
     !!{RST
     A dust screen whose attenuation depends on galaxy stellar mass and redshift.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: delta0                       , deltaMass          , &
          &                                                deltaRedshift                , deltaMassRedshift  , &
          &                                                redshiftPivot                , massStellarMinimum , &
          &                                                massStellarMaximum           , wavelengthReference
   contains
     !![
     <methods docformat="rst">
       <method method="attenuationMean" description="Return the mean attenuation, :math:`\bar{A}` (in magnitudes, and not limited to be non-negative), at the reference wavelength for the given node."/>
     </methods>
     !!]
     final     ::                      stellarMassRedshiftDestructor
     procedure :: depthOpticalV     => stellarMassRedshiftDepthOpticalV
     procedure :: attenuationMean   => stellarMassRedshiftAttenuationMean
     procedure :: supportsComponent => stellarMassRedshiftSupportsComponent
     procedure :: request           => stellarMassRedshiftRequest
  end type dustAttenuationStellarMassRedshift

  interface dustAttenuationStellarMassRedshift
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationStellarMassRedshift` dust attenuation class.
     !!}
     module procedure stellarMassRedshiftConstructorParameters
     module procedure stellarMassRedshiftConstructorInternal
  end interface dustAttenuationStellarMassRedshift

contains

  function stellarMassRedshiftConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationStellarMassRedshift` dust attenuation class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustAttenuationStellarMassRedshift)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (dustExtinctionCurveClass          ), pointer       :: dustExtinctionCurve_
    class           (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    double precision                                                    :: delta0              , deltaMass          , &
         &                                                                 deltaRedshift       , deltaMassRedshift  , &
         &                                                                 redshiftPivot       , massStellarMinimum , &
         &                                                                 massStellarMaximum  , wavelengthReference

    !![
    <inputParameter docformat="rst">
      <name>delta0</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The offset, :math:`\delta_0`, in magnitudes, of the attenuation from the relation of :cite:t:`garn_predicting_2010`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>deltaMass</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The change, :math:`\delta_\mathrm{M}`, in the slope of the attenuation with :math:`\log_{10}` stellar mass, in magnitudes
      per dex.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>deltaRedshift</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The slope, :math:`\delta_z`, of the attenuation with :math:`u = \ln[(1+z)/(1+z_\mathrm{p})]`, in magnitudes.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>deltaMassRedshift</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The coefficient, :math:`\delta_{\mathrm{M}z}`, of the product :math:`X u` of log stellar mass and redshift, in magnitudes.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshiftPivot</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The pivot redshift, :math:`z_\mathrm{p}`, of the redshift dependence.
      </description>
      <source>parameters</source>
      <minimum>0.0d0</minimum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massStellarMinimum</name>
      <defaultValue>1.0d8</defaultValue>
      <description>
      The stellar mass (in :math:`\mathrm{M}_\odot`) below which the attenuation is evaluated at this mass.
      </description>
      <source>parameters</source>
      <minimum inclusive="false">0.0d0</minimum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massStellarMaximum</name>
      <defaultValue>1.0d11</defaultValue>
      <description>
      The stellar mass (in :math:`\mathrm{M}_\odot`) above which the attenuation is evaluated at this mass.
      </description>
      <source>parameters</source>
      <minimum inclusive="false">0.0d0</minimum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavelengthReference</name>
      <defaultValue>6564.61d0</defaultValue>
      <description>
      The wavelength (in Å) at which the attenuation is specified. The default is that of H\ :math:`\alpha` (in vacuum).
      </description>
      <source>parameters</source>
      <minimum inclusive="false">0.0d0</minimum>
    </inputParameter>
    <objectBuilder class="dustExtinctionCurve" name="dustExtinctionCurve_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=dustAttenuationStellarMassRedshift(delta0,deltaMass,deltaRedshift,deltaMassRedshift,redshiftPivot,massStellarMinimum,massStellarMaximum,wavelengthReference,dustExtinctionCurve_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustExtinctionCurve_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function stellarMassRedshiftConstructorParameters

  function stellarMassRedshiftConstructorInternal(delta0,deltaMass,deltaRedshift,deltaMassRedshift,redshiftPivot,massStellarMinimum,massStellarMaximum,wavelengthReference,dustExtinctionCurve_,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationStellarMassRedshift` dust attenuation class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (dustAttenuationStellarMassRedshift)                        :: self
    double precision                                    , intent(in   )         :: delta0              , deltaMass          , &
         &                                                                         deltaRedshift       , deltaMassRedshift  , &
         &                                                                         redshiftPivot       , massStellarMinimum , &
         &                                                                         massStellarMaximum  , wavelengthReference
    class           (dustExtinctionCurveClass          ), intent(in   ), target :: dustExtinctionCurve_
    class           (cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="delta0, deltaMass, deltaRedshift, deltaMassRedshift, redshiftPivot, massStellarMinimum, massStellarMaximum, wavelengthReference, *dustExtinctionCurve_, *cosmologyFunctions_"/>
    !!]

    if (massStellarMinimum >= massStellarMaximum) call Error_Report('`massStellarMinimum` < `massStellarMaximum` is required'//{introspection:location})
    return
  end function stellarMassRedshiftConstructorInternal

  subroutine stellarMassRedshiftDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationStellarMassRedshift` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationStellarMassRedshift), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine stellarMassRedshiftDestructor

  double precision function stellarMassRedshiftAttenuationMean(self,node) result(attenuation)
    !!{RST
    Return the mean attenuation, :math:`\bar{A}` (in magnitudes, and not limited to be non-negative), at the reference
    wavelength for the given node.
    !!}
    use :: Galactic_Structure_Options, only : massTypeStellar
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (dustAttenuationStellarMassRedshift), intent(inout)         :: self
    type            (treeNode                          ), intent(inout), target :: node
    class           (massDistributionClass             ), pointer               :: massDistribution_
    class           (nodeComponentBasic                ), pointer               :: basic
    double precision                                                            :: massStellar      , massLogarithmic, &
         &                                                                         redshift         , redshiftTerm

    ! Find the stellar mass, limited to the range over which the relation is applied.
    massDistribution_ => node             %massDistribution(massType=massTypeStellar)
    massStellar       =  massDistribution_%massTotal       (                        )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    massStellar    =min(max(massStellar,self%massStellarMinimum),self%massStellarMaximum)
    massLogarithmic=log10(massStellar/1.0d10)
    ! Find the redshift of the node.
    basic        => node%basic()
    redshift     =  self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(basic%time()))
    redshiftTerm =  log((1.0d0+redshift)/(1.0d0+self%redshiftPivot))
    ! Evaluate the mean attenuation.
    attenuation=+0.91d0                                              &
         &      +0.77d0*massLogarithmic                              &
         &      +0.11d0*massLogarithmic**2                           &
         &      -0.09d0*massLogarithmic**3                           &
         &      +self%delta0                                         &
         &      +self%deltaMass        *massLogarithmic              &
         &      +self%deltaRedshift                    *redshiftTerm &
         &      +self%deltaMassRedshift*massLogarithmic*redshiftTerm
    return
  end function stellarMassRedshiftAttenuationMean

  double precision function stellarMassRedshiftDepthOpticalV(self,node,componentType) result(depthOpticalV)
    !!{RST
    Return the :math:`V`-band optical depth of the screen: that which gives an attenuation of :math:`\max(\bar{A},0)`
    magnitudes at the reference wavelength, given the extinction curve.
    !!}
    implicit none
    class(dustAttenuationStellarMassRedshift), intent(inout)         :: self
    type (treeNode                          ), intent(inout), target :: node
    type (enumerationComponentTypeType      ), intent(in   )         :: componentType
    !$GLC attributes unused :: componentType

    ! An attenuation of A magnitudes is an optical depth of A/(2.5 log10 e) = 0.4 A ln 10.
    depthOpticalV=+max(self%attenuationMean(node),0.0d0)                                   &
         &        *0.4d0                                                                   &
         &        *log(10.0d0)                                                             &
         &        /self%dustExtinctionCurve_%attenuationRelative(self%wavelengthReference)
    return
  end function stellarMassRedshiftDepthOpticalV

  logical function stellarMassRedshiftSupportsComponent(self,componentType)
    !!{RST
    The attenuation depends only on the galaxy as a whole, so may be applied to any component, and to the sum of all
    components.
    !!}
    implicit none
    class(dustAttenuationStellarMassRedshift), intent(inout) :: self
    type (enumerationComponentTypeType      ), intent(in   ) :: componentType
    !$GLC attributes unused :: self, componentType

    stellarMassRedshiftSupportsComponent=.true.
    return
  end function stellarMassRedshiftSupportsComponent

  function stellarMassRedshiftRequest(self) result(request)
    !!{RST
    The attenuation depends only on the galaxy as a whole, so emission need not be decomposed by component, age,
    metallicity, or radius.
    !!}
    implicit none
    type (decompositionRequest              )                :: request
    class(dustAttenuationStellarMassRedshift), intent(inout) :: self
    !$GLC attributes unused :: self

    request%resolveComponents =.false.
    request%resolveMetallicity=.false.
    request%resolveRadius     =.false.
    return
  end function stellarMassRedshiftRequest
