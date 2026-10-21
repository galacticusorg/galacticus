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
Implements an N-body dark matter halo mass error class which models errors arising from merger tree construction.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <nbodyHaloMassError name="nbodyHaloMassErrorTreeConstruction" docformat="rst">
   <description>
   An N-body dark matter halo mass error class which models the uncertainty in the masses of progenitor halos arising from errors in the construction of merger trees---for example, halos being incorrectly linked into a tree as a result of close passages with other halos, or other halo finder failures. Such errors accumulate as a halo is traced further back through a tree, so the error is defined relative to a reference halo (typically the parent halo) from which the halo was reached. The fractional mass error is :math:`\sigma(M) = \sigma_\mathrm{t} (M/10^{12}\mathrm{M}_\odot)^\gamma \Delta^{\beta/2}`, where :math:`\Delta = \ln(a_\mathrm{reference}/a) = \ln[(1+z)/(1+z_\mathrm{reference})]` is the time lag between the halo and the reference halo, :math:`\sigma_\mathrm{t}=`\ ``[normalization]``, :math:`\gamma=`\ ``[exponentMass]``, and :math:`\beta=`\ ``[exponentTimeLag]``. For :math:`\beta=1` the variance grows linearly with the number of snapshots traversed (for snapshots uniformly spaced in :math:`\ln a`), as expected if each link in the tree contributes an independent error. The error is zero if no reference halo is given, or if the halo is not earlier than the reference halo. Errors are assumed to be uncorrelated between halos. Tree construction errors may also bias the masses of progenitor halos---for example, if poorly-resolved progenitors are systematically assigned too little mass, or the main branch is lost to a less massive halo. The mean measured mass is then :math:`(1+\delta) M`, with :math:`\delta = \delta_\mathrm{max} \{ 1 - \exp[-(N/N_0)^\alpha \Delta^{\beta_\delta}] \}`, where :math:`N=M/m_\mathrm{p}` is the number of particles in the halo, :math:`m_\mathrm{p}=`\ ``[massParticle]``, :math:`\delta_\mathrm{max}=`\ ``[biasMaximum]``, :math:`N_0=`\ ``[biasParticleNumber]``, :math:`\alpha=`\ ``[biasExponentParticleNumber]``, and :math:`\beta_\delta=`\ ``[biasExponentTimeLag]``. For :math:`\alpha &lt; 0` the bias saturates at :math:`\delta_\mathrm{max}` for poorly-resolved halos and vanishes for well-resolved ones, with the transition moving to larger :math:`N` as the time lag increases. The bias is zero if no reference halo is given, if the halo is not earlier than the reference halo, or if :math:`\delta_\mathrm{max}=0` (the default). This class is intended to be combined with a model for particle sampling errors (e.g. :galacticus-class:`nbodyHaloMassErrorPowerLaw`) using :galacticus-class:`nbodyHaloMassErrorSummation`.
   </description>
  </nbodyHaloMassError>
  !!]
  type, extends(nbodyHaloMassErrorClass) :: nbodyHaloMassErrorTreeConstruction
     !!{RST
     An N-body halo mass error class which models errors arising from merger tree construction.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: normalization                , exponentMass              , &
          &                                                exponentTimeLag              , biasMaximum               , &
          &                                                biasParticleNumber           , biasExponentParticleNumber, &
          &                                                biasExponentTimeLag          , massParticle
   contains
     final     ::                    treeConstructionDestructor
     procedure :: errorFractional => treeConstructionErrorFractional
     procedure :: biasFractional  => treeConstructionBiasFractional
     procedure :: correlation     => treeConstructionCorrelation
  end type nbodyHaloMassErrorTreeConstruction

  interface nbodyHaloMassErrorTreeConstruction
     !!{RST
     Constructors for the :galacticus-class:`nbodyHaloMassErrorTreeConstruction` N-body halo mass error class.
     !!}
     module procedure treeConstructionConstructorParameters
     module procedure treeConstructionConstructorInternal
  end interface nbodyHaloMassErrorTreeConstruction

contains

  function treeConstructionConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nbodyHaloMassErrorTreeConstruction` N-body halo mass error class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyHaloMassErrorTreeConstruction)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    double precision                                                    :: normalization      , exponentMass              , &
         &                                                                 exponentTimeLag    , biasMaximum               , &
         &                                                                 biasParticleNumber , biasExponentParticleNumber, &
         &                                                                 biasExponentTimeLag, massParticle

    ! Check and read parameters.
    !![
    <inputParameter docformat="rst">
      <name>normalization</name>
      <source>parameters</source>
      <description>
      Parameter :math:`\sigma_\mathrm{t}` appearing in the model for merger tree construction errors: the fractional mass error for a halo of mass :math:`10^{12}\mathrm{M}_\odot` at a time lag of :math:`\Delta=1` from the reference halo.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentMass</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      Parameter :math:`\gamma` appearing in the model for merger tree construction errors: the exponent of halo mass in the fractional mass error.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentTimeLag</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>
      Parameter :math:`\beta` appearing in the model for merger tree construction errors: the exponent of the time lag, :math:`\Delta=\ln(a_\mathrm{reference}/a)`, in the variance of the fractional mass error.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>biasMaximum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      Parameter :math:`\delta_\mathrm{max}` appearing in the model for merger tree construction errors: the maximum fractional bias in the mass of a progenitor halo (negative values correspond to masses which are underestimated). Must be greater than :math:`-1`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>biasParticleNumber</name>
      <source>parameters</source>
      <defaultValue>1000.0d0</defaultValue>
      <description>
      Parameter :math:`N_0` appearing in the model for merger tree construction errors: the particle number scale of the fractional mass bias.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>biasExponentParticleNumber</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>
      Parameter :math:`\alpha` appearing in the model for merger tree construction errors: the exponent of particle number in the fractional mass bias.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>biasExponentTimeLag</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>
      Parameter :math:`\beta_\delta` appearing in the model for merger tree construction errors: the exponent of the time lag, :math:`\Delta=\ln(a_\mathrm{reference}/a)`, in the fractional mass bias.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massParticle</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The mass of a particle, :math:`m_\mathrm{p}`, in the N-body simulation, used to compute the number of particles in a halo for the fractional mass bias. Required if :math:`\delta_\mathrm{max} \ne 0`.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nbodyHaloMassErrorTreeConstruction(normalization,exponentMass,exponentTimeLag,biasMaximum,biasParticleNumber,biasExponentParticleNumber,biasExponentTimeLag,massParticle,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function treeConstructionConstructorParameters

  function treeConstructionConstructorInternal(normalization,exponentMass,exponentTimeLag,biasMaximum,biasParticleNumber,biasExponentParticleNumber,biasExponentTimeLag,massParticle,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nbodyHaloMassErrorTreeConstruction` N-body halo mass error class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nbodyHaloMassErrorTreeConstruction)                        :: self
    double precision                                    , intent(in   )         :: normalization      , exponentMass              , &
         &                                                                         exponentTimeLag    , biasMaximum               , &
         &                                                                         biasParticleNumber , biasExponentParticleNumber, &
         &                                                                         biasExponentTimeLag, massParticle
    class           (cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="normalization, exponentMass, exponentTimeLag, biasMaximum, biasParticleNumber, biasExponentParticleNumber, biasExponentTimeLag, massParticle, *cosmologyFunctions_"/>
    !!]

    if (normalization      <   0.0d0                            ) call Error_Report('normalization must be non-negative'                      //{introspection:location})
    if (biasMaximum        <= -1.0d0                            ) call Error_Report('biasMaximum must be greater than -1'                     //{introspection:location})
    if (biasParticleNumber <=  0.0d0                            ) call Error_Report('biasParticleNumber must be positive'                     //{introspection:location})
    if (biasMaximum        /=  0.0d0 .and. massParticle <= 0.0d0) call Error_Report('massParticle must be positive if biasMaximum is non-zero'//{introspection:location})
    return
  end function treeConstructionConstructorInternal

  subroutine treeConstructionDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nbodyHaloMassErrorTreeConstruction` N-body halo mass error class.
    !!}
    implicit none
    type(nbodyHaloMassErrorTreeConstruction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine treeConstructionDestructor

  double precision function treeConstructionErrorFractional(self,node,nodeReference)
    !!{RST
    Return the fractional error on the mass of an N-body halo in the merger tree construction error model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nbodyHaloMassErrorTreeConstruction), intent(inout)           :: self
    type            (treeNode                          ), intent(inout)           :: node
    type            (treeNode                          ), intent(inout), optional :: nodeReference
    class           (nodeComponentBasic                ), pointer                 :: basic
    class           (nodeComponentBasic                ), pointer                 :: basicReference
    double precision                                    , parameter               :: massNormalization=1.0d12
    double precision                                                              :: timeLag

    ! Errors arise only for halos reached by tracing back through a tree from a later reference halo.
    treeConstructionErrorFractional=0.0d0
    if (.not.present(nodeReference)) return
    basic          => node         %basic()
    basicReference => nodeReference%basic()
    if (basic%time() >= basicReference%time()) return
    ! Evaluate the time lag between the halo and the reference halo.
    timeLag                        =+log(                                                                  &
         &                               +self%cosmologyFunctions_%expansionFactor(basicReference%time())  &
         &                               /self%cosmologyFunctions_%expansionFactor(basic         %time())  &
         &                              )
    treeConstructionErrorFractional=+self%normalization                                             &
         &                          *(basic%mass()/massNormalization)**       self%exponentMass     &
         &                          * timeLag                        **(0.5d0*self%exponentTimeLag)
    return
  end function treeConstructionErrorFractional

  double precision function treeConstructionBiasFractional(self,node,nodeReference)
    !!{RST
    Return the fractional bias in the mass of an N-body halo in the merger tree construction error model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nbodyHaloMassErrorTreeConstruction), intent(inout)           :: self
    type            (treeNode                          ), intent(inout)           :: node
    type            (treeNode                          ), intent(inout), optional :: nodeReference
    class           (nodeComponentBasic                ), pointer                 :: basic
    class           (nodeComponentBasic                ), pointer                 :: basicReference
    double precision                                                              :: timeLag

    ! Bias arises only for halos reached by tracing back through a tree from a later reference halo.
    treeConstructionBiasFractional=0.0d0
    if (self%biasMaximum == 0.0d0 .or. .not.present(nodeReference)) return
    basic          => node         %basic()
    basicReference => nodeReference%basic()
    if (basic%time() >= basicReference%time()) return
    ! Evaluate the time lag between the halo and the reference halo.
    timeLag                       =+log(                                                                 &
         &                              +self%cosmologyFunctions_%expansionFactor(basicReference%time()) &
         &                              /self%cosmologyFunctions_%expansionFactor(basic         %time()) &
         &                             )
    treeConstructionBiasFractional=+self%biasMaximum                                                &
         &                         *(                                                               &
         &                           +1.0d0                                                         &
         &                           -exp(                                                          &
         &                                -(basic%mass()/self%massParticle/self%biasParticleNumber) &
         &                                 **self%biasExponentParticleNumber                        &
         &                                *timeLag                                                  &
         &                                 **self%biasExponentTimeLag                               &
         &                               )                                                          &
         &                          )
    return
  end function treeConstructionBiasFractional

  double precision function treeConstructionCorrelation(self,node1,node2)
    !!{RST
    Return the correlation of the masses of a pair of N-body halos. Merger tree construction errors are assumed to be uncorrelated
    between halos.
    !!}
    implicit none
    class(nbodyHaloMassErrorTreeConstruction), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node1, node2
    !$GLC attributes unused :: self, node1, node2

    treeConstructionCorrelation=0.0d0
    return
  end function treeConstructionCorrelation
