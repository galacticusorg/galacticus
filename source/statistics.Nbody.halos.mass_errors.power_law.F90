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
Implements an N-body dark matter halo mass error class in which errors are a power-law in halo mass.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <nbodyHaloMassError name="nbodyHaloMassErrorPowerLaw" docformat="rst">
   <description>
   An N-body dark matter halo mass error class in which the fractional mass uncertainty is modeled as a power-law function of halo mass, useful for parameterizing simulation resolution effects. The model is characterized by the normalization parameters :math:`\sigma_{12}` and :math:`\sigma_\infty`, and the power-law exponent with respect to mass. Optionally, correlations between mass errors of pairs of halos may be modeled as a power-law in the ratio of halo masses and expansion factors: :math:`C_{12} = C_0 [M_2/M_1]^\alpha [(1+z_2)/(1+z_1)]^\beta`. If this option is disabled (the default), a trivial correlation model is used in which the correlation is unity for halos with identical mass and time, and zero otherwise.
   </description>
  </nbodyHaloMassError>
  !!]
  type, extends(nbodyHaloMassErrorClass) :: nbodyHaloMassErrorPowerLaw
     !!{RST
     An N-body halo mass error class in which errors are a power-law in halo mass.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_            => null()
     double precision                                   :: normalizationSquared                    , exponent                   , &
          &                                                fractionalErrorHighMassSquared          , normalization              , &
          &                                                fractionalErrorHighMass                 , correlationNormalization   , &
          &                                                correlationMassExponent                 , correlationRedshiftExponent
     logical                                            :: correlationModelTrivial
   contains
     final     ::                    powerLawDestructor
     procedure :: errorFractional => powerLawErrorFractional
     procedure :: correlation     => powerLawCorrelation
  end type nbodyHaloMassErrorPowerLaw

  interface nbodyHaloMassErrorPowerLaw
     !!{RST
     Constructors for the :galacticus-class:`nbodyHaloMassErrorPowerLaw` N-body halo mass error class.
     !!}
     module procedure nbodyHaloMassErrorPowerLawConstructorParameters
     module procedure nbodyHaloMassErrorPowerLawConstructorInternal
  end interface nbodyHaloMassErrorPowerLaw

  ! Reference mass used in the error model.
  double precision :: massReference=1.0d12

contains

  function nbodyHaloMassErrorPowerLawConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nbodyHaloMassErrorPowerLaw` N-body halo mass error class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyHaloMassErrorPowerLaw)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    double precision                                            :: normalization              , fractionalErrorHighMass    , &
         &                                                         exponent                   , correlationNormalization   , &
         &                                                         correlationMassExponent    , correlationRedshiftExponent
    logical                                                     :: correlationModelTrivial

    ! Check and read parameters.
    !![
    <inputParameter docformat="rst">
      <name>normalization</name>
      <source>parameters</source>
      <description>
      Parameter :math:`\sigma_{12}` appearing in model for random errors in the halo mass function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>fractionalErrorHighMass</name>
      <source>parameters</source>
      <description>
      Parameter :math:`\sigma_\infty` appearing in model for random errors in the halo mass function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponent</name>
      <source>parameters</source>
      <description>
      Parameter :math:`\gamma` appearing in model for random errors in the halo mass function. Specifically, the fractional error is given by

      .. math::

         \sigma(M) = \left[ \sigma^2_{12} \left({M_\mathrm{halo} \over 10^{12}\mathrm{M}_\odot}\right)^{2\gamma} + \sigma^2_\infty \right]^{1/2},

      where :math:`\sigma_{12}=`\ ``[normalization]``, and :math:`\gamma=`\ ``[exponent]``.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>correlationModelTrivial</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, the correlation between mass errors of pairs of halos is unity for halos with identical mass and time, and zero otherwise. If false, a power-law correlation model in mass ratio and expansion factor ratio is used instead.
      </description>
    </inputParameter>
    !!]
    if (correlationModelTrivial) then
       self=nbodyHaloMassErrorPowerLaw(normalization,exponent,fractionalErrorHighMass,correlationModelTrivial)
    else
       !![
       <inputParameter docformat="rst">
         <name>correlationNormalization</name>
         <source>parameters</source>
         <defaultValue>0.0d0</defaultValue>
         <description>
         Variable :math:`C_0` in the model for the correlation between halo mass errors: :math:`C_{12} = C_0 [M_2/M_1]^\alpha [(1+z_2)/(1+z_1)]^\beta`.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>correlationMassExponent</name>
         <source>parameters</source>
         <defaultValue>0.0d0</defaultValue>
         <description>
         Variable :math:`\alpha` in the model for the correlation between halo mass errors: :math:`C_{12} = C_0 [M_2/M_1]^\alpha [(1+z_2)/(1+z_1)]^\beta`.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>correlationRedshiftExponent</name>
         <source>parameters</source>
         <defaultValue>0.0d0</defaultValue>
         <description>
         Variable :math:`\beta` in the model for the correlation between halo mass errors: :math:`C_{12} = C_0 [M_2/M_1]^\alpha [(1+z_2)/(1+z_1)]^\beta`.
         </description>
       </inputParameter>
       <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
       !!]
       self=nbodyHaloMassErrorPowerLaw(normalization,exponent,fractionalErrorHighMass,correlationModelTrivial,correlationNormalization,correlationMassExponent,correlationRedshiftExponent,cosmologyFunctions_)
       !![
       <objectDestructor name="cosmologyFunctions_"/>
       !!]
    end if
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nbodyHaloMassErrorPowerLawConstructorParameters

  function nbodyHaloMassErrorPowerLawConstructorInternal(normalization,exponent,fractionalErrorHighMass,correlationModelTrivial,correlationNormalization,correlationMassExponent,correlationRedshiftExponent,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nbodyHaloMassErrorPowerLaw` N-body halo mass error class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nbodyHaloMassErrorPowerLaw)                                  :: self
    double precision                            , intent(in   )                   :: normalization              , exponent                   , &
         &                                                                           fractionalErrorHighMass
    logical                                     , intent(in   )                   :: correlationModelTrivial
    double precision                            , intent(in   ), optional         :: correlationNormalization   , correlationMassExponent    , &
         &                                                                           correlationRedshiftExponent
    class           (cosmologyFunctionsClass   ), intent(in   ), optional, target :: cosmologyFunctions_
    !![
    <constructorAssign variables="normalization, exponent, fractionalErrorHighMass, correlationModelTrivial"/>
    !!]

    self%normalizationSquared          =normalization          **2
    self%fractionalErrorHighMassSquared=fractionalErrorHighMass**2
    if (correlationModelTrivial) then
       self%correlationNormalization    =  0.0d0
       self%correlationMassExponent     =  0.0d0
       self%correlationRedshiftExponent =  0.0d0
       self%cosmologyFunctions_         => null()
    else
       if (.not.(present(correlationNormalization).and.present(correlationMassExponent).and.present(correlationRedshiftExponent))) &
            & call Error_Report('all parameters of correlation model must be provided'//{introspection:location})
       if (.not.present(cosmologyFunctions_)) &
            & call Error_Report('cosmology functions must be provided for correlation model'//{introspection:location})
       self%correlationNormalization    =  correlationNormalization
       self%correlationMassExponent     =  correlationMassExponent
       self%correlationRedshiftExponent =  correlationRedshiftExponent
       self%cosmologyFunctions_         => cosmologyFunctions_
       !![
       <referenceCountIncrement owner="self" isResult="yes" object="cosmologyFunctions_"/>
       !!]
    end if
    return
  end function nbodyHaloMassErrorPowerLawConstructorInternal

  subroutine powerLawDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nbodyHaloMassErrorPowerLaw` N-body halo mass error class.
    !!}
    implicit none
    type(nbodyHaloMassErrorPowerLaw), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine powerLawDestructor

  double precision function powerLawErrorFractional(self,node)
    !!{RST
    Return the fractional error on the mass of an N-body halo in the power-law error model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nbodyHaloMassErrorPowerLaw), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    class           (nodeComponentBasic        ), pointer       :: basic

    basic                   =>  node%basic        ()
    powerLawErrorFractional =   sqrt(                                     &
         &                           +self%normalizationSquared           &
         &                           *(                                   &
         &                             +basic%mass         ()             &
         &                             /      massReference               &
         &                            )**(2.0d0*self%exponent)            &
         &                           +self%fractionalErrorHighMassSquared &
         &                          )
    return
  end function powerLawErrorFractional

  double precision function powerLawCorrelation(self,node1,node2)
    !!{RST
    Return the correlation of the masses of a pair of N-body halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nbodyHaloMassErrorPowerLaw), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node1    , node2
    class           (nodeComponentBasic        ), pointer       :: basic1   , basic2
    double precision                                            :: massRatio, expansionFactorRatio

    basic1 => node1%basic()
    basic2 => node2%basic()
    if (self%correlationModelTrivial) then
       if     (                                &
            &   basic1%mass() == basic2%mass() &
            &  .and.                           &
            &   basic1%time() == basic2%time() &
            & ) then
          powerLawCorrelation=1.0d0
       else
          powerLawCorrelation=0.0d0
       end if
    else
       ! Extract mass and expansion factor ratios.
       massRatio            =  +                                         basic2%mass ()  &
            &                  /                                         basic1%mass ()
       expansionFactorRatio =  +self%cosmologyFunctions_%expansionFactor(basic2%time ()) &
            &                  /self%cosmologyFunctions_%expansionFactor(basic1%time ())
       ! Ensure ratios are below unity, invert otherwise.
       if (           massRatio > 1.0d0)            massRatio=1.0d0/           massRatio
       if (expansionFactorRatio > 1.0d0) expansionFactorRatio=1.0d0/expansionFactorRatio
       ! Evaluate the correlation.
       powerLawCorrelation  =  +self%correlationNormalization                                   &
            &                  *                    massRatio**self%correlationMassExponent     &
            &                  *         expansionFactorRatio**self%correlationRedshiftExponent
    end if
    return
  end function powerLawCorrelation
