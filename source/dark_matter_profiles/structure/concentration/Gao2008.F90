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
  An implementation of dark matter halo profile concentrations using the :cite:t:`gao_redshift_2008` algorithm.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOEinasto
  use :: Virial_Density_Contrast   , only : virialDensityContrastFixed

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationGao2008" docformat="rst">
   <description>
   A dark matter profile concentration class in which the concentration is computed using a fitting function from :cite:t:`gao_redshift_2008`:

   .. math::

      \log_{10} c = A \log_{10} M_\mathrm{halo} + B.

   The parameters are a function of expansion factor, :math:`a`. We use the following fits to the :cite:t:`gao_redshift_2008` results:

   .. math::

      A &amp; = -0.140 \exp\left[-\left(\left\{\log_{10}a+0.05\right\}/0.35\right)^2\right], \\
      B &amp; = 2.646 \exp\left[-\left(\log_{10}a/0.50\right)^2\right].

   For halo masses of :math:`10^{11}` to :math:`10^{15}h^{-1}\mathrm{M}_\odot`, these fits reproduce the concentrations given in
   Table 1 of :cite:t:`gao_redshift_2008` to about 3% at :math:`z \le 0.5`, but are about 10--19% too high at :math:`z=1`--2,
   and about 10--15% too low at :math:`z=3`.

   Following :cite:t:`gao_redshift_2008`, whose concentrations are :math:`c=r_{200}/r_{-2}` from fits of Einasto profiles with
   the shape parameter fixed by :math:`\alpha(\nu)` (see :galacticus-class:`darkMatterProfileShapeGao2008`), the density
   profile which defines these concentrations is Einasto, with that :math:`\alpha(\nu)`, and the density contrast is 200
   times the critical density.
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationGao2008
     !!{RST
     A dark matter halo profile concentration class implementing the algorithm of :cite:t:`gao_redshift_2008`.
     !!}
     private
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_             => null()
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_              => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_             => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_        => null()
     type            (virialDensityContrastFixed   ), pointer :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMOEinasto  ), pointer :: darkMatterProfileDMODefinition_  => null()
     double precision                                      :: scatter
   contains
     final     ::                                   gao2008Destructor
     procedure :: concentration                  => gao2008Concentration
     procedure :: concentrationMean              => gao2008ConcentrationMean
     procedure :: densityContrastDefinition      => gao2008DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => gao2008DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationGao2008

  interface darkMatterProfileConcentrationGao2008
     !!{RST
     Constructors for the :galacticus-class:`darkMatterProfileConcentrationGao2008` dark matter halo profile concentration class.
     !!}
     module procedure gao2008ConstructorParameters
     module procedure gao2008ConstructorInternal
  end interface darkMatterProfileConcentrationGao2008

contains

  function gao2008ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterProfileConcentrationGao2008` dark matter halo profile concentration class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileConcentrationGao2008)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    class           (criticalOverdensityClass             ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass        ), pointer       :: cosmologicalMassVariance_
    double precision                                                       :: scatter

    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <inputParameter docformat="rst">
      <name>scatter</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The scatter (in dex) to assume in the halo concentration distribution at fixed mass.
      </description>
      <minimum>0.0</minimum>
    </inputParameter>
    !!]
    self=darkMatterProfileConcentrationGao2008(scatter,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function gao2008ConstructorParameters

  function gao2008ConstructorInternal(scatter,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterProfileConcentrationGao2008` dark matter halo profile concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Dark_Matter_Profiles_Shape, only : darkMatterProfileShapeGao2008
    use :: Virial_Density_Contrast   , only : fixedDensityTypeCritical
    implicit none
    type            (darkMatterProfileConcentrationGao2008             )                         :: self
    class           (cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    class           (criticalOverdensityClass                          ), intent(in   ), target  :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                     ), intent(in   ), target  :: cosmologicalMassVariance_
    double precision                                                    , intent(in   )          :: scatter
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    type            (darkMatterProfileShapeGao2008                     )               , pointer :: darkMatterProfileShapeDefinition_
    !![
    <constructorAssign variables="scatter, *cosmologyParameters_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    allocate(self%darkMatterProfileDMODefinition_ )
    allocate(     darkMatterHaloScaleDefinition_  )
    allocate(     darkMatterProfileShapeDefinition_)
    allocate(self%virialDensityContrastDefinition_)
    !![
    <referenceConstruct owner="self" isResult="yes" object="virialDensityContrastDefinition_">
     <constructor>
      virialDensityContrastFixed                        (                                                                            &amp;
       &amp;                                             densityContrastValue                =200.0d0                              , &amp;
       &amp;                                             densityType                         =fixedDensityTypeCritical             , &amp;
       &amp;                                             turnAroundOverVirialRadius          =2.0d0                                , &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_               &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct                             object="darkMatterHaloScaleDefinition_"  >
     <constructor>
      darkMatterHaloScaleVirialDensityContrastDefinition(                                                                            &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_             , &amp;
       &amp;                                             virialDensityContrast_              =self%virialDensityContrastDefinition_  &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct                             object="darkMatterProfileShapeDefinition_">
     <constructor>
      darkMatterProfileShapeGao2008                     (                                                                            &amp;
       &amp;                                             criticalOverdensity_                =self%criticalOverdensity_            , &amp;
       &amp;                                             cosmologicalMassVariance_           =self%cosmologicalMassVariance_         &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct owner="self" isResult="yes" object="darkMatterProfileDMODefinition_" >
     <constructor>
      darkMatterProfileDMOEinasto                       (                                                                            &amp;
       &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleDefinition_       , &amp;
       &amp;                                             darkMatterProfileShape_             =darkMatterProfileShapeDefinition_      &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <objectDestructor                               name  ="darkMatterHaloScaleDefinition_"   />
    <objectDestructor                               name  ="darkMatterProfileShapeDefinition_"/>
    !!]
    return
  end function gao2008ConstructorInternal

  subroutine gao2008Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`darkMatterProfileConcentrationGao2008` dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationGao2008), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%criticalOverdensity_"            />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine gao2008Destructor

  double precision function gao2008Concentration(self,node) result(concentration)
    !!{RST
    Return the concentration of the dark matter halo profile of ``node`` using the :cite:t:`gao_redshift_2008` algorithm.
    !!}
    implicit none
    class(darkMatterProfileConcentrationGao2008), intent(inout), target :: self
    type (treeNode                             ), intent(inout), target :: node

    ! Get the mean concentration.
    concentration=self%concentrationMean(node)
    ! Add scatter if necessary.
    if (self%scatter > 0.0d0)                                                                    &
         &  concentration=+concentration                                                         &
         &                *10.0d0**(                                                             &
         &                          +self%scatter                                                &
         &                          *node%hostTree%randomNumberGenerator_%standardNormalSample() &
         &                         )
    return
  end function gao2008Concentration

  double precision function gao2008ConcentrationMean(self,node) result(concentration)
    !!{RST
    Return the mean concentration of the dark matter halo profile of ``node`` using the :cite:t:`gao_redshift_2008` algorithm.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileConcentrationGao2008), intent(inout)          :: self
    type            (treeNode                             ), intent(inout), target  :: node
    class           (nodeComponentBasic                   )               , pointer :: basic
    double precision                                       , parameter              :: littleHubbleConstantGao2008=0.73d0
    double precision                                                                :: logarithmExpansionFactor          , logarithmHaloMass, &
         &                                                                             parameterA                        , parameterB

    ! Get the basic component.
    basic                   => node%basic()
    ! Compute the concentration.
    logarithmHaloMass       =log10(littleHubbleConstantGao2008             *basic%mass())
    logarithmExpansionFactor=log10(self%cosmologyFunctions_%expansionFactor(basic%time()))
    parameterA              =- 0.140d0*exp(-((logarithmExpansionFactor+0.05d0)/0.35d0)**2)
    parameterB              =+ 2.646d0*exp(-((logarithmExpansionFactor+0.00d0)/0.50d0)**2)
    concentration          =+10.0d0**(parameterA*logarithmHaloMass+parameterB)
    return
  end function gao2008ConcentrationMean

  function gao2008DensityContrastDefinition(self)
    !!{RST
    Return a virial density contrast object defining that used in the definition of concentration in the :cite:t:`gao_redshift_2008` algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass           ), pointer       :: gao2008DensityContrastDefinition
    class(darkMatterProfileConcentrationGao2008), intent(inout) :: self

    gao2008DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function gao2008DensityContrastDefinition

  function gao2008DarkMatterProfileDefinition(self)
    !!{RST
    Return a dark matter density profile object defining that used in the definition of concentration in the :cite:t:`gao_redshift_2008` algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass            ), pointer       :: gao2008DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationGao2008), intent(inout) :: self

    gao2008DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function gao2008DarkMatterProfileDefinition

