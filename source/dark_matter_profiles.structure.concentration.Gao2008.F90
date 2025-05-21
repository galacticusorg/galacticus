!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  An implementation of dark matter halo profile concentrations using the \cite{gao_redshift_2008} algorithm.
  !!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMONFW
  use :: Virial_Density_Contrast , only : virialDensityContrastFixed

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationGao2008">
   <description>
    A dark matter profile concentration class in which the concentration is computed using a fitting function from
    \cite{gao_redshift_2008}:
    \begin{equation}
    \log_{10} c = A \log_{10} M_\mathrm{halo} + B.
    \end{equation}
    The parameters are a function of expansion factor, $a$. We use the following fits to the \cite{gao_redshift_2008} results:
    \begin{eqnarray}
    A &amp;=&amp; -0.140 \exp\left[-\left(\left\{\log_{10}a+0.05\right\}/0.35\right)^2\right], \\
    B &amp;=&amp;  2.646 \exp\left[-\left(\log_{10}a/0.50\right)^2\right].
    \end{eqnarray}
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
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{gao_redshift_2008}.
     !!}
     private
     class           (cosmologyParametersClass  ), pointer :: cosmologyParameters_             => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_              => null()
     type            (virialDensityContrastFixed), pointer :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMONFW   ), pointer :: darkMatterProfileDMODefinition_  => null()
     double precision                                      :: scatter
   contains
     final     ::                                   gao2008Destructor
     procedure :: concentration                  => gao2008Concentration
     procedure :: concentrationMean              => gao2008ConcentrationMean
     procedure :: densityContrastDefinition      => gao2008DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => gao2008DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationGao2008

  interface darkMatterProfileConcentrationGao2008
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationGao2008} dark matter halo profile concentration class.
     !!}
     module procedure gao2008ConstructorParameters
     module procedure gao2008ConstructorInternal
  end interface darkMatterProfileConcentrationGao2008

contains

  function gao2008ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationGao2008} dark matter halo profile concentration class which takes a parameter
    list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileConcentrationGao2008)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    double precision                                                       :: scatter

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <inputParameter>
      <name>scatter</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The scatter (in dex) to assume in the halo concentration distribution at fixed mass.</description>
    </inputParameter>
    !!]
    self=darkMatterProfileConcentrationGao2008(scatter,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function gao2008ConstructorParameters

  function gao2008ConstructorInternal(scatter,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileConcentrationGao2008} dark matter halo profile concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    implicit none
    type            (darkMatterProfileConcentrationGao2008             )                         :: self
    class           (cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    double precision                                                    , intent(in   )          :: scatter
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="scatter, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    allocate(self%darkMatterProfileDMODefinition_ )
    allocate(     darkMatterHaloScaleDefinition_  )
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
    <referenceConstruct owner="self" isResult="yes" object="darkMatterProfileDMODefinition_" >
     <constructor>
      darkMatterProfileDMONFW                           (                                                                            &amp;
       &amp;                                             velocityDispersionUseSeriesExpansion=.true.                               , &amp;
       &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleDefinition_         &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <objectDestructor                               name  ="darkMatterHaloScaleDefinition_" />
    !!]
    return
  end function gao2008ConstructorInternal

  subroutine gao2008Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationGao2008} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationGao2008), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine gao2008Destructor

  double precision function gao2008Concentration(self,node) result(concentration)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite{gao_redshift_2008} algorithm.
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
    !!{
    Return the mean concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite{gao_redshift_2008} algorithm.
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
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{gao_redshift_2008} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass           ), pointer       :: gao2008DensityContrastDefinition
    class(darkMatterProfileConcentrationGao2008), intent(inout) :: self

    gao2008DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function gao2008DensityContrastDefinition

  function gao2008DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{gao_redshift_2008} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass            ), pointer       :: gao2008DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationGao2008), intent(inout) :: self

    gao2008DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function gao2008DarkMatterProfileDefinition

