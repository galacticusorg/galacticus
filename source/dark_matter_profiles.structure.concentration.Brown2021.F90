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
  An implementation of dark matter halo profile concentrations using the
  \cite{brown_towards_2022} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOEinasto
  use :: Virial_Density_Contrast   , only : virialDensityContrastFixed

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationBrown2021">
   <description>Dark matter halo concentrations are computed using the algorithm of \cite{brown_towards_2022}.</description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationBrown2021
     !!{
     A dark matter halo profile concentration class implementing the algorithm of
     \cite[][eqn. 20]{brown_towards_2022}. Specifically the concentration is given by
     \begin{equation}
       c_\mathrm{200c} = 4.39 \nu_\mathrm{c}^{-0.87},
     \end{equation}
     where $\nu_\mathrm{c} = \delta_\mathrm{c}/\sigma_\mathrm{c}(M)$ is the peak height.
     
     This implementation accepts any \refClass{cosmologicalMassVarianceClass} object for use in computing for computing
     $\sigma_\mathrm{c}(M)$. \emph{However}, \cite{brown_towards_2022} recommend using $\sigma_\mathrm{c}(M)$ computed using a
     generalized top-hat window function (\refClass{powerSpectrumWindowFunctionTopHatGeneralized}) with $\mu_\mathrm{g}=0.2138$.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer     :: cosmologyFunctions_              => null()
     class           (cosmologyParametersClass     ), pointer     :: cosmologyParameters_             => null()
     class           (criticalOverdensityClass     ), pointer     :: criticalOverdensity_             => null()
     class           (cosmologicalMassVarianceClass), pointer     :: cosmologicalMassVariance_        => null()
     type            (virialDensityContrastFixed   ), pointer     :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMOEinasto  ), pointer     :: darkMatterProfileDMODefinition_  => null()
   contains
     final     ::                                   brown2021Destructor
     procedure :: concentration                  => brown2021Concentration
     procedure :: densityContrastDefinition      => brown2021DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => brown2021DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationBrown2021

  interface darkMatterProfileConcentrationBrown2021
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationBrown2021} dark matter halo profile concentration class.
     !!}
     module procedure brown2021ConstructorParameters
     module procedure brown2021ConstructorInternal
  end interface darkMatterProfileConcentrationBrown2021

contains

  function brown2021ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily brown2021} dark matter halo
    profile concentration class.
    !!}
    implicit none
    type (darkMatterProfileConcentrationBrown2021)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class(criticalOverdensityClass               ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass          ), pointer       :: cosmologicalMassVariance_

    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
     self=darkMatterProfileConcentrationBrown2021(cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function brown2021ConstructorParameters

  function brown2021ConstructorInternal(cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationBrown2021} dark matter halo profile
    concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    implicit none
    type (darkMatterProfileConcentrationBrown2021           )                         :: self
    class(cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    class(cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class(criticalOverdensityClass                          ), intent(in   ), target  :: criticalOverdensity_
    class(cosmologicalMassVarianceClass                     ), intent(in   ), target  :: cosmologicalMassVariance_
    type (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    allocate(     darkMatterHaloScaleDefinition_  )
    allocate(self%virialDensityContrastDefinition_)
    allocate(self%darkMatterProfileDMODefinition_ )
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
      darkMatterProfileDMOEinasto                       (                                                                            &amp;
       &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleDefinition_         &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <objectDestructor                               name  ="darkMatterHaloScaleDefinition_" />
    !!]
    return
  end function brown2021ConstructorInternal

  subroutine brown2021Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationBrown2021} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationBrown2021), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%criticalOverdensity_"            />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine brown2021Destructor

  double precision function brown2021Concentration(self,node)
    !!{
    Return the mean concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the \cite{brown_towards_2022} algorithm.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic
    implicit none
    class           (darkMatterProfileConcentrationBrown2021), intent(inout), target  :: self
    type            (treeNode                               ), intent(inout), target  :: node
    class           (nodeComponentBasic                     )               , pointer :: basic
    double precision                                                                  :: peakHeight

    basic      =>  node%basic()
    peakHeight =  +self%criticalOverdensity_     %value                          (time=basic%time(),mass=basic%mass(),node=node) &
         &        /self%cosmologicalMassVariance_%rootVariance                   (time=basic%time(),mass=basic%mass()          )
    ! Evaluate the concentration using equation (20) of Brown et al. (2021).
    brown2021Concentration=+4.39d0             &
         &                 /peakHeight**0.87d0
    return
  end function brown2021Concentration

  function brown2021DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{brown_towards_2022} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass             ), pointer       :: brown2021DensityContrastDefinition
    class(darkMatterProfileConcentrationBrown2021), intent(inout) :: self

    brown2021DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function brown2021DensityContrastDefinition

  function brown2021DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{brown_towards_2022} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass              ), pointer       :: brown2021DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationBrown2021), intent(inout) :: self

    brown2021DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function brown2021DarkMatterProfileDefinition
