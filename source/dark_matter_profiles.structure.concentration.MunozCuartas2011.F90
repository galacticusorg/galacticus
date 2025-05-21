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
  An implementation of dark matter halo profile concentrations using the \cite{munoz-cuartas_redshift_2011} algorithm.
  !!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMONFW
  use :: Virial_Density_Contrast , only : virialDensityContrastBryanNorman1998

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationMunozCuartas2011">
   <description>
    A dark matter profile concentration class in which the concentration is computed using a fitting function from
    \cite{munoz-cuartas_redshift_2011}:
    \begin{equation}
    \log_{10} c = a \log_{10} \left( {M_\mathrm{halo} \over h^{-1}M_\odot} \right) + b.
    \end{equation}
    The parameters are a function of redshift, $z$, given by
    \begin{eqnarray}
    a &amp;=&amp; wz-m, \\
    b &amp;=&amp; {\alpha \over (z+\gamma)} + {\beta \over (z+\gamma)^2},
    \end{eqnarray}
    where $w=0.029$, $m=0.097$, $\alpha=-110.001$, $\beta=2469.720$, $\gamma=16.885$.
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationMunozCuartas2011
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{munoz-cuartas_redshift_2011}.
     !!}
     private
     class(cosmologyFunctionsClass             ), pointer :: cosmologyFunctions_              => null()
     class(cosmologyParametersClass            ), pointer :: cosmologyParameters_             => null()
     type (virialDensityContrastBryanNorman1998), pointer :: virialDensityContrastDefinition_ => null()
     type (darkMatterProfileDMONFW             ), pointer :: darkMatterProfileDMODefinition_  => null()
   contains
     final     ::                                   munozCuartas2011Destructor
     procedure :: concentration                  => munozCuartas2011Concentration
     procedure :: densityContrastDefinition      => munozCuartas2011DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => munozCuartas2011DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationMunozCuartas2011

  interface darkMatterProfileConcentrationMunozCuartas2011
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationMunozCuartas2011} dark matter halo profile concentration class.
     !!}
     module procedure munozCuartas2011ConstructorParameters
     module procedure munozCuartas2011ConstructorInternal
  end interface darkMatterProfileConcentrationMunozCuartas2011

contains

  function munozCuartas2011ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationMunozCuartas2011} dark matter halo profile concentration class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileConcentrationMunozCuartas2011)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                       ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass                      ), pointer       :: cosmologyParameters_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationMunozCuartas2011(cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function munozCuartas2011ConstructorParameters

  function munozCuartas2011ConstructorInternal(cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileConcentrationMunozCuartas2011} dark matter halo profile concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    implicit none
    type (darkMatterProfileConcentrationMunozCuartas2011    )                         :: self
    class(cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class(cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    type (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    allocate(self%virialDensityContrastDefinition_)
    allocate(self%darkMatterProfileDMODefinition_ )
    allocate(     darkMatterHaloScaleDefinition_  )
    !![
    <referenceConstruct owner="self" object="virialDensityContrastDefinition_">
     <constructor>
      virialDensityContrastBryanNorman1998              (                                                                            &amp;
       &amp;                                             allowUnsupportedCosmology           =     .true.                          , &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_               &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct              object="darkMatterHaloScaleDefinition_"  >
     <constructor>
      darkMatterHaloScaleVirialDensityContrastDefinition(                                                                            &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_             , &amp;
       &amp;                                             virialDensityContrast_              =self%virialDensityContrastDefinition_  &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct owner="self" object="darkMatterProfileDMODefinition_" >
     <constructor>
      darkMatterProfileDMONFW                           (                                                                            &amp;
       &amp;                                             velocityDispersionUseSeriesExpansion=.true.                               , &amp;
       &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleDefinition_         &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <objectDestructor                name  ="darkMatterHaloScaleDefinition_" />
    !!]
    return
  end function munozCuartas2011ConstructorInternal

  subroutine munozCuartas2011Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationMunozCuartas2011} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationMunozCuartas2011), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine munozCuartas2011Destructor

  double precision function munozCuartas2011Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite{munoz-cuartas_redshift_2011} algorithm.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Galacticus_Nodes    , only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileConcentrationMunozCuartas2011), intent(inout), target  :: self
    type            (treeNode                                      ), intent(inout), target  :: node
    class           (nodeComponentBasic                            )               , pointer :: basic
    double precision                                                , parameter              :: alpha                   =-110.001d0, beta               =2469.720d0, &
         &                                                                                      gamma                   =  16.885d0, m                  =0.097d0   , &
         &                                                                                      w                       =   0.029d0
    double precision                                                                         :: a                                  , b                             , &
         &                                                                                      concentrationLogarithmic           , haloMassLogarithmic           , &
         &                                                                                      redshift
    !$GLC attributes unused :: self

    ! Get required objects.
    basic => node%basic()
    ! Compute the concentration.
    redshift                     =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(basic%time()))
    a                            =w*redshift-m
    b                            =alpha/(redshift+gamma)+beta/(redshift+gamma)**2
    haloMassLogarithmic          =log10(basic%mass()*self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH))
    concentrationLogarithmic     =a*haloMassLogarithmic+b
    munozCuartas2011Concentration=10.0d0**concentrationLogarithmic
    return
  end function munozCuartas2011Concentration

  function munozCuartas2011DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{munoz-cuartas_redshift_2011} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass                    ), pointer       :: munozCuartas2011DensityContrastDefinition
    class(darkMatterProfileConcentrationMunozCuartas2011), intent(inout) :: self

    munozCuartas2011DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function munozCuartas2011DensityContrastDefinition

  function munozCuartas2011DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{munoz-cuartas_redshift_2011} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass                     ), pointer       :: munozCuartas2011DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationMunozCuartas2011), intent(inout) :: self

    munozCuartas2011DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function munozCuartas2011DarkMatterProfileDefinition
