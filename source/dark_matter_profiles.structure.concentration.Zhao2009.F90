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

  !!{
  An implementation of dark matter halo profile concentrations using the
  \cite{zhao_accurate_2009} algorithm.
  !!}

  use :: Cosmology_Functions                      , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                     , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistory  , darkMatterHaloMassAccretionHistoryClass
  use :: Dark_Matter_Profiles_DMO                 , only : darkMatterProfileDMONFW
  use :: Virial_Density_Contrast                  , only : virialDensityContrastBryanNorman1998

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationZhao2009">
   <description>
    A dark matter profile concentration class in which the concentration is computed using a fitting function from
    \cite{zhao_accurate_2009}:
    \begin{equation}
     c = 4 \left(1 + \left[ {t  \over 3.75 t_\mathrm{form}}\right]^{8.4}\right)^{1/8},
    \end{equation}
    where $t$ is the time for the halo and $t_\mathrm{form}$ is a formation time defined by \cite{zhao_accurate_2009} as the
    time at which the main branch progenitor of the halo had a mass equal to $0.04$ of the current halo mass. This formation
    time is computed directly from the merger tree branch associated with each halo. If the no branch exists or does not extend
    to the formation time then the formation time is computed by extrapolating the mass of the earliest resolved main branch
    progenitor to earlier times using the selected \refClass{darkMatterHaloMassAccretionHistoryClass}.
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationZhao2009
     !!{
     A dark matter halo profile concentration class implementing the algorithm of
     \cite{zhao_accurate_2009}.
     !!}
     private
     class(cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_                 => null()
     class(cosmologyParametersClass               ), pointer :: cosmologyParameters_                => null()
     class(darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     type (virialDensityContrastBryanNorman1998   ), pointer :: virialDensityContrastDefinition_    => null()
     type (darkMatterProfileDMONFW                ), pointer :: darkMatterProfileDMODefinition_     => null()
   contains
     final     ::                                   zhao2009Destructor
     procedure :: concentration                  => zhao2009Concentration
     procedure :: densityContrastDefinition      => zhao2009DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => zhao2009DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationZhao2009

  interface darkMatterProfileConcentrationZhao2009
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationZhao2009} dark matter halo profile
     concentration class.
     !!}
     module procedure zhao2009ConstructorParameters
     module procedure zhao2009ConstructorInternal
  end interface darkMatterProfileConcentrationZhao2009

contains

  function zhao2009ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationZhao2009} dark matter halo profile
    concentration class which takes an input parameter list.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileConcentrationZhao2009 )                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class(darkMatterHaloMassAccretionHistoryClass), pointer       :: darkMatterHaloMassAccretionHistory_

    !![
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"/>
    <objectBuilder class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"/>
    <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationZhao2009(cosmologyFunctions_,cosmologyParameters_,darkMatterHaloMassAccretionHistory_)
    !![
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="cosmologyParameters_"               />
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    !!]
    return
    !![
    <inputParametersValidate source="parameters"/>
    !!]
  end function zhao2009ConstructorParameters

  function zhao2009ConstructorInternal(cosmologyFunctions_,cosmologyParameters_,darkMatterHaloMassAccretionHistory_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileConcentrationZhao2009} dark matter halo profile
    concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    implicit none
    type (darkMatterProfileConcentrationZhao2009            )                         :: self
    class(cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    class(cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class(darkMatterHaloMassAccretionHistoryClass           ), intent(in   ), target  :: darkMatterHaloMassAccretionHistory_
    type (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_, *darkMatterHaloMassAccretionHistory_"/>
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
  end function zhao2009ConstructorInternal

  subroutine zhao2009Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationZhao2009} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationZhao2009), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%cosmologyParameters_"               />
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="self%virialDensityContrastDefinition_"   />
    <objectDestructor name="self%darkMatterProfileDMODefinition_"    />
    !!]
    return
  end subroutine zhao2009Destructor

  double precision function zhao2009Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the \cite{zhao_accurate_2009} algorithm.
    !!}
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time
    use :: Galacticus_Nodes                , only : nodeComponentBasic             , treeNode
    implicit none
    class           (darkMatterProfileConcentrationZhao2009), intent(inout), target  :: self
    type            (treeNode                              ), intent(inout), target  :: node
    class           (nodeComponentBasic                    )               , pointer :: basic
    double precision                                        , parameter              :: concentrationMinimum =4.00d0
    double precision                                        , parameter              :: formationMassFraction=0.04d0
    double precision                                                                 :: timeFormation               , timeNode
    !$GLC attributes unused :: self

    ! Get the basic component.
    basic => node%basic()
    ! Compute the concentration.
    timeNode     =basic%time()
    timeFormation=Dark_Matter_Halo_Formation_Time(node,formationMassFraction,self%darkMatterHaloMassAccretionHistory_)
    ! Compute the concentration from the formation time using the Zhao et al. (2009) fitting formula.
    zhao2009Concentration=concentrationMinimum*(1.0d0+(timeNode/3.75d0/timeFormation)**8.4d0)**0.125d0
   return
  end function zhao2009Concentration

  function zhao2009DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of
    concentration in the \cite{zhao_accurate_2009} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass            ), pointer       :: zhao2009DensityContrastDefinition
    class(darkMatterProfileConcentrationZhao2009), intent(inout) :: self

    zhao2009DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function zhao2009DensityContrastDefinition

  function zhao2009DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of
    concentration in the \cite{zhao_accurate_2009} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass             ), pointer       :: zhao2009DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationZhao2009), intent(inout) :: self

    zhao2009DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function zhao2009DarkMatterProfileDefinition

