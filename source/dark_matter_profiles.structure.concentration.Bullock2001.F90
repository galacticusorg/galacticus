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
  An implementation of dark matter halo profile concentrations using the \cite{navarro_structure_1996} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Virial_Density_Contrast   , only : virialDensityContrastClass   , virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationBullock2001">
   <description>Dark matter halo concentrations are computed using the algorithm of \cite{bullock_profiles_2001}.</description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationBullock2001
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{bullock_profiles_2001}.
     !!}
     private
     class           (cosmologyParametersClass                                      ), pointer :: cosmologyParameters_             => null()
     class           (cosmologyFunctionsClass                                       ), pointer :: cosmologyFunctions_              => null()
     class           (criticalOverdensityClass                                      ), pointer :: criticalOverdensity_             => null()
     class           (cosmologicalMassVarianceClass                                 ), pointer :: cosmologicalMassVariance_        => null()
     class           (virialDensityContrastClass                                    ), pointer :: virialDensityContrast_           => null()
     type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMONFW                                       ), pointer :: darkMatterProfileDMODefinition_  => null()
     double precision                                                                          :: F                                         , K
   contains
     final     ::                                   bullock2001Destructor
     procedure :: concentration                  => bullock2001Concentration
     procedure :: densityContrastDefinition      => bullock2001DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => bullock2001DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationBullock2001

  interface darkMatterProfileConcentrationBullock2001
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationBullock2001} dark matter halo profile concentration class.
     !!}
     module procedure bullock2001ConstructorParameters
     module procedure bullock2001ConstructorInternal
  end interface darkMatterProfileConcentrationBullock2001

contains

  function bullock2001ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily bullock2001} dark matter halo profile
    concentration class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileConcentrationBullock2001)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass                 ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass            ), pointer       :: cosmologicalMassVariance_
    class           (virialDensityContrastClass               ), pointer       :: virialDensityContrast_
    double precision                                                           :: F                        , K

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>F</name>
      <source>parameters</source>
      <defaultValue>0.01d0</defaultValue>
      <defaultSource>\cite{bullock_profiles_2001}</defaultSource>
      <description>The parameter $F$ appearing in the halo concentration algorithm of \cite{bullock_profiles_2001}.</description>
    </inputParameter>
    <inputParameter>
      <name>K</name>
      <source>parameters</source>
      <defaultValue>4.0d0</defaultValue>
      <defaultSource>\cite{bullock_profiles_2001}</defaultSource>
      <description>The parameter $K$ appearing in the halo concentration algorithm of \cite{bullock_profiles_2001}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="virialDensityContrast"    name="virialDensityContrast_"    source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationBullock2001(F,K,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="virialDensityContrast_"   />
    !!]
    return
  end function bullock2001ConstructorParameters

  function bullock2001ConstructorInternal(F,K,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationBullock2001} dark matter halo profile
    concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    implicit none
    type            (darkMatterProfileConcentrationBullock2001         )                         :: self
    double precision                                                    , intent(in   )          :: F                              , K
    class           (cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    class           (criticalOverdensityClass                          ), intent(in   ), target  :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                     ), intent(in   ), target  :: cosmologicalMassVariance_
    class           (virialDensityContrastClass                        ), intent(in   ), target  :: virialDensityContrast_
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="F,K,*cosmologyParameters_,*cosmologyFunctions_,*criticalOverdensity_,*cosmologicalMassVariance_,*virialDensityContrast_"/>
    !!]

    allocate(self%darkMatterProfileDMODefinition_ )
    allocate(     darkMatterHaloScaleDefinition_  )
    allocate(self%virialDensityContrastDefinition_)
    !![
    <referenceConstruct owner="self" object="virialDensityContrastDefinition_">
     <constructor>
      virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                            &amp;
       &amp;                                                         tableStore                          =.true.                               , &amp;
       &amp;                                                         cosmologyFunctions_                 =self%cosmologyFunctions_               &amp;
       &amp;                                                        )
     </constructor>
    </referenceConstruct>
    <referenceConstruct              object="darkMatterHaloScaleDefinition_"  >
     <constructor>
      darkMatterHaloScaleVirialDensityContrastDefinition            (                                                                            &amp;
       &amp;                                                         cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                                         cosmologyFunctions_                 =self%cosmologyFunctions_             , &amp;
       &amp;                                                         virialDensityContrast_              =self%virialDensityContrastDefinition_  &amp;
       &amp;                                                        )
     </constructor>
    </referenceConstruct>
    <referenceConstruct owner="self" object="darkMatterProfileDMODefinition_" >
     <constructor>
      darkMatterProfileDMONFW                                       (                                                                            &amp;
       &amp;                                                         velocityDispersionUseSeriesExpansion=.true.                               , &amp;
       &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScaleDefinition_         &amp;
       &amp;                                                        )
     </constructor>
    </referenceConstruct>
    <objectDestructor                name  ="darkMatterHaloScaleDefinition_" />
    !!]
    return
  end function bullock2001ConstructorInternal

  subroutine bullock2001Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationBullock2001} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationBullock2001), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%criticalOverdensity_"            />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine bullock2001Destructor

  double precision function bullock2001Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the \cite{bullock_profiles_2001} algorithm.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Virial_Density_Contrast             , only : virialDensityContrastClass
    implicit none
    class           (darkMatterProfileConcentrationBullock2001), intent(inout), target  :: self
    type            (treeNode                                 ), intent(inout), target  :: node
    class           (virialDensityContrastClass               ), pointer                :: virialDensityContrast_
    class           (nodeComponentBasic                       )               , pointer :: basic
    double precision                                                                    :: massHalo              , massHaloFormation       , &
         &                                                                                 sigmaFormation        , timeFormation           , &
         &                                                                                 expansionFactor       , expansionFactorFormation, &
         &                                                                                 densityContrast

    ! Compute the characteristic mass at formation time.
    basic                  =>  node                               %basic                    ()
    virialDensityContrast_ =>  self                               %densityContrastDefinition()
    densityContrast        =   virialDensityContrast_             %densityContrast          (                                                     &
         &                                                                                                          basic%mass            ()    , &
         &                                                                                                          basic%timeLastIsolated()      &
         &                                                                                  )
    massHalo               =   Dark_Matter_Profile_Mass_Definition                          (                                                     &
         &                                                                                   node                                               , &
         &                                                                                   densityContrast                                    , &
         &                                                                                   cosmologyParameters_  =self %cosmologyParameters_  , &
         &                                                                                   cosmologyFunctions_   =self %cosmologyFunctions_   , &
         &                                                                                   virialDensityContrast_=self %virialDensityContrast_  &
         &                                                                                  )
    massHaloFormation      =  +self%F*massHalo
    ! Compute the corresponding rms fluctuation in the density field (i.e. σ(M)).
    sigmaFormation          =+self%cosmologicalMassVariance_%rootVariance   (                    massHaloFormation,     self%cosmologyFunctions_%cosmicTime(1.0d0))
    ! Get the time at which this equals the critical overdensity for collapse.
    timeFormation           =+self%criticalOverdensity_     %timeOfCollapse (criticalOverdensity=sigmaFormation   ,mass=massHalo                                  )
    ! Get the corresponding expansion factors.
    expansionFactorFormation=+self%cosmologyFunctions_      %expansionFactor(                    timeFormation                                                    )
    expansionFactor         =+self%cosmologyFunctions_      %expansionFactor(                    basic%time()                                                     )
    ! Compute the concentration.
    bullock2001Concentration=+self%K                   &
         &                   *expansionFactor          &
         &                   /expansionFactorFormation
    return
  end function bullock2001Concentration

  function bullock2001DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{bullock_profiles_2001} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass               ), pointer       :: bullock2001DensityContrastDefinition
    class(darkMatterProfileConcentrationBullock2001), intent(inout) :: self

    bullock2001DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function bullock2001DensityContrastDefinition

  function bullock2001DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{bullock_profiles_2001} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass                ), pointer       :: bullock2001DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationBullock2001), intent(inout) :: self

    bullock2001DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function bullock2001DarkMatterProfileDefinition

