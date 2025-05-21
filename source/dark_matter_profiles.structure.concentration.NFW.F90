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
  An implementation of dark matter halo profile concentrations using the \cite{navarro_structure_1996} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMONFW
  use :: Virial_Density_Contrast   , only : virialDensityContrast        , virialDensityContrastClass, virialDensityContrastFixed
  use :: Root_Finder               , only : rootFinder

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationNFW1996">
   <description>
    A dark matter profile concentration class in which the concentration is computed using the algorithm from
    \cite{navarro_structure_1996}. In this algorithm, for a given halo of mass $M$ at time $t_0$, a formation time is defined as
    the epoch at which there is a 50\% probability (according to extended Press-Schechter theory) for a progenitor halo to have
    a mass greater than $fM$, where $f=${\normalfont \ttfamily [f]} is a parameter of the algorithm. This implies formation
    when the critical overdensity for collapse is
    \begin{equation}
     \delta_\mathrm{crit}(t_\mathrm{form}) = \left[ 2 \nu_{1/2}^2 \left\{\sigma(fM)^22-\sigma(M)^2\right\}
     \right]^{1/2}+\delta_\mathrm{crit}(t_0),
    \end{equation}
    where $\nu_{1/2} = [\hbox{erfc}^{-1}(1/2)]^{1/2}$. \cite{navarro_structure_1996} then assume an overdensity at collapse of
    \begin{equation}
     \Delta(t_\mathrm{form}) = C  \left[ {a(t_0) \over a(t_\mathrm{form})} \right]^3
    \end{equation}
    where $C=${\normalfont \ttfamily [C]} is a parameter of the algorithm. The concentration is then determined by solving
    \begin{equation}
     {\Delta(t_\mathrm{form}) \over \Delta_\mathrm{virial}(t_0)} = {c^3 \over 3 [\ln(1+c)-c/(1+c)]}.
    \end{equation}
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationNFW1996
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{navarro_structure_1996}.
     !!}
     private
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_             => null()
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_              => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_             => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_        => null()
     class           (virialDensityContrastClass   ), pointer :: virialDensityContrast_           => null()
     type            (virialDensityContrastFixed   ), pointer :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMONFW      ), pointer :: darkMatterProfileDMODefinition_  => null()
     type            (rootFinder                   )          :: finder
     double precision                                         :: f                                         , C
   contains
     final     ::                                   nfw1996Destructor
     procedure :: concentration                  => nfw1996Concentration
     procedure :: densityContrastDefinition      => nfw1996DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => nfw1996DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationNFW1996

  interface darkMatterProfileConcentrationNFW1996
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationNFW1996} dark matter halo profile concentration class.
     !!}
     module procedure nfw1996ConstructorParameters
     module procedure nfw1996ConstructorInternal
  end interface darkMatterProfileConcentrationNFW1996

  ! Target value used in concentration root finder.
  double precision :: rootTarget
  !$omp threadprivate(rootTarget)

contains

  function nfw1996ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily nfw1996} dark matter halo profile
    concentration class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileConcentrationNFW1996)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass             ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass        ), pointer       :: cosmologicalMassVariance_
    class           (virialDensityContrastClass           ), pointer       :: virialDensityContrast_
    double precision                                                       :: f                        , C

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>f</name>
      <source>parameters</source>
      <variable>f</variable>
      <defaultValue>0.01d0</defaultValue>
      <defaultSource>\cite{navarro_structure_1996}</defaultSource>
      <description>The parameter $f$ appearing in the halo concentration algorithm of \cite{navarro_structure_1996}.</description>
    </inputParameter>
    <inputParameter>
      <name>C</name>
      <source>parameters</source>
      <variable>C</variable>
      <defaultValue>2000.0d0</defaultValue>
      <defaultSource>\cite{navarro_structure_1996}</defaultSource>
      <description>The parameter $C$ appearing in the halo concentration algorithm of \cite{navarro_structure_1996}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="virialDensityContrast"    name="virialDensityContrast_"    source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationNFW1996(f,C,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="virialDensityContrast_"   />
    !!]
    return
  end function nfw1996ConstructorParameters

  function nfw1996ConstructorInternal(f,C,cosmologyParameters_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,virialDensityContrast_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationNFW1996} dark matter halo profile
    concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    use :: Root_Finder            , only : rangeExpandMultiplicative
    implicit none
    type            (darkMatterProfileConcentrationNFW1996             )                         :: self
    double precision                                                    , intent(in   )          :: f                                   , C
    class           (cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    class           (criticalOverdensityClass                          ), intent(in   ), target  :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                     ), intent(in   ), target  :: cosmologicalMassVariance_
    class           (virialDensityContrastClass                        ), intent(in   ), target  :: virialDensityContrast_
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    double precision                                                    , parameter              :: toleranceAbsolute             =0.0d0, toleranceRelative=1.0d-6
    !![
    <constructorAssign variables="f, C, *cosmologyParameters_, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *virialDensityContrast_"/>
    !!]
    
    self%finder=rootFinder(                                               &
         &                 rootFunction       =nfw1996ConcentrationRoot , &
         &                 rangeExpandDownward=0.5d0                    , &
         &                 rangeExpandUpward  =2.0d0                    , &
         &                 rangeExpandType    =rangeExpandMultiplicative, &
         &                 toleranceAbsolute  =toleranceAbsolute        , &
         &                 toleranceRelative  =toleranceRelative          &
         &                )
    allocate(self%virialDensityContrastDefinition_)
    allocate(     darkMatterHaloScaleDefinition_  )
    allocate(self%darkMatterProfileDMODefinition_ )
    !![
    <referenceConstruct owner="self" object="virialDensityContrastDefinition_">
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
  end function nfw1996ConstructorInternal

  subroutine nfw1996Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationNFW1996} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationNFW1996), intent(inout) :: self

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
  end subroutine nfw1996Destructor

  double precision function nfw1996Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the \cite{navarro_structure_1996} algorithm.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileConcentrationNFW1996), intent(inout), target  :: self
    type            (treeNode                             ), intent(inout), target  :: node
    double precision                                       , parameter              :: fitParameterNuHalf         =0.47693628d0
    class           (nodeComponentBasic                   )               , pointer :: basic
    double precision                                                                :: collapseCriticalOverdensity             , collapseExpansionFactor       , &
         &                                                                             collapseMass                            , collapseOverdensity           , &
         &                                                                             collapseTime                            , expansionFactor               , &
         &                                                                             nodeMass                                , nodeTime                      , &
         &                                                                             nowTime

    ! Get the basic component.
    basic                      =>  node                     %basic          (        )
    ! Get the properties of the node.
    nodeMass                   =   basic                    %mass           (        )
    nodeTime                   =   basic                    %time           (        )
    expansionFactor            =   self %cosmologyFunctions_%expansionFactor(nodeTime)
    ! Compute the mass of a progenitor as defined by NFW.
    collapseMass               =  +self%F   &
         &                        *nodeMass
    ! Find the time of collapse for this progenitor. The critical overdensity at collapse is scaled by a factor σ(M,t₀)/σ(M,t)
    ! because the definition of critical overdensity in Galacticus does not include the 1/D(t) factor that is included in the
    ! definition used by Navarro et al. (1996).
    nowTime                    =  +self%cosmologyFunctions_%cosmicTime      (   1.0d0)
    collapseCriticalOverdensity=+(                                                                                         &
         &                        +sqrt(                                                                                   &
         &                              +2.0d0                                                                             &
         &                              *fitParameterNuHalf**2                                                             &
         &                              *(                                                                                 &
         &                                +self%cosmologicalMassVariance_%rootVariance(time=nodeTime,mass=collapseMass)**2 &
         &                                -self%cosmologicalMassVariance_%rootVariance(time=nodeTime,mass=    nodeMass)**2 &
         &                               )                                                                                 &
         &                             )                                                                                   &
         &                        +        self%criticalOverdensity_     %value       (time=nodeTime,mass=    nodeMass)    &
         &                       )                                                                                         &
         &                      *          self%cosmologicalMassVariance_%rootVariance(time= nowTime,mass=    nodeMass)    &
         &                      /          self%cosmologicalMassVariance_%rootVariance(time=nodeTime,mass=    nodeMass)
    collapseTime               =self%criticalOverdensity_%timeOfCollapse (collapseCriticalOverdensity,mass=nodeMass)
    collapseExpansionFactor    =self%cosmologyFunctions_ %expansionFactor(collapseTime                             )
    ! Compute the overdensity of the progenitor at collapse using the scaling given by NFW.
    collapseOverdensity        =self%C*(expansionFactor/collapseExpansionFactor)**3
    ! Find the ratio of this overdensity to that at for the present node.
    rootTarget                 =collapseOverdensity/self%virialDensityContrast_%densityContrast(nodeMass,nodeTime)
     ! Find the concentration.
    nfw1996Concentration       =self%finder%find(rootRange=[1.0d0,20.0d0])
    return
  end function nfw1996Concentration

  double precision function nfw1996ConcentrationRoot(concentration)
    !!{
    Root function used in finding concentrations in the \cite{navarro_structure_1996} method.
    !!}
    implicit none
    double precision, intent(in   ) :: concentration

    nfw1996ConcentrationRoot=concentration**3/(log(1.0d0+concentration)-concentration/(1.0d0+concentration))/3.0d0-rootTarget
    return
  end function nfw1996ConcentrationRoot

  function nfw1996DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of
    concentration in the \cite{navarro_structure_1996} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass           ), pointer       :: nfw1996DensityContrastDefinition
    class(darkMatterProfileConcentrationNfw1996), intent(inout) :: self

    nfw1996DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function nfw1996DensityContrastDefinition

  function nfw1996DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{navarro_structure_1996} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass            ), pointer       :: nfw1996DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationNFW1996), intent(inout) :: self

    nfw1996DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function nfw1996DarkMatterProfileDefinition

