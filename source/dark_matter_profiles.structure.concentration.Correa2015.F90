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
  An implementation of dark matter halo profile concentrations using the \cite{correa_accretion_2015} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMONFW
  use :: Linear_Growth             , only : linearGrowthClass
  use :: Virial_Density_Contrast   , only : virialDensityContrastFixed

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationCorrea2015">
   <description>Dark matter halo concentrations are computed using the algorithm of \cite{correa_accretion_2015}.</description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationCorrea2015
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{correa_accretion_2015}.
     !!}
     private
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_             => null()
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_              => null()
     class           (linearGrowthClass            ), pointer :: linearGrowth_                    => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_        => null()
     type            (virialDensityContrastFixed   ), pointer :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMONFW      ), pointer :: darkMatterProfileDMODefinition_  => null()
     double precision                                         :: A
   contains
     final     ::                                   correa2015Destructor
     procedure :: concentration                  => correa2015Concentration
     procedure :: densityContrastDefinition      => correa2015DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => correa2015DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationCorrea2015

  interface darkMatterProfileConcentrationCorrea2015
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationCorrea2015} dark matter halo profile concentration class.
     !!}
     module procedure correa2015ConstructorParameters
     module procedure correa2015ConstructorInternal
  end interface darkMatterProfileConcentrationCorrea2015

contains

  function correa2015ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily correa2015} dark matter halo profile concentration class.
    !!}
    implicit none
    type            (darkMatterProfileConcentrationCorrea2015)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (cosmologyParametersClass                ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass                       ), pointer       :: linearGrowth_
    class           (cosmologicalMassVarianceClass           ), pointer       :: cosmologicalMassVariance_
    double precision                                                          :: A

    !![
    <inputParameter>
      <name>A</name>
      <source>parameters</source>
      <variable>A</variable>
      <defaultValue>887.0d0</defaultValue>
      <defaultSource>\cite{correa_accretion_2015}</defaultSource>
      <description>The parameter $A$ appearing in eqn.~(17) of \cite{correa_accretion_2015}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationCorrea2015(A,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function correa2015ConstructorParameters

  function correa2015ConstructorInternal(A,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,cosmologicalMassVariance_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationCorrea2015} dark matter halo profile concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    implicit none
    type            (darkMatterProfileConcentrationCorrea2015          )                         :: self
    double precision                                                    , intent(in   )          :: A
    class           (cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class           (cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    class           (linearGrowthClass                                 ), intent(in   ), target  :: linearGrowth_
    class           (cosmologicalMassVarianceClass                     ), intent(in   ), target  :: cosmologicalMassVariance_
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="A, *cosmologyParameters_, *cosmologyFunctions_, *linearGrowth_, *cosmologicalMassVariance_"/>
    !!]

    allocate(self%virialDensityContrastDefinition_)
    allocate(self%darkMatterProfileDMODefinition_ )
    allocate(     darkMatterHaloScaleDefinition_  )
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
  end function correa2015ConstructorInternal

  subroutine correa2015Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationCorrea2015} dark matter profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationCorrea2015), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%linearGrowth_"                   />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine correa2015Destructor

  double precision function correa2015Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite{correa_accretion_2015} algorithm.
    !!}
    use :: Dark_Matter_Halos_Correa2015, only : Dark_Matter_Halo_Correa2015_Fit_Parameters
    use :: Galacticus_Nodes            , only : nodeComponentBasic                        , treeNode
    use :: Root_Finder                 , only : rangeExpandMultiplicative                 , rootFinder
    implicit none
    class           (darkMatterProfileConcentrationCorrea2015), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    class           (nodeComponentBasic                      )               , pointer :: basic
    double precision                                          , parameter              :: toleranceRelative=1.0d-6
    double precision                                          , parameter              :: toleranceAbsolute=0.0d+0
    type            (rootFinder                              ), save                   :: finder
    logical                                                   , save                   :: finderConstructed=.false.
    !$omp threadprivate(finder)
    double precision                                                                   :: mass                     , time           , &
         &                                                                                aTilde                   , bTilde         , &
         &                                                                                redshift                 , expansionFactor

    ! Get properties of the node.
    basic => node %basic()
    mass  =  basic%mass ()
    time  =  basic%time ()
    ! Determine the base redshift.
    expansionFactor=self%cosmologyFunctions_%expansionFactor            (time           )
    redshift       =self%cosmologyFunctions_%redshiftFromExpansionFactor(expansionFactor)
    ! Find the a~ and b~ parameters.
    call Dark_Matter_Halo_Correa2015_Fit_Parameters(mass,expansionFactor,self%cosmologyFunctions_,self%linearGrowth_,self%cosmologicalMassVariance_,aTilde,bTilde)
    ! Solve for the redshift corresponding to this mass.
    if (.not.finderConstructed) then
       finder=rootFinder(                                               &
            &            rootFunction       =concentrationSolver      , &
            &            toleranceAbsolute  =toleranceAbsolute        , &
            &            toleranceRelative  =toleranceRelative        , &
            &            rangeExpandDownward=0.5d0                    , &
            &            rangeExpandUpward  =2.0d0                    , &
            &            rangeExpandType    =rangeExpandMultiplicative  &
            &           )
       finderConstructed=.true.
    end if
    correa2015Concentration=finder%find(rootGuess=6.0d0)
    return

  contains

    double precision function concentrationSolver(concentration)
      !!{
      Solver used in finding halo concentration using the \cite{correa_accretion_2015} algorithm.
      !!}
      implicit none
      double precision, intent(in   ) :: concentration
      double precision                :: Y1           , Yc               , &
           &                             f1           , f2               , &
           &                             densityInner , redshiftFormation

      ! Evaluate left-hand side of eqn. 18 of Correa et al. (2015).
      Y1                 =log(2.0d0              )-0.5d0
      Yc                 =log(1.0d0+concentration)-concentration/(1.0d0+concentration)
      f1                 =+Y1 &
           &              /Yc
      ! Evaluate mean inner density.
      densityInner       =+200.0d0          &
           &              *concentration**3 &
           &              *Y1               &
           &              /Yc
      ! Evaluate eqn. 17 of Correa et al. (2015), rearranged to give formation redshift.
      redshiftFormation  =+(                                                              &
           &                +(                                                            &
           &                  +(                                                          &
           &                    +1.0d0                                                    &
           &                    +redshift                                                 &
           &                   )                                          **3             &
           &                  +self%cosmologyParameters_%omegaDarkEnergy()                &
           &                  /self%cosmologyParameters_%omegaMatter    ()                &
           &                 )                                                            &
           &                *densityInner                                                 &
           &                /  self%A                                                     &
           &                -  self%cosmologyParameters_%omegaDarkEnergy()                &
           &                /  self%cosmologyParameters_%omegaMatter    ()                &
           &               )                                              **(1.0d0/3.0d0) &
           &              -1.0d0
      ! Evaluate right-hand side of eqn. 19 of Correa et al. (2015).
      f2                 =+(                        &
           &                +1.0d0                  &
           &                +redshiftFormation      &
           &                -redshift               &
           &               )**aTilde                &
           &              *exp(                     &
           &                   +(                   &
           &                     +redshiftFormation &
           &                     -redshift          &
           &                    )                   &
           &                   *bTilde              &
           &                  )
      ! Left- and right-hand sides will be equal for correct concentration.
      concentrationSolver=+f1 &
           &              -f2
      return
    end function concentrationSolver

  end function correa2015Concentration

  function correa2015DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{correa_accretion_2015} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass              ), pointer       :: correa2015DensityContrastDefinition
    class(darkMatterProfileConcentrationCorrea2015), intent(inout) :: self
    !$GLC attributes unused :: self

   correa2015DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function correa2015DensityContrastDefinition

  function correa2015DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{correa_accretion_2015} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass               ), pointer       :: correa2015DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationCorrea2015), intent(inout) :: self

    correa2015DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function correa2015DarkMatterProfileDefinition
