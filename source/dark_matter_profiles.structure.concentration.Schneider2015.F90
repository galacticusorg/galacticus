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
  An implementation of halo profile concentrations using the algorithm of \cite{schneider_structure_2015}.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Root_Finder               , only : rootFinder

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationSchneider2015">
   <description>
    A dark matter profile concentration class in which the concentration using the algorithm of
    \cite{schneider_structure_2015}. Specifically, a reference model for concentrations in defined in a specific cosmological
    model. The concentration for a halo of given mass, $M$, and redshift, $z_0$, is then found by finding the redshift of
    collapse, $z_\mathrm{c}$ for the halo by solving:
    \begin{equation}
      \delta_\mathrm{c}(z_\mathrm{c}) = \left( {\pi \over 2} \left[ \sigma^2(f M) - \sigma^2(M) \right]
      \right)^{1/2}+\delta_\mathrm{c})(z_0),
    \end{equation}
    where $\delta_\mathrm{c}(z)$ is the critical overdensity for collapse at redshift $z$, and $f$ is the fraction of a halo's
    mass assembled at formation time (given by the {\normalfont \ttfamily [massFractionFormation]} parameter. From this, the
    mass of a halo in the reference model with the same redshift of collapse is found, and the reference model is used to
    compute the concentration of a halo of that mass.
   </description>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationSchneider2015
     !!{
     A dark matter halo profile concentration class implementing the algorithm of \cite{schneider_structure_2015}.
     !!}
     private
     class           (darkMatterProfileConcentrationClass), pointer :: referenceConcentration            => null()
     class           (cosmologyFunctionsClass            ), pointer :: referenceCosmologyFunctions       => null(), cosmologyFunctions_       => null()
     class           (criticalOverdensityClass           ), pointer :: referenceCriticalOverdensity      => null(), criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass      ), pointer :: referenceCosmologicalMassVariance => null(), cosmologicalMassVariance_ => null()
     type            (rootFinder                         )          :: finder
     double precision                                               :: massFractionFormation
   contains
     !![
     <methods>
     <method method="concentrationCompute" description="Compute the concentration for the given {\normalfont \ttfamily node}"/>
     </methods>
     !!]
    final     ::                                   schneider2015Destructor
    procedure :: concentration                  => schneider2015Concentration
    procedure :: concentrationMean              => schneider2015ConcentrationMean
    procedure :: concentrationCompute           => schneider2015ConcentrationCompute
    procedure :: densityContrastDefinition      => schneider2015DensityContrastDefinition
    procedure :: darkMatterProfileDMODefinition => schneider2015DarkMatterProfileDefinition
 end type darkMatterProfileConcentrationSchneider2015

  interface darkMatterProfileConcentrationSchneider2015
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationSchneider2015} dark matter halo profile concentration class.
     !!}
     module procedure schneider2015ConstructorParameters
     module procedure schneider2015ConstructorInternal
  end interface darkMatterProfileConcentrationSchneider2015

  ! Module-scope variables for root finding.
  class           (darkMatterProfileConcentrationSchneider2015), pointer  :: self_
  double precision                                                        :: massReferencePrevious        , timeCollapseReference                  , &
       &                                                                     time_                        , referenceCollapseMassRootPrevious      , &
       &                                                                     timeNowReference
  !$omp threadprivate(self_,massReferencePrevious,timeCollapseReference,time_,timeNowReference,referenceCollapseMassRootPrevious)

  ! Upper limit to the reference mass used during root finding.
  double precision                                             , parameter :: massReferenceMinimum=1.0d-30, massReferenceMaximum             =1.0d30

contains

 function schneider2015ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily schneider2015} dark matter halo profile concentration class.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileConcentrationSchneider2015)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    type            (inputParameters                            )                :: referenceParameters
    class           (darkMatterProfileConcentrationClass        ), pointer       :: referenceConcentration
    class           (criticalOverdensityClass                   ), pointer       :: referenceCriticalOverdensity     , criticalOverdensity_
    class           (cosmologicalMassVarianceClass              ), pointer       :: referenceCosmologicalMassVariance, cosmologicalMassvariance_
    class           (cosmologyFunctionsClass                    ), pointer       :: referenceCosmologyFunctions      , cosmologyFunctions_
    double precision                                                             :: massFractionFormation

    if (.not.parameters%isPresent('reference',requireValue=.false.)) call Error_Report('parameters must contain a "reference" section'//{introspection:location})
    referenceParameters=parameters%subParameters('reference',requireValue=.false.)
    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massFractionFormation</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <description>The fraction of a halo's mass assembled at ``formation'' in the halo concentration algorithm of \cite{schneider_structure_2015}.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileConcentration" name="referenceConcentration"             source="referenceParameters"/>
    <objectBuilder class="criticalOverdensity"            name="referenceCriticalOverdensity"       source="referenceParameters"/>
    <objectBuilder class="criticalOverdensity"            name="         criticalOverdensity_"      source=         "parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="referenceCosmologicalMassVariance"  source="referenceParameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="         cosmologicalMassVariance_" source=         "parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="referenceCosmologyFunctions"        source="referenceParameters"/>
    <objectBuilder class="cosmologyFunctions"             name="         cosmologyFunctions_"       source=         "parameters"/>
    !!]
    self=darkMatterProfileConcentrationSchneider2015(massFractionFormation,referenceConcentration,referenceCriticalOverdensity,referenceCosmologicalMassVariance,referenceCosmologyFunctions,criticalOverdensity_,cosmologicalMassvariance_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"         />
    <inputParametersValidate source="referenceParameters"/>
    <objectDestructor name="referenceConcentration"           />
    <objectDestructor name="referenceCriticalOverdensity"     />
    <objectDestructor name="criticalOverdensity_"             />
    <objectDestructor name="referenceCosmologicalMassVariance"/>
    <objectDestructor name="cosmologicalMassVariance_"        />
    <objectDestructor name="referenceCosmologyFunctions"      />
    <objectDestructor name="cosmologyFunctions_"              />
    !!]
    return
  end function schneider2015ConstructorParameters

  function schneider2015ConstructorInternal(massFractionFormation,referenceConcentration,referenceCriticalOverdensity,referenceCosmologicalMassVariance,referenceCosmologyFunctions,criticalOverdensity_,cosmologicalMassvariance_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileConcentrationSchneider2015} dark matter halo concentration class.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative
    implicit none
    type            (darkMatterProfileConcentrationSchneider2015)                        :: self
    class           (darkMatterProfileConcentrationClass        ), intent(in   ), target :: referenceConcentration
    class           (criticalOverdensityClass                   ), intent(in   ), target :: referenceCriticalOverdensity           , criticalOverdensity_
    class           (cosmologicalMassVarianceClass              ), intent(in   ), target :: referenceCosmologicalMassVariance      , cosmologicalMassvariance_
    class           (cosmologyFunctionsClass                    ), intent(in   ), target :: referenceCosmologyFunctions            , cosmologyFunctions_
    double precision                                             , intent(in   )         :: massFractionFormation
    double precision                                             , parameter             :: toleranceAbsolute                =0.0d0, toleranceRelative        =1.0d-6
    !![
    <constructorAssign variables="massFractionFormation, *referenceConcentration, *referenceCriticalOverdensity, *referenceCosmologicalMassVariance, *referenceCosmologyFunctions, *criticalOverdensity_, *cosmologicalMassvariance_, *cosmologyFunctions_"/>
    !!]

    self%finder=rootFinder(                                                            &
         &                 rootFunction       =schneider2015ReferenceCollapseMassRoot, &
         &                 toleranceAbsolute  =toleranceAbsolute                     , &
         &                 toleranceRelative  =toleranceRelative                     , &
         &                 rangeExpandUpward  =2.000d0                               , &
         &                 rangeExpandDownward=0.999d0                               , &
         &                 rangeExpandType    =rangeExpandMultiplicative             , &
         &                 rangeUpwardLimit   =massReferenceMaximum                  , &
         &                 rangeDownwardLimit =massReferenceMinimum                    &
         &                )
    return
  end function schneider2015ConstructorInternal

  subroutine schneider2015Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationSchneider2015} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationSchneider2015), intent(inout) :: self

    !![
    <objectDestructor name="self%referenceConcentration"           />
    <objectDestructor name="self%referenceCriticalOverdensity"     />
    <objectDestructor name="self%criticalOverdensity_"             />
    <objectDestructor name="self%referenceCosmologicalMassVariance"/>
    <objectDestructor name="self%cosmologicalMassVariance_"        />
    <objectDestructor name="self%referenceCosmologyFunctions"      />
    <objectDestructor name="self%cosmologyFunctions_"              />
    !!]
    return
  end subroutine schneider2015Destructor

  double precision function schneider2015Concentration(self,node) result(concentration)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the algorithm of
    \cite{schneider_structure_2015}.
    !!}
    implicit none
    class(darkMatterProfileConcentrationSchneider2015), intent(inout), target :: self
    type (treeNode                                   ), intent(inout), target :: node

    concentration=self%concentrationCompute(node,mean=.false.)
    return
  end function schneider2015Concentration
  
  double precision function schneider2015ConcentrationMean(self,node) result(concentration)
    !!{
    Return the mean concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the algorithm of
    \cite{schneider_structure_2015}.
    !!}
    implicit none
    class(darkMatterProfileConcentrationSchneider2015), intent(inout)         :: self
    type (treeNode                                   ), intent(inout), target :: node

    concentration=self%concentrationCompute(node,mean=.true.)
    return
  end function schneider2015ConcentrationMean
  
  double precision function schneider2015ConcentrationCompute(self,node,mean) result(concentration)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the algorithm of
    \cite{schneider_structure_2015}.
    !!}
    use :: Error                   , only : Error_Report      , errorStatusSuccess
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileConcentrationSchneider2015), intent(inout), target  :: self
    type            (treeNode                                   ), intent(inout), target  :: node
    logical                                                      , intent(in   )          :: mean
    class           (nodeComponentBasic                         )               , pointer :: basic
    double precision                                             , parameter              :: toleranceTimeRelative      =1.0d-6
    integer                                                                               :: status
    double precision                                                                      :: mass                                            , &
         &                                                                                   collapseCriticalOverdensity       , timeCollapse, &
         &                                                                                   massReference                     , variance    , &
         &                                                                                   timeNow

    ! Get the basic component and the halo mass and time.
    basic => node %basic()
    mass  =  basic%mass ()
    time_ =  basic%time ()
    ! Find critical overdensity at collapse for this node. The critical overdensity at collapse is scaled by a factor
    ! σ(M,t₀)/σ(M,t), as the "timeOfCollapse" function expects a critical overdensity divided by the linear growth factor. (This
    ! is because the definition of critical overdensity in Galacticus does not include the 1/D(t) factor that is included in the
    ! definition used by Schneider et al. (2015).)
    timeNow                    =     self%cosmologyFunctions_        %cosmicTime  (1.0d0                                )
    timeNowReference           =     self%referenceCosmologyFunctions%cosmicTime  (1.0d0                                )
    variance                   =max(                                                                                          &
         &                          +0.0d0                                                                                  , &
         &                          +self%cosmologicalMassVariance_  %rootVariance(mass*self%massFractionFormation,time_)**2  &
         &                          -self%cosmologicalMassVariance_  %rootVariance(mass                           ,time_)**2  &
         &                         )
    collapseCriticalOverdensity=+(                                                                         &
         &                        +sqrt(                                                                   &
         &                              +Pi                                                                &
         &                              /2.0d0                                                             &
         &                              *variance                                                          &
         &                             )                                                                   &
         &                        +  self%criticalOverdensity_       %value       (time=time_  ,mass=mass) &
         &                       )                                                                         &
         &                      *    self%cosmologicalMassVariance_  %rootVariance(time=timeNow,mass=mass) &
         &                      /    self%cosmologicalMassVariance_  %rootVariance(time=time_  ,mass=mass)

    ! Compute the corresponding epoch of collapse.
    timeCollapse         =self%criticalOverdensity_%timeOfCollapse(collapseCriticalOverdensity,mass)
    ! Compute time of collapse in the reference model, assuming same redshift of collapse in both models.
    timeCollapseReference=self%referenceCosmologyFunctions%cosmicTime(self%cosmologyFunctions_%expansionFactor(timeCollapse))
    ! Check for a collapse time very close to the present time.
    if (timeCollapseReference > time_*(1.0d0-toleranceTimeRelative)) then
       ! Collapse is happening at the current time - which will correspond to arbitrarily massive halos in the reference model.
       massReference=massReferenceMaximum
    else
       ! Find the mass of a halo collapsing at the same time in the reference model.
       self_                 => self
       massReferencePrevious =  -1.0d0
       if (schneider2015ReferenceCollapseMassRoot(massReferenceMaximum) > -time_*toleranceTimeRelative) then
          ! No solution can be found even at the maximum allowed mass. Simply set the reference mass to the maximum allowed mass -
          ! the choice should not matter too much as the abundances of such halos should be hugely suppressed.
          massReference        =massReferenceMaximum
       else
          massReference        =self%finder%find(rootGuess=mass,status=status)
          if (status /= errorStatusSuccess) call Error_Report('failed to determine reference mass'//{introspection:location})
       end if
   end if
    ! Compute the concentration of a node of this mass in the reference model.
    call basic%massSet(massReference)
    if (mean) then
       concentration=self%referenceConcentration%concentrationMean(node)
    else
       concentration=self%referenceConcentration%concentration    (node)
    end if
    call basic%massSet(mass         )
    return
  end function schneider2015ConcentrationCompute

  double precision function schneider2015ReferenceCollapseMassRoot(massReference)
    !!{
    Root function used to find the mass collapsing at given time in dark matter halo concentration algorithm of
    \cite{schneider_structure_2015}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: massReference
    double precision                :: variance     , collapseCriticalOverdensity
    
    if (massReference /= massReferencePrevious) then
       massReferencePrevious            =+massReference
       variance                         =max(                                                                                                           &
            &                                +0.0d0                                                                                                   , &
            &                                +self_%referenceCosmologicalMassVariance%rootVariance(massReference*self_%massFractionFormation,time_)**2  &
            &                                -self_%referenceCosmologicalMassVariance%rootVariance(massReference                            ,time_)**2  &
            &                               )
       ! Find critical overdensity at collapse for this node. The critical overdensity at collapse is scaled by a factor
       ! σ(M,t₀)/σ(M,t), as the "timeOfCollapse" function expects a critical overdensity divided by the linear growth factor. (This
       ! is because the definition of critical overdensity in Galacticus does not include the 1/D(t) factor that is included in the
       ! definition used by Schneider et al. (2015).)
       collapseCriticalOverdensity=+(                                                                                                &
            &                        +sqrt(                                                                                          &
            &                              +Pi                                                                                       &
            &                              /2.0d0                                                                                    &
            &                              *variance                                                                                 &
            &                             )                                                                                          &
            &                        +self_%referenceCriticalOverdensity     %value       (time=time_           ,mass=massReference) &
            &                       )                                                                                                &
            &                      *  self_%referenceCosmologicalMassVariance%rootVariance(time=timeNowReference,mass=massReference) &
            &                      /  self_%referenceCosmologicalMassVariance%rootVariance(time=time_           ,mass=massReference)
       ! Construct the root to get collapse at the required collapse time in the reference model.
       referenceCollapseMassRootPrevious=-self_%referenceCriticalOverdensity%timeOfCollapse       (collapseCriticalOverdensity,massReference) &
            &                            +                                   timeCollapseReference
    end if
    schneider2015ReferenceCollapseMassRoot=referenceCollapseMassRootPrevious    
    return
  end function schneider2015ReferenceCollapseMassRoot

  function schneider2015DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{schneider_structure_2015} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass                 ), pointer       :: schneider2015DensityContrastDefinition
    class(darkMatterProfileConcentrationSchneider2015), intent(inout) :: self

    schneider2015DensityContrastDefinition => self%referenceConcentration%densityContrastDefinition()
    return
  end function schneider2015DensityContrastDefinition

  function schneider2015DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{schneider_structure_2015} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass                  ), pointer       :: schneider2015DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationSchneider2015), intent(inout) :: self

    schneider2015DarkMatterProfileDefinition => self%referenceConcentration%darkMatterProfileDMODefinition()
    return
  end function schneider2015DarkMatterProfileDefinition
