!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% An implementation of halo profile concentrations using the algorithm of \cite{schneider_structure_2015}.

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Root_Finder               , only : rootFinder

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationSchneider2015">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{schneider_structure_2015}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationSchneider2015
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{schneider_structure_2015}.
     private
     class           (darkMatterProfileConcentrationClass), pointer :: referenceConcentration            => null()
     class           (cosmologyFunctionsClass            ), pointer :: referenceCosmologyFunctions       => null(), cosmologyFunctions_       => null()
     class           (criticalOverdensityClass           ), pointer :: referenceCriticalOverdensity      => null(), criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass      ), pointer :: referenceCosmologicalMassVariance => null(), cosmologicalMassVariance_ => null()
     type            (rootFinder                         )          :: finder
     double precision                                               :: massFractionFormation
  contains
    final     ::                                   schneider2015Destructor
    procedure :: concentration                  => schneider2015Concentration
    procedure :: densityContrastDefinition      => schneider2015DensityContrastDefinition
    procedure :: darkMatterProfileDMODefinition => schneider2015DarkMatterProfileDefinition
 end type darkMatterProfileConcentrationSchneider2015

  interface darkMatterProfileConcentrationSchneider2015
     !% Constructors for the {\normalfont \ttfamily Schneider2015} dark matter halo profile concentration class.
     module procedure schneider2015ConstructorParameters
     module procedure schneider2015ConstructorInternal
  end interface darkMatterProfileConcentrationSchneider2015

  ! Module-scope variables for root finding.
  class           (darkMatterProfileConcentrationSchneider2015), pointer :: schneider2015Self
  double precision                                                       :: schneider2015MassReferencePrevious, schneider2015TimeCollapseReference            , &
       &                                                                    schneider2015Time                 , schneider2015ReferenceCollapseMassRootPrevious
  !$omp threadprivate(schneider2015Self,schneider2015MassReferencePrevious,schneider2015TimeCollapseReference,schneider2015Time,schneider2015ReferenceCollapseMassRootPrevious)

contains

 function schneider2015ConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily schneider2015} dark matter halo profile concentration class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type            (darkMatterProfileConcentrationSchneider2015)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    type            (inputParameters                            )                :: referenceParameters
    class           (darkMatterProfileConcentrationClass        ), pointer       :: referenceConcentration
    class           (criticalOverdensityClass                   ), pointer       :: referenceCriticalOverdensity     , criticalOverdensity_
    class           (cosmologicalMassVarianceClass              ), pointer       :: referenceCosmologicalMassVariance, cosmologicalMassvariance_
    class           (cosmologyFunctionsClass                    ), pointer       :: referenceCosmologyFunctions      , cosmologyFunctions_
    double precision                                                             :: massFractionFormation

    if (.not.parameters%isPresent('reference',requireValue=.false.)) call Galacticus_Error_Report('parameters must contain a "reference" section'//{introspection:location})
    referenceParameters=parameters%subParameters('reference',requireValue=.false.)
    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>massFractionFormation</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.05d0</defaultValue>
    !#   <description>The fraction of a halo's mass assembled at ``formation'' in the halo concentration algorithm of \cite{schneider_structure_2015}.</description>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterProfileConcentration" name="referenceConcentration"             source="referenceParameters"/>
    !# <objectBuilder class="criticalOverdensity"            name="referenceCriticalOverdensity"       source="referenceParameters"/>
    !# <objectBuilder class="criticalOverdensity"            name="         criticalOverdensity_"      source=         "parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance"       name="referenceCosmologicalMassVariance"  source="referenceParameters"/>
    !# <objectBuilder class="cosmologicalMassVariance"       name="         cosmologicalMassVariance_" source=         "parameters"/>
    !# <objectBuilder class="cosmologyFunctions"             name="referenceCosmologyFunctions"        source="referenceParameters"/>
    !# <objectBuilder class="cosmologyFunctions"             name="         cosmologyFunctions_"       source=         "parameters"/>
    self=darkMatterProfileConcentrationSchneider2015(massFractionFormation,referenceConcentration,referenceCriticalOverdensity,referenceCosmologicalMassVariance,referenceCosmologyFunctions,criticalOverdensity_,cosmologicalMassvariance_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"         />
    !# <inputParametersValidate source="referenceParameters"/>
    !# <objectDestructor name="referenceConcentration"           />
    !# <objectDestructor name="referenceCriticalOverdensity"     />
    !# <objectDestructor name="criticalOverdensity_"             />
    !# <objectDestructor name="referenceCosmologicalMassVariance"/>
    !# <objectDestructor name="cosmologicalMassVariance_"        />
    !# <objectDestructor name="referenceCosmologyFunctions"      />
    !# <objectDestructor name="cosmologyFunctions_"              />
    return
  end function schneider2015ConstructorParameters

  function schneider2015ConstructorInternal(massFractionFormation,referenceConcentration,referenceCriticalOverdensity,referenceCosmologicalMassVariance,referenceCosmologyFunctions,criticalOverdensity_,cosmologicalMassvariance_,cosmologyFunctions_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily schneider2015} dark matter halo concentration class.
    implicit none
    type            (darkMatterProfileConcentrationSchneider2015)                        :: self
    class           (darkMatterProfileConcentrationClass        ), intent(in   ), target :: referenceConcentration
    class           (criticalOverdensityClass                   ), intent(in   ), target :: referenceCriticalOverdensity     , criticalOverdensity_
    class           (cosmologicalMassVarianceClass              ), intent(in   ), target :: referenceCosmologicalMassVariance, cosmologicalMassvariance_
    class           (cosmologyFunctionsClass                    ), intent(in   ), target :: referenceCosmologyFunctions      , cosmologyFunctions_
    double precision                                             , intent(in   )         :: massFractionFormation
    !# <constructorAssign variables="massFractionFormation, *referenceConcentration, *referenceCriticalOverdensity, *referenceCosmologicalMassVariance, *referenceCosmologyFunctions, *criticalOverdensity_, *cosmologicalMassvariance_, *cosmologyFunctions_"/>

    self%finder=rootFinder()
    return
  end function schneider2015ConstructorInternal

  subroutine schneider2015Destructor(self)
    !% Destructor for the {\normalfont \ttfamily schneider2015} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationSchneider2015), intent(inout) :: self

    !# <objectDestructor name="self%referenceConcentration"           />
    !# <objectDestructor name="self%referenceCriticalOverdensity"     />
    !# <objectDestructor name="self%criticalOverdensity_"             />
    !# <objectDestructor name="self%referenceCosmologicalMassVariance"/>
    !# <objectDestructor name="self%cosmologicalMassVariance_"        />
    !# <objectDestructor name="self%referenceCosmologyFunctions"      />
    !# <objectDestructor name="self%cosmologyFunctions_"              />
    return
  end subroutine schneider2015Destructor

  double precision function schneider2015Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the algorithm of
    !% \cite{schneider_structure_2015}.
    use :: Galacticus_Nodes        , only : nodeComponentBasic       , treeNode
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rangeExpandMultiplicative
    implicit none
    class           (darkMatterProfileConcentrationSchneider2015), intent(inout), target  :: self
    type            (treeNode                                   ), intent(inout), target  :: node
    class           (nodeComponentBasic                         )               , pointer :: basic
    double precision                                             , parameter              :: toleranceAbsolute          =0.0d00, toleranceRelative=1.0d-6, &
         &                                                                                   massReferenceMaximum       =1.0d20
    double precision                                                                      :: mass                                                        , &
         &                                                                                   collapseCriticalOverdensity       , timeCollapse            , &
         &                                                                                   massReference                     , variance                , &
         &                                                                                   timeNow

    ! Get the basic component and the halo mass and time.
    basic             => node %basic()
    mass              =  basic%mass ()
    schneider2015Time =  basic%time ()
    ! Find critical overdensity at collapse for this node. The critical overdensity at collapse is scaled by a factor
    ! σ(M,t₀)/σ(M,t) because the definition of critical overdensity in Galacticus does not include the 1/D(t) factor that is
    ! included in the definition used by Schneider et al. (2015).
    timeNow =self%cosmologyFunctions_%cosmicTime(1.0d0)
    variance=max(                                                                                                    &
         &       +0.0d0                                                                                            , &
         &       +self%cosmologicalMassVariance_%rootVariance(mass*self%massFractionFormation,schneider2015Time)**2  &
         &       -self%cosmologicalMassVariance_%rootVariance(mass                           ,schneider2015Time)**2  &
         &      )
    collapseCriticalOverdensity=+(                                                                               &
         &                        +sqrt(                                                                         &
         &                              +Pi                                                                      &
         &                              /2.0d0                                                                   &
         &                              *variance                                                                &
         &                             )                                                                         &
         &                        +self%criticalOverdensity_     %value       (time=schneider2015Time,mass=mass) &
         &                       )                                                                               &
         &                      *  self%cosmologicalMassVariance_%rootVariance(time=timeNow          ,mass=mass) &
         &                      /  self%cosmologicalMassVariance_%rootVariance(time=schneider2015Time,mass=mass)
    ! Compute the corresponding epoch of collapse.
    timeCollapse                      =self%criticalOverdensity_%timeOfCollapse(collapseCriticalOverdensity,mass)
    ! Compute time of collapse in the reference model, assuming same redshift of collapse in both models.
    schneider2015TimeCollapseReference=self%referenceCosmologyFunctions%cosmicTime(self%cosmologyFunctions_%expansionFactor(timeCollapse))
    ! Find the mass of a halo collapsing at the same time in the reference model.
    if (.not.self%finder%isInitialized()) then
       call self%finder%rootFunction(schneider2015ReferenceCollapseMassRoot)
       call self%finder%tolerance   (toleranceAbsolute,toleranceRelative)
       call self%finder%rangeExpand (                                                   &
            &                        rangeExpandUpward  =2.000d0                      , &
            &                        rangeExpandDownward=0.999d0                      , &
            &                        rangeExpandType    =rangeExpandMultiplicative    , &
            &                        rangeUpwardLimit   =massReferenceMaximum           &
            &                       )
    end if
    schneider2015Self                  => self
    schneider2015MassReferencePrevious =  -1.0d0
    if (schneider2015ReferenceCollapseMassRoot(massReferenceMaximum) > 0.0d0) then
       ! No solution can be found even at the maximum allowed mass. Simply set the reference mass to the maximum allowed mass -
       ! the choice shouldn't matter too much as the abundances of such halos should be hugely suppressed.
       massReference        =massReferenceMaximum
    else
       massReference        =self%finder%find(rootGuess=mass)
    end if
    ! Compute the concentration of a node of this mass in the reference model.
    call basic%massSet(massReference)
    schneider2015Concentration=self%referenceConcentration%concentration(node)
    call basic%massSet(mass         )
    return
  end function schneider2015Concentration

  double precision function schneider2015ReferenceCollapseMassRoot(massReference)
    !% Root function used to find the mass collapsing at given time in dark matter halo concentration algorithm of
    !% \cite{schneider_structure_2015}.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: massReference
    double precision                :: variance

    if (massReference /= schneider2015MassReferencePrevious) then
       schneider2015MassReferencePrevious            =+massReference
       variance                                      =max(                                                                                                                                               &
            &                                             +0.0d0                                                                                                                                       , &
            &                                             +schneider2015Self%referenceCosmologicalMassVariance%rootVariance(massReference*schneider2015Self%massFractionFormation,schneider2015Time)**2  &
            &                                             -schneider2015Self%referenceCosmologicalMassVariance%rootVariance(massReference                                        ,schneider2015Time)**2  &
            &                                            )
       schneider2015ReferenceCollapseMassRootPrevious=+sqrt(                                                                                                                                             &
            &                                               +Pi                                                                                                                                          &
            &                                               /2.0d0                                                                                                                                       &
            &                                               *variance                                                                                                                                    &
            &                                              )                                                                                                                                             &
            &                                         +schneider2015Self%referenceCriticalOverdensity%value(time=schneider2015Time                 ,mass=massReference)                                  &
            &                                         -schneider2015Self%referenceCriticalOverdensity%value(time=schneider2015TimeCollapseReference,mass=massReference)
    end if
    schneider2015ReferenceCollapseMassRoot=schneider2015ReferenceCollapseMassRootPrevious
    return
  end function schneider2015ReferenceCollapseMassRoot

  function schneider2015DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the
    !% \cite{schneider_structure_2015} algorithm.
    implicit none
    class(virialDensityContrastClass                 ), pointer       :: schneider2015DensityContrastDefinition
    class(darkMatterProfileConcentrationSchneider2015), intent(inout) :: self

    schneider2015DensityContrastDefinition => self%referenceConcentration%densityContrastDefinition()
    return
  end function schneider2015DensityContrastDefinition

  function schneider2015DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{schneider_structure_2015} algorithm.
    implicit none
    class(darkMatterProfileDMOClass                  ), pointer       :: schneider2015DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationSchneider2015), intent(inout) :: self

    schneider2015DarkMatterProfileDefinition => self%referenceConcentration%darkMatterProfileDMODefinition()
    return
  end function schneider2015DarkMatterProfileDefinition
