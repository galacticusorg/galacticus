!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  use Cosmological_Mass_Variance
  use Critical_Overdensities
  use Cosmology_Functions
  use Root_Finder
  
  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationSchneider2015">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{schneider_structure_2015}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationSchneider2015
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{schneider_structure_2015}.
     private
     class           (darkMatterProfileConcentrationClass), pointer :: referenceConcentration
     class           (cosmologyFunctionsClass            ), pointer :: referenceCosmologyFunctions      , cosmologyFunctions_
     class           (criticalOverdensityClass           ), pointer :: referenceCriticalOverdensity     , criticalOverdensity_
     class           (cosmologicalMassVarianceClass      ), pointer :: referenceCosmologicalMassVariance, cosmologicalMassVariance_
     type            (rootFinder                         )          :: finder
     double precision                                               :: massFractionFormation
  contains
     procedure :: concentration => schneider2015Concentration
  end type darkMatterProfileConcentrationSchneider2015

  interface darkMatterProfileConcentrationSchneider2015
     !% Constructors for the {\normalfont \ttfamily Schneider2015} dark matter halo profile concentration class.
     module procedure schneider2015ConstructorParameters
     module procedure schneider2015ConstructorInternal
  end interface darkMatterProfileConcentrationSchneider2015

contains

  function schneider2015ConstructorParameters(parameters)
    !% Default constructor for the {\normalfont \ttfamily schneider2015} dark matter halo profile concentration class.
    use Input_Parameters2
    use Galacticus_Error
    implicit none
    type (darkMatterProfileConcentrationSchneider2015)                :: schneider2015ConstructorParameters
    type(inputParameters                             ), intent(inout) :: parameters
    type(inputParameters                             )                :: referenceParameters
    !# <inputParameterList label="allowedParameterNames"          source="parameters"          />
    !# <inputParameterList label="allowedReferenceParameterNames" source="referenceParameters" />

    if (.not.parameters%isPresent('reference',requireValue=.false.)) call Galacticus_Error_Report('schneider2015ConstructorParameters','parameters must contain a "reference" section')
    referenceParameters=parameters%subParameters('reference',requireValue=.false.)
    ! Check and read parameters.
    call          parameters%checkParameters(allowedParameterNames         )    
    call referenceParameters%checkParameters(allowedReferenceParameterNames)    
    !# <inputParameter>
    !#   <name>massFractionFormation</name>
    !#   <source>parameters</source>
    !#   <variable>schneider2015ConstructorParameters%massFractionFormation</variable>    
    !#   <defaultValue>0.05d0</defaultValue>
    !#   <description>The fraction of a halo's mass assembled at ``formation'' in the halo concentration algorithm of \cite{schneider_structure_2015}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>    
    !# <objectBuilder class="darkMatterProfileConcentration" name="schneider2015ConstructorParameters%referenceConcentration"             source="referenceParameters"/>
    !# <objectBuilder class="criticalOverdensity"            name="schneider2015ConstructorParameters%referenceCriticalOverdensity"       source="referenceParameters"/>
    !# <objectBuilder class="criticalOverdensity"            name="schneider2015ConstructorParameters%         criticalOverdensity_"      source=         "parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance"       name="schneider2015ConstructorParameters%referenceCosmologicalMassVariance"  source="referenceParameters"/>
    !# <objectBuilder class="cosmologicalMassVariance"       name="schneider2015ConstructorParameters%         cosmologicalMassVariance_" source=         "parameters"/>
    !# <objectBuilder class="cosmologyFunctions"             name="schneider2015ConstructorParameters%referenceCosmologyFunctions"        source="referenceParameters"/>
    !# <objectBuilder class="cosmologyFunctions"             name="schneider2015ConstructorParameters%         cosmologyFunctions_"       source=         "parameters"/>
    return
  end function schneider2015ConstructorParameters

  function schneider2015ConstructorInternal(referenceConcentration,referenceCriticalOverdensity,referenceCosmologicalMassVariance,referenceCosmologyFunctions,criticalOverdensity_,cosmologicalMassvariance_,cosmologyFunctions_)
    !% Generic constructor for the {\normalfont \ttfamily schneider2015} dark matter halo concentration class.
    implicit none
    type (darkMatterProfileConcentrationSchneider2015)                        :: schneider2015ConstructorInternal
    class(darkMatterProfileConcentrationClass        ), intent(in   ), target :: referenceConcentration
    class(criticalOverdensityClass                   ), intent(in   ), target :: referenceCriticalOverdensity     , criticalOverdensity_
    class(cosmologicalMassVarianceClass              ), intent(in   ), target :: referenceCosmologicalMassVariance, cosmologicalMassvariance_
    class(cosmologyFunctionsClass                    ), intent(in   ), target :: referenceCosmologyFunctions      , cosmologyFunctions_

    ! Construct the object.
    schneider2015ConstructorInternal%referenceConcentration             => referenceConcentration
    schneider2015ConstructorInternal%referenceCriticalOverdensity       => referenceCriticalOverdensity
    schneider2015ConstructorInternal%         criticalOverdensity_      =>          criticalOverdensity_
    schneider2015ConstructorInternal%referenceCosmologicalMassVariance  => referenceCosmologicalMassVariance
    schneider2015ConstructorInternal%         cosmologicalMassVariance_ =>          cosmologicalMassVariance_
    schneider2015ConstructorInternal%referenceCosmologyFunctions        => referenceCosmologyFunctions
    schneider2015ConstructorInternal%         cosmologyFunctions_       =>          cosmologyFunctions_
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
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileConcentrationSchneider2015), intent(inout)          :: self
    type            (treeNode                                   ), intent(inout), pointer :: node
    class           (nodeComponentBasic                         )               , pointer :: basic
    double precision                                             , parameter              :: toleranceAbsolute          =0.0d0, toleranceRelative    =1.0d-6
    double precision                                                                      :: mass                             , time                        , &
         &                                                                                   collapseCriticalOverdensity      , timeCollapse                , &
         &                                                                                   massReference                    , timeCollapseReference
      
    ! Get the basic component and the halo mass and time.
    basic => node %basic()
    mass  =  basic%mass ()
    time  =  basic%time ()
    ! Find critical overdensity at collapse for this node.
    collapseCriticalOverdensity=+sqrt(                                                                                   &
         &                            +Pi                                                                                &
         &                            /2.0d0                                                                             &
         &                            *(                                                                                 &
         &                              +self%cosmologicalMassVariance_%rootVariance(mass*self%massFractionFormation)**2 &
         &                              -self%cosmologicalMassVariance_%rootVariance(mass                           )**2 &
         &                             )                                                                                 &
         &                           )                                                                                   &
         &                      +self%criticalOverdensity_%value(time=time,mass=mass)
    ! Compute the corresponding epoch of collapse.
    timeCollapse=self%criticalOverdensity_%timeOfCollapse(collapseCriticalOverdensity,mass)
    ! Compute time of collapse in the reference model, assuming same redshift of collapse in both models.
    timeCollapseReference=self%referenceCosmologyFunctions%cosmicTime(self%cosmologyFunctions_%expansionFactor(timeCollapse))
    ! Find the mass of a halo collapsing at the same time in the reference model.
    if (.not.self%finder%isInitialized()) then
       call self%finder%rootFunction(referenceCollapseMassRoot                  )
       call self%finder%tolerance   (toleranceAbsolute        ,toleranceRelative)
       call self%finder%rangeExpand (                                                   &
            &                        rangeExpandUpward  =2.0d0                        , &
            &                        rangeExpandDownward=0.5d0                        , &
            &                        rangeExpandType    =rangeExpandMultiplicative      &
            &                       )
    end if
    massReference=self%finder%find(rootGuess=mass)
    ! Compute the concentration of a node of this mass in the reference model.
    call basic%massSet(massReference)
    schneider2015Concentration=self%referenceConcentration%concentration(node)
    call basic%massSet(mass         )
    return

  contains

    double precision function referenceCollapseMassRoot(massReference)
      !% Root function used to find the mass collapsing at given time in dark matter halo concentration algorithm of
      !% \cite{schneider_structure_2015}.
      implicit none
      double precision, intent(in   ) :: massReference

      referenceCollapseMassRoot=+sqrt(                                                                                                    &
         &                            +Pi                                                                                                 &
         &                            /2.0d0                                                                                              &
         &                            *(                                                                                                  &
         &                              +self%referenceCosmologicalMassVariance%rootVariance(massReference*self%massFractionFormation)**2 &
         &                              -self%referenceCosmologicalMassVariance%rootVariance(massReference                           )**2 &
         &                             )                                                                                                  &
         &                           )                                                                                                    &
         &                      +self%referenceCriticalOverdensity%value(time=time                 ,mass=massReference)                   &
         &                      -self%referenceCriticalOverdensity%value(time=timeCollapseReference,mass=massReference)
      return
    end function referenceCollapseMassRoot
    
  end function schneider2015Concentration

