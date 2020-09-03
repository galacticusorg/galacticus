!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module implementing a 

module Node_Component_Scale_Virial_Theorem
  !% Implements a
  use Galacticus_Nodes                  , only : treeNode                         , nodeComponentDarkMatterProfile  
  use Dark_Matter_Profile_Scales        , only : darkMatterProfileScaleRadiusClass
  use Dark_Matter_Profiles              , only : darkMatterProfileClass
  use Virial_Orbits                     , only : virialOrbitClass
  use Dark_Matter_Halo_Scales           , only : darkMatterHaloScaleClass
  use Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Vrl_Thrm_Rate_Compute     , Node_Component_Dark_Matter_Profile_Vrl_Thrm_Tree_Initialize    , &
       &    Node_Component_Dark_Matter_Profile_Vrl_Thrm_Initialize       , Node_Component_Dark_Matter_Profile_Vrl_Thrm_Scale_Set          , &
       &    Node_Component_Dark_Matter_Profile_Vrl_Thrm_Thread_Initialize, Node_Component_Dark_Matter_Profile_Vrl_Thrm_Thread_Uninitialize, &
       &    Node_Component_Dark_Matter_Profile_Vrl_Thrm_Plausibility

  !# <component>
  !#  <class>darkMatterProfile</class>
  !#  <name>virialTheorem</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>scale</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="megaParsec" comment="Scale radius of the dark matter profile [Mpc]."/>
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>scaleGrowthRate</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>scaleIsLimited</name>
  !#     <type>logical</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>.true.</classDefault>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Module-scope variables used in solving for the scale radius as a function of energy.
  type            (treeNode                         ), pointer :: nodeActive
  class           (nodeComponentDarkMatterProfile   ), pointer :: darkMatterProfileActive
  double precision                                             :: energyTotal
  !$omp threadprivate(nodeActive,darkMatterProfileActive,energyTotal)
  
  ! Objects used by this component.
  class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_
  class           (darkMatterProfileClass           ), pointer :: darkMatterProfile_
  class           (virialOrbitClass                 ), pointer :: virialOrbit_
  class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_
  class           (mergerTreeMassResolutionClass    ), pointer :: mergerTreeMassResolution_
  !$omp threadprivate(darkMatterProfileScaleRadius_,darkMatterProfile_,virialOrbit_,darkMatterHaloScale_,mergerTreeMassResolution_)
  
  ! Parameter controlling scaling of orbital energy with mass ratio.
  double precision :: darkMatterProfileScaleVirialTheoremMassExponent, darkMatterProfileScaleVirialTheoremEnergyBoost, darkMatterProfileScaleVirialTheoremUnresolvedEnergy
  
contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Vrl_Thrm_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Initialize(globalParameters_)
    !% Initializes the ``virialTheorem'' implementation of the dark matter profile component.
    use Input_Parameters
    implicit none
    type(inputParameters), intent(inout) :: globalParameters_

    !# <inputParameter>
    !#   <name>darkMatterProfileScaleVirialTheoremEnergyBoost</name>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <source>globalParameters_</source>
    !#   <description>A boost to the energy</description>
    !#   <type>double</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>darkMatterProfileScaleVirialTheoremMassExponent</name>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <source>globalParameters_</source>
    !#   <description>The exponent of mass ratio appearing in the orbital energy term in the ``virial theorem'' dark matter profile scale model.</description>
    !#   <type>double</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>darkMatterProfileScaleVirialTheoremUnresolvedEnergy</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <source>globalParameters_</source>
    !#   <description>Factor multiplying the estimate of the internal energy of unresolved accretion in the ``virial theorem'' dark matter profile scale model.</description>
    !#   <type>double</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
   return
  end subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Initialize
  
  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Vrl_Thrm_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Thread_Initialize(parameters)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks                      , only : nodePromotionEvent               , openMPThreadBindingAtLevel
    use :: Input_Parameters                  , only : inputParameters                  , inputParameter
    use :: Galacticus_Nodes                  , only : defaultDarkMatterProfileComponent
    use :: Dark_Matter_Profile_Scales        , only : darkMatterProfileScaleRadius
    use :: Dark_Matter_Profiles              , only : darkMatterProfile
    use :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScale
    use :: Virial_Orbits                     , only : virialOrbit
    use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolution
    implicit none
    type(inputParameters), intent(inout) :: parameters

    if (defaultDarkMatterProfileComponent%virialTheoremIsActive()) then
       !# <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
       !# <objectBuilder class="darkMatterProfile"            name="darkMatterProfile_"            source="parameters"/>
       !# <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
       !# <objectBuilder class="virialOrbit"                  name="virialOrbit_"                  source="parameters"/>
       !# <objectBuilder class="mergerTreeMassResolution"     name="mergerTreeMassResolution_"     source="parameters"/>
       call nodePromotionEvent%attach(defaultDarkMatterProfileComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentDarkMatterProfileVirialTheorem')
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Vrl_Thrm_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none

    if (defaultDarkMatterProfileComponent%virialTheoremIsActive()) then
       !# <objectDestructor name="darkMatterProfileScaleRadius_"/>
       !# <objectDestructor name="darkMatterProfile_"           />
       !# <objectDestructor name="darkMatterHaloScale_"         />
       !# <objectDestructor name="virialOrbit_"                 />
       !# <objectDestructor name="mergerTreeMassResolution_"    />
        call nodePromotionEvent%detach(defaultDarkMatterProfileComponent,nodePromotion)
   end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Thread_Uninitialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Vrl_Thrm_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Compute the rate of change of the scale radius.
    use Galacticus_Nodes, only : treeNode, nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileVirialTheorem, propertyTypeInactive
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (                              ), intent(inout), pointer :: interruptProcedure
    integer                                         , intent(in   )          :: propertyType
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: growthRate
    !$GLC attributes unused :: interrupt, interruptProcedure
    
    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileVirialTheorem)
       growthRate=darkMatterProfile%scaleGrowthRate()
       call darkMatterProfile%scaleRate(growthRate)
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Rate_Compute

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Vrl_Thrm_Plausibility</unitName>
  !#  <after>Node_Component_Basic_Standard_Plausibility</after>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Plausibility(node)
    !% Determines whether the dark matter profile is physically plausible for radius solving tasks.
    use Galacticus_Nodes, only : treeNode, nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileVirialTheorem
    implicit none
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    ! Return immediately if already non-plausible.
    if (.not.(node%isPhysicallyPlausible.and.node%isSolvable)) return
    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileVirialTheorem)
       if (darkMatterProfile%scale() <= 0.0d0) node%isPhysicallyPlausible=.false.
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Plausibility

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Vrl_Thrm_Tree_Initialize</unitName>
  !#  <sortName>darkMatterProfile</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Tree_Initialize(node)
    !% Initialize the scale radius of {\normalfont \ttfamily node}.
    use Galacticus_Nodes                , only : nodeComponentBasic                         , nodeComponentDarkMatterProfile     , nodeComponentSatellite       , defaultDarkMatterProfileComponent, &
         &                                       nodeComponentDarkMatterProfileVirialTheorem
    use Root_Finder                     , only : rootFinder                                 , rangeExpandMultiplicative          , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    use Kepler_Orbits                   , only : keplerOrbit
    use Merger_Tree_Walkers             , only : mergerTreeWalkerAllNodes
    use Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use Beta_Functions                  , only : Beta_Function                              , Beta_Function_Incomplete_Normalized
    use Hypergeometric_Functions
    implicit none
    type            (treeNode                      ), intent(inout), pointer   :: node
    type            (treeNode                      )               , pointer   :: nodeChild                                      , nodeSibling                     , &
         &                                                                        nodeWork                                       , nodeUnresolved
    class           (nodeComponentBasic            )               , pointer   :: basicChild                                     , basicSibling                    , &
         &                                                                        basic                                          , basicParent                     , &
         &                                                                        basicUnresolved
    class           (nodeComponentDarkMatterProfile)               , pointer   :: darkMatterProfile                              , darkMatterProfileSibling        , &
         &                                                                        darkMatterProfileChild                         , darkMatterProfileParent         , &
         &                                                                        darkMatterProfileUnresolved
    class           (nodeComponentSatellite        )               , pointer   :: satelliteSibling
    double precision                                               , parameter :: massFunctionSlopeLogarithmic            =1.90d0
    double precision                                               , parameter :: energyInternalFormFactorSlopeLogarithmic=0.03d0
    double precision                                                           :: energyOrbital                                  , massRatio                       , &
         &                                                                        radiusScaleChild                               , radiusScale                     , &
         &                                                                        massUnresolved                                 , deltaTime                       , &
         &                                                                        radiusScaleUnresolved                          , massResolution                  , &
         &                                                                        energyPotentialSubresolutionFactor             , energyKineticSubresolutionFactor, &
         &                                                                        a                                              , b                               , &
         &                                                                        energyPotential                                , energyKinetic                   , &
         &                                                                        energyInternalSubresolutionFactor
    type            (rootFinder                    )                           :: finder
    type            (keplerOrbit                   )                           :: orbit
    type            (mergerTreeWalkerAllNodes      )                           :: treeWalker
    
    ! Check if we are the default method.
    if (defaultDarkMatterProfileComponent%virialTheoremIsActive()) then
       ! Get the darkMatterProfile component.
       darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
       if (darkMatterProfile%scale() <= 0.0d0) then
          ! Get the resolution of the tree.
          massResolution=mergerTreeMassResolution_%resolution(node%hostTree)
          ! Perform our own depth-first tree walk to set scales in all nodes of the tree. This is necessary as we require access
          ! to the parent scale to set scale growth rates, but must initialize scales in a strictly depth-first manner as some
          ! algorithms rely on knowing the progenitor structure of the tree to compute scale radii.
          treeWalker=mergerTreeWalkerAllNodes(node%hostTree,spanForest=.false.)
          do while (treeWalker%next(nodeWork))
             darkMatterProfile => nodeWork%darkMatterProfile(autoCreate=.true.)
             ! If this node has no children, draw its dark matter profile scale radius from a distribution.
             if (.not.associated(nodeWork%firstChild)) then
                radiusScale=darkMatterProfileScaleRadius_%radius(nodeWork)
             else
                basic                   =>  nodeWork              %basic            ()
                nodeChild               =>  nodeWork              %firstChild
                darkMatterProfileChild  =>  nodeChild             %darkMatterProfile()
                basicChild              =>  nodeChild             %basic            ()
                radiusScaleChild        =   darkMatterProfileChild%scale            ()
                massUnresolved          =  +basic                 %mass             () &
                     &                     -basicChild            %mass             ()
                ! Iterate over progenitors and sum their energies.
                nodeSibling =>                           nodeChild
                energyTotal =  darkMatterProfile_%energy(nodeSibling)
                do while (associated(nodeSibling%sibling))
                   nodeSibling              =>  nodeSibling       %sibling
                   basicSibling             =>  nodeSibling       %basic            (                      )
                   darkMatterProfileSibling =>  nodeSibling       %darkMatterProfile(                      )
                   satelliteSibling         =>  nodeSibling       %satellite        (autoCreate=.true.     )
                   orbit                    =   satelliteSibling  %virialOrbit      (                      )
                   massRatio                =  +basicSibling      %mass             (                      ) &
                        &                      /basicChild        %mass             (                      )
                   energyOrbital            =  +orbit             %energy           (                      ) &
                        &                      *basicSibling      %mass             (                      ) &
                        &                      /(                                                            &
                        &                        +1.0d0                                                      &
                        &                        +massRatio                                                  &
                        &                       )**darkMatterProfileScaleVirialTheoremMassExponent
                   massUnresolved           =  +massUnresolved                                               &
                        &                      -basicSibling      %mass             (                      )
                   ! Add orbital energy of this sibling.
                   energyTotal              = +energyTotal                                                   &
                        &                     +energyOrbital                                                 &
                        &                     *(                                                             &
                        &                       +1.0d0                                                       &
                        &                       +darkMatterProfileScaleVirialTheoremEnergyBoost              &
                        &                      )
                   ! Add the internal energy of the sibling.
                   energyTotal              =  +energyTotal                                                  &
                        &                      +darkMatterProfile_%energy           (           nodeSibling) &
                        &                      /(                                                            &
                        &                        +1.0d0                                                      &
                        &                        +massRatio                                                  &
                        &                       )**darkMatterProfileScaleVirialTheoremMassExponent           &
                        &                      *(                                                            &
                        &                        +1.0d0                                                      &
                        &                        +darkMatterProfileScaleVirialTheoremEnergyBoost             &
                        &                       )
                end do
                ! Account for unresolved accretion. We assume that unresolved halos are accreted with the mean orbital energy of
                ! the virial orbital parameter distribution, plus an internal energy corresponding to that of a halo with mass
                ! equal to the total unresolved mass scaled by some correction factor (to account for the fact that the unresolved
                ! accretion will not in fact be in a single halo).
                nodeUnresolved              => treeNode                                       (                         )
                basicUnresolved             => nodeUnresolved               %basic            (autoCreate=.true.        )
                darkMatterProfileUnresolved => nodeUnresolved               %darkMatterProfile(autoCreate=.true.        )
                call basicUnresolved            %massSet            (min(massResolution   ,massUnresolved       ))
                call basicUnresolved            %timeSet            (    basicChild%time()                       )
                call basicUnresolved            %timeLastIsolatedSet(    basicChild%time()                       )
                radiusScaleUnresolved       =  darkMatterProfileScaleRadius_%radius           (           nodeUnresolved)
                call darkMatterProfileUnresolved%scaleSet           (                      radiusScaleUnresolved )
                ! Compute a correction factor to the orbital energy which takes into account the mass dependence of the 1/(1+m/M)ᵅ
                ! term that is applied to the orbital energy. Averaging this over a power-law mass function gives the result
                ! below. In the case that α=0 the result is identically 1 - in this case we avoid computing beta functions.
                massRatio                       =+basicUnresolved%mass()                                               &
                     &                           /basicChild     %mass()
                a                               =+2.0d0                                                                &
                     &                           -massFunctionSlopeLogarithmic
                b                               =+massFunctionSlopeLogarithmic                                         &
                     &                           +darkMatterProfileScaleVirialTheoremMassExponent                      &
                     &                           -1.0d0
                energyKineticSubresolutionFactor=+Beta_Function_Incomplete_Normalized(a,b,massRatio/(1.0d0+massRatio)) &
                     &                           *Beta_Function                      (a,b                            ) &
                     &                           *           (2.0d0-massFunctionSlopeLogarithmic)                      &
                     &                           /massRatio**(2.0d0-massFunctionSlopeLogarithmic)
                if (darkMatterProfileScaleVirialTheoremMassExponent == 0.0d0) then
                   energyPotentialSubresolutionFactor=+1.0d0
                else
                   b                                 =+massFunctionSlopeLogarithmic                                         &
                        &                             +darkMatterProfileScaleVirialTheoremMassExponent                      &
                        &                             -2.0d0
                   energyPotentialSubresolutionFactor=+Beta_Function_Incomplete_Normalized(a,b,massRatio/(1.0d0+massRatio)) &
                        &                             *Beta_Function                      (a,b                            ) &
                        &                             *           (2.0d0-massFunctionSlopeLogarithmic)                      &
                        &                             /massRatio**(2.0d0-massFunctionSlopeLogarithmic)
                end if
                ! Compute a correction factor to the internal energy which takes into account the mass dependence of the 1/(1+m/M)ᵅ
                ! term that is applied to the orbital energy. Averaging this over a power-law mass function gives the result
                ! below.
                energyInternalSubresolutionFactor=+(                                                                                                                                                          &
                     &                              +2.0d0                                                                                                                                                    &
                     &                              -(massFunctionSlopeLogarithmic+energyInternalFormFactorSlopeLogarithmic)                                                                                  &
                     &                             )                                                                                                                                                          &
                     &                            *Hypergeometric_2F1(                                                                                                                                        &
                     &                                                [darkMatterProfileScaleVirialTheoremMassExponent, 8.0d0/3.0d0-(massFunctionSlopeLogarithmic+energyInternalFormFactorSlopeLogarithmic)], &
                     &                                                [                                                11.0d0/3.0d0-(massFunctionSlopeLogarithmic+energyInternalFormFactorSlopeLogarithmic)], &
                     &                                                -1.0d0/massRatio                                                                                                                        &
                     &                                               )                                                                                                                                        &
                     &                            /(                                                                                                                                                          &
                     &                                                                                                  8.0d0/3.0d0-(massFunctionSlopeLogarithmic+energyInternalFormFactorSlopeLogarithmic)   &
                     &                              )
                ! Determine the orbital and internal energies.
                energyKinetic  =+0.5d0                                                                  &
                     &          *virialOrbit_%velocityTotalRootMeanSquared(nodeUnresolved,nodeChild)**2 &
                     &          /(1.0d0+massRatio)
                energyPotential=+virialOrbit_%energyMean                  (nodeUnresolved,nodeChild)    &
                     &          -energyKinetic
                energyOrbital  =+energyPotential                                       &
                     &          *energyPotentialSubresolutionFactor                    &
                     &          +energyKinetic                                         &
                     &          *energyKineticSubresolutionFactor
                energyTotal    =+energyTotal                                           &
                     &          +massUnresolved                                        &
                     &          *darkMatterProfileScaleVirialTheoremUnresolvedEnergy   &
                     &          *(                                                     &
                     &            +energyOrbital                                       &
                     &            +energyInternalSubresolutionFactor                   &
                     &            *darkMatterProfile_%energy(nodeUnresolved)           &
                     &            /massResolution                                      &
                     &           )                                                     &
                     &          *(                                                     &
                     &            +1.0d0                                               &
                     &            +darkMatterProfileScaleVirialTheoremEnergyBoost      &
                     &           )
                ! Add mutual gravitational binding energy of any sibling halo and any unresolved mass.
                if (associated(nodeChild%sibling))                                      &
                     & energyTotal=+energyTotal                                         &
                     &             -gravitationalconstantGalacticus                     &
                     &             *basicSibling        %mass        (         )        &
                     &             *massUnresolved                                      &
                     &             /0.5d0                                               &
                     &             /darkMatterHaloScale_%virialRadius(nodeChild)        &
                     &             *darkMatterProfileScaleVirialTheoremUnresolvedEnergy &
                     &             *(                                                   &
                     &               +1.0d0                                             &
                     &               +darkMatterProfileScaleVirialTheoremEnergyBoost    &
                     &              )
                call nodeUnresolved%destroy()
                deallocate(nodeUnresolved)
                ! Convert energy back to scale radius.
                nodeActive              => nodeWork
                darkMatterProfileActive => darkMatterProfile
                if (.not.finder%isInitialized()) then
                   call finder%rootFunction(radiusScaleRoot                                  )
                   call finder%tolerance   (toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
                   call finder%rangeExpand (                                                             &
                        &                   rangeExpandDownward          =0.5d0                        , &
                        &                   rangeExpandUpward            =2.0d0                        , &
                        &                   rangeExpandType              =rangeExpandMultiplicative    , &
                        &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
                        &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative  &
                        &                  )
                end if
                radiusScale=finder%find(rootGuess=radiusScaleChild)
             end if
             call darkMatterProfile%scaleSet(radiusScale)             
          end do
          ! Walk the tree again to set growth rates.
          treeWalker=mergerTreeWalkerAllNodes(node%hostTree,spanForest=.false.)
          do while (treeWalker%next(nodeWork))
             darkMatterProfile => nodeWork%darkMatterProfile()
             ! Check if this node is the primary progenitor.
             if (nodeWork%isPrimaryProgenitor()) then
                ! It is, so compute the scale radius growth rate.
                ! Now compute the growth rate.
                basic       =>  nodeWork          %basic()
                basicParent =>  nodeWork   %parent%basic()
                deltaTime   =  +basicParent       %time () &
                     &         -basic             %time ()
                if (deltaTime > 0.0d0) then
                   darkMatterProfileParent => nodeWork%parent%darkMatterProfile()
                   call darkMatterProfile%scaleGrowthRateSet(                                  &
                        &                                    (                                 &
                        &                                     +darkMatterProfileParent%scale() &
                        &                                     -darkMatterProfile      %scale() &
                        &                                    )                                 &
                        &                                    /deltaTime                        &
                        &                                   )
                else
                   call darkMatterProfile%scaleGrowthRateSet(0.0d0)
                end if
             else
                ! It is not, so set scale radius growth rate to zero.
                call    darkMatterProfile%scaleGrowthRateSet(0.0d0)
             end if
          end do
       end if
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Tree_Initialize

  double precision function radiusScaleRoot(radiusScale)
    !% Function used in root-finding to compute the scale radius of a dark matter profile as a given energy.
    implicit none
    double precision, intent(in   ) :: radiusScale

    call darkMatterProfileActive%scaleSet(radiusScale)
    radiusScaleRoot=+energyTotal                           &
         &          -darkMatterProfile_%energy(nodeActive)
    return
  end function radiusScaleRoot
  
  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the growth rate of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use Galacticus_Error, only : Galacticus_Error_Report
    use Galacticus_Nodes, only : treeNode               , nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileVirialTheorem, nodeComponentBasic
    implicit none
    class(*                             ), intent(inout)          :: self
    type (treeNode                      ), intent(inout), target  :: node
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfileParent, darkMatterProfile
    class(nodeComponentBasic            )               , pointer :: basicParent            , basic
    !$GLC attributes unused :: self

    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileVirialTheorem)
       darkMatterProfileParent => node%parent%darkMatterProfile()
       basic                   => node       %basic            ()
       basicParent             => node%parent%basic            ()
       if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
       ! Adjust the scale radius to that of the parent node.
       call darkMatterProfile%scaleSet          (darkMatterProfileParent%scale          ())
       ! Adjust the growth rate to that of the parent node.
       call darkMatterProfile%scaleGrowthRateSet(darkMatterProfileParent%scaleGrowthRate())
    end select
    return
  end subroutine nodePromotion

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Vrl_Thrm_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use Galacticus_Nodes, only : treeNode, nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileVirialTheorem
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile

    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (darkMatterProfile)
       class is (nodeComponentDarkMatterProfileVirialTheorem)
       ! Set scale for the scale radius.
       call darkMatterProfile%scaleScale(darkMatterProfile%scale())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Vrl_Thrm_Scale_Set

end module Node_Component_Scale_Virial_Theorem
