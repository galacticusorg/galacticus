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
  Implements a node operator class that implements starvation of subhalos by removal of their \gls{cgm}.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <nodeOperator name="nodeOperatorCGMStarvation">
   <description>
    A node operator class that implements starvation of subhalos by removal of their \gls{cgm}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMStarvation
     !!{
     A node operator class that implements starvation of subhalos by removal of their \gls{cgm}.
     !!}
     private
     class  (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     logical                                    :: starveOutflowsOnly            , fractionBaryonLimitInNodeMerger
   contains
     final     ::                              cgmStarvationDestructor
     procedure :: nodesMerge                => cgmStarvationNodesMerge
     procedure :: differentialEvolutionPost => cgmStarvationDifferentialEvolutionPost
  end type nodeOperatorCGMStarvation
  
  interface nodeOperatorCGMStarvation
     !!{
     Constructors for the \refClass{nodeOperatorCGMStarvation} node operator class.
     !!}
     module procedure cgmStarvationConstructorParameters
     module procedure cgmStarvationConstructorInternal
  end interface nodeOperatorCGMStarvation
  
contains
  
  function cgmStarvationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMStarvation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorCGMStarvation)                :: self
    type   (inputParameters          ), intent(inout) :: parameters
    class  (cosmologyParametersClass ), pointer       :: cosmologyParameters_
    logical                                           :: starveOutflowsOnly  , fractionBaryonLimitInNodeMerger
    
    !![
    <inputParameter>
      <name>starveOutflowsOnly</name>
      <defaultValue>.false.</defaultValue>
      <description>
	If true, only starve outflowing gas by transferring it to the \gls{cgm} of the host halo.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionBaryonLimitInNodeMerger</name>
      <defaultValue>.false.</defaultValue>
      <description>
	Controls whether the \gls{cgm} gas content of nodes should be limited to not exceed the universal baryon fraction at node
        merger events. If set to {\normalfont \ttfamily true}, \gls{cgm} hot gas (and angular momentum, abundances, and chemicals
        proportionally) will be removed from the merged halo to the unaccreted gas reservoir to limit the baryonic mass to the
        universal baryon fraction where possible.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=nodeOperatorCGMStarvation(starveOutflowsOnly,fractionBaryonLimitInNodeMerger,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>  
    !!]
    return
  end function cgmStarvationConstructorParameters

  function cgmStarvationConstructorInternal(starveOutflowsOnly,fractionBaryonLimitInNodeMerger,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMStarvation} node operator class.
    !!}
    implicit none
    type   (nodeOperatorCGMStarvation)                        :: self
    class  (cosmologyParametersClass ), intent(in   ), target :: cosmologyParameters_
    logical                           , intent(in   )         :: starveOutflowsOnly  , fractionBaryonLimitInNodeMerger
    !![
    <constructorAssign variables="starveOutflowsOnly, fractionBaryonLimitInNodeMerger, *cosmologyParameters_"/>
    !!]

    return
  end function cgmStarvationConstructorInternal

  subroutine cgmStarvationDestructor(self)    
    !!{
    Destructor for the \refClass{nodeOperatorCGMStarvation} class.
    !!}
    implicit none
    type(nodeOperatorCGMStarvation), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine cgmStarvationDestructor

  subroutine cgmStarvationNodesMerge(self,node)
    !!{
    Update the \gls{cgm} content of a node as a result of a merger.
    !!}
    use :: Abundances_Structure         , only : zeroAbundances
    use :: Accretion_Halos              , only : accretionModeHot      , accretionModeTotal
    use :: Chemical_Abundances_Structure, only : zeroChemicalAbundances
    use :: Galactic_Structure_Options   , only : componentTypeAll      , massTypeBaryonic
    use :: Galacticus_Nodes             , only : nodeComponentBasic    , nodeComponentHotHalo, nodeComponentSpin
    use :: Mass_Distributions           , only : massDistributionClass
    use :: Error                        , only : Error_Report 
    implicit none
    class           (nodeOperatorCGMStarvation), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: node
    type            (treeNode                 ), pointer       :: nodeParent
    class           (nodeComponentHotHalo     ), pointer       :: hotHaloParent      , hotHalo
    class           (nodeComponentBasic       ), pointer       :: basicParent        , basic
    class           (nodeComponentSpin        ), pointer       :: spinParent
    class           (massDistributionClass    ), pointer       :: massDistribution_
    double precision                                           :: baryonicMassCurrent, baryonicMassMaximum, &
         &                                                        fractionRemove        

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       ! Find the parent node and its hot halo component.
       nodeParent    => node      %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       ! Get the basic components.
       basic       => node      %basic()
       basicParent => nodeParent%basic()
       spinParent  => nodeParent%spin ()
       if (                                                     .not.self%starveOutflowsOnly) then
          call        hotHaloParent%                    massSet(                                          &
               &                                                 hotHaloParent%mass                    () &
               &                                                +hotHalo      %mass                    () &
               &                                               )
          if (hotHaloParent%angularMomentumIsSettable())                                                  &
               & call hotHaloParent%         angularMomentumSet(                                          &
               &                                                 hotHaloParent%angularMomentum         () &
               &                                                +hotHalo      %mass                    () &
               &                                                *spinParent   %angularMomentum         () &
               &                                                /basicParent  %mass                    () &
               &                                               )
       end if
       call           hotHaloParent%           outflowedMassSet(                                          &
            &                                                    hotHaloParent%outflowedMass           () &
            &                                                   +hotHalo      %outflowedMass           () &
            &                                                  )
       if (hotHaloParent%outflowedAngularMomentumIsSettable()                                     )       &
            & call    hotHaloParent%outflowedAngularMomentumSet(                                          &
            &                                                    hotHaloParent%outflowedAngularMomentum() &
            &                                                   +hotHalo      %outflowedMass           () &
            &                                                   *spinParent   %angularMomentum         () &
            &                                                   /basicParent  %mass                    () &
            &                                                  )
       if (                                                     .not.self%starveOutflowsOnly) then
          call        hotHalo      %                    massSet(                                          &
               &                                                 0.0d0                                    &
               &                                               )
          if (hotHalo%            angularMomentumIsSettable()                                    )        &
               & call hotHalo      %         angularMomentumSet(                                          &
               &                                                 0.0d0                                    &
               &                                               )
       end if
       call           hotHalo      %           outflowedMassSet(                                          &
            &                                                    0.0d0                                    &
            &                                                  )
       if (hotHalo% outflowedAngularMomentumIsSettable()                                         )        &
            & call    hotHalo      %outflowedAngularMomentumSet(                                          &
            &                                                    0.0d0                                    &
            &                                                  )
       if (                                                     .not.self%starveOutflowsOnly) then
          call hotHaloParent%              abundancesSet(                                          &
               &                                          hotHaloParent%abundances              () &
               &                                         +hotHalo      %abundances              () &
               &                                        )
          call hotHalo      %              abundancesSet(                                          &
               &                                          zeroAbundances                           &
               &                                        )
       end if
       if (hotHaloParent%outflowedAbundancesIsSettable()                                    ) then
          call hotHaloParent%     outflowedAbundancesSet(                                          &
               &                                          hotHaloParent%outflowedAbundances     () &
               &                                         +hotHalo      %outflowedAbundances     () &
               &                                        )
          call hotHalo      %     outflowedAbundancesSet(                                          &
               &                                          zeroAbundances                           &
               &                                        )
       end if
       if (hotHaloParent%           chemicalsIsSettable() .and. .not.self%starveOutflowsOnly) then
          call hotHaloParent%               chemicalsSet(                                          &
               &                                          hotHaloParent%chemicals               () &
               &                                         +hotHalo      %chemicals               () &
               &                                        )
          call hotHalo      %               chemicalsSet(                                          &
               &                                          zeroChemicalAbundances                   &
               &                                        )
       end if
       if (hotHaloParent% outflowedChemicalsIsSettable()                                    ) then
          call hotHaloParent%      outflowedChemicalsSet(                                          &
               &                                          hotHaloParent%outflowedChemicals      () &
               &                                         +hotHalo      %outflowedChemicals      () &
               &                                        )
          call hotHalo      %      outflowedChemicalsSet(                                          &
               &                                          zeroChemicalAbundances                   &
               &                                        )
       end if
       ! Check if the baryon fraction in the parent hot halo exceeds the universal value. If it does, mitigate this by moving
       ! some of the mass to the failed accretion reservoir.
       if (self%fractionBaryonLimitInNodeMerger) then
          ! Get the default cosmology.
          massDistribution_   =>  nodeParent                            %massDistribution(massType=massTypeBaryonic)
          baryonicMassMaximum =  +basicParent                           %mass            (                         ) &
               &                 *self             %cosmologyParameters_%omegaBaryon     (                         ) &
               &                 /self             %cosmologyParameters_%omegaMatter     (                         )
          baryonicMassCurrent =  +massDistribution_                     %massTotal       (                         )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
          if (baryonicMassCurrent > baryonicMassMaximum .and. hotHaloParent%mass() > 0.0d0) then
             fractionRemove=min((baryonicMassCurrent-baryonicMassMaximum)/hotHaloParent%massTotal(),1.0d0)
             call hotHaloParent%      unaccretedMassSet(                                                             &
                  &                                      hotHaloParent%unaccretedMass      ()                        &
                  &                                     +hotHaloParent%mass                ()*       fractionRemove  &
                  &                                    )
             call hotHaloParent%unaccretedAbundancesSet(                                                             &
                  &                                      hotHaloParent%unaccretedAbundances()                        &
                  &                                     +hotHaloParent%abundances          ()*       fractionRemove  &
                  &                                    )
             call hotHaloParent%                massSet( hotHaloParent%mass                ()*(1.0d0-fractionRemove))
             call hotHaloParent%     angularMomentumSet( hotHaloParent%angularMomentum     ()*(1.0d0-fractionRemove))
             call hotHaloParent%          abundancesSet( hotHaloParent%abundances          ()*(1.0d0-fractionRemove))
             call hotHaloParent%          chemicalsSet ( hotHaloParent%chemicals           ()*(1.0d0-fractionRemove))
          end if
       end if
    end select
    return
  end subroutine cgmStarvationNodesMerge

  subroutine cgmStarvationDifferentialEvolutionPost(self,node)
    !!{
    Starve subhalos at the end of each differential evolution step.
    !!}
    use :: Abundances_Structure                 , only : zeroAbundances
    use :: Chemical_Abundances_Structure        , only : zeroChemicalAbundances 
    use :: Galacticus_Nodes                     , only : nodeComponentHotHalo  , nodeComponentSpin, nodeComponentBasic
    implicit none
    class(nodeOperatorCGMStarvation), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    type (treeNode                 ), pointer       :: nodeParent
    class(nodeComponentBasic       ), pointer       :: basicParent
    class(nodeComponentHotHalo     ), pointer       :: hotHaloParent, hotHalo
    class(nodeComponentSpin        ), pointer       :: spinParent

    if (.not.node%isSatellite()) return
    hotHalo => node%hotHalo(autoCreate=.true.)
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       ! Transfer any outflowed gas to the hot halo of the parent node.
       nodeParent => node%parent
       do while (nodeParent%isSatellite())
          nodeParent => nodeParent%parent
       end do
       basicParent   => nodeParent%basic  (                 )
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       spinParent    => nodeParent%spin   (                 )
       if (hotHaloParent%outflowedAngularMomentumIsSettable()) then
          call hotHaloParent%outflowedAngularMomentumSet(                                           &
               &                                          hotHaloParent %outflowedAngularMomentum() &
               &                                         +hotHalo       %outflowedMass           () &
               &                                         *spinParent    %angularMomentum         () &
               &                                         /basicParent   %mass                    () &
               &                                        )
          call hotHalo      %outflowedAngularMomentumSet(                                           &
               &                                          0.0d0                                     &
               &                                        )
       end if
       if (hotHaloParent%           outflowedMassIsSettable()) then
          call hotHaloParent%outflowedMassSet           (                                           &
               &                                          hotHaloParent %outflowedMass           () &
               &                                         +hotHalo       %outflowedMass           () &
               &                                        )
          call hotHalo      %outflowedMassSet           (                                           &
               &                                          0.0d0                                     &
               &                                        )
       end if
       if (hotHaloParent%     outflowedAbundancesIsSettable()) then
          call hotHaloParent%outflowedAbundancesSet     (                                           &
               &                                          hotHaloParent %outflowedAbundances     () &
               &                                         +hotHalo       %outflowedAbundances     () &
               &                                        )
          call hotHalo      %outflowedAbundancesSet     (                                           &
               &                                          zeroAbundances                            &
               &                                        )
       end if
       if (hotHaloParent%      outflowedChemicalsIsSettable()) then
          call hotHaloParent%outflowedChemicalsSet      (                                           &
               &                                          hotHaloParent %outflowedChemicals      () &
               &                                         +hotHalo       %outflowedChemicals      () &
               &                                        )
          call hotHalo      %outflowedChemicalsSet      (                                           &
               &                                          zeroChemicalAbundances                    &
               &                                        )
       end if
    end select
    return
  end subroutine cgmStarvationDifferentialEvolutionPost
