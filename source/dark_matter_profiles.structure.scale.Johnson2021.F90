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
  Implements a dark matter profile scale radius class using the energy-based model of \cite{johnson_random_2021}.
  !!}
  
  use :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMOClass
  use :: Galacticus_Nodes                  , only : nodeComponentDarkMatterProfile
  use :: Virial_Orbits                     , only : virialOrbitClass
  use :: Merger_Trees_Build_Mass_Resolution, only : mergerTreeMassResolutionClass

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusJohnson2021">
    <description>
      A dark matter profile scale radius class that computes scale radii based on the energy conservation approach of
      \cite{johnson_random_2021}.
    </description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusJohnson2021
     !!{
     A dark matter profile scale radius class that assigns dark matter profile scale radii using the energy-based model of
     \cite{johnson_random_2021}.
     !!}
     private
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_          => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_         => null()
     class           (virialOrbitClass                 ), pointer :: virialOrbit_                  => null()
     class           (mergerTreeMassResolutionClass    ), pointer :: mergerTreeMassResolution_     => null()
     double precision                                             :: massExponent                           , energyBoost, &
          &                                                          unresolvedEnergy
   contains
     final     ::           darkMatterProfileScaleJohnson2021Destructor
     procedure :: radius => darkMatterProfileScaleJohnson2021Radius
  end type darkMatterProfileScaleRadiusJohnson2021
  
  interface darkMatterProfileScaleRadiusJohnson2021
     !!{
     Constructors for the {\normalfont \ttfamily darkMatterProfileScaleJohnson2021} node operator class.
     !!}
     module procedure darkMatterProfileScaleJohnson2021ConstructorParameters
     module procedure darkMatterProfileScaleJohnson2021ConstructorInternal
  end interface darkMatterProfileScaleRadiusJohnson2021

  ! Sub-module-scope variables used in root finding.
  double precision                                                   :: energyTotal
  class           (darkMatterProfileScaleRadiusJohnson2021), pointer :: self_
  class           (nodeComponentDarkMatterProfile         ), pointer :: darkMatterProfile_
  type            (treeNode                               ), pointer :: node_
  !$omp threadprivate(energyTotal,self_,node_,darkMatterProfile_)

contains
  
  function darkMatterProfileScaleJohnson2021ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily randomWalk} dark matter profile scale radius class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileScaleRadiusJohnson2021)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (darkMatterProfileScaleRadiusClass      ), pointer       :: darkMatterProfileScaleRadius_
    class           (darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass              ), pointer       :: darkMatterProfileDMO_
    class           (virialOrbitClass                       ), pointer       :: virialOrbit_
    class           (mergerTreeMassResolutionClass          ), pointer       :: mergerTreeMassResolution_
    double precision                                                         :: massExponent                 , energyBoost, &
         &                                                                      unresolvedEnergy
    
    !![
    <inputParameter>
      <name>energyBoost</name>
      <defaultValue>0.797d0</defaultValue>
      <defaultSource>\citep{johnson_random_2021}</defaultSource>
      <source>parameters</source>
      <description>A boost to the energy.</description>
    </inputParameter>
    <inputParameter>
      <name>massExponent</name>
      <defaultValue>2.168d0</defaultValue>
      <defaultSource>\citep{johnson_random_2021}</defaultSource>
      <source>parameters</source>
      <description>The exponent of mass ratio appearing in the orbital energy term.</description>
    </inputParameter>
    <inputParameter>
      <name>unresolvedEnergy</name>
      <defaultValue>0.550d0</defaultValue>
      <defaultSource>\citep{johnson_random_2021}</defaultSource>
      <source>parameters</source>
      <description>Factor multiplying the estimate of the internal energy of unresolved accretion.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="virialOrbit"                  name="virialOrbit_"                  source="parameters"/>
    <objectBuilder class="mergerTreeMassResolution"     name="mergerTreeMassResolution_"     source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusJohnson2021(massExponent,energyBoost,unresolvedEnergy,darkMatterProfileScaleRadius_,darkMatterHaloScale_,darkMatterProfileDMO_,virialOrbit_,mergerTreeMassResolution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="virialOrbit_"                 />
    <objectDestructor name="mergerTreeMassResolution_"    />
    <objectDestructor name="darkMatterProfileDMO_"        />
    !!]
    return
  end function darkMatterProfileScaleJohnson2021ConstructorParameters

  function darkMatterProfileScaleJohnson2021ConstructorInternal(massExponent,energyBoost,unresolvedEnergy,darkMatterProfileScaleRadius_,darkMatterHaloScale_,darkMatterProfileDMO_,virialOrbit_,mergerTreeMassResolution_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily randomWalk} dark matter profile scale radius class.
    !!}
    implicit none
    type            (darkMatterProfileScaleRadiusJohnson2021)                        :: self
    class           (darkMatterProfileScaleRadiusClass      ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    class           (virialOrbitClass                       ), intent(in   ), target :: virialOrbit_
    class           (mergerTreeMassResolutionClass          ), intent(in   ), target :: mergerTreeMassResolution_
    class           (darkMatterProfileDMOClass              ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                         , intent(in   )         :: massExponent                 , energyBoost, &
         &                                                                              unresolvedEnergy
    !![
    <constructorAssign variables="massExponent, energyBoost, unresolvedEnergy, *darkMatterProfileScaleRadius_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *virialOrbit_, *mergerTreeMassResolution_"/>
    !!]

    return
  end function darkMatterProfileScaleJohnson2021ConstructorInternal

  subroutine darkMatterProfileScaleJohnson2021Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily randomWalk} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusJohnson2021), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    <objectDestructor name="self%darkMatterHaloScale_"         />
    <objectDestructor name="self%virialOrbit_"                 />
    <objectDestructor name="self%mergerTreeMassResolution_"    />
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    !!]
    return
  end subroutine darkMatterProfileScaleJohnson2021Destructor

  double precision function darkMatterProfileScaleJohnson2021Radius(self,node) result(radiusScale)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    use :: Calculations_Resets             , only : Calculations_Reset
    use :: Galacticus_Nodes                , only : nodeComponentBasic             , nodeComponentDarkMatterProfile     , nodeComponentSatellite
    use :: Root_Finder                     , only : rootFinder                     , rangeExpandMultiplicative          , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    use :: Kepler_Orbits                   , only : keplerOrbit
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Beta_Functions                  , only : Beta_Function                  , Beta_Function_Incomplete_Normalized
    use :: Hypergeometric_Functions        , only : Hypergeometric_2F1
    use :: Mass_Distributions              , only : massDistributionClass
    implicit none
    class           (darkMatterProfileScaleRadiusJohnson2021), intent(inout), target    :: self
    type            (treeNode                               ), intent(inout), target    :: node
    type            (treeNode                               )               , pointer   :: nodeChild                                      , nodeSibling                      , &
         &                                                                                 nodeUnresolved
    class           (nodeComponentBasic                     )               , pointer   :: basicChild                                     , basicSibling                     , &
         &                                                                                 basic                                          , basicUnresolved
    class           (nodeComponentDarkMatterProfile         )               , pointer   :: darkMatterProfile                              , darkMatterProfileSibling         , &
         &                                                                                 darkMatterProfileChild                         , darkMatterProfileUnresolved
    class           (nodeComponentSatellite                 )               , pointer   :: satelliteSibling
    class           (massDistributionClass                  )               , pointer   :: massDistribution_
    double precision                                                        , parameter :: massFunctionSlopeLogarithmic            =1.90d0
    double precision                                                        , parameter :: energyInternalFormFactorSlopeLogarithmic=0.03d0
    double precision                                                                    :: energyOrbital                                  , massRatio                        , &
         &                                                                                 radiusScaleChild                               , radiusScaleOriginal              , &
         &                                                                                 massUnresolved                                 , energyInternalSubresolutionFactor, &
         &                                                                                 radiusScaleUnresolved                          , massResolution                   , &
         &                                                                                 energyPotentialSubresolutionFactor             , energyKineticSubresolutionFactor , &
         &                                                                                 a                                              , b                                , &
         &                                                                                 energyPotential                                , energyKinetic                    , &
         &                                                                                 radiusVirial
    type            (rootFinder                             )                           :: finder
    type            (keplerOrbit                            )                           :: orbit
    
    ! Get the resolution of the tree.
    massResolution=self%mergerTreeMassResolution_%resolution(node%hostTree)
    ! Build a root finder to find scale radii.
    finder=rootFinder(                                                             &
         &            rootFunction                 =radiusScaleRoot              , &
         &            toleranceAbsolute            =1.0d-6                       , &
         &            toleranceRelative            =1.0d-3                       , &
         &            rangeExpandDownward          =0.5d+0                       , &
         &            rangeExpandUpward            =2.0d+0                       , &
         &            rangeExpandType              =rangeExpandMultiplicative    , &
         &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
         &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative  &
         &           )
    ! Create the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    ! If this node has no children, set its dark matter profile scale radius using the fallback method.
    if (.not.associated(node%firstChild)) then
       radiusScale=self%darkMatterProfileScaleRadius_%radius(node)
    else
       basic                   =>  node                  %basic            ()
       nodeChild               =>  node                  %firstChild
       darkMatterProfileChild  =>  nodeChild             %darkMatterProfile()
       basicChild              =>  nodeChild             %basic            ()
       radiusScaleChild        =   darkMatterProfileChild%scale            ()
       massUnresolved          =  +basic                 %mass             () &
            &                     -basicChild            %mass             ()       
       ! Iterate over progenitors and sum their energies.
       nodeSibling       =>                                                       nodeChild
       massDistribution_ => self             %darkMatterProfileDMO_ %get         (nodeSibling                   )
       radiusVirial      =  self             %darkMatterHaloScale_  %radiusVirial(nodeSibling                   )
       energyTotal       =  massDistribution_%energy                             (radiusVirial,massDistribution_)
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
       do while (associated(nodeSibling%sibling))
          nodeSibling              =>  nodeSibling       %sibling
          basicSibling             =>  nodeSibling       %basic                             (                                        )
          darkMatterProfileSibling =>  nodeSibling       %darkMatterProfile                 (                                        )
          satelliteSibling         =>  nodeSibling       %satellite                         (autoCreate=.true.                       )
          massDistribution_        =>  self              %darkMatterProfileDMO_%get         (nodeSibling                             )
          radiusVirial             =   self              %darkMatterHaloScale_ %radiusVirial(nodeSibling                             )
          orbit                    =   satelliteSibling  %virialOrbit                       (                                        )
          massRatio                =  +basicSibling      %mass                              (                                        ) &
               &                      /basicChild        %mass                              (                                        )
          energyOrbital            =  +orbit             %energy                            (                                        ) &
               &                      *basicSibling      %mass                              (                                        ) &
               &                      /(                                                                                               &
               &                        +1.0d0                                                                                         &
               &                        +massRatio                                                                                     &
               &                       )**self%massExponent
          massUnresolved           =  +massUnresolved                                                                                  &
               &                      -basicSibling      %mass                              (                                        )
          ! Add orbital energy of this sibling.
          energyTotal              = +energyTotal                                                                                      &
               &                     +energyOrbital                                                                                    &
               &                     *(                                                                                                &
               &                       +1.0d0                                                                                          &
               &                       +self%energyBoost                                                                               &
               &                      )
          ! Add the internal energy of the sibling.
          energyTotal              =  +energyTotal                                                                                     &
               &                      +massDistribution_ %energy                            (radiusVirial,massDistribution_          ) &
               &                      /(                                                                                               &
               &                        +1.0d0                                                                                         &
               &                        +massRatio                                                                                     &
               &                       )**self%massExponent                                                                            &
               &                      *(                                                                                               &
               &                        +1.0d0                                                                                         &
               &                        +self%energyBoost                                                                              &
               &                       )
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]
       end do
       ! Account for unresolved accretion. We assume that unresolved halos are accreted with the mean orbital energy of
       ! the virial orbital parameter distribution, plus an internal energy corresponding to that of a halo with mass
       ! equal to the total unresolved mass scaled by some correction factor (to account for the fact that the unresolved
       ! accretion will not in fact be in a single halo).
       if (massUnresolved > 0.0d0) then
          nodeUnresolved                       => treeNode                                                      (                         )
          nodeUnresolved             %hostTree => node                                        %hostTree
          basicUnresolved                      => nodeUnresolved                              %basic            (autoCreate=.true.        )
          darkMatterProfileUnresolved          => nodeUnresolved                              %darkMatterProfile(autoCreate=.true.        )
          call basicUnresolved            %massSet            (min(massResolution   ,massUnresolved       ))
          call basicUnresolved            %timeSet            (    basicChild%time()                       )
          call basicUnresolved            %timeLastIsolatedSet(    basicChild%time()                       )
          radiusScaleUnresolved                =  self          %darkMatterProfileScaleRadius_%radius           (           nodeUnresolved)
          call darkMatterProfileUnresolved%scaleSet           (                      radiusScaleUnresolved )
          ! Compute a correction factor to the orbital energy which takes into account the mass dependence of the 1/(1+m/M)ᵅ
          ! term that is applied to the orbital energy. Averaging this over a power-law mass function gives the result
          ! below. In the case that α=0 the result is identically 1 - in this case we avoid computing β functions.
          massRatio                       =+basicUnresolved%mass()                                               &
               &                           /basicChild     %mass()
          a                               =+2.0d0                                                                &
               &                           -massFunctionSlopeLogarithmic
          b                               =+massFunctionSlopeLogarithmic                                         &
               &                           +self%massExponent                                                    &
               &                           -1.0d0
          energyKineticSubresolutionFactor=+Beta_Function_Incomplete_Normalized(a,b,massRatio/(1.0d0+massRatio)) &
               &                           *Beta_Function                      (a,b                            ) &
               &                           *           (2.0d0-massFunctionSlopeLogarithmic)                      &
               &                           /massRatio**(2.0d0-massFunctionSlopeLogarithmic)
          if (self%massExponent == 0.0d0) then
             energyPotentialSubresolutionFactor=+1.0d0
          else
             b                                 =+massFunctionSlopeLogarithmic                                         &
                  &                             +self%massExponent                                                    &
                  &                             -2.0d0
             energyPotentialSubresolutionFactor=+Beta_Function_Incomplete_Normalized(a,b,massRatio/(1.0d0+massRatio)) &
                  &                             *Beta_Function                      (a,b                            ) &
                  &                             *           (2.0d0-massFunctionSlopeLogarithmic)                      &
                  &                             /massRatio**(2.0d0-massFunctionSlopeLogarithmic)
          end if
          ! Compute a correction factor to the internal energy which takes into account the mass dependence of the 1/(1+m/M)ᵅ
          ! term that is applied to the orbital energy. Averaging this over a power-law mass function gives the result
          ! below.
          energyInternalSubresolutionFactor=+(                                                                                                                            &
               &                              +2.0d0                                                                                                                      &
               &                              -(massFunctionSlopeLogarithmic+energyInternalFormFactorSlopeLogarithmic)                                                    &
               &                             )                                                                                                                            &
               &                            *Hypergeometric_2F1(                                                                                                          &
               &                                                [self%massExponent, 8.0d0/3.0d0-(massFunctionSlopeLogarithmic+energyInternalFormFactorSlopeLogarithmic)], &
               &                                                [                  11.0d0/3.0d0-(massFunctionSlopeLogarithmic+energyInternalFormFactorSlopeLogarithmic)], &
               &                                                -1.0d0/massRatio                                                                                          &
               &                                               )                                                                                                          &
               &                            /(                                                                                                                            &
               &                                                                    8.0d0/3.0d0-(massFunctionSlopeLogarithmic+energyInternalFormFactorSlopeLogarithmic)   &
               &                              )
          ! Determine the orbital and internal energies.
          massDistribution_ => self%darkMatterProfileDMO_%get         (nodeUnresolved)
          radiusVirial      =  self%darkMatterHaloScale_ %radiusVirial(nodeUnresolved)
          energyKinetic     =  +0.5d0                                                                       &
               &               *self%virialOrbit_%velocityTotalRootMeanSquared(nodeUnresolved,nodeChild)**2 &
               &               /(1.0d0+massRatio)
          energyPotential   =  +self%virialOrbit_%energyMean                  (nodeUnresolved,nodeChild)    &
               &               -                  energyKinetic
          energyOrbital     =  +energyPotential                                            &
               &               *energyPotentialSubresolutionFactor                         &
               &               +energyKinetic                                              &
               &               *energyKineticSubresolutionFactor
          energyTotal       =  +energyTotal                                                &
               &               +massUnresolved                                             &
               &               *self%unresolvedEnergy                                      &
               &               *(                                                          &
               &                 +energyOrbital                                            &
               &                 +energyInternalSubresolutionFactor                        &
               &                 *massDistribution_%energy(radiusVirial,massDistribution_) &
               &                 /massResolution                                           &
               &                )                                                          &
               &               *(                                                          &
               &                 +1.0d0                                                    &
               &                 +self%energyBoost                                         &
               &                )
          call nodeUnresolved%destroy()
          deallocate(nodeUnresolved)
          !![
          <objectDestructor name="massDistribution_"/>
	  !!]
       end if
       ! Add mutual gravitational binding energy of any sibling halo and any unresolved mass.
       if (associated(nodeChild%sibling))                                      &
            & energyTotal=+energyTotal                                         &
            &             -gravitationalConstant_internal                      &
            &             *basicSibling             %mass        (         )   &
            &             *massUnresolved                                      &
            &             /0.5d0                                               &
            &             /self%darkMatterHaloScale_%radiusVirial(nodeChild)   &
            &             *self%unresolvedEnergy                               &
            &             *(                                                   &
            &               +1.0d0                                             &
            &               +self%energyBoost                                  &
            &              )
       ! Convert energy back to scale radius.
       self_               => self
       node_               => node
       darkMatterProfile_  => darkMatterProfile
       radiusScaleOriginal =  darkMatterProfile%scale(                          )
       radiusScale         =  finder           %find (rootGuess=radiusScaleChild)
       call darkMatterProfile%scaleSet(radiusScaleOriginal)       
       call Calculations_Reset(node)
    end if
    return
  end function darkMatterProfileScaleJohnson2021Radius

  double precision function radiusScaleRoot(radiusScale)
    !!{
    Function used in root-finding to compute the scale radius of a dark matter profile as a given energy.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Mass_Distributions , only : massDistributionClass
    implicit none
    double precision                       , intent(in   ) :: radiusScale
    class           (massDistributionClass), pointer       :: massDistribution_
    
    call darkMatterProfile_%scaleSet(radiusScale)
    call Calculations_Reset(node_)
    massDistribution_ =>  self_%darkMatterProfileDMO_%get        (                                        node_                   )
    radiusScaleRoot    =  +                           energyTotal                                                                   &
         &                -     massDistribution_    %energy     (self_%darkMatterHaloScale_%radiusVirial(node_),massDistribution_)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function radiusScaleRoot
