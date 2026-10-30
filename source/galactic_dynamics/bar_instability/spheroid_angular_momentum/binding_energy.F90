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

!+    Contributions to this file made by: Andrew Robertson, Codex, Claude.

!!{RST
Implements a binding-energy-conserving pseudo-angular momentum prescription for spheroids formed by bar instabilities.
!!}

  !![
  <barInstabilitySpheroidAngularMomentum name="barInstabilitySpheroidAngularMomentumBindingEnergy" docformat="rst">
   <description>
   Assigns pseudo-angular momentum to a spheroid by blending the angular-momentum-transfer prescription for a nascent spheroid into a binding-energy-conserving prescription for an established spheroid. The binding-energy proxy is

   .. math::

      \mathcal{B} = c_\mathrm{d}\frac{G M_\mathrm{d}^2}{R_\mathrm{d}}
      + c_\mathrm{s}\frac{G M_\mathrm{s}^2}{R_\mathrm{s}}
      + f_\mathrm{int}\frac{G M_\mathrm{d}M_\mathrm{s}}{R_\mathrm{d}+R_\mathrm{s}},

   where radii are baryonic half-mass radii. Total derivatives of this proxy include the response of the equilibrium structure solver to changes in component masses and pseudo-angular momenta. The interpolation coordinate is the baryonic spheroid fraction

   .. math::

      \mu_\mathrm{s}=\frac{M_\mathrm{s}}{M_\mathrm{d}+M_\mathrm{s}}.

   The angular-momentum-transfer prescription is used below ``fractionBaryonicSpheroidTransitionMinimum`` and the binding-energy
   prescription above ``fractionBaryonicSpheroidTransitionMaximum``. Between the two limits, their rates are combined using a
   cubic smoothstep. The default limits are 0.003 and 0.010. Equal limits give a sharp transition; explicitly setting both limits
   to zero recovers the unregularized energy prescription for non-zero spheroids where the energy derivative is well defined. This
   is not recommended for production runs because the derivative becomes poorly conditioned for nascent spheroids. The derivatives
   are only as accurate as the equilibrium radii: at the default settings the prescription does not conserve the proxy, so
   ``stepFiniteDifferenceRelative`` and the structure solver's ``solutionTolerance`` should both be reduced, for example to
   :math:`10^{-3.5}` and :math:`10^{-7}`. See :doc:`/manuals/user-guide/bar-instability-spheroid-sizes` for the method and its
   limitations.
   </description>
  </barInstabilitySpheroidAngularMomentum>
  !!]
  type, extends(barInstabilitySpheroidAngularMomentumClass) :: barInstabilitySpheroidAngularMomentumBindingEnergy
     !!{RST
     Binding-energy-conserving pseudo-angular momentum prescription for a bar-built spheroid.
     !!}
     private
     class           (galacticStructureSolverClass), pointer :: galacticStructureSolver_                  => null()
     double precision                                        :: stepFiniteDifferenceRelative                       , formFactorDisk                           , &
          &                                                     formFactorSpheroid                                 , formFactorInteraction                    , &
          &                                                     fractionBaryonicSpheroidTransitionMaximum          , fractionBaryonicSpheroidTransitionMinimum
   contains
     final     ::         bindingEnergyDestructor
     procedure :: rate => bindingEnergyRate
  end type barInstabilitySpheroidAngularMomentumBindingEnergy

  interface barInstabilitySpheroidAngularMomentumBindingEnergy
     module procedure bindingEnergyConstructorParameters
     module procedure bindingEnergyConstructorInternal
  end interface barInstabilitySpheroidAngularMomentumBindingEnergy

contains

  function bindingEnergyConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the binding-energy prescription from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (barInstabilitySpheroidAngularMomentumBindingEnergy)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    class           (galacticStructureSolverClass                      ), pointer       :: galacticStructureSolver_
    double precision                                                                    :: stepFiniteDifferenceRelative             , &
         &                                                                                 formFactorDisk                           , &
         &                                                                                 formFactorSpheroid                       , &
         &                                                                                 formFactorInteraction                    , &
         &                                                                                 fractionBaryonicSpheroidTransitionMaximum, &
         &                                                                                 fractionBaryonicSpheroidTransitionMinimum

    !![
    <inputParameter docformat="rst">
      <name>stepFiniteDifferenceRelative</name>
      <defaultValue>1.0d-2</defaultValue>
      <description>
      The fractional step used to measure the response of the binding energy to bar-driven changes in masses and pseudo-angular momenta.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>formFactorDisk</name>
      <defaultValue>5.0d-1</defaultValue>
      <description>
      The coefficient :math:`c_\mathrm{d}` multiplying the disk self-binding term.
      </description>
      <source>parameters</source>
      <minimum>0.0</minimum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>formFactorSpheroid</name>
      <defaultValue>5.0d-1</defaultValue>
      <description>
      The coefficient :math:`c_\mathrm{s}` multiplying the spheroid self-binding term.
      </description>
      <source>parameters</source>
      <minimum>0.0</minimum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>formFactorInteraction</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The coefficient :math:`f_\mathrm{int}` multiplying the disk--spheroid interaction term.
      </description>
      <source>parameters</source>
      <minimum>0.0</minimum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>fractionBaryonicSpheroidTransitionMinimum</name>
      <defaultValue>3.0d-3</defaultValue>
      <description>
      The baryonic spheroid fraction below which the angular-momentum-transfer prescription is used. The baryonic spheroid fraction includes both gas and stars in the disk and spheroid.
      </description>
      <source>parameters</source>
      <minimum>0.0</minimum>
      <maximum>1.0</maximum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>fractionBaryonicSpheroidTransitionMaximum</name>
      <defaultValue>1.0d-2</defaultValue>
      <description>
      The baryonic spheroid fraction above which the binding-energy prescription is used. Between the minimum and maximum transition fractions, the two angular-momentum rates are blended with a cubic smoothstep. Setting this equal to the minimum gives a sharp transition.
      </description>
      <source>parameters</source>
      <minimum>0.0</minimum>
      <maximum>1.0</maximum>
    </inputParameter>
    <objectBuilder class="galacticStructureSolver" name="galacticStructureSolver_" source="parameters"/>
    !!]
    self=barInstabilitySpheroidAngularMomentumBindingEnergy(                                           &
         &                                                  stepFiniteDifferenceRelative             , &
         &                                                  formFactorDisk                           , &
         &                                                  formFactorSpheroid                       , &
         &                                                  formFactorInteraction                    , &
         &                                                  fractionBaryonicSpheroidTransitionMinimum, &
         &                                                  fractionBaryonicSpheroidTransitionMaximum, &
         &                                                  galacticStructureSolver_                   &
         &                                                 )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructureSolver_"/>
    !!]
    return
  end function bindingEnergyConstructorParameters

  function bindingEnergyConstructorInternal(stepFiniteDifferenceRelative,formFactorDisk,formFactorSpheroid,formFactorInteraction,fractionBaryonicSpheroidTransitionMinimum,fractionBaryonicSpheroidTransitionMaximum,galacticStructureSolver_) result(self)
    !!{RST
    Internal constructor for the binding-energy prescription.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galactic_Structure_Solvers, only : galacticStructureSolverFixed, galacticStructureSolverNull
    implicit none
    type            (barInstabilitySpheroidAngularMomentumBindingEnergy)                        :: self
    double precision                                                    , intent(in   )         :: stepFiniteDifferenceRelative             , &
         &                                                                                         formFactorDisk                           , &
         &                                                                                         formFactorSpheroid                       , &
         &                                                                                         formFactorInteraction                    , &
         &                                                                                         fractionBaryonicSpheroidTransitionMaximum, &
         &                                                                                         fractionBaryonicSpheroidTransitionMinimum
    class           (galacticStructureSolverClass                      ), intent(in   ), target :: galacticStructureSolver_
    !![
    <constructorAssign variables="stepFiniteDifferenceRelative, formFactorDisk, formFactorSpheroid, formFactorInteraction, fractionBaryonicSpheroidTransitionMinimum, fractionBaryonicSpheroidTransitionMaximum, *galacticStructureSolver_"/>
    !!]

    if (self%fractionBaryonicSpheroidTransitionMaximum < self%fractionBaryonicSpheroidTransitionMinimum) &
         & call Error_Report('the maximum baryonic spheroid-fraction transition must be greater than or equal to the minimum'//{introspection:location})
    ! Reject configurations for which the binding-energy derivative vanishes or is ill-conditioned.
    if (self%stepFiniteDifferenceRelative <= 0.0d0) &
         & call Error_Report('[stepFiniteDifferenceRelative] must be positive'//{introspection:location})
    if (self%formFactorSpheroid <= 0.0d0 .and. self%formFactorInteraction <= 0.0d0) &
         & call Error_Report('at least one of [formFactorSpheroid] and [formFactorInteraction] must be positive - otherwise the binding energy depends on the spheroid angular momentum only weakly, through the disk radius, and its derivative is ill-conditioned'//{introspection:location})
    select type (galacticStructureSolver_)
    class is (galacticStructureSolverFixed)
       call Error_Report('the binding-energy prescription requires a galactic structure solver whose radii depend on angular momentum - the "fixed" solver is not supported'//{introspection:location})
    class is (galacticStructureSolverNull )
       call Error_Report('the binding-energy prescription requires a galactic structure solver whose radii depend on angular momentum - the "null" solver is not supported' //{introspection:location})
    end select
    !![
    <componentPropertyAssert class="disk"     properties="massGas massStellar angularMomentum" require="gettable settable"/>
    <componentPropertyAssert class="spheroid" properties="massGas massStellar angularMomentum" require="gettable settable"/>
    <componentPropertyAssert class="disk"     properties="halfMassRadius"                      require="gettable"         />
    <componentPropertyAssert class="spheroid" properties="halfMassRadius"                      require="gettable"         />
    !!]
    return
  end function bindingEnergyConstructorInternal

  subroutine bindingEnergyDestructor(self)
    !!{RST
    Destructor for the binding-energy prescription.
    !!}
    implicit none
    type(barInstabilitySpheroidAngularMomentumBindingEnergy), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticStructureSolver_"/>
    !!]
    return
  end subroutine bindingEnergyDestructor

  double precision function bindingEnergyRate(self,node,rateMassGasTransfer,rateMassStellarTransfer,rateAngularMomentumDisk,rateAngularMomentumTransfer,fractionAngularMomentumRetainedSpheroid)
    !!{RST
    Return the spheroid pseudo-angular momentum rate, blending from the angular-momentum-transfer seed to the rate that conserves the binding-energy proxy along the bar-driven direction in state space.
    !!}
    use, intrinsic :: IEEE_Arithmetic    , only : ieee_is_finite
    use            :: Calculations_Resets, only : Calculations_Reset
    use            :: Galacticus_Nodes   , only : nodeComponentDisk , nodeComponentSpheroid
    implicit none
    class           (barInstabilitySpheroidAngularMomentumBindingEnergy), intent(inout)         :: self
    type            (treeNode                                          ), intent(inout), target :: node
    double precision                                                    , intent(in   )         :: rateMassGasTransfer                           , &
         &                                                                                         rateMassStellarTransfer                       , &
         &                                                                                         rateAngularMomentumDisk                       , &
         &                                                                                         rateAngularMomentumTransfer                   , &
         &                                                                                         fractionAngularMomentumRetainedSpheroid
    class           (nodeComponentDisk                                 ), pointer               :: disk
    class           (nodeComponentSpheroid                             ), pointer               :: spheroid
    logical                                                                                     :: perturbationDefined
    double precision                                                                            :: angularMomentumDisk                           , &
         &                                                                                         angularMomentumSpheroid                       , &
         &                                                                                         baryonMassDisk                                , &
         &                                                                                         baryonMassSpheroid                            , &
         &                                                                                         baryonMassTotal                               , &
         &                                                                                         bindingEnergyInitial                          , &
         &                                                                                         bindingEnergyKnownDirection                   , &
         &                                                                                         bindingEnergySpheroidPerturbed                , &
         &                                                                                         derivativeBindingEnergyKnown                  , &
         &                                                                                         derivativeBindingEnergySpheroidAngularMomentum, &
         &                                                                                         massGasDisk                                   , &
         &                                                                                         massGasSpheroid                               , &
         &                                                                                         massStellarDisk                               , &
         &                                                                                         massStellarSpheroid                           , &
         &                                                                                         perturbationAngularMomentumSpheroid           , &
         &                                                                                         perturbationTime                              , &
         &                                                                                         rateBindingEnergy                             , &
         &                                                                                         rateRetained                                  , &
         &                                                                                         spheroidBaryonFraction                        , &
         &                                                                                         timescaleMinimum                              , &
         &                                                                                         transitionCoordinate                          , &
         &                                                                                         transitionWeight

    disk     => node%disk    ()
    spheroid => node%spheroid()
    angularMomentumDisk      =disk    %angularMomentum()
    angularMomentumSpheroid  =spheroid%angularMomentum()
    massGasDisk              =disk    %massGas        ()
    massGasSpheroid          =spheroid%massGas        ()
    massStellarDisk          =disk    %massStellar    ()
    massStellarSpheroid      =spheroid%massStellar    ()
    baryonMassDisk           =massGasDisk    +massStellarDisk
    baryonMassSpheroid       =massGasSpheroid+massStellarSpheroid
    baryonMassTotal          =baryonMassDisk+baryonMassSpheroid
    rateRetained             =fractionAngularMomentumRetainedSpheroid*rateAngularMomentumTransfer

    ! Seed nascent spheroids with the standard prescription. The binding-energy derivative is ill-conditioned while their
    ! self-binding term is negligible, and is not defined for a zero-mass or zero-angular-momentum spheroid.
    if     (                                      &
         &       baryonMassSpheroid      <= 0.0d0 &
         &  .or. baryonMassTotal         <= 0.0d0 &
         &  .or. angularMomentumSpheroid <= 0.0d0 &
         & ) then
       bindingEnergyRate=rateRetained
       return
    end if
    spheroidBaryonFraction=baryonMassSpheroid/baryonMassTotal
    if (spheroidBaryonFraction <= self%fractionBaryonicSpheroidTransitionMinimum) then
       bindingEnergyRate=rateRetained
       return
    end if

    bindingEnergyInitial=bindingEnergyProxy(self,node)
    if (.not.ieee_is_finite(bindingEnergyInitial) .or. bindingEnergyInitial <= 0.0d0) then
       bindingEnergyRate=rateRetained
       return
    end if

    ! Choose a perturbation corresponding to the requested fractional change in the most rapidly changing property.
    timescaleMinimum  =huge(timescaleMinimum)
    perturbationDefined=.false.
    if (rateMassGasTransfer > 0.0d0 .and. massGasDisk > 0.0d0) then
       timescaleMinimum  =min(timescaleMinimum,massGasDisk/rateMassGasTransfer)
       perturbationDefined=.true.
    end if
    if (rateMassStellarTransfer > 0.0d0 .and. massStellarDisk > 0.0d0) then
       timescaleMinimum  =min(timescaleMinimum,massStellarDisk/rateMassStellarTransfer)
       perturbationDefined=.true.
    end if
    if (abs(rateAngularMomentumDisk) > 0.0d0 .and. angularMomentumDisk > 0.0d0) then
       timescaleMinimum  =min(timescaleMinimum,angularMomentumDisk/abs(rateAngularMomentumDisk))
       perturbationDefined=.true.
    end if
    if (.not.perturbationDefined .or. .not.ieee_is_finite(timescaleMinimum)) then
       bindingEnergyRate=rateRetained
       return
    end if
    perturbationTime=self%stepFiniteDifferenceRelative*timescaleMinimum
    call disk    %massGasSet        (massGasDisk        -rateMassGasTransfer    *perturbationTime)
    call disk    %massStellarSet    (massStellarDisk    -rateMassStellarTransfer*perturbationTime)
    call disk    %angularMomentumSet(angularMomentumDisk+rateAngularMomentumDisk*perturbationTime)
    call spheroid%massGasSet        (massGasSpheroid    +rateMassGasTransfer    *perturbationTime)
    call spheroid%massStellarSet    (massStellarSpheroid+rateMassStellarTransfer*perturbationTime)
    call Calculations_Reset(node)
    call self%galacticStructureSolver_%revert(node)
    call self%galacticStructureSolver_%solve (node)
    bindingEnergyKnownDirection=bindingEnergyProxy(self,node)
    if (.not.ieee_is_finite(bindingEnergyKnownDirection) .or. bindingEnergyKnownDirection <= 0.0d0) then
       call restoreNode()
       bindingEnergyRate=rateRetained
       return
    end if

    ! Move directly to the spheroid-angular-momentum perturbation. The solver's revert state still holds the original radii, so
    ! restoring equilibrium between the two probes would perform a redundant structure solve.
    call restoreProperties()
    perturbationAngularMomentumSpheroid=self%stepFiniteDifferenceRelative*angularMomentumSpheroid
    call spheroid%angularMomentumSet(angularMomentumSpheroid+perturbationAngularMomentumSpheroid)
    call Calculations_Reset(node)
    call self%galacticStructureSolver_%revert(node)
    call self%galacticStructureSolver_%solve (node)
    bindingEnergySpheroidPerturbed=bindingEnergyProxy(self,node)
    call restoreNode()
    if (.not.ieee_is_finite(bindingEnergySpheroidPerturbed) .or. bindingEnergySpheroidPerturbed <= 0.0d0) then
       bindingEnergyRate=rateRetained
       return
    end if
    derivativeBindingEnergyKnown=(bindingEnergyKnownDirection-bindingEnergyInitial)/perturbationTime
    derivativeBindingEnergySpheroidAngularMomentum=(bindingEnergySpheroidPerturbed-bindingEnergyInitial)/perturbationAngularMomentumSpheroid
    if     (                                                                                          &
         &   .not.ieee_is_finite(derivativeBindingEnergyKnown)                                        &
         &  .or. .not.ieee_is_finite(derivativeBindingEnergySpheroidAngularMomentum)                  &
         &  .or. abs(derivativeBindingEnergySpheroidAngularMomentum) <= epsilon(bindingEnergyInitial) &
         &       *max(abs(bindingEnergyInitial),1.0d0)/max(abs(angularMomentumSpheroid),1.0d0)        &
         & ) then
       bindingEnergyRate=rateRetained
    else
       rateBindingEnergy=-derivativeBindingEnergyKnown/derivativeBindingEnergySpheroidAngularMomentum
       if (.not.ieee_is_finite(rateBindingEnergy)) then
          bindingEnergyRate=rateRetained
       else
          if (self%fractionBaryonicSpheroidTransitionMaximum <= self%fractionBaryonicSpheroidTransitionMinimum) then
             transitionWeight=1.0d0
          else
             transitionCoordinate=(spheroidBaryonFraction-self%fractionBaryonicSpheroidTransitionMinimum) &
                  &               /(self%fractionBaryonicSpheroidTransitionMaximum-self%fractionBaryonicSpheroidTransitionMinimum)
             transitionCoordinate=max(0.0d0,min(1.0d0,transitionCoordinate))
             transitionWeight=transitionCoordinate**2*(3.0d0-2.0d0*transitionCoordinate)
          end if
          bindingEnergyRate=(1.0d0-transitionWeight)*rateRetained+transitionWeight*rateBindingEnergy
       end if
    end if
    return

  contains

    subroutine restoreNode()
      !!{RST
      Restore the unperturbed component state and its corresponding equilibrium structure.
      !!}
      implicit none

      call restoreProperties()
      call Calculations_Reset(node)
      call self%galacticStructureSolver_%revert(node)
      call self%galacticStructureSolver_%solve (node)
      return
    end subroutine restoreNode

    subroutine restoreProperties()
      !!{RST
      Restore the unperturbed masses and pseudo-angular momenta without solving for structure.
      !!}
      implicit none

      call disk    %massGasSet        (massGasDisk            )
      call disk    %massStellarSet    (massStellarDisk        )
      call disk    %angularMomentumSet(angularMomentumDisk    )
      call spheroid%massGasSet        (massGasSpheroid        )
      call spheroid%massStellarSet    (massStellarSpheroid    )
      call spheroid%angularMomentumSet(angularMomentumSpheroid)
      return
    end subroutine restoreProperties

  end function bindingEnergyRate

  double precision function bindingEnergyProxy(self,node)
    !!{RST
    Return the magnitude of the disk--spheroid binding-energy proxy.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDisk             , nodeComponentSpheroid
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (barInstabilitySpheroidAngularMomentumBindingEnergy), intent(inout) :: self
    type            (treeNode                                          ), intent(inout) :: node
    class           (nodeComponentDisk                                 ), pointer       :: disk
    class           (nodeComponentSpheroid                             ), pointer       :: spheroid
    double precision                                                                    :: massDisk  , massSpheroid  , &
         &                                                                                 radiusDisk, radiusSpheroid

    disk           => node    %disk          ()
    spheroid       => node    %spheroid      ()
    massDisk       =  disk    %massGas       ()+disk    %massStellar()
    massSpheroid   =  spheroid%massGas       ()+spheroid%massStellar()
    radiusDisk     =  disk    %halfMassRadius()
    radiusSpheroid =  spheroid%halfMassRadius()
    if (massDisk < 0.0d0 .or. massSpheroid < 0.0d0 .or. radiusDisk <= 0.0d0 .or. radiusSpheroid <= 0.0d0) then
       bindingEnergyProxy=-1.0d0
    else
       bindingEnergyProxy=gravitationalConstant_internal*(                                                &
            &                                             + self%formFactorDisk                           &
            &                                             *      massDisk          **2                    &
            &                                             /      radiusDisk                               &
            &                                             + self%formFactorSpheroid                       &
            &                                                                         *massSpheroid  **2  &
            &                                                                         /radiusSpheroid     &
            &                                             + self%formFactorInteraction                    &
            &                                             *      massDisk                                 &
            &                                             *                            massSpheroid       &
            &                                             /(     radiusDisk           +radiusSpheroid   ) &
            &                                            )
    end if
    return
  end function bindingEnergyProxy
