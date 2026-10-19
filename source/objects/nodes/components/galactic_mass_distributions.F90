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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module which builds the mass distributions of the galactic (disk, spheroid, and nuclear star cluster) node components.
!!}

module Node_Components_Galactic_Mass_Distributions
  !!{RST
  Builds the mass distributions of the galactic (disk, spheroid, and nuclear star cluster) node components, by scaling their
  dimensionless stellar and gas mass distributions (optionally truncated, in the case of spherical distributions) to the size
  and masses of a particular component. The procedure for untruncated distributions is written once as a template and
  instantiated for cylindrical and spherical distributions.
  !!}
  implicit none
  private
  public :: Node_Component_Mass_Distribution_Scaled_Cylindrical, Node_Component_Mass_Distribution_Scaled_Spherical, Node_Component_Mass_Distribution_Truncated_Spherical

  !![
  <generic identifier="geometry">
   <instance label="Cylindrical" intrinsic="massDistributionCylindrical" scaler="massDistributionCylindricalScaler"/>
   <instance label="Spherical"   intrinsic="massDistributionSpherical"   scaler="massDistributionSphericalScaler"  />
  </generic>
  !!]

contains

  function Node_Component_Mass_Distribution_Scaled_{geometry¦label}(includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,scalerStellarPool,scalerGasPool) result(massDistribution_)
    !!{RST
    Return the mass distribution of a galactic component with the given stellar and gas masses and scale radius, built by
    scaling the component's dimensionless stellar and gas mass distributions. Scalers are drawn from the given per-thread pools,
    re-initialized if re-used or constructed if not. A null distribution is returned if the scale radius is non-positive, or
    if neither stars nor gas are to be included.
    !!}
    use :: Mass_Distributions, only : massDistributionClass, {geometry¦intrinsic}, {geometry¦scaler}, kinematicsDistributionLocal
    use :: Object_Pools      , only : objectPool
    implicit none
    class           (massDistributionClass      ), pointer                :: massDistribution_
    logical                                      , intent(in   )          :: includeStars              , includeGas
    double precision                             , intent(in   )          :: massStellar               , massGas             , &
         &                                                                   radiusScale
    class           (massDistributionClass      ), intent(in   ), pointer :: massDistributionStellar_  , massDistributionGas_
    type            (kinematicsDistributionLocal), intent(in   ), pointer :: kinematicDistribution_
    type            (objectPool                 ), intent(inout)          :: scalerStellarPool         , scalerGasPool
    class           (massDistributionClass      ), pointer                :: massDistributionStellar   , massDistributionGas
    logical                                                               :: reusedStellar             , reusedGas
    integer                                                               :: iStellar                  , iGas

    if (radiusScale <= 0.0d0 .or. .not.(includeGas .or. includeStars)) then
       ! The component has non-positive size, or no components matched. Return a null distribution.
       massDistribution_ => null()
    else
       ! Build the individual distributions.
       massDistributionStellar => null()
       massDistributionGas     => null()
       if (includeStars) then
          ! Acquire a stellar scaler from the per-thread pool: re-initialize its scaling factors on
          ! re-use, or construct a new one (wrapping the dimensionless stellar distribution) on a miss.
          call scalerStellarPool%acquire(iStellar,reusedStellar)
          if (.not.reusedStellar) allocate({geometry¦scaler} :: scalerStellarPool%slots(iStellar)%object_)
          select type (scaler_ => scalerStellarPool%slots(iStellar)%object_)
          type is ({geometry¦scaler})
             if (reusedStellar) then
                call scaler_%initialize(factorScalingLength=radiusScale,factorScalingMass=massStellar)
             else
                select type (massDistributionStellar_)
                class is ({geometry¦intrinsic})
                   !![
                   <referenceConstruct object="scaler_" constructor="{geometry¦scaler}(factorScalingLength=radiusScale,factorScalingMass=massStellar,massDistribution_=massDistributionStellar_)"/>
                   !!]
                end select
             end if
             call scaler_%setKinematicsDistribution(kinematicDistribution_)
             ! Add the local-handle reference (the pool retains its own); from here the function manages
             ! this handle exactly as it would a freshly-constructed scaler.
             !![
             <referenceCountIncrement object="scaler_"/>
             !!]
             massDistributionStellar => scaler_
          end select
       end if
       if (includeGas  ) then
          ! Acquire a gas scaler from the per-thread pool (see the stellar case above).
          call scalerGasPool%acquire(iGas,reusedGas)
          if (.not.reusedGas) allocate({geometry¦scaler} :: scalerGasPool%slots(iGas)%object_)
          select type (scaler_ => scalerGasPool%slots(iGas)%object_)
          type is ({geometry¦scaler})
             if (reusedGas) then
                call scaler_%initialize(factorScalingLength=radiusScale,factorScalingMass=massGas)
             else
                select type (massDistributionGas_)
                class is ({geometry¦intrinsic})
                   !![
                   <referenceConstruct object="scaler_" constructor="{geometry¦scaler}(factorScalingLength=radiusScale,factorScalingMass=massGas,massDistribution_=massDistributionGas_)"/>
                   !!]
                end select
             end if
             call scaler_%setKinematicsDistribution(kinematicDistribution_)
             !![
             <referenceCountIncrement object="scaler_"/>
             !!]
             massDistributionGas => scaler_
          end select
       end if
       ! Combine the distributions as necessary.
       massDistribution_ => Node_Component_Mass_Distribution_Combine(includeStars,includeGas,massDistributionStellar,massDistributionGas)
    end if
    return
  end function Node_Component_Mass_Distribution_Scaled_{geometry¦label}

  function Node_Component_Mass_Distribution_Truncated_Spherical(includeStars,includeGas,massStellar,massGas,radiusScale,massDistributionStellar_,massDistributionGas_,kinematicDistribution_,galacticComponentTruncation_,fractionMassRetained) result(massDistribution_)
    !!{RST
    Return the mass distribution of a spherical galactic component with the given stellar and gas masses and scale radius, built
    by truncating the component's dimensionless stellar and gas mass distributions (using ``galacticComponentTruncation_``) such
    that they retain a fraction ``fractionMassRetained`` of their mass, and then scaling them to the size and masses of the
    component. Since the truncated distributions differ from node to node they are constructed afresh, rather than being drawn
    from a pool. A null distribution is returned if the scale radius is non-positive, or if neither stars nor gas are to be
    included.
    !!}
    use :: Galactic_Component_Truncations, only : galacticComponentTruncationClass
    use :: Mass_Distributions            , only : massDistributionClass           , kinematicsDistributionLocal
    implicit none
    class           (massDistributionClass           ), pointer                :: massDistribution_
    logical                                           , intent(in   )          :: includeStars                , includeGas
    double precision                                  , intent(in   )          :: massStellar                 , massGas             , &
         &                                                                        radiusScale                 , fractionMassRetained
    class           (massDistributionClass           ), intent(in   ), pointer :: massDistributionStellar_    , massDistributionGas_
    type            (kinematicsDistributionLocal     ), intent(in   ), pointer :: kinematicDistribution_
    class           (galacticComponentTruncationClass), intent(inout)          :: galacticComponentTruncation_
    class           (massDistributionClass           ), pointer                :: massDistributionStellar     , massDistributionGas

    if (radiusScale <= 0.0d0 .or. .not.(includeGas .or. includeStars)) then
       ! The component has non-positive size, or no components matched. Return a null distribution.
       massDistribution_ => null()
    else
       ! Build the individual distributions.
       massDistributionStellar => null()
       massDistributionGas     => null()
       if (includeStars) massDistributionStellar => scaledTruncated(massDistributionStellar_,massStellar)
       if (includeGas  ) massDistributionGas     => scaledTruncated(massDistributionGas_    ,massGas    )
       ! Combine the distributions as necessary.
       massDistribution_ => Node_Component_Mass_Distribution_Combine(includeStars,includeGas,massDistributionStellar,massDistributionGas)
    end if
    return

  contains

    function scaledTruncated(massDistributionDimensionless_,mass) result(massDistributionScaled)
      !!{RST
      Truncate the given dimensionless mass distribution, and scale it to the size and the given mass of the component.
      !!}
      use :: Error             , only : Error_Report
      use :: Mass_Distributions, only : massDistributionSpherical, massDistributionSphericalScaler
      implicit none
      class           (massDistributionClass          ), pointer                :: massDistributionScaled
      class           (massDistributionClass          ), intent(in   ), pointer :: massDistributionDimensionless_
      double precision                                 , intent(in   )          :: mass
      class           (massDistributionSpherical      ), pointer                :: massDistributionTruncated
      type            (massDistributionSphericalScaler), pointer                :: massDistributionScaler

      select type (massDistributionDimensionless_)
      class is (massDistributionSpherical)
         massDistributionTruncated => galacticComponentTruncation_%massDistribution(massDistributionDimensionless_,fractionMassRetained)
      class default
         massDistributionTruncated => null()
         call Error_Report('truncation requires a spherically-symmetric mass distribution'//{introspection:location})
      end select
      ! Scale the truncated distribution such that its total mass is that of the component.
      allocate(massDistributionScaler)
      !![
      <referenceConstruct object="massDistributionScaler" constructor="massDistributionSphericalScaler(factorScalingLength=radiusScale,factorScalingMass=mass/massDistributionTruncated%massTotal(),massDistribution_=massDistributionTruncated)"/>
      <objectDestructor name="massDistributionTruncated"/>
      !!]
      call massDistributionScaler%setKinematicsDistribution(kinematicDistribution_)
      massDistributionScaled => massDistributionScaler
      return
    end function scaledTruncated

  end function Node_Component_Mass_Distribution_Truncated_Spherical

  function Node_Component_Mass_Distribution_Combine(includeStars,includeGas,massDistributionStellar,massDistributionGas) result(massDistribution_)
    !!{RST
    Combine the stellar and gas mass distributions of a galactic component, returning a composite of the two if both are
    included, or the one which is included otherwise. Ownership of the references held by ``massDistributionStellar`` and
    ``massDistributionGas`` passes to the returned distribution.
    !!}
    use :: Mass_Distributions, only : massDistributionClass, massDistributionComposite, massDistributionList, massDistributionListAcquire
    implicit none
    class  (massDistributionClass    ), pointer                :: massDistribution_
    logical                           , intent(in   )          :: includeStars              , includeGas
    class  (massDistributionClass    ), intent(inout), pointer :: massDistributionStellar   , massDistributionGas
    type   (massDistributionComposite), pointer                :: massDistributionTotal
    type   (massDistributionList     ), pointer                :: massDistributionComponents

    massDistribution_ => null()
    if      (includeStars .and. includeGas) then
       ! Composite the stellar and gas distributions.
       allocate(massDistributionTotal          )
       massDistributionComponents                        => massDistributionListAcquire()
       massDistributionComponents%next                   => massDistributionListAcquire()
       massDistributionComponents     %massDistribution_ => massDistributionStellar
       massDistributionComponents%next%massDistribution_ => massDistributionGas
       !![
       <referenceConstruct object="massDistributionTotal" constructor="massDistributionComposite(massDistributionComponents)"/>
       <objectDestructor name="massDistributionStellar"/>
       <objectDestructor name="massDistributionGas"    />
       !!]
       nullify(massDistributionComponents)
       massDistribution_ => massDistributionTotal
    else if (includeStars                 ) then
       ! Return just the stellar component.
       massDistribution_ => massDistributionStellar
    else if (                   includeGas) then
       ! Return just the gas component.
       massDistribution_ => massDistributionGas
    end if
    return
  end function Node_Component_Mass_Distribution_Combine

end module Node_Components_Galactic_Mass_Distributions
