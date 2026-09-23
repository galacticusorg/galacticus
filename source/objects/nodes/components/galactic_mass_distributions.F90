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
  dimensionless stellar and gas mass distributions to the size and masses of a particular component. The procedure is written
  once as a template and instantiated for cylindrical and spherical distributions.
  !!}
  implicit none
  private
  public :: Node_Component_Mass_Distribution_Scaled_Cylindrical, Node_Component_Mass_Distribution_Scaled_Spherical

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
    use :: Mass_Distributions, only : massDistributionClass, {geometry¦intrinsic}       , {geometry¦scaler}          , massDistributionComposite, &
         &                            massDistributionList , massDistributionListAcquire, kinematicsDistributionLocal
    use :: Object_Pools      , only : objectPool
    implicit none
    class           (massDistributionClass      ), pointer                :: massDistribution_
    logical                                      , intent(in   )          :: includeStars              , includeGas
    double precision                             , intent(in   )          :: massStellar               , massGas             , &
         &                                                                   radiusScale
    class           (massDistributionClass      ), intent(in   ), pointer :: massDistributionStellar_  , massDistributionGas_
    type            (kinematicsDistributionLocal), intent(in   ), pointer :: kinematicDistribution_
    type            (objectPool                 ), intent(inout)          :: scalerStellarPool         , scalerGasPool
    type            ({geometry¦scaler}          ), pointer                :: massDistributionStellar   , massDistributionGas
    type            (massDistributionComposite  ), pointer                :: massDistributionTotal
    type            (massDistributionList       ), pointer                :: massDistributionComponents
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
    end if
    return
  end function Node_Component_Mass_Distribution_Scaled_{geometry¦label}

end module Node_Components_Galactic_Mass_Distributions
