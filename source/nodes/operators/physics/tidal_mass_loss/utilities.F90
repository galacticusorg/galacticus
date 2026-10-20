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
Contains a module which applies tidal mass loss to galactic components.
!!}

module Tidal_Mass_Loss_Utilities
  !!{RST
  Applies tidal mass loss to galactic components. Since the disk, spheroid, and nuclear star cluster components share no parent
  type below :galacticus-type:`nodeComponent`, the procedure is written once as a template and instantiated for each component.
  !!}
  implicit none
  private
  public :: Tidal_Mass_Loss_Apply_Disk, Tidal_Mass_Loss_Apply_Spheroid, Tidal_Mass_Loss_Apply_NSC

  !![
  <generic identifier="tidalcomponent">
   <instance label="Disk"     intrinsic="nodeComponentDisk"     accessor="disk"     name="disk"                />
   <instance label="Spheroid" intrinsic="nodeComponentSpheroid" accessor="spheroid" name="spheroid"            />
   <instance label="NSC"      intrinsic="nodeComponentNSC"      accessor="NSC"      name="nuclear star cluster"/>
  </generic>
  !!]

contains

  subroutine Tidal_Mass_Loss_Apply_{tidalcomponent¦label}(node,tidalStripping_,tidalTruncationID)
    !!{RST
    Apply the rates of tidal mass loss from the {tidalcomponent¦name} of ``node``: gas is moved to the outflowed reservoir of the
    hot halo, stellar mass is removed, and angular momentum is lost at the specific angular momentum of the stripped material
    (as given by ``tidalStripping_``). The degree of tidal truncation of the {tidalcomponent¦name}, :math:`x=-\log f` where
    :math:`f` is the fraction of the mass of its untruncated mass distribution which it retains (stored in the meta-property
    with ID ``tidalTruncationID``), is increased at rate :math:`\dot{x} = \dot{M}_\mathrm{strip}/M`.
    !!}
    use :: Abundances_Structure          , only : operator(*)
    use :: Galacticus_Nodes              , only : {tidalcomponent¦intrinsic}, nodeComponentHotHalo, treeNode
    use :: Histories                     , only : operator(*)               , history
    use :: Stellar_Luminosities_Structure, only : operator(*)               , stellarLuminosities , zeroStellarLuminosities, max
    use :: Tidal_Stripping_Mass_Loss_Rate, only : tidalStrippingClass
    implicit none
    type            (treeNode                  ), intent(inout), target  :: node
    class           (tidalStrippingClass       ), intent(inout)          :: tidalStripping_
    integer                                     , intent(in   )          :: tidalTruncationID
    class           ({tidalcomponent¦intrinsic})               , pointer :: component
    class           (nodeComponentHotHalo      )               , pointer :: hotHalo
    type            (stellarLuminosities       ), save                   :: luminositiesTransferRate
    !$omp threadprivate(luminositiesTransferRate)
    double precision                                                     :: fractionGas             , fractionStellar             , &
         &                                                                  massLossRate            , ratioSpecificAngularMomentum
    type            (history                   )                         :: historyTransferRate

    ! Return if the component has no mass.
    component => node%{tidalcomponent¦accessor}()
    if (component%massGas()+component%massStellar() <= 0.0d0) return
    ! Return if the tidal mass loss rate is zero.
    massLossRate=tidalStripping_%rateMassLoss(component,ratioSpecificAngularMomentum)
    if (massLossRate                                <= 0.0d0) return
    ! Transfer stripped material from the component.
    !! Gas is moved to the hot halo component.
    hotHalo         => node%hotHalo()
    fractionGas     =  min(1.0d0,max(0.0d0,component%massGas()/(component%massGas()+component%massStellar())))
    fractionStellar =  1.0d0-fractionGas
    if (fractionGas     > 0.0d0 .and. component%massGas    () > 0.0d0) then
       call component       %                  massGasRate(-fractionGas    *massLossRate                                                                                                         )
       call component       %            abundancesGasRate(-fractionGas    *massLossRate*component%abundancesGas    ()/ component%massGas()                                                      )
       call   hotHalo       %           outflowingMassRate(+fractionGas    *massLossRate                                                                                                         )
       call   hotHalo       %     outflowingAbundancesRate(+fractionGas    *massLossRate*component%abundancesGas    ()/ component%massGas()                                                      )
       call   hotHalo       %outflowingAngularMomentumRate(+fractionGas    *massLossRate*component%angularMomentum  ()/(component%massGas()+component%massStellar())*ratioSpecificAngularMomentum)
    end if
    ! Stellar mass is simply removed.
    if (fractionStellar > 0.0d0 .and. component%massStellar() > 0.0d0) then
       ! If luminosities are being treated as inactive properties this is an error - they appear on the right-hand side
       ! of the following ODE terms so are not inactive. (An approach similar to what is used for transfer of
       ! luminosities to the spheroid by bar instabilities could work here.)
       !! Stellar mass and metals.
       call component       %              massStellarRate(-fractionStellar*massLossRate                                                                                                         )
       call component       %        abundancesStellarRate(-fractionStellar*massLossRate*component%abundancesStellar()/                     component%massStellar()                              )
       !! Stellar luminosities.
       luminositiesTransferRate=max(zeroStellarLuminosities,component%luminositiesStellar())
       call component       %      luminositiesStellarRate(-fractionStellar*massLossRate*luminositiesTransferRate     /                     component%massStellar()                              )
       !! Stellar properties history.
       historyTransferRate=component%stellarPropertiesHistory()
       if (historyTransferRate%exists()) &
            & call component% stellarPropertiesHistoryRate(-fractionStellar*massLossRate*historyTransferRate          /                     component%massStellar()                              )
       call historyTransferRate%destroy()
       !! Star formation history.
       historyTransferRate=component%starFormationHistory()
       if (historyTransferRate%exists()) &
            & call component% starFormationHistoryRate    (-fractionStellar*massLossRate*historyTransferRate          /                     component%massStellar()                              )
    end if
    ! Angular momentum is lost at the specific angular momentum of the stripped material.
    call component          % angularMomentumRate         (-                massLossRate*component%angularMomentum  ()/(component%massGas()+component%massStellar())*ratioSpecificAngularMomentum)
    ! The degree of truncation increases in proportion to the fractional rate of mass loss.
    call component          % floatRank0MetaPropertyRate  (tidalTruncationID,+massLossRate                            /(component%massGas()+component%massStellar())                             )
    return
  end subroutine Tidal_Mass_Loss_Apply_{tidalcomponent¦label}

end module Tidal_Mass_Loss_Utilities
