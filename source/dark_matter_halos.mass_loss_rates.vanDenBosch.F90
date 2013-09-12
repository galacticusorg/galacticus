!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a calculation of dark matter halo mass loss rates using the method of
!% \cite{van_den_bosch_mass_2005}.

module Dark_Matter_Halos_Mass_Loss_Rates_vanDenBosch
  !% Implements a calculation of dark matter halo mass loss rates using the method of \cite{van_den_bosch_mass_2005}.
  implicit none
  private
  public :: Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch_Initialize

  ! Parameters of the mass loss model.
  double precision, parameter :: massLossTimescaleNormalization=0.13d0 !    Mass loss timescale normalization [Gyr].
  double precision, parameter :: zeta                          =0.36d0 !    Mass loss scaling with halo mass.

  ! Pre-computed mass loss rate normalization.
  double precision            :: massLossRateNormalization

contains

  !# <darkMatterHaloMassLossRateMethod>
  !#  <unitName>Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch_Initialize</unitName>
  !# </darkMatterHaloMassLossRateMethod>
  subroutine Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch_Initialize(darkMatterHaloMassLossRateMethod,Dark_Matter_Halos_Mass_Loss_Rate_Get)
    !% Initializes the ``vanDenBosch2005'' dark matter halo mass loss rate method.
    use ISO_Varying_String
    use Virial_Density_Contrast
    use Cosmology_Functions
    implicit none
    type     (varying_string                              ), intent(in   )          :: darkMatterHaloMassLossRateMethod
    procedure(Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch), intent(inout), pointer :: Dark_Matter_Halos_Mass_Loss_Rate_Get
    class    (cosmologyFunctionsClass                     )               , pointer :: cosmologyFunctionsDefault

    if (darkMatterHaloMassLossRateMethod == 'vanDenBosch2005') then
       ! Set a pointer to our implementation.
       Dark_Matter_Halos_Mass_Loss_Rate_Get => Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch
       ! Get the default cosmology functions object.
       cosmologyFunctionsDefault => cosmologyFunctions()
       ! Pre-compute the mass loss rate normalization factor.
       massLossRateNormalization=massLossTimescaleNormalization*sqrt(Halo_Virial_Density_Contrast(cosmologyFunctionsDefault%cosmicTime(1.0d0)))
    end if
    return
  end subroutine Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch_Initialize

  double precision function Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch(thisNode)
    !% Returns the rate of mass loss from dark matter halos using the prescription of \cite{van_den_bosch_mass_2005}.
    use Galacticus_Nodes
    use Virial_Density_Contrast
    use Cosmology_Functions
    implicit none
    type            (treeNode               ), intent(inout), pointer :: thisNode
    class           (nodeComponentBasic     )               , pointer :: parentBasic              , thisBasic
    class           (nodeComponentSatellite )               , pointer :: thisSatellite
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctionsDefault
    double precision                                                  :: massLossTimescale        , satelliteBoundMass, &
         &                                                               satelliteHostMassRatio   , satelliteTime

    thisSatellite      => thisNode     %satellite()
    satelliteBoundMass =  thisSatellite%boundMass()
    if (satelliteBoundMass > 0.0d0) then
       ! Get the default cosmology functions object.
       cosmologyFunctionsDefault => cosmologyFunctions()
      thisBasic             => thisNode %basic()
       satelliteTime         =  thisBasic%time ()
       massLossTimescale     =   massLossRateNormalization                                             &
            &                   *     cosmologyFunctionsDefault%expansionFactor(satelliteTime) **1.5d0 &
            &                   /sqrt(Halo_Virial_Density_Contrast             (satelliteTime))
       parentBasic           => thisNode%parent%basic()
       satelliteHostMassRatio=  satelliteBoundMass/parentBasic%mass()
       Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch=-satelliteBoundMass*satelliteHostMassRatio**zeta/massLossTimescale
    else
       Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch=0.0d0
    end if
    return
  end function Dark_Matter_Halos_Mass_Loss_Rate_vanDenBosch

end module Dark_Matter_Halos_Mass_Loss_Rates_vanDenBosch
