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

!% Contains custom functions for the standard hot halo component.

subroutine Node_Component_Hot_Halo_Standard_Mass_Sink(self,setValue)
  !% Account for a sink of gaseous material in the standard hot halo hot gas.
  use Galacticus_Error
  implicit none
  class(nodeComponentHotHaloStandard), intent(inout) :: self
  double precision                   , intent(in   ) :: setValue
  
  ! Trap cases where an attempt is made to add gas via this sink function.
  if (setValue > 0.0d0) call Galacticus_Error_Report('Node_Component_Hot_Halo_Standard_Mass_Sink','attempt to add mass via sink in hot halo')
  
  ! Proportionally adjust the rates of all components of the hot gas reservoir.
  call Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate(self,setValue)
  
  return
end subroutine Node_Component_Hot_Halo_Standard_Mass_Sink

subroutine Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate(self,gasMassRate)
  !% Adjusts the rates of all components of the hot gas reservoir under the assumption of uniformly distributed properties
  !% (e.g. fully-mixed metals).
  implicit none
  class(nodeComponentHotHaloStandard), intent(inout) :: self
  double precision                   , intent(in   ) :: gasMassRate
  double precision                                   :: gasMass
  
  ! Exit immediately for zero rate.
  if (gasMassRate == 0.0d0) return
  ! Get the gas mass present.
  gasMass=self%mass()
  ! If gas is present, adjust the rates.
  if (gasMass > 0.0d0) then
     ! Mass.
     call self%           massRate(                       gasMassRate        )
     ! Angular momentum.
     call self%angularMomentumRate(self%angularMomentum()*gasMassRate/gasMass)
     ! Metal abundances.
     call self%     abundancesRate(self%abundances     ()*gasMassRate/gasMass)
     ! Chemical abundances.
     call self%      chemicalsRate(self%chemicals      ()*gasMassRate/gasMass)
  end if
  return
end subroutine Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate
