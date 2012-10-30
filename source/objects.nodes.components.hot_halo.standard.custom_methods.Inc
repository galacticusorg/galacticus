!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
