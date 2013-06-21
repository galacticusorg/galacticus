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

!% Contains a module which calculates the specific angular momentum of cooling gas assuming all gas has the mean specific angular
!% momentum of the hot gas halo.

module Cooling_Specific_Angular_Momenta_Mean
  !% Calculates the specific angular momentum of cooling gas assuming all gas has the mean specific angular momentum of the hot gas
  !% halo.
  implicit none
  private
  public :: Cooling_Specific_AM_Mean_Initialize

contains

  !# <coolingSpecificAngularMomentumMethod>
  !#  <unitName>Cooling_Specific_AM_Mean_Initialize</unitName>
  !# </coolingSpecificAngularMomentumMethod>
  subroutine Cooling_Specific_AM_Mean_Initialize(coolingSpecificAngularMomentumMethod,Cooling_Specific_Angular_Momentum_Get)
    !% Initializes the ``mean'' specific angular momentum of cooling gas module.
    use Galacticus_Nodes
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type     (varying_string                        ), intent(in   )          :: coolingSpecificAngularMomentumMethod   
    procedure(Cooling_Specific_Angular_Momentum_Mean), intent(inout), pointer :: Cooling_Specific_Angular_Momentum_Get  
                                                                                                                     
    if (coolingSpecificAngularMomentumMethod == 'mean') then
       Cooling_Specific_Angular_Momentum_Get => Cooling_Specific_Angular_Momentum_Mean
       ! Check that the required properties are gettable.
       if     (                                                               &
            &  .not.(                                                         &
            &        defaultHotHaloComponent%           massIsGettable().and. &
            &        defaultHotHaloComponent%angularMomentumIsGettable()      &
            &       )                                                         &
            & ) call Galacticus_Error_Report('Cooling_Specific_AM_Mean_Initialize','this method requires that the "mass" and "angularMomentum" properties of the hot halo be gettable')
    end if
    return
  end subroutine Cooling_Specific_AM_Mean_Initialize

  double precision function Cooling_Specific_Angular_Momentum_Mean(thisNode,radius)
    !% Return the specific angular momentum of cooling gas in the mean model.
    use Galacticus_Nodes
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode              
    double precision                      , intent(in   )          :: radius                
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent  
    
    ! Compute mean specific angular momentum from the hot halo component.                                                                                     
    thisHotHaloComponent => thisNode%hotHalo()
    Cooling_Specific_Angular_Momentum_Mean= thisHotHaloComponent%angularMomentum() &
         &                                 /thisHotHaloComponent%mass           ()
    return
  end function Cooling_Specific_Angular_Momentum_Mean
  
end module Cooling_Specific_Angular_Momenta_Mean
