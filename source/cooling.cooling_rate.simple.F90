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

!% Contains a module which implements a simple cooling rate calculation in which the cooling rate equals the mass of hot gas
!% divided by a fixed timescale.

module Cooling_Rates_Simple
  !% Implements a simple cooling rate calculation in which the cooling rate equals the mass of hot gas
  !% divided by a fixed timescale.
  use Galacticus_Nodes
  implicit none
  private
  public :: Cooling_Rate_Simple_Initialize

  ! The fixed timescale for cooling.
  double precision :: coolingRateSimpleTimescale

contains

  !# <coolingRateMethod>
  !#  <unitName>Cooling_Rate_Simple_Initialize</unitName>
  !# </coolingRateMethod>
  subroutine Cooling_Rate_Simple_Initialize(coolingRateMethod,Cooling_Rate_Get)
    !% Initializes the ``simple'' cooling rate module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type     (varying_string     ), intent(in   )          :: coolingRateMethod
    procedure(Cooling_Rate_Simple), intent(inout), pointer :: Cooling_Rate_Get

    if (coolingRateMethod == 'simple') then
       Cooling_Rate_Get => Cooling_Rate_Simple

       ! Get cooling rate parameters.
       !@ <inputParameter>
       !@   <name>coolingRateSimpleTimescale</name>
       !@   <defaultValue>1 Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The timescale (in Gyr) for cooling in the simple cooling rate model.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingRateSimpleTimescale',coolingRateSimpleTimescale,defaultValue=1.0d0)

       ! Check that the properties we need are gettable.
       if (.not.defaultHotHaloComponent%massIsGettable())                                 &
            & call Galacticus_Error_Report(                                               &
            &                              'Cooling_Rate_Simple_Initialize'             , &
            &                              'hot halo component must have gettable mass'   &
            &                             )
    end if
    return
  end subroutine Cooling_Rate_Simple_Initialize

  double precision function Cooling_Rate_Simple(thisNode)
    !% Computes the mass cooling rate in a hot gas halo assuming a fixed timescale for cooling.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentHotHalo)               , pointer :: thisHotHaloComponent

    thisHotHaloComponent => thisNode%hotHalo()
    Cooling_Rate_Simple=thisHotHaloComponent%mass()/coolingRateSimpleTimescale
    return
  end function Cooling_Rate_Simple

end module Cooling_Rates_Simple
