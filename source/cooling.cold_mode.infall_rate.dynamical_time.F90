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

!% Contains a module which implements a calculation of cold mode infall rates assuming infall on a dynamical timescale.

module Cooling_Cold_Mode_Infall_Rates_Dynamical_Time
  !% Implements a calculation of cold mode infall rates assuming infall on a dynamical timescale.
  use Galacticus_Nodes
  implicit none
  private
  public :: Cooling_Cold_Mode_Infall_Rate_Dynamical_Time_Initialize

  ! Infall rate.
  double precision :: coldModeInfallRateDynamicalTime

contains

  !# <coldModeInfallRateMethod>
  !#  <unitName>Cooling_Cold_Mode_Infall_Rate_Dynamical_Time_Initialize</unitName>
  !# </coldModeInfallRateMethod>
  subroutine Cooling_Cold_Mode_Infall_Rate_Dynamical_Time_Initialize(coldModeInfallRateMethod,Cooling_Cold_Mode_Infall_Rate_Get)
    !% Initializes the ``dynamical time'' cold mode infall rate module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type     (varying_string                              ), intent(in   )          :: coldModeInfallRateMethod
    procedure(Cooling_Cold_Mode_Infall_Rate_Dynamical_Time), intent(inout), pointer :: Cooling_Cold_Mode_Infall_Rate_Get

    if (coldModeInfallRateMethod == 'dynamicalTime') then
       Cooling_Cold_Mode_Infall_Rate_Get => Cooling_Cold_Mode_Infall_Rate_Dynamical_Time

       ! Get cooling rate parameters.
       !@ <inputParameter>
       !@   <name>coldModeInfallRateDynamicalTime</name>
       !@   <defaultValue>2.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The timescale (in units of the halo dynamical time) for infall of the cold mode component.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coldModeInfallRateDynamicalTime',coldModeInfallRateDynamicalTime,defaultValue=2.0d0)

       ! Check that the properties we need are gettable.
       if (.not.defaultHotHaloComponent%massColdIsGettable())                                         &
            & call Galacticus_Error_Report(                                                           &
            &                              'Cooling_Cold_Mode_Infall_Rate_Dynamical_Time_Initialize', &
            &                              'hot halo component must have gettable cold mass'          &
            &                             )
    end if
    return
  end subroutine Cooling_Cold_Mode_Infall_Rate_Dynamical_Time_Initialize

  double precision function Cooling_Cold_Mode_Infall_Rate_Dynamical_Time(thisNode)
    !% Computes the mass cooling rate in a hot gas halo assuming a fixed timescale for cooling.
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode                ), intent(inout), pointer :: thisNode
    class(nodeComponentHotHalo    )               , pointer :: thisHotHalo
    class(darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_

    thisHotHalo                                  => thisNode%hotHalo   ()
    darkMatterHaloScale_                         => darkMatterHaloScale()
    Cooling_Cold_Mode_Infall_Rate_Dynamical_Time =  coldModeInfallRateDynamicalTime*thisHotHalo%massCold()/darkMatterHaloScale_%dynamicalTimescale(thisNode)
    return
  end function Cooling_Cold_Mode_Infall_Rate_Dynamical_Time

end module Cooling_Cold_Mode_Infall_Rates_Dynamical_Time
