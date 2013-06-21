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

!% Contains a module that implements calculations of the cooling time.

module Cooling_Times
  !% Implements calculations of the cooling function.
  use Abundances_Structure
  use Chemical_Abundances_Structure
  use Radiation_Structure
  use ISO_Varying_String
  !# <include directive="coolingTimeMethod" type="moduleUse">
  include 'cooling.cooling_time.modules.inc'
  !# </include>
  implicit none
  private
  public :: Cooling_Time, Cooling_Time_Density_Log_Slope, Cooling_Time_Temperature_Log_Slope

  ! Flag to indicate if this module has been initialized.
  logical                                       :: coolingTimeInitialized                =.false.

  ! Name of cooling time available method used.
  type     (varying_string           )          :: coolingTimeMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Cooling_Time_Get_Template), pointer :: Cooling_Time_Get                      =>null()
  procedure(Cooling_Time_Get_Template), pointer :: Cooling_Time_Density_Log_Slope_Get    =>null()
  procedure(Cooling_Time_Get_Template), pointer :: Cooling_Time_Temperature_Log_Slope_Get=>null()
  abstract interface
     double precision function Cooling_Time_Get_Template(temperature,density,gasAbundances,chemicalDensities,radiation)
       import abundances, radiationStructure, chemicalAbundances
       double precision                    , intent(in   ) :: density          , temperature
       type            (abundances        ), intent(in   ) :: gasAbundances
       type            (chemicalAbundances), intent(in   ) :: chemicalDensities
       type            (radiationStructure), intent(in   ) :: radiation
     end function Cooling_Time_Get_Template
  end interface

contains

  subroutine Cooling_Time_Initialize
    !% Initialize the cooling time module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.coolingTimeInitialized) then
       !$omp critical(Cooling_Time_Initialization)
       if (.not.coolingTimeInitialized) then
          ! Get the cooling time method parameter.
          !@ <inputParameter>
          !@   <name>coolingTimeMethod</name>
          !@   <defaultValue>simple</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be use for computing cooling times.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('coolingTimeMethod',coolingTimeMethod,defaultValue='simple')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="coolingTimeMethod" type="functionCall" functionType="void">
          !#  <functionArgs>coolingTimeMethod,Cooling_Time_Get,Cooling_Time_Density_Log_Slope_Get,Cooling_Time_Temperature_Log_Slope_Get</functionArgs>
          include 'cooling.cooling_time.inc'
          !# </include>
          if (.not.(associated(Cooling_Time_Get).and.associated(Cooling_Time_Density_Log_Slope_Get) &
               & .and.associated(Cooling_Time_Temperature_Log_Slope_Get))) call&
               & Galacticus_Error_Report('Cooling_Time','method ' //char(coolingTimeMethod)//' is unrecognized')
          coolingTimeInitialized=.true.
       end if
       !$omp end critical(Cooling_Time_Initialization)
    end if
    return
  end subroutine Cooling_Time_Initialize

  double precision function Cooling_Time(temperature,density,gasAbundances,chemicalDensities,radiation)
    !% Return the cooling time at the given temperature and density for the specified set of abundances and radiation
    !% field. Units of the returned cooling time are the Gyr.
    implicit none
    double precision                    , intent(in   ) :: density          , temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities
    type            (radiationStructure), intent(in   ) :: radiation

    ! Initialize the module if necessary.
    call Cooling_Time_Initialize

    ! Get the cooling time using the selected method.
    Cooling_Time=Cooling_Time_Get(temperature,density,gasAbundances,chemicalDensities,radiation)

    return
  end function Cooling_Time

  double precision function Cooling_Time_Density_Log_Slope(temperature,density,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic slope of the cooling time-density relation.
    implicit none
    double precision                    , intent(in   ) :: density          , temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities
    type            (radiationStructure), intent(in   ) :: radiation

    ! Initialize the module if necessary.
    call Cooling_Time_Initialize

    ! Get the cooling time using the selected method.
    Cooling_Time_Density_Log_Slope=Cooling_Time_Density_Log_Slope_Get(temperature,density,gasAbundances,chemicalDensities,radiation)

    return
  end function Cooling_Time_Density_Log_Slope

  double precision function Cooling_Time_Temperature_Log_Slope(temperature,density,gasAbundances,chemicalDensities,radiation)
    !% Return the logarithmic slope of the cooling time-temperature relation.
    implicit none
    double precision                    , intent(in   ) :: density          , temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances), intent(in   ) :: chemicalDensities
    type            (radiationStructure), intent(in   ) :: radiation

    ! Initialize the module if necessary.
    call Cooling_Time_Initialize

    ! Get the cooling time using the selected method.
    Cooling_Time_Temperature_Log_Slope=Cooling_Time_Temperature_Log_Slope_Get(temperature,density,gasAbundances,chemicalDensities,radiation)

    return
  end function Cooling_Time_Temperature_Log_Slope

end module Cooling_Times
