!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements a simple cooling radius calculation (finds the radius at which the time available for
!% cooling equals the cooling time).

module Cooling_Radii_Simple
  !% Implements a simple cooling radius calculation (finds the radius at which the time available for cooling equals the cooling
  !% time).
  use, intrinsic :: ISO_C_Binding
  use Tree_Nodes
  use Kind_Numbers
  private
  public :: Cooling_Radius_Simple_Initialize, Cooling_Radius_Simple_Reset

  ! Module global variable that stores the time available for cooling.
  double precision :: coolingTimeAvailable
  !$omp threadprivate(coolingTimeAvailable)

  ! Module global pointer to the active node.
  type(treeNode), pointer :: activeNode
  !$omp threadprivate(activeNode)

  ! Internal record of the number of abundance properties.
  integer :: abundancesCount

  ! Record of unique ID of node which we last computed results for.
  integer(kind=kind_int8) :: lastUniqueID=-1
  !$omp threadprivate(lastUniqueID)

  ! Record of whether or not cooling radius has already been computed for this node.
  logical :: coolingRadiusComputed=.false.,coolingRadiusGrowthRateComputed=.false.
  !$omp threadprivate(coolingRadiusComputed,coolingRadiusGrowthRateComputed)

  ! Stored values of cooling radius.
  double precision :: coolingRadiusStored,coolingRadiusGrowthRateStored
  !$omp threadprivate(coolingRadiusStored,coolingRadiusGrowthRateStored)

contains

  !# <coolingRadiusMethod>
  !#  <unitName>Cooling_Radius_Simple_Initialize</unitName>
  !# </coolingRadiusMethod>
  subroutine Cooling_Radius_Simple_Initialize(coolingRadiusMethod,Cooling_Radius_Get,Cooling_Radius_Growth_Rate_Get)
    !% Initializes the ``simple'' cooling radius module.
    use ISO_Varying_String
    use Abundances_Structure
    implicit none
    type(varying_string),          intent(in)    :: coolingRadiusMethod
    procedure(),          pointer, intent(inout) :: Cooling_Radius_Get,Cooling_Radius_Growth_Rate_Get
    
    if (coolingRadiusMethod == 'simple') then
       Cooling_Radius_Get => Cooling_Radius_Simple
       Cooling_Radius_Growth_Rate_Get => Cooling_Radius_Growth_Rate_Simple
       ! Get a count of the number of abundances properties.
       abundancesCount=Abundances_Property_Count()
    end if
    return
  end subroutine Cooling_Radius_Simple_Initialize

  !# <calculationResetTask>
  !# <unitName>Cooling_Radius_Simple_Reset</unitName>
  !# </calculationResetTask>
  subroutine Cooling_Radius_Simple_Reset(thisNode)
    !% Reset the cooling radius calculation.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    coolingRadiusComputed          =.false.
    coolingRadiusGrowthRateComputed=.false.
    lastUniqueID                   =thisNode%uniqueID()
    return
  end subroutine Cooling_Radius_Simple_Reset

  double precision function Cooling_Radius_Growth_Rate_Simple(thisNode)
    !% Return the growth rate of the cooling radius in the ``simple'' model in Mpc/Gyr.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Hot_Halo_Temperature_Profile
    use Hot_Halo_Density_Profile
    use Cooling_Times
    use Abundances_Structure
    use Radiation_Structure
    use Cooling_Times_Available
    implicit none
    type(treeNode),            intent(inout), pointer     :: thisNode
    double precision,          dimension(abundancesCount) :: abundancesMassFraction
    double precision                                      :: virialRadius,coolingRadius,coolingTimeAvailable &
         &,coolingTimeAvailableIncreaseRate,densityLogSlope,temperatureLogSlope,density,temperature,coolingTimeDensityLogSlope &
         &,coolingTimeTemperatureLogSlope
    type(abundancesStructure), save                       :: abundances
    !$omp threadprivate(abundances)
    type(radiationStructure)                              :: radiation

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Cooling_Radius_Simple_Reset(thisNode)

    ! Check if cooling radius growth rate is already computed.
    if (.not.coolingRadiusGrowthRateComputed) then
       ! Flag that cooling radius is now computed.
       coolingRadiusGrowthRateComputed=.true.
       
       ! Get the virial radius.
       virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
       
       ! Get the cooling radius.
       coolingRadius=Cooling_Radius_Simple(thisNode)
       
       ! Check if cooling radius has reached virial radius.
       if (coolingRadius >= virialRadius) then
          coolingRadiusGrowthRateStored=0.0d0
       else
          ! Get the time available for cooling in thisNode.
          coolingTimeAvailable=Cooling_Time_Available(thisNode)
          
          ! Get the rate of increase of the time available for cooling.
          coolingTimeAvailableIncreaseRate=Cooling_Time_Available_Increase_Rate(thisNode)
          
          ! Logarithmic slope of density profile.
          densityLogSlope=Hot_Halo_Density_Log_Slope(thisNode,coolingRadius)
          
          ! Logarithmic slope of density profile.
          temperatureLogSlope=Hot_Halo_Temperature_Logarithmic_Slope(thisNode,coolingRadius)
          
          ! Get cooling density, temperature and metallicity.
          density=Hot_Halo_Density(activeNode,coolingRadius)
          temperature=Hot_Halo_Temperature(activeNode,coolingRadius)
          
          ! Set the radiation field.
          call radiation%setCMB(Tree_Node_Time(thisNode))
          
          ! Get the abundances for this node.
          call Tree_Node_Hot_Halo_Abundances(thisNode,abundancesMassFraction)
          call abundances%pack(abundancesMassFraction)
          call abundances%massToMassFraction(Tree_Node_Hot_Halo_Mass(thisNode))

          ! Logarithmic slope of the cooling time-density relation.
          coolingTimeDensityLogSlope=Cooling_Time_Density_Log_Slope(temperature,density,abundances,radiation)
          
          ! Logarithmic slope of the cooling time-temperature relation.
          coolingTimeTemperatureLogSlope=Cooling_Time_Temperature_Log_Slope(temperature,density,abundances,radiation)
          
          ! Compute rate at which cooling radius grows.
          if (coolingRadius > 0.0d0) then
             coolingRadiusGrowthRateStored=(coolingRadius/coolingTimeAvailable)*coolingTimeAvailableIncreaseRate&
                  &/(densityLogSlope*coolingTimeDensityLogSlope+temperatureLogSlope*coolingTimeTemperatureLogSlope)
          else
             coolingRadiusGrowthRateStored=0.0d0
          end if
       end if
    end if
    ! Return the stored value.
    Cooling_Radius_Growth_Rate_Simple=coolingRadiusGrowthRateStored       
    return
  end function Cooling_Radius_Growth_Rate_Simple

  double precision function Cooling_Radius_Simple(thisNode)
    !% Return the cooling radius in the simple model.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Cooling_Times_Available
    use Root_Finder
    use FGSL
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    double precision,        parameter              :: zeroRadius=0.0d0
    type(fgsl_function),     save                   :: rootFunction
    type(fgsl_root_fsolver), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision,        parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6
    type(c_ptr)                                     :: parameterPointer
    double precision                                :: virialRadius 

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Cooling_Radius_Simple_Reset(thisNode)

    ! Check if cooling radius is already computed.
    if (.not.coolingRadiusComputed) then
       ! Flag that cooling radius is now computed.
       coolingRadiusComputed=.true.

       ! Get the time available for cooling in thisNode.
       coolingTimeAvailable=Cooling_Time_Available(thisNode)
       
       ! Make the node available to the root finding routine.
       activeNode => thisNode
       
       ! Check if cooling time at halo center is reached.
       if (Cooling_Radius_Root(zeroRadius,parameterPointer) > 0.0d0) then
          ! Cooling time at halo center exceeds the time available, return zero radius.
          coolingRadiusStored=zeroRadius
          Cooling_Radius_Simple=coolingRadiusStored
         return
       end if
       
       ! Check if cooling time at halo virial radius is reached.
       virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
       if (Cooling_Radius_Root(virialRadius,parameterPointer) < 0.0d0) then
          ! Cooling time available exceeds cooling time at virial radius, return virial radius.
          coolingRadiusStored=virialRadius
          Cooling_Radius_Simple=coolingRadiusStored
          return
       end if
       
       ! Cooling radius is between zero and virial radii. Search for the virial radius.
       coolingRadiusStored=Root_Find(zeroRadius,virialRadius,Cooling_Radius_Root,parameterPointer,rootFunction,rootFunctionSolver &
            &,toleranceAbsolute,toleranceRelative)
       Cooling_Radius_Simple=coolingRadiusStored
    else
       Cooling_Radius_Simple=coolingRadiusStored
    end if
    return
  end function Cooling_Radius_Simple
  
  function Cooling_Radius_Root(radius,parameterPointer) bind(c)
    !% Root function which evaluates the difference between the cooling time at {\tt radius} and the time available for cooling.
    use Cooling_Times
    use Abundances_Structure
    use Radiation_Structure
    use Hot_Halo_Density_Profile
    use Hot_Halo_Temperature_Profile
    implicit none
    real(c_double)                                        :: Cooling_Radius_Root
    real(c_double),            value                      :: radius
    type(c_ptr),               value                      :: parameterPointer
    double precision,          dimension(abundancesCount) :: abundancesMassFraction
    double precision                                      :: coolingTime,density,temperature
    type(abundancesStructure), save                       :: abundances
    !$omp threadprivate(abundances)
    type(radiationStructure)                              :: radiation

    ! Compute density, temperature and abundances.
    density    =Hot_Halo_Density    (activeNode,radius)
    temperature=Hot_Halo_Temperature(activeNode,radius)
 
    ! Get the abundances for this node.
    call Tree_Node_Hot_Halo_Abundances(activeNode,abundancesMassFraction)
    call abundances%pack(abundancesMassFraction)
    call abundances%massToMassFraction(Tree_Node_Hot_Halo_Mass(activeNode))
    ! Set the radiation field.
    call radiation%setCMB(Tree_Node_Time(activeNode))
    ! Compute the cooling time at the specified radius.
    coolingTime=Cooling_Time(temperature,density,abundances,radiation)
    ! Return the difference between cooling time and time available.
    Cooling_Radius_Root=coolingTime-coolingTimeAvailable

    return
  end function Cooling_Radius_Root

end module Cooling_Radii_Simple
