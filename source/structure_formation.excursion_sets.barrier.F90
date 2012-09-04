!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations of barriers for excursion set calculations.

module Excursion_Sets_Barriers
  !% Implements calculations of barriers for excursion set calculations.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Excursion_Sets_Barrier, Excursion_Sets_Barrier_Gradient, Excursion_Sets_Barrier_Name

  ! Flag to indicate if this module has been initialized.  
  logical                                         :: barrierModuleInitalized=.false.

  ! Name of method to use for barrier function.
  type(varying_string)                            :: excursionSetBarrierMethod
  type(varying_string), allocatable, dimension(:) :: excursionSetBarrierRemapMethods,excursionSetBarrierRatesRemapMethods

  ! Fully qualified name describing this barrier.
  type(varying_string)                                :: barrierName,barrierRatesName

  ! Pointer to the function that actually does the calculation.
  procedure(Excursion_Sets_Barrier         ), pointer :: Excursion_Sets_Barrier_Get          => null()
  procedure(Excursion_Sets_Barrier_Gradient), pointer :: Excursion_Sets_Barrier_Gradient_Get => null()
  
contains

  subroutine Excursion_Sets_Barrier_Initialize
    !% Initialize the excursion sets barrier module.
    use Input_Parameters
    use Galacticus_Error
    use String_Handling
    use Memory_Management
    !# <include directive="excursionSetBarrierMethod" type="moduleUse">
    include 'structure_formation.excursion_sets.barrier.moduleUse.inc'
    !# </include>
    !# <include directive="excursionSetBarrierRemapInitialize" type="moduleUse">
    include 'structure_formation.excursion_sets.barrier.remap.moduleUse.inc'
    !# </include>
    implicit none
    integer                                         :: remapCount,matchedCount
    logical                                         :: rateCalculation
    type(varying_string)                            :: name
    type(varying_string), allocatable, dimension(:) :: methods

    ! Initialize if necessary.
    if (.not.barrierModuleInitalized) then
       !$omp critical(Excursion_Sets_Barrier_Initialization) 
       if (.not.barrierModuleInitalized) then
          ! Get the barrier and remap method parameters.
          !@ <inputParameter>
          !@   <name>excursionSetBarrierMethod</name>
          !@   <defaultValue>criticalOverdensity</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of excursion set barriers.
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('excursionSetBarrierMethod',excursionSetBarrierMethod,defaultValue='criticalOverdensity')
          ! Include file that makes calls to all available method initialization routines.
          name=""
          !# <include directive="excursionSetBarrierMethod" type="code" action="subroutine">
          !#  <subroutineArgs>excursionSetBarrierMethod,Excursion_Sets_Barrier_Get,Excursion_Sets_Barrier_Gradient_Get,name</subroutineArgs>
          include 'structure_formation.excursion_sets.barrier.inc'
          !# </include>          
          if     (.not.(                                                              &
               &        associated(Excursion_Sets_Barrier_Get         ).and.          &
               &        associated(Excursion_Sets_Barrier_Gradient_Get)               &
               &       )                                                              &
               & )                                                                    &
               & call Galacticus_Error_Report(                                        &
               &       'Excursion_Sets_Barrier_Initialize'                            &
               &      ,'method '//char(excursionSetBarrierMethod)//' is unrecognized' &
               &                             )
          ! Record the barrier name.
          barrierName     =name
          barrierRatesName=name
          ! Determine how many remappings are to be applied.
          remapCount=Get_Input_Parameter_Array_Size('excursionSetBarrierRemapMethods')
          ! Allocate remap array and read remap names.
          !@ <inputParameter>
          !@   <name>excursionSetBarrierRemapMethods</name>
          !@   <defaultValue>null</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for remapping excursion set barriers.
          !@   </description>
          !@ </inputParameter>
          if (remapCount > 0) then
             allocate(excursionSetBarrierRemapMethods(remapCount))
             call Memory_Usage_Record(sizeof(excursionSetBarrierRemapMethods))
             call Get_Input_Parameter('excursionSetBarrierRemapMethods',excursionSetBarrierRemapMethods)
          else
             remapCount=1
             allocate(excursionSetBarrierRemapMethods(remapCount))
             call Memory_Usage_Record(sizeof(excursionSetBarrierRemapMethods))
             call Get_Input_Parameter('excursionSetBarrierRemapMethods',excursionSetBarrierRemapMethods,defaultValue=['null'])
          end if
          ! Include file that makes calls to all available method initialization routines.
          matchedCount   =0
          rateCalculation=.false.
          name           =""
          allocate(methods(size(excursionSetBarrierRemapMethods)))
          methods        =excursionSetBarrierRemapMethods
          !# <include directive="excursionSetBarrierRemapInitialize" type="code" action="subroutine">
          !#  <subroutineArgs>methods,name,rateCalculation,matchedCount</subroutineArgs>
          include 'structure_formation.excursion_sets.barrier.remap.initialize.inc'
          !# </include>
          deallocate(methods)
          if (matchedCount /= remapCount)                                                                                       &
               & call Galacticus_Error_Report(                                                                                  &
               &       'Excursion_Sets_Barrier_Initialize'                                                                      &
               &      ,'one or more of methods ['//char(String_Join(excursionSetBarrierRemapMethods,", "))//'] is unrecognized' &
               &                             )
          ! Record the barrier name.
          barrierName=barrierName//name
          ! Determine how many remappings are to be applied.
          remapCount=Get_Input_Parameter_Array_Size('excursionSetBarrierRatesRemapMethods')
          ! Allocate remap array and read remap names.
          !@ <inputParameter>
          !@   <name>excursionSetBarrierRatesRemapMethods</name>
          !@   <defaultValue>null</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for remapping excursion set barriers for rate calculations.
          !@   </description>
          !@ </inputParameter>
          if (remapCount > 0) then
             allocate(excursionSetBarrierRatesRemapMethods(remapCount))
             call Memory_Usage_Record(sizeof(excursionSetBarrierRatesRemapMethods))
             call Get_Input_Parameter('excursionSetBarrierRatesRemapMethods',excursionSetBarrierRatesRemapMethods)
          else
             remapCount=1
             allocate(excursionSetBarrierRatesRemapMethods(remapCount))
             call Memory_Usage_Record(sizeof(excursionSetBarrierRatesRemapMethods))
             call Get_Input_Parameter('excursionSetBarrierRatesRemapMethods',excursionSetBarrierRatesRemapMethods,defaultValue=['null'])
          end if
          ! Include file that makes calls to all available method initialization routines.
          matchedCount   =0
          rateCalculation=.true.
          name           =""
          allocate(methods(size(excursionSetBarrierRatesRemapMethods)))
          methods        =excursionSetBarrierRatesRemapMethods
          include 'structure_formation.excursion_sets.barrier.remap.initialize.inc'
          deallocate(methods)
          if (matchedCount /= remapCount)                                                                                            &
               & call Galacticus_Error_Report(                                                                                       &
               &       'Excursion_Sets_Barrier_Initialize'                                                                           &
               &      ,'one or more of methods ['//char(String_Join(excursionSetBarrierRatesRemapMethods,", "))//'] is unrecognized' &
               &                             )
          ! Record the barrier name.
          barrierRatesName=barrierRatesName//name
          ! Mark that the module is initalized.
          barrierModuleInitalized=.true.
       end if
       !$omp end critical(Excursion_Sets_Barrier_Initialization) 
    end if
    return
  end subroutine Excursion_Sets_Barrier_Initialize

  double precision function Excursion_Sets_Barrier(variance,time,ratesCalculation)
    !% Return the barrier for excursion sets at the given {\tt variance} and {\tt time}.
    include 'structure_formation.excursion_sets.barrier.remap.moduleUse.inc'
    implicit none
    double precision, intent(in)           :: variance,time
    logical         , intent(in), optional :: ratesCalculation
    logical                                :: ratesCalculationActual
    integer                                :: iRemap,remapCount

    ! Initialize the module if necessary.
    call Excursion_Sets_Barrier_Initialize

    ! Compute the original barrier.
    Excursion_Sets_Barrier=Excursion_Sets_Barrier_Get(variance,time)
    ! Remap the barrier.
    ratesCalculationActual=.false.
    if (present(ratesCalculation)) ratesCalculationActual=ratesCalculation
    if (ratesCalculationActual) then
       remapCount=size(excursionSetBarrierRatesRemapMethods)
    else
       remapCount=size(excursionSetBarrierRemapMethods     )
    end if

    do iRemap=1,remapCount
       !# <include directive="excursionSetBarrierRemap" type="code" action="subroutine">
       !#  <subroutineArgs>Excursion_Sets_Barrier,variance,time,ratesCalculationActual,iRemap</subroutineArgs>
       include 'structure_formation.excursion_sets.barrier.remap.inc'
       !# </include>
    end do
    return
  end function Excursion_Sets_Barrier

  double precision function Excursion_Sets_Barrier_Gradient(variance,time,ratesCalculation)
    !% Return the gradient (with respect to mass) of the barrier for excursion sets at the given {\tt variance} and {\tt time}.
    include 'structure_formation.excursion_sets.barrier.remap.moduleUse.inc'
    implicit none
    double precision, intent(in)           :: variance,time
    logical         , intent(in), optional :: ratesCalculation
    double precision                       :: barrier
    logical                                :: ratesCalculationActual
    integer                                :: iRemap,remapCount

    ! Initialize the module if necessary.
    call Excursion_Sets_Barrier_Initialize

    ! Compute the original barrier and its gradient.
    barrier                        =Excursion_Sets_Barrier_Get         (variance,time)
    Excursion_Sets_Barrier_Gradient=Excursion_Sets_Barrier_Gradient_Get(variance,time)
    ! Remap the barrier.
    ratesCalculationActual=.false.
    if (present(ratesCalculation)) ratesCalculationActual=ratesCalculation
    if (ratesCalculationActual) then
       remapCount=size(excursionSetBarrierRatesRemapMethods)
    else
       remapCount=size(excursionSetBarrierRemapMethods     )
    end if
    do iRemap=1,remapCount
       !# <include directive="excursionSetBarrierRemapGradient" type="code" action="subroutine">
       !#  <subroutineArgs>barrier,Excursion_Sets_Barrier_Gradient,variance,time,ratesCalculationActual,iRemap</subroutineArgs>
       include 'structure_formation.excursion_sets.barrier.remap.gradient.inc'
       !# </include>
    end do
    return
  end function Excursion_Sets_Barrier_Gradient

  function Excursion_Sets_Barrier_Name()
    !% Return the fully-qualified name of the selected excursion set barrier.
    type(varying_string) :: Excursion_Sets_Barrier_Name

    Excursion_Sets_Barrier_Name=barrierName
    return
  end function Excursion_Sets_Barrier_Name
  
end module Excursion_Sets_Barriers

