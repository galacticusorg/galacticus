!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a linear barrier for excursion set calculations of dark matter halo formation.

module Excursion_Sets_Barriers_Linear
  !% Implements a linear barrier for excursion set calculations of dark matter halo formation.
  private
  public :: Excursion_Sets_Barriers_Linear_Initialize

  ! Parameters controlling the barrier.
  double precision :: excursionSetBarrierConstantCoefficient,excursionSetBarrierLinearCoefficient

contains

  !# <excursionSetBarrierMethod>
  !#  <unitName>Excursion_Sets_Barriers_Linear_Initialize</unitName>
  !# </excursionSetBarrierMethod>
  subroutine Excursion_Sets_Barriers_Linear_Initialize(excursionSetBarrierMethod,Excursion_Sets_Barrier_Get,Excursion_Sets_Barrier_Gradient_Get,barrierName)
    !% Initialize the linear excursion set barrier module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: excursionSetBarrierMethod
    procedure(double precision), pointer, intent(inout) :: Excursion_Sets_Barrier_Get,Excursion_Sets_Barrier_Gradient_Get
    type(varying_string),                 intent(inout) :: barrierName
    character(len=10)                                   :: label

    if (excursionSetBarrierMethod == 'linear') then
       Excursion_Sets_Barrier_Get          => Excursion_Sets_Barrier_Linear
       Excursion_Sets_Barrier_Gradient_Get => Excursion_Sets_Barrier_Gradient_Linear
       !@ <inputParameter>
       !@   <name>excursionSetBarrierConstantCoefficient</name>
       !@   <defaultValue>1.67</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The constant term in the excursion set barrier.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('excursionSetBarrierConstantCoefficient',excursionSetBarrierConstantCoefficient,defaultValue=1.67d0)
       !@ <inputParameter>
       !@   <name>excursionSetBarrierLinearCoefficient</name>
       !@   <defaultValue>0.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The coefficient of the linear term in the excursion set barrier.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('excursionSetBarrierLinearCoefficient',excursionSetBarrierLinearCoefficient,defaultValue=0.0d0)
       ! Construct a name for this barrier.
       write (label,'(e10.4)') excursionSetBarrierConstantCoefficient
       barrierName=barrierName//":barrierLinear:constantCoefficient:"//label
       write (label,'(e10.4)') excursionSetBarrierLinearCoefficient
       barrierName=barrierName//":linearCoefficient"//label
    end if
    return
  end subroutine Excursion_Sets_Barriers_Linear_Initialize

  double precision function Excursion_Sets_Barrier_Linear(variance,time)
    !% Return a linear barrier for excursion set calculations at the given {\tt variance}.
    implicit none
    double precision, intent(in) :: variance,time

    Excursion_Sets_Barrier_Linear=excursionSetBarrierConstantCoefficient+excursionSetBarrierLinearCoefficient*variance
    return
  end function Excursion_Sets_Barrier_Linear

  double precision function Excursion_Sets_Barrier_Gradient_Linear(variance,time)
    !% Return the gradient of a linear barrier for excursion set calculations at the given {\tt variance}.
    implicit none
    double precision, intent(in) :: variance,time

    Excursion_Sets_Barrier_Gradient_Linear=excursionSetBarrierLinearCoefficient
    return
  end function Excursion_Sets_Barrier_Gradient_Linear
  
end module Excursion_Sets_Barriers_Linear
