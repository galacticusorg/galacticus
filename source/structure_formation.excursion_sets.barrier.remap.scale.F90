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


!% Contains a module which implements a scale remapping of excursion set barriers.

module Excursion_Sets_Barriers_Remap_Scale
  !% Implements a scaling remapping of excursion set barriers.
  private
  public :: Excursion_Sets_Barriers_Remap_Scale_Initialize, Excursion_Sets_Barrier_Remap_Scale,&
       & Excursion_Sets_Barrier_Gradient_Remap_Scale

  ! Record of the position of this remapping in the list of those to be applied.
  integer          :: methodPosition=-1,methodRatesPosition=-1
  
  ! Factor by which the barrier should be scaled.
  double precision :: excursionSetBarrierRemapScalingFactor
  logical          :: parametersInitialized=.false.
  
contains

  !# <excursionSetBarrierRemapInitialize>
  !#  <unitName>Excursion_Sets_Barriers_Remap_Scale_Initialize</unitName>
  !# </excursionSetBarrierRemapInitialize>
  subroutine Excursion_Sets_Barriers_Remap_Scale_Initialize(excursionSetBarrierRemapMethods,barrierName,ratesCalculation,matchedCount)
    !% Initialize the scale excursion set barrier remapping module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string), intent(in   ), dimension(:) :: excursionSetBarrierRemapMethods
    type(varying_string), intent(inout)               :: barrierName
    logical             , intent(in   )               :: ratesCalculation
    integer             , intent(inout)               :: matchedCount
    integer                                           :: i,position

    if (any(excursionSetBarrierRemapMethods == 'scale')) then
       ! Locate the position of the scale method in the list.
       position=-1
       do i=1,size(excursionSetBarrierRemapMethods)
          if (excursionSetBarrierRemapMethods(i) == 'scale') then
             position=i
             exit
          end if
       end do
       ! Record that our method is active.
       if (ratesCalculation) then
          methodRatesPosition=position
       else
          methodPosition     =position
       end if
       ! Increment the count of matched methods.
       matchedCount=matchedCount+1
       ! Construct a name for this barrier.
       barrierName=barrierName//":barrierRemapScale"
       ! Get the scaling factor.
       if (.not.parametersInitialized) then
          !@ <inputParameter>
          !@   <name>excursionSetBarrierRemapScalingFactor</name>
          !@   <defaultValue>$1$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The factor by which the excursion set barrier should be rescaled if the {\tt scale} remapping method is active.
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('excursionSetBarrierRemapScalingFactor',excursionSetBarrierRemapScalingFactor,defaultValue=1.0d0)
          parametersInitialized=.true.
       end if
    end if
    return
  end subroutine Excursion_Sets_Barriers_Remap_Scale_Initialize

  !# <excursionSetBarrierRemap>
  !#  <unitName>Excursion_Sets_Barrier_Remap_Scale</unitName>
  !# </excursionSetBarrierRemap>
  subroutine Excursion_Sets_Barrier_Remap_Scale(barrier,variance,time,ratesCalculation,iRemap)
    !% Return the barrier for excursion set calculations unmodified.
    implicit none
    double precision, intent(inout) :: barrier
    double precision, intent(in   ) :: variance,time
    logical         , intent(in   ) :: ratesCalculation
    integer         , intent(in   ) :: iRemap

    if ((ratesCalculation.and.iRemap == methodRatesPosition).or.(.not.ratesCalculation.and.iRemap == methodPosition)) barrier=barrier*excursionSetBarrierRemapScalingFactor
    return
  end subroutine Excursion_Sets_Barrier_Remap_Scale

  !# <excursionSetBarrierRemapGradient>
  !#  <unitName>Excursion_Sets_Barrier_Gradient_Remap_Scale</unitName>
  !# </excursionSetBarrierRemapGradient>
  subroutine Excursion_Sets_Barrier_Gradient_Remap_Scale(barrier,barrierGradient,variance,time,ratesCalculation,iRemap)
    !% Return the gradient of the barrier for excursion set calculations unmodified.
    implicit none
    double precision, intent(inout) :: barrierGradient
    double precision, intent(in   ) :: barrier,variance,time
    logical         , intent(in   ) :: ratesCalculation
    integer         , intent(in   ) :: iRemap

    if ((ratesCalculation.and.iRemap == methodRatesPosition).or.(.not.ratesCalculation.and.iRemap == methodPosition)) barrierGradient=barrierGradient*excursionSetBarrierRemapScalingFactor
    return
  end subroutine Excursion_Sets_Barrier_Gradient_Remap_Scale
  
end module Excursion_Sets_Barriers_Remap_Scale
