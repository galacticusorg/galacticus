!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a null remapping of excursion set barriers.

module Excursion_Sets_Barriers_Remap_Null
  !% Implements a null remapping of excursion set barriers.
  private
  public :: Excursion_Sets_Barriers_Remap_Null_Initialize, Excursion_Sets_Barrier_Remap_Null,&
       & Excursion_Sets_Barrier_Gradient_Remap_Null

  ! Record of whether this remap method is active.
  logical :: methodIsActive=.false., methodRatesIsActive=.false.

contains

  !# <excursionSetBarrierRemapInitialize>
  !#  <unitName>Excursion_Sets_Barriers_Remap_Null_Initialize</unitName>
  !# </excursionSetBarrierRemapInitialize>
  subroutine Excursion_Sets_Barriers_Remap_Null_Initialize(excursionSetBarrierRemapMethods,barrierName,ratesCalculation,matchedCount)
    !% Initialize the null excursion set barrier remapping module.
    use ISO_Varying_String
    implicit none
    type   (varying_string), dimension(:), intent(in   ) :: excursionSetBarrierRemapMethods
    type   (varying_string)              , intent(inout) :: barrierName
    logical                              , intent(in   ) :: ratesCalculation
    integer                              , intent(inout) :: matchedCount

    if (any(excursionSetBarrierRemapMethods == 'null')) then
       ! Record that our method is active.
       if (ratesCalculation) then
          methodRatesIsActive=.true.
       else
          methodIsActive     =.true.
       end if
       ! Increment the count of matched methods.
       matchedCount=matchedCount+1
       ! Construct a name for this barrier.
       barrierName=barrierName//":barrierRemapNull"
    end if
    return
  end subroutine Excursion_Sets_Barriers_Remap_Null_Initialize

  !# <excursionSetBarrierRemap>
  !#  <unitName>Excursion_Sets_Barrier_Remap_Null</unitName>
  !# </excursionSetBarrierRemap>
  subroutine Excursion_Sets_Barrier_Remap_Null(barrier,variance,time,ratesCalculation,iRemap)
    !% Return the barrier for excursion set calculations unmodified.
    implicit none
    double precision, intent(inout) :: barrier
    double precision, intent(in   ) :: time            , variance
    logical         , intent(in   ) :: ratesCalculation
    integer         , intent(in   ) :: iRemap

    return
  end subroutine Excursion_Sets_Barrier_Remap_Null

  !# <excursionSetBarrierRemapGradient>
  !#  <unitName>Excursion_Sets_Barrier_Gradient_Remap_Null</unitName>
  !# </excursionSetBarrierRemapGradient>
  subroutine Excursion_Sets_Barrier_Gradient_Remap_Null(barrier,barrierGradient,variance,time,ratesCalculation,iRemap)
    !% Return the gradient of the barrier for excursion set calculations unmodified.
    implicit none
    double precision, intent(inout) :: barrierGradient
    double precision, intent(in   ) :: barrier         , time, variance
    logical         , intent(in   ) :: ratesCalculation
    integer         , intent(in   ) :: iRemap

    return
  end subroutine Excursion_Sets_Barrier_Gradient_Remap_Null

end module Excursion_Sets_Barriers_Remap_Null
