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

!% Contains a module which implements a \cite{sheth_ellipsoidal_2001} remapping of excursion set barriers.

module Excursion_Sets_Barriers_Remap_SMT
  !% Implements a \cite{sheth_ellipsoidal_2001} remapping of excursion set barriers.
  private
  public :: Excursion_Sets_Barriers_Remap_SMT_Initialize, Excursion_Sets_Barrier_Remap_SMT,&
       & Excursion_Sets_Barrier_Gradient_Remap_SMT

  ! Parameters of the remapping function
  double precision, parameter :: smtFitParameterA=0.707d0
  double precision, parameter :: smtFitParameterB=0.500d0
  double precision, parameter :: smtFitParameterC=0.600d0

  ! Record of the position of this remapping in the list of those to be applied.
  integer                     :: methodPosition  =-1     , methodRatesPosition=-1

contains

  !# <excursionSetBarrierRemapInitialize>
  !#  <unitName>Excursion_Sets_Barriers_Remap_SMT_Initialize</unitName>
  !# </excursionSetBarrierRemapInitialize>
  subroutine Excursion_Sets_Barriers_Remap_SMT_Initialize(excursionSetBarrierRemapMethods,barrierName,ratesCalculation,matchedCount)
    !% Initialize the \cite{sheth_ellipsoidal_2001} excursion set barrier remapping module.
    use ISO_Varying_String
    implicit none
    type   (varying_string), dimension(:), intent(in   ) :: excursionSetBarrierRemapMethods
    type   (varying_string)              , intent(inout) :: barrierName
    logical                              , intent(in   ) :: ratesCalculation
    integer                              , intent(inout) :: matchedCount
    integer                                              :: i                              , position

    if (any(excursionSetBarrierRemapMethods == 'Sheth-Mo-Tormen')) then
       ! Locate the position of the scale method in the list.
       position=-1
       do i=1,size(excursionSetBarrierRemapMethods)
          if (excursionSetBarrierRemapMethods(i) == 'Sheth-Mo-Tormen') then
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
       barrierName=barrierName//":barrierRemapShethMoTormen"
    end if
    return
  end subroutine Excursion_Sets_Barriers_Remap_SMT_Initialize

  !# <excursionSetBarrierRemap>
  !#  <unitName>Excursion_Sets_Barrier_Remap_SMT</unitName>
  !# </excursionSetBarrierRemap>
  subroutine Excursion_Sets_Barrier_Remap_SMT(barrier,variance,time,ratesCalculation,iRemap)
    !% Return the barrier for excursion set calculations remapped according to \cite{sheth_ellipsoidal_2001}.
    implicit none
    double precision, intent(inout) :: barrier
    double precision, intent(in   ) :: time            , variance
    logical         , intent(in   ) :: ratesCalculation
    integer         , intent(in   ) :: iRemap

    if ((ratesCalculation.and.iRemap == methodRatesPosition).or.(.not.ratesCalculation.and.iRemap == methodPosition)) barrier=sqrt(smtFitParameterA)*barrier*(1.0d0+smtFitParameterB*(variance/smtFitParameterA/barrier**2)**smtFitParameterC)
    return
  end subroutine Excursion_Sets_Barrier_Remap_SMT

  !# <excursionSetBarrierRemapGradient>
  !#  <unitName>Excursion_Sets_Barrier_Gradient_Remap_SMT</unitName>
  !# </excursionSetBarrierRemapGradient>
  subroutine Excursion_Sets_Barrier_Gradient_Remap_SMT(barrier,barrierGradient,variance,time,ratesCalculation,iRemap)
    !% Return the gradient of the barrier for excursion set calculations remapped according to \cite{sheth_ellipsoidal_2001}.
    implicit none
    double precision, intent(inout) :: barrierGradient
    double precision, intent(in   ) :: barrier         , time, variance
    logical         , intent(in   ) :: ratesCalculation
    integer         , intent(in   ) :: iRemap

    if ((ratesCalculation.and.iRemap == methodRatesPosition).or.(.not.ratesCalculation.and.iRemap == methodPosition)) then
       if (variance <= 0.0d0) then
          barrierGradient=0.0d0
       else
          barrierGradient=                                                                                                       &
               &           sqrt(smtFitParameterA)*barrierGradient*(                                                             &
               &                                                      1.0d0                                                      &
               &                                                     +smtFitParameterB                                           &
               &                                                     * (variance/smtFitParameterA/barrier**2)**smtFitParameterC) &
               &           +sqrt(smtFitParameterA)*barrier        *(                                                            &
               &                                                      smtFitParameterB                                           &
               &                                                     *smtFitParameterC                                           &
               &                                                     *((variance/smtFitParameterA/barrier**2)**smtFitParameterC) &
               &                                                     *(1.0d0-2.0d0*variance*barrierGradient/barrier)             &
               &                                                     /variance                                                   &
               &                                                    )
       end if
    end if
    return
  end subroutine Excursion_Sets_Barrier_Gradient_Remap_SMT

end module Excursion_Sets_Barriers_Remap_SMT
