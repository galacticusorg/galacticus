!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!!{
Contains a module which implements timers.
!!}

module Timers
  !!{
  Implements timers.
  !!}
  implicit none
  private
  public :: timer

  type :: timer
     !!{
     Type used to perform timing.
     !!}
     double precision :: timeStart, timeStop
   contains
     !![
     <methods>
       <method description="Start the timer."                                      method="start"     />
       <method description="Stop the timer."                                       method="stop"      />
       <method description="Report the time recorded as a double precision value." method="report"    />
       <method description="Report the time recorded as a character value."        method="reportText"/>
     </methods>
     !!]
     procedure :: start      => timerStart
     procedure :: stop       => timerStop
     procedure :: report     => timerReport
     procedure :: reportText => timerReportText
  end type timer

  interface timer
     !!{
     Constructors for the {\normalfont \ttfamily timer} class.
     !!}
     module procedure timerConstructorInternal
  end interface timer
  
contains

  function timerConstructorInternal() result(self)
    !!{
    Constructor for the {\normalfont \ttfamily timer} class.
    !!}
    type(timer) :: self

    self%timeStart=-huge(0.0d0)
    self%timeStop =-huge(0.0d0)
    return
  end function timerConstructorInternal

  subroutine timerStart(self)
    !!{
    Start the timer.
    !!}
    use :: OMP_Lib, only : OMP_Get_wTime
    implicit none
    class(timer), intent(inout) :: self
    
    self%timeStart=OMP_Get_wTime()
    return
  end subroutine timerStart
  
  subroutine timerStop(self)
    !!{
    Stop the timer.
    !!}
    use :: OMP_Lib, only : OMP_Get_wTime
    implicit none
    class(timer), intent(inout) :: self
    
    self%timeStop=OMP_Get_wTime()
    return
  end subroutine timerStop

  double precision function timerReport(self)
    !!{
    Report the time recorded by the time as a double precision value.
    !!}
    implicit none
    class(timer), intent(inout) :: self

    timerReport=+self%timeStop  &
         &      -self%timeStart
    return
  end function timerReport
  
  function timerReportText(self)
    !!{
    Report the time recorded by the time as a text value.
    !!}
    use :: Numerical_Constants_Prefixes, only : siFormat
    implicit none
    character(len=16)                :: timerReportText
    character(len=64)                :: timerReportText_
    class    (timer ), intent(inout) :: self

    timerReportText_=siFormat(self%report(),'f12.6')
    write (timerReportText,'(a,a1)') trim(timerReportText_),'s'
    return
  end function timerReportText

end module Timers
