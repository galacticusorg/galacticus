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
Contains a module which implements computation of formatted dates and times.
!!}

module Dates_and_Times
  !!{
  Implements computation of formatted dates and times.
  !!}
  implicit none
  private
  public :: Formatted_Date_and_Time

  ! Names of days of the week and months of the year
  character(len=3), dimension(7),  parameter :: daysOfWeek=['Sun','Mon','Tue','Wed','Thu','Fri','Sat']
  character(len=3), dimension(12), parameter :: months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

  ! Tables used to determine what day of the week each month starts on.
  integer, dimension(12) :: monthTableLeapYear   =[6,2,3,6,1,4,6,2,5,0,3,5]
  integer, dimension(12) :: monthTableNonLeapYear=[0,3,3,6,1,4,6,2,5,0,3,5]

contains

  function Formatted_Date_and_Time()
    !!{
    Return a formatted date and time.
    !!}
    use :: ISO_Varying_String, only : assignment(=), varying_string, operator(//), extract, &
         &                            len
    use :: String_Handling   , only : operator(//)
    implicit none
    type     (varying_string) :: Formatted_Date_and_Time
    integer                   :: c                         , century      , &
         &                       dayOfWeek                 , m            , &
         &                       timeValues             (8), y            , &
         &                       year
    logical                   :: isLeapYear
    character(len=2         ) :: label
    character(len=2         ) :: hourOffset                , minuteOffset
    character(len=1         ) :: signOffset

    ! Get the date and time in numerical form.
    call Date_and_Time(values=timeValues)

    ! Compute the day of the week (following the algorithm at:
    ! http://en.wikipedia.org/wiki/Calculating_the_day_of_the_week#An_algorithm_to_calculate_the_day_of_the_week
    century=timeValues(1)/100
    c=2*(3-mod(century,4))
    year=timeValues(1)-100*century
    y=year+(year/4)
    isLeapYear=(mod(year,4)==0).and.(century/=0.or.mod(century,4)==0)
    select case (isLeapYear)
    case (.true.)
       m=monthTableLeapYear   (timeValues(2))
    case (.false.)
       m=monthTableNonLeapYear(timeValues(2))
    end select
    dayOfWeek=mod(c+y+m+timeValues(3),7)+1

    ! Construct a formatted string.
    Formatted_Date_and_Time=daysOfWeek(dayOfWeek)//" "//months(timeValues(2))//" "
    Formatted_Date_and_Time=Formatted_Date_and_Time//timeValues(3)//" "
    write (label,'(i2)') timeValues(5)
    Formatted_Date_and_Time=Formatted_Date_and_Time//label//":"
    write (label,'(i2.2)') timeValues(6)
    Formatted_Date_and_Time=Formatted_Date_and_Time//label//":"
    write (label,'(i2.2)') timeValues(7)
    Formatted_Date_and_Time=Formatted_Date_and_Time//label//" "
    write (hourOffset,'(i2.2)') abs(timeValues(4))/60
    write (minuteOffset,'(i2.2)') mod(abs(timeValues(4)),60)
    if (timeValues(4) > 0) then
       signOffset="+"
    else
       signOffset="-"
    end if
    Formatted_Date_and_Time=Formatted_Date_and_Time//"(GMT "//signOffset//hourOffset//":"//minuteOffset//") "
    Formatted_Date_and_Time=Formatted_Date_and_Time//timeValues(1)
    return
  end function Formatted_Date_and_Time

end module Dates_and_Times
