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
Contains a module which reports on timestepping criteria.
!!}

module Evolve_To_Time_Reports
  !!{
  Contains functions which report on timestepping criteria.
  !!}
  implicit none
  private
  public :: Evolve_To_Time_Report

contains

  subroutine Evolve_To_Time_Report(message,time,index)
    !!{
    Display a report on evolution timestep criteria.
    !!}
    use :: Display           , only : displayMessage
    use :: ISO_Varying_String, only : assignment(=) , operator(//), varying_string
    use :: Kind_Numbers      , only : kind_int8
    use :: String_Handling   , only : operator(//)
    implicit none
    character       (len=*         ), intent(in   )           :: message
    double precision                , intent(in   )           :: time
    integer         (kind=kind_int8), intent(in   ), optional :: index
    type            (varying_string)                          :: vMessage
    character       (len=12        )                          :: label
    character       (len=32        )                          :: paddedMessage

    write (paddedMessage,'(a32  )') message
    write (label        ,'(e12.6)') time
    vMessage=paddedMessage//label
    if (present(index)) vMessage=vMessage//" ["//index//"]"
    call displayMessage(vMessage)
    return
  end subroutine Evolve_To_Time_Report

end module Evolve_To_Time_Reports
