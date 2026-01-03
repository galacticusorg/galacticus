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
Contains a module which implements conversion of integers to Roman numerals.
!!}

module Numerical_Roman_Numerals
  !!{
  Implements conversion of integers to Roman numerals.
  !!}
  implicit none
  private
  public :: Roman_Numerals

  ! Paired list of smallest decimal value and the corresponding Roman numeral expression.
  integer         , dimension(13) :: decimal=[   1,   4,   5,   9,  10,  40,  50,  90, 100, 400, 500, 900,1000]
  character(len=2), dimension(13) :: roman  =[' I','IV',' V','IX',' X','XL',' L','XC',' C','CD',' D','CM',' M']
  
contains

  function Roman_Numerals(i) result(numeral)
    !!{
    Convert the given integer to Roman numerals.
    !!}
    use :: ISO_Varying_String, only : varying_string, var_str, operator(//)
    implicit none
    type   (varying_string)                :: numeral, numeral_
    integer                , intent(in   ) :: i
    integer                                :: i_     , j
    
    i_     =i
    numeral=var_str('')
    do j=13,1,-1       
       do while (i_ >= decimal(j))
          numeral_=numeral //trim(adjustl(roman  (j)))
          numeral =numeral_
          i_      =i_     -               decimal(j)
       end do
    end do
    return
  end function Roman_Numerals

end module Numerical_Roman_Numerals
