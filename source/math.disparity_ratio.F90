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
Contains a module which implements calculations of disparity ratios, i.e. $r(a,b)=\hbox{max}(|a|/|b|,|b|/|a|)$.
!!}

module Disparity_Ratios
  !!{
  Implements calculations of disparity ratios, i.e. $r(a,b)=\hbox{max}(|a|/|b|,|b|/|a|)$.
  !!}
  implicit none
  private
  public :: Disparity_Ratio

contains

  double precision function Disparity_Ratio(a,b)
    !!{
    Computes the disparity ratio, $r(a,b)=\hbox{max}(|a|/|b|,|b|/|a|)$.
    !!}
    implicit none
    double precision, intent(in   ) :: a       , b
    double precision                :: xMinimum, xMaximum
    
    xMinimum=min(abs(a),abs(b))
    xMaximum=max(abs(a),abs(b))
    if (xMinimum == 0.0d0 .and. xMaximum == 0.0d0) then
       Disparity_Ratio=1.0d0
    else if (exponent(xMaximum)-exponent(xMinimum) >= maxexponent(xMaximum) .or. xMinimum == 0.0d0) then
       Disparity_Ratio=huge(0.0d0)
    else
       Disparity_Ratio=abs(xMaximum/xMinimum)
    end if    
    return
  end function Disparity_Ratio
  
end module Disparity_Ratios
