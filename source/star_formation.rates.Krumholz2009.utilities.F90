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

  !+    Contributions to this file made by: Matías Liempi

  !!{  
  Contains a module that implements various useful utility functions for calculations for the \cite{krumholz_star_2009}
  star formation surface density rate law.
  !!}

module Star_Formation_Rate_Krumholz2009_Utilities
  !!{  
  Implements various useful utility functions for calculations for the \cite{krumholz_star_2009} star formation surface
  density rate law.
  !!}
  implicit none
  private
  public :: krumholz2009MolecularFractionSlow

contains

  double precision function krumholz2009MolecularFractionSlow(s) result(fraction)
    !!{
    Slow (but more accurate at low molecular fraction) fitting function from \cite{krumholz_star_2009} for the molecular
    hydrogen fraction.
    !!}
    implicit none
    double precision, intent(in   ) :: s
    double precision, parameter     :: sTiny        =1.000000d-06
    double precision, parameter     :: sHuge        =1.000000d+10
    double precision, parameter     :: deltaInfinity=0.214008d+00 ! The value of δ for s → ∞.
    double precision, parameter     :: sMaximum     =1.000000d+01
    double precision                :: delta

    if      (s <  sTiny   ) then
       ! Series expansion for very small s.
       fraction=1.0d0-0.75d0*s
    else if (s >= sHuge   ) then
       ! Truncate to zero for extremely large s.
       fraction=0.0d0
    else if (s >= sMaximum) then
       ! Simplified form for very large s.
       fraction=1.0d0/(0.75d0/(1.0d0+deltaInfinity))**5/5.0d0/s**5
    else
       ! Full expression.
       delta   =0.0712d0/((0.1d0/s+0.675d0)**2.8d0)
       fraction=1.0d0-1.0d0/((1.0d0+(((1.0d0+delta)/0.75d0/s)**5))**0.2d0)
    end if
    return
  end function krumholz2009MolecularFractionSlow

end module Star_Formation_Rate_Krumholz2009_Utilities
