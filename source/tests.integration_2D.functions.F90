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

!+    Contributions to this file made by: Niusha Ahvazi

!!{RST
Contains a module that provides simple analytic test integrand functions used by the two-dimensional integration unit tests.
!!}

module Test_Integration2D_Functions
  !!{RST
  Provides simple analytic test integrand functions used by the two-dimensional integration unit tests.
  !!}
  implicit none
  public

contains
  
  double precision function integrand1(x, y)
    !!{RST
    Integrand function used in tests.
    !!}
    implicit none
    double precision, intent(in) :: x, y
    
    integrand1=sin(x**2)*cos(y)
    return
  end function integrand1
  
  double precision function integrand2(x, y)
    !!{RST
    Integrand function used in tests.
    !!}
    implicit none
    double precision, intent(in) :: x, y
    
    integrand2=x**2*cos(y)
    return
  end function integrand2
  
  double precision function integrand3(x, y)
    !!{RST
    Integrand function used in tests.
    !!}
    implicit none
    double precision, intent(in) :: x, y
    
    integrand3=sqrt(x)*y**2
    return
  end function integrand3
  
  double precision function integrand4(x, y)
    !!{RST
    Integrand function used in tests.
    !!}
    implicit none
    double precision, intent(in) :: x, y
    
    integrand4=y/sqrt(x)
    return
  end function integrand4
  
end module Test_Integration2D_Functions

