!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of error functions.

module Error_Functions
  !% Implements calculations of error functions.
  implicit none
  private
  public :: Error_Function,Error_Function_Complementary

contains

  double precision function Error_Function(argument)
    !% Computes the error function.
    use FGSL
    implicit none
    double precision, intent(in   ) :: argument

    Error_Function=FGSL_SF_Erf(argument)
    return
  end function Error_Function

  double precision function Error_Function_Complementary(argument)
    !% Computes the complementary error function.
    use FGSL
    implicit none
    double precision, intent(in   ) :: argument

    Error_Function_Complementary=FGSL_SF_ErfC(argument)
    return
  end function Error_Function_Complementary

end module Error_Functions
