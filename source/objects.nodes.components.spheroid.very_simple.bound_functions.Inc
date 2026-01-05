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
Contains custom functions for the very simple spheroid component.
!!}

double precision function Node_Component_Spheroid_Very_Simple_Half_Mass_Radius(self)
  !!{
  Return the half-mass radius of the very simple spheroid.
  !!}
  implicit none
  class(nodeComponentSpheroidVerySimple), intent(inout) :: self

  Node_Component_Spheroid_Very_Simple_Half_Mass_Radius=self%radius()
  return
end function Node_Component_Spheroid_Very_Simple_Half_Mass_Radius

double precision function Node_Component_Spheroid_Very_Simple_Mass_Baryonic(self) result(massBaryonic)
  !!{
  Return the baryonic mass for the very simple spheroid component.
  !!}
  implicit none
  class(nodeComponentSpheroidVerySimple), intent(inout) :: self

  massBaryonic=+max(0.0d0,self%massStellar()) &
       &       +max(0.0d0,self%massGas    ())
  return
end function Node_Component_Spheroid_Very_Simple_Mass_Baryonic
