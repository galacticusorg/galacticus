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
Contains a module which stores data for the standard black hole node component.
!!}

module Node_Component_Black_Hole_Standard_Data
  !!{
  Stores data for the standard black hole node component.
  !!}
  use :: Object_Pools, only : objectPool
  implicit none
  public

  ! A per-thread pool of black hole mass distributions, re-used across nodes to avoid allocating and
  ! destroying a new object on every call to the component's mass distribution function. The pool is
  ! released in the component's thread uninitialization task.
  type (objectPool) :: massDistributionPool
  !$omp threadprivate(massDistributionPool)

end module Node_Component_Black_Hole_Standard_Data
