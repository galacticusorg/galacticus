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


!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which records whether the disk component truncates its mass distribution.
!!}

module Galactic_Component_Truncation_Disk_Status
  !!{RST
  Records whether the disk component truncates its mass distribution (i.e. whether the
  :galacticus-class:`galacticComponentTruncationClass` object that it uses is active). This allows classes which must use the
  same truncation as the disk (and which cannot depend on the disk component) to check that they do.
  !!}
  implicit none
  private
  public :: diskTruncationIsActive

  ! Set by the initialization task of the disk component.
  logical :: diskTruncationIsActive=.false.

end module Galactic_Component_Truncation_Disk_Status
