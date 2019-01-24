!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Contains a module which provides a class to store N-body simulation data.

module NBody_Simulation_Data
  !% Provides a class to store N-body simulation data.
  use IO_HDF5
  use Kind_Numbers
  implicit none
  private
  public :: nBodyData
  
  type :: nBodyData
     !% A class to store N-body simulation data.
     type            (hdf5Object)                              :: analysis
     double precision            , allocatable, dimension(:,:) :: position           , velocity
     integer                     , allocatable, dimension(  :) :: identifier
     double precision                                          :: lengthSoftening    , massParticle
     integer(kind=kind_int8)     , allocatable, dimension(  :) :: particleIDs        , particleIDsPrevious
     integer                     , allocatable, dimension(:,:) :: boundStatusPrevious, sampleWeightPrevious
  end type nBodyData

  interface nBodyData
     module procedure nBodyDataConstructor
  end interface nBodyData

contains

  function nBodyDataConstructor() result (self)
    !% A default constructor for the {\normalfont \ttfamily nBodyData} class.
    implicit none
    type(nBodyData) :: self

    return
  end function nBodyDataConstructor

end module NBody_Simulation_Data
