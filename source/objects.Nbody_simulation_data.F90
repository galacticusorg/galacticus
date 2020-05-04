!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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
  use :: IO_HDF5           , only : hdf5Object
  use :: ISO_Varying_String, only : varying_string
  use :: Kind_Numbers      , only : kind_int8
  use :: Hashes            , only : rank1IntegerSizeTHash, rank2IntegerSizeTHash, rank1DoubleHash, rank2DoubleHash
  implicit none
  private
  public :: nBodyData, nBodyDataPropertyType

  type :: nBodyData
     !% A class to store N-body simulation data.
     type            (varying_string       )                              :: label
     type            (hdf5Object           )                              :: analysis
     double precision                       , allocatable, dimension(:,:) :: position              , velocity
     integer                                , allocatable, dimension(  :) :: identifier
     double precision                                                     :: lengthSoftening       , massParticle
     integer         (kind_int8            ), allocatable, dimension(  :) :: particleIDs
     type            (rank1IntegerSizeTHash)                              :: propertiesInteger
     type            (rank1DoubleHash      )                              :: propertiesReal
     type            (rank2IntegerSizeTHash)                              :: propertiesIntegerRank1
     type            (rank2DoubleHash      )                              :: propertiesRealRank1
  end type nBodyData

  interface nBodyData
     module procedure nBodyDataConstructor
  end interface nBodyData

  !# <enumeration>
  !#  <name>propertyType</name>
  !#  <description>Enumeration of property types for N-body data properties.</description>
  !#  <visibility>public</visibility>
  !#  <entry label="unknown"/>
  !#  <entry label="integer"/>
  !#  <entry label="real"   />
  !# </enumeration>
  
contains

  function nBodyDataConstructor() result (self)
    !% A default constructor for the {\normalfont \ttfamily nBodyData} class.
    implicit none
    type(nBodyData) :: self

    return
  end function nBodyDataConstructor

  integer function nBodyDataPropertyType(propertyName)
    !% Returns the type of the named property.
    implicit none
    character(len=*), intent(in   ) :: propertyName

    select case (propertyName)
    case('particleID'               )
       nBodyDataPropertyType=propertyTypeInteger
    case('descendentID'             )
       nBodyDataPropertyType=propertyTypeInteger
    case('progenitorCount'          )
       nBodyDataPropertyType=propertyTypeInteger
    case('hostID'                   )
       nBodyDataPropertyType=propertyTypeInteger
    case('hostRootID'               )
       nBodyDataPropertyType=propertyTypeInteger
    case('descendentHostID'         )
       nBodyDataPropertyType=propertyTypeInteger
    case('isPhantom'                )
       nBodyDataPropertyType=propertyTypeInteger
    case('alwaysIsolated'           )
       nBodyDataPropertyType=propertyTypeInteger
    case('expansionFactor'          )
       nBodyDataPropertyType=propertyTypeReal
    case('descendentExpansionFactor')
       nBodyDataPropertyType=propertyTypeReal
    case('massVirial'               )
       nBodyDataPropertyType=propertyTypeReal
    case default
       nBodyDataPropertyType=propertyTypeUnknown
    end select
    return
  end function nBodyDataPropertyType

end module NBody_Simulation_Data
