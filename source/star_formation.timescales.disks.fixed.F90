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

  !% Implementation of a fixed timescale for star formation in galactic disks.
  
  !# <starFormationTimescaleDisks name="starFormationTimescaleDisksFixed">
  !#  <description>A fixed timescale for star formation in galactic disks.</description>
  !# </starFormationTimescaleDisks>
  type, extends(starFormationTimescaleDisksClass) :: starFormationTimescaleDisksFixed
     !% Implementation of a fixed timescale for star formation in galactic disks.
     private
     double precision :: timescaleValue
   contains
     procedure :: timescale => fixedTimescale
  end type starFormationTimescaleDisksFixed

  interface starFormationTimescaleDisksFixed
     !% Constructors for the {\normalfont \ttfamily fixed} timescale for star formation in disks class.
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface starFormationTimescaleDisksFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily fixed} timescale for star formation in disks class which takes a
    !% parameter set as input.
    use Input_Parameters
    implicit none
    type            (starFormationTimescaleDisksFixed)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    double precision                                                  :: timescale

    !# <inputParameter>
    !#   <name>timescale</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The timescale for star formation in the fixed timescale model for disks.</description>
    !#   <group>starFormation</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=starFormationTimescaleDisksFixed(timescale)
    !# <inputParametersValidate source="parameters"/>
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(timescale) result(self)
    !% Internal constructor for the {\normalfont \ttfamily fixed} timescale for star formation in disks class.
    implicit none
    type            (starFormationTimescaleDisksFixed)                :: self
    double precision                                  , intent(in   ) :: timescale

    self%timescaleValue=timescale
    return
  end function fixedConstructorInternal

  double precision function fixedTimescale(self,node)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\normalfont \ttfamily node}, assuming a fixed
    !% timecale.
    implicit none
    class(starFormationTimescaleDisksFixed), intent(inout), target :: self
    type (treeNode                        ), intent(inout), target :: node
    !GCC$ attributes unused :: node

    fixedTimescale=self%timescaleValue
    return
  end function fixedTimescale
