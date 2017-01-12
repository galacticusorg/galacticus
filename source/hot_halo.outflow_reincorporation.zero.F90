!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% An implementation of the hot halo outflow reincorporation class which gives zero reincorporation rate.
  
  !# <hotHaloOutflowReincorporation name="hotHaloOutflowReincorporationZero">
  !#  <description>An implementation of the hot halo outflow reincorporation class which gives zero reincorporation rate.</description>
  !# </hotHaloOutflowReincorporation>
  type, extends(hotHaloOutflowReincorporationClass) :: hotHaloOutflowReincorporationZero
     !% An implementation of the hot halo outflow reincorporation class which gives zero reincorporation rate.
     private
   contains
     procedure :: rate => zeroRate
  end type hotHaloOutflowReincorporationZero

contains

  double precision function zeroRate(self,node)
    !% Return the rate of mass reincorporation for outflowed gas in the hot halo.
    use Cosmology_Functions
    use Dark_Matter_Profiles
    implicit none
    class(hotHaloOutflowReincorporationZero), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    !GCC$ attributes unused :: self, node
    
    zeroRate=0.0d0
    return
  end function zeroRate
