!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% An implementation of the hot halo temperature class which uses an isothermal virial temperature.
  
  !# <hotHaloTemperatureProfile name="hotHaloTemperatureProfileVirial">
  !#  <description>Provides an implementation of the hot halo temperature profile class which uses an isothermal virial temperature.</description>
  !# </hotHaloTemperatureProfile>
  type, extends(hotHaloTemperatureProfileClass) :: hotHaloTemperatureProfileVirial
     !% An implementation of the hot halo temperature profile class which uses an isothermal virial temperature.
     private
   contains
     procedure :: temperature         => virialTemperature
     procedure :: temperatureLogSlope => virialTemperatureLogSlope
  end type hotHaloTemperatureProfileVirial

  interface hotHaloTemperatureProfileVirial
     !% Constructors for the {\normalfont \ttfamily virial} hot halo temperature profile class.
     module procedure virialDefaultConstructor
  end interface hotHaloTemperatureProfileVirial
  
contains

  function virialDefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily virial} hot halo temperature profile class.
    implicit none
    type(hotHaloTemperatureProfileVirial) :: virialDefaultConstructor

    return
  end function virialDefaultConstructor
  
  double precision function virialTemperature(self,node,radius)
    !% Return the density in a {\normalfont \ttfamily virial} hot halo mass distribution.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (hotHaloTemperatureProfileVirial), intent(inout)          :: self
    type            (treeNode                       ), intent(inout), pointer :: node
    double precision                                 , intent(in   )          :: radius
    class           (darkMatterHaloScaleClass       )               , pointer :: darkMatterHaloScale_ 
    
    darkMatterHaloScale_ => darkMatterHaloScale                   (    )
    virialTemperature    =  darkMatterHaloScale_%virialTemperature(node)
    return
  end function virialTemperature
  
  double precision function virialTemperatureLogSlope(self,node,radius)
    !% Return the logarithmic slope of the density profile in a {\normalfont \ttfamily virial} hot halo mass
    !% distribution.
    implicit none
    class           (hotHaloTemperatureProfileVirial), intent(inout)          :: self
    type            (treeNode                       ), intent(inout), pointer :: node
    double precision                                 , intent(in   )          :: radius

    virialTemperatureLogSlope=0.0d0
    return
  end function virialTemperatureLogSlope
  
