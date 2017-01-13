!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% An implementation of the hot halo temperature class which uses the ``hydrostatic'' solution from the Enzo code.
  
  !# <hotHaloTemperatureProfile name="hotHaloTemperatureProfileEnzoHydrostatic">
  !#  <description>Provides an implementation of the hot halo temperature profile class which uses the ``hydrostatic'' solution from the Enzo code.</description>
  !# </hotHaloTemperatureProfile>
  type, extends(hotHaloTemperatureProfileClass) :: hotHaloTemperatureProfileEnzoHydrostatic
     !% An implementation of the hot halo temperature profile class which uses a ``hydrostatic'' profile as used by the Enzo code.
     private
   contains
     procedure :: temperature         => enzoHydrostaticTemperature
     procedure :: temperatureLogSlope => enzoHydrostaticTemperatureLogSlope
  end type hotHaloTemperatureProfileEnzoHydrostatic

  double precision, parameter :: enzoHydrostaticTemperatureMinimum=1.0d2
  
contains

  double precision function enzoHydrostaticTemperature(self,node,radius)
    !% Return the density in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Atomic
    use Numerical_Constants_Prefixes
    use Dark_Matter_Profiles
    implicit none
    class           (hotHaloTemperatureProfileEnzoHydrostatic), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: radius
    class           (darkMatterProfileClass                  ), pointer       :: darkMatterProfile_
    double precision                                                          :: enclosedMass
    !GCC$ attributes unused :: self
    
    if (radius == 0.0d0) then
       enzoHydrostaticTemperature=enzoHydrostaticTemperatureMinimum
    else
       darkMatterProfile_ => darkMatterProfile              (        &
            &                                               )
       enclosedMass       =  darkMatterProfile_%enclosedMass(        &
            &                                                node  , &
            &                                                radius  &
            &                                               )
       enzoHydrostaticTemperature=max(                                     &
            &                         +kilo                           **2  &
            &                         *gravitationalConstantGalacticus     &
            &                         *enclosedMass                        &
            &                         *meanAtomicMassPrimordial            &
            &                         *massHydrogenAtom                    &
            &                         /3.0d0                               &
            &                         /boltzmannsConstant                  &
            &                         /radius                            , &
            &                         enzoHydrostaticTemperatureMinimum    &
            &                        )
    end if
    return
  end function enzoHydrostaticTemperature
  
  double precision function enzoHydrostaticTemperatureLogSlope(self,node,radius)
    !% Return the logarithmic slope of the density profile in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass
    !% distribution.
    use Numerical_Constants_Math
    use Dark_Matter_Profiles
    implicit none
    class           (hotHaloTemperatureProfileEnzoHydrostatic), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: radius
    class           (darkMatterProfileClass                  ), pointer       :: darkMatterProfile_
    double precision                                                          :: enclosedMass      , density
    !GCC$ attributes unused :: self
    
    if (self%temperature(node,radius) <= enzoHydrostaticTemperatureMinimum) then
       enzoHydrostaticTemperatureLogSlope=0.0d0
    else
       darkMatterProfile_                => darkMatterProfile              (        &
            &                                                              )
       enclosedMass                      =  darkMatterProfile_%enclosedMass(        &
            &                                                               node  , &
            &                                                               radius  &
            &                                                              )
       density                           =  darkMatterProfile_%density     (        &
            &                                                               node  , &
            &                                                               radius  &
            &                                                              )
       enzoHydrostaticTemperatureLogSlope= +4.0d0           &
            &                               *Pi              &
            &                               *radius      **3 &
            &                               *density         &
            &                               /enclosedMass    &
            &                               -1.0d0
    end if
    return
  end function enzoHydrostaticTemperatureLogSlope
  
