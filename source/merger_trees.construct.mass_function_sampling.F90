!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a class for sampling the halo mass function when constructing merger trees.

module Merger_Trees_Mass_Function_Sampling
  !% Implements a class for sampling the halo mass function when constructing merger trees.
  private
  
  !# <functionClass>
  !#  <name>mergerTreeHaloMassFunctionSampling</name>
  !#  <descriptiveName>Merger Tree Halo Mass Function Sampling</descriptiveName>
  !#  <description>Class providing methods for sampling from the halo mass function when building merger trees.</description>
  !#  <default>haloMassFunction</default>
  !#  <method name="sample" >
  !#   <description>Returns the sampling rate for merger trees of the given {\normalfont \ttfamily mass}, per decade of halo mass.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: mass,time,massMinimum,massMaximum</argument>
  !#  </method>
  !# </functionClass>

end module Merger_Trees_Mass_Function_Sampling
