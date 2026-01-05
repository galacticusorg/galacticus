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
Contains a module which implements a dark matter halo bias class.
!!}

module Dark_Matter_Halo_Biases
  !!{
  Implements a dark matter halo bias class.
  !!}
  use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>darkMatterHaloBias</name>
   <descriptiveName>Dark matter halo biases.</descriptiveName>
   <description>
    Class providing models of the bias of dark matter halos.
   </description>
   <default>tinker2010</default>
   <method name="biasByMass" >
    <description>Returns the bias of a halo specified by a mass (in $M_\odot$) and time (in Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: mass, time</argument>
    <argument>double precision, intent(in   ), optional :: radius</argument>
   </method>
   <method name="biasByNode" >
    <description>Returns the bias of the halo in the supplied \gls{node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)           :: node</argument>
    <argument>double precision          , intent(in   ), optional :: radius</argument>
    <code>
     class(nodeComponentBasic), pointer :: basic
     basic                        => node%basic     (                                )
     darkMatterHaloBiasBiasByNode =  self%biasByMass(basic%mass(),basic%time(),radius)
    </code>
   </method>
   <generic name="bias">
    <method>biasByMass</method>
    <method>biasByNode</method>
   </generic>
  </functionClass>
  !!]

end module Dark_Matter_Halo_Biases
