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
Contains a module which provides a class that implements ram pressure stripping.
!!}

module Ram_Pressure_Stripping_Mass_Loss_Rate
  !!{
  Provides a class that implements calculations of ram pressure stripping.
  !!}
  use :: Galacticus_Nodes, only : nodeComponent
  private

  !![
  <functionClass>
   <name>ramPressureStripping</name>
   <descriptiveName>Ram pressure stripping</descriptiveName>
   <description>
    Class providing models of ram pressure stripping-induced rates of mass loss.
   </description>
   <default>simpleCylindrical</default>
   <method name="rateMassLoss" >
    <description>Returns the rate of mass loss (in $M_\odot$~Gyr$^{-1}$) due to ram pressure stripping of the given {\normalfont \ttfamily component}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(nodeComponent), intent(inout) :: component</argument>
   </method>
  </functionClass>
  !!]

end module Ram_Pressure_Stripping_Mass_Loss_Rate
