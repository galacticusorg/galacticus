!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains custom functions for the very simple satellite component.
!!}

double precision function Node_Component_Satellite_Very_Simple_Time_Until_Merging(self)
  !!{
  Return the {\normalfont \ttfamily timeUntilMerging} property of the {\normalfont \ttfamily satelliteVerySimple} component.
  !!}
  implicit none
  class(nodeComponentSatelliteVerySimple), intent(inout) :: self
  type (treeNode                        ), pointer       :: selfNode

  selfNode => self%host()
  if (selfNode%isSatellite()) then
     Node_Component_Satellite_Very_Simple_Time_Until_Merging=self%timeUntilMergingData
  else
     Node_Component_Satellite_Very_Simple_Time_Until_Merging=-1.0d0
  end if
  return
end function Node_Component_Satellite_Very_Simple_Time_Until_Merging

double precision function Node_Component_Satellite_Very_Simple_Time_Of_Merging(self)
  !!{
  Return the {\normalfont \ttfamily timeOfMerging} property of the {\normalfont \ttfamily satelliteVerySimple} component.
  !!}
  implicit none
  class(nodeComponentSatelliteVerysimple), intent(inout) :: self
  type (treeNode                        ), pointer       :: selfNode
  class(nodeComponentBasic              ), pointer       :: selfBasicComponent


  selfNode => self%host()
  if (selfNode%isSatellite()) then
     selfBasicComponent => selfNode%basic()
     Node_Component_Satellite_Very_Simple_Time_Of_Merging=self%timeUntilMergingData+selfBasicComponent%time()
  else
     Node_Component_Satellite_Very_Simple_Time_Of_Merging=-1.0d0
  end if
  return
end function Node_Component_Satellite_Very_Simple_Time_Of_Merging
