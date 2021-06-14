!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains custom functions for the Cole2000 formation times component.

double precision function FormationTimeCole2000FormationTime(self)
  !% Return the {\normalfont \ttfamily formationTime} property of the {\normalfont \ttfamily formationTime} component class.
  implicit none
  class(nodeComponentFormationTimeCole2000), intent(inout) :: self
  type (treeNode                          ), pointer       :: formationNode         , selfNode
  class(nodeComponentFormationTime        ), pointer       :: formationTimeComponent

  selfNode                          => self                  %host         ()
  if (associated(selfNode%formationNode)) then
     formationNode                     => selfNode              %formationNode
  else
     formationNode => selfNode
  end if
  formationTimeComponent            => formationNode         %formationTime()
  select type (formationTimeComponent)
  class is (nodeComponentFormationTimeCole2000)
     FormationTimeCole2000FormationTime=  formationTimeComponent%formationTimeData
  class default
     FormationTimeCole2000FormationTime=0.0d0
     call Galacticus_Error_Report('incorrect class'//{introspection:location})
  end select
  return
end function FormationTimeCole2000FormationTime
