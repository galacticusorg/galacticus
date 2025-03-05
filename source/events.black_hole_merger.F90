!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which performs tasks associated with black hole merger events.
!!}

module Events_Black_Hole_Merger
  !!{
  Performs tasks associated with black hole merger events.
  !!}
  implicit none
  private
  public :: Event_Black_Hole_Merger

contains

  subroutine Event_Black_Hole_Merger(blackHole1,blackHole2,blackHoleMerged)
    !!{
    Perform tasks associated with a merger beween {\normalfont \ttfamily blackHole1} and {\normalfont \ttfamily blackHole2}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class(nodeComponentBlackHole), intent(inout) :: blackHole1     , blackHole2, &
         &                                          blackHoleMerged

    !![
    <eventHook name="blackHoleMerger">
     <import>
      <module name="Galacticus_Nodes" symbols="nodeComponentBlackHole"/>
     </import>
     <interface>
       class(nodeComponentBlackHole), intent(inout) :: blackHole1, blackHole2, blackHoleMerged
     </interface>
     <callWith>blackHole1,blackHole2,blackHoleMerged</callWith>
    </eventHook>
    !!]
    return
  end subroutine Event_Black_Hole_Merger

end module Events_Black_Hole_Merger
