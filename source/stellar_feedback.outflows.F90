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
Contains a module which provides a class that implements outflows due to stellar feedback.
!!}

module Stellar_Feedback_Outflows
  !!{
  Provides a class that implements ejective stellar feedback.
  !!}
  use :: Galacticus_Nodes, only : nodeComponent
  private

  !![
  <functionClass>
   <name>stellarFeedbackOutflows</name>
   <descriptiveName>Stellar feedback.</descriptiveName>
   <description>
    Class providing models of outflows due to stellar feedback.
   </description>
   <default>powerLaw</default>
   <method name="outflowRate" >
    <description>Returns the outflow rates (both ejective and expulsive) due to stellar feedback in the given {\normalfont \ttfamily component} in units of $M_\odot/$Gyr.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class           (nodeComponent), intent(inout) :: component</argument>
    <argument>double precision               , intent(in   ) :: rateStarFormation  , rateEnergyInput</argument>
    <argument>double precision               , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive</argument>
   </method>
  </functionClass>
  !!]

end module Stellar_Feedback_Outflows
