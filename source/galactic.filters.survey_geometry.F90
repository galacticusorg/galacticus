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

!!{
Contains a module which implements a filter which passes only nodes that lie within a survey geometry.
!!}

  use :: Geometry_Surveys, only : surveyGeometryClass

  !![
  <galacticFilter name="galacticFilterSurveyGeometry">
   <description>A filter which passes only nodes that lie within a survey geometry.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterSurveyGeometry
     !!{
     A galactic filter class which passes only nodes that lie within a survey geometry.
     !!}
     private
     class(surveyGeometryClass), pointer :: surveyGeometry_ => null()
   contains
     final     ::           surveyGeometryDestructor
     procedure :: passes => surveyGeometryPasses
  end type galacticFilterSurveyGeometry

  interface galacticFilterSurveyGeometry
     !!{
     Constructors for the ``surveyGeometry'' galactic filter class.
     !!}
     module procedure surveyGeometryConstructorParameters
     module procedure surveyGeometryConstructorInternal
  end interface galacticFilterSurveyGeometry

contains

  function surveyGeometryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``surveyGeometry'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Geometry_Surveys, only : surveyGeometry, surveyGeometryClass
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (galacticFilterSurveyGeometry)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(surveyGeometryClass         ), pointer       :: surveyGeometry_

    !![
    <objectBuilder class="surveyGeometry" name="surveyGeometry_" source="parameters"/>
    !!]
    self=galacticFilterSurveyGeometry(surveyGeometry_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="surveyGeometry_"/>
    !!]
    return
  end function surveyGeometryConstructorParameters

  function surveyGeometryConstructorInternal(surveyGeometry_) result(self)
    !!{
    Internal constructor for the ``surveyGeometry'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterSurveyGeometry)                        :: self
    class(surveyGeometryClass        ), intent(in   ), target :: surveyGeometry_
    !![
    <constructorAssign variables="*surveyGeometry_"/>
    !!]

    return
  end function surveyGeometryConstructorInternal

  subroutine surveyGeometryDestructor(self)
    !!{
    Destructor for the ``surveyGeometry'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterSurveyGeometry), intent(inout) :: self

    !![
    <objectDestructor name="self%surveyGeometry_"/>
    !!]
    return
  end subroutine surveyGeometryDestructor

  logical function surveyGeometryPasses(self,node)
    !!{
    Implement a galactic filter which passes only nodes with a survey geometry.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentPosition, nodeComponentSpheroid, treeNode
    implicit none
    class           (galacticFilterSurveyGeometry), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    class           (nodeComponentDisk           ), pointer               :: disk
    class           (nodeComponentSpheroid       ), pointer               :: spheroid
    class           (nodeComponentPosition       ), pointer               :: position
    double precision                                                      :: massStellar

    position    => node    %position   ()
    disk        => node    %disk       ()
    spheroid    => node    %spheroid   ()
    massStellar = +disk    %massStellar() &
         &        +spheroid%massStellar()
    surveyGeometryPasses=self%surveyGeometry_%pointIncluded(position%position(),massStellar)
    return
  end function surveyGeometryPasses
