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
Implements a filter which passes only nodes that lie within a survey geometry.
!!}

  use :: Geometry_Surveys        , only : surveyGeometryClass
  use :: Node_Property_Extractors, only : nodePropertyExtractorPositionOrbital
  
  !![
  <enumeration>
   <name>positionType</name>
   <description>Enumeration of position types.</description>
   <visibility>public</visibility>
   <encodeFunction>yes</encodeFunction>
   <entry label="position"/>
   <entry label="orbital" />
  </enumeration>
  !!]

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
     class(surveyGeometryClass                 ), pointer :: surveyGeometry_  => null()
     type (enumerationPositionTypeType         )          :: positionType
     type (nodePropertyExtractorPositionOrbital)          :: positionOrbital_
   contains
     final     ::           surveyGeometryDestructor
     procedure :: passes => surveyGeometryPasses
  end type galacticFilterSurveyGeometry

  interface galacticFilterSurveyGeometry
     !!{
     Constructors for the \refClass{galacticFilterSurveyGeometry} galactic filter class.
     !!}
     module procedure surveyGeometryConstructorParameters
     module procedure surveyGeometryConstructorInternal
  end interface galacticFilterSurveyGeometry

contains

  function surveyGeometryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterSurveyGeometry} galactic filter class which takes a parameter set as input.
    !!}
    use :: Geometry_Surveys, only : surveyGeometry, surveyGeometryClass
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (galacticFilterSurveyGeometry)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(surveyGeometryClass         ), pointer       :: surveyGeometry_
    type (varying_string              )                :: positionType
    
    !![
    <inputParameter>
      <name>positionType</name>
      <source>parameters</source>
      <defaultValue>var_str('position')</defaultValue>
      <description>The type of position to use in survey geometry filters.</description>
    </inputParameter>
    <objectBuilder class="surveyGeometry" name="surveyGeometry_" source="parameters"/>
    !!]
    self=galacticFilterSurveyGeometry(enumerationPositionTypeEncode(positionType,includesPrefix=.false.),surveyGeometry_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="surveyGeometry_"/>
    !!]
    return
  end function surveyGeometryConstructorParameters

  function surveyGeometryConstructorInternal(positionType,surveyGeometry_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterSurveyGeometry} galactic filter class.
    !!}
    implicit none
    type (galacticFilterSurveyGeometry)                        :: self
    class(surveyGeometryClass         ), intent(in   ), target :: surveyGeometry_
    type (enumerationPositionTypeType ), intent(in   )         :: positionType
    !![
    <constructorAssign variables="positionType, *surveyGeometry_"/>
    !!]

    self%positionOrbital_=nodePropertyExtractorPositionOrbital()
    return
  end function surveyGeometryConstructorInternal

  subroutine surveyGeometryDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterSurveyGeometry} galactic filter class.
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
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk, nodeComponentPosition, nodeComponentSpheroid
    implicit none
    class           (galacticFilterSurveyGeometry), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    class           (nodeComponentBasic          ), pointer               :: basic
    class           (nodeComponentDisk           ), pointer               :: disk
    class           (nodeComponentSpheroid       ), pointer               :: spheroid
    class           (nodeComponentPosition       ), pointer               :: position
    double precision                              , dimension(3)          :: position_
    double precision                                                      :: massStellar

    disk        => node    %disk       ()
    spheroid    => node    %spheroid   ()
    massStellar = +disk    %massStellar() &
         &        +spheroid%massStellar()
    select case (self%positionType%ID)
    case (positionTypePosition%ID)
       position  => node                     %position(                 )
       position_ =  position                 %position(                 )
    case (positionTypeOrbital %ID)
       basic     => node                     %basic   (                 )
       position_ =  self    %positionOrbital_%extract (node,basic%time())
    end select
    surveyGeometryPasses=self%surveyGeometry_%pointIncluded(position_,massStellar)
    return
  end function surveyGeometryPasses
