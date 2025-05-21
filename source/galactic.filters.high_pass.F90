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
Implements a high-pass filter on any node property.
!!}

  use :: Node_Property_Extractors, only : nodePropertyExtractorScalar

  !![
  <galacticFilter name="galacticFilterHighPass">
   <description>A high-pass filter on any node property.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterHighPass
     !!{
     A high-pass galactic filter class on any node property.
     !!}
     private
     class           (nodePropertyExtractorScalar), pointer :: nodePropertyExtractor_ => null()
     double precision                                       :: threshold
   contains
     final     ::           highPassDestructor
     procedure :: passes => highPassPasses
  end type galacticFilterHighPass

  interface galacticFilterHighPass
     !!{
     Constructors for the \refClass{galacticFilterHighPass} galactic filter class.
     !!}
     module procedure highPassConstructorParameters
     module procedure highPassConstructorInternal
  end interface galacticFilterHighPass

contains
  
  function highPassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterHighPass} galactic filter class which takes a parameter set as input.
    !!}
    use :: Error                   , only : Error_Report
    use :: Input_Parameters        , only : inputParameter            , inputParameters
    use :: Node_Property_Extractors, only : nodePropertyExtractorClass
    implicit none
    type            (galacticFilterHighPass    )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (nodePropertyExtractorClass), pointer       :: nodePropertyExtractor_
    double precision                                            :: threshold

    !![
    <inputParameter>
      <name>threshold</name>
      <source>parameters</source>
      <description>The threshold value above which to pass.</description>
    </inputParameter>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       self=galacticFilterHighPass(threshold,nodePropertyExtractor_)
     class default
       call Error_Report('extracted property must be a real scalar'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function highPassConstructorParameters

  function highPassConstructorInternal(threshold,nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterHighPass} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterHighPass     )                        :: self
    class           (nodePropertyExtractorScalar), intent(in   ), target :: nodePropertyExtractor_
    double precision                             , intent(in   )         :: threshold
    !![
    <constructorAssign variables="threshold, *nodePropertyExtractor_"/>
    !!]

    return
  end function highPassConstructorInternal

  subroutine highPassDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterHighPass} galactic filter class.
    !!}
    implicit none
    type(galacticFilterHighPass), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine highPassDestructor

  logical function highPassPasses(self,node)
    !!{
    Implement a high-pass galactic filter.
    !!}
    implicit none
    class(galacticFilterHighPass), intent(inout)         :: self
    type (treeNode              ), intent(inout), target :: node
   
    highPassPasses=self%nodePropertyExtractor_%extract(node) >= self%threshold
    return
  end function highPassPasses
