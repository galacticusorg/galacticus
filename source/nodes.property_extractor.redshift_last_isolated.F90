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
Implements a property extractor class that extracts the redshift at which a \gls{node} was last isolated.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRedshiftLastIsolated">
   <description>A node property extractor class which extracts the redshift at which a \gls{node} was last isolated---named ``{\normalfont \ttfamily redshiftLastIsolated}''.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRedshiftLastIsolated
     !!{
     A redshiftLastIsolated property extractor class.
     !!}
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
   contains
     final     ::                redshiftLastIsolatedDestructor
     procedure :: extract     => redshiftLastIsolatedExtract
     procedure :: name        => redshiftLastIsolatedName
     procedure :: description => redshiftLastIsolatedDescription
     procedure :: unitsInSI   => redshiftLastIsolatedUnitsInSI
  end type nodePropertyExtractorRedshiftLastIsolated

  interface nodePropertyExtractorRedshiftLastIsolated
     !!{
     Constructors for the {\normalfont \ttfamily redshiftLastIsolated} output analysis class.
     !!}
     module procedure redshiftLastIsolatedConstructorParameters
     module procedure redshiftLastIsolatedConstructorInternal
  end interface nodePropertyExtractorRedshiftLastIsolated

contains

  function redshiftLastIsolatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily redshiftLastIsolated} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRedshiftLastIsolated)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRedshiftLastIsolated(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function redshiftLastIsolatedConstructorParameters

  function redshiftLastIsolatedConstructorInternal(cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily redshiftLastIsolated} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRedshiftLastIsolated)                        :: self
    class(cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    return
  end function redshiftLastIsolatedConstructorInternal

  subroutine redshiftLastIsolatedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily redshiftLastIsolated} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRedshiftLastIsolated), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine redshiftLastIsolatedDestructor

  double precision function redshiftLastIsolatedExtract(self,node,instance)
    !!{
    Implement a last isolated redshift output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorRedshiftLastIsolated), intent(inout), target   :: self
    type (treeNode                                 ), intent(inout), target   :: node
    type (multiCounter                             ), intent(inout), optional :: instance
    class(nodeComponentBasic                       ), pointer                 :: basic
    !$GLC attributes unused :: self, instance

    basic                       => node              %basic()
    redshiftLastIsolatedExtract =  self %cosmologyFunctions_%redshiftFromExpansionFactor(    &
         &                         self %cosmologyFunctions_%expansionFactor             (   &
         &                         basic                    %timeLastIsolated             () &
         &                                                                               )   &
         &                                                                              )
    return
  end function redshiftLastIsolatedExtract

  function redshiftLastIsolatedName(self)
    !!{
    Return the name of the last isolated redshift property.
    !!}
    implicit none
    type (varying_string                           )                :: redshiftLastIsolatedName
    class(nodePropertyExtractorRedshiftLastIsolated), intent(inout) :: self
    !$GLC attributes unused :: self

    redshiftLastIsolatedName=var_str('redshiftLastIsolated')
    return
  end function redshiftLastIsolatedName

  function redshiftLastIsolatedDescription(self)
    !!{
    Return a description of the redshiftLastIsolated property.
    !!}
    implicit none
    type (varying_string                           )                :: redshiftLastIsolatedDescription
    class(nodePropertyExtractorRedshiftLastIsolated), intent(inout) :: self
    !$GLC attributes unused :: self

    redshiftLastIsolatedDescription=var_str('The redshift of the epoch at which the galaxy was last isolated.')
    return
  end function redshiftLastIsolatedDescription

  double precision function redshiftLastIsolatedUnitsInSI(self)
    !!{
    Return the units of the last isolated redshift property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorRedshiftLastIsolated), intent(inout) :: self
    !$GLC attributes unused :: self

    redshiftLastIsolatedUnitsInSI=0.0d0
    return
  end function redshiftLastIsolatedUnitsInSI

