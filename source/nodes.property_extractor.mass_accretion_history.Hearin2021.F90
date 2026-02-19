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
Implements a node property extractor class for parameters of the \cite{hearin_differentiable_2021} mass accretion history model.
!!}!

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassAccretionHistoryHearin2021">
   <description>A node property extractor class for parameters of the \cite{hearin_differentiable_2021} mass accretion history model.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorMassAccretionHistoryHearin2021
     !!{
     A property extractor class for parameters of the \cite{hearin_differentiable_2021} mass accretion history model.
     !!}
     private
   contains
     procedure :: elementCount => massAccretionHistoryHearin2021ElementCount
     procedure :: extract      => massAccretionHistoryHearin2021Extract
     procedure :: names        => massAccretionHistoryHearin2021Names
     procedure :: descriptions => massAccretionHistoryHearin2021Descriptions
     procedure :: unitsInSI    => massAccretionHistoryHearin2021UnitsInSI
  end type nodePropertyExtractorMassAccretionHistoryHearin2021

  interface nodePropertyExtractorMassAccretionHistoryHearin2021
     !!{
     Constructors for the \refClass{nodePropertyExtractorMassAccretionHistoryHearin2021} output analysis class.
     !!}
     module procedure massAccretionHistoryHearin2021ConstructorParameters
  end interface nodePropertyExtractorMassAccretionHistoryHearin2021

contains

  function massAccretionHistoryHearin2021ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMassAccretionHistoryHearin2021} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorMassAccretionHistoryHearin2021)                :: self
    type(inputParameters                                    ), intent(inout) :: parameters

    self=nodePropertyExtractorMassAccretionHistoryHearin2021()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massAccretionHistoryHearin2021ConstructorParameters

  integer function massAccretionHistoryHearin2021ElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily massAccretionHistoryHearin2021} property extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorMassAccretionHistoryHearin2021), intent(inout) :: self
    double precision                                                     , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    massAccretionHistoryHearin2021ElementCount=3
    return
  end function massAccretionHistoryHearin2021ElementCount

  function massAccretionHistoryHearin2021Extract(self,node,time,instance)
    !!{
    Implement extraction of halo environment properties.
    !!}
    implicit none
    double precision                                                     , dimension(:) , allocatable :: massAccretionHistoryHearin2021Extract
    class           (nodePropertyExtractorMassAccretionHistoryHearin2021), intent(inout), target      :: self
    type            (treeNode                                           ), intent(inout), target      :: node
    double precision                                                     , intent(in   )              :: time
    type            (multiCounter                                       ), intent(inout), optional    :: instance
    double precision                                                                                  :: powerLawIndexEarly                   , powerLawIndexLate, &
         &                                                                                               log10TimeZero
    !$GLC attributes unused :: time, instance

    if (node%hostTree%properties%exists('treeMAHHearin2021PowerLawIndexEarly')) then
       powerLawIndexEarly=node%hostTree%properties%value('treeMAHHearin2021PowerLawIndexEarly')
    else
       ! Set an unphysical value if the property does not exist.
       powerLawIndexEarly=-huge(0.0d0)
    end if
    if (node%hostTree%properties%exists('treeMAHHearin2021PowerLawIndexLate' )) then
       powerLawIndexLate =node%hostTree%properties%value('treeMAHHearin2021PowerLawIndexLate' )
    else
       ! Set an unphysical value if the property does not exist.
       powerLawIndexLate =-huge(0.0d0)
    end if
    if (node%hostTree%properties%exists('treeMAHHearin2021Log10TimeZero'     )) then
       log10TimeZero     =node%hostTree%properties%value('treeMAHHearin2021Log10TimeZero'     )
    else
       ! Set an unphysical value if the property does not exist.
       log10TimeZero     =-huge(0.0d0)
    end if
    allocate(massAccretionHistoryHearin2021Extract(3))
    massAccretionHistoryHearin2021Extract=[                    &
         &                                 powerLawIndexEarly, &
         &                                 powerLawIndexLate , &
         &                                 log10TimeZero       &
         &                                ]
    return
  end function massAccretionHistoryHearin2021Extract

  subroutine massAccretionHistoryHearin2021Names(self,time,names)
    !!{
    Return the name of the {\normalfont \ttfamily massAccretionHistoryHearin2021} property.
    !!}
    implicit none
    class           (nodePropertyExtractorMassAccretionHistoryHearin2021), intent(inout)                             :: self
    double precision                                                     , intent(in   )                             :: time
    type            (varying_string                                     ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)='treeMAHHearin2021PowerLawIndexEarly'
    names(2)='treeMAHHearin2021PowerLawIndexLate'
    names(3)='treeMAHHearin2021Log10TimeZero'
    return
  end subroutine massAccretionHistoryHearin2021Names

  subroutine massAccretionHistoryHearin2021Descriptions(self,time,descriptions)
    !!{
    Return a description of the {\normalfont \ttfamily massAccretionHistoryHearin2021} property.
    !!}
    implicit none
    class           (nodePropertyExtractorMassAccretionHistoryHearin2021), intent(inout)                             :: self
    double precision                                                     , intent(in   )                             :: time
    type            (varying_string                                     ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)='Early-time power-law index in the Hearin et al. (2021) mass accretion history for this tree.'
    descriptions(2)='Late-time power-law index in the Hearin et al. (2021) mass accretion history for this tree.'
    descriptions(3)='Parameter log₁₀t₀ in the Hearin et al. (2021) mass accretion history for this tree.'
    return
  end subroutine massAccretionHistoryHearin2021Descriptions

  function massAccretionHistoryHearin2021UnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily massAccretionHistoryHearin2021} property in the SI system.
    !!}
    implicit none
    double precision                                                     , allocatable  , dimension(:) :: massAccretionHistoryHearin2021UnitsInSI
    class           (nodePropertyExtractorMassAccretionHistoryHearin2021), intent(inout)               :: self
    double precision                                                     , intent(in   )               :: time
    !$GLC attributes unused :: self, time

    allocate(massAccretionHistoryHearin2021UnitsInSI(3))
    massAccretionHistoryHearin2021UnitsInSI=0.0d0
    return
  end function massAccretionHistoryHearin2021UnitsInSI

