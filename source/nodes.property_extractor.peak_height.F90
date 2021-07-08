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

  use :: Cosmological_Density_Field, only : cosmologicalMassVariance, cosmologicalMassVarianceClass, criticalOverdensity, criticalOverdensityClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorPeakHeight">
   <description>A property extractor class for peak height of the halo.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorPeakHeight
     !!{
     A property extractor which extracts peakHeight properties.
     !!}
     private
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
   contains
     final     ::                 peakHeightDestructor
     procedure :: elementCount => peakHeightElementCount
     procedure :: extract      => peakHeightExtract
     procedure :: names        => peakHeightNames
     procedure :: descriptions => peakHeightDescriptions
     procedure :: unitsInSI    => peakHeightUnitsInSI
     procedure :: type         => peakHeightType
  end type nodePropertyExtractorPeakHeight

  interface nodePropertyExtractorPeakHeight
     !!{
     Constructors for the ``peakHeight'' output extractor class.
     !!}
     module procedure peakHeightConstructorParameters
     module procedure peakHeightConstructorInternal
  end interface nodePropertyExtractorPeakHeight

contains

  function peakHeightConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``peakHeight'' property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorPeakHeight)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(criticalOverdensityClass       ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=nodePropertyExtractorPeakHeight(criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function peakHeightConstructorParameters

  function peakHeightConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the ``peakHeight'' output extractor property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorPeakHeight)                        :: self
    class(criticalOverdensityClass       ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="*criticalOverdensity_,*cosmologicalMassVariance_"/>
    !!]

    return
  end function peakHeightConstructorInternal

  subroutine peakHeightDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily peakHeight} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorPeakHeight), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine peakHeightDestructor

  integer function peakHeightElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily peakHeight} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorPeakHeight), intent(inout) :: self
    double precision                                 , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    peakHeightElementCount=3
    return
  end function peakHeightElementCount

  function peakHeightExtract(self,node,time,instance)
    !!{
    Implement a peakHeight output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    double precision                                 , dimension(:) , allocatable :: peakHeightExtract
    class           (nodePropertyExtractorPeakHeight), intent(inout), target      :: self
    type            (treeNode                       ), intent(inout), target      :: node
    double precision                                 , intent(in   )              :: time
    type            (multiCounter                   ), intent(inout), optional    :: instance
    class           (nodeComponentBasic             ), pointer                    :: basic
    double precision                                                              :: criticalOverdensityLastIsolated, densityFieldRootVariance, &
         &                                                                           peakHeightNu
    !$GLC attributes unused :: time, instance

    basic => node%basic()
    criticalOverdensityLastIsolated=self%criticalOverdensity_     %value       (mass=basic%mass(),time=basic%timeLastIsolated())
    densityFieldRootVariance       =self%cosmologicalMassVariance_%rootVariance(mass=basic%mass(),time=basic%timeLastIsolated())
    peakHeightNu                   =criticalOverdensityLastIsolated/densityFieldRootVariance
    allocate(peakHeightExtract(3))
    peakHeightExtract=[                                 &
         &             criticalOverdensityLastIsolated, &
         &             densityFieldRootVariance       , &
         &             peakHeightNu                     &
         &            ]
    return
  end function peakHeightExtract

  function peakHeightNames(self,time)
    !!{
    Return the names of the {\normalfont \ttfamily peakHeight} properties.
    !!}
    implicit none
    type            (varying_string                 ), dimension(:) , allocatable :: peakHeightNames
    class           (nodePropertyExtractorPeakHeight), intent(inout)              :: self
    double precision                                 , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(peakHeightNames(3))
    peakHeightNames=[                                              &
         &           var_str('criticalOverdensityLastIsolated'  ), &
         &           var_str('haloSigma'                        ), &
         &           var_str('haloPeakHeightNu'                 )  &
         &          ]
    return
  end function peakHeightNames

  function peakHeightDescriptions(self,time)
    !!{
    Return the descriptions of the {\normalfont \ttfamily peakHeight} properties.
    !!}
    implicit none
    type            (varying_string                 ), dimension(:) , allocatable :: peakHeightDescriptions
    class           (nodePropertyExtractorPeakHeight), intent(inout)              :: self
    double precision                                 , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(peakHeightDescriptions(3))
    peakHeightDescriptions=[                                                                              &
         &                  var_str('Critical overdensity at the time when the halo was last isolated.'), &
         &                  var_str('The mass fluctuation on the scale of the halo.'                   ), &
         &                  var_str('The peak height, ν, of the halo.'                                 )  &
         &                 ]
    return
  end function peakHeightDescriptions

  function peakHeightUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily peakHeight} properties in the SI system.
    !!}
    implicit none
    double precision                                 , dimension(:) , allocatable :: peakHeightUnitsInSI
    class           (nodePropertyExtractorPeakHeight), intent(inout)              :: self
    double precision                                 , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(peakHeightUnitsInSI(3))
    peakHeightUnitsInSI=[0.0d0,0.0d0,0.0d0]
    return
  end function peakHeightUnitsInSI

  integer function peakHeightType(self)
    !!{
    Return the type of the peakHeight property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorPeakHeight), intent(inout) :: self
    !$GLC attributes unused :: self

    peakHeightType=outputAnalysisPropertyTypeLinear
    return
  end function peakHeightType
