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

  use :: Cosmological_Density_Field, only : cosmologicalMassVariance, cosmologicalMassVarianceClass, criticalOverdensity, criticalOverdensityClass
  use :: Linear_Growth             , only : linearGrowthClass

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
     class(linearGrowthClass)            , pointer :: linearGrowth_             => null()
   contains
     final     ::                 peakHeightDestructor
     procedure :: elementCount => peakHeightElementCount
     procedure :: extract      => peakHeightExtract
     procedure :: names        => peakHeightNames
     procedure :: descriptions => peakHeightDescriptions
     procedure :: unitsInSI    => peakHeightUnitsInSI
  end type nodePropertyExtractorPeakHeight

  interface nodePropertyExtractorPeakHeight
     !!{
     Constructors for the {\normalfont \ttfamily peakHeight} output extractor class.
     !!}
     module procedure peakHeightConstructorParameters
     module procedure peakHeightConstructorInternal
  end interface nodePropertyExtractorPeakHeight

contains

  function peakHeightConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily peakHeight} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorPeakHeight)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(criticalOverdensityClass       ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_
    class(linearGrowthClass              ), pointer       :: linearGrowth_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !!]
    self=nodePropertyExtractorPeakHeight(criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    !!]
    return
  end function peakHeightConstructorParameters

  function peakHeightConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily peakHeight} output extractor property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorPeakHeight)                        :: self
    class(criticalOverdensityClass       ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_
    class(linearGrowthClass              ), intent(in   ), target :: linearGrowth_ 
    !![
    <constructorAssign variables="*criticalOverdensity_, *cosmologicalMassVariance_, *linearGrowth_"/>
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
    <objectDestructor name="self%linearGrowth_"            />
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

    peakHeightElementCount=4
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
         &                                                                           peakHeightNu                   , linearGrowthFactor
    !$GLC attributes unused :: time, instance

    basic                           =>  node                          %basic       (                                                         )
    criticalOverdensityLastIsolated =   self%criticalOverdensity_     %value       (mass=basic%mass(),time=basic%timeLastIsolated(),node=node)
    densityFieldRootVariance        =   self%cosmologicalMassVariance_%rootVariance(mass=basic%mass(),time=basic%timeLastIsolated()          )
    peakHeightNu                    =  +criticalOverdensityLastIsolated                                                                        &
         &                             /densityFieldRootVariance
    linearGrowthFactor              =   self%LinearGrowth_            %value        (                 time=basic%timeLastIsolated()          )
    allocate(peakHeightExtract(4))
    peakHeightExtract=[                                 &
         &             criticalOverdensityLastIsolated, &
         &             densityFieldRootVariance       , &
         &             peakHeightNu                   , &
         &             linearGrowthFactor               &
         &            ]
    return
  end function peakHeightExtract

  subroutine peakHeightNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily peakHeight} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorPeakHeight), intent(inout)                             :: self
    double precision                                 , intent(in   )                             :: time
    type            (varying_string                 ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(4))
    names(1)=var_str('criticalOverdensityLastIsolated'  )
    names(2)=var_str('haloSigma'                        )
    names(3)=var_str('haloPeakHeightNu'                 )
    names(4)=var_str('linearGrowth'                     )
    return
  end subroutine peakHeightNames

  subroutine peakHeightDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily peakHeight} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorPeakHeight), intent(inout)                             :: self
    double precision                                 , intent(in   )                             :: time
    type            (varying_string                 ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(4))
    descriptions(1)=var_str('Critical overdensity at the time when the halo was last isolated.'          )
    descriptions(2)=var_str('The mass fluctuation on the scale of the halo.'                             )
    descriptions(3)=var_str('The peak height, ν, of the halo.'                                           )
    descriptions(4)=var_str('The cosmological growth factor at the time when the halo was last isolated.')
    return
  end subroutine peakHeightDescriptions

  function peakHeightUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily peakHeight} properties in the SI system.
    !!}
    implicit none
    double precision                                 , dimension(:) , allocatable :: peakHeightUnitsInSI
    class           (nodePropertyExtractorPeakHeight), intent(inout)              :: self
    double precision                                 , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(peakHeightUnitsInSI(4))
    peakHeightUnitsInSI=[0.0d0,0.0d0,0.0d0,0.0d0]
    return
  end function peakHeightUnitsInSI
