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
  Implementation of a merger tree masses class which samples masses from a distribution uniformly.
  !!}

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesSampledDistributionUniform">
   <description>A merger tree masses class which samples masses from a distribution uniformly.</description>
  </mergerTreeBuildMasses>
  !!]
  type, extends(mergerTreeBuildMassesSampledDistribution) :: mergerTreeBuildMassesSampledDistributionUniform
     !!{
     Implementation of a merger tree masses class which samples masses from a distribution with uniform sampling.
     !!}
     private
   contains
     final     ::              sampledDistributionUniformDestructor
     procedure :: sampleCMF => sampledDistributionUniformSampleCMF
  end type mergerTreeBuildMassesSampledDistributionUniform

  interface mergerTreeBuildMassesSampledDistributionUniform
     module procedure sampledDistributionUniformConstructorParameters
     module procedure sampledDistributionUniformConstructorInternal
  end interface mergerTreeBuildMassesSampledDistributionUniform

contains

  function sampledDistributionUniformConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassesSampledDistributionUniform} merger tree masses class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeBuildMassesSampledDistributionUniform)                :: self
    type(inputParameters                                ), intent(inout) :: parameters

    self%mergerTreeBuildMassesSampledDistribution=mergerTreeBuildMassesSampledDistribution(parameters)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sampledDistributionUniformConstructorParameters

  function sampledDistributionUniformConstructorInternal(massTreeMinimum,massTreeMaximum,treesPerDecade,mergerTreeBuildMassDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassesSampledDistributionUniform} merger tree masses class.
    !!}
    implicit none
    type            (mergerTreeBuildMassesSampledDistributionUniform)                        :: self
    class           (mergerTreeBuildMassDistributionClass           ), intent(in   ), target :: mergerTreeBuildMassDistribution_
    double precision                                                 , intent(in   )         :: massTreeMinimum                 , massTreeMaximum, &
         &                                                                                      treesPerDecade
    !![
    <constructorAssign variables="massTreeMinimum, massTreeMaximum, treesPerDecade, *mergerTreeBuildMassDistribution_"/>
    !!]

    return
  end function sampledDistributionUniformConstructorInternal

  subroutine sampledDistributionUniformDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildMassesSampledDistributionUniform} merger tree masses class.
    !!}
    implicit none
    type(mergerTreeBuildMassesSampledDistributionUniform), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBuildMassDistribution_"/>
    !!]
    return
  end subroutine sampledDistributionUniformDestructor

  subroutine sampledDistributionUniformSampleCMF(self,x)
    !!{
    Generate a uniform sample of points from the merger tree mass distribution.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLinear
    implicit none
    class           (mergerTreeBuildMassesSampledDistributionUniform), intent(inout)               :: self
    double precision                                                 , intent(  out), dimension(:) :: x
    !$GLC attributes unused :: self

    x=Make_Range(0.0d0,1.0d0,size(x),rangeType=rangeTypeLinear,rangeBinned=.true.)
    return
  end subroutine sampledDistributionUniformSampleCMF
