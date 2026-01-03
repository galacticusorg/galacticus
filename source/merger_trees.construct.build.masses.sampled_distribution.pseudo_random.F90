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
  Implementation of a merger tree masses class which samples masses from a distribution using pseudo-random sampling.
  !!}

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesSampledDistributionPseudoRandom">
   <description>A merger tree masses class which samples masses from a distribution using pseudo-random sampling.</description>
  </mergerTreeBuildMasses>
  !!]
  type, extends(mergerTreeBuildMassesSampledDistribution) :: mergerTreeBuildMassesSampledDistributionPseudoRandom
     !!{
     Implementation of a merger tree masses class which samples masses from a distribution with pseudo-random sampling.
     !!}
     private
     class(randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::              sampledDistributionPseudoRandomDestructor
     procedure :: sampleCMF => sampledDistributionPseudoRandomSampleCMF
  end type mergerTreeBuildMassesSampledDistributionPseudoRandom

  interface mergerTreeBuildMassesSampledDistributionPseudoRandom
     module procedure sampledDistributionPseudoRandomConstructorParameters
     module procedure sampledDistributionPseudoRandomConstructorInternal
  end interface mergerTreeBuildMassesSampledDistributionPseudoRandom

contains

  function sampledDistributionPseudoRandomConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassesSampledDistributionPseudoRandom} merger tree masses class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeBuildMassesSampledDistributionPseudoRandom)                :: self
    type(inputParameters                                     ), intent(inout) :: parameters

    self%mergerTreeBuildMassesSampledDistribution=mergerTreeBuildMassesSampledDistribution(parameters)
    !![
    <objectBuilder class="randomNumberGenerator" name="self%randomNumberGenerator_" source="parameters"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sampledDistributionPseudoRandomConstructorParameters

  function sampledDistributionPseudoRandomConstructorInternal(massTreeMinimum,massTreeMaximum,treesPerDecade,mergerTreeBuildMassDistribution_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassesSampledDistributionPseudoRandom} merger tree masses class.
    !!}
    implicit none
    type            (mergerTreeBuildMassesSampledDistributionPseudoRandom)                        :: self
    class           (mergerTreeBuildMassDistributionClass                ), intent(in   ), target :: mergerTreeBuildMassDistribution_
    class           (randomNumberGeneratorClass                          ), intent(in   ), target :: randomNumberGenerator_
    double precision                                                      , intent(in   )         :: massTreeMinimum                 , massTreeMaximum, &
         &                                                                                           treesPerDecade
    !![
    <constructorAssign variables="massTreeMinimum, massTreeMaximum, treesPerDecade, *mergerTreeBuildMassDistribution_, *randomNumberGenerator_"/>
    !!]

    return
  end function sampledDistributionPseudoRandomConstructorInternal

  subroutine sampledDistributionPseudoRandomDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildMassesSampledDistributionPseudoRandom} merger tree masses class.
    !!}
    implicit none
    type(mergerTreeBuildMassesSampledDistributionPseudoRandom), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBuildMassDistribution_"/>
    <objectDestructor name="self%randomNumberGenerator_"          />
    !!]
    return
  end subroutine sampledDistributionPseudoRandomDestructor

  subroutine sampledDistributionPseudoRandomSampleCMF(self,x)
    !!{
    Generate a pseudo-random sample of points from the merger tree mass distribution.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (mergerTreeBuildMassesSampledDistributionPseudoRandom), intent(inout)               :: self
    double precision                                                      , intent(  out), dimension(:) :: x
    integer         (c_size_t                                            )                              :: iTree

    do iTree=1,size(x)
       x(iTree)=self%randomNumberGenerator_%uniformSample()
    end do
    return
  end subroutine sampledDistributionPseudoRandomSampleCMF
