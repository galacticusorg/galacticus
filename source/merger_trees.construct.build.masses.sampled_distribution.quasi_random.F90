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
  Implementation of a merger tree masses class which samples masses from a distribution using quasi-random sampling.
  !!}

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesSampledDistributionQuasiRandom">
   <description>A merger tree masses class which samples masses from a distribution using quasi-random sampling.</description>
  </mergerTreeBuildMasses>
  !!]
  type, extends(mergerTreeBuildMassesSampledDistribution) :: mergerTreeBuildMassesSampledDistributionQuasiRandom
     !!{
     Implementation of a merger tree masses class which samples masses from a distribution with pseudo-random sampling.
     !!}
     private
   contains
     final     ::              sampledDistributionQuasiRandomDestructor
     procedure :: sampleCMF => sampledDistributionQuasiRandomSampleCMF
  end type mergerTreeBuildMassesSampledDistributionQuasiRandom

  interface mergerTreeBuildMassesSampledDistributionQuasiRandom
     module procedure sampledDistributionQuasiRandomConstructorParameters
     module procedure sampledDistributionQuasiRandomConstructorInternal
  end interface mergerTreeBuildMassesSampledDistributionQuasiRandom

contains

  function sampledDistributionQuasiRandomConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassesSampledDistributionQuasiRandom} merger tree masses class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeBuildMassesSampledDistributionQuasiRandom)                :: self
    type(inputParameters                                    ), intent(inout) :: parameters

    self%mergerTreeBuildMassesSampledDistribution=mergerTreeBuildMassesSampledDistribution(parameters)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sampledDistributionQuasiRandomConstructorParameters

  function sampledDistributionQuasiRandomConstructorInternal(massTreeMinimum,massTreeMaximum,treesPerDecade,mergerTreeBuildMassDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassesSampledDistributionQuasiRandom} merger tree masses class.
    !!}
    implicit none
    type            (mergerTreeBuildMassesSampledDistributionQuasiRandom)                        :: self
    class           (mergerTreeBuildMassDistributionClass               ), intent(in   ), target :: mergerTreeBuildMassDistribution_
    double precision                                                     , intent(in   )         :: massTreeMinimum                 , massTreeMaximum, &
         &                                                                                          treesPerDecade
    !![
    <constructorAssign variables="massTreeMinimum, massTreeMaximum, treesPerDecade, *mergerTreeBuildMassDistribution_"/>
    !!]

    return
  end function sampledDistributionQuasiRandomConstructorInternal

  subroutine sampledDistributionQuasiRandomDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildMassesSampledDistributionQuasiRandom} merger tree masses class.
    !!}
    implicit none
    type(mergerTreeBuildMassesSampledDistributionQuasiRandom), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBuildMassDistribution_"/>
    !!]
    return
  end subroutine sampledDistributionQuasiRandomDestructor

  subroutine sampledDistributionQuasiRandomSampleCMF(self,x)
    !!{
    Generate a quasiRandom sample of points from the merger tree mass distribution.
    !!}
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Quasi_Random_Sequences, only : quasiRandomNumberGenerator, gsl_qrng_sobol
    implicit none
    class           (mergerTreeBuildMassesSampledDistributionQuasiRandom), intent(inout)               :: self
    double precision                                                     , intent(  out), dimension(:) :: x
    type            (quasiRandomNumberGenerator                         )                              :: quasiRandomSequence
    integer         (c_size_t                                           )                              :: iTree
    !$GLC attributes unused :: self

    quasiRandomSequence=quasiRandomNumberGenerator(gsl_qrng_sobol)
    do iTree=1,size(x)
       x(iTree)=quasiRandomSequence%get()
    end do
    return
  end subroutine sampledDistributionQuasiRandomSampleCMF
