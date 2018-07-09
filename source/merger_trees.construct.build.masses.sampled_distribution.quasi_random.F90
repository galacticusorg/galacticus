!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implementation of a merger tree masses class which samples masses from a distribution using quasi-random sampling.

  !# <mergerTreeBuildMasses name="mergerTreeBuildMassesSampledDistributionQuasiRandom" defaultThreadPrivate="yes">
  !#  <description>A merger tree masses class which samples masses from a distribution using quasi-random sampling.</description>
  !# </mergerTreeBuildMasses>
  type, extends(mergerTreeBuildMassesSampledDistribution) :: mergerTreeBuildMassesSampledDistributionQuasiRandom
     !% Implementation of a merger tree masses class which samples masses from a distribution with pseudi-random sampling.
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
    !% Constructor for the {\normalfont \ttfamily sampledDistributionQuasiRandom} merger tree masses class which takes a parameter set
    !% as input.
    use Input_Parameters
    implicit none
    type(mergerTreeBuildMassesSampledDistributionQuasiRandom)                :: self
    type(inputParameters                                    ), intent(inout) :: parameters

    self%mergerTreeBuildMassesSampledDistribution=mergerTreeBuildMassesSampledDistribution(parameters)
    return
  end function sampledDistributionQuasiRandomConstructorParameters

  function sampledDistributionQuasiRandomConstructorInternal(massTreeMinimum,massTreeMaximum,treesPerDecade,mergerTreeBuildMassDistribution_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sampledDistributionQuasiRandom} merger tree masses class.
    implicit none
    type            (mergerTreeBuildMassesSampledDistributionQuasiRandom)                        :: self
    class           (mergerTreeBuildMassDistributionClass               ), intent(in   ), target :: mergerTreeBuildMassDistribution_
    double precision                                                     , intent(in   )         :: massTreeMinimum                 , massTreeMaximum, &
         &                                                                                          treesPerDecade
    !# <constructorAssign variables="massTreeMinimum, massTreeMaximum, treesPerDecade, *mergerTreeBuildMassDistribution_"/>

    return
  end function sampledDistributionQuasiRandomConstructorInternal

  subroutine sampledDistributionQuasiRandomDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sampledDistributionQuasiRandom} merger tree masses class.
    implicit none
    type(mergerTreeBuildMassesSampledDistributionQuasiRandom), intent(inout) :: self
    
    !# <objectDestructor name="self%mergerTreeBuildMassDistribution_"/>
    return
  end subroutine sampledDistributionQuasiRandomDestructor
  
  subroutine sampledDistributionQuasiRandomSampleCMF(self,x)
    !% Generate a quasiRandom sample of points from the merger tree mass distribution.
    use, intrinsic :: ISO_C_Binding
    use               FGSL
    use               Quasi_Random
    implicit none
    class           (mergerTreeBuildMassesSampledDistributionQuasiRandom), intent(inout)               :: self
    double precision                                                     , intent(  out), dimension(:) :: x
    type            (fgsl_qrng                                          )                              :: quasiSequenceObject
    logical                                                                                            :: quasiSequenceReset =.true.
    integer         (c_size_t                                           )                              :: iTree
    !GCC$ attributes unused :: self

    do iTree=1,size(x)
       x(iTree)=Quasi_Random_Get(quasiSequenceObject,reset=quasiSequenceReset)
    end do
    call Quasi_Random_Free(quasiSequenceObject)
    return
  end subroutine sampledDistributionQuasiRandomSampleCMF
