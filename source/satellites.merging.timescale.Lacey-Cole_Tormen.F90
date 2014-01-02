!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Implements calculations of satellite merging times using the \cite{lacey_merger_1993} method with a parameterization of
  !% orbital parameters designed to fit the results of \cite{tormen_rise_1997} as described by \cite{cole_hierarchical_2000}.
 
  !# <satelliteMergingTimescales name="satelliteMergingTimescalesLaceyCole1993Tormen" />
  use FGSL

  type, extends(satelliteMergingTimescalesLaceyCole1993) :: satelliteMergingTimescalesLaceyCole1993Tormen
     !% A class implementing the \cite{cole_hierarchical_2000} method for satellite merging timescales.
     private
     type   (fgsl_rng) :: clonedPseudoSequenceObject, randomSequenceObject
     logical           :: resetRandomSequence       , resetRandomSequenceSnapshot
   contains
     final     ::                     laceyCole1993TormenDestructor
     procedure :: stateStore       => laceyCole1993TormenStateStore
     procedure :: stateRestore     => laceyCole1993TormenStateRestore
     procedure :: stateSnapshot    => laceyCole1993TormenStateSnapshot
     procedure :: timeUntilMerging => laceyCole1993TormenTimeUntilMerging
  end type satelliteMergingTimescalesLaceyCole1993Tormen

  interface satelliteMergingTimescalesLaceyCole1993Tormen
     !% Constructors for the \cite{cole_hierarchical_2000} merging timescale class.
     module procedure laceyCole1993TormenDefaultConstructor
  end interface satelliteMergingTimescalesLaceyCole1993Tormen

contains

  function laceyCole1993TormenDefaultConstructor()
    !% Default constructor for the \cite{cole_hierarchical_2000} merging timescale class.
    implicit none
    type(satelliteMergingTimescalesLaceyCole1993Tormen) :: laceyCole1993TormenDefaultConstructor

    laceyCole1993TormenDefaultConstructor%resetRandomSequence=.true.
    return
  end function laceyCole1993TormenDefaultConstructor

  subroutine laceyCole1993TormenDestructor(self)
    !% Destructor for the \cite{cole_hierarchical_2000} merging timescale class.
    use Gaussian_Random
    implicit none
    type(satelliteMergingTimescalesLaceyCole1993Tormen), intent(inout) :: self

    ! Destroy the random number object.
    if (self%resetRandomSequence        ) call Gaussian_Random_Free(self%randomSequenceObject      )
    if (self%resetRandomSequenceSnapshot) call Gaussian_Random_Free(self%clonedPseudoSequenceObject)
    return
  end subroutine laceyCole1993TormenDestructor

  double precision function laceyCole1993TormenTimeUntilMerging(self,thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{lacey_merger_1993} method with a parameterization of orbital
    !% parameters designed to fit the results of \cite{tormen_rise_1997} as described by \cite{cole_hierarchical_2000}.
    use Galacticus_Nodes
    use Kepler_Orbits
    use Gaussian_Random
    implicit none
    class           (satelliteMergingTimescalesLaceyCole1993Tormen)           , intent(inout)          :: self
    type            (treeNode                                     )           , intent(inout), pointer :: thisNode
    type            (keplerOrbit                                  )           , intent(inout)          :: thisOrbit
    type            (treeNode                                     )                          , pointer :: hostNode
    double precision                                               , parameter                         :: orbitalFactorDistributionSigma=0.26d0                          !   Cole et al. (2000).
    double precision                                               , parameter                         :: orbitalFactorDistributionMean =-0.14d0                         !   Cole et al. (2000).
    double precision                                                                                   :: log10OrbitalFactor                          , randomDeviate, &
         &                                                                                                orbitalFactor

    ! Find the host node.
    hostNode => thisNode%parent
    ! Compute the orbital factor - selected at random from a lognormal distribution.
    randomDeviate=Gaussian_Random_Get(self%randomSequenceObject,orbitalFactorDistributionSigma,self%resetRandomSequence)
    log10OrbitalFactor=orbitalFactorDistributionMean+randomDeviate
    orbitalFactor=10.0d0**log10OrbitalFactor
    ! Compute the timescale.
    laceyCole1993TormenTimeUntilMerging=orbitalFactor*self%timeUntilMergingMassDependence(thisNode)
    return
  end function laceyCole1993TormenTimeUntilMerging

  subroutine laceyCole1993TormenStateSnapshot(self)
    !% Store a snapshot of the random number generator internal state.
    implicit none
    class(satelliteMergingTimescalesLaceyCole1993Tormen), intent(inout) :: self

    if (.not.self%resetRandomSequence) self%clonedPseudoSequenceObject=FGSL_Rng_Clone(self%randomSequenceObject)
    self%resetRandomSequenceSnapshot=self%resetRandomSequence
    return
  end subroutine laceyCole1993TormenStateSnapshot

  subroutine laceyCole1993TormenStateStore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    class  (satelliteMergingTimescalesLaceyCole1993Tormen), intent(inout) :: self
    integer                                               , intent(in   ) :: stateFile
    type   (fgsl_file                                    ), intent(in   ) :: fgslStateFile
    
    write (stateFile) self%resetRandomSequenceSnapshot
    if (.not.self%resetRandomSequenceSnapshot) call Pseudo_Random_Store(self%clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine laceyCole1993TormenStateStore
  
  subroutine laceyCole1993TormenStateRestore(self,stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use Pseudo_Random
    implicit none
    class  (satelliteMergingTimescalesLaceyCole1993Tormen), intent(inout) :: self
    integer                                               , intent(in   ) :: stateFile
    type   (fgsl_file                                    ), intent(in   ) :: fgslStateFile
    
    read (stateFile) self%resetRandomSequence
    if (.not.self%resetRandomSequence) call Pseudo_Random_Retrieve(self%randomSequenceObject,fgslStateFile)
    return
  end subroutine laceyCole1993TormenStateRestore
  
