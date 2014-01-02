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

!+    Contributions to this file made by:  Martin White.

  !% Implements calculations of satellite merging times using the \cite{wetzel_what_2010} method.

  !# <satelliteMergingTimescales name="satelliteMergingTimescalesWetzelWhite2010" />

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesWetzelWhite2010
     !% A class implementing the \cite{wetzel_what_2010} method for satellite merging timescales.
     private
   contains
     final     ::                     wetzelWhite2010Destructor
     procedure :: timeUntilMerging => wetzelWhite2010TimeUntilMerging
  end type satelliteMergingTimescalesWetzelWhite2010

  interface satelliteMergingTimescalesWetzelWhite2010
     !% Constructors for the \cite{wetzel_what_2010} merging timescale class.
     module procedure wetzelWhite2010DefaultConstructor
  end interface satelliteMergingTimescalesWetzelWhite2010

contains

  function wetzelWhite2010DefaultConstructor()
    !% Default constructor for the \cite{wetzel_what_2010} merging timescale class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type(satelliteMergingTimescalesWetzelWhite2010) :: wetzelWhite2010DefaultConstructor

    return
  end function wetzelWhite2010DefaultConstructor

  elemental subroutine wetzelWhite2010Destructor(self)
    !% Default constructor for the \cite{wetzel_what_2010} merging timescale class.
    implicit none
    type(satelliteMergingTimescalesWetzelWhite2010), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine wetzelWhite2010Destructor

  double precision function wetzelWhite2010TimeUntilMerging(self,thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{wetzel_what_2010} method.
    use Galacticus_Nodes
    use Dynamical_Friction_Timescale_Utilities
    use Cosmology_Functions
    use Kepler_Orbits
    implicit none
    class           (satelliteMergingTimescalesWetzelWhite2010), intent(inout)          :: self
    type            (treeNode                                 ), intent(inout), pointer :: thisNode
    type            (keplerOrbit                              ), intent(inout)          :: thisOrbit
    type            (treeNode                                 )               , pointer :: hostNode
    class           (nodeComponentBasic                       )               , pointer :: hostBasic                , thisBasic
    class           (cosmologyFunctionsClass                  )               , pointer :: cosmologyFunctionsDefault
    double precision                                           , parameter              :: timeScaleNormalization   =0.2d0      !   C_dyn from Wetzel & White (2010).
    double precision                                                                    :: massRatio

    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Find the host node.
    hostNode => thisNode%parent
    ! Compute mass ratio.
    thisBasic => thisNode%basic()
    hostBasic => hostNode%basic()
    massRatio=hostBasic%mass()/thisBasic%mass()
    ! Compute dynamical friction timescale using eqn. (2) from Wetzel & White (2010).
    wetzelWhite2010TimeUntilMerging= Dynamical_Friction_Timescale_Multiplier()   &
         &                          *timeScaleNormalization                      &
         &                          /cosmologyFunctionsDefault%expansionRate(    &
         &                            cosmologyFunctionsDefault%expansionFactor( &
         &                             thisBasic%time()                          &
         &                                                                     ) &
         &                                                                  )    &
         &                          *           massRatio                        &
         &                          /log(1.0d0+massRatio)
    return
  end function wetzelWhite2010TimeUntilMerging
