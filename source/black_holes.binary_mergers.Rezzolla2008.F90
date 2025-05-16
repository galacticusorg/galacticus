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
  Implements a black hole binary merger in which the black hole mass and spin resulting from binary mergers utilizing the
  approximations of \cite{rezzolla_final_2008}.
  !!}

  !![
  <blackHoleBinaryMerger name="blackHoleBinaryMergerRezzolla2008">
   <description>
    A black hole binary merger class that uses the fitting function of \cite{rezzolla_final_2008} to compute the spin of the
    black hole resulting from a binary merger. The mass of the resulting black hole is assumed to equal the sum of the mass of
    the initial black holes (i.e. there is negligible energy loss through gravitational waves).
   </description>
  </blackHoleBinaryMerger>
  !!]
  type, extends(blackHoleBinaryMergerClass) :: blackHoleBinaryMergerRezzolla2008
     !!{
     A black hole binary merger class in which the black hole mass and spin resulting from binary mergers utilizing the approximations of \cite{rezzolla_final_2008}.
     !!}
     private
   contains
     procedure :: merge => rezzolla2008Merge
  end type blackHoleBinaryMergerRezzolla2008

  interface blackHoleBinaryMergerRezzolla2008
     !!{
     Constructors for the \refClass{blackHoleBinaryMergerRezzolla2008} black hole binary merger class.
     !!}
     module procedure rezzolla2008ConstructorParameters
  end interface blackHoleBinaryMergerRezzolla2008

contains

  function rezzolla2008ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleBinaryMergerRezzolla2008} black hole binary merger class which takes a parameter list as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(blackHoleBinaryMergerRezzolla2008)                :: self
    type(inputParameters                  ), intent(inout) :: parameters

    self=blackHoleBinaryMergerRezzolla2008()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function rezzolla2008ConstructorParameters

  subroutine rezzolla2008Merge(self,massBlackHoleA,massBlackHoleB,spinBlackHoleA,spinBlackHoleB,massBlackHoleFinal,spinBlackHoleFinal)
    !!{
    Computes the mass and spin of a black hole resulting from a binary merger utilizing the approximations of
    \cite{rezzolla_final_2008}.
    !!}
    implicit none
    class           (blackHoleBinaryMergerRezzolla2008), intent(inout) :: self
    double precision                                   , intent(in   ) :: massBlackHoleA             , massBlackHoleB             , spinBlackHoleA                 , &
         &                                                                spinBlackHoleB
    double precision                                   , intent(  out) :: massBlackHoleFinal         , spinBlackHoleFinal
    ! Parameters of the fitting functions used by Rezzolla et al.
    double precision                                   , parameter     :: s4                =-0.129d0, s5                =-0.384d0, t0                    =-2.686d0, &
         &                                                                t2                =-3.454d0, t3                =2.353d0
    ! Assumed fixed values for the various alignment angles.
    double precision                                   , parameter     :: cosinePhi         =1.0d0   , cosineTheta       =1.0d0   , cosineXi              =1.0d0
    double precision                                                   :: massBlackHole1             , massBlackHole2             , spinBlackHole1                 , &
         &                                                                spinBlackHole2             , massRatio                  , orbitalAngularMomentum         , &
         &                                                                symmetricMassRatio         , argumentSqrt
    !$GLC attributes unused :: self

    ! Check for case of two zero-mass black holes.
    if (massBlackHoleA <= 0.0d0 .and. massBlackHoleB <= 0.0d0) then
       massBlackHoleFinal=0.0d0
       spinBlackHoleFinal=0.0d0
       return
    end if
    ! Find which is the more massive of the two black holes.
    if (massBlackHoleA < massBlackHoleB) then
       massBlackHole1=massBlackHoleB
       massBlackHole2=massBlackHoleA
       spinBlackHole1=spinBlackHoleB
       spinBlackHole2=spinBlackHoleA
    else
       massBlackHole1=massBlackHoleA
       massBlackHole2=massBlackHoleB
       spinBlackHole1=spinBlackHoleA
       spinBlackHole2=spinBlackHoleB
    end if
    ! Assume negligible mass loss (Rezzolla et al. report 5-7% mass loss via gravitational radiation).
    massBlackHoleFinal=massBlackHole1+massBlackHole2
    ! Compute the black hole mass ratio.
    massRatio=massBlackHole2/massBlackHole1
    symmetricMassRatio=massRatio/(1.0d0+massRatio)**2
    ! Evaluate the fitting formula for the orbital angular momentum (which accounts for losses to gravitational waves).
    orbitalAngularMomentum=s4*(spinBlackHole1**2+spinBlackHole2**2*massRatio**4+2.0d0*spinBlackHole1*spinBlackHole2*massRatio**2&
         &*cosinePhi)/(1.0d0+massRatio**2)**2+((s5*symmetricMassRatio+t0+2.0d0)/(1.0d0+massRatio**2))*(spinBlackHole1*cosineTheta&
         &+spinBlackHole2*massRatio**2*cosinePhi)+2.0d0*sqrt(3.0d0)+t2*symmetricMassRatio+t3*symmetricMassRatio**2
    ! Compute the spin of the final black hole by evaluating the fitting formula.
    argumentSqrt      =+max(                                                                  &
         &                  +0.0d0                                                          , &
         &                  +      spinBlackHole1               **2                           &
         &                  +                     spinBlackHole2**2*massRatio**4              &
         &                  +2.0d0*spinBlackHole1*spinBlackHole2   *massRatio**2*cosinePhi    &
         &                  +2.0d0                                                            &
         &                  *(                                                                &
         &                    +    spinBlackHole1                               *cosineTheta  &
         &                    +                   spinBlackHole2   *massRatio**2*cosineXi     &
         &                   )                                                                &
         &                  * orbitalAngularMomentum*massRatio                                &
         &                  +(orbitalAngularMomentum*massRatio)**2                            &
         &                  )
    spinBlackHoleFinal=+sqrt(argumentSqrt)                                                    &
         &             /(                                                                     &
         &               +1.0d0                                                               &
         &               +massRatio                                                           &
         &              )**2
    return
  end subroutine rezzolla2008Merge
