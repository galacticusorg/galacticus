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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

  !!{
  Implements a black hole binary recoil velocity class which follows the formulae in \cite{campanelli_large_2007}, derived from
  post-Newtonian evaluations.
  !!}

  !![
  <blackHoleBinaryRecoil name="blackHoleBinaryRecoilCampanelli2007">
   <description>
    A black hole binary recoil class that computes the recoil velocity during a black hole binary merger due to the emission of
    gravitational waves, following the formulae derived in \cite {campanelli_large_2007}
    \begin{equation}
    V_\mathrm{recoil}=V_\mathrm{m}\mathbf{\hat{e_1}}+V_\perp(\cos\xi\mathbf{\hat{e_1}}+\sin\xi\mathbf{\hat{e_2}})+V_\parallel\mathbf{\hat{e_2}}
    \end{equation}
    with:
    \begin{equation}
    V_\mathrm{m}=A\frac{q^2(1-q)}{(1+q)^5}(1+B\frac{q}{(1+q)^2})
    \end{equation}
    \begin{equation}
    V_\perp=H\frac{q^2}{(1+q)^5}(\alpha^\parallel_2-q\alpha^\parallel_1)
    \end{equation}
    \begin{equation}
    V_\parallel=K\cos(\theta-\theta_0)\frac{q^2}{(1+q)^5}(\alpha^\perp_2-q\alpha^\perp_1)
    \end{equation}
    where $\theta$ is defined as the angle between the inplane \gls{component} of $\Delta$ and the infall direction at
    merger. $q$ is the mass ratio of the black holes as $q=M_{\bullet,1}/M_{\bullet,2}$ and
    $\alpha_i=\mathbf{S_i}/M_{\bullet,i}$ depends of the spin and mass of the black hole $\xi$ measures the angle between the
    unequal mass and the spin contribution to the recoil velocity in the orbital plane. $\mathbf{\hat{e_1}} ,
    \mathbf{\hat{e_2}}$ are orthogonal unit vectors in the orbital plane. Our method assumes the spin of the second black hole
    is randomly generated, while that of the first is aligned with the angular momentum of the system. The constants used are
    retrieved from the articles by: \cite{koppitz_recoil_2007} for $H=(7.3\pm 0.3)10^3$~km/s, \cite{gonzalez_maximum_2007} for
    $A=1.2 \times 10^4$~km/s $B=-0.93$, \cite{gonzalez_supermassive_2007} for $K\cos(\delta\theta)=(6,-5.3)10^4$~km/s and
    $K=(6.0\pm 0.1)10^4$~km/s.
   </description>
  </blackHoleBinaryRecoil>
  !!]
  type, extends(blackHoleBinaryRecoilClass) :: blackHoleBinaryRecoilCampanelli2007
     !!{
     A black hole binary recoil class which follows the formulae in \cite{campanelli_large_2007}, derived from post-Newtonian
     evaluations.
     !!}
     private
   contains
     procedure :: velocity => campanelli2007Velocity
  end type blackHoleBinaryRecoilCampanelli2007

  interface blackHoleBinaryRecoilCampanelli2007
     !!{
     Constructors for the \refClass{blackHoleBinaryRecoilCampanelli2007} black hole binary recoil class.
     !!}
     module procedure campanelli2007ConstructorParameters
  end interface blackHoleBinaryRecoilCampanelli2007

contains

  function campanelli2007ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleBinaryRecoilCampanelli2007} black hole binary recoiled class which takes a parameter list as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(blackHoleBinaryRecoilCampanelli2007)                :: self
    type(inputParameters                    ), intent(inout) :: parameters

    self=blackHoleBinaryRecoilCampanelli2007()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function campanelli2007ConstructorParameters

  double precision function campanelli2007Velocity(self,blackHole1,blackHole2)
    !!{
    Returns the recoil velocity of a black hole binary, accounting for the parallel and perpendicular velocities, plus that of
    the binary's center of mass. Constants used are retrieved from the articles by: \cite{koppitz_recoil_2007} for $H=(7.3\pm
    0.3)10^3$~km/s, \cite{gonzalez_maximum_2007} for $A=1.2 \times 10^4$~km/s $B=-0.93$, \cite{gonzalez_supermassive_2007} for $K
    \cos(\delta\theta)=(6,-5.3)10^4$~km/s and $K=(6.0\pm 0.1)10^4$~km/s.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (blackHoleBinaryRecoilCampanelli2007), intent(inout) :: self
    class           (nodeComponentBlackHole             ), intent(inout) :: blackHole1            , blackHole2
    double precision                                     , parameter     :: A               =1.2d4, B             =-0.93d0, H                 =7.3d3, &
         &                                                                  K               =6.0d4
    double precision                                                     :: q                     , velocityMass          , velocityOrthogonal      , &
         &                                                                  velocityParallel
    double precision                                                     :: alpha1Orthogonal      , alpha1Parallel        , alpha2Orthogonal        , &
         &                                                                  alpha2Parallel
    double precision                                                     :: phi                   , theta
    !$GLC attributes unused :: self

    ! If either black hole has non-positive mass, recoil velocity must be zero.
    if (blackHole1%mass() <= 0.0d0 .or. blackHole2%mass() <= 0.0d0) then
       campanelli2007Velocity=0.0d0
       return
    end if
    ! Select angular coordinates at random for the spin of the black hole.
    phi  =     2.0d0*Pi*blackHole1%hostNode%hostTree%randomNumberGenerator_%uniformSample()
    theta=acos(2.0d0*   blackHole1%hostNode%hostTree%randomNumberGenerator_%uniformSample()-1.0d0)
    ! Compute the mass ratio of the two black holes.
    q=blackHole1%mass()/blackHole2%mass()
    ! Compute α (and components), the angular momentum of the black hole per unit mass. This is equal to the spin scalar
    ! (since we don't track the direction of spin).
    alpha1Orthogonal=0.0d0
    alpha1Parallel  =blackHole1%spin()
    alpha2Orthogonal=blackHole2%spin()*sin(theta)*cos(phi)
    alpha2Parallel  =blackHole2%spin()*cos(theta)
    ! Compute velocity kicks in each direction.
    velocityOrthogonal=H*(q**2      /(1.0d0+q)**5)*(alpha2Parallel  -q*alpha1Parallel  )
    velocityParallel  =K*(q**2      /(1.0d0+q)**5)*(alpha2Orthogonal-q*alpha1Orthogonal)
    velocityMass      =A*(q**2*(1-q)/(1.0d0+q)**5)*(1.0d0+B*(q/(1.0d0+q)**2))
    ! Compute the net recoil velocity.
    campanelli2007Velocity=sqrt(velocityParallel**2+velocityMass**2+velocityOrthogonal**2)
    return
  end function campanelli2007Velocity

