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
  Implementation of a kinematic distribution class for the soliton-NFW mass distribution.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionSolitonNFW">
    <description>
      A kinematic distribution class for the soliton-NFW mass distribution. In the NFW region, the velocity dispersion is computed
      using the analytic solution of \cite{lokas_properties_2001}, plus a correction term which accounts for the difference in
      mass inside the soliton radius between the soliton and the NFW profile. Specifically, we solve the Jeans equation for a
      point mass (since outside of the soliton radius, the mass within it acts as a point mass) and the NFW density profile:
      \begin{equation}
      \frac{\mathrm{d}(\rho\sigma^2}{\mathrm{d}r} = - \frac{\mathrm{G} \Delta M}{r^2} \rho, 
      \end{equation}
      where $\Delta M = M_\mathrm{soliton}(r_\mathrm{soliton}) - M_\mathrm{NFW}(r_\mathrm{soliton})$, and which has the solution:
      \begin{equation}
      \rho \sigma^2 = - \frac{\mathrm{G} \Delta M \rho_\mathrm{s}}{r_\mathrm{s}} \left( \frac{4}{r/r_\mathrm{s}} - \frac{1}{(r/r_\mathrm{s})^2} + \frac{2}{1+r/r_\mathrm{s}} + 6 \log \left[\frac{r}{r+r_\mathrm{s}}\right] \right).
      \end{equation}
      Inside the soliton, the NFW solution is used as a boundary condition at $r_\mathrm{soliton}$, and the Jeans equation is then
      solved using the soliton density profile, which results in a solution
      \begin{equation}
       \frac{\pi \mathrm{G} r_\mathrm{c}^2 \rho_\mathrm{c}^2}{20038287360 a^{3/2}}                                                
       \left[                                                  
       -\sqrt{a} (169995 + 631540 y + 1200199 y^2 + 1317888 y^3 + 849849 y^4 + 300300 y^5 + 45045 y^6) (28672 + 169995 y + 631540 y^2 + 1200199 y^3 + 1317888 y^4 + 849849 y^5 + 300300 y^6 + 45045 y^7)/(1+y)^{14}
       +45045 \tan^{-1}(\sqrt{a} x) \left\{ -2 ( 14336 + 169995 y + 631540 y^2 + 1200199 y^3 + 1317888 y^4 + 849849 y^5 + 300300 y^6 + 45045 y^7)/[x (1+y)^7]      
       -45045 \sqrt{a} \tan^{-1}(\sqrt{a} x) \right\}
       \right],
      \end{equation}
      where $x = r/r_\mathrm{c}$ and $y = a x^2$.
    </description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionCollisionlessTabulated) :: kinematicsDistributionSolitonNFW
     !!{
     A kinematics distribution for the Soliton-NFW mass distribution.
     !!}
   contains
     procedure :: velocityDispersion1D => solitonNFWKinematicsVelocityDispersion1D
  end type kinematicsDistributionSolitonNFW

  interface kinematicsDistributionSolitonNFW
     !!{
     Constructors for the \refClass{kinematicsDistributionSolitonNFW} kinematic distribution class.
     !!}
     module procedure solitonNFWKinematicsConstructorParameters
     module procedure solitonNFWKinematicsConstructorInternal
     module procedure solitonNFWKinematicsConstructorDecorated
  end interface kinematicsDistributionSolitonNFW

  ! Coefficient of the dimensionless radius in the soliton profile.
   double precision, parameter :: coefficientCore=0.091d0 ! Schive et al. (2014; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S; equation 3).

contains

  function solitonNFWKinematicsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionSolitonNFW} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionSolitonNFW)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    double precision                                                  :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum

    !![
    <inputParameter>
      <name>toleranceRelativeVelocityDispersion</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersionMaximum</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    !!]
    self=kinematicsDistributionSolitonNFW(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function solitonNFWKinematicsConstructorParameters

  function solitonNFWKinematicsConstructorInternal(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionSolitonNFW} kinematic distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionSolitonNFW)                          :: self
    double precision                                  , intent(in   ), optional :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum"/>
    !!]

    return
  end function solitonNFWKinematicsConstructorInternal
  
  function solitonNFWKinematicsConstructorDecorated(kinematicsDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionSolitonNFW} kinematic distribution class.
    !!}
    implicit none
    type (kinematicsDistributionSolitonNFW)                :: self
    class(kinematicsDistributionClass     ), intent(in   ) :: kinematicsDistribution_

    self%toleranceRelativeVelocityDispersion       =kinematicsDistribution_%toleranceRelativeVelocityDispersion
    self%toleranceRelativeVelocityDispersionMaximum=kinematicsDistribution_%toleranceRelativeVelocityDispersionMaximum
    return
  end function solitonNFWKinematicsConstructorDecorated
  
  logical function solitonNFWKinematicsIsCollisional(self)
    !!{
    Return false indicating that the soliton-NFW distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionSolitonNFW), intent(inout) :: self
    !$GLC attributes unused :: self
    
    solitonNFWKinematicsIsCollisional=.false.
    return
  end function solitonNFWKinematicsIsCollisional

  double precision function solitonNFWKinematicsVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in a soliton-NFW kinematic distribution.
    !!}
    use :: Error                           , only : Error_Report
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (kinematicsDistributionSolitonNFW), intent(inout)          :: self
    class           (coordinate                      ), intent(in   )          :: coordinates
    class           (massDistributionClass           ), intent(inout), target  :: massDistribution_ , massDistributionEmbedding
    class           (massDistributionClass           )               , pointer :: massDistribution__
    double precision                                                           :: radius

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating NFW distribution we have an analytic solution for the velocity dispersion.
       select type (massDistributionEmbedding)
       class is (massDistributionSolitonNFW)
          ! Determine if the radius is within the soliton.
          radius=+coordinates%rSpherical()
          if (radius < massDistributionEmbedding%radiusSoliton) then
             ! In the soliton region.
             !! Evaluate ρσ² at the soliton radius as our boundary condition.
             velocityDispersion=+velocityDispersionSquareNFW(massDistributionEmbedding%radiusSoliton) &
                  &             *massDistributionEmbedding%densityNormalizationNFW                    &
                  &             /        massDistributionEmbedding%radiusSolitonScaleFree             &
                  &             /(+1.0d0+massDistributionEmbedding%radiusSolitonScaleFree)**2
             ! Add on the solution for ρσ² for the soliton from the boundary to the current radius.
             velocityDispersion=+velocityDispersion                                            &
                  &             +jeansIntegralSoliton(massDistributionEmbedding%radiusSoliton) &
                  &             -jeansIntegralSoliton(                          radius       )
             ! Convert to velocity dispersion.
             velocityDispersion=+sqrt(                                                                              &
                  &                   +velocityDispersion                                                           &
                  &                   /massDistributionEmbedding%densitySolitonCentral                              &
                  &                   *(+1.0d0+coefficientCore*(radius/massDistributionEmbedding%radiusCore)**2)**8 &
                  &                  )
          else
             ! In the NFW region.
             velocityDispersion=sqrt(velocityDispersionSquareNFW(radius))
          end if
       class default
          velocityDispersion=+0.0d0
          call Error_Report('expecting a soliton-NFW mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our tabulated distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion=+self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return

  contains

    double precision function jeansIntegralSoliton(radius)
      !!{
      Compute the Jeans integral in the soliton region of a soliton-NFW profile.
      !!}
      implicit none
      double precision, intent(in   ) :: radius
      double precision                :: radiusScaleFree, radiusFactor

      select type (massDistributionEmbedding)
      class is (massDistributionSolitonNFW)
         radiusScaleFree     =+                          radius                   &
              &               /massDistributionEmbedding%radiusCore
         radiusFactor        =+coefficientCore                                    &
              &               *radiusScaleFree**2
         jeansIntegralSoliton=+Pi                                                 &
              &               *gravitationalConstant_internal                     &
              &               *massDistributionEmbedding%radiusCore           **2 &
              &               *massDistributionEmbedding%densitySolitonCentral**2 &
              &               /2.003828736d10                                     &
              &               /coefficientCore**1.5d0                             &
              &               *(                                                  &
              &                 -sqrt(coefficientCore)                            &
              &                 *(                                                &
              &                   + 169995.0d0                                    &
              &                   + 631540.0d0*radiusFactor                       &
              &                   +1200199.0d0*radiusFactor**2                    &
              &                   +1317888.0d0*radiusFactor**3                    &
              &                   + 849849.0d0*radiusFactor**4                    &
              &                   + 300300.0d0*radiusFactor**5                    &
              &                   +  45045.0d0*radiusFactor**6                    &
              &                  )                                                &
              &                 *(                                                &
              &                   +  28672.0d0                                    &
              &                   + 169995.0d0*radiusFactor                       &
              &                   + 631540.0d0*radiusFactor**2                    &
              &                   +1200199.0d0*radiusFactor**3                    &
              &                   +1317888.0d0*radiusFactor**4                    &
              &                   + 849849.0d0*radiusFactor**5                    &
              &                   + 300300.0d0*radiusFactor**6                    &
              &                   +  45045.0d0*radiusFactor**7                    &
              &                  )                                                &
              &                 /(1.0d0+radiusFactor)**14                         &
              &                 +45045.0d0                                        &
              &                 *atan(sqrt(coefficientCore)*radiusScaleFree)      &
              &                 *(                                                &
              &                   -2.0d0                                          &
              &                   *(                                              &
              &                     +  14336.0d0                                  &
              &                     + 169995.0d0*radiusFactor                     &
              &                     + 631540.0d0*radiusFactor**2                  &
              &                     +1200199.0d0*radiusFactor**3                  &
              &                     +1317888.0d0*radiusFactor**4                  &
              &                     + 849849.0d0*radiusFactor**5                  &
              &                     + 300300.0d0*radiusFactor**6                  &
              &                     +  45045.0d0*radiusFactor**7                  &
              &                    )                                              &
              &                   /(radiusScaleFree*(1.0d0+radiusFactor)**7)      &
              &                   -45045.0d0                                      &
              &                   *     sqrt(coefficientCore)                     &
              &                   *atan(sqrt(coefficientCore)*radiusScaleFree)    &
              &                  )                                                &
              &                )
      class default
         jeansIntegralSoliton=0.0d0
         call Error_Report('expecting a soliton-NFW mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
      end select
      return
    end function jeansIntegralSoliton
         
    double precision function velocityDispersionSquareNFW(radius)
      !!{
      Compute the square of the velocity dispersion in the NFW region of a soliton-NFW profile.
      !!}
      use :: Dilogarithms, only : Dilogarithm
      implicit none
      double precision, intent(in   ) :: radius
      double precision, parameter     :: nfwMinimumRadiusForExactSolution=1.0d-2
      double precision, parameter     :: nfwMaximumRadiusForExactSolution=1.0d+2    
      double precision, parameter     :: nfwNormalizationFactorUnitRadius=-8.5d0+Pi**2-6.0d0*log(2.0d0)+6.0d0*log(2.0d0)**2 ! Precomputed NFW normalization factor for unit radius.
      double precision                :: radiusScaleFree                                                                   , logRadiusScaleFree, &
           &                             onePlusRadiusScaleFree                                                            , logOnePlusRadius

      select type (massDistributionEmbedding)
      class is (massDistributionSolitonNFW)
         radiusScaleFree=+                          radius      &
              &          /massDistributionEmbedding%radiusScale
         if (radiusScaleFree == 1.0d0) then
            velocityDispersionSquareNFW=nfwNormalizationFactorUnitRadius
         else if (radiusScaleFree >= nfwMaximumRadiusForExactSolution) then
            logRadiusScaleFree         =+log(radiusScaleFree)
            velocityDispersionSquareNFW=+(-   3.0d0+   4.0d0*logRadiusScaleFree)/(    16.0d0*radiusScaleFree   ) &
                 &                      +(   69.0d0+  20.0d0*logRadiusScaleFree)/(   200.0d0*radiusScaleFree**2) &
                 &                      +(-  97.0d0-  60.0d0*logRadiusScaleFree)/(  1200.0d0*radiusScaleFree**3) &
                 &                      +(   71.0d0+ 105.0d0*logRadiusScaleFree)/(  3675.0d0*radiusScaleFree**4) &
                 &                      +(-   1.0d0-  56.0d0*logRadiusScaleFree)/(  3136.0d0*radiusScaleFree**5) &
                 &                      +(-1271.0d0+2520.0d0*logRadiusScaleFree)/(211680.0d0*radiusScaleFree**6)
         else if (radiusScaleFree >= nfwMinimumRadiusForExactSolution) then
            onePlusRadiusScaleFree     =+     1.0d0+radiusScaleFree
            logRadiusScaleFree         =+log(       radiusScaleFree)
            logOnePlusRadius           =+log(onePlusRadiusScaleFree)
            velocityDispersionSquareNFW=+0.5d0                               &
                 &                      *       radiusScaleFree              &
                 &                      *onePlusRadiusScaleFree**2           &
                 &                      *(                                   &
                 &                        +Pi**2                             &
                 &                        -logRadiusScaleFree                &
                 &                        -1.0d0/       radiusScaleFree      &
                 &                        -1.0d0/onePlusRadiusScaleFree**2   &
                 &                        -6.0d0/onePlusRadiusScaleFree      &
                 &                        +(                                 &
                 &                          +1.0d0+ 1.0d0/radiusScaleFree**2 &
                 &                                - 4.0d0/radiusScaleFree    &
                 &                          -2.0d0/onePlusRadiusScaleFree    &
                 &                         )                                 &
                 &                        *logOnePlusRadius                  &
                 &                        +3.0d0                             &
                 &                        *logOnePlusRadius**2               &
                 &                        +6.0d0                             &
                 &                        *Dilogarithm(-radiusScaleFree)     &
                 &                       )
         else if (radiusScaleFree > 0.0d0) then
            logRadiusScaleFree         =+log(radiusScaleFree)
            velocityDispersionSquareNFW=+ 1.0d0/   4.0d0*(-23.0d0      + 2.0d0*Pi**2- 2.0d0*logRadiusScaleFree)*radiusScaleFree    &
                 &                      +                (-59.0d0/6.0d0+       Pi**2-       logRadiusScaleFree)*radiusScaleFree**2 &
                 &                      + 1.0d0/  24.0d0*(-101.0d0     +12.0d0*Pi**2-12.0d0*logRadiusScaleFree)*radiusScaleFree**3 &
                 &                      +11.0d0/  60.0d0                                                       *radiusScaleFree**4 &
                 &                      -13.0d0/ 240.0d0                                                       *radiusScaleFree**5 &
                 &                      +37.0d0/1400.0d0                                                       *radiusScaleFree**6
         else
            velocityDispersionSquareNFW=0.0d0
         end if
         velocityDispersionSquareNFW=+velocityDispersionSquareNFW                          &
              &                      *4.0d0                                                &
              &                      *Pi                                                   &
              &                      *gravitationalConstant_internal                       &
              &                      *massDistributionEmbedding%densityNormalizationNFW    &
              &                      *massDistributionEmbedding%radiusScale            **2
         ! Add on the correction term for the difference in the mass of the soliton and the mass of the NFW profile within the soliton radius.
         if (radiusScaleFree > nfwMaximumRadiusForExactSolution) then
            ! For very large radii use a series solution.
            velocityDispersionSquareNFW=+velocityDispersionSquareNFW                        &
                 &                      -gravitationalConstant_internal                     &
                 &                      *(                                                  &
                 &                        -massDistributionEmbedding%massSolitionTransition &
                 &                        +massDistributionEmbedding%massNFWTransition      &
                 &                       )                                                  &
                 &                      /  massDistributionEmbedding%radiusScale            &
                 &                      *(                                                  &
                 &                        +1.0d0/ 4.0d0/radiusScaleFree                     &
                 &                        +1.0d0/10.0d0/radiusScaleFree**2                  &
                 &                        +1.0d0/20.0d0/radiusScaleFree**3                  &
                 &                        +1.0d0/35.0d0/radiusScaleFree**4                  &
                 &                      )
         else
            ! For other radii use the exact solution.
            velocityDispersionSquareNFW=+velocityDispersionSquareNFW                        &
                 &                      +gravitationalConstant_internal                     &
                 &                      *(                                                  &
                 &                        -massDistributionEmbedding%massSolitionTransition &
                 &                        +massDistributionEmbedding%massNFWTransition      &
                 &                       )                                                  &
                 &                      /  massDistributionEmbedding%radiusScale            &
                 &                      /2.0d0                                              &
                 &                      *(                                                  &
                 &                        +4.0d0/             radiusScaleFree               &
                 &                        -1.0d0/             radiusScaleFree**2            &
                 &                        +2.0d0/     (+1.0d0+radiusScaleFree   )           &
                 &                        +6.0d0*log(                                       &
                 &                                   +        radiusScaleFree               &
                 &                                   /(+1.0d0+radiusScaleFree   )           &
                 &                                  )                                       &
                 &                       )                                                  &
                 &                      *(                                                  &
                 &                        +                   radiusScaleFree               &
                 &                        *           (+1.0d0+radiusScaleFree   )**2        &
                 &                       )
         end if
      class default
         velocityDispersionSquareNFW=0.0d0
         call Error_Report('expecting a soliton-NFW mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
      end select
      return
    end function velocityDispersionSquareNFW
    
  end function solitonNFWKinematicsVelocityDispersion1D
