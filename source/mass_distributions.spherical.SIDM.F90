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
  Implements an abstract spherical mass distribution class for SIDM models.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass
  
  !![
  <massDistribution name="massDistributionSphericalSIDM" abstract="yes">
   <description>
     An abstract mass distribution class for spherical SIDM models. Provides a method to compute interaction radii.
   </description>
  </massDistribution>
  !!]
  type, abstract, extends(massDistributionSphericalDecorator) :: massDistributionSphericalSIDM
     !!{
     Implementation of a spherical mass distribution for SIDM models.
     !!}
     private
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
     double precision                                   :: timeAge
  contains
     !![
     <methods>
       <method method="radiusInteraction" description="Computes the characteristic interaction radius of the halo."/>
     </methods>
     !!]
     procedure :: radiusInteraction => sidmRadiusInteraction
  end type massDistributionSphericalSIDM

  ! Submodule-scope variables used in root finding.
  class           (massDistributionSphericalSIDM), pointer :: self_
  class           (kinematicsDistributionClass  ), pointer :: kinematicsDistribution_
  double precision                                         :: crossSection_
  !$omp threadprivate(self_,kinematicsDistribution_,crossSection_)

contains

  double precision function sidmRadiusInteraction(self) result(radiusInteraction)
    !!{
    Returns the characteristic interaction radius (in Mpc) of the self-interacting dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Dark_Matter_Particles           , only : darkMatterParticleSelfInteractingDarkMatter
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Prefixes    , only : centi                                      , kilo
    use :: Numerical_Constants_Astronomical, only : megaParsec                                 , massSolar
    use :: Root_Finder                     , only : rootFinder                                 , rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: ISO_Varying_String              , only : char
    implicit none
    class           (massDistributionSphericalSIDM), intent(inout), target :: self
    type            (rootFinder                   ), save                  :: finder
    logical                                        , save                  :: finderInitialized=.false.
    double precision                               , parameter             :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-3
    !$omp threadprivate(finder,finderInitialized)

    self_ => self
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       crossSection_=+darkMatterParticle_%crossSectionSelfInteraction() &
            &        *centi     **2                                     &
            &        /megaParsec**2                                     &
            &        *kilo                                              &
            &        *massSolar
    class default
       call Error_Report('expected self-interacting dark matter particle but found type "'//char(darkMatterParticle_%objectType())//'"'//{introspection:location})
    end select
    if (.not.finderInitialized) then
       finder=rootFinder(                                             &
            &            rootFunction     =sidmRadiusInteractionRoot, &
            &            toleranceAbsolute=toleranceAbsolute        , &
            &            toleranceRelative=toleranceRelative          &
            &           )
       call finder%rangeExpand(                                                             &
            &                  rangeExpandUpward            =2.0d0                        , &
            &                  rangeExpandDownward          =0.5d0                        , &
            &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                  rangeExpandType              =rangeExpandMultiplicative      &
            &                 )
       finderInitialized=.true.
    end if
    kinematicsDistribution_ => self  %massDistribution_%kinematicsDistribution(               )
    radiusInteraction       =  finder                  %find                  (rootGuess=1.0d0)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function sidmRadiusInteraction

  double precision function sidmRadiusInteractionRoot(radius) result(residual)
    !!{
    Root function used in seeking the characteristic interaction radius in self-interacting dark matter profiles.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    use :: Coordinates                     , only : coordinateSpherical, assignment(=)
    implicit none
    double precision                     , intent(in   ) :: radius
    type            (coordinateSpherical)                :: coordinates
    
    coordinates= [radius,0.0d0,0.0d0]
    residual   =+4.0d0                                                                                                           &
         &      /sqrt(Pi)                                                                                                        &
         &      /MpcPerKmPerSToGyr                                                                                               &
         &      *self_%massDistribution_      %density             (coordinates                                                ) &
         &      *      kinematicsDistribution_%velocityDispersion1D(coordinates,self_%massDistribution_,self_%massDistribution_) &
         &      *crossSection_                                                                                                   &
         &      -1.0d0                                                                                                           &
         &      /self_%timeAge
    return
  end function sidmRadiusInteractionRoot
