!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  An abstract dark matter halo profile class for SIDM profiles.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <darkMatterProfile name="darkMatterProfileSIDM" abstract="yes">
    <description>
      Abstract dark matter halo profile for self-interacting dark matter particles. Provides a method to compute the
      characteristic radius for interactions, $r_1$, defined as (e.g. Jiang et al. 2020):
      \begin{equation}
      \frac{4}{\sqrt{\pi}} \rho_\mathrm{dm}(r_1) v_\mathrm{rms}(r_1) \frac{\sigma}{m} = \frac{1}{t_\mathrm{age}},
      \end{equation}      
      where the left-hand side is the scattering rate per particle, with $\rho_\mathrm{dm}(r)$ being the dark matter density
      profile, $v_\mathrm{rms}$ the average relative velocity between DM particles (which is approximated by the 1D velocity
      dispersion), and $\sigma/m$ is the self-interaction cross-section per unit mass.
    </description>
  </darkMatterProfile>
  !!]
  type, abstract, extends(darkMatterProfileClass) :: darkMatterProfileSIDM
     !!{     
     An abstract dark matter halo profile class for SIDM profiles.
     !!}
     private
     class(darkMatterProfileClass ), pointer :: darkMatterProfile_  => null()
     class(darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
   contains
     !![
     <methods>
       <method method="radiusInteraction" description="Computes the characteristic interaction radius of the halo."/>
     </methods>
     !!]
     procedure :: radiusInteraction => sidmRadiusInteraction
  end type darkMatterProfileSIDM

  ! Submodule-scope variables used in root finding.
  class           (darkMatterProfileSIDM), pointer :: self_
  type            (treeNode             ), pointer :: node_
  double precision                                 :: timeAge_, crossSection_
  !$omp threadprivate(self_,node_,timeAge_,crossSection_)
  
contains
  
  double precision function sidmRadiusInteraction(self,node,timeAge)
    !!{
    Returns the characteristic interaction radius (in Mpc) of the self-interacting dark matter profile of {\normalfont \ttfamily node}.
    !!}
    use :: Dark_Matter_Particles           , only : darkMatterParticleSelfInteractingDarkMatter
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: Numerical_Constants_Prefixes    , only : centi                                      , kilo
    use :: Numerical_Constants_Astronomical, only : megaParsec                                 , massSolar
    use :: Root_Finder                     , only : rootFinder                                 , rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    class           (darkMatterProfileSIDM), intent(inout), target   :: self
    type            (treeNode             ), intent(inout), target   :: node
    double precision                       , intent(in   ), optional :: timeAge
    class           (nodeComponentBasic   )               , pointer  :: basic
    type            (rootFinder           )                          :: finder
    double precision                       , parameter               :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    
    self_         =>  self
    node_         =>  node
    if (present(timeAge)) then
       timeAge_   =         timeAge
    else
       basic      =>  node %basic  ()
       timeAge_   =   basic%time   ()
    end if
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       crossSection_=+darkMatterParticle_%crossSectionSelfInteraction() &
            &        *centi     **2                                     &
            &        /megaParsec**2                                     &
            &        *kilo                                              &
            &        *massSolar
    class default
       call Error_Report('expected self-interacting dark matter particle'//{introspection:location})
    end select
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
    sidmRadiusInteraction=finder%find(rootGuess=self%darkMatterHaloScale_%radiusVirial(node))
    return
  end function sidmRadiusInteraction

  double precision function sidmRadiusInteractionRoot(radius)
    !!{
    Root function used in seeking the characteristic interaction radius in self-interacting dark matter profiles.
    !!}
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    sidmRadiusInteractionRoot=+4.0d0                                                           &
         &                    /sqrt(Pi)                                                        &
         &                    /Mpc_per_km_per_s_To_Gyr                                         &
         &                    *self_%darkMatterProfile_%density                 (node_,radius) &
         &                    *self_%darkMatterProfile_%radialVelocityDispersion(node_,radius) &
         &                    *crossSection_                                                   &
         &                    -1.0d0                                                           &
         &                    /timeAge_
    return
  end function sidmRadiusInteractionRoot
  
