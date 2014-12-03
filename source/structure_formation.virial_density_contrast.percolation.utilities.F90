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

!% Contains a module of utilities needed by the {\tt percolation} virial density contrast class.

module Virial_Density_Contrast_Percolation_Utilities
  !% Provides utilities needed by the {\tt percolation} virial density contrast class.
  private
  public :: Virial_Density_Contrast_Percolation_Solver

contains

  !# <functionGlobal>
  !#  <unitName>Virial_Density_Contrast_Percolation_Solver</unitName>
  !#  <type>double precision</type>
  !#  <arguments>double precision , intent(in   )          :: mass, time, linkingLength</arguments>
  !#  <arguments>double precision , intent(in   ), pointer :: densityContrastCurrent</arguments>
  !# </functionGlobal>
  double precision function Virial_Density_Contrast_Percolation_Solver(mass,time,linkingLength,densityContrastCurrent)
    !% Return the virial density contrast at the given epoch, based on the percolation algorithm of \cite{more_overdensity_2011}.
    use Dark_Matter_Profiles
    use Cosmology_Parameters
    use Cosmology_Functions
    use Dark_Matter_Profiles_Concentration
    use Root_Finder
    use Galacticus_Nodes
    use Galacticus_Error
    use Numerical_Constants_Math
    use Galacticus_Calculations_Resets
    implicit none
    double precision                                     , intent(in   )          :: mass                           , time     , &
         &                                                                           linkingLength
    double precision                                     , intent(in   ), pointer :: densityContrastCurrent
    double precision                                     , parameter              :: percolationThreshold=0.652960d0
    class           (cosmologyParametersClass           ), pointer                :: cosmologyParameters_
    class           (cosmologyFunctionsClass            ), pointer                :: cosmologyFunctions_
    class           (darkMatterProfileConcentrationClass), pointer                :: darkMatterProfileConcentration_
    class           (darkMatterProfileClass             ), pointer                :: darkMatterProfile_
    type            (treeNode                           ), pointer                :: workNode
    class           (nodeComponentBasic                 ), pointer                :: workBasic
    class           (nodeComponentDarkMatterProfile     ), pointer                :: workDarkMatterProfile
    double precision                                                              :: boundingDensity                , radiusHalo, &
         &                                                                           timeActual
    logical                                                                       :: collapsingActual
    type            (rootFinder                         )                         :: finder

    ! Get required objects.
    cosmologyFunctions_             => cosmologyFunctions            ()
    cosmologyParameters_            => cosmologyParameters           ()
    darkMatterProfileConcentration_ => darkMatterProfileConcentration()
    darkMatterProfile_              => darkMatterProfile             ()
    ! Compute the bounding density, based on percolation theory (eq. 5 of More et al.).
    boundingDensity=+cosmologyFunctions_%matterDensityEpochal(time) &
         &          *percolationThreshold                           &
         &          /linkingLength                         **3
    ! Create a node and set the mass and time.
    workNode              => treeNode                  (                 )
    workBasic             => workNode%basic            (autoCreate=.true.)
    workDarkMatterProfile => workNode%darkMatterProfile(autoCreate=.true.)
    call workBasic%massSet            (mass)
    call workBasic%timeSet            (time)
    call workBasic%timeLastIsolatedSet(time)
    call Galacticus_Calculations_Reset(workNode)
    ! Make an initial guess at the halo radius.
    radiusHalo=(mass/4.0d0/Pi/boundingDensity)**(1.0d0/3.0d0)
    ! Find the corresponding halo radius.
    call finder   %tolerance          (                                               &
         &                             toleranceRelative  =1.0d-3                     &
         &                            )
    call finder   %rangeExpand        (                                               &
         &                             rangeExpandUpward  =2.0d0                    , &
         &                             rangeExpandDownward=0.5d0                    , &
         &                             rangeExpandType    =rangeExpandMultiplicative  &
         &                            )
    call finder   %rootFunction       (                                               &
         &                                                 haloRadiusRootFunction     &
         &                            )
    radiusHalo=finder%find(rootGuess=radiusHalo)
    call workNode%destroy()
    deallocate(workNode)
    ! Compute the corresponding density contrast.
    Virial_Density_Contrast_Percolation_Solver=+3.0d0                                             &
         &                                     *mass                                              &
         &                                     /4.0d0                                             &
         &                                     /Pi                                                &
         &                                     /radiusHalo                                    **3 &
         &                                     /cosmologyFunctions_%matterDensityEpochal(time)
    return

  contains
    
    double precision function haloRadiusRootFunction(haloRadiusTrial)
      !% Root function used to find the radius of a halo giving the correct bounding density.
      use Dark_Matter_Profile_Scales
      implicit none
      double precision, intent(in   ) :: haloRadiusTrial
      double precision                :: scaleRadius         , densityHaloRadius

      ! Construct the current density contrast.
      densityContrastCurrent=+3.0d0                                             &
           &                 *mass                                              &
           &                 /4.0d0                                             &
           &                 /Pi                                                &
           &                 /haloRadiusTrial                               **3 &
           &                 /cosmologyFunctions_%matterDensityEpochal(time)
      ! Find scale radius of the halo.
      scaleRadius=Dark_Matter_Profile_Scale(workNode,darkMatterProfileConcentration_)
      call workDarkMatterProfile%scaleSet(scaleRadius)
      call Galacticus_Calculations_Reset(workNode)
      ! Compute density at the halo radius.
      densityHaloRadius=darkMatterProfile_%density(workNode,haloRadiusTrial)
      ! Find difference from target density.
      haloRadiusRootFunction=boundingDensity-densityHaloRadius
      return
    end function haloRadiusRootFunction

  end function Virial_Density_Contrast_Percolation_Solver

end module Virial_Density_Contrast_Percolation_Utilities
