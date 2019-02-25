!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module of utilities needed by the {\normalfont \ttfamily percolation} virial density contrast class.

module Virial_Density_Contrast_Percolation_Utilities
  !% Provides utilities needed by the {\normalfont \ttfamily percolation} virial density contrast class.
  use Galacticus_Nodes          , only : treeNode                                 , nodeComponentDarkMatterProfile
  use Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius             , darkMatterProfileScaleRadiusClass, &
       &                                 darkMatterProfileScaleRadiusConcentration
  use Dark_Matter_Profiles      , only : darkMatterProfile                        , darkMatterProfileClass
  private
  public :: Virial_Density_Contrast_Percolation_Solver

  ! Module-scope variables used in root finding.
  type            (treeNode                                 ), pointer :: workNode
  class           (nodeComponentDarkMatterProfile           ), pointer :: workDarkMatterProfile
  class           (darkMatterProfileClass                   ), pointer :: darkMatterProfile_
  type            (darkMatterProfileScaleRadiusConcentration), pointer :: darkMatterProfileScaleRadius_
  double precision                                                     :: boundingDensity              , densityMatterMean, &
       &                                                                  massHalo
  double precision                                           , pointer :: densityContrast
  !$omp threadprivate(workNode,workDarkMatterProfile,boundingDensity,densityMatterMean,massHalo,densityContrast,darkMatterProfileScaleRadius_,darkMatterProfile_)

contains

  !# <functionGlobal>
  !#  <unitName>Virial_Density_Contrast_Percolation_Solver</unitName>
  !#  <type>double precision</type>
  !#  <arguments>double precision , intent(in   )         :: mass, time, linkingLength</arguments>
  !#  <arguments>double precision , intent(in   ), target :: densityContrastCurrent</arguments>
  !# </functionGlobal>
  double precision function Virial_Density_Contrast_Percolation_Solver(mass,time,linkingLength,densityContrastCurrent)
    !% Return the virial density contrast at the given epoch, based on the percolation algorithm of \cite{more_overdensity_2011}.
    use Cosmology_Parameters
    use Cosmology_Functions
    use Root_Finder
    use Galacticus_Error
    use Numerical_Constants_Math
    use Galacticus_Calculations_Resets
    use Dark_Matter_Halo_Scales           , only : darkMatterHaloScale           , darkMatterHaloScaleClass
    use Dark_Matter_Profiles_Concentration, only : darkMatterProfileConcentration, darkMatterProfileConcentrationClass
    use Virial_Density_Contrast           , only : virialDensityContrast         , virialDensityContrastClass
    use Galacticus_Nodes                  , only : nodeComponentBasic
    implicit none
    double precision                                     , intent(in   )         :: mass                                      , time, &
         &                                                                          linkingLength
    double precision                                     , intent(in   ), target :: densityContrastCurrent
    double precision                                     , parameter             :: percolationThreshold           =0.652960d0
    class           (cosmologyParametersClass           ), pointer               :: cosmologyParameters_
    class           (cosmologyFunctionsClass            ), pointer               :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass           ), pointer               :: darkMatterHaloScale_
    class           (darkMatterProfileConcentrationClass), pointer               :: darkMatterProfileConcentration_
    class           (virialDensityContrastClass         ), pointer               :: virialDensityContrast_
    class           (nodeComponentBasic                 ), pointer               :: workBasic
    type            (rootFinder                         )                        :: finder
    double precision                                                             :: radiusHalo

    ! Initialize module-scope variables.
    massHalo        =  mass
    densityContrast => densityContrastCurrent
    ! Get required objects.
    cosmologyFunctions_             => cosmologyFunctions            ()
    cosmologyParameters_            => cosmologyParameters           ()
    darkMatterProfile_              => darkMatterProfile             ()
    darkMatterHaloScale_            => darkMatterHaloScale           ()
    darkMatterProfileConcentration_ => darkMatterProfileConcentration()
    virialDensityContrast_          => virialDensityContrast         ()
    ! Build a scale radius object.
    allocate(darkMatterProfileScaleRadius_)
    !# <referenceConstruct object="darkMatterProfileScaleRadius_">
    !#  <constructor>
    !#   darkMatterProfileScaleRadiusConcentration(                                                                   &amp;
    !#    &amp;                                    correctForConcentrationDefinition=.true.                         , &amp;
    !#    &amp;                                    useMeanConcentration             =.true.                         , &amp;
    !#    &amp;                                    cosmologyParameters_             =cosmologyParameters_           , &amp;
    !#    &amp;                                    cosmologyFunctions_              =cosmologyFunctions_            , &amp;
    !#    &amp;                                    darkMatterHaloScale_             =darkMatterHaloScale_           , &amp;
    !#    &amp;                                    darkMatterProfile_               =darkMatterProfile_             , &amp;
    !#    &amp;                                    virialDensityContrast_           =virialDensityContrast_         , &amp;
    !#    &amp;                                    darkMatterProfileConcentration_  =darkMatterProfileConcentration_  &amp;
    !#    &amp;                                   )
    !#  </constructor>
    !# </referenceConstruct>
    ! Compute the bounding density, based on percolation theory (eq. 5 of More et al.).
    densityMatterMean=cosmologyFunctions_%matterDensityEpochal(time)
    boundingDensity=+densityMatterMean       &
         &          *percolationThreshold    &
         &          /linkingLength       **3
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
    !# <objectDestructor name="darkMatterProfileScaleRadius_"/>
    ! Compute the corresponding density contrast.    
    Virial_Density_Contrast_Percolation_Solver=+3.0d0                &
         &                                     *mass                 &
         &                                     /4.0d0                &
         &                                     /Pi                   &
         &                                     /radiusHalo       **3 &
         &                                     /densityMatterMean
    return
  end function Virial_Density_Contrast_Percolation_Solver

  double precision function haloRadiusRootFunction(haloRadiusTrial)
    !% Root function used to find the radius of a halo giving the correct bounding density.
    use Dark_Matter_Profiles_Shape
    use Numerical_Constants_Math
    use Galacticus_Calculations_Resets
    implicit none
    double precision                             , intent(in   ) :: haloRadiusTrial
    double precision                                             :: scaleRadius            , densityHaloRadius
    class           (darkMatterProfileShapeClass), pointer       :: darkMatterProfileShape_

    ! Construct the current density contrast.
    densityContrast=+3.0d0                &
         &          *massHalo             &
         &          /4.0d0                &
         &          /Pi                   &
         &          /haloRadiusTrial  **3 &
         &          /densityMatterMean
    ! Find scale radius of the halo.
    scaleRadius=darkMatterProfileScaleRadius_%radius(workNode)
    call workDarkMatterProfile%scaleSet(scaleRadius)
    if (workDarkMatterProfile%shapeIsSettable()) then
       darkMatterProfileShape_ => darkMatterProfileShape()
       call workDarkMatterProfile%shapeSet(darkMatterProfileShape_%shape(workNode))
    end if
    call Galacticus_Calculations_Reset(workNode)
    ! Compute density at the halo radius.
    densityHaloRadius=darkMatterProfile_%density(workNode,haloRadiusTrial)
    ! Find difference from target density.
    haloRadiusRootFunction=boundingDensity-densityHaloRadius
    return
  end function haloRadiusRootFunction

end module Virial_Density_Contrast_Percolation_Utilities
