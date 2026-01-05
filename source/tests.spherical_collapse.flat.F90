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
Contains a program which tests spherical collapse calculations for a flat Universe.
!!}

program Tests_Spherical_Collapse_Flat
  !!{
  Tests spherical collapse calculations for a flat Universe. Compares results to the fitting formula of
  \cite{bryan_statistical_1998}.
  !!}
  use :: Cosmology_Functions     , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters    , only : cosmologyParametersSimple
  use :: Display                 , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                                                        , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Virial_Density_Contrast , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  implicit none
  double precision                                                               , dimension(7) :: redshift              =[0.0d0,1.0d0,3.0d0,7.0d0,15.0d0,31.0d0,63.0d0]
  type           (cosmologyParametersSimple                                     )               :: cosmologyParameters_
  type           (cosmologyFunctionsMatterLambda                                )               :: cosmologyFunctions_
  type           (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)               :: virialDensityContrast_
  double precision                                                               , parameter    :: mass                  =1.0d12
  character       (len=1024                                                     )               :: message
  integer                                                                                       :: iExpansion
  double precision                                                                              :: age                                                                  , bryanNormanFit , &
       &                                                                                           densityContrast                                                      , expansionFactor, &
       &                                                                                           x

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: flat cosmology")
  ! Test spherical collapse in a flat universe.
  !![
  <referenceConstruct object="cosmologyParameters_"  >
   <constructor>
    cosmologyParametersSimple                                     (                                           &amp;
     &amp;                                                         OmegaMatter         = 0.1d0              , &amp;
     &amp;                                                         OmegaBaryon         = 0.0d0              , &amp;
     &amp;                                                         OmegaDarkEnergy     = 0.9d0              , &amp;
     &amp;                                                         temperatureCMB      = 2.7d0              , &amp;
     &amp;                                                         HubbleConstant      =70.0d0                &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"   >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                           &amp;
     &amp;                                                         cosmologyParameters_=cosmologyParameters_  &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_">
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                           &amp;
     &amp;                                                         cosmologyFunctions_ =cosmologyFunctions_ , &amp;
     &amp;                                                         tableStore          =.true.                &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshift(iExpansion))
     age            =cosmologyFunctions_%cosmicTime(expansionFactor)
     densityContrast=virialDensityContrast_%densityContrast(mass,age)
     x              =cosmologyFunctions_%omegaMatterEpochal(age)-1.0d0
     bryanNormanFit =(18.0d0*Pi**2+82.0d0*x-39.0d0*x**2)/cosmologyFunctions_%omegaMatterEpochal(age)
     write (message,'(a,f6.1,a,f6.4,a)') "virial density contrast [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctions_%omegaMatterEpochal(age),"]"
     call Assert(trim(message),densityContrast,bryanNormanFit,relTol=1.1d-2)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Spherical_Collapse_Flat
