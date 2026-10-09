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

!+    Contributions to this file made by: Andrew Robertson, Codex.

program Test_Survey_Geometry_Martin2010_ALFALFA
  !!{RST
  Test the treatment of non-positive HI masses by the Martin et al. (2010) ALFALFA survey geometry.
  !!}
  use :: Cosmology_Parameters, only : cosmologyParametersSimple
  use :: Geometry_Surveys    , only : surveyGeometryMartin2010ALFALFA
  use :: Unit_Tests          , only : Assert, Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type(cosmologyParametersSimple       )         :: cosmologyParameters
  type(surveyGeometryMartin2010ALFALFA)         :: surveyGeometry

  call Unit_Tests_Begin_Group("Martin et al. (2010) ALFALFA survey geometry")
  cosmologyParameters=cosmologyParametersSimple(                 &
       &                                         OmegaMatter    =0.3d0    , &
       &                                         OmegaBaryon    =0.0455d0 , &
       &                                         OmegaDarkEnergy=0.7d0    , &
       &                                         temperatureCMB =2.72548d0, &
       &                                         HubbleConstant =70.0d0     &
       &                                        )
  surveyGeometry=surveyGeometryMartin2010ALFALFA(cosmologyParameters)

  call Assert("zero HI mass has zero survey depth"    ,surveyGeometry%distanceMaximum(mass= 0.0d0),0.0d0,absTol=0.0d0)
  call Assert("negative HI mass has zero survey depth",surveyGeometry%distanceMaximum(mass=-1.0d0),0.0d0,absTol=0.0d0)
  call Assert("positive HI mass has positive depth"   ,surveyGeometry%distanceMaximum(mass= 1.0d9) > 0.0d0,.true.)

  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Survey_Geometry_Martin2010_ALFALFA
