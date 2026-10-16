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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests the computational domain volume integrators.
!!}

program Test_Computational_Domain_Volume_Integrators
  !!{RST
  Tests the computational domain volume integrators. Each integrator is checked against analytically-known integrals over its own
  coordinate system, for a smooth integrand and for the discontinuous integrand produced by a constant density cloud whose surface
  cuts through the cell---the case which arises in radiative transfer models.
  !!}
  use :: Computational_Domain_Volume_Integrators                , only : computationalDomainVolumeIntegratorCartesian3D, computationalDomainVolumeIntegratorCylindrical, &
       &                                                                 computationalDomainVolumeIntegratorSpherical
  use :: Display                                                , only : displayMessage                                , displayVerbositySet                           , &
       &                                                                 verbosityLevelStandard
  use :: ISO_Varying_String                                     , only : varying_string                                , var_str                                       , &
       &                                                                 operator(//)                                  , assignment(=)
  use :: Numerical_Constants_Math                               , only : Pi
  use :: String_Handling                                        , only : operator(//)
  use :: Test_Computational_Domain_Volume_Integrators_Functions , only : countEvaluations                              , countEvaluationsReset                         , &
       &                                                                 integrandCartesian                            , integrandCloud                                , &
       &                                                                 integrandCylindrical                          , integrandSpherical                            , &
       &                                                                 integrandUnity                                , radiusCloud
  use :: Unit_Tests                                             , only : Assert                                        , Unit_Tests_Begin_Group                        , &
       &                                                                 Unit_Tests_End_Group                          , Unit_Tests_Finish
  implicit none
  type            (computationalDomainVolumeIntegratorCartesian3D)                 :: integratorCartesian3D
  type            (computationalDomainVolumeIntegratorSpherical  )                 :: integratorSpherical
  type            (computationalDomainVolumeIntegratorCylindrical)                 :: integratorCylindrical
  double precision                                                , dimension(3,2) :: boundariesCartesian3D
  double precision                                                , dimension(2,2) :: boundariesCylindrical
  double precision                                                , dimension(2  ) :: boundariesSpherical
  double precision                                                                 :: integral             , integralAbsolute

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Computational domain volume integrators")

  ! Cartesian 3D cells.
  call Unit_Tests_Begin_Group("Cartesian 3D")
  !! A unit integrand must integrate to the volume of the cell.
  boundariesCartesian3D(1,:)=[0.0d0,1.0d0]
  boundariesCartesian3D(2,:)=[0.0d0,1.0d0]
  boundariesCartesian3D(3,:)=[0.0d0,1.0d0]
  integratorCartesian3D     =computationalDomainVolumeIntegratorCartesian3D(boundariesCartesian3D)
  integral                  =reportEvaluations(integratorCartesian3D%integrate(integrandUnity),"unit integrand")
  call Assert("∫dV over 0<x<1, 0<y<1, 0<z<1 equals the cell volume",integral,integratorCartesian3D%volume(),relTol=1.0d-3)
  !! A smooth integrand.
  integral=reportEvaluations(integratorCartesian3D%integrate(integrandCartesian),"smooth integrand")
  call Assert("∫(x²+y²+z²)dV over 0<x<1, 0<y<1, 0<z<1",integral,1.0d0,relTol=1.0d-3)
  !! A discontinuous integrand, with the discontinuity cutting through the cell. The cell contains one octant of the cloud.
  integral=reportEvaluations(integratorCartesian3D%integrate(integrandCloud),"discontinuous integrand (one octant of cloud in cell)")
  call Assert("∫dV over the part of a unit-radius cloud lying in 0<x<1, 0<y<1, 0<z<1",integral,Pi/6.0d0,relTol=2.0d-2)
  !! A discontinuous integrand, with the cloud inscribed in the cell.
  boundariesCartesian3D(1,:)=[-1.0d0,+1.0d0]
  boundariesCartesian3D(2,:)=[-1.0d0,+1.0d0]
  boundariesCartesian3D(3,:)=[-1.0d0,+1.0d0]
  integratorCartesian3D     =computationalDomainVolumeIntegratorCartesian3D(boundariesCartesian3D)
  integral                  =reportEvaluations(integratorCartesian3D%integrate(integrandCloud),"discontinuous integrand (cloud inscribed in cell)")
  call Assert("∫dV over a unit-radius cloud inscribed in -1<x<1, -1<y<1, -1<z<1",integral,4.0d0*Pi/3.0d0,relTol=2.0d-2)
  !! The same integral, but with an absolute tolerance set to 1% of the cell volume. Since the integrand here is a density, this
  !! asks for the mean density in the cell to 1%, which is the accuracy actually required---rather than 1% of whatever the integral
  !! happens to be, which becomes arbitrarily expensive for cells containing only a sliver of the cloud.
  integratorCartesian3D     =computationalDomainVolumeIntegratorCartesian3D(boundariesCartesian3D,toleranceAbsolute=1.0d-2*8.0d0)
  integral                  =reportEvaluations(integratorCartesian3D%integrate(integrandCloud),"discontinuous integrand (cloud inscribed in cell; absolute tolerance)")
  call Assert("∫dV over a unit-radius cloud inscribed in -1<x<1, -1<y<1, -1<z<1, to an absolute tolerance",integral,4.0d0*Pi/3.0d0,absTol=1.0d-2*8.0d0)
  !! A cell containing only a thin sliver of the cloud---here a spherical cap of height 0.01, of volume πh²(3R-h)/3. This is the
  !! case which arises for the ~24% of cells straddling the surface of the cloud in radiative transfer models, and the case for
  !! which a purely relative tolerance is pathological: it asks for 1% of a number which is small compared to the cell volume, and
  !! so compared to the accuracy actually needed.
  boundariesCartesian3D(1,:)=[+0.99d0,+1.00d0]
  boundariesCartesian3D(2,:)=[-1.00d0,+1.00d0]
  boundariesCartesian3D(3,:)=[-1.00d0,+1.00d0]
  integratorCartesian3D     =computationalDomainVolumeIntegratorCartesian3D(boundariesCartesian3D)
  integral                  =reportEvaluations(integratorCartesian3D%integrate(integrandCloud),"discontinuous integrand (sliver of cloud in cell)")
  call Assert("∫dV over the sliver of a unit-radius cloud lying in 0.99<x<1, -1<y<1, -1<z<1",integral,Pi*1.0d-4*(3.0d0-0.01d0)/3.0d0,relTol=2.0d-2)
  !! The same integral, but with an absolute tolerance set to 1% of the cell volume.
  integratorCartesian3D     =computationalDomainVolumeIntegratorCartesian3D(boundariesCartesian3D,toleranceAbsolute=1.0d-2*0.04d0)
  integral                  =reportEvaluations(integratorCartesian3D%integrate(integrandCloud),"discontinuous integrand (sliver of cloud in cell; absolute tolerance)")
  call Assert("∫dV over the sliver of a unit-radius cloud lying in 0.99<x<1, -1<y<1, -1<z<1, to an absolute tolerance",integral,Pi*1.0d-4*(3.0d0-0.01d0)/3.0d0,absTol=1.0d-2*0.04d0)
  !! A cell with the geometry actually arising in radiative transfer models: the cell is much smaller than the cloud, so the
  !! surface of the cloud cuts through it as an almost flat sheet. This is a far easier discontinuity than the cases above, in
  !! which the cell and the cloud are comparable in size, and is the regime in which the cost of populating a domain is set.
  radiusCloud               =6.0d-6
  boundariesCartesian3D(1,:)=[+5.7d-6,+6.3d-6]
  boundariesCartesian3D(2,:)=[+0.0d+0,+6.0d-7]
  boundariesCartesian3D(3,:)=[+0.0d+0,+6.0d-7]
  integratorCartesian3D     =computationalDomainVolumeIntegratorCartesian3D(boundariesCartesian3D)
  integral                  =reportEvaluations(integratorCartesian3D%integrate(integrandCloud),"discontinuous integrand (model-like cell; relative tolerance)")
  integratorCartesian3D     =computationalDomainVolumeIntegratorCartesian3D(boundariesCartesian3D,toleranceAbsolute=1.0d-2*2.16d-19)
  integralAbsolute          =reportEvaluations(integratorCartesian3D%integrate(integrandCloud),"discontinuous integrand (model-like cell; absolute tolerance)")
  call Assert("∫dV over a model-like cell straddling the cloud surface agrees between tolerance specifications",integralAbsolute,integral,absTol=1.0d-2*2.16d-19)
  radiusCloud               =1.0d0
  call Unit_Tests_End_Group()

  ! Spherical cells.
  call Unit_Tests_Begin_Group("spherical")
  !! A unit integrand must integrate to the volume of the cell.
  boundariesSpherical  =[1.0d0,2.0d0]
  integratorSpherical  =computationalDomainVolumeIntegratorSpherical(boundariesSpherical)
  integral             =reportEvaluations(integratorSpherical%integrate(integrandUnity),"unit integrand")
  call Assert("∫dV over 1<r<2 equals the cell volume",integral,integratorSpherical%volume(),relTol=1.0d-3)
  !! A smooth integrand.
  integral=reportEvaluations(integratorSpherical%integrate(integrandSpherical),"smooth integrand")
  call Assert("∫cos²(θ)dV over 1<r<2",integral,28.0d0*Pi/9.0d0,relTol=1.0d-3)
  !! A discontinuous integrand, with the discontinuity cutting through the cell.
  boundariesSpherical  =[0.0d0,2.0d0]
  integratorSpherical  =computationalDomainVolumeIntegratorSpherical(boundariesSpherical)
  integral             =reportEvaluations(integratorSpherical%integrate(integrandCloud),"discontinuous integrand")
  call Assert("∫dV over a unit-radius cloud contained in 0<r<2",integral,4.0d0*Pi/3.0d0,relTol=2.0d-2)
  call Unit_Tests_End_Group()

  ! Cylindrical cells.
  call Unit_Tests_Begin_Group("cylindrical")
  !! A unit integrand must integrate to the volume of the cell.
  boundariesCylindrical(1,:)=[0.0d0,1.0d0]
  boundariesCylindrical(2,:)=[0.0d0,2.0d0]
  integratorCylindrical     =computationalDomainVolumeIntegratorCylindrical(boundariesCylindrical)
  integral                  =reportEvaluations(integratorCylindrical%integrate(integrandUnity),"unit integrand")
  call Assert("∫dV over 0<r<1, 0<z<2 equals the cell volume",integral,integratorCylindrical%volume(),relTol=1.0d-3)
  !! A smooth integrand.
  integral=reportEvaluations(integratorCylindrical%integrate(integrandCylindrical),"smooth integrand")
  call Assert("∫rz dV over 0<r<1, 0<z<2",integral,4.0d0*Pi/3.0d0,relTol=1.0d-3)
  !! A discontinuous integrand, with the discontinuity cutting through the cell.
  boundariesCylindrical(1,:)=[+0.0d0,+2.0d0]
  boundariesCylindrical(2,:)=[-2.0d0,+2.0d0]
  integratorCylindrical     =computationalDomainVolumeIntegratorCylindrical(boundariesCylindrical)
  integral                  =reportEvaluations(integratorCylindrical%integrate(integrandCloud),"discontinuous integrand")
  call Assert("∫dV over a unit-radius cloud contained in 0<r<2, -2<z<2",integral,4.0d0*Pi/3.0d0,relTol=2.0d-2)
  !! The same integral, but with an absolute tolerance set to 1% of the cell volume.
  integratorCylindrical     =computationalDomainVolumeIntegratorCylindrical(boundariesCylindrical,toleranceAbsolute=1.0d-2*16.0d0*Pi)
  integral                  =reportEvaluations(integratorCylindrical%integrate(integrandCloud),"discontinuous integrand (absolute tolerance)")
  call Assert("∫dV over a unit-radius cloud contained in 0<r<2, -2<z<2, to an absolute tolerance",integral,4.0d0*Pi/3.0d0,absTol=1.0d-2*16.0d0*Pi)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  double precision function reportEvaluations(integral,label)
    !!{RST
    Report the number of integrand evaluations used to compute  integral, and reset the counter ready for the next integration.
    !!}
    implicit none
    double precision                , intent(in   ) :: integral
    character       (len=*         ), intent(in   ) :: label
    type            (varying_string)                :: message

    message  =var_str('  ')//label//': '//countEvaluations//' evaluations'
    call displayMessage(message,verbosityLevelStandard)
    call countEvaluationsReset()
    reportEvaluations=integral
    return
  end function reportEvaluations

end program Test_Computational_Domain_Volume_Integrators
