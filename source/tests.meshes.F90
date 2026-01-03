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
Contains a program to test the array functions.
!!}

program Test_Meshes
  !!{
  Test mesh functions.
  !!}
  use            :: Display      , only : displayVerbositySet, verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding, only : c_double_complex
  use            :: Meshes       , only : Meshes_Apply_Point , cloudTypeCubic        , cloudTypePoint      , cloudTypeTriangular
  use            :: Unit_Tests   , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  complex(c_double_complex), dimension(10,10,10) :: mesh,meshExpected

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Mesh functions")

  ! Apply a point to the mesh.
  mesh        =0.0d0
  meshExpected=0.0d0
  call Meshes_Apply_Point(mesh,10.0d0,[3.4d0,4.3d0,8.1d0],pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypePoint)
  meshExpected(4,5,9)=cmplx(1.0d0,0.0d0,kind=c_double_complex)
  call Assert('Apply point to mesh',reshape(real(mesh),[1000]),reshape(real(meshExpected),[1000]),absTol=1.0d-2)

  ! Apply a cubic cloud to the mesh.
  mesh        =0.0d0
  meshExpected=0.0d0
  call Meshes_Apply_Point(mesh,10.0d0,[3.4d0,4.3d0,8.1d0],pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeCubic)
  meshExpected(4,5,9)=cmplx(0.9d0*0.8d0*0.6d0,0.0d0,kind=c_double_complex)
  meshExpected(3,5,9)=cmplx(0.1d0*0.8d0*0.6d0,0.0d0,kind=c_double_complex)
  meshExpected(4,4,9)=cmplx(0.9d0*0.2d0*0.6d0,0.0d0,kind=c_double_complex)
  meshExpected(4,5,8)=cmplx(0.9d0*0.8d0*0.4d0,0.0d0,kind=c_double_complex)
  meshExpected(3,4,9)=cmplx(0.1d0*0.2d0*0.6d0,0.0d0,kind=c_double_complex)
  meshExpected(3,5,8)=cmplx(0.1d0*0.8d0*0.4d0,0.0d0,kind=c_double_complex)
  meshExpected(4,4,8)=cmplx(0.9d0*0.2d0*0.4d0,0.0d0,kind=c_double_complex)
  meshExpected(3,4,8)=cmplx(0.1d0*0.2d0*0.4d0,0.0d0,kind=c_double_complex)
  call Assert('Apply cubic cloud to mesh',reshape(real(mesh),[1000]),reshape(real(meshExpected),[1000]),absTol=1.0d-2)

  ! Apply a triangular cloud to the mesh.
  mesh        =0.0d0
  meshExpected=0.0d0
  call Meshes_Apply_Point(mesh,10.0d0,[3.4d0,4.3d0,8.1d0],pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
  meshExpected(4,5,9)=cmplx(0.98d0*0.92d0*0.68d0,0.0d0,kind=c_double_complex)
  meshExpected(3,5,9)=cmplx(0.02d0*0.92d0*0.68d0,0.0d0,kind=c_double_complex)
  meshExpected(4,4,9)=cmplx(0.98d0*0.08d0*0.68d0,0.0d0,kind=c_double_complex)
  meshExpected(4,5,8)=cmplx(0.98d0*0.92d0*0.32d0,0.0d0,kind=c_double_complex)
  meshExpected(3,4,9)=cmplx(0.02d0*0.08d0*0.68d0,0.0d0,kind=c_double_complex)
  meshExpected(3,5,8)=cmplx(0.02d0*0.92d0*0.32d0,0.0d0,kind=c_double_complex)
  meshExpected(4,4,8)=cmplx(0.98d0*0.08d0*0.32d0,0.0d0,kind=c_double_complex)
  meshExpected(3,4,8)=cmplx(0.02d0*0.08d0*0.32d0,0.0d0,kind=c_double_complex)
  call Assert('Apply triangular cloud to mesh',reshape(real(mesh),[1000]),reshape(real(meshExpected),[1000]),absTol=1.0d-2)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Meshes
