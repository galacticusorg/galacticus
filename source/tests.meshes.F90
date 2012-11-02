!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a program to test the array functions.

program Test_Meshes
  !% Test mesh functions.
  use, intrinsic :: ISO_C_Binding
  use Unit_Tests
  use Meshes
  implicit none
  complex(c_double_complex), dimension(10,10,10) :: mesh,meshExpected

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
