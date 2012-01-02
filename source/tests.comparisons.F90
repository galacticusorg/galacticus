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


!% Contains a program to test numerical comparison functions.

program Test_Comparison
  !% Tests that numerical comparison functions work.
  use Unit_Tests
  use Numerical_Comparison
  implicit none
 
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("numerical comparison")

  ! Check results.
  call Assert('double values agree' ,[                                                           &
       &                              Values_Agree ( 1.0d0, 1.000d0                            ), &
       &                              Values_Agree ( 1.0d0, 2.000d0                            ), &
       &                              Values_Agree ( 1.0d0, 1.005d0,absTol=0.01d0              ), &
       &                              Values_Agree ( 1.0d0, 1.015d0,absTol=0.01d0              ), &
       &                              Values_Agree (10.0d0,10.050d0,              relTol=0.01d0), &
       &                              Values_Agree (10.0d0,10.150d0,              relTol=0.01d0), &
       &                              Values_Agree (10.0d0,10.050d0,absTol=0.10d0,relTol=0.01d0), &
       &                              Values_Agree (10.0d0,10.150d0,absTol=0.20d0,relTol=0.01d0), &
       &                              Values_Agree (10.0d0,10.150d0,absTol=0.10d0,relTol=0.02d0), &
       &                              Values_Agree (10.0d0,10.250d0,absTol=0.10d0,relTol=0.02d0)  &
       &                             ],                                                           &
       &                             [                                                            &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .false.                                                     &
       &                             ]                                                            &
       &     )
  call Assert('double values differ',[                                                            &
       &                              Values_Differ( 1.0d0, 1.000d0                            ), &
       &                              Values_Differ( 1.0d0, 2.000d0                            ), &
       &                              Values_Differ( 1.0d0, 1.005d0,absTol=0.01d0              ), &
       &                              Values_Differ( 1.0d0, 1.015d0,absTol=0.01d0              ), &
       &                              Values_Differ(10.0d0,10.050d0,              relTol=0.01d0), &
       &                              Values_Differ(10.0d0,10.150d0,              relTol=0.01d0), &
       &                              Values_Differ(10.0d0,10.050d0,absTol=0.10d0,relTol=0.01d0), &
       &                              Values_Differ(10.0d0,10.150d0,absTol=0.20d0,relTol=0.01d0), &
       &                              Values_Differ(10.0d0,10.150d0,absTol=0.10d0,relTol=0.02d0), &
       &                              Values_Differ(10.0d0,10.250d0,absTol=0.10d0,relTol=0.02d0)  &
       &                             ],                                                           &
       &                             [                                                            &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .true.                                                      &
       &                              ]                                                           &
       &     )
  call Assert('real values agree'   ,[                                                           &
       &                              Values_Agree ( 1.0e0, 1.000e0                            ), &
       &                              Values_Agree ( 1.0e0, 2.000e0                            ), &
       &                              Values_Agree ( 1.0e0, 1.005e0,absTol=0.01e0              ), &
       &                              Values_Agree ( 1.0e0, 1.015e0,absTol=0.01e0              ), &
       &                              Values_Agree (10.0e0,10.050e0,              relTol=0.01e0), &
       &                              Values_Agree (10.0e0,10.150e0,              relTol=0.01e0), &
       &                              Values_Agree (10.0e0,10.050e0,absTol=0.10e0,relTol=0.01e0), &
       &                              Values_Agree (10.0e0,10.150e0,absTol=0.20e0,relTol=0.01e0), &
       &                              Values_Agree (10.0e0,10.150e0,absTol=0.10e0,relTol=0.02e0), &
       &                              Values_Agree (10.0e0,10.250e0,absTol=0.10e0,relTol=0.02e0)  &
       &                             ],                                                           &
       &                             [                                                            &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .false.                                                     &
       &                             ]                                                            &
       &     )
  call Assert('real values differ'  ,[                                                            &
       &                              Values_Differ( 1.0e0, 1.000e0                            ), &
       &                              Values_Differ( 1.0e0, 2.000e0                            ), &
       &                              Values_Differ( 1.0e0, 1.005e0,absTol=0.01e0              ), &
       &                              Values_Differ( 1.0e0, 1.015e0,absTol=0.01e0              ), &
       &                              Values_Differ(10.0e0,10.050e0,              relTol=0.01e0), &
       &                              Values_Differ(10.0e0,10.150e0,              relTol=0.01e0), &
       &                              Values_Differ(10.0e0,10.050e0,absTol=0.10e0,relTol=0.01e0), &
       &                              Values_Differ(10.0e0,10.150e0,absTol=0.20e0,relTol=0.01e0), &
       &                              Values_Differ(10.0e0,10.150e0,absTol=0.10e0,relTol=0.02e0), &
       &                              Values_Differ(10.0e0,10.250e0,absTol=0.10e0,relTol=0.02e0)  &
       &                             ],                                                           &
       &                             [                                                            &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .true.                                                      &
       &                              ]                                                           &
       &     )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Comparison
