!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a program to test array search functions.

program Test_Search
  !% Tests that array search functions work.
  use Unit_Tests
  use Arrays_Search
  use ISO_Varying_String
  implicit none
  double precision,     dimension(10) :: myArray =[0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0]
  double precision,     dimension(10) :: mySearch=[3.4d0,9.0d0,4.2d0,-1.0d0,10.0d0,5.5d0,5.999999d0,6.000001d0,1.1d0,7.5d0]
  integer,              dimension(10) :: myIndices
  type(varying_string), dimension(26) :: stringArray
  integer                         :: i 

  ! Define an array of varying strings.
  stringArray=[            &
       &       'Alpha   ', &
       &       'Bravo   ', &
       &       'Charlie ', &
       &       'Delta   ', &
       &       'Echo    ', &
       &       'Foxtrot ', &
       &       'Golf    ', &
       &       'Hotel   ', &
       &       'India   ', &
       &       'Juliet  ', &
       &       'Kilo    ', &
       &       'Lima    ', &
       &       'Mike    ', &
       &       'November', &
       &       'Oscar   ', &
       &       'Papa    ', &
       &       'Quebec  ', &
       &       'Romeo   ', &
       &       'Sierra  ', &
       &       'Tango   ', &
       &       'Uniform ', &
       &       'Victor  ', &
       &       'Whiskey ', &
       &       'X-ray   ', &
       &       'Yankee  ', &
       &       'Zulu    '  &
       &      ]


  ! Begin unit tests.
  call Unit_Tests_Begin_Group("array search")

  ! Test searching of double arrays.
  do i=1,size(mySearch)
     myIndices(i)=Search_Array(myArray,mySearch(i))
  end do
  call Assert('search double array',myIndices,max(min(int(mySearch)+1,9),1))
  do i=1,size(mySearch)
     myIndices(i)=Search_Array_For_Closest(myArray,mySearch(i))
  end do
  call Assert('search double array for closest match',myIndices,max(min(nint(mySearch)+1,10),1))

  ! Test searching of varying string arrays.
  call Assert('search string array',[                                               &
       &                             Search_Array(stringArray,var_str('Monkey'  )), &
       &                             Search_Array(stringArray,var_str('Unicorn' )), &
       &                             Search_Array(stringArray,var_str('Beaver'  )), &
       &                             Search_Array(stringArray,var_str('Fox'     )), &
       &                             Search_Array(stringArray,var_str('Oscar'   )), &
       &                             Search_Array(stringArray,var_str('Shoggoth'))  &
       &                            ],                                              &
       &                            [                                               &
       &                             13,                                            &
       &                             20,                                            &
       &                              1,                                            &
       &                              5,                                            &
       &                             15,                                            &
       &                             18                                             &
       &                            ]                                               &
       &     )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Search
