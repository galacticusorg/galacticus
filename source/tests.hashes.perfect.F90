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


!% Contains a program to test perfect hashing algorithms.

program Test_Perfect_Hashes
  !% Tests perfect hashing algorithms.
  use Unit_Tests
  use Hashes_Perfect
  use Kind_Numbers
  use Memory_Management
  implicit none
  integer,                parameter                        :: keyCount=11
  integer(kind=kind_int8)                                  :: i
  integer(kind=kind_int8),             dimension(keyCount) :: keys,values,retrievedValues
  integer,                allocatable, dimension(:)        :: bucketCount
  logical,                             dimension(keyCount) :: keyPresent
  type(hashPerfect)                                        :: hash

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.hashes.perfect.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Perfect hashes")

  ! Create a list of keys.
  keys  =[6,28,65,14,88,184,38,44,523,12,98]

  ! Create a list of values.
  values=[99,88,77,66,55,44,33,22,11,86,34]

  ! Create the hash function.
  call hash%create(keys,values)

  ! Allocate arrays.
  call Alloc_Array(bucketCount,int([hash%size()]),lowerBounds=[0])

  ! Look up indices and presence.
  bucketCount=0
  do i=1,keyCount
     keyPresent(i)=hash%isPresent(keys(i))
     bucketCount(hash%index(keys(i)))=bucketCount(hash%index(keys(i)))+1
     retrievedValues(i)=hash%value(keys(i))
  end do

  ! Make assertions.
  call Assert('All keys present in hash',all(keyPresent                             ),.true.)
  call Assert('No hash collisions'      ,all(bucketCount >= 0 .and. bucketCount <= 1),.true.)
  call Assert('All values correct'      ,retrievedValues                             ,values)

  ! Destroy the hash.
  call hash%destroy()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Perfect_Hashes
