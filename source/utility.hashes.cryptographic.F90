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

! Specify an explicit dependence on the utility.hashes.cryptographic.md5.o object file.
!: ./work/build/utility.hashes.cryptographic.md5.o

module Hashes_Cryptographic
  use, intrinsic :: ISO_C_Binding
  private
  public :: Hash_MD5

  interface
     subroutine md5(textLength,text,hash) bind(c,name='md5')
       !% Template for a C function that returns the MD5 of the input.
       import
       integer  (kind=c_int ), value :: textLength
       character(kind=c_char)        :: hash      (35        )
       character(kind=c_char)        :: text      (textLength)
     end subroutine md5
  end interface

contains

  function Hash_MD5(text)
    use ISO_Varying_String
    use String_Handling
    implicit none
    type     (varying_string)                               :: Hash_MD5
    type     (varying_string), intent(in   )                :: text
    character(kind=c_char   ), allocatable  , dimension(:)  :: textC
    character(kind=c_char   )               , dimension(35) :: hash
    integer  (kind=c_int    )                               :: textLen
    integer                                                 :: i

    textLen=len(text)+1
    allocate(textC(textLen))
    do i=1,textLen-1
       textC(i)=extract(text,i,i)
    end do
    textC(textLen)=char(0)
    call md5(textLen,textC,hash)
    deallocate(textC)
    do i=13,34
       if (hash(i) == "/") hash(i)="@"
    end do
    Hash_MD5=String_C_to_Fortran(hash(13:34))
    return
  end function Hash_MD5

end module Hashes_Cryptographic
