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

! Specify an explicit dependence on the md5.o and sha256.o object files.
!: $(BUILDPATH)/external/MD5/md5.o
!: $(BUILDPATH)/external/SHA256/sha256.o

module Hashes_Cryptographic
  use, intrinsic :: ISO_C_Binding, only : c_char, c_int
  private
  public :: Hash_MD5, Hash_SHA256_File

  interface
     subroutine md5Hash(textLength,input,result) bind(c,name='md5Hash')
       !!{RST
       Template for a C function that returns the MD5 hash of the input.
       !!}
       import
       integer  (kind=c_int ), value :: textLength
       character(kind=c_char)        :: input     (textLength)
       character(kind=c_char)        :: result    (        33)
     end subroutine md5Hash
  end interface

  interface
     function sha256File(fileName,result) bind(c,name='sha256File')
       !!{RST
       Template for a C function that returns the SHA-256 hash of the contents of the named file. The function result is zero on
       success, and non-zero if the file could not be read.
       !!}
       import
       integer  (kind=c_int ) :: sha256File
       character(kind=c_char) :: fileName  (*)
       character(kind=c_char) :: result    (65)
     end function sha256File
  end interface

contains

  function Hash_MD5(text)
    use :: ISO_Varying_String, only : assignment(=)      , extract, len, operator(//), &
          &                           varying_string
    use :: String_Handling   , only : String_C_to_Fortran, char
    implicit none
    type     (varying_string)                               :: Hash_MD5
    type     (varying_string), intent(in   )                :: text
    character(kind=c_char   ), allocatable  , dimension(: ) :: textC
    character(kind=c_char   )               , dimension(33) :: hash
    integer  (kind=c_int    )                               :: textLen
    integer                                                 :: i

    textLen=len(text)+1
    allocate(textC(textLen))
    do i=1,textLen-1
       textC(i)=extract(text,i,i)
    end do
    textC(textLen)=char(0)
    call md5Hash(textLen,textC,hash)
    deallocate(textC)  
    Hash_MD5=String_C_to_Fortran(hash)
    return
  end function Hash_MD5

  function Hash_SHA256_File(fileName,status)
    !!{RST
    Return the SHA-256 hash of the contents of the file ``fileName``, as a string of 64 lowercase hexadecimal digits.

    If the file can not be read then an empty string is returned, and ``status`` (if present) is set non-zero. If ``status`` is
    not present, failure to read the file is fatal---a caller verifying the integrity of a file must not be allowed to silently
    treat an unreadable file as a match.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=)      , extract, len, varying_string, &
         &                            operator(//)
    use :: String_Handling   , only : String_C_to_Fortran, char
    implicit none
    type     (varying_string)                               :: Hash_SHA256_File
    type     (varying_string), intent(in   )                :: fileName
    integer                  , intent(  out), optional      :: status
    character(kind=c_char   ), allocatable  , dimension(: ) :: fileNameC
    character(kind=c_char   )               , dimension(65) :: hash
    integer  (kind=c_int    )                               :: status_
    integer                                                 :: fileNameLength  , i

    fileNameLength=len(fileName)+1
    allocate(fileNameC(fileNameLength))
    do i=1,fileNameLength-1
       fileNameC(i)=extract(fileName,i,i)
    end do
    fileNameC(fileNameLength)=char(0)
    status_=sha256File(fileNameC,hash)
    deallocate(fileNameC)
    if (present(status)) status=status_
    if (status_ == 0) then
       Hash_SHA256_File=String_C_to_Fortran(hash)
    else
       Hash_SHA256_File=''
       if (.not.present(status)) call Error_Report('unable to read file "'//fileName//'"'//{introspection:location})
    end if
    return
  end function Hash_SHA256_File

end module Hashes_Cryptographic
