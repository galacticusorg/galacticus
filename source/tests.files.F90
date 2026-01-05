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
Contains a program to test file functions.
!!}

program Test_Files
  !!{
  Tests that file functions work.
  !!}
  use :: Display           , only : displayVerbositySet, verbosityLevelStandard
  use :: File_Utilities    , only : File_Exists        , File_Rename           , File_Remove         , File_Name_Expand
  use :: Input_Paths       , only : inputPath          , pathTypeExec
  use :: ISO_Varying_String, only : operator(//)       , var_str
  use :: Unit_Tests        , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer :: i
  
  call displayVerbositySet   (verbosityLevelStandard                                                                       )
  call Unit_Tests_Begin_Group("File utilities"                                                                             )
  ! File existence.
  call Assert                ('file exists'        ,File_Exists(inputPath(pathTypeExec)//'source/tests.files.F90' ),.true. )
  call Assert                ('file does not exist',File_Exists(inputPath(pathTypeExec)//'source/tests.bork.crump'),.false.)
  ! File renaming.
  open(newunit=i,file='tmp.file',status='unknown',form='formatted')
  write (i,*) "test file"
  close(i)
  call File_Rename(var_str('tmp.file'),var_str('mvd.file'))
  call Assert                ('file rename' ,File_Exists('mvd.file') .and. .not.File_Exists('tmp.file'),.true. )
  ! File removal.
  call File_Remove('mvd.file')
  call Assert                ('file removal',                                   File_Exists('tmp.file'),.false.)
  ! File name expansion.
  call Assert                ('file name expansion',File_Name_Expand('/opt%PARAM_PATH%/myParameters.xml'),var_str('/opt/home/galacticus/super/secret/stuff/myParameters.xml'))
  ! Finish.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Files
