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
Contains a program to test SI prefix functions.
!!}

program Test_SI_Prefixes
  !!{
  Tests that Roman numeral conversion functions work.
  !!}
  use :: Display                     , only : displayVerbositySet, verbosityLevelStandard
  use :: ISO_Varying_String          , only : var_str
  use :: Numerical_Constants_Prefixes, only : siFormat
  use :: Unit_Tests                  , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none

  call displayVerbositySet(verbosityLevelStandard                                                    )
  call Unit_Tests_Begin_Group        ("SI prefix formatting"                                           )
  call Assert                        ("1.23d-32 →   '0.01 q'"       ,trim(adjustl(siFormat(1.23d-32,'f7.2,1x'))),  '0.01 q')
  call Assert                        ("1.23d-05 →  '12.30 μ'"       ,trim(adjustl(siFormat(1.23d-05,'f7.2,1x'))), '12.30 μ')
  call Assert                        ("1.23d-02 →  '12.30 m'"       ,trim(adjustl(siFormat(1.23d-02,'f7.2,1x'))), '12.30 m')
  call Assert                        ("1.23d+01 →  '12.30'  "       ,trim(adjustl(siFormat(1.23d+01,'f7.2,1x'))), '12.30  '  )
  call Assert                        ("1.23d+04 →  '12.30 k'"       ,trim(adjustl(siFormat(1.23d+04,'f7.2,1x'))), '12.30 k')
  call Assert                        ("1.23d+07 →  '12.30 M'"       ,trim(adjustl(siFormat(1.23d+07,'f7.2,1x'))), '12.30 M')
  call Assert                        ("1.23d+32 → '123.00 Q'"       ,trim(adjustl(siFormat(1.23d+32,'f7.2,1x'))),'123.00 Q')
  call Unit_Tests_End_Group          (                                                                     )
  call Unit_Tests_Finish             (                                                                     )
end program Test_SI_Prefixes
