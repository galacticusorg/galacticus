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
Contains a program to test the \cite{inoue_updated_2014} algorithm for IGM absorption.
!!}

program Test_Inoue2014
  !!{
  Tests the \cite{inoue_updated_2014} algorithm for IGM absorption.
  !!}
  use :: Display                               , only : displayVerbositySet                           , verbosityLevelStandard
  use :: Stellar_Population_Spectra_Postprocess, only : stellarPopulationSpectraPostprocessorInoue2014
  use :: Unit_Tests                            , only : Assert                                        , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
          &                                             compareLessThanOrEqual
  implicit none
  type            (stellarPopulationSpectraPostprocessorInoue2014)               :: postprocessor
  integer                                                                        :: inoueUnit        , ioStatus        , i
  double precision                                                               :: wavelength       , relativeError   , relativeErrorMaximum
  double precision                                                , dimension(7) :: transmissionInoue, transmissionSelf
  double precision                                                , dimension(7) :: redshift=[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0]

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Inoue et al. (2014) IGM attenuation model")
  ! Compare our calculation of attenuation with that from Inoue's own code and record the maximum error. The file read below was
  ! extracted from the tarball containing Inoue's code (where it was originally called "test") downloaded from:
  !   http://www.las.osaka-sandai.ac.jp/~inoue/ANAIGM/ANAIGM.tar.gz
  postprocessor=stellarPopulationSpectraPostprocessorInoue2014()
  relativeErrorMaximum=0.0d0
  open(newunit=inoueUnit,file="testSuite/data/inoue2014igmAttenuationModel.txt",status='old',form='formatted',iostat=ioStatus)
  do while (ioStatus == 0)
     read (inoueUnit,*,iostat=ioStatus) wavelength,transmissionInoue
     transmissionSelf=1.0d0
     do i=1,7
        transmissionSelf(i)=postprocessor%multiplier(wavelength,0.0d0,redshift(i))
        relativeError=abs(transmissionSelf(i)-transmissionInoue(i))/abs(transmissionInoue(i))
        relativeErrorMaximum=max(relativeError,relativeErrorMaximum)
     end do
  end do
  close(inoueUnit)
  call Assert('Maximum error less than 1%',relativeErrorMaximum,1.0d-2,compareLessThanOrEqual)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Test_Inoue2014
