!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of a spectrum postprocessor that applies the \cite{inoue_updated_2014} calculation of the attenuation of spectra by the intergalactic medium.
  !!}

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorInoue2014">
   <description>
    A stellar population postprocessing class that postprocesses spectra through absorption by the \gls{igm} using the results
    of \cite{inoue_updated_2014}.
   </description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorInoue2014
     !!{
     A spectrum postprocessor applying the \cite{inoue_updated_2014} calculation of the attenuation of spectra by the intergalactic medium.
     !!}
     private
   contains
     procedure :: multiplier => inoue2014Multiplier
  end type stellarPopulationSpectraPostprocessorInoue2014

  interface stellarPopulationSpectraPostprocessorInoue2014
     !!{
     Constructors for the \refClass{stellarPopulationSpectraPostprocessorInoue2014} stellar population spectra postprocessor class.
     !!}
     module procedure inoue2014ConstructorParameters
  end interface stellarPopulationSpectraPostprocessorInoue2014

  ! Fitting function coefficients.
  double precision, dimension(3,2:40) :: aLAF=reshape(                                        &
       &                                              [                                       &
       &                                               1.68976d-02, 2.35379d-03, 1.02611d-04, &
       &                                               4.69229d-03, 6.53625d-04, 2.84940d-05, &
       &                                               2.23898d-03, 3.11884d-04, 1.35962d-05, &
       &                                               1.31901d-03, 1.83735d-04, 8.00974d-06, &
       &                                               8.70656d-04, 1.21280d-04, 5.28707d-06, &
       &                                               6.17843d-04, 8.60640d-05, 3.75186d-06, &
       &                                               4.60924d-04, 6.42055d-05, 2.79897d-06, &
       &                                               3.56887d-04, 4.97135d-05, 2.16720d-06, &
       &                                               2.84278d-04, 3.95992d-05, 1.72628d-06, &
       &                                               2.31771d-04, 3.22851d-05, 1.40743d-06, &
       &                                               1.92348d-04, 2.67936d-05, 1.16804d-06, &
       &                                               1.62155d-04, 2.25878d-05, 9.84689d-07, &
       &                                               1.38498d-04, 1.92925d-05, 8.41033d-07, &
       &                                               1.19611d-04, 1.66615d-05, 7.26340d-07, &
       &                                               1.04314d-04, 1.45306d-05, 6.33446d-07, &
       &                                               9.17397d-05, 1.27791d-05, 5.57091d-07, &
       &                                               8.12784d-05, 1.13219d-05, 4.93564d-07, &
       &                                               7.25069d-05, 1.01000d-05, 4.40299d-07, &
       &                                               6.50549d-05, 9.06198d-06, 3.95047d-07, &
       &                                               5.86816d-05, 8.17421d-06, 3.56345d-07, &
       &                                               5.31918d-05, 7.40949d-06, 3.23008d-07, &
       &                                               4.84261d-05, 6.74563d-06, 2.94068d-07, &
       &                                               4.42740d-05, 6.16726d-06, 2.68854d-07, &
       &                                               4.06311d-05, 5.65981d-06, 2.46733d-07, &
       &                                               3.73821d-05, 5.20723d-06, 2.27003d-07, &
       &                                               3.45377d-05, 4.81102d-06, 2.09731d-07, &
       &                                               3.19891d-05, 4.45601d-06, 1.94255d-07, &
       &                                               2.97110d-05, 4.13867d-06, 1.80421d-07, &
       &                                               2.76635d-05, 3.85346d-06, 1.67987d-07, &
       &                                               2.58178d-05, 3.59636d-06, 1.56779d-07, &
       &                                               2.41479d-05, 3.36374d-06, 1.46638d-07, &
       &                                               2.26347d-05, 3.15296d-06, 1.37450d-07, &
       &                                               2.12567d-05, 2.96100d-06, 1.29081d-07, &
       &                                               1.99967d-05, 2.78549d-06, 1.21430d-07, &
       &                                               1.88476d-05, 2.62543d-06, 1.14452d-07, &
       &                                               1.77928d-05, 2.47850d-06, 1.08047d-07, &
       &                                               1.68222d-05, 2.34330d-06, 1.02153d-07, &
       &                                               1.59286d-05, 2.21882d-06, 9.67268d-08, &
       &                                               1.50996d-05, 2.10334d-06, 9.16925d-08  &
       &                                              ]                                     , &
       &                                              [3,39]                                  &
       &                                             )
  double precision, dimension(2,2:40) :: aDLA=reshape(                                        &
       &                                              [                                       &
       &                                               1.61698d-04, 5.38995d-05,              &
       &                                               1.54539d-04, 5.15129d-05,              &
       &                                               1.49767d-04, 4.99222d-05,              &
       &                                               1.46031d-04, 4.86769d-05,              &
       &                                               1.42893d-04, 4.76312d-05,              &
       &                                               1.40159d-04, 4.67196d-05,              &
       &                                               1.37714d-04, 4.59048d-05,              &
       &                                               1.35495d-04, 4.51650d-05,              &
       &                                               1.33452d-04, 4.44841d-05,              &
       &                                               1.31561d-04, 4.38536d-05,              &
       &                                               1.29785d-04, 4.32617d-05,              &
       &                                               1.28117d-04, 4.27056d-05,              &
       &                                               1.26540d-04, 4.21799d-05,              &
       &                                               1.25041d-04, 4.16804d-05,              &
       &                                               1.23614d-04, 4.12046d-05,              &
       &                                               1.22248d-04, 4.07494d-05,              &
       &                                               1.20938d-04, 4.03127d-05,              &
       &                                               1.19681d-04, 3.98938d-05,              &
       &                                               1.18469d-04, 3.94896d-05,              &
       &                                               1.17298d-04, 3.90995d-05,              &
       &                                               1.16167d-04, 3.87225d-05,              &
       &                                               1.15071d-04, 3.83572d-05,              &
       &                                               1.14011d-04, 3.80037d-05,              &
       &                                               1.12983d-04, 3.76609d-05,              &
       &                                               1.11972d-04, 3.73241d-05,              &
       &                                               1.11002d-04, 3.70005d-05,              &
       &                                               1.10051d-04, 3.66836d-05,              &
       &                                               1.09125d-04, 3.63749d-05,              &
       &                                               1.08220d-04, 3.60734d-05,              &
       &                                               1.07337d-04, 3.57789d-05,              &
       &                                               1.06473d-04, 3.54909d-05,              &
       &                                               1.05629d-04, 3.52096d-05,              &
       &                                               1.04802d-04, 3.49340d-05,              &
       &                                               1.03991d-04, 3.46636d-05,              &
       &                                               1.03198d-04, 3.43994d-05,              &
       &                                               1.02420d-04, 3.41402d-05,              &
       &                                               1.01657d-04, 3.38856d-05,              &
       &                                               1.00908d-04, 3.36359d-05,              &
       &                                               1.00168d-04, 3.33895d-05               &
       &                                              ]                                     , &
       &                                              [2,39]                                  &
       &                                             )

contains

  function inoue2014ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSpectraPostprocessorInoue2014} stellar population spectra postprocessor class which takes a
    parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(stellarPopulationSpectraPostprocessorInoue2014)                :: self
    type(inputParameters                               ), intent(inout) :: parameters

    self=stellarPopulationSpectraPostprocessorInoue2014()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function inoue2014ConstructorParameters

  double precision function inoue2014Multiplier(self,wavelength,age,redshift)
    !!{
    Apply the \cite{inoue_updated_2014} calculation of the attenuation of spectra by the intergalactic medium.
    !!}
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    implicit none
    class           (stellarPopulationSpectraPostprocessorInoue2014), intent(inout) :: self
    double precision                                                , intent(in   ) :: age                      , redshift                        , &
         &                                                                             wavelength
    double precision                                                , parameter     :: redshiftZero       =0.0d0
    integer                                                                         :: i
    double precision                                                                :: opticalDepth             , wavelengthObservedLymanContinuum, &
         &                                                                             wavelengthLymanLine      , wavelengthScaled
    !$GLC attributes unused :: self, age

    ! Return if this is a zero redshift case.
    inoue2014Multiplier=1.0d0
    if (redshift <= 0.0d0) return
    ! Initialize optical depth to zero.
    opticalDepth=0.0d0
    ! Line absorption.
    do i=2,40
       ! Lyman-α forest.
       wavelengthLymanLine=lymanSeriesLimitWavelengthHydrogen_atomic/(1.0d0-1.0d0/dble(i**2))
       wavelengthScaled   =wavelength*(1.0d0+redshift)/wavelengthLymanLine
       if (wavelengthScaled < 1.0d0+redshiftZero .or. wavelengthScaled > 1.0d0+redshift) cycle
       if      (wavelengthScaled < 2.2d0) then
          opticalDepth=opticalDepth+aLAF(1,i)*(wavelengthScaled)**1.2d0
       else if (wavelengthScaled < 5.7d0) then
          opticalDepth=opticalDepth+aLAF(2,i)*(wavelengthScaled)**3.7d0
       else
          opticalDepth=opticalDepth+aLAF(3,i)*(wavelengthScaled)**5.5d0
       end if
       ! DLAs.
       if      (wavelengthScaled < 3.0d0) then
          opticalDepth=opticalDepth+aDLA(1,i)*(wavelengthScaled)**2.0d0
       else
          opticalDepth=opticalDepth+aDLA(2,i)*(wavelengthScaled)**3.0d0
       end if
    end do
    ! Compute the observed wavelength in units of the Lyman-continuum wavelength.
    wavelengthObservedLymanContinuum=wavelength*(1.0d0+redshift)/lymanSeriesLimitWavelengthHydrogen_atomic
    ! Add continuum absorption is wavelength is sufficiently short.
    if (wavelengthObservedLymanContinuum < 1.0d0+redshift) then
       ! Lyman-α forest continuum absorption.
       if      (redshift < 1.2d0) then
          if      (wavelengthObservedLymanContinuum < 1.0d0+redshift)                                         &
               & opticalDepth=+opticalDepth                                                                   &
               &              +0.3248d+0                                                                      &
               &              *(                                                                              &
               &                +                                   (wavelengthObservedLymanContinuum)**1.2d0 &
               &                -        (1.0d0+redshift)**(-0.9d0)*(wavelengthObservedLymanContinuum)**2.1d0 &
               &               )
       else if (redshift < 4.7d0) then
          if      (wavelengthObservedLymanContinuum < 2.2d0         ) then
             opticalDepth    =+opticalDepth                                                                   &
                  &           +2.5450d-2*(1.0d0+redshift)**(+1.6d0)*(wavelengthObservedLymanContinuum)**2.1d0 &
                  &           +0.3248d+0                           *(wavelengthObservedLymanContinuum)**1.2d0 &
                  &           -0.2496d+0                           *(wavelengthObservedLymanContinuum)**2.1d0
          else if (wavelengthObservedLymanContinuum< 1.0d0+redshift) then
             opticalDepth    =+opticalDepth                                                                   &
                  &           +2.5450d-2                                                                      &
                  &           *(                                                                              &
                  &             +        (1.0d0+redshift)**(+1.6d0)*(wavelengthObservedLymanContinuum)**2.1d0 &
                  &             -                                   (wavelengthObservedLymanContinuum)**3.7d0 &
                  &            )
          end if
       else
          if      (wavelengthObservedLymanContinuum< 2.2d0         ) then
             opticalDepth    =+opticalDepth                                                                   &
                  &           +5.2210d-4*(1.0d0+redshift)**(+3.4d0)*(wavelengthObservedLymanContinuum)**2.1d0 &
                  &           +0.3248d+0                           *(wavelengthObservedLymanContinuum)**1.2d0 &
                  &           -3.1400d-2                           *(wavelengthObservedLymanContinuum)**2.1d0
          else if (wavelengthObservedLymanContinuum< 5.7d0         ) then
             opticalDepth    =+opticalDepth                                                                   &
                  &           +5.2210d-4*(1.0d0+redshift)**(+3.4d0)*(wavelengthObservedLymanContinuum)**2.1d0 &
                  &           +0.2182d+0                           *(wavelengthObservedLymanContinuum)**2.1d0 &
                  &           -2.5450d-2                           *(wavelengthObservedLymanContinuum)**3.7d0
          else if (wavelengthObservedLymanContinuum< 1.0d0+redshift) then
             opticalDepth    =+opticalDepth                                                                   &
                  &           +5.2210d-4                                                                      &
                  &           *(                                                                              &
                  &             +        (1.0d0+redshift)**(+3.4d0)*(wavelengthObservedLymanContinuum)**2.1d0 &
                  &             -                                   (wavelengthObservedLymanContinuum)**5.5d0 &
                  &            )
          end if
       end if
       ! DLA continuum absorption.
       if (redshift < 2.0d0) then
          if      (wavelengthObservedLymanContinuum< 1.0d0+redshift)                                          &
               & opticalDepth=+opticalDepth                                                                   &
               &              +0.2113d+0*(1.0d0+redshift)**(+2.0d0)                                           &
               &              -7.6610d-2*(1.0d0+redshift)**(+2.3d0)/(wavelengthObservedLymanContinuum)**0.3d0 &
               &              -0.1347d+0                           *(wavelengthObservedLymanContinuum)**2.0d0
       else
          if      (wavelengthObservedLymanContinuum< 3.0d0         ) then
             opticalDepth    =+opticalDepth                                                                   &
                  &           +0.6340d+0                                                                      &
                  &           +4.6960d-2*(1.0d0+redshift)**(+3.0d0)                                           &
                  &           -1.7790d-2*(1.0d0+redshift)**(+3.3d0)/(wavelengthObservedLymanContinuum)**0.3d0 &
                  &           -0.1347d+0                           *(wavelengthObservedLymanContinuum)**2.0d0 &
                  &           -0.2905d+0                           /(wavelengthObservedLymanContinuum)**0.3d0
          else if (wavelengthObservedLymanContinuum< 1.0d0+redshift) then
             opticalDepth    =+opticalDepth                                                                   &
                  &           +4.6960d-2*(1.0d0+redshift)**(+3.0d0)                                           &
                  &           -1.7790d-2*(1.0d0+redshift)**(+3.3d0)/(wavelengthObservedLymanContinuum)**0.3d0 &
                  &           -2.9160d-2                           *(wavelengthObservedLymanContinuum)**3.0d0
          end if
       end if
    end if
    ! Compute attenuation from optical depth.
    inoue2014Multiplier=exp(-opticalDepth)
    return
  end function inoue2014Multiplier

