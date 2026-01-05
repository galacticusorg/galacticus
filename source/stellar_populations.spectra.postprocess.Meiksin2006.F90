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
  An implementation of a spectrum postprocessor that applies the \cite{meiksin_colour_2006} calculation of the attenuation of spectra by the intergalactic medium.
  !!}

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorMeiksin2006">
   <description>
    A stellar population postprocessor class that postprocesses spectra through absorption by the \gls{igm} using the results
    of \cite{meiksin_colour_2006}.
   </description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorMeiksin2006
     !!{
     An spectrum postprocessor implementing the \cite{meiksin_colour_2006} calculation of the attenuation of spectra by the intergalactic medium.
     !!}
     private
   contains
     procedure :: multiplier => meiksin2006Multiplier
  end type stellarPopulationSpectraPostprocessorMeiksin2006

  interface stellarPopulationSpectraPostprocessorMeiksin2006
     !!{
     Constructors for the \refClass{stellarPopulationSpectraPostprocessorMeiksin2006} stellar population spectra postprocessor class.
     !!}
     module procedure meiksin2006ConstructorParameters
  end interface stellarPopulationSpectraPostprocessorMeiksin2006

contains

  function meiksin2006ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSpectraPostprocessorMeiksin2006} stellar population spectra postprocessor class which takes a
    parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(stellarPopulationSpectraPostprocessorMeiksin2006)                :: self
    type(inputParameters                                 ), intent(inout) :: parameters

    self=stellarPopulationSpectraPostprocessorMeiksin2006()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function meiksin2006ConstructorParameters

  double precision function meiksin2006Multiplier(self,wavelength,age,redshift)
    !!{
    Suppress the Lyman continuum in a spectrum.
    !!}
    use :: Factorials                , only : Factorial
    use :: Gamma_Functions           , only : Gamma_Function_Logarithmic
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    implicit none
    class           (stellarPopulationSpectraPostprocessorMeiksin2006), intent(inout) :: self
    double precision                                                  , intent(in   ) :: age                                    , redshift           , &
         &                                                                               wavelength
    ! Parameters of the Lyman-limit system distribution.
    double precision                                                  , parameter     :: N0                              =0.25d0
    double precision                                                  , parameter     :: beta                            =1.50d0
    double precision                                                  , parameter     :: gamma                           =1.50d0
    double precision                                                  , dimension(31) :: opticalDepthLymanLines                 , redshiftLymanLines
    integer                                                                           :: iLine
    double precision                                                                  :: nFactorial                             , opticalDepth       , &
         &                                                                               seriesSolutionTermA                    , seriesSolutionTermB, &
         &                                                                               wavelengthObservedLymanContinuum
    !$GLC attributes unused :: self, age

    ! Check if this is a zero redshift case.
    if (redshift <= 0.0d0) then
       ! It is, so return no attenuation modification.
       meiksin2006Multiplier=1.0d0
       return
    else
       ! Compute the observed wavelength in units of the Lyman-continuum wavelength.
       wavelengthObservedLymanContinuum=wavelength*(1.0d0+redshift)/lymanSeriesLimitWavelengthHydrogen_atomic
       ! Evaluate redshifts of various Lyman-series lines.
       forall (iLine=3:9)
          redshiftLymanLines(iLine)=wavelengthObservedLymanContinuum*(1.0d0-1.0d0/dble(iLine**2))-1.0d0
       end forall
       ! Evaluate optical depths relative to Lyman-α.
       opticalDepthLymanLines(2)=1.0d0 ! By definition.
       if (redshiftLymanLines(3) < 3.0d0) then
          opticalDepthLymanLines(3)=0.348d0*(0.25d0*(1.0+redshiftLymanLines(3)))**0.3333d0
       else
          opticalDepthLymanLines(3)=0.348d0*(0.25d0*(1.0+redshiftLymanLines(3)))**0.1667d0
       end if
       if (redshiftLymanLines(4) < 3.0d0) then
          opticalDepthLymanLines(4)=0.179d0*(0.25d0*(1.0d0+redshiftLymanLines(4)))**0.3333d0
       else
          opticalDepthLymanLines(4)=0.179d0*(0.25d0*(1.0d0+redshiftLymanLines(4)))**0.1667d0
       end if
       if (redshiftLymanLines(5) < 3.0d0) then
          opticalDepthLymanLines(5)=0.109d0*(0.25d0*(1.0d0+redshiftLymanLines(5)))**0.3333d0
       else
          opticalDepthLymanLines(5)=0.109d0*(0.25d0*(1.0d0+redshiftLymanLines(5)))**0.1667d0
       end if
       opticalDepthLymanLines(6)=0.0722d0*(0.25d0*(1.0d0+redshiftLymanLines(6)))**0.3333d0
       opticalDepthLymanLines(7)=0.0508d0*(0.25d0*(1.0d0+redshiftLymanLines(7)))**0.3333d0
       opticalDepthLymanLines(8)=0.0373d0*(0.25d0*(1.0d0+redshiftLymanLines(8)))**0.3333d0
       opticalDepthLymanLines(9)=0.0283d0*(0.25d0*(1.0d0+redshiftLymanLines(9)))**0.3333d0
       forall (iLine=10:31)
          opticalDepthLymanLines(iLine)=opticalDepthLymanLines(9)*720.0d0/dble(iLine)/dble(iLine**2-1)
       end forall
       ! Scale optical depths by Lyman-α optical depth.
       if (redshift <= 4.0d0) then
          forall (iLine=2:31)
             opticalDepthLymanLines(iLine)=opticalDepthLymanLines(iLine)*0.00211d0*(wavelengthObservedLymanContinuum*(1.0d0-1.0d0/dble(iLine**2)))**3.70d0
          end forall
       else
          forall (iLine=2:31)
             opticalDepthLymanLines(iLine)=opticalDepthLymanLines(iLine)*0.00058d0*(wavelengthObservedLymanContinuum*(1.0d0-1.0d0/dble(iLine**2)))**4.50d0
          end forall
       end if
       ! Accumulate optical depths if line falls within the required redshift range.
       opticalDepth=0.0d0
       do iLine=2,31
          if (wavelengthObservedLymanContinuum < (1.0d0+redshift)/(1.0d0-1.0d0/(dble(iLine)**2))) opticalDepth=opticalDepth+opticalDepthLymanLines(iLine)
       end do
       if (wavelengthObservedLymanContinuum < (1.0d0+redshift)) then
          ! Add in photoelectric absorption contributions.
          seriesSolutionTermA=0.0d0
          seriesSolutionTermB=0.0d0
          do iLine=0,9
             nFactorial=Factorial(iLine)
             seriesSolutionTermA=seriesSolutionTermA+dble(-1**iLine)*(beta-1.0d0)/(dble(iLine)+1.0d0-beta)/nFactorial
             seriesSolutionTermB=seriesSolutionTermB+dble(-1**iLine)*(beta-1.0d0)*(((1.0d0+redshift)**(gamma+1.0d0-3.0d0 &
                  &*dble(iLine))*(wavelengthObservedLymanContinuum**(3.0d0*dble(iLine))) -(wavelengthObservedLymanContinuum**(gamma &
                  &+1.0d0))))/(dble(iLine)+1.0d0-beta)/(3.0d0*dble(iLine)-gamma-1.0d0)/nFactorial
          end do
          ! Add contribution due to Lyman-limit systems.
          opticalDepth=opticalDepth+N0*(exp(Gamma_Function_Logarithmic(2.0d0-beta))-exp(-1.0d0)-seriesSolutionTermA)*(((1.0d0+redshift)**(&
               &-3.0d0*(beta-1.0d0)+gamma+1.0d0))*(wavelengthObservedLymanContinuum**(3.0d0*(beta-1.0d0)))&
               &-(wavelengthObservedLymanContinuum**(gamma+1.0d0)))/(4.0d0+gamma-3.0d0*beta)-N0*seriesSolutionTermB
          ! Add contribution due to optically thin systems.
          opticalDepth=opticalDepth+0.805d0*(wavelengthObservedLymanContinuum**3)*(1.0d0/wavelengthObservedLymanContinuum-1.0d0 &
               &/(1.0d0+redshift))
       end if
       ! Compute attenuation from optical depth.
       meiksin2006Multiplier=exp(-opticalDepth)
    end if
    return
  end function meiksin2006Multiplier

