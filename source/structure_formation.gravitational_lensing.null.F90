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
Implements a null gravitational lensing distribution.
!!}

  !![
  <gravitationalLensing name="gravitationalLensingNull">
   <description>Implements a null gravitational lensing distribution.</description>
  </gravitationalLensing>
  !!]
  type, extends(gravitationalLensingClass) :: gravitationalLensingNull
     private
   contains
     procedure :: magnificationPDF => nullMagnificationPDF
     procedure :: magnificationCDF => nullMagnificationCDF
   end type gravitationalLensingNull

  interface gravitationalLensingNull
     !!{
     Constructors for the null gravitational lensing class.
     !!}
     module procedure nullConstructorParameters
  end interface gravitationalLensingNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{takahashi_probability_2011} gravitational lensing class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(gravitationalLensingNull)                :: self
    type(inputParameters         ), intent(inout) :: parameters
    
    self=gravitationalLensingNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  double precision function nullMagnificationPDF(self,magnification,redshift,scaleSource)
    !!{
    Compute the magnification probability density function for a null lensing case.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (gravitationalLensingNull), intent(inout) :: self
    double precision                          , intent(in   ) :: magnification, redshift, &
         &                                                       scaleSource

    nullMagnificationPDF=0.0d0
    call Error_Report('the PDF is a Dirac delta function, so is not implemented'//{introspection:location})
    return
  end function nullMagnificationPDF

  double precision function nullMagnificationCDF(self,magnification,redshift,scaleSource)
    !!{
    Compute the magnification probability density function at the given {\normalfont \ttfamily magnification} and {\normalfont
    \ttfamily redshift} for a null lensing case.
    !!}
    implicit none
    class           (gravitationalLensingNull), intent(inout) :: self
    double precision                          , intent(in   ) :: magnification, redshift, &
         &                                                       scaleSource

    if (magnification < 1.0d0) then
       nullMagnificationCDF=0.0d0
    else
       nullMagnificationCDF=1.0d0
    end if
    return
  end function nullMagnificationCDF
