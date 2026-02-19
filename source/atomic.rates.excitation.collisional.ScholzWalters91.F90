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

!+ Contributions to this file made by: Daniel McAndrew.

  !!{
  An implementation of atomic collisional excitation using the fitting functions of \cite{scholz_collisional_1991}.
  !!}

  !![
  <atomicExcitationRateCollisional name="atomicExcitationRateCollisionalScholzWalters1991">
   <description>Atomic collisional excitation using the fitting functions of \cite{scholz_collisional_1991}.</description>
  </atomicExcitationRateCollisional>
  !!]
  type, extends(atomicExcitationRateCollisionalClass) :: atomicExcitationRateCollisionalScholzWalters1991
     !!{
     An atomic collisional excitation class using the fitting functions of \cite{scholz_collisional_1991}.
     !!}
     private
   contains
     procedure :: coolingRate => scholzWalters1991CoolingRate
  end type atomicExcitationRateCollisionalScholzWalters1991

  interface atomicExcitationRateCollisionalScholzWalters1991
     !!{
     Constructors for the \refClass{atomicExcitationRateCollisionalScholzWalters1991} atomic collisional excitation class.
     !!}
     module procedure scholzWalters1991ConstructorParameters
  end interface atomicExcitationRateCollisionalScholzWalters1991

contains

  function scholzWalters1991ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{atomicExcitationRateCollisionalScholzWalters1991} atomic collisional excitation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(atomicExcitationRateCollisionalScholzWalters1991)                :: self
    type(inputParameters                                 ), intent(inout) :: parameters

    self=atomicExcitationRateCollisionalScholzWalters1991()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function scholzWalters1991ConstructorParameters

  double precision function scholzWalters1991CoolingRate(self,atomicNumber,electronNumber,temperature)
    !!{
    Return collisional excitation cooling rates, in units of J m$^3$ s$^{-1}$, for ion {\normalfont \ttfamily Ion} at
    temperature {\normalfont \ttfamily T} (in Kelvin) using the fitting functions of \cite{scholz_collisional_1991}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (atomicExcitationRateCollisionalScholzWalters1991), intent(inout) :: self
    double precision                                                  , intent(in   ) :: temperature
    integer                                                           , intent(in   ) :: atomicNumber          , electronNumber
    double precision                                                                  :: temperatureLogarithmic
    !$GLC attributes unused :: self

    ! Zero cooling rate by default in case we can't identify the ionic species.
    scholzWalters1991CoolingRate=0.0d0
    ! Return immediately on unphysical temperature.
    if (temperature <= 0.0d0) return
    ! Compute cooling rate for HI atoms.
    if (atomicNumber == 1 .and. electronNumber == 1) then
       ! Fitting formula from Scholz & Walters, 1991, ApJ, 380, 302.
       temperatureLogarithmic = log(temperature)
       if (temperature <= 2.0d3) then
          scholzWalters1991CoolingRate=0.0d0
       else if (temperature > 2.0d3 .and. temperature < 1.0d5) then
          scholzWalters1991CoolingRate=+1.0d-33                                     &
               &                       *exp(                                        &
               &                            +2.1379130d+2                           &
               &                            -1.1394920d+2*temperatureLogarithmic    &
               &                            +2.5060620d+1*temperatureLogarithmic**2 &
               &                            -2.7627550d+0*temperatureLogarithmic**3 &
               &                            +1.5153520d-1*temperatureLogarithmic**4 &
               &                            -3.2903820d-3*temperatureLogarithmic**5 &
               &                            -118348.00d+0/temperature               &
               &                           )
       else
          scholzWalters1991CoolingRate=+1.0d-33                                     &
               &                       *exp(                                        &
               &                            +2.7125446d+2                           &
               &                            -9.8019455d+1*temperatureLogarithmic    &
               &                            +1.4007276d+1*temperatureLogarithmic**2 &
               &                            -9.7808421d-1*temperatureLogarithmic**3 &
               &                            +3.3562890d-2*temperatureLogarithmic**4 &
               &                            -4.5533231d-4*temperatureLogarithmic**5 &
               &                            -118348.00d+0/temperature               &
               &                           )
       end if
    end if
    ! Compute cooling rate for HeII ions.
    if (atomicNumber == 2 .and. electronNumber == 1) then
       if (temperature < 10.0d0) then
          scholzWalters1991CoolingRate=0.0d0
       else
          ! Fitting formula based on results from Mappings III.
          temperatureLogarithmic=log(temperature)
          if (temperature <= 2.0d3) then
             scholzWalters1991CoolingRate=0.0d0
          else if (temperature > 2.0d3 .and. temperature < 1.9952623149687d6) then
             scholzWalters1991CoolingRate=+1.0d-31                                      &
                  &                       *exp(                                         &
                  &                            +4.26241200d+0                           &
                  &                            -6.16369360d-2*temperatureLogarithmic    &
                  &                            -4.28615360d-2*temperatureLogarithmic**2 &
                  &                            -3.47087772d-3*temperatureLogarithmic**3 &
                  &                            +8.48727210d-5*temperatureLogarithmic**4 &
                  &                            +2.29204450d-5*temperatureLogarithmic**5 &
                  &                            +2.25993770d-6*temperatureLogarithmic**6 &
                  &                            -1.73417580d-7*temperatureLogarithmic**7 &
                  &                            -473638.000d+0/temperature               &
                  &                           )
          else
             scholzWalters1991CoolingRate=+1.0d-31                                     &
                  &                       *exp(                                        &
                  &                            +5.5213270d+0                           &
                  &                            -0.3369025d+0*temperatureLogarithmic    &
                  &                            +5.7712696d-3*temperatureLogarithmic**2 &
                  &                            -1.1971748d-3*temperatureLogarithmic**3 &
                  &                            +4.4052693d-5*temperatureLogarithmic**4 &
                  &                            -4.4549989d-7*temperatureLogarithmic**5 &
                  &                            -473638.00d+0/temperature               &
                  &                           )
          end if
       end if
    end if
    return
  end function scholzWalters1991CoolingRate
