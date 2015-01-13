!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which computes cooling rates due to collisional excitation using the fitting functions of
!% \cite{scholz_collisional_1991}.

module Atomic_Rates_Excitation_Collisional_ScholzWalters91
  !% Contains a function which computes cooling rates due to collisional excitation using the fitting functions of
  !% \cite{scholz_collisional_1991}.
  private
  public :: Collisional_Excitation_Cooling_Rate_ScholzWalters91_Initialize
  
contains

  !# <collisionalExcitationMethod>
  !#  <unitName>Collisional_Excitation_Cooling_Rate_ScholzWalters91_Initialize</unitName>
  !# </collisionalExcitationMethod>
  subroutine Collisional_Excitation_Cooling_Rate_ScholzWalters91_Initialize(collisionalExcitationMethod,Collisional_Excitation_Rate_Get)
    !% Initializes the ``ScholzWalters91'' collisional excitation rate module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                     ), intent(in   )          :: collisionalExcitationMethod
    procedure(Collisional_Excitation_Cooling_Rate), intent(inout), pointer :: Collisional_Excitation_Rate_Get

    ! Check if this collisional excitation rate method has been selected.
    if (collisionalExcitationMethod == 'ScholzWalters91') Collisional_Excitation_Rate_Get => Collisional_Excitation_Cooling_Rate
    return
  end subroutine Collisional_Excitation_Cooling_Rate_ScholzWalters91_Initialize
  
  double precision function Collisional_Excitation_Cooling_Rate(atomicNumber, electronNumber, temperature)
    !% Return collisional excitation cooling rates, in units of J/m$^3$/s, for ion {\normalfont \ttfamily Ion} at temperature {\normalfont \ttfamily T} (in Kelvin)
    !% using the fitting functions of \cite{scholz_collisional_1991}.
    use Galacticus_Error
    implicit none
    double precision,  intent(in   ) :: temperature                     
    integer         ,  intent(in   ) :: atomicNumber          , electronNumber
    double precision                 :: temperatureLogarithmic

    ! Zero cooling rate by default in case we can't identify the ionic species.
    Collisional_Excitation_Cooling_Rate=0.0d0
    ! Return immediately on unphysical temperature.
    if (temperature <= 0.0d0) return
    ! Compute cooling rate for HI atoms.
    if (atomicNumber == 1 .and. electronNumber == 1) then
       ! Fitting formula from Scholz & Walters, 1991, ApJ, 380, 302.
       temperatureLogarithmic = log(temperature)
       if (temperature <= 2.0d3) then
          Collisional_Excitation_Cooling_Rate=0.0d0
       else if (temperature > 2.0d3 .and. temperature < 1.0d5) then
          Collisional_Excitation_Cooling_Rate=+1.0d-33                                     &
               &                              *exp(                                        &
               &                                   +2.1379130d+2                           &
               &                                   -1.1394920d+2*temperatureLogarithmic    &
               &                                   +2.5060620d+1*temperatureLogarithmic**2 &
               &                                   -2.7627550d+0*temperatureLogarithmic**3 &
               &                                   +1.5153520d-1*temperatureLogarithmic**4 &
               &                                   -3.2903820d-3*temperatureLogarithmic**5 &
               &                                   -118348.00d+0/temperature               &
               &                                  )
       else
          Collisional_Excitation_Cooling_Rate=+1.0d-33                                     &
               &                              *exp(                                        &
               &                                   +2.7125446d+2                           &
               &                                   -9.8019455d+1*temperatureLogarithmic    &
               &                                   +1.4007276d+1*temperatureLogarithmic**2 &
               &                                   -9.7808421d-1*temperatureLogarithmic**3 &
               &                                   +3.3562890d-2*temperatureLogarithmic**4 &
               &                                   -4.5533231d-4*temperatureLogarithmic**5 &
               &                                   -118348.00d+0/temperature               &
               &                              )           
       end if
       return
    end if
    ! Compute cooling rate for HeII ions.
    if (atomicNumber == 2 .and. electronNumber == 1) then
       if (temperature < 10.0d0) then
          Collisional_Excitation_Cooling_Rate=0.0d0
       else
          ! Fitting formula based on results from Mappings III.
          temperatureLogarithmic=log(temperature)
          if (temperature <= 2.0d3) then
             Collisional_Excitation_Cooling_Rate=0.0d0           
          else if (temperature > 2.0d3 .and. temperature < 1.9952623149687d6) then
             Collisional_Excitation_Cooling_Rate=+1.0d-31                                      &
                  &                              *exp(                                         &
                  &                                   +4.26241200d+0                           &
                  &                                   -6.16369360d-2*temperatureLogarithmic    &
                  &                                   -4.28615360d-2*temperatureLogarithmic**2 &
                  &                                   -3.47087772d-3*temperatureLogarithmic**3 &
                  &                                   +8.48727210d-5*temperatureLogarithmic**4 &
                  &                                   +2.29204450d-5*temperatureLogarithmic**5 &
                  &                                   +2.25993770d-6*temperatureLogarithmic**6 &
                  &                                   -1.73417580d-7*temperatureLogarithmic**7 &
                  &                                   -473638.000d+0/temperature               &
                  &                                  )
          else
             Collisional_Excitation_Cooling_Rate=+1.0d-31                                     &
                  &                              *exp(                                        &
                  &                                   +5.5213270d+0                           &
                  &                                   -0.3369025d+0*temperatureLogarithmic    &
                  &                                   +5.7712696d-3*temperatureLogarithmic**2 &
                  &                                   -1.1971748d-3*temperatureLogarithmic**3 &
                  &                                   +4.4052693d-5*temperatureLogarithmic**4 &
                  &                                   -4.4549989d-7*temperatureLogarithmic**5 &
                  &                                   -473638.00d+0/temperature               &
                  &                                  )
          end if
       end if
       return
    end if
    return
  end function Collisional_Excitation_Cooling_Rate

end module Atomic_Rates_Excitation_Collisional_ScholzWalters91
