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

!+ Contributions to this file made by: Daniel McAndrew.

!% Contains a module which computes rates of dielectronic recombination.

module Atomic_Rates_Recombination_Dielectronic_Arnaud1985
  !% Contains a function which computes rates of dielectronic recombination.
  private
  public :: Dielectronic_Recom_Arnaud1985_Initialize

  ! Array to hold the coefficients of the fitting function.
  double precision, dimension(4,28,28) :: coefficients
  logical                              :: coefficientsSet=.false.

contains

  !# <dielectronicRecombinationMethod>
  !#  <unitName>Dielectronic_Recom_Arnaud1985_Initialize</unitName>
  !# </dielectronicRecombinationMethod>
  subroutine Dielectronic_Recom_Arnaud1985_Initialize(dielectronicRecombinationMethod,Dielectronic_Recombination_Rate_Get)
    !% Initializes the ``Arnaud1985'' dielectronic recombination rate module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                            ), intent(in   )          :: dielectronicRecombinationMethod
    procedure(Dielectronic_Recombination_Rate_Arnaud1985), intent(inout), pointer :: Dielectronic_Recombination_Rate_Get

    ! Check if this dielectronic recombination rate method has been selected.
    if (dielectronicRecombinationMethod == 'Arnaud85') Dielectronic_Recombination_Rate_Get => Dielectronic_Recombination_Rate_Arnaud1985
    return
  end subroutine Dielectronic_Recom_Arnaud1985_Initialize

  double precision function Dielectronic_Recombination_Rate_Arnaud1985(atomicNumber,electronNumber,temperature)
    !% This function calculates rates of dielectric recombination for all ionization stages of all elements from H to Ni ($Z=28$)
    !% by use of the fits from \cite{aldrovandi_radiative_1973}, \cite{shull_ionization_1982} and \cite{arnaud_updated_1985}.
    !% Input parameters: {\tt atomicNumber}: atomic number; {\tt electronNumber}: number of electrons; {\tt temperature}: temperature [K].
    !% Output parameter: rate coefficient [cm$^3$ s$^{-1}$].
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: temperature                     
    integer         , intent(in   ) :: atomicNumber, electronNumber

    ! Initialize data arrays.
    if (.not.coefficientsSet) call Set_Coefficients()
    ! Set zero rate by default.
    Dielectronic_Recombination_Rate_Arnaud1985=0.0d0
    ! Return on unphysical temperature.
    if (temperature <= 0.0) return
    ! Abort on unphysical conditions.
    if (atomicNumber   < 1 .or. atomicNumber   >           28) call Galacticus_Error_Report('Dielectronic_Recombination_Rate_Arnaud1985()','Fatal: atomic number is unphysical or too large')
    if (electronNumber < 1 .or. electronNumber > atomicNumber) call Galacticus_Error_Report('Dielectronic_Recombination_Rate_Arnaud1985()','Fatal: electron number is unphysical or too large')
    ! Evaluate the fitting function.
    Dielectronic_Recombination_Rate_Arnaud1985=                             &
         & +       coefficients(1,atomicNumber,electronNumber)              &
         & *exp  (-coefficients(3,atomicNumber,electronNumber)/temperature) &
         & *(                                                               &
         &   +1.0d0                                                         &
         &   +     coefficients(2,atomicNumber,electronNumber)              &
         &   *exp(-coefficients(4,atomicNumber,electronNumber)/temperature) &
         & )                                                                &
         & /temperature**1.5d0
    return
  end function Dielectronic_Recombination_Rate_Arnaud1985

  subroutine  Set_coefficients
    !% Populate the array of fitting function coefficients for the {\tt Arnaud1985} dielectronic recombination rates.
    implicit none

    coefficientsSet=.true.
    coefficients(:, 2, 2)=[1.90d-03, 3.00d-01, 4.70d+05, 9.40d+04]
    coefficients(:, 6, 6)=[2.54d-03, 4.42d-02, 1.57d+05, 3.74d+05]
    coefficients(:, 6, 5)=[6.15d-03, 5.88d-02, 1.41d+05, 1.41d+05]
    coefficients(:, 6, 4)=[1.62d-03, 3.43d-01, 8.19d+04, 1.59d+05]
    coefficients(:, 6, 3)=[4.78d-02, 3.62d-01, 3.44d+06, 5.87d+05]
    coefficients(:, 6, 2)=[3.22d-02, 3.15d-01, 4.06d+06, 8.31d+05]
    coefficients(:, 7, 7)=[2.98d-03, 0.00d+00, 2.20d+05, 1.00d+05]
    coefficients(:, 7, 6)=[7.41d-03, 7.64d-02, 2.01d+05, 7.37d+04]
    coefficients(:, 7, 5)=[1.13d-02, 1.64d-01, 1.72d+05, 2.25d+05]
    coefficients(:, 7, 4)=[2.62d-03, 2.43d-01, 1.02d+05, 1.25d+05]
    coefficients(:, 7, 3)=[7.50d-02, 3.50d-01, 4.75d+06, 8.35d+05]
    coefficients(:, 7, 2)=[4.61d-02, 3.09d-01, 5.44d+06, 1.14d+06]
    coefficients(:, 8, 8)=[1.11d-03, 9.25d-02, 1.75d+05, 1.45d+05]
    coefficients(:, 8, 7)=[5.07d-03, 1.81d-01, 1.98d+05, 3.35d+05]
    coefficients(:, 8, 6)=[1.48d-02, 3.05d-01, 2.41d+05, 2.83d+05]
    coefficients(:, 8, 5)=[1.84d-02, 1.00d-01, 2.12d+05, 2.83d+05]
    coefficients(:, 8, 4)=[4.13d-03, 1.62d-01, 1.25d+05, 2.27d+05]
    coefficients(:, 8, 3)=[1.06d-01, 3.40d-01, 6.25d+06, 1.12d+06]
    coefficients(:, 8, 2)=[4.72d-02, 0.00d+00, 6.00d+06, 1.00d+05]
    coefficients(:,10,10)=[9.77d-04, 7.30d-02, 3.11d+05, 2.06d+05]
    coefficients(:,10, 9)=[2.65d-03, 2.42d-01, 2.84d+05, 3.07d+05]
    coefficients(:,10, 8)=[3.69d-03, 1.01d+00, 2.24d+05, 2.94d+05]
    coefficients(:,10, 7)=[1.18d-02, 3.91d-01, 2.70d+05, 5.50d+05]
    coefficients(:,10, 6)=[2.44d-02, 2.52d+00, 3.09d+05, 9.91d+05]
    coefficients(:,10, 5)=[3.02d-02, 4.45d-01, 2.83d+05, 1.73d+06]
    coefficients(:,10, 4)=[6.10d-03, 2.54d-01, 1.68d+05, 6.13d+05]
    coefficients(:,10, 3)=[2.52d-01, 3.04d-01, 1.40d+07, 1.80d+06]
    coefficients(:,10, 2)=[7.14d-02, 2.96d-01, 1.10d+07, 2.24d+06]
    coefficients(:,12,12)=[4.49d-04, 2.10d-02, 5.01d+04, 2.81d+04]
    coefficients(:,12,11)=[1.95d-03, 7.40d-02, 6.06d+05, 1.44d+06]
    coefficients(:,12,10)=[5.12d-03, 3.23d-01, 4.69d+05, 7.55d+05]
    coefficients(:,12, 9)=[7.74d-03, 6.36d-01, 3.74d+05, 7.88d+05]
    coefficients(:,12, 8)=[1.17d-02, 8.07d-01, 3.28d+05, 1.02d+06]
    coefficients(:,12, 7)=[3.69d-02, 3.51d-01, 4.80d+05, 9.73d+05]
    coefficients(:,12, 6)=[3.63d-02, 5.48d-01, 3.88d+05, 7.36d+05]
    coefficients(:,12, 5)=[4.15d-02, 2.33d-01, 3.39d+05, 3.82d+05]
    coefficients(:,12, 4)=[8.86d-03, 3.18d-01, 2.11d+05, 1.54d+06]
    coefficients(:,12, 3)=[2.52d-01, 3.15d-01, 1.40d+07, 2.64d+06]
    coefficients(:,12, 2)=[9.28d-02, 0.00d+00, 1.45d+07, 1.00d+05]
    coefficients(:,13,13)=[2.00d-03, 1.34d+00, 4.93d+04, 1.01d+05]
    coefficients(:,13,12)=[2.35d-03, 1.25d-01, 6.72d+04, 2.45d+04]
    coefficients(:,13,11)=[2.97d-03, 3.49d-01, 7.90d+05, 9.20d+05]
    coefficients(:,13,10)=[6.97d-03, 1.12d-01, 8.42d+05, 3.34d+05]
    coefficients(:,13, 9)=[1.22d-02, 1.29d+00, 4.29d+05, 8.42d+05]
    coefficients(:,13, 8)=[1.84d-02, 1.36d+00, 3.72d+05, 1.45d+06]
    coefficients(:,13, 7)=[3.32d-02, 1.24d+00, 4.17d+05, 1.55d+06]
    coefficients(:,13, 6)=[3.78d-02, 8.02d-01, 3.66d+05, 1.41d+06]
    coefficients(:,13, 5)=[5.12d-02, 1.72d-01, 3.59d+05, 3.39d+05]
    coefficients(:,13, 4)=[1.15d-02, 9.43d-01, 2.31d+05, 2.97d+06]
    coefficients(:,13, 3)=[2.85d-01, 3.12d-01, 1.61d+07, 3.11d+06]
    coefficients(:,13, 2)=[1.02d-01, 1.25d-01, 1.71d+07, 1.91d+06]
    coefficients(:,14,14)=[1.10d-03, 0.00d+00, 7.70d+04, 1.00d+05]
    coefficients(:,14,13)=[5.87d-03, 7.53d-01, 9.63d+04, 6.46d+04]
    coefficients(:,14,12)=[5.03d-03, 1.88d-01, 8.75d+04, 4.71d+04]
    coefficients(:,14,11)=[5.43d-03, 4.50d-01, 1.05d+06, 7.98d+05]
    coefficients(:,14,10)=[8.86d-03, 0.00d+00, 1.14d+06, 1.00d+05]
    coefficients(:,14, 9)=[1.68d-02, 1.80d+00, 4.85d+05, 1.03d+06]
    coefficients(:,14, 8)=[2.49d-02, 1.88d+00, 4.15d+05, 1.91d+06]
    coefficients(:,14, 7)=[3.13d-02, 2.01d+00, 3.66d+05, 2.11d+06]
    coefficients(:,14, 6)=[4.25d-02, 1.22d+00, 3.63d+05, 2.14d+06]
    coefficients(:,14, 5)=[6.18d-02, 3.03d-01, 3.88d+05, 1.12d+06]
    coefficients(:,14, 4)=[1.38d-02, 1.42d+00, 2.51d+05, 3.93d+06]
    coefficients(:,14, 3)=[3.27d-01, 3.06d-01, 1.88d+07, 3.60d+06]
    coefficients(:,14, 2)=[1.13d-01, 2.86d-01, 1.99d+07, 4.14d+06]
    coefficients(:,16,16)=[1.62d-03, 0.00d+00, 1.25d+05, 1.00d+05]
    coefficients(:,16,15)=[1.09d-02, 1.20d-02, 1.92d+05, 1.80d+04]
    coefficients(:,16,14)=[3.35d-02, 6.59d-02, 1.89d+05, 1.59d+05]
    coefficients(:,16,13)=[3.14d-02, 6.89d-02, 1.68d+05, 8.04d+04]
    coefficients(:,16,12)=[1.27d-02, 1.87d-01, 1.38d+05, 1.71d+05]
    coefficients(:,16,11)=[1.47d-02, 1.29d-01, 1.80d+06, 1.75d+06]
    coefficients(:,16,10)=[1.34d-02, 1.04d+00, 6.90d+05, 2.15d+06]
    coefficients(:,16, 9)=[2.38d-02, 1.12d+00, 5.84d+05, 2.59d+06]
    coefficients(:,16, 8)=[3.19d-02, 1.40d+00, 5.17d+05, 2.91d+06]
    coefficients(:,16, 7)=[7.13d-02, 1.00d+00, 6.66d+05, 2.32d+06]
    coefficients(:,16, 6)=[8.00d-02, 5.55d-01, 6.00d+05, 2.41d+06]
    coefficients(:,16, 5)=[7.96d-02, 1.63d+00, 5.09d+05, 6.37d+06]
    coefficients(:,16, 4)=[1.34d-02, 3.04d-01, 2.91d+05, 1.04d+06]
    coefficients(:,16, 3)=[4.02d-01, 2.98d-01, 2.41d+07, 4.67d+06]
    coefficients(:,16, 2)=[1.45d-01, 2.81d-01, 2.54d+07, 5.30d+06]
    coefficients(:,18,18)=[1.00d-03, 5.00d-03, 3.20d+05, 3.10d+05]
    coefficients(:,18,17)=[1.10d-02, 4.50d-02, 2.90d+05, 5.50d+05]
    coefficients(:,18,16)=[3.40d-02, 5.70d-02, 2.39d+05, 6.00d+05]
    coefficients(:,18,15)=[6.85d-02, 8.70d-02, 2.56d+05, 3.81d+05]
    coefficients(:,18,14)=[9.00d-02, 7.69d-02, 2.50d+05, 3.30d+05]
    coefficients(:,18,13)=[6.35d-02, 1.40d-01, 2.10d+05, 2.15d+05]
    coefficients(:,18,12)=[2.60d-02, 1.20d-01, 1.80d+05, 2.15d+05]
    coefficients(:,18,11)=[1.70d-02, 1.00d-01, 2.70d+06, 3.30d+06]
    coefficients(:,18,10)=[2.10d-02, 1.92d+00, 8.30d+05, 3.50d+06]
    coefficients(:,18, 9)=[3.50d-02, 1.66d+00, 6.95d+05, 3.60d+06]
    coefficients(:,18, 8)=[4.30d-02, 1.67d+00, 6.05d+05, 3.80d+06]
    coefficients(:,18, 7)=[7.13d-02, 1.40d+00, 6.68d+05, 2.90d+06]
    coefficients(:,18, 6)=[9.60d-02, 1.31d+00, 6.50d+05, 3.60d+06]
    coefficients(:,18, 5)=[8.50d-02, 1.02d+00, 5.30d+05, 2.80d+06]
    coefficients(:,18, 4)=[1.70d-02, 2.45d-01, 3.55d+05, 1.10d+06]
    coefficients(:,18, 3)=[4.76d-01, 2.94d-01, 3.01d+07, 6.05d+06]
    coefficients(:,18, 2)=[2.97d-01, 2.77d-01, 3.13d+07, 6.54d+06]
    coefficients(:,20,20)=[3.28d-04, 9.07d-02, 3.46d+04, 1.64d+04]
    coefficients(:,20,19)=[5.84d-02, 1.10d-01, 3.85d+05, 2.45d+05]
    coefficients(:,20,18)=[1.12d-01, 1.74d-02, 4.08d+05, 4.27d+05]
    coefficients(:,20,17)=[1.32d-01, 1.32d-01, 3.82d+05, 6.92d+05]
    coefficients(:,20,16)=[1.33d-01, 1.14d-01, 3.53d+05, 8.78d+05]
    coefficients(:,20,15)=[1.26d-01, 1.62d-01, 3.19d+05, 7.43d+05]
    coefficients(:,20,14)=[1.39d-01, 8.78d-02, 3.22d+05, 6.99d+05]
    coefficients(:,20,13)=[9.55d-02, 2.63d-01, 2.47d+05, 4.43d+05]
    coefficients(:,20,12)=[4.02d-02, 6.27d-02, 2.29d+05, 2.81d+05]
    coefficients(:,20,11)=[4.19d-02, 6.16d-02, 3.73d+06, 5.84d+06]
    coefficients(:,20,10)=[2.57d-02, 2.77d+00, 9.26d+05, 4.89d+06]
    coefficients(:,20, 9)=[4.45d-02, 2.23d+00, 7.96d+05, 4.62d+06]
    coefficients(:,20, 8)=[5.48d-02, 2.00d+00, 6.90d+05, 4.52d+06]
    coefficients(:,20, 7)=[7.13d-02, 1.82d+00, 6.70d+05, 3.32d+06]
    coefficients(:,20, 6)=[1.09d-01, 1.74d+00, 7.00d+05, 4.93d+06]
    coefficients(:,20, 5)=[1.10d-01, 2.43d-01, 5.67d+05, 4.41d+06]
    coefficients(:,20, 4)=[2.05d-02, 1.85d-01, 4.21d+05, 2.27d+06]
    coefficients(:,20, 3)=[5.49d-01, 2.92d-01, 3.65d+07, 7.25d+06]
    coefficients(:,20, 2)=[2.68d-01, 0.00d+00, 3.74d+07, 1.00d+05]
    coefficients(:,26,26)=[1.58d-03, 4.56d-01, 6.00d+04, 8.97d+04]
    coefficients(:,26,25)=[8.38d-03, 3.23d-01, 1.94d+05, 1.71d+05]
    coefficients(:,26,24)=[1.54d-02, 3.10d-01, 3.31d+05, 2.73d+05]
    coefficients(:,26,23)=[3.75d-02, 4.11d-01, 4.32d+05, 3.49d+05]
    coefficients(:,26,22)=[1.17d-01, 3.59d-01, 6.28d+05, 5.29d+05]
    coefficients(:,26,21)=[2.54d-01, 9.75d-02, 7.50d+05, 4.69d+05]
    coefficients(:,26,20)=[2.91d-01, 2.29d-01, 7.73d+05, 6.54d+05]
    coefficients(:,26,19)=[1.50d-01, 4.20d+00, 2.62d+05, 1.32d+06]
    coefficients(:,26,18)=[1.40d-01, 3.30d+00, 2.50d+05, 1.33d+06]
    coefficients(:,26,17)=[1.00d-01, 5.30d+00, 2.57d+05, 1.41d+06]
    coefficients(:,26,16)=[2.00d-01, 1.50d+00, 2.84d+05, 1.52d+06]
    coefficients(:,26,15)=[2.40d-01, 7.00d-01, 8.69d+05, 1.51d+06]
    coefficients(:,26,14)=[2.60d-01, 6.00d-01, 4.21d+05, 1.82d+06]
    coefficients(:,26,13)=[1.90d-01, 5.00d-01, 4.57d+05, 1.84d+06]
    coefficients(:,26,12)=[1.20d-01, 1.00d+00, 2.85d+05, 2.31d+06]
    coefficients(:,26,11)=[3.50d-01, 0.00d+00, 8.18d+06, 1.00d+05]
    coefficients(:,26,10)=[6.60d-02, 7.80d+00, 1.51d+06, 9.98d+06]
    coefficients(:,26, 9)=[1.00d-01, 6.30d+00, 1.30d+06, 9.98d+06]
    coefficients(:,26, 8)=[1.30d-01, 5.50d+00, 1.19d+06, 1.00d+07]
    coefficients(:,26, 7)=[2.30d-01, 3.60d+00, 1.09d+06, 1.10d+07]
    coefficients(:,26, 6)=[1.40d-01, 4.90d+00, 9.62d+05, 8.34d+06]
    coefficients(:,26, 5)=[1.10d-01, 1.60d+00, 7.23d+05, 1.01d+07]
    coefficients(:,26, 4)=[4.10d-02, 4.20d+00, 4.23d+05, 1.07d+07]
    coefficients(:,26, 3)=[7.47d-01, 2.84d-01, 5.84d+07, 1.17d+07]
    coefficients(:,26, 2)=[3.69d-01, 0.00d+00, 6.00d+07, 1.00d+05]
    coefficients(:,28,28)=[1.41d-03, 4.69d-01, 9.82d+04, 1.01d+05]
    coefficients(:,28,27)=[5.20d-03, 3.57d-01, 2.01d+05, 1.91d+05]
    coefficients(:,28,26)=[1.38d-02, 2.81d-01, 3.05d+05, 2.32d+05]
    coefficients(:,28,25)=[2.30d-02, 1.28d-01, 4.20d+05, 3.18d+05]
    coefficients(:,28,24)=[4.19d-02, 4.17d-02, 5.56d+05, 4.55d+05]
    coefficients(:,28,23)=[6.83d-02, 5.58d-02, 6.72d+05, 5.51d+05]
    coefficients(:,28,22)=[1.22d-01, 3.46d-02, 7.93d+05, 5.28d+05]
    coefficients(:,28,21)=[3.00d-01, 0.00d+00, 9.00d+05, 1.00d+05]
    coefficients(:,28,20)=[1.50d-01, 1.90d+00, 1.00d+06, 5.50d+05]
    coefficients(:,28,19)=[6.97d-01, 2.77d-01, 7.81d+05, 8.87d+05]
    coefficients(:,28,18)=[7.09d-01, 1.35d-01, 7.64d+05, 1.80d+06]
    coefficients(:,28,17)=[6.44d-01, 1.34d-01, 7.44d+05, 1.25d+06]
    coefficients(:,28,16)=[5.25d-01, 1.92d-01, 6.65d+05, 1.89d+06]
    coefficients(:,28,15)=[4.46d-01, 3.32d-01, 5.97d+05, 8.84d+05]
    coefficients(:,28,14)=[3.63d-01, 3.37d-01, 5.24d+05, 1.29d+06]
    coefficients(:,28,13)=[3.02d-01, 1.21d-01, 4.96d+05, 6.24d+05]
    coefficients(:,28,12)=[1.02d-01, 5.14d-02, 4.46d+05, 1.59d+06]
    coefficients(:,28,11)=[2.70d-01, 1.83d-01, 8.49d+06, 8.01d+06]
    coefficients(:,28,10)=[4.67d-02, 7.56d+00, 1.36d+06, 9.32d+06]
    coefficients(:,28, 9)=[8.35d-02, 4.55d+00, 1.23d+06, 9.45d+06]
    coefficients(:,28, 8)=[9.96d-02, 4.87d+00, 1.06d+06, 9.45d+06]
    coefficients(:,28, 7)=[1.99d-01, 2.19d+00, 1.25d+06, 8.01d+06]
    coefficients(:,28, 6)=[2.40d-01, 1.15d+00, 1.23d+06, 7.57d+06]
    coefficients(:,28, 5)=[1.15d-01, 1.23d+00, 3.32d+05, 2.64d+06]
    coefficients(:,28, 4)=[3.16d-02, 1.32d-01, 6.45d+05, 1.93d+06]
    coefficients(:,28, 3)=[8.03d-01, 2.89d-01, 6.65d+07, 1.19d+07]
    coefficients(:,28, 2)=[5.75d-01, 2.86d-01, 6.81d+07, 9.08d+06]
  end subroutine Set_Coefficients

end module Atomic_Rates_Recombination_Dielectronic_Arnaud1985
