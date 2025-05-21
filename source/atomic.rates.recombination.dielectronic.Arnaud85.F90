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

  !+ Contributions to this file made by: Andrew Benson, Daniel McAndrew.

  !!{
  Implements an atomic dielectronic recombination class which uses the fits from \cite{aldrovandi_radiative_1973},
  \cite{shull_ionization_1982} and \cite{arnaud_updated_1985}.
  !!}

  !![
  <atomicRecombinationRateDielectronic name="atomicRecombinationRateDielectronicArnaud1985">
   <description>
    Implements an atomic dielectronic recombination class which uses the fits from \cite{aldrovandi_radiative_1973}, \cite{shull_ionization_1982} and \cite{arnaud_updated_1985}.
   </description>
  </atomicRecombinationRateDielectronic>
  !!]
  type, extends(atomicRecombinationRateDielectronicClass) :: atomicRecombinationRateDielectronicArnaud1985
     !!{
     Implements an atomic dielectronic recombination class which uses the fits from \cite{aldrovandi_radiative_1973},
     \cite{shull_ionization_1982} and \cite{arnaud_updated_1985}.
     !!}
     private
   contains
     procedure :: rate => arnaud1985Rate
  end type atomicRecombinationRateDielectronicArnaud1985

  interface atomicRecombinationRateDielectronicArnaud1985
     !!{
     Constructors for the \refClass{atomicRecombinationRateDielectronicArnaud1985} atomic dielectronic recombination rate class.
     !!}
     module procedure arnaud1985ConstructorParameters
  end interface atomicRecombinationRateDielectronicArnaud1985

  ! Set coefficients for fitting functions.
  double precision, dimension(4,28,28) :: coefficients
  integer                              :: i
  data(coefficients(i, 2, 2),i=1,4) /1.90d-03, 3.00d-01, 4.70d+05, 9.40d+04/
  data(coefficients(i, 6, 6),i=1,4) /2.54d-03, 4.42d-02, 1.57d+05, 3.74d+05/
  data(coefficients(i, 6, 5),i=1,4) /6.15d-03, 5.88d-02, 1.41d+05, 1.41d+05/
  data(coefficients(i, 6, 4),i=1,4) /1.62d-03, 3.43d-01, 8.19d+04, 1.59d+05/
  data(coefficients(i, 6, 3),i=1,4) /4.78d-02, 3.62d-01, 3.44d+06, 5.87d+05/
  data(coefficients(i, 6, 2),i=1,4) /3.22d-02, 3.15d-01, 4.06d+06, 8.31d+05/
  data(coefficients(i, 7, 7),i=1,4) /2.98d-03, 0.00d+00, 2.20d+05, 1.00d+05/
  data(coefficients(i, 7, 6),i=1,4) /7.41d-03, 7.64d-02, 2.01d+05, 7.37d+04/
  data(coefficients(i, 7, 5),i=1,4) /1.13d-02, 1.64d-01, 1.72d+05, 2.25d+05/
  data(coefficients(i, 7, 4),i=1,4) /2.62d-03, 2.43d-01, 1.02d+05, 1.25d+05/
  data(coefficients(i, 7, 3),i=1,4) /7.50d-02, 3.50d-01, 4.75d+06, 8.35d+05/
  data(coefficients(i, 7, 2),i=1,4) /4.61d-02, 3.09d-01, 5.44d+06, 1.14d+06/
  data(coefficients(i, 8, 8),i=1,4) /1.11d-03, 9.25d-02, 1.75d+05, 1.45d+05/
  data(coefficients(i, 8, 7),i=1,4) /5.07d-03, 1.81d-01, 1.98d+05, 3.35d+05/
  data(coefficients(i, 8, 6),i=1,4) /1.48d-02, 3.05d-01, 2.41d+05, 2.83d+05/
  data(coefficients(i, 8, 5),i=1,4) /1.84d-02, 1.00d-01, 2.12d+05, 2.83d+05/
  data(coefficients(i, 8, 4),i=1,4) /4.13d-03, 1.62d-01, 1.25d+05, 2.27d+05/
  data(coefficients(i, 8, 3),i=1,4) /1.06d-01, 3.40d-01, 6.25d+06, 1.12d+06/
  data(coefficients(i, 8, 2),i=1,4) /4.72d-02, 0.00d+00, 6.00d+06, 1.00d+05/
  data(coefficients(i,10,10),i=1,4) /9.77d-04, 7.30d-02, 3.11d+05, 2.06d+05/
  data(coefficients(i,10, 9),i=1,4) /2.65d-03, 2.42d-01, 2.84d+05, 3.07d+05/
  data(coefficients(i,10, 8),i=1,4) /3.69d-03, 1.01d+00, 2.24d+05, 2.94d+05/
  data(coefficients(i,10, 7),i=1,4) /1.18d-02, 3.91d-01, 2.70d+05, 5.50d+05/
  data(coefficients(i,10, 6),i=1,4) /2.44d-02, 2.52d+00, 3.09d+05, 9.91d+05/
  data(coefficients(i,10, 5),i=1,4) /3.02d-02, 4.45d-01, 2.83d+05, 1.73d+06/
  data(coefficients(i,10, 4),i=1,4) /6.10d-03, 2.54d-01, 1.68d+05, 6.13d+05/
  data(coefficients(i,10, 3),i=1,4) /2.52d-01, 3.04d-01, 1.40d+07, 1.80d+06/
  data(coefficients(i,10, 2),i=1,4) /7.14d-02, 2.96d-01, 1.10d+07, 2.24d+06/
  data(coefficients(i,12,12),i=1,4) /4.49d-04, 2.10d-02, 5.01d+04, 2.81d+04/
  data(coefficients(i,12,11),i=1,4) /1.95d-03, 7.40d-02, 6.06d+05, 1.44d+06/
  data(coefficients(i,12,10),i=1,4) /5.12d-03, 3.23d-01, 4.69d+05, 7.55d+05/
  data(coefficients(i,12, 9),i=1,4) /7.74d-03, 6.36d-01, 3.74d+05, 7.88d+05/
  data(coefficients(i,12, 8),i=1,4) /1.17d-02, 8.07d-01, 3.28d+05, 1.02d+06/
  data(coefficients(i,12, 7),i=1,4) /3.69d-02, 3.51d-01, 4.80d+05, 9.73d+05/
  data(coefficients(i,12, 6),i=1,4) /3.63d-02, 5.48d-01, 3.88d+05, 7.36d+05/
  data(coefficients(i,12, 5),i=1,4) /4.15d-02, 2.33d-01, 3.39d+05, 3.82d+05/
  data(coefficients(i,12, 4),i=1,4) /8.86d-03, 3.18d-01, 2.11d+05, 1.54d+06/
  data(coefficients(i,12, 3),i=1,4) /2.52d-01, 3.15d-01, 1.40d+07, 2.64d+06/
  data(coefficients(i,12, 2),i=1,4) /9.28d-02, 0.00d+00, 1.45d+07, 1.00d+05/
  data(coefficients(i,13,13),i=1,4) /2.00d-03, 1.34d+00, 4.93d+04, 1.01d+05/
  data(coefficients(i,13,12),i=1,4) /2.35d-03, 1.25d-01, 6.72d+04, 2.45d+04/
  data(coefficients(i,13,11),i=1,4) /2.97d-03, 3.49d-01, 7.90d+05, 9.20d+05/
  data(coefficients(i,13,10),i=1,4) /6.97d-03, 1.12d-01, 8.42d+05, 3.34d+05/
  data(coefficients(i,13, 9),i=1,4) /1.22d-02, 1.29d+00, 4.29d+05, 8.42d+05/
  data(coefficients(i,13, 8),i=1,4) /1.84d-02, 1.36d+00, 3.72d+05, 1.45d+06/
  data(coefficients(i,13, 7),i=1,4) /3.32d-02, 1.24d+00, 4.17d+05, 1.55d+06/
  data(coefficients(i,13, 6),i=1,4) /3.78d-02, 8.02d-01, 3.66d+05, 1.41d+06/
  data(coefficients(i,13, 5),i=1,4) /5.12d-02, 1.72d-01, 3.59d+05, 3.39d+05/
  data(coefficients(i,13, 4),i=1,4) /1.15d-02, 9.43d-01, 2.31d+05, 2.97d+06/
  data(coefficients(i,13, 3),i=1,4) /2.85d-01, 3.12d-01, 1.61d+07, 3.11d+06/
  data(coefficients(i,13, 2),i=1,4) /1.02d-01, 1.25d-01, 1.71d+07, 1.91d+06/
  data(coefficients(i,14,14),i=1,4) /1.10d-03, 0.00d+00, 7.70d+04, 1.00d+05/
  data(coefficients(i,14,13),i=1,4) /5.87d-03, 7.53d-01, 9.63d+04, 6.46d+04/
  data(coefficients(i,14,12),i=1,4) /5.03d-03, 1.88d-01, 8.75d+04, 4.71d+04/
  data(coefficients(i,14,11),i=1,4) /5.43d-03, 4.50d-01, 1.05d+06, 7.98d+05/
  data(coefficients(i,14,10),i=1,4) /8.86d-03, 0.00d+00, 1.14d+06, 1.00d+05/
  data(coefficients(i,14, 9),i=1,4) /1.68d-02, 1.80d+00, 4.85d+05, 1.03d+06/
  data(coefficients(i,14, 8),i=1,4) /2.49d-02, 1.88d+00, 4.15d+05, 1.91d+06/
  data(coefficients(i,14, 7),i=1,4) /3.13d-02, 2.01d+00, 3.66d+05, 2.11d+06/
  data(coefficients(i,14, 6),i=1,4) /4.25d-02, 1.22d+00, 3.63d+05, 2.14d+06/
  data(coefficients(i,14, 5),i=1,4) /6.18d-02, 3.03d-01, 3.88d+05, 1.12d+06/
  data(coefficients(i,14, 4),i=1,4) /1.38d-02, 1.42d+00, 2.51d+05, 3.93d+06/
  data(coefficients(i,14, 3),i=1,4) /3.27d-01, 3.06d-01, 1.88d+07, 3.60d+06/
  data(coefficients(i,14, 2),i=1,4) /1.13d-01, 2.86d-01, 1.99d+07, 4.14d+06/
  data(coefficients(i,16,16),i=1,4) /1.62d-03, 0.00d+00, 1.25d+05, 1.00d+05/
  data(coefficients(i,16,15),i=1,4) /1.09d-02, 1.20d-02, 1.92d+05, 1.80d+04/
  data(coefficients(i,16,14),i=1,4) /3.35d-02, 6.59d-02, 1.89d+05, 1.59d+05/
  data(coefficients(i,16,13),i=1,4) /3.14d-02, 6.89d-02, 1.68d+05, 8.04d+04/
  data(coefficients(i,16,12),i=1,4) /1.27d-02, 1.87d-01, 1.38d+05, 1.71d+05/
  data(coefficients(i,16,11),i=1,4) /1.47d-02, 1.29d-01, 1.80d+06, 1.75d+06/
  data(coefficients(i,16,10),i=1,4) /1.34d-02, 1.04d+00, 6.90d+05, 2.15d+06/
  data(coefficients(i,16, 9),i=1,4) /2.38d-02, 1.12d+00, 5.84d+05, 2.59d+06/
  data(coefficients(i,16, 8),i=1,4) /3.19d-02, 1.40d+00, 5.17d+05, 2.91d+06/
  data(coefficients(i,16, 7),i=1,4) /7.13d-02, 1.00d+00, 6.66d+05, 2.32d+06/
  data(coefficients(i,16, 6),i=1,4) /8.00d-02, 5.55d-01, 6.00d+05, 2.41d+06/
  data(coefficients(i,16, 5),i=1,4) /7.96d-02, 1.63d+00, 5.09d+05, 6.37d+06/
  data(coefficients(i,16, 4),i=1,4) /1.34d-02, 3.04d-01, 2.91d+05, 1.04d+06/
  data(coefficients(i,16, 3),i=1,4) /4.02d-01, 2.98d-01, 2.41d+07, 4.67d+06/
  data(coefficients(i,16, 2),i=1,4) /1.45d-01, 2.81d-01, 2.54d+07, 5.30d+06/
  data(coefficients(i,18,18),i=1,4) /1.00d-03, 5.00d-03, 3.20d+05, 3.10d+05/
  data(coefficients(i,18,17),i=1,4) /1.10d-02, 4.50d-02, 2.90d+05, 5.50d+05/
  data(coefficients(i,18,16),i=1,4) /3.40d-02, 5.70d-02, 2.39d+05, 6.00d+05/
  data(coefficients(i,18,15),i=1,4) /6.85d-02, 8.70d-02, 2.56d+05, 3.81d+05/
  data(coefficients(i,18,14),i=1,4) /9.00d-02, 7.69d-02, 2.50d+05, 3.30d+05/
  data(coefficients(i,18,13),i=1,4) /6.35d-02, 1.40d-01, 2.10d+05, 2.15d+05/
  data(coefficients(i,18,12),i=1,4) /2.60d-02, 1.20d-01, 1.80d+05, 2.15d+05/
  data(coefficients(i,18,11),i=1,4) /1.70d-02, 1.00d-01, 2.70d+06, 3.30d+06/
  data(coefficients(i,18,10),i=1,4) /2.10d-02, 1.92d+00, 8.30d+05, 3.50d+06/
  data(coefficients(i,18, 9),i=1,4) /3.50d-02, 1.66d+00, 6.95d+05, 3.60d+06/
  data(coefficients(i,18, 8),i=1,4) /4.30d-02, 1.67d+00, 6.05d+05, 3.80d+06/
  data(coefficients(i,18, 7),i=1,4) /7.13d-02, 1.40d+00, 6.68d+05, 2.90d+06/
  data(coefficients(i,18, 6),i=1,4) /9.60d-02, 1.31d+00, 6.50d+05, 3.60d+06/
  data(coefficients(i,18, 5),i=1,4) /8.50d-02, 1.02d+00, 5.30d+05, 2.80d+06/
  data(coefficients(i,18, 4),i=1,4) /1.70d-02, 2.45d-01, 3.55d+05, 1.10d+06/
  data(coefficients(i,18, 3),i=1,4) /4.76d-01, 2.94d-01, 3.01d+07, 6.05d+06/
  data(coefficients(i,18, 2),i=1,4) /2.97d-01, 2.77d-01, 3.13d+07, 6.54d+06/
  data(coefficients(i,20,20),i=1,4) /3.28d-04, 9.07d-02, 3.46d+04, 1.64d+04/
  data(coefficients(i,20,19),i=1,4) /5.84d-02, 1.10d-01, 3.85d+05, 2.45d+05/
  data(coefficients(i,20,18),i=1,4) /1.12d-01, 1.74d-02, 4.08d+05, 4.27d+05/
  data(coefficients(i,20,17),i=1,4) /1.32d-01, 1.32d-01, 3.82d+05, 6.92d+05/
  data(coefficients(i,20,16),i=1,4) /1.33d-01, 1.14d-01, 3.53d+05, 8.78d+05/
  data(coefficients(i,20,15),i=1,4) /1.26d-01, 1.62d-01, 3.19d+05, 7.43d+05/
  data(coefficients(i,20,14),i=1,4) /1.39d-01, 8.78d-02, 3.22d+05, 6.99d+05/
  data(coefficients(i,20,13),i=1,4) /9.55d-02, 2.63d-01, 2.47d+05, 4.43d+05/
  data(coefficients(i,20,12),i=1,4) /4.02d-02, 6.27d-02, 2.29d+05, 2.81d+05/
  data(coefficients(i,20,11),i=1,4) /4.19d-02, 6.16d-02, 3.73d+06, 5.84d+06/
  data(coefficients(i,20,10),i=1,4) /2.57d-02, 2.77d+00, 9.26d+05, 4.89d+06/
  data(coefficients(i,20, 9),i=1,4) /4.45d-02, 2.23d+00, 7.96d+05, 4.62d+06/
  data(coefficients(i,20, 8),i=1,4) /5.48d-02, 2.00d+00, 6.90d+05, 4.52d+06/
  data(coefficients(i,20, 7),i=1,4) /7.13d-02, 1.82d+00, 6.70d+05, 3.32d+06/
  data(coefficients(i,20, 6),i=1,4) /1.09d-01, 1.74d+00, 7.00d+05, 4.93d+06/
  data(coefficients(i,20, 5),i=1,4) /1.10d-01, 2.43d-01, 5.67d+05, 4.41d+06/
  data(coefficients(i,20, 4),i=1,4) /2.05d-02, 1.85d-01, 4.21d+05, 2.27d+06/
  data(coefficients(i,20, 3),i=1,4) /5.49d-01, 2.92d-01, 3.65d+07, 7.25d+06/
  data(coefficients(i,20, 2),i=1,4) /2.68d-01, 0.00d+00, 3.74d+07, 1.00d+05/
  data(coefficients(i,26,26),i=1,4) /1.58d-03, 4.56d-01, 6.00d+04, 8.97d+04/
  data(coefficients(i,26,25),i=1,4) /8.38d-03, 3.23d-01, 1.94d+05, 1.71d+05/
  data(coefficients(i,26,24),i=1,4) /1.54d-02, 3.10d-01, 3.31d+05, 2.73d+05/
  data(coefficients(i,26,23),i=1,4) /3.75d-02, 4.11d-01, 4.32d+05, 3.49d+05/
  data(coefficients(i,26,22),i=1,4) /1.17d-01, 3.59d-01, 6.28d+05, 5.29d+05/
  data(coefficients(i,26,21),i=1,4) /2.54d-01, 9.75d-02, 7.50d+05, 4.69d+05/
  data(coefficients(i,26,20),i=1,4) /2.91d-01, 2.29d-01, 7.73d+05, 6.54d+05/
  data(coefficients(i,26,19),i=1,4) /1.50d-01, 4.20d+00, 2.62d+05, 1.32d+06/
  data(coefficients(i,26,18),i=1,4) /1.40d-01, 3.30d+00, 2.50d+05, 1.33d+06/
  data(coefficients(i,26,17),i=1,4) /1.00d-01, 5.30d+00, 2.57d+05, 1.41d+06/
  data(coefficients(i,26,16),i=1,4) /2.00d-01, 1.50d+00, 2.84d+05, 1.52d+06/
  data(coefficients(i,26,15),i=1,4) /2.40d-01, 7.00d-01, 8.69d+05, 1.51d+06/
  data(coefficients(i,26,14),i=1,4) /2.60d-01, 6.00d-01, 4.21d+05, 1.82d+06/
  data(coefficients(i,26,13),i=1,4) /1.90d-01, 5.00d-01, 4.57d+05, 1.84d+06/
  data(coefficients(i,26,12),i=1,4) /1.20d-01, 1.00d+00, 2.85d+05, 2.31d+06/
  data(coefficients(i,26,11),i=1,4) /3.50d-01, 0.00d+00, 8.18d+06, 1.00d+05/
  data(coefficients(i,26,10),i=1,4) /6.60d-02, 7.80d+00, 1.51d+06, 9.98d+06/
  data(coefficients(i,26, 9),i=1,4) /1.00d-01, 6.30d+00, 1.30d+06, 9.98d+06/
  data(coefficients(i,26, 8),i=1,4) /1.30d-01, 5.50d+00, 1.19d+06, 1.00d+07/
  data(coefficients(i,26, 7),i=1,4) /2.30d-01, 3.60d+00, 1.09d+06, 1.10d+07/
  data(coefficients(i,26, 6),i=1,4) /1.40d-01, 4.90d+00, 9.62d+05, 8.34d+06/
  data(coefficients(i,26, 5),i=1,4) /1.10d-01, 1.60d+00, 7.23d+05, 1.01d+07/
  data(coefficients(i,26, 4),i=1,4) /4.10d-02, 4.20d+00, 4.23d+05, 1.07d+07/
  data(coefficients(i,26, 3),i=1,4) /7.47d-01, 2.84d-01, 5.84d+07, 1.17d+07/
  data(coefficients(i,26, 2),i=1,4) /3.69d-01, 0.00d+00, 6.00d+07, 1.00d+05/
  data(coefficients(i,28,28),i=1,4) /1.41d-03, 4.69d-01, 9.82d+04, 1.01d+05/
  data(coefficients(i,28,27),i=1,4) /5.20d-03, 3.57d-01, 2.01d+05, 1.91d+05/
  data(coefficients(i,28,26),i=1,4) /1.38d-02, 2.81d-01, 3.05d+05, 2.32d+05/
  data(coefficients(i,28,25),i=1,4) /2.30d-02, 1.28d-01, 4.20d+05, 3.18d+05/
  data(coefficients(i,28,24),i=1,4) /4.19d-02, 4.17d-02, 5.56d+05, 4.55d+05/
  data(coefficients(i,28,23),i=1,4) /6.83d-02, 5.58d-02, 6.72d+05, 5.51d+05/
  data(coefficients(i,28,22),i=1,4) /1.22d-01, 3.46d-02, 7.93d+05, 5.28d+05/
  data(coefficients(i,28,21),i=1,4) /3.00d-01, 0.00d+00, 9.00d+05, 1.00d+05/
  data(coefficients(i,28,20),i=1,4) /1.50d-01, 1.90d+00, 1.00d+06, 5.50d+05/
  data(coefficients(i,28,19),i=1,4) /6.97d-01, 2.77d-01, 7.81d+05, 8.87d+05/
  data(coefficients(i,28,18),i=1,4) /7.09d-01, 1.35d-01, 7.64d+05, 1.80d+06/
  data(coefficients(i,28,17),i=1,4) /6.44d-01, 1.34d-01, 7.44d+05, 1.25d+06/
  data(coefficients(i,28,16),i=1,4) /5.25d-01, 1.92d-01, 6.65d+05, 1.89d+06/
  data(coefficients(i,28,15),i=1,4) /4.46d-01, 3.32d-01, 5.97d+05, 8.84d+05/
  data(coefficients(i,28,14),i=1,4) /3.63d-01, 3.37d-01, 5.24d+05, 1.29d+06/
  data(coefficients(i,28,13),i=1,4) /3.02d-01, 1.21d-01, 4.96d+05, 6.24d+05/
  data(coefficients(i,28,12),i=1,4) /1.02d-01, 5.14d-02, 4.46d+05, 1.59d+06/
  data(coefficients(i,28,11),i=1,4) /2.70d-01, 1.83d-01, 8.49d+06, 8.01d+06/
  data(coefficients(i,28,10),i=1,4) /4.67d-02, 7.56d+00, 1.36d+06, 9.32d+06/
  data(coefficients(i,28, 9),i=1,4) /8.35d-02, 4.55d+00, 1.23d+06, 9.45d+06/
  data(coefficients(i,28, 8),i=1,4) /9.96d-02, 4.87d+00, 1.06d+06, 9.45d+06/
  data(coefficients(i,28, 7),i=1,4) /1.99d-01, 2.19d+00, 1.25d+06, 8.01d+06/
  data(coefficients(i,28, 6),i=1,4) /2.40d-01, 1.15d+00, 1.23d+06, 7.57d+06/
  data(coefficients(i,28, 5),i=1,4) /1.15d-01, 1.23d+00, 3.32d+05, 2.64d+06/
  data(coefficients(i,28, 4),i=1,4) /3.16d-02, 1.32d-01, 6.45d+05, 1.93d+06/
  data(coefficients(i,28, 3),i=1,4) /8.03d-01, 2.89d-01, 6.65d+07, 1.19d+07/
  data(coefficients(i,28, 2),i=1,4) /5.75d-01, 2.86d-01, 6.81d+07, 9.08d+06/

contains

  function arnaud1985ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{atomicRecombinationRateDielectronicArnaud1985} atomic ionization potential class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(atomicRecombinationRateDielectronicArnaud1985)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=atomicRecombinationRateDielectronicArnaud1985()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function arnaud1985ConstructorParameters

  double precision function arnaud1985Rate(self,atomicNumber,electronNumber,temperature)
    !!{
    Calculates rates of dielectric recombination for all ionization stages of all elements from H to Ni ($Z=28$) by use of the
    fits from \cite{aldrovandi_radiative_1973}, \cite{shull_ionization_1982} and \cite{arnaud_updated_1985}.  Input parameters:
    {\normalfont \ttfamily atomicNumber}: atomic number; {\normalfont \ttfamily electronNumber}: number of electrons;
    {\normalfont \ttfamily temperature}: temperature [K].  Output parameter: rate coefficient [cm$^3$ s$^{-1}$].
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (atomicRecombinationRateDielectronicArnaud1985), intent(inout) :: self
    double precision                                               , intent(in   ) :: temperature
    integer                                                        , intent(in   ) :: atomicNumber, electronNumber
    !$GLC attributes unused :: self

    ! Set zero rate by default.
    arnaud1985Rate=0.0d0
    ! Return on unphysical temperature.
    if (temperature <= 0.0) return
    ! Abort on unphysical conditions.
    if (atomicNumber   < 1 .or. atomicNumber   >           28) call Error_Report('atomic number is unphysical or too large'  //{introspection:location})
    if (electronNumber < 1 .or. electronNumber > atomicNumber) call Error_Report('electron number is unphysical or too large'//{introspection:location})
    ! Evaluate the fitting function.
    arnaud1985Rate=+       coefficients(1,atomicNumber,electronNumber)              &
         &         *exp  (-coefficients(3,atomicNumber,electronNumber)/temperature) &
         &         *(                                                               &
         &           +1.0d0                                                         &
         &           +     coefficients(2,atomicNumber,electronNumber)              &
         &           *exp(-coefficients(4,atomicNumber,electronNumber)/temperature) &
         &         )                                                                &
         &         /temperature**1.5d0
    return
  end function arnaud1985Rate
