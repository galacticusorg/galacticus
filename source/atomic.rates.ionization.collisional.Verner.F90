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
  An implementation of atomic collisional ionization rates based on the \href{https://web.archive.org/web/20220313133801/https://www.pa.uky.edu/~verner/dima/col/cfit.f}{code}
  originally written by Dima Verner.
  !!}

  !![
  <atomicIonizationRateCollisional name="atomicIonizationRateCollisionalVerner1996">
   <description>Atomic collisional ionization rates are computed based on the \href{https://web.archive.org/web/20220313133801/https://www.pa.uky.edu/~verner/dima/col/cfit.f}{code} originally written by Dima Verner.</description>
  </atomicIonizationRateCollisional>
  !!]
  type, extends(atomicIonizationRateCollisionalClass) :: atomicIonizationRateCollisionalVerner1996
     !!{
     A collisional ionization rate class based on the \href{https://web.archive.org/web/20220313133801/https://www.pa.uky.edu/~verner/dima/col/cfit.f}{code} originally written by Dima Verner.
     !!}
     private
   contains
     procedure :: rate => verner1996Rate
  end type atomicIonizationRateCollisionalVerner1996

  interface atomicIonizationRateCollisionalVerner1996
     !!{
     Constructors for the \refClass{atomicIonizationRateCollisionalVerner1996} atomic collisional ionization class.
     !!}
     module procedure verner1996ConstructorParameters
  end interface atomicIonizationRateCollisionalVerner1996

  ! Arrays to store coefficients of collisional ionization rate fitting functions.
  double precision :: fitCoefficient(5,28,28)
  integer          :: i

  ! Initialize the fitting coefficients.
  data(fitCoefficient(i, 1, 1),i=1,5) /   13.6d0,0.0d0,2.91d-08,0.2320d0,0.39d0/
  data(fitCoefficient(i, 2, 2),i=1,5) /   24.6d0,0.0d0,1.75d-08,0.1800d0,0.35d0/
  data(fitCoefficient(i, 2, 1),i=1,5) /   54.4d0,1.0d0,2.05d-09,0.2650d0,0.25d0/
  data(fitCoefficient(i, 3, 3),i=1,5) /    5.4d0,0.0d0,1.39d-07,0.4380d0,0.41d0/
  data(fitCoefficient(i, 3, 2),i=1,5) /   75.6d0,1.0d0,2.01d-09,0.2090d0,0.23d0/
  data(fitCoefficient(i, 3, 1),i=1,5) /  122.4d0,1.0d0,9.60d-10,0.5820d0,0.17d0/
  data(fitCoefficient(i, 4, 4),i=1,5) /    9.3d0,0.0d0,1.02d-07,0.3750d0,0.27d0/
  data(fitCoefficient(i, 4, 3),i=1,5) /   18.2d0,1.0d0,2.08d-08,0.4390d0,0.21d0/
  data(fitCoefficient(i, 4, 2),i=1,5) /  153.9d0,0.0d0,2.67d-09,0.6120d0,0.27d0/
  data(fitCoefficient(i, 4, 1),i=1,5) /  217.7d0,1.0d0,4.27d-10,0.6580d0,0.15d0/
  data(fitCoefficient(i, 5, 5),i=1,5) /    8.3d0,0.0d0,6.49d-08,0.2000d0,0.26d0/
  data(fitCoefficient(i, 5, 4),i=1,5) /   25.2d0,1.0d0,1.24d-08,0.2670d0,0.22d0/
  data(fitCoefficient(i, 5, 3),i=1,5) /   37.9d0,1.0d0,3.27d-09,0.2950d0,0.23d0/
  data(fitCoefficient(i, 5, 2),i=1,5) /  259.4d0,1.0d0,4.95d-10,0.4890d0,0.09d0/
  data(fitCoefficient(i, 5, 1),i=1,5) /  340.2d0,1.0d0,2.19d-10,0.6570d0,0.15d0/
  data(fitCoefficient(i, 6, 6),i=1,5) /   11.3d0,0.0d0,6.85d-08,0.1930d0,0.25d0/
  data(fitCoefficient(i, 6, 5),i=1,5) /   24.4d0,1.0d0,1.86d-08,0.2860d0,0.24d0/
  data(fitCoefficient(i, 6, 4),i=1,5) /   47.9d0,1.0d0,6.35d-09,0.4270d0,0.21d0/
  data(fitCoefficient(i, 6, 3),i=1,5) /   64.5d0,1.0d0,1.50d-09,0.4160d0,0.13d0/
  data(fitCoefficient(i, 6, 2),i=1,5) /  392.1d0,1.0d0,2.99d-10,0.6660d0,0.02d0/
  data(fitCoefficient(i, 6, 1),i=1,5) /  490.0d0,1.0d0,1.23d-10,0.6200d0,0.16d0/
  data(fitCoefficient(i, 7, 7),i=1,5) /   14.5d0,0.0d0,4.82d-08,0.0652d0,0.42d0/
  data(fitCoefficient(i, 7, 6),i=1,5) /   29.6d0,0.0d0,2.98d-08,0.3100d0,0.30d0/
  data(fitCoefficient(i, 7, 5),i=1,5) /   47.5d0,1.0d0,8.10d-09,0.3500d0,0.24d0/
  data(fitCoefficient(i, 7, 4),i=1,5) /   77.5d0,1.0d0,3.71d-09,0.5490d0,0.18d0/
  data(fitCoefficient(i, 7, 3),i=1,5) /   97.9d0,0.0d0,1.51d-09,0.0167d0,0.74d0/
  data(fitCoefficient(i, 7, 2),i=1,5) /  552.1d0,0.0d0,3.71d-10,0.5460d0,0.29d0/
  data(fitCoefficient(i, 7, 1),i=1,5) /  667.0d0,1.0d0,7.77d-11,0.6240d0,0.16d0/
  data(fitCoefficient(i, 8, 8),i=1,5) /   13.6d0,0.0d0,3.59d-08,0.0730d0,0.34d0/
  data(fitCoefficient(i, 8, 7),i=1,5) /   35.1d0,1.0d0,1.39d-08,0.2120d0,0.22d0/
  data(fitCoefficient(i, 8, 6),i=1,5) /   54.9d0,1.0d0,9.31d-09,0.2700d0,0.27d0/
  data(fitCoefficient(i, 8, 5),i=1,5) /   77.4d0,0.0d0,1.02d-08,0.6140d0,0.27d0/
  data(fitCoefficient(i, 8, 4),i=1,5) /  113.9d0,1.0d0,2.19d-09,0.6300d0,0.17d0/
  data(fitCoefficient(i, 8, 3),i=1,5) /  138.1d0,0.0d0,1.95d-09,0.3600d0,0.54d0/
  data(fitCoefficient(i, 8, 2),i=1,5) /  739.3d0,0.0d0,2.12d-10,0.3960d0,0.35d0/
  data(fitCoefficient(i, 8, 1),i=1,5) /  871.4d0,1.0d0,5.21d-11,0.6290d0,0.16d0/
  data(fitCoefficient(i, 9, 9),i=1,5) /   17.4d0,1.0d0,7.00d-08,0.1780d0,0.29d0/
  data(fitCoefficient(i, 9, 8),i=1,5) /   35.0d0,0.0d0,5.41d-08,0.5710d0,0.27d0/
  data(fitCoefficient(i, 9, 7),i=1,5) /   62.7d0,1.0d0,9.37d-09,0.3190d0,0.20d0/
  data(fitCoefficient(i, 9, 6),i=1,5) /   87.1d0,1.0d0,4.92d-09,0.3230d0,0.24d0/
  data(fitCoefficient(i, 9, 5),i=1,5) /  114.2d0,0.0d0,7.06d-09,0.6840d0,0.27d0/
  data(fitCoefficient(i, 9, 4),i=1,5) /  157.2d0,1.0d0,1.28d-09,0.6480d0,0.16d0/
  data(fitCoefficient(i, 9, 3),i=1,5) /  185.2d0,1.0d0,5.61d-10,0.7380d0,0.16d0/
  data(fitCoefficient(i, 9, 2),i=1,5) /  953.9d0,0.0d0,1.66d-10,0.5420d0,0.29d0/
  data(fitCoefficient(i, 9, 1),i=1,5) / 1103.1d0,1.0d0,3.74d-11,0.6590d0,0.15d0/
  data(fitCoefficient(i,10,10),i=1,5) /   21.6d0,1.0d0,1.50d-08,0.0329d0,0.43d0/
  data(fitCoefficient(i,10, 9),i=1,5) /   41.0d0,0.0d0,1.98d-08,0.2950d0,0.20d0/
  data(fitCoefficient(i,10, 8),i=1,5) /   63.5d0,1.0d0,7.03d-09,0.0677d0,0.39d0/
  data(fitCoefficient(i,10, 7),i=1,5) /   97.1d0,1.0d0,4.24d-09,0.0482d0,0.58d0/
  data(fitCoefficient(i,10, 6),i=1,5) /  126.2d0,1.0d0,2.79d-09,0.3050d0,0.25d0/
  data(fitCoefficient(i,10, 5),i=1,5) /  157.9d0,0.0d0,3.45d-09,0.5810d0,0.28d0/
  data(fitCoefficient(i,10, 4),i=1,5) /  207.3d0,1.0d0,9.56d-10,0.7490d0,0.14d0/
  data(fitCoefficient(i,10, 3),i=1,5) /  239.1d0,1.0d0,4.73d-10,0.9920d0,0.04d0/
  data(fitCoefficient(i,10, 2),i=1,5) / 1196.0d0,1.0d0,3.92d-11,0.2620d0,0.20d0/
  data(fitCoefficient(i,10, 1),i=1,5) / 1360.6d0,1.0d0,2.77d-11,0.6610d0,0.13d0/
  data(fitCoefficient(i,11,11),i=1,5) /    5.1d0,1.0d0,1.01d-07,0.2750d0,0.23d0/
  data(fitCoefficient(i,11,10),i=1,5) /   47.3d0,1.0d0,7.35d-09,0.0560d0,0.35d0/
  data(fitCoefficient(i,11, 9),i=1,5) /   71.6d0,1.0d0,8.10d-09,0.1480d0,0.32d0/
  data(fitCoefficient(i,11, 8),i=1,5) /   98.9d0,0.0d0,1.14d-08,0.5530d0,0.28d0/
  data(fitCoefficient(i,11, 7),i=1,5) /  138.4d0,1.0d0,2.63d-09,0.2300d0,0.29d0/
  data(fitCoefficient(i,11, 6),i=1,5) /  172.2d0,1.0d0,1.85d-09,0.3630d0,0.22d0/
  data(fitCoefficient(i,11, 5),i=1,5) /  208.5d0,0.0d0,2.82d-09,0.6740d0,0.27d0/
  data(fitCoefficient(i,11, 4),i=1,5) /  264.2d0,1.0d0,6.72d-10,0.7520d0,0.14d0/
  data(fitCoefficient(i,11, 3),i=1,5) /  299.9d0,1.0d0,2.80d-10,0.7810d0,0.15d0/
  data(fitCoefficient(i,11, 2),i=1,5) / 1465.1d0,1.0d0,4.63d-11,0.5580d0,0.16d0/
  data(fitCoefficient(i,11, 1),i=1,5) / 1648.7d0,1.0d0,2.16d-11,0.7430d0,0.13d0/
  data(fitCoefficient(i,12,12),i=1,5) /    7.6d0,0.0d0,6.21d-07,0.5920d0,0.39d0/
  data(fitCoefficient(i,12,11),i=1,5) /   15.2d0,0.0d0,1.92d-08,0.0027d0,0.85d0/
  data(fitCoefficient(i,12,10),i=1,5) /   80.1d0,1.0d0,5.56d-09,0.1070d0,0.30d0/
  data(fitCoefficient(i,12, 9),i=1,5) /  109.3d0,1.0d0,4.35d-09,0.1590d0,0.31d0/
  data(fitCoefficient(i,12, 8),i=1,5) /  141.3d0,0.0d0,7.10d-09,0.6580d0,0.25d0/
  data(fitCoefficient(i,12, 7),i=1,5) /  186.5d0,1.0d0,1.70d-09,0.2420d0,0.28d0/
  data(fitCoefficient(i,12, 6),i=1,5) /  224.9d0,1.0d0,1.22d-09,0.3430d0,0.23d0/
  data(fitCoefficient(i,12, 5),i=1,5) /  266.0d0,0.0d0,2.20d-09,0.8970d0,0.22d0/
  data(fitCoefficient(i,12, 4),i=1,5) /  328.2d0,1.0d0,4.86d-10,0.7510d0,0.14d0/
  data(fitCoefficient(i,12, 3),i=1,5) /  367.5d0,1.0d0,2.35d-10,1.0300d0,0.10d0/
  data(fitCoefficient(i,12, 2),i=1,5) / 1761.8d0,1.0d0,2.06d-11,0.1960d0,0.25d0/
  data(fitCoefficient(i,12, 1),i=1,5) / 1962.7d0,1.0d0,1.75d-11,0.8350d0,0.11d0/
  data(fitCoefficient(i,13,13),i=1,5) /    6.0d0,1.0d0,2.28d-07,0.3870d0,0.25d0/
  data(fitCoefficient(i,13,12),i=1,5) /   18.8d0,0.0d0,1.18d-07,2.2100d0,0.25d0/
  data(fitCoefficient(i,13,11),i=1,5) /   28.5d0,1.0d0,4.40d-09,0.1060d0,0.24d0/
  data(fitCoefficient(i,13,10),i=1,5) /  120.0d0,0.0d0,1.75d-08,0.8720d0,0.22d0/
  data(fitCoefficient(i,13, 9),i=1,5) /  153.8d0,1.0d0,2.61d-09,0.1590d0,0.31d0/
  data(fitCoefficient(i,13, 8),i=1,5) /  198.5d0,1.0d0,1.85d-09,0.1520d0,0.36d0/
  data(fitCoefficient(i,13, 7),i=1,5) /  241.4d0,1.0d0,1.14d-09,0.2280d0,0.29d0/
  data(fitCoefficient(i,13, 6),i=1,5) /  284.6d0,1.0d0,8.00d-10,0.4170d0,0.16d0/
  data(fitCoefficient(i,13, 5),i=1,5) /  390.2d0,1.0d0,5.83d-10,0.4970d0,0.23d0/
  data(fitCoefficient(i,13, 4),i=1,5) /  399.4d0,0.0d0,4.93d-10,0.7060d0,0.16d0/
  data(fitCoefficient(i,13, 3),i=1,5) /  442.0d0,1.0d0,9.77d-11,0.2780d0,0.17d0/
  data(fitCoefficient(i,13, 2),i=1,5) / 2086.6d0,0.0d0,3.94d-11,0.2860d0,0.36d0/
  data(fitCoefficient(i,13, 1),i=1,5) / 2304.1d0,1.0d0,1.38d-11,0.8350d0,0.11d0/
  data(fitCoefficient(i,14,14),i=1,5) /    8.2d0,1.0d0,1.88d-07,0.3760d0,0.25d0/
  data(fitCoefficient(i,14,13),i=1,5) /   16.4d0,1.0d0,6.43d-08,0.6320d0,0.20d0/
  data(fitCoefficient(i,14,12),i=1,5) /   33.5d0,1.0d0,2.01d-08,0.4730d0,0.22d0/
  data(fitCoefficient(i,14,11),i=1,5) /   54.0d0,1.0d0,4.94d-09,0.1720d0,0.23d0/
  data(fitCoefficient(i,14,10),i=1,5) /  166.8d0,1.0d0,1.76d-09,0.1020d0,0.31d0/
  data(fitCoefficient(i,14, 9),i=1,5) /  205.3d0,1.0d0,1.74d-09,0.1800d0,0.29d0/
  data(fitCoefficient(i,14, 8),i=1,5) /  246.5d0,1.0d0,1.23d-09,0.5180d0,0.07d0/
  data(fitCoefficient(i,14, 7),i=1,5) /  303.5d0,1.0d0,8.27d-10,0.2390d0,0.28d0/
  data(fitCoefficient(i,14, 6),i=1,5) /  351.1d0,1.0d0,6.01d-10,0.3050d0,0.25d0/
  data(fitCoefficient(i,14, 5),i=1,5) /  401.4d0,1.0d0,4.65d-10,0.6660d0,0.04d0/
  data(fitCoefficient(i,14, 4),i=1,5) /  476.4d0,1.0d0,2.63d-10,0.6660d0,0.16d0/
  data(fitCoefficient(i,14, 3),i=1,5) /  523.5d0,1.0d0,1.18d-10,0.7340d0,0.16d0/
  data(fitCoefficient(i,14, 2),i=1,5) / 2437.7d0,0.0d0,3.36d-11,0.3360d0,0.37d0/
  data(fitCoefficient(i,14, 1),i=1,5) / 2673.2d0,1.0d0,1.19d-11,0.9890d0,0.08d0/
  data(fitCoefficient(i,15,15),i=1,5) /   10.5d0,1.0d0,1.99d-07,0.5350d0,0.24d0/
  data(fitCoefficient(i,15,14),i=1,5) /   19.8d0,1.0d0,5.88d-08,0.5370d0,0.21d0/
  data(fitCoefficient(i,15,13),i=1,5) /   30.2d0,1.0d0,2.96d-08,0.8650d0,0.16d0/
  data(fitCoefficient(i,15,12),i=1,5) /   51.4d0,1.0d0,1.01d-08,0.5460d0,0.20d0/
  data(fitCoefficient(i,15,11),i=1,5) /   65.0d0,1.0d0,2.36d-09,0.1920d0,0.17d0/
  data(fitCoefficient(i,15,10),i=1,5) /  220.4d0,0.0d0,6.66d-09,1.0000d0,0.18d0/
  data(fitCoefficient(i,15, 9),i=1,5) /  263.2d0,1.0d0,1.24d-09,0.2150d0,0.26d0/
  data(fitCoefficient(i,15, 8),i=1,5) /  309.4d0,0.0d0,2.27d-09,0.7340d0,0.23d0/
  data(fitCoefficient(i,15, 7),i=1,5) /  371.7d0,1.0d0,6.14d-10,0.2560d0,0.27d0/
  data(fitCoefficient(i,15, 6),i=1,5) /  424.5d0,1.0d0,4.69d-10,0.3420d0,0.23d0/
  data(fitCoefficient(i,15, 5),i=1,5) /  479.6d0,0.0d0,6.14d-10,0.3340d0,0.39d0/
  data(fitCoefficient(i,15, 4),i=1,5) /  560.4d0,0.0d0,3.22d-10,0.8500d0,0.12d0/
  data(fitCoefficient(i,15, 3),i=1,5) /  611.9d0,1.0d0,9.32d-11,0.7340d0,0.16d0/
  data(fitCoefficient(i,15, 2),i=1,5) / 2816.9d0,0.0d0,3.79d-11,0.8050d0,0.22d0/
  data(fitCoefficient(i,15, 1),i=1,5) / 3069.9d0,1.0d0,9.73d-12,0.9910d0,0.08d0/
  data(fitCoefficient(i,16,16),i=1,5) /   10.4d0,1.0d0,5.49d-08,0.1000d0,0.25d0/
  data(fitCoefficient(i,16,15),i=1,5) /   23.3d0,1.0d0,6.81d-08,0.6930d0,0.21d0/
  data(fitCoefficient(i,16,14),i=1,5) /   34.8d0,1.0d0,2.14d-08,0.3530d0,0.24d0/
  data(fitCoefficient(i,16,13),i=1,5) /   47.3d0,1.0d0,1.66d-08,1.0300d0,0.14d0/
  data(fitCoefficient(i,16,12),i=1,5) /   72.6d0,1.0d0,6.12d-09,0.5800d0,0.19d0/
  data(fitCoefficient(i,16,11),i=1,5) /   88.1d0,1.0d0,1.33d-09,0.0688d0,0.35d0/
  data(fitCoefficient(i,16,10),i=1,5) /  280.9d0,0.0d0,4.93d-09,1.1300d0,0.16d0/
  data(fitCoefficient(i,16, 9),i=1,5) /  328.2d0,1.0d0,8.73d-10,0.1930d0,0.28d0/
  data(fitCoefficient(i,16, 8),i=1,5) /  379.1d0,0.0d0,1.35d-09,0.4310d0,0.32d0/
  data(fitCoefficient(i,16, 7),i=1,5) /  447.1d0,1.0d0,4.59d-10,0.2420d0,0.28d0/
  data(fitCoefficient(i,16, 6),i=1,5) /  504.8d0,1.0d0,3.49d-10,0.3050d0,0.25d0/
  data(fitCoefficient(i,16, 5),i=1,5) /  564.7d0,0.0d0,5.23d-10,0.4280d0,0.35d0/
  data(fitCoefficient(i,16, 4),i=1,5) /  651.6d0,0.0d0,2.59d-10,0.8540d0,0.12d0/
  data(fitCoefficient(i,16, 3),i=1,5) /  707.2d0,1.0d0,7.50d-11,0.7340d0,0.16d0/
  data(fitCoefficient(i,16, 2),i=1,5) / 3223.9d0,0.0d0,2.67d-11,0.5720d0,0.28d0/
  data(fitCoefficient(i,16, 1),i=1,5) / 3494.2d0,1.0d0,6.32d-12,0.5850d0,0.17d0/
  data(fitCoefficient(i,17,17),i=1,5) /   13.0d0,1.0d0,1.69d-07,0.4300d0,0.24d0/
  data(fitCoefficient(i,17,16),i=1,5) /   23.8d0,1.0d0,6.96d-08,0.6700d0,0.20d0/
  data(fitCoefficient(i,17,15),i=1,5) /   39.6d0,1.0d0,3.40d-08,0.8650d0,0.18d0/
  data(fitCoefficient(i,17,14),i=1,5) /   53.5d0,1.0d0,1.10d-08,0.3280d0,0.25d0/
  data(fitCoefficient(i,17,13),i=1,5) /   67.8d0,1.0d0,1.11d-08,1.3700d0,0.10d0/
  data(fitCoefficient(i,17,12),i=1,5) /   97.0d0,1.0d0,3.17d-09,0.3300d0,0.24d0/
  data(fitCoefficient(i,17,11),i=1,5) /  114.2d0,1.0d0,1.01d-09,0.1960d0,0.16d0/
  data(fitCoefficient(i,17,10),i=1,5) /  348.3d0,0.0d0,2.11d-09,0.3130d0,0.37d0/
  data(fitCoefficient(i,17, 9),i=1,5) /  400.1d0,1.0d0,6.32d-10,0.1730d0,0.30d0/
  data(fitCoefficient(i,17, 8),i=1,5) /  455.6d0,0.0d0,9.48d-10,0.3440d0,0.36d0/
  data(fitCoefficient(i,17, 7),i=1,5) /  529.3d0,1.0d0,3.69d-10,0.2730d0,0.26d0/
  data(fitCoefficient(i,17, 6),i=1,5) /  592.0d0,1.0d0,2.85d-10,0.3430d0,0.23d0/
  data(fitCoefficient(i,17, 5),i=1,5) /  656.7d0,0.0d0,4.81d-10,0.6580d0,0.27d0/
  data(fitCoefficient(i,17, 4),i=1,5) /  749.8d0,1.0d0,1.31d-10,0.6230d0,0.16d0/
  data(fitCoefficient(i,17, 3),i=1,5) /  809.4d0,1.0d0,6.13d-11,0.7360d0,0.16d0/
  data(fitCoefficient(i,17, 2),i=1,5) / 3658.4d0,0.0d0,1.90d-11,0.3790d0,0.36d0/
  data(fitCoefficient(i,17, 1),i=1,5) / 3946.3d0,1.0d0,5.14d-12,0.5530d0,0.18d0/
  data(fitCoefficient(i,18,18),i=1,5) /   15.8d0,1.0d0,5.99d-08,0.1360d0,0.26d0/
  data(fitCoefficient(i,18,17),i=1,5) /   27.6d0,1.0d0,6.07d-08,0.5440d0,0.21d0/
  data(fitCoefficient(i,18,16),i=1,5) /   40.9d0,1.0d0,3.43d-08,0.8340d0,0.17d0/
  data(fitCoefficient(i,18,15),i=1,5) /   52.3d0,0.0d0,3.00d-08,1.0300d0,0.25d0/
  data(fitCoefficient(i,18,14),i=1,5) /   75.0d0,1.0d0,8.73d-09,0.3660d0,0.31d0/
  data(fitCoefficient(i,18,13),i=1,5) /   91.0d0,1.0d0,5.78d-09,0.3140d0,0.34d0/
  data(fitCoefficient(i,18,12),i=1,5) /  124.3d0,1.0d0,2.98d-09,0.7030d0,0.16d0/
  data(fitCoefficient(i,18,11),i=1,5) /  143.5d0,1.0d0,7.25d-10,0.2070d0,0.15d0/
  data(fitCoefficient(i,18,10),i=1,5) /  422.4d0,1.0d0,1.40d-09,0.6960d0,0.13d0/
  data(fitCoefficient(i,18, 9),i=1,5) /  478.7d0,1.0d0,4.78d-10,0.1640d0,0.31d0/
  data(fitCoefficient(i,18, 8),i=1,5) /  539.0d0,0.0d0,8.02d-10,0.4390d0,0.32d0/
  data(fitCoefficient(i,18, 7),i=1,5) /  618.3d0,1.0d0,2.88d-10,0.2590d0,0.27d0/
  data(fitCoefficient(i,18, 6),i=1,5) /  686.1d0,1.0d0,2.32d-10,0.3620d0,0.22d0/
  data(fitCoefficient(i,18, 5),i=1,5) /  755.7d0,0.0d0,3.33d-10,0.4120d0,0.36d0/
  data(fitCoefficient(i,18, 4),i=1,5) /  854.8d0,1.0d0,1.27d-10,0.9100d0,0.13d0/
  data(fitCoefficient(i,18, 3),i=1,5) /  918.0d0,1.0d0,5.21d-11,0.7810d0,0.15d0/
  data(fitCoefficient(i,18, 2),i=1,5) / 4120.7d0,0.0d0,1.66d-11,0.4350d0,0.33d0/
  data(fitCoefficient(i,18, 1),i=1,5) / 4426.2d0,1.0d0,4.32d-12,0.5540d0,0.18d0/
  data(fitCoefficient(i,19,19),i=1,5) /    4.3d0,1.0d0,2.02d-07,0.2720d0,0.31d0/
  data(fitCoefficient(i,19,18),i=1,5) /   31.6d0,1.0d0,4.01d-08,0.3710d0,0.22d0/
  data(fitCoefficient(i,19,17),i=1,5) /   45.8d0,1.0d0,1.50d-08,0.4330d0,0.21d0/
  data(fitCoefficient(i,19,16),i=1,5) /   60.9d0,1.0d0,1.94d-08,0.8890d0,0.16d0/
  data(fitCoefficient(i,19,15),i=1,5) /   82.7d0,1.0d0,6.95d-09,0.4940d0,0.18d0/
  data(fitCoefficient(i,19,14),i=1,5) /   99.4d0,1.0d0,4.11d-09,0.5400d0,0.17d0/
  data(fitCoefficient(i,19,13),i=1,5) /  117.6d0,1.0d0,2.23d-09,0.5190d0,0.16d0/
  data(fitCoefficient(i,19,12),i=1,5) /  154.7d0,1.0d0,2.15d-09,0.8280d0,0.14d0/
  data(fitCoefficient(i,19,11),i=1,5) /  175.8d0,0.0d0,1.61d-09,0.6420d0,0.13d0/
  data(fitCoefficient(i,19,10),i=1,5) /  504.0d0,1.0d0,1.07d-09,0.6950d0,0.13d0/
  data(fitCoefficient(i,19, 9),i=1,5) /  564.7d0,1.0d0,3.78d-10,0.1730d0,0.30d0/
  data(fitCoefficient(i,19, 8),i=1,5) /  629.4d0,0.0d0,6.24d-10,0.4180d0,0.33d0/
  data(fitCoefficient(i,19, 7),i=1,5) /  714.6d0,1.0d0,2.29d-10,0.2450d0,0.28d0/
  data(fitCoefficient(i,19, 6),i=1,5) /  786.6d0,1.0d0,1.86d-10,0.3440d0,0.23d0/
  data(fitCoefficient(i,19, 5),i=1,5) /  861.1d0,0.0d0,2.69d-10,0.3960d0,0.37d0/
  data(fitCoefficient(i,19, 4),i=1,5) /  968.0d0,1.0d0,1.06d-10,0.9120d0,0.13d0/
  data(fitCoefficient(i,19, 3),i=1,5) / 1053.4d0,1.0d0,4.24d-11,0.7370d0,0.16d0/
  data(fitCoefficient(i,19, 2),i=1,5) / 4610.9d0,0.0d0,1.38d-11,0.4160d0,0.34d0/
  data(fitCoefficient(i,19, 1),i=1,5) / 4934.1d0,1.0d0,3.67d-12,0.5550d0,0.18d0/
  data(fitCoefficient(i,20,20),i=1,5) /    6.1d0,0.0d0,4.40d-07,0.8480d0,0.33d0/
  data(fitCoefficient(i,20,19),i=1,5) /   11.9d0,0.0d0,5.22d-08,0.1510d0,0.34d0/
  data(fitCoefficient(i,20,18),i=1,5) /   50.9d0,1.0d0,2.06d-08,0.4180d0,0.20d0/
  data(fitCoefficient(i,20,17),i=1,5) /   67.3d0,1.0d0,1.72d-08,0.6380d0,0.19d0/
  data(fitCoefficient(i,20,16),i=1,5) /   84.5d0,1.0d0,1.26d-08,1.0100d0,0.14d0/
  data(fitCoefficient(i,20,15),i=1,5) /  108.8d0,1.0d0,4.72d-09,0.5260d0,0.17d0/
  data(fitCoefficient(i,20,14),i=1,5) /  127.2d0,1.0d0,2.89d-09,0.5480d0,0.17d0/
  data(fitCoefficient(i,20,13),i=1,5) /  147.2d0,1.0d0,1.64d-09,0.5520d0,0.15d0/
  data(fitCoefficient(i,20,12),i=1,5) /  188.3d0,1.0d0,1.57d-09,0.7990d0,0.14d0/
  data(fitCoefficient(i,20,11),i=1,5) /  211.3d0,1.0d0,4.32d-10,0.2320d0,0.14d0/
  data(fitCoefficient(i,20,10),i=1,5) /  591.9d0,0.0d0,9.47d-10,0.3110d0,0.38d0/
  data(fitCoefficient(i,20, 9),i=1,5) /  657.2d0,1.0d0,2.98d-10,0.1630d0,0.31d0/
  data(fitCoefficient(i,20, 8),i=1,5) /  726.6d0,0.0d0,4.78d-10,0.3590d0,0.36d0/
  data(fitCoefficient(i,20, 7),i=1,5) /  817.6d0,1.0d0,1.86d-10,0.2440d0,0.28d0/
  data(fitCoefficient(i,20, 6),i=1,5) /  894.5d0,1.0d0,1.56d-10,0.3640d0,0.22d0/
  data(fitCoefficient(i,20, 5),i=1,5) /  974.0d0,0.0d0,2.16d-10,0.3570d0,0.39d0/
  data(fitCoefficient(i,20, 4),i=1,5) / 1087.0d0,1.0d0,7.70d-11,0.6550d0,0.15d0/
  data(fitCoefficient(i,20, 3),i=1,5) / 1157.0d0,1.0d0,3.58d-11,0.7360d0,0.16d0/
  data(fitCoefficient(i,20, 2),i=1,5) / 5128.9d0,0.0d0,1.28d-11,0.5200d0,0.30d0/
  data(fitCoefficient(i,20, 1),i=1,5) / 5469.9d0,1.0d0,3.08d-12,0.5280d0,0.19d0/
  data(fitCoefficient(i,21,21),i=1,5) /    6.6d0,1.0d0,3.16d-07,0.2040d0,0.28d0/
  data(fitCoefficient(i,21,20),i=1,5) /   12.8d0,1.0d0,8.61d-08,0.1810d0,0.25d0/
  data(fitCoefficient(i,21,19),i=1,5) /   24.8d0,1.0d0,5.08d-08,0.3570d0,0.24d0/
  data(fitCoefficient(i,21,18),i=1,5) /   73.5d0,1.0d0,1.00d-08,0.4530d0,0.15d0/
  data(fitCoefficient(i,21,17),i=1,5) /   91.9d0,1.0d0,6.76d-09,0.4600d0,0.15d0/
  data(fitCoefficient(i,21,16),i=1,5) /  110.7d0,1.0d0,5.27d-09,0.5610d0,0.17d0/
  data(fitCoefficient(i,21,15),i=1,5) /  138.0d0,1.0d0,3.40d-09,0.5600d0,0.16d0/
  data(fitCoefficient(i,21,14),i=1,5) /  158.1d0,1.0d0,2.18d-09,0.6120d0,0.15d0/
  data(fitCoefficient(i,21,13),i=1,5) /  180.0d0,1.0d0,1.26d-09,0.6100d0,0.14d0/
  data(fitCoefficient(i,21,12),i=1,5) /  225.1d0,1.0d0,1.24d-09,0.8520d0,0.13d0/
  data(fitCoefficient(i,21,11),i=1,5) /  249.8d0,1.0d0,3.62d-10,0.3490d0,0.05d0/
  data(fitCoefficient(i,21,10),i=1,5) /  687.4d0,1.0d0,5.52d-10,0.3750d0,0.28d0/
  data(fitCoefficient(i,21, 9),i=1,5) /  756.7d0,1.0d0,5.64d-10,0.8730d0,0.15d0/
  data(fitCoefficient(i,21, 8),i=1,5) /  830.8d0,1.0d0,4.50d-10,1.0500d0,0.13d0/
  data(fitCoefficient(i,21, 7),i=1,5) /  927.5d0,1.0d0,2.73d-10,0.8660d0,0.15d0/
  data(fitCoefficient(i,21, 6),i=1,5) / 1009.0d0,1.0d0,1.56d-10,0.7150d0,0.17d0/
  data(fitCoefficient(i,21, 5),i=1,5) / 1094.0d0,0.0d0,1.81d-10,1.1400d0,0.36d0/
  data(fitCoefficient(i,21, 4),i=1,5) / 1213.0d0,1.0d0,4.29d-11,0.7840d0,0.15d0/
  data(fitCoefficient(i,21, 3),i=1,5) / 1288.0d0,0.0d0,2.21d-11,0.0270d0,0.82d0/
  data(fitCoefficient(i,21, 2),i=1,5) / 5674.9d0,1.0d0,4.51d-12,0.9180d0,0.04d0/
  data(fitCoefficient(i,21, 1),i=1,5) / 6033.8d0,0.0d0,2.03d-12,0.0170d0,0.70d0/
  data(fitCoefficient(i,22,22),i=1,5) /    6.8d0,1.0d0,1.60d-07,0.3600d0,0.28d0/
  data(fitCoefficient(i,22,21),i=1,5) /   13.6d0,0.0d0,2.14d-07,0.8800d0,0.28d0/
  data(fitCoefficient(i,22,20),i=1,5) /   27.5d0,1.0d0,2.85d-08,0.2270d0,0.21d0/
  data(fitCoefficient(i,22,19),i=1,5) /   43.3d0,1.0d0,3.48d-08,0.3900d0,0.23d0/
  data(fitCoefficient(i,22,18),i=1,5) /   99.3d0,1.0d0,1.00d-08,0.5790d0,0.18d0/
  data(fitCoefficient(i,22,17),i=1,5) /  119.5d0,1.0d0,7.01d-09,0.6380d0,0.17d0/
  data(fitCoefficient(i,22,16),i=1,5) /  140.8d0,1.0d0,4.95d-09,0.7170d0,0.16d0/
  data(fitCoefficient(i,22,15),i=1,5) /  170.4d0,1.0d0,2.99d-09,0.6930d0,0.17d0/
  data(fitCoefficient(i,22,14),i=1,5) /  192.1d0,1.0d0,2.10d-09,0.7220d0,0.16d0/
  data(fitCoefficient(i,22,13),i=1,5) /  215.9d0,1.0d0,1.62d-09,0.7650d0,0.14d0/
  data(fitCoefficient(i,22,12),i=1,5) /  265.0d0,1.0d0,1.11d-09,0.8850d0,0.12d0/
  data(fitCoefficient(i,22,11),i=1,5) /  291.5d0,0.0d0,9.09d-10,0.9720d0,0.06d0/
  data(fitCoefficient(i,22,10),i=1,5) /  787.8d0,1.0d0,4.41d-10,0.3590d0,0.29d0/
  data(fitCoefficient(i,22, 9),i=1,5) /  863.1d0,1.0d0,4.39d-10,0.7810d0,0.17d0/
  data(fitCoefficient(i,22, 8),i=1,5) /  941.9d0,1.0d0,3.73d-10,1.0500d0,0.13d0/
  data(fitCoefficient(i,22, 7),i=1,5) / 1044.0d0,1.0d0,2.28d-10,0.8580d0,0.15d0/
  data(fitCoefficient(i,22, 6),i=1,5) / 1131.0d0,1.0d0,1.34d-10,0.7570d0,0.16d0/
  data(fitCoefficient(i,22, 5),i=1,5) / 1221.0d0,0.0d0,1.55d-10,1.1500d0,0.36d0/
  data(fitCoefficient(i,22, 4),i=1,5) / 1346.0d0,1.0d0,3.80d-11,0.8350d0,0.14d0/
  data(fitCoefficient(i,22, 3),i=1,5) / 1426.0d0,0.0d0,1.89d-11,0.0280d0,0.82d0/
  data(fitCoefficient(i,22, 2),i=1,5) / 6249.1d0,1.0d0,4.01d-12,0.9680d0,0.03d0/
  data(fitCoefficient(i,22, 1),i=1,5) / 6625.0d0,1.0d0,1.62d-12,0.6570d0,0.14d0/
  data(fitCoefficient(i,23,23),i=1,5) /    6.7d0,0.0d0,8.82d-07,0.3590d0,0.32d0/
  data(fitCoefficient(i,23,22),i=1,5) /   14.7d0,0.0d0,3.11d-07,0.4320d0,0.29d0/
  data(fitCoefficient(i,23,21),i=1,5) /   29.3d0,1.0d0,3.50d-08,0.2470d0,0.25d0/
  data(fitCoefficient(i,23,20),i=1,5) /   46.7d0,0.0d0,5.32d-08,1.1100d0,0.16d0/
  data(fitCoefficient(i,23,19),i=1,5) /   65.3d0,1.0d0,8.98d-09,0.1400d0,0.37d0/
  data(fitCoefficient(i,23,18),i=1,5) /  128.1d0,1.0d0,5.87d-09,0.5170d0,0.17d0/
  data(fitCoefficient(i,23,17),i=1,5) /  150.6d0,1.0d0,5.11d-09,0.6790d0,0.16d0/
  data(fitCoefficient(i,23,16),i=1,5) /  173.4d0,1.0d0,3.71d-09,0.7610d0,0.15d0/
  data(fitCoefficient(i,23,15),i=1,5) /  205.8d0,1.0d0,2.24d-09,0.7110d0,0.17d0/
  data(fitCoefficient(i,23,14),i=1,5) /  230.5d0,1.0d0,1.65d-09,0.7640d0,0.15d0/
  data(fitCoefficient(i,23,13),i=1,5) /  256.0d0,1.0d0,1.26d-09,0.7620d0,0.14d0/
  data(fitCoefficient(i,23,12),i=1,5) /  308.0d0,1.0d0,8.86d-10,0.8860d0,0.12d0/
  data(fitCoefficient(i,23,11),i=1,5) /  336.3d0,0.0d0,3.89d-10,0.1420d0,0.39d0/
  data(fitCoefficient(i,23,10),i=1,5) /  896.0d0,1.0d0,3.80d-10,0.4090d0,0.27d0/
  data(fitCoefficient(i,23, 9),i=1,5) /  976.0d0,0.0d0,4.84d-10,0.1730d0,0.57d0/
  data(fitCoefficient(i,23, 8),i=1,5) / 1060.0d0,1.0d0,2.49d-10,0.6500d0,0.14d0/
  data(fitCoefficient(i,23, 7),i=1,5) / 1168.0d0,0.0d0,5.91d-10,1.6100d0,0.18d0/
  data(fitCoefficient(i,23, 6),i=1,5) / 1260.0d0,0.0d0,5.02d-10,2.1200d0,0.15d0/
  data(fitCoefficient(i,23, 5),i=1,5) / 1355.0d0,1.0d0,5.38d-11,0.1370d0,0.40d0/
  data(fitCoefficient(i,23, 4),i=1,5) / 1486.0d0,1.0d0,5.56d-11,0.7080d0,0.10d0/
  data(fitCoefficient(i,23, 3),i=1,5) / 1571.0d0,0.0d0,2.84d-11,0.0240d0,0.79d0/
  data(fitCoefficient(i,23, 2),i=1,5) / 6851.3d0,0.0d0,2.54d-11,2.9200d0,0.09d0/
  data(fitCoefficient(i,23, 1),i=1,5) / 7246.1d0,0.0d0,1.32d-11,3.5100d0,0.07d0/
  data(fitCoefficient(i,24,24),i=1,5) /    6.8d0,1.0d0,1.03d-07,0.2170d0,0.27d0/
  data(fitCoefficient(i,24,23),i=1,5) /   16.5d0,0.0d0,2.45d-07,0.3810d0,0.32d0/
  data(fitCoefficient(i,24,22),i=1,5) /   31.0d0,0.0d0,1.09d-07,0.5180d0,0.27d0/
  data(fitCoefficient(i,24,21),i=1,5) /   49.1d0,1.0d0,1.52d-08,0.1820d0,0.30d0/
  data(fitCoefficient(i,24,20),i=1,5) /   69.5d0,0.0d0,3.25d-08,1.3600d0,0.13d0/
  data(fitCoefficient(i,24,19),i=1,5) /   90.6d0,1.0d0,5.50d-09,0.1430d0,0.37d0/
  data(fitCoefficient(i,24,18),i=1,5) /  160.2d0,1.0d0,5.13d-09,0.6570d0,0.16d0/
  data(fitCoefficient(i,24,17),i=1,5) /  184.7d0,1.0d0,3.85d-09,0.7220d0,0.15d0/
  data(fitCoefficient(i,24,16),i=1,5) /  209.3d0,1.0d0,2.81d-09,0.7590d0,0.15d0/
  data(fitCoefficient(i,24,15),i=1,5) /  244.4d0,1.0d0,1.76d-09,0.7320d0,0.16d0/
  data(fitCoefficient(i,24,14),i=1,5) /  271.0d0,1.0d0,1.30d-09,0.7640d0,0.15d0/
  data(fitCoefficient(i,24,13),i=1,5) /  298.0d0,1.0d0,1.02d-09,0.8100d0,0.13d0/
  data(fitCoefficient(i,24,12),i=1,5) /  354.8d0,1.0d0,7.19d-10,0.8870d0,0.12d0/
  data(fitCoefficient(i,24,11),i=1,5) /  384.2d0,1.0d0,1.61d-10,0.1500d0,0.22d0/
  data(fitCoefficient(i,24,10),i=1,5) / 1011.0d0,1.0d0,4.64d-10,0.9710d0,0.12d0/
  data(fitCoefficient(i,24, 9),i=1,5) / 1097.0d0,1.0d0,3.31d-10,0.9240d0,0.14d0/
  data(fitCoefficient(i,24, 8),i=1,5) / 1185.0d0,1.0d0,2.49d-10,0.9310d0,0.15d0/
  data(fitCoefficient(i,24, 7),i=1,5) / 1299.0d0,1.0d0,1.68d-10,0.9100d0,0.14d0/
  data(fitCoefficient(i,24, 6),i=1,5) / 1396.0d0,1.0d0,1.01d-10,0.8050d0,0.15d0/
  data(fitCoefficient(i,24, 5),i=1,5) / 1496.0d0,0.0d0,1.17d-10,1.2100d0,0.35d0/
  data(fitCoefficient(i,24, 4),i=1,5) / 1634.0d0,1.0d0,2.91d-11,0.8840d0,0.13d0/
  data(fitCoefficient(i,24, 3),i=1,5) / 1721.0d0,0.0d0,1.45d-11,0.0350d0,0.80d0/
  data(fitCoefficient(i,24, 2),i=1,5) / 7482.0d0,1.0d0,3.07d-12,0.9670d0,0.03d0/
  data(fitCoefficient(i,24, 1),i=1,5) / 7894.8d0,1.0d0,1.46d-12,0.1830d0,0.39d0/
  data(fitCoefficient(i,25,25),i=1,5) /    7.4d0,1.0d0,8.56d-08,0.1320d0,0.26d0/
  data(fitCoefficient(i,25,24),i=1,5) /   15.6d0,0.0d0,1.18d-07,0.3590d0,0.19d0/
  data(fitCoefficient(i,25,23),i=1,5) /   33.7d0,0.0d0,8.54d-08,0.3970d0,0.32d0/
  data(fitCoefficient(i,25,22),i=1,5) /   51.2d0,1.0d0,1.80d-08,0.2720d0,0.18d0/
  data(fitCoefficient(i,25,21),i=1,5) /   72.4d0,1.0d0,8.22d-09,0.1610d0,0.32d0/
  data(fitCoefficient(i,25,20),i=1,5) /   95.0d0,0.0d0,2.15d-08,1.5400d0,0.11d0/
  data(fitCoefficient(i,25,19),i=1,5) /  119.3d0,1.0d0,3.65d-09,0.1470d0,0.37d0/
  data(fitCoefficient(i,25,18),i=1,5) /  194.5d0,1.0d0,3.91d-09,0.6990d0,0.15d0/
  data(fitCoefficient(i,25,17),i=1,5) /  221.8d0,1.0d0,2.92d-09,0.7190d0,0.15d0/
  data(fitCoefficient(i,25,16),i=1,5) /  248.3d0,1.0d0,2.23d-09,0.8060d0,0.14d0/
  data(fitCoefficient(i,25,15),i=1,5) /  286.0d0,1.0d0,1.39d-09,0.7350d0,0.16d0/
  data(fitCoefficient(i,25,14),i=1,5) /  314.4d0,1.0d0,1.04d-09,0.7610d0,0.15d0/
  data(fitCoefficient(i,25,13),i=1,5) /  343.6d0,1.0d0,8.28d-10,0.8090d0,0.13d0/
  data(fitCoefficient(i,25,12),i=1,5) /  403.0d0,1.0d0,5.60d-10,0.7870d0,0.14d0/
  data(fitCoefficient(i,25,11),i=1,5) /  435.2d0,1.0d0,1.52d-10,0.2990d0,0.08d0/
  data(fitCoefficient(i,25,10),i=1,5) / 1133.0d0,1.0d0,4.03d-10,1.0400d0,0.11d0/
  data(fitCoefficient(i,25, 9),i=1,5) / 1244.0d0,1.0d0,2.74d-10,0.9230d0,0.14d0/
  data(fitCoefficient(i,25, 8),i=1,5) / 1317.0d0,1.0d0,2.18d-10,0.9900d0,0.14d0/
  data(fitCoefficient(i,25, 7),i=1,5) / 1437.0d0,1.0d0,1.49d-10,0.9680d0,0.13d0/
  data(fitCoefficient(i,25, 6),i=1,5) / 1539.0d0,1.0d0,8.70d-11,0.8020d0,0.15d0/
  data(fitCoefficient(i,25, 5),i=1,5) / 1644.0d0,0.0d0,1.02d-10,1.2200d0,0.35d0/
  data(fitCoefficient(i,25, 4),i=1,5) / 1788.0d0,1.0d0,2.54d-11,0.8830d0,0.13d0/
  data(fitCoefficient(i,25, 3),i=1,5) / 1880.0d0,0.0d0,1.28d-11,0.0330d0,0.81d0/
  data(fitCoefficient(i,25, 2),i=1,5) / 8141.0d0,1.0d0,2.77d-12,1.0100d0,0.02d0/
  data(fitCoefficient(i,25, 1),i=1,5) / 8571.9d0,1.0d0,1.32d-12,0.2190d0,0.37d0/
  data(fitCoefficient(i,26,26),i=1,5) /    7.9d0,0.0d0,2.52d-07,0.7010d0,0.25d0/
  data(fitCoefficient(i,26,25),i=1,5) /   16.2d0,1.0d0,2.21d-08,0.0330d0,0.45d0/
  data(fitCoefficient(i,26,24),i=1,5) /   30.6d0,0.0d0,4.10d-08,0.3660d0,0.17d0/
  data(fitCoefficient(i,26,23),i=1,5) /   54.8d0,0.0d0,3.53d-08,0.2430d0,0.39d0/
  data(fitCoefficient(i,26,22),i=1,5) /   75.0d0,1.0d0,1.04d-08,0.2850d0,0.17d0/
  data(fitCoefficient(i,26,21),i=1,5) /   99.0d0,1.0d0,1.23d-08,0.4110d0,0.21d0/
  data(fitCoefficient(i,26,20),i=1,5) /  125.0d0,1.0d0,9.47d-09,0.4580d0,0.21d0/
  data(fitCoefficient(i,26,19),i=1,5) /  151.1d0,1.0d0,4.71d-09,0.2800d0,0.28d0/
  data(fitCoefficient(i,26,18),i=1,5) /  233.6d0,1.0d0,3.02d-09,0.6970d0,0.15d0/
  data(fitCoefficient(i,26,17),i=1,5) /  262.1d0,1.0d0,2.34d-09,0.7640d0,0.14d0/
  data(fitCoefficient(i,26,16),i=1,5) /  290.0d0,1.0d0,1.76d-09,0.8050d0,0.14d0/
  data(fitCoefficient(i,26,15),i=1,5) /  331.0d0,1.0d0,1.14d-09,0.7730d0,0.15d0/
  data(fitCoefficient(i,26,14),i=1,5) /  361.0d0,1.0d0,8.66d-10,0.8050d0,0.14d0/
  data(fitCoefficient(i,26,13),i=1,5) /  392.0d0,1.0d0,6.61d-10,0.7620d0,0.14d0/
  data(fitCoefficient(i,26,12),i=1,5) /  457.0d0,1.0d0,4.41d-10,0.6980d0,0.16d0/
  data(fitCoefficient(i,26,11),i=1,5) /  489.3d0,1.0d0,1.18d-10,0.2110d0,0.15d0/
  data(fitCoefficient(i,26,10),i=1,5) / 1262.0d0,1.0d0,3.61d-10,1.1600d0,0.09d0/
  data(fitCoefficient(i,26, 9),i=1,5) / 1360.0d0,1.0d0,2.45d-10,0.9780d0,0.13d0/
  data(fitCoefficient(i,26, 8),i=1,5) / 1470.0d0,1.0d0,1.87d-10,0.9880d0,0.14d0/
  data(fitCoefficient(i,26, 7),i=1,5) / 1582.0d0,1.0d0,1.33d-10,1.0300d0,0.12d0/
  data(fitCoefficient(i,26, 6),i=1,5) / 1690.0d0,1.0d0,7.84d-11,0.8480d0,0.14d0/
  data(fitCoefficient(i,26, 5),i=1,5) / 1800.0d0,0.0d0,8.90d-11,1.2000d0,0.35d0/
  data(fitCoefficient(i,26, 4),i=1,5) / 1960.0d0,1.0d0,2.29d-11,0.9360d0,0.12d0/
  data(fitCoefficient(i,26, 3),i=1,5) / 2046.0d0,0.0d0,1.12d-11,0.0340d0,0.81d0/
  data(fitCoefficient(i,26, 2),i=1,5) / 8828.0d0,1.0d0,2.46d-12,1.0200d0,0.02d0/
  data(fitCoefficient(i,26, 1),i=1,5) / 9277.7d0,1.0d0,9.79d-13,0.6640d0,0.14d0/
  data(fitCoefficient(i,27,27),i=1,5) /    7.9d0,1.0d0,8.89d-08,0.1270d0,0.24d0/
  data(fitCoefficient(i,27,26),i=1,5) /   17.1d0,1.0d0,5.65d-08,0.1940d0,0.23d0/
  data(fitCoefficient(i,27,25),i=1,5) /   33.5d0,1.0d0,3.06d-08,0.2010d0,0.22d0/
  data(fitCoefficient(i,27,24),i=1,5) /   51.3d0,0.0d0,2.27d-08,0.5740d0,0.10d0/
  data(fitCoefficient(i,27,23),i=1,5) /   79.5d0,0.0d0,1.93d-08,0.1950d0,0.42d0/
  data(fitCoefficient(i,27,22),i=1,5) /  102.0d0,0.0d0,1.27d-08,0.1260d0,0.47d0/
  data(fitCoefficient(i,27,21),i=1,5) /  129.0d0,1.0d0,3.58d-09,0.1940d0,0.29d0/
  data(fitCoefficient(i,27,20),i=1,5) /  158.0d0,0.0d0,1.17d-08,1.9800d0,0.07d0/
  data(fitCoefficient(i,27,19),i=1,5) /  186.1d0,1.0d0,1.78d-09,0.1120d0,0.42d0/
  data(fitCoefficient(i,27,18),i=1,5) /  275.0d0,1.0d0,2.41d-09,0.7390d0,0.14d0/
  data(fitCoefficient(i,27,17),i=1,5) /  305.0d0,1.0d0,1.86d-09,0.7620d0,0.14d0/
  data(fitCoefficient(i,27,16),i=1,5) /  336.0d0,1.0d0,1.41d-09,0.8040d0,0.14d0/
  data(fitCoefficient(i,27,15),i=1,5) /  379.0d0,1.0d0,9.54d-10,0.8130d0,0.14d0/
  data(fitCoefficient(i,27,14),i=1,5) /  411.0d0,1.0d0,7.12d-10,0.8030d0,0.14d0/
  data(fitCoefficient(i,27,13),i=1,5) /  444.0d0,1.0d0,5.34d-10,0.7180d0,0.15d0/
  data(fitCoefficient(i,27,12),i=1,5) /  512.0d0,1.0d0,3.62d-10,0.6580d0,0.17d0/
  data(fitCoefficient(i,27,11),i=1,5) /  546.6d0,1.0d0,1.05d-10,0.2520d0,0.12d0/
  data(fitCoefficient(i,27,10),i=1,5) / 1397.0d0,1.0d0,3.10d-10,1.1700d0,0.09d0/
  data(fitCoefficient(i,27, 9),i=1,5) / 1486.0d0,1.0d0,1.56d-10,0.5720d0,0.15d0/
  data(fitCoefficient(i,27, 8),i=1,5) / 1603.0d0,1.0d0,1.32d-10,0.6820d0,0.13d0/
  data(fitCoefficient(i,27, 7),i=1,5) / 1735.0d0,1.0d0,9.08d-11,0.5110d0,0.17d0/
  data(fitCoefficient(i,27, 6),i=1,5) / 1846.0d0,0.0d0,3.45d-10,2.8400d0,0.11d0/
  data(fitCoefficient(i,27, 5),i=1,5) / 1962.0d0,1.0d0,5.01d-11,0.7140d0,0.11d0/
  data(fitCoefficient(i,27, 4),i=1,5) / 2119.0d0,1.0d0,1.92d-11,0.1170d0,0.42d0/
  data(fitCoefficient(i,27, 3),i=1,5) / 2219.0d0,1.0d0,1.95d-11,1.2000d0,0.00d0/
  data(fitCoefficient(i,27, 2),i=1,5) / 9544.0d0,0.0d0,1.68d-11,3.5200d0,0.05d0/
  data(fitCoefficient(i,27, 1),i=1,5) /10012.1d0,1.0d0,1.45d-12,0.6350d0,0.15d0/
  data(fitCoefficient(i,28,28),i=1,5) /    7.6d0,0.0d0,1.65d-07,0.4520d0,0.28d0/
  data(fitCoefficient(i,28,27),i=1,5) /   18.2d0,0.0d0,8.42d-08,0.6190d0,0.16d0/
  data(fitCoefficient(i,28,26),i=1,5) /   35.3d0,1.0d0,1.89d-08,0.2200d0,0.21d0/
  data(fitCoefficient(i,28,25),i=1,5) /   54.9d0,1.0d0,1.48d-08,0.2160d0,0.21d0/
  data(fitCoefficient(i,28,24),i=1,5) /   76.0d0,0.0d0,1.13d-08,0.5180d0,0.09d0/
  data(fitCoefficient(i,28,23),i=1,5) /  108.0d0,0.0d0,1.16d-08,0.1850d0,0.44d0/
  data(fitCoefficient(i,28,22),i=1,5) /  133.0d0,0.0d0,8.68d-09,0.1380d0,0.46d0/
  data(fitCoefficient(i,28,21),i=1,5) /  162.0d0,1.0d0,2.45d-09,0.1630d0,0.32d0/
  data(fitCoefficient(i,28,20),i=1,5) /  193.0d0,0.0d0,9.24d-09,2.2500d0,0.05d0/
  data(fitCoefficient(i,28,19),i=1,5) /  225.0d0,0.0d0,2.41d-09,0.0270d0,0.79d0/
  data(fitCoefficient(i,28,18),i=1,5) /  321.0d0,1.0d0,1.92d-09,0.7380d0,0.14d0/
  data(fitCoefficient(i,28,17),i=1,5) /  352.0d0,1.0d0,1.50d-09,0.7610d0,0.14d0/
  data(fitCoefficient(i,28,16),i=1,5) /  384.0d0,1.0d0,1.16d-09,0.8030d0,0.14d0/
  data(fitCoefficient(i,28,15),i=1,5) /  430.0d0,1.0d0,8.08d-10,0.8560d0,0.13d0/
  data(fitCoefficient(i,28,14),i=1,5) /  464.0d0,1.0d0,6.09d-10,0.8500d0,0.13d0/
  data(fitCoefficient(i,28,13),i=1,5) /  499.0d0,1.0d0,4.48d-10,0.7180d0,0.15d0/
  data(fitCoefficient(i,28,12),i=1,5) /  571.3d0,1.0d0,3.00d-10,0.6220d0,0.18d0/
  data(fitCoefficient(i,28,11),i=1,5) /  607.0d0,1.0d0,7.90d-11,0.1600d0,0.19d0/
  data(fitCoefficient(i,28,10),i=1,5) / 1541.0d0,1.0d0,2.78d-10,1.2500d0,0.08d0/
  data(fitCoefficient(i,28, 9),i=1,5) / 1648.0d0,1.0d0,1.92d-10,1.0400d0,0.12d0/
  data(fitCoefficient(i,28, 8),i=1,5) / 1756.0d0,1.0d0,1.51d-10,1.1100d0,0.12d0/
  data(fitCoefficient(i,28, 7),i=1,5) / 1894.0d0,1.0d0,1.05d-10,1.0900d0,0.11d0/
  data(fitCoefficient(i,28, 6),i=1,5) / 2011.0d0,1.0d0,6.04d-11,0.8490d0,0.14d0/
  data(fitCoefficient(i,28, 5),i=1,5) / 2131.0d0,0.0d0,6.91d-11,1.2100d0,0.35d0/
  data(fitCoefficient(i,28, 4),i=1,5) / 2295.0d0,1.0d0,1.84d-11,0.9910d0,0.11d0/
  data(fitCoefficient(i,28, 3),i=1,5) / 2399.0d0,0.0d0,9.03d-12,0.0420d0,0.79d0/
  data(fitCoefficient(i,28, 2),i=1,5) /10290.0d0,1.0d0,2.61d-12,0.5680d0,0.16d0/
  data(fitCoefficient(i,28, 1),i=1,5) /10775.3d0,1.0d0,1.39d-12,0.7290d0,0.13d0/

contains

  function verner1996ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{atomicIonizationRateCollisionalVerner1996} atomic collisional ionization class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(atomicIonizationRateCollisionalVerner1996)                :: self
    type(inputParameters                          ), intent(inout) :: parameters

    self=atomicIonizationRateCollisionalVerner1996()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function verner1996ConstructorParameters

  double precision function verner1996Rate(self,atomicNumber,ionizationState,temperature)
    !!{
    Computes the rate coefficient of direct collisional ionization (in units of cm$^3$ s$^{-1}$) at the specified {\normalfont \ttfamily
    temperature} for all ions of atoms with $Z<28$ by use of the fits from
    \citeauthor{voronov_practical_1997}~(\citeyear{voronov_practical_1997}; Version 2, March 24, 1997). Based on the
    \href{https://web.archive.org/web/20220313133801/https://www.pa.uky.edu/~verner/dima/col/cfit.f}{code} originally written by Dima Verner. The ionization state passed to
    this function should be that of the atom/ion prior to ionization.
    !!}
    use :: Numerical_Constants_Physical, only : boltzmannsConstant
    use :: Numerical_Constants_Units   , only : electronVolt
    implicit none
    class           (atomicIonizationRateCollisionalVerner1996), intent(inout) :: self
    integer                                                    , intent(in   ) :: atomicNumber  , ionizationState
    double precision                                           , intent(in   ) :: temperature
    integer                                                                    :: electronNumber
    double precision                                                           :: energyScaled  , temperatureScaled
    !$GLC attributes unused :: self

    ! Set a default rate coefficient of zero.
    verner1996Rate=0.0d0
    ! Compute number of electrons.
    electronNumber=atomicNumber-ionizationState+1
    ! Return immediately for unphysical conditions or out of range cases.
    if (temperature    <= 0.0d0                               ) return
    if (atomicNumber   <  1 .or. atomicNumber   > 28          ) return
    if (electronNumber <  1 .or. electronNumber > atomicNumber) return
    ! Compute scaled energy.
    temperatureScaled=+temperature        &
         &            *boltzmannsConstant &
         &            /electronVolt
    energyScaled     =+fitCoefficient(1,atomicNumber,electronNumber) &
         &            /temperatureScaled
    ! Return if the scaled energy is out of range.
    if (energyScaled > 80.0d0) return
    ! Compute the rate coefficient
    verner1996Rate=+                      fitCoefficient(3,atomicNumber,electronNumber) &
         &         *(                                                                   &
         &           +1.0d0                                                             &
         &           +                    fitCoefficient(2,atomicNumber,electronNumber) &
         &           *sqrt(+energyScaled)                                               &
         &          )                                                                   &
         &         /(                                                                   &
         &           +                    fitCoefficient(4,atomicNumber,electronNumber) &
         &           +      energyScaled                                                &
         &          )                                                                   &
         &         *        energyScaled**fitCoefficient(5,atomicNumber,electronNumber) &
         &         *exp   (-energyScaled)
    return
  end function verner1996Rate
