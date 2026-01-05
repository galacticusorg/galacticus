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
  An implementation of atomic radiative recombination rates based on the
  \href{https://web.archive.org/web/20220313133758/https://www.pa.uky.edu/~verner/dima/rec/rrfit.f}{code} originally written by Dima Verner.
  !!}

  !![
  <atomicRecombinationRateRadiative name="atomicRecombinationRateRadiativeVerner1996">
   <description>Atomic radiative recombination rates are computed based on the \href{https://web.archive.org/web/20220313133758/https://www.pa.uky.edu/~verner/dima/rec/rrfit.f}{code} originally written by Dima Verner.</description>
  </atomicRecombinationRateRadiative>
  !!]
  type, extends(atomicRecombinationRateRadiativeClass) :: atomicRecombinationRateRadiativeVerner1996
     !!{
     A radiative recombination rate class based on the \href{https://web.archive.org/web/20220313133758/https://www.pa.uky.edu/~verner/dima/rec/rrfit.f}{code} originally written by Dima Verner.
     !!}
     private
   contains
     procedure :: rate => verner1996Rate
  end type atomicRecombinationRateRadiativeVerner1996

  interface atomicRecombinationRateRadiativeVerner1996
     !!{
     Constructors for the \refClass{atomicRecombinationRateRadiativeVerner1996} atomic radiative recombination class.
     !!}
     module procedure verner1996ConstructorParameters
  end interface atomicRecombinationRateRadiativeVerner1996

  ! Arrays to hold coefficients of fitting functions.
  double precision :: coefficientsIron  (3,10   ), fitCoefficients(2,30,30), &
       &              fitCoefficientsNew(4,30,30)
  integer          :: i

  ! Set the fitting coefficients.
  data(fitCoefficients   (i, 4, 4),i=1, 2) /4.500d-13,0.6480d0/
  data(fitCoefficients   (i, 5, 4),i=1, 2) /2.680d-12,0.7250d0/
  data(fitCoefficients   (i, 6, 4),i=1, 2) /4.900d-12,0.8030d0/
  data(fitCoefficients   (i, 7, 4),i=1, 2) /9.400d-12,0.7650d0/
  data(fitCoefficients   (i, 8, 4),i=1, 2) /1.590d-11,0.7590d0/
  data(fitCoefficients   (i, 9, 4),i=1, 2) /1.919d-11,0.7849d0/
  data(fitCoefficients   (i,10, 4),i=1, 2) /2.800d-11,0.7710d0/
  data(fitCoefficients   (i,11, 4),i=1, 2) /4.030d-11,0.7873d0/
  data(fitCoefficients   (i,12, 4),i=1, 2) /6.830d-11,0.7650d0/
  data(fitCoefficients   (i,13, 4),i=1, 2) /8.343d-11,0.8291d0/
  data(fitCoefficients   (i,14, 4),i=1, 2) /1.430d-10,0.8230d0/
  data(fitCoefficients   (i,15, 4),i=1, 2) /1.460d-10,0.8391d0/
  data(fitCoefficients   (i,16, 4),i=1, 2) /2.000d-10,0.8060d0/
  data(fitCoefficients   (i,17, 4),i=1, 2) /2.091d-10,0.8691d0/
  data(fitCoefficients   (i,18, 4),i=1, 2) /3.070d-10,0.8190d0/
  data(fitCoefficients   (i,19, 4),i=1, 2) /3.210d-10,0.8934d0/
  data(fitCoefficients   (i,20, 4),i=1, 2) /4.610d-10,0.8330d0/
  data(fitCoefficients   (i,21, 4),i=1, 2) /4.870d-10,0.8060d0/
  data(fitCoefficients   (i,22, 4),i=1, 2) /5.139d-10,0.7781d0/
  data(fitCoefficients   (i,23, 4),i=1, 2) /5.850d-10,0.7570d0/
  data(fitCoefficients   (i,24, 4),i=1, 2) /6.556d-10,0.7359d0/
  data(fitCoefficients   (i,25, 4),i=1, 2) /7.238d-10,0.7242d0/
  data(fitCoefficients   (i,27, 4),i=1, 2) /8.404d-10,0.7167d0/
  data(fitCoefficients   (i,28, 4),i=1, 2) /1.360d-09,0.8420d0/
  data(fitCoefficients   (i,29, 4),i=1, 2) /1.880d-09,0.8420d0/
  data(fitCoefficients   (i,30, 4),i=1, 2) /2.400d-09,0.8420d0/
  data(fitCoefficients   (i, 5, 5),i=1, 2) /4.600d-13,0.6360d0/
  data(fitCoefficients   (i, 6, 5),i=1, 2) /2.300d-12,0.6450d0/
  data(fitCoefficients   (i, 7, 5),i=1, 2) /5.000d-12,0.6760d0/
  data(fitCoefficients   (i, 8, 5),i=1, 2) /9.600d-12,0.6700d0/
  data(fitCoefficients   (i, 9, 5),i=1, 2) /1.558d-11,0.6816d0/
  data(fitCoefficients   (i,10, 5),i=1, 2) /2.300d-11,0.7040d0/
  data(fitCoefficients   (i,11, 5),i=1, 2) /3.253d-11,0.7075d0/
  data(fitCoefficients   (i,12, 5),i=1, 2) /4.600d-11,0.7110d0/
  data(fitCoefficients   (i,13, 5),i=1, 2) /5.951d-11,0.7125d0/
  data(fitCoefficients   (i,14, 5),i=1, 2) /7.700d-11,0.7140d0/
  data(fitCoefficients   (i,15, 5),i=1, 2) /1.042d-10,0.7330d0/
  data(fitCoefficients   (i,16, 5),i=1, 2) /1.400d-10,0.7550d0/
  data(fitCoefficients   (i,17, 5),i=1, 2) /1.760d-10,0.7682d0/
  data(fitCoefficients   (i,18, 5),i=1, 2) /2.140d-10,0.7740d0/
  data(fitCoefficients   (i,19, 5),i=1, 2) /2.629d-10,0.7772d0/
  data(fitCoefficients   (i,20, 5),i=1, 2) /3.240d-10,0.7800d0/
  data(fitCoefficients   (i,21, 5),i=1, 2) /3.970d-10,0.7840d0/
  data(fitCoefficients   (i,22, 5),i=1, 2) /4.700d-10,0.7873d0/
  data(fitCoefficients   (i,23, 5),i=1, 2) /5.500d-10,0.7910d0/
  data(fitCoefficients   (i,24, 5),i=1, 2) /6.301d-10,0.7952d0/
  data(fitCoefficients   (i,25, 5),i=1, 2) /7.058d-10,0.7982d0/
  data(fitCoefficients   (i,27, 5),i=1, 2) /8.252d-10,0.8004d0/
  data(fitCoefficients   (i,28, 5),i=1, 2) /8.710d-10,0.8000d0/
  data(fitCoefficients   (i,29, 5),i=1, 2) /9.170d-10,0.8000d0/
  data(fitCoefficients   (i,30, 5),i=1, 2) /9.630d-10,0.7990d0/
  data(fitCoefficients   (i, 6, 6),i=1, 2) /4.700d-13,0.6240d0/
  data(fitCoefficients   (i, 7, 6),i=1, 2) /2.200d-12,0.6390d0/
  data(fitCoefficients   (i, 8, 6),i=1, 2) /5.100d-12,0.6600d0/
  data(fitCoefficients   (i, 9, 6),i=1, 2) /9.171d-12,0.6757d0/
  data(fitCoefficients   (i,10, 6),i=1, 2) /1.500d-11,0.6840d0/
  data(fitCoefficients   (i,11, 6),i=1, 2) /2.191d-11,0.6875d0/
  data(fitCoefficients   (i,12, 6),i=1, 2) /3.200d-11,0.6910d0/
  data(fitCoefficients   (i,13, 6),i=1, 2) /4.308d-11,0.6970d0/
  data(fitCoefficients   (i,14, 6),i=1, 2) /5.800d-11,0.7030d0/
  data(fitCoefficients   (i,15, 6),i=1, 2) /7.316d-11,0.7027d0/
  data(fitCoefficients   (i,16, 6),i=1, 2) /9.200d-11,0.7140d0/
  data(fitCoefficients   (i,17, 6),i=1, 2) /1.198d-10,0.7508d0/
  data(fitCoefficients   (i,18, 6),i=1, 2) /1.580d-10,0.7900d0/
  data(fitCoefficients   (i,19, 6),i=1, 2) /2.048d-10,0.8032d0/
  data(fitCoefficients   (i,20, 6),i=1, 2) /2.600d-10,0.8000d0/
  data(fitCoefficients   (i,21, 6),i=1, 2) /3.280d-10,0.7990d0/
  data(fitCoefficients   (i,22, 6),i=1, 2) /3.966d-10,0.7973d0/
  data(fitCoefficients   (i,23, 6),i=1, 2) /4.760d-10,0.8000d0/
  data(fitCoefficients   (i,24, 6),i=1, 2) /5.547d-10,0.8027d0/
  data(fitCoefficients   (i,25, 6),i=1, 2) /6.313d-10,0.8058d0/
  data(fitCoefficients   (i,27, 6),i=1, 2) /7.503d-10,0.8085d0/
  data(fitCoefficients   (i,28, 6),i=1, 2) /7.940d-10,0.8080d0/
  data(fitCoefficients   (i,29, 6),i=1, 2) /8.380d-10,0.8080d0/
  data(fitCoefficients   (i,30, 6),i=1, 2) /8.810d-10,0.8070d0/
  data(fitCoefficients   (i, 7, 7),i=1, 2) /4.100d-13,0.6080d0/
  data(fitCoefficients   (i, 8, 7),i=1, 2) /2.000d-12,0.6460d0/
  data(fitCoefficients   (i, 9, 7),i=1, 2) /5.231d-12,0.6615d0/
  data(fitCoefficients   (i,10, 7),i=1, 2) /9.100d-12,0.6680d0/
  data(fitCoefficients   (i,11, 7),i=1, 2) /1.447d-11,0.6814d0/
  data(fitCoefficients   (i,12, 7),i=1, 2) /2.300d-11,0.6950d0/
  data(fitCoefficients   (i,13, 7),i=1, 2) /3.145d-11,0.6915d0/
  data(fitCoefficients   (i,14, 7),i=1, 2) /4.300d-11,0.6880d0/
  data(fitCoefficients   (i,15, 7),i=1, 2) /5.659d-11,0.7023d0/
  data(fitCoefficients   (i,16, 7),i=1, 2) /7.400d-11,0.7160d0/
  data(fitCoefficients   (i,17, 7),i=1, 2) /9.561d-11,0.7102d0/
  data(fitCoefficients   (i,18, 7),i=1, 2) /1.230d-10,0.7020d0/
  data(fitCoefficients   (i,19, 7),i=1, 2) /1.587d-10,0.7105d0/
  data(fitCoefficients   (i,20, 7),i=1, 2) /2.040d-10,0.7300d0/
  data(fitCoefficients   (i,21, 7),i=1, 2) /2.630d-10,0.7490d0/
  data(fitCoefficients   (i,22, 7),i=1, 2) /3.220d-10,0.7683d0/
  data(fitCoefficients   (i,23, 7),i=1, 2) /3.950d-10,0.7830d0/
  data(fitCoefficients   (i,24, 7),i=1, 2) /4.671d-10,0.7967d0/
  data(fitCoefficients   (i,25, 7),i=1, 2) /5.407d-10,0.8058d0/
  data(fitCoefficients   (i,27, 7),i=1, 2) /6.611d-10,0.8121d0/
  data(fitCoefficients   (i,28, 7),i=1, 2) /7.080d-10,0.8110d0/
  data(fitCoefficients   (i,29, 7),i=1, 2) /7.550d-10,0.8100d0/
  data(fitCoefficients   (i,30, 7),i=1, 2) /8.020d-10,0.8090d0/
  data(fitCoefficients   (i, 8, 8),i=1, 2) /3.100d-13,0.6780d0/
  data(fitCoefficients   (i, 9, 8),i=1, 2) /1.344d-12,0.6708d0/
  data(fitCoefficients   (i,10, 8),i=1, 2) /4.400d-12,0.6750d0/
  data(fitCoefficients   (i,11, 8),i=1, 2) /7.849d-12,0.6952d0/
  data(fitCoefficients   (i,12, 8),i=1, 2) /1.400d-11,0.7160d0/
  data(fitCoefficients   (i,13, 8),i=1, 2) /2.049d-11,0.7090d0/
  data(fitCoefficients   (i,14, 8),i=1, 2) /3.000d-11,0.7020d0/
  data(fitCoefficients   (i,15, 8),i=1, 2) /4.125d-11,0.6965d0/
  data(fitCoefficients   (i,16, 8),i=1, 2) /5.500d-11,0.7110d0/
  data(fitCoefficients   (i,17, 8),i=1, 2) /7.280d-11,0.7518d0/
  data(fitCoefficients   (i,18, 8),i=1, 2) /9.550d-11,0.7930d0/
  data(fitCoefficients   (i,19, 8),i=1, 2) /1.235d-10,0.8052d0/
  data(fitCoefficients   (i,20, 8),i=1, 2) /1.580d-10,0.8000d0/
  data(fitCoefficients   (i,21, 8),i=1, 2) /2.060d-10,0.7990d0/
  data(fitCoefficients   (i,22, 8),i=1, 2) /2.537d-10,0.7977d0/
  data(fitCoefficients   (i,23, 8),i=1, 2) /3.190d-10,0.8020d0/
  data(fitCoefficients   (i,24, 8),i=1, 2) /3.844d-10,0.8071d0/
  data(fitCoefficients   (i,25, 8),i=1, 2) /4.564d-10,0.8124d0/
  data(fitCoefficients   (i,27, 8),i=1, 2) /5.842d-10,0.8168d0/
  data(fitCoefficients   (i,28, 8),i=1, 2) /6.380d-10,0.8160d0/
  data(fitCoefficients   (i,29, 8),i=1, 2) /6.920d-10,0.8150d0/
  data(fitCoefficients   (i,30, 8),i=1, 2) /7.460d-10,0.8140d0/
  data(fitCoefficients   (i, 9, 9),i=1, 2) /6.273d-13,0.6798d0/
  data(fitCoefficients   (i,10, 9),i=1, 2) /1.500d-12,0.6930d0/
  data(fitCoefficients   (i,11, 9),i=1, 2) /3.399d-12,0.7054d0/
  data(fitCoefficients   (i,12, 9),i=1, 2) /7.700d-12,0.7180d0/
  data(fitCoefficients   (i,13, 9),i=1, 2) /1.275d-11,0.7170d0/
  data(fitCoefficients   (i,14, 9),i=1, 2) /2.110d-11,0.7160d0/
  data(fitCoefficients   (i,15, 9),i=1, 2) /2.975d-11,0.6945d0/
  data(fitCoefficients   (i,16, 9),i=1, 2) /4.000d-11,0.6960d0/
  data(fitCoefficients   (i,17, 9),i=1, 2) /5.281d-11,0.7491d0/
  data(fitCoefficients   (i,18, 9),i=1, 2) /6.920d-11,0.8110d0/
  data(fitCoefficients   (i,19, 9),i=1, 2) /9.044d-11,0.8251d0/
  data(fitCoefficients   (i,20, 9),i=1, 2) /1.180d-10,0.8100d0/
  data(fitCoefficients   (i,21, 9),i=1, 2) /1.580d-10,0.8040d0/
  data(fitCoefficients   (i,22, 9),i=1, 2) /1.983d-10,0.7980d0/
  data(fitCoefficients   (i,23, 9),i=1, 2) /2.570d-10,0.8040d0/
  data(fitCoefficients   (i,24, 9),i=1, 2) /3.154d-10,0.8101d0/
  data(fitCoefficients   (i,25, 9),i=1, 2) /3.837d-10,0.8183d0/
  data(fitCoefficients   (i,27, 9),i=1, 2) /5.147d-10,0.8253d0/
  data(fitCoefficients   (i,28, 9),i=1, 2) /5.750d-10,0.8240d0/
  data(fitCoefficients   (i,29, 9),i=1, 2) /6.350d-10,0.8230d0/
  data(fitCoefficients   (i,30, 9),i=1, 2) /6.960d-10,0.8210d0/
  data(fitCoefficients   (i,10,10),i=1, 2) /2.200d-13,0.7590d0/
  data(fitCoefficients   (i,11,10),i=1, 2) /8.775d-13,0.7467d0/
  data(fitCoefficients   (i,12,10),i=1, 2) /3.500d-12,0.7340d0/
  data(fitCoefficients   (i,13,10),i=1, 2) /6.481d-12,0.7345d0/
  data(fitCoefficients   (i,14,10),i=1, 2) /1.200d-11,0.7350d0/
  data(fitCoefficients   (i,15,10),i=1, 2) /1.834d-11,0.7285d0/
  data(fitCoefficients   (i,16,10),i=1, 2) /2.700d-11,0.7330d0/
  data(fitCoefficients   (i,17,10),i=1, 2) /3.711d-11,0.7641d0/
  data(fitCoefficients   (i,18,10),i=1, 2) /4.900d-11,0.8010d0/
  data(fitCoefficients   (i,19,10),i=1, 2) /6.444d-11,0.8175d0/
  data(fitCoefficients   (i,20,10),i=1, 2) /8.510d-11,0.8200d0/
  data(fitCoefficients   (i,21,10),i=1, 2) /1.170d-10,0.8220d0/
  data(fitCoefficients   (i,22,10),i=1, 2) /1.494d-10,0.8242d0/
  data(fitCoefficients   (i,23,10),i=1, 2) /2.010d-10,0.8280d0/
  data(fitCoefficients   (i,24,10),i=1, 2) /2.525d-10,0.8311d0/
  data(fitCoefficients   (i,25,10),i=1, 2) /3.177d-10,0.8341d0/
  data(fitCoefficients   (i,27,10),i=1, 2) /4.552d-10,0.8364d0/
  data(fitCoefficients   (i,28,10),i=1, 2) /5.250d-10,0.8360d0/
  data(fitCoefficients   (i,29,10),i=1, 2) /5.950d-10,0.8360d0/
  data(fitCoefficients   (i,30,10),i=1, 2) /6.650d-10,0.8350d0/
  data(fitCoefficients   (i,12,12),i=1, 2) /1.400d-13,0.8550d0/
  data(fitCoefficients   (i,13,12),i=1, 2) /7.197d-13,0.7697d0/
  data(fitCoefficients   (i,14,12),i=1, 2) /3.700d-12,0.6930d0/
  data(fitCoefficients   (i,15,12),i=1, 2) /7.980d-12,0.6829d0/
  data(fitCoefficients   (i,16,12),i=1, 2) /1.200d-11,0.7010d0/
  data(fitCoefficients   (i,17,12),i=1, 2) /1.800d-11,0.7232d0/
  data(fitCoefficients   (i,18,12),i=1, 2) /2.690d-11,0.7440d0/
  data(fitCoefficients   (i,19,12),i=1, 2) /3.748d-11,0.7628d0/
  data(fitCoefficients   (i,20,12),i=1, 2) /5.040d-11,0.7800d0/
  data(fitCoefficients   (i,21,12),i=1, 2) /7.240d-11,0.7950d0/
  data(fitCoefficients   (i,22,12),i=1, 2) /9.440d-11,0.8107d0/
  data(fitCoefficients   (i,23,12),i=1, 2) /1.350d-10,0.8220d0/
  data(fitCoefficients   (i,24,12),i=1, 2) /1.751d-10,0.8340d0/
  data(fitCoefficients   (i,25,12),i=1, 2) /2.298d-10,0.8417d0/
  data(fitCoefficients   (i,27,12),i=1, 2) /3.461d-10,0.8469d0/
  data(fitCoefficients   (i,28,12),i=1, 2) /4.030d-10,0.8460d0/
  data(fitCoefficients   (i,29,12),i=1, 2) /4.600d-10,0.8450d0/
  data(fitCoefficients   (i,30,12),i=1, 2) /5.170d-10,0.8440d0/
  data(fitCoefficients   (i,13,13),i=1, 2) /3.980d-13,0.8019d0/
  data(fitCoefficients   (i,14,13),i=1, 2) /1.000d-12,0.7860d0/
  data(fitCoefficients   (i,15,13),i=1, 2) /2.558d-12,0.7629d0/
  data(fitCoefficients   (i,16,13),i=1, 2) /5.700d-12,0.7550d0/
  data(fitCoefficients   (i,17,13),i=1, 2) /1.011d-11,0.7703d0/
  data(fitCoefficients   (i,18,13),i=1, 2) /1.580d-11,0.7930d0/
  data(fitCoefficients   (i,19,13),i=1, 2) /2.448d-11,0.8052d0/
  data(fitCoefficients   (i,20,13),i=1, 2) /3.760d-11,0.8100d0/
  data(fitCoefficients   (i,21,13),i=1, 2) /5.870d-11,0.8150d0/
  data(fitCoefficients   (i,22,13),i=1, 2) /7.972d-11,0.8206d0/
  data(fitCoefficients   (i,23,13),i=1, 2) /1.130d-10,0.8270d0/
  data(fitCoefficients   (i,24,13),i=1, 2) /1.470d-10,0.8325d0/
  data(fitCoefficients   (i,25,13),i=1, 2) /1.908d-10,0.8372d0/
  data(fitCoefficients   (i,27,13),i=1, 2) /2.976d-10,0.8406d0/
  data(fitCoefficients   (i,28,13),i=1, 2) /3.630d-10,0.8400d0/
  data(fitCoefficients   (i,29,13),i=1, 2) /4.280d-10,0.8390d0/
  data(fitCoefficients   (i,30,13),i=1, 2) /4.940d-10,0.8390d0/
  data(fitCoefficients   (i,14,14),i=1, 2) /5.900d-13,0.6010d0/
  data(fitCoefficients   (i,15,14),i=1, 2) /1.294d-12,0.6766d0/
  data(fitCoefficients   (i,16,14),i=1, 2) /2.700d-12,0.7450d0/
  data(fitCoefficients   (i,17,14),i=1, 2) /5.165d-12,0.7893d0/
  data(fitCoefficients   (i,18,14),i=1, 2) /9.120d-12,0.8110d0/
  data(fitCoefficients   (i,19,14),i=1, 2) /1.513d-11,0.8186d0/
  data(fitCoefficients   (i,20,14),i=1, 2) /2.400d-11,0.8200d0/
  data(fitCoefficients   (i,21,14),i=1, 2) /3.960d-11,0.8220d0/
  data(fitCoefficients   (i,22,14),i=1, 2) /5.518d-11,0.8245d0/
  data(fitCoefficients   (i,23,14),i=1, 2) /8.370d-11,0.8280d0/
  data(fitCoefficients   (i,24,14),i=1, 2) /1.123d-10,0.8313d0/
  data(fitCoefficients   (i,25,14),i=1, 2) /1.525d-10,0.8342d0/
  data(fitCoefficients   (i,26,14),i=1, 2) /2.000d-10,0.8360d0/
  data(fitCoefficients   (i,27,14),i=1, 2) /2.537d-10,0.8364d0/
  data(fitCoefficients   (i,28,14),i=1, 2) /3.160d-10,0.8360d0/
  data(fitCoefficients   (i,29,14),i=1, 2) /3.780d-10,0.8360d0/
  data(fitCoefficients   (i,30,14),i=1, 2) /4.410d-10,0.8350d0/
  data(fitCoefficients   (i,15,15),i=1, 2) /9.761d-13,0.6209d0/
  data(fitCoefficients   (i,16,15),i=1, 2) /1.800d-12,0.6860d0/
  data(fitCoefficients   (i,17,15),i=1, 2) /3.320d-12,0.7579d0/
  data(fitCoefficients   (i,18,15),i=1, 2) /6.030d-12,0.8120d0/
  data(fitCoefficients   (i,19,15),i=1, 2) /1.063d-11,0.8269d0/
  data(fitCoefficients   (i,20,15),i=1, 2) /1.800d-11,0.8200d0/
  data(fitCoefficients   (i,21,15),i=1, 2) /3.130d-11,0.8180d0/
  data(fitCoefficients   (i,22,15),i=1, 2) /4.451d-11,0.8153d0/
  data(fitCoefficients   (i,23,15),i=1, 2) /6.830d-11,0.8200d0/
  data(fitCoefficients   (i,24,15),i=1, 2) /9.206d-11,0.8246d0/
  data(fitCoefficients   (i,25,15),i=1, 2) /1.250d-10,0.8302d0/
  data(fitCoefficients   (i,26,15),i=1, 2) /1.640d-10,0.8340d0/
  data(fitCoefficients   (i,27,15),i=1, 2) /2.093d-10,0.8349d0/
  data(fitCoefficients   (i,28,15),i=1, 2) /2.630d-10,0.8340d0/
  data(fitCoefficients   (i,29,15),i=1, 2) /3.170d-10,0.8330d0/
  data(fitCoefficients   (i,30,15),i=1, 2) /3.700d-10,0.8320d0/
  data(fitCoefficients   (i,16,16),i=1, 2) /4.100d-13,0.6300d0/
  data(fitCoefficients   (i,17,16),i=1, 2) /1.248d-12,0.7663d0/
  data(fitCoefficients   (i,18,16),i=1, 2) /3.230d-12,0.8690d0/
  data(fitCoefficients   (i,19,16),i=1, 2) /6.384d-12,0.8790d0/
  data(fitCoefficients   (i,20,16),i=1, 2) /1.070d-11,0.8400d0/
  data(fitCoefficients   (i,21,16),i=1, 2) /1.920d-11,0.8210d0/
  data(fitCoefficients   (i,22,16),i=1, 2) /2.765d-11,0.8012d0/
  data(fitCoefficients   (i,23,16),i=1, 2) /4.650d-11,0.8060d0/
  data(fitCoefficients   (i,24,16),i=1, 2) /6.539d-11,0.8099d0/
  data(fitCoefficients   (i,25,16),i=1, 2) /9.539d-11,0.8202d0/
  data(fitCoefficients   (i,26,16),i=1, 2) /1.330d-10,0.8280d0/
  data(fitCoefficients   (i,27,16),i=1, 2) /1.769d-10,0.8299d0/
  data(fitCoefficients   (i,28,16),i=1, 2) /2.290d-10,0.8280d0/
  data(fitCoefficients   (i,29,16),i=1, 2) /2.810d-10,0.8260d0/
  data(fitCoefficients   (i,30,16),i=1, 2) /3.330d-10,0.8240d0/
  data(fitCoefficients   (i,17,17),i=1, 2) /1.010d-12,0.7380d0/
  data(fitCoefficients   (i,18,17),i=1, 2) /1.950d-12,0.7520d0/
  data(fitCoefficients   (i,19,17),i=1, 2) /3.766d-12,0.7662d0/
  data(fitCoefficients   (i,20,17),i=1, 2) /7.080d-12,0.7800d0/
  data(fitCoefficients   (i,21,17),i=1, 2) /1.430d-11,0.7920d0/
  data(fitCoefficients   (i,22,17),i=1, 2) /2.152d-11,0.8038d0/
  data(fitCoefficients   (i,23,17),i=1, 2) /3.740d-11,0.8120d0/
  data(fitCoefficients   (i,24,17),i=1, 2) /5.335d-11,0.8207d0/
  data(fitCoefficients   (i,25,17),i=1, 2) /7.807d-11,0.8260d0/
  data(fitCoefficients   (i,26,17),i=1, 2) /1.090d-10,0.8290d0/
  data(fitCoefficients   (i,27,17),i=1, 2) /1.459d-10,0.8296d0/
  data(fitCoefficients   (i,28,17),i=1, 2) /1.910d-10,0.8290d0/
  data(fitCoefficients   (i,29,17),i=1, 2) /2.360d-10,0.8280d0/
  data(fitCoefficients   (i,30,17),i=1, 2) /2.810d-10,0.8280d0/
  data(fitCoefficients   (i,18,18),i=1, 2) /3.770d-13,0.6510d0/
  data(fitCoefficients   (i,19,18),i=1, 2) /1.304d-12,0.6753d0/
  data(fitCoefficients   (i,20,18),i=1, 2) /3.960d-12,0.7000d0/
  data(fitCoefficients   (i,21,18),i=1, 2) /1.130d-11,0.7240d0/
  data(fitCoefficients   (i,22,18),i=1, 2) /1.857d-11,0.7484d0/
  data(fitCoefficients   (i,23,18),i=1, 2) /3.170d-11,0.7680d0/
  data(fitCoefficients   (i,24,18),i=1, 2) /4.479d-11,0.7883d0/
  data(fitCoefficients   (i,25,18),i=1, 2) /6.106d-11,0.8020d0/
  data(fitCoefficients   (i,26,18),i=1, 2) /8.130d-11,0.8100d0/
  data(fitCoefficients   (i,27,18),i=1, 2) /1.098d-10,0.8118d0/
  data(fitCoefficients   (i,28,18),i=1, 2) /1.500d-10,0.8100d0/
  data(fitCoefficients   (i,29,18),i=1, 2) /1.900d-10,0.8080d0/
  data(fitCoefficients   (i,30,18),i=1, 2) /2.300d-10,0.8060d0/
  data(fitCoefficients   (i,19,19),i=1, 2) /2.762d-13,0.8023d0/
  data(fitCoefficients   (i,20,19),i=1, 2) /6.780d-13,0.8000d0/
  data(fitCoefficients   (i,21,19),i=1, 2) /2.330d-12,0.7980d0/
  data(fitCoefficients   (i,22,19),i=1, 2) /3.983d-12,0.7955d0/
  data(fitCoefficients   (i,23,19),i=1, 2) /1.150d-11,0.7940d0/
  data(fitCoefficients   (i,24,19),i=1, 2) /1.906d-11,0.7919d0/
  data(fitCoefficients   (i,25,19),i=1, 2) /3.620d-11,0.7907d0/
  data(fitCoefficients   (i,26,19),i=1, 2) /6.050d-11,0.7900d0/
  data(fitCoefficients   (i,27,19),i=1, 2) /8.818d-11,0.7898d0/
  data(fitCoefficients   (i,28,19),i=1, 2) /1.190d-10,0.7900d0/
  data(fitCoefficients   (i,29,19),i=1, 2) /1.500d-10,0.7900d0/
  data(fitCoefficients   (i,30,19),i=1, 2) /1.810d-10,0.7900d0/
  data(fitCoefficients   (i,20,20),i=1, 2) /1.120d-13,0.9000d0/
  data(fitCoefficients   (i,21,20),i=1, 2) /6.540d-13,0.8670d0/
  data(fitCoefficients   (i,22,20),i=1, 2) /1.196d-12,0.8344d0/
  data(fitCoefficients   (i,23,20),i=1, 2) /5.330d-12,0.8090d0/
  data(fitCoefficients   (i,24,20),i=1, 2) /9.471d-12,0.7846d0/
  data(fitCoefficients   (i,25,20),i=1, 2) /2.169d-11,0.7683d0/
  data(fitCoefficients   (i,26,20),i=1, 2) /4.120d-11,0.7590d0/
  data(fitCoefficients   (i,27,20),i=1, 2) /6.409d-11,0.7570d0/
  data(fitCoefficients   (i,28,20),i=1, 2) /8.910d-11,0.7590d0/
  data(fitCoefficients   (i,29,20),i=1, 2) /1.140d-10,0.7610d0/
  data(fitCoefficients   (i,30,20),i=1, 2) /1.390d-10,0.7630d0/
  data(fitCoefficients   (i,21,21),i=1, 2) /1.170d-13,0.8980d0/
  data(fitCoefficients   (i,22,21),i=1, 2) /5.330d-12,0.8640d0/
  data(fitCoefficients   (i,23,21),i=1, 2) /1.060d-11,0.8300d0/
  data(fitCoefficients   (i,24,21),i=1, 2) /1.580d-11,0.7960d0/
  data(fitCoefficients   (i,25,21),i=1, 2) /2.100d-11,0.7620d0/
  data(fitCoefficients   (i,26,21),i=1, 2) /2.620d-11,0.7280d0/
  data(fitCoefficients   (i,27,21),i=1, 2) /2.822d-11,0.7280d0/
  data(fitCoefficients   (i,28,21),i=1, 2) /3.040d-11,0.7280d0/
  data(fitCoefficients   (i,29,21),i=1, 2) /3.260d-11,0.7280d0/
  data(fitCoefficients   (i,30,21),i=1, 2) /3.480d-11,0.7280d0/
  data(fitCoefficients   (i,22,22),i=1, 2) /1.220d-13,0.8970d0/
  data(fitCoefficients   (i,23,22),i=1, 2) /3.870d-12,0.8480d0/
  data(fitCoefficients   (i,24,22),i=1, 2) /7.610d-12,0.7980d0/
  data(fitCoefficients   (i,25,22),i=1, 2) /1.140d-11,0.7480d0/
  data(fitCoefficients   (i,26,22),i=1, 2) /1.510d-11,0.6990d0/
  data(fitCoefficients   (i,27,22),i=1, 2) /1.626d-11,0.6990d0/
  data(fitCoefficients   (i,28,22),i=1, 2) /1.750d-11,0.6990d0/
  data(fitCoefficients   (i,29,22),i=1, 2) /1.870d-11,0.6990d0/
  data(fitCoefficients   (i,30,22),i=1, 2) /2.000d-11,0.6990d0/
  data(fitCoefficients   (i,23,23),i=1, 2) /1.270d-13,0.8950d0/
  data(fitCoefficients   (i,24,23),i=1, 2) /2.680d-12,0.8240d0/
  data(fitCoefficients   (i,25,23),i=1, 2) /5.240d-12,0.7530d0/
  data(fitCoefficients   (i,26,23),i=1, 2) /7.800d-12,0.6820d0/
  data(fitCoefficients   (i,27,23),i=1, 2) /8.402d-12,0.6820d0/
  data(fitCoefficients   (i,28,23),i=1, 2) /9.050d-12,0.6820d0/
  data(fitCoefficients   (i,29,23),i=1, 2) /9.700d-12,0.6820d0/
  data(fitCoefficients   (i,30,23),i=1, 2) /1.030d-11,0.6820d0/
  data(fitCoefficients   (i,24,24),i=1, 2) /1.320d-13,0.8940d0/
  data(fitCoefficients   (i,25,24),i=1, 2) /1.730d-12,0.8200d0/
  data(fitCoefficients   (i,26,24),i=1, 2) /3.320d-12,0.7460d0/
  data(fitCoefficients   (i,27,24),i=1, 2) /3.575d-12,0.7460d0/
  data(fitCoefficients   (i,28,24),i=1, 2) /3.580d-12,0.7460d0/
  data(fitCoefficients   (i,29,24),i=1, 2) /3.580d-12,0.7460d0/
  data(fitCoefficients   (i,30,24),i=1, 2) /3.590d-12,0.7460d0/
  data(fitCoefficients   (i,25,25),i=1, 2) /1.370d-13,0.8920d0/
  data(fitCoefficients   (i,26,25),i=1, 2) /1.020d-12,0.8430d0/
  data(fitCoefficients   (i,27,25),i=1, 2) /1.278d-12,0.7682d0/
  data(fitCoefficients   (i,28,25),i=1, 2) /1.600d-12,0.7000d0/
  data(fitCoefficients   (i,29,25),i=1, 2) /1.920d-12,0.7000d0/
  data(fitCoefficients   (i,30,25),i=1, 2) /2.240d-12,0.7000d0/
  data(fitCoefficients   (i,26,26),i=1, 2) /1.420d-13,0.8910d0/
  data(fitCoefficients   (i,27,26),i=1, 2) /4.459d-13,0.7897d0/
  data(fitCoefficients   (i,28,26),i=1, 2) /1.400d-12,0.7000d0/
  data(fitCoefficients   (i,29,26),i=1, 2) /1.500d-12,0.7000d0/
  data(fitCoefficients   (i,30,26),i=1, 2) /1.600d-12,0.7000d0/
  data(fitCoefficients   (i,27,27),i=1, 2) /2.510d-13,0.7950d0/
  data(fitCoefficients   (i,28,27),i=1, 2) /1.000d-12,0.7000d0/
  data(fitCoefficients   (i,29,27),i=1, 2) /1.000d-12,0.7000d0/
  data(fitCoefficients   (i,30,27),i=1, 2) /1.100d-12,0.7000d0/
  data(fitCoefficients   (i,28,28),i=1, 2) /3.600d-13,0.7000d0/
  data(fitCoefficients   (i,29,28),i=1, 2) /3.600d-13,0.7000d0/
  data(fitCoefficients   (i,30,28),i=1, 2) /3.600d-13,0.7000d0/
  data(fitCoefficients   (i,29,29),i=1, 2) /3.600d-13,0.7000d0/
  data(fitCoefficients   (i,30,29),i=1, 2) /3.600d-13,0.7000d0/
  data(fitCoefficients   (i,30,30),i=1, 2) /3.600d-13,0.7000d0/
  data(fitCoefficientsNew(i, 1, 1),i=1, 4) /7.982d-11,0.7480d0,3.148d+00,7.036d+05/
  data(fitCoefficientsNew(i, 2, 1),i=1, 4) /1.891d-10,0.7524d0,9.370d+00,2.774d+06/
  data(fitCoefficientsNew(i, 3, 1),i=1, 4) /3.039d-10,0.7539d0,1.871d+01,6.209d+06/
  data(fitCoefficientsNew(i, 4, 1),i=1, 4) /4.290d-10,0.7557d0,3.000d+01,1.093d+07/
  data(fitCoefficientsNew(i, 5, 1),i=1, 4) /5.437d-10,0.7560d0,4.576d+01,1.706d+07/
  data(fitCoefficientsNew(i, 6, 1),i=1, 4) /6.556d-10,0.7567d0,6.523d+01,2.446d+07/
  data(fitCoefficientsNew(i, 7, 1),i=1, 4) /7.586d-10,0.7563d0,9.015d+01,3.338d+07/
  data(fitCoefficientsNew(i, 8, 1),i=1, 4) /8.616d-10,0.7563d0,1.191d+02,4.352d+07/
  data(fitCoefficientsNew(i, 9, 1),i=1, 4) /9.712d-10,0.7566d0,1.499d+02,5.498d+07/
  data(fitCoefficientsNew(i,10, 1),i=1, 4) /1.085d-09,0.7570d0,1.834d+02,6.776d+07/
  data(fitCoefficientsNew(i,11, 1),i=1, 4) /1.163d-09,0.7558d0,2.328d+02,8.262d+07/
  data(fitCoefficientsNew(i,12, 1),i=1, 4) /1.317d-09,0.7574d0,2.585d+02,9.769d+07/
  data(fitCoefficientsNew(i,13, 1),i=1, 4) /1.419d-09,0.7578d0,3.057d+02,1.143d+08/
  data(fitCoefficientsNew(i,14, 1),i=1, 4) /1.517d-09,0.7574d0,3.601d+02,1.329d+08/
  data(fitCoefficientsNew(i,15, 1),i=1, 4) /1.586d-09,0.7560d0,4.327d+02,1.534d+08/
  data(fitCoefficientsNew(i,16, 1),i=1, 4) /1.729d-09,0.7568d0,4.725d+02,1.746d+08/
  data(fitCoefficientsNew(i,17, 1),i=1, 4) /1.791d-09,0.7565d0,5.591d+02,1.972d+08/
  data(fitCoefficientsNew(i,18, 1),i=1, 4) /1.913d-09,0.7567d0,6.175d+02,2.212d+08/
  data(fitCoefficientsNew(i,19, 1),i=1, 4) /2.033d-09,0.7569d0,6.797d+02,2.463d+08/
  data(fitCoefficientsNew(i,20, 1),i=1, 4) /2.129d-09,0.7570d0,7.591d+02,2.739d+08/
  data(fitCoefficientsNew(i,21, 1),i=1, 4) /2.262d-09,0.7578d0,8.186d+02,3.000d+08/
  data(fitCoefficientsNew(i,22, 1),i=1, 4) /2.370d-09,0.7574d0,9.002d+02,3.307d+08/
  data(fitCoefficientsNew(i,23, 1),i=1, 4) /2.415d-09,0.7565d0,1.032d+03,3.635d+08/
  data(fitCoefficientsNew(i,24, 1),i=1, 4) /2.537d-09,0.7571d0,1.108d+03,3.954d+08/
  data(fitCoefficientsNew(i,25, 1),i=1, 4) /2.618d-09,0.7565d0,1.225d+03,4.307d+08/
  data(fitCoefficientsNew(i,26, 1),i=1, 4) /2.735d-09,0.7568d0,1.314d+03,4.659d+08/
  data(fitCoefficientsNew(i,27, 1),i=1, 4) /2.809d-09,0.7565d0,1.444d+03,5.042d+08/
  data(fitCoefficientsNew(i,28, 1),i=1, 4) /3.002d-09,0.7581d0,1.467d+03,5.409d+08/
  data(fitCoefficientsNew(i,29, 1),i=1, 4) /3.022d-09,0.7564d0,1.666d+03,5.855d+08/
  data(fitCoefficientsNew(i,30, 1),i=1, 4) /3.127d-09,0.7567d0,1.779d+03,6.246d+08/
  data(fitCoefficientsNew(i, 2, 2),i=1, 4) /9.356d-10,0.7892d0,4.266d-02,4.677d+06/
  data(fitCoefficientsNew(i, 3, 2),i=1, 4) /1.112d-10,0.6926d0,2.437d+01,8.323d+06/
  data(fitCoefficientsNew(i, 4, 2),i=1, 4) /1.317d-10,0.6691d0,8.473d+01,1.412d+07/
  data(fitCoefficientsNew(i, 5, 2),i=1, 4) /1.922d-10,0.6717d0,1.272d+02,1.975d+07/
  data(fitCoefficientsNew(i, 6, 2),i=1, 4) /2.765d-10,0.6858d0,1.535d+02,2.556d+07/
  data(fitCoefficientsNew(i, 7, 2),i=1, 4) /3.910d-10,0.6988d0,1.611d+02,3.271d+07/
  data(fitCoefficientsNew(i, 8, 2),i=1, 4) /4.897d-10,0.7048d0,1.906d+02,4.093d+07/
  data(fitCoefficientsNew(i, 9, 2),i=1, 4) /5.602d-10,0.7052d0,2.476d+02,5.077d+07/
  data(fitCoefficientsNew(i,10, 2),i=1, 4) /6.161d-10,0.7029d0,3.274d+02,6.243d+07/
  data(fitCoefficientsNew(i,11, 2),i=1, 4) /6.833d-10,0.7018d0,4.060d+02,7.491d+07/
  data(fitCoefficientsNew(i,12, 2),i=1, 4) /7.510d-10,0.7020d0,4.921d+02,8.643d+07/
  data(fitCoefficientsNew(i,13, 2),i=1, 4) /8.182d-10,0.7008d0,5.875d+02,1.007d+08/
  data(fitCoefficientsNew(i,14, 2),i=1, 4) /8.722d-10,0.6996d0,7.098d+02,1.155d+08/
  data(fitCoefficientsNew(i,15, 2),i=1, 4) /9.142d-10,0.6961d0,8.682d+02,1.335d+08/
  data(fitCoefficientsNew(i,16, 2),i=1, 4) /9.692d-10,0.6945d0,1.017d+03,1.517d+08/
  data(fitCoefficientsNew(i,17, 2),i=1, 4) /1.021d-09,0.6932d0,1.184d+03,1.695d+08/
  data(fitCoefficientsNew(i,18, 2),i=1, 4) /1.087d-09,0.6936d0,1.329d+03,1.880d+08/
  data(fitCoefficientsNew(i,19, 2),i=1, 4) /1.145d-09,0.6921d0,1.503d+03,2.098d+08/
  data(fitCoefficientsNew(i,20, 2),i=1, 4) /1.179d-09,0.6893d0,1.757d+03,2.344d+08/
  data(fitCoefficientsNew(i,21, 2),i=1, 4) /1.265d-09,0.6902d0,1.877d+03,2.555d+08/
  data(fitCoefficientsNew(i,22, 2),i=1, 4) /1.322d-09,0.6885d0,2.092d+03,2.829d+08/
  data(fitCoefficientsNew(i,23, 2),i=1, 4) /1.375d-09,0.6885d0,2.321d+03,3.056d+08/
  data(fitCoefficientsNew(i,24, 2),i=1, 4) /1.422d-09,0.6874d0,2.589d+03,3.336d+08/
  data(fitCoefficientsNew(i,25, 2),i=1, 4) /1.488d-09,0.6867d0,2.802d+03,3.623d+08/
  data(fitCoefficientsNew(i,26, 2),i=1, 4) /1.542d-09,0.6859d0,3.073d+03,3.926d+08/
  data(fitCoefficientsNew(i,27, 2),i=1, 4) /1.589d-09,0.6846d0,3.373d+03,4.267d+08/
  data(fitCoefficientsNew(i,28, 2),i=1, 4) /1.676d-09,0.6861d0,3.530d+03,4.538d+08/
  data(fitCoefficientsNew(i,29, 2),i=1, 4) /1.686d-09,0.6824d0,4.031d+03,4.948d+08/
  data(fitCoefficientsNew(i,30, 2),i=1, 4) /1.758d-09,0.6834d0,4.254d+03,5.258d+08/
  data(fitCoefficientsNew(i, 3, 3),i=1, 4) /1.036d-11,0.3880d0,1.077d+02,1.177d+07/
  data(fitCoefficientsNew(i, 4, 3),i=1, 4) /2.338d-11,0.4211d0,3.647d+02,1.215d+07/
  data(fitCoefficientsNew(i, 5, 3),i=1, 4) /4.487d-11,0.4644d0,5.371d+02,1.465d+07/
  data(fitCoefficientsNew(i, 6, 3),i=1, 4) /8.540d-11,0.5247d0,5.014d+02,1.479d+07/
  data(fitCoefficientsNew(i, 7, 3),i=1, 4) /1.169d-10,0.5470d0,6.793d+02,1.650d+07/
  data(fitCoefficientsNew(i, 8, 3),i=1, 4) /2.053d-10,0.6019d0,4.772d+02,1.711d+07/
  data(fitCoefficientsNew(i, 9, 3),i=1, 4) /2.739d-10,0.6188d0,5.033d+02,2.064d+07/
  data(fitCoefficientsNew(i,10, 3),i=1, 4) /3.200d-10,0.6198d0,6.329d+02,2.616d+07/
  data(fitCoefficientsNew(i,11, 3),i=1, 4) /3.873d-10,0.6295d0,7.000d+02,2.989d+07/
  data(fitCoefficientsNew(i,12, 3),i=1, 4) /4.284d-10,0.6287d0,8.748d+02,3.586d+07/
  data(fitCoefficientsNew(i,13, 3),i=1, 4) /4.881d-10,0.6326d0,9.941d+02,4.085d+07/
  data(fitCoefficientsNew(i,14, 3),i=1, 4) /5.373d-10,0.6337d0,1.164d+03,4.677d+07/
  data(fitCoefficientsNew(i,15, 3),i=1, 4) /5.876d-10,0.6354d0,1.341d+03,5.292d+07/
  data(fitCoefficientsNew(i,16, 3),i=1, 4) /6.571d-10,0.6400d0,1.452d+03,5.796d+07/
  data(fitCoefficientsNew(i,17, 3),i=1, 4) /7.076d-10,0.6397d0,1.653d+03,6.555d+07/
  data(fitCoefficientsNew(i,18, 3),i=1, 4) /7.538d-10,0.6388d0,1.889d+03,7.306d+07/
  data(fitCoefficientsNew(i,19, 3),i=1, 4) /8.182d-10,0.6411d0,2.044d+03,8.057d+07/
  data(fitCoefficientsNew(i,20, 3),i=1, 4) /8.577d-10,0.6403d0,2.334d+03,8.850d+07/
  data(fitCoefficientsNew(i,21, 3),i=1, 4) /9.162d-10,0.6413d0,2.543d+03,9.690d+07/
  data(fitCoefficientsNew(i,22, 3),i=1, 4) /9.844d-10,0.6440d0,2.708d+03,1.044d+08/
  data(fitCoefficientsNew(i,23, 3),i=1, 4) /1.020d-09,0.6427d0,3.057d+03,1.140d+08/
  data(fitCoefficientsNew(i,24, 3),i=1, 4) /1.091d-09,0.6445d0,3.225d+03,1.229d+08/
  data(fitCoefficientsNew(i,25, 3),i=1, 4) /1.151d-09,0.6451d0,3.461d+03,1.334d+08/
  data(fitCoefficientsNew(i,26, 3),i=1, 4) /1.198d-09,0.6443d0,3.789d+03,1.437d+08/
  data(fitCoefficientsNew(i,27, 3),i=1, 4) /1.211d-09,0.6406d0,4.357d+03,1.572d+08/
  data(fitCoefficientsNew(i,28, 3),i=1, 4) /1.288d-09,0.6440d0,4.506d+03,1.651d+08/
  data(fitCoefficientsNew(i,29, 3),i=1, 4) /1.372d-09,0.6472d0,4.627d+03,1.740d+08/
  data(fitCoefficientsNew(i,30, 3),i=1, 4) /1.412d-09,0.6454d0,5.053d+03,1.891d+08/
  data(fitCoefficientsNew(i,11,11),i=1, 4) /5.641d-12,0.1749d0,3.077d+02,2.617d+06/
  data(fitCoefficientsNew(i,12,11),i=1, 4) /1.920d-11,0.3028d0,4.849d+02,5.890d+06/
  data(fitCoefficientsNew(i,13,11),i=1, 4) /3.753d-11,0.3585d0,6.848d+02,9.035d+06/
  data(fitCoefficientsNew(i,14,11),i=1, 4) /5.942d-11,0.3930d0,8.962d+02,1.213d+07/
  data(fitCoefficientsNew(i,15,11),i=1, 4) /1.721d-10,0.5429d0,2.848d+02,3.975d+07/
  data(fitCoefficientsNew(i,16,11),i=1, 4) /3.502d-10,0.6266d0,1.532d+02,1.755d+07/
  data(fitCoefficientsNew(i,17,11),i=1, 4) /2.502d-10,0.5580d0,5.303d+02,4.558d+07/
  data(fitCoefficientsNew(i,18,11),i=1, 4) /2.862d-10,0.5621d0,7.002d+02,4.885d+07/
  data(fitCoefficientsNew(i,19,11),i=1, 4) /2.757d-10,0.5364d0,1.204d+03,7.013d+07/
  data(fitCoefficientsNew(i,20,11),i=1, 4) /5.273d-10,0.6281d0,5.329d+02,3.188d+07/
  data(fitCoefficientsNew(i,21,11),i=1, 4) /3.890d-10,0.5645d0,1.391d+03,6.295d+07/
  data(fitCoefficientsNew(i,22,11),i=1, 4) /4.207d-10,0.5646d0,1.688d+03,6.872d+07/
  data(fitCoefficientsNew(i,23,11),i=1, 4) /4.605d-10,0.5659d0,1.949d+03,7.419d+07/
  data(fitCoefficientsNew(i,24,11),i=1, 4) /4.975d-10,0.5655d0,2.257d+03,8.072d+07/
  data(fitCoefficientsNew(i,25,11),i=1, 4) /5.349d-10,0.5658d0,2.577d+03,8.710d+07/
  data(fitCoefficientsNew(i,26,11),i=1, 4) /7.688d-10,0.6173d0,1.653d+03,6.161d+07/
  data(fitCoefficientsNew(i,27,11),i=1, 4) /5.850d-10,0.5598d0,3.538d+03,1.052d+08/
  data(fitCoefficientsNew(i,28,11),i=1, 4) /6.347d-10,0.5631d0,3.780d+03,1.116d+08/
  data(fitCoefficientsNew(i,29,11),i=1, 4) /6.619d-10,0.5602d0,4.322d+03,1.210d+08/
  data(fitCoefficientsNew(i,30,11),i=1, 4) /7.002d-10,0.5612d0,4.726d+03,1.287d+08/
  data(fitCoefficientsNew(i, 6, 4),i=1, 4) /2.020d-09,0.7798d0,6.690d-01,2.425d+06/
  data(fitCoefficientsNew(i, 7, 4),i=1, 4) /1.595d-11,0.3529d0,9.870d+03,2.584d+07/
  data(fitCoefficientsNew(i, 8, 4),i=1, 4) /2.008d-11,0.3567d0,1.520d+04,3.843d+07/
  data(fitCoefficientsNew(i,10, 4),i=1, 4) /2.793d-11,0.3533d0,3.017d+04,7.872d+07/
  data(fitCoefficientsNew(i, 6, 5),i=1, 4) /8.577d-10,0.7837d0,7.286d-01,1.140d+07/
  data(fitCoefficientsNew(i, 7, 5),i=1, 4) /7.039d-10,0.8607d0,2.203d+00,3.029d+06/
  data(fitCoefficientsNew(i, 8, 5),i=1, 4) /1.542d-10,0.6712d0,1.775d+02,1.535d+08/
  data(fitCoefficientsNew(i,10, 5),i=1, 4) /2.515d-10,0.7011d0,3.028d+02,8.903d+06/
  data(fitCoefficientsNew(i, 6, 6),i=1, 4) /7.651d-09,0.8027d0,1.193d-03,9.334d+12/
  data(fitCoefficientsNew(i, 7, 6),i=1, 4) /5.989d-10,0.7560d0,1.377d+00,1.517d+09/
  data(fitCoefficientsNew(i, 8, 6),i=1, 4) /3.672d-09,0.7676d0,2.900d-01,8.521d+07/
  data(fitCoefficientsNew(i,10, 6),i=1, 4) /9.353d-11,0.6270d0,9.031d+02,2.387d+07/
  data(fitCoefficientsNew(i, 7, 7),i=1, 4) /1.243d-08,0.7022d0,1.136d-03,1.015d+13/
  data(fitCoefficientsNew(i, 8, 7),i=1, 4) /1.816d-08,0.7170d0,6.717d-03,3.286d+12/
  data(fitCoefficientsNew(i,10, 7),i=1, 4) /4.227d-11,0.5395d0,1.687d+03,1.491d+17/
  data(fitCoefficientsNew(i, 8, 8),i=1, 4) /1.341d-10,0.6159d0,1.673d+00,6.366d+16/
  data(fitCoefficientsNew(i,10, 8),i=1, 4) /9.563d-12,0.3067d0,9.768d+03,4.851d+17/
  data(fitCoefficientsNew(i,10, 9),i=1, 4) /5.417d-08,0.6930d0,1.179d-03,1.060d+07/
  data(fitCoefficientsNew(i,10,10),i=1, 4) /5.023d-12,0.2420d0,3.181d+02,1.450d+18/
  data(fitCoefficientsNew(i,26,12),i=1, 4) /8.110d-10,0.5442d0,1.755d+03,6.799d+07/
  data(fitCoefficientsNew(i,26,13),i=1, 4) /6.967d-10,0.5602d0,1.727d+03,5.618d+07/
  data(fitCoefficientsNew(i,26,14),i=1, 4) /3.236d-08,0.3247d0,2.338d+01,8.337d+12/
  data(fitCoefficientsNew(i,26,15),i=1, 4) /2.664d-08,0.3285d0,2.297d+01,6.672d+12/
  data(fitCoefficientsNew(i,26,16),i=1, 4) /2.165d-08,0.3403d0,2.195d+01,6.383d+12/
  data(fitCoefficientsNew(i,26,17),i=1, 4) /2.460d-08,0.3387d0,1.487d+01,5.228d+12/
  data(fitCoefficientsNew(i,26,18),i=1, 4) /1.907d-08,0.3768d0,1.216d+01,5.431d+12/
  data(fitCoefficientsNew(i,26,19),i=1, 4) /1.439d-08,0.4170d0,1.006d+01,4.898d+12/
  data(fitCoefficientsNew(i,26,20),i=1, 4) /1.184d-08,0.4798d0,5.883d+00,2.582d+12/
  data(fitCoefficientsNew(i,26,21),i=1, 4) /1.036d-08,0.5428d0,2.743d+00,1.014d+12/
  data(fitCoefficientsNew(i,26,22),i=1, 4) /8.288d-09,0.6012d0,1.216d+00,1.182d+12/
  data(fitCoefficientsNew(i,26,23),i=1, 4) /6.330d-09,0.6355d0,5.457d-01,8.545d+11/
  data(fitCoefficientsNew(i,26,24),i=1, 4) /5.422d-09,0.5067d0,4.998d-01,8.079d+11/
  data(fitCoefficientsNew(i,26,25),i=1, 4) /6.076d-09,0.3112d0,3.401d-01,1.960d+12/
  data(fitCoefficientsNew(i,26,26),i=1, 4) /8.945d-09,0.2156d0,4.184d-02,5.353d+13/
  data(coefficientsIron  (1,i    ),i=1,10) /4.33d-10,3.91d-10,3.49d-10,3.16d-10,2.96d-10,2.59d-10,2.24d-10,1.91d-10,1.68d-10,1.46d-10/
  data(coefficientsIron  (2,i    ),i=1,10) /0.531d0,0.523d0,0.521d0,0.534d0,0.557d0,0.567d0,0.579d0,0.601d0,0.602d0,0.597d0/
  data(coefficientsIron  (3,i    ),i=1,10) /5.77d-02,6.15d-02,6.22d-02,6.02d-02,5.79d-02,5.65d-02,5.49d-02,5.10d-02,5.07d-02,5.22d-02/

contains

  function verner1996ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{atomicRecombinationRateRadiativeVerner1996} atomic radiative recombination class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(atomicRecombinationRateRadiativeVerner1996)                :: self
    type(inputParameters                           ), intent(inout) :: parameters

    self=atomicRecombinationRateRadiativeVerner1996()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function verner1996ConstructorParameters

  double precision function verner1996Rate(self,atomicNumber,ionizationState,temperature,level)
    !!{
    Computes the rate coefficient of radiative recombination (in units of cm$^3$ s$^{-1}$) at the specified {\normalfont \ttfamily temperature} for all ions
    of all elements from H through Zn (selected by the {\normalfont \ttfamily atomicNumber} and the {\normalfont \ttfamily ionizationState} \emph{of the recombined
    ion}) use of the following fits:
    \begin{itemize}
    \item H-like, He-like, Li-like, Na-like: \citep{verner_atomic_1996};
    \item Other ions of C, N, O, Ne: \citep{pequignot_total_1991}, refitted by Verner \& Ferland formula to ensure correct asymptotes;
    \item Fe XVII-XXIII: \citep{arnaud_iron_1992};
    \item Fe I-XV: refitted by Verner \& Ferland formula to ensure correct asymptotes;
    \item Other ions of Mg, Si, S, Ar, Ca, Fe, Ni: \citep{shull_ionization_1982};
    \item Other ions of Na, Al: \citep{landini_x-uv_1990};
    \item Other ions of F, P, Cl, K, Ti, Cr, Mn, Co (excluding Ti I-II, Cr I-IV, Mn I-V, Co I): \citep{landini_ion_1991};
    \item All other species: interpolations of the power-law fits.
    \end{itemize}
    Based on the \href{https://web.archive.org/web/20220313133758/https://www.pa.uky.edu/~verner/dima/rec/rrfit.f}{code} originally written by Dima Verner. The ionization state
    passed to this function should be that of the atom/ion post recombination.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (atomicRecombinationRateRadiativeVerner1996), intent(inout)           :: self
    integer                                                     , intent(in   )           :: atomicNumber     , ionizationState
    double precision                                            , intent(in   )           :: temperature
    type            (enumerationRecombinationCaseType          ), intent(in   ), optional :: level
    ! Coefficients in the fitting function of Fergusen & Ferland (1997).
    double precision                                            , parameter               :: ffFitA=-9.97652d0
    double precision                                            , parameter               :: ffFitB= 0.03506d0
    double precision                                            , parameter               :: ffFitC= 0.15861d0
    double precision                                            , parameter               :: ffFitD=-0.03762d0
    double precision                                            , parameter               :: ffFitE= 0.30113d0
    double precision                                            , parameter               :: ffFitF= 0.00762d0
    double precision                                            , parameter               :: ffFitG=-0.06397d0
    double precision                                            , parameter               :: ffFitH=-0.00023d0
    double precision                                            , parameter               :: ffFitI= 0.00127d0
    integer                                                                               :: electronNumber
    double precision                                                                      :: temperatureScaled, fitFactor       , &
         &                                                                                   logTemperature
    !$GLC attributes unused :: self
    !![
    <optionalArgument name="level" defaultsTo="recombinationCaseA" />
    !!]

    ! Set zero rate by default.
    verner1996Rate=0.0d0
    ! If temperature is unphysical, return.
    if (temperature <= 0.0d0) return
    ! Determine type of recombination coefficient required.
    select case (level_%ID)
    case (recombinationCaseA%ID)
       ! Ensure atomic number is in range.
       if (atomicNumber  < 1 .or. atomicNumber    > 30          ) call Error_Report('atomic number is out of range'  //{introspection:location})
       ! Compute number of electrons.
       electronNumber=atomicNumber-ionizationState+1
       ! Ensure electron number is in range.
       if (electronNumber < 1 .or. electronNumber > atomicNumber) call Error_Report('electron number is out of range'//{introspection:location})
       ! Compute rate using the relevant fitting function.
       if     (                                                  &
            &    electronNumber <=  3                            &
            &  .or.                                              &
            &    electronNumber == 11                            &
            &  .or.                                              &
            &   (atomicNumber   >   5 .and. atomicNumber   <  9) &
            &  .or.                                              &
            &    atomicNumber   == 10                            &
            &  .or.                                              &
            &   (atomicNumber   == 26 .and. electronNumber > 11) &
            & ) then
          temperatureScaled=+sqrt(                                                       &
               &                  +    temperature                                       &
               &                  /    fitCoefficientsNew(3,atomicNumber,electronNumber) &
               &                 )
          verner1996Rate   =+          fitCoefficientsNew(1,atomicNumber,electronNumber) &
               &            /(                                                           &
               &              +        temperatureScaled                                 &
               &              *(                                                         &
               &                +      temperatureScaled                                 &
               &                +1.0d0                                                   &
               &               )**(                                                      &
               &                   +1.0d0                                                &
               &                   -   fitCoefficientsNew(2,atomicNumber,electronNumber) &
               &                 )                                                       &
               &              *(                                                         &
               &                +1.0d0                                                   &
               &                +sqrt(                                                   &
               &                      +temperature                                       &
               &                      /fitCoefficientsNew(4,atomicNumber,electronNumber) &
               &                     )                                                   &
               &               )**(                                                      &
               &                   +1.0d0                                                &
               &                   +   fitCoefficientsNew(2,atomicNumber,electronNumber) &
               &                  )                                                      &
               &             )
       else
          temperatureScaled=temperature/1.0d4
          if (atomicNumber == 26 .and. electronNumber <= 13) then
             verner1996Rate=+                     coefficientsIron(1,electronNumber) &
                  &         /temperatureScaled**(                                    &
                  &                              +coefficientsIron(2,electronNumber) &
                  &                              +coefficientsIron(3,electronNumber) &
                  &                              *log10(temperatureScaled)           &
                  &                             )
          else
             verner1996Rate=fitCoefficients(1,atomicNumber,electronNumber)&
                  &/temperatureScaled **fitCoefficients(2,atomicNumber,electronNumber)
          end if
       end if
    case (recombinationCaseB%ID)
       ! Case B recombination coefficient was requested.
       select case (atomicNumber)
       case (1) ! Hydrogen
          select case (ionizationState)
          case (1) ! H⁺
             ! Fergusen & Ferland (1997).
             logTemperature=+log10(temperature)
             fitFactor     =+(                      &
                  &           +      ffFitA         &
                  &           +      logTemperature &
                  &           *(                    &
                  &             +    ffFitC         &
                  &             +    logTemperature &
                  &             *(                  &
                  &               +  ffFitE         &
                  &               +  logTemperature &
                  &               *(                &
                  &                 +ffFitG         &
                  &                 +logTemperature &
                  &                 *ffFitI         &
                  &                )                &
                  &              )                  &
                  &            )                    &
                  &          )                      &
                  &         /(                      &
                  &           +1.0d0                &
                  &           +      logTemperature &
                  &           *(                    &
                  &             +    ffFitB         &
                  &             +    logTemperature &
                  &             *(                  &
                  &               +  ffFitD         &
                  &               +  logTemperature &
                  &               *(                &
                  &                 +ffFitF         &
                  &                 +logTemperature &
                  &                 *ffFitH         &
                  &                )                &
                  &              )                  &
                  &            )                    &
                  &          )
             if (fitFactor < 0.0d0) then
                verner1996Rate=+(10.0d0**fitFactor) &
                     &         /temperature
             else
                verner1996Rate=0.0d0
             end if
          case default
             call Error_Report('ionization state invalid for hydrogen'//{introspection:location})
          end select
       case (2) ! Helium
          select case (ionizationState)
          case (1) ! HeII
             ! Fit from Aparna Venkatesan.
             verner1996Rate=+4.3d-13       &
                  &         /(             &
                  &           +temperature &
                  &           /1.0d4       &
                  &         )**0.672d0
          case (2) ! HeIII
             ! Fit from Aparna Venkatesan, via Mike Shull. This is a fit to the data in Table 1 of Storey & Hummer (1995; MNRAS; 272; 41) for Z=2 and Ne=100.
             if (temperature < 1.0d6) then
                logTemperature=log10(temperature)
                verner1996Rate=+2.06d-11            &
                     &         *4.00d+00            &
                     &         /sqrt(temperature)   &
                     &         *(                   &
                     &           +6.572345200d0     &
                     &           -1.361055000d0     &
                     &           *logTemperature    &
                     &           +0.045686398d0     &
                     &           *logTemperature**2 &
                     &          )
             else
                ! Fit behaves badly above 10⁶K, so we fix the logarithmic slope of the polynomial part at its value at 10⁵K to
                ! avoid negative recombination rates at higher temperatures.
                logTemperature=+5.00d+00
                verner1996Rate=+2.06d-11            &
                     &         *4.00d+00            &
                     &         /sqrt(temperature)   &
                     &         *(                   &
                     &           +6.572345200d0     &
                     &           -1.361055000d0     &
                     &           *logTemperature    &
                     &           +0.045686398d0     &
                     &           *logTemperature**2 &
                     &          )                   &
                     &         /(                   &
                     &           +temperature       &
                     &           /1.0d5             &
                     &          )**0.392685d0
             end if
          case default
             call Error_Report('ionization state invalid for helium'//{introspection:location})
          end select
       case default
          call Error_Report('case B coefficients unavailable for requested atomic number'//{introspection:location})
       end select
    case default
       ! Recombination coefficient for an individual level was requested. We can not compute it, so report an error.
       call Error_Report('coefficients for individual levels are not available'//{introspection:location})
    end select
    return
  end function verner1996Rate
