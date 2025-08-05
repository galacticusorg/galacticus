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
Implements an atomic ionization potential class, which provides potentials for all ionization stages of all atoms from H to Zn
using data taken from Dima Verner's \href{https://web.archive.org/web/20220313133801/https://www.pa.uky.edu/~verner/dima/col/cfit.f}{code}.
!!}

  !![
  <atomicIonizationPotential name="atomicIonizationPotentialVerner">
   <description>
    Implements an atomic ionization potential class, which provides potentials for all ionization stages of all atoms from H to Zn using data taken from Dima Verner's \href{https://web.archive.org/web/20220313133801/https://www.pa.uky.edu/~verner/dima/col/cfit.f}{code}.
   </description>
  </atomicIonizationPotential>
  !!]
  type, extends(atomicIonizationPotentialClass) :: atomicIonizationPotentialVerner
     !!{
     Implements an atomic ionization potential class, which provides potentials for all ionization stages of all atoms from H
     to Zn using data taken from Dima Verner's \href{https://web.archive.org/web/20220313133801/https://www.pa.uky.edu/~verner/dima/col/cfit.f}{code}.
     !!}
     private
   contains
     procedure :: potential => vernerPotential
  end type atomicIonizationPotentialVerner

  interface atomicIonizationPotentialVerner
     !!{
     Constructors for the \refClass{atomicIonizationPotentialVerner} atomic ionization potential class.
     !!}
     module procedure vernerConstructorParameters
  end interface atomicIonizationPotentialVerner

  ! Ionization potential data.
  double precision, dimension(83,83) :: potentialIonization
  data potentialIonization( 1, 1) /   13.59900d0/
  data potentialIonization( 2, 2) /   24.58800d0/
  data potentialIonization( 2, 1) /   54.41800d0/
  data potentialIonization( 3, 3) /    5.39200d0/
  data potentialIonization( 3, 2) /   75.64100d0/
  data potentialIonization( 3, 1) /  122.45500d0/
  data potentialIonization( 4, 4) /    9.32300d0/
  data potentialIonization( 4, 3) /   18.21100d0/
  data potentialIonization( 4, 2) /  153.89600d0/
  data potentialIonization( 4, 1) /  217.72000d0/
  data potentialIonization( 5, 5) /    8.29800d0/
  data potentialIonization( 5, 4) /   25.15500d0/
  data potentialIonization( 5, 3) /   37.93100d0/
  data potentialIonization( 5, 2) /  259.37700d0/
  data potentialIonization( 5, 1) /  340.22900d0/
  data potentialIonization( 6, 6) /   11.26000d0/
  data potentialIonization( 6, 5) /   24.38400d0/
  data potentialIonization( 6, 4) /   47.88800d0/
  data potentialIonization( 6, 3) /   64.49400d0/
  data potentialIonization( 6, 2) /  392.09000d0/
  data potentialIonization( 6, 1) /  489.99700d0/
  data potentialIonization( 7, 7) /   14.53400d0/
  data potentialIonization( 7, 6) /   29.60200d0/
  data potentialIonization( 7, 5) /   47.45000d0/
  data potentialIonization( 7, 4) /   77.47400d0/
  data potentialIonization( 7, 3) /   97.89100d0/
  data potentialIonization( 7, 2) /  552.07400d0/
  data potentialIonization( 7, 1) /  667.05100d0/
  data potentialIonization( 8, 8) /   13.61800d0/
  data potentialIonization( 8, 7) /   35.11800d0/
  data potentialIonization( 8, 6) /   54.93600d0/
  data potentialIonization( 8, 5) /   77.41400d0/
  data potentialIonization( 8, 4) /  113.90000d0/
  data potentialIonization( 8, 3) /  138.12100d0/
  data potentialIonization( 8, 2) /  739.33800d0/
  data potentialIonization( 8, 1) /  871.41700d0/
  data potentialIonization( 9, 9) /   17.42300d0/
  data potentialIonization( 9, 8) /   34.97100d0/
  data potentialIonization( 9, 7) /   62.70900d0/
  data potentialIonization( 9, 6) /   87.14100d0/
  data potentialIonization( 9, 5) /  114.24000d0/
  data potentialIonization( 9, 4) /  157.16600d0/
  data potentialIonization( 9, 3) /  185.18800d0/
  data potentialIonization( 9, 2) /  953.91500d0/
  data potentialIonization( 9, 1) / 1103.12600d0/
  data potentialIonization(10,10) /   21.56500d0/
  data potentialIonization(10, 9) /   40.96400d0/
  data potentialIonization(10, 8) /   63.46000d0/
  data potentialIonization(10, 7) /   97.12000d0/
  data potentialIonization(10, 6) /  126.22000d0/
  data potentialIonization(10, 5) /  157.93000d0/
  data potentialIonization(10, 4) /  207.28000d0/
  data potentialIonization(10, 3) /  239.10000d0/
  data potentialIonization(10, 2) / 1195.83000d0/
  data potentialIonization(10, 1) / 1362.20000d0/
  data potentialIonization(11,11) /    5.13900d0/
  data potentialIonization(11,10) /   47.28700d0/
  data potentialIonization(11, 9) /   71.62000d0/
  data potentialIonization(11, 8) /   98.92000d0/
  data potentialIonization(11, 7) /  138.39000d0/
  data potentialIonization(11, 6) /  172.15000d0/
  data potentialIonization(11, 5) /  208.48000d0/
  data potentialIonization(11, 4) /  264.19000d0/
  data potentialIonization(11, 3) /  299.88000d0/
  data potentialIonization(11, 2) / 1465.13000d0/
  data potentialIonization(11, 1) / 1648.71000d0/
  data potentialIonization(12,12) /    7.64600d0/
  data potentialIonization(12,11) /   15.03500d0/
  data potentialIonization(12,10) /   80.14400d0/
  data potentialIonization(12, 9) /  109.25000d0/
  data potentialIonization(12, 8) /  141.27000d0/
  data potentialIonization(12, 7) /  186.51000d0/
  data potentialIonization(12, 6) /  224.95000d0/
  data potentialIonization(12, 5) /  265.96000d0/
  data potentialIonization(12, 4) /  328.24000d0/
  data potentialIonization(12, 3) /  367.54000d0/
  data potentialIonization(12, 2) / 1761.85000d0/
  data potentialIonization(12, 1) / 1962.68000d0/
  data potentialIonization(13,13) /    5.98600d0/
  data potentialIonization(13,12) /   18.82900d0/
  data potentialIonization(13,11) /   28.44800d0/
  data potentialIonization(13,10) /  119.99000d0/
  data potentialIonization(13, 9) /  153.82600d0/
  data potentialIonization(13, 8) /  190.48000d0/
  data potentialIonization(13, 7) /  241.44000d0/
  data potentialIonization(13, 6) /  284.60000d0/
  data potentialIonization(13, 5) /  330.11000d0/
  data potentialIonization(13, 4) /  399.37000d0/
  data potentialIonization(13, 3) /  442.08000d0/
  data potentialIonization(13, 2) / 2086.04000d0/
  data potentialIonization(13, 1) / 2304.16000d0/
  data potentialIonization(14,14) /    8.15200d0/
  data potentialIonization(14,13) /   16.34600d0/
  data potentialIonization(14,12) /   33.49300d0/
  data potentialIonization(14,11) /   45.14200d0/
  data potentialIonization(14,10) /  166.77000d0/
  data potentialIonization(14, 9) /  205.06000d0/
  data potentialIonization(14, 8) /  246.53000d0/
  data potentialIonization(14, 7) /  303.18000d0/
  data potentialIonization(14, 6) /  351.11000d0/
  data potentialIonization(14, 5) /  401.38000d0/
  data potentialIonization(14, 4) /  476.08000d0/
  data potentialIonization(14, 3) /  523.52000d0/
  data potentialIonization(14, 2) / 2437.74000d0/
  data potentialIonization(14, 1) / 2673.20000d0/
  data potentialIonization(15,15) /   10.49000d0/
  data potentialIonization(15,14) /   19.73000d0/
  data potentialIonization(15,13) /   30.20400d0/
  data potentialIonization(15,12) /   51.44400d0/
  data potentialIonization(15,11) /   65.02600d0/
  data potentialIonization(15,10) /  220.43400d0/
  data potentialIonization(15, 9) /  263.23000d0/
  data potentialIonization(15, 8) /  309.42000d0/
  data potentialIonization(15, 7) /  371.74000d0/
  data potentialIonization(15, 6) /  424.51000d0/
  data potentialIonization(15, 5) /  479.59000d0/
  data potentialIonization(15, 4) /  560.43000d0/
  data potentialIonization(15, 3) /  611.90000d0/
  data potentialIonization(15, 2) / 2817.02000d0/
  data potentialIonization(15, 1) / 3069.86000d0/
  data potentialIonization(16,16) /   10.36000d0/
  data potentialIonization(16,15) /   23.33000d0/
  data potentialIonization(16,14) /   34.83000d0/
  data potentialIonization(16,13) /   47.30500d0/
  data potentialIonization(16,12) /   72.68000d0/
  data potentialIonization(16,11) /   88.05400d0/
  data potentialIonization(16,10) /  280.94000d0/
  data potentialIonization(16, 9) /  328.24000d0/
  data potentialIonization(16, 8) /  379.11000d0/
  data potentialIonization(16, 7) /  447.10000d0/
  data potentialIonization(16, 6) /  504.79000d0/
  data potentialIonization(16, 5) /  564.67000d0/
  data potentialIonization(16, 4) /  651.65000d0/
  data potentialIonization(16, 3) /  707.16000d0/
  data potentialIonization(16, 2) / 3223.92000d0/
  data potentialIonization(16, 1) / 3494.21000d0/
  data potentialIonization(17,17) /   12.96800d0/
  data potentialIonization(17,16) /   23.81400d0/
  data potentialIonization(17,15) /   39.61000d0/
  data potentialIonization(17,14) /   53.46600d0/
  data potentialIonization(17,13) /   67.80000d0/
  data potentialIonization(17,12) /   97.03000d0/
  data potentialIonization(17,11) /  114.19700d0/
  data potentialIonization(17,10) /  348.29000d0/
  data potentialIonization(17, 9) /  400.06000d0/
  data potentialIonization(17, 8) /  455.63000d0/
  data potentialIonization(17, 7) /  529.28000d0/
  data potentialIonization(17, 6) /  591.99000d0/
  data potentialIonization(17, 5) /  656.71000d0/
  data potentialIonization(17, 4) /  749.76000d0/
  data potentialIonization(17, 3) /  809.41000d0/
  data potentialIonization(17, 2) / 3658.51000d0/
  data potentialIonization(17, 1) / 3946.32000d0/
  data potentialIonization(18,18) /   15.76000d0/
  data potentialIonization(18,17) /   27.63000d0/
  data potentialIonization(18,16) /   40.74000d0/
  data potentialIonization(18,15) /   59.81000d0/
  data potentialIonization(18,14) /   75.02000d0/
  data potentialIonization(18,13) /   91.01000d0/
  data potentialIonization(18,12) /  124.32400d0/
  data potentialIonization(18,11) /  143.45800d0/
  data potentialIonization(18,10) /  422.45000d0/
  data potentialIonization(18, 9) /  478.69000d0/
  data potentialIonization(18, 8) /  538.96000d0/
  data potentialIonization(18, 7) /  618.26000d0/
  data potentialIonization(18, 6) /  686.11000d0/
  data potentialIonization(18, 5) /  755.75000d0/
  data potentialIonization(18, 4) /  854.78000d0/
  data potentialIonization(18, 3) /  918.04000d0/
  data potentialIonization(18, 2) / 4120.87000d0/
  data potentialIonization(18, 1) / 4426.24000d0/
  data potentialIonization(19,19) /    4.34069d0/
  data potentialIonization(19,18) /   31.63000d0/
  data potentialIonization(19,17) /   45.80600d0/
  data potentialIonization(19,16) /   60.91000d0/
  data potentialIonization(19,15) /   82.66000d0/
  data potentialIonization(19,14) /   99.40000d0/
  data potentialIonization(19,13) /  117.56000d0/
  data potentialIonization(19,12) /  154.71000d0/
  data potentialIonization(19,11) /  175.81800d0/
  data potentialIonization(19,10) /  503.80000d0/
  data potentialIonization(19, 9) /  564.70000d0/
  data potentialIonization(19, 8) /  629.40000d0/
  data potentialIonization(19, 7) /  714.60000d0/
  data potentialIonization(19, 6) /  786.60000d0/
  data potentialIonization(19, 5) /  861.60000d0/
  data potentialIonization(19, 4) /  968.00000d0/
  data potentialIonization(19, 3) / 1033.40000d0/
  data potentialIonization(19, 2) / 4611.06000d0/
  data potentialIonization(19, 1) / 4934.07000d0/
  data potentialIonization(20,20) /    6.11321d0/
  data potentialIonization(20,19) /   11.87181d0/
  data potentialIonization(20,18) /   50.91350d0/
  data potentialIonization(20,17) /   67.27000d0/
  data potentialIonization(20,16) /   84.50000d0/
  data potentialIonization(20,15) /  108.78000d0/
  data potentialIonization(20,14) /  127.20000d0/
  data potentialIonization(20,13) /  147.24000d0/
  data potentialIonization(20,12) /  188.30000d0/
  data potentialIonization(20,11) /  211.27700d0/
  data potentialIonization(20,10) /  591.90000d0/
  data potentialIonization(20, 9) /  657.20000d0/
  data potentialIonization(20, 8) /  726.60000d0/
  data potentialIonization(20, 7) /  817.60000d0/
  data potentialIonization(20, 6) /  894.50000d0/
  data potentialIonization(20, 5) /  974.00000d0/
  data potentialIonization(20, 4) / 1087.00000d0/
  data potentialIonization(20, 3) / 1157.00000d0/
  data potentialIonization(20, 2) / 5129.16000d0/
  data potentialIonization(20, 1) / 5469.88000d0/
  data potentialIonization(21,21) /    6.56154d0/
  data potentialIonization(21,20) /   12.79987d0/
  data potentialIonization(21,19) /   24.75704d0/
  data potentialIonization(21,18) /   73.49000d0/
  data potentialIonization(21,17) /   91.90000d0/
  data potentialIonization(21,15) /  138.00000d0/
  data potentialIonization(21,14) /  158.10000d0/
  data potentialIonization(21,13) /  180.02000d0/
  data potentialIonization(21,12) /  225.10000d0/
  data potentialIonization(21,11) /  249.83700d0/
  data potentialIonization(21,10) /  687.36000d0/
  data potentialIonization(21, 9) /  756.70000d0/
  data potentialIonization(21, 8) /  830.80000d0/
  data potentialIonization(21, 7) /  927.50000d0/
  data potentialIonization(21, 6) / 1009.00000d0/
  data potentialIonization(21, 5) / 1094.00000d0/
  data potentialIonization(21, 4) / 1213.00000d0/
  data potentialIonization(21, 3) / 1288.00000d0/
  data potentialIonization(21, 2) / 5675.26000d0/
  data potentialIonization(21, 1) / 6033.76900d0/
  data potentialIonization(22,22) /    6.82000d0/
  data potentialIonization(22,21) /   13.58000d0/
  data potentialIonization(22,20) /   27.49190d0/
  data potentialIonization(22,19) /   43.26750d0/
  data potentialIonization(22,18) /   99.30000d0/
  data potentialIonization(22,17) /  119.53000d0/
  data potentialIonization(22,15) /  170.40000d0/
  data potentialIonization(22,14) /  192.10000d0/
  data potentialIonization(22,13) /  215.90000d0/
  data potentialIonization(22,12) /  265.00000d0/
  data potentialIonization(22,11) /  291.50200d0/
  data potentialIonization(22,10) /  787.84000d0/
  data potentialIonization(22, 9) /  863.10000d0/
  data potentialIonization(22, 8) /  941.90000d0/
  data potentialIonization(22, 7) / 1044.00000d0/
  data potentialIonization(22, 6) / 1131.00000d0/
  data potentialIonization(22, 5) / 1221.00000d0/
  data potentialIonization(22, 4) / 1346.00000d0/
  data potentialIonization(22, 3) / 1425.30000d0/
  data potentialIonization(22, 2) / 6249.20000d0/
  data potentialIonization(22, 1) / 6625.82000d0/
  data potentialIonization(23,23) /    6.74000d0/
  data potentialIonization(23,22) /   14.66000d0/
  data potentialIonization(23,21) /   29.31100d0/
  data potentialIonization(23,20) /   46.70900d0/
  data potentialIonization(23,19) /   65.28220d0/
  data potentialIonization(23,18) /  128.13000d0/
  data potentialIonization(23,17) /  150.60000d0/
  data potentialIonization(23,16) /  173.40000d0/
  data potentialIonization(23,15) /  205.80000d0/
  data potentialIonization(23,14) /  230.50000d0/
  data potentialIonization(23,13) /  255.70000d0/
  data potentialIonization(23,12) /  308.10000d0/
  data potentialIonization(23,11) /  336.28000d0/
  data potentialIonization(23,10) /  896.00000d0/
  data potentialIonization(23, 9) /  976.00000d0/
  data potentialIonization(23, 8) / 1060.00000d0/
  data potentialIonization(23, 7) / 1168.00000d0/
  data potentialIonization(23, 6) / 1260.00000d0/
  data potentialIonization(23, 5) / 1355.00000d0/
  data potentialIonization(23, 4) / 1486.00000d0/
  data potentialIonization(23, 3) / 1569.60000d0/
  data potentialIonization(23, 2) / 6851.77000d0/
  data potentialIonization(23, 1) / 7246.13000d0/
  data potentialIonization(24,24) /    6.76660d0/
  data potentialIonization(24,23) /   16.49750d0/
  data potentialIonization(24,22) /   30.96000d0/
  data potentialIonization(24,20) /   69.46000d0/
  data potentialIonization(24,19) /   90.63600d0/
  data potentialIonization(24,18) /  160.18000d0/
  data potentialIonization(24,17) /  184.70000d0/
  data potentialIonization(24,16) /  209.30000d0/
  data potentialIonization(24,15) /  244.40000d0/
  data potentialIonization(24,14) /  270.80000d0/
  data potentialIonization(24,13) /  298.00000d0/
  data potentialIonization(24,12) /  354.80000d0/
  data potentialIonization(24,11) /  384.17000d0/
  data potentialIonization(24,10) / 1010.60000d0/
  data potentialIonization(24, 9) / 1097.00000d0/
  data potentialIonization(24, 8) / 1185.00000d0/
  data potentialIonization(24, 7) / 1299.00000d0/
  data potentialIonization(24, 6) / 1396.00000d0/
  data potentialIonization(24, 5) / 1496.00000d0/
  data potentialIonization(24, 4) / 1634.00000d0/
  data potentialIonization(24, 3) / 1721.60000d0/
  data potentialIonization(24, 2) / 7482.40000d0/
  data potentialIonization(24, 1) / 7894.79000d0/
  data potentialIonization(25,25) /    7.43410d0/
  data potentialIonization(25,24) /   15.64011d0/
  data potentialIonization(25,20) /   95.74800d0/
  data potentialIonization(25,19) /  119.27000d0/
  data potentialIonization(25,18) /  194.50000d0/
  data potentialIonization(25,17) /  221.80000d0/
  data potentialIonization(25,16) /  248.30000d0/
  data potentialIonization(25,15) /  286.00000d0/
  data potentialIonization(25,14) /  314.40000d0/
  data potentialIonization(25,13) /  343.60000d0/
  data potentialIonization(25,12) /  403.00000d0/
  data potentialIonization(25,11) /  435.16000d0/
  data potentialIonization(25,10) / 1133.10000d0/
  data potentialIonization(25, 9) / 1224.00000d0/
  data potentialIonization(25, 8) / 1317.00000d0/
  data potentialIonization(25, 7) / 1437.00000d0/
  data potentialIonization(25, 6) / 1539.00000d0/
  data potentialIonization(25, 5) / 1644.00000d0/
  data potentialIonization(25, 4) / 1788.00000d0/
  data potentialIonization(25, 3) / 1879.90000d0/
  data potentialIonization(25, 2) / 8141.35000d0/
  data potentialIonization(25, 1) / 8571.94000d0/
  data potentialIonization(26,26) /    7.87000d0/
  data potentialIonization(26,25) /   16.18800d0/
  data potentialIonization(26,24) /   30.65200d0/
  data potentialIonization(26,23) /   54.80000d0/
  data potentialIonization(26,22) /   75.00000d0/
  data potentialIonization(26,21) /   99.10000d0/
  data potentialIonization(26,20) /  125.00000d0/
  data potentialIonization(26,19) /  151.06000d0/
  data potentialIonization(26,18) /  233.60000d0/
  data potentialIonization(26,17) /  262.10000d0/
  data potentialIonization(26,16) /  290.30000d0/
  data potentialIonization(26,15) /  330.80000d0/
  data potentialIonization(26,14) /  361.00000d0/
  data potentialIonization(26,13) /  392.20000d0/
  data potentialIonization(26,12) /  457.00000d0/
  data potentialIonization(26,11) /  489.27000d0/
  data potentialIonization(26,10) / 1262.00000d0/
  data potentialIonization(26, 9) / 1358.00000d0/
  data potentialIonization(26, 8) / 1456.00000d0/
  data potentialIonization(26, 7) / 1582.00000d0/
  data potentialIonization(26, 6) / 1689.00000d0/
  data potentialIonization(26, 5) / 1799.00000d0/
  data potentialIonization(26, 4) / 1950.00000d0/
  data potentialIonization(26, 3) / 2046.00000d0/
  data potentialIonization(26, 2) / 8828.80000d0/
  data potentialIonization(26, 1) / 9277.65000d0/
  data potentialIonization(27,27) /    7.86400d0/
  data potentialIonization(27,26) /   17.08300d0/
  data potentialIonization(27,25) /   33.50000d0/
  data potentialIonization(27,20) /  157.80000d0/
  data potentialIonization(27,19) /  186.13000d0/
  data potentialIonization(27,18) /  275.40000d0/
  data potentialIonization(27,17) /  305.30000d0/
  data potentialIonization(27,16) /  336.00000d0/
  data potentialIonization(27,15) /  379.00000d0/
  data potentialIonization(27,13) /  444.00000d0/
  data potentialIonization(27,12) /  511.96000d0/
  data potentialIonization(27,11) /  546.59000d0/
  data potentialIonization(27,10) / 1397.20000d0/
  data potentialIonization(27, 9) / 1504.60000d0/
  data potentialIonization(27, 8) / 1603.00000d0/
  data potentialIonization(27, 7) / 1735.00000d0/
  data potentialIonization(27, 6) / 1846.00000d0/
  data potentialIonization(27, 5) / 1962.00000d0/
  data potentialIonization(27, 4) / 2119.00000d0/
  data potentialIonization(27, 3) / 2219.00000d0/
  data potentialIonization(27, 2) / 9544.83000d0/
  data potentialIonization(27, 1) /10012.07000d0/
  data potentialIonization(28,28) /    7.63800d0/
  data potentialIonization(28,27) /   18.16898d0/
  data potentialIonization(28,20) /  193.00000d0/
  data potentialIonization(28,18) /  321.00000d0/
  data potentialIonization(28,17) /  352.00000d0/
  data potentialIonization(28,16) /  384.00000d0/
  data potentialIonization(28,15) /  430.00000d0/
  data potentialIonization(28,14) /  464.00000d0/
  data potentialIonization(28,13) /  499.00000d0/
  data potentialIonization(28,12) /  571.30000d0/
  data potentialIonization(28,11) /  607.06000d0/
  data potentialIonization(28,10) / 1541.00000d0/
  data potentialIonization(28, 9) / 1648.00000d0/
  data potentialIonization(28, 8) / 1756.00000d0/
  data potentialIonization(28, 7) / 1894.00000d0/
  data potentialIonization(28, 6) / 2011.00000d0/
  data potentialIonization(28, 5) / 2131.00000d0/
  data potentialIonization(28, 4) / 2295.00000d0/
  data potentialIonization(28, 3) / 2399.20000d0/
  data potentialIonization(28, 2) /10289.50000d0/
  data potentialIonization(28, 1) /10775.30000d0/
  data potentialIonization(29,29) /    7.72640d0/
  data potentialIonization(29,28) /   20.29210d0/
  data potentialIonization(29, 9) / 1793.00000d0/
  data potentialIonization(29, 8) / 1905.00000d0/
  data potentialIonization(29, 4) / 2459.00000d0/
  data potentialIonization(29, 3) / 2585.00000d0/
  data potentialIonization(29, 2) /11063.10000d0/
  data potentialIonization(29, 1) /11567.50000d0/
  data potentialIonization(30,30) /    9.39430d0/
  data potentialIonization(30,29) /   17.96450d0/
  data potentialIonization(30, 9) / 1953.00000d0/
  data potentialIonization(30, 8) / 2070.00000d0/
  data potentialIonization(30, 4) / 2647.00000d0/
  data potentialIonization(30, 3) / 2780.00000d0/
  data potentialIonization(30, 2) /11865.60000d0/
  data potentialIonization(30, 1) /12388.80000d0/
  data potentialIonization(31,31) /    5.99900d0/
  data potentialIonization(31,30) /   20.51000d0/
  data potentialIonization(31,29) /   30.71000d0/
  data potentialIonization(31, 9) / 2120.00000d0/
  data potentialIonization(31, 8) / 2242.00000d0/
  data potentialIonization(31, 3) / 2982.00000d0/
  data potentialIonization(31, 2) /12697.00000d0/
  data potentialIonization(31, 1) /13239.00000d0/
  data potentialIonization(32,32) /    7.89900d0/
  data potentialIonization(32,31) /   15.93400d0/
  data potentialIonization(32,30) /   34.22000d0/
  data potentialIonization(32,29) /   45.71000d0/
  data potentialIonization(32, 9) / 2294.00000d0/
  data potentialIonization(32, 8) / 2421.00000d0/
  data potentialIonization(32, 3) / 3192.00000d0/
  data potentialIonization(32, 2) /13558.00000d0/
  data potentialIonization(32, 1) /14119.00000d0/
  data potentialIonization(33,33) /    9.81000d0/
  data potentialIonization(33, 9) / 2474.00000d0/
  data potentialIonization(33, 8) / 2606.00000d0/
  data potentialIonization(33, 3) / 3406.00000d0/
  data potentialIonization(33, 2) /14448.00000d0/
  data potentialIonization(33, 1) /15029.00000d0/
  data potentialIonization(34, 9) / 2661.00000d0/
  data potentialIonization(34, 8) / 2798.00000d0/
  data potentialIonization(34, 3) / 3633.00000d0/
  data potentialIonization(34, 2) /14448.00000d0/
  data potentialIonization(34, 1) /15029.00000d0/
  data potentialIonization(35,35) /   11.81400d0/
  data potentialIonization(35, 9) / 2855.00000d0/
  data potentialIonization(35, 4) / 3684.00000d0/
  data potentialIonization(35, 3) / 3865.00000d0/
  data potentialIonization(35, 2) /16317.00000d0/
  data potentialIonization(35, 1) /16937.00000d0/
  data potentialIonization(36, 9) / 3056.00000d0/
  data potentialIonization(36, 4) / 3912.00000d0/
  data potentialIonization(36, 3) / 4105.00000d0/
  data potentialIonization(36, 2) /17296.00000d0/
  data potentialIonization(36, 1) /17936.00000d0/
  data potentialIonization(37,37) /    4.17700d0/
  data potentialIonization(38,38) /    5.69500d0/
  data potentialIonization(38,37) /   11.03000d0/
  data potentialIonization(39,39) /    6.38000d0/
  data potentialIonization(39,38) /   12.24000d0/
  data potentialIonization(42,42) /    7.09900d0/
  data potentialIonization(45,45) /    7.46000d0/
  data potentialIonization(47,47) /    7.57600d0/
  data potentialIonization(48,48) /    8.99300d0/
  data potentialIonization(48,47) /   16.90800d0/
  data potentialIonization(49,49) /    5.78600d0/
  data potentialIonization(50,50) /    7.34400d0/
  data potentialIonization(53,53) /   10.45700d0/
  data potentialIonization(54,54) /   12.13000d0/
  data potentialIonization(55,55) /    3.89400d0/
  data potentialIonization(56,56) /    5.21200d0/
  data potentialIonization(56,55) /   10.00400d0/
  data potentialIonization(60,59) /   10.72000d0/
  data potentialIonization(63,63) /    5.67000d0/
  data potentialIonization(69,69) /    6.18000d0/
  data potentialIonization(70,70) /    6.25400d0/
  data potentialIonization(70,69) /   12.17000d0/
  data potentialIonization(77,77) /    9.10000d0/
  data potentialIonization(79,79) /    9.22500d0/
  data potentialIonization(80,80) /   10.43700d0/
  data potentialIonization(81,81) /    6.10800d0/
  data potentialIonization(82,82) /    7.41600d0/
  data potentialIonization(83,83) /    7.28900d0/

contains

  function vernerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{atomicIonizationPotentialVerner} atomic ionization potential class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(atomicIonizationPotentialVerner)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=atomicIonizationPotentialVerner()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function vernerConstructorParameters

  double precision function vernerPotential(self,atomicNumber,electronNumber)
    !!{
    Return the ionization potential (in units of electron volts) for the ion with given {\normalfont \ttfamily atomicNumber} and
    {\normalfont \ttfamily electronNumber} using data taken from Dima Verner's
    \href{https://web.archive.org/web/20220313133801/https://www.pa.uky.edu/~verner/dima/col/cfit.f}{code}.
    !!}
    implicit none
    class  (atomicIonizationPotentialVerner), intent(inout) :: self
    integer                                 , intent(in   ) :: atomicNumber, electronNumber
    !$GLC attributes unused :: self

    vernerPotential=potentialIonization(atomicNumber,electronNumber)
    return
  end function vernerPotential
