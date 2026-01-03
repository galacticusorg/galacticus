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
Contains a module of useful mathematical constants.
!!}

module Numerical_Constants_Math
  !!{
  Contains various useful mathematical constants.
  !!}
  use :: Kind_Numbers, only : kind_quad
  implicit none

  !![
  <constant variable="e" gslSymbol="M_E" gslHeader="gsl_math" symbol="\mathrm{e}" units="dimensionless" unitsInSI="1.0" description="Euler's number---the base of natural logarithms." externalDescription="https://en.wikipedia.org/wiki/E_(mathematical_constant)" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/math.html#mathematical-constants" group="math"/>
  <constant variable="Pi" gslSymbol="M_PI" gslHeader="gsl_math" symbol="\pi" units="dimensionless" unitsInSI="1.0" description="The ratio of a circle's perimeter to its diameter." externalDescription="https://en.wikipedia.org/wiki/Pi" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/math.html#mathematical-constants" group="math"/>
  <constant variable="PiQuadPrecision" value="3.141592653589793238462643383279502884197_kind_quad" type="real(kind_quad)" symbol="\pi" units="dimensionless" unitsInSI="1.0" description="The ratio of a circle's perimeter to its diameter (to quadruple precision)." externalDescription="https://en.wikipedia.org/wiki/Pi" reference="OEIS: A002117" referenceURL="https://oeis.org/A002117"  group="math"/>
  <constant variable="ln2" gslSymbol="M_LN2" gslHeader="gsl_math" symbol="\log 2" units="dimensionless" unitsInSI="1.0" description="The natural logarithm of 2." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/math.html#mathematical-constants" group="math"/>
  <constant variable="ln10" gslSymbol="M_LN10" gslHeader="gsl_math" symbol="\log 10" units="dimensionless" unitsInSI="1.0" description="The natural logarithm of 10." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/math.html#mathematical-constants" group="math"/>
  <constant variable="eulersConstant" gslSymbol="M_EULER" gslHeader="gsl_math" symbol="\gamma" units="dimensionless" unitsInSI="1.0" description="Euler's constant." externalDescription="https://en.wikipedia.org/wiki/Euler%27s_constant" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/math.html#mathematical-constants" group="math"/>
  <constant variable="riemannZeta3" value="1.20205690315959428539973816151144999076d0" symbol="\zeta(3)" units="dimensionless" unitsInSI="1.0" description="Riemann zeta function evaluated at $s=3$." externalDescription="https://en.wikipedia.org/wiki/Riemann_zeta_function" reference="OEIS: A002117" referenceURL="https://oeis.org/A002117" group="math"/>
  <constant variable="catalan" value="0.91596559417721901505460351493238411077d0" symbol="G" units="dimensionless" unitsInSI="1.0" description="Catalan's constant." externalDescription="https://en.wikipedia.org/wiki/Catalan%27s_constant" reference="OEIS: A006752" referenceURL="https://oeis.org/A006752" group="math"/>
  !!]

end module Numerical_Constants_Math
