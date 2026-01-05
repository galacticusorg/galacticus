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
Contains a module of useful astronomical constants.
!!}

module Numerical_Constants_Astronomical
  !!{
  Contains various useful astronomical constants.
  !!}
  use :: Numerical_Constants_Atomic  , only : atomicMassHelium     , atomicMassHydrogen, atomicMassLithium7
  use :: Numerical_Constants_Math    , only : Pi
  use :: Numerical_Constants_Physical, only : gravitationalConstant
  use :: Numerical_Constants_Prefixes, only : giga                 , hecto             , kilo              , mega
  use :: Numerical_Constants_Units   , only : ergs
  implicit none
  public

  ! Physical properties of the Sun.
  !![
  <constant variable="massSolar" group="astrophysical" gslSymbol="GSL_CONST_MKSA_SOLAR_MASS" gslHeader="gsl_const_mksa" symbol="\mathrm{M}_\odot" description="Solar mass." units="kg" unitsInSI="1.0" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_SOLAR_MASS"/>
  <constant variable="radiusSolar" group="astrophysical" value="6.95508d8" symbol="\mathrm{R}_\odot" description="Solar radius." units="m" unitsInSI="1.0" reference="Allen's Astrophysical Quantities, page 340" referenceURL="https://link.springer.com/book/10.1007/978-1-4612-1186-0"/>
  <constant variable="luminositySolar" group="astrophysical" value="3.845d26" symbol="\mathrm{L}_\odot" description="Solar luminosity." units="W" unitsInSI="1.0" reference="Allen's Astrophysical Quantities, page 340" referenceURL="https://link.springer.com/book/10.1007/978-1-4612-1186-0"/>
  !!]

  ! Solar composition.
  !![
  <constant variable="hydrogenByMassSolar" value="0.7070d0"  group="astrophysical" symbol="\mathrm{X}_\odot" description="Solar hydrogen fraction by mass." units="dimensionless" unitsInSI="1.0" reference="Allen's Astrophysical Quantities, page 28" referenceURL="https://link.springer.com/book/10.1007/978-1-4612-1186-0"/>
  <constant variable="heliumByMassSolar"   value="0.2740d0"  group="astrophysical" symbol="\mathrm{Y}_\odot" description="Solar helium fraction by mass."   units="dimensionless" unitsInSI="1.0" reference="Allen's Astrophysical Quantities, page 28" referenceURL="https://link.springer.com/book/10.1007/978-1-4612-1186-0"/>
  <constant variable="metallicitySolar"    value="0.0188d0"  group="astrophysical" symbol="\mathrm{Z}_\odot" description="Solar metal fraction by mass."    units="dimensionless" unitsInSI="1.0" reference="Allen's Astrophysical Quantities, page 28" referenceURL="https://link.springer.com/book/10.1007/978-1-4612-1186-0"/>
  <constant variable="heliumToHydrogenAbundanceSolar" value="(heliumByMassSolar/atomicMassHelium)/(hydrogenByMassSolar/atomicMassHydrogen)"  group="astrophysical" symbol="\mathrm{He/H}_\odot" description="Solar helium to hydrogen ratio by number." units="dimensionless" unitsInSI="1.0" reference="Derived."/>
  !!]

  ! Primordial composition.
  !![
  <constant variable="hydrogenByMassPrimordial" value="0.7514d+00"  group="astrophysical" symbol="\mathrm{X}_\mathrm{p}" description="Primordial hydrogen fraction by mass." units="dimensionless" unitsInSI="1.0" reference="\cite{cyburt_update_2008}"/>
  <constant variable="heliumByMassPrimordial"   value="0.2486d+00"  group="astrophysical" symbol="\mathrm{Y}_\mathrm{p}" description="Primordial helium fraction by mass. (Theoretical expectation based on WMAP results.)" units="dimensionless" unitsInSI="1.0" reference="\cite{cyburt_update_2008}"/>
  <constant variable="lithiumToHydrogenAbundancePrimordial" value="5.2400d-10"  group="astrophysical" symbol="\mathrm{Z}_\mathrm{Li,p}" description="Primordial lithium fraction by mass." units="dimensionless" unitsInSI="1.0" reference="\cite{cyburt_update_2008}"/>
  <constant variable="heliumToHydrogenAbundancePrimordial" value="(heliumByMassPrimordial/atomicMassHelium)/(hydrogenByMassPrimordial/atomicMassHydrogen)"  group="astrophysical" symbol="\mathrm{He/H}_\mathrm{p}" description="Primordial helium to hydrogen ratio by number." units="dimensionless" unitsInSI="1.0" reference="Derived."/>
  <constant variable="metallicityPrimordial" value="lithiumToHydrogenAbundancePrimordial*atomicMassLithium7/atomicMassHydrogen"  group="astrophysical" symbol="\mathrm{Z}_\mathrm{p}" description="Primordial metallicity." units="dimensionless" unitsInSI="1.0" reference="Derived."/>
  <constant variable="meanAtomicMassPrimordial" value="1.0d0/(2.0d0*hydrogenByMassPrimordial/atomicMassHydrogen+3.0d0*heliumByMassPrimordial/atomicMassHelium)"  group="astrophysical" symbol="\mu_\mathrm{p}" description="Primordial mean atomic mass." units="dimensionless" unitsInSI="1.0" reference="Derived."/>
  !!]

  ! Astrophysical scales.
  !![
  <constant variable="parsec" gslSymbol="GSL_CONST_MKSA_PARSEC" gslHeader="gsl_const_mksa"  group="units:astrophysical" symbol="\mathrm{pc}" description="Parsec---the distance at which 1 AU subtends an angle of one arcsecond." externalDescription="https://en.wikipedia.org/wiki/Parsec" units="m" unitsInSI="1.0" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_PARSEC"/>
  <constant variable="kiloParsec" value="kilo*parsec"  group="units:astrophysical" symbol="\mathrm{kpc}" description="Kiloparsec---$10^3$ parsecs." externalDescription="https://en.wikipedia.org/wiki/Parsec" units="m" unitsInSI="1.0" reference="Derived."/>
  <constant variable="megaParsec" value="mega*parsec"  group="units:astrophysical" symbol="\mathrm{Mpc}" description="Megaparsec---$10^6$ parsecs." externalDescription="https://en.wikipedia.org/wiki/Parsec" units="m" unitsInSI="1.0" reference="Derived."/>
  <constant variable="year" value="3.15581497635456d7"  group="units" symbol="\mathrm{yr}" description="Sidereal year." externalDescription="https://en.wikipedia.org/wiki/Sidereal_year" units="s" unitsInSI="1.0" reference="Earth observation center." referenceURL="https://hpiers.obspm.fr/eop-pc/models/constants.html"/>
  <constant variable="gigayear" value="giga*year"  group="units:astrophysical" symbol="\mathrm{Gyr}" description="Gigayear---$10^9$ years." externalDescription="https://en.wikipedia.org/wiki/Sidereal_year" units="s" unitsInSI="1.0" reference="Derived."/>
  !!]

  ! Newton's gravitational constant (in Galacticus' M_Solar, Mpc, km/s unit system).
  !![
  <constant variable="gravitationalConstant_internal" value="gravitationalConstant*massSolar/(kilo**2)/megaParsec" group="units:physical" description="Newton's gravitational constant in Galacticus' $\mathrm{M}_\odot$, Mpc, km/s unit system." units="km$^2$s$^{-2}$Mpc/$\mathrm{M}_\odot$" unitsInSI="6.456507" reference="Derived."/>
  !!]

  ! Conversion from Mpc/(km/s) to Gyr.
  !![
  <constant variable="MpcPerKmPerSToGyr" value="megaParsec/kilo/gigaYear" group="units:astrophysical" description="The conversion from Mpc/(km s$^{-1})$ to Gyr." units="Mpc/km s$^{-1}$/Gyr" unitsInSI="1.022712e-3" reference="Derived."/>
  !!]

  ! Conversion from magnitudes to optical depth.
  !![
  <constant variable="opticalDepthToMagnitudes" value="2.5d0/log(10.0d0)" group="astrophysical" description="The conversion factor between magnitudes of extinction and optical depth." units="dimensionless" unitsInSI="1.0" reference="Derived."/>
  !!]

  ! AB magnitude system.
  !![
  <constant variable="offsetAB" value="48.57d0" group="astrophysical" description="The zero point offset in the AB magnitude system: $m = -2.5\log_{10}(F_\nu/\mathrm{[ergs/s/cm}^2\mathrm{/Hz]})-48.57$. This zero point is chosen such that AB and Vega magnitude systems agree in the V-band. The original definition of the AB system zero point used a value of $48.60$ \citep{oke_secondary_1983} to achieve this. However, more \href{http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Vega}{recent measurements} give Vega a magnitude of $+0.03$ in the Vega magnitude system, which is why the AB magnitude zero point is now also shifted by $0.03$ relative to the original definition." externalDescription="https://en.wikipedia.org/wiki/AB_magnitude" units="dimensionless" unitsInSI="1.0" reference="\cite{oke_secondary_1983}"/>
  <constant variable="luminosityZeroPointAB" value="(10.0d0**(-offsetAB/2.5d0))*4.0d0*Pi*((10.0d0*parsec*hecto)**2)*ergs" group="astrophysical" description="The luminosity zero point offset in the AB magnitude system, found by computing the flux of a zeroth magnitude source at 10pc." externalDescription="https://en.wikipedia.org/wiki/AB_magnitude" units="W/Hz" unitsInSI="1.0" reference="\cite{oke_secondary_1983}"/>
  !!]

  ! Angular conversions.
  !![
  <constant variable="arcminutesToDegrees" value="1.0d0/60.0d0" group="units:astrophysical" description="Conversion factor from arcminutes to degrees." units="$^\circ/^\prime$" unitsInSI="1.0" reference="Defined."/>
  <constant variable="arcsecondsToDegrees" value="1.0d0/3600.0d0" group="units:astrophysical" description="Conversion factor from arcseconds to degrees." units="$^\circ/^{\prime\prime}$" unitsInSI="1.0" reference="Defined."/>
  <constant variable="hoursToDegrees" value="360.0d0/24.0d0" group="units:astrophysical" description="Conversion factor from hours to degrees." units="$^\circ/^{\prime}$" unitsInSI="1.0" reference="Defined."/>
  <constant variable="minutesToDegrees" value="360.0d0/24.0d0/60.0d0" group="units:astrophysical" description="Conversion factor from minutes to degrees." units="$^\circ/^{\prime\prime}$" unitsInSI="1.0" reference="Defined."/>
  <constant variable="secondsToDegrees" value="360.0d0/24.0d0/3600.0d0" group="units:astrophysical" description="Conversion factor from seconds to degrees." units="$^\circ/^{\prime\prime}$" unitsInSI="1.0" reference="Defined."/>
  <constant variable="hoursToRadians" value="2.0d0*Pi/24.0d0" group="units:astrophysical" description="Conversion factor from hours to radians." units="$\mathrm{rad}/^{\prime}$" unitsInSI="1.0" reference="Defined."/>
  <constant variable="degreesToRadians" value="2.0d0*Pi/360.0d0" group="units:astrophysical" description="Conversion factor from degrees to radians." units="$\mathrm{rad}/^\circ$" unitsInSI="1.0" reference="Defined."/>
  !!]

end module Numerical_Constants_Astronomical
