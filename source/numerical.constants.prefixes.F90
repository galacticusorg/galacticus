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
Contains a module of useful numerical prefixes.
!!}

module Numerical_Constants_Prefixes
  !!{
  Contains useful numerical prefixes.
  !!}
  implicit none
  public
  
  !![
  <constant variable="quetta" value="1.0d+30" symbol="\mathrm{Q}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="ronna"  value="1.0d+27" symbol="\mathrm{R}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="yotta"  value="1.0d+24" symbol="\mathrm{Y}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="zetta"  value="1.0d+21" symbol="\mathrm{Z}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="exa"    value="1.0d+18" symbol="\mathrm{E}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="peta"   value="1.0d+15" symbol="\mathrm{P}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="tera"   value="1.0d+12" symbol="\mathrm{T}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="giga"   value="1.0d+09" symbol="\mathrm{G}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="mega"   value="1.0d+06" symbol="\mathrm{M}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="kilo"   value="1.0d+03" symbol="\mathrm{k}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="hecto"  value="1.0d+02" symbol="\mathrm{h}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="deca"   value="1.0d+01" symbol="\mathrm{da}" description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="deci"   value="1.0d-01" symbol="\mathrm{d}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="centi"  value="1.0d-02" symbol="\mathrm{c}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="milli"  value="1.0d-03" symbol="\mathrm{m}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="micro"  value="1.0d-06" symbol="\mu"         description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="nano"   value="1.0d-09" symbol="\mathrm{n}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="pico"   value="1.0d-12" symbol="\mathrm{p}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="femto"  value="1.0d-15" symbol="\mathrm{f}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="atto"   value="1.0d-18" symbol="\mathrm{a}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="zepto"  value="1.0d-21" symbol="\mathrm{z}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="yocto"  value="1.0d-24" symbol="\mathrm{y}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="ronto"  value="1.0d-27" symbol="\mathrm{r}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  <constant variable="quecto" value="1.0d-30" symbol="\mathrm{q}"  description="SI prefix." units="dimensionless" unitsInSI="1.0" reference="NIST" referenceURL="https://www.nist.gov/pml/owm/metric-si-prefixes" group="prefixes"/>
  !!]
  
  character(len=2), parameter, dimension(-10:10) :: siPrefix=['q ','r ','y ','z ','a ','f ','p ','n ','μ','m ',' ','k ','M ','G ','T ','P ','E ','Z ','Y ','R ','Q ']

contains

  function siFormat(value_,format_)
    !!{
    Format a value using SI prefixes.
    !!}   
    implicit none
    character       (len=64            )                :: siFormat
    double precision                    , intent(in   ) :: value_
    character       (len=*             ), intent(in   ) :: format_
    integer                                             :: iPrefix
    character       (len=len(format_)+5)                :: formatExtended
    
    if (value_ == 0.0d0) then
       iPrefix=0
    else
       iPrefix=min(max(floor(dble(int(log10(abs(value_)))/3.0d0)),lbound(siPrefix,dim=1)),ubound(siPrefix,dim=1))
    end if
    formatExtended="("//trim(format_)//",a2)"
    write (siFormat,formatExtended) value_/10.0d0**(3*iPrefix),siPrefix(iPrefix)
    return
  end function siFormat
  
end module Numerical_Constants_Prefixes
