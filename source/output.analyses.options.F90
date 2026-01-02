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
Contains a module which provides options and enumerations for on-the-fly analyses.
!!}

module Output_Analyses_Options
  !!{
  Provides options and enumerations for on-the-fly analyses.
  !!}
  public

  !![
  <enumeration>
   <name>outputAnalysisPropertyType</name>
   <description>Property types.</description>
   <entry label="linear"    />
   <entry label="log10"     />
   <entry label="magnitude" />
   <entry label="unknown"   />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>outputAnalysisPropertyQuantity</name>
   <description>Property quantities.</description>
   <entry label="unknown"          />
   <entry label="mass"             />
   <entry label="starFormationRate"/>
   <entry label="luminosity"       />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>outputAnalysisCovarianceModel</name>
   <description>Output analyses covariance models.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="poisson" />
   <entry label="binomial"/>
  </enumeration>
  !!]

end module Output_Analyses_Options
