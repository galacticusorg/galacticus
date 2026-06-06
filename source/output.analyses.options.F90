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
   <description>Enumeration of the scaling types used for output analysis properties, distinguishing linear, logarithmic (base-10), magnitude, and unknown scalings.</description>
   <entry label="linear"    />
   <entry label="log10"     />
   <entry label="magnitude" />
   <entry label="unknown"   />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>outputAnalysisPropertyQuantity</name>
   <description>Enumeration of the physical quantity types represented by output analysis properties, such as mass, star formation rate, luminosity, or unknown.</description>
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

  !![
  <enumeration>
   <name>outputAnalysisState</name>
   <description>Output analyses states.</description>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <entry label="unknown" />
   <entry label="zero"    />
   <entry label="positive"/>
   <entry label="negative"/>
  </enumeration>
  !!]

contains

  function outputAnalysisState(distribution) result(state)
    !!{
    Determine the state of an output analysis distribution.
    !!}
    implicit none
    type            (enumerationOutputAnalysisStateType)               :: state
    double precision                                    , dimension(:) :: distribution

    state=outputAnalysisStateUnknown
    if      (all(distribution == 0.0d0)) then
       state=outputAnalysisStateZero
    else if (all(distribution >= 0.0d0)) then
       state=outputAnalysisStatePositive
    else if (all(distribution <= 0.0d0)) then
       state=outputAnalysisStateNegative
    end if
    return
  end function outputAnalysisState
  
end module Output_Analyses_Options
