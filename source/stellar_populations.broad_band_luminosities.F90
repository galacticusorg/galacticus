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
Contains a module which provides a class that implements broad band luminosities of stellar populations.
!!}

module Stellar_Population_Broad_Band_Luminosities
  !!{
  Provides a class that implements broad band luminosities of stellar populations.
  !!}
  use :: Stellar_Populations                   , only : stellarPopulationClass
  use :: Stellar_Population_Spectra_Postprocess, only : stellarPopulationSpectraPostprocessorList
  use :: Abundances_Structure                  , only : abundances
  private

  !![
  <functionClass>
   <name>stellarPopulationBroadBandLuminosities</name>
   <descriptiveName>Stellar Population Broad Band Luminosities</descriptiveName>
   <description>Class providing broad band luminosities of stellar populations.</description>
   <default>standard</default>
   <method name="luminosities" >
    <description>Returns the luminosity for a $1 M_\odot$ simple {\normalfont \ttfamily stellarPopulation\_} of given {\normalfont \ttfamily abundances} and {\normalfont \ttfamily age} and observed through the filter specified by {\normalfont \ttfamily filterIndex}.</description>
    <type>double precision, dimension(size(luminosityIndex))</type>
    <pass>yes</pass>
    <argument>integer                                                    , intent(in   ), dimension(:  )              :: luminosityIndex                       , filterIndex </argument>
    <argument>type            (stellarPopulationSpectraPostprocessorList), intent(inout), dimension(:  )              :: stellarPopulationSpectraPostprocessor_              </argument>
    <argument>class           (stellarPopulationClass                   ), intent(inout)                              :: stellarPopulation_                                  </argument>
    <argument>type            (abundances                               ), intent(in   )                              :: abundancesStellar                                   </argument>
    <argument>double precision                                           , intent(in   ), dimension(:  )              :: age                                   , redshift    </argument>
   </method>
   <method name="luminosityTracks" >
    <description>Returns the luminosity for a $1 M_\odot$ simple stellar population of given {\normalfont \ttfamily abundances} drawn from the given {\normalfont \ttfamily stellarPopulation} and observed through the filter specified by {\normalfont \ttfamily filterIndex}, for all available ages.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer                                                    , intent(in   ), dimension(:  )              :: luminosityIndex                       , filterIndex</argument>
    <argument>type            (stellarPopulationSpectraPostprocessorList), intent(inout), dimension(:  )              :: stellarPopulationSpectraPostprocessor_             </argument>
    <argument>class           (stellarPopulationClass                   ), intent(inout)                              :: stellarPopulation_                                 </argument>
    <argument>type            (abundances                               ), intent(in   )                              :: abundancesStellar                                  </argument>
    <argument>double precision                                           , intent(in   ), dimension(:  )              :: redshift                                           </argument>
    <argument>double precision                                           , intent(  out), dimension(:  ), allocatable :: ages                                               </argument>
    <argument>double precision                                           , intent(  out), dimension(:,:), allocatable :: luminosities                                       </argument>
   </method>
  </functionClass>
  !!]

end module Stellar_Population_Broad_Band_Luminosities
