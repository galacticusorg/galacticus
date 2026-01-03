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
Contains a module which implements a class for stellar populations.
!!}

module Stellar_Populations
  !!{
  Implements a class for stellar populations.
  !!}
  use            :: Abundances_Structure      , only : abundances
  use            :: Hashes                    , only : integerSizeTHash
  use, intrinsic :: ISO_C_Binding             , only : c_size_t
  use            :: Stellar_Population_Spectra, only : stellarPopulationSpectraClass
  implicit none
  private

  !![
  <functionClass>
   <name>stellarPopulation</name>
   <descriptiveName>Stellar Populations</descriptiveName>
   <description>Class providing stellar populations.</description>
   <default>standard</default>
   <data>integer(c_size_t) :: uniqueID_=-1_c_size_t</data>
   <method name="rateRecycling" >
    <description>Return the rate of recycling from this population.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (abundances), intent(in   ) :: abundances_            </argument>
    <argument>double precision            , intent(in   ) :: ageMinimum , ageMaximum</argument>
   </method>
   <method name="uniqueID" >
    <description>Return the uniqueID corresponding to this population.</description>
    <type>integer(c_size_t)</type>
    <pass>yes</pass>
    <code>
     if (self%uniqueID_ &lt; 0_c_size_t) call stellarPopulationUniqueIDAssign(self)
     stellarPopulationUniqueID=self%uniqueID_
    </code>
   </method>
   <method name="rateYield" >
    <description>Return the rate of element yield from this population.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (abundances), intent(in   )           :: abundances_             </argument>
    <argument>double precision            , intent(in   )           :: ageMinimum  , ageMaximum</argument>
    <argument>integer                     , intent(in   ), optional :: elementIndex            </argument>
   </method>
   <method name="rateEnergy" >
    <description>Return the rate of energy input from this population.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (abundances), intent(in   ) :: abundances_            </argument>
    <argument>double precision            , intent(in   ) :: ageMinimum , ageMaximum</argument>
   </method>
   <method name="recycledFractionInstantaneous" >
    <description>Return the recycled fraction from this population in the instantaneous approximation.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="yieldInstantaneous" >
    <description>Return the metal yield from this population in the instantaneous approximation.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="spectra" >
    <description>Return at set of stellar spectra for this population.</description>
    <type>class(stellarPopulationSpectraClass)</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

  ! Dictionary of unique IDs by descriptor.
  type   (integerSizeTHash) :: descriptors
  logical                   :: descriptorsInitialized=.false.
  integer(c_size_t        ) :: uniqueID_             =0_c_size_t

contains

  subroutine stellarPopulationUniqueIDAssign(self)
    !!{
    Assign a unique ID to a stellar population. Populations are distinguished based on the hash of their descriptor.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(stellarPopulationClass), intent(inout) :: self
    type (varying_string        )                :: hashedDescriptor

    hashedDescriptor=self%hashedDescriptor()
    !$omp critical(stellarPopulationUniqueIDAssign)
    if (.not.descriptorsInitialized) then
       call descriptors%initialize()
       descriptorsInitialized=.true.
    end if
    if (.not.descriptors%exists(hashedDescriptor)) then
       uniqueID_=uniqueID_+1_c_size_t
       if (uniqueID_ < 0_c_size_t) call Error_Report('ran out of unique IDs for stellar populations'//{introspection:location})
       call descriptors%set(hashedDescriptor,uniqueID_)
    end if
    self%uniqueID_=descriptors%value(hashedDescriptor)
    !$omp end critical(stellarPopulationUniqueIDAssign)
    return
  end subroutine stellarPopulationUniqueIDAssign

end module Stellar_Populations
