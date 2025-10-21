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

!+    Contributions to this file made by:  Alex Merson.

!!{
Contains a module which initializes data for the stellar luminosities class.
!!}

module Stellar_Luminosities_Initialization
  !!{
  Initializes data for the stellar luminosities class.
  !!}
  implicit none
  private
  public :: Stellar_Luminosities_Initialize
  
contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Stellar_Luminosities_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Stellar_Luminosities_Initialize(parameters)
    !!{
    Extract and store a list of output redshifts.
    !!}
    use, intrinsic :: ISO_C_Binding            , only : c_size_t
    use            :: Input_Parameters         , only : inputParameter, inputParameters
    use            :: Output_Times             , only : outputTimes   , outputTimesClass
    use            :: Stellar_Luminosities_Data, only : outputCount   , outputRedshifts
    implicit none
    type   (inputParameters ), intent(inout) :: parameters
    class  (outputTimesClass), pointer       :: outputTimes_
    integer(c_size_t        )                :: i
    
    !![
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    outputCount=outputTimes_%count()
    if (allocated(outputRedshifts)) deallocate(outputRedshifts)
    allocate(outputRedshifts(outputCount))
    do i=1,outputCount
       outputRedshifts(i)=outputTimes_%redshift(i)
    end do
    !![
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end subroutine Stellar_Luminosities_Initialize

end module Stellar_Luminosities_Initialization
