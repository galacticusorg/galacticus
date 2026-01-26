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
Contains a module which implements a class for stellar spectra postprocessors.
!!}

module Stellar_Population_Spectra_Postprocess
  !!{
  Implements a class for stellar spectra postprocessors.
  !!}
  private
  public :: stellarPopulationSpectraPostprocessorList

  !![
  <functionClass>
   <name>stellarPopulationSpectraPostprocessor</name>
   <descriptiveName>Postprocessors for stellar population spectra</descriptiveName>
   <description>
    Class providing postprocessors for stellar population spectra. Postprocessors apply effects such as absorption by the \gls{igm}.
   </description>
   <default>inoue2014</default>
   <method name="multiplier" >
    <description>Return the multiplicative modification to the spectrum.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavelength, age, redshift</argument>
   </method>
   <method name="isRedshiftDependent" >
    <description>Return true if the postprocessor is redshift dependent.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      stellarPopulationSpectraPostprocessorIsRedshiftDependent=.true.
    </code>
   </method>
  </functionClass>
  !!]

  type :: stellarPopulationSpectraPostprocessorList
     !!{
     Type used to build linked list of stellar population spectra postprocessors.
     !!}
     class(stellarPopulationSpectraPostprocessorClass), pointer :: stellarPopulationSpectraPostprocessor_ => null()
   contains
     final     ::                  stellarPopulationSpectraPostprocessorListDestructor
     procedure ::                  stellarPopulationSpectraPostprocessorListAssign
     generic   :: assignment(=) => stellarPopulationSpectraPostprocessorListAssign
  end type stellarPopulationSpectraPostprocessorList

  !![
  <functionClass>
   <name>stellarPopulationSpectraPostprocessorBuilder</name>
   <descriptiveName>Builder for postprocessors for stellar population spectra</descriptiveName>
   <description>
    Class providing builders for postprocessors for stellar population spectra. These act as a factory for {\normalfont
    \ttfamily stellarPopulationSpectraPostprocessor} objects. Different postprocessors can be applied to different filters. The
    {\normalfont \ttfamily [luminosityPostprocessSet]} parameter specifies, for each filter, a descriptor which is passed to
    the builder object, which then uses that descriptor to build a postprocessor. (If this parameter is not present then
    ``{\normalfont \ttfamily default}'' is assumed for all filters.)
   </description>
   <default>lookup</default>
   <method name="build" >
    <description>Build and return a postprocessor.</description>
    <type>class(stellarPopulationSpectraPostprocessorClass)</type>
    <pass>yes</pass>
    <argument>type(varying_string), intent(in   ) :: descriptor</argument>
   </method>
  </functionClass>
  !!]

contains

  subroutine stellarPopulationSpectraPostprocessorListDestructor(self)
    !!{
    Destructor for elements of stellar population spectra postprocessor lists.
    !!}
    implicit none
    type(stellarPopulationSpectraPostprocessorList), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarPopulationSpectraPostprocessor_"/>
    !!]
    return
  end subroutine stellarPopulationSpectraPostprocessorListDestructor

  recursive subroutine stellarPopulationSpectraPostprocessorListAssign(self,from)
    !!{
    Perform assignment for the \refClass{stellarPopulationSpectraPostprocessorList} class.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorList), intent(  out) :: self
    class(stellarPopulationSpectraPostprocessorList), intent(in   ) :: from

    self%stellarPopulationSpectraPostprocessor_ => from%stellarPopulationSpectraPostprocessor_
    !![
    <referenceCountIncrement owner="self" object="stellarPopulationSpectraPostprocessor_"/>
    !!]
    return
  end subroutine stellarPopulationSpectraPostprocessorListAssign

end module Stellar_Population_Spectra_Postprocess
