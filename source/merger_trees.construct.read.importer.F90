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
Contains a module which provides an object that implements importing of merger trees from file.
!!}

module Merger_Tree_Read_Importers
  !!{
  Provides an object that implements importing of merger trees from file.
  !!}
  use            :: Galacticus_Nodes  , only : treeNode
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: ISO_Varying_String, only : varying_string
  use            :: Kind_Numbers      , only : kind_int8
  private
  public :: nodeData, nodeDataMinimal
  !![
  <workaround type="gfortran" PR="88632" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=88632">
   <description>
    importerUnitConvert is used by submodules, so must be exported to the object file. gfortran currently does not do this if
    the symbol is private, so we mark it as public.
   </description>
  </workaround>
  !!]
  public :: importerUnitConvert

  ! Type used to specify units.
  type :: importerUnits
     logical          :: status
     double precision :: unitsInSI
     integer          :: hubbleExponent, scaleFactorExponent
   contains
     !![
     <methods>
       <method description="Multiply by another {\normalfont \ttfamily importerUnits} object." method="operator(*)" />
       <method description="Raise to the given integer power."                                 method="operator(**)"/>
       <method description="Return true if the provided units are equal."                      method="operator(==)"/>
       <method description="Return true if the provided units are not equal."                  method="operator(/=)"/>
     </methods>
     !!]
     procedure :: multiply     => importerUnitsMultiply
     procedure :: exponentiate => importerUnitsExponentiate
     procedure :: areEqual     => importerUnitsAreEqual
     procedure :: areNotEqual  => importerUnitsAreNotEqual
     generic   :: operator(*)  => multiply
     generic   :: operator(**) => exponentiate
     generic   :: operator(==) => areEqual
     generic   :: operator(/=) => areNotEqual
  end type importerUnits

  ! Types used to store raw data.
  type :: nodeDataMinimal
     !!{
     Structure used to store minimal raw data read from merger tree files.
     !!}
     integer         (kind=kind_int8)               :: descendantIndex, hostIndex, &
          &                                            nodeIndex
     double precision                               :: nodeTime       , nodeMass
  end type nodeDataMinimal

  ! Type used to store raw data.
  type, extends(nodeDataMinimal) :: nodeData
     !!{
     Structure used to store default raw data read from merger tree files.
     !!}
     integer         (kind_int8)                            :: isolatedNodeIndex            , mergesWithIndex                    , &
          &                                                    particleCount                , primaryIsolatedNodeIndex
     double precision                                       :: angularMomentum              , halfMassRadius                     , &
          &                                                    velocityDispersion           , spin                               , &
          &                                                    scaleRadius                  , velocityMaximum
     double precision           , dimension(3)              :: position                     , velocity                           , &
          &                                                    angularMomentum3D            , spin3D
     double precision           , dimension(:), allocatable :: reals
     integer         (kind_int8), dimension(:), allocatable :: integers
     logical                                                :: childIsSubhalo               , isSubhalo
     class           (nodeData ), pointer                   :: descendant         => null() , host                      => null(), &
          &                                                    parent             => null()
     type            (treeNode ), pointer                   :: node               => null()
  end type nodeData

  interface importerUnitConvert
     !!{
     Unit converters for merger tree importers.
     !!}
     module procedure importerUnitConvertScalar
     module procedure importerUnitConvert1D
     module procedure importerUnitConvert2D
  end interface importerUnitConvert

  !![
  <functionClass>
   <name>mergerTreeImporter</name>
   <description>
    Class providing functions for importing merger trees. When merger trees are to be read from file, a number of different
    file formats are supported. This ``importer'' class is used to read these files and place the contents into internal data
    structures that \glc\ can then manipulate.
   </description>
   <descriptiveName>Merger Tree Importer</descriptiveName>
   <default>galacticus</default>
   <method name="open" >
    <description>Opens the file.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(varying_string), intent(in   ) :: fileName</argument>
   </method>
   <method name="treesHaveSubhalos" >
    <description>Returns a Boolean integer specifying whether or not the trees have subhalos.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="massesIncludeSubhalos" >
    <description>Returns a Boolean specifying whether halo masses include the contribution from their subhalos.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="angularMomentaIncludeSubhalos" >
    <description>Returns a Boolean specifying whether halo angular momenta (or spins) include the contribution from their subhalos.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="treesAreSelfContained" >
    <description>Returns a Boolean integer specifying whether trees are self-contained.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="velocitiesIncludeHubbleFlow" >
    <description>Returns a Boolean integer specifying whether velocities include the Hubble flow.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="positionsArePeriodic" >
    <description>Returns a Boolean integer specifying whether positions are periodic.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="canReadSubsets" >
    <description>Returns true if arbitrary subsets of halos from a forest can be read.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="cubeLength" >
    <description>Returns the length of the simulation cube.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )           :: time</argument>
    <argument>integer         , intent(  out), optional :: status</argument>
   </method>
   <method name="treeCount" >
    <description>Returns a count of the number of trees available.</description>
    <type>integer(kind=c_size_t)</type>
    <pass>yes</pass>
   </method>
   <method name="treeIndex" >
    <description>Returns the index of the $i^\mathrm{th}$ tree.</description>
    <type>integer(kind=kind_int8)</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: i</argument>
   </method>
   <method name="nodeCount" >
    <description>Returns the number of nodes in the $i^\mathrm{th}$ tree.</description>
    <type>integer(kind=c_size_t)</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: i</argument>
   </method>
   <method name="treeWeight" >
    <description>Returns the weight to assign to the $i^\mathrm{th}$ tree.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: i</argument>
   </method>
   <method name="positionsAvailable" >
    <description>Return true if positions and/or velocities are available.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>logical, intent(in   ) :: positions, velocities</argument>
   </method>
   <method name="scaleRadiiAvailable" >
    <description>Return true if scale radii are available.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="particleCountAvailable" >
    <description>Return true if particle counts are available.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="velocityMaximumAvailable" >
    <description>Return true if rotation curve velocity maxima are available.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="velocityDispersionAvailable" >
    <description>Return true if halo velocity dispersions are available.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="angularMomentaAvailable" >
    <description>Return true if angular momenta (magnitudes) are available.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="angularMomenta3DAvailable" >
    <description>Return true if angular momenta (vectors) are available.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="spinAvailable" >
    <description>Return true if spin (magnitudes) are available.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="spin3DAvailable" >
    <description>Return true if spin (vectors) are available.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="import" >
    <description>Imports the $i^\mathrm{th}$ tree.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer                 , intent(in   )                            :: i</argument>
    <argument>class  (nodeDataMinimal), intent(  out), allocatable, dimension(:) :: nodes</argument>
    <argument>integer(c_size_t       ), intent(in   ), optional   , dimension(:) :: nodeSubset</argument>
    <argument>logical                 , intent(in   ), optional                  :: requireScaleRadii, requireAngularMomenta, requireAngularMomenta3D, requireSpin, requireSpin3D, requirePositions, structureOnly</argument>
    <argument>type   (varying_string ), intent(in   ), optional   , dimension(:) :: requireNamedReals, requireNamedIntegers</argument>
   </method>
   <method name="subhaloTrace" >
    <description>Supplies epochs, positions, and velocities for traced subhalos.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class           (nodeData), intent(in   )                 :: node</argument>
    <argument>double precision          , intent(  out), dimension(:  ) :: time</argument>
    <argument>double precision          , intent(  out), dimension(:,:) :: position, velocity</argument>
   </method>
   <method name="subhaloTraceCount" >
    <description>Returns the length of a node's subhalo trace.</description>
    <type>integer(kind=c_size_t)</type>
    <pass>yes</pass>
    <argument>class(nodeData), intent(in   ) :: node</argument>
   </method>
  </functionClass>
  !!]

contains

  function importerUnitsMultiply(units1,units2)
    !!{
    Multiply to {\normalfont \ttfamily importerUnits} objects.
    !!}
    implicit none
    type (importerUnits)                :: importerUnitsMultiply
    class(importerUnits), intent(in   ) :: units1
    type (importerUnits), intent(in   ) :: units2

    importerUnitsMultiply%status             =units1%status             .and.units2%status
    importerUnitsMultiply%unitsInSI          =units1%unitsInSI            *  units2%unitsInSI
    importerUnitsMultiply%scaleFactorExponent=units1%scaleFactorExponent  +  units2%scaleFactorExponent
    importerUnitsMultiply%hubbleExponent     =units1%hubbleExponent       +  units2%hubbleExponent
    return
  end function importerUnitsMultiply

  function importerUnitsExponentiate(units1,exponent)
    !!{
    Exponentiate {\normalfont \ttfamily importerUnits} objects.
    !!}
    implicit none
    type   (importerUnits)                :: importerUnitsExponentiate
    class  (importerUnits), intent(in   ) :: units1
    integer               , intent(in   ) :: exponent

    importerUnitsExponentiate%status             =units1%status
    importerUnitsExponentiate%unitsInSI          =units1%unitsInSI          **exponent
    importerUnitsExponentiate%scaleFactorExponent=units1%scaleFactorExponent* exponent
    importerUnitsExponentiate%hubbleExponent     =units1%hubbleExponent     * exponent
    return
  end function importerUnitsExponentiate

  logical function importerUnitsAreEqual(units1,units2)
    !!{
    Test whether two {\normalfont \ttfamily importerUnits} objects are equal.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(importerUnits), intent(in   ) :: units1
    type (importerUnits), intent(in   ) :: units2

    if (.not.(units1%status.and.units2%status)) call Error_Report('units not defined'//{introspection:location})
    importerUnitsAreEqual= units1%unitsInSI           ==  units2%unitsInSI           &
         &                .and.                                                      &
         &                 units1%scaleFactorExponent ==  units2%scaleFactorExponent &
         &                .and.                                                      &
         &                 units1%hubbleExponent      ==  units2%hubbleExponent
    return
  end function importerUnitsAreEqual

  logical function importerUnitsAreNotEqual(units1,units2)
    !!{
    Test whether two {\normalfont \ttfamily importerUnits} objects are not equal.
    !!}
    implicit none
    class(importerUnits), intent(in   ) :: units1
    type (importerUnits), intent(in   ) :: units2

    importerUnitsAreNotEqual=.not.(units1 == units2)
    return
  end function importerUnitsAreNotEqual

  function importerUnitConvertScalar(values,times,units,requiredUnits,cosmologyParameters_,cosmologyFunctions_)
    !!{
    Convert a set of values for \glc\ internal units.
    !!}
    use :: Cosmology_Functions , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters, only : cosmologyParametersClass, hubbleUnitsLittleH
    use :: Error               , only : Error_Report
    implicit none
    double precision                          , intent(in   ) :: values                    , times
    type            (importerUnits           ), intent(in   ) :: units
    double precision                          , intent(in   ) :: requiredUnits
    class           (cosmologyParametersClass), intent(inout) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(inout) :: cosmologyFunctions_
    double precision                                          :: importerUnitConvertScalar

    if (.not.units%status) call Error_Report('units are not defined'//{introspection:location})
    importerUnitConvertScalar=values*(units%unitsInSI/requiredUnits)*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**units%hubbleExponent
    if (units%scaleFactorExponent /= 0) then
       importerUnitConvertScalar=importerUnitConvertScalar*cosmologyFunctions_%expansionFactor(times)**units%scaleFactorExponent
    end if
    return
  end function importerUnitConvertScalar

  function importerUnitConvert1D(values,times,units,requiredUnits,cosmologyParameters_,cosmologyFunctions_)
    !!{
    Convert a set of values for \glc\ internal units.
    !!}
    use :: Cosmology_Functions , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters, only : cosmologyParametersClass, hubbleUnitsLittleH
    use :: Error               , only : Error_Report
    implicit none
    double precision                          , intent(in   ), dimension(           :) :: values               , times
    type            (importerUnits           ), intent(in   )                          :: units
    double precision                          , intent(in   )                          :: requiredUnits
    class           (cosmologyParametersClass), intent(inout)                          :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(inout)                          :: cosmologyFunctions_
    double precision                                         , dimension(size(values)) :: importerUnitConvert1D
    integer                                                                            :: i

    if (.not.units%status) call Error_Report('units are not defined'//{introspection:location})
    importerUnitConvert1D=values*(units%unitsInSI/requiredUnits)*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**units%hubbleExponent
    if (units%scaleFactorExponent /= 0) then
       do i=1,size(values)
          importerUnitConvert1D(i)=importerUnitConvert1D(i)*cosmologyFunctions_%expansionFactor(times(i))**units%scaleFactorExponent
       end do
    end if
    return
  end function importerUnitConvert1D

  function importerUnitConvert2D(values,times,units,requiredUnits,cosmologyParameters_,cosmologyFunctions_)
    !!{
    Convert a set of values for \glc\ internal units.
    !!}
    use :: Cosmology_Functions , only : cosmologyFunctionsClass
    use :: Cosmology_Parameters, only : cosmologyParametersClass, hubbleUnitsLittleH
    use :: Error               , only : Error_Report
    implicit none
    double precision                          , intent(in   ), dimension(                 :,                 :) :: values
    double precision                          , intent(in   ), dimension(                                    :) :: times
    type            (importerUnits           ), intent(in   )                                                   :: units
    double precision                          , intent(in   )                                                   :: requiredUnits
    class           (cosmologyParametersClass), intent(inout)                                                   :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(inout)                                                   :: cosmologyFunctions_
    double precision                                         , dimension(size(values,dim=1),size(values,dim=2)) :: importerUnitConvert2D
    integer                                                                                                     :: i

    if (.not.units%status) call Error_Report('units are not defined'//{introspection:location})
    importerUnitConvert2D=values*(units%unitsInSI/requiredUnits)*cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**units%hubbleExponent
    if (units%scaleFactorExponent /= 0) then
       do i=1,size(values,dim=2)
          importerUnitConvert2D(:,i)=importerUnitConvert2D(:,i)*cosmologyFunctions_%expansionFactor(times(i))**units%scaleFactorExponent
       end do
    end if
    return
  end function importerUnitConvert2D

end module Merger_Tree_Read_Importers
