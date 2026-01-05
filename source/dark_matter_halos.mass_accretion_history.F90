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
Contains a module which provides a class for calculations of dark matter halo mass accretion histories.
!!}

module Dark_Matter_Halo_Mass_Accretion_Histories
  !!{
  Provides a class for calculations of dark matter halo mass accretion histories.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>darkMatterHaloMassAccretionHistory</name>
   <descriptiveName>Dark Matter Halo Mass Accretion Histories</descriptiveName>
   <description>
    Class providing dark matter halo mass accretion histories.
   </description>
   <default>wechsler2002</default>
   <method name="time">
    <description>Returns the time at which the given halo mass was reached.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: mass</argument>
    <modules>Root_Finder Galacticus_Nodes</modules>
    <code>
     class  (nodeComponentBasic), pointer :: basic
     type   (rootFinder        ), save    :: finder
     logical                    , save    :: finderInitialized=.false.
     !$omp threadprivate(finder,finderInitialized)
     if (.not.finderInitialized) then
       finder=  rootFinder(                                                             &amp;
        &amp;              rootFunction                 =timeRoot                     , &amp;
        &amp;              toleranceRelative            =1.0d-3                       , &amp;
        &amp;              rangeExpandDownward          =0.5d+0                       , &amp;
        &amp;              rangeExpandUpward            =2.0d+0                       , &amp;
        &amp;              rangeExpandType              =rangeExpandMultiplicative    , &amp;
        &amp;              rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
        &amp;              rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &amp;
        &amp;             )
       finderInitialized=.true.
     end if
     node_                                  => node
     self_                                  => self
     mass_                                  =  mass
     basic                                  => node  %basic(                      )
     darkMatterHaloMassAccretionHistoryTime =  finder%find (rootGuess=basic%time())
    </code>
   </method>
   <method name="mass">
    <description>Returns the halo mass at the given halo time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target :: node</argument>
    <argument>double precision          , intent(in   )         :: time</argument>
    <modules>Root_Finder Galacticus_Nodes</modules>
    <code>
     class  (nodeComponentBasic), pointer :: basic
     type   (rootFinder        ), save    :: finder
     logical                    , save    :: finderInitialized=.false.
     !$omp threadprivate(finder,finderInitialized)
     if (.not.finderInitialized) then
       finder=  rootFinder(                                                             &amp;
        &amp;              rootFunction                 =massRoot                     , &amp;
        &amp;              toleranceRelative            =1.0d-3                       , &amp;
        &amp;              rangeExpandDownward          =0.5d+0                       , &amp;
        &amp;              rangeExpandUpward            =2.0d+0                       , &amp;
        &amp;              rangeExpandType              =rangeExpandMultiplicative    , &amp;
        &amp;              rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
        &amp;              rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &amp;
        &amp;             )
       finderInitialized=.true.
     end if
     node_                                  => node
     self_                                  => self
     time_                                  =  time
     basic                                  => node  %basic(                      )
     darkMatterHaloMassAccretionHistoryMass =  finder%find (rootGuess=basic%mass())
    </code>
   </method>
   <method name="massAccretionRate">
    <description>Returns the mass accretion rate at the specified time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: time</argument>
   </method>
  </functionClass>
  !!]

  ! Module-scope variables used in root-finding.
  class           (darkMatterHaloMassAccretionHistoryClass), pointer :: self_
  type            (treeNode                               ), pointer :: node_
  double precision                                                   :: time_, mass_
  !$omp threadprivate(self_,node_,time_,mass_)
  
contains

  double precision function timeRoot(time)
    !!{
    Root function used for solving for the time at a given mass.
    !!}
    implicit none
    double precision, intent(in   ) :: time

    timeRoot=self_%mass(node_,time)-mass_
    return
  end function timeRoot
  
  double precision function massRoot(mass)
    !!{
    Root function used for solving for the mass at a given time.
    !!}
    implicit none
    double precision, intent(in   ) :: mass

    massRoot=self_%time(node_,mass)-time_
    return
  end function massRoot
  
end module Dark_Matter_Halo_Mass_Accretion_Histories
