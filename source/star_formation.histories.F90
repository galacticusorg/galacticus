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

!!{
Contains a module which implements a class for computation and output of star formation histories for galaxies.
!!}

module Star_Formation_Histories
  !!{
  Implements a class for computation and output of star formation histories for galaxies.
  !!}
  use            :: Abundances_Structure      , only : abundances
  use            :: Galacticus_Nodes          , only : treeNode
  use            :: Galactic_Structure_Options, only : enumerationComponentTypeType
  use            :: Histories                 , only : history
  use, intrinsic :: ISO_C_Binding             , only : c_size_t
  use            :: Kind_Numbers              , only : kind_int8
  use            :: Locks                     , only : ompLock
  implicit none
  private

  ! Enumeration of possible star formation history age bin structures.
  !![
  <enumeration>
   <name>starFormationHistoryAges</name>
   <description>Used to specify distribution of age bins in star formation histories.</description>
   <visibility>public</visibility>
   <entry label="arbitrary"      description="Ages are arbitrary and may vary between galaxies/components/outputs."/>
   <entry label="fixedPerOutput" description="Ages are fixed per output."                                          />
   <entry label="fixed"          description="Ages are fixed globally."                                          />
  </enumeration>
  !!]

  !![
  <functionClass>
   <name>starFormationHistory</name>
   <descriptiveName>Star Formation Histories</descriptiveName>
   <description>
     Class providing recording and output of star formation histories.
   </description>
   <default>null</default>
   <method name="create" >
    <description>Create the star formation history object.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout), target   :: node</argument>
    <argument>type            (history ), intent(inout)           :: historyStarFormation</argument>
    <argument>double precision          , intent(in   )           :: timeBegin</argument>
    <argument>double precision          , intent(in   ), optional :: timeEnd</argument>
   </method>
   <method name="scales" >
    <description>Set ODE solver absolute scales for a star formation history object.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (history   ), intent(inout) :: historyStarFormation</argument>
    <argument>type            (treeNode  ), intent(inout) :: node</argument>
    <argument>double precision            , intent(in   ) :: massStellar, massGas</argument>
    <argument>type            (abundances), intent(in   ) :: abundancesStellar</argument>
   </method>
   <method name="rate" >
    <description>Record the rate of star formation in this history.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (treeNode  ), intent(inout) :: node</argument>
    <argument>type            (history   ), intent(inout) :: historyStarFormation</argument>
    <argument>type            (abundances), intent(in   ) :: abundancesFuel</argument>
    <argument>double precision            , intent(in   ) :: rateStarFormation</argument>
   </method>
   <method name="metallicityBoundaries" >
    <description>Return a (zero-indexed) array of metallicity boundaries for this history.</description>
    <type>double precision, allocatable, dimension(:)</type>
    <pass>yes</pass>
   </method>
   <method name="times" >
    <description>Return an array of times for this history \emph{if} the tabulation in time is static per output.</description>
    <type>double precision, allocatable, dimension(:)</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout), optional :: node                </argument>
    <argument>integer         (c_size_t), intent(in   ), optional :: indexOutput         </argument>
    <argument>type            (history ), intent(in   ), optional :: starFormationHistory</argument>
    <argument>logical                   , intent(in   ), optional :: allowTruncation     </argument>
    <argument>double precision          , intent(  out), optional :: timeStart           </argument>
    <modules>Error</modules>
    <code>
      !$GLC attributes unused :: allowTruncation
      if (     present(node).and.     present(indexOutput         )) call Error_Report('only one of `node` and `indexOutput` may be provided'        //{introspection:location})
      if (.not.present(node).and..not.present(indexOutput         )) call Error_Report('one of `node` and `indexOutput` must be provided'            //{introspection:location})
      if (     present(node).and..not.present(starFormationHistory)) call Error_Report('`node` requires that `starFormationHistory` must be provided'//{introspection:location})
      if (present(indexOutput)) then
        allocate(starFormationHistoryTimes(0))
        call Error_Report('times are not static'//{introspection:location})
      else if (present(node)) then
        allocate(starFormationHistoryTimes(size(starFormationHistory%time)))
        starFormationHistoryTimes=starFormationHistory%time
      end if
      if (present(timeStart)) timeStart=0.0d0
    </code>
   </method>
   <method name="timeNext">
     <description></description>
     <type>double precision</type>
     <pass>yes</pass>
     <argument>type(treeNode), intent(inout) :: node                </argument>
     <argument>type(history ), intent(in   ) :: starFormationHistory</argument>
     <modules>Galacticus_Nodes Arrays_Search</modules>
     <code>
       class  (nodeComponentBasic), pointer :: basic
       integer(c_size_t          )          :: i
       basic                        => node%basic()
       i                            =  searchArray(starFormationHistory%time,basic%time())
       starFormationHistoryTimeNext =  starFormationHistory%time(i+1)
     </code>
   </method>
   <method name="masses" >
    <description>Return an array of masses of stars formed for this history.</description>
    <type>double precision, allocatable, dimension(:,:)</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout)           :: node                </argument>
    <argument>type   (history ), intent(in   )           :: starFormationHistory</argument>
    <argument>logical          , intent(in   ), optional :: allowTruncation     </argument>
    <code>
      !$GLC attributes unused :: allowTruncation
      allocate(starFormationHistoryMasses(size(starFormationHistory%data,dim=1),size(starFormationHistory%data,dim=2)))
      starFormationHistoryMasses=starFormationHistory%data
    </code>
   </method>
   <method name="ageDistribution" >
    <description>Return an enumeration member indicating what may be assumed about the distribution of ages in the star formation histories.</description>
    <type>type(enumerationStarFormationHistoryAgesType)</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     starFormationHistoryAgeDistribution=starFormationHistoryAgesArbitrary
    </code>
   </method>
   <method name="update">
     <description>Update the star formation history after an output time is reached.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>type   (treeNode), intent(inout), target :: node                </argument>
     <argument>integer(c_size_t), intent(in   )         :: indexOutput         </argument>
     <argument>type   (history ), intent(inout)         :: historyStarFormation</argument>
     <code>
       !$GLC attributes unused :: self, node, indexOutput, historyStarFormation
     </code>
   </method>
   <method name="rangeIsSufficient">
     <description>Return true if the star formation history spans a sufficient range of times.</description>
     <type>logical</type>
     <pass>yes</pass>
     <argument>type(history), intent(in   ) :: starFormationHistory, rangeHistory</argument>
     <code>
       !$GLC attributes unused :: self, starFormationHistory, rangeHistory
       starFormationHistoryRangeIsSufficient=.true.
     </code>
   </method>
   <method name="extend">
     <description>Extend a star formation history to span a sufficient range of times.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>type            (history), intent(inout)               :: starFormationHistory</argument>
     <argument>double precision         , intent(in   ), dimension(:) :: times               </argument>
     <modules>Error</modules>
     <code>
       !$GLC attributes unused :: self, starFormationHistory, times
       call Error_Report("unexpected attempt to extend star formation history"//{introspection:location})
     </code>
   </method>
   <method name="move">
     <description>Move one star formation history into another.</description>
     <type>void</type>
     <pass>yes</pass>
     <argument>type(treeNode), intent(inout) :: node1                , node2                </argument>
     <argument>type(history ), intent(inout) :: starFormationHistory1, starFormationHistory2</argument>
     <modules>Error</modules>
     <code>
       !$GLC attributes unused :: self, node1, node2
       call starFormationHistory1%increment(starFormationHistory2,autoExtend=.true.)
       call starFormationHistory2%reset    (                                       )
     </code>
   </method>
  </functionClass>
  !!]

end module Star_Formation_Histories
