!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements a class for computation and output of star formation histories for galaxies.

module Star_Formation_Histories
  !% Implements a class for computation and output of star formation histories for galaxies.
  use            :: Abundances_Structure, only : abundances
  use            :: Galacticus_Nodes    , only : treeNode
  use            :: Histories           , only : history
  use, intrinsic :: ISO_C_Binding       , only : c_size_t
  use            :: Kind_Numbers        , only : kind_int8
  implicit none
  private

  !# <functionClass>
  !#  <name>starFormationHistory</name>
  !#  <descriptiveName>Star Formation Histories</descriptiveName>
  !#  <description>
  !#    Class providing recording and output of star formation histories.
  !#
  !#    Conventionally, star formation histories are output as follows:
  !#    \begin{verbatim}
  !#    HDF5 "galacticus.hdf5" {
  !#    GROUP "starFormationHistories" {
  !#       COMMENT "Star formation history data."
  !#       GROUP "Output1" {
  !#          COMMENT "Star formation histories for all trees at each out"
  !#          GROUP "mergerTree1" {
  !#             COMMENT "Star formation histories for each tree."
  !#             DATASET "diskSFH&lt;nodeID&gt;" {
  !#             COMMENT "Star formation history stellar masses of the disk "
  !#                DATATYPE  H5T_IEEE_F64LE
  !#                DATASPACE  SIMPLE { }
  !#             }
  !#             DATASET "diskTime&lt;nodeID&gt;" {
  !#             COMMENT "Star formation history times of the disk component"
  !#                DATATYPE  H5T_IEEE_F64LE
  !#                DATASPACE  SIMPLE { }
  !#             }
  !#             DATASET "spheroidSFH&lt;nodeID&gt;" {
  !#             COMMENT "Star formation history stellar masses of the spher"
  !#                DATATYPE  H5T_IEEE_F64LE
  !#                DATASPACE  SIMPLE { }
  !#             }
  !#             DATASET "spheroidTime&lt;nodeID&gt;" {
  !#             COMMENT "Star formation history times of the spheroid compo"
  !#                DATATYPE  H5T_IEEE_F64LE
  !#                DATASPACE  SIMPLE { }
  !#             }
  !#          }
  !#          GROUP "mergerTree2" {
  !#          .
  !#          .
  !#          .
  !#          }
  !#       }
  !#       GROUP "Output1" {
  !#       .
  !#       .
  !#       .
  !#       }
  !#    }
  !#    }
  !#    \end{verbatim}
  !#    where {\normalfont \ttfamily nodeID} is the index of the relevant node. The specifics of each dataset will depend on the selected star formation history method.
  !#  </description>
  !#  <default>null</default>
  !#  <method name="create" >
  !#   <description>Create the star formation history object.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>type            (history ), intent(inout) :: historyStarFormation</argument>
  !#   <argument>double precision          , intent(in   ) :: timeBegin</argument>
  !#  </method>
  !#  <method name="scales" >
  !#   <description>Set ODE solver absolute scales for a star formation history object.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (history   ), intent(inout) :: historyStarFormation</argument>
  !#   <argument>double precision            , intent(in   ) :: massStellar</argument>
  !#   <argument>type            (abundances), intent(in   ) :: abundancesStellar</argument>
  !#  </method>
  !#  <method name="rate" >
  !#   <description>Record the rate of star formation in this history.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode  ), intent(inout) :: node</argument>
  !#   <argument>type            (history   ), intent(inout) :: historyStarFormation</argument>
  !#   <argument>type            (abundances), intent(in   ) :: abundancesFuel</argument>
  !#   <argument>double precision            , intent(in   ) :: rateStarFormation</argument>
  !#  </method>
  !#  <method name="output" >
  !#   <description>Output the star formation history.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type     (treeNode      ), intent(inout), target :: node</argument>
  !#   <argument>logical                  , intent(in   )         :: nodePassesFilter</argument>
  !#   <argument>type     (history       ), intent(inout)         :: historyStarFormation</argument>
  !#   <argument>integer  (c_size_t      ), intent(in   )         :: indexOutput</argument>
  !#   <argument>integer  (kind=kind_int8), intent(in   )         :: indexTree</argument>
  !#   <argument>character(len=*         ), intent(in   )         :: labelComponent</argument>
  !#  </method>
  !# </functionClass>

end module Star_Formation_Histories
