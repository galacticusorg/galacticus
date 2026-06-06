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
Contains a module which provides a class that implements on-the-fly analyses.
!!}

module Output_Analyses
  !!{
  Provides a class that implements on-the-fly analyses.
  !!}
  use            :: Galacticus_Nodes, only : treeNode, mergerTree
  use, intrinsic :: ISO_C_Binding   , only : c_size_t
  private

  !![
  <functionClass>
   <name>outputAnalysis</name>
   <descriptiveName>Output Analysis</descriptiveName>
   <description>Class providing on-the-fly analysis of \glc\ model outputs, computing statistics such as stellar mass
    functions, luminosity functions, size--mass relations, and other observational comparisons during the simulation
    run rather than as a post-processing step. Analyses are finalized (e.g. accumulated across trees and MPI tasks)
    after all trees have been processed, and results are written to the output \gls{hdf5} file. A log-likelihood
    method enables use within posterior sampling frameworks.</description>
   <default>null</default>
   <method name="analyze" >
    <description>Extract properties from the given \mono{node} at output index \mono{iOutput} and accumulate them into the analysis statistic (e.g. histogram bin counts or likelihood contributions).</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout) :: node</argument>
    <argument>integer(c_size_t), intent(in   ) :: iOutput</argument>
   </method>
   <method name="newTree" >
    <description>Perform any initialization or bookkeeping required at the start of analyzing a new merger tree at output index \mono{iOutput}, for example resetting per-tree accumulators.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type   (mergerTree), intent(inout) :: tree</argument>
    <argument>integer(c_size_t  ), intent(in   ) :: iOutput</argument>
    <code>
     !$GLC attributes unused :: self, tree, iOutput
    </code>
   </method>
   <method name="finalize" >
    <description>Finalize the analysis after all trees and MPI tasks have been processed, computing derived statistics and writing results to the output HDF5 group specified by \mono{groupName}.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(varying_string), intent(in   ), optional :: groupName</argument>
   </method>
   <method name="reduce" >
    <description>Reduce (aggregate) this analysis object's accumulated statistics onto the \mono{reduced} object, typically used for MPI reduction to combine results from multiple processes.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(outputAnalysisClass), intent(inout) :: reduced</argument>
    <code>
     !$GLC attributes unused :: self, reduced
    </code>
   </method>
   <method name="logLikelihood" >
    <description>Return the log-likelihood of the model predictions compared to the target observational data, summed over all bins or data points included in this analysis.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Output_Analyses
