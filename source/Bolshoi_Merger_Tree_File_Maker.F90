!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+ Contributions to this file made by: Stephanie Dörschner.

!% Contains a driver program for reading ASCII files of Rockstar merger trees from the Bolshoi Simulation and converting to
!% \glc's HDF5 format merger tree files.

program Bolshoi_Merger_Tree_File_Maker
  !% Driver program for reading CSV files of Rockstar merger trees from the Bolshoi Simulation and converting to \glc's HDF5
  !% format merger tree files. Run {\normalfont \ttfamily scripts/aux/Bolshoi\_Trees\_Grab.pl} on raw data files before converting to HDF5.
  use Merger_Tree_Data_Structure
  use Merger_Trees_Bolshoi
  use Command_Arguments
  use Memory_Management
  use HDF5
  implicit none
  integer                   :: hdfChunkSize=1024, hdfCompressionLevel=9
  character(len=1024      ) :: nodesFile    , outputFile       , outputFormat
  type     (mergerTreeData) :: mergerTrees

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Bolshoi_Merger_Tree_File_Maker.size')

  ! Get the name of the input and output files.
  if (Command_Argument_Count() < 3 .or. Command_Argument_Count() > 5) stop "Usage: Bolshoi_Merger_Tree_File_Maker.exe <nodesFile> <outputFile> <outputFormat> [<hdfChunkSize> [hdfCompressionLevel]]"
  call Get_Argument(1,nodesFile     )
  call Get_Argument(2,outputFile    )
  call Get_Argument(3,outputFormat  )
  if (Command_Argument_Count() >= 4) call Get_Argument(4,hdfChunkSize       )
  if (Command_Argument_Count() == 5) call Get_Argument(5,hdfCompressionLevel)

  ! Process file.
  call Merger_Trees_Bolshoi_Process(nodesFile,mergerTrees)

  ! Output HDF5 file.
  call mergerTrees%export(outputFile,enumerationMergerTreeFormatEncode(trim(outputFormat)),int(hdfChunkSize,kind=hsize_t),hdfCompressionLevel)

end program Bolshoi_Merger_Tree_File_Maker
