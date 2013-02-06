!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a driver program for reading CSV files of merger trees from the Millennium Simulation database and converting to
!% \glc's HDF5 format merger tree files.

program Millennium_Merger_Tree_File_Maker
  !% Driver program for reading CSV files of merger trees from the Millennium Simulation database and converting to \glc's HDF5
  !% format merger tree files.
  use Merger_Tree_Data_Structure
  use Merger_Trees_Millennium
  use Command_Arguments
  use Memory_Management
  implicit none
  integer              :: hdfChunkSize=1024, hdfCompressionLevel=9, generation
  character(len=1024)  :: nodesFile,particlesFile,outputFile,outputFormat
  type(mergerTreeData) :: mergerTrees

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Millennium_Merger_Tree_File_Maker.size')

  ! Get the name of the input and output files.
  if (Command_Argument_Count() < 5 .or. Command_Argument_Count() > 7) stop "Usage: Merger_Tree_File_Maker.exe <nodesFile> <particlesFile> <outputFile> <outputFormat> <generation> [<hdfChunkSize> [hdfCompressionLevel]]"
  call Get_Argument(1,nodesFile     )
  call Get_Argument(2,particlesFile )
  call Get_Argument(3,outputFile    )
  call Get_Argument(4,outputFormat  )
  call Get_Argument(5,generation    )
  if (Command_Argument_Count() >= 6) call Get_Argument(6,hdfChunkSize       )
  if (Command_Argument_Count() == 7) call Get_Argument(7,hdfCompressionLevel)

  ! Process file.
  call Merger_Trees_Millennium_Process(nodesFile,particlesFile,mergerTrees,generation)

  ! Output HDF5 file.
  call mergerTrees%export(outputFile,outputFormat,hdfChunkSize,hdfCompressionLevel)

end program Millennium_Merger_Tree_File_Maker
