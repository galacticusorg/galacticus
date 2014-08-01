!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program which computes the projected correlation function for galaxies in a given
!% mass range based on a halo occupation model.

program Projected_Correlation_Function
  !% Computes the projected correlation function for galaxies in a given mass range based on a
  !% halo occupation model.
  use, intrinsic :: ISO_C_Binding
  use               FGSL
  use               Command_Arguments
  use               ISO_Varying_String
  use               Memory_Management
  use               Input_Parameters
  use               Galacticus_Error
  use               Geometry_Surveys
  use               Galacticus_Nodes
  use               Cosmology_Functions
  use               Conditional_Mass_Functions
  use               Numerical_Integration
  use               Dark_Matter_Profiles_Concentration
  use               Node_Component_Dark_Matter_Profile_Scale
  use               Numerical_Ranges
  use               Numerical_Constants_Math
  use               Power_Spectra
  use               FFTLogs
  use               Tables
  use               IO_HDF5
  use               Halo_Model_Projected_Correlations
  implicit none
  class           (conditionalMassFunctionClass       ), pointer                   :: conditionalMassFunction_
  double precision                                     , allocatable, dimension(:) :: projectedSeparationBinned                    , projectedCorrelationBinned
  type            (varying_string                     )                            :: parameterFile                                , projectedCorrelationFunctionOutputFileName
  type            (hdf5Object                         )                            :: outputFile
  double precision                                                                 :: projectedCorrelationFunctionSeparationMinimum, projectedCorrelationFunctionSeparationMaximum, &
       &                                                                              projectedCorrelationFunctionMassMinimum      , projectedCorrelationFunctionMassMaximum      , &
       &                                                                              projectedCorrelationFunctionHaloMassMinimum  , projectedCorrelationFunctionHaloMassMaximum  , &
       &                                                                              projectedCorrelationFunctionLineOfSightDepth
  integer                                                                          :: projectedCorrelationFunctionSeparationCount

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Projected_Correlation_Function.size')
  ! Check that correct number of arguments have been supplied.
  if (Command_Argument_Count() /= 1) call Galacticus_Error_Report(message="Usage: Projected_Correlation_Function.exe <parameterFile>")
  ! Get the name of the parameter file from the first command line argument.
  call Get_Argument              (1,parameterFile)
  ! Open the parameter file.
  call Input_Parameters_File_Open(  parameterFile)
  ! Read parameters controlling the calculation.
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionSeparationMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The minimum separation at which to compute the projected correlation function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionSeparationMinimum',projectedCorrelationFunctionSeparationMinimum)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionSeparationMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The maximum separation at which to compute the projected correlation function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionSeparationMaximum',projectedCorrelationFunctionSeparationMaximum)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionSeparationCount</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The number of separations at which to compute the projected correlation function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionSeparationCount',projectedCorrelationFunctionSeparationCount)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionLineOfSightDepth</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The maximum line of sight depth to which to integrate when computing the projected correlation function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionLineOfSightDepth',projectedCorrelationFunctionLineOfSightDepth)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionOutputFileName</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The name of the file to which to output projected correlation function data.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionOutputFileName',projectedCorrelationFunctionOutputFileName)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionMassMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^8M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum mass of galaxies to include in the projected correlation function calculation.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionMassMinimum',projectedCorrelationFunctionMassMinimum,defaultValue=1.0d8)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionMassMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^{12}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum mass of galaxies to include in the projected correlation function calculation.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionMassMaximum',projectedCorrelationFunctionMassMaximum,defaultValue=1.0d12)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionHaloMassMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^6M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum halo mass to use when integrating over the halo mass function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionHaloMassMinimum',projectedCorrelationFunctionHaloMassMinimum,defaultValue=1.0d6)
  !@ <inputParameter>
  !@   <name>projectedCorrelationFunctionHaloMassMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum halo mass to use when integrating over the halo mass function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('projectedCorrelationFunctionHaloMassMaximum',projectedCorrelationFunctionHaloMassMaximum,defaultValue=1.0d16)  
  ! Ensure the nodes objects are initialized.
  call Galacticus_Nodes_Initialize                        ()
  call Node_Component_Dark_Matter_Profile_Scale_Initialize()
  ! Get the default conditional mass function object.
  conditionalMassFunction_  => conditionalMassFunction()
  ! Generate binned separations.
  call Alloc_Array(projectedSeparationBinned ,[projectedCorrelationFunctionSeparationCount])
  call Alloc_Array(projectedCorrelationBinned,[projectedCorrelationFunctionSeparationCount])
  projectedSeparationBinned=Make_Range(projectedCorrelationFunctionSeparationMinimum,projectedCorrelationFunctionSeparationMaximum,projectedCorrelationFunctionSeparationCount,rangeTypeLogarithmic)
  ! Compute projected correlation function.
  call Halo_Model_Projected_Correlation(                                              &
       &                                conditionalMassFunction_                    , &
       &                                projectedSeparationBinned                   , &
       &                                projectedCorrelationFunctionMassMinimum     , &
       &                                projectedCorrelationFunctionMassMaximum     , &
       &                                projectedCorrelationFunctionHaloMassMinimum , &
       &                                projectedCorrelationFunctionHaloMassMaximum , &
       &                                projectedCorrelationFunctionLineOfSightDepth, &
       &                                projectedCorrelationBinned                    &
       &                               )
  ! Write the data to file.
  call outputFile%openFile(char(projectedCorrelationFunctionOutputFileName))
  call outputFile%writeDataset(projectedSeparationBinned ,"separation"          ,commentText="Projected separation in units of Mpc." )
  call outputFile%writeDataset(projectedCorrelationBinned,"projectedCorrelation",commentText="Projected correlation in units of Mpc.")
  call outputFile%close()
  ! Close the parameter file.
  call Input_Parameters_File_Close

end program Projected_Correlation_Function
