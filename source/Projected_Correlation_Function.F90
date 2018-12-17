!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a program which computes the projected correlation function for galaxies in a given
!% mass range based on a halo occupation model.

program Projected_Correlation_Function
  !% Computes the projected correlation function for galaxies in a given mass range based on a
  !% halo occupation model.
  use, intrinsic :: ISO_C_Binding
  use               Command_Arguments
  use               ISO_Varying_String
  use               Memory_Management
  use               Input_Parameters
  use               Galacticus_Error
  use               Geometry_Surveys
  use               Galacticus_Nodes                        , only : nodeClassHierarchyInitialize
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
  type            (inputParameters                    )                            :: parameters
  double precision                                                                 :: projectedCorrelationFunctionSeparationMinimum, projectedCorrelationFunctionSeparationMaximum, &
       &                                                                              projectedCorrelationFunctionMassMinimum      , projectedCorrelationFunctionMassMaximum      , &
       &                                                                              projectedCorrelationFunctionHaloMassMinimum  , projectedCorrelationFunctionHaloMassMaximum  , &
       &                                                                              projectedCorrelationFunctionLineOfSightDepth
  integer                                                                          :: projectedCorrelationFunctionSeparationCount
  logical                                                                          :: projectedCorrelationFunctionHalfIntegral

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Projected_Correlation_Function.size')
  ! Check that correct number of arguments have been supplied.
  if (Command_Argument_Count() /= 1) call Galacticus_Error_Report(message="Usage: Projected_Correlation_Function.exe <parameterFile>")
  ! Get the name of the parameter file from the first command line argument.
  call Get_Argument              (1,parameterFile)
  ! Open the parameter file.
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Read parameters controlling the calculation.
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionSeparationMinimum</name>
  !#   <cardinality>1</cardinality>
  !#   <description>The minimum separation at which to compute the projected correlation function.</description>
  !#   <source>parameters</source>
  !#   <type>string</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionSeparationMaximum</name>
  !#   <cardinality>1</cardinality>
  !#   <description>The maximum separation at which to compute the projected correlation function.</description>
  !#   <source>parameters</source>
  !#   <type>string</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionSeparationCount</name>
  !#   <cardinality>1</cardinality>
  !#   <description>The number of separations at which to compute the projected correlation function.</description>
  !#   <source>parameters</source>
  !#   <type>string</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionLineOfSightDepth</name>
  !#   <cardinality>1</cardinality>
  !#   <description>The maximum line of sight depth to which to integrate when computing the projected correlation function.</description>
  !#   <source>parameters</source>
  !#   <type>string</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionHalfIntegral</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>.false.</defaultValue>
  !#   <description>Set to {\normalfont \ttfamily true} if the projected correlation function is computed as $w_\mathrm{p}(r_\mathrm{p})=\int_0^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$, instead of the usual $w_\mathrm{p}(r_\mathrm{p})=\int_{-\pi_\mathrm{max}}^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$.</description>
  !#   <source>parameters</source>
  !#   <type>string</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionOutputFileName</name>
  !#   <cardinality>1</cardinality>
  !#   <description>The name of the file to which to output projected correlation function data.</description>
  !#   <source>parameters</source>
  !#   <type>string</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionMassMinimum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>1.0d8</defaultValue>
  !#   <description>The minimum mass of galaxies to include in the projected correlation function calculation.</description>
  !#   <source>parameters</source>
  !#   <type>string</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionMassMaximum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>1.0d12</defaultValue>
  !#   <description>The maximum mass of galaxies to include in the projected correlation function calculation.</description>
  !#   <source>parameters</source>
  !#   <type>string</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionHaloMassMinimum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>1.0d6</defaultValue>
  !#   <description>The minimum halo mass to use when integrating over the halo mass function.</description>
  !#   <source>parameters</source>
  !#   <type>real</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>projectedCorrelationFunctionHaloMassMaximum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>1.0d16</defaultValue>
  !#   <description>The maximum halo mass to use when integrating over the halo mass function.</description>
  !#   <source>parameters</source>
  !#   <type>real</type>
  !# </inputParameter>
  ! Ensure the nodes objects are initialized.
  call nodeClassHierarchyInitialize                       (parameters)
  call Node_Component_Dark_Matter_Profile_Scale_Initialize(parameters)
  ! Get the default conditional mass function object.
  conditionalMassFunction_  => conditionalMassFunction()
  ! Generate binned separations.
  call allocateArray(projectedSeparationBinned ,[projectedCorrelationFunctionSeparationCount])
  call allocateArray(projectedCorrelationBinned,[projectedCorrelationFunctionSeparationCount])
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
       &                                projectedCorrelationFunctionHalfIntegral    , &
       &                                projectedCorrelationBinned                    &
       &                               )
  ! Write the data to file.
  call outputFile%openFile(char(projectedCorrelationFunctionOutputFileName))
  call outputFile%writeDataset(projectedSeparationBinned ,"separation"          ,commentText="Projected separation in units of Mpc." )
  call outputFile%writeDataset(projectedCorrelationBinned,"projectedCorrelation",commentText="Projected correlation in units of Mpc.")
  call outputFile%close()
end program Projected_Correlation_Function
