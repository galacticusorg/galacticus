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

!% Compute covariance matrices for mass function estimates.

program Mass_Function_Covariance
  !% Compute covariance matrices for mass function estimates.
  use IO_HDF5
  use ISO_Varying_String
  use Memory_Management
  use Galacticus_Error
  use Input_Parameters
  use Statistics_Mass_Function_Covariance
  use Numerical_Constants_Astronomical
  implicit none
  integer                             , parameter                   :: fileNameLengthMaximum=1024
  double precision                    , allocatable, dimension(:  ) :: mass,massFunction,massFunctionObserved
  double precision                    , allocatable, dimension(:,:) :: covariance,covariancePoisson,covarianceHalo,covarianceLSS&
       &,correlation
  logical                                                           :: massFunctionCovarianceIncludePoisson,massFunctionCovarianceIncludeHalo,massFunctionCovarianceIncludeLSS
  character(len=fileNameLengthMaximum)                              :: parameterFileCharacter
  type     (varying_string           )                              :: parameterFile,massFunctionCovarianceOutputFileName
  type     (hdf5Object               )                              :: outputFile,thisDataset
  integer                                                           :: massFunctionCovarianceBinCount
  double precision                                                  :: massFunctionCovarianceRedshiftMinimum,massFunctionCovarianceRedshiftMaximum,&
       &massFunctionCovarianceMassMinimum ,massFunctionCovarianceMassMaximum

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Mass_Function_Covariance.size')

  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  if (len_trim(parameterFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Mass_Function_Covariance.exe <parameterFile>")
  parameterFile=parameterFileCharacter

  ! Open the parameter file.
  call Input_Parameters_File_Open(parameterFile,outputFile,allowedParametersFile='Mass_Function_Covariance.parameters.xml')

  ! Get the name of the output file.
  !@ <inputParameter>
  !@   <name>massFunctionCovarianceOutputFileName</name>
  !@   <defaultValue>covariance.hdf5</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The name of the file to which the covariance matrix should be written.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceOutputFileName',massFunctionCovarianceOutputFileName,defaultValue='covariance.hdf5',writeOutput=.false.)

  ! Open the output file.
  call outputFile%openFile(char(massFunctionCovarianceOutputFileName),overWrite=.false.)
  ! Read parameters controlling the calculation.
  !@ <inputParameter>
  !@   <name>massFunctionCovarianceRedshiftMinimum</name>
  !@   <defaultValue>0.0</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The minimum redshift at which calculations of the mass function covariance should be carried out.
  !@   </description>
  !@   <type>boolean</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceRedshiftMinimum',massFunctionCovarianceRedshiftMinimum,defaultValue=0.0d0)
  !@ <inputParameter>
  !@   <name>massFunctionCovarianceRedshiftMaximum</name>
  !@   <defaultValue>0.0</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The maximum redshift at which calculations of the mass function covariance should be carried out.
  !@   </description>
  !@   <type>boolean</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceRedshiftMaximum',massFunctionCovarianceRedshiftMaximum,defaultValue=0.1d0)
 !@ <inputParameter>
  !@   <name>massFunctionCovarianceBinCount</name>
  !@   <defaultValue>0.0</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The number of bins in the mass function for covariance calculations.
  !@   </description>
  !@   <type>boolean</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceBinCount',massFunctionCovarianceBinCount,defaultValue=10)
  !@ <inputParameter>
  !@   <name>massFunctionCovarianceMassMinimum</name>
  !@   <defaultValue>$10^8M_\odot$</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The minimum mass in the mass function for covariance calculations.
  !@   </description>
  !@   <type>boolean</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceMassMinimum',massFunctionCovarianceMassMinimum,defaultValue=1.0d08)
  !@ <inputParameter>
  !@   <name>massFunctionCovarianceMassMaximum</name>
  !@   <defaultValue>$10^8M_\odot$</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The maximum mass in the mass function for covariance calculations.
  !@   </description>
  !@   <type>boolean</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceMassMaximum',massFunctionCovarianceMassMaximum,defaultValue=1.0d13)
  !@ <inputParameter>
  !@   <name>massFunctionCovarianceIncludePoisson</name>
  !@   <defaultValue>true</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     Specifies whether or not to include the Poisson contribution to mass function covariance matrices.
  !@   </description>
  !@   <type>boolean</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceIncludePoisson',massFunctionCovarianceIncludePoisson,defaultValue=.true.)
  !@ <inputParameter>
  !@   <name>massFunctionCovarianceIncludeHalo</name>
  !@   <defaultValue>true</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     Specifies whether or not to include the halo contribution to mass function covariance matrices.
  !@   </description>
  !@   <type>boolean</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceIncludeHalo',massFunctionCovarianceIncludeHalo,defaultValue=.true.)
  !@ <inputParameter>
  !@   <name>massFunctionCovarianceIncludeLSS</name>
  !@   <defaultValue>true</defaultValue>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     Specifies whether or not to include the large-scale structure contribution to mass function covariance matrices.
  !@   </description>
  !@   <type>boolean</type>
  !@   <cardinality>1</cardinality>
  !@   <group>output</group>
  !@ </inputParameter>
  call Get_Input_Parameter('massFunctionCovarianceIncludeLSS',massFunctionCovarianceIncludeLSS,defaultValue=.true.)

  ! Read the observed mass function if available.
  if (outputFile%hasDataset("massFunctionObserved")) call outputFile%readDataset("massFunctionObserved",massFunctionObserved)

  ! Compute the covariance matrix.
  call Mass_Function_Covariance_Matrix(massFunctionCovarianceRedshiftMinimum,massFunctionCovarianceRedshiftMaximum&
       &,massFunctionCovarianceBinCount ,massFunctionCovarianceMassMinimum,massFunctionCovarianceMassMaximum,massFunctionObserved&
       & ,massFunctionCovarianceIncludePoisson ,massFunctionCovarianceIncludeHalo,massFunctionCovarianceIncludeLSS,mass &
       &,massFunction,covariance,covariancePoisson ,covarianceHalo,covarianceLSS,correlation)

  ! Write out the covariance matrix.
  call outputFile %writeDataset  (mass               ,"mass"             ,"Mass; M [M☉]"                                         ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(massSolar          ,"unitsInSI"                                    )
  call thisDataset%close         (                                                                   )
  call outputFile %writeDataset  (massFunction       ,"massFunction"     ,"Mass function; dn/dln(M) [Mpc⁻^]"                     ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
  call thisDataset%close         (                                                                   )
  call outputFile %writeDataset  (covariance         ,"covariance"       ,"Covariance of mass function; [Mpc⁻⁶]"            ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
  call thisDataset%close         (                                                                   )
  call outputFile %writeDataset  (covariancePoisson  ,"covariancePoisson","Covariance due to Poisson noise; [Mpc⁻⁶]"        ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
  call thisDataset%close         (                                                                   )
  call outputFile %writeDataset  (covarianceHalo     ,"covarianceHalo"   ,"Covariance due to halo effect; [Mpc⁻⁶]"          ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
  call thisDataset%close         (                                                                   )
  call outputFile %writeDataset  (covarianceLSS      ,"covarianceLSS"    ,"Covariance due to large scale structure; [Mpc⁻⁶]",datasetReturned=thisDataset)
  call thisDataset%writeAttribute(1.0d0/megaParsec**3,"unitsInSI"                                    )
  call thisDataset%close         (                                                                   )
  call outputFile %writeDataset  (correlation        ,"correlation"      ,"Correlation matrix for stellar mass function; []"                            )
  
  ! Close the parameter and output files.
  call Close_Parameters_Group     ()
  call outputFile%close           ()
  call Input_Parameters_File_Close()
  
  ! All done, finish.
end program Mass_Function_Covariance
