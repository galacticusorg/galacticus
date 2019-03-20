!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Compute covariance matrices for mass function estimates.

program Mass_Function_Covariance
  !% Compute covariance matrices for mass function estimates.
  use IO_HDF5
  use ISO_Varying_String
  use Galacticus_Error
  use Input_Parameters
  use Statistics_Mass_Function_Covariance
  use Numerical_Constants_Astronomical
  implicit none
  integer                             , parameter                   :: fileNameLengthMaximum=1024
  double precision                    , allocatable, dimension(:  ) :: mass,massFunction,massFunctionObserved, completenessObserved,numberObserved,massObserved,massWidthObserved
  double precision                    , allocatable, dimension(:,:) :: covariance,covariancePoisson,covarianceHalo,covarianceLSS&
       &,correlation
  logical                                                           :: massFunctionCovarianceIncludePoisson,massFunctionCovarianceIncludeHalo,massFunctionCovarianceIncludeLSS
  character(len=fileNameLengthMaximum)                              :: parameterFileCharacter
  type     (varying_string           )                              :: parameterFile,massFunctionCovarianceOutputFileName
  type     (hdf5Object               )                              :: outputFile,thisDataset
  type     (inputParameters          )                              :: parameters
  integer                                                           :: massFunctionCovarianceBinCount
  double precision                                                  :: massFunctionCovarianceRedshiftMinimum,massFunctionCovarianceRedshiftMaximum,&
       &massFunctionCovarianceMassMinimum ,massFunctionCovarianceMassMaximum, completenessErrorObserved

  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  if (len_trim(parameterFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Mass_Function_Covariance.exe <parameterFile>")
  parameterFile=parameterFileCharacter

  ! Open the parameter file.
  parameters=inputParameters(parameterFile,outputParametersGroup=outputFile)
  call parameters%markGlobal()

  ! Get the name of the output file.
  !# <inputParameter>
  !#   <name>massFunctionCovarianceOutputFileName</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>var_str('covariance.hdf5')</defaultValue>
  !#   <description>The name of the file to which the covariance matrix should be written.</description>
  !#   <group>output</group>
  !#   <source>globalParameters</source>
  !#   <type>string</type>
  !#   <writeOutput>yes</writeOutput>
  !# </inputParameter>

  ! Open the output file.
  call outputFile%openFile(char(massFunctionCovarianceOutputFileName),overWrite=.false.)
  ! Read parameters controlling the calculation.
  !# <inputParameter>
  !#   <name>massFunctionCovarianceRedshiftMinimum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>0.0d0</defaultValue>
  !#   <description>The minimum redshift at which calculations of the mass function covariance should be carried out.</description>
  !#   <group>output</group>
  !#   <source>globalParameters</source>
  !#   <type>boolean</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>massFunctionCovarianceRedshiftMaximum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>0.1d0</defaultValue>
  !#   <description>The maximum redshift at which calculations of the mass function covariance should be carried out.</description>
  !#   <group>output</group>
  !#   <source>globalParameters</source>
  !#   <type>boolean</type>
  !# </inputParameter>
 !# <inputParameter>
 !#   <name>massFunctionCovarianceBinCount</name>
 !#   <cardinality>1</cardinality>
 !#   <defaultValue>10</defaultValue>
 !#   <description>The number of bins in the mass function for covariance calculations.</description>
 !#   <group>output</group>
 !#   <source>globalParameters</source>
 !#   <type>boolean</type>
 !# </inputParameter>
  !# <inputParameter>
  !#   <name>massFunctionCovarianceMassMinimum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>1.0d08</defaultValue>
  !#   <description>The minimum mass in the mass function for covariance calculations.</description>
  !#   <group>output</group>
  !#   <source>globalParameters</source>
  !#   <type>boolean</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>massFunctionCovarianceMassMaximum</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>1.0d13</defaultValue>
  !#   <description>The maximum mass in the mass function for covariance calculations.</description>
  !#   <group>output</group>
  !#   <source>globalParameters</source>
  !#   <type>boolean</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>massFunctionCovarianceIncludePoisson</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>.true.</defaultValue>
  !#   <description>Specifies whether or not to include the Poisson contribution to mass function covariance matrices.</description>
  !#   <group>output</group>
  !#   <source>globalParameters</source>
  !#   <type>boolean</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>massFunctionCovarianceIncludeHalo</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>.true.</defaultValue>
  !#   <description>Specifies whether or not to include the halo contribution to mass function covariance matrices.</description>
  !#   <group>output</group>
  !#   <source>globalParameters</source>
  !#   <type>boolean</type>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>massFunctionCovarianceIncludeLSS</name>
  !#   <cardinality>1</cardinality>
  !#   <defaultValue>.true.</defaultValue>
  !#   <description>Specifies whether or not to include the large-scale structure contribution to mass function covariance matrices.</description>
  !#   <group>output</group>
  !#   <source>globalParameters</source>
  !#   <type>boolean</type>
  !# </inputParameter>

  ! Read the observed mass function if available.
  completenessErrorObserved=0.0d0
  if (outputFile%hasDataset  ("massFunctionObserved")) call outputFile%readDataset  ("massFunctionObserved",massFunctionObserved     )
  if (outputFile%hasDataset  ("completenessObserved")) call outputFile%readDataset  ("completenessObserved",completenessObserved     )
  if (outputFile%hasDataset  ("numberObserved"      )) call outputFile%readDataset  ("numberObserved"      ,numberObserved           )
  if (outputFile%hasAttribute("completenessError"   )) call outputFile%readAttribute("completenessError"   ,completenessErrorObserved)
  if (outputFile%hasDataset  ("massObserved"        )) call outputFile%readDataset  ("massObserved"        ,massObserved             )
  if (outputFile%hasDataset  ("massWidthObserved"   )) call outputFile%readDataset  ("massWidthObserved"   ,massWidthObserved        )

  ! Compute the covariance matrix.
  call Mass_Function_Covariance_Matrix(massFunctionCovarianceRedshiftMinimum,massFunctionCovarianceRedshiftMaximum&
       &,massFunctionCovarianceBinCount ,massFunctionCovarianceMassMinimum,massFunctionCovarianceMassMaximum,massObserved,massWidthObserved,massFunctionObserved,completenessObserved,numberObserved,completenessErrorObserved&
       & ,massFunctionCovarianceIncludePoisson ,massFunctionCovarianceIncludeHalo,massFunctionCovarianceIncludeLSS,mass &
       &,massFunction,covariance,covariancePoisson ,covarianceHalo,covarianceLSS,correlation)

  ! Write out the covariance matrix.
  call outputFile %writeDataset  (mass               ,"mass"             ,"Mass; M [M☉]"                                    ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(massSolar          ,"unitsInSI"                                    )
  call thisDataset%close         (                                                                   )
  call outputFile %writeDataset  (massFunction       ,"massFunction"     ,"Mass function; dn/dln(M) [Mpc⁻^]"                ,datasetReturned=thisDataset)
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
  
  ! Close the output file.
  call parameters%destroy()
  call outputFile%close  ()

  ! All done, finish.
end program Mass_Function_Covariance
