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

!% Contains a code which computes dark matter halo spin distributions.

program haloSpinDistributions
  !% Computes dark matter halo spin distributions.
  use Galacticus_Error
  use Galacticus_Nodes          , only : treeNode                    , nodeComponentBasic, nodeComponentSpin, nodeComponentDarkMatterProfile, &
       &                                 nodeClassHierarchyInitialize
  use Input_Parameters
  use Functions_Global_Utilities
  use IO_HDF5
  use ISO_Varying_String
  use Cosmology_Functions
  use Halo_Spin_Distributions
  use Node_Components
  use Memory_Management
  implicit none
  integer                                         , parameter                   :: fileNameLengthMaximum=1024
  class           (cosmologyFunctionsClass       ), pointer                     :: cosmologyFunctions_
  class           (haloSpinDistributionClass     ), pointer                     :: haloSpinDistribution_
  type            (treeNode                      ), pointer                     :: node
  class           (nodeComponentBasic            ), pointer                     :: nodeBasic
  class           (nodeComponentSpin             ), pointer                     :: nodeSpin
  class           (nodeComponentDarkMatterProfile), pointer                     :: nodeDarkMatterProfile
  double precision                                , allocatable, dimension(:  ) :: outputRedshifts           , spin
  double precision                                , allocatable, dimension(:,:) :: spinDistribution
  character       (len=fileNameLengthMaximum     )                              :: fileCharacter
  type            (varying_string                )                              :: outputFileName            , parameterFileName
  type            (hdf5Object                    )                              :: outputFile
  type            (inputParameters               )                              :: parameters
  integer                                                                       :: iOutput                   , spinCount          , &
       &                                                                           iSpin
  double precision                                                              :: spinMinimum               , spinMaximum        , &
       &                                                                           spinPointsPerDecade       , haloMassMinimum
  
  ! Get the name of the parameter file from the first command line argument.
  if (Command_Argument_Count() /= 2) call Galacticus_Error_Report(message="Usage: haloSpinDistribution.exe <parameterFile> <outputFile>")
  call Get_Command_Argument(1,fileCharacter)
  parameterFileName=fileCharacter
  call Get_Command_Argument(2,fileCharacter)
  outputFileName   =fileCharacter
  ! Initialize global functions.
  call Functions_Global_Set()
  ! Open the output file.
  call outputFile%openFile(char(outputFileName),overWrite=.true.,objectsOverwritable=.false.)
  ! Open the parameter file.
  parameters=inputParameters(parameterFileName,outputParametersGroup=outputFile)
  call parameters%markGlobal()
  ! Initialize nodes and components.
  call nodeClassHierarchyInitialize(parameters)
  call Node_Components_Initialize  (parameters)
  ! Get required objects.
  cosmologyFunctions_   => cosmologyFunctions  ()
  haloSpinDistribution_ => haloSpinDistribution()
  ! Read parameters controlling spin distribution tabulation.
  !# <inputParameter>
  !#   <name>spinMinimum</name>
  !#   <variable>spinMinimum</variable>
  !#   <defaultValue>3.0d-4</defaultValue>
  !#   <description>Minimum spin for which the distribution function should be calculated.</description>
  !#   <source>parameters</source>
  !#   <type>real</type>
  !#   <cardinality>0..1</cardinality>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>spinMaximum</name>
  !#   <variable>spinMaximum</variable>
  !#   <defaultValue>0.5d0</defaultValue>
  !#   <description>Maximum spin for which the distribution function should be calculated.</description>
  !#   <source>parameters</source>
  !#   <type>real</type>
  !#   <cardinality>0..1</cardinality>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>spinPointsPerDecade</name>
  !#   <variable>spinPointsPerDecade</variable>
  !#   <defaultValue>10.0d0</defaultValue>
  !#   <description>Number of points per decade of spin at which to calculate the distribution.</description>
  !#   <source>parameters</source>
  !#   <type>real</type>
  !#   <cardinality>0..1</cardinality>
  !# </inputParameter>
  !# <inputParameter>
  !#   <name>haloMassMinimum</name>
  !#   <variable>haloMassMinimum</variable>
  !#   <defaultValue>0.0d0</defaultValue>
  !#   <description>Minimum halo mass above which spin distribution should be averaged.</description>
  !#   <source>parameters</source>
  !#   <type>real</type>
  !#   <cardinality>0..1</cardinality>
  !# </inputParameter>
  call allocateArray(outputRedshifts,[max(globalParameters%count('outputRedshifts',zeroIfNotPresent=.true.),1)])
  !# <inputParameter>
  !#   <name>outputRedshifts</name>
  !#   <variable>outputRedshifts</variable>
  !#   <defaultValue>[0.0d0]</defaultValue>
  !#   <description>Redshifts for which the spin distribution should be computed.</description>
  !#   <source>parameters</source>
  !#   <type>real</type>
  !#   <cardinality>0..1</cardinality>
  !# </inputParameter>
  ! Create a tree node.
  node                  => treeNode                  (                 )
  nodeBasic             => node    %basic            (autoCreate=.true.)
  nodeSpin              => node    %spin             (autoCreate=.true.)
  nodeDarkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
  ! Compute number of spin points to tabulate.
  spinCount=int(log10(spinMaximum/spinMinimum)*spinPointsPerDecade)+1
  call allocateArray(spin            ,[                      spinCount])
  call allocateArray(spinDistribution,[size(outputRedshifts),spinCount])
  ! Iterate over output redshifts.
  do iOutput=1,size(outputRedshifts)
     ! Set the epoch for the node.
     call   nodeBasic            %timeSet                    (                           &
          &  cosmologyFunctions_ %cosmicTime                  (                          &
          &   cosmologyFunctions_%expansionFactorFromRedshift  (                         &
          &                                                     outputRedshifts(iOutput) &
          &                                                    )                         &
          &                                                   )                          &
          &                                                  )
     ! Iterate over spins.
     do iSpin=1,spinCount
        spin(iSpin)=exp(log(spinMinimum)+log(spinMaximum/spinMinimum)*dble(iSpin-1)/dble(spinCount-1))
        call nodeSpin%spinSet(spin(iSpin))
        ! Evaluate the distribution.
        if (haloMassMinimum <= 0.0d0) then
           ! No minimum halo mass specified - simply evaluate the spin distribution.
           spinDistribution(iOutput,iSpin)=haloSpinDistribution_%distribution(node)
        else
           ! A minimum halo mass was specified. Evaluate the spin distribution averaged over halo mass, if the distribution class
           ! supports this.
           select type (haloSpinDistribution_)
           class is (haloSpinDistributionNbodyErrors)
              spinDistribution(iOutput,iSpin)=haloSpinDistribution_%distributionAveraged(node,haloMassMinimum)
           class default
              call Galacticus_Error_Report('halo spin distribution class does not support averaging over halo mass'//{introspection:location})
           end select
        end if
     end do
  end do
  ! Store the distribution, redshifts, and spins.
  call outputFile%writeDataset(outputRedshifts ,'redshift'    ,'Redshifts at which the spin distribution is tabulated.')
  call outputFile%writeDataset(spin            ,'spin'        ,'Spins at which the spin distribution is tabulated.'    )
  call outputFile%writeDataset(spinDistribution,'distribution','Spin parameter distribution.'                          )
  ! Close the output file.
  call parameters%destroy()
  call outputFile%close  ()
  ! Uninitialize node components
  call Node_Components_Uninitialize()
end program haloSpinDistributions
