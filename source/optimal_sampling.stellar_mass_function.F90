!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program which computes the optimal sampling of the halo mass function to minimize errors on the stellar mass
!% function.

program Optimal_Sampling_SMF
  !% Compute the optimal number of trees to run of each mass.
  use, intrinsic :: ISO_C_Binding
  use Memory_Management
  use Numerical_Ranges
  use Numerical_Constants_Astronomical
  use Input_Parameters
  use Galacticus_Error
  use ISO_Varying_String
  use FGSL
  use Numerical_Integration
  use Halo_Mass_Function
  use Galacticus_Meta_Compute_Times
  use Cosmology_Functions
  use Merger_Trees_Mass_Function_Sampling
  use IO_HDF5
  implicit none
  integer                                     , parameter                 :: fileNameLengthMaximum=1024
  character       (len=fileNameLengthMaximum )                            :: parameterFileCharacter
  type            (varying_string            )                            :: parameterFile,optimalSamplingDensityOutputFileName
  integer                                     , parameter                 :: stellarMassPointsPerDecade=100
  integer                                     , parameter                 :: haloMassPointsPerDecade   =100
  double precision                            , parameter                 :: stellarMassMinimum=1.0d8 , stellarMassMaximum=1.0d13
  double precision                            , parameter                 :: haloMassMinimum   =1.0d10, haloMassMaximum   =1.0d15
  double precision                            , parameter                 :: toleranceAbsolute =1.0d-16,  toleranceRelative =1.0d-3
  double precision                            , allocatable, dimension(:) :: stellarMassTableMass,stellarMassTableMassFunction,haloMassTableMass,haloMassTableXi,treeTiming,haloMassFunction,samplingDensity
  integer                                                                 :: stellarMassTableCount,haloMassTableCount,iMass,iUnit
  double precision                                                        :: stellarMass,time,optimalSamplingLogarithmicBinWidth,haloMass
  class           (cosmologyFunctionsClass   ), pointer                   :: cosmologyFunctionsDefault
  type            (c_ptr                     )                            :: parameterPointer
  type            (fgsl_function             )                            :: integrandFunction
  type            (fgsl_integration_workspace)                            :: integrationWorkspace
  logical                                                                 :: integrationReset=.true.
  type            (hdf5Object                )                            :: outputFile,thisDataset

  ! Read in basic code memory usage.
  call Code_Memory_Usage('optimal_sampling.stellar_mass_function.size')

  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,parameterFileCharacter)
  if (len_trim(parameterFileCharacter) == 0) call Galacticus_Error_Report(message="Usage: optimal_sampling.stellar_mass_function.exe <parameterFile>")
  parameterFile=parameterFileCharacter

  ! Open the parameter file.
  call Input_Parameters_File_Open(parameterFile)

  ! Get the name of the output file.
  !@ <inputParameter>
  !@   <name>optimalSamplingDensityOutputFileName</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The name of a file to which the optimal tree mass sampling density function will be written.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('optimalSamplingDensityOutputFileName',optimalSamplingDensityOutputFileName)
  !@ <inputParameter>
  !@   <name>optimalSamplingLogarithmicBinWidth</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The logarithmic width of bins in the stellar mass function to use when computing the optimal tree mass sampling density function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('optimalSamplingLogarithmicBinWidth',optimalSamplingLogarithmicBinWidth)

  ! Compute number of points to tabulate.
  stellarMassTableCount=int(log10(stellarMassMaximum/stellarMassMinimum)*dble(stellarMassPointsPerDecade)+1.0d0)
  haloMassTableCount   =int(log10(   haloMassMaximum/   haloMassMinimum)*dble(   haloMassPointsPerDecade)+1.0d0)

  ! Allocate arrays.
  call Alloc_Array(stellarMassTableMass        ,[stellarMassTableCount])
  call Alloc_Array(stellarMassTableMassFunction,[stellarMassTableCount])
  call Alloc_Array(   haloMassTableMass        ,[   haloMassTableCount])
  call Alloc_Array(   haloMassTableXi          ,[   haloMassTableCount])
  call Alloc_Array(   haloMassFunction         ,[   haloMassTableCount])
  call Alloc_Array(samplingDensity             ,[   haloMassTableCount])
  call Alloc_Array(treeTiming                  ,[   haloMassTableCount])

  ! Create mass tabulations.
  stellarMassTableMass=Make_Range(stellarMassMinimum,stellarMassMaximum,stellarMassTableCount,rangeType=rangeTypeLogarithmic)
  haloMassTableMass   =Make_Range(   haloMassMinimum,   haloMassMaximum,   haloMassTableCount,rangeType=rangeTypeLogarithmic)
  ! Get the default cosmology functions object.
  cosmologyFunctionsDefault => cosmologyFunctions()
  ! Get the cosmic time at the present day.
  time=cosmologyFunctionsDefault%cosmicTime(0.93457d0)

  ! Open an HDF5 file to write the results to.
  call outputFile%openFile(char(optimalSamplingDensityOutputFileName),overWrite=.true.)

  ! Loop over stellar masses.
  do iMass=1,stellarMassTableCount

     ! Set the stellar mass to use.
     stellarMass=stellarMassTableMass(iMass)

     ! Integrate over the mass function weighting by the stellar conditional mass function.
     stellarMassTableMassFunction(iMass)=Integrate(haloMassMinimum,haloMassMaximum,Stellar_Mass_Function_Integrand,parameterPointer&
          &,integrandFunction ,integrationWorkspace,toleranceAbsolute=toleranceAbsolute,toleranceRelative=toleranceRelative,reset&
          &=integrationReset)

  end do
  call Integrate_Done(integrandFunction,integrationWorkspace)
  integrationReset=.true.

  ! Store stellar mass function data to file.
  call outputFile%writeDataset(stellarMassTableMass        ,'stellarMass'        ,'Stellar mass'         ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(massSolar          ,'unitsInSI')
  call thisDataset%close()
  call outputFile%writeDataset(stellarMassTableMassFunction,'stellarMassFunction','Stellar mass function',datasetReturned=thisDataset)
  call thisDataset%writeAttribute(1.0d0/megaParsec**3,'unitsInSI')
  call thisDataset%close()

   ! Open a file for the tree weighting function data.
   open(newunit=iUnit,file='treeSampling.data',status='unknown',form='formatted')

   ! Loop over halo masses.
   do iMass=1,haloMassTableCount

      ! Set the halo mass to use.
      haloMass=haloMassTableMass(iMass)

      ! Get the sampling density function.
      samplingDensity(iMass)=Merger_Tree_Construct_Mass_Function_Sampling(haloMass,time,haloMassMinimum,haloMassMaximum)

      ! Get the halo mass function.
      haloMassFunction(iMass)=Halo_Mass_Function_Differential(time,haloMass)*haloMass

      ! Get the tree processing time.
      treeTiming(iMass)=Galacticus_Time_Per_Tree(haloMass)

   end do
   close(iUnit)

  ! Store tree mass sampling function data to file.
  call outputFile%writeDataset(haloMassTableMass,'haloMass'        ,'Dark matter halo mass'         ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(massSolar          ,'unitsInSI')
  call thisDataset%close()
  call outputFile%writeDataset(samplingDensity  ,'samplingDensity' ,'Halo mass sampling density')
  call outputFile%writeDataset(haloMassFunction ,'haloMassFunction','Dark matter halo mass function',datasetReturned=thisDataset)
  call thisDataset%writeAttribute(1.0d0/megaParsec**3,'unitsInSI')
  call thisDataset%close()
  call outputFile%writeDataset(treeTiming       ,'treeTiming'      ,'Time to process tree'          ,datasetReturned=thisDataset)
  call thisDataset%writeAttribute(1.0d0              ,'unitsInSI')
  call thisDataset%close()

   ! Close the output file.
   call outputFile%close()

  ! Close the parameter file.
  call Input_Parameters_File_Close

contains

  function Stellar_Mass_Function_Integrand(mass,parameterPointer) bind(c)
    !% The integrand (as a function of halo mass) giving the stellar mass function.
    use Conditional_Mass_Functions
    implicit none
    real            (kind=c_double)            :: Stellar_Mass_Function_Integrand
    real            (kind=c_double), value     :: mass
    type            (c_ptr        ), value     :: parameterPointer
    double precision               , parameter :: deltaLogMass=0.097d0
    double precision                           :: conditionalStellarMassFunction

    conditionalStellarMassFunction=                                                                                       &
         & (                                                                                                              &
         &   Cumulative_Conditional_Mass_Function(mass,stellarMass*(10.0d0**(-0.5d0*optimalSamplingLogarithmicBinWidth))) &
         &  -Cumulative_Conditional_Mass_Function(mass,stellarMass*(10.0d0**(+0.5d0*optimalSamplingLogarithmicBinWidth))) &
         & )                                                                                                              &
         & /(10.0d0**(+0.5d0*optimalSamplingLogarithmicBinWidth)-10.0d0**(-0.5d0*optimalSamplingLogarithmicBinWidth))

    Stellar_Mass_Function_Integrand=Halo_Mass_Function_Differential(time,mass)*conditionalStellarMassFunction
    return
  end function Stellar_Mass_Function_Integrand

end program Optimal_Sampling_SMF
