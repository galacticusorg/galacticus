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

!% Contains a module which implements calculations of halo mass functions and related properties for output.

module Halo_Mass_Function_Tasks
  !% Implements calculations of halo mass functions and related properties for output.
  use IO_HDF5
  implicit none
  private
  public :: Halo_Mass_Function_Compute, Halo_Mass_Function_Open_File, Halo_Mass_Function_Close_File, Halo_Mass_Function_Output

  ! HDF5 object for the output file.
  type            (hdf5Object), public                      :: haloMassFunctionOutputFile

  ! Arrays of halo mass function data.
  double precision            , allocatable, dimension(:,:) :: haloMassFunction_Mass                 , haloMassFunction_bias             , &
       &                                                       haloMassFunction_cumulative           , haloMassFunction_dndM             , &
       &                                                       haloMassFunction_dndlnM               , haloMassFunction_massFraction     , &
       &                                                       haloMassFunction_nu                   , haloMassFunction_sigma            , &
       &                                                       haloMassFunction_virialRadius         , haloMassFunction_virialTemperature, &
       &                                                       haloMassFunction_virialVelocity

  ! Arrays of output time data.
  double precision            , allocatable, dimension(:  ) :: outputCharacteristicMass              , outputCriticalOverdensities       , &
       &                                                       outputExpansionFactors                , outputGrowthFactors               , &
       &                                                       outputRedshifts                       , outputTimes                       , &
       &                                                       outputVirialDensityContrast

  ! The upper limit to halo mass used when computing cumulative mass functions.
  double precision            , parameter                   :: haloMassEffectiveInfinity      =1.0d16

contains

  subroutine Halo_Mass_Function_Open_File(outputFileName)
    !% Open the output file for halo mass function data.
    use ISO_Varying_String
    use HDF5
    implicit none
    type(varying_string), intent(in   ) :: outputFileName

    ! Open the output file.
    call haloMassFunctionOutputFile%openFile(char(outputFileName),overWrite=.true.,objectsOverwritable=.false.)

    ! Set default chunking and compression levels.
    call IO_HDF5_Set_Defaults(chunkSize=int(128,kind=hsize_t),compressionLevel=9)

    return
  end subroutine Halo_Mass_Function_Open_File

  subroutine Halo_Mass_Function_Close_File
    !% Close the output file for halo mass function data.
    implicit none

    call haloMassFunctionOutputFile%close()
    return
  end subroutine Halo_Mass_Function_Close_File

  subroutine Halo_Mass_Function_Compute
    !% Computes mass functions and related properties for output.
    use Galacticus_Nodes
    use Halo_Mass_Function
    use Dark_Matter_Halo_Biases
    use Memory_Management
    use Numerical_Ranges
    use Input_Parameters
    use Critical_Overdensity
    use Power_Spectra
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions
    use Linear_Growth
    use Virial_Density_Contrast
    use Galacticus_Display
    use Galacticus_Calculations_Resets
    implicit none
    class           (nodeComponentBasic     ), pointer :: thisBasicComponent
    integer                                            :: haloMassFunctionsCount      , haloMassFunctionsPointsPerDecade, &
         &                                                iMass                       , iOutput                         , &
         &                                                outputCount                 , verbosityLevel
    double precision                                   :: haloMassFunctionsMassMaximum, haloMassFunctionsMassMinimum
    type            (treeNode               ), pointer :: thisNode
    class           (cosmologyFunctionsClass), pointer :: cosmologyFunctionsDefault

    ! Get the verbosity level parameter.
    !@ <inputParameter>
    !@   <name>verbosityLevel</name>
    !@   <defaultValue>1</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The level of verbosity for \glc\ (higher values give more verbosity).
    !@   </description>
    !@   <type>integer</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('verbosityLevel',verbosityLevel,1)
    call Galacticus_Verbosity_Level_Set(verbosityLevel)

    ! Get the verbosity level parameter.
    !@ <inputParameter>
    !@   <name>verbosityLevel</name>
    !@   <defaultValue>1</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The level of verbosity for \glc\ (higher values give more verbosity).
    !@   </description>
    !@   <type>integer</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('verbosityLevel',verbosityLevel,1)
    call Galacticus_Verbosity_Level_Set(verbosityLevel)

    ! Get the requested output redshifts.
    outputCount=max(Get_Input_Parameter_Array_Size('outputRedshifts'),1)
    call Alloc_Array(outputTimes                ,[outputCount])
    call Alloc_Array(outputRedshifts            ,[outputCount])
    call Alloc_Array(outputExpansionFactors     ,[outputCount])
    call Alloc_Array(outputGrowthFactors        ,[outputCount])
    call Alloc_Array(outputCriticalOverdensities,[outputCount])
    call Alloc_Array(outputVirialDensityContrast,[outputCount])
    call Alloc_Array(outputCharacteristicMass   ,[outputCount])
    !@ <inputParameter>
    !@   <name>outputRedshifts</name>
    !@   <defaultValue>0</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     A list of redshifts at which halo mass functions should be computed.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    if (outputCount == 1) then
       ! If only one (or zero) output redshifts present, make redshift zero the default.
       call Get_Input_Parameter('outputRedshifts',outputRedshifts,defaultValue=[0.0d0])
    else
       call Get_Input_Parameter('outputRedshifts',outputRedshifts                     )
    end if
    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Compute output time properties.
    do iOutput=1,outputCount
       outputExpansionFactors     (iOutput)=cosmologyFunctionsDefault%expansionFactorFromRedshift      (outputRedshifts       (iOutput))
       outputTimes                (iOutput)=cosmologyFunctionsDefault%cosmicTime                       (outputExpansionFactors(iOutput))
       outputGrowthFactors        (iOutput)=Linear_Growth_Factor                (outputTimes           (iOutput))
       outputCriticalOverdensities(iOutput)=Critical_Overdensity_for_Collapse   (outputTimes           (iOutput))
       outputVirialDensityContrast(iOutput)=Halo_Virial_Density_Contrast        (outputTimes           (iOutput))
       outputCharacteristicMass   (iOutput)=Critical_Overdensity_Collapsing_Mass(outputTimes           (iOutput))
    end do

    ! Find the mass range and increment size.
    !@ <inputParameter>
    !@   <name>haloMassFunctionsMassMinimum</name>
    !@   <defaultValue>$10^{10}$</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The minimum mass at which to tabulate halo mass functions.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('haloMassFunctionsMassMinimum',haloMassFunctionsMassMinimum,defaultValue=1.0d10)
    !@ <inputParameter>
    !@   <name>haloMassFunctionsMassMaximum</name>
    !@   <defaultValue>$10^{15}$</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The maximum mass at which to tabulate halo mass functions.
    !@   </description>
    !@   <type>real</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('haloMassFunctionsMassMaximum',haloMassFunctionsMassMaximum,defaultValue=1.0d15)
    !@ <inputParameter>
    !@   <name>haloMassFunctionsPointsPerDecade</name>
    !@   <defaultValue>10</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The number of points per decade of halo mass at which to tabulate halo mass functions.
    !@   </description>
    !@   <type>integer</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('haloMassFunctionsPointsPerDecade',haloMassFunctionsPointsPerDecade,defaultValue=10)

    ! Compute number of tabulation points.
    haloMassFunctionsCount=int(log10(haloMassFunctionsMassMaximum/haloMassFunctionsMassMinimum)*dble(haloMassFunctionsPointsPerDecade))+1

    ! Allocate arrays for halo mass functions.
    call Alloc_Array(haloMassFunction_Mass             ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_dndM             ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_dndlnM           ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_cumulative       ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_massFraction     ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_bias             ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_sigma            ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_nu               ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_virialVelocity   ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_virialTemperature,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_virialRadius     ,[haloMassFunctionsCount,outputCount])

    ! Create a node object.
    thisNode => treeNode()

    ! Get the basic component.
    thisBasicComponent => thisNode%basic(autoCreate=.true.)

    ! Loop over all output times.
    do iOutput=1,outputCount
       ! Set the time in the node.
       call thisBasicComponent%timeSet(outputTimes(iOutput))

       ! Build a range of halo masses.
       haloMassFunction_Mass(:,iOutput)=Make_Range(haloMassFunctionsMassMinimum,haloMassFunctionsMassMaximum,haloMassFunctionsCount,rangeTypeLogarithmic)

       ! Loop over all halo masses.
       do iMass=1,haloMassFunctionsCount
          ! Reset calculations.
          call Galacticus_Calculations_Reset(thisNode)
          ! Set the mass in the node.
          call thisBasicComponent%massSet(haloMassFunction_Mass(iMass,iOutput))
          ! Compute halo properties.
          haloMassFunction_dndM             (iMass,iOutput)=Halo_Mass_Function_Differential(outputTimes(iOutput),haloMassFunction_Mass(iMass,iOutput))
          haloMassFunction_dndlnM           (iMass,iOutput)=haloMassFunction_dndM(iMass,iOutput)*haloMassFunction_Mass(iMass,iOutput)
          ! haloMassFunction_cumulative       (iMass,iOutput)=Halo_Mass_Function_Integrated(outputTimes(iOutput)&
          !      & haloMassFunction_Mass(iMass,iOutput),haloMassEffectiveInfinity)
          ! haloMassFunction_massFraction     (iMass,iOutput)=Halo_Mass_Fraction_Integrated(outputTimes(iOutput)&
          !      &,haloMassFunction_Mass(iMass,iOutput),haloMassEffectiveInfinity)
          haloMassFunction_sigma            (iMass,iOutput)=Cosmological_Mass_Root_Variance(haloMassFunction_Mass(iMass,iOutput))
          haloMassFunction_nu               (iMass,iOutput)=outputCriticalOverdensities(iOutput)/haloMassFunction_sigma(iMass,iOutput)
          haloMassFunction_bias             (iMass,iOutput)=Dark_Matter_Halo_Bias              (thisNode)
          haloMassFunction_virialVelocity   (iMass,iOutput)=Dark_Matter_Halo_Virial_Velocity   (thisNode)
          haloMassFunction_virialTemperature(iMass,iOutput)=Dark_Matter_Halo_Virial_Temperature(thisNode)
          haloMassFunction_virialRadius     (iMass,iOutput)=Dark_Matter_Halo_Virial_Radius     (thisNode)
       end do

    end do

    return
  end subroutine Halo_Mass_Function_Compute

  subroutine Halo_Mass_Function_Output
    !% Outputs halo mass function data.
    use Numerical_Constants_Astronomical
    implicit none
    type(hdf5Object) :: massFunctionGroup, outputsGroup, thisDataset

    ! Open the group for output time information.
    outputsGroup=haloMassFunctionOutputFile%openGroup('Outputs','Group containing datasets relating to output times.')

    ! Write output time data.
    call outputsGroup%writeDataset(outputTimes            ,'outputTime'              ,'The time corresponding to each output.'                 &
         & ,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(gigaYear ,'unitsInSI')
    call thisDataset%close()
    call outputsGroup%writeDataset(outputCharacteristicMass,'outputCharacteristicMass','The characteristic mass corresponding to each output.' &
         & ,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(massSolar,'unitsInSI')
    call thisDataset%close()
    call outputsGroup%writeDataset(outputRedshifts            ,'outputRedshift'             ,'The redshift corresponding to each output.'               )
    call outputsGroup%writeDataset(outputExpansionFactors     ,'outputExpansionFactor'      ,'The expansion factor corresponding to each output.'       )
    call outputsGroup%writeDataset(outputGrowthFactors        ,'outputGrowthFactors'        ,'The growth factor corresponding to each output.'          )
    call outputsGroup%writeDataset(outputCriticalOverdensities,'outputCriticalOverdensities','The critical overdensities corresponding to each output.' )
    call outputsGroup%writeDataset(outputVirialDensityContrast,'outputVirialDensityContrast','The virial density contrast corresponding to each output.')

    ! Close the outputs group.
    call outputsGroup%close()

    ! Write mass function datasets.
    massFunctionGroup=haloMassFunctionOutputFile%openGroup('haloMassFunctions','Group containing datasets relating to&
         & halo mass functions.')

    ! Write the halo mass function data.
    call massFunctionGroup%writeDataset(haloMassFunction_Mass,'haloMass','The mass of the halo.',datasetReturned=thisDataset)
    call thisDataset%writeAttribute(massSolar,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_dndM,'haloMassFunctionM','The halo mass function (per unit halo mass).'&
         &,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(1.0d0/megaParsec**3/massSolar,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_dndlnM,'haloMassFunctionLnM','The halo mass function (per logarithmic   &
         &       halo mass).',datasetReturned=thisDataset)
    call thisDataset%writeAttribute(1.0d0/megaParsec**3,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_cumulative,'haloMassFunctionCumulative','The halo cumulative mass&
         & function.',datasetReturned=thisDataset)
    call thisDataset%writeAttribute(1.0d0/megaParsec**3,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_virialVelocity,'haloVirialVelocity','The virial velocity of halos.' &
         &,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(kilo,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_virialTemperature,'haloVirialTemperature','The virial temperature of halos.' &
         &,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(1.0d0,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_virialRadius,'haloVirialRadius','The virial radius of halos.'&
         &,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(megaParsec,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_bias,'haloBias','The large scale linear bias of halos.')
    call massFunctionGroup%writeDataset(haloMassFunction_sigma,'haloSigma','The mass fluctuation on the scale of the halo.')
    call massFunctionGroup%writeDataset(haloMassFunction_nu,'haloNu','The peak height of the halo.')
    call massFunctionGroup%writeDataset(haloMassFunction_massFraction,'haloMassFractionCumulative','The halo cumulative mass fraction.')

    ! Close the datasets group.
    call massFunctionGroup%close()

    return
  end subroutine Halo_Mass_Function_Output

end module Halo_Mass_Function_Tasks
