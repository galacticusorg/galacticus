!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu

!% Contains a module which implements calculations of halo mass functions and related properties for output.

module Halo_Mass_Function_Tasks
  !% Implements calculations of halo mass functions and related properties for output.
  use IO_HDF5
  private
  public :: Halo_Mass_Function_Compute, Halo_Mass_Function_Open_File, Halo_Mass_Function_Close_File, Halo_Mass_Function_Output
  
  ! HDF5 object for the output file.
  type(hdf5Object), public :: haloMassFunctionOutputFile

  ! Arrays of halo mass function data.
  double precision, allocatable, dimension(:,:) :: haloMassFunction_Mass,haloMassFunction_dndM,haloMassFunction_dndlnM &
       &,haloMassFunction_bias,haloMassFunction_sigma,haloMassFunction_nu,haloMassFunction_virialTemperature&
       &,haloMassFunction_virialVelocity ,haloMassFunction_virialRadius,haloMassFunction_cumulative,haloMassFunction_massFraction
  
  ! Arrays of output time data.
  double precision, allocatable, dimension(:  ) :: outputRedshifts,outputExpansionFactors,outputTimes,outputGrowthFactors &
       &,outputCriticalOverdensities,outputVirialDensityContrast,outputCharacteristicMass

  ! The upper limit to halo mass used when computing cumulative mass functions.
  double precision, parameter                   :: haloMassEffectiveInfinity=1.0d16

contains
  
  subroutine Halo_Mass_Function_Open_File(outputFileName)
    !% Open the output file for halo mass function data.
    use ISO_Varying_String
    use HDF5
    implicit none
    type(varying_string), intent(in) :: outputFileName

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
    use Tree_Nodes
    use Merger_Trees
    use Halo_Mass_Function
    use Dark_Matter_Halo_Biases
    use Memory_Management
    use Numerical_Ranges
    use Input_Parameters
    use Critical_Overdensity
    use CDM_Power_Spectrum
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions
    use Linear_Growth
    use Virial_Density_Contrast
    use Galacticus_Calculations_Resets
    implicit none
    integer          :: haloMassFunctionsPointsPerDecade,haloMassFunctionsCount,iMass,outputCount,iOutput
    double precision :: haloMassFunctionsMassMinimum,haloMassFunctionsMassMaximum,criticalOverdensity,characteristicMass
    type(mergerTree) :: thisTree

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
    !@ </inputParameter>
    if (outputCount == 1) then
       ! If only one (or zero) output redshifts present, make redshift zero the default.
       call Get_Input_Parameter('outputRedshifts',outputRedshifts,defaultValue=[0.0d0])
    else
       call Get_Input_Parameter('outputRedshifts',outputRedshifts                     )
    end if
    
    ! Compute output time properties.
    do iOutput=1,outputCount
       outputExpansionFactors     (iOutput)=Expansion_Factor_from_Redshift      (outputRedshifts       (iOutput))
       outputTimes                (iOutput)=Cosmology_Age                       (outputExpansionFactors(iOutput))
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
    !@ </inputParameter>
    call Get_Input_Parameter('haloMassFunctionsMassMinimum',haloMassFunctionsMassMinimum,defaultValue=1.0d10)
    !@ <inputParameter>
    !@   <name>haloMassFunctionsMassMaximum</name>
    !@   <defaultValue>$10^{15}$</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The maximum mass at which to tabulate halo mass functions.
    !@   </description>
    !@ </inputParameter>
    call Get_Input_Parameter('haloMassFunctionsMassMaximum',haloMassFunctionsMassMaximum,defaultValue=1.0d15)
    !@ <inputParameter>
    !@   <name>haloMassFunctionsPointsPerDecade</name>
    !@   <defaultValue>10</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The number of points per decade of halo mass at which to tabulate halo mass functions.
    !@   </description>
    !@ </inputParameter>
    call Get_Input_Parameter('haloMassFunctionsPointsPerDecade',haloMassFunctionsPointsPerDecade,defaultValue=10)

    ! Compute number of tabulation points.
    haloMassFunctionsCount=int(dlog10(haloMassFunctionsMassMaximum/haloMassFunctionsMassMinimum)*dble(haloMassFunctionsPointsPerDecade))+1

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
    call thisTree%createNode(thisTree%baseNode)

    ! Loop over all output times.
    do iOutput=1,outputCount
       ! Set the time in the node.
       call Tree_Node_Time_Set(thisTree%baseNode,outputTimes(iOutput))
       
       ! Build a range of halo masses.
       haloMassFunction_Mass(:,iOutput)=Make_Range(haloMassFunctionsMassMinimum,haloMassFunctionsMassMaximum,haloMassFunctionsCount,rangeTypeLogarithmic)
       
       ! Loop over all halo masses.
       do iMass=1,haloMassFunctionsCount
          ! Reset calculations.
          call Galacticus_Calculations_Reset(thisTree%baseNode)
          ! Set the mass in the node.
          call Tree_Node_Mass_Set(thisTree%baseNode,haloMassFunction_Mass(iMass,iOutput))
          ! Compute halo properties.
          haloMassFunction_dndM             (iMass,iOutput)=Halo_Mass_Function_Differential(outputTimes(iOutput),haloMassFunction_Mass(iMass,iOutput))
          haloMassFunction_dndlnM           (iMass,iOutput)=haloMassFunction_dndM(iMass,iOutput)*haloMassFunction_Mass(iMass,iOutput)
          haloMassFunction_cumulative       (iMass,iOutput)=Halo_Mass_Function_Integrated(outputTimes(iOutput)&
               &,haloMassFunction_Mass(iMass,iOutput),haloMassEffectiveInfinity)
          haloMassFunction_massFraction     (iMass,iOutput)=Halo_Mass_Fraction_Integrated(outputTimes(iOutput)&
               &,haloMassFunction_Mass(iMass,iOutput),haloMassEffectiveInfinity)
          haloMassFunction_sigma            (iMass,iOutput)=sigma_CDM(haloMassFunction_Mass(iMass,iOutput))
          haloMassFunction_nu               (iMass,iOutput)=outputCriticalOverdensities(iOutput)/haloMassFunction_sigma(iMass,iOutput)
          haloMassFunction_bias             (iMass,iOutput)=Dark_Matter_Halo_Bias              (thisTree%baseNode)
          haloMassFunction_virialVelocity   (iMass,iOutput)=Dark_Matter_Halo_Virial_Velocity   (thisTree%baseNode)
          haloMassFunction_virialTemperature(iMass,iOutput)=Dark_Matter_Halo_Virial_Temperature(thisTree%baseNode)
          haloMassFunction_virialRadius     (iMass,iOutput)=Dark_Matter_Halo_Virial_Radius     (thisTree%baseNode)
       end do
       
    end do

    return
  end subroutine Halo_Mass_Function_Compute

  subroutine Halo_Mass_Function_Output
    !% Outputs halo mass function data.
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    implicit none
    type(hdf5Object) :: outputsGroup,thisDataset,massFunctionGroup
    
    ! Open the group for output time information.
    outputsGroup=IO_HDF5_Open_Group(haloMassFunctionOutputFile,'Outputs','Group containing datasets relating to output times.')
    
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
    massFunctionGroup=IO_HDF5_Open_Group(haloMassFunctionOutputFile,'haloMassFunctions','Group containing datasets relating to&
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
