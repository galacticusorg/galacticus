!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements calculations of halo mass functions and related properties for output.

module Halo_Mass_Functions_Tasks
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
       &                                                       haloMassFunction_virialVelocity       , haloMassFunction_scaleRadius      , &
       &                                                       haloMassFunction_velocityMaximum      , haloMassFunction_nuFnu            , &
       &                                                       haloMassFunction_alpha                , haloMassFunction_cumulativeSubhalo

  ! Arrays of output time data.
  double precision            , allocatable, dimension(:  ) :: outputCharacteristicMass              , outputCriticalOverdensities       , &
       &                                                       outputExpansionFactors                , outputGrowthFactors               , &
       &                                                       outputRedshifts                       , outputTimes                       , &
       &                                                       outputVirialDensityContrast

  ! The upper limit to halo mass used when computing cumulative mass functions.
  double precision            , parameter                   :: haloMassEffectiveInfinity      =1.0d16

  ! Largest subhalo mass (in units of host mass) for which we expect significant unevolved subhalo mass function.
  double precision            , parameter                   :: subhaloMassMaximum             =1.0d02

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
    use Node_Components
    use Halo_Mass_Functions
    use Unevolved_Subhalo_Mass_Functions
    use Dark_Matter_Halo_Biases
    use Dark_Matter_Profiles
    use Dark_Matter_Profile_Scales
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Integration
    use Input_Parameters
    use Critical_Overdensities
    use Cosmological_Mass_Variance
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions
    use Linear_Growth
    use Virial_Density_Contrast
    use Galacticus_Display
    use Galacticus_Calculations_Resets
    use ISO_Varying_String
    use Cosmology_Parameters
    implicit none
    class           (nodeComponentBasic               ), pointer :: thisBasic
    class           (nodeComponentDarkMatterProfile   ), pointer :: thisDarkMatterProfile
    type            (treeNode                         ), pointer :: thisNode
    class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_
    class           (cosmologyParametersClass      ), pointer :: cosmologyParameters_
    class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_
    class           (darkMatterProfileClass           ), pointer :: darkMatterProfile_
    class           (virialDensityContrastClass       ), pointer :: virialDensityContrast_
    class           (criticalOverdensityClass         ), pointer :: criticalOverdensity_
    class           (linearGrowthClass                ), pointer :: linearGrowth_
    class           (cosmologicalMassVarianceClass    ), pointer :: cosmologicalMassVariance_
    class           (haloMassFunctionClass            ), pointer :: haloMassFunction_
    class           (unevolvedSubhaloMassFunctionClass), pointer :: unevolvedSubhaloMassFunction_
    type            (fgsl_function                    )          :: integrandFunction
    type            (fgsl_integration_workspace       )          :: integrationWorkspace
    integer                                                      :: haloMassFunctionsCount       , haloMassFunctionsPointsPerDecade, &
         &                                                          iMass                        , iOutput                         , &
         &                                                          outputCount                  , verbosityLevel
    double precision                                             :: haloMassFunctionsMassMaximum , haloMassFunctionsMassMinimum

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
    
    ! Initialize nodes and components.
    call Galacticus_Nodes_Initialize()
    call Node_Components_Initialize ()

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
    ! Get required objects.
    cosmologyParameters_   => cosmologyParameters  ()
    cosmologyFunctions_           => cosmologyFunctions          ()
    virialDensityContrast_        => virialDensityContrast       ()
    darkMatterProfile_            => darkMatterProfile           ()
    criticalOverdensity_          => criticalOverdensity         ()
    linearGrowth_                 => linearGrowth                ()
    haloMassFunction_             => haloMassFunction            ()
    unevolvedSubhaloMassFunction_ => unevolvedSubhaloMassFunction()
    
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

    ! Compute output time properties.
    do iOutput=1,outputCount
       outputExpansionFactors     (iOutput)=cosmologyFunctions_   %expansionFactorFromRedshift(                             outputRedshifts       (iOutput))
       outputTimes                (iOutput)=cosmologyFunctions_   %cosmicTime                 (                             outputExpansionFactors(iOutput))
       outputGrowthFactors        (iOutput)=linearGrowth_         %value                      (                             outputTimes           (iOutput))
       outputCriticalOverdensities(iOutput)=criticalOverdensity_  %value                      (                             outputTimes           (iOutput))
       outputVirialDensityContrast(iOutput)=virialDensityContrast_%densityContrast            (haloMassFunctionsMassMinimum,outputTimes           (iOutput))
       outputCharacteristicMass   (iOutput)=criticalOverdensity_  %collapsingMass             (                             outputTimes           (iOutput))
    end do

    ! Allocate arrays for halo mass functions.
    call Alloc_Array(haloMassFunction_Mass             ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_dndM             ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_dndlnM           ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_cumulative       ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_cumulativeSubhalo,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_massFraction     ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_bias             ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_sigma            ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_alpha            ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_nu               ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_nuFnu            ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_virialVelocity   ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_virialTemperature,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_virialRadius     ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_scaleRadius      ,[haloMassFunctionsCount,outputCount])
    call Alloc_Array(haloMassFunction_velocityMaximum  ,[haloMassFunctionsCount,outputCount])

    ! Get required objects.
    darkMatterHaloScale_      => darkMatterHaloScale     ()
    virialDensityContrast_    => virialDensityContrast   ()
    cosmologicalMassVariance_ => cosmologicalMassVariance()

    ! Create a node object.
    thisNode => treeNode()

    ! Get the basic and dark matter profile components.
    thisBasic             => thisNode%basic            (autoCreate=.true.)
    thisDarkMatterProfile => thisNode%darkMatterProfile(autoCreate=.true.)

    ! Loop over all output times.
    do iOutput=1,outputCount

       ! Set the time in the node.
       call thisBasic%timeSet(outputTimes(iOutput))

       ! Build a range of halo masses.
       haloMassFunction_Mass(:,iOutput)=Make_Range(haloMassFunctionsMassMinimum,haloMassFunctionsMassMaximum,haloMassFunctionsCount,rangeTypeLogarithmic)
       
       ! Loop over all halo masses.
       do iMass=1,haloMassFunctionsCount

          ! Reset calculations.
          call Galacticus_Calculations_Reset(thisNode)
          ! Set the mass in the node.
          call thisBasic            %massSet (haloMassFunction_Mass(iMass,iOutput))
          ! Set the node scale radius.
          call thisDarkMatterProfile%scaleSet(Dark_Matter_Profile_Scale(thisNode))          
          ! Compute halo properties.
          haloMassFunction_dndM             (iMass,iOutput)=haloMassFunction_%differential(outputTimes(iOutput),haloMassFunction_Mass(iMass,iOutput))
          haloMassFunction_dndlnM           (iMass,iOutput)=haloMassFunction_dndM(iMass,iOutput)*haloMassFunction_Mass(iMass,iOutput)
          haloMassFunction_cumulative       (iMass,iOutput)=haloMassFunction_%integrated  (outputTimes(iOutput),haloMassFunction_Mass(iMass,iOutput),haloMassEffectiveInfinity)
          haloMassFunction_massFraction     (iMass,iOutput)=haloMassFunction_%massFraction(outputTimes(iOutput),haloMassFunction_Mass(iMass,iOutput),haloMassEffectiveInfinity)
          haloMassFunction_sigma            (iMass,iOutput)=cosmologicalMassVariance_%rootVariance                   (haloMassFunction_Mass(iMass,iOutput))
          haloMassFunction_alpha            (iMass,iOutput)=cosmologicalMassVariance_%rootVarianceLogarithmicGradient(haloMassFunction_Mass(iMass,iOutput))
          haloMassFunction_nu               (iMass,iOutput)=outputCriticalOverdensities(iOutput)/haloMassFunction_sigma(iMass,iOutput)
          haloMassFunction_nuFnu            (iMass,iOutput)=+haloMassFunction_Mass                                                              (iMass,iOutput)**2 &
               &                                            *haloMassFunction_dndM                                                              (iMass,iOutput)    &
               &                                            /cosmologyParameters_         %densityCritical                ()                                       &
               &                                            /cosmologyParameters_         %OmegaMatter                    ()                                       &
               &                                            /abs(                                                                                                  &
               &                                                 cosmologicalMassVariance_%rootVarianceLogarithmicGradient(                                        &
               &                                                                                                           haloMassFunction_Mass(iMass,iOutput)    &
               &                                                                                                          )                                        &
               &                                                )
          haloMassFunction_bias             (iMass,iOutput)=Dark_Matter_Halo_Bias                        (thisNode)
          haloMassFunction_virialVelocity   (iMass,iOutput)=darkMatterHaloScale_ %virialVelocity         (thisNode)
          haloMassFunction_virialTemperature(iMass,iOutput)=darkMatterHaloScale_ %virialTemperature      (thisNode)
          haloMassFunction_virialRadius     (iMass,iOutput)=darkMatterHaloScale_ %virialRadius           (thisNode)
          haloMassFunction_scaleRadius      (iMass,iOutput)=thisDarkMatterProfile%scale                  (        )
          haloMassFunction_velocityMaximum  (iMass,iOutput)=darkMatterProfile_   %circularVelocityMaximum(thisNode)
          ! Integrate the unevolved subhalo mass function over the halo mass function to get the total subhalo mass function.
          haloMassFunction_cumulativeSubhalo(iMass,iOutput)=Integrate(                                          &
               &                                                      log(haloMassFunction_Mass(1,iOutput)/subhaloMassMaximum       ), &
               &                                                      log(                                 haloMassEffectiveInfinity), &
               &                                                      subhaloMassFunctionIntegrand            , &
               &                                                      integrandFunction                       , &
               &                                                      integrationWorkspace                    , &
               &                                                      toleranceAbsolute    =0.0d+0            , &
               &                                                      toleranceRelative    =1.0d-4            , &
               &                                                      integrationRule      =FGSL_Integ_Gauss15  &
               &                                                     )
          call Integrate_Done(integrandFunction,integrationWorkspace)
       end do

    end do

    return

  contains
    
    double precision function subhaloMassFunctionIntegrand(logMass)
      !% Integrand function used to find the cumulative subhalo mass function.
      implicit none
      double precision, intent(in   ) :: logMass
      double precision                :: mass

      ! Extract integrand parameters.
      mass=exp(logMass)
      ! Return the differential halo mass function multiplied by the integrated unevolved subhalo mass function in such hosts.
      subhaloMassFunctionIntegrand=+                                                                                                                               mass  &
           &                       *haloMassFunction_            %differential(outputTimes(iOutput)                                                               ,mass) &
           &                       *unevolvedSubhaloMassFunction_%integrated  (outputTimes(iOutput),haloMassFunction_Mass(iMass,iOutput),haloMassEffectiveInfinity,mass)
      return
    end function subhaloMassFunctionIntegrand
  
  end subroutine Halo_Mass_Function_Compute

  subroutine Halo_Mass_Function_Output
    !% Outputs halo mass function data.
    use Cosmology_Parameters
    use Numerical_Constants_Astronomical
    implicit none
    class(cosmologyParametersClass), pointer :: cosmologyParameters_
    type (hdf5Object              )          :: massFunctionGroup   , outputsGroup, thisDataset

    ! Get required objects.
    cosmologyParameters_ => cosmologyParameters()
    
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

    ! Store other usual information.
    call massFunctionGroup%writeAttribute(cosmologyParameters_%densityCritical(),'densityCritical')
    
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
    call massFunctionGroup%writeDataset(haloMassFunction_cumulativeSubhalo,'subhaloMassFunctionCumulative','The subhalo cumulative mass&
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
    call massFunctionGroup%writeDataset(haloMassFunction_scaleRadius,'haloScaleRadius','The scale radius of halos.'&
         &,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(megaParsec,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_velocityMaximum,'haloVelocityMaximum','The maximum circular velocity of halos.' &
         &,datasetReturned=thisDataset)
    call thisDataset%writeAttribute(kilo,'unitsInSI')
    call thisDataset%close()
    call massFunctionGroup%writeDataset(haloMassFunction_bias,'haloBias','The large scale linear bias of halos.')
    call massFunctionGroup%writeDataset(haloMassFunction_sigma,'haloSigma','The mass fluctuation on the scale of the halo.')
    call massFunctionGroup%writeDataset(haloMassFunction_alpha,'haloAlpha','dlog(sigma)/dlog(m).')
    call massFunctionGroup%writeDataset(haloMassFunction_nu,'haloNu','The peak height of the halo.')
    call massFunctionGroup%writeDataset(haloMassFunction_massFraction,'haloMassFractionCumulative','The halo cumulative mass fraction.')
    call massFunctionGroup%writeDataset(haloMassFunction_nuFnu,'haloMassFunctionNuFNu','The halo mass fraction function as a function of the nu parameter.')


    ! Close the datasets group.
    call massFunctionGroup%close()

    return
  end subroutine Halo_Mass_Function_Output

end module Halo_Mass_Functions_Tasks
