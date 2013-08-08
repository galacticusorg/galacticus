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

!% Contains a program which computes the conditional stellar mass function in bins of stellar mass for a fixed halo mass for use
!% in calculation of constraints.

program Constraint_CSMF
  !% Computes the conditional stellar mass function in bins of stellar mass for a fixed halo mass for use in calculation of constraints.
  use, intrinsic :: ISO_C_Binding
  use FGSL
  use IO_HDF5
  use Memory_Management
  use Numerical_Ranges
  use Numerical_Integration
  use Input_Parameters
  use ISO_Varying_String
  use Conditional_Stellar_Mass_Functions
  use Galacticus_Error
  use Command_Arguments
  use Cosmology_Functions
  use Geometry_Surveys
  implicit none
  double precision                     , allocatable, dimension(:) :: massStellar,conditionalStellarMassFunction
  type     (varying_string            )                            :: parameterFile,conditionalStellarMassFunctionOutputFileName
  character(len=32                    )                            :: conditionalStellarMassFunctionHaloMassText
  integer                                                          :: conditionalStellarMassFunctionMassCount,iMass
  logical                                                          :: conditionalStellarMassFunctionUseSurveyLimits&
       &,integrateOverHaloMassFunction,integrationReset=.true.,integrationResetNormalization=.true.
  double precision                                                 :: conditionalStellarMassFunctionHaloMass &
       &,conditionalStellarMassFunctionMassMinimum,conditionalStellarMassFunctionMassMaximum &
       &,conditionalStellarMassFunctionRedshiftMinimum,conditionalStellarMassFunctionRedshiftMaximum &
       &,conditionalStellarMassFunctionHaloMassMinimum ,conditionalStellarMassFunctionHaloMassMaximum,massStellarLogarithmDelta &
       &,massBinMinimum,massBinMaximum,timeMinimum,timeMaximum,binTimeMinimum,binTimeMaximum,time,logHaloMassLower&
       &,logHaloMassUpper,distanceMaximum
  type     (fgsl_function             )                            :: integrandFunction,integrandFunctionNormalization
  type     (fgsl_integration_workspace)                            :: integrationWorkspace,integrationWorkspaceNormalization
  type     (c_ptr                     )                            :: parameterPointer
  type     (hdf5Object                )                            :: outputFile

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Conditional_Stellar_Mass_Function.size')

  ! Check that correct number of arguments have been supplied.
  if (Command_Argument_Count() /= 1) call Galacticus_Error_Report(message="Usage: Conditional_Stellar_Mass_Function.exe <parameterFile>")

  ! Get the name of the parameter file from the first command line argument.
  call Get_Argument(1,parameterFile   )

  ! Open the parameter file.
  call Input_Parameters_File_Open(parameterFile)

  ! Read parameters controlling the calculation.
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionOutputFileName</name>
  !@   <attachedTo>module</attachedTo>
  !@   <description>
  !@     The name of the file to which the computed conditional stellar mass function should be output.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionOutputFileName',conditionalStellarMassFunctionOutputFileName)
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionHaloMass</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>all</defaultValue>
  !@   <description>
  !@     The halo mass for which to compute the conditional stellar mass function. A value of ``all'' will cause the conditional stellar mass function to be integrated over the halo mass function, giving the stellar mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionHaloMass',conditionalStellarMassFunctionHaloMassText,defaultValue="all")
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionMassMinimum</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>$10^8M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum stellar mass for which to compute the conditional stellar mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionMassMinimum',conditionalStellarMassFunctionMassMinimum,defaultValue=1.0d8)
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionMassMaximum</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>$10^{12}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum stellar mass for which to compute the conditional stellar mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionMassMaximum',conditionalStellarMassFunctionMassMaximum,defaultValue=1.0d12)
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionMassCount</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>21</defaultValue>
  !@   <description>
  !@     The number of bins for which to compute the conditional stellar mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionMassCount',conditionalStellarMassFunctionMassCount,defaultValue=21)
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionRedshiftMinimum</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>0</defaultValue>
  !@   <description>
  !@     The minimum redshift for which to compute the conditional stellar mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionRedshiftMinimum',conditionalStellarMassFunctionRedshiftMinimum,defaultValue=0.0d0)
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionRedshiftMaximum</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>0</defaultValue>
  !@   <description>
  !@     The maximum redshift for which to compute the conditional stellar mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionRedshiftMaximum',conditionalStellarMassFunctionRedshiftMaximum,defaultValue=0.0d0)
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionUseSurveyLimits</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>false</defaultValue>
  !@   <description>
  !@     Specifies whether the limiting redshifts for integrating over the halo mass function should be limited by those of a galaxy survey.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionUseSurveyLimits',conditionalStellarMassFunctionUseSurveyLimits,defaultValue=.false.)
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionHaloMassMinimum</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>$10^6M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum halo mass to use when integrating over the halo mass function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionHaloMassMinimum',conditionalStellarMassFunctionHaloMassMinimum,defaultValue=1.0d6)
  !@ <inputParameter>
  !@   <name>conditionalStellarMassFunctionHaloMassMaximum</name>
  !@   <attachedTo>module</attachedTo>
  !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum halo mass to use when integrating over the halo mass function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalStellarMassFunctionHaloMassMaximum',conditionalStellarMassFunctionHaloMassMaximum,defaultValue=1.0d16)

  ! Decode the halo mass parameter.
  integrateOverHaloMassFunction=(conditionalStellarMassFunctionHaloMassText == "all")
  if (.not.integrateOverHaloMassFunction) read (conditionalStellarMassFunctionHaloMassText,*) conditionalStellarMassFunctionHaloMass

  ! Compute the time corresponding to the specified redshift.
  timeMinimum=Cosmology_Age(Expansion_Factor_from_Redshift(conditionalStellarMassFunctionRedshiftMaximum))
  timeMaximum=Cosmology_Age(Expansion_Factor_from_Redshift(conditionalStellarMassFunctionRedshiftMinimum))

  ! Find logarithmic limits for halo mass in integrations.
  logHaloMassLower=log10(conditionalStellarMassFunctionHaloMassMinimum)
  logHaloMassUpper=log10(conditionalStellarMassFunctionHaloMassMaximum)

  ! Compute a range of stellar masses.
  call Alloc_Array(massStellar                   ,[conditionalStellarMassFunctionMassCount])
  call Alloc_Array(conditionalStellarMassFunction,[conditionalStellarMassFunctionMassCount])
  massStellar              =Make_Range(conditionalStellarMassFunctionMassMinimum,conditionalStellarMassFunctionMassMaximum,conditionalStellarMassFunctionMassCount,rangeType=rangeTypeLogarithmic)
  massStellarLogarithmDelta=log(massStellar(2)/massStellar(1))

  ! Compute the conditional stellar mass function and output to file.
  do iMass=1,conditionalStellarMassFunctionMassCount
     massBinMinimum=exp(log(massStellar(iMass))-0.5d0*massStellarLogarithmDelta)
     massBinMaximum=exp(log(massStellar(iMass))+0.5d0*massStellarLogarithmDelta)
     ! Branch on whether the conditional mass function is to be integrated over the halo mass function.
     if (integrateOverHaloMassFunction) then
        ! Branch on whether a range of redshifts was given.
        if (conditionalStellarMassFunctionRedshiftMaximum <= conditionalStellarMassFunctionRedshiftMinimum) then
           ! No range of redshifts given. Compute the mass function at the minimum redshift.
           time=timeMaximum
           conditionalStellarMassFunction(iMass)=Integrate(                                            &
                &                                          logHaloMassLower                          , &
                &                                          logHaloMassUpper                          , &
                &                                          Stellar_Mass_Function_Halo_Mass_Integrand , &
                &                                          parameterPointer                          , &
                &                                          integrandFunction                         , &
                &                                          integrationWorkspace                      , &
                &                                          toleranceRelative=1.0d-3                  , &
                &                                          reset=integrationReset                      &
                &                                         )
        else
           if (conditionalStellarMassFunctionUseSurveyLimits) then
              ! A survey geometry is imposed. Find the maximum distance at which a galaxy of the present
              ! mass can be detected in this survey.
              distanceMaximum=Geometry_Survey_Distance_Maximum(sqrt(massBinMinimum*massBinMaximum))
              ! Set integration limits appropriately.
              binTimeMinimum=max(timeMinimum,Time_From_Comoving_Distance(distanceMaximum))
              binTimeMaximum=timeMaximum
           else
              ! No survey geometry is imposed, so use the full range of specified redshifts.
              binTimeMinimum=timeMinimum
              binTimeMaximum=timeMaximum
           end if
           ! Range of redshifts was given, integrate the mass function over this time interval.
           conditionalStellarMassFunction(iMass)= Integrate(                                                     &
                &                                           binTimeMinimum                                     , &
                &                                           binTimeMaximum                                     , &
                &                                           Stellar_Mass_Function_Time_Integrand               , &
                &                                           parameterPointer                                   , &
                &                                           integrandFunction                                  , &
                &                                           integrationWorkspace                               , &
                &                                           toleranceRelative=1.0d-3                           , &
                &                                           reset=integrationReset                               &
                &                                          )                                                     &
                &                                   /Integrate(                                                  &
                &                                           binTimeMinimum                                     , &
                &                                           binTimeMaximum                                     , &
                &                                           Stellar_Mass_Function_Time_Normalization_Integrand , &
                &                                           parameterPointer                                   , &
                &                                           integrandFunctionNormalization                     , &
                &                                           integrationWorkspaceNormalization                  , &
                &                                           toleranceRelative=1.0d-3                           , &
                &                                           reset=integrationResetNormalization                  &
                &                                          )
        end if
     else
        conditionalStellarMassFunction(iMass)=                                                                                       &
             &                                (                                                                                      &
             &                                  Cumulative_Conditional_Stellar_Mass_Function(                                        &
             &                                                                               conditionalStellarMassFunctionHaloMass, &
             &                                                                               massBinMinimum                          &
             &                                                                              )                                        &
             &                                 -Cumulative_Conditional_Stellar_Mass_Function(                                        &
             &                                                                               conditionalStellarMassFunctionHaloMass, &
             &                                                                               massBinMaximum                          &
             &                                                                              )                                        &
             &                                )
     end if
  end do
  conditionalStellarMassFunction=conditionalStellarMassFunction/massStellarLogarithmDelta

  ! Write the data to file.
  call outputFile%openFile(char(conditionalStellarMassFunctionOutputFileName))
  call outputFile%writeDataset(massStellar,"stellarMass" ,commentText="stellar mass in units of M☉")
  if (integrateOverHaloMassFunction) then
     call outputFile%writeDataset(conditionalStellarMassFunction,"massFunction",commentText="stellar mass function in units of Mpc⁻³ per lo(mass)")
  else
     call outputFile%writeDataset(conditionalStellarMassFunction,"massFunction",commentText="conditional stellar mass function in units of per log(mass)")
  end if
  call outputFile%close()

  ! Close the parameter file.
  call Input_Parameters_File_Close

contains  
  
  function Stellar_Mass_Function_Time_Integrand(timePrime,parameterPointer) bind(c)
    !% Integral over time.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double                  )        :: Stellar_Mass_Function_Time_Integrand
    real(c_double                  ), value :: timePrime
    type(c_ptr                     ), value :: parameterPointer
    type(fgsl_function             ), save  :: integrandFunctionTime
    type(fgsl_integration_workspace), save  :: integrationWorkspaceTime
    type(c_ptr                     )        :: parameterPointerTime
    logical                         , save  :: integrationResetTime=.true.

    time=timePrime
    Stellar_Mass_Function_Time_Integrand= Integrate(                                             &
         &                                           logHaloMassLower                          , &
         &                                           logHaloMassUpper                          , &
         &                                           Stellar_Mass_Function_Halo_Mass_Integrand , &
         &                                           parameterPointerTime                      , &
         &                                           integrandFunctionTime                     , &
         &                                           integrationWorkspaceTime                  , &
         &                                           toleranceRelative=1.0d-3                  , &
         &                                           reset=integrationResetTime                  &
         &                                          )                                            &
         &                               *Comoving_Volume_Element_Time(timePrime)
    return
  end function Stellar_Mass_Function_Time_Integrand

  function Stellar_Mass_Function_Time_Normalization_Integrand(timePrime,parameterPointer) bind(c)
    !% Normalization integral over time.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)        :: Stellar_Mass_Function_Time_Normalization_Integrand
    real(c_double), value :: timePrime
    type(c_ptr),    value :: parameterPointer

    Stellar_Mass_Function_Time_Normalization_Integrand=Comoving_Volume_Element_Time(timePrime)
    return
  end function Stellar_Mass_Function_Time_Normalization_Integrand

  function Stellar_Mass_Function_Halo_Mass_Integrand(logMass,parameterPointer) bind(c)
    !% Integral over halo mass function.
    use, intrinsic :: ISO_C_Binding
    use Halo_Mass_Function
    use Conditional_Stellar_Mass_Functions
    implicit none
    real(c_double)          :: Stellar_Mass_Function_Halo_Mass_Integrand
    real(c_double)  , value :: logMass
    type(c_ptr   )  , value :: parameterPointer
    double precision        :: mass

    mass=10.0d0**logMass
    Stellar_Mass_Function_Halo_Mass_Integrand= Halo_Mass_Function_Differential(time,mass)                          &
         &                                    *                                     mass                           &
         &                                    *log(10.0d0)                                                         &
         &                                    *(                                                                   &
         &                                       Cumulative_Conditional_Stellar_Mass_Function(mass,massBinMinimum) &
         &                                      -Cumulative_Conditional_Stellar_Mass_Function(mass,massBinMaximum) &
         &                                     )
    return
  end function Stellar_Mass_Function_Halo_Mass_Integrand

end program Constraint_CSMF
