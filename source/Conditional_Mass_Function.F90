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

!% Contains a program which computes the conditional mass function in bins of mass for a fixed halo mass for use
!% in calculation of constraints.

program Conditional_Mass_Function
  !% Computes the conditional mass function in bins of mass for a fixed halo mass for use in calculation of constraints.
  use, intrinsic :: ISO_C_Binding
  use FGSL
  use IO_HDF5
  use Memory_Management
  use Numerical_Ranges
  use Numerical_Integration
  use Input_Parameters
  use ISO_Varying_String
  use Conditional_Mass_Functions
  use Galacticus_Error
  use Command_Arguments
  use Cosmology_Functions
  use Geometry_Surveys
  implicit none
  double precision                              , allocatable, dimension(:) :: mass,thisConditionalMassFunction
  class           (cosmologyFunctionsClass     ), pointer                   :: cosmologyFunctions_
  class           (conditionalMassFunctionClass), pointer                   :: conditionalMassFunction_
  class           (surveyGeometryClass         ), pointer                   :: surveyGeometry_
  type            (varying_string              )                            :: parameterFile                         , conditionalMassFunctionOutputFileName
  character       (len=32                      )                            :: conditionalMassFunctionHaloMassText
  integer                                                                   :: conditionalMassFunctionMassCount      , iMass                                 , &
       &                                                                       fieldCount                            , iField
  logical                                                                   :: conditionalMassFunctionUseSurveyLimits, integrateOverHaloMassFunction         , &
       &                                                                       integrationReset=.true.               , integrationResetNormalization=.true.
  double precision                                                          :: conditionalMassFunctionHaloMass       , conditionalMassFunctionMassMinimum    , &
       &                                                                       conditionalMassFunctionMassMaximum    , conditionalMassFunctionRedshiftMinimum, &
       &                                                                       conditionalMassFunctionRedshiftMaximum, conditionalMassFunctionHaloMassMinimum, &
       &                                                                       conditionalMassFunctionHaloMassMaximum, massLogarithmDelta                    , &
       &                                                                       massBinMinimum                        , massBinMaximum                        , &
       &                                                                       timeMinimum                           , timeMaximum                           , &
       &                                                                       binTimeMinimum                        , binTimeMaximum                        , &
       &                                                                       time                                  , logHaloMassLower                      , &
       &                                                                       logHaloMassUpper                      , distanceMaximum                       , &
       &                                                                       massFunctionIntegrand                 , volumeIntegrand
  type     (fgsl_function             )                                     :: integrandFunction                     , integrandFunctionNormalization
  type     (fgsl_integration_workspace)                                     :: integrationWorkspace                  , integrationWorkspaceNormalization
  type     (c_ptr                     )                                     :: parameterPointer
  type     (hdf5Object                )                                     :: outputFile

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Conditional_Mass_Function.size')

  ! Check that correct number of arguments have been supplied.
  if (Command_Argument_Count() /= 1) call Galacticus_Error_Report(message="Usage: Conditional_Mass_Function.exe <parameterFile>")

  ! Get the name of the parameter file from the first command line argument.
  call Get_Argument(1,parameterFile   )

  ! Open the parameter file.
  call Input_Parameters_File_Open(parameterFile)

  ! Read parameters controlling the calculation.
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionOutputFileName</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The name of the file to which the computed conditional mass function should be output.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionOutputFileName',conditionalMassFunctionOutputFileName)
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionHaloMass</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>all</defaultValue>
  !@   <description>
  !@     The halo mass for which to compute the conditional mass function. A value of ``all'' will cause the conditional mass function to be integrated over the halo mass function, giving the mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionHaloMass',conditionalMassFunctionHaloMassText,defaultValue="all")
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionMassMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^8M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum mass for which to compute the conditional mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionMassMinimum',conditionalMassFunctionMassMinimum,defaultValue=1.0d8)
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionMassMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^{12}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum mass for which to compute the conditional mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionMassMaximum',conditionalMassFunctionMassMaximum,defaultValue=1.0d12)
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionMassCount</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>21</defaultValue>
  !@   <description>
  !@     The number of bins for which to compute the conditional mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionMassCount',conditionalMassFunctionMassCount,defaultValue=21)
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionRedshiftMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>0</defaultValue>
  !@   <description>
  !@     The minimum redshift for which to compute the conditional mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionRedshiftMinimum',conditionalMassFunctionRedshiftMinimum,defaultValue=0.0d0)
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionRedshiftMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>0</defaultValue>
  !@   <description>
  !@     The maximum redshift for which to compute the conditional mass function.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionRedshiftMaximum',conditionalMassFunctionRedshiftMaximum,defaultValue=0.0d0)
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionUseSurveyLimits</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>false</defaultValue>
  !@   <description>
  !@     Specifies whether the limiting redshifts for integrating over the halo mass function should be limited by those of a galaxy survey.
  !@   </description>
  !@   <type>string</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionUseSurveyLimits',conditionalMassFunctionUseSurveyLimits,defaultValue=.false.)
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionHaloMassMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^6M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum halo mass to use when integrating over the halo mass function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionHaloMassMinimum',conditionalMassFunctionHaloMassMinimum,defaultValue=1.0d6)
  !@ <inputParameter>
  !@   <name>conditionalMassFunctionHaloMassMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
  !@   <description>
  !@     The maximum halo mass to use when integrating over the halo mass function.
  !@   </description>
  !@   <type>real</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('conditionalMassFunctionHaloMassMaximum',conditionalMassFunctionHaloMassMaximum,defaultValue=1.0d16)
  
  ! Get the default cosmology functions, conditional mass function, and survey geometry objects.
  cosmologyFunctions_      => cosmologyFunctions     ()     
  conditionalMassFunction_ => conditionalMassFunction()
  surveyGeometry_          => surveyGeometry         ()

  ! Decode the halo mass parameter.
  integrateOverHaloMassFunction=(conditionalMassFunctionHaloMassText == "all")
  if (.not.integrateOverHaloMassFunction) read (conditionalMassFunctionHaloMassText,*) conditionalMassFunctionHaloMass

  ! Compute the time corresponding to the specified redshift.
  timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(conditionalMassFunctionRedshiftMaximum))
  timeMaximum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(conditionalMassFunctionRedshiftMinimum))

  ! Find logarithmic limits for halo mass in integrations.
  logHaloMassLower=log10(conditionalMassFunctionHaloMassMinimum)
  logHaloMassUpper=log10(conditionalMassFunctionHaloMassMaximum)

  ! Compute a range of masses.
  call Alloc_Array(mass                       ,[conditionalMassFunctionMassCount])
  call Alloc_Array(thisConditionalMassFunction,[conditionalMassFunctionMassCount])
  mass              =Make_Range(conditionalMassFunctionMassMinimum,conditionalMassFunctionMassMaximum,conditionalMassFunctionMassCount,rangeType=rangeTypeLogarithmic)
  massLogarithmDelta=log(mass(2)/mass(1))

  ! Compute the conditional mass function and output to file.
  do iMass=1,conditionalMassFunctionMassCount
     massBinMinimum=exp(log(mass(iMass))-0.5d0*massLogarithmDelta)
     massBinMaximum=exp(log(mass(iMass))+0.5d0*massLogarithmDelta)
     ! Branch on whether the conditional mass function is to be integrated over the halo mass function.
     if (integrateOverHaloMassFunction) then
        ! Branch on whether a range of redshifts was given.
        if (conditionalMassFunctionRedshiftMaximum <= conditionalMassFunctionRedshiftMinimum) then
           ! No range of redshifts given. Compute the mass function at the minimum redshift.
           time=timeMaximum
           thisConditionalMassFunction(iMass)=Integrate(                                    &
                &                                       logHaloMassLower                  , &
                &                                       logHaloMassUpper                  , &
                &                                       Mass_Function_Halo_Mass_Integrand , &
                &                                       parameterPointer                  , &
                &                                       integrandFunction                 , &
                &                                       integrationWorkspace              , &
                &                                       toleranceRelative=1.0d-3          , &
                &                                       reset=integrationReset              &
                &                                      )
        else
           ! Determine number of fields to integrate over.
           if (conditionalMassFunctionUseSurveyLimits) then
              fieldCount=surveyGeometry_%fieldCount()
           else
              fieldCount=1
           end if
           massFunctionIntegrand=0.0d0
           volumeIntegrand      =0.0d0
           do iField=1,fieldCount
              if (conditionalMassFunctionUseSurveyLimits) then
                 ! A survey geometry is imposed. Find the maximum distance at which a galaxy of the present
                 ! mass can be detected in this survey.
                 distanceMaximum=surveyGeometry_%distanceMaximum(sqrt(massBinMinimum*massBinMaximum),iField)
                 ! Set integration limits appropriately.
                 binTimeMinimum=max(timeMinimum,cosmologyFunctions_%timeAtDistanceComoving(distanceMaximum))
                 binTimeMaximum=timeMaximum
              else
                 ! No survey geometry is imposed, so use the full range of specified redshifts.
                 binTimeMinimum=timeMinimum
                 binTimeMaximum=timeMaximum
              end if
              ! Range of redshifts was given, integrate the mass function over this time interval.
              massFunctionIntegrand=                                                        &
                   &                +massFunctionIntegrand                                  &
                   &                +surveyGeometry_%solidAngle(iField)                     &
                   &                *Integrate(                                             &
                   &                           binTimeMinimum                             , &
                   &                           binTimeMaximum                             , &
                   &                           Mass_Function_Time_Integrand               , &
                   &                           parameterPointer                           , &
                   &                           integrandFunction                          , &
                   &                           integrationWorkspace                       , &
                   &                           toleranceRelative=1.0d-3                   , &
                   &                           reset=integrationReset                       &
                   &                          )
              volumeIntegrand      =                                                        &
                   &                +volumeIntegrand                                        &
                   &                +surveyGeometry_%solidAngle(iField)                     &
                   &                *Integrate(                                             &
                   &                           binTimeMinimum                             , &
                   &                           binTimeMaximum                             , &
                   &                           Mass_Function_Time_Normalization_Integrand , &
                   &                           parameterPointer                           , &
                   &                           integrandFunctionNormalization             , &
                   &                           integrationWorkspaceNormalization          , &
                   &                           toleranceRelative=1.0d-3                   , &
                   &                           reset=integrationResetNormalization          &
                   &                          )
           end do
           thisConditionalMassFunction(iMass)=massFunctionIntegrand/volumeIntegrand
        end if
     else
        thisConditionalMassFunction(iMass)=                                                                         &
             &                             (                                                                        &
             &                               conditionalMassFunction_%massFunction(                                 &
             &                                                                     conditionalMassFunctionHaloMass, &
             &                                                                     massBinMinimum                   &
             &                                                                    )                                 &
             &                              -conditionalMassFunction_%massFunction(                                 &
             &                                                                     conditionalMassFunctionHaloMass, &
             &                                                                     massBinMaximum                   &
             &                                                                    )                                 &
             &                             )
     end if
  end do
  thisConditionalMassFunction=thisConditionalMassFunction/massLogarithmDelta

  ! Write the data to file.
  call outputFile%openFile(char(conditionalMassFunctionOutputFileName))
  call outputFile%writeDataset(mass,"mass" ,commentText="mass in units of M☉")
  if (integrateOverHaloMassFunction) then
     call outputFile%writeDataset(thisConditionalMassFunction,"massFunction",commentText="mass function in units of Mpc⁻³ per log(mass)")
  else
     call outputFile%writeDataset(thisConditionalMassFunction,"massFunction",commentText="conditional mass function in units of per log(mass)")
  end if
  call outputFile%close()

  ! Close the parameter file.
  call Input_Parameters_File_Close

contains  
  
  function Mass_Function_Time_Integrand(timePrime,parameterPointer) bind(c)
    !% Integral over time.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double                  )        :: Mass_Function_Time_Integrand
    real(c_double                  ), value :: timePrime
    type(c_ptr                     ), value :: parameterPointer
    type(fgsl_function             ), save  :: integrandFunctionTime
    type(fgsl_integration_workspace), save  :: integrationWorkspaceTime
    type(c_ptr                     )        :: parameterPointerTime
    logical                         , save  :: integrationResetTime=.true.

    time=timePrime
    Mass_Function_Time_Integrand= Integrate(                                             &
         &                                           logHaloMassLower                  , &
         &                                           logHaloMassUpper                  , &
         &                                           Mass_Function_Halo_Mass_Integrand , &
         &                                           parameterPointerTime              , &
         &                                           integrandFunctionTime             , &
         &                                           integrationWorkspaceTime          , &
         &                                           toleranceRelative=1.0d-3          , &
         &                                           reset=integrationResetTime          &
         &                                          )                                    &
         &                               *cosmologyFunctions_%comovingVolumeElementTime(timePrime)
    return
  end function Mass_Function_Time_Integrand

  function Mass_Function_Time_Normalization_Integrand(timePrime,parameterPointer) bind(c)
    !% Normalization integral over time.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)        :: Mass_Function_Time_Normalization_Integrand
    real(c_double), value :: timePrime
    type(c_ptr),    value :: parameterPointer

    Mass_Function_Time_Normalization_Integrand=cosmologyFunctions_%comovingVolumeElementTime(timePrime)
    return
  end function Mass_Function_Time_Normalization_Integrand

  function Mass_Function_Halo_Mass_Integrand(logMass,parameterPointer) bind(c)
    !% Integral over halo mass function.
    use, intrinsic :: ISO_C_Binding
    use Halo_Mass_Function
    use Conditional_Mass_Functions
    implicit none
    real(c_double)          :: Mass_Function_Halo_Mass_Integrand
    real(c_double)  , value :: logMass
    type(c_ptr   )  , value :: parameterPointer
    class            (conditionalMassFunctionClass), pointer :: conditionalMassFunction_
    double precision        :: mass

    conditionalMassFunction_ => conditionalMassFunction()
    mass=10.0d0**logMass
    Mass_Function_Halo_Mass_Integrand= Halo_Mass_Function_Differential(time,mass)                           &
         &                                    *                             mass                            &
         &                                    *log(10.0d0)                                                  &
         &                                    *(                                                            &
         &                                      +conditionalMassFunction_%massFunction(mass,massBinMinimum) &
         &                                      -conditionalMassFunction_%massFunction(mass,massBinMaximum) &
         &                                     )
    return
  end function Mass_Function_Halo_Mass_Integrand

end program Conditional_Mass_Function
