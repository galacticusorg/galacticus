!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines methods for a mass function likelihood function.

function likelihoodMassFunctionConstructor(haloMassMinimum,haloMassMaximum,redshiftMinimum,redshiftMaximum,useSurveyLimits,massFunctionFileName,modelSurfaceBrightness,surfaceBrightnessLimit)
  !% Constructor for the mass function likelihood class.
  use Galacticus_Input_Paths
  use IO_HDF5
  use Memory_Management
  implicit none
  type            (likelihoodMassFunction)                              :: likelihoodMassFunctionConstructor
  double precision                        , intent(in   )               :: redshiftMinimum                  , redshiftMaximum       , &
       &                                                                   haloMassMinimum                  , haloMassMaximum       , &
       &                                                                   surfaceBrightnessLimit
  logical                                 , intent(in   )               :: useSurveyLimits                  , modelSurfaceBrightness
  character       (len=*                 ), intent(in   )               :: massFunctionFileName
  double precision                        , allocatable  , dimension(:) :: massBinWidth
  type            (hdf5Object            )                              :: massFunctionFile
  integer                                                               :: i

  likelihoodMassFunctionConstructor%useSurveyLimits=useSurveyLimits
  likelihoodMassFunctionConstructor%       haloMassMinimum=      haloMassMinimum
  likelihoodMassFunctionConstructor%       haloMassMaximum=      haloMassMaximum
  likelihoodMassFunctionConstructor%    logHaloMassMinimum=log10(haloMassMinimum)
  likelihoodMassFunctionConstructor%    logHaloMassMaximum=log10(haloMassMaximum)
  likelihoodMassFunctionConstructor%       redshiftMinimum=      redshiftMinimum
  likelihoodMassFunctionConstructor%       redshiftMaximum=      redshiftMaximum
  likelihoodMassFunctionConstructor%modelSurfaceBrightness=      modelSurfaceBrightness
  likelihoodMassFunctionConstructor%surfaceBrightnessLimit=      surfaceBrightnessLimit
  ! Read the mass function file.
  !$omp critical(HDF5_Access)
  call massFunctionFile%openFile(char(Galacticus_Input_Path())//massFunctionFileName,readOnly=.true.)
  call massFunctionFile%readDataset("mass"                ,likelihoodMassFunctionConstructor%mass                )
  call massFunctionFile%readDataset("massFunctionObserved",likelihoodMassFunctionConstructor%massFunctionObserved)
  call massFunctionFile%readDataset("covariance"          ,likelihoodMassFunctionConstructor%covarianceMatrix    )
  ! Compute mass bin limits.
  call Alloc_Array(likelihoodMassFunctionConstructor%massFunction,shape(likelihoodMassFunctionConstructor%mass))
  call Alloc_Array(likelihoodMassFunctionConstructor%massMinimum ,shape(likelihoodMassFunctionConstructor%mass))
  call Alloc_Array(likelihoodMassFunctionConstructor%massMaximum ,shape(likelihoodMassFunctionConstructor%mass))
  if (massFunctionFile%hasDataset("massWidthObserved")) then
     call massFunctionFile%readDataset("massWidthObserved",massBinWidth)
     do i=1,size(likelihoodMassFunctionConstructor%mass)
        likelihoodMassFunctionConstructor%massMinimum(i)=likelihoodMassFunctionConstructor%mass(i)/sqrt(massBinWidth(i))
        likelihoodMassFunctionConstructor%massMaximum(i)=likelihoodMassFunctionConstructor%mass(i)*sqrt(massBinWidth(i))
     end do
     call Dealloc_Array(massBinWidth)
  else
     do i=1,size(likelihoodMassFunctionConstructor%mass)
        if (i ==                                            1) then
           likelihoodMassFunctionConstructor              %massMinimum(i  ) &
                & =      likelihoodMassFunctionConstructor%mass       (i  ) &
                & *sqrt(                                                    &
                &        likelihoodMassFunctionConstructor%mass       (i  ) &
                &       /likelihoodMassFunctionConstructor%mass       (i+1) &
                &      )
        else
           likelihoodMassFunctionConstructor              %massMinimum(i  ) &
                & =sqrt(                                                    &
                &        likelihoodMassFunctionConstructor%mass       (i-1) &
                &       *likelihoodMassFunctionConstructor%mass       (i  ) &
                &      )
        end if
        if (i == size(likelihoodMassFunctionConstructor%mass)) then
           likelihoodMassFunctionConstructor              %massMaximum(i  ) &
                & =      likelihoodMassFunctionConstructor%mass       (i  ) &
                & *sqrt(                                                    &
                &        likelihoodMassFunctionConstructor%mass       (i  ) &
                &       /likelihoodMassFunctionConstructor%mass       (i-1) &
                &      )
        else
           likelihoodMassFunctionConstructor              %massMaximum(i  ) &
                & =sqrt(                                                    &
                &        likelihoodMassFunctionConstructor%mass       (i  ) &
                &       *likelihoodMassFunctionConstructor%mass       (i+1) &
                &      )
        end if
     end do
  end if
  call massFunctionFile%close()
  !$omp end critical(HDF5_Access)
  ! Find the inverse covariance matrix.
  likelihoodMassFunctionConstructor%covariance       =likelihoodMassFunctionConstructor%covarianceMatrix
  likelihoodMassFunctionConstructor%inverseCovariance=likelihoodMassFunctionConstructor%covariance      %invert()
  return
end function likelihoodMassFunctionConstructor

double precision function likelihoodMassFunctionEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return the log-likelihood for the mass function likelihood function.
  use, intrinsic :: ISO_C_Binding
  use Constraints_State
  use Constraints_Constants
  use Constraints_Convergence
  use Conditional_Mass_Functions
  use Cosmology_Functions
  use Geometry_Surveys
  use Numerical_Integration
  use Galacticus_Error
  use Mass_Function_Incompletenesses
  implicit none
  class           (likelihoodMassFunction                     ), intent(inout)               :: self
  class           (state                                      ), intent(in   )               :: simulationState
  type            (mappingList                                ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence                                ), intent(inout)               :: simulationConvergence
  double precision                                             , intent(in   )               :: temperature                    , logLikelihoodCurrent              , &
       &                                                                                        logPriorCurrent                , logPriorProposed
  double precision                                             , allocatable  , dimension(:) :: stateVector
  class           (cosmologyFunctionsClass                    ), pointer                     :: cosmologyFunctions_
  class           (surveyGeometryClass                        ), pointer                     :: surveyGeometry_
  type            (conditionalMassFunctionBehroozi2010        )                              :: thisConditionalMassFunction
  type            (massFunctionIncompletenessSurfaceBrightness)                              :: massFunctionIncompletenessModel
  type            (fgsl_function                              )                              :: integrandFunction               , integrandFunctionNormalization
  type            (fgsl_integration_workspace                 )                              :: integrationWorkspace            , integrationWorkspaceNormalization
  type            (c_ptr                                      )                              :: parameterPointer
  type            (vector                                     )                              :: difference
  integer                                                                                    :: i                               , j
  double precision                                                                           :: binTimeMinimum                  , binTimeMaximum                    , &
       &                                                                                        timeMinimum                     , timeMaximum                       , &
       &                                                                                        distanceMaximum                 , time                              , &
       &                                                                                        massFunction                    , normalization
  logical                                                                                    :: integrationReset                , integrationResetNormalization

  ! Do not evaluate if the proposed prior is impossible.
  if (logPriorProposed <= logImpossible) then
     likelihoodMassFunctionEvaluate=0.0d0
     return
  end if
  ! Get the default cosmology functions object.
  cosmologyFunctions_ => cosmologyFunctions()     
  ! Get the default survey geometry object.
  surveyGeometry_     => surveyGeometry    ()     
  ! Construct the conditional mass function object.
  stateVector=simulationState%get()
  if     (                                                                  &
       &   (.not.self%modelSurfaceBrightness .and. size(stateVector) /= 11) &
       &  .or.                                                              &
       &   (     self%modelSurfaceBrightness .and. size(stateVector) /= 14) &
       & )                                                                  &
       & call Galacticus_Error_Report('likelihoodMassFunctionEvaluate','11 or 14 parameters are required for this likelihood function')
  do i=1,size(stateVector)
     stateVector(i)=parameterMappings(i)%thisMapping%unmap(stateVector(i))
  end do
  thisConditionalMassFunction                                  &
       & =conditionalMassFunctionBehroozi2010(                 &
       &                                      stateVector( 1), &
       &                                      stateVector( 2), &
       &                                      stateVector( 3), &
       &                                      stateVector( 4), &
       &                                      stateVector( 5), &
       &                                      stateVector( 6), &
       &                                      stateVector( 7), &
       &                                      stateVector( 8), &
       &                                      stateVector( 9), &
       &                                      stateVector(10), &
       &                                      stateVector(11)  &
       &                                     )
  ! Extract surface brightness parameters and get an incompleteness model.
  if (self%modelSurfaceBrightness)                                                  &
       & massFunctionIncompletenessModel                                            &
       &  =massFunctionIncompletenessSurfaceBrightness(                             &
       &                                               self%surfaceBrightnessLimit, &
       &                                               1.0d0                      , &
       &                                               stateVector(12)            , &
       &                                               stateVector(13)            , &
       &                                               stateVector(14)              &
       &                                               )
  ! Compute the time corresponding to the specified redshifts.
  timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftMaximum))
  timeMaximum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftMinimum))
  ! Compute the mass function.
  integrationReset             =.true.
  integrationResetNormalization=.true.
  do i=1,size(self%mass)
     ! A survey geometry is imposed - sum over fields.
     massFunction =0.0d0
     normalization=0.0d0
     do j=1,surveyGeometry_%fieldCount()
        ! Find the maximum distance at which a galaxy of the present mass can be detected in this survey.
        distanceMaximum=surveyGeometry_%distanceMaximum(sqrt(self%massMinimum(i)*self%massMaximum(i)),j)
        ! Set integration limits appropriately.
        binTimeMaximum=                                                       timeMaximum
        binTimeMinimum=min(                                                                 &
             &                                                                timeMaximum , &
             &             max(                                                             &
             &                                                                timeMinimum , &
             &                 cosmologyFunctions_%timeAtDistanceComoving(distanceMaximum)  &
             &                )                                                             &
             &            )
        ! Integrate the mass function over this time interval.
        massFunction=                                                               &
             &       +massFunction                                                  &
             &       +surveyGeometry_%solidAngle(j)                                 &
             &       *Integrate(                                                    &
             &                   binTimeMinimum                                   , &
             &                   binTimeMaximum                                   , &
             &                   likelihoodMassFunctionTimeIntegrand              , &
             &                   parameterPointer                                 , &
             &                   integrandFunction                                , &
             &                   integrationWorkspace                             , &
             &                   toleranceRelative=1.0d-3                         , &
             &                   reset=integrationReset                             &
             &                  )
        normalization=                                                              &
             &        +normalization                                                &
             &        +surveyGeometry_%solidAngle(j)                                &
             &        *Integrate(                                                   &
             &                   binTimeMinimum                                   , &
             &                   binTimeMaximum                                   , &
             &                   likelihoodMassFunctionTimeNormalizationIntegrand , &
             &                   parameterPointer                                 , &
             &                   integrandFunctionNormalization                   , &
             &                   integrationWorkspaceNormalization                , &
             &                   toleranceRelative=1.0d-3                         , &
             &                   reset=integrationResetNormalization                &
             &                  )
     end do
     self%massFunction(i)=                          &
          &               +massFunction             &
          &               /normalization            &
          &               /log(                     &
          &                     self%massMaximum(i) &
          &                    /self%massMinimum(i) &
          &                   )
     ! Account for surface brightness incompleteness if required.
     if (self%modelSurfaceBrightness) self%massFunction(i)=self%massFunction(i)*massFunctionIncompletenessModel%completeness(self%mass(i))
  end do
  call Integrate_Done(integrandFunction             ,integrationWorkspace             )
  call Integrate_Done(integrandFunctionNormalization,integrationWorkspaceNormalization)
  deallocate(stateVector)
  ! Evaluate the log-likelihood.
  difference                    =self%massFunction-self%massFunctionObserved
  likelihoodMassFunctionEvaluate=-0.5d0*(difference*(self%inverseCovariance*difference))
  return
  
contains
  
  function likelihoodMassFunctionTimeIntegrand(timePrime,parameterPointer) bind(c)
    !% Integral over time.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    use String_Handling
    implicit none
    real     (c_double                  )        :: likelihoodMassFunctionTimeIntegrand
    real     (c_double                  ), value :: timePrime
    type     (c_ptr                     ), value :: parameterPointer
    type     (fgsl_function             )        :: integrandFunctionTime
    type     (fgsl_integration_workspace)        :: integrationWorkspaceTime
    type     (c_ptr                     )        :: parameterPointerTime
    integer                                      :: errorStatus
    logical                                      :: integrationResetTime
    type     (varying_string            )        :: message
    character(len=14                    )        :: label

    time                =timePrime
    integrationResetTime=.true.
    likelihoodMassFunctionTimeIntegrand= Integrate(                                            &
         &                                           self%logHaloMassMinimum                 , &
         &                                           self%logHaloMassMaximum                 , &
         &                                           likelihoodMassFunctionHaloMassIntegrand , &
         &                                           parameterPointerTime                    , &
         &                                           integrandFunctionTime                   , &
         &                                           integrationWorkspaceTime                , &
         &                                           toleranceRelative=1.0d-3                , &
         &                                           toleranceAbsolute=1.0d-9                , &
         &                                           reset      =integrationResetTime        , &
         &                                           errorStatus=errorStatus                   &
         &                                          )                                          &
         &                               *cosmologyFunctions_%comovingVolumeElementTime(timePrime)
    call Integrate_Done(integrandFunctionTime,integrationWorkspaceTime)
    if (errorStatus /= errorStatusSuccess) then
       message='integration failed - state vector follows'
       do i=1,size(stateVector)
          write (label,'(e12.6)') stateVector(i)
          message=message//char(10)//" state ["//i//"] = "//trim(label)
        end do       
       call Galacticus_Error_Report('likelihoodMassFunctionTimeIntegrand',message)
    end if
    return
  end function likelihoodMassFunctionTimeIntegrand
  
  function likelihoodMassFunctionTimeNormalizationIntegrand(timePrime,parameterPointer) bind(c)
    !% Normalization integral over time.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)        :: likelihoodMassFunctionTimeNormalizationIntegrand
    real(c_double), value :: timePrime
    type(c_ptr),    value :: parameterPointer
    
    likelihoodMassFunctionTimeNormalizationIntegrand=cosmologyFunctions_%comovingVolumeElementTime(timePrime)
    return
  end function likelihoodMassFunctionTimeNormalizationIntegrand
  
  function likelihoodMassFunctionHaloMassIntegrand(logMass,parameterPointer) bind(c)
    !% Integral over halo mass function.
    use, intrinsic :: ISO_C_Binding
    use Halo_Mass_Function
    implicit none
    real(c_double)          :: likelihoodMassFunctionHaloMassIntegrand
    real(c_double)  , value :: logMass
    type(c_ptr   )  , value :: parameterPointer
    double precision        :: mass

    mass=10.0d0**logMass
    likelihoodMassFunctionHaloMassIntegrand= Halo_Mass_Function_Differential(time,mass)                                 &
         &                                    *                                   mass                                  &
         &                                    *log(10.0d0)                                                              &
         &                                    *max(                                                                     &
         &                                         +0.0d0                                                             , &
         &                                         +thisConditionalMassFunction%massFunction(mass,self%massMinimum(i))  &
         &                                         -thisConditionalMassFunction%massFunction(mass,self%massMaximum(i))  &
         &                                        )
    return
  end function likelihoodMassFunctionHaloMassIntegrand
    
end function likelihoodMassFunctionEvaluate

subroutine likelihoodMassFunctionFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodMassFunction), intent(inout) :: self
  
  return
end subroutine likelihoodMassFunctionFunctionChanged
