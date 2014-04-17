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

function likelihoodMassFunctionConstructor(haloMassMinimum,haloMassMaximum,redshiftMinimum,redshiftMaximum,useSurveyLimits,massFunctionFileName)
  !% Constructor for the mass function likelihood class.
  use Galacticus_Input_Paths
  use IO_HDF5
  use Memory_Management
  implicit none
  type            (likelihoodMassFunction)                :: likelihoodMassFunctionConstructor
  double precision                        , intent(in   ) :: redshiftMinimum                  , redshiftMaximum, &
       &                                                     haloMassMinimum                  , haloMassMaximum
  logical                                 , intent(in   ) :: useSurveyLimits
  character       (len=*                 ), intent(in   ) :: massFunctionFileName
  type            (hdf5Object            )                :: massFunctionFile
  integer                                                 :: i

  likelihoodMassFunctionConstructor%useSurveyLimits=useSurveyLimits
  likelihoodMassFunctionConstructor%   haloMassMinimum=      haloMassMinimum
  likelihoodMassFunctionConstructor%   haloMassMaximum=      haloMassMaximum
  likelihoodMassFunctionConstructor%logHaloMassMinimum=log10(haloMassMinimum)
  likelihoodMassFunctionConstructor%logHaloMassMaximum=log10(haloMassMaximum)
  likelihoodMassFunctionConstructor%   redshiftMinimum=      redshiftMinimum
  likelihoodMassFunctionConstructor%   redshiftMaximum=      redshiftMaximum
  ! Read the mass function file.
  !$omp critical(HDF5_Access)
  call massFunctionFile%openFile(char(Galacticus_Input_Path())//massFunctionFileName,readOnly=.true.)
  call massFunctionFile%readDataset("mass"                ,likelihoodMassFunctionConstructor%mass                )
  call massFunctionFile%readDataset("massFunctionObserved",likelihoodMassFunctionConstructor%massFunctionObserved)
  call massFunctionFile%readDataset("covariance"          ,likelihoodMassFunctionConstructor%covarianceMatrix    )
  call massFunctionFile%close()
  !$omp end critical(HDF5_Access)
  ! Find the inverse covariance matrix.
  likelihoodMassFunctionConstructor%covariance       =likelihoodMassFunctionConstructor%covarianceMatrix
  likelihoodMassFunctionConstructor%inverseCovariance=likelihoodMassFunctionConstructor%covariance      %invert()
  ! Compute mass bin limits.
  call Alloc_Array(likelihoodMassFunctionConstructor%massFunction,shape(likelihoodMassFunctionConstructor%mass))
  call Alloc_Array(likelihoodMassFunctionConstructor%massMinimum ,shape(likelihoodMassFunctionConstructor%mass))
  call Alloc_Array(likelihoodMassFunctionConstructor%massMaximum ,shape(likelihoodMassFunctionConstructor%mass))
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
  return
end function likelihoodMassFunctionConstructor

double precision function likelihoodMassFunctionEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return the log-likelihood for the mass function likelihood function.
  use, intrinsic :: ISO_C_Binding
  use Constraints_State
  use Constraints_Convergence
  use Conditional_Mass_Functions
  use Cosmology_Functions
  use Geometry_Surveys
  use Numerical_Integration
  use Galacticus_Error
use mpi_utilities
  implicit none
  class           (likelihoodMassFunction             ), intent(inout)               :: self
  class           (state                              ), intent(in   )               :: simulationState
  type            (mappingList                        ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence                        ), intent(inout)               :: simulationConvergence
  double precision                                     , intent(in   )               :: temperature                , logLikelihoodCurrent              , &
       &                                                                                logPriorCurrent            , logPriorProposed
  double precision                                     , allocatable  , dimension(:) :: stateVector
  class           (cosmologyFunctionsClass            ), pointer                     :: cosmologyFunctions_
  type            (conditionalMassFunctionBehroozi2010)                              :: thisConditionalMassFunction
  type            (fgsl_function                      )                              :: integrandFunction           , integrandFunctionNormalization
  type            (fgsl_integration_workspace         )                              :: integrationWorkspace        , integrationWorkspaceNormalization
  type            (c_ptr                              )                              :: parameterPointer
  type            (vector                             )                              :: difference
  integer                                                                            :: i
  double precision                                                                   :: binTimeMinimum             , binTimeMaximum                    , &
       &                                                                                timeMinimum                , timeMaximum                       , &
       &                                                                                distanceMaximum            , time
  logical                                                                            :: integrationReset           , integrationResetNormalization

  ! Do not evaluate if the proposed prior is impossible.
  if (logPriorProposed <= logImpossible) then
     likelihoodMassFunctionEvaluate=0.0d0
     return
  end if
  ! Get the default cosmology functions object.
  cosmologyFunctions_ => cosmologyFunctions()     
  ! Construct the conditional mass function object.
  stateVector=simulationState%get()
  if (size(stateVector) /= 11) call Galacticus_Error_Report('likelihoodMassFunctionEvaluate','11 parameters are required for this likelihood function')
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
  deallocate(stateVector)
  ! Compute the time corresponding to the specified redshifts.
  timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftMaximum))
  timeMaximum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftMinimum))
  ! Compute the mass function.
  integrationReset             =.true.
  integrationResetNormalization=.true.
  do i=1,size(self%mass)
     ! A survey geometry is imposed. Find the maximum distance at which a galaxy of the present
     ! mass can be detected in this survey.
     distanceMaximum=Geometry_Survey_Distance_Maximum(sqrt(self%massMinimum(i)*self%massMaximum(i)))
     ! Set integration limits appropriately.
     binTimeMinimum=max(timeMinimum,cosmologyFunctions_%timeAtDistanceComoving(distanceMaximum))
     binTimeMaximum=    timeMaximum
     ! Integrate the mass function over this time interval.
     self%massFunction(i)= Integrate(                                                    &
          &                          binTimeMinimum                                   , &
          &                          binTimeMaximum                                   , &
          &                          likelihoodMassFunctionTimeIntegrand              , &
          &                          parameterPointer                                 , &
          &                          integrandFunction                                , &
          &                          integrationWorkspace                             , &
          &                          toleranceRelative=1.0d-3                         , &
          &                          reset=integrationReset                             &
          &                         )                                                   &
          &               /Integrate(                                                   &
          &                          binTimeMinimum                                   , &
          &                          binTimeMaximum                                   , &
          &                          likelihoodMassFunctionTimeNormalizationIntegrand , &
          &                          parameterPointer                                 , &
          &                          integrandFunctionNormalization                   , &
          &                          integrationWorkspaceNormalization                , &
          &                          toleranceRelative=1.0d-3                         , &
          &                          reset=integrationResetNormalization                &
          &                         )                                                   &
          &               /log      (                                                   &
          &                           self%massMaximum(i)                               &
          &                          /self%massMinimum(i)                               &
          &                         )
 end do
 call Integrate_Done(integrandFunction             ,integrationWorkspace             )
 call Integrate_Done(integrandFunctionNormalization,integrationWorkspaceNormalization)
 ! Evaluate the log-likelihood.
 difference                    =self%massFunction-self%massFunctionObserved
 likelihoodMassFunctionEvaluate=-0.5d0*(difference*(self%inverseCovariance*difference))
 return
 
contains
  
  function likelihoodMassFunctionTimeIntegrand(timePrime,parameterPointer) bind(c)
    !% Integral over time.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double                  )        :: likelihoodMassFunctionTimeIntegrand
    real(c_double                  ), value :: timePrime
    type(c_ptr                     ), value :: parameterPointer
    type(fgsl_function             )        :: integrandFunctionTime
    type(fgsl_integration_workspace)        :: integrationWorkspaceTime
    type(c_ptr                     )        :: parameterPointerTime
    logical                                 :: integrationResetTime

    time=timePrime
    integrationResetTime=.true.
    likelihoodMassFunctionTimeIntegrand= Integrate(                                             &
         &                                           self%logHaloMassMinimum                 , &
         &                                           self%logHaloMassMaximum                 , &
         &                                           likelihoodMassFunctionHaloMassIntegrand , &
         &                                           parameterPointerTime                    , &
         &                                           integrandFunctionTime                   , &
         &                                           integrationWorkspaceTime                , &
         &                                           toleranceRelative=1.0d-3                , &
         &                                           reset=integrationResetTime                &
         &                                          )                                          &
         &                               *cosmologyFunctions_%comovingVolumeElementTime(timePrime)
    call Integrate_Done(integrandFunctionTime,integrationWorkspaceTime)
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
    likelihoodMassFunctionHaloMassIntegrand= Halo_Mass_Function_Differential(time,mass)                             &
         &                                    *                                   mass                              &
         &                                    *log(10.0d0)                                                          &
         &                                    *(                                                                    &
         &                                       thisConditionalMassFunction%massFunction(mass,self%massMinimum(i)) &
         &                                      -thisConditionalMassFunction%massFunction(mass,self%massMaximum(i)) &
         &                                     )
    return
  end function likelihoodMassFunctionHaloMassIntegrand
    
end function likelihoodMassFunctionEvaluate

subroutine likelihoodMassFunctionFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodMassFunction), intent(inout) :: self
  
  return
end subroutine likelihoodMassFunctionFunctionChanged
