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

!% Defines methods for a projected correlation function likelihood function.

function likelihoodProjectedCorrelationFunctionConstructor(haloMassMinimum,haloMassMaximum,lineOfSightDepth,projectedCorrelationFunctionFileName)
  !% Constructor for the projected correlation function likelihood class.
  use Galacticus_Input_Paths
  use IO_HDF5
  use Memory_Management
  use Node_Component_Dark_Matter_Profile_Scale
  use Galacticus_Nodes
  implicit none
  type            (likelihoodProjectedCorrelationFunction)                              :: likelihoodProjectedCorrelationFunctionConstructor
  double precision                                        , intent(in   )               :: haloMassMinimum                                  , haloMassMaximum, &
       &                                                                                   lineOfSightDepth
  character       (len=*                 )                , intent(in   )               :: projectedCorrelationFunctionFileName
  double precision                                        , allocatable  , dimension(:) :: separationBinWidth
  type            (hdf5Object            )                                              :: projectedCorrelationFunctionFile
  integer                                                                               :: i

  likelihoodProjectedCorrelationFunctionConstructor% haloMassMinimum= haloMassMinimum
  likelihoodProjectedCorrelationFunctionConstructor% haloMassMaximum= haloMassMaximum
  likelihoodProjectedCorrelationFunctionConstructor%lineOfSightDepth=lineOfSightDepth
  ! Read the projected correlation function file.
  !$omp critical(HDF5_Access)
  call projectedCorrelationFunctionFile%openFile(char(Galacticus_Input_Path())//projectedCorrelationFunctionFileName,readOnly=.true.)
  call projectedCorrelationFunctionFile%readDataset("separation"                          ,likelihoodProjectedCorrelationFunctionConstructor%separation                          )
  call projectedCorrelationFunctionFile%readDataset("projectedCorrelationFunctionObserved",likelihoodProjectedCorrelationFunctionConstructor%projectedCorrelationFunctionObserved)
  call projectedCorrelationFunctionFile%readDataset("covariance"                          ,likelihoodProjectedCorrelationFunctionConstructor%covarianceMatrix                    )
  call projectedCorrelationFunctionFile%readDataset("massMinimum"                         ,likelihoodProjectedCorrelationFunctionConstructor%massMinimum                         )
  call projectedCorrelationFunctionFile%readDataset("massMaximum"                         ,likelihoodProjectedCorrelationFunctionConstructor%massMaximum                         )
  call projectedCorrelationFunctionFile%close()
  !$omp end critical(HDF5_Access)
  ! Allocate storage for the model projected correlation function.
  call Alloc_Array(likelihoodProjectedCorrelationFunctionConstructor%projectedCorrelationFunction,[size(likelihoodProjectedCorrelationFunctionConstructor%separation),size(likelihoodProjectedCorrelationFunctionConstructor%massMinimum)])
  ! Find the inverse covariance matrix.
  likelihoodProjectedCorrelationFunctionConstructor%covariance       =likelihoodProjectedCorrelationFunctionConstructor%covarianceMatrix
  likelihoodProjectedCorrelationFunctionConstructor%inverseCovariance=likelihoodProjectedCorrelationFunctionConstructor%covariance      %invert()
  ! Ensure the nodes objects are initialized.
  call Galacticus_Nodes_Initialize                        ()
  call Node_Component_Dark_Matter_Profile_Scale_Initialize()
  return
end function likelihoodProjectedCorrelationFunctionConstructor

double precision function likelihoodProjectedCorrelationFunctionEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return the log-likelihood for the projected correlation function likelihood function.
  use Constraints_State
  use Constraints_Convergence
  use Conditional_Mass_Functions
  use Cosmology_Functions
  use Galacticus_Error
  use Halo_Model_Projected_Correlations
  implicit none
  class           (likelihoodProjectedCorrelationFunction), intent(inout)               :: self
  class           (state                                 ), intent(in   )               :: simulationState
  type            (mappingList                           ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence                           ), intent(inout)               :: simulationConvergence
  double precision                                        , intent(in   )               :: temperature                , logLikelihoodCurrent, &
       &                                                                                   logPriorCurrent            , logPriorProposed
  double precision                                        , allocatable  , dimension(:) :: stateVector
  class           (cosmologyFunctionsClass               ), pointer                     :: cosmologyFunctions_
  type            (conditionalMassFunctionBehroozi2010   )                              :: thisConditionalMassFunction
  type            (vector                                )                              :: difference
  integer                                                                               :: i

  ! Do not evaluate if the proposed prior is impossible.
  if (logPriorProposed <= logImpossible) then
     likelihoodProjectedCorrelationFunctionEvaluate=0.0d0
     return
  end if
  ! Get the default cosmology functions object.
  cosmologyFunctions_ => cosmologyFunctions()     
  ! Construct the conditional mass function object.
  stateVector=simulationState%get()
  if (size(stateVector) /= 11) call Galacticus_Error_Report('likelihoodProjectedCorrelationFunctionEvaluate','11 parameters are required for this likelihood function')
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
  ! Compute the projected correlation function.
  do i=1,size(self%massMinimum)
     call Halo_Model_Projected_Correlation(                                        &
          &                                thisConditionalMassFunction           , &
          &                                self%separation                       , &
          &                                self%massMinimum                 (  i), &
          &                                self%massMaximum                 (  i), &
          &                                self%haloMassMinimum                  , &
          &                                self%haloMassMaximum                  , &
          &                                self%lineOfSightDepth                 , &
          &                                self%projectedCorrelationFunction(:,i)  &
          &                               )
  end do
  ! Evaluate the log-likelihood.
  difference                                    =reshape(                                                         &
       &                                                  +     self%projectedCorrelationFunction                 &
       &                                                  -     self%projectedCorrelationFunctionObserved       , &
       &                                                 [                                                        &
       &                                                   size(self%projectedCorrelationFunction        ,dim=1)  &
       &                                                  *size(self%projectedCorrelationFunction        ,dim=2)  &
       &                                                 ]                                                        &
       &                                                )

  likelihoodProjectedCorrelationFunctionEvaluate=-0.5d0*(difference*(self%inverseCovariance*difference))

write (0,*) likelihoodProjectedCorrelationFunctionEvaluate

  return
end function likelihoodProjectedCorrelationFunctionEvaluate

subroutine likelihoodProjectedCorrelationFunctionFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodProjectedCorrelationFunction), intent(inout) :: self
  
  return
end subroutine likelihoodProjectedCorrelationFunctionFunctionChanged
