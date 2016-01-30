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

function likelihoodHaloMassFunctionConstructor(fileName,redshift,massRangeMinimum)
  !% Constructor for the halo mass function likelihood class.
  use IO_HDF5
  use Memory_Management
  use Cosmology_Functions
  use Galacticus_Display
  implicit none
  type            (likelihoodHaloMassFunction)                                :: likelihoodHaloMassFunctionConstructor
  character       (len=*                     ), intent(in   )                 :: fileName
  double precision                            , intent(in   )                 :: redshift                             , massRangeMinimum
  class           (cosmologyFunctionsClass   ), pointer                       :: cosmologyFunctions_
  double precision                            , allocatable  , dimension(:  ) :: eigenValueArray                      , massOriginal     , &
       &                                                                         massFunctionOriginal
  double precision                            , allocatable  , dimension(:,:) :: massFunctionCovarianceOriginal
  character       (len=6                     )                                :: redshiftLabel
  type            (hdf5Object                )                                :: massFunctionFile                     , massFunctionGroup, &
       &                                                                         analysisGroup
  integer                                                                     :: i                                    , j                , &
       &                                                                         ii                                   , jj               , &
       &                                                                         massCountReduced
  double precision                                                            :: massIntervalLogarithmic
  type            (matrix                    )                                :: eigenVectors
  type            (vector                    )                                :: eigenValues

  ! Convert redshift to time.
  cosmologyFunctions_ => cosmologyFunctions()
  likelihoodHaloMassFunctionConstructor%time=cosmologyFunctions_ %cosmicTime                 (          &
       &                                      cosmologyFunctions_%expansionFactorFromRedshift (         &
       &                                                                                       redshift &
       &                                                                                      )         &
       &                                                                                     )
  ! Read the halo mass function file.
  write (redshiftLabel,'(f5.2)') redshift
  !$omp critical(HDF5_Access)
  call massFunctionFile %openFile(trim(fileName),readOnly=.true.)
  analysisGroup    =massFunctionFile%openGroup('analysis'                                       )
  massFunctionGroup=analysisGroup   %openGroup('haloMassFunctionZ'//trim(adjustl(redshiftLabel)))
  call massFunctionGroup%readDataset("mass"                  ,massOriginal                  )
  call massFunctionGroup%readDataset("massFunction"          ,massFunctionOriginal          )
  call massFunctionGroup%readDataset("massFunctionCovariance",massFunctionCovarianceOriginal)
  call massFunctionGroup%close()
  call analysisGroup    %close()
  call massFunctionFile %close()
  !$omp end critical(HDF5_Access)
  ! Find a reduced mass function excluding any empty bins.
  massCountReduced=count(massFunctionOriginal > 0.0d0 .and. massOriginal > massRangeMinimum)
  call Alloc_Array(likelihoodHaloMassFunctionConstructor%mass            ,[massCountReduced                 ])
  call Alloc_Array(likelihoodHaloMassFunctionConstructor%massFunction    ,[massCountReduced                 ])
  call Alloc_Array(likelihoodHaloMassFunctionConstructor%covarianceMatrix,[massCountReduced,massCountReduced])
  ii=0
  do i=1,size(massOriginal)
     if (massFunctionOriginal(i) > 0.0d0 .and. massOriginal(i) > massRangeMinimum) then
        ii=ii+1
        likelihoodHaloMassFunctionConstructor%mass        (ii)=massOriginal        (i)
        likelihoodHaloMassFunctionConstructor%massFunction(ii)=massFunctionOriginal(i)
        jj=0
        do j=1,size(massOriginal)
           if (massFunctionOriginal(j) > 0.0d0 .and. massOriginal(j) > massRangeMinimum) then
              jj=jj+1
              likelihoodHaloMassFunctionConstructor%covarianceMatrix(ii,jj)=massFunctionCovarianceOriginal(i,j)
           end if
        end do
     end if
  end do
  ! Compute mass ranges for bins.
  massIntervalLogarithmic=+log(                                  &
       &                       +massOriginal(size(massOriginal)) &
       &                       /massOriginal(                 1) &
       &                      )                                  &
       &                  /dble(                                 &
       &                        +size(massOriginal)              &
       &                        -1                               &
       &                       )
  call Alloc_Array(likelihoodHaloMassFunctionConstructor%massMinimum,shape(likelihoodHaloMassFunctionConstructor%mass))
  call Alloc_Array(likelihoodHaloMassFunctionConstructor%massMaximum,shape(likelihoodHaloMassFunctionConstructor%mass))
  do i=1,size(likelihoodHaloMassFunctionConstructor%mass)
     likelihoodHaloMassFunctionConstructor%massMinimum(i)=likelihoodHaloMassFunctionConstructor%mass(i)*exp(-0.5d0*massIntervalLogarithmic)
     likelihoodHaloMassFunctionConstructor%massMaximum(i)=likelihoodHaloMassFunctionConstructor%mass(i)*exp(+0.5d0*massIntervalLogarithmic)
  end do
  ! Find the inverse covariance matrix.
  likelihoodHaloMassFunctionConstructor%covariance       =likelihoodHaloMassFunctionConstructor%covarianceMatrix
  likelihoodHaloMassFunctionConstructor%inverseCovariance=likelihoodHaloMassFunctionConstructor%covariance      %invert()
  ! Symmetrize the inverse covariance matrix.
  call likelihoodHaloMassFunctionConstructor%inverseCovariance%symmetrize()
  ! Get eigenvalues and vectors of the inverse covariance matrix.
  allocate(eigenValueArray(size(likelihoodHaloMassFunctionConstructor%mass)))
  call likelihoodHaloMassFunctionConstructor%inverseCovariance%eigenSystem(eigenVectors,eigenValues)
  eigenValueArray=eigenValues
  if (any(eigenValueArray < 0.0d0)) call Galacticus_Display_Message('WARNING: inverse covariance matrix is not semi-positive definite')
  deallocate(eigenValueArray               )
  deallocate(massOriginal                  )
  deallocate(massFunctionOriginal          )
  deallocate(massFunctionCovarianceOriginal)
  return
end function likelihoodHaloMassFunctionConstructor

double precision function likelihoodHaloMassFunctionEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
  !% Return the log-likelihood for the halo mass function likelihood function.
  use Constraints_State
  use Constraints_Constants
  use Constraints_Convergence
  use Galacticus_Error
  use Halo_Mass_Functions
  use Cosmology_Parameters
  use Critical_Overdensities
  use Cosmological_Mass_Variance
  implicit none
  class           (likelihoodHaloMassFunction   ), intent(inout)               :: self
  class           (state                        ), intent(in   )               :: simulationState
  type            (mappingList                  ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence                  ), intent(inout)               :: simulationConvergence
  double precision                               , intent(in   )               :: temperature              , logLikelihoodCurrent, &
       &                                                                          logPriorCurrent          , logPriorProposed
  real                                           , intent(inout)               :: timeEvaluate
  double precision                               , intent(  out), optional     :: logLikelihoodVariance
  double precision                               , allocatable  , dimension(:) :: stateVector              , massFunction
  class           (cosmologyParametersClass     ), pointer                     :: cosmologyParameters_
  class           (cosmologicalMassVarianceClass), pointer                     :: cosmologicalMassVariance_
  class           (criticalOverdensityClass     ), pointer                     :: criticalOverdensity_
  type            (haloMassFunctionShethTormen  )                              :: haloMassFunction_
  type            (vector                       )                              :: difference
  integer                                                                      :: i

  ! There is no variance in our likelihood estimate.
  if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
  ! Do not evaluate if the proposed prior is impossible.
  if (logPriorProposed <= logImpossible) then
     likelihoodHaloMassFunctionEvaluate=0.0d0
     return
  end if
  ! Get required objects.
  cosmologyParameters_      => cosmologyParameters     ()
  cosmologicalMassVariance_ => cosmologicalMassVariance()
  criticalOverdensity_      => criticalOverdensity     ()
  ! Build the halo mass function object.
  stateVector=simulationState%get()
  if (size(stateVector) /= 3) call Galacticus_Error_Report('likelihoodHaloMassFunctionEvaluate','3 parameters are required for this likelihood function')
  do i=1,size(stateVector)
     stateVector(i)=parameterMappings(i)%thisMapping%unmap(stateVector(i))
  end do
  haloMassFunction_=haloMassFunctionShethTormen(                              &
       &                                        cosmologyParameters_        , &
       &                                        cosmologicalMassVariance_   , &
       &                                        criticalOverdensity_        , &
       &                                        stateVector              (1), &
       &                                        stateVector              (2), &
       &                                        stateVector              (3)  &
       &                                       )
  ! Compute the mass function.
  allocate(massFunction(size(self%mass)))
  do i=1,size(self%mass)
     massFunction(i)=+haloMassFunction_%integrated(                     &
          &                                        self%time          , &
          &                                        self%massMinimum(i), &
          &                                        self%massMaximum(i)  &
          &                                       )                     &
          &          /log(                                              &
          &                                       +self%massMaximum(i)  &
          &                                       /self%massMinimum(i)  &
          &              )     
  end do
  ! Evaluate the log-likelihood.
  difference                        =massFunction-self%massFunction
  likelihoodHaloMassFunctionEvaluate=-0.5d0*(difference*(self%inverseCovariance*difference))
  ! Clean up.
  deallocate(stateVector )
  deallocate(massFunction)
  return    
end function likelihoodHaloMassFunctionEvaluate

subroutine likelihoodHaloMassFunctionFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodHaloMassFunction), intent(inout) :: self
  
  return
end subroutine likelihoodHaloMassFunctionFunctionChanged
