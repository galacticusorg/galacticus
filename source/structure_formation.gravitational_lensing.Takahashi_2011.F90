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

!% Implements the gravitational lensing distributions of \cite{takahashi_probability_2011}.

  !# <gravitationalLensing name="gravitationalLensingTakahashi2011">
  !#  <description>Implements the gravitational lensing distributions of \cite{takahashi_probability_2011}.</description>
  !# </gravitationalLensing>

  type, extends(gravitationalLensingClass) :: gravitationalLensingTakahashi2011
     logical                                    :: tableInitialized         , cdfInitialized
     type            (table1DGeneric          ) :: convergencePDF
     type            (table1DLogarithmicLinear) :: magnificationCDFTable
     double precision                           :: redshiftPrevious         , convergenceEmptyBeam                , &
          &                                        convergenceVariance      , convergenceScale                    , &
          &                                        convergenceVarianceScaled, convergenceDistributionNormalization, &
          &                                        aConvergence             , omegaConvergence                    , &
          &                                        scaleSourcePrevious
   contains
     !@ <objectMethods>
     !@   <object>gravitationalLensingTakahashi2011</object>
     !@   <objectMethod>
     !@     <method>lensingDistributionConstruct</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ redshift\argin, \doublezero\ scaleSource\argin</arguments>
     !@     <description>Construct the gravitational lensing distribution functions for the specified redshift and source scale.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>convergenceDistribution</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ convergence\argin</arguments>
     !@     <description>Returns the gravitatoinal lensing convergence probability density function at the given convergence.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: magnificationPDF             => takahashi2011MagnificationPDF
     procedure :: magnificationCDF             => takahashi2011MagnificationCDF
     procedure :: convergenceDistribution      => takahashi2011ConvergenceDistribution
     procedure :: lensingDistributionConstruct => takahashi2011LensingDistributionConstruct
   end type gravitationalLensingTakahashi2011

  interface gravitationalLensingTakahashi2011
     !% Constructors for the \cite{takahashi_probability_2011} gravitational lensing class.
     module procedure takahashi2011DefaultConstructor
  end interface gravitationalLensingTakahashi2011

  ! Smallest variance for which calculations are stable.
  double precision, parameter :: takahashi2011ConvergenceVarianceSmall=1.0d-5

  ! Smallest redshift for which to compute lensing.
  double precision, parameter :: takahashi2011RedshiftTiny            =1.0d-2

contains

  function takahashi2011DefaultConstructor()
    !% Default constructor for the \cite{takahashi_probability_2011} gravitational lensing class.
    implicit none
    type(gravitationalLensingTakahashi2011) :: takahashi2011DefaultConstructor

    takahashi2011DefaultConstructor%tableInitialized=.false.
    takahashi2011DefaultConstructor%cdfInitialized  =.false.
    takahashi2011DefaultConstructor%redshiftPrevious=-2.0d0
    return
  end function takahashi2011DefaultConstructor

  double precision function takahashi2011MagnificationPDF(self,magnification,redshift,scaleSource)
    !% Compute the magnification probability density function at the given {\tt magnification} and {\tt redshift} using the
    !% \cite{takahashi_probability_2011} formalism.
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout) :: self
    double precision                                   , intent(in   ) :: magnification, redshift, &
         &                                                                scaleSource

    ! Handle redshift zero case.
    if (redshift <= takahashi2011RedshiftTiny) then
       if (magnification == 1.0d0) then
          takahashi2011MagnificationPDF=1.0d0
       else
          takahashi2011MagnificationPDF=0.0d0
       end if
       return
    else
       ! Construct the distribution.
       call self%lensingDistributionConstruct(redshift,scaleSource)
       ! Approximate a delta function for small redshifts.
       if (self%convergenceVariance < takahashi2011ConvergenceVarianceSmall) then
          if (magnification == 1.0d0) then
             takahashi2011MagnificationPDF=0.0d0
          else
             takahashi2011MagnificationPDF=1.0d0
          end if
          return
       end if
       ! Evaluate the magnification PDF (eqn. 11 of Takahashi et al.).
       takahashi2011MagnificationPDF=takahashi2011MagnificationDistribution(self,magnification)
    end if
    return
  end function takahashi2011MagnificationPDF
  
  double precision function takahashi2011MagnificationCDF(self,magnification,redshift,scaleSource)
    !% Compute the magnification probability density function at the given {\tt magnification} and {\tt redshift} using the
    !% \cite{takahashi_probability_2011} formalism.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout) :: self
    double precision                                   , intent(in   ) :: magnification, redshift       , &
         &                                                                scaleSource
    double precision                                   , parameter     :: magnificationMinimum=   1.0d-2
    double precision                                   , parameter     :: magnificationMaximum=1000.0d0
    integer                                            , parameter     :: magnificationCount   =1000
    integer                                                            :: i
    type            (fgsl_function                    )                :: integrandFunction
    type            (fgsl_integration_workspace       )                :: integrationWorkspace
    type            (c_ptr                            )                :: parameterPointer
    logical                                                            :: integrationReset
    double precision                                                   :: magnificationLower             , magnificationUpper, &
         &                                                                cdf                            , cdfPrevious

    ! Handle redshift zero case.
    if (redshift <= takahashi2011RedshiftTiny) then
       if (magnification < 1.0d0) then
          takahashi2011MagnificationCDF=0.0d0
       else
          takahashi2011MagnificationCDF=1.0d0
       end if
       return
    else
       ! Construct the distribution.
       call self%lensingDistributionConstruct(redshift,scaleSource)
       ! Approximate a delta function for small redshifts.
       if (self%convergenceVariance < takahashi2011ConvergenceVarianceSmall) then
          if (magnification < 1.0d0) then
             takahashi2011MagnificationCDF=0.0d0
          else
             takahashi2011MagnificationCDF=1.0d0
          end if
          return
       end if
       ! Tabulate the cumulative distribution function if table does not yet exist.
       if (.not.self%cdfInitialized) then
          call self%magnificationCDFTable%create(magnificationMinimum,magnificationMaximum,magnificationCount,tableCount=1)
          do i=1,magnificationCount
             integrationReset=.true.
             if (i == 1 ) then
                magnificationLower=0.0d0
                cdfPrevious       =0.0d0
             else
                magnificationLower=self%magnificationCDFTable%x(i-1)
                cdfPrevious       =self%magnificationCDFTable%y(i-1)
             end if
             magnificationUpper   =self%magnificationCDFTable%x(i  )
             cdf=Integrate(                           &
                  &        magnificationLower       , &
                  &        magnificationUpper       , &
                  &        magnificationPDFIntegrand, &
                  &        parameterPointer         , &
                  &        integrandFunction        , &
                  &        integrationWorkspace     , &
                  &        toleranceRelative=1.0d-6 , &
                  &        reset=integrationReset     &
                  &       )
             call self%magnificationCDFTable%populate(cdf+cdfPrevious,i)
             call Integrate_Done(integrandFunction,integrationWorkspace)
          end do
          self%cdfInitialized=.true.
       end if
       ! Interpolate in the tabulated cumulative distribution function.
       if      (magnification < magnificationMinimum) then
          takahashi2011MagnificationCDF=0.0d0
       else if (magnification > magnificationMaximum) then
          takahashi2011MagnificationCDF=1.0d0
       else
          takahashi2011MagnificationCDF=self%magnificationCDFTable%interpolate(magnification)/self%magnificationCDFTable%y(-1)
       end if
    end if
    return

  contains

    function magnificationPDFIntegrand(magnification,parameterPointer) bind(c)
      !% Integral for the magnification probability distribution function.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real            (c_double)        :: magnificationPDFIntegrand
      real            (c_double), value :: magnification
      type            (c_ptr   ), value :: parameterPointer
      
      magnificationPDFIntegrand=takahashi2011MagnificationDistribution(self,magnification)
     return
    end function magnificationPDFIntegrand
    
  end function takahashi2011MagnificationCDF
  
  double precision function takahashi2011MagnificationDistribution(self,magnification)
    !% The gravitational lensing magnification distribution from \cite[][eq.~11]{takahashi_probability_2011}.
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout) :: self
    double precision                                   , intent(in   ) :: magnification
    double precision                                   , parameter     :: magnificationZeroPoint=3.0d0
    double precision                                   , parameter     :: convergenceZeroPoint  =1.0d0-1.0d0/sqrt(magnificationZeroPoint)
    double precision                                                   :: convergence
    
    if (magnification <= 0.0d0) then
       takahashi2011MagnificationDistribution=0.0d0
    else
       convergence                           =          &
            & +1.0d0                                    &
            & -1.0d0                                    &
            & /sqrt(magnification)
       takahashi2011MagnificationDistribution=          &
            & +0.5d0                                    &
            & *(                                        &
            &   +1.0d0                                  &
            &   -convergence                            &
            & )**3                                      &
            & *self%convergenceDistribution(convergence)
       if (magnification > 1.0d0)                                  &
            & takahashi2011MagnificationDistribution=              &
            &  +takahashi2011MagnificationDistribution             &
            &  +exp(                                               &
            &       -0.25d0                                        &
            &       /(                                             &
            &         +magnification                               &
            &         -1.0d0                                       &
            &        )**4                                          &
            &      )                                               &
            &  *0.5d0                                              &
            &  *(                                                  &
            &    +1.0d0                                            &
            &    -convergenceZeroPoint                             &
            &   )**3                                               &
            &  *self%convergenceDistribution(convergenceZeroPoint) &
            &  /(                                                  &
            &     magnification                                    &
            &    /magnificationZeroPoint                           &
            &   )**3
    end if
    return
  end function takahashi2011MagnificationDistribution

  double precision function takahashi2011ConvergenceDistribution(self,convergence)
    !% The distribution function for gravitational lensing convergence \citep[][eqn.~8]{takahashi_probability_2011}.
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout) :: self
    double precision                                   , intent(in   ) :: convergence
    double precision                                                   :: scaledConvergence
    
    scaledConvergence=convergence/self%convergenceScale
    if (scaledConvergence <= -1.0d0) then
       takahashi2011ConvergenceDistribution=+0.0d0
    else
       takahashi2011ConvergenceDistribution=                                           &
            &                               +self%convergenceDistributionNormalization &
            &                               *exp(                                      &
            &                                    -0.5d0                                &
            &                                    /self%omegaConvergence**2             &
            &                                    *(                                    &
            &                                      +log(                               &
            &                                           +1.0d0                         &
            &                                           +scaledConvergence             &
            &                                          )                               &
            &                                      +self%omegaConvergence**2           &
            &                                      /2.0d0                              &
            &                                     )**2                                 &
            &                                    *(                                    &
            &                                      +1.0d0                              &
            &                                      +self%aConvergence                  &
            &                                      /(                                  &
            &                                        +1.0d0                            &
            &                                        +scaledConvergence                &
            &                                       )                                  &
            &                                     )                                    &
            &                                   )                                      &
            &                               /(                                         &
            &                                 +1.0d0                                   &
            &                                 +scaledConvergence                       &
            &                                )                                         &
            &                               /self%convergenceScale
    end if
    return
  end function takahashi2011ConvergenceDistribution
  
  subroutine takahashi2011LensingDistributionConstruct(self,redshift,scaleSource)
    !% Construct the lensing distribution function for the \cite{takahashi_probability_2011} formalism.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    use Numerical_Ranges
    use Numerical_Comparison
    use Cosmology_Parameters
    use Cosmology_Functions
    use Root_Finder
    use Galacticus_Error
    use Galacticus_Input_Paths
    use IO_HDF5
    use File_Utilities
    implicit none
    class           (gravitationalLensingTakahashi2011), intent(inout)               :: self
    double precision                                   , intent(in   )               :: redshift                                 , scaleSource
    class           (cosmologyParametersClass         ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer                     :: cosmologyFunctions_
    double precision                                   , parameter                   :: redshiftZero                          =   0.0d0
    double precision                                   , parameter                   :: convergenceParametersToleranceAbsolute=   0.0d0
    double precision                                   , parameter                   :: convergenceParametersToleranceRelative=   1.0d-6
    double precision                                   , parameter                   :: convergenceMaximumFactor              =  10.0d0
    double precision                                   , parameter                   :: omegaKappaMinimum                     =   1.0d-2
    double precision                                   , parameter                   :: omegaKappaMaximum                     =   1.0d0
    double precision                                   , parameter                   :: magnificationMinimum                  =   0.0d0
    double precision                                   , parameter                   :: magnificationMaximum                  =1000.0d0
    integer                                            , parameter                   :: omegaKappaCount                       =1000
    double precision                                   , allocatable  , dimension(:) :: tableOmegaKappa                           , tableAKappa, &
         &                                                                              tableNKappa                               , tableConvergenceVariance
    integer                                                                          :: i
    type            (rootFinder                       )                              :: finder
    logical                                                                          :: integrationReset
    type            (fgsl_function                    )                              :: integrandFunction
    type            (fgsl_integration_workspace       )                              :: integrationWorkspace
    type            (c_ptr                            )                              :: parameterPointer
    double precision                                                                 :: distanceComovingSource                    , timeLens             , &
         &                                                                              convergencePdfMoment0                     , convergencePdfMoment1, &
         &                                                                              convergencePdfMoment2                     , convergenceMinimum   , &
         &                                                                              convergenceMaximum                        , convergence          , &
         &                                                                              magnificationPdfMoment0
    type            (hdf5Object                       )                              :: parametersFile

    ! Construct tabulation of convergence distribution function parameters if necessary.
    if (.not.self%tableInitialized) then
       ! Determine the paramters, A_kappa and omega_kappa, of the convergence distribution (eq. 8
       ! of Takahashi et al.). To do this, we use a look-up table of precomputed values.
       ! Check if a precomputed file exists.
       !$omp critical(takahashi2011Tabulate)
       if (File_Exists(Galacticus_Input_Path()//"data/largeScaleStructure/gravitationalLensingConvergenceTakahashi2011.hdf5")) then
          ! Read the results from file.
          !$omp critical(HDF5_Access)
          call parametersFile%openFile(char(Galacticus_Input_Path()//"data/largeScaleStructure/gravitationalLensingConvergenceTakahashi2011.hdf5"),readOnly=.true.)
          call parametersFile%readDataset("convergenceVariance",tableConvergenceVariance)
          call parametersFile%readDataset(             "NKappa",tableNKappa             )
          call parametersFile%readDataset(             "AKappa",tableAKappa             )
          call parametersFile%readDataset(         "omegaKappa",tableOmegaKappa         )
          call parametersFile%close      (                                              )
          !$omp end critical(HDF5_Access)
       else
          ! Create a root-finder to solve the parameters.
          call finder%rootFunction(convergencePdfParameterSolver)
          call finder%tolerance   (                                        &
               &                   convergenceParametersToleranceAbsolute, &
               &                   convergenceParametersToleranceRelative  &
               &                  )
          call finder%rangeExpand (                                                             &
               &                   rangeExpandUpward            =2.0d0                        , &
               &                   rangeExpandDownward          =0.5d0                        , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          ! Construct a range of omega_kappa values.
          tableOmegaKappa=Make_Range(                      &
               &                     omegaKappaMinimum   , &
               &                     omegaKappaMaximum   , &
               &                     omegaKappaCount     , &
               &                     rangeTypeLogarithmic  &
               &                    )
          allocate(tableAKappa             (omegaKappaCount))
          allocate(tableNKappa             (omegaKappaCount))
          allocate(tableConvergenceVariance(omegaKappaCount))
          ! Iterate over values of omega_k
          do i=1,omegaKappaCount
             ! Set omega_kappa parameter.
             self%omegaConvergence                    =tableOmegaKappa(i)
             ! Set convergence scale to unity for table building. This implies a minimum convegence of -1.
             self%convergenceScale                    =+1.0d0
             convergenceMinimum                       =-self%convergenceScale
             convergenceMaximum                       =+convergenceMaximumFactor*sqrt(self%omegaConvergence)
             ! Set normalization to unity.
             self%convergenceDistributionNormalization=1.0d0
             ! Solve for a value of A_k that gives consistent first and second moments of the
             ! convergence distribution (eq. 9 of Takahashi et al.).
             self%aConvergence                        =finder%find(rootGuess=1.0d0)
             ! Evaluate zeroth and second moments of the convergence distribution.
             integrationReset     =.true.
             convergencePdfMoment0=Integrate(                                         &
                  &                             convergenceMinimum                     , &
                  &                             convergenceMaximum                     , &
                  &                             convergenceDistributionMoment0Integrand, &
                  &                             parameterPointer                       , &
                  &                             integrandFunction                      , &
                  &                             integrationWorkspace                   , &
                  &                             toleranceRelative=1.0d-4               , &
                  &                             reset=integrationReset                   &
                  &                            )
             call Integrate_Done(integrandFunction,integrationWorkspace)
             integrationReset     =.true.
             convergencePdfMoment1=Integrate(                                         &
                  &                             convergenceMinimum                     , &
                  &                             convergenceMaximum                     , &
                  &                             convergenceDistributionMoment1Integrand, &
                  &                             parameterPointer                       , &
                  &                             integrandFunction                      , &
                  &                             integrationWorkspace                   , &
                  &                             toleranceRelative=1.0d-4               , &
                  &                             reset=integrationReset                   &
                  &                            )
             call Integrate_Done(integrandFunction,integrationWorkspace)
             integrationReset     =.true.
             convergencePdfMoment2=Integrate(                                        &
                  &                         convergenceMinimum                     , &
                  &                         convergenceMaximum                     , &
                  &                         convergenceDistributionMoment2Integrand, &
                  &                         parameterPointer                       , &
                  &                         integrandFunction                      , &
                  &                         integrationWorkspace                   , &
                  &                         toleranceRelative=1.0d-4               , &
                  &                         reset=integrationReset                   &
                  &                        )
             call Integrate_Done(integrandFunction,integrationWorkspace)
             ! Store the A parameter.
             tableAKappa             (i)=self%aConvergence
             ! Compute the normalization of the convergence distribution.
             tableNKappa             (i)=1.0d0                /convergencePdfMoment0
             ! Determine what convergence variance this parameter combination corresponds to.
             tableConvergenceVariance(i)=convergencePdfMoment2/convergencePdfMoment0
             ! Check that the ratio of first to second moments equals -2.
             if (Values_Differ(convergencePdfMoment1/convergencePdfMoment2,-2.0d0,absTol=2.0d-3)) &
                  & call Galacticus_Error_Report('takahashi2011MagnificationPDF','convergence PDF does not satisfy consistency criterion')
          end do
          ! Store the results to file.
          !$omp critical(HDF5_Access)
          call parametersFile%openFile(char(Galacticus_Input_Path()//"data/largeScaleStructure/gravitationalLensingConvergenceTakahashi2011.hdf5"))
          call parametersFile%writeDataset(tableConvergenceVariance,"convergenceVariance","Dimensionless variance of lensing convergence"     )
          call parametersFile%writeDataset(tableNKappa             ,"NKappa"             ,"Parameter N_kappa from Takahashi et al. (2011)"    )
          call parametersFile%writeDataset(tableAKappa             ,"AKappa"             ,"Parameter A_kappa from Takahashi et al. (2011)"    )
          call parametersFile%writeDataset(tableOmegaKappa         ,"omegaKappa"         ,"Parameter omega_kappa from Takahashi et al. (2011)")
          call parametersFile%close()
          !$omp end critical(HDF5_Access)
       end if
       !$omp end critical(takahashi2011Tabulate)
       ! Create a table.
       call self%convergencePDF%create  (tableConvergenceVariance,tableCount=3)
       call self%convergencePDF%populate(tableNKappa             ,           1)
       call self%convergencePDF%populate(tableAKappa             ,           2)
       call self%convergencePDF%populate(tableOmegaKappa         ,           3)
       deallocate(tableConvergenceVariance)
       deallocate(tableNKappa             )
       deallocate(tableAKappa             )
       deallocate(tableOmegaKappa         )
       ! Record that the table is now initialized.
       self%tableInitialized=.true.
    end if
    ! Check if redshift has changed since previous call.
    if (redshift /= self%redshiftPrevious .or. scaleSource /= self%scaleSourcePrevious) then
       ! Update convergences.
       self%redshiftPrevious   =redshift
       self%scaleSourcePrevious=scaleSource
       if (self%cdfInitialized) then
          self%cdfInitialized  =.false.
          call self%magnificationCDFTable%destroy()
       end if
       ! Get cosmology functions.
       cosmologyFunctions_  => cosmologyFunctions ()
       cosmologyParameters_ => cosmologyParameters()
       ! Find the comoving distance to the source.
       distanceComovingSource=cosmologyFunctions_  %distanceComoving           (           &
            &                  cosmologyFunctions_ %cosmicTime                  (          &
            &                   cosmologyFunctions_%expansionFactorFromRedshift  (         &
            &                                                                     redshift &
            &                                                                    )         &
            &                                                                   )          &
            &                                                                  )
       ! Find the convergence of an empty beam.
       integrationReset         =.true.
       self%convergenceEmptyBeam=Integrate(                               &
            &                              redshiftZero                 , &
            &                              redshift                     , &
            &                              emptyBeamConvergenceIntegrand, &
            &                              parameterPointer             , &
            &                              integrandFunction            , &
            &                              integrationWorkspace         , &
            &                              toleranceRelative=1.0d-3     , &
            &                              reset=integrationReset         &
            &                            )
       call Integrate_Done(integrandFunction,integrationWorkspace)
       ! Find the variance of the convergence.
       integrationReset         =.true.
       self%convergenceVariance =Integrate(                               &
            &                              redshiftZero                 , &
            &                              redshift                     , &
            &                              convergenceVarianceIntegrand , &
            &                              parameterPointer             , &
            &                              integrandFunction            , &
            &                              integrationWorkspace         , &
            &                              toleranceRelative=1.0d-3     , &
            &                              reset=integrationReset         &
            &                             )
       call Integrate_Done(integrandFunction,integrationWorkspace)
       ! Determine the parameters of the convergence distribution.
       self%convergenceScale                    =abs(self%convergenceEmptyBeam)
       self%convergenceVarianceScaled           =self%convergenceVariance/self%convergenceScale**2
       self%convergenceDistributionNormalization=self%convergencePDF%interpolate(self%convergenceVarianceScaled,1)
       self%aConvergence                        =self%convergencePDF%interpolate(self%convergenceVarianceScaled,2)
       self%omegaConvergence                    =self%convergencePDF%interpolate(self%convergenceVarianceScaled,3)
       ! Integrate the modified magnification distribution in order to find the normalization.
       integrationReset        =.true.
       magnificationPdfMoment0=Integrate(                             &
            &                              magnificationMinimum     , &
            &                              magnificationMaximum     , &
            &                              magnificationPDFIntegrand, &
            &                              parameterPointer         , &
            &                              integrandFunction        , &
            &                              integrationWorkspace     , &
            &                              toleranceRelative=1.0d-3 , &
            &                              reset=integrationReset     &
            &                             )
       call Integrate_Done(integrandFunction,integrationWorkspace)
       self%convergenceDistributionNormalization         &
            & =self%convergenceDistributionNormalization &
            & /magnificationPdfMoment0
    end if
    return
    
  contains

    function magnificationPDFIntegrand(magnification,parameterPointer) bind(c)
      !% Integral for the magnification probability distribution function.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real            (c_double)        :: magnificationPDFIntegrand
      real            (c_double), value :: magnification
      type            (c_ptr   ), value :: parameterPointer

      magnificationPDFIntegrand=takahashi2011MagnificationDistribution(self,magnification)
      return
    end function magnificationPDFIntegrand
    
    double precision function convergencePdfParameterSolver(a)
      !% Root function used in finding equivalent circular orbits.
      implicit none
      double precision, intent(in   ) :: a
      
      self%aConvergence    =a
      integrationReset     =.true.
      convergencePdfMoment1=Integrate(                                        &
           &                         convergenceMinimum                     , &
           &                         convergenceMaximum                     , &
           &                         convergenceDistributionMoment1Integrand, &
           &                         parameterPointer                       , &
           &                         integrandFunction                      , &
           &                         integrationWorkspace                   , &
           &                         toleranceRelative=1.0d-3               , &
           &                         reset=integrationReset                   &
           &                        )
      call Integrate_Done(integrandFunction,integrationWorkspace)
      integrationReset     =.true.
      convergencePdfMoment2=Integrate(                                        &
           &                         convergenceMinimum                     , &
           &                         convergenceMaximum                     , &
           &                         convergenceDistributionMoment2Integrand, &
           &                         parameterPointer                       , &
           &                         integrandFunction                      , &
           &                         integrationWorkspace                   , &
           &                         toleranceRelative=1.0d-3               , &
           &                         reset=integrationReset                   &
           &                        )
      call Integrate_Done(integrandFunction,integrationWorkspace)
      if (convergencePdfMoment2 <= 0.0d0) then
         convergencePdfParameterSolver=0.0d0
      else
         convergencePdfParameterSolver=2.0d0+convergencePdfMoment1/convergencePdfMoment2*self%convergenceScale
      end if
    return
    end function convergencePdfParameterSolver
    
    function emptyBeamConvergenceIntegrand(redshiftLens,parameterPointer) bind(c)
      !% Integral for gravitational lensing convergence in an empty beam.
      use, intrinsic :: ISO_C_Binding
      use Numerical_Constants_Physical
      use Numerical_Constants_Prefixes
      implicit none
      real            (c_double)        :: emptyBeamConvergenceIntegrand
      real            (c_double), value :: redshiftLens
      type            (c_ptr   ), value :: parameterPointer
      double precision                  :: time                         , distanceComovingLens

      ! Find cosmic time at this redshift.
      timeLens            =cosmologyFunctions_ %cosmicTime                 (              &
           &                cosmologyFunctions_%expansionFactorFromRedshift (             &
           &                                                                 redshiftLens &
           &                                                                )             &
           &                                                               )
      ! Find comoving distance to the lens.
      distanceComovingLens=cosmologyFunctions_ %distanceComoving           (              &
           &                                                                 timeLens     &
           &                                                               )
      ! Evaluate the integrand.
      emptyBeamConvergenceIntegrand=                                                          &
           &                        -1.5d0                                                    &
           &                        /speedLight                                               &
           &                        *kilo                                                     &
           &                        *cosmologyParameters_%HubbleConstant        (        )**2 &
           &                        /cosmologyFunctions_ %hubbleParameterEpochal(timeLens)    &
           &                        *cosmologyParameters_%OmegaMatter           (        )    &
           &                        *(                                                        &
           &                          +1.0d0                                                  &
           &                          +redshiftLens                                           &
           &                         )                                                        &
           &                        *  distanceComovingLens                                   &
           &                        *(                                                        &
           &                          +distanceComovingSource                                 &
           &                          -distanceComovingLens                                   &
           &                        )                                                         &
           &                        /  distanceComovingSource
      return
    end function emptyBeamConvergenceIntegrand
    
    function convergenceVarianceIntegrand(redshiftLens,parameterPointer) bind(c)
      !% Integral for variance in the gravitational lensing convergence.
      use, intrinsic :: ISO_C_Binding
      use Numerical_Constants_Physical
      use Numerical_Constants_Prefixes
      use Numerical_Constants_Math
      implicit none
      real            (c_double                  )            :: convergenceVarianceIntegrand
      real            (c_double                  ), value     :: redshiftLens
      type            (c_ptr                     ), value     :: parameterPointer
      double precision                            , parameter :: wavenumberMinimum            =0.0d0
      logical                                                 :: integrationReset
      double precision                                        :: distanceComovingLens               , wavenumberMaximum, &
           &                                                     lensingPower
      type            (fgsl_function             )            :: integrandFunction
      type            (fgsl_integration_workspace)            :: integrationWorkspace
      
      ! Find cosmic time at this redshift.
      timeLens            =cosmologyFunctions_ %cosmicTime                 (              &
           &                cosmologyFunctions_%expansionFactorFromRedshift (             &
           &                                                                 redshiftLens &
           &                                                                )             &
           &                                                               )
      ! Find comoving distance to the lens.
      distanceComovingLens=cosmologyFunctions_ %distanceComoving           (              &
           &                                                                 timeLens     &
           &                                                               )
      ! Integrate over the power spectrum.
      wavenumberMaximum   =1.0d0/scaleSource
      integrationReset    =.true.
      lensingPower        =Integrate(                                           &
           &                         wavenumberMinimum                        , &
           &                         wavenumberMaximum                        , &
           &                         convergenceVariancePowerSpectrumIntegrand, &
           &                         parameterPointer                         , &
           &                         integrandFunction                        , &
           &                         integrationWorkspace                     , &
           &                         toleranceRelative=1.0d-3                 , &
           &                         reset=integrationReset                     &
           &                        )
      call Integrate_Done(integrandFunction,integrationWorkspace)
      ! Evaluate the integrand.
      convergenceVarianceIntegrand =                                                          &
           &                        +9.0d0                                                    &
           &                        /8.0d0                                                    &
           &                        /Pi                                                       &
           &                        /speedLight                                           **3 &
           &                        *kilo                                                 **3 &
           &                        *cosmologyParameters_%HubbleConstant        (        )**4 &
           &                        /cosmologyFunctions_ %hubbleParameterEpochal(timeLens)    &
           &                        *cosmologyParameters_%OmegaMatter           (        )**2 &
           &                        *(                                                        &
           &                          +1.0d0                                                  &
           &                          +redshiftLens                                           &
           &                         )                                                    **2 &
           &                        *(                                                        &
           &                          +  distanceComovingLens                                 &
           &                          *(                                                      &
           &                            +distanceComovingSource                               &
           &                            -distanceComovingLens                                 &
           &                           )                                                      &
           &                          / distanceComovingSource                                &
           &                         )                                                    **2 &
           &                        *lensingPower
      return
    end function convergenceVarianceIntegrand
    
    function convergenceVariancePowerSpectrumIntegrand(wavenumber,parameterPointer) bind(c)
      !% Integral over power spectrum used in computing the variance in the gravitational lensing convergence.
      use, intrinsic :: ISO_C_Binding
      use Power_Spectra_Nonlinear
      implicit none
      real            (c_double)        :: convergenceVariancePowerSpectrumIntegrand
      real            (c_double), value :: wavenumber
      type            (c_ptr   ), value :: parameterPointer

      convergenceVariancePowerSpectrumIntegrand=wavenumber*Power_Spectrum_Nonlinear(wavenumber,timeLens)
      return
    end function convergenceVariancePowerSpectrumIntegrand

    function convergenceDistributionMoment0Integrand(scaledConvergence,parameterPointer) bind(c)
      !% Integral over scaled convergence distribution.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real            (c_double)        :: convergenceDistributionMoment0Integrand
      real            (c_double), value :: scaledConvergence
      type            (c_ptr   ), value :: parameterPointer

      convergenceDistributionMoment0Integrand=self%convergenceDistribution(scaledConvergence)
      return
    end function convergenceDistributionMoment0Integrand

    function convergenceDistributionMoment1Integrand(scaledConvergence,parameterPointer) bind(c)
      !% Integral of first moment over scaled convergence distribution.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real            (c_double)        :: convergenceDistributionMoment1Integrand
      real            (c_double), value :: scaledConvergence
      type            (c_ptr   ), value :: parameterPointer

      convergenceDistributionMoment1Integrand=scaledConvergence*self%convergenceDistribution(scaledConvergence)
      return
    end function convergenceDistributionMoment1Integrand

    function convergenceDistributionMoment2Integrand(scaledConvergence,parameterPointer) bind(c)
      !% Integral of second moment over scaled convergence distribution.
      use, intrinsic :: ISO_C_Binding
      implicit none
      real            (c_double)        :: convergenceDistributionMoment2Integrand
      real            (c_double), value :: scaledConvergence
      type            (c_ptr   ), value :: parameterPointer

      convergenceDistributionMoment2Integrand=scaledConvergence**2*self%convergenceDistribution(scaledConvergence)
      return
    end function convergenceDistributionMoment2Integrand

  end subroutine takahashi2011LensingDistributionConstruct
