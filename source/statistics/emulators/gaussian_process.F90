!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a Gaussian process emulator read from an emulator file.
  !!}

  type :: emulatorGaussianProcessComponent
     !!{RST
     A Gaussian process emulating one principal component coefficient.
     !!}
     double precision                              :: amplitude
     double precision, allocatable, dimension(:  ) :: lengthScales    , alpha
     ! The transpose of the lower-triangular Cholesky factor, L, of the covariance matrix of the training data (so that
     ! factorTransposed(:,i) holds row i of L, and the triangular solve accesses memory contiguously).
     double precision, allocatable, dimension(:,:) :: factorTransposed
  end type emulatorGaussianProcessComponent

  !![
  <emulator name="emulatorGaussianProcess" docformat="rst">
   <description>
   An emulator read from an emulator file (see :ref:`manual-sec-EmulatorFileFormat`), whose observable is compressed by
   principal components analysis, with a Gaussian process (with a Matérn :math:`\nu=5/2` covariance function with a separate
   length scale for each input) emulating each standardized principal component coefficient. The emulator is that for the
   observable ``[label]`` in the file ``[fileName]``. Predictions follow exactly the definition in the emulator file format:
   the mean and variance of each Gaussian process are transformed back to principal component coefficients, and projected
   onto the bins of the observable, neglecting covariance between components.
   </description>
  </emulator>
  !!]
  type, extends(emulatorClass) :: emulatorGaussianProcess
     !!{RST
     A Gaussian process emulator read from an emulator file.
     !!}
     private
     type            (varying_string                  )                              :: fileName         , label
     type            (varying_string                  ), allocatable, dimension(:  ) :: inputNames_
     double precision                                  , allocatable, dimension(:,:) :: inputs           , pcaComponents   , &
          &                                                                             covarianceTarget_
     double precision                                  , allocatable, dimension(:  ) :: binMean          , binScale        , &
          &                                                                             coefficientMean  , coefficientScale, &
          &                                                                             x_               , yTarget_
     type            (emulatorGaussianProcessComponent), allocatable, dimension(:  ) :: components
     double precision                                                                :: jitter
   contains
     procedure :: countInputs  => gaussianProcessCountInputs
     procedure :: inputNames   => gaussianProcessInputNames
     procedure :: countOutputs => gaussianProcessCountOutputs
     procedure :: outputs      => gaussianProcessOutputs
     procedure :: target       => gaussianProcessTarget
     procedure :: predict      => gaussianProcessPredict
  end type emulatorGaussianProcess

  interface emulatorGaussianProcess
     !!{RST
     Constructors for the :galacticus-class:`emulatorGaussianProcess` emulator class.
     !!}
     module procedure gaussianProcessConstructorParameters
     module procedure gaussianProcessConstructorInternal
  end interface emulatorGaussianProcess

contains

  function gaussianProcessConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`emulatorGaussianProcess` emulator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(emulatorGaussianProcess)                :: self
    type(inputParameters        ), intent(inout) :: parameters
    type(varying_string         )                :: fileName  , label

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <description>
      The name of the emulator file.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>label</name>
      <description>
      The label of the emulated observable within the emulator file.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=emulatorGaussianProcess(fileName,label)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gaussianProcessConstructorParameters

  function gaussianProcessConstructorInternal(fileName,label) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`emulatorGaussianProcess` emulator class: read the emulator from file.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File     , hdf5Group
    use :: ISO_Varying_String, only : char         , operator(//), operator(/=)   , var_str
    use :: Linear_Algebra    , only : assignment(=), matrix      , matrixCholesky
    use :: String_Handling   , only : operator(//)
    implicit none
    type            (emulatorGaussianProcess)                              :: self
    type            (varying_string         ), intent(in   )               :: fileName       , label
    double precision                         , allocatable, dimension(:,:) :: covariance     , factor
    double precision                         , allocatable, dimension(:  ) :: logLengthScales, noiseVariance
    type            (varying_string         )                              :: formatName     , kernel
    integer                                                                :: formatVersion  , countComponents, &
         &                                                                    k              , i              , &
         &                                                                    j
    double precision                                                       :: logAmplitude
    logical                                                                :: hasFactor
    !![
    <constructorAssign variables="fileName, label"/>
    !!]

    !$ call hdf5Access%set()
    block
      type(hdf5File ) :: file
      type(hdf5Group) :: emulatorGroup, trainingSetGroup
      file=hdf5File(char(fileName),readOnly=.true.)
      call file%readAttribute('format'       ,formatName   )
      call file%readAttribute('formatVersion',formatVersion)
      if (formatName /= 'galacticusEmulator') call Error_Report("'"//fileName//"' is not a Galacticus emulator file"//{introspection:location})
      if (formatVersion /= 1) call Error_Report(var_str("'")//fileName//"' has format version "//formatVersion//"; only version 1 is supported"//{introspection:location})
      if (.not.file%hasGroup('emulators/'//char(label))) call Error_Report("no emulator for '"//label//"' in '"//fileName//"'"//{introspection:location})
      emulatorGroup   =file%openGroup('emulators/'   //char(label))
      trainingSetGroup=file%openGroup('trainingSets/'//char(label))
      call emulatorGroup%readAttribute('kernel'         ,kernel         )
      call emulatorGroup%readAttribute('jitter'         ,self%jitter    )
      call emulatorGroup%readAttribute('countComponents',countComponents)
      if (kernel /= 'matern52ARD') call Error_Report("unsupported kernel '"//kernel//"'"//{introspection:location})
      call emulatorGroup   %readDataset('inputNames'      ,self%inputNames_      )
      call emulatorGroup   %readDataset('inputs'          ,self%inputs           )
      call emulatorGroup   %readDataset('binMean'         ,self%binMean          )
      call emulatorGroup   %readDataset('binScale'        ,self%binScale         )
      call emulatorGroup   %readDataset('pcaComponents'   ,self%pcaComponents    )
      call emulatorGroup   %readDataset('coefficientMean' ,self%coefficientMean  )
      call emulatorGroup   %readDataset('coefficientScale',self%coefficientScale )
      call trainingSetGroup%readDataset('x'               ,self%x_               )
      call trainingSetGroup%readDataset('yTarget'         ,self%yTarget_         )
      call trainingSetGroup%readDataset('covarianceTarget',self%covarianceTarget_)
      allocate(self%components(countComponents))
      do k=1,countComponents
         block
           type(hdf5Group) :: componentGroup_
           componentGroup_=emulatorGroup%openGroup(char(var_str('component')//k))
           call componentGroup_%readAttribute('logAmplitude'   ,logAmplitude                )
           call componentGroup_%readDataset  ('logLengthScales',logLengthScales             )
           call componentGroup_%readDataset  ('alpha'          ,self%components(k)%alpha    )
           call componentGroup_%readDataset  ('noiseVariance'  ,noiseVariance               )
           hasFactor=componentGroup_%hasDataset('choleskyFactor')
           if (hasFactor) call componentGroup_%readDataset('choleskyFactor',self%components(k)%factorTransposed)
         end block
         self%components(k)%amplitude   =exp(logAmplitude   )
         self%components(k)%lengthScales=exp(logLengthScales)
         if (.not.hasFactor) then
            ! No stored factor - compute it from the covariance matrix of the training data.
            allocate(covariance(size(self%inputs,dim=2),size(self%inputs,dim=2)))
            allocate(factor    (size(self%inputs,dim=2),size(self%inputs,dim=2)))
            do i=1,size(self%inputs,dim=2)
               do j=1,i
                  covariance(i,j)=gaussianProcessKernel(self%inputs(:,i),self%inputs(:,j),self%components(k)%amplitude,self%components(k)%lengthScales)
                  covariance(j,i)=covariance(i,j)
               end do
               covariance(i,i)=covariance(i,i)+noiseVariance(i)+self%jitter
            end do
            block
              type(matrixCholesky) :: decomposition
              decomposition=matrixCholesky(matrix(covariance))
              factor       =decomposition
            end block
            self%components(k)%factorTransposed=transpose(factor)
            deallocate(covariance,factor)
         end if
      end do
    end block
    !$ call hdf5Access%unset()
    ! Validate shapes.
    if (size(self%inputNames_) /= size(self%inputs,dim=1)) call Error_Report('inconsistent number of inputs'//{introspection:location})
    if (size(self%pcaComponents,dim=1) /= size(self%binMean) .or. size(self%pcaComponents,dim=2) /= countComponents) &
         & call Error_Report('inconsistent shape of principal components'//{introspection:location})
    return
  end function gaussianProcessConstructorInternal

  double precision function gaussianProcessKernel(input1,input2,amplitude,lengthScales) result(kernel)
    !!{RST
    The Matérn :math:`\nu=5/2` covariance function with a separate length scale for each input,
    :math:`k = A (1 + \sqrt{5} r + 5 r^2/3) \exp(-\sqrt{5} r)`, with :math:`r^2 = \sum_i [(u_i-u^\prime_i)/\ell_i]^2`.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: input1      , input2, &
         &                                           lengthScales
    double precision, intent(in   )               :: amplitude
    double precision                              :: radius

    radius=sqrt(sum(((input1-input2)/lengthScales)**2))
    kernel=+amplitude                                        &
         & *(1.0d0+sqrt(5.0d0)*radius+5.0d0*radius**2/3.0d0) &
         & *exp(-sqrt(5.0d0)*radius)
    return
  end function gaussianProcessKernel

  integer function gaussianProcessCountInputs(self)
    !!{RST
    Return the number of inputs to the emulator.
    !!}
    implicit none
    class(emulatorGaussianProcess), intent(inout) :: self

    gaussianProcessCountInputs=size(self%inputNames_)
    return
  end function gaussianProcessCountInputs

  subroutine gaussianProcessInputNames(self,names)
    !!{RST
    Return the names of the inputs to the emulator.
    !!}
    implicit none
    class  (emulatorGaussianProcess), intent(inout)                            :: self
    type   (varying_string         ), intent(  out), allocatable, dimension(:) :: names
    integer                                                                    :: i

    ! An array of varying strings is not allocated on assignment, so allocate it, and copy each element.
    allocate(names(size(self%inputNames_)))
    do i=1,size(self%inputNames_)
       names(i)=self%inputNames_(i)
    end do
    return
  end subroutine gaussianProcessInputNames

  integer function gaussianProcessCountOutputs(self)
    !!{RST
    Return the number of outputs of the emulator.
    !!}
    implicit none
    class(emulatorGaussianProcess), intent(inout) :: self

    gaussianProcessCountOutputs=size(self%binMean)
    return
  end function gaussianProcessCountOutputs

  function gaussianProcessOutputs(self) result(outputs)
    !!{RST
    Return the abscissae of the outputs of the emulator.
    !!}
    implicit none
    double precision                         , allocatable  , dimension(:) :: outputs
    class           (emulatorGaussianProcess), intent(inout)               :: self

    outputs=self%x_
    return
  end function gaussianProcessOutputs

  subroutine gaussianProcessTarget(self,yTarget,covarianceTarget)
    !!{RST
    Return the target data of the emulated observable.
    !!}
    implicit none
    class           (emulatorGaussianProcess), intent(inout)                              :: self
    double precision                         , intent(  out), allocatable, dimension(:  ) :: yTarget
    double precision                         , intent(  out), allocatable, dimension(:,:) :: covarianceTarget

    yTarget         =self%yTarget_
    covarianceTarget=self%covarianceTarget_
    return
  end subroutine gaussianProcessTarget

  subroutine gaussianProcessPredict(self,quantiles,mean,variance)
    !!{RST
    Predict the mean and variance of the emulated observable. For each component, :math:`k`, the Gaussian process gives a mean
    :math:`m_k = \mathbf{k}_*^\mathrm{T} \alpha_k` and variance :math:`v_k = A_k - |\mathsf{L}_k^{-1} \mathbf{k}_*|^2`
    (limited to be non-negative), which are transformed back to the coefficients of the principal components and projected
    onto the bins.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (emulatorGaussianProcess), intent(inout)               :: self
    double precision                         , intent(in   ), dimension(:) :: quantiles
    double precision                         , intent(  out), dimension(:) :: mean           , variance
    double precision                         , allocatable  , dimension(:) :: covarianceCross, solved
    double precision                                                       :: meanComponent  , varianceComponent  , &
         &                                                                    coefficient    , varianceCoefficient
    integer                                                                :: k              , i                  , &
         &                                                                    countPoints

    if (size(quantiles) /= size(self%inputNames_)                                          ) &
         & call Error_Report('incorrect number of inputs' //{introspection:location})
    if (size(mean     ) /= size(self%binMean    ) .or. size(variance) /= size(self%binMean)) &
         & call Error_Report('incorrect number of outputs'//{introspection:location})
    countPoints=size(self%inputs,dim=2)
    allocate(covarianceCross(countPoints))
    allocate(solved         (countPoints))
    mean    =self%binMean
    variance=0.0d0
    do k=1,size(self%components)
       ! The covariance between the input point and each training point.
       do i=1,countPoints
          covarianceCross(i)=gaussianProcessKernel(quantiles,self%inputs(:,i),self%components(k)%amplitude,self%components(k)%lengthScales)
       end do
       meanComponent=dot_product(covarianceCross,self%components(k)%alpha)
       ! Solve L v = k* by forward substitution; the variance is A - |v|².
       do i=1,countPoints
          solved(i)=+(                                                                         &
               &      +covarianceCross(i)                                                      &
               &      -dot_product(self%components(k)%factorTransposed(1:i-1,i),solved(1:i-1)) &
               &     )                                                                         &
               &    /self%components(k)%factorTransposed(i,i)
       end do
       varianceComponent  =max(self%components(k)%amplitude-dot_product(solved,solved),0.0d0)
       ! Transform back to the coefficient of this principal component, and project onto the bins.
       coefficient        =self%coefficientMean (k)+self%coefficientScale(k)   *meanComponent
       varianceCoefficient=                         self%coefficientScale(k)**2*varianceComponent
       mean               =mean    +self%binScale   *coefficient        *self%pcaComponents(:,k)
       variance           =variance+self%binScale**2*varianceCoefficient*self%pcaComponents(:,k)**2
    end do
    return
  end subroutine gaussianProcessPredict
